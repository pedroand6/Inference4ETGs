"""Run the MultiNest spectral fits of the exported jobs, one single-core fit per process.

Each job is a file <jobs-dir>/<job_id>.npz written by the reference notebook (sm.export_job). Each finished
fit is <output-dir>/<job_id>.h5, so a run can be stopped and started again: finished jobs are skipped.

1. Environment (the PC or cluster needs the same software as the notebook)
   - Python 3.10 with numpy, pandas, h5py, scipy, astropy, extinction, bagpipes and pymultinest.
   - The Bagpipes model grids, installed as for the notebook.
   - libmultinest (the MultiNest C library, built with its Fortran compiler), with its folder on
     LD_LIBRARY_PATH, for example: export LD_LIBRARY_PATH=$HOME/MultiNest/lib:$LD_LIBRARY_PATH

2. Files to copy to the machine, keeping this layout (into one folder, from which you run the commands)
   - synthetic_model.py, run_spectral_fits.py, slurm_fit.sh
   - splus-filters/
   - spectral_fits/jobs/  (all the .npz files)

3. Run on one machine, one fit per physical core (hardware threads of one core are not counted).
   By default the number of workers is the physical cores this process may use (Linux /proc/cpuinfo, within
   the CPU affinity mask) minus one, to leave a core for you; where the topology is not readable it is half
   the usable CPUs, with a printed note. The chosen number is printed at start.
       python run_spectral_fits.py --status
       python run_spectral_fits.py --workers 19

4. Run on a Slurm cluster, one fit per array task (edit the array size and paths in slurm_fit.sh)
       sbatch slurm_fit.sh
   Each array task writes its own log file, <output-dir>/fit_logs/<job_id>.csv, so tasks never write to the
   same file on a shared filesystem.

5. Copy the results back: every <job_id>.h5 in the output folder (spectral_fits/) goes into the
   notebook's spectral_fits/ folder. Then run the notebook's analysis cells again.
   fit_log.csv (local runs) and fit_logs/*.csv (Slurm tasks) list the run time, iterations and convergence
   of every fit made here; copy them back too (the notebook reads both).
"""

# One core per fit: these must be set before numpy is imported, or the math libraries start many threads
import os

SINGLE_THREAD = "1"
for thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[thread_variable] = SINGLE_THREAD

import argparse
import csv
import itertools
import multiprocessing
import socket
import sys
import time

try:
    import pymultinest    # noqa: F401  (only checked here; the fits import it themselves)
except (ImportError, OSError) as import_error:
    sys.exit(f"pymultinest could not be imported ({import_error}).\n"
             "It needs the MultiNest C library: put the folder holding libmultinest.so on LD_LIBRARY_PATH, e.g.\n"
             "  export LD_LIBRARY_PATH=$HOME/MultiNest/lib:$LD_LIBRARY_PATH")

import synthetic_model as sm

DEFAULT_JOBS_DIRECTORY = "spectral_fits/jobs"
DEFAULT_OUTPUT_DIRECTORY = "spectral_fits"
LOG_FILE_NAME = "fit_log.csv"             # local runs: one file, written by the parent process only
TASK_LOG_DIRECTORY = "fit_logs"           # Slurm tasks: one <job_id>.csv each, so no two tasks share a file
CPUINFO_FILE = "/proc/cpuinfo"
LOG_COLUMNS = ["job_id", "status", "converged", "minutes", "n_iterations", "error", "host"]
MAX_BASENAME_LENGTH = 100         # [characters] MultiNest limit on the output basename (the path with its "/")
SECONDS_PER_MINUTE = 60.0
SLURM_TASK_VARIABLE = "SLURM_ARRAY_TASK_ID"


# =============================================================================
# 1. Jobs and their order
# =============================================================================
def list_job_ids(jobs_directory):
    """Job ids found in the jobs folder, in sorted order.

    Inputs: jobs_directory: folder with the <job_id>.npz files (str).
    Returns: sorted list of job ids (str).
    """
    return sorted(name[:-len(".npz")] for name in os.listdir(jobs_directory) if name.endswith(".npz"))


def is_done(job_id, output_directory):
    """Whether the fit of a job has its result file.

    Inputs: job_id: job identifier (str); output_directory: folder of the .h5 files (str).
    Returns: True if <output_directory>/<job_id>.h5 exists.
    """
    return os.path.exists(os.path.join(output_directory, f"{job_id}.h5"))


def interleave_calibration_and_galaxies(job_ids):
    """Alternate calibration mocks and real galaxies, so a partial run has some of both.

    Inputs: job_ids: job identifiers (list of str).
    Returns: the same ids, calibration mocks (prefix sm.CALIBRATION_ID_PREFIX) and galaxies alternating,
             with the longer kind's remainder at the end.
    """
    mocks = [job_id for job_id in job_ids if job_id.startswith(sm.CALIBRATION_ID_PREFIX)]
    galaxies = [job_id for job_id in job_ids if not job_id.startswith(sm.CALIBRATION_ID_PREFIX)]
    order = []
    for pair in itertools.zip_longest(mocks, galaxies):
        order.extend(job_id for job_id in pair if job_id is not None)
    return order


# =============================================================================
# 2. Fitting and logging
# =============================================================================
def fit_one(arguments):
    """Fit one job file and describe the outcome. Runs in a worker process and never raises.

    Inputs: arguments: tuple (job_path, output_directory, max_iterations, n_live_points).
    Returns: dict with the fit_and_save summary (job_id, status, converged, run_time, n_iterations, error)
             plus minutes (wall time of this call) and host (machine name).
    """
    job_path, output_directory, max_iterations, n_live_points = arguments
    job_id = os.path.basename(job_path)[:-len(".npz")]
    start = time.time()
    try:
        summary = sm.fit_and_save(sm.load_job(job_path), output_directory, max_iterations, n_live_points)
    except Exception as error:    # an unreadable job file
        summary = {"job_id": job_id, "status": "failed", "converged": None,
                   "run_time": None, "n_iterations": None, "error": repr(error)}
    summary["minutes"] = (time.time() - start) / SECONDS_PER_MINUTE
    summary["host"] = socket.gethostname()
    return summary


def append_to_log(log_path, summary):
    """Append one finished job to a log file, writing the header first if the file is new.

    Inputs: log_path: the log file (str); summary: dict from fit_one.
    Returns: None.
    """
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
    is_new = not os.path.exists(log_path)
    with open(log_path, "a", newline="") as log_file:
        writer = csv.DictWriter(log_file, fieldnames=LOG_COLUMNS, extrasaction="ignore")
        if is_new:
            writer.writeheader()
        writer.writerow(summary)


def run_jobs(job_paths, output_directory, workers, max_iterations, n_live_points, log_path):
    """Fit the jobs with a pool of single-core worker processes, logging each as it finishes.

    A fresh process per fit (maxtasksperchild=1) so that no memory or sampler state is carried over.
    Inputs: job_paths: job files (list of str); output_directory: folder for results and log (str);
            workers: number of processes (count); max_iterations: MultiNest iteration limit (count);
            n_live_points: live points (count); log_path: file that receives one line per finished job (str).
    Returns: None.
    """
    tasks = [(path, output_directory, max_iterations, n_live_points) for path in job_paths]
    with multiprocessing.Pool(workers, maxtasksperchild=1) as pool:
        for finished_count, summary in enumerate(pool.imap_unordered(fit_one, tasks), start=1):
            append_to_log(log_path, summary)
            print(f"[{finished_count}/{len(tasks)}] {summary['job_id']}: {summary['status']}, "
                  f"converged={summary['converged']}, {summary['minutes']:.1f} min, "
                  f"iterations={summary['n_iterations']}" + (f", error={summary['error']}" if summary["error"] else ""),
                  flush=True)


# =============================================================================
# 3. Command line
# =============================================================================
def usable_cpu_ids():
    """Ids of the logical CPUs this process may use (respects Slurm and taskset where Linux reports it).

    Returns: set of processor numbers (int); every CPU if the affinity mask is not available.
    """
    if hasattr(os, "sched_getaffinity"):
        return set(os.sched_getaffinity(0))
    return set(range(os.cpu_count() or 1))


def count_physical_cores(cpuinfo_text, allowed_cpus):
    """Count the physical cores behind some logical CPUs, so hardware threads of one core count once.

    A core is a unique (physical id, core id) pair in the Linux /proc/cpuinfo text.
    Inputs: cpuinfo_text: contents of /proc/cpuinfo (str); allowed_cpus: processor numbers to count (set of int).
    Returns: number of physical cores (count), or None if the text has no usable topology fields.
    """
    cores = set()
    for block in cpuinfo_text.split("\n\n"):
        fields = dict(line.split(":", 1) for line in block.splitlines() if ":" in line)
        fields = {name.strip(): value.strip() for name, value in fields.items()}
        if {"processor", "physical id", "core id"} <= fields.keys() and int(fields["processor"]) in allowed_cpus:
            cores.add((fields["physical id"], fields["core id"]))
    return len(cores) or None


def default_workers():
    """Default number of parallel fits: physical cores usable by this process minus one, at least 1.

    Falls back to half the usable logical CPUs (assuming two threads per core) when /proc/cpuinfo cannot
    be read or has no topology.
    Returns: tuple (workers (count), note describing how it was chosen (str)).
    """
    allowed_cpus = usable_cpu_ids()
    try:
        with open(CPUINFO_FILE) as cpuinfo_file:
            physical_cores = count_physical_cores(cpuinfo_file.read(), allowed_cpus)
    except OSError:
        physical_cores = None
    if physical_cores is not None:
        return max(1, physical_cores - 1), f"{physical_cores} physical cores usable, minus one"
    return (max(1, len(allowed_cpus) // 2 - 1),
            f"NOTE: core topology unreadable; assuming 2 threads per core on {len(allowed_cpus)} usable CPUs, "
            "minus one (set --workers to override)")


def parse_arguments():
    """Read the command-line options (see the module docstring for how to use them).

    Returns: argparse.Namespace.
    """
    parser = argparse.ArgumentParser(description="Run the MultiNest spectral fits of the exported jobs.")
    parser.add_argument("--jobs-dir", default=DEFAULT_JOBS_DIRECTORY, help="folder with the job .npz files")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIRECTORY,
                        help="folder for the .h5 results and fit_log.csv (keep the path short: MultiNest limit)")
    parser.add_argument("--workers", type=int, default=None,
                        help="parallel fits, one physical core each (default: physical cores usable by this "
                             "process minus one; half the usable CPUs minus one if the topology is unreadable)")
    parser.add_argument("--max-iterations", type=int, default=sm.MAX_ITERATIONS, help="MultiNest iteration limit")
    parser.add_argument("--n-live", type=int, default=sm.N_LIVE_POINTS, help="MultiNest live points")
    parser.add_argument("--task-index", type=int,
                        default=int(os.environ[SLURM_TASK_VARIABLE]) if SLURM_TASK_VARIABLE in os.environ else None,
                        help=f"run only this job, by position in the sorted job list (default: ${SLURM_TASK_VARIABLE})")
    parser.add_argument("--status", action="store_true", help="print how many jobs are done and pending, then exit")
    return parser.parse_args()


def main():
    """Fit the pending jobs (or one array task) as the command line asks, or print the job status.

    Returns: None; exits with a message on a missing job folder, a too-long output path or a bad task index.
    """
    arguments = parse_arguments()
    job_ids = list_job_ids(arguments.jobs_dir)
    if not job_ids:
        sys.exit(f"No .npz job files in {arguments.jobs_dir}")

    # 1. Status: counts only
    pending_ids = [job_id for job_id in job_ids if not is_done(job_id, arguments.output_dir)]
    if arguments.status:
        print(f"{len(job_ids)} jobs: {len(job_ids) - len(pending_ids)} done, {len(pending_ids)} pending")
        return

    # 2. MultiNest cannot take output basenames of 100 characters or more
    longest_basename = max(len(os.path.join(arguments.output_dir, job_id) + "/") for job_id in job_ids)
    if longest_basename >= MAX_BASENAME_LENGTH:
        sys.exit(f"The MultiNest output path would be {longest_basename} characters (limit "
                 f"{MAX_BASENAME_LENGTH - 1}). Use a shorter relative --output-dir, e.g. 'fits'.")
    os.makedirs(arguments.output_dir, exist_ok=True)

    # 3. Which jobs to run: one array task, or every pending job
    if arguments.task_index is not None:
        if not 0 <= arguments.task_index < len(job_ids):
            sys.exit(f"--task-index {arguments.task_index} is outside 0..{len(job_ids) - 1} ({len(job_ids)} jobs)")
        selected_ids = [job_ids[arguments.task_index]]
        workers = 1
        workers_note = "Slurm array task"
        log_path = os.path.join(arguments.output_dir, TASK_LOG_DIRECTORY, f"{selected_ids[0]}.csv")
    else:
        selected_ids = interleave_calibration_and_galaxies(pending_ids)
        if arguments.workers is None:
            requested_workers, workers_note = default_workers()
        else:
            requested_workers, workers_note = arguments.workers, "set with --workers"
        workers = min(requested_workers, max(1, len(selected_ids)))
        log_path = os.path.join(arguments.output_dir, LOG_FILE_NAME)
    print(f"Running {len(selected_ids)} job(s) on {workers} worker(s) ({workers_note}), "
          f"max {arguments.max_iterations} iterations, {arguments.n_live} live points", flush=True)

    job_paths = [os.path.join(arguments.jobs_dir, f"{job_id}.npz") for job_id in selected_ids]
    run_jobs(job_paths, arguments.output_dir, workers, arguments.max_iterations, arguments.n_live, log_path)


if __name__ == "__main__":
    main()
