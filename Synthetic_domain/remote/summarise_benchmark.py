"""Summarise one benchmark trial: seconds per iteration of each fit, and the memory used by the job.

Used by benchmark_job.sh as:  python summarise_benchmark.py <fits-folder> <memory-samples.csv> <stage-name>
Run from the folder that holds synthetic_model.py (the bundle folder).
"""
import csv
import glob
import os
import sys

import synthetic_model as sm

BYTES_PER_MB = 1e6    # [bytes/MB]


def read_fit_times(fits_directory):
    """Run time and iterations of every fit in a folder.

    Inputs: fits_directory: folder of <job_id>.h5 files (str).
    Returns: list of (job_id, run_time [s], n_iterations) tuples.
    """
    fits = []
    for path in sorted(glob.glob(os.path.join(fits_directory, "*.h5"))):
        _, metadata = sm.load_fit(path)
        job_id = os.path.basename(path)[:-len(".h5")]
        fits.append((job_id, metadata["run_time"], metadata["n_iterations"]))
    return fits


def peak_memory_of_stage(samples_path, stage_name):
    """Highest sampled memory of the job during one stage.

    Inputs: samples_path: CSV written by start_memory_sampler (str); stage_name: stage to look at (str).
    Returns: tuple (highest memory.current [bytes], highest anonymous memory [bytes]); zeros if no samples.
    """
    peak_current = 0
    peak_anon = 0
    with open(samples_path) as samples_file:
        for row in csv.DictReader(samples_file):
            if row["stage"] == stage_name:
                peak_current = max(peak_current, int(row["memory_current_bytes"]))
                peak_anon = max(peak_anon, int(row["memory_anon_bytes"]))
    return peak_current, peak_anon


def main():
    """Print the trial summary."""
    fits_directory, samples_path, stage_name = sys.argv[1:4]
    fits = read_fit_times(fits_directory)
    print(f"{stage_name}: {len(fits)} fits")
    seconds_per_iteration = []
    for job_id, run_time, n_iterations in fits:
        seconds_per_iteration.append(run_time / n_iterations)
        print(f"  {job_id}: {run_time:.0f} s for {n_iterations} iterations "
              f"= {run_time / n_iterations:.3f} s/iteration")
    if seconds_per_iteration:
        mean_value = sum(seconds_per_iteration) / len(seconds_per_iteration)
        print(f"  mean {mean_value:.3f} s/iteration, slowest {max(seconds_per_iteration):.3f} s/iteration")
    peak_current, peak_anon = peak_memory_of_stage(samples_path, stage_name)
    print(f"  memory during the trial: peak memory.current {peak_current / BYTES_PER_MB:.0f} MB, "
          f"peak anonymous {peak_anon / BYTES_PER_MB:.0f} MB")


if __name__ == "__main__":
    main()
