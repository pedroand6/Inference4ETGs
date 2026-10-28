#!/usr/bin/env bash
# Runs the spectral or photometric fits (FITS_SUBDIR), INSIDE an hljob job (started by submit_fits.sh as
# `bash inputs/fit_job.sh`).
#
# The job's working folder is /work. inputs/ holds this file and bundle.tar.gz; whatever is written to out/
# comes back to the laptop, also when the job is killed by its time limit (so results are copied to out/
# every few minutes, not only at the end).
#
# Environment variables (set by submit_fits.sh with hljob -e):
#   WORKERS          parallel fits, one core each        [default 14]
#   MAX_ITERATIONS   MultiNest iteration limit per fit   [default 50000]
#   SYNC_SECONDS     seconds between copies of the results to out/   [default 600]
#   FITS_SUBDIR      folder of the jobs and results: spectral_fits, photometric_fits or photometric_fits_dblplaw
#                    [default spectral_fits]
set -euo pipefail

# 1. Settings
WORKERS="${WORKERS:-14}"
MAX_ITERATIONS="${MAX_ITERATIONS:-50000}"
SYNC_SECONDS="${SYNC_SECONDS:-600}"
JOB_DIR="/work"
export OUT_DIR="$JOB_DIR/out"
RUN_DIR="$JOB_DIR/run"        # short path: MultiNest limits the length of its output paths
FIT_DIR="${FITS_SUBDIR:-spectral_fits}"   # relative to RUN_DIR, as run_spectral_fits.py expects

mkdir -p "$OUT_DIR" "$RUN_DIR"
# fetch_fits.sh reads this, so the results of photometric fits can never land among the spectral ones
echo "$FIT_DIR" > "$OUT_DIR/fits_subdir.txt"

# 2. Unpack the bundle
echo "== 2. Unpacking the bundle"
tar -xzf "$JOB_DIR/inputs/bundle.tar.gz" -C "$RUN_DIR"
cd "$RUN_DIR"
source job_common.sh

# 3. Environment (Python 3.10, packages, MultiNest, grid check); stops the job if the grids differ
echo "== 3. Building the environment"
bash build_env.sh
source env/activate.sh

# 4. Put the fits finished in earlier jobs where the runner looks for them, so it skips them
echo "== 4. Restoring finished fits"
mkdir -p "$FIT_DIR"
if compgen -G "finished_fits/*.h5" > /dev/null; then
    cp finished_fits/*.h5 "$FIT_DIR/"
fi
python run_spectral_fits.py --jobs-dir "$FIT_DIR/jobs" --output-dir "$FIT_DIR" --status

copy_results_to_out() {
    # Copies the results to out/. Whole files are copied every time (they are small) so that a copy caught
    # in the middle of a write is replaced by the next one.
    if compgen -G "$FIT_DIR/*.h5" > /dev/null; then
        cp "$FIT_DIR"/*.h5 "$OUT_DIR/"
    fi
    if [[ -f "$FIT_DIR/fit_log.csv" ]]; then
        cp "$FIT_DIR/fit_log.csv" "$OUT_DIR/"
    fi
    if [[ -d "$FIT_DIR/fit_logs" ]]; then
        cp -r "$FIT_DIR/fit_logs" "$OUT_DIR/"
    fi
    echo "$(read_memory_peak_bytes)" > "$OUT_DIR/memory_peak_bytes.txt"
}

sync_forever() {
    while true; do
        sleep "$SYNC_SECONDS"
        copy_results_to_out
    done
}

# 5. Copy the results periodically in the background, and once more however the script ends
echo "== 5. Starting the fits: $WORKERS workers, up to $MAX_ITERATIONS iterations each"
start_memory_sampler "$OUT_DIR/memory_samples.csv"
sync_forever &
sync_pid=$!
trap 'kill "$sync_pid" "$MEMORY_SAMPLER_PID" 2>/dev/null || true; copy_results_to_out' EXIT

# 6. The fits. A failing runner must not skip the final copy, so its exit code is kept instead of stopping the script
set_stage "fits" "$OUT_DIR/memory_samples.csv"
runner_exit_code=0
python run_spectral_fits.py --jobs-dir "$FIT_DIR/jobs" --output-dir "$FIT_DIR" \
    --workers "$WORKERS" --max-iterations "$MAX_ITERATIONS" || runner_exit_code=$?

echo "== 6. Runner finished with exit code $runner_exit_code"
python run_spectral_fits.py --jobs-dir "$FIT_DIR/jobs" --output-dir "$FIT_DIR" --status
exit "$runner_exit_code"
