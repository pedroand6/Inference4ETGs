#!/usr/bin/env bash
# Runs LOCALLY: packs everything the fits need into one bundle and submits it to the job server with hljob.
# The job is detached, so this prints the job id and returns at once; follow it with the commands printed at the end.
#
# Fits already finished in Synthetic_domain/spectral_fits/*.h5 travel in the bundle, so a second job
# continues where the previous one stopped instead of repeating fits.
#
# Settings are environment variables (defaults for the full run); extra arguments go to `hljob run` as they are
# (for example  -e SYNC_SECONDS=300).
#   CORES            vCPUs of the job (the server's maximum is 16)           [default 16]
#   MEMORY_MB        memory limit in MB; --pause-llm frees the AI models, so up to about 7000   [default 7000]
#   TIME_LIMIT       wall-clock limit HH:MM:SS (the server's maximum is 24 h)   [default 24:00:00]
#   WORKERS          parallel fits, one core each                              [default 14]
#   MAX_ITERATIONS   MultiNest iteration limit per fit                         [default 50000]
#   JOB_SCRIPT       script run inside the job: fit_job.sh or benchmark_job.sh   [default fit_job.sh]
#   JOB_NAME         name shown in `hljob ls`                                  [default spectral_fits]
#   MULTINEST_SOURCE the MultiNest source folder                               [default ~/Programs/MultiNest]
set -euo pipefail

# 1. Settings and paths
CORES="${CORES:-16}"
MEMORY_MB="${MEMORY_MB:-7000}"
TIME_LIMIT="${TIME_LIMIT:-24:00:00}"
WORKERS="${WORKERS:-14}"
MAX_ITERATIONS="${MAX_ITERATIONS:-50000}"
JOB_SCRIPT="${JOB_SCRIPT:-fit_job.sh}"
JOB_NAME="${JOB_NAME:-spectral_fits}"
MULTINEST_SOURCE="${MULTINEST_SOURCE:-$HOME/Programs/MultiNest}"

REMOTE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYNTHETIC_DIR="$(dirname "$REMOTE_DIR")"
export PATH="$HOME/.local/bin:$PATH"       # where hljob is installed

# 2. Collect the bundle contents in a scratch folder, so the tarball has the layout the job expects
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
bundle_dir="$scratch_dir/bundle_contents"
mkdir -p "$bundle_dir/spectral_fits" "$bundle_dir/finished_fits"

cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/run_spectral_fits.py" "$bundle_dir/"
cp -rL "$SYNTHETIC_DIR/splus-filters" "$bundle_dir/"
cp -rL "$SYNTHETIC_DIR/spectral_fits/jobs" "$bundle_dir/spectral_fits/jobs"
# Jobs listed in skipped_fit_jobs.txt are deliberately left unfitted, so they are not sent
skipped_count=0
if [[ -f "$SYNTHETIC_DIR/skipped_fit_jobs.txt" ]]; then
    while read -r skipped_job_id; do
        if [[ -n "$skipped_job_id" && -f "$bundle_dir/spectral_fits/jobs/$skipped_job_id.npz" ]]; then
            rm "$bundle_dir/spectral_fits/jobs/$skipped_job_id.npz"
            skipped_count=$((skipped_count + 1))
        fi
    done < <(sed 's/#.*//; s/[[:space:]]//g' "$SYNTHETIC_DIR/skipped_fit_jobs.txt" | grep -v '^$' || true)
fi
cp "$REMOTE_DIR/grid_md5.txt" "$REMOTE_DIR/build_env.sh" "$REMOTE_DIR/job_common.sh" \
   "$REMOTE_DIR/fit_job.sh" "$REMOTE_DIR/benchmark_job.sh" "$REMOTE_DIR/summarise_benchmark.py" "$bundle_dir/"

# 3. MultiNest source without the local build products (they are rebuilt in the job)
mkdir -p "$bundle_dir/multinest"
tar -C "$MULTINEST_SOURCE" --exclude=./build --exclude=./lib --exclude=./bin -cf - . | tar -C "$bundle_dir/multinest" -xf -

# 4. Fits finished so far (top level of spectral_fits/ only; pilot/ and superseded/ are not needed)
finished_count=0
for finished_fit in "$SYNTHETIC_DIR"/spectral_fits/*.h5; do
    if [[ -f "$finished_fit" ]]; then
        cp "$finished_fit" "$bundle_dir/finished_fits/"
        finished_count=$((finished_count + 1))
    fi
done

# 5. The tarball
tar -C "$bundle_dir" -czf "$scratch_dir/bundle.tar.gz" .
echo "Bundle: $(du -h "$scratch_dir/bundle.tar.gz" | cut -f1), with $finished_count finished fit(s), $skipped_count skipped job(s) left out"

# 6. Submit. Both files are uploaded to inputs/; the job script unpacks the bundle itself.
job_id="$(hljob run --detach --pause-llm -n "$JOB_NAME" -c "$CORES" -m "$MEMORY_MB" -t "$TIME_LIMIT" \
    -e "WORKERS=$WORKERS" -e "MAX_ITERATIONS=$MAX_ITERATIONS" "$@" \
    -f "$scratch_dir/bundle.tar.gz" -f "$REMOTE_DIR/$JOB_SCRIPT" -- bash "inputs/$JOB_SCRIPT")"

echo "Submitted job $job_id ($CORES cores, $MEMORY_MB MB, $TIME_LIMIT, $WORKERS workers)"
echo "  follow:  hljob logs $job_id -f        (or --tail 30)"
echo "  wait:    hljob wait $job_id"
echo "  fetch:   $REMOTE_DIR/fetch_fits.sh $job_id"
echo "  resume:  run this script again after fetching (finished fits travel in the bundle)"
