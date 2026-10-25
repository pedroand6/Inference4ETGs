#!/usr/bin/env bash
# Runs LOCALLY: packs everything the training needs into one bundle and submits it to the job server (GPU) with
# hljob. The job is detached, so this prints the job id and returns at once; follow it with the commands printed
# at the end. The trained runs come back in out/ml_runs/<embedding>/.
#
# Settings are environment variables (defaults for the full run); extra arguments go to `hljob run` as they are.
#   CORES            vCPUs of the job                                          [default 8]
#   MEMORY_MB        memory limit in MB                                        [default 7000]
#   TIME_LIMIT       wall-clock limit HH:MM:SS                                 [default 08:00:00]
#   EMBEDDINGS       embeddings to train, one after the other                  [default "mlp cnn"]
#   ALLOW_CPU_TRAINING set to 1 to let the job train on the CPU if no GPU is usable    [default: the job fails]
#   JOB_NAME         name shown in `hljob ls`                                  [default ml_training]
#   MULTINEST_SOURCE the MultiNest source folder (build_env.sh builds it)      [default ~/Programs/MultiNest]
#   BUNDLE_ONLY      set to a folder path to only write bundle.tar.gz there and NOT submit   [default: submit]
#   LOCAL_PYTHON     Python with pandas, used to shrink the catalogue          [default ../../.venv/bin/python]
set -euo pipefail

# 1. Settings and paths
CORES="${CORES:-8}"
MEMORY_MB="${MEMORY_MB:-7000}"
TIME_LIMIT="${TIME_LIMIT:-08:00:00}"
EMBEDDINGS="${EMBEDDINGS:-mlp cnn}"
ALLOW_CPU_TRAINING="${ALLOW_CPU_TRAINING:-0}"
JOB_NAME="${JOB_NAME:-ml_training}"
MULTINEST_SOURCE="${MULTINEST_SOURCE:-$HOME/Programs/MultiNest}"
BUNDLE_ONLY="${BUNDLE_ONLY:-}"

REMOTE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYNTHETIC_DIR="$(dirname "$REMOTE_DIR")"
LOCAL_PYTHON="${LOCAL_PYTHON:-$(dirname "$SYNTHETIC_DIR")/.venv/bin/python}"
export PATH="$HOME/.local/bin:$PATH"       # where hljob is installed

# 2. Collect the bundle contents in a scratch folder, so the tarball has the layout the job expects
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
bundle_dir="$scratch_dir/bundle_contents"
mkdir -p "$bundle_dir"

cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/ml_data.py" "$SYNTHETIC_DIR/ml_networks.py" \
   "$SYNTHETIC_DIR/train_ml_model.py" "$SYNTHETIC_DIR/training_set.parquet" \
   "$SYNTHETIC_DIR/synthetic_validation.parquet" "$SYNTHETIC_DIR/splus_noise_model.csv" "$bundle_dir/"
cp -rL "$SYNTHETIC_DIR/splus-filters" "$bundle_dir/"
cp "$REMOTE_DIR/grid_md5.txt" "$REMOTE_DIR/build_env.sh" "$REMOTE_DIR/job_common.sh" "$bundle_dir/"

# 3. The DR3 catalogue is large; the training only draws magnitudes and redshifts from it, so only those columns
# (and the star/galaxy probability used to select galaxies) go in, written under the same name
"$LOCAL_PYTHON" - "$SYNTHETIC_DIR/dr3_corrected.csv" "$bundle_dir/dr3_corrected.csv" <<'SMALL_CATALOGUE'
import sys
import pandas as pd
pd.read_csv(sys.argv[1], usecols=["PROB_GAL", "r_auto", "zml"]).to_csv(sys.argv[2], index=False)
SMALL_CATALOGUE

# 4. MultiNest source without the local build products (they are rebuilt in the job)
mkdir -p "$bundle_dir/multinest"
tar -C "$MULTINEST_SOURCE" --exclude=./build --exclude=./lib --exclude=./bin -cf - . | tar -C "$bundle_dir/multinest" -xf -

# 5. The tarball
tar -C "$bundle_dir" -czf "$scratch_dir/bundle.tar.gz" .
echo "Bundle: $(du -h "$scratch_dir/bundle.tar.gz" | cut -f1)"
if [[ -n "$BUNDLE_ONLY" ]]; then
    mkdir -p "$BUNDLE_ONLY"
    cp "$scratch_dir/bundle.tar.gz" "$BUNDLE_ONLY/"
    echo "BUNDLE_ONLY set: bundle written to $BUNDLE_ONLY/bundle.tar.gz, nothing submitted"
    exit 0
fi

# 6. Submit (one GPU job that trains every embedding). Both files are uploaded to inputs/.
job_id="$(hljob run --detach --gpu -n "$JOB_NAME" -c "$CORES" -m "$MEMORY_MB" -t "$TIME_LIMIT" \
    -e EMBEDDINGS="$EMBEDDINGS" -e ALLOW_CPU_TRAINING="$ALLOW_CPU_TRAINING" "$@" \
    -f "$scratch_dir/bundle.tar.gz" -f "$REMOTE_DIR/train_job.sh" -- bash inputs/train_job.sh)"

echo "Submitted job $job_id (GPU, $CORES cores, $MEMORY_MB MB, $TIME_LIMIT, embeddings: $EMBEDDINGS)"
echo "  follow:  hljob logs $job_id -f        (or --tail 30)"
echo "  wait:    hljob wait $job_id"
echo "  fetch:   hljob fetch $job_id -o /tmp/ml_job   then copy /tmp/ml_job/out/ml_runs/* to Synthetic_domain/ml_runs/"
