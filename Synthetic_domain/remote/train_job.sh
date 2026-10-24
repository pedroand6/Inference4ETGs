#!/usr/bin/env bash
# Trains the SFH posterior estimators on the GPU, INSIDE an hljob job (started by submit_training.sh as
# `bash inputs/train_job.sh`).
#
# The job's working folder is /work. inputs/ holds this file and bundle.tar.gz; whatever is written to out/
# comes back to the laptop, also when the job is killed by its time limit. Each run is written straight to
# out/ml_runs/<embedding>/, so a run that finished before the limit is never lost.
#
# Environment variables (set by submit_training.sh with hljob -e):
#   EMBEDDINGS   embeddings to train, one after the other   [default "mlp cnn"]
set -euo pipefail

# 1. Settings
EMBEDDINGS="${EMBEDDINGS:-mlp cnn}"
JOB_DIR="/work"
export OUT_DIR="$JOB_DIR/out"
RUN_DIR="$JOB_DIR/run"

mkdir -p "$OUT_DIR/ml_runs" "$RUN_DIR"

# 2. Unpack the bundle
echo "== 2. Unpacking the bundle"
tar -xzf "$JOB_DIR/inputs/bundle.tar.gz" -C "$RUN_DIR"
cd "$RUN_DIR"
source job_common.sh

# 3. Environment: the usual build plus PyTorch (GPU build) and zuko
echo "== 3. Building the environment"
INSTALL_TORCH=1 bash build_env.sh
source env/activate.sh
nvidia-smi || echo "WARNING: nvidia-smi not found; the training will run on the CPU"

# 4. Memory samples, written to out/ for the record; the sampler is stopped however the script ends
start_memory_sampler "$OUT_DIR/memory_samples.csv"
trap 'kill "$MEMORY_SAMPLER_PID" 2>/dev/null || true; echo "$(read_memory_peak_bytes)" > "$OUT_DIR/memory_peak_bytes.txt"' EXIT

# 5. One training run per embedding. A failing run must not stop the next one, so exit codes are kept.
failed_embeddings=""
for embedding in $EMBEDDINGS; do
    echo "== 5. Training the '$embedding' embedding"
    set_stage "train_$embedding" "$OUT_DIR/memory_samples.csv"
    training_exit_code=0
    python train_ml_model.py --embedding "$embedding" --output-dir "$OUT_DIR/ml_runs" \
        >> "$OUT_DIR/train_$embedding.log" 2>&1 || training_exit_code=$?
    echo "'$embedding' finished with exit code $training_exit_code"
    if [[ "$training_exit_code" -ne 0 ]]; then
        failed_embeddings="$failed_embeddings $embedding"
    fi
done

# 6. The runs are already in out/ml_runs/ (written there directly); report and set the exit code
echo "== 6. Runs in out/ml_runs:"
ls "$OUT_DIR/ml_runs"
if [[ -n "$failed_embeddings" ]]; then
    echo "ERROR: training failed for:$failed_embeddings" >&2
    exit 1
fi
