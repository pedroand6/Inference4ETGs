#!/usr/bin/env bash
# Trains the SFH posterior estimators on the GPU, INSIDE an hljob job (started by submit_training.sh as
# `bash inputs/train_job.sh`).
#
# The job's working folder is /work. inputs/ holds this file and bundle.tar.gz; whatever is written to out/
# comes back to the laptop, also when the job is killed by its time limit. Each run is written straight to
# out/ml_runs/<embedding>/ (with the model's suffix, e.g. mlp_dblplaw), so a run that finished before the
# limit is never lost.
#
# Environment variables (set by submit_training.sh with hljob -e):
#   EMBEDDINGS   embeddings to train, one after the other   [default "mlp"]
#   SFH_MODEL    SFH model to train: delayed or dblplaw          [default delayed]
#   ALLOW_CPU_TRAINING  set to 1 to train on the CPU when no GPU is usable   [default: fail without a GPU]
set -euo pipefail

# 1. Settings
EMBEDDINGS="${EMBEDDINGS:-mlp}"
SFH_MODEL="${SFH_MODEL:-delayed}"
ALLOW_CPU_TRAINING="${ALLOW_CPU_TRAINING:-0}"
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
INSTALL_TORCH=1 SKIP_MULTINEST=1 bash build_env.sh
source env/activate.sh

# 3b. GPU check. Without it `--device auto` would silently fall back to the CPU and burn the whole time limit
# (up to 8 h), so fail fast unless CPU training was explicitly allowed. Both the driver (nvidia-smi) and the
# torch build in the venv (a CPU-only wheel cannot see the GPU) must work.
if [[ "$ALLOW_CPU_TRAINING" == "1" ]]; then
    echo "ALLOW_CPU_TRAINING=1: training on the CPU, no GPU check"
    DEVICE="cpu"
else
    if ! nvidia-smi; then
        echo "ERROR: nvidia-smi failed: no GPU visible to the job. Set ALLOW_CPU_TRAINING=1 to train on the CPU." >&2
        exit 1
    fi
    if ! python -c "import torch; assert torch.cuda.is_available()"; then
        echo "ERROR: torch cannot use CUDA (CPU-only build or unseen GPU). Set ALLOW_CPU_TRAINING=1 to train on the CPU." >&2
        exit 1
    fi
    DEVICE="cuda"
fi
echo "Training device: $DEVICE"

# 4. Memory samples, written to out/ for the record; the sampler is stopped however the script ends
start_memory_sampler "$OUT_DIR/memory_samples.csv"
trap 'kill "$MEMORY_SAMPLER_PID" 2>/dev/null || true; echo "$(read_memory_peak_bytes)" > "$OUT_DIR/memory_peak_bytes.txt"' EXIT

# 5. One training run per embedding. A failing run must not stop the next one, so exit codes are kept.
failed_embeddings=""
for embedding in $EMBEDDINGS; do
    echo "== 5. Training the '$embedding' embedding of the '$SFH_MODEL' SFH model"
    set_stage "train_${embedding}_${SFH_MODEL}" "$OUT_DIR/memory_samples.csv"
    training_exit_code=0
    python train_ml_model.py --sfh-model "$SFH_MODEL" --embedding "$embedding" \
        --output-dir "$OUT_DIR/ml_runs" --device "$DEVICE" \
        >> "$OUT_DIR/train_${embedding}_${SFH_MODEL}.log" 2>&1 || training_exit_code=$?
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
