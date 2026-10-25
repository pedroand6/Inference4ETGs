#!/usr/bin/env bash
# Runs LOCALLY: packs a notebook and the files it reads, and submits a job that executes it on the job server,
# so heavy notebooks never run on the laptop. The job is detached; this prints the job id and returns at once.
# Bring the results back with fetch_notebook.sh <job-id>.
#
# Usage:  submit_notebook.sh <notebook> <input> [<input> ...]
#   <notebook>  notebook to execute, relative to Synthetic_domain (e.g. ml_model.ipynb)
#   <input>     files or folders it reads, relative to Synthetic_domain (the Python modules, splus-filters/ and
#               the job scripts are always included)
#
# Settings are environment variables:
#   GPU              set to 1 to run on the GPU (and fail if it is not usable)      [default 0]
#   NEEDS_MULTINEST  set to 1 to build MultiNest in the job (CPU jobs only)         [default 0]
#   CORES            vCPUs of the job                                                [default 8]
#   MEMORY_MB        memory limit in MB                                              [default 16000]
#   TIME_LIMIT       wall-clock limit HH:MM:SS                                       [default 04:00:00]
#   JOB_NAME         name shown in `hljob ls`                                        [default: the notebook name]
#   MULTINEST_SOURCE the MultiNest source folder                                     [default ~/Programs/MultiNest]
#   FLOAT32_FILES    parquet inputs sent with float32 instead of float64 columns, to halve the upload (the job
#                    server buffers the whole upload in a small memory limit)                       [default none]
#   HLJOB_OPTIONS    extra `hljob run` options, e.g. "--pause-llm" to free the memory the LLM holds   [default none]
set -euo pipefail

# 1. Arguments, settings and paths
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <notebook> <input> [<input> ...]" >&2
    exit 1
fi
notebook="$1"
shift
GPU="${GPU:-0}"
NEEDS_MULTINEST="${NEEDS_MULTINEST:-0}"
CORES="${CORES:-8}"
MEMORY_MB="${MEMORY_MB:-16000}"
TIME_LIMIT="${TIME_LIMIT:-04:00:00}"
JOB_NAME="${JOB_NAME:-$(basename "$notebook" .ipynb)}"
MULTINEST_SOURCE="${MULTINEST_SOURCE:-$HOME/Programs/MultiNest}"
read -r -a float32_files <<< "${FLOAT32_FILES:-}"
read -r -a hljob_options <<< "${HLJOB_OPTIONS:-}"

REMOTE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYNTHETIC_DIR="$(dirname "$REMOTE_DIR")"
export PATH="$HOME/.local/bin:$PATH"       # where hljob is installed

# 2. Bundle contents: the notebook, its inputs, the modules and the job helpers, with their Synthetic_domain paths.
# Credentials (.env) are never packed.
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
bundle_dir="$scratch_dir/bundle_contents"
mkdir -p "$bundle_dir"
(
    cd "$SYNTHETIC_DIR"
    tar --exclude="__pycache__" --exclude=".env" -cf - "$notebook" ./*.py splus-filters "$@"
) | tar -C "$bundle_dir" -xf -
cp "$REMOTE_DIR/grid_md5.txt" "$REMOTE_DIR/build_env.sh" "$REMOTE_DIR/job_common.sh" "$bundle_dir/"

# 2b. Large parquet inputs as float32 (about 1e-7 relative precision, far below any noise in the data)
for float32_file in "${float32_files[@]}"; do
    "$SYNTHETIC_DIR/../.venv/bin/python" - "$bundle_dir/$float32_file" <<'TO_FLOAT32'
import sys
import pandas as pd
table = pd.read_parquet(sys.argv[1])
float64_columns = table.select_dtypes("float64").columns
table[float64_columns] = table[float64_columns].astype("float32")
table.to_parquet(sys.argv[1], index=False)
TO_FLOAT32
    echo "Sent as float32: $float32_file ($(du -h "$bundle_dir/$float32_file" | cut -f1))"
done

# 3. MultiNest source without the local build products, only when the notebook runs fits
if [[ "$NEEDS_MULTINEST" == "1" ]]; then
    mkdir -p "$bundle_dir/multinest"
    tar -C "$MULTINEST_SOURCE" --exclude=./build --exclude=./lib --exclude=./bin -cf - . | tar -C "$bundle_dir/multinest" -xf -
fi

# 4. The tarball
tar -C "$bundle_dir" -czf "$scratch_dir/bundle.tar.gz" .
echo "Bundle: $(du -h "$scratch_dir/bundle.tar.gz" | cut -f1)"

# 5. Submit
gpu_option=()
if [[ "$GPU" == "1" ]]; then
    gpu_option=(--gpu)
fi
job_id="$(hljob run --detach "${gpu_option[@]}" "${hljob_options[@]}" -n "$JOB_NAME" -c "$CORES" -m "$MEMORY_MB" -t "$TIME_LIMIT" \
    -e NOTEBOOK="$notebook" -e NEEDS_GPU="$GPU" -e NEEDS_MULTINEST="$NEEDS_MULTINEST" \
    -f "$scratch_dir/bundle.tar.gz" -f "$REMOTE_DIR/notebook_job.sh" -- bash inputs/notebook_job.sh)"

echo "Submitted job $job_id ($notebook, GPU=$GPU, $CORES cores, $MEMORY_MB MB, $TIME_LIMIT)"
echo "  follow:  hljob logs $job_id --tail 30"
echo "  fetch:   remote/fetch_notebook.sh $job_id"
