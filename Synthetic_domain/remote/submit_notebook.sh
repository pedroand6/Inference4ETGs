#!/usr/bin/env bash
# Runs LOCALLY: packs a notebook and the files it reads, and submits a job that executes it on the job server,
# so heavy notebooks never run on the laptop. The job is detached; this prints the job id and returns at once.
# Bring the results back with fetch_notebook.sh <job-id>.
#
# Usage:  submit_notebook.sh <notebook> <input> [<input> ...]
#   <notebook>  notebook to execute, relative to Synthetic_domain (e.g. ml_model.ipynb)
#   <input>     files or folders it reads, relative to Synthetic_domain (the Python modules, splus-filters/ and
#               the job scripts are always included). An input starting with ../ lies elsewhere in the repository
#               (e.g. ../Data/sdss_eboss_firefly-dr16.fits, ../dustMaps/csfd): it is packed at its repository path.
#
# The bundle keeps the repository layout (Synthetic_domain/..., Data/..., dustMaps/...) and the notebook runs from
# Synthetic_domain/, so its ../ paths work in the job as they do here. Credentials (.env) are never packed.
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
#   BUNDLE_ONLY      set to a folder path to only write bundle.tar.gz there and NOT submit           [default: submit]
#   EXCLUDE_NAMES    files or folders never packed from the Synthetic_domain inputs, matched by name at any
#                    depth: results of an earlier model (superseded/) and the pilot fits (pilot/), which no
#                    notebook reads and which swell the upload; empty packs everything [default "superseded pilot"]
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
read -r -a exclude_names <<< "${EXCLUDE_NAMES-superseded pilot}"

BUNDLE_ONLY="${BUNDLE_ONLY:-}"

REMOTE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYNTHETIC_DIR="$(dirname "$REMOTE_DIR")"
REPO_DIR="$(dirname "$SYNTHETIC_DIR")"
export PATH="$HOME/.local/bin:$PATH"       # where hljob is installed

# 2. Inputs of Synthetic_domain and inputs elsewhere in the repository (given with a leading ../)
synthetic_inputs=()
repository_inputs=()
for input in "$@"; do
    if [[ "$input" == ../* ]]; then
        repository_input="${input#../}"
        if [[ "$repository_input" == ../* || "$repository_input" == /* ]]; then
            echo "Input $input lies outside the repository; only ../<path inside the repository> is allowed" >&2
            exit 1
        fi
        repository_inputs+=("$repository_input")
    else
        synthetic_inputs+=("$input")
    fi
done

# 2a. Bundle contents with their repository paths: the notebook, its inputs and the modules in Synthetic_domain/,
# the repository-level inputs at their own paths, and the job helpers at the top. Credentials (.env) are never
# packed, nor the EXCLUDE_NAMES (spectral_fits/ alone carries ~330 MB of superseded 16-bin fits)
exclude_options=()
for exclude_name in "${exclude_names[@]}"; do
    exclude_options+=(--exclude="$exclude_name")
done
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
bundle_dir="$scratch_dir/bundle_contents"
mkdir -p "$bundle_dir/Synthetic_domain"
(
    cd "$SYNTHETIC_DIR"
    tar --exclude="__pycache__" --exclude=".env" "${exclude_options[@]}" \
        -cf - "$notebook" ./*.py splus-filters "${synthetic_inputs[@]}"
) | tar -C "$bundle_dir/Synthetic_domain" -xf -
if [[ ${#repository_inputs[@]} -gt 0 ]]; then
    (
        cd "$REPO_DIR"
        tar --exclude="__pycache__" --exclude=".env" -cf - "${repository_inputs[@]}"
    ) | tar -C "$bundle_dir" -xf -
fi
cp "$REMOTE_DIR/grid_md5.txt" "$REMOTE_DIR/build_env.sh" "$REMOTE_DIR/job_common.sh" "$bundle_dir/"

# 2b. Large parquet inputs as float32 (about 1e-7 relative precision, far below any noise in the data)
for float32_file in "${float32_files[@]}"; do
    "$SYNTHETIC_DIR/../.venv/bin/python" - "$bundle_dir/Synthetic_domain/$float32_file" <<'TO_FLOAT32'
import sys
import pandas as pd
table = pd.read_parquet(sys.argv[1])
float64_columns = table.select_dtypes("float64").columns
table[float64_columns] = table[float64_columns].astype("float32")
table.to_parquet(sys.argv[1], index=False)
TO_FLOAT32
    echo "Sent as float32: $float32_file ($(du -h "$bundle_dir/Synthetic_domain/$float32_file" | cut -f1))"
done

# 3. MultiNest source without the local build products, only when the notebook runs fits
if [[ "$NEEDS_MULTINEST" == "1" ]]; then
    mkdir -p "$bundle_dir/multinest"
    tar -C "$MULTINEST_SOURCE" --exclude=./build --exclude=./lib --exclude=./bin -cf - . | tar -C "$bundle_dir/multinest" -xf -
fi

# 4. The tarball
tar -C "$bundle_dir" -czf "$scratch_dir/bundle.tar.gz" .
echo "Bundle: $(du -h "$scratch_dir/bundle.tar.gz" | cut -f1)"
if [[ -n "$BUNDLE_ONLY" ]]; then
    mkdir -p "$BUNDLE_ONLY"
    cp "$scratch_dir/bundle.tar.gz" "$BUNDLE_ONLY/"
    echo "BUNDLE_ONLY set: bundle written to $BUNDLE_ONLY/bundle.tar.gz, nothing submitted"
    exit 0
fi

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
