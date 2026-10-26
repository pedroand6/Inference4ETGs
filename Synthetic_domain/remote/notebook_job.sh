#!/usr/bin/env bash
# Executes one notebook top to bottom, INSIDE an hljob job (started by submit_notebook.sh as
# `bash inputs/notebook_job.sh`).
#
# The job's working folder is /work. inputs/ holds this file and bundle.tar.gz (a copy of the repository files the
# notebook reads, in the repository layout: Synthetic_domain/..., Data/..., dustMaps/...). The notebook runs from
# Synthetic_domain/. Everything it writes or changes — the executed notebook itself, tables, figures, also outside
# Synthetic_domain/ — is copied to out/files/ with its repository path, so fetch_notebook.sh can put it back in place.
#
# Environment variables (set by submit_notebook.sh with hljob -e):
#   NOTEBOOK        notebook to execute, relative to Synthetic_domain                      [required]
#   NEEDS_GPU       set to 1 to fail at once when torch cannot use the GPU                  [default 0]
#   NEEDS_MULTINEST set to 1 to build MultiNest (only notebooks that run fits need it)     [default 0]
set -euo pipefail

# 1. Settings
NOTEBOOK="${NOTEBOOK:?NOTEBOOK must be set}"
NEEDS_GPU="${NEEDS_GPU:-0}"
NEEDS_MULTINEST="${NEEDS_MULTINEST:-0}"
JOB_DIR="/work"
export OUT_DIR="$JOB_DIR/out"
RUN_DIR="$JOB_DIR/run"
NOTEBOOK_DIR="$RUN_DIR/Synthetic_domain"   # the notebooks' own folder: their relative paths start here
KERNEL_NAME="inference4etgs"          # the kernel name stored in the notebooks' metadata
# Jupyter, to execute the notebook, and the packages only spectroscopic_reference.ipynb imports (SDSS queries and
# the CSFD dust map); versions as in the local environment
NOTEBOOK_PACKAGES=(nbconvert==7.17.2 nbclient==0.11.0 ipykernel==7.2.0 astroquery==0.4.11 dustmaps==1.0.14)

mkdir -p "$OUT_DIR" "$RUN_DIR"

# 2. Unpack the bundle
echo "== 2. Unpacking the bundle"
tar -xzf "$JOB_DIR/inputs/bundle.tar.gz" -C "$RUN_DIR"
cd "$RUN_DIR"
source job_common.sh

# 3. Environment: the usual build, with PyTorch and zuko (ml_* modules), plus NOTEBOOK_PACKAGES.
# The kernel is registered inside the venv under the name the notebook asks for.
echo "== 3. Building the environment"
if [[ "$NEEDS_MULTINEST" == "1" ]]; then
    skip_multinest=0
else
    skip_multinest=1
fi
INSTALL_TORCH=1 SKIP_MULTINEST="$skip_multinest" bash build_env.sh
source env/activate.sh
export PATH="$HOME/.local/bin:$PATH"     # uv, installed by build_env.sh
uv pip install --quiet --python "$(command -v python)" "${NOTEBOOK_PACKAGES[@]}"
python -m ipykernel install --sys-prefix --name "$KERNEL_NAME"

# 3b. GPU check: without it the notebook would silently sample on the CPU for hours
if [[ "$NEEDS_GPU" == "1" ]]; then
    if ! python -c "import torch; assert torch.cuda.is_available()"; then
        echo "ERROR: torch cannot use CUDA; the notebook would run on the CPU." >&2
        exit 1
    fi
    nvidia-smi --query-gpu=name,memory.total --format=csv
fi

# 4. Memory samples for the record; the sampler is stopped however the script ends
start_memory_sampler "$OUT_DIR/memory_samples.csv"
trap 'kill "$MEMORY_SAMPLER_PID" 2>/dev/null || true; echo "$(read_memory_peak_bytes)" > "$OUT_DIR/memory_peak_bytes.txt"' EXIT

# 5. Execute the notebook in place. A marker file taken before the run tells which files the notebook wrote.
echo "== 5. Executing $NOTEBOOK"
set_stage "notebook" "$OUT_DIR/memory_samples.csv"
touch "$RUN_DIR/.before_notebook"
notebook_exit_code=0
cd "$NOTEBOOK_DIR"
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 \
    --ExecutePreprocessor.kernel_name="$KERNEL_NAME" "$NOTEBOOK" \
    >> "$OUT_DIR/notebook.log" 2>&1 || notebook_exit_code=$?
cd "$RUN_DIR"
echo "Notebook finished with exit code $notebook_exit_code"

# 6. Copy back what the notebook wrote (also after a failure: the outputs up to the failing cell are useful).
# The environment and caches are not results. Paths are relative to the repository root; files_layout.txt tells
# fetch_notebook.sh so (jobs made before it put Synthetic_domain paths in out/files/).
echo "== 6. Files written by the notebook"
find . -type f -newer "$RUN_DIR/.before_notebook" \
    -not -path "./env/*" -not -path "./pipes/*" -not -path "./Synthetic_domain/pipes/*" \
    -not -path "*/__pycache__/*" -not -path "*/.ipynb_checkpoints/*" \
    -not -name ".before_notebook" -print | tee "$OUT_DIR/written_files.txt"
echo "repository" > "$OUT_DIR/files_layout.txt"
while read -r written_file; do
    mkdir -p "$OUT_DIR/files/$(dirname "$written_file")"
    cp -p "$written_file" "$OUT_DIR/files/$written_file"
done < "$OUT_DIR/written_files.txt"
exit "$notebook_exit_code"
