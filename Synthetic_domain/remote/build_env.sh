#!/usr/bin/env bash
# Builds, INSIDE a job, the software the spectral fits need: Python 3.10 with the pinned packages, and MultiNest.
#
# Why so much care: every job starts from a fresh container, and the fits must come from the same model as the
# local ones, so versions are pinned and the Bagpipes model grids are checked against the local md5 sums.
#
# Run from the folder where the bundle was unpacked (it holds multinest/, grid_md5.txt). It writes
# env/activate.sh; run `source env/activate.sh` afterwards to get the Python and the MultiNest library.
#
# Environment variables (all optional):
#   OUT_DIR   folder of the job results; the build log and the build time are kept there   [default /work/out]
set -euo pipefail

# 1. Settings
BUNDLE_DIR="$PWD"
ENV_DIR="$BUNDLE_DIR/env"
VENV_DIR="$ENV_DIR/venv"
MULTINEST_PREFIX="$ENV_DIR/multinest"            # where libmultinest.so ends up
OUT_DIR="${OUT_DIR:-/work/out}"
PYTHON_VERSION="3.10"
BUILD_JOBS="$(nproc)"
# Same versions as the notebook's environment
PINNED_PACKAGES=(
    numpy==1.26.4 scipy==1.15.3 astropy==6.1.7 pandas==2.3.3 bagpipes==1.3.5 pymultinest==2.12
    h5py==3.16.0 joblib==1.5.3 pyarrow==25.0.1 spectres==2.2.2 extinction==0.4.8
    matplotlib==3.10.8 corner==2.2.3 nautilus-sampler==1.0.6
)

mkdir -p "$OUT_DIR" "$ENV_DIR"
exec > >(tee -a "$OUT_DIR/build_env.log") 2>&1     # keep a log in out/, so it comes back with the results
build_start_seconds=$SECONDS

# 2. uv, which can install a Python version the container does not have
echo "== 2. Installing uv"
python3 -m pip install --quiet --user uv 2>/dev/null \
    || python3 -m pip install --quiet --user --break-system-packages uv   # the container's Python may be "externally managed"
export PATH="$HOME/.local/bin:$PATH"
uv --version

# 3. Python 3.10 and a virtual environment
echo "== 3. Python $PYTHON_VERSION and the virtual environment"
uv venv --quiet --python "$PYTHON_VERSION" "$VENV_DIR"

# 4. The pinned packages
echo "== 4. Installing the pinned packages"
uv pip install --quiet --python "$VENV_DIR/bin/python" "${PINNED_PACKAGES[@]}"

# 5. MultiNest, built from the shipped source. Its CMake script writes into the source folder,
# so the build happens in a copy; only the shared library is needed (the MPI one is not used).
# The container has the FlexiBLAS runtime (BLAS and LAPACK in one library) but not the development symlink that
# CMake searches for, so CMake is told the library file directly. The job is not root, so nothing can be installed.
echo "== 5. Building MultiNest"
blas_lapack_library="$(ldconfig -p | awk '/libflexiblas\.so\.3 / {print $NF; exit}')"
if [[ -z "$blas_lapack_library" ]]; then
    echo "ERROR: no BLAS/LAPACK runtime (libflexiblas.so.3) found in the container." >&2
    exit 1
fi
echo "BLAS and LAPACK: $blas_lapack_library"
multinest_build_source="$ENV_DIR/multinest_source"
cp -r "$BUNDLE_DIR/multinest" "$multinest_build_source"
mkdir -p "$multinest_build_source/build"
(
    cd "$multinest_build_source/build"
    cmake .. -DBLAS_LIBRARIES="$blas_lapack_library" -DLAPACK_LIBRARIES="$blas_lapack_library"
    make -j"$BUILD_JOBS" multinest_shared
)
mkdir -p "$MULTINEST_PREFIX/lib"
cp "$multinest_build_source"/lib/libmultinest.so* "$MULTINEST_PREFIX/lib/"

# 6. A file that sets the environment for the fits (the library folder must be on LD_LIBRARY_PATH for PyMultiNest)
cat > "$ENV_DIR/activate.sh" <<ACTIVATE
export PATH="$VENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="$MULTINEST_PREFIX/lib:\${LD_LIBRARY_PATH:-}"
ACTIVATE
source "$ENV_DIR/activate.sh"

# 7. Check that the Bagpipes grids are the ones the local fits used; a different grid would silently change the model
echo "== 7. Checking the Bagpipes grid md5 sums"
grid_directory="$(python - <<'FIND_GRIDS'
import importlib.util, os
bagpipes_folder = importlib.util.find_spec("bagpipes").submodule_search_locations[0]
print(os.path.join(bagpipes_folder, "models", "grids"))
FIND_GRIDS
)"
if ! (cd "$grid_directory" && md5sum --check "$BUNDLE_DIR/grid_md5.txt"); then
    echo "ERROR: the installed Bagpipes grids differ from the local ones (see the md5 list above)." >&2
    exit 1
fi

# 8. Versions, and a test that PyMultiNest finds the library
echo "== 8. Versions"
python --version
python - <<'PRINT_VERSIONS'
import importlib.metadata
for package in ["numpy", "scipy", "astropy", "pandas", "bagpipes", "pymultinest", "h5py", "joblib",
                "pyarrow", "spectres", "extinction", "matplotlib", "corner", "nautilus-sampler"]:
    print(f"{package}=={importlib.metadata.version(package)}")
import pymultinest    # fails here, not in the middle of a fit, if libmultinest.so is not found
print("pymultinest imports and finds libmultinest.so")
PRINT_VERSIONS

build_seconds=$((SECONDS - build_start_seconds))
echo "$build_seconds" > "$OUT_DIR/env_build_seconds.txt"
echo "Environment ready in $build_seconds s"
