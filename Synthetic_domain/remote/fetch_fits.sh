#!/usr/bin/env bash
# Runs LOCALLY: downloads a job's results and copies the fits and logs into Synthetic_domain/spectral_fits/,
# without overwriting a local file that is newer than the downloaded one.
# Works on a running job too (it gets the files copied to out/ so far) and after a timeout.
#
# Usage:  fetch_fits.sh <job-id>
set -euo pipefail

# 1. Arguments and paths
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <job-id>" >&2
    exit 1
fi
job_id="$1"
REMOTE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FITS_DIR="$(dirname "$REMOTE_DIR")/spectral_fits"
export PATH="$HOME/.local/bin:$PATH"

# 2. Download into a scratch folder, never straight into the repository
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
hljob fetch "$job_id" -o "$scratch_dir/job"
downloaded_out="$scratch_dir/job/out"
mkdir -p "$FITS_DIR" "$FITS_DIR/job_logs"

# 3. Fits: cp -u copies only when the downloaded file is newer than the local one (or the local one is missing)
if compgen -G "$downloaded_out/*.h5" > /dev/null; then
    cp -p -u "$downloaded_out"/*.h5 "$FITS_DIR/"
    echo "Fits copied: $(ls "$downloaded_out"/*.h5 | wc -l) file(s) in the job output"
else
    echo "No .h5 files in the job output yet"
fi

# 4. Run logs. fit_log.csv of separate jobs must not overwrite one another, so each is appended to the local log
# without repeating a job_id it already lists.
if [[ -f "$downloaded_out/fit_log.csv" ]]; then
    local_log="$FITS_DIR/fit_log.csv"
    if [[ ! -f "$local_log" ]]; then
        cp "$downloaded_out/fit_log.csv" "$local_log"
    else
        # Skip the header (line 1 of the new file) and the rows whose job_id (field 1) is already in the local log
        awk -F, 'NR == FNR {known[$1] = 1; next} FNR > 1 && !($1 in known)' "$local_log" "$downloaded_out/fit_log.csv" >> "$local_log"
    fi
fi
if [[ -d "$downloaded_out/fit_logs" ]]; then
    mkdir -p "$FITS_DIR/fit_logs"
    cp -p -u "$downloaded_out"/fit_logs/*.csv "$FITS_DIR/fit_logs/"
fi

# 5. Job diagnostics (build log, memory) kept apart, one folder per job
for diagnostic in build_env.log env_build_seconds.txt memory_peak_bytes.txt memory_samples.csv; do
    if [[ -f "$downloaded_out/$diagnostic" ]]; then
        mkdir -p "$FITS_DIR/job_logs/$job_id"
        cp "$downloaded_out/$diagnostic" "$FITS_DIR/job_logs/$job_id/"
    fi
done
cp -p "$scratch_dir/job/stdout.log" "$scratch_dir/job/stderr.log" "$scratch_dir/job/job.json" "$FITS_DIR/job_logs/$job_id/" 2>/dev/null || true

echo "Done. Next: python run_spectral_fits.py --status   (from Synthetic_domain/)"
