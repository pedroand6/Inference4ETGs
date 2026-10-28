#!/usr/bin/env bash
# Runs LOCALLY: downloads the results of a notebook job (submit_notebook.sh) and puts every file the notebook wrote
# back at its place in the repository — the executed notebook, its tables and figures, also files outside
# Synthetic_domain. The job log and memory record go to Synthetic_domain/notebook_job_logs/<job-id>/.
#
# Usage:  fetch_notebook.sh <job-id>
set -euo pipefail

# 1. Arguments and paths
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <job-id>" >&2
    exit 1
fi
job_id="$1"
REMOTE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYNTHETIC_DIR="$(dirname "$REMOTE_DIR")"
REPO_DIR="$(dirname "$SYNTHETIC_DIR")"
LOG_DIR="$SYNTHETIC_DIR/notebook_job_logs/$job_id"
export PATH="$HOME/.local/bin:$PATH"

# 2. Download into a scratch folder, never straight into the repository
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
hljob fetch "$job_id" -o "$scratch_dir/job"
downloaded_out="$scratch_dir/job/out"

# 3. The files the notebook wrote. Jobs that record files_layout.txt list them by repository path; older jobs
# listed them by their path inside Synthetic_domain
destination_dir="$SYNTHETIC_DIR"
if [[ -f "$downloaded_out/files_layout.txt" ]]; then
    destination_dir="$REPO_DIR"
fi
if [[ -d "$downloaded_out/files" ]]; then
    # 3b. Calibration and cross-fit mocks are never replaced: the stored fits were made on the local mocks'
    # truths, so a job that drew new mocks (its bundle lacked the fit folder) must not overwrite them
    kept_mocks=0
    while IFS= read -r -d '' fetched_file; do
        local_file="$destination_dir/${fetched_file#"$downloaded_out/files/"}"
        if [[ -f "$local_file" ]] && ! cmp -s "$fetched_file" "$local_file"; then
            echo "WARNING: kept the local ${local_file#"$REPO_DIR/"}; the job's differing copy is not used" >&2
            rm "$fetched_file"
            kept_mocks=$((kept_mocks + 1))
        fi
    done < <(find "$downloaded_out/files" -type f \( -name calibration_mocks.joblib -o -name cross_fit_mocks.joblib \
        -o -path "*/jobs/calib_*.npz" -o -path "*/jobs/cross_calib_*.npz" \) -print0)
    if [[ $kept_mocks -gt 0 ]]; then
        echo "WARNING: $kept_mocks mock file(s) differed from the local ones: check that the job's bundle held the" \
            "fit folders" >&2
    fi
    cp -rp "$downloaded_out/files/." "$destination_dir/"
    echo "Files put back in $destination_dir:"
    sed 's/^/  /' "$downloaded_out/written_files.txt"
else
    echo "The job wrote no files (see the log below)"
fi

# 4. Logs
mkdir -p "$LOG_DIR"
for log_file in notebook.log build_env.log memory_peak_bytes.txt memory_samples.csv written_files.txt \
    files_layout.txt; do
    if [[ -f "$downloaded_out/$log_file" ]]; then
        cp "$downloaded_out/$log_file" "$LOG_DIR/"
    fi
done
cp -p "$scratch_dir/job/stdout.log" "$scratch_dir/job/stderr.log" "$scratch_dir/job/job.json" "$LOG_DIR/" 2>/dev/null || true
echo "Logs in $LOG_DIR"
