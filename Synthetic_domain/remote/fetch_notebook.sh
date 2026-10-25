#!/usr/bin/env bash
# Runs LOCALLY: downloads the results of a notebook job (submit_notebook.sh) and puts every file the notebook wrote
# back at its place in Synthetic_domain — the executed notebook, its tables and figures. The job log and memory
# record go to notebook_job_logs/<job-id>/.
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
LOG_DIR="$SYNTHETIC_DIR/notebook_job_logs/$job_id"
export PATH="$HOME/.local/bin:$PATH"

# 2. Download into a scratch folder, never straight into the repository
scratch_dir="$(mktemp -d)"
trap 'rm -rf "$scratch_dir"' EXIT
hljob fetch "$job_id" -o "$scratch_dir/job"
downloaded_out="$scratch_dir/job/out"

# 3. The files the notebook wrote, at their Synthetic_domain paths
if [[ -d "$downloaded_out/files" ]]; then
    cp -rp "$downloaded_out/files/." "$SYNTHETIC_DIR/"
    echo "Files put back in Synthetic_domain:"
    sed 's/^/  /' "$downloaded_out/written_files.txt"
else
    echo "The job wrote no files (see the log below)"
fi

# 4. Logs
mkdir -p "$LOG_DIR"
for log_file in notebook.log build_env.log memory_peak_bytes.txt memory_samples.csv written_files.txt; do
    if [[ -f "$downloaded_out/$log_file" ]]; then
        cp "$downloaded_out/$log_file" "$LOG_DIR/"
    fi
done
cp -p "$scratch_dir/job/stdout.log" "$scratch_dir/job/stderr.log" "$scratch_dir/job/job.json" "$LOG_DIR/" 2>/dev/null || true
echo "Logs in $LOG_DIR"
