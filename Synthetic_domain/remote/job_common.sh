#!/usr/bin/env bash
# Helpers shared by fit_job.sh and benchmark_job.sh (source this file; it only defines functions).
#
# The memory numbers come from the job's cgroup (the limit the job is killed at), not from `ps`, because the
# job is killed when the cgroup passes its limit, whatever the sum of the process sizes says.

CGROUP_DIRECTORY="/sys/fs/cgroup"
MEMORY_SAMPLE_SECONDS=5

read_memory_peak_bytes() {
    # Prints the highest memory use of the whole job so far [bytes], or "unavailable".
    if [[ -r "$CGROUP_DIRECTORY/memory.peak" ]]; then
        cat "$CGROUP_DIRECTORY/memory.peak"
    else
        echo "unavailable"
    fi
}

start_memory_sampler() {
    # Starts a background loop that writes "stage, seconds, memory.current, anon" to a CSV every few seconds.
    # memory.peak alone cannot tell the fits from the environment build (compilers peak too), and "anon" is the
    # memory that cannot be reclaimed, which is what decides an out-of-memory kill. The stage name is read from
    # the file $1.stage on every sample, so the caller can label the phases.
    # Inputs: $1: path of the CSV (the stage file is "$1.stage"). Sets MEMORY_SAMPLER_PID.
    local csv_path="$1"
    echo "stage,seconds,memory_current_bytes,memory_anon_bytes" > "$csv_path"
    echo "setup" > "$csv_path.stage"
    (
        while true; do
            if [[ -r "$CGROUP_DIRECTORY/memory.current" ]]; then
                current_bytes="$(cat "$CGROUP_DIRECTORY/memory.current")"
                anon_bytes="$(awk '$1 == "anon" {print $2}' "$CGROUP_DIRECTORY/memory.stat")"
                echo "$(cat "$csv_path.stage"),$SECONDS,$current_bytes,$anon_bytes" >> "$csv_path"
            fi
            sleep "$MEMORY_SAMPLE_SECONDS"
        done
    ) &
    MEMORY_SAMPLER_PID=$!
}

set_stage() {
    # Labels the memory samples that follow. Inputs: $1: stage name, $2: path of the sampler CSV.
    echo "$1" > "$2.stage"
}
