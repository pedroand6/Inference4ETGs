#!/usr/bin/env bash
# Benchmark, INSIDE an hljob job: how long does the environment take to build, and how fast and how large are
# the fits with 10 and with 14 workers? It decides how many workers the full run should use.
#
# Each trial fits real jobs for a few hundred iterations only (--max-iterations), one job per worker, on copies
# of job files so that nothing real is touched. Results go to out/benchmark/.
#
# Environment variables: BENCHMARK_ITERATIONS [default 400]; WORKER_TRIALS [default "10 14"];
#   MEMORY_LIMIT_FOR_NEXT_TRIAL_BYTES [default 6000000000]: a larger trial runs only if the previous one
#   stayed below this memory peak.
set -euo pipefail

# 1. Settings
BENCHMARK_ITERATIONS="${BENCHMARK_ITERATIONS:-400}"
WORKER_TRIALS="${WORKER_TRIALS:-10 14}"
MEMORY_LIMIT_FOR_NEXT_TRIAL_BYTES="${MEMORY_LIMIT_FOR_NEXT_TRIAL_BYTES:-6000000000}"
JOB_DIR="/work"
export OUT_DIR="$JOB_DIR/out"
RUN_DIR="$JOB_DIR/run"
BENCHMARK_OUT="$OUT_DIR/benchmark"
MEMORY_CSV="$OUT_DIR/memory_samples.csv"

mkdir -p "$OUT_DIR" "$RUN_DIR" "$BENCHMARK_OUT"

# 2. Unpack and build, timing the build
tar -xzf "$JOB_DIR/inputs/bundle.tar.gz" -C "$RUN_DIR"
cd "$RUN_DIR"
source job_common.sh
start_memory_sampler "$MEMORY_CSV"
trap 'kill "$MEMORY_SAMPLER_PID" 2>/dev/null || true' EXIT
set_stage "build" "$MEMORY_CSV"
bash build_env.sh
source env/activate.sh
echo "Container: $(nproc) CPUs, memory limit $(cat /sys/fs/cgroup/memory.max 2>/dev/null || echo unknown) bytes" \
    | tee "$BENCHMARK_OUT/container.txt"

# 3. One trial per worker count
previous_peak_bytes=0
for workers in $WORKER_TRIALS; do
    if (( previous_peak_bytes > MEMORY_LIMIT_FOR_NEXT_TRIAL_BYTES )); then
        echo "Skipping $workers workers: the previous trial used $previous_peak_bytes bytes"
        continue
    fi
    trial_name="workers_$workers"
    trial_jobs="benchmark_jobs_$workers"       # short relative paths: MultiNest limits path length
    trial_fits="bench_$workers"
    mkdir -p "$trial_jobs" "$trial_fits"

    # 3a. One real job per worker, taken at even spacing so mocks and galaxies are both in. Copies, not links.
    all_jobs=(spectral_fits/jobs/*.npz)
    spacing=$(( ${#all_jobs[@]} / workers ))
    for (( index = 0; index < workers; index++ )); do
        cp "${all_jobs[$(( index * spacing ))]}" "$trial_jobs/"
    done

    # 3b. The trial, with the memory samples labelled by its name
    echo "== Trial: $workers workers, $BENCHMARK_ITERATIONS iterations"
    set_stage "$trial_name" "$MEMORY_CSV"
    trial_start_seconds=$SECONDS
    python run_spectral_fits.py --jobs-dir "$trial_jobs" --output-dir "$trial_fits" \
        --workers "$workers" --max-iterations "$BENCHMARK_ITERATIONS" \
        | tee "$BENCHMARK_OUT/$trial_name.stdout.txt"
    echo "$(( SECONDS - trial_start_seconds ))" > "$BENCHMARK_OUT/$trial_name.wall_seconds.txt"
    set_stage "between" "$MEMORY_CSV"

    # 3c. Keep the per-fit logs and summarise this trial's memory
    cp "$trial_fits/fit_log.csv" "$BENCHMARK_OUT/$trial_name.fit_log.csv"
    python summarise_benchmark.py "$trial_fits" "$MEMORY_CSV" "$trial_name" | tee "$BENCHMARK_OUT/$trial_name.summary.txt"
    previous_peak_bytes="$(awk -F, -v stage="$trial_name" '$1 == stage && $3 > peak {peak = $3} END {print peak + 0}' "$MEMORY_CSV")"
done

echo "memory.peak of the whole job: $(read_memory_peak_bytes)" | tee "$BENCHMARK_OUT/job_memory_peak.txt"
