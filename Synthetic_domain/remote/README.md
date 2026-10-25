# Running the spectral fits on the job server

The 93 fits (33 calibration mocks and 60 galaxies, `spectral_fits/jobs/`) run on the job server through the
`hljob` client (installed in `~/.local/bin`, logged in with `hljob login`; the token lives in
`~/.config/hljob/config`, never commit it). Every job starts from a fresh container, so each job builds its
own Python 3.10 environment and MultiNest (about 40 s) before fitting.

| File | Runs | Does |
|---|---|---|
| `submit_fits.sh` | locally | packs the bundle, submits a detached job, prints the job id |
| `fetch_fits.sh <ID>` | locally | downloads a job and copies fits and logs into `spectral_fits/` |
| `fit_job.sh` | in the job | unpacks, builds the environment, runs `run_spectral_fits.py`, copies results to `out/` every 10 min |
| `build_env.sh` | in the job | uv, Python 3.10, pinned packages, MultiNest, md5 check of the Bagpipes grids |
| `benchmark_job.sh`, `summarise_benchmark.py` | in the job | worker-count benchmark (see below) |
| `job_common.sh` | in the job | memory sampling from the job's cgroup |
| `grid_md5.txt` | both | md5 sums of the three Bagpipes grids of the local install; a mismatch stops the job |

## Submit

From any folder:

    Synthetic_domain/remote/submit_fits.sh

Defaults are the full run: 16 cores, 7000 MB, 24 h, `WORKERS=14`, `MAX_ITERATIONS=50000`. Change them with
environment variables, for example `WORKERS=12 TIME_LIMIT=12:00:00 submit_fits.sh`. Extra arguments go to
`hljob run` (for example `-e SYNC_SECONDS=300`). The script prints the job id.

## Follow

    hljob logs ID -f          # live log (Ctrl-C only stops watching)
    hljob logs ID --tail 30
    hljob ls                  # state of all jobs; only one job runs at a time
    hljob wait ID             # blocks until it ends
    hljob cancel ID

## Fetch

    Synthetic_domain/remote/fetch_fits.sh ID

Works on a running job (it gets the files copied to `out/` so far), after a normal end, and after a timeout.
It copies `*.h5` into `spectral_fits/` (a local file that is newer is kept), appends new rows to
`spectral_fits/fit_log.csv` and puts the build log, memory samples and job logs in `spectral_fits/job_logs/ID/`.
Then `python run_spectral_fits.py --status` (from `Synthetic_domain/`) tells how many fits are done.

## Resume

After a timeout, a cancel or a failure, fetch the job and submit again. The finished `*.h5` files of
`spectral_fits/` travel in the bundle (`finished_fits/`) and the runner skips them. A fit that was running when the
job stopped starts again from the beginning.

## Server facts (probe and benchmark)

- Container: Fedora 43, Python 3.14 (not used), no root, 16 vCPUs, 7000 MB limit with `--pause-llm`.
- `uv` installs Python 3.10 from the network; the pinned package versions are those of the local environment.
- The container has no BLAS development symlink, so `build_env.sh` points CMake at `libflexiblas.so.3`.
- Benchmark numbers (400 iterations per fit) are in the Task 12 report.

# Training the ML posterior estimator (GPU)

One GPU job trains both embeddings (MLP and CNN) one after the other with `train_ml_model.py`.

| File | Runs | Does |
|---|---|---|
| `submit_training.sh` | locally | packs the bundle (code, training and validation sets, noise model, a 3-column copy of the DR3 catalogue, MultiNest source), submits a detached GPU job, prints the job id |
| `train_job.sh` | in the job | unpacks, runs `build_env.sh` with `INSTALL_TORCH=1` (PyTorch 2.5.1 CUDA 12.1 and zuko 1.2.0), trains each embedding, writes `out/ml_runs/<embedding>/` and `out/train_<embedding>.log` |

## Submit

    Synthetic_domain/remote/submit_training.sh

Defaults: GPU, 8 cores, 7000 MB, 08:00:00, `EMBEDDINGS="mlp cnn"`. Change them with environment variables, for
example `EMBEDDINGS=mlp TIME_LIMIT=04:00:00 submit_training.sh`; extra arguments go to `hljob run`. Only one job runs
at a time, so submit it after the spectral fits have finished. `BUNDLE_ONLY=<folder> submit_training.sh` writes the
bundle there and submits nothing (to check what travels).

## Follow

    hljob logs ID -f          # the job log, including the environment build
    hljob wait ID

The training log of each embedding (one line per epoch) is `out/train_<embedding>.log`; it comes back with the fetch.

## Fetch

    hljob fetch ID -o /tmp/ml_job
    mkdir -p Synthetic_domain/ml_runs
    cp -r /tmp/ml_job/out/ml_runs/* Synthetic_domain/ml_runs/

The runs land in `Synthetic_domain/ml_runs/<embedding>/` as `model.pt`, `config.json` and `history.csv`. The weights
(`model.pt`) are ignored by git; `config.json` and `history.csv` are committed. A run that finished before a timeout
is already in `out/ml_runs/`, so it can be fetched after a timeout too.

## Test locally (CPU, a few minutes)

    cd Synthetic_domain
    python train_ml_model.py --embedding mlp --output-dir /tmp/ml_test --epochs 2 \
        --max-train-rows 3000 --max-validation-rows 500 --device cpu
