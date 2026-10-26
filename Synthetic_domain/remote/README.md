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

**Photometric fits** (spec amendment A: Bagpipes fits of the 12-band aperture photometry, exported by
`spectroscopic_reference.ipynb` to `photometric_fits/jobs/`) use the same scripts with `FITS_SUBDIR`:

    FITS_SUBDIR=photometric_fits Synthetic_domain/remote/submit_fits.sh

They share job ids (specobjids, `calib_<k>`) with the spectral fits, so they never share a folder: the
job records its folder in `out/fits_subdir.txt`, `fetch_fits.sh` copies into that folder (and stops if a
`FITS_SUBDIR` you set differs from it), and `sm.fit_and_save` never replaces a fit of the other kind.
`skipped_fit_jobs.txt` applies to the spectral fits only (`SKIP_LIST` overrides it). `BUNDLE_ONLY=<folder>`
writes the bundle and submits nothing. `JOB_SCRIPT=benchmark_job.sh` benchmarks the jobs of `FITS_SUBDIR` too.

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

# Running a notebook on the job server

`submit_notebook.sh <notebook> <input> ...` packs the notebook, every module, `splus-filters/` and the
listed inputs (paths relative to `Synthetic_domain/`; a path starting with `../` is elsewhere in the
repository, e.g. `../dustMaps/csfd`) into a bundle with the repository layout, and the job executes the
notebook from `Synthetic_domain/`. `.env` is never packed. `fetch_notebook.sh ID` puts every file the notebook
wrote back at its repository path. `BUNDLE_ONLY=<folder>` writes the bundle and submits nothing. The job
installs Jupyter and the packages only `spectroscopic_reference.ipynb` needs (astroquery, dustmaps). That
notebook logs in to the S-PLUS Cloud only for galaxies missing from `splus_aperture_photometry.csv`, so with a
complete cache it needs no `.env`.

## Notebook 2 (`spectroscopic_reference.ipynb`)

The job server buffers the whole upload in a small memory limit, so keep the bundle small:

- **Firefly:** the full `../Data/sdss_eboss_firefly-dr16.fits` is 2.4 GB. Make the slim copy first (only the
  columns cell 3 reads and the spectra that can match an S-PLUS galaxy; the notebook gives the same result):

      cd Synthetic_domain && python make_slim_firefly.py

  Cell 1 reads `../Data/sdss_eboss_firefly-dr16_slim.fits` whenever it exists. Remake it when
  `dr3_corrected.csv` changes.
- **CSFD dust map:** `../dustMaps/csfd` is 289 MB (`csfd_ebv.fits` 201 MB, `mask.fits` 101 MB). It cannot be
  cut down, so if the upload fails, try `HLJOB_OPTIONS=--pause-llm` first.
- **Training set:** only `redshift` and the bands are read, so send it as float32 (`FLOAT32_FILES`).

**Before every notebook-2 job, supersede locally** (file moves only, seconds). Cell 10 moves outdated fits,
mock files, calibration jobs and run logs to `superseded/`, but in the job it moves only the job's copies: a move
keeps the old mtime, and `notebook_job.sh` brings back only files newer than its start. Without the local step, the
old `.h5` files would be bundled by `submit_fits.sh` as finished fits (and reported failed), and the fetched
`calibration_mocks.joblib` would overwrite the old one, its only copy:

    cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -c "import synthetic_model as sm; \
        [print(d, len(sm.supersede_outdated_fits(d))) for d in ('spectral_fits', 'photometric_fits')]"

First run (exports `photometric_fits/jobs/` and `photometric_fits/calibration_mocks.joblib`):

    cd Synthetic_domain
    FLOAT32_FILES=training_set.parquet remote/submit_notebook.sh spectroscopic_reference.ipynb \
        dr3_corrected.csv ../Data/sdss_eboss_firefly-dr16_slim.fits ../dustMaps/csfd training_set.parquet \
        spectra spectral_fits skipped_fit_jobs.txt splus_aperture_photometry.csv splus_noise_model.csv \
        spectroscopic_reference.csv literature
    remote/fetch_notebook.sh ID

Then run and fetch the photometric fits (`FITS_SUBDIR=photometric_fits remote/submit_fits.sh`,
`remote/fetch_fits.sh ID`). **The second run, which reads the fits, must also bundle `photometric_fits`** (the
`.h5` files and `calibration_mocks.joblib`): add `photometric_fits` to the inputs above. Without the mock file,
cell 18 would draw new mocks whose truths differ from the fitted ones, so it stops when `calib_*.h5` fits exist
but the mock file is missing.
