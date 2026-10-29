# Running the spectral fits on the job server

The 201 spectral fits (33 calibration mocks and 168 galaxies, `spectral_fits/jobs/`) run on the job server through the
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

**Photometric fits** (Bagpipes fits of the 12-band aperture photometry under the same delayed-tau prior, exported by
`spectroscopic_reference.ipynb` to `photometric_fits/jobs/`: the 87 validation galaxies and 33 calibration mocks)
use the same scripts with `FITS_SUBDIR`:

    FITS_SUBDIR=photometric_fits Synthetic_domain/remote/submit_fits.sh

The double-power-law photometric fits (section 6c of notebook 2: the DPL calibration mocks, the validation galaxies
and the delayed-tau mocks' cross-fits `cross_calib_<k>`) have their own folder:

    FITS_SUBDIR=photometric_fits_dblplaw Synthetic_domain/remote/submit_fits.sh

The DPL mocks' cross-fits (fitted with delayed-tau) are exported to `photometric_fits/jobs/`, so a
`FITS_SUBDIR=photometric_fits` run after it fits only those: the finished delayed-tau fits travel in its bundle.

The double-power-law spectral fits (section 6d of notebook 2: every galaxy with a prepared spectrum, refitted
with the DPL prior, and 33 DPL calibration mocks, about 1400 fits) have their own folder too:

    FITS_SUBDIR=spectral_fits_dblplaw Synthetic_domain/remote/submit_fits.sh

They take about 30 h, so two 24 h jobs: fetch the first (`fetch_fits.sh ID`) and submit again; the second fits
only the rest, because the finished fits travel in its bundle.

They share job ids (specobjids, `calib_<k>`) with the spectral fits, so they never share a folder: the
job records its folder in `out/fits_subdir.txt`, `fetch_fits.sh` copies into that folder (and stops if a
`FITS_SUBDIR` you set differs from it), and `sm.fit_and_save` never replaces a fit of the other kind.
`skipped_fit_jobs.txt` applies to the delayed-tau spectral fits and `skipped_fit_jobs_dblplaw.txt` to the DPL
spectral fits (`SKIP_LIST` overrides either; the photometric folders have none). Each folder needs its own
list, because every folder's mocks reuse the `calib_<k>` ids. Both are empty: the ids `skipped_fit_jobs.txt`
listed named 16-bin mocks. `BUNDLE_ONLY=<folder>`
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

**Fits of the 16-bin model must never be at the top level of the fits folder.** The runner skips by file presence, so
a 16-bin `.h5` there would be bundled as a finished fit and never refitted. Every fit of the current model records
`sfh_model = "delayed"` in its metadata; `sm.fit_and_save` refuses a file without it, and
`sm.supersede_outdated_fits` moves such files (and the 16-bin mocks, calibration jobs and run logs) to
`superseded/`. The local command is in "Notebook 2" below. `fetch_fits.sh` appends a row to `fit_log.csv` only for
job ids not already in it, which is why the old log is superseded too.

## Server facts (probe and benchmark)

- Container: Fedora 43, Python 3.14 (not used), no root, 16 vCPUs, 7000 MB limit with `--pause-llm`.
- `uv` installs Python 3.10 from the network; the pinned package versions are those of the local environment.
- The container has no BLAS development symlink, so `build_env.sh` points CMake at `libflexiblas.so.3`.
- Benchmark numbers (400 iterations per fit) are in the Task 12 report.

# Training the ML posterior estimator (GPU)

One GPU job trains the MLP embedding of one SFH model with `train_ml_model.py --sfh-model <model>`:
- `SFH_MODEL=delayed`, the default, trains the delayed-τ flow (9 targets);
- `SFH_MODEL=dblplaw` trains the double-power-law flow (10 targets).

The CNN lost clearly in the 16-bin comparison and is not trained. The code still accepts `--embedding cnn`;
list it in `EMBEDDINGS` to train both.

`config.json` records each run's `sfh_model` and `target_names` (`ml_data.target_names`).
`ml_networks.load_run(directory, sfh_model)` refuses:
- a run of the other model;
- a run trained for other targets, such as the old 16-bin runs in `ml_runs/superseded/` (`mlp_16bin`, `cnn_16bin`).

A run saved without `sfh_model`, today's `ml_runs/mlp`, counts as delayed-τ.

| File | Runs | Does |
|---|---|---|
| `submit_training.sh` | locally | packs the bundle (code, the model's training and validation sets, noise model, a 3-column copy of the DR3 catalogue, MultiNest source), submits a detached GPU job, prints the job id |
| `train_job.sh` | in the job | unpacks, runs `build_env.sh` with `INSTALL_TORCH=1` (PyTorch 2.5.1 CUDA 12.1 and zuko 1.2.0), trains each embedding in `EMBEDDINGS` for `SFH_MODEL`, writes `out/ml_runs/<embedding>[_dblplaw]/` and `out/train_<embedding>_<model>.log` |

The bundle packs the model's sets that `training_set.ipynb` writes:
- `training_set.parquet` and `synthetic_validation.parquet` for delayed-τ;
- the same names with `_dblplaw` for the double power law.

Run that notebook first (next section).

## Submit

    Synthetic_domain/remote/submit_training.sh                     # delayed-τ
    SFH_MODEL=dblplaw Synthetic_domain/remote/submit_training.sh   # double power law

Defaults: GPU, 8 cores, 7000 MB, 08:00:00, `EMBEDDINGS=mlp`, `SFH_MODEL=delayed`.
- Change them with environment variables, for example `TIME_LIMIT=04:00:00 submit_training.sh`; extra arguments go to `hljob run`.
- Only one job runs at a time, so the jobs queue.
- `BUNDLE_ONLY=<folder> submit_training.sh` writes the bundle there and submits nothing, to check what travels.
- New modules must be added by hand to the `cp` list of `submit_training.sh` and to the one in `submit_fits.sh`. `submit_notebook.sh` packs every `./*.py` by itself.

## Follow

    hljob logs ID -f          # the job log, including the environment build
    hljob wait ID

The training log of each embedding (one line per epoch) is `out/train_<embedding>_<model>.log`. It comes back with the fetch.

## Fetch

    hljob fetch ID -o /tmp/ml_job
    mkdir -p Synthetic_domain/ml_runs
    cp -r /tmp/ml_job/out/ml_runs/* Synthetic_domain/ml_runs/

The runs land in `Synthetic_domain/ml_runs/<embedding>/` for delayed-τ and in `ml_runs/<embedding>_dblplaw/` for DPL. Each holds `model.pt`, `config.json` and `history.csv`.
- The weights (`model.pt`) are ignored by git; `config.json` and `history.csv` are committed.
- A run that finished before a timeout is already in `out/ml_runs/`, so it can be fetched after a timeout too.

## Test locally (CPU, a few minutes)

    cd Synthetic_domain
    python train_ml_model.py --sfh-model delayed --embedding mlp --output-dir /tmp/ml_test --epochs 2 \
        --max-train-rows 3000 --max-validation-rows 500 --device cpu

# Running a notebook on the job server

`submit_notebook.sh <notebook> <input> ...` packs the notebook, every module, `splus-filters/` and the
listed inputs (paths relative to `Synthetic_domain/`; a path starting with `../` is elsewhere in the
repository, e.g. `../dustMaps/csfd`) into a bundle with the repository layout, and the job executes the
notebook from `Synthetic_domain/`. `.env` is never packed, nor any file or folder named in `EXCLUDE_NAMES`
(default `superseded pilot`: results of an earlier model and the pilot fits, which no notebook reads).
`fetch_notebook.sh ID` puts every file the notebook wrote back at its repository path. `BUNDLE_ONLY=<folder>`
writes the bundle and submits nothing. The job
installs Jupyter and the packages only `spectroscopic_reference.ipynb` needs (astroquery, dustmaps). That
notebook logs in to the S-PLUS Cloud only for galaxies missing from `splus_aperture_photometry.csv`, so with a
complete cache it needs no `.env`.

## Notebook 1 (`training_set.ipynb`): the training and synthetic sets of each SFH model

CPU job. It reads only the DR3 catalogue (the filter curves always travel); the noise model is refitted from the
catalogue and rewritten, so the `splus_noise_model.csv` sent below is only overwritten. The delayed-tau run
(commit 0d6316d) used all 16 cores with 15 workers:

    cd Synthetic_domain
    CORES=16 MEMORY_MB=7100 TIME_LIMIT=06:00:00 HLJOB_OPTIONS="--pause-llm -e N_JOBS=15" \
        remote/submit_notebook.sh training_set.ipynb dr3_corrected.csv splus_noise_model.csv
    remote/fetch_notebook.sh ID

`SFH_MODELS` in the parameter cell lists the models to build. A model whose three set files exist in the job is
skipped (`REBUILD_SETS = False`), so a model's stored sets travel only to be kept. The run that made the
double-power-law sets (`*_dblplaw.parquet`, commit 623cd5c) added `training_set.parquet
synthetic_validation.parquet synthetic_test.parquet` to the inputs above, and the delayed-tau sets were kept.

The delayed-tau run took 12 min (2^19 + 2 x 2^14 Bagpipes models and the coverage test) with a 6.2 GB memory peak.
The notebook takes its worker count from the environment variable `N_JOBS` (default 7, for 8 cores); on 16 cores
pass it through hljob. The server refuses more than 7168 MB per job, so run it with
`HLJOB_OPTIONS="--pause-llm -e N_JOBS=15" CORES=16 MEMORY_MB=7100`. (The DPL run, job 2af348, took 16 min with
the default 7 workers.)
It writes `training_set.parquet` (125 MB), `synthetic_validation.parquet` and
`synthetic_test.parquet` (6 MB each), `splus_noise_model.csv` and the two figures, and `fetch_notebook.sh` brings
back every file it wrote. The parquet files are git-ignored and reproducible from the seeds; the CSV and the
figures are committed.

## Notebook 2 (`spectroscopic_reference.ipynb`)

The job server buffers the whole upload in a small memory limit, so keep the bundle small:

- **Firefly:** the full `../Data/sdss_eboss_firefly-dr16.fits` is 2.4 GB. Make the slim copy first (only the
  columns cell 3 reads and the spectra that can match an S-PLUS galaxy; the notebook gives the same result):

      cd Synthetic_domain && python make_slim_firefly.py

  Cell 1 reads `../Data/sdss_eboss_firefly-dr16_slim.fits` whenever it exists. Remake it when
  `dr3_corrected.csv` changes.
- **CSFD dust map:** `../dustMaps/csfd` is 289 MB (`csfd_ebv.fits` 201 MB, `mask.fits` 101 MB). It cannot be
  cut down, so if the upload fails, try `HLJOB_OPTIONS=--pause-llm` first.
- **Training set:** the coverage test reads only `redshift` and the bands, so send it as float32
  (`FLOAT32_FILES`). It must be the delayed-tau set from notebook 1, because the coverage test, and with it the
  validation/test split, is computed against it.

**Before every notebook-2 job, supersede locally** (file moves only, seconds). Cell 10 moves outdated fits,
mock files, calibration jobs and run logs to `superseded/`, but in the job it moves only the job's copies: a move
keeps the old mtime, and `notebook_job.sh` brings back only files newer than its start. Without the local step, the
old `.h5` files would be bundled by `submit_fits.sh` as finished fits (and reported failed), and the fetched
`calibration_mocks.joblib` would overwrite the old one, its only copy:

    cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -c "import synthetic_model as sm; \
        [print(d, len(sm.supersede_outdated_fits(d, m))) for d, m in sm.FIT_FOLDER_MODELS.items()]"

The notebook runs three times, in this order (the full run order of the delayed-tau model is: notebook 1, notebook 2,
spectral fits, MLP training, notebook 2, photometric fits, notebook 2, then `ml_model.ipynb`):

1. **Export.** With no spectral fits yet, it exports the 168 galaxy jobs and the 33 calibration mocks to
   `spectral_fits/jobs/` (its own superseding finds nothing left to move, because the local step ran first).
   Section 6b draws no photometric mocks and exports no photometric job, because it needs the validation half,
   and that comes from the spectral fits.
2. **After the spectral fits** are fetched: it reads them, makes the split and exports `photometric_fits/jobs/`
   (the validation galaxies and 33 mocks) with `photometric_fits/calibration_mocks.joblib`. Before it,
   `photometric_fits/jobs/` must be empty: `submit_fits.sh` fits every job file there, and a galaxy job left from
   an earlier split (the 87 of the 16-bin model were moved to `photometric_fits/superseded/jobs/` on 2026-10-26)
   would be fitted but never read. After the fetch, it holds exactly this run's export.
3. **After the photometric fits** are fetched: it reads them and writes the `phot_*` columns.

The command is the same each time. It lists whole folders, so `submit_notebook.sh` leaves out their
`superseded/` and `pilot/` subfolders (`EXCLUDE_NAMES`): the notebook reads only the top-level fits, and these
hold about 410 MB on disk (the export bundle that still carried them was 619 MB). Without them the run-3 bundle
is 390 MB, mostly the CSFD map and the float32 training set. Check the size with `BUNDLE_ONLY=<folder>` first:

    cd Synthetic_domain
    FLOAT32_FILES=training_set.parquet remote/submit_notebook.sh spectroscopic_reference.ipynb \
        dr3_corrected.csv ../Data/sdss_eboss_firefly-dr16_slim.fits ../dustMaps/csfd training_set.parquet \
        spectra spectral_fits skipped_fit_jobs.txt splus_aperture_photometry.csv splus_noise_model.csv \
        spectroscopic_reference.csv literature
    remote/fetch_notebook.sh ID

Run and fetch the spectral fits (`remote/submit_fits.sh`, `remote/fetch_fits.sh ID`) between runs 1 and 2, and the
photometric fits (`FITS_SUBDIR=photometric_fits remote/submit_fits.sh`, `remote/fetch_fits.sh ID`) between runs 2 and
3. **Run 3, which reads the photometric fits, must also bundle `photometric_fits`** (the
`.h5` files and `calibration_mocks.joblib`): add `photometric_fits` to the inputs above. Without the mock file,
cell 18 would draw new mocks whose truths differ from the fitted ones, so it stops when `calib_*` fits (or their
MultiNest folders) exist but the mock file is missing.

**DPL runs.** The notebook-2 run that exports the DPL jobs needs `photometric_fits` in its inputs (the fits,
`calibration_mocks.joblib` and, from then on, `cross_fit_mocks.joblib`); the run that reads the DPL fits also needs
`photometric_fits_dblplaw`. Between the two, fit `photometric_fits_dblplaw` and then `photometric_fits` (above), and
fetch both. The reading run writes the `phot_dblplaw_*` columns, `phot_log_evidence` and `photometric_cross_fits.csv`;
it recomputes the evidence error of the older delayed-tau fits from their samples (a few minutes).

**Locked test.** `EXPORT_LOCKED_TEST_JOBS = True` in notebook 2's parameter cell is set once, with `RUN_LOCKED_TEST` of
`ml_model.ipynb`, after the validation results have been inspected. The export run then also writes the test-half
galaxy jobs of both models, the two fits runs above fit them, and the reading run fills the test galaxies'
`phot_*` and `phot_dblplaw_*` columns. Leave it True afterwards: with False, the next reading run leaves those
columns empty for the test half again.
