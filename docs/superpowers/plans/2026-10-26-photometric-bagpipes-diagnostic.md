# Photometric Bagpipes diagnostic (step A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fit each validation-half reference galaxy's 12-band 3″ aperture photometry with Bagpipes + MultiNest under the training prior and the simulator's likelihood. Then test whether the ML posterior agrees with these fits.

**Architecture:**
- The photometric fit lives in `synthetic_model.py` section 6, next to `fit_spectrum`. It shares `components()`, `posterior_labels`, `save_fit`/`load_fit` and the job files.
- `run_spectral_fits.py` and the remote scripts learn a job *kind* and an output folder, so the same machinery runs `photometric_fits/` on the PC.
- Notebook 2 exports the photometric jobs and reads back their posteriors into `spectroscopic_reference.csv` (`phot_<label>_p16/_median/_p84`).
- `ml_model.ipynb` compares the ML against them with the spec's "works" test.

**Tech Stack:** Python 3.10, Bagpipes 1.3.5 (bc03_miles), PyMultiNest, numpy/pandas, pytest; the hljob job service on the PC.

**Spec:** `docs/superpowers/specs/2026-10-23-ml-sfh-posterior-design.md`, section "Amendment (2026-10-26): diagnosing the gap between the ML and the spectral reference (A → B → C)". It works within the code rules of `docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md`.

## Global Constraints

- **Code style:** the ten rules of the parts 1–2 spec:
  - English names, numbered banners and why-comments;
  - no compact code;
  - docstrings with inputs, units and returns;
  - named constants;
  - logic in modules, orchestration in notebooks;
  - parallel cases in one pattern;
  - concise, not prolix;
  - lines ≤ ~110 characters.
- **Prior and simulator:** they come only from `synthetic_model.py`. Never re-implement the prior or the noise model.
- **Likelihood:** Gaussian in magnitudes with the catalogue errors and no ERROR_FLOOR, over detected bands only (`sm.detected_bands`). The magnitude offset is free with a flat prior and marginalised analytically.
- **Locked data:** only the `validation` half is fitted or compared. The `test` half and `synthetic_test.parquet` stay untouched.
- **Compute:** NEVER execute a notebook, a fit or a posterior-sampling run on the laptop. Fits and notebooks run as `hljob` jobs on the PC (`remote/submit_fits.sh`, `remote/submit_notebook.sh`). Locally only unit tests and quick checks run (seconds).
- **Secrets and safety:** never print `.env`, and never pack it into a bundle. Do not touch `Synthetic_domain/pipes/` or the repo-root PDF. Never push.

## Review Focus

1. **Undetected bands:** a galaxy with some undetected bands (99 or NaN magnitudes, or error ≥ 1 mag). Those bands drop out of the likelihood, and the fit still runs on the detected ones (≥ `MIN_DETECTED_BANDS`).
2. **Job ids:** the spectral and photometric fits of one galaxy share a specobjid. The photometric `.h5` files must never overwrite or be mistaken for spectral ones (separate folders, kind stored in the metadata).
3. **The offset:** the analytic offset marginal must equal a brute-force numerical integral over the offset on a toy case. Mass labels must shift by −0.4 Δm per sample.
4. **Notebook 2 without credentials:** on the PC, notebook 2 must run with the aperture-photometry cache complete and no `.env`. The S-PLUS login may only happen when a galaxy is missing from the cache.
5. **Spectral fits unchanged:** existing spectral job files and `.h5` files still load and run unchanged (the default kind is "spectrum").

---

### Task 1: Photometric likelihood and fit

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (section 6)
- Test: `Synthetic_domain/tests/test_synthetic_model.py`

**Interfaces:**
- Produces:
  - `PHOTOMETRIC_DIMS: list[str]`, which is `DIMS` without `"z"`.
  - `photometric_components(cube, redshift) -> (bagpipes_components, drawn)`, mapping a PHOTOMETRIC_DIMS unit-cube point through `components(..., fixed_redshift=redshift)`.
  - `offset_marginal_log_likelihood(model_magnitudes, observed_magnitudes, errors, detected) -> (log_likelihood, best_offset, offset_sigma)`. The magnitudes are 1-D over the 12 bands. The offset is observed minus model, in mag.
  - `fit_photometry(magnitudes, errors, redshift, output_basename, max_iterations, n_live_points=N_LIVE_POINTS) -> dict` with the same keys as `fit_spectrum`.
  - `photometric_posterior_labels(unit_cube_samples, redshift, magnitudes, errors, rng) -> DataFrame`. These are the labels per sample, with the mass labels rescaled by an offset drawn from N(best_offset, offset_sigma) for that sample's model.

- [ ] **Step 1:** Write the failing tests.
  - The offset marginal equals `scipy.integrate.quad` over the offset for a 3-band toy (to 1e-6). An undetected band does not change it.
  - `photometric_components` gives the same components as `components()` with the z dimension removed.
  - `photometric_posterior_labels` shifts `log_mstar` by −0.4 × offset for a known offset.
- [ ] **Step 2:** Run them and see them fail: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py -q -k photometr`.
- [ ] **Step 3:** Implement.
  - Build the Bagpipes model once with `filt_list` = the S-PLUS filters (as `posterior_labels` does) and update it in place.
  - The likelihood uses `offset_marginal_log_likelihood`.
  - Wrap the likelihood in the same `LOG_ZERO` guard as `safe_spectral_log_likelihood`, so an exception never stops MultiNest.
- [ ] **Step 4:** Run the tests and see them pass.
- [ ] **Step 5:** Add a smoke test of `fit_photometry` on a noiseless mock with tiny settings (n_live 25, max_iterations 200). Mark it `@pytest.mark.slow` if it takes more than about 30 s, and skip it if pymultinest cannot import.
- [ ] **Step 6:** Commit: `feat: photometric Bagpipes fit under the training prior (offset marginalised analytically)`.

### Task 2: Job kind, output folder and remote plumbing

**Files:**
- Modify:
  - `Synthetic_domain/synthetic_model.py`: `export_job`, `load_job`, `fit_and_save`, `save_fit` metadata.
  - `Synthetic_domain/run_spectral_fits.py`.
  - `Synthetic_domain/remote/submit_fits.sh`, `fit_job.sh` and `fetch_fits.sh`.
  - `Synthetic_domain/remote/submit_notebook.sh`, so inputs outside `Synthetic_domain/` can be bundled.
  - `Synthetic_domain/spectroscopic_reference.ipynb` cell 8: the S-PLUS login happens only for galaxies missing from the cache.
- Test: `Synthetic_domain/tests/test_synthetic_model.py`, `Synthetic_domain/tests/test_run_spectral_fits.py` (if it exists; otherwise add the tests to the former).

**Interfaces:**
- **Job dicts:** they gain `"kind"`, either `"spectrum"` (the default when absent) or `"photometry"`. A photometry job holds `magnitudes`, `errors` (12) and `redshift`.
- **Dispatch:** `fit_and_save` dispatches on the kind and stores `kind` in the `.h5` metadata.
- **Fit scripts:**
  - `submit_fits.sh` / `fit_job.sh` / `fetch_fits.sh` take `FITS_SUBDIR` (default `spectral_fits`). It names the local folder for jobs and results and the job's output folder.
  - `fetch_fits.sh <job-id>` copies into `$FITS_SUBDIR`.
- **`submit_notebook.sh`:** it packs inputs relative to the repo root when they start with `../`. The bundle keeps the repo layout (`Synthetic_domain/...`, `Data/...`, `dustMaps/...`), and the notebook runs from `Synthetic_domain/`. `notebook_job.sh` and `fetch_notebook.sh` are adjusted to match. Files the notebook writes outside `Synthetic_domain/` are brought back at their repo paths too.

- [ ] **Step 1:** Write tests.
  - An export/load round trip works for a photometry job.
  - An old job file without `kind` loads as a spectrum job.
  - `fit_and_save` routes by kind (monkeypatch both fitters).
- [ ] **Step 2:** Implement. Check the shell scripts with `bash -n`. Make a `BUNDLE_ONLY`-style dry run of `submit_notebook.sh`, or list the tarball, to show the repo layout.
- [ ] **Step 3:** Make notebook 2 cell 8 lazy: connect to S-PLUS Cloud only if some galaxy lacks cached aperture photometry. Do not execute the notebook.
- [ ] **Step 4:** Commit: `feat: photometric fit jobs share the fit machinery (job kind, FITS_SUBDIR); notebook jobs can bundle repo-level inputs`.

### Task 3: Notebook 2 — export and read the photometric fits

**Files:**
- Modify: `Synthetic_domain/spectroscopic_reference.ipynb` (a new subsection after the spectral fits, before section 7, plus the CSV columns).

**What to add:**
1. **Constants:** `PHOTOMETRIC_FITS_DIRECTORY = "photometric_fits"` and `PHOTOMETRIC_JOBS_DIRECTORY`.
2. **Export:** one photometry job per fitted galaxy of the `validation` half, from its 3″ aperture magnitudes and errors and the spectroscopic z. Also one per calibration mock: its noiseless model photometry from `sm.posterior_labels`/`model_galaxy` at the mock's cube and redshift, noised once with `sm.add_splus_noise` and a fixed seed. Job ids are the specobjid, and `calib_<k>` for the mocks.
3. **Read back:** when `photometric_fits/<id>.h5` exists, add `phot_<label>_p16`, `phot_<label>_median` and `phot_<label>_p84` columns for mw_age, Z_fin, f_old, log_mstar and log_sfr. Galaxies without a fit get NaN.
4. **Calibration of the photometric fitter:** `sm.calibration_coverage` on the photometric calibration mocks, printed in the same format as the spectral one.
5. **Markdown:** say why the subsection exists (spec amendment A).

- [ ] **Step 1:** Edit the notebook programmatically, keeping the outputs of the other cells. Validate it with nbformat and compile every code cell.
- [ ] **Step 2:** Commit: `feat: notebook 2 exports photometric fit jobs and reads their posteriors into the reference CSV`.

### Task 4: `ml_model.ipynb` — the "works" test against photometric Bagpipes

**Files:**
- Modify: `Synthetic_domain/ml_model.ipynb` (section 3).
- Modify: `Synthetic_domain/ml_evaluation.py` (add `agreement_verdict`).
- Test: `Synthetic_domain/tests/test_ml_evaluation.py`.

**Interfaces:**
- Produces `agreement_verdict(comparison_rows, labels, min_inside_95=0.80, max_abs_median_offset=1.0) -> DataFrame`. It gives one row per label with the fraction inside 95%, the median offset and `agrees`. A final verdict row says whether all labels agree.

- [ ] **Step 1:** Write failing tests for `agreement_verdict`:
  - a perfectly agreeing set gives True;
  - a set shifted by 3σ gives False;
  - a label with no galaxies gives NaN and False.
- [ ] **Step 2:** Implement it, and add these to the notebook:
  - a cell comparing the ML against the `phot_*` labels (validation half) with `ev.compare_with_spectral`;
  - the verdict for both reference sets, photometric and spectral;
  - the decision printed per the spec: "A works" or "A does not work → step B".
  - The existing spectral comparison stays.
- [ ] **Step 3:** Run the unit tests. Validate the notebook JSON. Do not execute it.
- [ ] **Step 4:** Commit: `feat: ML agreement verdict against photometric and spectral Bagpipes references`.

### Operations (controller, after the tasks)

1. Wait for the running spectral-fit job of the +108 galaxies and fetch it with `remote/fetch_fits.sh`.
2. Run notebook 2 on the PC. It refreshes the split, the calibration and the validation layers, and exports the photometric jobs. Fetch it.
3. Submit the photometric fits: `FITS_SUBDIR=photometric_fits remote/submit_fits.sh`. Fetch them, then run notebook 2 again on the PC to read them.
4. Run `ml_model.ipynb` on the PC GPU (`RUN_LOCKED_TEST=False`), fetch it, and read the verdict.
5. If A does not work, write the B plan (joint fits). If B does not work, write the C plan (reweighted prior).
