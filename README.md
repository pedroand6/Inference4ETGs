# Inference4ETGs — star-formation histories from S-PLUS colours

This project predicts the **star-formation history (SFH)** of a galaxy from its **S-PLUS 12-band photometry** (as colours) and its redshift. It returns a full posterior, not a single best value. The predictions are meant to agree with what **Bagpipes** (BC03-MILES, `bc03_miles`) would infer under one fixed prior. Two SFH models run side by side, each with **one metallicity** for all stars:

- a single Bagpipes **delayed-τ** component, SFR(t) ∝ (t − t₀) e^(−(t − t₀)/τ);
- Bagpipes' **double power law** (DPL, `dblplaw`), SFR(t) ∝ [(t/τ)^α + (t/τ)^−β]^−1, with t counted from the Big Bang.

Each model has its own training sets, photometric Bagpipes fits and flow, and a last section compares the two. The posterior covers:

- the SFH parameters: for delayed-τ `age` (time since star formation began), `tau` and `metallicity`; for DPL `turnover_lookback` (time from the turnover τ to the observation), `alpha`, `beta` and `metallicity`;
- the dust and nebular parameters `av_ism`, `tau_bc` and `logU`, the formed mass `log_mass_formed` and the redshift;
- labels derived from them: `log_mstar`, `log_sfr`, `log_ssfr`, `mw_age`, `f_young`, `f_mid`, `f_old`.

The predictor is a **neural posterior estimator (NPE)**: an embedding network followed by a `zuko` normalizing flow. It is trained on synthetic Bagpipes galaxies. The project tries hard to **validate honestly**: there are held-out synthetic sets, a set of real galaxies whose labels come from Bagpipes fits of their SDSS spectra (and of their S-PLUS photometry), a locked test set, and independent physical checks that do not depend on full-spectrum fitting.

> Status in one line (2026-10-26): the pipeline was rebuilt around the delayed-τ model (it replaces the earlier 16-bin non-parametric SFH, which is kept in git as the tag `nonparametric-16bin`). The code and its unit tests are done, and the synthetic sets are generated (commit `0d6316d`; the Part 1 numbers are in [§8](#8-current-status-and-results)). **The fits, the training and the evaluation (Parts 2 and 3) are still pending the PC runs**, so §8 has no results for them yet.
> 2026-10-28: the double power law was added next to delayed-τ ([spec](docs/superpowers/specs/2026-10-28-double-power-law-sfh-design.md)). Its sets are generated (commit `623cd5c`); its fits, its training and the model comparison are pending the PC runs.

---

## Contents

1. [The pipeline at a glance](#1-the-pipeline-at-a-glance)
2. [Repository map](#2-repository-map)
3. [Reading order](#3-reading-order)
4. [What each part does](#4-what-each-part-does)
5. [Key design choices and why](#5-key-design-choices-and-why)
6. [Data flow: inputs, outputs, what git ignores](#6-data-flow-inputs-outputs-what-git-ignores)
7. [How to run things](#7-how-to-run-things)
8. [Current status and results](#8-current-status-and-results)
9. [Known limitations and open decisions](#9-known-limitations-and-open-decisions)
10. [Where the history and the rationale live](#10-where-the-history-and-the-rationale-live)

---

## 1. The pipeline at a glance

```
             sfh_models.py       (the shared prior; one section per SFH model, delayed-τ and DPL; MODELS)
             synthetic_model.py  (labels, noise model, coverage test, fits; re-exports the prior)
                    │                        │                                  │
                    ▼                        ▼                                  ▼
 PART 1  training_set.ipynb        PART 2  spectroscopic_reference.ipynb   PART 3  ml_model.ipynb
 ───────────────────────────       ─────────────────────────────────────   ──────────────────────────
 Latin-hypercube draws of each     S-PLUS DR3 × Firefly DR16 sample        loads the trained MLP runs
 prior → Bagpipes → S-PLUS mags    → SDSS spectra + 3″ S-PLUS apertures    → synthetic calibration
 2^19 noise-free training models   → export fit jobs (.npz)                → "works" test against photometric
 + 2 noisy 2^14 held-out sets            │                                    and spectral Bagpipes posteriors
 + S-PLUS noise model                    ▼  (on the PC, hljob)               (validation half; test LOCKED)
 + coverage test v2                run_spectral_fits.py: MultiNest         → DR3 catalogues; model comparison
        │                          under exactly the training prior                    ▲
        │                                │                                             │
        │                                ▼                                             │
        │                          calibration test on 33 prior mocks,                 │
        │                          labels, coverage, hash split,                       │
        │                          validation rulers (indices, Hα, Gallazzi,           │
        │                          scaling relations) → spectroscopic_reference.csv ───┤
        │                          + §6b/6c photometric fits → phot_* columns ─────────┤
        │                                                                              │
        └──► train_ml_model.py (on the PC GPU, hljob): NPE = MLP + zuko flow ──────────┘
             ml_data.py · ml_networks.py · ml_evaluation.py
```

Everything lives in [`Synthetic_domain/`](Synthetic_domain/). The three numbered notebooks are the entry points. Logic lives in the Python modules, which have unit tests.

---

## 2. Repository map

Only the files that matter for the current pipeline are described in detail. Legacy material is grouped at the end.

```
Inference4ETGs/
├── README.md                         this file
├── requirements.txt                  Python dependencies (torch 2.5.1, zuko 1.2.0 pinned; bagpipes, pymultinest, …)
├── extinctions.py                    Galactic-extinction correction of S-PLUS magnitudes (CSFD E(B−V), CCM89, R_V = 3.1)
├── dustMaps/csfd/                    CSFD dust map (FITS, git-ignored)
├── docs/superpowers/
│   ├── specs/                        design documents: WHAT and WHY (with dated amendments)
│   └── plans/                        step-by-step implementation plans (HOW)
├── Data/                             mostly legacy (TNG work); holds the Firefly DR16 catalogue used by part 2
│   └── sdss_eboss_firefly-dr16.fits  Firefly DR16 catalogue (git-ignored, must be present locally)
└── Synthetic_domain/                 ← the pipeline
    ├── sfh_models.py                 the shared prior (§1), one section per SFH model (§2 delayed-τ, §3 DPL: unit cube,
    │                                 components, labels) and the registry MODELS of SfhModel entries (§4)
    ├── synthetic_model.py            shared core: re-exports the prior, labels, noise model, coverage test, spectra,
    │                                 spectral likelihood, MultiNest fits (spectral and photometric), job files and
    │                                 job kinds, calibration mocks and cross-fits, superseding of old fits
    ├── training_set.ipynb            PART 1 entry point: training set, held-out sets, noise model, coverage test
    ├── spectroscopic_reference.ipynb PART 2 entry point: real-galaxy reference set from Bagpipes spectral fits
    ├── ml_model.ipynb                PART 3 entry point: model choice, calibration, agreement, DR3 catalogue
    ├── run_spectral_fits.py          runs the exported MultiNest fit jobs (spectral or photometric kind),
    │                                 one single-core process per fit
    ├── make_slim_firefly.py          writes a small Firefly copy (needed columns, ≤ 5″ from S-PLUS) for PC uploads
    ├── slurm_fit.sh                  example Slurm job array for run_spectral_fits.py (not used on the PC)
    ├── spectral_validation.py        independent rulers: spectral indices, light-weighted age/Z, Speagle MS, bands
    ├── fetch_literature.py           downloads MPA-JHU and Gallazzi (2005) tables into literature/
    ├── ml_data.py                    NPE inputs (colours), targets (9 delayed-τ, 10 DPL), decoding to every label,
    │                                 per-epoch simulated observations
    ├── ml_networks.py                MLP (and an unused 1-D CNN) embedding, zuko spline flow, training loop,
    │                                 checkpoints that record their targets, sampling
    ├── ml_evaluation.py              coverage curves, ranks, bias/NMAD, width ratios, reference comparison,
    │                                 the `agreement_verdict` of the "works" test, model comparison and evidence
    ├── train_ml_model.py             training CLI (--sfh-model delayed|dblplaw; --embedding mlp, the default),
    │                                 run on the PC GPU
    ├── tests/                        pytest suite (266 tests) for every module above
    ├── pytest.ini                    test paths and warning filters
    ├── remote/                       scripts that run heavy work as hljob jobs on the PC (see §7)
    ├── splus-filters/                the 12 S-PLUS filter curves
    ├── literature/                   MPA-JHU / Gallazzi reference tables and relation bands (+ SOURCES.md)
    ├── ml_runs/mlp/                  the delayed-τ run (after training): config.json, history.csv (committed); model.pt (ignored).
    │                                 ml_runs/mlp_dblplaw/ is the DPL run, same layout.
    │                                 The 16-bin runs are in ml_runs/superseded/{mlp,cnn}_16bin; load_run refuses them
    ├── dr3_corrected.csv             S-PLUS DR3 catalogue, extinction-corrected (noise model, coverage, DR3 run)
    ├── splus_aperture_photometry.csv cache of raw 3″ DR4 aperture magnitudes of all 1,362 S-PLUS×Firefly matches
    ├── splus_aperture_photometry.md  how that cache was queried (tables, columns)
    ├── splus_noise_model.csv         fitted S-PLUS noise model (part 1 output)
    ├── spectroscopic_reference.csv   the reference set: labels, coverage, split, validation columns (part 2 output)
    ├── validation_summary.csv        summary of the three validation layers (part 2 output)
    ├── ml_evaluation_summary.csv     part 3 metrics of delayed-τ; ml_evaluation_summary_dblplaw.csv for DPL
    ├── ml_model_comparison.csv       part 3 delayed-τ against DPL: shared labels, spectra, evidence
    ├── photometric_cross_fits.csv    ln Z of every calibration mock under both models (part 2 §6c output)
    ├── skipped_fit_jobs.txt          spectral fit jobs deliberately left unfitted (empty for the delayed-τ model)
    ├── spectral_fits/                spectral fit jobs, selection.csv and .h5 posteriors (git-ignored);
    │                                 superseded/ holds the 16-bin fits, mocks, jobs and logs
    ├── photometric_fits/             photometric fit jobs and .h5 posteriors (notebook 2 §6b; git-ignored), same layout
    ├── photometric_fits_dblplaw/     the DPL photometric fits and their mocks (notebook 2 §6c; git-ignored)
    ├── color_color.png, coverage_criterion.png              part 1 figures (`_dblplaw` copies for DPL)
    └── validation_{indices,halpha_gallazzi,relations}.png   part 2 figures
```

**Legacy or exploratory material. The current pipeline does not use it.**

| Path | What it is |
|---|---|
| [`Data/`](Data/) notebooks (`catalogue_estimation.ipynb`, `data_estimation.ipynb`, `tng_correction.ipynb`), `*.pkl`, `*.csv`, `outputs/` | Earlier IllustrisTNG-based ETG work (comments in Portuguese). Part 2 uses only `Data/sdss_eboss_firefly-dr16.fits` (or its slim copy `sdss_eboss_firefly-dr16_slim.fits`, when present). |
| `Synthetic_domain/synthetic_domain.ipynb`, `grid.ipynb`, `simple_sampling.ipynb` | First explorations of the synthetic colour domain. `training_set.ipynb` replaced them. (The regression fixture that pinned `components` to their prior is gone: the delayed-τ prior-map tests pin the new prior.) |
| `Synthetic_domain/grid_models.csv`, `output.csv`, `Realistic_spectra.csv` | Outputs and inputs of those legacy notebooks. (`dr3_data.csv`, by contrast, is the raw DR3 table that `extinctions.py` turns into `dr3_corrected.csv`.) |
| `Bapgipes-models/` (sic) | BPASS v2.2.1 grids and `conversion.py`, which converts them to the Bagpipes format. Not used: the pipeline uses `bc03_miles`. |
| `output.png`, `.kilo/` | Old figure; editor tool folder. |
| `Synthetic_domain/pipes/`, `Data/pipes/` | Bagpipes' own working folders (`pipes/` is created by Bagpipes). Leave them alone. |

---

## 3. Reading order

To understand the project in about an hour:

1. **This README**, then the delayed-τ design spec
   [`delayed-tau-sfh-design.md`](docs/superpowers/specs/2026-10-26-delayed-tau-sfh-design.md) (with its amendment at the end), which defines the current model, and then the two earlier specs it modifies:
   [`training-and-reference-sets-design.md`](docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md) (parts 1–2) and
   [`ml-sfh-posterior-design.md`](docs/superpowers/specs/2026-10-23-ml-sfh-posterior-design.md) (part 3). Their 16-bin SFH, labels and targets are superseded; their code style, locked data, report-only evaluation and validation rulers still hold. Read their dated **amendments** at the end: they record every change of plan and the reason for it.
   The double-power-law spec [`double-power-law-sfh-design.md`](docs/superpowers/specs/2026-10-28-double-power-law-sfh-design.md) (with its amendment) adds the second model and the comparison.
2. **The opening markdown cell of each notebook**. Each one gives the goal, inputs, outputs and run order: [`training_set.ipynb`](Synthetic_domain/training_set.ipynb) → [`spectroscopic_reference.ipynb`](Synthetic_domain/spectroscopic_reference.ipynb) → [`ml_model.ipynb`](Synthetic_domain/ml_model.ipynb).
3. **[`sfh_models.py`](Synthetic_domain/sfh_models.py)**, then **[`synthetic_model.py`](Synthetic_domain/synthetic_model.py)**, section by section. Every module has numbered banners (`# 1. Prior …`), and every function has a docstring with inputs, units and returns.

**Where each concept lives.** The section numbers are the banners in [`sfh_models.py`](Synthetic_domain/sfh_models.py), written `sfh §n`, and in [`synthetic_model.py`](Synthetic_domain/synthetic_model.py), written `sm §n`.

| Concept | Where |
|---|---|
| Shared prior ranges (redshift, Z, dust, nebular) | sfh §1 (constants at the top) |
| Oldest allowed age (or turnover lookback) at a redshift | sfh §1, `maximum_age` (age of the Universe at z minus `UNIVERSE_AGE_MARGIN`) |
| Each model's ranges and unit cube → Bagpipes `delayed` or `dblplaw` component and drawn parameters | sfh §2 `delayed_components`, sfh §3 `dblplaw_components`; called as `MODELS[name].components(u, fixed_redshift=None)` |
| The model registry: `MODELS[name].dims`, `photometric_dims`, `drawn_labels`, `derived_labels`, `training_labels` (parquet and fit column order), `target_names`, `calibration_labels`, file suffix (`model_file`) | sfh §4, `SfhModel`, `MODELS` |
| Labels (`mw_age`, `f_young/mid/old`, `log_sfr`, `log_ssfr`) | delayed-τ analytic: sfh §2, `delayed_tau_labels`; DPL summed on Bagpipes' age grid: sfh §3, `dblplaw_labels`; per model: sm §2, `summary_labels` |
| AB magnitudes, model simulation, mass rescaling | sm §2, `ab_magnitudes`, `simulate_models`, `rescale_mass_labels` |
| Detection rule, S-PLUS noise model, brightness borrowing | sm §3, `detected_bands`, `fit_noise_model`, `add_splus_noise`, `scale_to_observed_r` |
| Coverage test v2 | sm §4, `coverage_test` |
| Spectrum preparation, calibration polynomial, spectral likelihood, hash split | sm §5 |
| MultiNest fit, posterior labels, mocks, `.h5` fit files | sm §6 |
| Photometric fit (offset marginalised analytically) and its evidence | sm §6, `offset_marginal_log_likelihood`, `fit_photometry`, `photometric_posterior_labels`, `fit_evidence` |
| Job kinds (`"spectrum"` / `"photometry"`) | sm §6, `job_kind`, `JOB_FITTERS`, `fit_and_save` |
| Fit-job files, calibration mocks, cross-fits and the calibration test | sm §7, `draw_photometric_calibration_mocks`, `cross_fit_jobs` |
| Superseding the fits, mocks and logs of an earlier model (each folder judged by its own `sfh_model`) | sm §8, `FIT_FOLDER_MODELS`, `supersede_outdated_fits`; `sm.fit_is_current` |
| Spectral indices, light-weighted age/Z, Speagle main sequence | [`spectral_validation.py`](Synthetic_domain/spectral_validation.py) |
| NPE inputs/targets, encoding/decoding, per-epoch simulation | [`ml_data.py`](Synthetic_domain/ml_data.py) |
| Networks, flow, training loop, `posterior_samples` | [`ml_networks.py`](Synthetic_domain/ml_networks.py) |
| Metrics (coverage, ranks, bias, NMAD, width ratios, agreement, `agreement_verdict`, model comparison, evidence) | [`ml_evaluation.py`](Synthetic_domain/ml_evaluation.py) |
| Running anything heavy on the PC | [`Synthetic_domain/remote/`](Synthetic_domain/remote/) and its [README](Synthetic_domain/remote/README.md) |

---

## 4. What each part does

### Part 1 — synthetic training set ([`training_set.ipynb`](Synthetic_domain/training_set.ipynb))

- Builds the sets of each model in `SFH_MODELS` (`["delayed", "dblplaw"]`). A model whose three set files exist is skipped unless `REBUILD_SETS = True`, so the stored delayed-τ sets are never rewritten.
- Draws **2¹⁹ = 524,288** points of the model's unit cube by Latin hypercube (`scipy.stats.qmc`, seed 0; 7 dimensions for delayed-τ, 8 for DPL). `sm.MODELS[name].components` maps each point to one Bagpipes `delayed` or `dblplaw` component by inverse CDFs. Bagpipes then computes the 12 S-PLUS magnitudes and `log_mstar` for each model. The training set is **noise-free**: noise is added during training, freshly at each epoch.
- **The delayed-τ prior** (fixed before looking at the data; `MODELS["delayed"].dims` = z, log_age_fraction, log_tau, log_metallicity, av_ism, log_tau_bc, logU):

  | parameter | prior |
  |---|---|
  | `redshift` | U(0.002, 0.7) |
  | `age` (time since star formation began) | log-uniform, 0.05 Gyr to t_max(z) = age of the Universe at z minus 0.1 Gyr. The cube holds a fraction of this log range, so one cube works at every redshift |
  | `tau` | log-uniform, 0.01–20 Gyr |
  | `metallicity` | log-uniform, 0.005–2.5 Z☉ (the `bc03_miles` grid range), one value for all stars |
  | `av_ism` | U(0.01, 2) mag (CF00 slope 0.7, birth-cloud age 0.01 Gyr) |
  | `tau_bc` (extra birth-cloud A_V < 10 Myr = 1.086 τ_BC) | log-uniform, 0.5–4 |
  | `logU` (nebular) | U(−4, −1) |
  | `log_mass_formed` | fixed at 10 for the models, then rescaled to the observed r magnitude |

- **The DPL prior** (`MODELS["dblplaw"].dims` = z, log_turnover_lookback_fraction, log_alpha, log_beta, log_metallicity, av_ism, log_tau_bc, logU):

  | parameter | prior |
  |---|---|
  | `redshift` | U(0.002, 0.7) |
  | `turnover_lookback` (time from the turnover τ to the observation) | log-uniform, 0.05 Gyr to t_max(z) = age of the Universe at z minus 0.1 Gyr; drawn as a fraction of the log range, like the delayed-τ `age` |
  | `alpha` (falling slope after the turnover) | log-uniform, 0.01–1000 (Carnall et al. 2018) |
  | `beta` (rising slope before the turnover) | log-uniform, 0.01–1000 (Carnall et al. 2018) |
  | `metallicity`, `av_ism`, `tau_bc`, `logU` | as delayed-τ |
  | `log_mass_formed` | fixed at 10, then rescaled to the observed r magnitude as for delayed-τ |

  Bagpipes' τ is the age of the Universe at z minus `turnover_lookback` (Gyr since the Big Bang), so it is always positive. Star formation starts at the Big Bang, so DPL has no age parameter.

- **Delayed-τ labels per model** (`sm.MODELS["delayed"].training_labels`, in this order):
  - the drawn parameters `redshift, age, tau, metallicity, av_ism, tau_bc, logU, log_mass_formed` (the redshift is never called `"z"`, which is the z-band column);
  - `log_age, log_tau, log_metallicity`;
  - `log_mstar`, the surviving stellar mass, from Bagpipes;
  - analytic from the delayed-τ form: `log_sfr` (mean SFR over the last 100 Myr) and `log_ssfr = log_sfr − log_mstar`, then `mw_age` (mass-weighted age, Gyr ago) and `f_young`, `f_mid`, `f_old` (formed-mass fractions younger than 0.1 Gyr, between 0.1 and 1 Gyr, older than 1 Gyr).
  - About 20% of the prior has no star formation in the last 100 Myr. For those, `log_sfr` sits exactly on a floor (10⁻¹⁴ yr⁻¹ × M*), so it stays finite. Bagpipes' own `ssfr` is never used (it is NaN there).
  - A unit test checks the analytic labels against Bagpipes' own `sfh` arrays to 1%.
- **DPL labels** (`sm.MODELS["dblplaw"].training_labels`): the drawn `redshift, turnover_lookback, tau, alpha, beta, metallicity, av_ism, tau_bc, logU, log_mass_formed`, then `log_turnover_lookback, log_alpha, log_beta, log_metallicity`, `log_mstar`, and `log_sfr, log_ssfr, mw_age, f_young, f_mid, f_old`. These last six have no closed form for DPL: they are summed on Bagpipes' own age grid (`sfh_models.dblplaw_labels`), so they match Bagpipes' discretisation of sharp peaks (α or β near 1000). The same function serves the training labels and the ML decoding.
- **S-PLUS noise model** (`splus_noise_model.csv`), fitted to `dr3_corrected.csv` (PROB_GAL > 0.8). For each band and each 0.25 mag bin, it stores 21 quantiles of log error. A band is detected when its error is below 1 mag.
- **Two noisy held-out sets** of 2¹⁴ models each, with their own seeds: `synthetic_validation.parquet` and `synthetic_test.parquet`. DPL files carry the suffix `_dblplaw` (`training_set_dblplaw.parquet`, `synthetic_validation_dblplaw.parquet`, `synthetic_test_dblplaw.parquet`), with the same sizes, seeds and noise.
- **Coverage test v2** on DR3 (§5): which real galaxies the training domain can explain.

### Part 2 — spectroscopic reference set ([`spectroscopic_reference.ipynb`](Synthetic_domain/spectroscopic_reference.ipynb))

1. **Sample.** S-PLUS DR3 is matched to Firefly DR16 within 1.5″ (1,362 matches). The selection is a 4 redshift × 3 Firefly-age grid, with the highest S/N first in each cell, 14 galaxies per cell (168). The selection file `spectral_fits/selection.csv` is append-only.
2. **Spectra.** SDSS spectra are downloaded and dereddened with the same map and law as the photometry. They are converted to air wavelengths, masked (bad pixels, 5577/6300 Å sky lines, outside 3800–9200 Å) and rebinned to 5 Å.
3. **3″ aperture photometry** from the S-PLUS Cloud, matched to the SDSS fibre and cached. It sets the mass normalisation and feeds the coverage test. It **never enters the spectral likelihood**, so the reference labels stay independent of the ML inputs.
4. **Superseding.** Cell 10 (and `sm.supersede_outdated_fits`) moves everything made under an earlier model to `superseded/` inside every folder of `sm.FIT_FOLDER_MODELS` (`spectral_fits/` and `photometric_fits/` are delayed-τ, `photometric_fits_dblplaw/` is DPL): the fits that fail `sm.fit_is_current` (a stored `sfh_model` other than the folder's, which every new fit writes, or fewer iterations than `sm.MAX_ITERATIONS`) with their MultiNest folders, the old calibration mocks and their jobs, and the run logs. The new mocks reuse the `calib_<k>` ids, so `skipped_fit_jobs.txt` is emptied.
5. **Fits.** The notebook exports one `.npz` job per galaxy and per calibration mock. The fits themselves run outside the notebook, on the PC, through [`run_spectral_fits.py`](Synthetic_domain/run_spectral_fits.py).
   - Prior: the training prior with redshift fixed (6 dimensions: log age fraction, log τ, log Z, av_ism, log τ_BC, logU), plus `veldisp` U(40, 400) km/s and a log-uniform noise scale of 1–5 (8 dimensions in total). The age upper bound uses the fixed redshift.
   - Model: the Bagpipes spectrum × an order-2 calibration polynomial, with a Gaussian likelihood.
   - Sampler: MultiNest with 300 live points, sampling efficiency 0.8 and at most 50,000 iterations.
6. **Calibration test.** 33 mocks are drawn from the full prior. The fitter counts as calibrated when, for `log_age`, `log_tau`, `log_metallicity` and `log_mstar`, at least 85% of the truths lie inside the 95% interval. The 68% coverage is reported too. Age, τ and Z are compared in log because the prior spans several decades of each.
7. **Labels, coverage and split.** Each posterior sample is re-evaluated to give all delayed-τ training labels, and the mass labels are rescaled by the aperture r offset. Galaxies that **converged and are covered** are split 50/50 into `validation` / `test` by an md5 hash of `specobjid`, which no galaxy can leave.
   **6b. Photometric Bagpipes fits (cells 17–20).** Needs the split, so it runs after the spectral fits are read. Exports one `"photometry"`-kind job per calibration mock and per validation-half galaxy into `photometric_fits/jobs/`; the prior is the training prior with redshift fixed (6 dimensions) and the brightness offset is marginalised analytically. It prints a calibration verdict in the same format as section 5. Once the fits exist, it adds `phot_<label>_p16/_median/_p84` (`sm.photometric_summary_labels("delayed")`: the calibration labels, `mw_age`, `log_sfr` and the shared labels `av_ism`, `log_ssfr`, `f_young`, `f_mid`, `f_old`), `phot_converged`, `phot_fitter_calibrated`, `phot_log_evidence` (ln Z) and `phot_log_evidence_error` (MultiNest's √(H / n_live), recomputed by `sm.fit_evidence` for fits that stored none) to `spectroscopic_reference.csv`.
   **6c. DPL photometric fits and cross-fits (cells 21–25).** The spectral reference stays delayed-τ; DPL is fitted to the aperture photometry only, for the same validation galaxies.
   - 33 DPL calibration mocks are drawn from the DPL prior, each with the redshift and brightness of a validation galaxy, and noised once. The test uses `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity` and `log_mstar`.
   - Cross-fits (`cross_calib_<k>`): the delayed-τ mocks are fitted with DPL in `photometric_fits_dblplaw/`, and the DPL mocks with delayed-τ in `photometric_fits/`. Their ln Z go to `photometric_cross_fits.csv`: since the truth is known, they show whether the colours can tell the models apart.
   - Once the fits exist, it adds `phot_dblplaw_<label>_p16/_median/_p84`, `phot_dblplaw_converged`, `phot_dblplaw_fitter_calibrated`, `phot_dblplaw_log_evidence` and `phot_dblplaw_log_evidence_error`.
   - `EXPORT_LOCKED_TEST_JOBS = True` (set once, with `RUN_LOCKED_TEST`) also exports the test-half jobs of both models.
8. **Validation of the fits** (report-only, three layers; §5). The Gallazzi and scaling-relation comparisons use the single `metallicity`. `fetch_literature.py` was re-run so MPA-JHU and Gallazzi cover all 168 galaxies.

### Part 3 — ML posterior ([`ml_model.ipynb`](Synthetic_domain/ml_model.ipynb), trained by [`train_ml_model.py`](Synthetic_domain/train_ml_model.py))

- **Inputs (37 numbers, unchanged):**
  - 11 colours (band − r), their errors and their 11 detection flags;
  - r and σ_r, used only for the absolute masses;
  - z_obs and σ_z (σ_z = 0 for a spectroscopic redshift).
- **Targets (`MODELS[name].target_names`, read through `ml_data.target_names(sfh_model)`, in this order; the redshift is last):**
  - delayed-τ, 9: `log_age`, `log_tau`, `log_metallicity`, `av_ism`, `log_tau_bc`, `logU`, `log_mass_formed`, `log_mstar`, `redshift`;
  - DPL, 10: `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity`, `av_ism`, `log_tau_bc`, `logU`, `log_mass_formed`, `log_mstar`, `redshift`.

  `log_sfr` is not a target: its floor is a point mass that a flow cannot learn, so it is derived.
  Every other label (for delayed-τ `age`, `tau`, `metallicity`, `tau_bc`, `log_sfr`, `log_ssfr`, `mw_age`, `f_young`, `f_mid`, `f_old`; for DPL the same with `turnover_lookback`, `tau`, `alpha`, `beta` in place of `age` and `tau`) is decoded from each posterior sample by `ml_data.decode_samples(theta, sfh_model, n_jobs)`, with the same functions as the training labels, at that sample's own redshift. A decoded sample outside the prior is clipped into it (age or turnover lookback into [0.05, t_max(z)]) and left unchanged when inside. DPL decoding sums over Bagpipes' grid (about 0.1 ms per sample), so it runs in parallel blocks.
- **Simulation at each epoch:** a random DR3 r magnitude at similar z (masses shift by −0.4 Δm), S-PLUS noise and non-detections, and the redshift. With probability 0.5 the redshift is spectroscopic (σ_z = 0). Otherwise it is photometric: σ_z = 0.025 (1 + z), and z_obs is drawn from N(z, σ_z).
- **Network:**
  - MLP embedding: 3 × 256 SiLU → a 64-number summary, feeding a `zuko` neural spline flow (8 transforms, 2 × 256);
  - training: AdamW, learning rate 10⁻³ with cosine decay, batch 4,096, ≤ 200 epochs, early stopping with patience 15;
  - only the MLP is trained (the CNN lost clearly in the 16-bin comparison);
  - one run per model, `train_ml_model.py --sfh-model delayed|dblplaw`: `ml_runs/mlp/` and `ml_runs/mlp_dblplaw/`;
  - a checkpoint records its `sfh_model` and `target_names`, and `ml_networks.load_run(directory, sfh_model)` raises if either differs from the requested model, so a run of the other model or a 16-bin run cannot be loaded by mistake.
- **Evaluation (report-only, sections 1–4 for each model in `SFH_MODELS`):**
  - synthetic calibration (coverage curves, rank histograms), bias and NMAD, and information per parameter (posterior/prior width), for the 9 targets and the derived labels;
  - the "works" test (`agreement_verdict`, validation half only): a label agrees when ≥ 80% of the reference medians lie inside the ML 95% interval and |median offset| ≤ 1 combined σ. The judged labels are `log_age`, `log_tau`, `log_metallicity` and `log_mstar`, with `mw_age` reported. It is run against both references: photometric Bagpipes (does the ML reproduce Bagpipes on the colours?) and spectral Bagpipes (does it agree with the spectra?). Both verdicts are printed side by side, with each fitter's calibration flag; the notebook draws no further conclusion;
  - for DPL the "works" test judges its calibration labels (`log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity`, `log_mstar`); both models also report the other shared labels;
  - a DR3 catalogue per model (`ml_sfh_catalogue.parquet`, `ml_sfh_catalogue_dblplaw.parquet`) for covered galaxies, with the `out_of_domain` and `no_r_band` flags; the metrics go to `ml_evaluation_summary.csv` and `ml_evaluation_summary_dblplaw.csv`.
- **Model comparison (section 5, report-only, to `ml_model_comparison.csv`)** on the shared labels `log_mstar`, `mw_age`, `log_metallicity`, `av_ism`, `log_ssfr`, `f_young`, `f_mid`, `f_old`:
  1. each model's ML against its own photometric Bagpipes, side by side;
  2. delayed-τ against DPL (ML on the validation galaxies and DR3, photometric Bagpipes on the validation galaxies): median difference, scatter, and the fraction of overlapping 68% intervals;
  3. four estimators (ML and photometric Bagpipes of each model) against the spectral Bagpipes, split by passive and star-forming galaxies. The spectral reference assumes delayed-τ, so its `mw_age` and mass fractions carry that SFH shape;
  4. the Bayesian evidence ln B = `phot_log_evidence` − `phot_dblplaw_log_evidence` per galaxy (positive favours delayed-τ), its categories (|ln B| < 1, 1–3, > 3), the MultiNest error √(H / n_live), and how often the 66 cross-fitted mocks prefer their true model.
- Cells marked **LOCKED** run only with `RUN_LOCKED_TEST = True`.

---

## 5. Key design choices and why

| Choice | Why |
|---|---|
| **One shared module for the prior** (`sfh_models.py`, re-exported by `synthetic_model.py`) | The training set, the spectral fits and the ML must use exactly the same prior, labels and noise model. Nothing re-implements them. Every model-dependent function takes the model by name and reads it from `MODELS`. The prior-map tests (bounds, log-uniform spacing, age ≤ t_max(z), stratified Latin hypercube) pin the mapping. |
| **Delayed-τ with one metallicity replaces the 16 bins** (user decision, 2026-10-26; the 16-bin code is the tag `nonparametric-16bin`) | The 16-bin model had 15 populated bin fractions and a prior with three SFH branches, an enrichment law and an age capped at about 10.2 Gyr by the bin edges, and its ML disagreed with the Bagpipes fits of real galaxies (§8). The working interpretation, not a tested result, is that 12 colours cannot constrain that many numbers, so the prior shaped much of the answer. The delayed-τ SFH is the parametric model Bagpipes was built around. It has three numbers (age, τ, Z), its wide log-uniform prior covers the whole parameter space, and its labels are analytic, so the training labels, the fits and the ML decoding share one set of formulas. It is replaced, not kept alongside, so nothing in the pipeline mixes the two models. The cost is a restrictive SFH shape: bursts, quenching and rejuvenation cannot be represented (out of scope). |
| **DPL added next to delayed-τ** (user decision, 2026-10-28) | Compare SFH assumptions on the same colours. Each model has its own sets, photometric fits and flow; jobs, fits and checkpoints record their model, so the two are never mixed. Delayed-τ results are kept, not recomputed, and the spectral reference stays delayed-τ. |
| **Latin hypercube over a unit cube** | It covers the 7-D prior evenly with 2¹⁹ models. The same unit cube parametrises the MultiNest fits, so "the same prior" is literally the same mapping. The age is a fraction of the log range up to t_max(z), so one cube is valid at every redshift. |
| **Age bounded by the Universe at the draw's redshift (minus 0.1 Gyr)** | Bagpipes flags a model that forms stars before the Big Bang as unphysical and returns NaN magnitudes. Fits use the spectroscopic redshift for the bound, not the unused z slot of the cube. |
| **Analytic labels from (age, τ)** | `mw_age` and the fractions are closed-form incomplete-gamma expressions, written to stay finite at the prior corners (τ = 0.01 Gyr at the oldest age, τ = 20 Gyr at 0.05 Gyr), with fractions in [0, 1] summing to 1. The ML decodes with the same functions. They are tested against Bagpipes' `sfh` arrays to 1%. |
| **`log_sfr` is derived, not a flow target** | About 20% of prior draws have no star formation in the last 100 Myr, so `log_sfr` has a point mass on its floor, which a flow cannot learn. It is computed from the decoded parameters instead. The same point masses make `f_young` and `f_old` pile up at 0 and 1: read their calibration as report-only. |
| **Age, τ and Z are compared in log** | The prior spans several decades of each, so the calibration test and the "works" test use `log_age`, `log_tau`, `log_metallicity` and `log_mstar` (the log columns sit next to the linear ones). |
| **Colours, not magnitudes, as ML inputs** | The SFH information is in the colours. Brightness only sets the mass. r enters separately for the absolute mass and SFR. |
| **Coverage test v2 with a free magnitude offset** | The old geometric test accepted 99% of galaxies because its tolerance was too wide (median 1.19 mag), and it ignored redshift. The new test uses a redshift window, χ² over the detected bands with the offset solved analytically, and covered = p > 0.01 for n_detected − 1 degrees of freedom. The free offset is the same as a colour fit with the correct error covariance. A 0.05 mag floor stands in for model error. |
| **Spectral fits under exactly the training prior** | They answer "does the ML reproduce Bagpipes on real galaxies?" with labels the ML can in principle match. No catalogue from another code (Firefly, Starlight, pPXF) can be the label, because none reproduces the Bagpipes delayed-τ model. Priors are never tuned per galaxy. |
| **Spectrum-only likelihood with a calibration polynomial** | It keeps the labels independent of the S-PLUS inputs, and it absorbs flux-calibration and aperture errors. The consequence (seen in the 16-bin run, §8) is that the fit reads the absorption features, not the continuum shape. |
| **Photometric Bagpipes fits as a second reference** | They fit the 3″ photometry the ML sees, under the simulator's likelihood. Agreement with them says whether the ML reproduces Bagpipes on colours; agreement with the spectral fits says whether colours recover what the spectra say. The two answers are different and are printed side by side. |
| **Old results are superseded by file, not read** | Every new fit records `sfh_model = "delayed"`. A fit without it is moved to `superseded/` and never read. A checkpoint without matching `target_names` raises. Tests that read the parquet sets skip with a reason while those are 16-bin sets; the PC regenerated them (commit 0d6316d), so they now run. |
| **Calibration mocks (33, from the full prior) instead of 6 fixed mocks** | The first gate (≥ 5 of 6 mocks) failed about 15% of the time even for a calibrated fitter. Debugging found no bug: the misses were prior-volume effects of draws from the prior tails. A coverage test on 33 prior draws measures calibration directly. |
| **Hash-based validation/test split** (md5 of `specobjid`) | Adding galaxies later never moves an existing one, so the test half stays untouched across passes. |
| **Locked test sets** (`synthetic_test.parquet`, the reference `test` half) | They are used once, for the final evaluation. No architecture, epoch or threshold choice may look at them. **Not yet used.** |
| **Report-only evaluation** | Every metric is printed, and the researcher judges. No automatic pass/fail gates decide the science. The one exception is the explicit "works" test (section 3b of `ml_model.ipynb`). |
| **Firefly dropped as a validation** | Firefly is another full-spectrum fitting code, so disagreeing with it cannot tell which fit is wrong. It now only provides the sample (spectrum list, redshifts, age tertiles for the selection grid). |
| **Independent rulers instead** (report-only; no galaxy is filtered) | Layer 1: model-independent spectral indices (D4000ₙ, Hδ_A, Mgb, Fe5270, Fe5335, [MgFe]′) against each fit's posterior-predictive indices. Layer 2: MPA-JHU Hα fibre SFR (BPT classes 1–2) and Gallazzi et al. (2005) light-weighted age and Z, which are built on BC03 like Bagpipes. Layer 3: scaling relations (Gallazzi mass–metallicity, Speagle 2014 main sequence, age–σ of passive galaxies). All sources are Kroupa IMF, so no IMF offset is applied. |
| **NPE with a `zuko` spline flow** | It gives full, correlated posteriors (the age–τ–Z degeneracies matter). A mixture-density network would be less flexible. The flow models the logs of age, τ, Z and τ_BC, which are log-uniform in the prior. |
| **MLP only** | The CNN was a controlled check, and it lost clearly in the 16-bin comparison (validation loss −44.61 against −65.98 for the MLP, chosen on validation data only). The code still accepts `--embedding cnn`. |
| **Redshift as the last target** | The oldest allowed age depends on z, so z is decoded first. For photometric redshifts, each sample is decoded at its own sampled z, and its age is clipped into [0.05, t_max(z)] only when it lies outside that range. |
| **Heavy compute only on the PC** | User rule: fits, notebook executions, training and posterior sampling run as `hljob` jobs on the PC (16-vCPU containers, RTX 3090). The laptop runs only unit tests and quick checks. |

---

## 6. Data flow: inputs, outputs, what git ignores

| Step | Reads | Writes | In git? |
|---|---|---|---|
| Extinction correction (`extinctions.py`, run once from the repo root) | `Synthetic_domain/dr3_data.csv`, `dustMaps/csfd/` | `Synthetic_domain/dr3_corrected.csv` (adds `*_ext` columns, corrects the `*_auto` magnitudes) | yes |
| Part 1 `training_set.ipynb` (PC CPU via `remote/submit_notebook.sh`) | `splus-filters/`, `dr3_corrected.csv` | `training_set.parquet`, `synthetic_validation.parquet`, `synthetic_test.parquet` (18 label columns and the 12 magnitudes); for DPL `training_set_dblplaw.parquet`, `synthetic_validation_dblplaw.parquet`, `synthetic_test_dblplaw.parquet` (21 label columns) | **no** (`*.parquet`, reproducible from the seeds) |
| | | `splus_noise_model.csv`, `color_color.png`, `coverage_criterion.png` (and their `_dblplaw` copies) | yes |
| `fetch_literature.py` (run once) | SkyServer SQL, MPA-JHU DR4 pages | `literature/*.csv`, `literature/relations/*.csv`, `SOURCES.md` | yes (`literature/raw/` ignored) |
| Part 2 `spectroscopic_reference.ipynb` | `dr3_corrected.csv`, `Data/sdss_eboss_firefly-dr16.fits`, `dustMaps/`, `literature/`, `training_set.parquet`, `splus_aperture_photometry.csv` (S-PLUS Cloud only for galaxies missing from it; credentials in the repo-root `.env`) | `spectra/`, `spectral_fits/jobs/*.npz`, `selection.csv`, `failures.csv`, `calibration_mocks.joblib`; moves old results to `superseded/` | **no** (`spectra/`, `spectral_fits/`) |
| Superseding before every notebook-2 job (run **locally**, `sm.supersede_outdated_fits`; see [`remote/README.md`](Synthetic_domain/remote/README.md)) | `spectral_fits/`, `photometric_fits/`, `photometric_fits_dblplaw/` | their `superseded/` subfolders | **no** |
| Fits on the PC (`run_spectral_fits.py` via `remote/submit_fits.sh`) | `spectral_fits/jobs/` | `spectral_fits/<job_id>.h5`, `fit_log.csv`, `job_logs/` | **no** |
| Part 2 analysis (same notebook, re-run) | the `.h5` fits | `spectroscopic_reference.csv`, `validation_summary.csv`, `validation_*.png` | yes |
| Photometric export (notebook 2 §6b, after the spectral fits are read) | the reference set, the spectral `calibration_mocks.joblib` | `photometric_fits/jobs/*.npz`, `photometric_fits/calibration_mocks.joblib` | **no** (`photometric_fits/`) |
| Photometric fits on the PC (`FITS_SUBDIR=photometric_fits remote/submit_fits.sh`) | `photometric_fits/jobs/` | `photometric_fits/<job_id>.h5`, `job_logs/` | **no** |
| DPL export (notebook 2 §6c) | the reference set, `photometric_fits/` | `photometric_fits_dblplaw/jobs/*.npz`, `photometric_fits_dblplaw/calibration_mocks.joblib`, `cross_fit_mocks.joblib` in both photometric folders | **no** |
| DPL photometric fits and cross-fits on the PC (`FITS_SUBDIR=photometric_fits_dblplaw`, then `photometric_fits`) | the two `jobs/` folders | `photometric_fits_dblplaw/<job_id>.h5`, `photometric_fits/cross_calib_<k>.h5` | **no** |
| Reading the DPL fits (notebook 2 §6c) | both photometric folders | `phot_dblplaw_*` and `phot_log_evidence` columns of `spectroscopic_reference.csv`, `photometric_cross_fits.csv` | yes |
| `make_slim_firefly.py` (before a notebook-2 PC run) | `Data/sdss_eboss_firefly-dr16.fits`, `dr3_corrected.csv` | `Data/sdss_eboss_firefly-dr16_slim.fits` | no (`*.fits`) |
| Training on the PC GPU (`train_ml_model.py` via `remote/submit_training.sh`, `SFH_MODEL=delayed` or `dblplaw`) | the model's `training_set*.parquet` and `synthetic_validation*.parquet`, `splus_noise_model.csv`, r/zml of `dr3_corrected.csv` | `ml_runs/<embedding>[_dblplaw]/config.json` (with `sfh_model` and `target_names`), `history.csv` | yes |
| | | `ml_runs/<embedding>[_dblplaw]/model.pt` | **no** |
| Part 3 `ml_model.ipynb` (PC GPU via `remote/submit_notebook.sh`) | `ml_runs/`, the parquet sets, `spectroscopic_reference.csv`, `photometric_cross_fits.csv`, `dr3_corrected.csv` | `ml_evaluation_summary.csv`, `ml_evaluation_summary_dblplaw.csv`, `ml_model_comparison.csv` | yes |
| | | `ml_sfh_catalogue.parquet`, `ml_sfh_catalogue_dblplaw.parquet` | **no** |
| | | `notebook_job_logs/<job-id>/` | **no** |

Also ignored ([`.gitignore`](.gitignore)): `.venv`, `*.env` (S-PLUS credentials; never print or commit them), `*.fits`, `*.pkl`, `Synthetic_domain/dr3_alldata.csv` and `.superpowers` (development ledgers, local only).

Note that `spectroscopic_reference.csv` stores `specobjid` as text. Read it with `dtype={"specobjid": str}`, because eBOSS ids exceed int64.

---

## 7. How to run things

### Local environment (the laptop: tests and light work only)

- `.venv/` holds Python 3.10 with [`requirements.txt`](requirements.txt) (Bagpipes 1.3.5, PyMultiNest, torch 2.5.1, zuko 1.2.0, …). The MultiNest shared library is built into `.venv/lib`.
- The Jupyter kernel **"Inference4ETGs (MultiNest)"** (`inference4etgs`) runs `.venv/bin/python` with `LD_LIBRARY_PATH=.../.venv/lib`, so PyMultiNest loads. Select it when you open a notebook.
- Run the tests (266 tests, CPU, about 4 to 6 minutes; without `LD_LIBRARY_PATH` the MultiNest tests are skipped; the two tests that read the generated parquet sets run now that the delayed-tau sets exist, and skip with a reason if 16-bin sets are put back):

  ```bash
  cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q
  ```

- **Do not execute notebooks, MultiNest fits, training or posterior sampling on the laptop.** These run on the PC.

### Heavy work: `hljob` jobs on the PC ([`Synthetic_domain/remote/`](Synthetic_domain/remote/))

Each job runs in a fresh container on the PC. Inside the job, `build_env.sh` rebuilds Python 3.10 with the pinned versions and MultiNest (skipped for training jobs), and checks the Bagpipes grids against `grid_md5.txt`. Only one job runs at a time. Follow a job with `hljob logs ID -f`, `hljob ls` and `hljob wait ID`.

| Task | Submit (runs locally) | Fetch the results |
|---|---|---|
| MultiNest fits of exported jobs | `Synthetic_domain/remote/submit_fits.sh` (defaults: 16 cores, 7000 MB, 24 h, `WORKERS=14`, `MAX_ITERATIONS=50000`, `FITS_SUBDIR=spectral_fits` (set `photometric_fits` for the photometric fits, `photometric_fits_dblplaw` for the DPL ones); `SKIP_LIST` defaults to `skipped_fit_jobs.txt` for `spectral_fits` only; finished `.h5` files travel along, so a resubmission resumes; `BUNDLE_ONLY=<folder>` packs without submitting) | `remote/fetch_fits.sh ID` → the folder the job recorded in `out/fits_subdir.txt` (works on running jobs too) |
| NPE training (GPU) | `Synthetic_domain/remote/submit_training.sh` (defaults: `EMBEDDINGS=mlp`, `SFH_MODEL=delayed` (set `dblplaw` for the DPL flow), 8 h; fails fast without a usable GPU) | `hljob fetch ID -o /tmp/ml_job`, then copy `out/ml_runs/*` into `Synthetic_domain/ml_runs/` |
| Execute a notebook | `Synthetic_domain/remote/submit_notebook.sh <notebook> <inputs…>` (defaults: 8 cores, 16000 MB, 4 h; `GPU=1`, `NEEDS_MULTINEST=1`, `FLOAT32_FILES` to halve big parquet uploads, `BUNDLE_ONLY`; inputs written `../path` keep their repo path; `.env` is never packed; the exact notebook-1 and notebook-2 commands are in the remote README) | `remote/fetch_notebook.sh ID` → puts the executed notebook and every file it wrote back in place |

`python run_spectral_fits.py --status` (from `Synthetic_domain/`) counts finished fits. The details, server facts and benchmark are in [`Synthetic_domain/remote/README.md`](Synthetic_domain/remote/README.md).

### Run order of the delayed-τ model

One job at a time, all on the PC except step 2:

1. `training_set.ipynb` (CPU notebook job): the training, validation and test sets and the noise model.
2. **On the laptop, before every notebook-2 job** (file moves only, seconds): `sm.supersede_outdated_fits` moves the 16-bin fits, mocks, jobs and logs to `superseded/`. A notebook job brings back only files it wrote, not its moves, so the job's own superseding never reaches the laptop. Skipping this step leaves 16-bin `.h5` files at the top level, which `submit_fits.sh` would send as finished fits, and lets the fetched `calibration_mocks.joblib` overwrite the only copy of the 16-bin mocks. The command is in the remote README.
3. `spectroscopic_reference.ipynb` (CPU): exports the spectral jobs and mocks; then the spectral fits with `remote/submit_fits.sh` (14 workers).
4. `remote/submit_training.sh` (GPU): the MLP.
5. The notebook again (after step 2 again) reads the spectral fits, makes the split and exports the photometric jobs; then the photometric fits with `FITS_SUBDIR=photometric_fits remote/submit_fits.sh`.
6. The notebook once more reads the photometric fits; then `ml_model.ipynb` (GPU) evaluates.

### Run order of the DPL model

The delayed-τ results above are kept. After the local supersede step before each notebook-2 job:

1. `training_set.ipynb` (CPU): the DPL sets (the delayed-τ sets exist, so they are skipped).
2. `SFH_MODEL=dblplaw remote/submit_training.sh` (GPU, `MEMORY_MB=7100`): the DPL MLP.
3. `spectroscopic_reference.ipynb` (CPU, a small bundle with only the fitted spectra, plus `spectral_fits`, `photometric_fits` and, once it exists, `photometric_fits_dblplaw`): exports the DPL jobs, the DPL mocks and the cross-fits. The fit folders must travel: without their mock files the notebook would draw new mocks, and fetching those would pair the stored fits with truths they were never fitted on (`fetch_notebook.sh` refuses to overwrite a differing mock file).
4. `FITS_SUBDIR=photometric_fits_dblplaw remote/submit_fits.sh` (16 cores, 15 workers), then `FITS_SUBDIR=photometric_fits` for the 33 cross-fits of the DPL mocks.
5. The notebook again, with the same inputs and `photometric_fits_dblplaw`: reads all fits.
6. `ml_model.ipynb` (`GPU=1 CORES=16 MEMORY_MB=7100 TIME_LIMIT=12:00:00 HLJOB_OPTIONS=--pause-llm`; inputs as its Inputs line, including `photometric_cross_fits.csv`): both models and the comparison.

### Typical cycle

1. Edit the code and run `pytest` locally.
2. Submit the notebook or the fits to the PC.
3. Fetch the results.
4. Read the saved outputs in VS Code.
5. Commit the small artifacts (CSVs, figures, configs). Large data stays ignored.

---

## 8. Current status and results

As of 2026-10-26. The code is rebuilt and unit-tested (159 tests), and **Part 1 has been run on the PC** (commit `0d6316d`: the executed `training_set.ipynb` and the regenerated `color_color.png` and `coverage_criterion.png`). **There are no results yet for Parts 2 and 3**: the fits, the training and the evaluation are pending the PC runs listed in §7. The committed outputs of notebooks 2 and 3, their CSVs and the part-2 figures still describe the 16-bin model and are replaced as those runs finish; the 16-bin ML runs are kept in `ml_runs/superseded/` (`mlp_16bin`, `cnn_16bin`).

### Part 1 results (delayed-τ, from the committed outputs of `training_set.ipynb`)

- **Sets:** 524,288 training models, and 16,384 each in the validation and test sets (`synthetic_validation.parquet`, `synthetic_test.parquet`), none without a magnitude offset.
- **Detected bands:** the held-out sets have 9.2 (validation) and 9.1 (test) detected bands on average, against 8.7 for DR3 (the mean over all galaxies with `PROB_GAL` above the cut, as the notebook prints it). The 51,456 DR3 galaxies kept for the coverage test have 8.9 on average; 4,632 of them have `zml` outside the prior range.
- **Coverage of DR3 by the training domain** (coverage v2, p > 0.01): **91.2% of the 51,456 DR3 galaxies**, and 93.0% of the 46,824 with `zml` inside the redshift range. The 16-bin domain covered 92.1%. The old geometric test (K = 2) gives 99.1%.
- **Synthetic validation galaxies covered:** 96.2% (the expected value is about 99%).
- **Where coverage is lowest:** `zml` in (0.7, 1.0] (72.2%, outside the prior range) and the bright end (r ≤ 17: 77.4%; 17–18: 84.1%), then the faintest bin (21.0–21.5: 86.7%). It is 89–96% for `zml` up to 0.7 and 90–94% for r between 18 and 21, and it does not depend on the number of detected bands (90–92%).

### What exists for the delayed-τ model

- The prior, components and analytic labels (`synthetic_model.py` §1–2), the 8-dimensional spectral and 6-dimensional photometric fit priors, `sfh_model` metadata and the superseding of old fits (§6–8).
- `training_set.ipynb`, `spectroscopic_reference.ipynb` and `ml_model.ipynb` rewritten for the new labels and targets, and checked only for valid JSON and compiling code (notebook 1 has since been executed on the PC, see Part 1 results above; notebooks 2 and 3 are pending).
- `ml_data.py` / `ml_networks.py` / `train_ml_model.py` with 9 targets, analytic decoding and `target_names` checkpoints.
- Literature tables (MPA-JHU, Gallazzi) refetched for all 168 galaxies.

### What the PC runs will answer

- Part 2: whether the spectral and photometric fitters pass the calibration test (95% coverage ≥ 0.85 for `log_age`, `log_tau`, `log_metallicity`, `log_mstar`, with 68% reported), the new validation/test split, and the three validation layers.
- Part 3: synthetic calibration of the MLP, the two "works" verdicts (against photometric and spectral Bagpipes) on the validation half, and the DR3 catalogue. The locked test runs once, after the validation results are inspected.

### The 16-bin model, for reference (tag `nonparametric-16bin`)

These results are superseded; they explain why the model changed.

- 168 reference galaxies (validation 87 / test 78 / excluded 3). The spectral fitter passed the calibration test on 30 of 33 mocks (3 skipped as too slow) and the photometric fitter on 33 mocks.
- The MLP (validation loss −65.98; CNN −44.61) was calibrated on synthetic data, but its posteriors **disagreed with the Bagpipes spectral fits** of real galaxies: for the 32 first-pass validation galaxies the fraction of spectral medians inside the ML 95% interval was 0.22 (`mw_age`), 0.31 (`Z_fin`), 0.28 (`log_mstar`) and 0.38 (`f_old`).
- The inputs were not the cause (units agree, DR3 total colours give the same result, the galaxies lie inside the training colour domain). The working interpretation was that the spectral fits, with their calibration polynomial, read the absorption features while 12 colours constrain the continuum shape.
- The photometric Bagpipes fits (the emulator check) agreed with the ML on `Z_fin`, `f_old` and `log_mstar`, but the ML `mw_age` interval was too narrow (0.49 of the references inside its 95%).
- The planned next step of that diagnostic was a joint spectrum + photometry reference. The user instead replaced the 16-bin model by the delayed-τ model (spec 2026-10-26), so that diagnostic chain is superseded.

---

## 9. Known limitations and open decisions

- **The spectral reference is delayed-τ only.** Its `mw_age` and mass fractions carry the delayed-τ SFH shape, so part of a DPL disagreement with the spectra is a model difference, not an ML error.
- **The delayed-τ family is restrictive.** One smooth, single-peaked SFH with one metallicity cannot represent bursts, quenching, rejuvenation or chemical evolution. A galaxy whose true SFH is outside the family gets the closest delayed-τ posterior, and the calibration tests cannot detect that.
- **Point masses in the derived labels.** About 20% of the prior has no star formation in the last 100 Myr, so `log_sfr` sits on its floor there, and `f_young` and `f_old` are often exactly 0 or 1. Coverage and rank statistics of these labels are inflated by ties; they are report-only.
- **Degeneracies** between age, τ and metallicity are likely with 12 bands, and the wide log-uniform prior may leave parts of the posterior prior-dominated. Whether the new reference comparisons pass is open until the PC runs.
- **DR3 vs DR4 photometry.** The noise model and the coverage reference use DR3. The reference-set apertures are DR4 (`dr4_dual`), whose g band is about 0.07–0.1 mag fainter than DR3 for the same objects. This is an unresolved decision: move everything to DR4, or fetch DR3 apertures.
- **Spectral-fit calibration.** In the 16-bin run the spectral fits were overconfident at 1σ (68% coverage ≈ 0.5–0.67). Whether the delayed-τ fits are is open.
- **[Mg/Fe] > 0** in passive galaxies cannot be modelled by solar-scaled BC03 (the Mgb excess seen in the 16-bin run).
- **BOSS galaxies** have 2″ fibres against 3″ apertures. They are reported separately and have no MPA-JHU or Gallazzi values.
- **Notebook outputs lag their source.** `training_set.ipynb` and its two figures (`color_color.png`, `coverage_criterion.png`) are delayed-τ. The saved outputs of notebooks 2 and 3, `spectroscopic_reference.csv`, `validation_summary.csv`, `ml_evaluation_summary.csv` and the part-2 figures still show the 16-bin pass until the PC re-runs replace them.
- **Notebook 2 on the PC** needs network access (SDSS queries) and a large upload (the CSFD map alone is 289 MB); see the remote README.

---

## 10. Where the history and the rationale live

- **Specs** ([`docs/superpowers/specs/`](docs/superpowers/specs/)) record what was decided and why. The delayed-τ model is defined by `2026-10-26-delayed-tau-sfh-design.md` and its amendment (`log_sfr` derived, log-space labels, checkpoint `target_names`); the double power law and the comparison by `2026-10-28-double-power-law-sfh-design.md` and its amendment. The two 2026-10-23 specs carry a superseded banner for their 16-bin parts; their amendments (external fitting and the calibration test; Firefly as a check, then dropped; the validation rulers; the photometric diagnostic) still explain the history.
- **Plans** ([`docs/superpowers/plans/`](docs/superpowers/plans/)) give task-by-task implementation with tests:
  - `2026-10-23-training-and-reference-sets.md` (parts 1–2);
  - `2026-10-24-ml-sfh-posterior.md` (part 3);
  - `2026-10-25-bagpipes-validation-rulers.md` (validation layers, Firefly removal);
  - `2026-10-26-photometric-bagpipes-diagnostic.md` (the photometric fits, now reused as the second reference);
  - `2026-10-26-delayed-tau-sfh.md` (the delayed-τ rewrite and the PC operations that follow it). The plans before it are historical and describe the 16-bin model;
  - `2026-10-28-double-power-law-sfh.md` (the DPL model, its fits and flow, and the comparison).
- **Ledgers** (`.superpowers/sdd/*/progress.md`, local and git-ignored) hold the run log of every task, review, ruling and result, with job ids and numbers.
- **Code style** (binding for all new code, from the parts 1–2 spec):
  - English names that say what they hold;
  - numbered section banners;
  - comments that explain *why*;
  - no compact code;
  - docstrings with inputs, units and returns;
  - UPPER_CASE constants with units;
  - logic in modules, orchestration in notebooks;
  - parallel cases in one pattern;
  - concise.
