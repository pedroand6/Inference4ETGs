# Inference4ETGs — star-formation histories from S-PLUS colours

This project predicts the **star-formation history (SFH)** of a galaxy from its **S-PLUS 12-band photometry** (as colours) and its redshift. It returns a full posterior, not a single best value. The predictions are meant to agree with what **Bagpipes** (BC03-MILES, `bc03_miles`) would infer under one fixed prior:

- the mass fraction and the SFR in each of **16 age bins** (only 15 are ever populated, see below);
- the metallicity of each bin (set by one number, `Z_fin`);
- summary labels: `log_mstar`, `log_sfr`, `mw_age`, `f_young`, `f_mid`, `f_old`.

The predictor is a **neural posterior estimator (NPE)**: an embedding network followed by a `zuko` normalizing flow. It is trained on synthetic Bagpipes galaxies. The project tries hard to **validate honestly**: there are held-out synthetic sets, a set of real galaxies whose labels come from Bagpipes fits of their SDSS spectra, a locked test set, and independent physical checks that do not depend on full-spectrum fitting.

> Status in one line (2026-10-26): on synthetic data the ML posterior is well calibrated, but it **disagrees with the Bagpipes spectral fits of real galaxies**. Step A of the diagnostic chain A → B → C (photometric Bagpipes fits) is coded but has not been run yet ([§8](#8-current-status-and-results)).

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
             synthetic_model.py  (the ONE prior, SFH branches, labels, noise model, coverage test, fits)
                    │                        │                                  │
                    ▼                        ▼                                  ▼
 PART 1  training_set.ipynb        PART 2  spectroscopic_reference.ipynb   PART 3  ml_model.ipynb
 ───────────────────────────       ─────────────────────────────────────   ──────────────────────────
 Latin-hypercube draws of the      S-PLUS DR3 × Firefly DR16 sample        loads trained runs (MLP, CNN)
 prior → Bagpipes → S-PLUS mags    → SDSS spectra + 3″ S-PLUS apertures    → model choice (validation only)
 2^19 noise-free training models   → export fit jobs (.npz)                → synthetic calibration
 + 2 noisy 2^14 held-out sets            │                                 → agreement with Bagpipes
 + S-PLUS noise model                    ▼  (on the PC, hljob)               spectral posteriors
 + coverage test v2                run_spectral_fits.py: MultiNest           (validation half; test LOCKED)
        │                          under exactly the training prior        → DR3 catalogue of posteriors
        │                                │                                          ▲
        │                                ▼                                          │
        │                          calibration test on 33 prior mocks,              │
        │                          labels, coverage, hash split,                    │
        │                          validation rulers (indices, Hα, Gallazzi,        │
        │                          scaling relations) → spectroscopic_reference.csv─┤
        │                          + §6b photometric fits (step A) → phot_* columns ─┤
        │                                                                           │
        └──► train_ml_model.py (on the PC GPU, hljob): NPE = embedding + zuko flow ─┘
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
    ├── synthetic_model.py            shared core: prior, SFH branches, labels, noise model, coverage test,
    │                                 spectra, spectral likelihood, MultiNest fits (spectral and photometric),
    │                                 job files and job kinds, calibration test
    ├── training_set.ipynb            PART 1 entry point: training set, held-out sets, noise model, coverage test
    ├── spectroscopic_reference.ipynb PART 2 entry point: real-galaxy reference set from Bagpipes spectral fits
    ├── ml_model.ipynb                PART 3 entry point: model choice, calibration, agreement, DR3 catalogue
    ├── run_spectral_fits.py          runs the exported MultiNest fit jobs (spectral or photometric kind),
    │                                 one single-core process per fit
    ├── make_slim_firefly.py          writes a small Firefly copy (needed columns, ≤ 5″ from S-PLUS) for PC uploads
    ├── slurm_fit.sh                  example Slurm job array for run_spectral_fits.py (not used on the PC)
    ├── spectral_validation.py        independent rulers: spectral indices, light-weighted age/Z, Speagle MS, bands
    ├── fetch_literature.py           downloads MPA-JHU and Gallazzi (2005) tables into literature/
    ├── ml_data.py                    NPE inputs (colours) and targets (19 numbers), per-epoch simulated observations
    ├── ml_networks.py                MLP / 1-D CNN embeddings, zuko spline flow, training loop, checkpoints, sampling
    ├── ml_evaluation.py              coverage curves, ranks, bias/NMAD, width ratios, reference comparison,
    │                                 the A→B→C `agreement_verdict`
    ├── train_ml_model.py             training CLI (--embedding mlp|cnn), run on the PC GPU
    ├── tests/                        pytest suite (125 tests) for every module above
    ├── pytest.ini                    test paths and warning filters
    ├── remote/                       scripts that run heavy work as hljob jobs on the PC (see §7)
    ├── splus-filters/                the 12 S-PLUS filter curves
    ├── literature/                   MPA-JHU / Gallazzi reference tables and relation bands (+ SOURCES.md)
    ├── ml_runs/{mlp,cnn}/            trained runs: config.json, history.csv (committed); model.pt (ignored)
    ├── dr3_corrected.csv             S-PLUS DR3 catalogue, extinction-corrected (noise model, coverage, DR3 run)
    ├── splus_aperture_photometry.csv cache of raw 3″ DR4 aperture magnitudes of all 1,362 S-PLUS×Firefly matches
    ├── splus_aperture_photometry.md  how that cache was queried (tables, columns)
    ├── splus_noise_model.csv         fitted S-PLUS noise model (part 1 output)
    ├── spectroscopic_reference.csv   the reference set: labels, coverage, split, validation columns (part 2 output)
    ├── validation_summary.csv        summary of the three validation layers (part 2 output)
    ├── ml_evaluation_summary.csv     part 3 metrics (part 3 output)
    ├── skipped_fit_jobs.txt          spectral fit jobs deliberately left unfitted (3 very slow calibration mocks)
    ├── spectral_fits/                spectral fit jobs, selection.csv and .h5 posteriors (git-ignored)
    ├── photometric_fits/             step-A fit jobs and .h5 posteriors (written by notebook 2 §6b; not created yet)
    ├── color_color.png, coverage_criterion.png              part 1 figures
    └── validation_{indices,halpha_gallazzi,relations}.png   part 2 figures
```

**Legacy or exploratory material. The current pipeline does not use it.**

| Path | What it is |
|---|---|
| [`Data/`](Data/) notebooks (`catalogue_estimation.ipynb`, `data_estimation.ipynb`, `tng_correction.ipynb`), `*.pkl`, `*.csv`, `outputs/` | Earlier IllustrisTNG-based ETG work (comments in Portuguese). Part 2 uses only `Data/sdss_eboss_firefly-dr16.fits` (or its slim copy `sdss_eboss_firefly-dr16_slim.fits`, when present). |
| `Synthetic_domain/synthetic_domain.ipynb`, `grid.ipynb`, `simple_sampling.ipynb` | First explorations of the synthetic colour domain. `training_set.ipynb` replaced them; the regression test in `tests/legacy_components.json` pins `synthetic_model.components` to the legacy prior. |
| `Synthetic_domain/grid_models.csv`, `output.csv`, `Realistic_spectra.csv` | Outputs and inputs of those legacy notebooks. (`dr3_data.csv`, by contrast, is the raw DR3 table that `extinctions.py` turns into `dr3_corrected.csv`.) |
| `Bapgipes-models/` (sic) | BPASS v2.2.1 grids and `conversion.py`, which converts them to the Bagpipes format. Not used: the pipeline uses `bc03_miles`. |
| `output.png`, `.kilo/` | Old figure; editor tool folder. |
| `Synthetic_domain/pipes/`, `Data/pipes/` | Bagpipes' own working folders (`pipes/` is created by Bagpipes). Leave them alone. |

---

## 3. Reading order

To understand the project in about an hour:

1. **This README**, then the two design specs, in this order:
   [`training-and-reference-sets-design.md`](docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md) (parts 1–2) and
   [`ml-sfh-posterior-design.md`](docs/superpowers/specs/2026-10-23-ml-sfh-posterior-design.md) (part 3). Read their dated **amendments** at the end: they record every change of plan and the reason for it.
2. **The opening markdown cell of each notebook**. Each one gives the goal, inputs, outputs and run order: [`training_set.ipynb`](Synthetic_domain/training_set.ipynb) → [`spectroscopic_reference.ipynb`](Synthetic_domain/spectroscopic_reference.ipynb) → [`ml_model.ipynb`](Synthetic_domain/ml_model.ipynb).
3. **[`synthetic_model.py`](Synthetic_domain/synthetic_model.py)**, section by section. Every module has numbered banners (`# 1. Prior …`), and every function has a docstring with inputs, units and returns.

**Where each concept lives.** The section numbers are the banners in [`synthetic_model.py`](Synthetic_domain/synthetic_model.py), written `sm §n`.

| Concept | Where |
|---|---|
| Prior ranges, age bins, unit-cube dimensions `DIMS` | sm §1 (constants at the top) |
| SFH branches (Dirichlet, skew-normal, skew-normal + burst) | sm §1, `dirichlet_fractions`, `skew_normal_fractions`, `skew_normal_burst_fractions`, `choose_sfh_branch` |
| Enrichment law | sm §1, `bin_metallicities` |
| Unit cube → Bagpipes components | sm §1, `components(u, fixed_redshift=None)` |
| Labels (16-bin and summary), AB magnitudes, mass rescaling | sm §2, `bin_labels`, `summary_labels`, `rescale_mass_labels` |
| Detection rule, S-PLUS noise model, brightness borrowing | sm §3, `detected_bands`, `fit_noise_model`, `add_splus_noise`, `scale_to_observed_r` |
| Coverage test v2 | sm §4, `coverage_test` |
| Spectrum preparation, calibration polynomial, spectral likelihood, hash split | sm §5 |
| MultiNest fit, posterior labels, mocks, `.h5` fit files | sm §6 |
| Photometric fit of step A (offset marginalised analytically) | sm §6, `offset_marginal_log_likelihood`, `fit_photometry`, `photometric_posterior_labels` |
| Job kinds (`"spectrum"` / `"photometry"`) | sm §6, `job_kind`, `JOB_FITTERS`, `fit_and_save` |
| Fit-job files, calibration mocks and the calibration test | sm §7 |
| Spectral indices, light-weighted age/Z, Speagle main sequence | [`spectral_validation.py`](Synthetic_domain/spectral_validation.py) |
| NPE inputs/targets, encoding/decoding, per-epoch simulation | [`ml_data.py`](Synthetic_domain/ml_data.py) |
| Networks, flow, training loop, `posterior_samples` | [`ml_networks.py`](Synthetic_domain/ml_networks.py) |
| Metrics (coverage, ranks, bias, NMAD, width ratios, agreement, `agreement_verdict`) | [`ml_evaluation.py`](Synthetic_domain/ml_evaluation.py) |
| Running anything heavy on the PC | [`Synthetic_domain/remote/`](Synthetic_domain/remote/) and its [README](Synthetic_domain/remote/README.md) |

---

## 4. What each part does

### Part 1 — synthetic training set ([`training_set.ipynb`](Synthetic_domain/training_set.ipynb))

- Draws **2¹⁹ = 524,288** points of a 28-dimensional unit cube by Latin hypercube (`scipy.stats.qmc`, seed 0). `sm.components` maps each point to Bagpipes components by inverse CDFs. Bagpipes then computes the 12 S-PLUS magnitudes and the labels for each model. The training set is **noise-free**: noise is added during training, freshly at each epoch.
- **The prior** (fixed before looking at the data):

  | parameter | prior |
  |---|---|
  | redshift | U(0.002, 0.7) |
  | SFH branch | ~50% Dirichlet(1) bin fractions; ~46% skew-normal (peak ~ N(12, 3) Gyr ago, width ~ N(1, 0.4) Gyr, shape U(−0.25, 0.75)); ~4% skew-normal + burst (10–500 Myr ago, 40–50 Myr wide, 1–10% of the mass) |
  | Z_fin | log-uniform 0.1–2.5 Z☉ |
  | ISM A_V | U(0.01, 2) mag (CF00 slope 0.7) |
  | extra birth-cloud A_V (< 10 Myr) | 1.086 τ_BC, τ_BC log-uniform 0.5–4 |
  | log U (nebular) | U(−4, −1) |

- **Labels per model:** `frac_bin_00…15`, `logsfr_bin_00…15` (for 10¹⁰ M☉ formed), `zmet_bin_00…15`, and the summary labels. Bins older than the Universe at that redshift get `frac = 0` and NaN SFR and metallicity.
- **S-PLUS noise model** (`splus_noise_model.csv`), fitted to `dr3_corrected.csv` (PROB_GAL > 0.8). For each band and each 0.25 mag bin, it stores 21 quantiles of log error. A band is detected when its error is below 1 mag.
- **Two noisy held-out sets** of 2¹⁴ models each, with their own seeds: `synthetic_validation.parquet` and `synthetic_test.parquet`.
- **Coverage test v2** on DR3 (§5): which real galaxies the training domain can explain.

### Part 2 — spectroscopic reference set ([`spectroscopic_reference.ipynb`](Synthetic_domain/spectroscopic_reference.ipynb))

1. **Sample.** S-PLUS DR3 is matched to Firefly DR16 within 1.5″ (1,362 matches). The selection is a 4 redshift × 3 Firefly-age grid, with the highest S/N first in each cell. It started at 5 galaxies per cell (60) and now uses 14 per cell (168). The selection file `spectral_fits/selection.csv` is append-only.
2. **Spectra.** SDSS spectra are downloaded and dereddened with the same map and law as the photometry. They are converted to air wavelengths, masked (bad pixels, 5577/6300 Å sky lines, outside 3800–9200 Å) and rebinned to 5 Å.
3. **3″ aperture photometry** from the S-PLUS Cloud, matched to the SDSS fibre and cached. It sets the mass normalisation and feeds the coverage test. It **never enters the spectral likelihood**, so the reference labels stay independent of the ML inputs.
4. **Fits.** The notebook exports one `.npz` job per galaxy and per calibration mock. The fits themselves run outside the notebook, on the PC, through [`run_spectral_fits.py`](Synthetic_domain/run_spectral_fits.py).
   - Prior: the training prior with redshift fixed (27 dimensions), plus `veldisp` U(40, 400) km/s and a log-uniform noise scale of 1–5 (29 dimensions in total).
   - Model: the Bagpipes spectrum × an order-2 calibration polynomial, with a Gaussian likelihood.
   - Sampler: MultiNest with 300 live points, sampling efficiency 0.8 and at most 50,000 iterations.
5. **Calibration test.** 33 mocks are drawn from the full prior. The fitter counts as calibrated when, for `mw_age`, `Z_fin` and `f_old`, at least 85% of the truths lie inside the 95% interval.
6. **Labels, coverage and split.** Each posterior sample is re-evaluated to give the labels, and the mass labels are rescaled by the aperture r offset. Galaxies that **converged and are covered** are split 50/50 into `validation` / `test` by an md5 hash of `specobjid`.
   **6b. Photometric Bagpipes fits (step A, cells 17–20).** Exports one `"photometry"`-kind job per calibration mock and per validation-half galaxy into `photometric_fits/jobs/`. It prints a calibration verdict in the same format as section 5. Once the fits exist, it adds `phot_<label>_p16/_median/_p84` (mw_age, Z_fin, f_old, log_mstar, log_sfr), `phot_converged` and `phot_fitter_calibrated` to `spectroscopic_reference.csv`.
7. **Validation of the fits** (report-only, three layers; §5).

### Part 3 — ML posterior ([`ml_model.ipynb`](Synthetic_domain/ml_model.ipynb), trained by [`train_ml_model.py`](Synthetic_domain/train_ml_model.py))

- **Inputs (37 numbers):**
  - 11 colours (band − r), their errors and their 11 detection flags;
  - r and σ_r, used only for the absolute masses;
  - z_obs and σ_z (σ_z = 0 for a spectroscopic redshift).
- **Targets (19 numbers):**
  - 14 additive log-ratios of the 15 populated bin fractions, with bin 13 as reference (it exists at every redshift of the prior);
  - log Z_fin;
  - ISM A_V;
  - log M_formed and log M_living;
  - redshift.
  Per-bin SFRs, metallicities, `mw_age` and `f_*` are derived from each posterior sample with the `synthetic_model` functions, at that sample's redshift.
- **Simulation at each epoch:** a random DR3 r magnitude at similar z (masses shift by −0.4 Δm), S-PLUS noise and non-detections, and the redshift. With probability 0.5 the redshift is spectroscopic (σ_z = 0). Otherwise it is photometric: σ_z = 0.025 (1 + z), and z_obs is drawn from N(z, σ_z).
- **Networks:**
  - MLP embedding: 3 × 256 SiLU → a 64-number summary;
  - 1-D CNN embedding over the bands ordered by wavelength;
  - both feed a `zuko` neural spline flow (8 transforms, 2 × 256);
  - training: AdamW, learning rate 10⁻³ with cosine decay, batch 4,096, ≤ 200 epochs, early stopping with patience 15.
- **Evaluation (report-only):**
  - synthetic calibration (coverage curves, rank histograms);
  - bias and NMAD;
  - information per bin (posterior/prior width);
  - agreement with the spectral posteriors of the reference set (section 3);
  - the A→B→C "works" verdict (section 3b, `agreement_verdict`) against the photometric and the spectral Bagpipes references, on the same ML draws. It reports "inconclusive" unless the photometric fitter is calibrated and every validation galaxy has a fit;
  - a DR3 catalogue (`ml_sfh_catalogue.parquet`) for covered galaxies, with the `out_of_domain` and `no_r_band` flags.
- Cells marked **LOCKED** run only with `RUN_LOCKED_TEST = True`.

---

## 5. Key design choices and why

| Choice | Why |
|---|---|
| **One shared module for the prior** (`synthetic_model.py`) | The training set, the spectral fits and the ML must use exactly the same prior, labels and noise model. Nothing re-implements them. A legacy regression test pins `components()` to the original notebook. |
| **Latin hypercube over a unit cube** | It covers the 28-D prior evenly with 2¹⁹ models. The same unit cube parametrises the MultiNest fits, so "the same prior" is literally the same mapping. |
| **Three SFH branches with one signature** | Dirichlet(1) bins give flexible, non-parametric SFHs. The skew-normal branch gives smooth, mostly old SFHs. The burst branch adds recent episodes. All three are interchangeable functions `(params, edges) → fractions`. |
| **16 fixed age bins (0.007–14.2 Gyr centres), truncated at the age of the Universe** | Bin labels are comparable across redshift. With this cosmology the 16th bin (14.2 Gyr) is always older than the Universe, so **only 15 bins are populated** and `*_bin_15` is padding. This is kept as in the legacy prior, by decision. |
| **Enrichment law** Z rises from Z_ini = 0.4 Z_fin to Z_fin with the mass already formed | A physical, monotonic chemical evolution with one free number. All 16 `zmet_bin_*` labels therefore carry only `Z_fin`. |
| **Colours, not magnitudes, as ML inputs** | The SFH information is in the colours. Brightness only sets the mass. r enters separately for the absolute mass and SFR. |
| **Coverage test v2 with a free magnitude offset** | The old geometric test accepted 99% of galaxies because its tolerance was too wide (median 1.19 mag), and it ignored redshift. The new test uses a redshift window, χ² over the detected bands with the offset solved analytically, and covered = p > 0.01 for n_detected − 1 degrees of freedom. The free offset is the same as a colour fit with the correct error covariance. A 0.05 mag floor stands in for model error. |
| **Spectral fits under exactly the training prior** | They answer "does the ML reproduce Bagpipes on real galaxies?" with labels the ML can in principle match. No catalogue from another code (Firefly, Starlight, pPXF) can be the label, because none reproduces the Bagpipes 16-bin model. Priors are never tuned per galaxy. |
| **Spectrum-only likelihood with a calibration polynomial** | It keeps the labels independent of the S-PLUS inputs, and it absorbs flux-calibration and aperture errors. The consequence (see §8) is that the fit reads the absorption features, not the continuum shape. |
| **Calibration mocks (33, from the full prior) instead of 6 fixed mocks** | The first gate (≥ 5 of 6 mocks) failed about 15% of the time even for a calibrated fitter. Debugging found no bug: the misses were prior-volume effects of draws from the prior tails. A coverage test on 33 prior draws measures calibration directly. |
| **Hash-based validation/test split** (md5 of `specobjid`) | Adding galaxies later never moves an existing one, so the test half stays untouched across passes. |
| **Locked test sets** (`synthetic_test.parquet`, the reference `test` half) | They are used once, for the final evaluation. No architecture, epoch or threshold choice may look at them. **Not yet used.** |
| **Report-only evaluation** | Every metric is printed, and the researcher judges. No automatic pass/fail gates decide the science. The one exception is the explicit "works" test of the A→B→C diagnostic (§8). |
| **Firefly dropped as a validation** | Firefly is another full-spectrum fitting code, so disagreeing with it cannot tell which fit is wrong. It now only provides the sample (spectrum list, redshifts, age tertiles for the selection grid). |
| **Independent rulers instead** (report-only; no galaxy is filtered) | Layer 1: model-independent spectral indices (D4000ₙ, Hδ_A, Mgb, Fe5270, Fe5335, [MgFe]′) against each fit's posterior-predictive indices. Layer 2: MPA-JHU Hα fibre SFR (BPT classes 1–2) and Gallazzi et al. (2005) light-weighted age and Z, which are built on BC03 like Bagpipes. Layer 3: scaling relations (Gallazzi mass–metallicity, Speagle 2014 main sequence, age–σ of passive galaxies). All sources are Kroupa IMF, so no IMF offset is applied. |
| **NPE with a `zuko` spline flow** | It gives full, correlated posteriors (bin-to-bin correlations matter). A mixture-density network would be less flexible, and a posterior over the 27-D unit cube is nuisance-heavy and discontinuous at the branch switches. |
| **MLP vs CNN, chosen on validation data only** | The CNN is a controlled check, not an expected improvement. The MLP won (validation loss −65.98 against −44.61). |
| **Redshift as the 19th target** | Decoding per-bin labels needs z (bin edges, bin 14 exists only below z ≈ 0.48). For photometric redshifts, each sample is decoded at its own sampled z. |
| **Heavy compute only on the PC** | User rule: fits, notebook executions, training and posterior sampling run as `hljob` jobs on the PC (16-vCPU containers, RTX 3090). The laptop runs only unit tests and quick checks. |

---

## 6. Data flow: inputs, outputs, what git ignores

| Step | Reads | Writes | In git? |
|---|---|---|---|
| Extinction correction (`extinctions.py`, run once from the repo root) | `Synthetic_domain/dr3_data.csv`, `dustMaps/csfd/` | `Synthetic_domain/dr3_corrected.csv` (adds `*_ext` columns, corrects the `*_auto` magnitudes) | yes |
| Part 1 `training_set.ipynb` | `splus-filters/`, `dr3_corrected.csv` | `training_set.parquet`, `synthetic_validation.parquet`, `synthetic_test.parquet` | **no** (`*.parquet`, reproducible from the seeds) |
| | | `splus_noise_model.csv`, `color_color.png`, `coverage_criterion.png` | yes |
| `fetch_literature.py` (run once) | SkyServer SQL, MPA-JHU DR4 pages | `literature/*.csv`, `literature/relations/*.csv`, `SOURCES.md` | yes (`literature/raw/` ignored) |
| Part 2 `spectroscopic_reference.ipynb` | `dr3_corrected.csv`, `Data/sdss_eboss_firefly-dr16.fits`, `dustMaps/`, `literature/`, `training_set.parquet`, `splus_aperture_photometry.csv` (S-PLUS Cloud only for galaxies missing from it; credentials in the repo-root `.env`) | `spectra/`, `spectral_fits/jobs/*.npz`, `selection.csv`, `failures.csv`, `calibration_mocks.joblib` | **no** (`spectra/`, `spectral_fits/`) |
| Fits on the PC (`run_spectral_fits.py` via `remote/submit_fits.sh`) | `spectral_fits/jobs/` | `spectral_fits/<job_id>.h5`, `fit_log.csv`, `job_logs/` | **no** |
| Part 2 analysis (same notebook, re-run) | the `.h5` fits | `spectroscopic_reference.csv`, `validation_summary.csv`, `validation_*.png` | yes |
| Step A export (notebook 2 §6b) | the reference set, the spectral `calibration_mocks.joblib` | `photometric_fits/jobs/*.npz`, `photometric_fits/calibration_mocks.joblib` | **not ignored yet** (see below) |
| Step A fits on the PC (`FITS_SUBDIR=photometric_fits remote/submit_fits.sh`) | `photometric_fits/jobs/` | `photometric_fits/<job_id>.h5`, `job_logs/` | **not ignored yet** |
| `make_slim_firefly.py` (before a notebook-2 PC run) | `Data/sdss_eboss_firefly-dr16.fits`, `dr3_corrected.csv` | `Data/sdss_eboss_firefly-dr16_slim.fits` | no (`*.fits`) |
| Training on the PC GPU (`train_ml_model.py` via `remote/submit_training.sh`) | `training_set.parquet`, `synthetic_validation.parquet`, `splus_noise_model.csv`, r/zml of `dr3_corrected.csv` | `ml_runs/<embedding>/config.json`, `history.csv` | yes |
| | | `ml_runs/<embedding>/model.pt` | **no** |
| Part 3 `ml_model.ipynb` (PC GPU via `remote/submit_notebook.sh`) | `ml_runs/`, the parquet sets, `spectroscopic_reference.csv`, `dr3_corrected.csv` | `ml_evaluation_summary.csv` | yes |
| | | `ml_sfh_catalogue.parquet` | **no** |
| | | `notebook_job_logs/<job-id>/` | **no** |

`.gitignore` lists `Synthetic_domain/spectral_fits/` but not `Synthetic_domain/photometric_fits/`: add it before step A runs, or a broad `git add` would pick those files up.

Also ignored ([`.gitignore`](.gitignore)): `.venv`, `*.env` (S-PLUS credentials; never print or commit them), `*.fits`, `*.pkl`, `Synthetic_domain/dr3_alldata.csv` and `.superpowers` (development ledgers, local only).

Note that `spectroscopic_reference.csv` stores `specobjid` as text. Read it with `dtype={"specobjid": str}`, because eBOSS ids exceed int64.

---

## 7. How to run things

### Local environment (the laptop: tests and light work only)

- `.venv/` holds Python 3.10 with [`requirements.txt`](requirements.txt) (Bagpipes 1.3.5, PyMultiNest, torch 2.5.1, zuko 1.2.0, …). The MultiNest shared library is built into `.venv/lib`.
- The Jupyter kernel **"Inference4ETGs (MultiNest)"** (`inference4etgs`) runs `.venv/bin/python` with `LD_LIBRARY_PATH=.../.venv/lib`, so PyMultiNest loads. Select it when you open a notebook.
- Run the tests (125 tests, CPU, about 5 minutes; without `LD_LIBRARY_PATH` the MultiNest tests are skipped):

  ```bash
  cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q
  ```

- **Do not execute notebooks, MultiNest fits, training or posterior sampling on the laptop.** These run on the PC.

### Heavy work: `hljob` jobs on the PC ([`Synthetic_domain/remote/`](Synthetic_domain/remote/))

Each job runs in a fresh container on the PC. Inside the job, `build_env.sh` rebuilds Python 3.10 with the pinned versions and MultiNest (skipped for training jobs), and checks the Bagpipes grids against `grid_md5.txt`. Only one job runs at a time. Follow a job with `hljob logs ID -f`, `hljob ls` and `hljob wait ID`.

| Task | Submit (runs locally) | Fetch the results |
|---|---|---|
| MultiNest fits of exported jobs | `Synthetic_domain/remote/submit_fits.sh` (defaults: 16 cores, 7000 MB, 24 h, `WORKERS=14`, `MAX_ITERATIONS=50000`, `FITS_SUBDIR=spectral_fits` (set `photometric_fits` for step A); `SKIP_LIST` defaults to `skipped_fit_jobs.txt` for `spectral_fits` only; finished `.h5` files travel along, so a resubmission resumes; `BUNDLE_ONLY=<folder>` packs without submitting) | `remote/fetch_fits.sh ID` → the folder the job recorded in `out/fits_subdir.txt` (works on running jobs too) |
| NPE training (GPU) | `Synthetic_domain/remote/submit_training.sh` (defaults: `EMBEDDINGS="mlp cnn"`, 8 h; fails fast without a usable GPU) | `hljob fetch ID -o /tmp/ml_job`, then copy `out/ml_runs/*` into `Synthetic_domain/ml_runs/` |
| Execute a notebook | `Synthetic_domain/remote/submit_notebook.sh <notebook> <inputs…>` (defaults: 8 cores, 16000 MB, 4 h; `GPU=1`, `NEEDS_MULTINEST=1`, `FLOAT32_FILES` to halve big parquet uploads, `BUNDLE_ONLY`; inputs written `../path` keep their repo path; `.env` is never packed; the exact notebook-2 command is in the remote README's "Notebook 2" section) | `remote/fetch_notebook.sh ID` → puts the executed notebook and every file it wrote back in place |

`python run_spectral_fits.py --status` (from `Synthetic_domain/`) counts finished fits. The details, server facts and benchmark are in [`Synthetic_domain/remote/README.md`](Synthetic_domain/remote/README.md).

### Typical cycle

1. Edit the code and run `pytest` locally.
2. Submit the notebook or the fits to the PC.
3. Fetch the results.
4. Read the saved outputs in VS Code.
5. Commit the small artifacts (CSVs, figures, configs). Large data stays ignored.

---

## 8. Current status and results

These numbers come from the committed outputs and the development ledgers, as of 2026-10-26.

### Part 1 (complete)

- 524,288 training models; the held-out sets have 8.6 detected bands on average, against 8.7 in DR3.
- Coverage test v2: 92.1% of 51,456 S-PLUS DR3 galaxies are covered (94.1% of those with zml inside the prior range). The old geometric test covered 99.1%. Synthetic validation coverage is 96.3%.
- Coverage is lowest for bright galaxies (78% at r < 17, 85% at 17–18) and above z = 0.7 (72%, outside the prior).

### Part 2 (60 galaxies analysed; +108 exported, their fits not yet fetched)

- Fits: all 60 galaxies of the first pass converged. Calibration mocks: 30 of 33 are fitted; 3 are skipped as too slow (`skipped_fit_jobs.txt`).
- Calibration verdict (ledger ruling Q9; the saved notebook outputs and the CSV's still-empty `fitter_calibrated` column predate it): **pass**. 95% coverage is Z_fin 0.900, f_old 0.933, mw_age 0.967.
  - Worst case, counting the 3 skipped mocks as misses: 0.818 / 0.848 / 0.879.
  - 68% coverage is low for Z_fin (0.50) and f_old (0.60), and 0.667 for mw_age: the intervals are too narrow at 1σ.
- Split of the 60 galaxies: **32 validation / 26 test / 2 excluded** (not covered).
- Validation layers (60 galaxies, [`validation_summary.csv`](Synthetic_domain/validation_summary.csv)):
  - **Layer 1, fraction of galaxies with |z| < 2:** D4000ₙ 0.80, Hδ_A 0.87, Fe5270 0.93, Fe5335 0.93, [MgFe]′ 0.82, Mgb 0.57.
    - Mgb is too strong in the data (median z +1.7) while Fe matches: the signature of [Mg/Fe] > 0, which solar-scaled BC03 cannot reproduce.
    - The index code agrees with MPA-JHU to ≤ 0.02 Å (median).
  - **Layer 2:** Gallazzi age offset −0.12 dex (scatter 0.12, n = 42) and metallicity offset +0.05 dex (scatter 0.17). Hα SFR offset −0.36 dex, but with n = 7 this is indicative only.
  - **Layer 3, fraction inside the 16–84% band:** mass–metallicity 0.59, main sequence 0.23 (the sample is not main-sequence selected), age–σ 0.34.
- Second pass: `GALAXIES_PER_CELL` went from 5 to 14 (168 galaxies). The 108 new jobs were exported (`spectral_fits/selection.csv` has 168 rows) and submitted as job `20261025-174838-6133cb`. Locally there are still only the 90 `.h5` files of the first pass (60 galaxies + 30 mocks). Once that job is fetched, notebook 2 must be re-run on the PC. The saved outputs of notebook 2 still show the 60-galaxy pass.

### Part 3 (validation-only run done; locked test NOT run)

- **Chosen model: MLP** (validation loss −65.98; CNN −44.61). The MLP trained for 200 epochs in 33.5 min on the GPU, with no overfitting.
- **Synthetic validation coverage (MLP):** mw_age, Z_fin, log_mstar and log_sfr are within about 0.02 of nominal (68% ≈ 0.66–0.68; 95% ≈ 0.94–0.95). f_young, f_mid and f_old are over-covered at 68% (≈ 0.80), because the fractions pile up at 0 and 1.
- **DR3 catalogue:** 49,603 of 54,150 galaxies predicted (8.4% out of domain, 23 without r), 17 galaxies/s on CUDA.
- **Agreement with the Bagpipes spectral posteriors: POOR** (validation half, 32 galaxies).

  | label | fraction of spectral medians inside ML 95% | median offset (combined σ) | in physical units (ML vs spectral) |
  |---|---|---|---|
  | mw_age | 0.22 | +9.8 | 9.2 vs 6.3 Gyr (+2.8 Gyr) |
  | Z_fin | 0.31 | −4.9 | 0.30 vs 1.63 Z☉ |
  | log_mstar | 0.28 | +1.9 | +0.24 dex |
  | f_old | 0.38 | +3.6 | — |
  | log_sfr | — | — | −1.1 dex |

- **The inputs are not the cause:**
  - the units agree;
  - DR3 total colours give the same result;
  - the galaxies lie inside the training colour domain (nearest-model χ²/dof ≈ 0.8);
  - the 30 training models nearest in colour are 9.2 Gyr old, as the ML says.

  The working interpretation: the spectral fits, with their calibration polynomial, read the absorption features, while 12-band colours constrain the continuum shape, and the two disagree.

### Open diagnostic: A → B → C (user decision, 2026-10-26)

The same "works" test applies at every step. It uses the `validation` half only; the `test` half stays locked.

- The labels are mw_age, Z_fin, f_old and log_mstar.
- A label agrees when both hold:
  - ≥ 80% of the reference medians lie inside the ML 95% interval;
  - |median offset| ≤ 1 combined σ.
- The step works when all four labels agree.

The steps run in order:

- **A — photometric Bagpipes fits (code done, not run yet).**
  - What is fitted: each validation galaxy's 3″ aperture photometry, the ML's own input, with MultiNest under the training prior at the spectroscopic z.
  - Likelihood: the same as the simulator. It is Gaussian in magnitudes with the catalogue errors and no error floor, over the detected bands only. The magnitude offset is marginalised analytically.
  - Calibration check: the 33 calibration mocks, noised with `add_splus_noise`.
  - Code: commits `76cbcab`, `b50b290`, `ca919bb`, `687b49b`, plus review fixes still uncommitted in the working tree. `fit_photometry` and helpers in `synthetic_model.py` §6 (prior `PHOTOMETRIC_DIMS` = `DIMS` without z, 28 dimensions); job kinds and `FITS_SUBDIR=photometric_fits` in the fit scripts; notebook 2 §6b (`phot_*` columns); `agreement_verdict` (`MIN_INSIDE_95 = 0.80`, `MAX_ABS_MEDIAN_OFFSET = 1.0`) in `ml_evaluation.py`, used by `ml_model.ipynb` §3b.
  - To run it (all on the PC): `python make_slim_firefly.py`; notebook 2 (exports the jobs); `FITS_SUBDIR=photometric_fits remote/submit_fits.sh`; notebook 2 again, bundling `photometric_fits/`; then `ml_model.ipynb`, which prints the A verdict.
  - If A works, the ML is a faithful Bagpipes emulator on colours, and the spectral gap is a measure of how much information 12 bands carry.
  - Plan: [`2026-10-26-photometric-bagpipes-diagnostic.md`](docs/superpowers/plans/2026-10-26-photometric-bagpipes-diagnostic.md).
- **B — joint spectrum + photometry fits**, used as the reference labels. B runs if A does not work.
- **C — a reweighted training prior.** C runs if B does not work.
  - It gives more weight to young and star-forming SFH branches and to high (super-solar) Z_fin.
  - The weights are chosen against the DR3 colour distributions and the coverage test.
  - The training set and the held-out sets are regenerated, and the MLP is retrained.

---

## 9. Known limitations and open decisions

- **mw_age is capped at ≈ 10.2 Gyr.** The oldest populated bin runs from 7.05 Gyr to the age of the Universe, and `*_bin_15` is always empty. Changing `AGE_BIN_CENTRES` would mean regenerating everything.
- **DR3 vs DR4 photometry.** The noise model and the coverage reference use DR3. The reference-set apertures are DR4 (`dr4_dual`), whose g band is about 0.07–0.1 mag fainter than DR3 for the same objects. This is an unresolved decision: move everything to DR4, or fetch DR3 apertures.
- **The spectral fits are overconfident at 1σ.** 68% coverage of the calibration mocks is ≈ 0.5–0.67.
- **[Mg/Fe] > 0** in passive galaxies cannot be modelled by solar-scaled BC03 (the Mgb excess).
- **BOSS galaxies** have 2″ fibres against 3″ apertures. They are reported separately and have no MPA-JHU or Gallazzi values.
- **The reference sets are pilot-sized**: 32 / 26 galaxies until the +108 fits are analysed.
- **Notebook outputs lag their source.** Notebook 2's saved outputs show the 60-galaxy pass without §6b; `ml_model.ipynb`'s predate section 3b and the single-sampling-pass change in section 3. Both need a PC re-run.
- **Notebook 2 on the PC** needs network access (SDSS queries) and a large upload (the CSFD map alone is 289 MB); see the remote README.

---

## 10. Where the history and the rationale live

- **Specs** ([`docs/superpowers/specs/`](docs/superpowers/specs/)) record what was decided and why. Read the dated amendments: external fitting and the calibration test; Firefly as a check, then dropped; the validation rulers; redshift as the 19th target; the A→B→C diagnostic.
- **Plans** ([`docs/superpowers/plans/`](docs/superpowers/plans/)) give task-by-task implementation with tests:
  - `2026-10-23-training-and-reference-sets.md` (parts 1–2);
  - `2026-10-24-ml-sfh-posterior.md` (part 3);
  - `2026-10-25-bagpipes-validation-rulers.md` (validation layers, Firefly removal);
  - `2026-10-26-photometric-bagpipes-diagnostic.md` (step A).
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
