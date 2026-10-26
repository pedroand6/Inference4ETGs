# Delayed-τ SFH with a single metallicity — design

**Date:** 2026-10-26
**Supersedes:**
- the 16-bin non-parametric SFH and the enrichment law of `2026-10-23-training-and-reference-sets-design.md` (parts 1–2);
- the 16-bin targets of `2026-10-23-ml-sfh-posterior-design.md` (part 3);
- the A → B → C diagnostic chain (ML spec amendment of 2026-10-26).

Every other rule of those specs still holds: code style, the locked data, report-only evaluation, and the independent validation rulers.

## Goal

The SFH model is Bagpipes' parametric delayed-τ SFH with one metallicity, the kind of model Bagpipes was built around. The training prior covers the whole delayed-τ parameter space. Everything downstream is rebuilt under this model:
- the training and synthetic validation/test sets;
- the Bagpipes fits of the reference set;
- the calibration tests;
- the ML posterior.

## Decisions (user, 2026-10-26)

1. **SFH and metallicity:** a single Bagpipes `delayed` component, SFR(t) ∝ (t − t₀) e^{−(t − t₀)/τ}, with one metallicity Z. There is no enrichment law and there are no age bins.
2. **Prior:** a Latin hypercube, log-uniform and wide:

| Parameter | Range | Scale |
|---|---|---|
| `redshift` | 0.002 – 0.7 | uniform |
| `age` (time since star formation began) | 0.05 Gyr – age of the Universe at the draw's z (minus `UNIVERSE_AGE_MARGIN` = 0.1 Gyr) | log-uniform; the upper bound depends on z |
| `tau` | 0.01 – 20 Gyr | log-uniform |
| `metallicity` Z | 0.005 – 2.5 Z☉ (the BC03 grid range) | log-uniform |
| `av_ism` | 0.01 – 2.0 mag | uniform |
| `tau_bc` (birth-cloud optical depth) | 0.5 – 4 | log-uniform |
| `logU` | −4 to −1 | uniform |
| `log_mass_formed` | fixed at 10 for the models; rescaled to the observed r magnitude as now (`scale_to_observed_r`, `rescale_mass_labels`) | — |

3. **Kept unchanged:**
   - dust: Charlot & Fall (2000), exponent 0.7, birth-cloud age 0.01 Gyr;
   - Bagpipes nebular emission;
   - the bc03_miles grids, the Kroupa IMF and the cosmology (H0 = 70, Ωm = 0.3);
   - the S-PLUS noise model;
   - the coverage test;
   - the spectral fit machinery: calibration polynomial, veldisp, noise_scale, MultiNest settings, `safe_*_log_likelihood`;
   - the photometric fit (offset marginalised analytically);
   - job kinds and `FITS_SUBDIR`.
4. **The 16-bin model is replaced, not kept alongside.** The commit before any change is tagged `nonparametric-16bin`, so every earlier result can be recovered.
5. **ML targets:** physical parameters plus derived labels. Only the MLP embedding is trained, because the CNN lost clearly in the 16-bin comparison.

## Labels

**Drawn parameters:**
- `redshift`, `age` [Gyr], `tau` [Gyr], `metallicity` [Z☉];
- `av_ism` [mag], `tau_bc`, `logU`;
- `log_mass_formed` [log M☉].

**Derived for each model or posterior sample:**
- **From Bagpipes for the model:**
  - `log_mstar` [log M☉], the surviving stellar mass;
  - `log_sfr` [log M☉ yr⁻¹], averaged over the last 100 Myr (as now), with the existing floor;
  - `log_ssfr`.
- **Analytic, from the delayed-τ form:**
  - `mw_age` [Gyr ago], the mass-weighted age;
  - `f_young`, `f_mid`, `f_old`, the formed-mass fractions younger than 0.1 Gyr, between 0.1 and 1 Gyr, and older than 1 Gyr (`LABEL_EDGES`).

  A unit test checks these against Bagpipes' own `sfh` arrays.

**Removed:** `frac_bin_*`, `logsfr_bin_*`, `zmet_bin_*`, `Z_fin`, `sfh_type` and the SFH branches (Dirichlet, skew-normal, burst).

**`TRAINING_LABELS` order (used by the parquet sets and the fits):** the drawn parameters, then the derived labels.

## Part 1: training set (`training_set.ipynb`, `synthetic_model.py` sections 1–2)

**Prior functions in `synthetic_model.py`:**
- `DIMS` becomes `["z", "log_age_fraction", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU"]`.
- The age draw is a unit-interval fraction of the log range [log 0.05, log t_max(z)], so a single unit cube works at every redshift.
- `components(u, fixed_redshift=None)` keeps its signature and returns the Bagpipes components dict with `delayed: {age, tau, massformed, metallicity}`.
- `summary_labels` is rewritten for the new labels.
- The functions for bin edges and bin labels are removed.

**Sets:** the training set, `synthetic_validation.parquet` and `synthetic_test.parquet` are regenerated with the same sizes and seeds as now. All of them are made by a notebook job on the PC (`training_set.ipynb` through `remote/submit_notebook.sh`).

## Part 2: reference set (`spectroscopic_reference.ipynb`)

**Kept:**
- the selection: 168 galaxies, 14 per redshift × age cell;
- `selection.csv`, the SDSS spectra and the 3″ aperture photometry;
- the hash-based split (validation 87, test 78, excluded 3), which no galaxy can leave.

**Moved to `spectral_fits/superseded/` and `photometric_fits/superseded/`:**
- the old fits;
- the old calibration mocks.

`skipped_fit_jobs.txt` is emptied: its entries name old mocks, and the new mocks reuse the `calib_<k>` ids.

The notebook already supersedes fits made under an outdated setup. The test for that now compares a stored `sfh_model = "delayed"` key in the fit metadata, which every new fit writes.

**New fits under the new prior (spectral fits use `SPECTRAL_DIMS` = DIMS without z, plus veldisp and noise_scale):**
- spectral fits of all 168 galaxies;
- 33 new calibration mocks;
- photometric fits of the 87 validation galaxies;
- 33 photometric calibration mocks.

**Calibration test:**
- labels are `age`, `tau`, `metallicity` and `log_mstar`;
- it passes when 95% coverage is at least 0.85 for each label;
- 68% coverage is reported too.

**Validation rulers:** layers 1–3 are unchanged. The Gallazzi and scaling-relation comparisons use the single `metallicity` where they used `Z_fin`. `fetch_literature.py` is re-run so MPA-JHU and Gallazzi cover all 168 galaxies, not only the original 60.

## Part 3: ML (`ml_data.py`, `ml_networks.py`, `ml_evaluation.py`, `train_ml_model.py`, `ml_model.ipynb`)

**Flow targets θ (10), in this order:**
- log age, log tau, log metallicity;
- av_ism, log tau_bc, logU;
- log_mass_formed, log_mstar, log_sfr;
- redshift.

`N_TARGETS` follows `TARGET_NAMES`.

**Decoding:** each posterior sample is decoded to the derived labels with the same analytic functions as part 1, at the sample's own redshift.

**Inputs:** unchanged (37 values: 11 colours, 11 colour errors, 11 flags, m_r, σ_r, z_obs, σ_z).

**Training:** MLP only, on the GPU, with the same flow, optimiser, early stopping and per-epoch noise simulation.

**Evaluation (report-only):**
- synthetic calibration and accuracy for every target and derived label;
- the "works" test (≥ 80% of reference medians inside the ML 95% interval and |median offset| ≤ 1 combined σ, validation half only), against both:
  - photometric Bagpipes, i.e. does the ML reproduce Bagpipes on the colours;
  - spectral Bagpipes, i.e. does it agree with the spectra;
  
  on `age`, `tau`, `metallicity` and `log_mstar`, with `mw_age` reported;
- the DR3 catalogue (median and 16/84th percentiles per label, plus the `out_of_domain` and `no_r_band` flags).

The locked test still runs once, after the user inspects the validation results.

## Compute (user rule, 2026-10-25)

The laptop only runs code changes, unit tests and quick checks. Everything else is an `hljob` job on the PC:
- set generation;
- notebook executions;
- fits;
- training;
- posterior sampling.

Order (one job at a time, so the jobs queue):
1. Training and synthetic sets: CPU notebook job.
2. Notebook 2: exports the fit jobs and moves the old fits to `superseded/`.
3. Spectral fits and mocks: CPU, 14 workers.
4. Literature refetch: a light job, or local if it is only network queries without credentials.
5. MLP training: GPU.
6. Photometric fits: CPU.
7. Notebook 2 run that reads all the fits.
8. `ml_model.ipynb`: GPU.

## Testing (CPU, fast)

- **Prior map:**
  - every parameter within bounds;
  - log-uniform spacing, i.e. a uniform cube gives uniform log values;
  - age never exceeds t_max(z);
  - the Latin hypercube is stratified.
- **Derived labels:**
  - analytic mw_age and f_young/f_mid/f_old match a numerical integral of Bagpipes' `sfh` for several (age, tau) values, including tau ≪ age and tau ≫ age, to 1%;
  - the mass fractions sum to 1.
- **Targets:** the encode/decode round trip is exact; decoded labels match `summary_labels` for simulated models.
- **Fits and jobs:** the existing tests are adapted to the new DIMS. A smoke fit runs on a noiseless mock. A fit without `sfh_model = "delayed"` is superseded.
- **Notebooks:** valid JSON and every code cell compiles. They are executed only on the PC.

## Out of scope

- Other parametric SFHs (exponential, log-normal, double power law) and non-parametric SFHs.
- Metallicity evolution.
- Changing the dust or nebular models.
- Hyper-parameter searches.
- The CNN.

## Amendment (2026-10-26, from the code inventory)

1. **`log_sfr` is derived, not a flow target.**
   - About 20% of delayed-τ prior draws have no star formation in the last 100 Myr. Their `log_sfr` sits exactly on the floor (`SFR_FLOOR × M*`), and a flow cannot learn a target with that point mass.
   - `log_sfr` is therefore computed analytically for training labels and posterior samples alike: SFR₁₀₀ = 10^log_mass_formed × F(<0.1 Gyr) / 10⁸ yr, with the floor.
   - The flow has 9 targets: log age, log tau, log metallicity, av_ism, log tau_bc, logU, log_mass_formed, log_mstar, redshift (redshift last).
2. **Age, τ and Z are compared in log.**
   - The calibration test and the agreement test use `log_age`, `log_tau`, `log_metallicity` and `log_mstar`, since the prior spans several decades.
   - These log columns are written next to the linear ones.
3. **Checkpoints record their `target_names`.** Loading a run whose names differ from `ml_data.TARGET_NAMES` raises an error. The old runs move to `ml_runs/superseded/`.
