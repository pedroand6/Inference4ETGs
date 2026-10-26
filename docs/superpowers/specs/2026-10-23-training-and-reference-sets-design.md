# Training set v2 and spectroscopic reference sets — design

> **Superseded in part (2026-10-26):** the 16-bin SFH, its labels, the enrichment law and the SFH branches (Part 1 "Labels" and "Prior", and the parts that use them) are replaced by `2026-10-26-delayed-tau-sfh-design.md`. Everything else here still holds. The 16-bin code is the git tag `nonparametric-16bin`.

**Date:** 2026-10-23
**Scope:** parts 1 and 2 of the Inference4ETGs ML pipeline. Part 3 (the ML model) gets its own spec after this one is done.

## Goal

Predict, from S-PLUS 12-band photometry and redshift, the star formation history of galaxies of **all types** in 16 age bins: mass fraction / SFR per bin and metallicity per bin, as Bagpipes with BC03-MILES (`bc03_miles`, the 2016 update of BC03) would infer them.

This spec builds the data the ML model needs:

1. **Training set v2** — the synthetic Bagpipes training set with 16-bin labels, a realistic S-PLUS noise model, synthetic validation/test sets and a stricter coverage test.
2. **Spectroscopic reference set** — real S-PLUS galaxies whose labels come from Bagpipes fits to their SDSS spectra, under exactly the training prior, checked against the Firefly catalogue and split into a validation set and a locked test set.

## Why this design

Two different questions validate the ML model, and they need different data:

- *Does it reproduce Bagpipes?* — answered by labels from the same SPS, SFH parametrisation and prior as the training set.
- *Is it physically right?* — answered by an independent literature catalogue (Firefly, SDSS DR16).

Fitting SDSS spectra with Bagpipes under the training prior answers the first question on real galaxies: the label is a Bagpipes posterior constrained by a spectrum, which carries far more information than 12 bands. Firefly then checks those labels against an independent code. No catalogue from another code (Starlight, Firefly, pPXF) can serve as the label itself, because none reproduces the Bagpipes 16-bin model.

Priors are **not** tuned per galaxy from the literature: that would make the Firefly check circular and would put the test labels under a different prior from the training set.

## Global constraints

- SPS: Bagpipes `bc03_miles` grids (the installed default). Cosmology: `FlatLambdaCDM(H0=70, Om0=0.3)`.
- Redshift range of the prior: `Z_RANGE = (0.002, 0.7)`.
- Training set size: `N_MODELS = 2**19`. Synthetic validation and test sets: `2**14` models each.
- Sampler for spectral fits: MultiNest through PyMultiNest. It is built in `.venv/lib`; the notebook kernel must have `LD_LIBRARY_PATH` pointing there.
- Spectral fitting budget: **at most 3 hours** wall time on 11 worker processes (the machine has 12 cores).
- SDSS spectra are corrected for Milky Way extinction with the same map and law as the S-PLUS photometry (`extinctions.py`): CSFD E(B−V), CCM89, R_V = 3.1.
- Spectral fits use the spectrum only. S-PLUS photometry is never part of the spectral likelihood, so reference labels stay independent of the ML inputs.
- Large downloaded or generated data (`spectra/`, `spectral_fits/`, `*.parquet`) is git-ignored; generated sets are reproducible from their seeds. The old tracked `training_set.csv` is removed from the repository (it stays in git history).

## Code style

These rules apply to every file this spec creates or rewrites.

1. English names that say what they hold (`model_colours`, `galaxy_errors`). No single letters except loop indices and standard maths (`z` for redshift).
2. One step per block, with numbered section banners (`# 1. Load the S-PLUS catalogue`).
3. Comments explain *why* a step exists or what physics it encodes, not what the Python does.
4. No compact code: no nested comprehensions, no chained one-liners, no lambdas, no clever indexing without a named intermediate. Prefer explicit `for` loops when they read better.
5. Every function has a docstring: purpose, inputs with units, return value.
6. Constants are UPPER_CASE at the top of the file, with units in a comment. No magic numbers inline.
7. Each notebook opens with a markdown cell giving its goal, inputs, outputs and run order. Each section starts with a short markdown cell explaining the method.
8. Notebook cells orchestrate and display. Logic shared by both notebooks, or logic that needs unit tests, lives in the module.
9. Parallel cases follow one pattern: e.g. the three SFH branches are three functions with the same signature, not an if-chain with inline maths.
10. Be concise: short comments, nothing that restates the code, markdown cells of a few sentences.

## Files

All in `Synthetic_domain/`:

| File | Role |
|---|---|
| `synthetic_model.py` (new) | Shared module: prior, SFH branches, enrichment law, S-PLUS filters and magnitudes, 16-bin labels, S-PLUS noise model, coverage test, spectral likelihood pieces. No plotting, no data paths. |
| `tests/test_synthetic_model.py` (new) | `pytest` tests for the module. |
| `training_set.ipynb` (rewritten) | Part 1. |
| `spectroscopic_reference.ipynb` (new) | Part 2. |
| `ml_model.ipynb` | Part 3. Name reserved; not created by this spec. |

### Module interface

Names below are the contract between the module, the notebooks and the tests.

- Constants: `AGE_BIN_CENTRES` (Gyr, 16 values, unchanged), `N_BINS = 16`, `Z_RANGE`, the prior ranges currently in the notebook, `DIMS` (unit-cube dimension names), `BANDS`, `FILTER_FILES`, `PIVOT_WAVELENGTHS`.
- `bin_edges(redshift) -> np.ndarray` — bin edges in Gyr ago, truncated 0.1 Gyr before the age of the Universe.
- SFH branches, same signature `(params: dict, edges: np.ndarray) -> np.ndarray` of mass fractions:
  `dirichlet_fractions`, `skew_normal_fractions`, `skew_normal_burst_fractions`.
- `choose_sfh_branch(sfh_type_draw: float) -> int` — 0 Dirichlet, 1 skew-normal, 2 skew-normal + burst.
- `bin_metallicities(fractions, z_final) -> np.ndarray` — Z/Z☉ per bin from the enrichment law.
- `components(u, fixed_redshift=None) -> (dict, dict)` — Bagpipes components and drawn parameters for one unit-cube vector. With `fixed_redshift`, the `z` coordinate of `u` is ignored.
- `bin_labels(fractions, metallicities, edges, log_mass_formed) -> dict` — `frac_bin_00…15`, `logsfr_bin_00…15`, `zmet_bin_00…15`, padded to 16.
- `summary_labels(bagpipes_sfh) -> dict` — `log_mstar`, `log_sfr`, `mw_age`, `f_young`, `f_mid`, `f_old`.
- `ab_magnitudes(flux_per_angstrom) -> np.ndarray`.
- `fit_noise_model(catalogue) -> pd.DataFrame` and `add_splus_noise(magnitudes, noise_table, rng) -> (noisy_magnitudes, errors, detected)`.
- `coverage_test(observed_magnitudes, observed_errors, detected, observed_z, observed_dz, model_magnitudes, model_z) -> pd.DataFrame`.
- `calibration_polynomial(model_flux, data_flux, data_ivar, order=2) -> np.ndarray` and `spectral_log_likelihood(...) -> float`.

## Part 1 — Training set v2 (`training_set.ipynb`)

### Labels

For each model, in addition to the drawn parameters (`z`, `sfh_type`, `Z_fin`, `Av_ism`, `Av_bc`, `logU`) and the summary labels (`log_mstar`, `log_sfr`, `mw_age`, `f_young`, `f_mid`, `f_old`):

- `frac_bin_00 … frac_bin_15`: fraction of the total mass formed in each age bin; sums to 1.
- `logsfr_bin_00 … logsfr_bin_15`: log of the mean SFR in each bin, for 10¹⁰ M☉ formed. Scales with the photometric offset Δm like `log_sfr`.
- `zmet_bin_00 … zmet_bin_15`: Z/Z☉ of each bin.

Bins older than the Universe at that redshift get `frac = 0`, `logsfr = NaN`, `zmet = NaN`.

The notebook states that `zmet_bin_*` is fully set by `Z_fin` and the SFH (Z_ini = 0.4 Z_fin, rising with accumulated mass), so the 16 metallicities carry one free number.

### Prior

Unchanged except `Z_RANGE = (0.002, 0.7)` and `N_MODELS = 2**19` (~22 min). The notebook markdown is updated to match the code.

### S-PLUS noise model

Fitted to `dr3_corrected.csv` (PROB_GAL > 0.8):

- per band, in 0.25 mag bins of measured magnitude: the empirical distribution of log errors, stored as 21 quantiles. A band counts as detected when its error is below 1 mag (the coverage test's rule), so the detection probability of a bin is the fraction of its errors below 1 mag;
- `add_splus_noise` first scales each model to the `r_auto` of a random DR3 galaxy with `|zml − z| < 0.05`; then, per band, draws an error from that bin's quantiles, marks the band detected if the error is below 1 mag, and perturbs the flux with a Gaussian of that error. Models fainter than the faintest populated bin are undetected.

### Synthetic validation and test sets

The training set is saved as `training_set.parquet`. Two independent Latin Hypercube draws (own seeds), `2**14` models each, with noise and detection masks applied: `synthetic_validation.parquet`, `synthetic_test.parquet`. The training set itself stays noise-free; noise augmentation during training belongs to the ML spec.

### Coverage test v2

Replaces the geometric test, which accepted 99% of galaxies mainly because its tolerance (quadrature sum of 11 colour errors, median 1.19 mag) was too wide, and which ignored redshift.

For each galaxy:

1. Keep only models with `|z_model − zml| < dz`, where `dz` follows the existing `obs_dz` rule (`Z_SIGMA`, `DZ_MIN`, `DZ_MAX`, `SIGMA_Z`).
2. For each kept model, χ² over detected bands with the magnitude offset solved analytically (weighted mean of `m_obs − m_model`). Non-detected bands have zero weight. Errors include the existing 0.05 mag floor, which stands in for model error.
3. The galaxy is **covered** if the minimum χ² has p-value > `P_MIN = 0.01` for `n_detected − 1` degrees of freedom.

Output columns: `chi2_min`, `n_dof`, `p_value`, `covered`, `best_model`. The notebook prints one comparison with the old test, then drops it. The colour-colour and criterion figures are kept, adapted to the new columns.

## Part 2 — Spectroscopic reference set (`spectroscopic_reference.ipynb`)

### 1. Sample

- Crossmatch S-PLUS DR3 (PROB_GAL > 0.8) with the Firefly DR16 catalogue (`Data/sdss_eboss_firefly-dr16.fits`) within 1.5″. Keep `CLASS_NOQSO == "GALAXY"` with `Z_NOQSO` inside `Z_RANGE`. About 1,400 galaxies today.
- Record `fiber_diameter`: 3″ for SDSS (`RUN2D == "26"`), 2″ for BOSS (`RUN2D == "v5_13_0"`).
- Select `N_REAL_GALAXIES = 60`: a grid of 4 redshift bins × 3 Firefly mass-weighted-age bins (quantiles of the matched sample), 5 galaxies per cell, highest median spectral S/N first. Cells with fewer than 5 galaxies give what they have.

### 2. Aperture photometry

Query the S-PLUS Cloud with `splusdata` for the 3″-diameter aperture magnitudes and errors of the selected galaxies. The first implementation step confirms the exact column names. Correct them for Galactic extinction with `extinctions.py`'s method. BOSS galaxies keep this aperture and are reported separately.

### 3. Spectra

- Download with `astroquery.sdss` by plate / MJD / fibre into `spectra/`.
- Correct for Milky Way extinction (CSFD E(B−V) at the galaxy position, CCM89, R_V = 3.1): flux × 10^(0.4 A_λ), inverse variance × 10^(−0.8 A_λ).
- Mask pixels with `ivar = 0`, ±10 Å around the 5577 Å and 6300 Å sky lines, and keep only 3800–9200 Å observed.
- Rebin to a regular 5 Å grid by inverse-variance-weighted mean of the unmasked pixels in each bin; empty bins get zero inverse variance (~30% faster model calls at a resolution that still carries the SFH information). `spectres` is not used because masked gaps distort its bin widths.

### 4. Spectral fit

- **Likelihood:** Bagpipes model spectrum from `synthetic_model.components(u, fixed_redshift=Z_NOQSO)` plus `veldisp`; multiplied by the least-squares order-2 calibration polynomial (its constant term absorbs the normalisation); Gaussian likelihood with errors multiplied by the noise-scaling factor.
- **Prior:** the training unit cube without `z`, plus `veldisp` uniform 40–400 km/s and `noise_scale` log-uniform 1–5.
- **Sampler:** PyMultiNest with `n_live_points = 300`, `sampling_efficiency = 0.8`, `multimodal = True`, and an iteration cap sized to ~25 min per fit. Fits that reach the cap are flagged `converged = False`.
- **Output:** per galaxy `spectral_fits/<specobjid>/` (MultiNest files) and `spectral_fits/<specobjid>.h5` (equal-weight posterior samples of the unit cube and of all labels, log-evidence, `converged`, run time). Existing outputs are skipped, so the notebook can be re-run with a larger `N_REAL_GALAXIES`.
- **Labels:** after sampling, each equal-weight posterior sample is re-evaluated once with a Bagpipes model that also computes S-PLUS photometry, giving `bin_labels`, `summary_labels` and the model `r` magnitude. Masses and SFRs are normalised with the aperture `r` magnitude offset, as in the training set. Median and 16th/84th percentiles go to `spectroscopic_reference.csv`.
- **Run plan within the 3 h budget:** 6 inject-and-recover mocks + 60 real galaxies = 66 fits, run as 11 parallel workers (six waves of ~25 min). The notebook first times one fit, prints the projected total, and lowers `N_REAL_GALAXIES` (saying so) if the projection exceeds 3 h.

### 5. Inject-and-recover

Before trusting real fits: 6 mock spectra from known `u` vectors (2 per SFH branch), at redshifts and S/N drawn from the selected sample, with SDSS-like noise. Pass: in at least 5 of 6 mocks the true `mw_age`, `Z_fin` and `f_old` lie inside the 95% interval. A failure stops the notebook with a clear message before the real fits.

### 6. Firefly check

Per galaxy, using `Chabrier_MILES_*_massW` and `Chabrier_MILES_total_mass` (−9999 → NaN):

- `delta_log_age`: Bagpipes `mw_age` vs Firefly mass-weighted age (years in Firefly).
- `delta_log_metallicity`: Bagpipes mass-weighted Z/Z☉ vs Firefly mass-weighted Z/Z☉ (units confirmed in the first implementation step).
- `delta_log_mass`: reported, not used for the flag (the spectral fit carries no normalisation; Bagpipes mass comes from aperture photometry).
- `firefly_consistent = |delta_log_age| < MAX_DELTA_LOG_AGE and |delta_log_metallicity| < MAX_DELTA_LOG_METALLICITY`, both 0.3 dex. NaN when Firefly values are missing.

The notebook reports the median offset and scatter for SDSS vs BOSS and for star-forming vs passive (Bagpipes `log_sfr − log_mstar` above/below −11). Failing galaxies are flagged, never deleted.

### 7. Coverage

`coverage_test` from the module, applied to the aperture photometry, at `Z_NOQSO` with `dz = DZ_MIN`.

### 8. Locked split

Galaxies with `converged`, `covered` and `firefly_consistent` all true are split 50/50 into `validation` and `test` by a hash of `specobjid`, so adding galaxies later never moves an existing one. The split is saved as a column. The notebook states that the test half is used once, for the final ML evaluation, and that ~30/30 is a pilot-size set.

## Testing and error handling

`pytest` tests (pytest added to the venv):

- **Refactor regression:** for fixed `u` vectors, `components` gives the same Bagpipes components as the current notebook code run with the new `Z_RANGE`.
- Bin fractions sum to 1; bins older than the Universe are masked; metallicity rises from the oldest bin to the youngest.
- `add_splus_noise` reproduces the DR3 detection rate and median error per band within 5% when applied to models scaled to DR3 magnitudes.
- `coverage_test` marks a noisy model as covered, and the same model with one colour shifted by 0.5 mag as uncovered.
- `calibration_polynomial` exactly recovers a known polynomial multiplied into a model spectrum.

Notebook check cells (`assert` with a clear message) at each hand-off: file sizes, no NaN magnitudes in the training set, every fitted galaxy present in the final CSV.

Failures: a failed download or fit is written to `spectral_fits/failures.csv` and skipped; it never stops the batch.

## Out of scope

- The ML model, its architecture and its training-time noise augmentation (next spec).
- Fornax3D / MUSE and any other external sample.
- Modelling the BOSS 2″ fibre vs 3″ aperture mismatch beyond flagging and reporting it separately.
- Running more than 60 real spectral fits in this pass.

## Amendment (2026-10-23): external fitting and a calibration gate

Approved by the user after the first full run stopped at the inject-and-recover gate.

**Why.**
- *No defect behind the stop.* Debugging found no fault in the model, likelihood or sampler. χ² at the true parameters matches the noise, and MultiNest returns the same posterior with sampling efficiency 0.3 or with 800 live points.
- *The misses are prior-volume effects.* Truths drawn from prior tails land where the data mildly prefer the region the prior favours.
- *The gate is too strict for a sound fitter.* "≥ 5 of 6 mocks" fails about 15% of the time even when the fitter is calibrated.
- *The time budget is too short.* Real galaxies need well over 8,500 iterations (more than 85 min each on this laptop under load), so the 3-hour budget cannot hold the sample.

**Changes (they supersede §4 "Run plan" and §5 "Inject-and-recover" of Part 2):**
- **Calibration test.**
  - `N_CALIBRATION_MOCKS = 33` mocks drawn from the full training prior, with `sfh_type` drawn too. Each uses the wavelengths, mask and errors of a random selected galaxy.
  - The fitter is calibrated when, for each of `mw_age`, `Z_fin` and `f_old`, the fraction of truths inside the 95% interval is at least `MIN_COVERAGE_95 = 0.85`.
  - The 68%-interval coverage is reported too.
  - The result is written into the reference CSV as `fitter_calibrated`; it does not stop the notebook.
- **Fitting outside the notebook.**
  - The notebook exports a portable job bundle: one `.npz` per job in `spectral_fits/jobs/`.
  - `Synthetic_domain/run_spectral_fits.py` runs one MultiNest fit per core as independent processes, with BLAS threads pinned to 1. It runs locally with `--workers N`, or as a Slurm job array with `--task-index` (example `slurm_fit.sh`).
  - It is resumable and writes the same `<job_id>.h5` files, plus a `fit_log.csv`.
  - The iteration cap is `MAX_ITERATIONS = 50000`, and there is no wall-clock budget.
- **Analysis.**
  - The notebook's analysis cells read whatever results exist (none yet is fine), report calibration coverage, and build the reference table and split as before.
  - Fits stored with a smaller iteration cap (the first run's) are moved to `spectral_fits/superseded/`, as are the first run's 6 fixed-branch mocks.

## Amendment (2026-10-25): the Firefly comparison is a check, not a filter

Decided by the user. The locked split takes every galaxy whose fit converged and that passes the coverage test. `firefly_consistent` stays in the CSV and the notebook reports the share of consistent fits as a sanity check of the spectral fits against the literature, but it no longer excludes galaxies from the validation and test sets. This supersedes the Firefly condition of Part 2 §8. The split stays a hash of `specobjid`, so the galaxies already assigned keep their half; the test half had not yet been used.

## Amendment (2026-10-25): validating the Bagpipes spectral fits against independent rulers; Firefly dropped as a validation

Decided by the user. Firefly is another full-spectrum fitting code, so disagreeing with it cannot tell which fit is wrong; it no longer validates anything. Firefly stays only as the catalogue the sample is drawn from (the SDSS spectrum list, redshifts and the age tertiles of the selection grid). Its comparison columns, figure and consistency flags are removed from notebook 2 and from `spectroscopic_reference.csv`. They neither filter nor report.

No catalogue gives the true stellar age or metallicity of an unresolved galaxy. The validation therefore uses three rulers that are independent of full-spectrum SED fitting, in order of strength.

### Layer 1 — model-independent spectral indices (primary; every fitted galaxy)
- **Indices and passbands** (air wavelengths, rest frame):
  - D4000ₙ (Balogh et al. 1999): the ratio of mean F_ν in 4000–4100 Å to 3850–3950 Å.
  - Hδ_A (Worthey & Ottaviani 1997): feature 4083.500–4122.250, blue 4041.600–4079.750, red 4128.500–4161.000.
  - Lick Mgb: 5160.125–5192.625 / 5142.625–5161.375 / 5191.375–5206.375.
  - Lick Fe5270: 5245.650–5285.650 / 5233.150–5248.150 / 5285.650–5318.150.
  - Lick Fe5335: 5312.125–5352.125 / 5304.625–5315.875 / 5353.375–5363.375.
  - The combination [MgFe]′ = √(Mgb · (0.72 Fe5270 + 0.28 Fe5335)) (Thomas et al. 2003).
  - Equivalent-width indices use the standard pseudo-continuum: a straight line between the mean fluxes of the blue and red bands.
- **Observed values:** measured on each galaxy's SDSS spectrum at native resolution, in air wavelengths, corrected for Milky Way extinction and shifted to the rest frame. Uncertainties are propagated from the inverse variance by Monte Carlo (200 noise draws). Pixels with ivar = 0 inside a band invalidate that index for that galaxy (NaN).
- **Model values (posterior predictive):** for 100 equal-weight posterior samples per galaxy, the Bagpipes model spectrum is computed on the galaxy's native observed wavelength grid. That means the same components, including `veldisp` and the nebular emission (which fills in Hδ in star-forming galaxies), and the calibration polynomial fitted to that galaxy's data. The same index code measures each sample.
- **Statistics per galaxy and index:**
  - z = (observed − model median) / √(σ_obs² + σ_model²);
  - the posterior-predictive p-value (the fraction of model samples beyond the observed value, two-sided).
- **Report:**
  - the fraction of galaxies with |z| < 2 for each index;
  - the median z (bias) and the z scatter, split by star-forming / passive and by SDSS / BOSS;
  - plots of observed against model for D4000ₙ–Hδ_A and for [MgFe]′.
- **Cross-check of the index code:** our measured indices against MPA-JHU `galSpecIndx` (d4000_n, lick_hd_a, lick_mgb, lick_fe5270, lick_fe5335) for SDSS-I/II spectra. The median difference and scatter are reported; a systematic difference larger than the MPA-JHU errors is a bug to investigate, not a result.

### Layer 2 — independent physical tracers (SDSS-I/II galaxies; BOSS galaxies are marked "not available")
- **SFR from Hα:** MPA-JHU `galSpecExtra` fibre SFR (`sfr_fib_p50`, Brinchmann et al. 2004; Balmer-decrement dust correction). It is compared only for galaxies whose MPA-JHU BPT class is star-forming. It traces the last ~10 Myr, so it is compared with Bagpipes' `logsfr_bin_00` (0–14 Myr) normalised to the 3″ aperture, with the Kroupa → Chabrier IMF offset stated in the notebook. Reported: the median offset and scatter.
- **Stellar ages and metallicities from absorption-line indices:** Gallazzi et al. (2005), r-band light-weighted log age and log Z/Z☉, built on BC03 like Bagpipes, so differences reflect the method and not the stellar models. Bagpipes light-weighted r-band age and metallicity are computed per posterior sample from per-bin model luminosities (each age bin's r-band luminosity from a Bagpipes model containing only that bin). Reported: offset, scatter, and the fraction inside the combined 1σ and 2σ.
- **Access:** MPA-JHU tables come from SkyServer SQL (astroquery, the way the S/N was queried). The Gallazzi catalogue comes from its public release; the first implementation step locates it (MPA-JHU DR4 release page or VizieR) and records the source. If it cannot be found, that sub-check is reported as unavailable, not replaced by an SED-fitting catalogue.

### Layer 3 — consolidated scaling relations (population level)
The fitted galaxies, overlaid on relations from large samples measured with independent methods:
- **Stellar mass–metallicity relation:** median and 16–84th percentiles per stellar-mass bin, built from the full Gallazzi et al. (2005) catalogue with MPA-JHU masses, against our light-weighted metallicities and aperture-normalised masses (the aperture is stated).
- **Star-forming main sequence:** Speagle et al. (2014), log SFR = (0.84 − 0.026 t) log M* − (6.51 − 0.11 t), with t the age of the Universe in Gyr at each galaxy's redshift, and its ±0.2 dex scatter, for the star-forming galaxies.
- **Age–velocity dispersion of passive galaxies:** median and percentiles of Gallazzi age against SDSS `velDisp` from the large sample, against our passive galaxies with their SDSS `velDisp`.
- Reported as figures with the fraction of our galaxies inside each relation's 16–84th band.

### Where it runs and what it outputs
A new tested module section (index definitions and measurement; light-weighted age and metallicity; Speagle relation) plus a new section of notebook 2 that computes all three layers for every fitted galaxy, including the coming +100. Outputs:
- per-galaxy columns in `spectroscopic_reference.csv`: the observed and model indices with their z, the Hα-SFR and Gallazzi comparisons, and the relation flags;
- a `validation_summary.csv`;
- the figures.

Like the Firefly check it replaces, it is a report, not a filter. The validation/test split remains "converged and covered".
