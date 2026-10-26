# Phase 3: ML posterior for 16-bin SFHs from S-PLUS colours — design

> **Superseded in part (2026-10-26):** the 16-bin targets (title, "Targets", the 19th-target amendment) and the A → B → C diagnostic amendments are replaced by `2026-10-26-delayed-tau-sfh-design.md` (9 flow targets, MLP only, the "works" test against photometric and spectral Bagpipes). Everything else here still holds. The 16-bin code is the git tag `nonparametric-16bin`.

**Date:** 2026-10-23
**Builds on:** `docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md` (training set v2, synthetic held-out sets, coverage test, spectroscopic reference set).

## Goal

A neural posterior estimator that maps S-PLUS photometry and redshift to a **full posterior** over each galaxy's SFH. The SFH is the mass fraction and SFR in each age bin, plus metallicity, consistent with Bagpipes / BC03-MILES under the training prior. It must work for galaxies with spectroscopic redshifts and for the general S-PLUS sample with photometric redshifts. Its outputs are judged by calibration and by agreement with the Bagpipes spectral posteriors of the reference set (report-only, no pass/fail gates).

## Decisions (agreed with the user)

- **Output:** full posteriors, using a conditional normalizing flow (simulation-based inference, neural posterior estimation).
- **Redshift:** one model takes the observed redshift and its uncertainty; σ_z = 0 for spectroscopic redshifts.
- **Inputs are colours, not magnitudes.**
  - 11 colours (band − r) carry the SFH information.
  - The r magnitude enters only to set the absolute mass and SFR.
  - The coverage test is already colour-based: its analytic free offset is equivalent to a colour fit with the correct error covariance, and it does not need r to be detected. The notebook text says so explicitly.
- **Targets:** 18 physical numbers (Approach A). The alternatives were rejected:
  - a posterior over the generative unit cube (27 nuisance-heavy dimensions, with discontinuities where the SFH branch switches);
  - a mixture-density network (less flexible, weaker bin-to-bin correlations).
- **Embedding networks:** an MLP and a 1-D CNN are both trained and compared on validation data only. The CNN is a controlled check, not an expected improvement.
- **Evaluation is report-only:** every metric is reported, and the user judges.
- **Compute:** training on the user's RTX 3090 through the `hljob` service; inference on CPU or GPU.

## Global constraints

- Code style: the ten rules of the parts 1–2 spec apply verbatim (English names, numbered banners, why-comments, no compact code, docstrings with inputs/units/returns, named constants, opening markdown per notebook, logic in modules / orchestration in notebooks, parallel cases in one pattern, concise).
- Libraries: `torch` (CUDA build on the GPU host) and `zuko` for the flows; no other ML framework.
- Prior and labels come only from `synthetic_model.py`. Phase 3 never re-implements the prior, the noise model or the label definitions.
- Locked data: `synthetic_test.parquet` and the reference set's `test` half are used once, in the final evaluation. No architecture, epoch or threshold choice may look at them.
- Generated weights are git-ignored (`ml_runs/*/model.pt`). Configs, training histories and the final catalogue summary are committed.
- With the current prior, 15 age bins are populated (bin index 15 is always padding), and bin 14 exists only below z ≈ 0.48.

## Inputs

For each galaxy (37 numbers before standardisation):
- 11 colours c_b = m_b − m_r for b ≠ r, in `BANDS` order;
- 11 colour errors σ_c,b = √(σ_b² + σ_r²);
- 11 detection flags (a colour is usable when band b and r are both detected; unusable colours are set to 0 after standardisation, with the flag 0);
- the r magnitude m_r and its error σ_r;
- z_obs and σ_z.

Galaxies without a detected r band are flagged `no_r_band` and not predicted (9 of 51k in DR3). Inputs are standardised with training-set means and standard deviations, which are stored in the model config.

## Targets

θ (18 numbers, standardised for the flow):
1. 14 additive log-ratios of the 15 bin mass fractions (bins 0–14), log10(f_i / f_13) for i ≠ 13. The reference is bin 13 (centre 5.455 Gyr), the oldest bin that exists at every redshift of the prior, so the transform is identical at every z. A bin that does not exist at the model's redshift (bin 14 above z ≈ 0.48) takes the prior floor `MIN_FRACTION` = 10⁻⁶ before the transform, so the dimension is fixed. Decoding exponentiates and renormalises to sum 1.
2. log10 Z_fin.
3. A_V (ISM).
4. log10 M_formed and log10 M_living [M_☉], absolute (see the brightness draw below).

**Derived for each posterior sample** by `synthetic_model` functions:
- log SFR per bin (M_formed · f_i / bin width in yr);
- metallicity per bin (enrichment law from Z_fin and the fractions);
- `mw_age`;
- `f_young`, `f_mid`, `f_old`.

For bins that do not exist at the galaxy's redshift, these labels are NaN or 0 exactly as in the training labels.

## Simulated observations (per epoch, on the fly)

From the noise-free `training_set.parquet` (2¹⁹ models):
1. **Brightness:** shift each model to the r magnitude of a random DR3 galaxy with |zml − z| < 0.05 (`scale_to_observed_r`). Mass targets shift by −0.4 Δm (`rescale_mass_labels`).
2. **Noise:** S-PLUS errors, flux noise and non-detections (`add_splus_noise` with `splus_noise_model.csv`).
3. **Redshift:**
   - with probability 0.5: σ_z = 0 and z_obs = z (spectroscopic);
   - otherwise σ_z = 0.025 (1 + z) and z_obs = z + N(0, σ_z) (photometric). This is the same σ_z rule as the coverage test's redshift window.

Validation uses `synthetic_validation.parquet` (fixed noise), in two fixed, seeded copies: spectroscopic z and noisy photometric z.

## Model

- **MLP embedding:** 37 inputs → 3 × 256 (SiLU) → 64-number summary.
- **CNN embedding:**
  - The 12 bands are ordered by pivot wavelength, with channels [m_b − m_r (0 for r), σ_b, detected, log pivot wavelength].
  - Then 3 × Conv1d (kernel 3, 64 channels, SiLU), global average pooling, and concatenation with [m_r, σ_r, z_obs, σ_z].
  - Then 2 × 256 (SiLU) → 64-number summary.
- **Flow:** `zuko` neural spline flow over the 18 standardised targets, with 8 transforms of 2 × 256 hidden units each, conditioned on the summary.
- **Loss:** the mean negative log-probability of θ.

## Training

- AdamW, learning rate 10⁻³ with cosine decay, batch size 4,096, up to 200 epochs.
- Early stopping on the validation loss (the mean over both redshift copies), with patience 15. The best weights are kept.
- Fixed seeds; the config stores the training-set commit hash, seeds, architecture, standardisation statistics and package versions.
- `train_ml_model.py --embedding mlp|cnn` writes `ml_runs/<embedding>/{model.pt, config.json, history.csv}`. It detects the device (CUDA or CPU).
- Both embeddings are trained with identical data, budget and seeds.
- **Model choice:** lower validation loss, with the synthetic-validation 68% and 95% coverage closest to nominal as a tie-breaker. Validation data only.

## Inference

`posterior_samples(model, galaxies, n_samples=1000, seed)` returns θ samples and derived labels per galaxy. Galaxies that fail the coverage test are flagged `out_of_domain` and not predicted.

## Evaluation (report-only, `ml_model.ipynb`)

During development the metrics use validation data. The final report repeats them once on the locked test data for the chosen model.
1. **Calibration** (synthetic, both redshift modes):
   - coverage curve (expected vs empirical central-interval coverage, levels 10–90%);
   - rank histograms of the truth among the posterior samples.
   Both cover `mw_age`, `Z_fin`, `f_young/f_mid/f_old`, `log_mstar`, `log_sfr`, and every `frac_bin` / `logsfr_bin`.
2. **Accuracy** (synthetic): median bias and NMAD of the posterior median against the truth per label, against r magnitude, redshift and SFH branch.
3. **Information per bin:** the posterior-to-prior ratio of 68% interval widths for each bin's fraction and SFR. Near 1 means prior-dominated; near 0 means well constrained.
4. **Agreement with the Bagpipes spectral posteriors** (headline): the ML posterior from 3″ aperture colours and spectroscopic z against each reference galaxy's spectral posterior.
   - Metrics: the median offset in units of the combined width, and the fraction of spectral medians inside the ML 68% and 95% intervals.
   - Splits: SF vs passive, SDSS vs BOSS. Development uses the `validation` half; the `test` half is used once.
   - Galaxies whose spectral fit is not converged, or whose `fitter_calibrated` is false, are reported separately.
5. **Firefly cross-check:** ML `mw_age` and mass-weighted metallicity against Firefly for the reference galaxies.
6. **Application:** photometric-z predictions for every covered DR3 galaxy, written to `ml_sfh_catalogue.parquet` (median, 16th and 84th percentiles per label, plus `out_of_domain` and `no_r_band` flags). The report includes the speed (galaxies per second on CPU and GPU) and the out-of-domain fraction.

## Files

All in `Synthetic_domain/`:

| File | Responsibility |
|---|---|
| `ml_data.py` | Colours and errors from magnitudes; per-epoch simulated observations; target encoding and decoding to labels; standardisation. |
| `ml_networks.py` | MLP and CNN embeddings, the `zuko` flow, the training loop with early stopping, checkpoints, `posterior_samples`. |
| `ml_evaluation.py` | Coverage curves, rank histograms, bias/NMAD, width ratios, comparison with spectral posteriors. |
| `train_ml_model.py` | Training CLI for the GPU host. |
| `ml_model.ipynb` | Part 3 story notebook: loads trained runs, model choice, evaluation, DR3 catalogue (it does not train). |
| `tests/test_ml_data.py`, `tests/test_ml_networks.py`, `tests/test_ml_evaluation.py` | CPU tests. |

## Compute host (`hljob`)

- Training runs as `hljob run --gpu …` jobs on the user's machine. Each job uploads a bundle (the modules, `train_ml_model.py`, `training_set.parquet`, `synthetic_validation.parquet`, `splus_noise_model.csv`, `dr3_corrected.csv` r/zml columns) into `inputs/` and writes `out/`.
- The first implementation step probes the job environment: Python, packages, CUDA, core and thread count, network access for `pip install`. It then fixes how dependencies are provided:
  - `synthetic_model.py` imports Bagpipes, so it must be installable there, or the training path must import only the Bagpipes-free parts.
- A job cannot exceed 24 h, 16 cores, or GPU memory above 8 GB unless ComfyUI is evicted.

## Testing (CPU, fast)

- **Encoding round trip:** fractions → log-ratios → fractions, including the floor for a missing bin. Labels decoded from θ match `bin_labels` / `summary_labels` for simulated models.
- **Simulated observations:** colours are invariant to the brightness shift; mass targets shift by −0.4 Δm; σ_z = 0 gives z_obs = z.
- **Flow recovery:** on a toy problem with an analytic Gaussian posterior, a small flow with each embedding reproduces the posterior mean and width within tolerance after a few hundred CPU steps.
- **Metrics on known cases:**
  - a calibrated toy posterior gives a diagonal coverage curve;
  - a truth equal to the median gives zero bias;
  - a posterior equal to the prior gives a width ratio of 1.
- **Checkpoints:** a save/load round trip gives identical samples for a fixed seed.

## Out of scope

- Hyper-parameter searches beyond the MLP-vs-CNN comparison, ensembles, and transformer embeddings (a later pass if the calibration plots call for it).
- Changing the prior, including the 16th-bin issue and any DR3/DR4 photometric re-baselining; these are open decisions recorded in the parts 1–2 ledger.
- Training on real galaxies (domain adaptation).

## Amendment (2026-10-24): the redshift is a 19th target

Decoding a posterior sample into per-bin labels needs the galaxy's redshift. The bin edges, the existence of bin 14 (z < 0.48) and the SFR per bin all depend on it. For photometric-redshift galaxies the true redshift is uncertain, so it becomes a learned target: θ gains `redshift` (19 numbers). Each posterior sample is decoded at its own sampled redshift, clipped to the prior range. For spectroscopic redshifts the posterior on z is simply narrow. Everything else in this spec is unchanged.

## Amendment (2026-10-26): diagnosing the gap between the ML and the spectral reference (A → B → C)

**Finding (validation run, 32 reference galaxies):**
- On synthetic data the MLP is calibrated: the 68% and 95% coverage of mw_age, Z_fin, log_mstar and log_sfr is within 0.02 of nominal.
- On real galaxies it disagrees with the Bagpipes spectral posteriors. ML minus spectral medians are:
  - mw_age +2.8 Gyr (9.2 against 6.3);
  - Z_fin 0.30 against 1.63 Z_sun;
  - log_sfr −1.1 dex;
  - log_mstar +0.24 dex.
- The inputs are not the cause:
  - the units agree;
  - DR3 total colours give the same result;
  - the galaxies lie inside the training colour domain (nearest-model χ²/dof ≈ 0.8);
  - the 30 training models nearest in colour are 9.2 Gyr old, as the ML says.
- The spectral fits use a calibration polynomial, so they read the absorption features, not the continuum shape.

**Decision (user, 2026-10-26):** try A. If A does not work, try B. If B does not work, try C.

**What "works" means.** It is the same test at each step: the ML posterior against a reference label set, on the `validation` half only. The `test` half stays locked.
- The labels are mw_age, Z_fin, f_old and log_mstar.
- A label agrees when both hold:
  - at least 80% of the reference medians lie inside the ML 95% interval;
  - the median of (ML median − reference median) / combined 1σ width is at most 1.0 in absolute value.
- The step works when all four labels agree.
- Every step also reports the spectral comparison, which stays the headline physical check.

### A. Bagpipes fits of the photometry alone (the ML's exact target)

- **What is fitted:** each reference galaxy's 12-band 3″ aperture photometry, the ML's own input.
  - The redshift is fixed at the spectroscopic z.
  - The prior is the training prior (`components(cube, fixed_redshift=z)`): DIMS without z, and no veldisp or noise_scale.
  - The sampler is MultiNest, with the same settings as the spectral fits.
- **Likelihood:** it must match the simulator that made the training data (`add_splus_noise`). It is Gaussian in magnitudes with the catalogue errors. No ERROR_FLOOR is added, because the simulator adds none.
  - Only detected bands are used (`detected_bands`).
  - The overall brightness is a free magnitude offset with a flat prior. It is marginalised analytically: χ²_min of the offset plus ½ ln(2π / Σw).
  - The mass labels of each posterior sample are rescaled by an offset drawn from its Gaussian conditional, through `rescale_mass_labels`, as in training.
- **Which galaxies and where it runs:**
  - Every fitted galaxy of the `validation` half, plus a calibration check: the 33 calibration mocks' noiseless photometry, noised with `add_splus_noise`. The check passes when 95% coverage is at least 0.85.
  - The fits are jobs on the PC, written to `photometric_fits/` (separate from `spectral_fits/`, because job ids are specobjids).
- **Reading the result:**
  - **A works**, i.e. the ML agrees with photometric Bagpipes: the ML is a faithful emulator of Bagpipes on these colours. The spectral gap is then about how much information 12 bands carry, a scientific result, and the ML evaluation reports both comparisons.
  - **A does not work:** go to B.

### B. Reference labels from joint spectrum + photometry fits

- **What is fitted:** spectrum and aperture photometry together, with the same prior.
  - The spectrum keeps its calibration polynomial.
  - The photometry fixes the continuum shape.
  - The likelihood is the sum of the spectral likelihood and the A photometric likelihood.
- **Galaxies:** the `validation` half, as jobs on the PC.
- **Reading the result:**
  - **B works**, i.e. the ML agrees with the joint labels: the joint fits become the ML reference labels.
  - **B does not work:** go to C.

### C. A training prior weighted towards young, star-forming, metal-rich galaxies

The aim is better coverage of S-PLUS galaxies:
- Raise the weight of the star-forming and young SFH branches.
- Raise the share of models with high Z_fin (super-solar).
- Choose the weights by comparing training-set colours with DR3 colour distributions, using the coverage test and colour histograms.

This lifts the "changing the prior" out-of-scope item for this step only.
- The training set is regenerated on the PC, and the MLP is retrained on the GPU.
- The agreement test is repeated against the B labels, or the spectral labels if B was not built.
- The synthetic validation and test sets are regenerated from the new prior too, so that synthetic calibration is still measured on the training distribution.

**Compute rule (user, 2026-10-25):** every fit, notebook execution, training run and posterior sampling runs as an `hljob` job on the PC, never on the laptop.
