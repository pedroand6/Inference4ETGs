# Double-power-law SFH alongside delayed-τ — design

**Date:** 2026-10-28
**Extends:** `2026-10-26-delayed-tau-sfh-design.md` and its amendment.

That spec still holds for the delayed-τ model. Its out-of-scope item "other parametric SFHs" is lifted for the double power law (DPL) only. This spec adds a second SFH model and the comparison of the two.

## Goal

The pipeline runs two SFH models side by side:
- delayed-τ: kept as it is, with its sets, fits and flow;
- Bagpipes' double power law (`dblplaw`): new.

For DPL, the training prior covers its whole parameter space in the same way as delayed-τ (a wide Latin hypercube, log-uniform where a parameter spans decades). DPL gets its own:
- training and synthetic validation/test sets;
- photometric Bagpipes fits and calibration mocks;
- ML flow (MLP embedding).

The final report compares, for the two models:
- ML against Bagpipes photometry;
- the two models against each other;
- both models against the spectra;
- the Bayesian evidence for each model, galaxy by galaxy.

## Decisions (user, 2026-10-28)

1. **Photometric fits only for DPL.** The spectral reference stays delayed-τ:
   - the 1359 spectra (755 fitted so far, 604 pending);
   - the 33 spectral mocks.

   DPL is fitted to the S-PLUS aperture photometry only.
2. **DPL prior.** The peak time is drawn the same way as the delayed-τ age: log-uniform in lookback time.
3. **The comparison reports all four parts:**
   - each model's ML against its own photometric Bagpipes;
   - delayed-τ against DPL on shared labels;
   - both models against the spectral reference;
   - the Bayesian evidence per galaxy.
4. **Delayed-τ results are kept, not recomputed.** This covers the sets, the fits, the mocks and `ml_runs/mlp`, under their current file names. The refactor must reproduce the stored delayed-τ sets exactly (see Testing).

## The DPL model

Bagpipes' `dblplaw` component:

SFR(t) ∝ [(t/τ)^α + (t/τ)^−β]^−1

- t is the time since the Big Bang.
- τ [Gyr] is the time of the peak, also counted from the Big Bang.
- α is the falling slope after the peak; β is the rising slope before it.
- Star formation starts at the Big Bang, so there is no age parameter.
- There is one metallicity, as in the delayed-τ model.

**Prior** (Latin hypercube, same sizes and seeds as delayed-τ):

| Parameter | Range | Scale |
|---|---|---|
| `redshift` | 0.002 – 0.7 | uniform |
| `peak_lookback` (time from the SFR peak to the observation) | 0.05 Gyr – t_max(z) = age of the Universe at z − `UNIVERSE_AGE_MARGIN` | log-uniform; drawn as a unit-interval fraction of the log range, like the delayed-τ `age` |
| `alpha` | 0.01 – 1000 | log-uniform (Carnall et al. 2018) |
| `beta` | 0.01 – 1000 | log-uniform (Carnall et al. 2018) |
| `metallicity`, `av_ism`, `tau_bc`, `logU` | as delayed-τ | as delayed-τ |
| `log_mass_formed` | fixed at 10, then rescaled to the observed r magnitude as for delayed-τ | — |

- **Bagpipes' τ:** τ = age of the Universe at z − `peak_lookback`. It is at least `UNIVERSE_AGE_MARGIN`, so it is always positive. It always lies below Bagpipes' "unphysical" limit of 1.5 × the age of the Universe.
- **Kept unchanged:** dust (CF00), nebular emission, grids, IMF, cosmology and the S-PLUS noise model are those of delayed-τ.

**`DIMS`:**
`["z", "log_peak_lookback_fraction", "log_alpha", "log_beta", "log_metallicity", "av_ism", "log_tau_bc", "logU"]`

`PHOTOMETRIC_DIMS` is the same list without z.

**Labels:**
- **Drawn:**
  - `redshift`, `peak_lookback` [Gyr], `tau` [Gyr since the Big Bang];
  - `alpha`, `beta`, `metallicity`;
  - `av_ism`, `tau_bc`, `logU`, `log_mass_formed`.
- **Derived:**
  - `log_peak_lookback`, `log_alpha`, `log_beta`, `log_metallicity`;
  - `log_mstar` (from Bagpipes);
  - `log_sfr` and `log_ssfr` (the 100-Myr mean, with the floor, as in the delayed-τ amendment);
  - `mw_age`, `f_young`, `f_mid`, `f_old`.

**Derived labels without a closed form.** For DPL, the mass fractions and `mw_age` have no analytic expression. They are summed numerically on Bagpipes' own age grid (its `ages` and `age_widths`), vectorised over samples:
- they match Bagpipes' `sfh` arrays, including its discretisation of sharp peaks (α or β near 1000);
- the same function serves the training labels and the ML decoding;
- `log_sfr` comes from the same sum (formed mass younger than 100 Myr).

## Code structure

**New module `sfh_models.py`.** It holds one section per model. Each section defines:
- its name (`"delayed"`, `"dblplaw"`);
- `dims` and `photometric_dims`;
- `components(u, fixed_redshift=None)`;
- drawn and derived label names, and the derived-label function;
- `target_names` (flow targets) and `calibration_labels`;
- the file suffix.

`MODELS = {"delayed": ..., "dblplaw": ...}` maps names to models. The delayed-τ section is the current code (`components`, `delayed_tau_labels` and the rest) moved out of `synthetic_model.py`, with the same behaviour.

**Model-dependent code takes the model by name.** This covers:
- `simulate_models`, `summary_labels`;
- the fitters and job export;
- the calibration mocks;
- `fit_is_current`, `supersede_outdated_fits`;
- `ml_data` encode/decode, `train_ml_model.py`.

The global `SFH_MODEL`, `DIMS` and `PHOTOMETRIC_DIMS` of `synthetic_model.py` go.

**Jobs and fits store their model:**
- Each fit job writes `sfh_model`, and the fitter builds the matching Bagpipes components from it.
- A job without that key is delayed-τ (every existing job).
- `fit_is_current` compares a fit's stored `sfh_model` with the model of its folder. A DPL fit is therefore never superseded as "outdated" in the DPL folder, and the reverse holds too.

**File names.** Delayed-τ keeps today's names. DPL adds the suffix `_dblplaw`:

| Product | Delayed-τ | DPL |
|---|---|---|
| Sets | `training_set.parquet`, `synthetic_validation.parquet`, `synthetic_test.parquet` | the same names with `_dblplaw` |
| Photometric fits, jobs, mocks | `photometric_fits/` | `photometric_fits_dblplaw/` (via `FITS_SUBDIR`) |
| Flow run | `ml_runs/mlp/` | `ml_runs/mlp_dblplaw/` |
| DR3 catalogue | `ml_sfh_catalogue.parquet` | `ml_sfh_catalogue_dblplaw.parquet` |

## Part 1: training sets (`training_set.ipynb`)

`SFH_MODELS = ["delayed", "dblplaw"]` in the parameter cell; the notebook builds the sets of each listed model.
- **Delayed-τ:** the notebook skips a model whose three parquet files exist, unless a rebuild flag is set. The stored delayed-τ sets are therefore not rewritten.
- **DPL:** the sets are new, with the same sizes, seeds and S-PLUS noise as delayed-τ.

## Part 2: DPL photometric fits (`spectroscopic_reference.ipynb`, section 6b)

**Galaxies.** DPL is fitted to the same galaxies as the delayed-τ photometric fits: the validation half. The test half follows, as for delayed-τ, only with the locked test.

**Calibration mocks.** There are no DPL spectral mocks to reuse, so:
- 33 DPL photometric mocks are drawn from the DPL prior;
- each takes the redshift and aperture brightness of a validation galaxy (a brightness donor, as now);
- each is noised once with `sm.add_splus_noise`.

The test passes when 95% coverage is at least 0.85 for `log_peak_lookback`, `log_alpha`, `log_beta`, `log_metallicity` and `log_mstar`; 68% coverage is reported too.

**Cross-fits for the evidence test:**
- the 33 delayed-τ photometric mocks are fitted with DPL;
- the 33 DPL mocks are fitted with delayed-τ.

The resulting 66 cheap photometric fits test whether the colours can tell the two models apart when the truth is known. Their ids are `cross_calib_<k>`, and each job records the true model.

**Results.** These columns go to `spectroscopic_reference.csv`:
- `phot_dblplaw_<label>_p16/_median/_p84` for the DPL calibration labels and the shared labels;
- `phot_dblplaw_fitter_calibrated`;
- `phot_dblplaw_log_evidence`, next to a new `phot_log_evidence` for delayed-τ.

**Running.** Jobs go to `photometric_fits_dblplaw/jobs/` and run on the PC (`FITS_SUBDIR=photometric_fits_dblplaw remote/submit_fits.sh`). The local supersede step runs before every notebook-2 job, as now.

## Part 3: DPL flow and evaluation (`ml_data.py`, `train_ml_model.py`, `ml_model.ipynb`)

**Flow targets (10, in this order).**
- log_peak_lookback, log_alpha, log_beta, log_metallicity;
- av_ism, log_tau_bc, logU;
- log_mass_formed, log_mstar;
- redshift.

**Inputs, embedding and training:**
- the inputs are unchanged;
- the embedding is the MLP;
- training is unchanged: on the GPU, with `train_ml_model.py --sfh-model dblplaw`.

**Checkpoints.** Each checkpoint records its `target_names` and `sfh_model`. Loading one into the wrong model raises an error.

**Decoding.** Each posterior sample is decoded with the DPL derived-label function at its own redshift.

**Per-model evaluation** (`ml_model.ipynb` loops over `SFH_MODELS`; report-only, the validation half):
- synthetic calibration and accuracy;
- the "works" test against the model's own photometric Bagpipes, on its calibration labels plus the shared labels;
- its DR3 catalogue.

## Part 4: model comparison (`ml_model.ipynb`, new last section)

**Shared, model-independent labels:**
- `log_mstar`, `mw_age`, `log_metallicity`;
- `av_ism`, `log_ssfr`;
- `f_young`, `f_mid`, `f_old`.

**1. ML against its own Bagpipes.** The two models' "works" tests are shown side by side: label × model, giving inside_95, the median offset in combined σ, and whether the label agrees.

**2. Delayed-τ against DPL on the shared labels.**
- Pairs compared:
  - the ML of each model, on the validation galaxies and on DR3;
  - the photometric Bagpipes of each model, on the validation galaxies.
- Statistics: the median difference, the scatter, and the fraction of galaxies whose 68% intervals overlap.

**3. Against the spectra.** The "works" test of four estimators against spectral Bagpipes, on the shared labels:
- ML delayed-τ, ML DPL;
- photometric Bagpipes delayed-τ, photometric Bagpipes DPL.

Results are split by passive and star-forming galaxies, as now. The notebook states the caveat that the spectral reference assumes delayed-τ, so its `mw_age` and mass fractions carry that SFH shape.

**4. Bayesian evidence.** ln B = `phot_log_evidence` − `phot_dblplaw_log_evidence` per validation galaxy (positive favours delayed-τ). The notebook reports:
- the distribution of ln B, overall and by passive/star-forming;
- the fraction of galaxies with |ln B| < 1 (inconclusive), 1–3 (positive) and > 3 (strong) for each model;
- the MultiNest evidence uncertainty, √(H / n_live);
- on the 66 cross-fitted mocks, how often the true model is preferred. This shows whether the colours can separate the models at all.

**Locked test.** It runs once, at the end, for both models. As before, this happens after the user has inspected the validation results.

## Compute (user rule, 2026-10-25)

The laptop runs only code changes, unit tests and quick checks. Everything else is an `hljob` job on the PC, one at a time:
1. `training_set.ipynb` (CPU): the DPL sets.
2. DPL MLP training (GPU, `MEMORY_MB=7100`).
3. Notebook 2 (CPU): exports the DPL photometric jobs, the DPL mocks and the cross-fits. It runs after the local supersede step and with a small bundle (only fitted spectra; hljobd has a 3 GB limit).
4. Photometric fits, DPL folder (CPU, 16 cores, 15 workers), then the 33 delayed-τ cross-fits in the delayed folder.
5. Notebook 2 (CPU): reads all fits.
6. `ml_model.ipynb` (GPU): both models and the comparison.

The 604 pending delayed-τ spectral fits are independent of this work. When they run, they enlarge the validation half, and both models' photometric fits then grow with it.

## Testing (CPU, fast)

**DPL prior map:**
- every parameter within bounds;
- uniform cube values give uniform log values for the log-uniform parameters;
- `peak_lookback` never exceeds t_max(z), and τ > 0;
- the Latin hypercube is stratified.

**DPL derived labels:**
- mass fractions, `mw_age` and `log_sfr` match a sum over Bagpipes' `sfh` arrays for several (α, β, τ) values, to 1%. The cases include a sharp peak (α = β = 1000), a flat SFH (α = β = 0.01), an early peak and a recent peak;
- the fractions sum to 1.

**Delayed-τ unchanged.** After the move to `sfh_models.py`:
- the delayed-τ path regenerates the first rows of the stored `training_set.parquet` and `synthetic_validation.parquet` exactly (to float tolerance);
- the existing delayed-τ tests pass unchanged except for import paths.

**Targets:** for both models, the encode/decode round trip is exact, and the decoded labels match `summary_labels` for simulated models.

**Jobs and fits:**
- a DPL job builds `dblplaw` components;
- a job without `sfh_model` is delayed-τ;
- `fit_is_current` keeps DPL fits in the DPL folder and supersedes them in the delayed-τ folder;
- a DPL photometric smoke fit runs on a noiseless mock.

**Evidence and comparison:**
- ln B and its categories are computed correctly on a small hand-made table;
- comparison statistics are computed on toy data with a known offset.

**Checkpoints:** a checkpoint whose `sfh_model` or `target_names` differs from the requested model raises an error.

**Notebooks:** valid JSON, and every code cell compiles. They are executed only on the PC.

## Out of scope

- DPL spectral fits and DPL spectral mocks.
- Other SFH models (exponential, log-normal, non-parametric).
- Model averaging or choosing one model per galaxy in the catalogue: the evidence is reported, not used.
- Changing any delayed-τ setting, prior or result.
- The CNN, and hyper-parameter searches.
