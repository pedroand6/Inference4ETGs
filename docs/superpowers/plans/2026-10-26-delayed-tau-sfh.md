# Delayed-τ SFH Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 16-bin non-parametric SFH and the enrichment law with a single Bagpipes `delayed` SFH and one metallicity, everywhere in the pipeline. The pipeline covers the prior, the labels, the training and synthetic sets, the reference fits and calibration, the validation rulers and the ML.

**Architecture:**
- `synthetic_model.py` stays the single source of the prior, the labels and the fits. Sections 1–2 are rewritten, and sections 5–7 are adapted to the new `DIMS`.
- New analytic delayed-τ functions in section 2 give `mw_age`, `f_young/f_mid/f_old` and `log_sfr`. Part 1 labels and ML decoding both use them, so labels and decoded samples agree exactly.
- The ML predicts 9 physical targets; everything else is derived.
- The notebooks are adapted but never executed on the laptop. They run as `hljob` jobs on the PC (see Operations).

**Tech Stack:** Python 3.10, Bagpipes 1.3.5 (bc03_miles), scipy (`gammainc`, `gammaincc`), PyMultiNest, torch 2.5.1 + zuko 1.2.0, pytest, `hljob`.

**Spec:** `docs/superpowers/specs/2026-10-26-delayed-tau-sfh-design.md`, including its "Amendment (2026-10-26, from the code inventory)".

**Line-level inventories:** these are mandatory reading for each task. They list every symbol, its line, its callers and the tests to change.
- `.superpowers/sdd/2026-10-26-delayed-tau-sfh/inventory-prior-labels.md`: `synthetic_model.py` sections 1–2, `training_set.ipynb`, the Bagpipes `delayed` component.
- `.superpowers/sdd/2026-10-26-delayed-tau-sfh/inventory-fits-reference.md`: sections 4–7, `run_spectral_fits.py`, `spectroscopic_reference.ipynb`, `spectral_validation.py`, `fetch_literature.py`.
- `.superpowers/sdd/2026-10-26-delayed-tau-sfh/inventory-ml.md`: `ml_data.py`, `ml_networks.py`, `train_ml_model.py`, `ml_evaluation.py`, `ml_model.ipynb`.
- `.superpowers/sdd/2026-10-26-delayed-tau-sfh/inventory-tests-scripts-docs.md`: every test (keep, rewrite or delete), the remote scripts, the README, `.gitignore`, the generated data.

## Global Constraints

- **Code style:** the ten rules of the parts 1–2 spec (`docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md`):
  - English names, numbered section banners and why-comments;
  - no compact code;
  - docstrings with inputs, units and returns;
  - named constants;
  - logic in modules, orchestration in notebooks;
  - parallel cases in one pattern;
  - concise, not prolix;
  - lines ≤ ~110 characters.
- **Prior (verbatim from the spec):**
  - `redshift` 0.002–0.7, uniform;
  - `age` 0.05 Gyr to t_max(z) = age of the Universe at z − 0.1 Gyr, log-uniform;
  - `tau` 0.01–20 Gyr, log-uniform;
  - `metallicity` 0.005–2.5 Z☉, log-uniform;
  - `av_ism` 0.01–2.0 mag, uniform;
  - `tau_bc` 0.5–4, log-uniform;
  - `logU` −4 to −1, uniform;
  - `log_mass_formed` 10, rescaled to the observed r.
- **Unchanged:** CF00 dust (n = 0.7, t_bc = 0.01 Gyr), nebular emission, bc03_miles, Kroupa IMF, cosmology (H0 = 70, Ωm = 0.3), the noise model, the coverage test, the spectral and photometric likelihoods, and the MultiNest settings.
- **Names:**
  - `DIMS = ["z", "log_age_fraction", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU"]`.
  - The drawn label keys are `redshift, age, tau, metallicity, av_ism, tau_bc, logU, log_mass_formed`. The redshift is never called `"z"`, which is the z-band column.
  - The derived labels are `log_age, log_tau, log_metallicity, log_mstar, log_sfr, log_ssfr, mw_age, f_young, f_mid, f_old`.
- **ML targets, in order:** `["log_age", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU", "log_mass_formed", "log_mstar", "redshift"]`. The redshift is last. `log_sfr` is derived.
- **Calibration and agreement labels:** `["log_age", "log_tau", "log_metallicity", "log_mstar"]`, with `mw_age` reported too. The thresholds are unchanged: 95% coverage ≥ 0.85; inside 95% ≥ 0.80 and |median offset| ≤ 1.0.
- **Compute:** NEVER execute a notebook, a set generation, a real fit, training or posterior sampling on the laptop. Unit tests and checks of a few seconds are allowed. The controller submits every PC job.
- **Safety:** never print `.env` and never pack it into a bundle. Do not touch `Synthetic_domain/pipes/` or the repo-root PDF. Never push.
- **Commits:** every message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Tests:** run with `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q`. They take about 5 minutes.

## Review Focus

1. **Extreme corners of the prior:** τ = 0.01 with age = t_max (age/τ ≈ 1360), and τ = 20 with age = 0.05 (age/τ = 0.0025). The analytic labels must be finite, without overflow or 0/0, and fractions must lie in [0, 1] and sum to 1.
2. **No star formation in the last 100 Myr:** `log_sfr` sits exactly at the floor and `log_ssfr = log_sfr − log_mstar` stays finite. Never use Bagpipes' `ssfr`, which is NaN there.
3. **Age at the redshift edge:** at z = 0.7 the drawn age never exceeds t_max(z). Bagpipes would otherwise flag the model unphysical and return NaN magnitudes. Also, a decoded posterior sample at its own sampled redshift is clipped into [0.05, t_max(z)] only when it lies outside that range.
4. **Old artefacts:**
   - an old `.h5` fit without `sfh_model = "delayed"` is superseded, never read;
   - an old ML checkpoint without matching `target_names` raises;
   - tests reading old parquets skip with a clear reason.
5. **Reused mock ids:** the `calib_<k>` ids of the new calibration mocks are not skipped (`skipped_fit_jobs.txt` is emptied), and old `calib_*` fits never mix into the new calibration.

---

### Task 1: Prior, analytic delayed-τ labels and simulation (`synthetic_model.py` sections 1–3)

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py`. Section 1 is lines 27–318 and section 2 is lines 321–498; see `inventory-prior-labels.md` A, B and F.
- Test: `Synthetic_domain/tests/test_synthetic_model.py`. Delete or rewrite the bin tests listed in `inventory-tests-scripts-docs.md` A7. Delete `tests/make_legacy_fixture.py` and its fixture (A8).

**Interfaces:**
- **Produces, constants:**
  - `AGE_MIN = 0.05`, `TAU_RANGE = (0.01, 20.0)`, `METALLICITY_RANGE = (0.005, 2.5)`;
  - `DIMS` as given in Global Constraints;
  - `DRAWN_LABELS = ["redshift", "age", "tau", "metallicity", "av_ism", "tau_bc", "logU", "log_mass_formed"]`;
  - `DERIVED_LABELS = ["log_age", "log_tau", "log_metallicity", "log_mstar", "log_sfr", "log_ssfr", "mw_age", "f_young", "f_mid", "f_old"]`;
  - `TRAINING_LABELS = DRAWN_LABELS + DERIVED_LABELS`;
  - `MASS_DEPENDENT_LABELS = ["log_mass_formed", "log_mstar", "log_sfr"]`.
- **Produces, functions:**
  - `maximum_age(redshift) -> float` [Gyr], equal to `COSMOLOGY.age(z) - UNIVERSE_AGE_MARGIN`. It accepts arrays too.
  - `components(u, fixed_redshift=None) -> (bagpipes_components, drawn)`. `bagpipes_components["delayed"] = {"age", "tau", "massformed": LOG_MASS_FORMED, "metallicity"}`; the dust and nebular entries are as before. `drawn` has the `DRAWN_LABELS` keys.
  - `delayed_mass_fraction_younger_than(age, tau, lookback) -> ndarray`, using the regularised-gamma form from the inventory (section D).
  - `delayed_mass_weighted_age(age, tau) -> ndarray` [Gyr ago].
  - `delayed_tau_labels(age, tau, log_mass_formed, log_mstar) -> dict`. It gives `mw_age`, `f_young`, `f_mid`, `f_old`, `log_sfr` (analytic 100 Myr average with the `SFR_FLOOR × 10**log_mstar` floor) and `log_ssfr`. It is vectorised over arrays.
  - `summary_labels(bagpipes_sfh, drawn) -> dict` with every `TRAINING_LABELS` key. `log_mstar` comes from Bagpipes and the rest from `delayed_tau_labels` and `drawn`.
  - `simulate_models(unit_vectors) -> DataFrame`. Its columns are `TRAINING_LABELS`, then `BANDS`, in that order.
  - `rescale_mass_labels(table, magnitude_offset)`, unchanged in behaviour, over the new `MASS_DEPENDENT_LABELS`.
- **Removes:**
  - the constants `AGE_BIN_CENTRES`, `N_BINS`, `LOG_ZFIN_RANGE`, `ZINI_FRACTION`, `P_SKEW`, `P_BURST`, `MIN_FRACTION`, the `SKEW_*` and `BURST_*` constants, and `DIRICHLET_DIMS`;
  - the functions `positive_normal`, `bin_edges`, `bin_edges_array`, the three branch functions, `SFH_BRANCHES`, `choose_sfh_branch`, `bin_metallicities`, `bin_metallicities_array`, `bins_from_components`, `bin_labels`, `BIN_LABEL_COLUMNS`, `SUMMARY_LABEL_COLUMNS` and `mass_weighted_metallicity`;
  - the unused imports `norm` and `skewnorm`.

- [ ] **Step 1: Write failing tests.** Example core tests:

```python
def test_prior_draws_stay_inside_bounds_and_age_below_universe_age():
    rng = np.random.default_rng(0)
    for u in rng.random((500, len(sm.DIMS))):
        _, drawn = sm.components(u)
        assert sm.Z_RANGE[0] <= drawn["redshift"] <= sm.Z_RANGE[1]
        assert sm.AGE_MIN <= drawn["age"] <= sm.maximum_age(drawn["redshift"]) + 1e-12
        assert sm.TAU_RANGE[0] <= drawn["tau"] <= sm.TAU_RANGE[1]
        assert sm.METALLICITY_RANGE[0] <= drawn["metallicity"] <= sm.METALLICITY_RANGE[1]


def test_uniform_cube_gives_log_uniform_tau_and_metallicity():
    cube = np.full(len(sm.DIMS), 0.5)
    _, drawn = sm.components(cube)
    assert np.isclose(np.log10(drawn["tau"]), np.mean(np.log10(sm.TAU_RANGE)))
    assert np.isclose(np.log10(drawn["metallicity"]), np.mean(np.log10(sm.METALLICITY_RANGE)))


@pytest.mark.parametrize("age, tau", [(13.0, 0.3), (5.0, 1.0), (2.0, 20.0), (8.0, 3.0), (0.5, 0.1), (0.05, 0.01)])
def test_analytic_labels_match_bagpipes_sfh(age, tau):
    component = {"age": age, "tau": tau, "massformed": 10.0, "metallicity": 1.0}
    sfh = pipes.models.star_formation_history({"redshift": 0.002, "delayed": component})
    mass = sfh.sfh * sfh.age_widths
    ages = sfh.ages / 1e9
    assert np.isclose(sm.delayed_mass_weighted_age(age, tau), sfh.mass_weighted_age, rtol=0.01)
    for limit in (0.1, 1.0):
        younger = mass[ages < limit].sum() / mass.sum()
        assert abs(sm.delayed_mass_fraction_younger_than(age, tau, limit) - younger) < 0.01


def test_analytic_labels_are_finite_in_the_prior_corners():
    t_max = sm.maximum_age(0.002)
    for age, tau in [(t_max, 0.01), (0.05, 20.0), (t_max, 20.0), (0.05, 0.01)]:
        labels = sm.delayed_tau_labels(np.array([age]), np.array([tau]), np.array([10.0]), np.array([9.8]))
        assert all(np.isfinite(value).all() for value in labels.values())
        fractions = labels["f_young"] + labels["f_mid"] + labels["f_old"]
        assert np.allclose(fractions, 1.0)


def test_no_recent_star_formation_sits_on_the_sfr_floor():
    labels = sm.delayed_tau_labels(np.array([13.0]), np.array([0.01]), np.array([10.0]), np.array([9.7]))
    assert np.isclose(labels["log_sfr"][0], np.log10(sm.SFR_FLOOR) + 9.7)
    assert np.isfinite(labels["log_ssfr"][0])


def test_simulated_columns_follow_training_labels_then_bands():
    table = sm.simulate_models(np.full((2, len(sm.DIMS)), 0.5))
    assert list(table.columns) == sm.TRAINING_LABELS + sm.BANDS
    assert np.isfinite(table[sm.BANDS].to_numpy()).all()
```

  Keep `test_simulated_redshift_and_z_band_are_distinct_columns`. Add a test that the analytic `log_sfr` matches Bagpipes `sfh.sfr` to 0.02 dex for (5, 1) and (2, 20).

- [ ] **Step 2: Run the tests and see them fail.** Run `pytest tests/test_synthetic_model.py -q -k "prior or analytic or floor or simulated"`.
- [ ] **Step 3: Implement.**
  - **Analytic formulas** (y = A/τ; P and Q are `scipy.special.gammainc` and `gammaincc`):

```python
def delayed_mass_fraction_younger_than(age, tau, lookback):
    """Share of the formed mass with lookback time below `lookback` for SFR(t) ∝ t e^{-t/τ}, t since the start.

    Why the Q form: Q(2, (A-L)/τ) - Q(2, A/τ) avoids the cancellation of 1 - P/P when τ << A.
    Inputs: age A, tau, lookback L [Gyr], scalars or arrays. Returns: fraction in [0, 1] (dimensionless).
    """
    age, tau, lookback = np.broadcast_arrays(np.asarray(age, float), np.asarray(tau, float), np.asarray(lookback, float))
    formed_before_now = gammainc(2, age / tau)
    older_start = np.clip(age - lookback, 0, None) / tau
    fraction = (gammaincc(2, older_start) - gammaincc(2, age / tau)) / formed_before_now
    return np.where(lookback >= age, 1.0, np.clip(fraction, 0.0, 1.0))


def delayed_mass_weighted_age(age, tau):
    """Formed-mass-weighted lookback age [Gyr]: A - 2τ P(3, A/τ) / P(2, A/τ)."""
    age, tau = np.asarray(age, float), np.asarray(tau, float)
    return age - 2 * tau * gammainc(3, age / tau) / gammainc(2, age / tau)
```

    - `log_sfr = log10(max(10**log_mass_formed * F(0.1) / 1e8, SFR_FLOOR * 10**log_mstar))`.
    - Write `GYR_TO_YEARS` usage only inside functions, or move the constant above section 2 (inventory F.5).
  - **`components`:**
    - `t_max = maximum_age(redshift)`;
    - `age = 10**uniform(params["log_age_fraction"], log10(AGE_MIN), log10(t_max))`;
    - `tau` and `metallicity` log-uniform over their ranges;
    - keep the dust and nebular code;
    - `drawn["tau_bc"] = 10**uniform(params["log_tau_bc"], *LOG_TAU_BC_RANGE)`; A_V,bc = TAU_TO_AV × tau_bc stays inside the dust `eta`;
    - `drawn["log_mass_formed"] = LOG_MASS_FORMED`.
  - **`simulate_models`:** cache one Bagpipes model (single component set). Build each row in `TRAINING_LABELS` order, then add the magnitudes.
- [ ] **Step 4: Run the whole `test_synthetic_model.py` and see it pass.** Run it with section 5–7 tests that fail only because of the not-yet-adapted fit code marked xfail with reason "Task 2". Task 2 removes those marks.
- [ ] **Step 5: Commit.** Message: `feat: delayed-tau prior with a single metallicity and analytic labels replace the 16-bin SFH`.

### Task 2: Fits, mocks, calibration and validation rulers (`synthetic_model.py` sections 5–7, `spectral_validation.py`)

**Files:**
- Modify:
  - `Synthetic_domain/synthetic_model.py` sections 5–7 (`inventory-fits-reference.md` B, C, D);
  - `Synthetic_domain/spectral_validation.py` (G);
  - `Synthetic_domain/run_spectral_fits.py`, but only if the inventory section E finds a coupling (it says none).
- Test: `Synthetic_domain/tests/test_synthetic_model.py` (sections 5–7 tests), `Synthetic_domain/tests/test_spectral_validation.py`.

**Interfaces:**
- **Consumes (Task 1):** `DIMS`, `components`, `summary_labels(bagpipes_sfh, drawn)`, `delayed_tau_labels`, `TRAINING_LABELS`, `MASS_DEPENDENT_LABELS`.
- **Produces:**
  - `SPECTRAL_DIMS` (8 entries) and `PHOTOMETRIC_DIMS` (6 entries), rebuilt automatically.
  - `SFH_MODEL = "delayed"`. `save_fit` stores `metadata["sfh_model"] = SFH_MODEL` on every fit, and `load_fit` returns it.
  - `fit_is_current(metadata) -> bool`: True only when `sfh_model == SFH_MODEL` and `max_iterations >= MAX_ITERATIONS`.
  - `posterior_labels(unit_cube_samples, redshift)` and `photometric_posterior_labels(...)` return the `TRAINING_LABELS` columns (plus `magnitude_offset` for the photometric version), through the rewritten `sample_label_row(model, bagpipes_components, drawn)`.
  - `CALIBRATION_LABELS = ["log_age", "log_tau", "log_metallicity", "log_mstar"]`.
  - `draw_calibration_cubes` and `draw_calibration_mocks` keep their signatures. Their truths carry `TRAINING_LABELS`.
  - `spectral_validation.light_weighted_age_and_metallicity(unit_cube_samples, redshift, n_samples)` keeps its output keys (`log_age_light`, `log_metallicity_light`).
    - **Age:** the r-band light-weighted log of the linear mean age (Gallazzi 2005). Weights = formed mass per Bagpipes age step × the SSP r-band luminosity per unit mass at that age and the galaxy's Z.
    - **Metallicity:** log10 Z exactly, since there is a single metallicity.

- [ ] **Step 1: Write failing tests.**
  - `fit_is_current` is False for metadata without `sfh_model` and True for `{"sfh_model": "delayed", "max_iterations": sm.MAX_ITERATIONS}`.
  - `save_fit`/`load_fit` round-trips `sfh_model`.
  - `posterior_labels` on 3 midpoint cubes returns the `TRAINING_LABELS` columns, all finite.
  - `calibration_coverage` runs on the new labels.
  - The light-weighted metallicity of a model equals log10 Z.
  - The light-weighted age of a τ = 0.01, age = 5 Gyr model lies within 0.05 dex of log10(5e9 − 2τ).
  - Rewrite `test_light_weighted_age_of_a_single_populated_bin_lies_in_that_bin`, and any other test the inventory marks.
- [ ] **Step 2: Run the tests and see them fail.**
- [ ] **Step 3: Implement.** Follow `inventory-fits-reference.md` B–D and G line by line. The spectral and photometric likelihoods stay unchanged. Remove the xfail marks added in Task 1.
- [ ] **Step 4: Run the whole suite** and confirm everything passes except tests in `ml_*` files (Task 3).
- [ ] **Step 5: Commit.** Message: `feat: fits, calibration mocks and validation rulers under the delayed-tau model (sfh_model recorded in every fit)`.

### Task 3: ML data, networks, training CLI and evaluation

**Files:**
- Modify:
  - `Synthetic_domain/ml_data.py`, `ml_networks.py`, `train_ml_model.py`, `ml_evaluation.py` (`inventory-ml.md` sections 1–4 and 7);
  - `Synthetic_domain/remote/submit_training.sh` and `remote/train_job.sh` (default `EMBEDDINGS="mlp"`).
- Test: `Synthetic_domain/tests/test_ml_data.py`, `test_ml_networks.py`, `test_train_ml_model.py`, `test_ml_evaluation.py`.

**Interfaces:**
- **Consumes:** `sm.TRAINING_LABELS`, `sm.delayed_tau_labels`, `sm.maximum_age`, `sm.AGE_MIN`, `sm.Z_RANGE`.
- **Produces:**
  - `ml_data.TARGET_NAMES` (the 9 targets in Global Constraints order) and `N_TARGETS = 9`;
  - `encode_targets(table) -> ndarray (n, 9)`: logs of age, tau, metallicity and tau_bc; the other targets as they are;
  - `decode_samples(theta) -> dict label -> ndarray`, with every `TRAINING_LABELS` key derived through `sm.delayed_tau_labels`. Redshift is clipped to `Z_RANGE`, then age to [AGE_MIN, maximum_age(z)] and τ and Z to their ranges. Clipping is a no-op for in-prior values.
  - `ml_networks.save_run` writes `config["target_names"]`. `load_run` raises `ValueError("run was trained for targets ...")` when they differ from `TARGET_NAMES`.

- [ ] **Step 1: Write failing tests.**
  - The encode → decode round trip is exact on 200 simulated rows (`sm.simulate_models` on a small LHS), for every `TRAINING_LABELS` column, to 1e-9.
  - Decoding a theta whose age exceeds t_max(z) clips to t_max(z).
  - `load_run` on a config without the right `target_names` raises.
  - The existing colour, noise and redshift tests keep passing.
  - Mark `test_validation_observations_modes` and `test_training_cli_writes_a_run` with `pytest.mark.skipif(<parquet lacks the 'age' column>, reason="synthetic sets not yet regenerated on the PC")`.
- [ ] **Step 2: Run the tests and see them fail.**
- [ ] **Step 3: Implement** per `inventory-ml.md`.
  - Remove `fraction_columns`, the bin log-ratio encoding, `mass_formed_younger_than` and every `N_BINS` use.
  - `ml_evaluation` keeps its generic metrics; delete or adapt only what the inventory marks as bin-specific.
- [ ] **Step 4: Run the full suite and see it pass.**
- [ ] **Step 5: Commit.** Message: `feat: ML targets are the delayed-tau parameters; labels derived analytically; checkpoints record target names`.

### Task 4: `training_set.ipynb`

**Files:** Modify `Synthetic_domain/training_set.ipynb` (`inventory-prior-labels.md` C).

- [ ] **Step 1:** Edit cells programmatically, keeping the cell ids. Do not execute anything.
  - **Cell 0 (rewrite):** the delayed-τ SFH, the prior table copied from the spec, the label list and the outputs.
  - **Cell 1:** `N_JOBS = int(os.environ.get("N_JOBS", 7))`. A PC job gets 8 cores, so 7 workers fit.
  - **Cell 2:** the cube has 7 dimensions, `len(sm.DIMS)`.
  - **Cell 3:** the display columns become `redshift, age, tau, metallicity, log_mstar, log_sfr, mw_age, f_old`.
  - **Every other cell:** keep, but fix any text that mentions bins.
  - Validate with `nbformat.validate` and compile every code cell.
- [ ] **Step 2: Commit.** Message: `feat: notebook 1 generates the delayed-tau training and synthetic sets`.

### Task 5: `spectroscopic_reference.ipynb`, the skip list and the literature for 168 galaxies

**Files:**
- Modify:
  - `Synthetic_domain/spectroscopic_reference.ipynb` (`inventory-fits-reference.md` F);
  - `Synthetic_domain/skipped_fit_jobs.txt` (keep the header comment and remove the ids);
  - `Synthetic_domain/fetch_literature.py`, only if needed for 168 galaxies (inventory H).

- [ ] **Step 1: Superseding (cell 10).**
  - Move every `.h5` file and MultiNest folder in `spectral_fits/` and `photometric_fits/` whose metadata fails `sm.fit_is_current` to `<folder>/superseded/`.
  - Do the same for `spectral_fits/calibration_mocks.joblib` and `photometric_fits/calibration_mocks.joblib` when their truths lack `log_age`, and for `photometric_fits/jobs`/`spectral_fits/jobs` entries of old calibration mocks.
  - Existing names in `superseded/` get a `_16bin` suffix rather than being skipped.
- [ ] **Step 2: Labels.**
  - Cell 14 and the photometric calibration cell use `sm.CALIBRATION_LABELS`.
  - Cell 16 `LABEL_COLUMNS` = the `TRAINING_LABELS` summaries (`<label>_p16/_median/_p84`).
  - Section 6b writes `phot_<label>_*` for `CALIBRATION_LABELS + ["mw_age", "log_sfr"]`.
  - The validation cells use `log_metallicity` where they used `Z_fin`.
  - Section 0's markdown describes the delayed-τ model.
- [ ] **Step 3: Literature.**
  - Run `fetch_literature.py` locally, but only if the inventory shows it is light: network queries only, no credentials, under a few minutes. Otherwise make it a step of the notebook-2 job.
  - It must cover all 168 selected galaxies. Commit the refreshed `literature/` tables.
- [ ] **Step 4:** Validate and compile both notebooks. Run the unit tests.
- [ ] **Step 5: Commit.** Message: `feat: notebook 2 supersedes 16-bin fits and reads delayed-tau fits; literature for all 168 galaxies`.

### Task 6: `ml_model.ipynb`

**Files:** Modify `Synthetic_domain/ml_model.ipynb` (`inventory-ml.md` 5 and 7).

- [ ] **Step 1:** Edit the notebook. Do not execute it.
  - `RUN_DIRECTORIES = {"mlp": "ml_runs/mlp"}`. Model choice becomes "MLP (the CNN lost in the 16-bin comparison)".
  - `KEY_LABELS` = the 9 targets plus `mw_age`, `f_young`, `f_mid`, `f_old`, `log_sfr`.
  - `REFERENCE_LABELS = sm.CALIBRATION_LABELS + ["mw_age"]`.
  - The information-per-bin figure becomes information per parameter (width ratio per target).
  - Remove the A → B → C text and prints. The agreement section prints both verdicts (photometric Bagpipes = emulator check; spectral Bagpipes = physical check) and no next-step decision.
  - The catalogue writes the new labels.
  - Validate and compile.
- [ ] **Step 2: Commit.** Message: `feat: ML notebook evaluates the delayed-tau MLP`.

### Task 7: Docs

**Files:** Modify `README.md` (`inventory-tests-scripts-docs.md` C) and `Synthetic_domain/remote/README.md`.

- [ ] **Step 1:** Update the model description, the labels, the targets, the design-choice table (add "delayed-τ replaces 16 bins, why"), the status (results pending the PC runs) and the file map.
- [ ] **Step 2: Commit.** Message: `docs: README describes the delayed-tau model`.

### Operations (controller, after review)

1. Move the old ML runs (`ml_runs/mlp`, `ml_runs/cnn`) to `ml_runs/superseded/`.
2. **PC, CPU:** run `training_set.ipynb` as a notebook job (inputs: `dr3_corrected.csv`, `splus_noise_model.csv`, `splus-filters`). Fetch its outputs. They are large; the notebook job brings back what it writes, `training_set.parquet` included, so check the download size or split the job.
3. **Laptop, before every notebook-2 job (required; file moves only, seconds):** supersede the old fits locally.
   ```bash
   cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -c "import synthetic_model as sm; [print(d, len(sm.supersede_outdated_fits(d))) for d in ('spectral_fits', 'photometric_fits')]"
   ```
   - Why: cell 10 also supersedes, but in the job it moves the job's own copies. `notebook_job.sh` returns only files newer than its start marker, and a move keeps the old mtime, so the moves never reach the laptop.
   - Without this step, the local 16-bin `.h5` files stay at the top level. `submit_fits.sh` would bundle them as `finished_fits/`, and `fit_and_save` would report every such id as failed. Also, the fetched new `calibration_mocks.joblib` would overwrite the local 16-bin one, which is the only copy (`spectral_fits/` is git-ignored).
   - Check: the first call moves about 233 items in `spectral_fits` and 155 in `photometric_fits` (dry run of Task 5); a repeat moves 0. Afterwards no 16-bin `.h5` (except an unreadable one, which `fit_and_save` refits) and no 16-bin `calibration_mocks.joblib` remain at the top level of either folder; the current delayed-tau fits and mocks stay there on later calls.
   - Bundle: `submit_notebook.sh` leaves out every `superseded/` and `pilot/` folder (`EXCLUDE_NAMES`), which no notebook reads; check the size with `BUNDLE_ONLY=<folder>` before submitting.
4. **PC, CPU:** run notebook 2, which exports the new spectral and photometric jobs (its own superseding then finds nothing to move). Then submit the spectral fits with `remote/submit_fits.sh` (14 workers).
5. **PC, GPU:** train the MLP with `remote/submit_training.sh` (`EMBEDDINGS=mlp`).
6. **PC, CPU:** first run notebook 2 again (after step 3) to read the spectral fits and export the photometric jobs; `photometric_fits/jobs/` must be empty before it, so it then holds only this export (the 87 jobs of the 16-bin split were moved to `photometric_fits/superseded/jobs/`). Then run the photometric fits with `FITS_SUBDIR=photometric_fits remote/submit_fits.sh`.
7. **PC:** run notebook 2 to read every fit (after step 3 again, with `photometric_fits` added to its inputs), then `ml_model.ipynb` on the GPU. Report.
