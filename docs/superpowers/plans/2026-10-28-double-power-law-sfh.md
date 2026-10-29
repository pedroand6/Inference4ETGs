# Double-Power-Law SFH Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Bagpipes' double-power-law (`dblplaw`) SFH as a second model next to delayed-τ. DPL gets its own:
- training and synthetic sets;
- photometric Bagpipes fits, calibration mocks and cross-fits;
- MLP flow.

The ML notebook then compares the two models: ML against its own Bagpipes, model against model, both against the spectra, and the Bayesian evidence.

**Architecture:**
- A new module `sfh_models.py` holds the shared prior constants and one section per SFH model. A small `SfhModel` registry (`MODELS`) gives each model its dims, components, labels, flow targets, calibration labels and file suffix.
- Every model-dependent function takes the model by name (`sfh_model`). Jobs, fits and checkpoints store that name, so files of the two models cannot be mixed up.
- Delayed-τ keeps today's file names and results. DPL files get the suffix `_dblplaw`.
- The notebooks loop over `SFH_MODELS`. They run only as `hljob` jobs on the PC.

**Tech Stack:** Python 3.10, Bagpipes 1.3.5 (bc03_miles, `dblplaw` and `delayed` components), PyMultiNest, scipy, joblib, torch 2.5.1 + zuko 1.2.0, pytest, `hljob`.

**Spec:** `docs/superpowers/specs/2026-10-28-double-power-law-sfh-design.md`, including its amendment (2026-10-28, from the plan review). It extends `docs/superpowers/specs/2026-10-26-delayed-tau-sfh-design.md`.

## Global Constraints

- **Code style (user rule):**
  - legible, didactic, not prolix, English names;
  - numbered `# ====` section banners;
  - docstrings with "Why:" when not obvious, "Inputs:" with units in brackets, and "Returns:";
  - named constants with `# [unit] meaning` comments;
  - logic in modules, orchestration in notebooks;
  - lines ≤ ~110 characters.
- **DPL prior (verbatim from the spec):**
  - `redshift` 0.002–0.7, uniform;
  - `turnover_lookback` 0.05 Gyr to t_max(z) = age of the Universe at z − 0.1 Gyr, log-uniform, drawn as a unit-interval fraction of the log range;
  - τ = age of the Universe at z − `turnover_lookback`;
  - `alpha` and `beta` 0.01–1000, log-uniform;
  - `metallicity`, `av_ism`, `tau_bc`, `logU` as delayed-τ;
  - `log_mass_formed` fixed at 10, then rescaled to the observed r magnitude.
- **Shared labels:** `log_mstar`, `mw_age`, `log_metallicity`, `av_ism`, `log_ssfr`, `f_young`, `f_mid`, `f_old`.
- **Delayed-τ results are kept, not recomputed:**
  - sets, fits, mocks, `ml_runs/mlp` and `ml_sfh_catalogue.parquet` keep their names and contents;
  - the refactored delayed-τ path must reproduce the stored sets exactly.
- **Existing tests:** they change only through the new `sfh_model` argument and the new metadata keys. Their delayed-τ assertions stay as they are.
- **Spectral fits stay delayed-τ only.** There are no DPL spectral fits or mocks.
- **Compute:** the laptop runs only code changes, unit tests and quick checks. Notebooks, fits and training are `hljob` jobs on the PC (bazzite):
  - one job at a time;
  - CPU jobs: `CORES=16`, `WORKERS=15`;
  - GPU jobs: `MEMORY_MB=7100`;
  - CPU notebook jobs run with `--pause-llm` when they need more than 7 GB.
- **Git:**
  - commit locally on branch `home`, and never merge or push;
  - every commit message ends with a blank line and then `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Never touch:**
  - `Synthetic_domain/pipes/`, `.env` or its values, the hljob token, or the repository-root PDF;
  - the PC's services (e.g. the hljobd memory limit) without asking.
- **Tests:** run from `Synthetic_domain`, as `../.venv/bin/python -m pytest tests/<file>.py -q`. The full suite is `../.venv/bin/python -m pytest -q`.

## Review Focus

These are the inputs and failure modes the spec implies but no single test guards. Each line names the task whose tests pin it.
1. **Mixing the two models' files.** A DPL fit, job or checkpoint read as delayed-τ, or the reverse, must raise or be skipped, never be silently used. This also covers a DPL fit landing in `photometric_fits/`, and a cross-fit counted as a calibration mock. Pinned by: Task 4 (`fit_is_current` per folder; the supersede guard); Task 8 (`load_run` guard); Task 5 (cross-fit ids distinct from `calib_`).
2. **Extreme DPL corners** (α = β = 1000 with the turnover at t_max, at z = 0.002 and z = 0.7; α = β = 0.01). Labels must stay finite, with no 0/0, and the fractions must sum to 1. Pinned by: Task 2 (corner tests); Task 7 (`test_decode_extreme_dblplaw_samples_stay_finite`).
3. **Flow samples outside the prior.** Decoding must clip each sample to the model's prior at the sample's own redshift. In-prior values must decode exactly. Pinned by: Task 7 (clipping and round-trip tests for both models).
4. **Delayed-τ silently changing in the refactor.** Covers draws, labels, file names and the superseding of existing delayed-τ fits. Pinned by: Task 1 (regression against the stored parquet rows); Task 4 (the delayed folders supersede exactly as before).
5. **Missing data in the comparison.** Covers galaxies without a converged fit of one model, labels absent from a reference, and test-half columns absent when the locked test was not exported. The comparison must report counts and skip, not crash or pass on empty rows. Pinned by: Task 9 (helpers on toy tables with NaNs); Task 10 (`available_labels`, and the converged-fit filters).

---
## Controller notes for executors

- **Expected test totals:** they were counted on dry runs of the plan's own code blocks. The binding gate is that every test passes, with skips only for the tests that skip themselves (the short MultiNest run in Task 4). If a total differs by a few tests from the number given, but every test passes, the step passes.
- **Line numbers:** they are approximate. Find each edit place by its content, as the steps say.
- **Bundle scripts:** Task 1 adds `sfh_models.py` to the bundles of `remote/submit_fits.sh` and `remote/submit_training.sh`. Later tasks edit those files on top of Task 1's version.
- **DPL decoding:** `ml_data.decode_samples(..., n_jobs)` decodes in threads (`joblib`, `prefer="threads"`). `ml_model.ipynb` keeps only one model's synthetic posteriors in memory at a time, because GPU jobs are capped at about 7.1 GB.

---

### Task 1: `sfh_models.py` with the shared prior and the delayed-τ model; `synthetic_model.py` re-exports it

**Files:**
- Create: `Synthetic_domain/sfh_models.py` (sections 1, 2 and 4; section 3 is left for Task 2).
- Modify: `Synthetic_domain/synthetic_model.py`:
  - imports, lines 19–24;
  - section 1, lines 29–115 (banner `# 1. Prior (fixed before looking at the data)` through the end of `components()`);
  - section 2 label lists, lines 152–156;
  - section 2 analytic functions, lines 179–233 (`delayed_mass_fraction_younger_than`, `delayed_mass_weighted_age`, `delayed_tau_labels`);
  - section 5, line 562 (`GYR_TO_YEARS`);
  - section 6, lines 816–825 (`SFH_MODEL`, `PHOTOMETRIC_DIMS`);
  - section 7, line 1375 (`CALIBRATION_LABELS`).
- Modify: `Synthetic_domain/remote/submit_fits.sh` line 58 and `Synthetic_domain/remote/submit_training.sh` lines 39–41 (the PC bundles must carry `sfh_models.py`, or every remote fit and training job fails on `import sfh_models`); `Synthetic_domain/slurm_fit.sh` line 3 (comment).
- Test: create `Synthetic_domain/tests/test_sfh_models.py`. `tests/test_synthetic_model.py`, `tests/test_spectral_validation.py` and `tests/test_ml_data.py` must pass unchanged.

**Interfaces:**
- Consumes: nothing new.
- Produces (`sfh_models`):
  - constants `COSMOLOGY, UNIVERSE_AGE_MARGIN, LABEL_EDGES, SFR_TIMESCALE, GYR_TO_YEARS, Z_RANGE, METALLICITY_RANGE, AV_ISM_RANGE, LOG_TAU_BC_RANGE, LOGU_RANGE, DUST_SLOPE, BIRTH_CLOUD_AGE, LOG_MASS_FORMED, TAU_TO_AV, SFR_FLOOR`, `AGE_MIN`, `TAU_RANGE`, `DELAYED_DIMS`;
  - `HISTORY_LABELS = ["mw_age", "f_young", "f_mid", "f_old", "log_sfr", "log_ssfr"]` (the keys every `analytic_labels` returns);
  - `SHARED_LABELS = ["log_mstar", "mw_age", "log_metallicity", "av_ism", "log_ssfr", "f_young", "f_mid", "f_old"]`;
  - `uniform(draw, low, high)`, `maximum_age(redshift)`;
  - `dust_and_nebular(params) -> (bagpipes_parts {"t_bc", "dust", "nebular"}, drawn {"av_ism", "tau_bc", "logU"})`;
  - `delayed_components(u, fixed_redshift=None) -> (bagpipes_components, drawn)`;
  - `delayed_mass_fraction_younger_than(age, tau, lookback)`, `delayed_mass_weighted_age(age, tau)`, `delayed_tau_labels(age, tau, log_mass_formed, log_mstar) -> dict`;
  - `delayed_analytic_labels(drawn_columns, log_mstar) -> dict` (the registry's `analytic_labels` for delayed-τ);
  - `SfhModel` (frozen dataclass, fields as in the contract; properties `training_labels`, `photometric_dims`), `DELAYED`, `MODELS = {"delayed": DELAYED}`;
  - `model_file(name, sfh_model) -> str`.
- Produces (`synthetic_model`): re-exports of all the above that callers use (`sm.Z_RANGE`, `sm.COSMOLOGY`, `sm.maximum_age`, `sm.delayed_tau_labels`, `sm.GYR_TO_YEARS`, `sm.MODELS`, `sm.model_file`, `sm.HISTORY_LABELS`, `sm.SHARED_LABELS`, ...). The delayed-only globals `DIMS, PHOTOMETRIC_DIMS, DRAWN_LABELS, DERIVED_LABELS, TRAINING_LABELS, CALIBRATION_LABELS, SFH_MODEL, components` become aliases of `MODELS["delayed"]`, each with the comment `# delayed-tau alias, removed in Task 11`.

- [ ] **Step 1: Write the failing test file.** Create `Synthetic_domain/tests/test_sfh_models.py`:

```python
"""Tests for sfh_models.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest tests/test_sfh_models.py -q"""

import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy.stats import qmc

import sfh_models as sfh
import synthetic_model as sm

# Sizes and seeds of the stored delayed-tau sets, as training_set.ipynb draws them
TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
N_TRAINING_MODELS = 2**19
N_HELD_OUT_MODELS = 2**14
TRAINING_SEED = 0
VALIDATION_SEED = 1
N_REGRESSION_ROWS = 4            # first rows of each stored set that are simulated again (Bagpipes: ~0.3 s each)
REGRESSION_RTOL = 1e-10          # "exactly", up to floating-point rounding
REGRESSION_ATOL = 1e-12


def parquet_has_column(path, column):
    """True when the parquet file exists and its schema holds the column (read without loading the data)."""
    try:
        return column in pq.read_schema(path).names
    except OSError:
        return False


def first_unit_vectors(seed, n_models):
    """First N_REGRESSION_ROWS Latin hypercube vectors of a delayed-tau set, drawn as training_set.ipynb does.

    Why the full draw: the Latin hypercube permutes all n_models strata, so its first rows depend on n_models.
    Inputs: seed: Latin hypercube seed; n_models: size of the stored set (count).
    Returns: array (N_REGRESSION_ROWS, len(DELAYED_DIMS)) of unit-cube values (dimensionless).
    """
    sampler = qmc.LatinHypercube(len(sfh.DELAYED_DIMS), seed=seed)
    return sampler.random(n_models)[:N_REGRESSION_ROWS]


# =============================================================================
# Delayed-tau unchanged by the move to sfh_models.py
# =============================================================================
@pytest.mark.skipif(not parquet_has_column(TRAINING_FILE, "age"), reason="no delayed-tau training set here")
def test_delayed_path_regenerates_the_first_rows_of_the_stored_training_set():
    simulated = sm.simulate_models(first_unit_vectors(TRAINING_SEED, N_TRAINING_MODELS))
    stored = pd.read_parquet(TRAINING_FILE).head(N_REGRESSION_ROWS)
    columns = sfh.DELAYED.training_labels + sm.BANDS

    assert list(stored.columns) == columns
    np.testing.assert_allclose(simulated[columns].to_numpy(), stored[columns].to_numpy(),
                               rtol=REGRESSION_RTOL, atol=REGRESSION_ATOL)


@pytest.mark.skipif(not parquet_has_column(VALIDATION_FILE, "age"), reason="no delayed-tau validation set here")
def test_delayed_path_regenerates_the_labels_of_the_stored_validation_set():
    """The stored labels are rescaled to an observed brightness: undo that with the stored magnitude offset."""
    simulated = sm.simulate_models(first_unit_vectors(VALIDATION_SEED, N_HELD_OUT_MODELS))
    stored = pd.read_parquet(VALIDATION_FILE).head(N_REGRESSION_ROWS)
    offsets = stored["magnitude_offset"].to_numpy()
    assert np.isfinite(offsets).all()
    unscaled = sm.rescale_mass_labels(stored, -offsets)
    columns = sfh.DELAYED.training_labels

    np.testing.assert_allclose(simulated[columns].to_numpy(), unscaled[columns].to_numpy(),
                               rtol=REGRESSION_RTOL, atol=REGRESSION_ATOL)


def test_synthetic_model_reexports_the_shared_prior_and_the_delayed_aliases():
    assert sm.Z_RANGE == sfh.Z_RANGE and sm.COSMOLOGY is sfh.COSMOLOGY
    assert sm.maximum_age is sfh.maximum_age and sm.delayed_tau_labels is sfh.delayed_tau_labels
    assert sm.GYR_TO_YEARS == sfh.GYR_TO_YEARS == 1e9
    delayed = sm.MODELS["delayed"]
    assert sm.DIMS == delayed.dims
    assert sm.PHOTOMETRIC_DIMS == delayed.photometric_dims
    assert sm.DRAWN_LABELS == delayed.drawn_labels
    assert sm.DERIVED_LABELS == delayed.derived_labels
    assert sm.TRAINING_LABELS == delayed.training_labels
    assert sm.CALIBRATION_LABELS == delayed.calibration_labels
    assert sm.SFH_MODEL == "delayed"
    assert sm.components is sfh.delayed_components


def test_sfh_models_does_not_import_synthetic_model():
    """synthetic_model imports sfh_models; the reverse would be a circular import."""
    command = "import sys, sfh_models; assert 'synthetic_model' not in sys.modules"
    subprocess.run([sys.executable, "-c", command], check=True)


# =============================================================================
# Shared prior and the registry
# =============================================================================
def test_dust_and_nebular_follow_the_shared_prior():
    params = {"av_ism": 0.5, "log_tau_bc": 0.5, "logU": 0.5}
    bagpipes_parts, drawn = sfh.dust_and_nebular(params)

    assert list(bagpipes_parts) == ["t_bc", "dust", "nebular"]
    assert np.isclose(drawn["av_ism"], np.mean(sfh.AV_ISM_RANGE))
    assert np.isclose(np.log10(drawn["tau_bc"]), np.mean(sfh.LOG_TAU_BC_RANGE))
    assert np.isclose(drawn["logU"], np.mean(sfh.LOGU_RANGE))
    dust = bagpipes_parts["dust"]
    assert dust["type"] == "CF00" and dust["n"] == sfh.DUST_SLOPE and dust["Av"] == drawn["av_ism"]
    assert np.isclose(dust["eta"], 1 + sfh.TAU_TO_AV * drawn["tau_bc"] / drawn["av_ism"])
    assert bagpipes_parts["nebular"] == {"logU": drawn["logU"]}
    assert bagpipes_parts["t_bc"] == sfh.BIRTH_CLOUD_AGE


@pytest.mark.parametrize("name", list(sfh.MODELS))
def test_registry_entries_are_consistent(name):
    model = sfh.MODELS[name]
    u = np.full(len(model.dims), 0.5)
    _, drawn = model.components(u)

    assert model.name == name
    assert model.dims[0] == "z"
    assert model.photometric_dims == model.dims[1:]
    assert list(drawn) == model.drawn_labels
    assert model.training_labels == model.drawn_labels + model.derived_labels
    for log_label, drawn_label in model.log_labels.items():
        assert log_label in model.derived_labels and drawn_label in model.drawn_labels
    assert set(sfh.SHARED_LABELS) <= set(model.training_labels)
    assert set(sfh.HISTORY_LABELS) <= set(model.derived_labels)
    assert model.target_names[-1] == "redshift"
    assert set(model.calibration_labels) <= set(model.training_labels)


def test_model_file_keeps_the_delayed_tau_names():
    for name in ("training_set.parquet", "photometric_fits", "ml_runs/mlp", "ml_sfh_catalogue.parquet"):
        assert sfh.model_file(name, "delayed") == name
```

  The two regression tests reproduce what `training_set.ipynb` did on the PC: `qmc.LatinHypercube(7, seed=0).random(2**19)` for the training set and `seed=1`, `2**14` for the validation set, split into blocks whose order `pd.concat` keeps. Only the first 4 vectors are simulated. The validation labels were rescaled to an observed brightness with `sm.rescale_mass_labels(models, offsets)`; the test undoes that with the stored `magnitude_offset` (the validation magnitudes carry noise and are not compared). Before the refactor these comparisons are exact (maximum difference 0.0).

- [ ] **Step 2: Run it and see it fail.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_sfh_models.py -q
```

  Expected: 1 error during collection, `ModuleNotFoundError: No module named 'sfh_models'`.

- [ ] **Step 3: Create `Synthetic_domain/sfh_models.py`.** Sections 1, 2 and 4 (section 3 comes in Task 2). The section-2 functions are moved verbatim from `synthetic_model.py`; `delayed_components` is today's `components` with the dust and nebular block moved into `dust_and_nebular` (same arithmetic, same dict key order, so the draws are bit-identical).

```python
"""
SFH models of the S-PLUS synthetic domain: the shared prior, one section per Bagpipes SFH, and the registry.

synthetic_model.py imports this module (never the reverse), so the training sets, the fits and the ML use
exactly the same prior and labels for each SFH model.
"""

import os
from dataclasses import dataclass
from typing import Callable

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.special import gammainc, gammaincc


# =============================================================================
# 1. Shared prior (fixed before looking at the data)
# =============================================================================
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3)
UNIVERSE_AGE_MARGIN = 0.1        # [Gyr] age upper bound = Universe age at z minus this
LABEL_EDGES = [0.1, 1.0]         # [Gyr] f_young < 0.1 < f_mid < 1.0 < f_old
SFR_TIMESCALE = 0.1              # [Gyr] log_sfr is the mean SFR over this lookback time, as Bagpipes' sfr
GYR_TO_YEARS = 1e9               # [yr/Gyr]

# Parameters every SFH model shares: redshift, one metallicity, dust, nebular emission and the formed mass
Z_RANGE = (0.002, 0.7)                               # redshift
METALLICITY_RANGE = (0.005, 2.5)                     # [Z/Z_sun] the bc03_miles grid range
AV_ISM_RANGE = (0.01, 2.0)                           # [mag]
LOG_TAU_BC_RANGE = (np.log10(0.5), np.log10(4.0))    # birth-cloud optical depth
LOGU_RANGE = (-4.0, -1.0)                            # nebular ionisation parameter
DUST_SLOPE = 0.7                                     # CF00 attenuation-law exponent
BIRTH_CLOUD_AGE = 0.01                               # [Gyr]
LOG_MASS_FORMED = 10.0                               # [log M_sun]
TAU_TO_AV = 1.086                                    # A_V = 1.086 * tau_V
SFR_FLOOR = 1e-14                # [yr^-1] sSFR floor for SFHs with no recent star formation

# Labels every model's analytic_labels returns, computed from the shape of its SFH
HISTORY_LABELS = ["mw_age", "f_young", "f_mid", "f_old", "log_sfr", "log_ssfr"]


def uniform(draw, low, high):
    """Map a unit-interval draw onto [low, high).

    Inputs: draw in [0, 1) (dimensionless); low, high: bounds in the units of the parameter.
    Returns: parameter value in [low, high), in the same units as low and high.
    """
    return low + (high - low) * draw


def maximum_age(redshift):
    """Oldest allowed age [Gyr]: the age of the Universe at the redshift minus UNIVERSE_AGE_MARGIN.

    Why the margin: Bagpipes flags a model forming stars before the Big Bang as unphysical (NaN magnitudes).
    Inputs: redshift: dimensionless z, scalar or array.
    Returns: maximum age [Gyr], a float for a scalar redshift or an array of the same shape.
    """
    return COSMOLOGY.age(redshift).value - UNIVERSE_AGE_MARGIN


def dust_and_nebular(params):
    """Bagpipes dust, birth-cloud and nebular entries of one model, the same for every SFH model.

    Inputs: params: dict of unit-cube values (dimensionless, in [0, 1)) holding "av_ism", "log_tau_bc" and "logU".
    Returns: (bagpipes_parts, drawn). bagpipes_parts holds "t_bc" [Gyr], "dust" (CF00) and "nebular";
             drawn holds av_ism [mag], tau_bc and logU (dimensionless).
    """
    av_ism = uniform(params["av_ism"], *AV_ISM_RANGE)
    tau_birth_cloud = 10 ** uniform(params["log_tau_bc"], *LOG_TAU_BC_RANGE)
    av_birth_cloud = TAU_TO_AV * tau_birth_cloud
    log_u = uniform(params["logU"], *LOGU_RANGE)

    bagpipes_parts = {
        "t_bc": BIRTH_CLOUD_AGE,
        "dust": {"type": "CF00", "Av": av_ism, "eta": 1 + av_birth_cloud / av_ism, "n": DUST_SLOPE},
        "nebular": {"logU": log_u},
    }
    drawn = {"av_ism": av_ism, "tau_bc": tau_birth_cloud, "logU": log_u}
    return bagpipes_parts, drawn


# =============================================================================
# 2. Delayed-tau SFH: SFR(t) ~ t exp(-t / tau) since star formation began
# =============================================================================
# One Bagpipes "delayed" component with a single metallicity. Ages and tau are drawn log-uniform,
# because they span several decades.
AGE_MIN = 0.05                   # [Gyr] youngest age (time since star formation began)
TAU_RANGE = (0.01, 20.0)         # [Gyr] e-folding time

# "log_age_fraction" places the age inside [log AGE_MIN, log t_max(z)], so one cube works at every z
DELAYED_DIMS = ["z", "log_age_fraction", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU"]


def delayed_components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one delayed-tau unit-cube vector u (order of DELAYED_DIMS).

    Inputs: u: unit-cube vector (dimensionless, in [0, 1)) ordered as DELAYED_DIMS;
            fixed_redshift: if given, replaces the redshift drawn from u (fits); dimensionless z.
    Returns: (bagpipes_components, drawn). drawn holds the delayed-tau drawn labels: redshift (dimensionless),
             age and tau [Gyr], metallicity [Z/Z_sun], av_ism [mag], tau_bc and logU (dimensionless),
             log_mass_formed [log M_sun].
    """
    params = dict(zip(DELAYED_DIMS, u))

    # 1. Redshift, which sets the oldest allowed age
    if fixed_redshift is None:
        redshift = uniform(params["z"], *Z_RANGE)
    else:
        redshift = fixed_redshift
    oldest_age = maximum_age(redshift)

    # 2. Star formation history and metallicity
    age = 10 ** uniform(params["log_age_fraction"], np.log10(AGE_MIN), np.log10(oldest_age))
    # Rounding in 10**log10 can push the top draw a hair above the bound
    age = min(age, oldest_age)
    tau = 10 ** uniform(params["log_tau"], *np.log10(TAU_RANGE))
    metallicity = 10 ** uniform(params["log_metallicity"], *np.log10(METALLICITY_RANGE))

    # 3. Dust and nebular emission
    dust_parts, dust_drawn = dust_and_nebular(params)

    bagpipes_components = {"redshift": redshift}
    bagpipes_components.update(dust_parts)
    bagpipes_components["delayed"] = {"age": age, "tau": tau, "massformed": LOG_MASS_FORMED,
                                      "metallicity": metallicity}
    drawn = {"redshift": redshift, "age": age, "tau": tau, "metallicity": metallicity,
             "av_ism": dust_drawn["av_ism"], "tau_bc": dust_drawn["tau_bc"], "logU": dust_drawn["logU"],
             "log_mass_formed": LOG_MASS_FORMED}
    return bagpipes_components, drawn


def delayed_mass_fraction_younger_than(age, tau, lookback):
    """Share of the formed mass with lookback time below `lookback` for SFR(t) ~ t exp(-t / tau).

    t runs from the start of star formation, `age` ago. With P and Q the regularised lower and upper
    incomplete gamma functions, the mass formed up to t is proportional to P(2, t / tau).
    Why the Q form: Q(2, (A - L) / tau) - Q(2, A / tau) avoids the cancellation of 1 - P / P when tau << A.
    Inputs: age A, tau, lookback L [Gyr], scalars or arrays that broadcast together.
    Returns: fraction in [0, 1] (dimensionless), an array of the broadcast shape; 1 where L >= A.
    """
    age, tau, lookback = np.broadcast_arrays(np.asarray(age, dtype=float), np.asarray(tau, dtype=float),
                                             np.asarray(lookback, dtype=float))
    formed_before_now = gammainc(2, age / tau)
    older_start = np.clip(age - lookback, 0, None) / tau
    fraction = (gammaincc(2, older_start) - gammaincc(2, age / tau)) / formed_before_now
    return np.where(lookback >= age, 1.0, np.clip(fraction, 0.0, 1.0))


def delayed_mass_weighted_age(age, tau):
    """Formed-mass-weighted lookback age of SFR(t) ~ t exp(-t / tau): A - 2 tau P(3, A / tau) / P(2, A / tau).

    Limits: A - 2 tau when tau << A, and A / 3 when tau >> A.
    Inputs: age A and tau [Gyr], scalars or arrays that broadcast together.
    Returns: mass-weighted age [Gyr ago], an array of the broadcast shape.
    """
    age = np.asarray(age, dtype=float)
    tau = np.asarray(tau, dtype=float)
    return age - 2 * tau * gammainc(3, age / tau) / gammainc(2, age / tau)


def delayed_tau_labels(age, tau, log_mass_formed, log_mstar):
    """Analytic labels of a delayed-tau SFH, shared by the training sets and the ML decoding.

    log_sfr is the mean SFR over the last SFR_TIMESCALE, floored at SFR_FLOOR * M_star: about a fifth
    of the prior has no star formation that recent, and log(0) would be -inf.
    Inputs: age, tau [Gyr]; log_mass_formed, log_mstar [log M_sun]; arrays of one shape (or scalars).
    Returns: dict of arrays: mw_age [Gyr ago], f_young, f_mid, f_old (dimensionless, sum to 1; split at
             LABEL_EDGES), log_sfr [log M_sun/yr] and log_ssfr [log yr^-1].
    """
    log_mass_formed = np.asarray(log_mass_formed, dtype=float)
    log_mstar = np.asarray(log_mstar, dtype=float)
    younger_than_young_edge = delayed_mass_fraction_younger_than(age, tau, LABEL_EDGES[0])
    younger_than_old_edge = delayed_mass_fraction_younger_than(age, tau, LABEL_EDGES[1])
    recent_fraction = delayed_mass_fraction_younger_than(age, tau, SFR_TIMESCALE)

    recent_sfr = 10**log_mass_formed * recent_fraction / (SFR_TIMESCALE * GYR_TO_YEARS)
    sfr_floor = SFR_FLOOR * 10**log_mstar
    log_sfr = np.log10(np.maximum(recent_sfr, sfr_floor))
    return {
        "mw_age": delayed_mass_weighted_age(age, tau),
        "f_young": younger_than_young_edge,
        "f_mid": younger_than_old_edge - younger_than_young_edge,
        "f_old": 1 - younger_than_old_edge,
        "log_sfr": log_sfr,
        "log_ssfr": log_sfr - log_mstar,
    }


def delayed_analytic_labels(drawn_columns, log_mstar):
    """The HISTORY_LABELS of delayed-tau models from their drawn labels (the registry's analytic_labels).

    Inputs: drawn_columns: dict of arrays holding at least "age", "tau" [Gyr] and "log_mass_formed"
            [log M_sun]; log_mstar: array [log M_sun] of the same shape.
    Returns: dict of arrays, as delayed_tau_labels.
    """
    return delayed_tau_labels(drawn_columns["age"], drawn_columns["tau"], drawn_columns["log_mass_formed"],
                              log_mstar)


# =============================================================================
# 4. The model registry: everything that differs between SFH models
# =============================================================================
# Labels with the same meaning under every SFH model, compared between the models
SHARED_LABELS = ["log_mstar", "mw_age", "log_metallicity", "av_ism", "log_ssfr", "f_young", "f_mid", "f_old"]


@dataclass(frozen=True)
class SfhModel:
    """One SFH model: its unit cube, Bagpipes components, labels, flow targets and file names."""

    name: str                    # "delayed" or "dblplaw", stored in every job, fit and checkpoint
    dims: list                   # unit-cube dimensions, "z" first
    components: Callable         # components(u, fixed_redshift=None) -> (bagpipes_components, drawn)
    drawn_labels: list           # keys of drawn, in this order
    derived_labels: list         # labels computed from drawn and from Bagpipes' stellar mass
    log_labels: dict             # derived log label -> drawn label, e.g. {"log_age": "age"}
    analytic_labels: Callable    # analytic_labels(drawn_columns, log_mstar) -> dict with the HISTORY_LABELS
    target_names: list           # flow targets, redshift last
    calibration_labels: list     # labels of the calibration (coverage) test of the photometric fits
    file_suffix: str             # "" for delayed-tau (today's file names), "_dblplaw" for the double power law

    @property
    def training_labels(self):
        """Label columns of the sets and the fits: the drawn labels, then the derived ones. Returns: list."""
        return self.drawn_labels + self.derived_labels

    @property
    def photometric_dims(self):
        """Unit-cube dimensions of a fit at fixed redshift: dims without "z". Returns: list."""
        photometric = []
        for dimension in self.dims:
            if dimension != "z":
                photometric.append(dimension)
        return photometric


DELAYED = SfhModel(
    name="delayed",
    dims=DELAYED_DIMS,
    components=delayed_components,
    drawn_labels=["redshift", "age", "tau", "metallicity", "av_ism", "tau_bc", "logU", "log_mass_formed"],
    derived_labels=["log_age", "log_tau", "log_metallicity", "log_mstar", "log_sfr", "log_ssfr",
                    "mw_age", "f_young", "f_mid", "f_old"],
    log_labels={"log_age": "age", "log_tau": "tau", "log_metallicity": "metallicity"},
    analytic_labels=delayed_analytic_labels,
    target_names=["log_age", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU", "log_mass_formed",
                  "log_mstar", "redshift"],
    calibration_labels=["log_age", "log_tau", "log_metallicity", "log_mstar"],
    file_suffix="",
)

MODELS = {"delayed": DELAYED}


def model_file(name, sfh_model):
    """File or folder name of one SFH model's product: the delayed-tau name with the model's suffix.

    Inputs: name: the delayed-tau file or folder name ("training_set.parquet", "photometric_fits", "ml_runs/mlp");
            sfh_model: a name in MODELS.
    Returns: the name with the suffix before a file's extension, or after a folder name
             ("training_set_dblplaw.parquet", "photometric_fits_dblplaw", "ml_runs/mlp_dblplaw");
             unchanged for delayed-tau.
    """
    suffix = MODELS[sfh_model].file_suffix
    stem, extension = os.path.splitext(name)
    return stem + suffix + extension
```

- [ ] **Step 4: Make `synthetic_model.py` import and re-export it.**

  4a. Imports (lines 19–24). Replace

```python
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simpson
from scipy.special import gammainc, gammaincc
from scipy.stats import chi2 as chi_squared

import bagpipes as pipes
```

  with

```python
from scipy.integrate import simpson
from scipy.stats import chi2 as chi_squared

import bagpipes as pipes
# The shared prior and the SFH models live in sfh_models.py. They are re-exported here, so sm.Z_RANGE,
# sm.maximum_age, sm.delayed_tau_labels and the rest keep working
from sfh_models import (AGE_MIN, AV_ISM_RANGE, BIRTH_CLOUD_AGE, COSMOLOGY, DUST_SLOPE, GYR_TO_YEARS,
                        HISTORY_LABELS, LABEL_EDGES, LOG_MASS_FORMED, LOG_TAU_BC_RANGE, LOGU_RANGE,
                        METALLICITY_RANGE, MODELS, SFR_FLOOR, SFR_TIMESCALE, SHARED_LABELS, TAU_RANGE,
                        TAU_TO_AV, UNIVERSE_AGE_MARGIN, Z_RANGE, delayed_mass_fraction_younger_than,
                        delayed_mass_weighted_age, delayed_tau_labels, maximum_age, model_file, uniform)
```

  (`FlatLambdaCDM`, `gammainc` and `gammaincc` have no other use in `synthetic_model.py`; check with `grep -n "gammainc\|FlatLambdaCDM" Synthetic_domain/synthetic_model.py`, which must print nothing afterwards.)

  4b. Section 1. Delete everything from the line `# =============================================================================` above `# 1. Prior (fixed before looking at the data)` (line 29) through `    return bagpipes_components, drawn` at the end of `components()` (line 115), and put in its place:

```python
# =============================================================================
# 1. Prior (fixed before looking at the data): see sfh_models.py
# =============================================================================
DIMS = MODELS["delayed"].dims                     # delayed-tau alias, removed in Task 11
components = MODELS["delayed"].components         # delayed-tau alias, removed in Task 11
```

  The two blank lines before the section-2 banner stay.

  4c. Section 2 label lists (lines 153–156). Replace

```python
DRAWN_LABELS = ["redshift", "age", "tau", "metallicity", "av_ism", "tau_bc", "logU", "log_mass_formed"]
DERIVED_LABELS = ["log_age", "log_tau", "log_metallicity", "log_mstar", "log_sfr", "log_ssfr",
                  "mw_age", "f_young", "f_mid", "f_old"]
TRAINING_LABELS = DRAWN_LABELS + DERIVED_LABELS
```

  with

```python
DRAWN_LABELS = MODELS["delayed"].drawn_labels           # delayed-tau alias, removed in Task 11
DERIVED_LABELS = MODELS["delayed"].derived_labels       # delayed-tau alias, removed in Task 11
TRAINING_LABELS = MODELS["delayed"].training_labels     # delayed-tau alias, removed in Task 11
```

  4d. Delete the three functions `delayed_mass_fraction_younger_than`, `delayed_mass_weighted_age` and `delayed_tau_labels` (lines 179–233, from `def delayed_mass_fraction_younger_than(` up to the two blank lines before `def summary_labels(`). `summary_labels` and `simulate_models` stay unchanged in this task: they now use the imported `delayed_tau_labels` and the aliases.

  4e. Section 5. Delete the line `GYR_TO_YEARS = 1e9                       # [yr/Gyr]` (line 562); it is imported now.

  4f. Section 6 (lines 818–825). Replace

```python
SFH_MODEL = "delayed"

# The photometric fit (spec amendment A) samples the training prior at the spectroscopic redshift and nothing
# else: its brightness offset is marginalised analytically, and the simulator has no error multiplier
PHOTOMETRIC_DIMS = []
for dimension in DIMS:
    if dimension != "z":
        PHOTOMETRIC_DIMS.append(dimension)
```

  with

```python
SFH_MODEL = MODELS["delayed"].name                      # delayed-tau alias, removed in Task 11

# The photometric fit (spec amendment A) samples the training prior at the spectroscopic redshift and nothing
# else: its brightness offset is marginalised analytically, and the simulator has no error multiplier
PHOTOMETRIC_DIMS = MODELS["delayed"].photometric_dims   # delayed-tau alias, removed in Task 11
```

  4g. Section 7 (line 1375). Replace

```python
CALIBRATION_LABELS = ["log_age", "log_tau", "log_metallicity", "log_mstar"]
```

  with

```python
CALIBRATION_LABELS = MODELS["delayed"].calibration_labels   # delayed-tau alias, removed in Task 11
```

  `SPECTRAL_DIMS` (section 5) keeps its loop over `DIMS`; Task 4 rebuilds it from `MODELS["delayed"]`.

- [ ] **Step 5: Ship `sfh_models.py` in the PC bundles.**
  - `Synthetic_domain/remote/submit_fits.sh` line 58: replace
    `cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/run_spectral_fits.py" "$bundle_dir/"`
    with
    `cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/sfh_models.py" "$SYNTHETIC_DIR/run_spectral_fits.py" "$bundle_dir/"`.
  - `Synthetic_domain/remote/submit_training.sh` lines 39–41: replace

```bash
cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/ml_data.py" "$SYNTHETIC_DIR/ml_networks.py" \
   "$SYNTHETIC_DIR/train_ml_model.py" "$SYNTHETIC_DIR/training_set.parquet" \
   "$SYNTHETIC_DIR/synthetic_validation.parquet" "$SYNTHETIC_DIR/splus_noise_model.csv" "$bundle_dir/"
```

    with

```bash
cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/sfh_models.py" "$SYNTHETIC_DIR/ml_data.py" \
   "$SYNTHETIC_DIR/ml_networks.py" "$SYNTHETIC_DIR/train_ml_model.py" "$SYNTHETIC_DIR/training_set.parquet" \
   "$SYNTHETIC_DIR/synthetic_validation.parquet" "$SYNTHETIC_DIR/splus_noise_model.csv" "$bundle_dir/"
```

  - `Synthetic_domain/slurm_fit.sh` line 3: `run_spectral_fits.py, synthetic_model.py, splus-filters/ and` becomes `run_spectral_fits.py, synthetic_model.py, sfh_models.py, splus-filters/ and`.
  - `remote/submit_notebook.sh` already packs `./*.py`, so it needs no change.
  - Check the syntax: `bash -n Synthetic_domain/remote/submit_fits.sh && bash -n Synthetic_domain/remote/submit_training.sh` prints nothing.

- [ ] **Step 6: Run the new tests and see them pass.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_sfh_models.py -q
```

  Expected: `7 passed` (about 5 s; the two regression tests skip only if the stored parquet sets are absent).

- [ ] **Step 7: Run the whole suite.** The existing tests must pass without any edit.

```bash
cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q
```

  Expected: `166 passed` (about 4 minutes).

- [ ] **Step 8: Commit.**

```bash
git add Synthetic_domain/sfh_models.py Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_sfh_models.py \
        Synthetic_domain/remote/submit_fits.sh Synthetic_domain/remote/submit_training.sh Synthetic_domain/slurm_fit.sh
git commit -m "$(cat <<'EOF'
refactor: sfh_models.py holds the shared prior and the delayed-tau model; synthetic_model re-exports it

The delayed-tau path regenerates the first rows of the stored training and validation sets exactly.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

### Task 2: Double-power-law prior, Bagpipes' age grid and `dblplaw_labels` (`sfh_models.py` section 3)

**Files:**
- Modify: `Synthetic_domain/sfh_models.py`:
  - imports (add `import bagpipes as pipes`);
  - new section 3, inserted between the end of section 2 (`delayed_analytic_labels`) and the section-4 banner;
  - section 4: the `DBLPLAW` entry and `MODELS`.
- Test: `Synthetic_domain/tests/test_sfh_models.py` (append).

**Interfaces:**
- Consumes (Task 1): `uniform`, `maximum_age`, `dust_and_nebular`, `COSMOLOGY`, `Z_RANGE`, `METALLICITY_RANGE`, `LOG_MASS_FORMED`, `LABEL_EDGES`, `SFR_TIMESCALE`, `SFR_FLOOR`, `GYR_TO_YEARS`, `SfhModel`, `HISTORY_LABELS`.
- Produces:
  - `TURNOVER_LOOKBACK_MIN = 0.05`, `ALPHA_RANGE = (0.01, 1000.0)`, `BETA_RANGE = (0.01, 1000.0)`, `DBLPLAW_DIMS`;
  - `BAGPIPES_LOG_AGE_START = 6.0`, `BAGPIPES_LOG_AGE_STEP = 0.0025`, `LOG_GYR_TO_YEARS = 9.0`, `DBLPLAW_CHUNK = 32`;
  - `YOUNG_EDGE_INDEX`, `OLD_EDGE_INDEX`, `RECENT_EDGE_INDEX`: the first grid index at or above `LABEL_EDGES[0]`, `LABEL_EDGES[1]` and `SFR_TIMESCALE` (800, 1200 and 800);
  - `dblplaw_components(u, fixed_redshift=None) -> (bagpipes_components with "dblplaw": {"tau", "alpha", "beta", "massformed", "metallicity"}, drawn)`;
  - `bagpipes_age_grid() -> (ages, widths)` [Gyr] and the module constants `BAGPIPES_AGES`, `BAGPIPES_AGE_WIDTHS` [Gyr];
  - `bagpipes_universe_age(redshift) -> age [Gyr]` (Bagpipes' interpolated age of the Universe; it decides which grid ages lie after the Big Bang, exactly as in Bagpipes);
  - `dblplaw_grid_sums(tau, alpha, beta, redshift) -> dict` of 1-D arrays `mw_age, f_young, f_mid, f_old, recent`. It is vectorised over blocks of `DBLPLAW_CHUNK` samples, with no Python loop over samples, and each block uses only the grid ages below the age of the Universe at the block's lowest redshift;
  - `dblplaw_labels(tau, alpha, beta, redshift, log_mass_formed, log_mstar) -> dict` with the `HISTORY_LABELS` keys, arrays of the broadcast shape;
  - `dblplaw_analytic_labels(drawn_columns, log_mstar) -> dict`;
  - `DBLPLAW` (the contract's registry entry), `MODELS = {"delayed": DELAYED, "dblplaw": DBLPLAW}`.

  Bagpipes facts this task relies on (checked in `.venv/lib/python3.10/site-packages/bagpipes/models/star_formation_history.py`):
  - The grid is built in `__init__`: `hubble_time = utils.age_at_z[utils.z_array == 0.]`, `ages = 10**arange(6, log10(hubble_time) + 9 + 2 * 0.0025, 0.0025)` [yr], `age_lhs = 10**make_bins(log ages, make_rhs=True)[0]` with the first edge set to 0 and the last to `1e9 * hubble_time`, `age_widths = diff(age_lhs)`. That is 1654 points.
  - `update` uses `age_of_universe = 1e9 * np.interp(redshift, utils.z_array, utils.age_at_z)`.
  - `dblplaw` sets `sfr[ages < age_of_universe] = ((t/tau)**alpha + (t/tau)**-beta)**-1` with `t = age_of_universe - ages`, and flags `tau > 1.5 * age_of_universe` as unphysical.
  - `mass_weighted_age = sum(sfh * age_widths * ages) / sum(sfh * age_widths)`.
  - `sfr` is the mean over `ages < 1e8` divided by the summed widths there. That is 0.1–0.3% from "mass younger than 100 Myr / 100 Myr", well inside 1%.

  **Speed.** Decoding the flow's posteriors calls `dblplaw_labels` on about 1e8 samples, so the grid sum is written for speed, measured on the laptop (one core, 20 000 random prior samples):
  - An earlier draft took 0.12 ms per sample. It used `log`, `logaddexp` and `exp`, boolean-mask sums and blocks of 2048 samples. Each (block × 1654) array was then 27 MB, which does not fit in the cache.
  - The version below takes 0.06–0.07 ms per sample. It needs three transcendental calls per grid point: `log`, plus two `exp`.
  - Its blocks hold 32 samples, so each array is 0.42 MB. Blocks of 32 and 64 samples gave 0.057 and 0.066 ms per sample; 128 gave 0.08 ms.
  - It uses only the grid ages that the block's lowest redshift can reach.
  - The label sums are slices of the sorted grid, because the grid ages increase.
  - It does no BLAS products. A threaded `formed_mass @ ages` made it 2–4× slower on the laptop.

  Task 7 adds the parallel decode, `decode_samples(theta, sfh_model, n_jobs)`. With 15 workers on the PC, 3e7 samples take a few minutes.

- [ ] **Step 1: Write the failing tests.** In `Synthetic_domain/tests/test_sfh_models.py`, add `import bagpipes as pipes` to the imports so they read

```python
import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy.stats import qmc

import bagpipes as pipes
import sfh_models as sfh
import synthetic_model as sm
```

  and append at the end of the file:

```python
# =============================================================================
# Double power law: prior
# =============================================================================
LOW_REDSHIFT = 0.002             # Bagpipes needs z > 0 to build a star formation history


def test_dblplaw_draws_stay_inside_bounds_with_the_turnover_after_the_big_bang():
    rng = np.random.default_rng(0)
    for u in rng.random((500, len(sfh.DBLPLAW_DIMS))):
        bagpipes_components, drawn = sfh.dblplaw_components(u)
        redshift = drawn["redshift"]
        assert sfh.Z_RANGE[0] <= redshift <= sfh.Z_RANGE[1]
        assert sfh.TURNOVER_LOOKBACK_MIN <= drawn["turnover_lookback"] <= sfh.maximum_age(redshift)
        assert np.isclose(drawn["tau"], sfh.COSMOLOGY.age(redshift).value - drawn["turnover_lookback"])
        assert drawn["tau"] > 0
        # Bagpipes flags tau > 1.5 x the age of the Universe as unphysical
        assert drawn["tau"] < 1.5 * sfh.bagpipes_universe_age(redshift)
        assert sfh.ALPHA_RANGE[0] <= drawn["alpha"] <= sfh.ALPHA_RANGE[1]
        assert sfh.BETA_RANGE[0] <= drawn["beta"] <= sfh.BETA_RANGE[1]
        assert sfh.METALLICITY_RANGE[0] <= drawn["metallicity"] <= sfh.METALLICITY_RANGE[1]
        assert sfh.AV_ISM_RANGE[0] <= drawn["av_ism"] <= sfh.AV_ISM_RANGE[1]
        assert 10 ** sfh.LOG_TAU_BC_RANGE[0] <= drawn["tau_bc"] <= 10 ** sfh.LOG_TAU_BC_RANGE[1]
        assert sfh.LOGU_RANGE[0] <= drawn["logU"] <= sfh.LOGU_RANGE[1]
        assert drawn["log_mass_formed"] == sfh.LOG_MASS_FORMED


def test_dblplaw_turnover_at_the_top_of_the_cube_never_precedes_the_big_bang():
    cube = np.full(len(sfh.DBLPLAW_DIMS), 1.0)
    for fixed_redshift in (None, sfh.Z_RANGE[1], 0.3):
        _, drawn = sfh.dblplaw_components(cube, fixed_redshift=fixed_redshift)
        assert drawn["turnover_lookback"] <= sfh.maximum_age(drawn["redshift"])
        assert drawn["tau"] >= sfh.UNIVERSE_AGE_MARGIN - 1e-9


def test_uniform_cube_gives_log_uniform_turnover_lookback_alpha_beta_and_metallicity():
    cube = np.full(len(sfh.DBLPLAW_DIMS), 0.5)
    _, drawn = sfh.dblplaw_components(cube)
    log_lookback_range = [np.log10(sfh.TURNOVER_LOOKBACK_MIN), np.log10(sfh.maximum_age(drawn["redshift"]))]
    assert np.isclose(np.log10(drawn["turnover_lookback"]), np.mean(log_lookback_range))
    assert np.isclose(np.log10(drawn["alpha"]), np.mean(np.log10(sfh.ALPHA_RANGE)))
    assert np.isclose(np.log10(drawn["beta"]), np.mean(np.log10(sfh.BETA_RANGE)))
    assert np.isclose(np.log10(drawn["metallicity"]), np.mean(np.log10(sfh.METALLICITY_RANGE)))


@pytest.mark.parametrize("label, value_range", [("alpha", sfh.ALPHA_RANGE), ("beta", sfh.BETA_RANGE)])
def test_latin_hypercube_is_stratified_in_log_alpha_and_log_beta(label, value_range):
    """Each of n equal slices of the log range holds exactly one of n Latin hypercube draws."""
    n_draws = 50
    cube = qmc.LatinHypercube(len(sfh.DBLPLAW_DIMS), seed=4).random(n_draws)
    log_values = np.zeros(n_draws)
    for row_index, u in enumerate(cube):
        _, drawn = sfh.dblplaw_components(u)
        log_values[row_index] = np.log10(drawn[label])
    log_range = np.log10(value_range)
    slices = np.floor((log_values - log_range[0]) / np.diff(log_range) * n_draws).astype(int)
    assert np.array_equal(np.sort(slices), np.arange(n_draws))


def test_dblplaw_components_hold_one_dblplaw_component_and_the_shared_dust():
    u = np.random.default_rng(5).random(len(sfh.DBLPLAW_DIMS))
    bagpipes_components, drawn = sfh.dblplaw_components(u)

    assert list(drawn) == sfh.DBLPLAW.drawn_labels
    assert bagpipes_components["dblplaw"] == {"tau": drawn["tau"], "alpha": drawn["alpha"], "beta": drawn["beta"],
                                              "massformed": sfh.LOG_MASS_FORMED,
                                              "metallicity": drawn["metallicity"]}
    params = dict(zip(sfh.DBLPLAW_DIMS, u))
    shared_parts, _ = sfh.dust_and_nebular(params)
    for key, value in shared_parts.items():
        assert bagpipes_components[key] == value
    assert "delayed" not in bagpipes_components


def test_dblplaw_fixed_redshift_overrides_the_drawn_redshift():
    bagpipes_components, drawn = sfh.dblplaw_components(np.full(len(sfh.DBLPLAW_DIMS), 0.5), fixed_redshift=0.123)
    assert bagpipes_components["redshift"] == 0.123 and drawn["redshift"] == 0.123
    assert np.isclose(drawn["tau"], sfh.COSMOLOGY.age(0.123).value - drawn["turnover_lookback"])


@pytest.mark.parametrize("cube_value", [0.0, 0.5, 0.999999])
def test_dblplaw_prior_corners_are_physical_for_bagpipes(cube_value):
    bagpipes_components, _ = sfh.dblplaw_components(np.full(len(sfh.DBLPLAW_DIMS), cube_value))
    history = pipes.models.star_formation_history(bagpipes_components)
    assert not history.unphysical
    assert np.isfinite(history.stellar_mass)


# =============================================================================
# Double power law: labels summed on Bagpipes' grid
# =============================================================================
def bagpipes_dblplaw_sfh(tau, alpha, beta, redshift):
    """Bagpipes SFH of one dblplaw component (10^10 M_sun formed, Z_sun).

    Inputs: tau [Gyr since the Big Bang]; alpha, beta, redshift (dimensionless).
    Returns: bagpipes star_formation_history object.
    """
    component = {"tau": tau, "alpha": alpha, "beta": beta, "massformed": 10.0, "metallicity": 1.0}
    return pipes.models.star_formation_history({"redshift": redshift, "dblplaw": component})


def test_bagpipes_grid_is_rebuilt_exactly():
    history = bagpipes_dblplaw_sfh(5.0, 2.0, 3.0, 0.1)
    np.testing.assert_allclose(sfh.BAGPIPES_AGES * 1e9, history.ages, rtol=1e-12, atol=0)
    np.testing.assert_allclose(sfh.BAGPIPES_AGE_WIDTHS * 1e9, history.age_widths, rtol=1e-12, atol=0)
    assert np.isclose(sfh.bagpipes_universe_age(0.1) * 1e9, history.age_of_universe, rtol=1e-12)


# (alpha, beta, redshift, where the turnover tau lies): a sharp peak (alpha = beta = 1000, where the SFR maximum is
# at tau), a flat SFH, early and recent turnovers, asymmetric slopes. "oldest" puts tau at UNIVERSE_AGE_MARGIN.
DBLPLAW_CASES = [
    (1000.0, 1000.0, LOW_REDSHIFT, "early"),
    (1000.0, 1000.0, LOW_REDSHIFT, "recent"),
    (1000.0, 1000.0, LOW_REDSHIFT, "oldest"),
    (1000.0, 1000.0, 0.7, "oldest"),
    (0.01, 0.01, 0.3, "middle"),
    (0.5, 10.0, LOW_REDSHIFT, "early"),
    (10.0, 0.5, 0.7, "recent"),
    (0.01, 1000.0, LOW_REDSHIFT, "recent"),
]
TURNOVER_LOOKBACKS = {"recent": 0.06, "middle": 3.0}     # [Gyr]; "early" is 1 Gyr after the Big Bang
FRACTION_TOLERANCE = 0.01        # absolute, on fractions
DEX_PER_PERCENT = 0.0043         # log10(1.01): a 1% difference in SFR [dex]


@pytest.mark.parametrize("alpha, beta, redshift, where", DBLPLAW_CASES)
def test_dblplaw_labels_match_bagpipes_sfh(alpha, beta, redshift, where):
    universe_age = sfh.COSMOLOGY.age(redshift).value
    if where == "early":
        turnover_lookback = universe_age - 1.0
    elif where == "oldest":
        turnover_lookback = sfh.maximum_age(redshift)
    else:
        turnover_lookback = TURNOVER_LOOKBACKS[where]
    tau = universe_age - turnover_lookback
    history = bagpipes_dblplaw_sfh(tau, alpha, beta, redshift)
    mass = history.sfh * history.age_widths
    ages = history.ages / 1e9
    labels = sfh.dblplaw_labels(tau, alpha, beta, redshift, 10.0, history.stellar_mass)

    assert np.isclose(labels["mw_age"], history.mass_weighted_age, rtol=0.01)
    assert abs(labels["f_young"] - mass[ages < 0.1].sum() / mass.sum()) < FRACTION_TOLERANCE
    assert abs(labels["f_mid"] - mass[(ages >= 0.1) & (ages < 1.0)].sum() / mass.sum()) < FRACTION_TOLERANCE
    assert abs(labels["f_old"] - mass[ages >= 1.0].sum() / mass.sum()) < FRACTION_TOLERANCE
    assert np.isclose(labels["f_young"] + labels["f_mid"] + labels["f_old"], 1.0)

    recent_sfr = mass[ages < sfh.SFR_TIMESCALE].sum() / (sfh.SFR_TIMESCALE * 1e9)
    sfr_floor = sfh.SFR_FLOOR * 10**history.stellar_mass
    assert abs(labels["log_sfr"] - np.log10(max(recent_sfr, sfr_floor))) < DEX_PER_PERCENT
    if history.sfr > sfr_floor:
        assert abs(labels["log_sfr"] - np.log10(history.sfr)) < DEX_PER_PERCENT
    assert np.isclose(labels["log_ssfr"], labels["log_sfr"] - history.stellar_mass)


def test_dblplaw_labels_are_finite_in_the_prior_corners():
    corners = []
    for redshift in (LOW_REDSHIFT, sfh.Z_RANGE[1]):
        for turnover_lookback in (sfh.TURNOVER_LOOKBACK_MIN, sfh.maximum_age(redshift)):
            for alpha in sfh.ALPHA_RANGE:
                for beta in sfh.BETA_RANGE:
                    tau = sfh.COSMOLOGY.age(redshift).value - turnover_lookback
                    corners.append((tau, alpha, beta, redshift))
    tau, alpha, beta, redshift = np.array(corners).T
    labels = sfh.dblplaw_labels(tau, alpha, beta, redshift, np.full(len(tau), 10.0), np.full(len(tau), 9.8))

    assert all(np.isfinite(value).all() for value in labels.values())
    for fraction in ("f_young", "f_mid", "f_old"):
        assert ((labels[fraction] >= 0) & (labels[fraction] <= 1)).all(), fraction
    assert np.allclose(labels["f_young"] + labels["f_mid"] + labels["f_old"], 1.0)
    assert ((labels["mw_age"] > 0) & (labels["mw_age"] <= sfh.bagpipes_universe_age(redshift))).all()


SHARP_TURNOVER_AGE_TOLERANCE = 0.08  # [Gyr] one step of Bagpipes' grid near a lookback time of 13 Gyr


@pytest.mark.parametrize("redshift", [LOW_REDSHIFT, sfh.Z_RANGE[1]])
def test_sharp_turnover_at_the_oldest_lookback_gives_finite_labels(redshift):
    """alpha = beta = 1000 with tau = 0.1 Gyr, the corner the ML decoding reaches when it clips to the prior.

    Why: the SFR is a spike narrower than the grid step there, so every grid SFR is tiny (about e^-300).
    Normalising by each sample's largest grid SFR keeps the sums finite: no 0 / 0.
    """
    turnover_lookback = sfh.maximum_age(redshift)
    tau = sfh.COSMOLOGY.age(redshift).value - turnover_lookback
    assert np.isclose(tau, sfh.UNIVERSE_AGE_MARGIN)
    labels = sfh.dblplaw_labels(tau, 1000.0, 1000.0, redshift, 10.0, 9.8)

    assert all(np.isfinite(value) for value in labels.values())
    assert np.isclose(labels["f_young"] + labels["f_mid"] + labels["f_old"], 1.0)
    # All the mass forms within a grid step of the turnover, more than 1 Gyr ago
    assert np.isclose(labels["f_old"], 1.0)
    assert abs(labels["mw_age"] - turnover_lookback) < SHARP_TURNOVER_AGE_TOLERANCE
    # No star formation in the last 100 Myr, so log_ssfr sits on the floor
    assert np.isclose(labels["log_ssfr"], np.log10(sfh.SFR_FLOOR))


def test_dblplaw_labels_do_not_depend_on_the_chunk_size(monkeypatch):
    """Blocks of 3 samples stop at their own oldest Universe, so they use fewer grid ages than one block of 10.

    This also checks that the grid ages each block skips hold no mass.
    """
    rng = np.random.default_rng(7)
    n_samples = 10
    redshift = rng.uniform(*sfh.Z_RANGE, n_samples)
    tau = sfh.COSMOLOGY.age(redshift).value - rng.uniform(0.1, 5.0, n_samples)
    alpha = 10 ** rng.uniform(-2, 3, n_samples)
    beta = 10 ** rng.uniform(-2, 3, n_samples)
    one_block = sfh.dblplaw_labels(tau, alpha, beta, redshift, np.full(n_samples, 10.0), np.full(n_samples, 9.8))
    monkeypatch.setattr(sfh, "DBLPLAW_CHUNK", 3)
    small_blocks = sfh.dblplaw_labels(tau, alpha, beta, redshift, np.full(n_samples, 10.0),
                                      np.full(n_samples, 9.8))
    for label in sfh.HISTORY_LABELS:
        np.testing.assert_allclose(small_blocks[label], one_block[label], rtol=1e-12, atol=1e-15)


def test_dblplaw_labels_keep_the_input_shape():
    scalar = sfh.dblplaw_labels(5.0, 2.0, 3.0, 0.1, 10.0, 9.8)
    assert all(np.shape(value) == () for value in scalar.values())
    grid = sfh.dblplaw_labels(np.full((2, 3), 5.0), 2.0, 3.0, 0.1, 10.0, 9.8)
    assert all(value.shape == (2, 3) for value in grid.values())


def test_dblplaw_analytic_labels_read_the_drawn_columns():
    _, drawn = sfh.dblplaw_components(np.full(len(sfh.DBLPLAW_DIMS), 0.5))
    drawn_columns = {}
    for label, value in drawn.items():
        drawn_columns[label] = np.array([value])
    from_registry = sfh.DBLPLAW.analytic_labels(drawn_columns, np.array([9.8]))
    direct = sfh.dblplaw_labels(drawn["tau"], drawn["alpha"], drawn["beta"], drawn["redshift"],
                                drawn["log_mass_formed"], 9.8)
    assert list(from_registry) == sfh.HISTORY_LABELS
    for label in sfh.HISTORY_LABELS:
        assert np.isclose(from_registry[label][0], direct[label])


def test_model_file_adds_the_dblplaw_suffix():
    assert sfh.model_file("training_set.parquet", "dblplaw") == "training_set_dblplaw.parquet"
    assert sfh.model_file("photometric_fits", "dblplaw") == "photometric_fits_dblplaw"
    assert sfh.model_file("ml_runs/mlp", "dblplaw") == "ml_runs/mlp_dblplaw"
    assert sfh.model_file("ml_sfh_catalogue.parquet", "dblplaw") == "ml_sfh_catalogue_dblplaw.parquet"
```

  `test_registry_entries_are_consistent` (Task 1) is parametrised over `sfh.MODELS`, so it now also checks `DBLPLAW`.

- [ ] **Step 2: Run them and see them fail.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_sfh_models.py -q
```

  Expected: 1 error during collection, `AttributeError: module 'sfh_models' has no attribute 'ALPHA_RANGE'` (raised by the `parametrize` decorator of the stratification test).

- [ ] **Step 3: Implement.**

  3a. Imports at the top of `sfh_models.py` become

```python
import os
from dataclasses import dataclass
from typing import Callable

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.special import gammainc, gammaincc

import bagpipes as pipes
```

  3b. Insert section 3 after `delayed_analytic_labels` (two blank lines before it and two after, so it sits just above the `# 4. The model registry` banner):

```python
# =============================================================================
# 3. Double power law: SFR(t) ~ [(t / tau)**alpha + (t / tau)**-beta]**-1, t since the Big Bang
# =============================================================================
# One Bagpipes "dblplaw" component with a single metallicity. Star formation starts at the Big Bang, so there is
# no age. tau is the turnover time since the Big Bang, alpha the falling slope after it and beta the rising slope
# before it. The SFR maximum lies at tau only when alpha = beta; in general it lies at t = tau (beta / alpha)**(1 /
# (alpha + beta)). turnover_lookback, the lookback time to tau, is drawn as the delayed-tau age, log-uniform in
# lookback time; alpha and beta are log-uniform (Carnall et al. 2018).
TURNOVER_LOOKBACK_MIN = 0.05     # [Gyr] most recent turnover, as AGE_MIN
ALPHA_RANGE = (0.01, 1000.0)     # falling slope after the turnover (dimensionless)
BETA_RANGE = (0.01, 1000.0)      # rising slope before the turnover (dimensionless)

# "log_turnover_lookback_fraction" places the turnover inside [log TURNOVER_LOOKBACK_MIN, log t_max(z)], as the
# delayed-tau age
DBLPLAW_DIMS = ["z", "log_turnover_lookback_fraction", "log_alpha", "log_beta", "log_metallicity", "av_ism",
                "log_tau_bc", "logU"]

# Bagpipes' own SFH age grid (star_formation_history.__init__, default log_sampling)
BAGPIPES_LOG_AGE_START = 6.0     # [log yr] youngest grid age
BAGPIPES_LOG_AGE_STEP = 0.0025   # [dex] grid step
LOG_GYR_TO_YEARS = 9.0           # [dex] log10 of GYR_TO_YEARS, which Bagpipes adds to the log Hubble time
DBLPLAW_CHUNK = 32               # samples per block of the grid sum: 32 x 1654 grid ages, 0.42 MB per array


def dblplaw_components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one DPL unit-cube vector u (order of DBLPLAW_DIMS).

    Why tau >= UNIVERSE_AGE_MARGIN: turnover_lookback is at most maximum_age(z), so the turnover always lies after
    the Big Bang, and tau stays below Bagpipes' "unphysical" limit of 1.5 times the age of the Universe.
    Inputs: u: unit-cube vector (dimensionless, in [0, 1)) ordered as DBLPLAW_DIMS;
            fixed_redshift: if given, replaces the redshift drawn from u (fits); dimensionless z.
    Returns: (bagpipes_components, drawn). drawn holds the DPL drawn labels: redshift (dimensionless),
             turnover_lookback [Gyr ago], tau [Gyr since the Big Bang], alpha and beta (dimensionless),
             metallicity [Z/Z_sun], av_ism [mag], tau_bc and logU (dimensionless), log_mass_formed [log M_sun].
    """
    params = dict(zip(DBLPLAW_DIMS, u))

    # 1. Redshift, which sets the oldest allowed turnover
    if fixed_redshift is None:
        redshift = uniform(params["z"], *Z_RANGE)
    else:
        redshift = fixed_redshift
    oldest_turnover = maximum_age(redshift)

    # 2. Star formation history and metallicity
    turnover_lookback = 10 ** uniform(params["log_turnover_lookback_fraction"], np.log10(TURNOVER_LOOKBACK_MIN),
                                      np.log10(oldest_turnover))
    # Rounding in 10**log10 can push the top draw a hair above the bound
    turnover_lookback = min(turnover_lookback, oldest_turnover)
    tau = COSMOLOGY.age(redshift).value - turnover_lookback
    alpha = 10 ** uniform(params["log_alpha"], *np.log10(ALPHA_RANGE))
    beta = 10 ** uniform(params["log_beta"], *np.log10(BETA_RANGE))
    metallicity = 10 ** uniform(params["log_metallicity"], *np.log10(METALLICITY_RANGE))

    # 3. Dust and nebular emission
    dust_parts, dust_drawn = dust_and_nebular(params)

    bagpipes_components = {"redshift": redshift}
    bagpipes_components.update(dust_parts)
    bagpipes_components["dblplaw"] = {"tau": tau, "alpha": alpha, "beta": beta, "massformed": LOG_MASS_FORMED,
                                      "metallicity": metallicity}
    drawn = {"redshift": redshift, "turnover_lookback": turnover_lookback, "tau": tau, "alpha": alpha,
             "beta": beta, "metallicity": metallicity, "av_ism": dust_drawn["av_ism"],
             "tau_bc": dust_drawn["tau_bc"], "logU": dust_drawn["logU"], "log_mass_formed": LOG_MASS_FORMED}
    return bagpipes_components, drawn


def bagpipes_age_grid():
    """Bagpipes' SFH age grid, rebuilt exactly as star_formation_history.__init__ builds it.

    Why: the DPL labels have no closed form. Summed on Bagpipes' own grid they match its sfh arrays, including
    its discretisation of sharp peaks (alpha or beta near 1000).
    Returns: (ages, widths) [Gyr], two arrays of the grid length: the lookback age of each grid point and the
             width of its bin. The first bin starts at 0 and the last ends at the Hubble time.
    """
    hubble_time = pipes.utils.age_at_z[pipes.utils.z_array == 0.0][0]     # [Gyr]
    log_age_max = np.log10(hubble_time) + LOG_GYR_TO_YEARS + 2 * BAGPIPES_LOG_AGE_STEP
    log_ages = np.arange(BAGPIPES_LOG_AGE_START, log_age_max, BAGPIPES_LOG_AGE_STEP)    # [log yr]
    log_edges = pipes.utils.make_bins(log_ages, make_rhs=True)[0]                       # [log yr]
    edges = 10**log_edges                                                                # [yr]
    edges[0] = 0.0
    edges[-1] = GYR_TO_YEARS * hubble_time
    ages = 10**log_ages / GYR_TO_YEARS
    widths = (edges[1:] - edges[:-1]) / GYR_TO_YEARS
    return ages, widths


BAGPIPES_AGES, BAGPIPES_AGE_WIDTHS = bagpipes_age_grid()      # [Gyr]

# The grid ages increase, so each label is a sum over one slice of the grid. Each index is the first grid age at
# or above the edge, which matches Bagpipes' "ages < edge"
YOUNG_EDGE_INDEX = int(np.searchsorted(BAGPIPES_AGES, LABEL_EDGES[0]))     # f_young: grid ages below 0.1 Gyr
OLD_EDGE_INDEX = int(np.searchsorted(BAGPIPES_AGES, LABEL_EDGES[1]))       # f_old: grid ages from 1 Gyr on
RECENT_EDGE_INDEX = int(np.searchsorted(BAGPIPES_AGES, SFR_TIMESCALE))     # log_sfr: grid ages below 0.1 Gyr


def bagpipes_universe_age(redshift):
    """Age of the Universe [Gyr] at the redshift, interpolated in Bagpipes' redshift table as Bagpipes does.

    Why not COSMOLOGY.age: the cosmology is the same, but Bagpipes interpolates; its value decides which grid
    ages lie after the Big Bang, exactly as in Bagpipes.
    Inputs: redshift (dimensionless), scalar or array.
    Returns: age [Gyr], same shape.
    """
    return np.interp(redshift, pipes.utils.z_array, pipes.utils.age_at_z)


def dblplaw_grid_sums(tau, alpha, beta, redshift):
    """Mass-weighted age and formed-mass fractions of double-power-law SFHs, summed on Bagpipes' age grid.

    The SFR is evaluated at every grid age younger than the Universe at the redshift (t = age of the Universe
    - lookback age), as Bagpipes' dblplaw does. Samples go DBLPLAW_CHUNK at a time, vectorised over the block,
    so the arrays stay in the CPU cache. A block uses only the grid ages below the age of the Universe at its
    lowest redshift, because older ones hold no stars.
    Why in logs: alpha and beta reach 1000, so (t / tau)**alpha overflows. With x = ln(t / tau), the exponents
    are alpha x and -beta x, the larger of the two is L >= 0, and
    SFR = 1 / (e^(alpha x) + e^(-beta x)) = e^-L / (1 + e^-((alpha + beta) |x|)).
    e^-L is divided by its largest value in each sample's row, which leaves the normalised SFH unchanged and
    keeps the sum finite when the SFR is a spike between two grid ages.
    Why no matrix product: a threaded BLAS product of these small blocks is slower than numpy's sums.
    Inputs: tau [Gyr since the Big Bang], alpha, beta, redshift (dimensionless): 1-D arrays of one length.
    Returns: dict of 1-D arrays: "mw_age" [Gyr ago]; "f_young", "f_mid", "f_old" (formed-mass fractions split
             at LABEL_EDGES, sum to 1); "recent" (formed-mass fraction younger than SFR_TIMESCALE).
    """
    n_samples = len(tau)
    sums = {}
    for key in ("mw_age", "f_young", "f_mid", "f_old", "recent"):
        sums[key] = np.zeros(n_samples)

    for start in range(0, n_samples, DBLPLAW_CHUNK):
        block = slice(start, start + DBLPLAW_CHUNK)
        block_tau = tau[block, np.newaxis]             # [Gyr], (block, 1)
        block_alpha = alpha[block, np.newaxis]
        block_beta = beta[block, np.newaxis]
        universe_age = bagpipes_universe_age(redshift[block])      # [Gyr], (block,)

        # 1. Grid ages that the oldest Universe of the block reaches
        n_ages = int(np.searchsorted(BAGPIPES_AGES, universe_age.max()))
        ages = BAGPIPES_AGES[:n_ages]
        widths = BAGPIPES_AGE_WIDTHS[:n_ages]

        # 2. x = ln(t / tau), with t the time since the Big Bang of each grid age, (block, n_ages). Grid ages
        #    before the Big Bang get t = tau (x = 0), so the log stays finite; step 4 gives them no mass
        log_time_ratio = np.subtract.outer(universe_age, ages)     # t [Gyr]
        before_big_bang = log_time_ratio <= 0
        np.copyto(log_time_ratio, block_tau, where=before_big_bang)
        log_time_ratio /= block_tau
        np.log(log_time_ratio, out=log_time_ratio)

        # 3. The two exponents of SFR = 1 / (e^(alpha x) + e^(-beta x))
        falling_exponent = block_alpha * log_time_ratio
        rising_exponent = -block_beta * log_time_ratio
        larger_exponent = np.maximum(falling_exponent, rising_exponent)          # L >= 0
        exponent_gap = np.abs(falling_exponent - rising_exponent)                # (alpha + beta) |x|

        # 4. Formed mass per grid bin, relative to the largest grid SFR of each sample
        larger_exponent[before_big_bang] = np.inf
        smallest_exponent = larger_exponent.min(axis=1, keepdims=True)
        formed_mass = np.exp(smallest_exponent - larger_exponent) / (1 + np.exp(-exponent_gap)) * widths

        # 5. The sums, as Bagpipes' mass_weighted_age and the label edges
        total_mass = formed_mass.sum(axis=1)
        sums["mw_age"][block] = (formed_mass * ages).sum(axis=1) / total_mass
        sums["f_young"][block] = formed_mass[:, :YOUNG_EDGE_INDEX].sum(axis=1) / total_mass
        sums["f_mid"][block] = formed_mass[:, YOUNG_EDGE_INDEX:OLD_EDGE_INDEX].sum(axis=1) / total_mass
        sums["f_old"][block] = formed_mass[:, OLD_EDGE_INDEX:].sum(axis=1) / total_mass
        sums["recent"][block] = formed_mass[:, :RECENT_EDGE_INDEX].sum(axis=1) / total_mass
    return sums


def dblplaw_labels(tau, alpha, beta, redshift, log_mass_formed, log_mstar):
    """Labels of a double-power-law SFH, shared by the training sets and the ML decoding.

    mw_age and the mass fractions have no closed form: they are summed on Bagpipes' age grid (dblplaw_grid_sums).
    log_sfr is the formed mass younger than SFR_TIMESCALE over SFR_TIMESCALE, floored at SFR_FLOOR * M_star,
    as for delayed-tau.
    Inputs: tau [Gyr since the Big Bang], alpha, beta, redshift (dimensionless), log_mass_formed and log_mstar
            [log M_sun]; arrays that broadcast together, or scalars.
    Returns: dict of arrays of the broadcast shape, the keys of delayed_tau_labels: mw_age [Gyr ago], f_young,
             f_mid, f_old (dimensionless, sum to 1), log_sfr [log M_sun/yr] and log_ssfr [log yr^-1].
    """
    tau, alpha, beta, redshift, log_mass_formed, log_mstar = np.broadcast_arrays(
        np.asarray(tau, dtype=float), np.asarray(alpha, dtype=float), np.asarray(beta, dtype=float),
        np.asarray(redshift, dtype=float), np.asarray(log_mass_formed, dtype=float),
        np.asarray(log_mstar, dtype=float))
    shape = tau.shape
    sums = dblplaw_grid_sums(tau.ravel(), alpha.ravel(), beta.ravel(), redshift.ravel())

    recent_sfr = 10**log_mass_formed.ravel() * sums["recent"] / (SFR_TIMESCALE * GYR_TO_YEARS)
    sfr_floor = SFR_FLOOR * 10**log_mstar.ravel()
    log_sfr = np.log10(np.maximum(recent_sfr, sfr_floor))
    labels = {
        "mw_age": sums["mw_age"],
        "f_young": sums["f_young"],
        "f_mid": sums["f_mid"],
        "f_old": sums["f_old"],
        "log_sfr": log_sfr,
        "log_ssfr": log_sfr - log_mstar.ravel(),
    }
    for label in labels:
        labels[label] = labels[label].reshape(shape)
    return labels


def dblplaw_analytic_labels(drawn_columns, log_mstar):
    """The HISTORY_LABELS of double-power-law models from their drawn labels (the registry's analytic_labels).

    Inputs: drawn_columns: dict of arrays holding at least "tau" [Gyr since the Big Bang], "alpha", "beta",
            "redshift" (dimensionless) and "log_mass_formed" [log M_sun]; log_mstar: array [log M_sun].
    Returns: dict of arrays, as dblplaw_labels.
    """
    return dblplaw_labels(drawn_columns["tau"], drawn_columns["alpha"], drawn_columns["beta"],
                          drawn_columns["redshift"], drawn_columns["log_mass_formed"], log_mstar)
```

  3c. In section 4, replace the line `MODELS = {"delayed": DELAYED}` with

```python
DBLPLAW = SfhModel(
    name="dblplaw",
    dims=DBLPLAW_DIMS,
    components=dblplaw_components,
    drawn_labels=["redshift", "turnover_lookback", "tau", "alpha", "beta", "metallicity", "av_ism", "tau_bc", "logU",
                  "log_mass_formed"],
    derived_labels=["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "log_mstar", "log_sfr",
                    "log_ssfr", "mw_age", "f_young", "f_mid", "f_old"],
    log_labels={"log_turnover_lookback": "turnover_lookback", "log_alpha": "alpha", "log_beta": "beta",
                "log_metallicity": "metallicity"},
    analytic_labels=dblplaw_analytic_labels,
    target_names=["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "av_ism", "log_tau_bc", "logU",
                  "log_mass_formed", "log_mstar", "redshift"],
    calibration_labels=["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "log_mstar"],
    file_suffix="_dblplaw",
)

MODELS = {"delayed": DELAYED, "dblplaw": DBLPLAW}
```

- [ ] **Step 4: Run the tests and see them pass.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_sfh_models.py -q
```

  Expected: `34 passed, 6 warnings` in about 8 s. The 6 warnings are `RuntimeWarning: overflow encountered in power`. They come from Bagpipes' own `dblplaw` (in the corner test at cube value 0.999999 and in five of the `test_dblplaw_labels_match_bagpipes_sfh` cases with α or β = 1000), which overflows to inf and so gives SFR 0, as intended. `dblplaw_labels` works in logs and raises none.

  The 34 tests are Task 1's 7 plus 27. The 27 are:
  - the `dblplaw` case of `test_registry_entries_are_consistent`;
  - 10 prior tests: 1 + 1 + 1 + 2 + 1 + 1, then 3 corner cases;
  - `test_bagpipes_grid_is_rebuilt_exactly`;
  - the 8 Bagpipes comparison cases;
  - the corner-finiteness test;
  - the 2 sharp-turnover cases;
  - the chunk-size, shape, analytic-labels and `model_file` tests.

- [ ] **Step 5: Check the speed of the grid sum.** This is a quick check (about 2 s), not a test:

```bash
cd Synthetic_domain && ../.venv/bin/python -c "
import time
import numpy as np
import sfh_models as sfh
rng = np.random.default_rng(0)
n_samples = 20000
redshift = rng.uniform(*sfh.Z_RANGE, n_samples)
log_lookback = rng.uniform(np.log10(sfh.TURNOVER_LOOKBACK_MIN), np.log10(sfh.maximum_age(redshift)))
tau = sfh.COSMOLOGY.age(redshift).value - 10**log_lookback
alpha = 10 ** rng.uniform(*np.log10(sfh.ALPHA_RANGE), n_samples)
beta = 10 ** rng.uniform(*np.log10(sfh.BETA_RANGE), n_samples)
start = time.perf_counter()
sfh.dblplaw_grid_sums(tau, alpha, beta, redshift)
print(f'{(time.perf_counter() - start) / n_samples * 1e3:.3f} ms per sample')
"
```

  Expected: between `0.05` and `0.08 ms per sample` on the laptop. Above 0.1 ms means the code has drifted from the version above, for example to larger blocks, a BLAS product (`@`) or boolean-mask sums. Fix it before going on: Task 7 decodes about 1e8 samples with this function.

- [ ] **Step 6: Run the whole suite.**

```bash
cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q
```

  Expected: `193 passed, 6 warnings` (the Bagpipes overflow warnings of Step 4). That is Task 1's 166 plus the 27 tests this task adds to `test_sfh_models.py`: 26 new functions and cases, plus the `dblplaw` case of `test_registry_entries_are_consistent`.

- [ ] **Step 7: Commit.**

```bash
git add Synthetic_domain/sfh_models.py Synthetic_domain/tests/test_sfh_models.py
git commit -m "$(cat <<'EOF'
feat: double-power-law prior and labels summed on Bagpipes' own age grid (sfh_models section 3)

The grid sum is vectorised over blocks of 32 samples and stops at each block's oldest Universe. It takes about
0.06 ms per sample on one laptop core. alpha = beta = 1000 with tau = 0.1 Gyr gives finite labels.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

### Task 3: `simulate_models` and `summary_labels` per SFH model; `training_set.ipynb` builds the sets of `SFH_MODELS`

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py`:
  - section 2: `summary_labels` (today lines 236–256, after Task 1 about lines 101–121) and `simulate_models` (today lines 271–296);
  - section 6: `sample_label_row`, its call `row = summary_labels(model.sfh, drawn)`.
- Modify: `Synthetic_domain/training_set.ipynb`. The cells, by id:
  - markdown `afa6d90a` (title), code `f136f1c3` (parameters);
  - markdown `5b628636` (§1), code `15413c23` (§1);
  - markdown `742607a1` (§4), code `9ff91c28` (§4);
  - markdown `33ef32cd` (§5), code `bd81d555` (§5);
  - markdown `fd06b7a0` (§6), code `64cba55b` (§6).
  
  Cells `1a99cd90`, `7be5189c`, `11040f58` and `19201d8e` (catalogue and noise model) are unchanged.
- Test:
  - `Synthetic_domain/tests/test_synthetic_model.py`: `test_summary_labels_give_every_training_label` and the three `test_simulated_*` tests, plus a new DPL summary test;
  - `Synthetic_domain/tests/test_ml_data.py`: the `simulated_models` fixture;
  - `Synthetic_domain/tests/test_sfh_models.py`: the two regression tests.

**Interfaces:**
- Consumes (Tasks 1–2): `MODELS`, `HISTORY_LABELS`, each model's `components`, `drawn_labels`, `log_labels`, `analytic_labels`, `training_labels`, `dims`; `model_file`.
- Produces:
  - `summary_labels(bagpipes_sfh, drawn, sfh_model) -> dict`, ordered as `MODELS[sfh_model].training_labels`;
  - `simulate_models(unit_vectors, sfh_model) -> DataFrame`, columns `MODELS[sfh_model].training_labels + BANDS`.

  Both take `sfh_model` as a required argument. `sample_label_row` passes the alias `SFH_MODEL` ("delayed") until Task 4 threads the job's model through the fits.
  - The notebook: `SFH_MODELS`, `REBUILD_SETS`, `set_files(sfh_model)`, `simulate_in_parallel(n_models, seed, sfh_model)`, `coverage_in_parallel(..., training_magnitudes, training_z)`, `colour_colour_figure(covered, training_magnitudes, output_file)`, `p_value_figure(splus_p_values, validation_p_values, output_file)`. The outputs per model come from `sm.model_file`: `training_set_dblplaw.parquet`, `synthetic_validation_dblplaw.parquet`, `synthetic_test_dblplaw.parquet`, `color_color_dblplaw.png`, `coverage_criterion_dblplaw.png`.

- [ ] **Step 1: Write the failing tests.**

  1a. `tests/test_synthetic_model.py`: add `import sfh_models` between `import bagpipes as pipes` and `import synthetic_model as sm`. Replace the whole `test_summary_labels_give_every_training_label` function (today lines 175–188) with:

```python
def test_summary_labels_give_every_training_label():
    u = np.random.default_rng(6).random(len(sm.MODELS["delayed"].dims))
    bagpipes_components, drawn = sm.MODELS["delayed"].components(u)
    sfh = pipes.models.star_formation_history(bagpipes_components)
    labels = sm.summary_labels(sfh, drawn, "delayed")

    assert list(labels) == sm.MODELS["delayed"].training_labels
    assert all(np.isfinite(value) for value in labels.values())
    assert labels["log_mstar"] == sfh.stellar_mass
    assert np.isclose(labels["log_age"], np.log10(drawn["age"]))
    assert np.isclose(labels["log_tau"], np.log10(drawn["tau"]))
    assert np.isclose(labels["log_metallicity"], np.log10(drawn["metallicity"]))
    assert np.isclose(labels["mw_age"], sm.delayed_mass_weighted_age(drawn["age"], drawn["tau"]))


def test_dblplaw_summary_labels_give_every_training_label():
    dblplaw = sm.MODELS["dblplaw"]
    u = np.random.default_rng(6).random(len(dblplaw.dims))
    bagpipes_components, drawn = dblplaw.components(u)
    sfh = pipes.models.star_formation_history(bagpipes_components)
    labels = sm.summary_labels(sfh, drawn, "dblplaw")

    assert list(labels) == dblplaw.training_labels
    assert all(np.isfinite(value) for value in labels.values())
    assert labels["log_mstar"] == sfh.stellar_mass
    for log_label, drawn_label in dblplaw.log_labels.items():
        assert np.isclose(labels[log_label], np.log10(drawn[drawn_label])), log_label
    expected = sfh_models.dblplaw_labels(drawn["tau"], drawn["alpha"], drawn["beta"], drawn["redshift"],
                                         drawn["log_mass_formed"], sfh.stellar_mass)
    for label in sm.HISTORY_LABELS:
        assert np.isclose(labels[label], expected[label]), label
    assert np.isclose(labels["mw_age"], sfh.mass_weighted_age, rtol=0.01)
```

  1b. Same file: replace `test_simulated_columns_follow_training_labels_then_bands`, `test_simulated_models_are_finite_and_their_fractions_sum_to_one` and `test_simulated_redshift_and_z_band_are_distinct_columns` (today lines 197–221, up to the `# Detection rule and noise model` banner) with:

```python
@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_simulated_columns_follow_training_labels_then_bands(sfh_model):
    model_definition = sm.MODELS[sfh_model]
    table = sm.simulate_models(np.full((2, len(model_definition.dims)), 0.5), sfh_model)
    assert list(table.columns) == model_definition.training_labels + sm.BANDS
    assert np.isfinite(table[sm.BANDS].to_numpy()).all()


@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_simulated_models_are_finite_and_their_fractions_sum_to_one(sfh_model):
    model_definition = sm.MODELS[sfh_model]
    unit_vectors = np.random.default_rng(2).random((4, len(model_definition.dims)))
    models = sm.simulate_models(unit_vectors, sfh_model)

    assert len(models) == 4
    assert np.isfinite(models[model_definition.training_labels + sm.BANDS].to_numpy()).all()
    assert np.allclose(models[["f_young", "f_mid", "f_old"]].sum(axis=1), 1.0)
    assert np.allclose(models["log_ssfr"], models["log_sfr"] - models["log_mstar"])


@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_simulated_redshift_and_z_band_are_distinct_columns(sfh_model):
    unit_vectors = np.random.default_rng(3).random((4, len(sm.MODELS[sfh_model].dims)))
    models = sm.simulate_models(unit_vectors, sfh_model)

    assert models["redshift"].between(*sm.Z_RANGE).all()
    # The z band is an AB magnitude; before the fix it held the redshift (below 1)
    assert models["z"].between(MIN_PLAUSIBLE_MAGNITUDE, MAX_PLAUSIBLE_MAGNITUDE).all()
```

  1c. `tests/test_ml_data.py`, the `simulated_models` fixture (lines 27–28): replace

```python
    unit_vectors = qmc.LatinHypercube(d=len(sm.DIMS), seed=30).random(200)
    return sm.simulate_models(unit_vectors)
```

  with

```python
    unit_vectors = qmc.LatinHypercube(d=len(sm.MODELS["delayed"].dims), seed=30).random(200)
    return sm.simulate_models(unit_vectors, "delayed")
```

  The old delayed-τ tests change in only two ways. They pass the new required `sfh_model` argument, and they use `sm.MODELS["delayed"]` in place of the aliases that Task 11 removes. Every delayed-τ assertion stays as it is: compare with the replaced lines before saving, and do not loosen any of them.

  1d. `tests/test_sfh_models.py`: in the two regression tests, `sm.simulate_models(first_unit_vectors(TRAINING_SEED, N_TRAINING_MODELS))` becomes `sm.simulate_models(first_unit_vectors(TRAINING_SEED, N_TRAINING_MODELS), "delayed")`, and `sm.simulate_models(first_unit_vectors(VALIDATION_SEED, N_HELD_OUT_MODELS))` becomes `sm.simulate_models(first_unit_vectors(VALIDATION_SEED, N_HELD_OUT_MODELS), "delayed")`.

- [ ] **Step 2: Run them and see them fail.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_synthetic_model.py tests/test_sfh_models.py -q -k "summary or simulated or delayed_path"
```

  Expected: 10 failed, with `TypeError: summary_labels() takes 2 positional arguments but 3 were given` and `TypeError: simulate_models() takes 1 positional argument but 2 were given`.

- [ ] **Step 3: Implement in `synthetic_model.py`.**

  3a. Replace the whole `summary_labels` function with:

```python
def summary_labels(bagpipes_sfh, drawn, sfh_model):
    """Every training label of one model of the given SFH model.

    Inputs: bagpipes_sfh: the model's sfh object (stellar_mass in log M_sun, living stars);
            drawn: dict from the model's components() with its drawn labels (extra keys are ignored);
            sfh_model: a name in MODELS ("delayed" or "dblplaw").
    Returns: dict ordered as MODELS[sfh_model].training_labels: the drawn values; the log labels (base-10 logs
             of drawn values); log_mstar from Bagpipes; the HISTORY_LABELS from the analytic_labels, as floats.
    """
    model_definition = MODELS[sfh_model]
    labels = {}
    for label in model_definition.drawn_labels:
        labels[label] = drawn[label]
    for log_label, drawn_label in model_definition.log_labels.items():
        labels[log_label] = np.log10(drawn[drawn_label])
    labels["log_mstar"] = bagpipes_sfh.stellar_mass

    # analytic_labels is vectorised: give it one-element columns
    drawn_columns = {}
    for label in model_definition.drawn_labels:
        drawn_columns[label] = np.array([drawn[label]])
    history_labels = model_definition.analytic_labels(drawn_columns, np.array([bagpipes_sfh.stellar_mass]))
    for label in HISTORY_LABELS:
        labels[label] = float(history_labels[label][0])

    ordered_labels = {}
    for label in model_definition.training_labels:
        ordered_labels[label] = labels[label]
    return ordered_labels
```

  3b. Replace the whole `simulate_models` function with:

```python
def simulate_models(unit_vectors, sfh_model):
    """Labels and AB magnitudes (10^10 M_sun formed) for a block of unit-cube vectors of one SFH model.

    Inputs: unit_vectors: iterable of unit-cube vectors (dimensionless, in [0, 1)), each ordered as
            MODELS[sfh_model].dims; sfh_model: a name in MODELS ("delayed" or "dblplaw").
    Returns: DataFrame with one row per model: the model's training_labels, then the AB magnitudes [mag] in
             BANDS order (NaN where a flux is not positive). The "z" column is the z-band magnitude, not the
             redshift.

    One Bagpipes model is built and then updated: update() is ~70x faster than building a new one, and safe
    because every model has the same component set (Bagpipes fixes the components when it builds a model).
    """
    model_definition = MODELS[sfh_model]
    model = None
    rows = []
    for u in unit_vectors:
        bagpipes_components, drawn = model_definition.components(u)
        if model is None:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
        else:
            model.update(bagpipes_components)

        row = summary_labels(model.sfh, drawn, sfh_model)
        magnitudes = ab_magnitudes(model.photometry)
        for band, magnitude in zip(BANDS, magnitudes):
            row[band] = magnitude
        rows.append(row)
    return pd.DataFrame(rows, columns=model_definition.training_labels + BANDS)
```

  3c. In `sample_label_row` (section 6), replace `    row = summary_labels(model.sfh, drawn)` with `    row = summary_labels(model.sfh, drawn, SFH_MODEL)`. The posterior labels of the delayed-τ fits are unchanged; Task 4 replaces `SFH_MODEL` with the job's model.

  Check that no caller is left with the old signatures. `grep -rn "summary_labels(\|simulate_models(" Synthetic_domain --include=*.py --include=*.ipynb | grep -v "/pipes/"` should list only:
  - the two definitions;
  - `sample_label_row`;
  - the tests above.

  The notebook calls `simulate_models` through joblib, as `delayed(sm.simulate_models)(block)`, which this grep does not match; Step 5 rewrites that call.

- [ ] **Step 4: Run the tests and see them pass.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_synthetic_model.py tests/test_sfh_models.py -q -k "summary or simulated or delayed_path"
```

  Expected: `10 passed`. The delayed-τ regression tests still pass, because the registry path gives bit-identical delayed-τ labels.

- [ ] **Step 5: Rewrite the notebook cells.** Do not execute the notebook. `$SCRATCH` below is your scratchpad directory. Write each source below to a file named after its cell id in a scratch folder (for example `$SCRATCH/training_set_cells/afa6d90a.md`, `.../f136f1c3.py`, ...; `.md` for markdown, `.py` for code). Then apply all of them with the script after the cells. It keeps the cell ids and clears the stale outputs and execution counts of the code cells it rewrites. NotebookEdit with the same cell ids is an equivalent route; the outputs of rewritten code cells must be cleared either way.

  Markdown `afa6d90a`, which replaces the cell starting `# Synthetic S-PLUS training set (part 1)`:

```markdown
# Synthetic S-PLUS training sets (part 1)

Bagpipes (BC03-MILES) models in the 12 S-PLUS bands, drawn by Latin Hypercube from wide priors fixed **before** looking at the data, for each SFH model listed in `SFH_MODELS`. All prior code lives in `sfh_models.py` (one section per SFH model) and `synthetic_model.py`, shared with `spectroscopic_reference.ipynb` and the ML.

**Inputs:** `splus-filters/`, `dr3_corrected.csv` (S-PLUS DR3 corrected for Galactic extinction).
**Outputs, per SFH model:** `training_set.parquet` (noise-free), `synthetic_validation.parquet` and `synthetic_test.parquet` (S-PLUS noise and non-detections), `color_color.png`, `coverage_criterion.png`. Delayed-τ keeps these names; the double power law adds `_dblplaw` (`training_set_dblplaw.parquet`, ...). Shared: `splus_noise_model.csv`.
**Skipping:** a model whose three set files exist is not rebuilt unless `REBUILD_SETS = True`, so the stored delayed-τ sets are never rewritten by accident. The coverage test and the figures run only for the models built in this run, so a kept model's coverage figures are not recomputed either. `splus_noise_model.csv` is written on every run, because the held-out sets need it; it is a deterministic fit of the same DR3 catalogue, so its content does not change.
**Run:** top to bottom with the *Inference4ETGs (MultiNest)* kernel, as a notebook job on the PC (8 cores, 7 workers). Never on the laptop.

**Shared by both SFH models:** the same sizes and seeds, a single metallicity $Z$, Charlot & Fall (2000) dust (exponent 0.7, birth-cloud age 0.01 Gyr), Bagpipes nebular emission, the Kroupa IMF and a flat cosmology with $H_0 = 70$, $\Omega_m = 0.3$.

| parameter | prior |
|---|---|
| `redshift` $z$ | $U(0.002, 0.7)$ |
| `metallicity` $Z$ | log-uniform, 0.005–2.5 $Z_\odot$ |
| `av_ism` (ISM $A_V$) | $U(0.01, 2)$ mag |
| `tau_bc` (birth-cloud optical depth; extra $A_V$ = $1.086\,\tau_{\rm BC}$ for stars < 10 Myr) | log-uniform, 0.5–4 |
| `logU` | $U(-4, -1)$ |
| `log_mass_formed` | fixed at 10 for the models; rescaled to the observed $r$ magnitude |

**Delayed-τ (`delayed`):** one Bagpipes `delayed` component, $\mathrm{SFR}(t) \propto (t - t_0)\,e^{-(t - t_0)/\tau}$. The age is the time since star formation began; it never exceeds the age of the Universe at the model's redshift minus 0.1 Gyr.

| parameter | prior |
|---|---|
| `age` (time since star formation began) | log-uniform, 0.05 Gyr to $t_{\max}(z)$ = age of the Universe at $z$ minus 0.1 Gyr |
| `tau` $\tau$ | log-uniform, 0.01–20 Gyr |

**Double power law (`dblplaw`):** one Bagpipes `dblplaw` component, $\mathrm{SFR}(t) \propto \left[(t/\tau)^{\alpha} + (t/\tau)^{-\beta}\right]^{-1}$, with $t$ the time since the Big Bang. Star formation starts at the Big Bang, so there is no age parameter.

| parameter | prior |
|---|---|
| `turnover_lookback` (lookback time to the turnover $\tau$) | log-uniform, 0.05 Gyr to $t_{\max}(z)$ |
| `tau` $\tau$ (turnover time since the Big Bang) | age of the Universe at $z$ minus `turnover_lookback`, so at least 0.1 Gyr |
| `alpha` $\alpha$ (falling slope, after the turnover) | log-uniform, 0.01–1000 |
| `beta` $\beta$ (rising slope, before the turnover) | log-uniform, 0.01–1000 |

$\tau$ is the turnover of the double power law. It is the SFR maximum only when $\alpha = \beta$. In general the maximum lies at $t = \tau\,(\beta/\alpha)^{1/(\alpha+\beta)}$, so for $\beta \ll \alpha$ the SFR may still be rising when the galaxy is observed.

**Drawn labels:**
- delayed-τ: `redshift, age, tau, metallicity, av_ism, tau_bc, logU, log_mass_formed`;
- DPL: `redshift, turnover_lookback, tau, alpha, beta, metallicity, av_ism, tau_bc, logU, log_mass_formed`.

The redshift is called `redshift`, never `z`, because `z` is the $z$-band column.

**Derived labels:**
- the base-10 logarithms of drawn values: `log_age`, `log_tau`, `log_metallicity` (delayed-τ); `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity` (DPL);
- `log_mstar` [$\log M_\odot$]: the surviving stellar mass, from Bagpipes;
- `log_sfr` [$\log M_\odot\,{\rm yr}^{-1}$]: the mean SFR over the last 100 Myr, with a floor of $10^{-14}\,{\rm yr}^{-1}$ times the stellar mass. About 20% of the delayed-τ prior has no star formation in that window and sits exactly on the floor;
- `log_ssfr` = `log_sfr` $-$ `log_mstar`;
- `mw_age` [Gyr ago]: the mass-weighted age;
- `f_young`, `f_mid`, `f_old`: the formed-mass fractions younger than 0.1 Gyr, between 0.1 and 1 Gyr, and older than 1 Gyr. They sum to 1.

`log_sfr`, `mw_age` and the fractions are analytic for delayed-τ. For DPL they have no closed form and are summed on Bagpipes' own age grid, so they match its SFH arrays.

For an observed galaxy, $\Delta m = m_{\rm obs} - m_{\rm model}$ rescales every mass-like label (`log_mass_formed`, `log_mstar`, `log_sfr`; `log_ssfr` does not move) by $-0.4\,\Delta m$ dex.
```

  Code `f136f1c3`, which replaces the parameter cell starting `import os`:

```python
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.stats import qmc

import synthetic_model as sm

# SFH models whose sets this run builds (names in sm.MODELS); a model whose three set files exist is skipped,
# unless REBUILD_SETS is True
SFH_MODELS = ["delayed", "dblplaw"]
REBUILD_SETS = False

# Sizes and seeds, the same for every SFH model
N_MODELS = 2**19                 # training set
N_HELD_OUT = 2**14               # each synthetic validation / test set
TRAINING_SEED = 0
VALIDATION_SEED = 1
TEST_SEED = 2
N_JOBS = int(os.environ.get("N_JOBS", 7))   # a PC job gets 8 cores, so 7 workers fit
BLOCKS_PER_JOB = 20              # small blocks keep the workers busy until the end

# Files: the delayed-tau names; sm.model_file adds each SFH model's suffix
SPLUS_FILE = "dr3_corrected.csv"
TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
TEST_FILE = "synthetic_test.parquet"
NOISE_MODEL_FILE = "splus_noise_model.csv"
COLOUR_FIGURE_FILE = "color_color.png"
P_VALUE_FIGURE_FILE = "coverage_criterion.png"

# Columns summarised after each training set is simulated
DISPLAY_COLUMNS = {
    "delayed": ["redshift", "age", "tau", "metallicity", "log_mstar", "log_sfr", "mw_age", "f_old"],
    "dblplaw": ["redshift", "turnover_lookback", "alpha", "beta", "metallicity", "log_mstar", "log_sfr", "mw_age",
                "f_old"],
}

MIN_PROB_GALAXY = 0.8
OLD_GEOMETRIC_COVERAGE = 0.991   # geometric test of commit c04268e (K = 2), kept for one comparison

# Coverage test and its summary tables
COVERAGE_BLOCKS_PER_JOB = 4                                       # galaxy blocks per worker
ZML_GROUP_EDGES = [-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 10]   # photometric redshift (dimensionless)
R_GROUP_EDGES = [0, 17, 18, 19, 20, 20.5, 21, 21.5, 30]           # r_auto [mag]
P_VALUE_FLOOR = 1e-8                                              # p-values are clipped here before taking log10
LOG_P_RANGE = (-8, 0)                                             # log10 p histogram range
N_LOG_P_BINS = 40
```

  Markdown `5b628636`, which replaces `## 1. Training set`:

```markdown
## 1. Training sets
Each model is a Latin Hypercube vector $u \in [0, 1)^{d}$, one dimension per entry of its SFH model's `dims` ($d = 7$ for delayed-τ, 8 for DPL). The model's `components` maps it to Bagpipes components by inverse CDFs, so the priors are exactly those of the tables. The delayed-τ age and the DPL turnover lookback time are fractions of the log range up to $t_{\max}(z)$, so one cube works at every redshift. An SFH model whose three set files exist is skipped unless `REBUILD_SETS = True`.
```

  Code `15413c23`, which replaces the cell starting `def simulate_in_parallel(n_models, seed):`:

```python
def set_files(sfh_model):
    """Training, synthetic validation and synthetic test file names of one SFH model.

    Inputs: sfh_model: a name in sm.MODELS ("delayed" or "dblplaw").
    Returns: dict "training", "validation", "test" -> file name (the delayed-tau names have no suffix).
    """
    return {
        "training": sm.model_file(TRAINING_FILE, sfh_model),
        "validation": sm.model_file(VALIDATION_FILE, sfh_model),
        "test": sm.model_file(TEST_FILE, sfh_model),
    }


def simulate_in_parallel(n_models, seed, sfh_model):
    """Simulate Bagpipes models of one SFH model drawn by Latin Hypercube, in parallel on N_JOBS workers.

    Inputs: n_models: number of models (count); seed: Latin Hypercube seed (integer);
            sfh_model: a name in sm.MODELS.
    Returns: DataFrame with one row per model: the model's training labels and the 12 magnitudes [mag].
    """
    sampler = qmc.LatinHypercube(len(sm.MODELS[sfh_model].dims), seed=seed)
    unit_vectors = sampler.random(n_models)
    blocks = np.array_split(unit_vectors, BLOCKS_PER_JOB * N_JOBS)
    tasks = []
    for block in blocks:
        tasks.append(delayed(sm.simulate_models)(block, sfh_model))
    tables = Parallel(n_jobs=N_JOBS, verbose=5)(tasks)
    return pd.concat(tables, ignore_index=True)


# 1. SFH models built in this run; the others keep their stored sets
models_to_build = []
for sfh_model in SFH_MODELS:
    all_sets_exist = True
    for file_name in set_files(sfh_model).values():
        if not os.path.exists(file_name):
            all_sets_exist = False
    if REBUILD_SETS or not all_sets_exist:
        models_to_build.append(sfh_model)
    else:
        print(f"{sfh_model}: the three set files exist, so they are kept (REBUILD_SETS = False)")

# 2. Noise-free training sets
for sfh_model in models_to_build:
    training_file = set_files(sfh_model)["training"]
    training_models = simulate_in_parallel(N_MODELS, TRAINING_SEED, sfh_model)
    training_models.to_parquet(training_file, index=False)

    # Check: Bagpipes returned positive fluxes in every band for every model
    finite_magnitudes = np.isfinite(training_models[sm.BANDS].to_numpy()).all(axis=1)
    print(f"{sfh_model}: {len(training_models)} models saved to {training_file}")
    assert finite_magnitudes.all(), f"{sfh_model}: {(~finite_magnitudes).sum()} models have non-finite magnitudes"
    display(training_models[DISPLAY_COLUMNS[sfh_model]].describe().T.round(3))
```

  Markdown `742607a1`, which replaces `## 4. Synthetic validation and test sets`:

```markdown
## 4. Synthetic validation and test sets
For each SFH model built in this run: independent draws with their own seeds. Each model takes the $r$ magnitude of a random DR3 galaxy at similar redshift, then gets S-PLUS errors, flux noise and non-detections. Mass-like labels are rescaled by the same magnitude offset.
```

  Code `9ff91c28`, which replaces the cell starting `catalogue_r = galaxy_catalogue["r_auto"].to_numpy()`:

```python
catalogue_r = galaxy_catalogue["r_auto"].to_numpy()
catalogue_z = galaxy_catalogue["zml"].to_numpy()
held_out_seeds = {"validation": VALIDATION_SEED, "test": TEST_SEED}

for sfh_model in models_to_build:
    for set_name, seed in held_out_seeds.items():
        output_file = set_files(sfh_model)[set_name]
        rng = np.random.default_rng(seed)
        models = simulate_in_parallel(N_HELD_OUT, seed, sfh_model)
        model_magnitudes = models[sm.BANDS].to_numpy()

        # 1. Realistic brightness at each redshift
        scaled_magnitudes, offsets = sm.scale_to_observed_r(model_magnitudes, models["redshift"].to_numpy(),
                                                            catalogue_r, catalogue_z, rng)
        models = sm.rescale_mass_labels(models, offsets)
        models["magnitude_offset"] = offsets

        # 2. S-PLUS noise and non-detections
        noisy_magnitudes, errors, detected_mask = sm.add_splus_noise(scaled_magnitudes, noise_table, rng)
        for band_index, band in enumerate(sm.BANDS):
            models[band] = noisy_magnitudes[:, band_index]
            models[f"e_{band}"] = errors[:, band_index]
        models.to_parquet(output_file, index=False)

        print(f"{output_file}: {len(models)} models, {detected_mask.sum(axis=1).mean():.1f} detected bands "
              f"on average (DR3: {all_detected.sum(axis=1).mean():.1f})")
        # No DR3 galaxy within the redshift window: these models have no realistic brightness
        n_without_offset = np.isnan(offsets).sum()
        print(f"{output_file}: {n_without_offset} models without a magnitude offset")
```

  Markdown `33ef32cd`, which replaces `## 5. Coverage test v2`:

```markdown
## 5. Coverage test v2
For each SFH model built in this run (`models_to_build`), against the training set it just wrote. A kept model keeps its stored coverage numbers and figures. For each galaxy, only models with $|z_{\rm model} - z_{\rm ml}| < \Delta z$ are compared. Over the detected bands, $\chi^2$ is minimised analytically over the magnitude offset (the model's mass). The galaxy is covered if the best $\chi^2$ has $p > 0.01$ for $n_{\rm detected} - 1$ degrees of freedom; the 0.05 mag error floor stands in for model error. Noisy synthetic galaxies should be covered ~99% of the time; that is the test's calibration check.
```

  Code `bd81d555`, which replaces the cell starting `training_magnitudes = training_models[sm.BANDS].to_numpy()`. The coverage runs only for the models built in this run, against the training set just written. A kept model, such as delayed-τ while its three sets exist, is not recomputed (spec decision 4). That saves a full coverage pass (2^19 models × about 9e4 galaxies) in the DPL job:

```python
def coverage_in_parallel(magnitudes, errors, detected_mask, z, dz, training_magnitudes, training_z):
    """Run sm.coverage_test against one training set, in parallel on N_JOBS workers.

    Inputs: magnitudes and errors [mag], arrays (n_galaxies, 12) in sm.BANDS order; detected_mask: bool array
            (n_galaxies, 12); z: redshifts (dimensionless), shape (n_galaxies,); dz: redshift half-windows
            (dimensionless), shape (n_galaxies,); training_magnitudes [mag], array (n_models, 12), and
            training_z (dimensionless), shape (n_models,): the noise-free training set.
    Returns: DataFrame with one row per galaxy (chi2_min, n_dof, p_value, covered, best_model).
    """
    blocks = np.array_split(np.arange(len(z)), COVERAGE_BLOCKS_PER_JOB * N_JOBS)
    tasks = []
    for block in blocks:
        task = delayed(sm.coverage_test)(magnitudes[block], errors[block], detected_mask[block], z[block],
                                         dz[block], training_magnitudes, training_z)
        tasks.append(task)
    tables = Parallel(n_jobs=N_JOBS)(tasks)
    return pd.concat(tables, ignore_index=True)


print(f"Old geometric test (K = 2): {OLD_GEOMETRIC_COVERAGE:.1%}")
if not models_to_build:
    print("No SFH model was built in this run, so the stored coverage numbers and figures stand")
# SFH model -> its S-PLUS table with the coverage columns, its validation coverage and its training magnitudes
coverage_results = {}
for sfh_model in models_to_build:
    files = set_files(sfh_model)
    training = pd.read_parquet(files["training"], columns=sm.BANDS + ["redshift"])
    training_magnitudes = training[sm.BANDS].to_numpy()
    training_z = training["redshift"].to_numpy()

    # 1. S-PLUS galaxies
    coverage = coverage_in_parallel(observed_magnitudes, observed_errors, detected, observed_z, observed_dz,
                                    training_magnitudes, training_z)
    model_splus = splus.copy()
    for column in coverage.columns:
        model_splus[column] = coverage[column].to_numpy()
    print(f"{sfh_model}: coverage v2 (p > {sm.P_MIN}): {model_splus['covered'].mean():.1%} of {len(model_splus)} "
          f"galaxies")
    # Galaxies outside the prior redshift range can only be uncovered, so also report the rest
    inside_prior = ~outside_prior
    covered_inside_prior = model_splus["covered"].to_numpy()[inside_prior]
    print(f"{sfh_model}: coverage v2 with zml inside Z_RANGE: {covered_inside_prior.mean():.1%} of "
          f"{inside_prior.sum()} galaxies")

    # 2. Calibration check on noisy synthetic galaxies
    validation = pd.read_parquet(files["validation"])
    validation_magnitudes = validation[sm.BANDS].to_numpy()
    validation_error_columns = []
    for band in sm.BANDS:
        validation_error_columns.append(f"e_{band}")
    validation_errors = validation[validation_error_columns].to_numpy()
    validation_detected = sm.detected_bands(validation_magnitudes, validation_errors)
    usable = validation_detected.sum(axis=1) >= sm.MIN_DETECTED_BANDS
    validation_z = validation["redshift"].to_numpy()[usable]
    validation_coverage = coverage_in_parallel(validation_magnitudes[usable], validation_errors[usable],
                                               validation_detected[usable], validation_z,
                                               sm.redshift_window(validation_z), training_magnitudes, training_z)
    print(f"{sfh_model}: synthetic validation galaxies covered: {validation_coverage['covered'].mean():.1%} "
          f"(expected about {1 - sm.P_MIN:.0%})")

    # 3. Coverage by redshift, magnitude and number of detected bands
    groupings = {
        "zml": pd.cut(model_splus["zml"], ZML_GROUP_EDGES),
        "r_auto": pd.cut(model_splus["r_auto"], R_GROUP_EDGES),
        "n_detected": model_splus["n_dof"] + 1,
    }
    for name, grouping in groupings.items():
        grouped = model_splus.groupby(grouping, observed=True)
        table = pd.DataFrame({
            "covered": grouped["covered"].mean(),
            "median p": grouped["p_value"].median(),
            "galaxies": grouped.size(),
        })
        print(f"{sfh_model}, by {name}:")
        display(table.round(3))

    coverage_results[sfh_model] = {"splus": model_splus, "validation": validation_coverage,
                                   "training_magnitudes": training_magnitudes}
```

  Markdown `fd06b7a0`, which replaces `## 6. Figures`:

```markdown
## 6. Figures
For each SFH model built in this run: colour-colour diagrams of the synthetic domain with covered and uncovered galaxies, and the p-value distributions of S-PLUS and synthetic galaxies (uniform if models and errors are right). The DPL figures carry the `_dblplaw` suffix.
```

  Code `64cba55b`, which replaces the cell starting `from matplotlib.colors import LinearSegmentedColormap`. The figure code is unchanged, wrapped in two functions and called once per model:

```python
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------------------
# Figure constants
# ---------------------------------------------------------------------------
POSTER = {"panel": "#1b1b3a", "grid": "#2a2a52", "edge": "#3d3d70", "text": "#d6d6ec",
          "muted": "#8a8ab8", "cyan": "#5fd0f5", "red": "#e0524f", "red_edge": "#8c2424",
          "covered": "#b8b8d0"}
DOMAIN_COLOURMAP = LinearSegmentedColormap.from_list("poster_blue", ["#22305f", "#2f5fa0", "#4aa6e0", POSTER["cyan"]])
POSTER_STYLE = {"figure.facecolor": POSTER["panel"], "axes.facecolor": POSTER["panel"],
                "axes.edgecolor": POSTER["edge"], "axes.labelcolor": POSTER["text"], "axes.labelsize": 20,
                "xtick.color": POSTER["muted"], "ytick.color": POSTER["muted"], "text.color": POSTER["text"],
                "font.size": 16, "legend.frameon": False}
PANELS = [("g", "i", "r", "z"), ("u", "g", "g", "r"), ("g", "r", "r", "i"), ("r", "i", "i", "z")]
N_GALAXIES_SHOWN = 1500                  # S-PLUS galaxies drawn per colour-colour panel (count)
DISPLAY_SEED = 3                         # random choice of the galaxies shown
PERCENTILE_RANGE = (0.5, 99.5)           # [%] of the shown galaxies' colours that sets the axis limits
X_MARGINS = np.array([-0.8, 1.7])        # [mag] added to the lower and upper x limit
Y_MARGINS = np.array([-0.8, 0.95])       # [mag] added to the lower and upper y limit
COLOUR_FIGURE_SIZE = (12, 10)            # [inch]
P_VALUE_FIGURE_SIZE = (10, 7)            # [inch]
LEGEND_RECT = (0, 0, 1, 0.96)            # figure area left to the panels, below the legend
COVERED_MARKER_SIZE = 8                  # [points^2]
UNCOVERED_MARKER_SIZE = 10               # [points^2]
LEGEND_MARKER_SIZE = 10                  # [points]
COVERED_ALPHA = 0.45
UNCOVERED_EDGE_WIDTH = 0.5               # [points]
HISTOGRAM_LINE_WIDTH = 2                 # [points]
FIGURE_DPI = 300
# ---------------------------------------------------------------------------


def colour(magnitudes, first_band, second_band):
    """Colour first_band - second_band of every object.

    Inputs: magnitudes [mag], array (n, 12) in sm.BANDS order; first_band and second_band: band names.
    Returns: colours [mag], shape (n,).
    """
    return magnitudes[:, sm.BANDS.index(first_band)] - magnitudes[:, sm.BANDS.index(second_band)]


def colour_colour_figure(covered, training_magnitudes, output_file):
    """Colour-colour diagrams of one training set with the covered and uncovered S-PLUS galaxies.

    Inputs: covered: bool array (n_galaxies,), the coverage of each S-PLUS galaxy (rows of observed_magnitudes);
            training_magnitudes [mag], array (n_models, 12) in sm.BANDS order; output_file: PNG file name.
    Returns: None; saves and shows the figure.
    """
    # 1. Random subsample of the S-PLUS galaxies, so the panels stay readable
    display_rng = np.random.default_rng(DISPLAY_SEED)
    shown = display_rng.random(len(covered)) < N_GALAXIES_SHOWN / len(covered)

    # 2. One panel per colour pair
    with plt.rc_context(POSTER_STYLE):
        figure, axes = plt.subplots(2, 2, figsize=COLOUR_FIGURE_SIZE)
        for axis, (band_a, band_b, band_c, band_d) in zip(axes.flat, PANELS):
            panel_detected = shown.copy()
            for band in (band_a, band_b, band_c, band_d):
                panel_detected = panel_detected & detected[:, sm.BANDS.index(band)]
            galaxy_x = colour(observed_magnitudes, band_a, band_b)[panel_detected]
            galaxy_y = colour(observed_magnitudes, band_c, band_d)[panel_detected]
            panel_covered = covered[panel_detected]

            x_percentiles = np.percentile(galaxy_x, PERCENTILE_RANGE)
            y_percentiles = np.percentile(galaxy_y, PERCENTILE_RANGE)
            x_limits = x_percentiles + X_MARGINS
            y_limits = y_percentiles + Y_MARGINS
            axis.hexbin(colour(training_magnitudes, band_a, band_b), colour(training_magnitudes, band_c, band_d),
                        extent=(*x_limits, *y_limits), bins="log", mincnt=1, cmap=DOMAIN_COLOURMAP, zorder=1)
            axis.scatter(galaxy_x[panel_covered], galaxy_y[panel_covered], s=COVERED_MARKER_SIZE,
                         c=POSTER["covered"], alpha=COVERED_ALPHA, lw=0, zorder=2)
            axis.scatter(galaxy_x[~panel_covered], galaxy_y[~panel_covered], s=UNCOVERED_MARKER_SIZE,
                         c=POSTER["red"], edgecolors=POSTER["red_edge"], lw=UNCOVERED_EDGE_WIDTH, zorder=3)
            axis.set(xlim=x_limits, ylim=y_limits, xlabel=f"${band_a} - {band_b}$",
                     ylabel=f"${band_c} - {band_d}$")

        handles = [
            plt.Line2D([], [], ls="", marker="h", ms=LEGEND_MARKER_SIZE, color=POSTER["cyan"],
                       label="Synthetic Domain"),
            plt.Line2D([], [], ls="", marker="o", ms=LEGEND_MARKER_SIZE, color=POSTER["covered"],
                       label="Covered Galaxies"),
            plt.Line2D([], [], ls="", marker="o", ms=LEGEND_MARKER_SIZE, color=POSTER["red"],
                       mec=POSTER["red_edge"], label="Uncovered Galaxies"),
        ]
        figure.legend(handles=handles, loc="upper center", ncol=3, labelcolor=POSTER["text"])
        figure.tight_layout(rect=LEGEND_RECT)
        figure.savefig(output_file, dpi=FIGURE_DPI, facecolor=figure.get_facecolor())
        plt.show()


def p_value_figure(splus_p_values, validation_p_values, output_file):
    """Distributions of log10 p of the S-PLUS and the synthetic validation galaxies against one training set.

    Inputs: splus_p_values, validation_p_values: coverage-test p-values (dimensionless; NaN rows are dropped);
            output_file: PNG file name.
    Returns: None; saves and shows the figure.
    """
    log_p_bins = np.linspace(*LOG_P_RANGE, N_LOG_P_BINS + 1)
    with plt.rc_context(POSTER_STYLE):
        figure, axis = plt.subplots(figsize=P_VALUE_FIGURE_SIZE)
        splus_p = np.clip(splus_p_values.dropna(), P_VALUE_FLOOR, 1)
        validation_p = np.clip(validation_p_values.dropna(), P_VALUE_FLOOR, 1)
        splus_log_p = np.log10(splus_p)
        validation_log_p = np.log10(validation_p)
        axis.hist(splus_log_p, bins=log_p_bins, density=True, histtype="step", lw=HISTOGRAM_LINE_WIDTH,
                  color=POSTER["covered"], label="S-PLUS DR3")
        axis.hist(validation_log_p, bins=log_p_bins, density=True, histtype="step", lw=HISTOGRAM_LINE_WIDTH,
                  color=POSTER["cyan"], label="Synthetic validation")
        axis.axvline(np.log10(sm.P_MIN), color=POSTER["red"], ls="--", lw=HISTOGRAM_LINE_WIDTH,
                     label=f"$p = {sm.P_MIN}$")
        axis.set(xlabel=r"$\log_{10} p$", ylabel="density")
        axis.grid(True, color=POSTER["grid"])
        axis.legend(labelcolor=POSTER["text"])
        figure.tight_layout()
        figure.savefig(output_file, dpi=FIGURE_DPI, facecolor=figure.get_facecolor())
        plt.show()


for sfh_model, results in coverage_results.items():
    print(f"{sfh_model}:")
    model_splus = results["splus"]
    colour_colour_figure(model_splus["covered"].to_numpy(), results["training_magnitudes"],
                         sm.model_file(COLOUR_FIGURE_FILE, sfh_model))
    p_value_figure(model_splus["p_value"], results["validation"]["p_value"],
                   sm.model_file(P_VALUE_FIGURE_FILE, sfh_model))
```

  The script, saved as `$SCRATCH/replace_cells.py`:

```python
"""Replace notebook cell sources by cell id; edited code cells lose their stale outputs.

Usage: python replace_cells.py <notebook.ipynb> <cells folder>; the folder holds <cell id>.py or <cell id>.md.
"""
import json
import os
import sys

notebook_path, cells_folder = sys.argv[1], sys.argv[2]
with open(notebook_path) as handle:
    notebook = json.load(handle)

new_sources = {}
for file_name in os.listdir(cells_folder):
    cell_id, _ = os.path.splitext(file_name)
    with open(os.path.join(cells_folder, file_name)) as handle:
        new_sources[cell_id] = handle.read().rstrip("\n")

replaced = set()
for cell in notebook["cells"]:
    if cell.get("id") in new_sources:
        cell["source"] = new_sources[cell["id"]].splitlines(keepends=True)
        if cell["cell_type"] == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
        replaced.add(cell["id"])
missing = set(new_sources) - replaced
assert not missing, f"cells not found: {missing}"

with open(notebook_path, "w") as handle:
    json.dump(notebook, handle, indent=1, ensure_ascii=False)
    handle.write("\n")
print(f"replaced {len(replaced)} cells")
```

  Run it:

```bash
.venv/bin/python $SCRATCH/replace_cells.py Synthetic_domain/training_set.ipynb $SCRATCH/training_set_cells
```

  Expected: `replaced 10 cells`.

  Notes on the cells:
  - joblib's `delayed` keeps its name. No variable is called `delayed`, and the SFH model is always the string `"delayed"`.
  - The delayed-τ sets exist, so with `REBUILD_SETS = False` only DPL is simulated. The delayed-τ held-out sets would be rebuilt with the same `np.random.default_rng(seed)` stream as before.
  - The coverage test (§5) and the figures (§6) loop over `models_to_build`; the figure cell loops over `coverage_results`, which §5 fills. A run that keeps the delayed-τ sets writes only `color_color_dblplaw.png` and `coverage_criterion_dblplaw.png`. `color_color.png` and `coverage_criterion.png` stay as they are.
  - Cell `19201d8e` is unchanged and still writes `splus_noise_model.csv`, which the held-out sets of every built model need. `sm.fit_noise_model` takes quantiles of the same DR3 catalogue and draws nothing at random, so the file's content does not change.

- [ ] **Step 6: Verify the notebook without running it.**

```bash
cd Synthetic_domain && ../.venv/bin/python -c "
import json
import nbformat
notebook = json.load(open('training_set.ipynb'))
for cell in notebook['cells']:
    if cell['cell_type'] == 'code':
        compile(''.join(cell['source']), cell.get('id', '?'), 'exec')
nbformat.validate(nbformat.read('training_set.ipynb', as_version=4))
print('all code cells compile; notebook valid')
"
grep -c "sm\.DIMS\|sm\.TRAINING_LABELS\|sm\.DRAWN_LABELS\|sm\.components\|sm\.PHOTOMETRIC_DIMS" training_set.ipynb
```

  Expected: `all code cells compile; notebook valid`, then `0`. `grep -c` exits with status 1 when it counts 0; that is the expected result.

- [ ] **Step 7: Run the whole suite.**

```bash
cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q
```

  Expected: `197 passed, 6 warnings` (about 4 minutes; the warnings are Task 2's Bagpipes overflow warnings). That is Task 2's 193 plus 4: `test_dblplaw_summary_labels_give_every_training_label`, and a `dblplaw` case for each of the three `test_simulated_*` tests. Later parts count from 197.

- [ ] **Step 8: Commit.**

```bash
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py \
        Synthetic_domain/tests/test_ml_data.py Synthetic_domain/tests/test_sfh_models.py Synthetic_domain/training_set.ipynb
git commit -m "$(cat <<'EOF'
feat: simulate_models and summary_labels take the SFH model; training_set.ipynb builds the sets of SFH_MODELS

Delayed-tau sets are kept when their three files exist (REBUILD_SETS = False); DPL outputs carry _dblplaw.
The coverage test and the figures run only for the models built in the run, so the delayed-tau coverage figures
are not recomputed. splus_noise_model.csv is rewritten with the same content, a deterministic fit of DR3.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
EOF
)"
```

  The DPL sets are generated later on the PC (`training_set.ipynb` as an `hljob` notebook job, CPU). Never generate them on the laptop.

### Task 4: Jobs and fits carry their SFH model (fitters, photometric dims per model, `fit_is_current` per folder, supersede, evidence error)

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py`. Line numbers below are today's (before Tasks 1–3); Tasks 1–3 move sections 1–2 out, so find every function by name.
  - imports at the top (the `from sfh_models import ...` statement added by Task 1);
  - section 5: `SPECTRAL_DIMS`, `training_unit_vector`, `spectral_components` (today ~738–785);
  - section 6: constants (~804–825), `fit_spectrum` docstring (~872–905), `run_multinest` (~907–953), `posterior_labels`, `thinned_samples`, `sample_label_row` (~956–1020), `photometric_components` (~1024–1033), `fit_photometry` (~1079–1120), `photometric_posterior_labels` (~1123–1155), new evidence helpers after it, `save_fit`, `load_fit`, `fit_is_current` (~1192–1248), new `job_model`, `fit_spectrum_job`, `fit_photometry_job`, `fit_and_save` (~1251–1351);
  - section 7: `export_job`, `load_job` (~1381–1423);
  - section 8: constants, `calibration_mocks_are_current`, `supersede_outdated_fits` (~1490–1594).
- Modify: `Synthetic_domain/run_spectral_fits.py` (module docstring only, lines 4–17). `slurm_fit.sh` already names `sfh_models.py` (Task 1), so it is not touched.
- Modify: `Synthetic_domain/spectroscopic_reference.ipynb`: the cell starting `# 1. Fits, calibration mocks, their job files and run logs` (section 4) and the four `sm.fit_is_current(metadata)` calls and the one `sm.photometric_components(...)` call (cells 14, 16, 18, 19, 20).
- Modify: `Synthetic_domain/remote/README.md` (the local supersede command, lines 182–183).
- Test: `Synthetic_domain/tests/test_synthetic_model.py` (sections "MultiNest fit and fit files", "Photometric fit", "Job files and the calibration test", "Superseding", "run_spectral_fits.py").

**Interfaces:**
- **Consumes (Tasks 1–3):** `sfh_models.MODELS`, `SfhModel.dims`, `.photometric_dims`, `.components(u, fixed_redshift=None)`, `.training_labels`, `.calibration_labels`; `SHARED_LABELS`; `model_file`; `sm.summary_labels(bagpipes_sfh, drawn, sfh_model)`.
- **Produces (constants):**
  - `SPECTRAL_MODEL = "delayed"`; `SPECTRAL_DIMS` = `MODELS["delayed"].photometric_dims + ["veldisp", "noise_scale"]` (the same 8 names as today);
  - `JOB_DEFAULT_SFH_MODEL = "delayed"`; `EVIDENCE_ERROR_SAMPLES = 200`;
  - `FIT_FOLDER_MODELS = {"spectral_fits": "delayed", "photometric_fits": "delayed", "photometric_fits_dblplaw": "dblplaw"}`.
- **Produces (functions, new or changed signature):**
  - `job_model(job) -> str`;
  - `training_unit_vector(cube_by_name, sfh_model) -> ndarray`;
  - `photometric_components(cube, redshift, sfh_model) -> (bagpipes_components, drawn)`;
  - `sample_label_row(model, bagpipes_components, drawn, sfh_model) -> dict`;
  - `thinned_samples(unit_cube_samples, max_samples=MAX_POSTERIOR_SAMPLES) -> ndarray`;
  - `fit_photometry(magnitudes, errors, redshift, sfh_model, output_basename, max_iterations, n_live_points=N_LIVE_POINTS) -> dict`;
  - `photometric_posterior_labels(unit_cube_samples, redshift, magnitudes, errors, rng, sfh_model) -> DataFrame`;
  - `run_multinest(...)` and so `fit_spectrum` / `fit_photometry` return `log_evidence_error` (pymultinest's `logZerr`) too;
  - `evidence_error(log_likelihoods, log_evidence, n_live_points) -> float`;
  - `photometric_log_likelihoods(unit_cube_samples, redshift, magnitudes, errors, sfh_model, max_samples=EVIDENCE_ERROR_SAMPLES) -> ndarray`;
  - `load_unit_cube_samples(path) -> ndarray`;
  - `fit_evidence(fit_file, job) -> (log_evidence, log_evidence_error)`;
  - `fit_is_current(metadata, sfh_model) -> bool`;
  - `calibration_mocks_are_current(mock_file, sfh_model) -> bool`;
  - `supersede_outdated_fits(fits_directory, sfh_model) -> list` (raises `ValueError` when the folder is in `FIT_FOLDER_MODELS` with another model).
- **Unchanged signatures, new behaviour:**
  - `export_job` writes `sfh_model` (always) and `truth_model` (when the job has one); `load_job` returns `sfh_model` (`"delayed"` for a file without it) and `truth_model` when stored;
  - `fit_and_save` fits with the job's model, refuses a same-kind fit of another model, and stores `sfh_model`, `log_evidence_error`, `n_live_points` and (cross-fits) `truth_model` in the metadata;
  - `save_fit` stores `metadata["sfh_model"]`, or `"delayed"` when the caller gives none;
  - `fit_spectrum_job` raises `ValueError` for a job whose model is not `"delayed"` (spectral fits are delayed-τ only);
  - `posterior_labels(unit_cube_samples, redshift)` and `spectral_components(cube, redshift)` keep their signatures and stay delayed-τ.

**Evidence uncertainty of the fits that already exist.** pymultinest's `solve` returns `logZerr` (the "Nested Sampling Global Log-Evidence error" of MultiNest, √(H / n_live)), so every new fit stores it. The 414 delayed-τ photometric fits already in `photometric_fits/` have no such key, and their MultiNest folders stayed in the PC job (only the `.h5` files are fetched). For those, `fit_evidence` recomputes the same quantity from the stored fit: H = mean of ln L over the stored equal-weight posterior samples − ln Z, where ln L is the fit's own likelihood, evaluated again at `EVIDENCE_ERROR_SAMPLES` evenly spaced samples. `n_live_points` is `sm.N_LIVE_POINTS` (300) for every old fit, because `submit_fits.sh` never set `--n-live`. Cost: about 3 ms per sample, so about 4 minutes for the 414 old fits, inside the notebook-2 job on the PC. A MultiNest-gated test checks that the recomputed error agrees with MultiNest's own to within a factor of 2.

**Note on `run_spectral_fits.py`:** the worker passes every job through unchanged. It calls `sm.load_job` and then `sm.fit_and_save`, and both now carry `sfh_model` and `truth_model`. `interleave_calibration_and_galaxies` treats `cross_calib_<k>` ids as galaxies, which only changes the order of the runs. Only the docstring changes, and a DPL job is tested end to end.

- [ ] **Step 1: Update the test helpers and the tests whose signatures change, and write the new failing tests**

**Rule for the existing tests (spec amendment 4).** An existing delayed-τ test changes in two ways only:
- it passes the new `sfh_model` argument (`"delayed"`);
- it expects the new metadata keys (`log_evidence_error`).

The delayed-only aliases of Task 1 are replaced by the registry entry they alias, which has the same value: `sm.PHOTOMETRIC_DIMS` becomes `sm.MODELS["delayed"].photometric_dims`, `sm.DIMS` becomes `sm.MODELS["delayed"].dims`, `sm.components` becomes `sm.MODELS["delayed"].components`, `sm.TRAINING_LABELS` becomes `sm.MODELS["delayed"].training_labels`, and `sm.SFH_MODEL` becomes `sm.MODELS["delayed"].name` or `"delayed"`. Every delayed-τ assertion stays as it is. The DPL cases go into new tests. Do not loosen or drop an old assertion.

In `tests/test_synthetic_model.py`, replace these helpers and tests, each with the function of the same name.

```python
def noiseless_photometry(cube, redshift, sfh_model="delayed"):
    """Model AB magnitudes [mag] of one photometric-fit point of sfh_model, and the unscaled stellar mass [log M_sun]."""
    import bagpipes as pipes

    bagpipes_components, _ = sm.photometric_components(cube, redshift, sfh_model)
    model = pipes.model_galaxy(bagpipes_components, filt_list=sm.FILTER_FILES)
    return sm.ab_magnitudes(model.photometry), model.sfh.stellar_mass


def test_photometric_components_are_the_training_prior_without_the_redshift():
    rng = np.random.default_rng(20)
    delayed = sm.MODELS["delayed"]
    cube = rng.random(len(delayed.photometric_dims))
    full_unit_vector = np.insert(cube, delayed.dims.index("z"), 0.9)   # any redshift draw: fixed_redshift replaces it
    expected_components, expected_drawn = delayed.components(full_unit_vector, fixed_redshift=0.12)

    bagpipes_components, drawn = sm.photometric_components(cube, 0.12, "delayed")

    expected_dims = list(delayed.dims)
    expected_dims.remove("z")
    assert delayed.photometric_dims == expected_dims
    assert bagpipes_components == expected_components
    assert drawn == expected_drawn
    assert "veldisp" not in bagpipes_components


def test_photometric_posterior_labels_shift_the_masses_by_the_offset():
    """A galaxy 1.5 mag fainter than its model (tiny errors) gets log_mstar lower by 0.4 x 1.5 = 0.6 dex."""
    cube = np.full(len(sm.MODELS["delayed"].photometric_dims), 0.5)
    model_magnitudes, unscaled_log_mass = noiseless_photometry(cube, 0.1)
    offset = 1.5
    errors = np.full(len(sm.BANDS), 1e-4)

    posterior = sm.photometric_posterior_labels(cube[np.newaxis, :], 0.1, model_magnitudes + offset, errors,
                                                np.random.default_rng(21), "delayed")

    assert np.isclose(posterior["magnitude_offset"].iloc[0], offset, atol=1e-3)
    assert np.isclose(posterior["log_mstar"].iloc[0], unscaled_log_mass - 0.4 * offset, atol=1e-3)
    assert np.isclose(posterior["log_mass_formed"].iloc[0], sm.LOG_MASS_FORMED - 0.4 * offset, atol=1e-3)
    for column in (*sm.MODELS["delayed"].training_labels, "model_r", "magnitude_offset"):
        assert column in posterior.columns, column


def test_photometric_posterior_labels_draw_the_offset_from_its_conditional():
    """Repeated copies of one point get offsets with the analytic mean and width (errors 0.1 mag, 12 bands)."""
    n_copies = 200
    cube = np.full(len(sm.MODELS["delayed"].photometric_dims), 0.5)
    model_magnitudes, _ = noiseless_photometry(cube, 0.1)
    errors = np.full(len(sm.BANDS), 0.1)
    samples = np.tile(cube, (n_copies, 1))

    posterior = sm.photometric_posterior_labels(samples, 0.1, model_magnitudes + 1.0, errors,
                                                np.random.default_rng(22), "delayed")

    offset_sigma = 0.1 / np.sqrt(len(sm.BANDS))
    assert len(posterior) == n_copies
    assert np.isclose(posterior["magnitude_offset"].mean(), 1.0, atol=3 * offset_sigma / np.sqrt(n_copies))
    assert np.isclose(posterior["magnitude_offset"].std(), offset_sigma, rtol=0.2)


def test_fit_photometry_needs_enough_detected_bands():
    magnitudes = np.full(len(sm.BANDS), 99.0)
    magnitudes[:sm.MIN_DETECTED_BANDS - 1] = 20.0
    errors = np.full(len(sm.BANDS), 0.05)
    with pytest.raises(ValueError, match="detected bands"):
        sm.fit_photometry(magnitudes, errors, 0.1, "delayed", "unused/", max_iterations=10, n_live_points=25)


@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_fit_photometry_smoke_run_on_a_noiseless_mock(tmp_path):
    """A short MultiNest run on noiseless model photometry returns unit-cube samples and run details."""
    cube = np.random.default_rng(23).random(len(sm.MODELS["delayed"].photometric_dims))
    magnitudes, _ = noiseless_photometry(cube, 0.1)
    errors = np.full(len(sm.BANDS), 0.05)
    magnitudes[0] = 99.0                     # one undetected band: the fit runs on the other eleven

    fit = sm.fit_photometry(magnitudes, errors, 0.1, "delayed", str(tmp_path) + "/", max_iterations=200,
                            n_live_points=25)

    assert set(fit) == {"samples", "log_evidence", "log_evidence_error", "n_iterations", "converged", "run_time"}
    assert fit["samples"].shape[1] == len(sm.MODELS["delayed"].photometric_dims)
    assert np.all((fit["samples"] >= 0) & (fit["samples"] <= 1))
    assert np.isfinite(fit["log_evidence"])
    assert np.isfinite(fit["log_evidence_error"]) and fit["log_evidence_error"] > 0
    assert fit["n_iterations"] <= 200


def photometry_job(job_id, rng, sfh_model="delayed"):
    """A photometry job of sfh_model at z = 0.1: noiseless magnitudes of a random point, 0.05 mag errors, u undetected.

    Inputs: job_id: identifier (str); rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: job dict with job_id, kind, sfh_model, magnitudes [mag], errors [mag] and redshift.
    """
    cube = rng.random(len(sm.MODELS[sfh_model].photometric_dims))
    magnitudes, _ = noiseless_photometry(cube, 0.1, sfh_model)
    magnitudes[0] = 99.0
    return {"job_id": job_id, "kind": "photometry", "sfh_model": sfh_model, "magnitudes": magnitudes,
            "errors": np.full(len(sm.BANDS), 0.05), "redshift": 0.1}


def fake_fit(n_dims):
    """What a fitter returns, without running MultiNest: two unit-cube samples of n_dims coordinates."""
    return {"samples": np.full((2, n_dims), 0.5), "log_evidence": -1.0, "log_evidence_error": 0.1,
            "n_iterations": 7, "converged": True, "run_time": 3.0}


@pytest.fixture
def fake_fitters(monkeypatch):
    """Replace both fitters and both labellers; returns the list of fitters called, in order."""
    called = []

    def fake_fit_spectrum(spectrum, redshift, output_basename, max_iterations, n_live_points):
        called.append("spectrum")
        return fake_fit(len(sm.SPECTRAL_DIMS))

    def fake_fit_photometry(magnitudes, errors, redshift, sfh_model, output_basename, max_iterations,
                            n_live_points):
        called.append("photometry")
        return fake_fit(len(sm.MODELS[sfh_model].photometric_dims))

    def fake_spectral_labels(unit_cube_samples, redshift):
        return pd.DataFrame({"mw_age": [1.0, 2.0]})

    def fake_photometric_labels(unit_cube_samples, redshift, magnitudes, errors, rng, sfh_model):
        return pd.DataFrame({"mw_age": [3.0, 4.0], "magnitude_offset": [1.0, 1.0]})

    monkeypatch.setattr(sm, "fit_spectrum", fake_fit_spectrum)
    monkeypatch.setattr(sm, "fit_photometry", fake_fit_photometry)
    monkeypatch.setattr(sm, "posterior_labels", fake_spectral_labels)
    monkeypatch.setattr(sm, "photometric_posterior_labels", fake_photometric_labels)
    return called


def write_fit(path, max_iterations=sm.MAX_ITERATIONS, sfh_model="delayed"):
    """A small fit file; sfh_model=None writes it as a 16-bin fit did, without the key."""
    metadata = {"converged": True, "run_time": 1.0, "n_iterations": 5, "max_iterations": max_iterations}
    sm.save_fit(str(path), pd.DataFrame({"log_age": [0.5]}), np.zeros((1, len(sm.SPECTRAL_DIMS))), metadata)
    with h5py.File(path, "a") as fit_file:
        if sfh_model is None:
            del fit_file.attrs["sfh_model"]
        else:
            fit_file.attrs["sfh_model"] = sfh_model


def test_save_and_load_fit_round_trip_the_sfh_model(tmp_path):
    """save_fit writes sfh_model on every fit, even when the caller's metadata lacks it."""
    path = str(tmp_path / "round_trip.h5")
    posterior = pd.DataFrame({"log_age": [0.1, 0.2]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "max_iterations": sm.MAX_ITERATIONS}
    sm.save_fit(path, posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), metadata)

    loaded_posterior, loaded_metadata = sm.load_fit(path)

    assert loaded_metadata["sfh_model"] == sm.MODELS["delayed"].name
    assert sm.fit_is_current(loaded_metadata, "delayed")
    assert list(loaded_posterior["log_age"]) == [0.1, 0.2]
    assert "sfh_model" not in metadata    # the caller's dict is left alone


def test_fit_is_current_needs_the_delayed_model_and_the_production_iteration_limit():
    """Fits without sfh_model (16-bin era), of another model or with a smaller limit are not current."""
    assert not sm.fit_is_current({"converged": True, "max_iterations": sm.MAX_ITERATIONS}, "delayed")
    assert sm.fit_is_current({"sfh_model": "delayed", "max_iterations": sm.MAX_ITERATIONS}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "nonparametric", "max_iterations": sm.MAX_ITERATIONS}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "delayed", "max_iterations": sm.MAX_ITERATIONS - 1}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "delayed"}, "delayed")
    assert sm.MODELS["delayed"].name == "delayed"


def test_calibration_mocks_without_truths_are_not_current(tmp_path):
    """A mock file holding no mock (e.g. every photometric mock lacked detected bands) is redrawn, not reused."""
    joblib.dump(([], {}), tmp_path / sm.CALIBRATION_MOCK_FILE)
    assert not sm.calibration_mocks_are_current(str(tmp_path / sm.CALIBRATION_MOCK_FILE), "delayed")
    write_mock_file(tmp_path / "current.joblib", sm.MODELS["delayed"].training_labels)
    assert sm.calibration_mocks_are_current(str(tmp_path / "current.joblib"), "delayed")
```

Next to these kept tests, add their DPL counterparts. Each goes at the place its comment names:

```python
# --- in "Photometric fit", after test_photometric_components_are_the_training_prior_without_the_redshift ---
def test_dblplaw_photometric_components_are_its_training_prior_without_the_redshift():
    """The DPL fit point maps to the DPL training prior at the fixed redshift, with dblplaw components only."""
    dblplaw = sm.MODELS["dblplaw"]
    cube = np.random.default_rng(46).random(len(dblplaw.photometric_dims))
    full_unit_vector = np.insert(cube, dblplaw.dims.index("z"), 0.9)   # any redshift draw: fixed_redshift replaces it
    expected_components, expected_drawn = dblplaw.components(full_unit_vector, fixed_redshift=0.12)

    bagpipes_components, drawn = sm.photometric_components(cube, 0.12, "dblplaw")

    expected_dims = list(dblplaw.dims)
    expected_dims.remove("z")
    assert dblplaw.photometric_dims == expected_dims
    assert bagpipes_components == expected_components
    assert drawn == expected_drawn
    assert "dblplaw" in bagpipes_components
    assert "delayed" not in bagpipes_components
    assert "veldisp" not in bagpipes_components


# --- in "MultiNest fit and fit files", after test_save_and_load_fit_round_trip_the_sfh_model ---
def test_save_fit_keeps_the_sfh_model_the_caller_gives(tmp_path):
    """A DPL fit is stored as dblplaw, so fit_is_current keeps it only in a DPL folder."""
    path = str(tmp_path / "dblplaw.h5")
    posterior = pd.DataFrame({"log_turnover_lookback": [0.1, 0.2]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "max_iterations": sm.MAX_ITERATIONS,
                "sfh_model": "dblplaw"}
    sm.save_fit(path, posterior, np.zeros((2, len(sm.MODELS["dblplaw"].photometric_dims))), metadata)

    _, loaded_metadata = sm.load_fit(path)

    assert loaded_metadata["sfh_model"] == "dblplaw"
    assert sm.fit_is_current(loaded_metadata, "dblplaw")
    assert not sm.fit_is_current(loaded_metadata, "delayed")


def test_fit_is_current_compares_the_fit_model_with_the_folder_model():
    """A fit is current only in a folder of its own model; FIT_FOLDER_MODELS names each folder's model."""
    full_run = sm.MAX_ITERATIONS
    assert not sm.fit_is_current({"sfh_model": "delayed", "max_iterations": full_run}, "dblplaw")
    assert sm.fit_is_current({"sfh_model": "dblplaw", "max_iterations": full_run}, "dblplaw")
    assert not sm.fit_is_current({"sfh_model": "dblplaw", "max_iterations": full_run}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "dblplaw", "max_iterations": full_run - 1}, "dblplaw")
    assert sm.FIT_FOLDER_MODELS == {"spectral_fits": "delayed", "photometric_fits": "delayed",
                                    "photometric_fits_dblplaw": "dblplaw"}


# --- in "Superseding", after test_calibration_mocks_without_truths_are_not_current ---
def test_calibration_mocks_are_current_only_with_their_model_labels(tmp_path):
    """A mock file is current only in a folder of the model whose calibration labels its truths carry."""
    write_mock_file(tmp_path / "delayed.joblib", sm.MODELS["delayed"].training_labels)
    write_mock_file(tmp_path / "dblplaw.joblib", sm.MODELS["dblplaw"].training_labels)
    assert not sm.calibration_mocks_are_current(str(tmp_path / "delayed.joblib"), "dblplaw")
    assert sm.calibration_mocks_are_current(str(tmp_path / "dblplaw.joblib"), "dblplaw")
    assert not sm.calibration_mocks_are_current(str(tmp_path / "dblplaw.joblib"), "delayed")
```

The existing supersede tests call `sm.supersede_outdated_fits` with one argument. Give them the folder's model:

```bash
cd Synthetic_domain
sed -i 's|sm.supersede_outdated_fits(str(tmp_path))|sm.supersede_outdated_fits(str(tmp_path), "delayed")|g; s|sm.supersede_outdated_fits(str(tmp_path / "missing"))|sm.supersede_outdated_fits(str(tmp_path / "missing"), "delayed")|' tests/test_synthetic_model.py
grep -c 'supersede_outdated_fits(str(tmp_path[^)]*), "delayed")' tests/test_synthetic_model.py   # expect 6
```

Then add the other new tests. Each goes at the place its comment names:

```python
# --- in "MultiNest fit and fit files", after test_save_and_load_fit_round_trip_the_sfh_model ---
def test_a_job_without_sfh_model_is_delayed():
    """Every job exported before the DPL model has no sfh_model: it is a delayed-tau job."""
    assert sm.job_model({"job_id": "123"}) == "delayed"
    assert sm.job_model({"job_id": "cross_calib_000", "sfh_model": "dblplaw"}) == "dblplaw"


def test_spectral_dims_stay_the_delayed_photometric_dims_plus_the_nuisances():
    """Spectral fits are delayed-tau only (DPL spec, decision 1); the DPL photometric fit has 7 dimensions."""
    assert sm.SPECTRAL_MODEL == "delayed"
    assert sm.SPECTRAL_DIMS == [*sm.MODELS["delayed"].photometric_dims, "veldisp", "noise_scale"]
    assert len(sm.MODELS["dblplaw"].photometric_dims) == 7


# --- in "Photometric fit", after test_fit_photometry_smoke_run_on_a_noiseless_mock ---
def test_evidence_error_is_the_square_root_of_the_information_over_the_live_points():
    """H = mean ln L - ln Z = -2 - (-5) = 3 nats, so the error is sqrt(3 / 300) = 0.1; no information gives 0."""
    assert np.isclose(sm.evidence_error(np.array([-1.0, -3.0]), -5.0, 300), 0.1)
    assert sm.evidence_error(np.array([-5.0, -5.0]), -5.0, 300) == 0.0


def test_photometric_log_likelihoods_are_the_fit_likelihood():
    """At each sample, the recomputed ln L equals the fit's own likelihood of that point (DPL model)."""
    samples = np.random.default_rng(45).random((3, len(sm.MODELS["dblplaw"].photometric_dims)))
    model_magnitudes, _ = noiseless_photometry(samples[0], 0.1, "dblplaw")
    observed = model_magnitudes + 1.0            # sample 0 has perfect colours
    errors = np.full(len(sm.BANDS), 0.05)
    detected = sm.detected_bands(observed, errors)

    log_likelihoods = sm.photometric_log_likelihoods(samples, 0.1, observed, errors, "dblplaw")

    assert len(log_likelihoods) == 3
    for sample_index, cube in enumerate(samples):
        sample_magnitudes, _ = noiseless_photometry(cube, 0.1, "dblplaw")
        expected = sm.safe_photometric_log_likelihood(sample_magnitudes, observed, errors, detected)
        assert np.isclose(log_likelihoods[sample_index], expected)
    perfect_fit = 0.5 * np.log(2 * np.pi / np.sum(1 / errors**2))    # chi2_min = 0
    assert np.isclose(log_likelihoods[0], perfect_fit)


def test_fit_evidence_reads_the_stored_error_or_recomputes_it(tmp_path, monkeypatch):
    """A fit with log_evidence_error returns it; an older fit gets sqrt(H / N_LIVE_POINTS) from its samples."""
    posterior = pd.DataFrame({"log_age": [0.5]})
    metadata = {"converged": True, "run_time": 1.0, "n_iterations": 5, "max_iterations": sm.MAX_ITERATIONS,
                "log_evidence": -5.0, "redshift": 0.1, "kind": "photometry"}
    samples = np.zeros((1, len(sm.MODELS["delayed"].photometric_dims)))
    sm.save_fit(str(tmp_path / "new.h5"), posterior, samples, dict(metadata, log_evidence_error=0.25))
    sm.save_fit(str(tmp_path / "old.h5"), posterior, samples, metadata)
    job = {"magnitudes": np.full(12, 20.0), "errors": np.full(12, 0.05), "redshift": 0.1}

    def fixed_log_likelihoods(unit_cube_samples, redshift, magnitudes, errors, sfh_model):
        return np.array([-1.0, -3.0])

    monkeypatch.setattr(sm, "photometric_log_likelihoods", fixed_log_likelihoods)

    assert sm.fit_evidence(str(tmp_path / "new.h5"), job) == (-5.0, 0.25)
    log_evidence, log_evidence_error = sm.fit_evidence(str(tmp_path / "old.h5"), job)
    assert log_evidence == -5.0
    assert np.isclose(log_evidence_error, np.sqrt(3 / sm.N_LIVE_POINTS))


# --- in "Job files and the calibration test", after test_job_file_written_before_kinds_loads_as_a_spectrum_job ---
def test_export_and_load_job_round_trip_the_models(tmp_path):
    """A cross-fit job keeps its model and true model; a job file written before the DPL model is delayed-tau."""
    job = {"job_id": "cross_calib_003", "kind": "photometry", "sfh_model": "dblplaw", "truth_model": "delayed",
           "magnitudes": np.full(12, 20.0), "errors": np.full(12, 0.05), "redshift": 0.1}
    loaded = sm.load_job(sm.export_job(job, str(tmp_path / "jobs")))
    assert loaded["sfh_model"] == "dblplaw"
    assert loaded["truth_model"] == "delayed"

    old_path = str(tmp_path / "123.npz")
    np.savez(old_path, job_id="123", kind="photometry", redshift=0.1, magnitudes=np.full(12, 20.0),
             errors=np.full(12, 0.05))
    old_job = sm.load_job(old_path)
    assert old_job["sfh_model"] == "delayed"
    assert "truth_model" not in old_job


def test_fit_and_save_fits_each_photometry_job_with_its_own_model(tmp_path, monkeypatch):
    """A DPL job is fitted and labelled as DPL; a job without sfh_model as delayed-tau; both store their model."""
    fitted_models = []
    labelled_models = []

    def capturing_fit_photometry(magnitudes, errors, redshift, sfh_model, output_basename, max_iterations,
                                 n_live_points):
        fitted_models.append(sfh_model)
        return fake_fit(len(sm.MODELS[sfh_model].photometric_dims))

    def capturing_labels(unit_cube_samples, redshift, magnitudes, errors, rng, sfh_model):
        labelled_models.append(sfh_model)
        return pd.DataFrame({"mw_age": [3.0, 4.0]})

    monkeypatch.setattr(sm, "fit_photometry", capturing_fit_photometry)
    monkeypatch.setattr(sm, "photometric_posterior_labels", capturing_labels)
    photometry = {"kind": "photometry", "magnitudes": np.full(12, 20.0), "errors": np.full(12, 0.05),
                  "redshift": 0.1}
    dblplaw_job = dict(photometry, job_id="cross_calib_000", sfh_model="dblplaw", truth_model="delayed")
    old_job = dict(photometry, job_id="123")      # no sfh_model, as every job exported before the DPL model

    for job in (dblplaw_job, old_job):
        summary = sm.fit_and_save(job, str(tmp_path), max_iterations=10, n_live_points=25)
        assert summary["status"] == "fitted", summary["error"]

    assert fitted_models == ["dblplaw", "delayed"]
    assert labelled_models == ["dblplaw", "delayed"]
    _, dblplaw_metadata = sm.load_fit(str(tmp_path / "cross_calib_000.h5"))
    assert dblplaw_metadata["sfh_model"] == "dblplaw"
    assert dblplaw_metadata["truth_model"] == "delayed"
    assert dblplaw_metadata["log_evidence_error"] == 0.1
    assert dblplaw_metadata["n_live_points"] == 25
    _, old_metadata = sm.load_fit(str(tmp_path / "123.h5"))
    assert old_metadata["sfh_model"] == "delayed"
    assert "truth_model" not in old_metadata


def test_fit_and_save_fails_a_spectral_job_of_the_dblplaw_model(tmp_path, fake_fitters):
    """Spectral fits are delayed-tau only (DPL spec, decision 1)."""
    job = {"job_id": "spectral_dpl", "kind": "spectrum", "sfh_model": "dblplaw", "redshift": 0.1,
           "spectrum": {"wavelength": np.arange(4000.0, 4010.0), "flux": np.ones(10), "ivar": np.ones(10)}}

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=10, n_live_points=25)

    assert summary["status"] == "failed"
    assert "delayed only" in summary["error"]
    assert fake_fitters == []


@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_fit_and_save_dblplaw_photometry_smoke_run(tmp_path):
    """A short MultiNest run of a DPL photometry job (noiseless mock) is saved with DPL labels and its model."""
    job = photometry_job("dpl_smoke", np.random.default_rng(38), "dblplaw")
    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=200, n_live_points=25)

    assert summary["status"] == "fitted", summary["error"]
    posterior, metadata = sm.load_fit(str(tmp_path / "dpl_smoke.h5"))
    assert metadata["sfh_model"] == "dblplaw"
    assert np.isfinite(metadata["log_evidence_error"])
    for column in (*sm.MODELS["dblplaw"].training_labels, "magnitude_offset"):
        assert column in posterior.columns, column
        assert np.isfinite(posterior[column]).all(), column


@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_recomputed_evidence_error_agrees_with_multinest(tmp_path):
    """sqrt(H / n_live) from the stored samples matches MultiNest's own error to a factor of 2."""
    job = photometry_job("evidence", np.random.default_rng(37))
    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=5000, n_live_points=50)
    assert summary["status"] == "fitted", summary["error"]
    if not summary["converged"]:
        pytest.skip("the short run did not converge, so its evidence is unfinished")
    fit_file = str(tmp_path / "evidence.h5")
    _, metadata = sm.load_fit(fit_file)

    samples = sm.load_unit_cube_samples(fit_file)
    log_likelihoods = sm.photometric_log_likelihoods(samples, 0.1, job["magnitudes"], job["errors"], "delayed")
    recomputed = sm.evidence_error(log_likelihoods, metadata["log_evidence"], 50)

    assert 0.5 < recomputed / metadata["log_evidence_error"] < 2.0


# --- in "Superseding", after test_supersede_moves_jobs_of_mocks_whose_file_is_missing ---
def test_supersede_judges_each_folder_by_its_own_model(tmp_path):
    """A DPL fit stays in photometric_fits_dblplaw and leaves photometric_fits; a delayed-tau fit the reverse."""
    dblplaw_folder = tmp_path / "photometric_fits_dblplaw"
    delayed_folder = tmp_path / "photometric_fits"
    for folder in (dblplaw_folder, delayed_folder):
        folder.mkdir()
        write_fit(folder / "dpl.h5", sfh_model="dblplaw")
        write_fit(folder / "delayed.h5", sfh_model="delayed")
    write_mock_file(dblplaw_folder / sm.CALIBRATION_MOCK_FILE, sm.MODELS["dblplaw"].training_labels)

    moved_dblplaw = sm.supersede_outdated_fits(str(dblplaw_folder), "dblplaw")
    moved_delayed = sm.supersede_outdated_fits(str(delayed_folder), "delayed")

    assert moved_dblplaw == ["delayed.h5"]
    assert moved_delayed == ["dpl.h5"]
    assert (dblplaw_folder / "dpl.h5").exists()
    assert (dblplaw_folder / sm.CALIBRATION_MOCK_FILE).exists()    # DPL truths are current in the DPL folder
    assert (delayed_folder / "delayed.h5").exists()


def test_supersede_refuses_a_model_other_than_the_folder_model(tmp_path):
    folder = tmp_path / "photometric_fits_dblplaw"
    folder.mkdir()
    write_fit(folder / "dpl.h5", sfh_model="dblplaw")

    with pytest.raises(ValueError, match="dblplaw"):
        sm.supersede_outdated_fits(str(folder), "delayed")
    assert (folder / "dpl.h5").exists()


# --- in "run_spectral_fits.py", after test_run_script_fits_an_exported_photometry_job ---
@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_run_script_fits_an_exported_dblplaw_photometry_job(tmp_path):
    """The worker passes the job's model and true model through to the fit file."""
    job = photometry_job("cross_calib_000", np.random.default_rng(39), "dblplaw")
    job["truth_model"] = "delayed"
    sm.export_job(job, str(tmp_path / "jobs"))

    result = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", str(tmp_path),
                        "--max-iterations", "200", "--n-live", "25", "--workers", "1")

    assert result.returncode == 0, result.stderr
    _, metadata = sm.load_fit(str(tmp_path / "cross_calib_000.h5"))
    assert metadata["sfh_model"] == "dblplaw"
    assert metadata["truth_model"] == "delayed"
```

- [ ] **Step 2: Run the tests and see them fail**

Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py -q -k "job_model or without_sfh_model or spectral_dims_stay or round_trip_the_models or own_model or dblplaw or fit_is_current or save_fit or save_and_load or photometric_components or photometric_posterior or fit_photometry or mocks_are_current or without_truths or supersede or evidence or log_likelihoods"`

Expected: FAIL. The failures include `AttributeError: module 'synthetic_model' has no attribute 'job_model'` (also `FIT_FOLDER_MODELS`, `SPECTRAL_MODEL`, `evidence_error`, `fit_evidence`), `TypeError: fit_is_current() takes 1 positional argument but 2 were given`, `TypeError: photometric_components() takes 2 positional arguments but 3 were given` and `TypeError: supersede_outdated_fits() takes 1 positional argument but 2 were given`.

- [ ] **Step 3: Implement in `synthetic_model.py`**

(a) **Imports.** Task 1 added a `from sfh_models import (...)` statement at the top. Make sure it also imports `MODELS`, `SHARED_LABELS` and `model_file`, adding any missing names to that statement.

(b) **Section 5.** Replace the `SPECTRAL_DIMS` loop, `training_unit_vector` and `spectral_components` with:

```python
# Spectral fits are delayed-tau only (DPL spec, decision 1): the DPL model is fitted to the photometry alone
SPECTRAL_MODEL = "delayed"

# Redshift is fixed to the spectroscopic value, so it leaves the list; veldisp and noise_scale join it
SPECTRAL_DIMS = list(MODELS[SPECTRAL_MODEL].photometric_dims)
SPECTRAL_DIMS.append("veldisp")
SPECTRAL_DIMS.append("noise_scale")


def training_unit_vector(cube_by_name, sfh_model):
    """Unit-cube vector in a model's dims order for a fit at fixed redshift, from the fit's own named coordinates.

    Inputs: cube_by_name: dict of unit-cube values (dimensionless, in [0, 1)) holding every dimension of the model
            but "z"; other names (the fit's nuisance parameters) are ignored; sfh_model: "delayed" or "dblplaw".
    Returns: array of len(MODELS[sfh_model].dims) values (dimensionless); the "z" slot holds IGNORED_REDSHIFT_DRAW.
    """
    dims = MODELS[sfh_model].dims
    unit_vector = np.zeros(len(dims))
    for dimension_index, dimension in enumerate(dims):
        if dimension == "z":
            # components() ignores this value, because fixed_redshift replaces the drawn redshift
            unit_vector[dimension_index] = IGNORED_REDSHIFT_DRAW
        else:
            unit_vector[dimension_index] = cube_by_name[dimension]
    return unit_vector


def spectral_components(cube, redshift):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit unit-cube point.

    The SFH, metallicity, dust and nebular part is the delayed-tau training prior (MODELS[SPECTRAL_MODEL]) at the
    fixed spectroscopic redshift, which also sets the oldest allowed age; veldisp and noise_scale are the
    spectral nuisance parameters.
    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as SPECTRAL_DIMS;
            redshift: spectroscopic redshift (dimensionless).
    Returns: (bagpipes_components, drawn, noise_scale). The components are the delayed-tau model's plus
             veldisp [km/s]; drawn holds the model's drawn labels plus veldisp [km/s]; noise_scale is the error
             multiplier (dimensionless, 1 to 5).
    """
    cube_by_name = dict(zip(SPECTRAL_DIMS, cube))
    unit_vector = training_unit_vector(cube_by_name, SPECTRAL_MODEL)
    bagpipes_components, drawn = MODELS[SPECTRAL_MODEL].components(unit_vector, fixed_redshift=redshift)
    veldisp = uniform(cube_by_name["veldisp"], *VELDISP_RANGE)
    bagpipes_components["veldisp"] = veldisp
    drawn["veldisp"] = veldisp
    # Uniform in log10 so that small and large error multipliers are equally likely per decade
    noise_scale = 10 ** uniform(cube_by_name["noise_scale"], *LOG_NOISE_SCALE_RANGE)
    return bagpipes_components, drawn, noise_scale
```

(c) **Section 6 constants.** After `PHOTOMETRY_JOB = ...`, add two constants. Leave the `SFH_MODEL` and `PHOTOMETRIC_DIMS` aliases from Task 1 in place; Task 11 removes them.

```python
JOB_DEFAULT_SFH_MODEL = "delayed"   # SFH model of a job or save_fit call naming none (every job before the DPL model)
EVIDENCE_ERROR_SAMPLES = 200        # posterior samples that recompute the ln Z error of a fit that stored none (count)
```

(d) **`fit_spectrum` docstring.** Replace its `Returns:` paragraph with:

```
    Returns: dict with samples (unit-cube points, shape (n_samples, len(SPECTRAL_DIMS))), log_evidence
             (dimensionless ln Z), log_evidence_error (MultiNest's uncertainty of ln Z, sqrt(H / n_live)),
             n_iterations (count, from multinest_iterations), converged (bool: the run certainly stopped before
             max_iterations, judged by the upper bound of the count) and run_time [s].
```

(e) **`run_multinest`.** Replace it whole:

```python
def run_multinest(log_likelihood, n_dims, output_basename, max_iterations, n_live_points):
    """One MultiNest run over the unit cube, with the settings shared by the spectral and photometric fits.

    The prior is the identity map: each fit's components function turns a unit-cube point into physical
    values, exactly as the training set does.
    Inputs: log_likelihood: function of one unit-cube point (dimensionless) giving a log-likelihood
            (dimensionless); it must never raise; n_dims: number of unit-cube dimensions (count);
            output_basename: directory prefix for the MultiNest files (str, under 100 characters);
            max_iterations: MultiNest iteration limit (count); n_live_points: live points (count).
    Returns: dict as fit_spectrum returns: samples (n_samples, n_dims), log_evidence, log_evidence_error,
             n_iterations, converged and run_time [s].
    """
    # pymultinest loads the MultiNest C library on import, so it needs LD_LIBRARY_PATH (set by the kernel)
    import pymultinest

    os.makedirs(output_basename, exist_ok=True)

    def unit_cube_prior(cube):
        """Identity map: the training prior is already uniform on the unit cube.

        Inputs: cube: unit-cube point (dimensionless). Returns: the same point (dimensionless).
        """
        return cube

    # 1. Run MultiNest, timing only the sampling itself
    start = time.time()
    result = pymultinest.solve(
        LogLikelihood=log_likelihood, Prior=unit_cube_prior, n_dims=n_dims,
        outputfiles_basename=output_basename, n_live_points=n_live_points,
        sampling_efficiency=SAMPLING_EFFICIENCY, multimodal=True, max_iter=max_iterations,
        # Standard nested-sampling evidence: the importance-sampled one underflows in short runs and
        # MultiNest then writes it in a form pymultinest cannot parse
        importance_nested_sampling=False,
        resume=False, verbose=False, use_MPI=False,
    )
    run_time = time.time() - start

    # 2. Hitting the iteration limit means the run stopped before it converged; the upper bound of the
    #    count is used so that a capped run is never reported as converged
    n_iterations, max_possible_iterations = multinest_iteration_bounds(output_basename)
    return {
        "samples": result["samples"],
        "log_evidence": result["logZ"],
        # MultiNest's "global log-evidence error", sqrt(H / n_live), read by pymultinest from stats.dat
        "log_evidence_error": result["logZerr"],
        "n_iterations": n_iterations,
        "converged": max_possible_iterations < max_iterations,
        "run_time": run_time,
    }
```

(f) **`posterior_labels`, `thinned_samples`, `sample_label_row`.** Replace all three:

```python
def posterior_labels(unit_cube_samples, redshift):
    """Labels of each spectral-fit posterior sample, from a Bagpipes model that also computes S-PLUS photometry.

    Keeps at most MAX_POSTERIOR_SAMPLES evenly spaced samples, because building labels is the slow part.
    model_r is the model r magnitude for 10^10 M_sun formed, used later to rescale mass-like labels.
    Inputs: unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, len(SPECTRAL_DIMS)));
            redshift: spectroscopic redshift (dimensionless).
    Returns: DataFrame with one row per kept sample: the delayed-tau training labels (mass-like ones for
             10^10 M_sun formed), then veldisp [km/s], noise_scale (dimensionless) and model_r [mag].
    """
    # 1. Thin the samples evenly: neighbouring MultiNest samples are nearly identical
    kept_samples = thinned_samples(unit_cube_samples)

    # 2. Labels per kept sample, with one photometry model updated in place
    model = None
    rows = []
    for cube in kept_samples:
        bagpipes_components, drawn, noise_scale = spectral_components(cube, redshift)
        model = updated_photometry_model(model, bagpipes_components)
        drawn_with_noise = dict(drawn)
        drawn_with_noise["noise_scale"] = noise_scale
        rows.append(sample_label_row(model, bagpipes_components, drawn_with_noise, SPECTRAL_MODEL))
    return pd.DataFrame(rows)


def thinned_samples(unit_cube_samples, max_samples=MAX_POSTERIOR_SAMPLES):
    """At most max_samples evenly spaced posterior samples, because computing anything per sample is the slow part.

    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, n_dims);
            max_samples: largest number of samples kept (count).
    Returns: every step-th row (dimensionless), step chosen so that no more than max_samples remain.
    """
    step = max(1, int(np.ceil(len(unit_cube_samples) / max_samples)))
    return unit_cube_samples[::step]


def sample_label_row(model, bagpipes_components, drawn, sfh_model):
    """Labels of one posterior sample, computed from its photometry model as for a training-set model.

    Inputs: model: photometry model_galaxy already at bagpipes_components; bagpipes_components: dict as from the
            model's components(), the point the model is at (not read: the labels come from drawn and the model);
            drawn: dict from components() with the model's drawn labels, plus any extra parameters of the fit
            (veldisp [km/s], noise_scale); sfh_model: "delayed" or "dblplaw".
    Returns: dict: the model's training labels in order (summary_labels), then the extra drawn parameters, then
             model_r [mag], the model r magnitude for 10^10 M_sun formed. Mass-like labels are for that mass.
    """
    row = summary_labels(model.sfh, drawn, sfh_model)
    for key, value in drawn.items():
        if key not in row:
            row[key] = value
    row["model_r"] = ab_magnitudes(model.photometry)[BANDS.index("r")]
    return row
```

(g) **`photometric_components`, `fit_photometry`, `photometric_posterior_labels`.** Replace all three:

```python
def photometric_components(cube, redshift, sfh_model):
    """Bagpipes components and drawn parameters for one photometric-fit unit-cube point of one SFH model.

    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as MODELS[sfh_model].photometric_dims;
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw".
    Returns: (bagpipes_components, drawn), exactly as MODELS[sfh_model].components(..., fixed_redshift=redshift)
             gives them.
    """
    cube_by_name = dict(zip(MODELS[sfh_model].photometric_dims, cube))
    unit_vector = training_unit_vector(cube_by_name, sfh_model)
    return MODELS[sfh_model].components(unit_vector, fixed_redshift=redshift)


def fit_photometry(magnitudes, errors, redshift, sfh_model, output_basename, max_iterations,
                   n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one galaxy's 12-band photometry under one SFH model's training prior at a fixed redshift.

    It fits the ML's own target: the model's training prior (its photometric_dims) and the likelihood of the
    simulator that made the training data (offset_marginal_log_likelihood), over the detected bands only.
    Inputs: magnitudes, errors: [mag], 12 values in BANDS order; a band is left out when detected_bands says
            it is not detected (99 or NaN magnitude, or error >= MAX_DETECTED_ERROR);
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw";
            output_basename, max_iterations and n_live_points: as in fit_spectrum.
    Returns: dict as fit_spectrum returns, with samples of shape (n_samples, len(MODELS[sfh_model].photometric_dims)).
    Raises: ValueError when fewer than MIN_DETECTED_BANDS bands are detected.
    """
    # 1. Bands that enter the likelihood
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)
    if detected.sum() < MIN_DETECTED_BANDS:
        raise ValueError(f"{detected.sum()} detected bands, fewer than MIN_DETECTED_BANDS = {MIN_DETECTED_BANDS}")

    # 2. Bagpipes model built once with the S-PLUS filters; every point has the same component set
    n_dims = len(MODELS[sfh_model].photometric_dims)
    midpoint_cube = np.full(n_dims, MIDPOINT_DRAW)
    starting_components, _ = photometric_components(midpoint_cube, redshift, sfh_model)
    model = updated_photometry_model(None, starting_components)

    # 3. Likelihood of one unit-cube point
    def log_likelihood(cube):
        """Log-likelihood of the photometry for one unit-cube point.

        Inputs: cube: unit-cube point (dimensionless, in the model's photometric_dims order).
        Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
        """
        bagpipes_components, _ = photometric_components(cube, redshift, sfh_model)
        # PyMultiNest calls sys.exit if the likelihood raises: a model Bagpipes cannot compute is impossible
        try:
            model.update(bagpipes_components)
            model_magnitudes = ab_magnitudes(model.photometry)
        except Exception:
            return LOG_ZERO
        return safe_photometric_log_likelihood(model_magnitudes, magnitudes, errors, detected)

    # 4. MultiNest over the unit cube
    return run_multinest(log_likelihood, n_dims, output_basename, max_iterations, n_live_points)


def photometric_posterior_labels(unit_cube_samples, redshift, magnitudes, errors, rng, sfh_model):
    """Labels of each photometric-fit posterior sample, at the brightness of the observed galaxy.

    Each kept sample gets its own magnitude offset, drawn from the offset's Gaussian given that sample's
    model (offset_marginal_log_likelihood), and its mass-like labels move by it (rescale_mass_labels), as the
    training labels move with the brightness offset of each model.
    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, len(MODELS[sfh_model].photometric_dims));
            redshift: spectroscopic redshift (dimensionless); magnitudes, errors: [mag] the fitted photometry,
            as given to fit_photometry; rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: DataFrame with one row per kept sample (at most MAX_POSTERIOR_SAMPLES): the model's training labels
             (MASS_DEPENDENT_LABELS already rescaled), model_r [mag, for 10^10 M_sun formed] and
             magnitude_offset [mag, observed minus model].
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)

    # 1. Labels per kept sample, each with an offset drawn given its own model
    kept_samples = thinned_samples(unit_cube_samples)
    model = None
    rows = []
    for cube in kept_samples:
        bagpipes_components, drawn = photometric_components(cube, redshift, sfh_model)
        model = updated_photometry_model(model, bagpipes_components)
        row = sample_label_row(model, bagpipes_components, drawn, sfh_model)
        model_magnitudes = ab_magnitudes(model.photometry)
        _, best_offset, offset_sigma = offset_marginal_log_likelihood(model_magnitudes, magnitudes, errors, detected)
        row["magnitude_offset"] = rng.normal(best_offset, offset_sigma)
        rows.append(row)

    # 2. Mass-like labels at the observed brightness
    labels = pd.DataFrame(rows)
    return rescale_mass_labels(labels, labels["magnitude_offset"].to_numpy())
```

(h) **Evidence helpers.** Add these right after `photometric_posterior_labels`:

```python
# -----------------------------------------------------------------------------
# Evidence uncertainty (DPL spec, part 4): stored by every new fit, recomputed for fits saved without it
# -----------------------------------------------------------------------------
def evidence_error(log_likelihoods, log_evidence, n_live_points):
    """MultiNest's uncertainty of ln Z, sqrt(H / n_live), with the information H taken from posterior samples.

    H = E_posterior[ln L] - ln Z is the information gained from prior to posterior [nats]; MultiNest reports
    sqrt(H / n_live) as the error of its ln Z.
    Inputs: log_likelihoods: ln L of equal-weight posterior samples (dimensionless), with the fit's own likelihood;
            log_evidence: the fit's ln Z (dimensionless); n_live_points: live points of the fit (count).
    Returns: the ln Z uncertainty (dimensionless float); 0 when the samples give H <= 0 (no information).
    """
    information = np.mean(np.asarray(log_likelihoods, dtype=float)) - log_evidence
    return float(np.sqrt(max(information, 0.0) / n_live_points))


def photometric_log_likelihoods(unit_cube_samples, redshift, magnitudes, errors, sfh_model,
                                max_samples=EVIDENCE_ERROR_SAMPLES):
    """The photometric fit's log-likelihood at its posterior samples, as fit_photometry computed it.

    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, len(MODELS[sfh_model].photometric_dims));
            redshift: spectroscopic redshift (dimensionless); magnitudes, errors: [mag] the fitted photometry;
            sfh_model: the fit's model; max_samples: at most this many evenly spaced samples are used (count).
    Returns: array of log-likelihoods (dimensionless), one per used sample.
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)
    model = None
    log_likelihoods = []
    for cube in thinned_samples(unit_cube_samples, max_samples):
        bagpipes_components, _ = photometric_components(cube, redshift, sfh_model)
        model = updated_photometry_model(model, bagpipes_components)
        model_magnitudes = ab_magnitudes(model.photometry)
        log_likelihoods.append(safe_photometric_log_likelihood(model_magnitudes, magnitudes, errors, detected))
    return np.array(log_likelihoods)


def load_unit_cube_samples(path):
    """Unit-cube posterior samples of a file written by save_fit.

    Inputs: path: file path (str). Returns: array (n_samples, n_dims), dimensionless.
    """
    with h5py.File(path, "r") as source:
        return source["unit_cube_samples"][()]


def fit_evidence(fit_file, job):
    """ln Z of a photometric fit and its uncertainty: the stored one, or sqrt(H / n_live) recomputed.

    Why recompute: fits saved before the DPL model stored no evidence error, and their MultiNest folders stayed
    in the PC job; the recomputed value is MultiNest's own estimate, from the stored posterior samples. A fit
    without n_live_points ran with N_LIVE_POINTS (submit_fits.sh never sets --n-live).
    Inputs: fit_file: path of a photometric fit (str); job: its photometry job (magnitudes and errors [mag],
            the photometry the fit used).
    Returns: (log_evidence, log_evidence_error), dimensionless floats.
    """
    _, metadata = load_fit(fit_file)
    if "log_evidence_error" in metadata:
        return float(metadata["log_evidence"]), float(metadata["log_evidence_error"])
    samples = load_unit_cube_samples(fit_file)
    log_likelihoods = photometric_log_likelihoods(samples, metadata["redshift"], job["magnitudes"], job["errors"],
                                                  metadata["sfh_model"])
    n_live_points = metadata.get("n_live_points", N_LIVE_POINTS)
    return float(metadata["log_evidence"]), evidence_error(log_likelihoods, metadata["log_evidence"], n_live_points)
```

(i) **`save_fit`, `load_fit`, `fit_is_current`, `job_model`.** Replace the first three and add `job_model` after `job_kind`:

```python
def save_fit(path, posterior, unit_cube_samples, metadata):
    """Write posterior labels, unit-cube samples and fit metadata to one HDF5 file, all or nothing.

    The file is written as path + ".tmp" and renamed to path only once complete, so an interrupted
    or failed write never leaves a partial file at path (which would count as a finished job). Every file
    records its sfh_model, so fits of an earlier or of the other SFH model can be told apart (fit_is_current).
    Inputs: path: output file path (str); posterior: DataFrame of labels, one column per label;
            unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, n_dims), n_dims that of the
            fit: len(SPECTRAL_DIMS) or the model's photometric_dims);
            metadata: dict of scalars (converged, run_time [s], n_iterations, max_iterations, n_live_points,
            log_evidence, log_evidence_error, redshift, kind, sfh_model and, for a cross-fit, truth_model; files
            written before kinds existed have no kind and are spectral). Without sfh_model, the fit is stored as
            JOB_DEFAULT_SFH_MODEL. The caller's dict is not changed.
    Returns: None. An existing file at path is replaced.
    """
    stored_metadata = dict(metadata)
    stored_metadata["sfh_model"] = metadata.get("sfh_model", JOB_DEFAULT_SFH_MODEL)
    temporary_path = path + ".tmp"
    with h5py.File(temporary_path, "w") as output:
        output.create_dataset("unit_cube_samples", data=unit_cube_samples)
        labels_group = output.create_group("posterior")
        for column in posterior.columns:
            labels_group.create_dataset(column, data=posterior[column].to_numpy(dtype=float))
        for key, value in stored_metadata.items():
            output.attrs[key] = value
    os.replace(temporary_path, path)


def load_fit(path):
    """Read a file written by save_fit.

    Inputs: path: file path (str).
    Returns: (posterior DataFrame with one column per label (in alphabetical order, as HDF5 stores them),
             metadata dict with the attributes save_fit wrote, sfh_model included; files written before the
             delayed-tau model have no sfh_model, and files written before the DPL model have no
             log_evidence_error and no n_live_points).
    """
    with h5py.File(path, "r") as source:
        columns = {}
        for column in source["posterior"]:
            columns[column] = source["posterior"][column][()]
        metadata = {}
        for key, value in source.attrs.items():
            # HDF5 returns numpy scalars; .item() gives back the plain Python value
            metadata[key] = value.item() if hasattr(value, "item") else value
    return pd.DataFrame(columns), metadata


def fit_is_current(metadata, sfh_model):
    """Whether a stored fit belongs to the current production run of a fits folder's SFH model.

    Fits that fail this (16-bin fits without sfh_model, fits of the other model, pilot runs with a smaller limit)
    are moved to superseded/ by supersede_outdated_fits (spectroscopic_reference.ipynb, section 4), never read
    as results. The folder's model comes from FIT_FOLDER_MODELS.
    Inputs: metadata: dict as load_fit returns it; sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: bool: True only when metadata["sfh_model"] == sfh_model and max_iterations >= MAX_ITERATIONS.
    """
    same_model = metadata.get("sfh_model") == sfh_model
    full_run = metadata.get("max_iterations", 0) >= MAX_ITERATIONS
    return same_model and full_run


def job_model(job):
    """SFH model a fit job is fitted with.

    Inputs: job: job dict (as fit_and_save takes).
    Returns: job["sfh_model"], or JOB_DEFAULT_SFH_MODEL ("delayed") for a job without one (every job exported
             before the DPL model).
    """
    return job.get("sfh_model", JOB_DEFAULT_SFH_MODEL)
```

(j) **`fit_spectrum_job`, `fit_photometry_job`, `fit_and_save`.** Replace all three:

```python
def fit_spectrum_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a spectrum job and label its posterior samples.

    Inputs: job: dict with spectrum (as from prepare_spectrum) and redshift (dimensionless); multinest_basename:
            directory prefix for the MultiNest files (str); max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_spectrum, posterior DataFrame from posterior_labels).
    Raises: ValueError for a job of a model other than SPECTRAL_MODEL (spectral fits are delayed-tau only).
    """
    if job_model(job) != SPECTRAL_MODEL:
        raise ValueError(f"spectral fits are {SPECTRAL_MODEL} only, not {job_model(job)} (DPL spec, decision 1)")
    fit = fit_spectrum(job["spectrum"], job["redshift"], multinest_basename, max_iterations, n_live_points)
    posterior = posterior_labels(fit["samples"], job["redshift"])
    return fit, posterior


def fit_photometry_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a photometry job with its SFH model and label its posterior samples at the galaxy's brightness.

    Inputs: job: dict with magnitudes and errors ([mag], 12 in BANDS order), redshift (dimensionless) and
            sfh_model (job_model); multinest_basename: directory prefix for the MultiNest files (str);
            max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_photometry, posterior DataFrame from photometric_posterior_labels).
    """
    sfh_model = job_model(job)
    fit = fit_photometry(job["magnitudes"], job["errors"], job["redshift"], sfh_model, multinest_basename,
                         max_iterations, n_live_points)
    # A fixed seed, so labelling the same samples again gives the same offsets
    rng = np.random.default_rng(OFFSET_DRAW_SEED)
    posterior = photometric_posterior_labels(fit["samples"], job["redshift"], job["magnitudes"], job["errors"],
                                             rng, sfh_model)
    return fit, posterior


# The fitter of each job kind: same signature, (job, multinest_basename, max_iterations, n_live_points)
JOB_FITTERS = {SPECTRUM_JOB: fit_spectrum_job, PHOTOMETRY_JOB: fit_photometry_job}


def fit_and_save(job, output_directory, max_iterations, n_live_points=N_LIVE_POINTS):
    """Fit one job and save <output_directory>/<job_id>.h5. Never raises.

    Skips jobs whose readable file of the same kind and SFH model already exists, so a batch can be re-run
    with more galaxies. An unreadable file does not count as finished: the job is fitted again and the file
    replaced. A readable file of the other kind is never replaced (job ids are specobjids, shared by the
    spectral and photometric fits of one galaxy), and neither is a fit of another SFH model (an earlier model's
    result, the only copy, to be moved to superseded/): the job fails instead. A failed fit comes back with
    status "failed" and the error text instead of an exception.
    Inputs: job: dict with job_id (str), redshift (dimensionless), kind (SPECTRUM_JOB when absent), sfh_model
            (JOB_DEFAULT_SFH_MODEL when absent) and, for a cross-fit, truth_model; a spectrum job holds spectrum
            (dict as from prepare_spectrum), a photometry job magnitudes and errors ([mag], 12 in BANDS order);
            output_directory: folder for the .h5 file and the MultiNest files (str); max_iterations: MultiNest
            iteration limit (count); n_live_points: live points (count).
    Returns: dict with job_id, status ("fitted", "skipped" or "failed"), converged (bool or None),
             run_time [s] (or None), n_iterations (count or None) and error (str or None).
    """
    # 1. Output path and a summary that the steps below fill in
    output_path = os.path.join(output_directory, f"{job['job_id']}.h5")
    kind = job_kind(job)
    sfh_model = job_model(job)
    summary = {"job_id": job["job_id"], "status": None, "converged": None,
               "run_time": None, "n_iterations": None, "error": None}

    # 2. Skip finished jobs; their run details come from the file. Any failure to read it (damaged,
    # truncated or missing fields) means the job is not finished, so it is fitted again below. A readable
    # fit of another SFH model is neither skipped (not a result of this model) nor overwritten (only copy)
    stored_kind = None
    stored_model = None
    if os.path.exists(output_path):
        try:
            _, metadata = load_fit(output_path)
            stored_details = {"converged": metadata["converged"], "run_time": metadata["run_time"],
                              "n_iterations": metadata["n_iterations"]}
            stored_kind = metadata.get("kind", SPECTRUM_JOB)
            stored_model = metadata.get("sfh_model")
        except Exception:
            pass
    if stored_kind == kind and stored_model != sfh_model:
        summary.update(status="failed", error=f"{output_path} holds a fit of another SFH model "
                                              f"({stored_model}, not {sfh_model}): move it to superseded/")
        return summary
    if stored_kind == kind:
        summary.update(status="skipped", **stored_details)
        return summary
    if stored_kind is not None:
        summary.update(status="failed", error=f"{output_path} holds a {stored_kind} fit and this job is a "
                                              f"{kind} fit: use a separate output folder for each kind")
        return summary

    # 3. Fit, label the posterior and save; any error is recorded rather than raised
    try:
        multinest_basename = os.path.join(output_directory, str(job["job_id"])) + "/"
        fit, posterior = JOB_FITTERS[kind](job, multinest_basename, max_iterations, n_live_points)
        metadata = {"converged": fit["converged"], "run_time": fit["run_time"],
                    "n_iterations": fit["n_iterations"], "max_iterations": max_iterations,
                    "n_live_points": n_live_points, "log_evidence": fit["log_evidence"],
                    "log_evidence_error": fit["log_evidence_error"], "redshift": job["redshift"], "kind": kind,
                    "sfh_model": sfh_model}
        if "truth_model" in job:
            metadata["truth_model"] = job["truth_model"]
        save_fit(output_path, posterior, fit["samples"], metadata)
        summary.update(status="fitted", converged=fit["converged"],
                       run_time=fit["run_time"], n_iterations=fit["n_iterations"])
    except Exception as error:
        summary.update(status="failed", error=repr(error))
    return summary
```

(k) **Section 7: `export_job` and `load_job`.** Replace both:

```python
def export_job(job, jobs_directory):
    """Write one fit job to <jobs_directory>/<job_id>.npz, a portable file for run_spectral_fits.py.

    Inputs: job: dict with job_id (str), redshift (dimensionless), kind (SPECTRUM_JOB when absent), sfh_model
            (JOB_DEFAULT_SFH_MODEL when absent) and, for a cross-fit, truth_model; a spectrum job holds spectrum
            (dict with wavelength [Angstrom, observed], flux and ivar), a photometry job magnitudes and errors
            ([mag], 12 in BANDS order); jobs_directory: folder for the file (str), created if missing.
    Returns: path of the written file (str).
    """
    os.makedirs(jobs_directory, exist_ok=True)
    path = os.path.join(jobs_directory, f"{job['job_id']}.npz")
    kind = job_kind(job)
    arrays = {"job_id": str(job["job_id"]), "redshift": job["redshift"], "kind": kind, "sfh_model": job_model(job)}
    if "truth_model" in job:
        arrays["truth_model"] = job["truth_model"]
    if kind == PHOTOMETRY_JOB:
        arrays["magnitudes"] = job["magnitudes"]
        arrays["errors"] = job["errors"]
    else:
        arrays["wavelength"] = job["spectrum"]["wavelength"]
        arrays["flux"] = job["spectrum"]["flux"]
        arrays["ivar"] = job["spectrum"]["ivar"]
    np.savez(path, **arrays)
    return path


def load_job(path):
    """Read a file written by export_job.

    Inputs: path: file path (str).
    Returns: job dict in the format fit_and_save takes: job_id (str), kind, sfh_model, redshift (dimensionless
             float), truth_model when the file has one, and either spectrum (wavelength, flux, ivar) or magnitudes
             and errors [mag]. A file without a kind (exported before the photometric fits existed) is a spectrum
             job; a file without sfh_model (exported before the DPL model) is a JOB_DEFAULT_SFH_MODEL job.
    """
    with np.load(path) as source:
        kind = SPECTRUM_JOB
        if "kind" in source.files:
            kind = str(source["kind"])
        sfh_model = JOB_DEFAULT_SFH_MODEL
        if "sfh_model" in source.files:
            sfh_model = str(source["sfh_model"])
        job = {"job_id": str(source["job_id"]), "kind": kind, "sfh_model": sfh_model,
               "redshift": float(source["redshift"])}
        if "truth_model" in source.files:
            job["truth_model"] = str(source["truth_model"])
        if kind == PHOTOMETRY_JOB:
            job["magnitudes"] = source["magnitudes"]
            job["errors"] = source["errors"]
        else:
            job["spectrum"] = {"wavelength": source["wavelength"], "flux": source["flux"], "ivar": source["ivar"]}
        return job
```

(l) **Section 8.** Add after `RUN_LOG_NAMES = ...`:

```python
# The SFH model of each fits folder: a fit of another model in it is superseded, never read
FIT_FOLDER_MODELS = {"spectral_fits": "delayed", "photometric_fits": "delayed", "photometric_fits_dblplaw": "dblplaw"}
```

Then replace `calibration_mocks_are_current` and `supersede_outdated_fits`:

```python
def calibration_mocks_are_current(mock_file, sfh_model):
    """Whether a stored calibration-mock file holds truths with the calibration labels of a folder's model.

    Inputs: mock_file: path of a joblib file with (jobs, truths), truths mapping job_id -> Series of the
            true labels; sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: bool: True when there is at least one truth and every truth has every
             MODELS[sfh_model].calibration_labels value (16-bin truths lack log_age, delayed-tau truths lack
             log_turnover_lookback). A file without mocks is not current, so it is redrawn, not reused.
    """
    _, truths = joblib.load(mock_file)
    if len(truths) == 0:
        return False
    for truth in truths.values():
        if not set(MODELS[sfh_model].calibration_labels) <= set(truth.index):
            return False
    return True


def supersede_outdated_fits(fits_directory, sfh_model):
    """Move every result of an earlier model or run out of a fits folder, into its SUPERSEDED_SUBDIRECTORY.

    Why: run_spectral_fits.py counts any top-level .h5 as done, and new calibration mocks reuse the calib_<k>
    ids, so an old file left in place would be skipped as finished or scored against a new truth. Each folder is
    judged against its own SFH model (FIT_FOLDER_MODELS): a DPL fit stays in photometric_fits_dblplaw and leaves
    photometric_fits, and the reverse.
    Moved, in this order:
      1. every readable .h5 fit that fails fit_is_current(metadata, sfh_model), with its MultiNest folder (the
         name without .h5); an unreadable .h5 stays, because fit_and_save fits that job again;
      2. CALIBRATION_MOCK_FILE when its truths lack the model's calibration labels (16-bin mocks) or it holds none;
      3. while no mock file is left, the jobs/calib_*.npz files, which then belong to no current mock;
      4. when any fit moved, the run logs RUN_LOG_NAMES (notebook 2 reads failures from them, and
         fetch_fits.sh drops new rows of ids already logged).
    Inputs: fits_directory: folder with the .h5 fits, jobs/ and the mock file (str); a missing folder is fine;
            sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: list of the moved names (str), relative to fits_directory; empty when everything is current.
    Raises: ValueError when the folder is one of FIT_FOLDER_MODELS and sfh_model is not its model.
    """
    folder_name = os.path.basename(os.path.normpath(fits_directory))
    if folder_name in FIT_FOLDER_MODELS and FIT_FOLDER_MODELS[folder_name] != sfh_model:
        raise ValueError(f"{fits_directory} holds {FIT_FOLDER_MODELS[folder_name]} fits, not {sfh_model}")
    if not os.path.isdir(fits_directory):
        return []
    moved = []

    # 1. Fits of an earlier model or run, or of the other model, with their MultiNest folders
    for name in sorted(os.listdir(fits_directory)):
        if not name.endswith(".h5"):
            continue
        try:
            _, metadata = load_fit(os.path.join(fits_directory, name))
        except (OSError, KeyError):
            continue
        if fit_is_current(metadata, sfh_model):
            continue
        for fit_name in (name, name[:-len(".h5")]):
            if move_to_superseded(fits_directory, fit_name) is not None:
                moved.append(fit_name)
    fits_were_moved = len(moved) > 0

    # 2. Calibration mocks of an earlier model
    mock_file = os.path.join(fits_directory, CALIBRATION_MOCK_FILE)
    if os.path.exists(mock_file) and not calibration_mocks_are_current(mock_file, sfh_model):
        move_to_superseded(fits_directory, CALIBRATION_MOCK_FILE)
        moved.append(CALIBRATION_MOCK_FILE)

    # 3. Job files of mocks that no stored mock file holds (the notebook stores mocks before exporting them)
    jobs_directory = os.path.join(fits_directory, JOBS_SUBDIRECTORY)
    if not os.path.exists(mock_file) and os.path.isdir(jobs_directory):
        for name in sorted(os.listdir(jobs_directory)):
            if name.startswith(CALIBRATION_ID_PREFIX):
                job_name = os.path.join(JOBS_SUBDIRECTORY, name)
                move_to_superseded(fits_directory, job_name)
                moved.append(job_name)

    # 4. Run logs of the moved fits
    if fits_were_moved:
        for name in RUN_LOG_NAMES:
            if move_to_superseded(fits_directory, name) is not None:
                moved.append(name)
    return moved
```

On the existing folders this behaves exactly as today. `spectral_fits` and `photometric_fits` are judged against `"delayed"`, which is what the removed `SFH_MODEL` check did, and their mock files hold delayed-τ truths, which have every delayed-τ calibration label.

- [ ] **Step 4: Docs of the fit worker**

In `run_spectral_fits.py`, edit the module docstring. Replace these lines:

```
the .h5 file. The two kinds share job ids (the specobjid), so they live in separate folders:
spectral_fits/ and photometric_fits/ (pass --jobs-dir photometric_fits/jobs --output-dir photometric_fits).
```

with:

```
the .h5 file, with the job's SFH model (sfh_model, "delayed" when the job file has none) and, for a cross-fit,
its true model (truth_model). The two kinds share job ids (the specobjid), so they live in separate folders:
spectral_fits/, photometric_fits/ (delayed-tau) and photometric_fits_dblplaw/ (double power law); pass e.g.
--jobs-dir photometric_fits_dblplaw/jobs --output-dir photometric_fits_dblplaw.
```

Also replace the docstring line 17, `   - synthetic_model.py, run_spectral_fits.py, slurm_fit.sh`, with `   - synthetic_model.py, sfh_models.py, run_spectral_fits.py, slurm_fit.sh` (Task 1 does not touch this docstring).

`slurm_fit.sh` and the bundle scripts need nothing here: Task 1 already lists `sfh_models.py` in `slurm_fit.sh` line 3 and copies it in `remote/submit_fits.sh` line 58 and `remote/submit_training.sh`.

- [ ] **Step 5: Notebook calls and the README supersede command**

Run this script from `Synthetic_domain/` with `../.venv/bin/python <script>`. Save it in your scratchpad, not in the repository.

```python
"""Task 4: section 4 supersedes every fits folder against its own model; fit_is_current and
photometric_components get the model argument."""
import json

NOTEBOOK = "spectroscopic_reference.ipynb"
with open(NOTEBOOK) as notebook_file:
    notebook = json.load(notebook_file)

SECTION_4_CELL = r'''# 1. Fits, calibration mocks, their job files and run logs of an earlier model or run (sm.supersede_outdated_fits),
# in every fits folder, each judged against its own SFH model (sm.FIT_FOLDER_MODELS). The photometric folders too:
# section 6b refuses calib_* fits whose mock file is missing
for fits_directory, sfh_model in sm.FIT_FOLDER_MODELS.items():
    superseded_directory = os.path.join(fits_directory, sm.SUPERSEDED_SUBDIRECTORY)
    moved = sm.supersede_outdated_fits(fits_directory, sfh_model)
    n_moved_fits = sum(name.endswith(".h5") for name in moved)
    if moved:
        print(f"{fits_directory}: moved {len(moved)} items ({n_moved_fits} fits) to {superseded_directory}/")
    else:
        print(f"{fits_directory}: nothing to supersede")

# 2. The first run's inject-and-recover mock file
if sm.move_to_superseded(FITS_DIRECTORY, os.path.basename(OLD_MOCK_FILE)) is not None:
    print(f"Moved {OLD_MOCK_FILE} to {os.path.join(FITS_DIRECTORY, sm.SUPERSEDED_SUBDIRECTORY)}/")'''

REPLACEMENTS = {
    "sm.fit_is_current(metadata)": 'sm.fit_is_current(metadata, "delayed")',
    'sm.photometric_components(photometric_cube, job["redshift"])':
        'sm.photometric_components(photometric_cube, job["redshift"], "delayed")',
}
EXPECTED_COUNTS = {"sm.fit_is_current(metadata)": 4,
                   'sm.photometric_components(photometric_cube, job["redshift"])': 1}

counts = {old: 0 for old in REPLACEMENTS}
section_4_cells = 0
for cell in notebook["cells"]:
    if cell["cell_type"] != "code":
        continue
    text = "".join(cell["source"])
    if text.startswith("# 1. Fits, calibration mocks, their job files and run logs"):
        cell["source"] = SECTION_4_CELL.splitlines(keepends=True)
        cell["outputs"] = []
        cell["execution_count"] = None
        section_4_cells += 1
        continue
    for old, new in REPLACEMENTS.items():
        counts[old] += text.count(old)
        text = text.replace(old, new)
    cell["source"] = text.splitlines(keepends=True)
assert section_4_cells == 1, section_4_cells
assert counts == EXPECTED_COUNTS, counts

with open(NOTEBOOK, "w") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print("notebook updated")
```

Verify, from `Synthetic_domain/`:

```bash
../.venv/bin/python -c "
import nbformat
nb = nbformat.read('spectroscopic_reference.ipynb', as_version=4); nbformat.validate(nb)
for i, c in enumerate(nb.cells):
    if c.cell_type == 'code': compile(c.source, f'cell {i}', 'exec')
print('ok', len(nb.cells))"
grep -c 'fit_is_current(metadata)' spectroscopic_reference.ipynb    # expect 0
```

Expected: `ok 31`, then `0`.

In `remote/README.md`, replace:

```
        [print(d, len(sm.supersede_outdated_fits(d))) for d in ('spectral_fits', 'photometric_fits')]"
```

with:

```
        [print(d, len(sm.supersede_outdated_fits(d, m))) for d, m in sm.FIT_FOLDER_MODELS.items()]"
```

- [ ] **Step 6: Run the tests and see them pass**

Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py -q -k "job_model or without_sfh_model or spectral_dims_stay or round_trip_the_models or own_model or dblplaw or fit_is_current or save_fit or save_and_load or photometric_components or photometric_posterior or fit_photometry or mocks_are_current or without_truths or supersede or evidence or log_likelihoods"`

Expected: PASS. `test_recomputed_evidence_error_agrees_with_multinest` may be SKIPPED if its short run does not converge. If it fails on the factor-2 bound, stop and report the two values; do not loosen the bound.

Then run the whole suite:

`cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q`

Expected: `211 passed` (about 5 minutes): Task 3's 194, plus the 17 new tests of this task. None of the existing tests is removed. If the short run of `test_recomputed_evidence_error_agrees_with_multinest` does not converge, the result is `210 passed, 1 skipped`.

- [ ] **Step 7: Commit**

```bash
git add Synthetic_domain/synthetic_model.py Synthetic_domain/run_spectral_fits.py \
        Synthetic_domain/spectroscopic_reference.ipynb Synthetic_domain/remote/README.md \
        Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: fit jobs and fit files carry their SFH model; each fits folder is judged by its own model

Photometric fits take the model from the job (sfh_model, delayed when absent), spectral fits stay
delayed-tau. New fits store MultiNest's ln Z error and their live points; fit_evidence recomputes
sqrt(H / n_live) from the stored samples of older fits. The existing delayed-tau tests
only gain the sfh_model argument and the new metadata keys; their assertions are unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Photometric calibration mocks from any SFH model's prior, and cross-fit jobs

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py`:
  - imports: add `from scipy.stats import qmc`;
  - section 7 (after `draw_calibration_mocks`, today ~1437–1455): new constants and functions;
  - section 8: `CROSS_FIT_MOCK_FILE`, the new `cross_fit_mocks_are_current`, and `supersede_outdated_fits` again (the Task 4 version, extended).
- Test: `Synthetic_domain/tests/test_synthetic_model.py` (sections "Job files and the calibration test" and "Superseding").

**Interfaces:**
- **Consumes (Task 4):** `photometric_components(cube, redshift, sfh_model)`, `sample_label_row(..., sfh_model)`, `updated_photometry_model`, `job_model`, `photometric_posterior_labels(..., sfh_model)`, `supersede_outdated_fits(fits_directory, sfh_model)`, `calibration_mocks_are_current(mock_file, sfh_model)`. **Consumes (sfh_models):** `MODELS`, `SHARED_LABELS`.
- **Produces (constants):**
  - `CROSS_CALIBRATION_ID_PREFIX = "cross_calib_"`;
  - `CROSS_FIT_MOCK_FILE = "cross_fit_mocks.joblib"`;
  - `PHOTOMETRIC_PREFIXES = {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}`;
  - `PHOTOMETRIC_EXTRA_LABELS = ["mw_age", "log_sfr"]`.
- **Produces (functions):**
  - `draw_photometric_cubes(n_mocks, sfh_model, seed) -> ndarray (n_mocks, len(photometric_dims))`: a seeded Latin hypercube;
  - `noiseless_photometric_mocks(cubes, redshifts, sfh_model) -> (magnitudes (n, 12), truths: list of Series)`;
  - `photometric_mock_jobs(job_ids, noiseless_magnitudes, redshifts, magnitude_offsets, truths, noise_table, rng, sfh_model) -> (jobs, rescaled_truths)`;
  - `draw_photometric_calibration_mocks(n_mocks, sfh_model, donor_redshifts, donor_r, noise_table, seed) -> (jobs, truths)`;
  - `cross_fit_id(mock_id) -> str`;
  - `cross_fit_jobs(mock_jobs, mock_truths, fitted_model) -> (jobs, truths)`;
  - `cross_fit_mocks_are_current(mock_file, sfh_model) -> bool`;
  - `photometric_summary_labels(sfh_model) -> list`;
  - `supersede_outdated_fits` also handles `CROSS_FIT_MOCK_FILE` and `jobs/cross_calib_*.npz`.

**Mock design (spec Part 2).** The delayed-τ photometric mocks keep their stored file. Their notebook code (Task 6) now calls the shared functions with the same rng order, `scale_to_observed_r` then `add_splus_noise` on `PHOTOMETRIC_MOCK_SEED`, so a redraw would give the same mocks.

The DPL mocks are built as follows:
- 33 points from a Latin hypercube over the 7 DPL photometric dimensions;
- each takes a random validation galaxy as donor: its spectroscopic redshift, and the shift to its aperture r;
- each is noised once with `add_splus_noise`;
- a mock with fewer than `MIN_DETECTED_BANDS` detected bands is left out, as for delayed-τ.

Cross-fit jobs copy a mock's photometry under the id `cross_calib_<k>`, with `sfh_model` set to the other model and `truth_model` to the mock's own. Each fits folder stores the cross-fits exported into it, with their truths, in `CROSS_FIT_MOCK_FILE`. That keeps the folder self-contained, as `CALIBRATION_MOCK_FILE` does for the mocks themselves.

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_synthetic_model.py`, at the end of the section "Job files and the calibration test":

```python
NOISE_MODEL_FILE = "splus_noise_model.csv"    # the stored S-PLUS error model (training_set.ipynb)


def test_photometric_cubes_are_a_seeded_latin_hypercube():
    cubes = sm.draw_photometric_cubes(33, "dblplaw", seed=2)

    assert cubes.shape == (33, len(sm.MODELS["dblplaw"].photometric_dims))
    assert np.array_equal(cubes, sm.draw_photometric_cubes(33, "dblplaw", seed=2))
    for column in cubes.T:    # one point in each of the 33 strata of every coordinate
        assert sorted(np.floor(column * 33).astype(int)) == list(range(33))


def test_noiseless_photometric_mocks_give_the_photometry_and_truth_of_each_point():
    cubes = np.random.default_rng(42).random((2, len(sm.MODELS["dblplaw"].photometric_dims)))
    redshifts = np.array([0.1, 0.2])

    magnitudes, truths = sm.noiseless_photometric_mocks(cubes, redshifts, "dblplaw")

    assert magnitudes.shape == (2, len(sm.BANDS))
    assert np.isfinite(magnitudes).all()
    for mock_index, truth in enumerate(truths):
        for label in (*sm.MODELS["dblplaw"].training_labels, "model_r"):
            assert label in truth.index, label
        assert truth["redshift"] == redshifts[mock_index]
        assert truth["model_r"] == magnitudes[mock_index, sm.BANDS.index("r")]
        _, drawn = sm.photometric_components(cubes[mock_index], redshifts[mock_index], "dblplaw")
        assert truth["alpha"] == drawn["alpha"]


def test_photometric_mock_jobs_move_the_truths_with_the_brightness_and_leave_out_undetected_mocks():
    noise_table = pd.read_csv(NOISE_MODEL_FILE)
    noiseless = np.full((2, len(sm.BANDS)), 18.0)
    truths = {}
    for job_id in ("m0", "m1"):
        truths[job_id] = pd.Series({"log_mass_formed": 10.0, "log_mstar": 9.8, "log_sfr": -1.0, "mw_age": 5.0})
    offsets = np.array([1.0, np.nan])     # m1 has no brightness donor: nothing is detected

    jobs, rescaled = sm.photometric_mock_jobs(["m0", "m1"], noiseless, np.array([0.1, 0.2]), offsets, truths,
                                              noise_table, np.random.default_rng(40), "dblplaw")
    repeat_jobs, _ = sm.photometric_mock_jobs(["m0", "m1"], noiseless, np.array([0.1, 0.2]), offsets, truths,
                                              noise_table, np.random.default_rng(40), "dblplaw")

    assert [job["job_id"] for job in jobs] == ["m0"]
    job = jobs[0]
    assert job["kind"] == sm.PHOTOMETRY_JOB and job["sfh_model"] == "dblplaw" and job["redshift"] == 0.1
    detected = sm.detected_bands(job["magnitudes"], job["errors"])
    assert detected.sum() >= sm.MIN_DETECTED_BANDS
    assert np.allclose(job["magnitudes"][detected], 19.0, atol=0.5)
    assert np.array_equal(job["magnitudes"], repeat_jobs[0]["magnitudes"], equal_nan=True)
    assert set(rescaled) == {"m0"}
    assert np.isclose(rescaled["m0"]["log_mstar"], 9.8 - sm.DEX_PER_MAGNITUDE)
    assert np.isclose(rescaled["m0"]["log_mass_formed"], 10.0 - sm.DEX_PER_MAGNITUDE)
    assert rescaled["m0"]["mw_age"] == 5.0


def test_draw_photometric_calibration_mocks_take_a_donor_redshift_and_brightness():
    noise_table = pd.read_csv(NOISE_MODEL_FILE)
    donor_redshifts = np.array([0.05, 0.1, 0.15])
    donor_r = np.array([16.0, 17.0, 99.0])    # the third galaxy is undetected in r: never a donor
    donor_r_by_redshift = {0.05: 16.0, 0.1: 17.0}

    jobs, truths = sm.draw_photometric_calibration_mocks(3, "dblplaw", donor_redshifts, donor_r, noise_table, 41)
    repeat_jobs, _ = sm.draw_photometric_calibration_mocks(3, "dblplaw", donor_redshifts, donor_r, noise_table, 41)

    assert len(jobs) >= 1
    assert {job["job_id"] for job in jobs} <= {"calib_000", "calib_001", "calib_002"}
    assert set(truths) == {job["job_id"] for job in jobs}
    for job, repeat in zip(jobs, repeat_jobs):
        assert job["redshift"] in donor_r_by_redshift
        assert job["sfh_model"] == "dblplaw"
        truth = truths[job["job_id"]]
        offset = donor_r_by_redshift[job["redshift"]] - truth["model_r"]
        assert np.isclose(truth["log_mass_formed"], sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * offset)
        for label in sm.MODELS["dblplaw"].calibration_labels:
            assert np.isfinite(truth[label]), label
        assert np.array_equal(job["magnitudes"], repeat["magnitudes"], equal_nan=True)


def test_calibration_coverage_runs_on_the_dblplaw_labels():
    """Truth = one posterior sample: every DPL calibration label has one usable mock and a finite coverage."""
    cubes = np.random.default_rng(43).random((3, len(sm.MODELS["dblplaw"].photometric_dims)))
    magnitudes, _ = sm.noiseless_photometric_mocks(cubes[:1], np.array([0.1]), "dblplaw")
    posterior = sm.photometric_posterior_labels(cubes, 0.1, magnitudes[0], np.full(len(sm.BANDS), 0.05),
                                                np.random.default_rng(44), "dblplaw")
    labels = sm.MODELS["dblplaw"].calibration_labels

    table = sm.calibration_coverage({"calib_000": posterior}, {"calib_000": posterior.iloc[1]}, labels)

    assert list(table["label"].unique()) == labels
    assert (table["n_jobs"] == 1).all()
    assert table["coverage"].notna().all()


def test_cross_fit_jobs_fit_one_model_s_mocks_with_the_other():
    mock_jobs = [{"job_id": "calib_000", "kind": "photometry", "magnitudes": np.full(12, 20.0),
                  "errors": np.full(12, 0.05), "redshift": 0.1}]        # stored before jobs named their model
    mock_truths = {"calib_000": pd.Series({"log_age": 0.5})}

    jobs, truths = sm.cross_fit_jobs(mock_jobs, mock_truths, "dblplaw")

    assert [job["job_id"] for job in jobs] == ["cross_calib_000"]
    assert jobs[0]["sfh_model"] == "dblplaw"
    assert jobs[0]["truth_model"] == "delayed"
    assert np.array_equal(jobs[0]["magnitudes"], mock_jobs[0]["magnitudes"])
    assert truths["cross_calib_000"] is mock_truths["calib_000"]
    assert "sfh_model" not in mock_jobs[0]     # the mock jobs are left alone
    with pytest.raises(ValueError, match="other model"):
        sm.cross_fit_jobs(mock_jobs, mock_truths, "delayed")
    assert sm.cross_fit_id("calib_007") == "cross_calib_007"
    with pytest.raises(ValueError):
        sm.cross_fit_id("10311134961714089984")


def test_photometric_summary_labels_and_prefixes():
    assert sm.photometric_summary_labels("delayed") == [
        "log_age", "log_tau", "log_metallicity", "log_mstar", "mw_age", "log_sfr",
        "av_ism", "log_ssfr", "f_young", "f_mid", "f_old"]
    assert sm.photometric_summary_labels("dblplaw") == [
        "log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "log_mstar", "mw_age", "log_sfr",
        "av_ism", "log_ssfr", "f_young", "f_mid", "f_old"]
    assert sm.PHOTOMETRIC_PREFIXES == {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}
```

At the end of the section "Superseding", add:

```python
def write_cross_mock_file(path, fitted_model, truth_model):
    """A cross-fit mock file as the notebook stores it: (jobs, truths) of one cross_calib_000 job."""
    job = {"job_id": "cross_calib_000", "kind": "photometry", "sfh_model": fitted_model, "truth_model": truth_model}
    labels = sm.MODELS[truth_model].training_labels
    joblib.dump(([job], {"cross_calib_000": pd.Series(np.zeros(len(labels)), index=labels)}), path)


def test_cross_fit_mocks_are_current_only_in_a_folder_of_their_fitted_model(tmp_path):
    path = tmp_path / sm.CROSS_FIT_MOCK_FILE
    write_cross_mock_file(path, "delayed", "dblplaw")
    assert sm.cross_fit_mocks_are_current(str(path), "delayed")
    assert not sm.cross_fit_mocks_are_current(str(path), "dblplaw")
    joblib.dump(([], {}), path)
    assert not sm.cross_fit_mocks_are_current(str(path), "delayed")


def test_supersede_moves_cross_jobs_only_without_a_current_cross_mock_file(tmp_path):
    """Current cross mocks, their jobs and fits stay; cross mocks made for the other folder go, with their jobs."""
    folder = tmp_path / "photometric_fits"
    jobs = folder / "jobs"
    jobs.mkdir(parents=True)
    (jobs / "cross_calib_000.npz").write_bytes(b"cross job")
    (jobs / "calib_000.npz").write_bytes(b"mock job")
    write_mock_file(folder / sm.CALIBRATION_MOCK_FILE, sm.MODELS["delayed"].training_labels)
    write_cross_mock_file(folder / sm.CROSS_FIT_MOCK_FILE, "delayed", "dblplaw")
    write_fit(folder / "cross_calib_000.h5")       # a delayed-tau fit of a DPL mock: current in this folder

    assert sm.supersede_outdated_fits(str(folder), "delayed") == []

    write_cross_mock_file(folder / sm.CROSS_FIT_MOCK_FILE, "dblplaw", "delayed")   # made for the other folder
    moved = sm.supersede_outdated_fits(str(folder), "delayed")

    assert sorted(moved) == sorted([sm.CROSS_FIT_MOCK_FILE, os.path.join("jobs", "cross_calib_000.npz")])
    assert (jobs / "calib_000.npz").exists()
    assert (folder / "cross_calib_000.h5").exists()
```

- [ ] **Step 2: Run the tests and see them fail**

Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py -q -k "photometric_cubes or noiseless_photometric or photometric_mock_jobs or draw_photometric or dblplaw_labels or cross or summary_labels"`

Expected: FAIL, with `AttributeError: module 'synthetic_model' has no attribute 'draw_photometric_cubes'` (and likewise `noiseless_photometric_mocks`, `photometric_mock_jobs`, `draw_photometric_calibration_mocks`, `cross_fit_jobs`, `cross_fit_id`, `photometric_summary_labels`, `CROSS_FIT_MOCK_FILE`, `cross_fit_mocks_are_current`).

- [ ] **Step 3: Implement**

(a) Imports: add `from scipy.stats import qmc` after `from scipy.stats import chi2 as chi_squared`.

(b) Section 7 constants. After `CALIBRATION_ID_PREFIX = ...`, add:

```python
CROSS_CALIBRATION_ID_PREFIX = "cross_calib_"   # cross-fit ids: cross_calib_007 is mock calib_007 fitted with the other model
# Reference-CSV columns of each model's photometric fits: <prefix><label>_p16, <prefix>log_evidence, ...
PHOTOMETRIC_PREFIXES = {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}
PHOTOMETRIC_EXTRA_LABELS = ["mw_age", "log_sfr"]   # summarised besides the calibration and the shared labels
```

(c) Section 7 functions. Add these after `draw_calibration_mocks`:

```python
# -----------------------------------------------------------------------------
# Photometric calibration mocks of any SFH model, and the cross-fits of the evidence test (DPL spec, part 2)
# -----------------------------------------------------------------------------
def draw_photometric_cubes(n_mocks, sfh_model, seed):
    """Unit-cube points of photometric calibration mocks: a Latin hypercube over the model's photometric dims.

    Why a Latin hypercube: as for the training sets, so that 33 mocks still cover every parameter range evenly.
    Inputs: n_mocks: number of points (count); sfh_model: "delayed" or "dblplaw"; seed: Latin hypercube seed (integer).
    Returns: array (n_mocks, len(MODELS[sfh_model].photometric_dims)), values in [0, 1) (dimensionless).
    """
    sampler = qmc.LatinHypercube(len(MODELS[sfh_model].photometric_dims), seed=seed)
    return sampler.random(n_mocks)


def noiseless_photometric_mocks(cubes, redshifts, sfh_model):
    """Noiseless S-PLUS photometry and true labels of known photometric-fit points, for 10^10 M_sun formed.

    Inputs: cubes: unit-cube points (dimensionless), shape (n_mocks, len(MODELS[sfh_model].photometric_dims));
            redshifts: one redshift per point (dimensionless); sfh_model: "delayed" or "dblplaw".
    Returns: (magnitudes [mag], shape (n_mocks, 12) in BANDS order, NaN where a flux is not positive;
             truths: list of Series, one per point, as sample_label_row gives them: the model's training labels,
             mass-like ones for 10^10 M_sun formed, and model_r [mag]).
    """
    magnitudes = np.zeros((len(cubes), len(BANDS)))
    truths = []
    model = None
    for mock_index, cube in enumerate(cubes):
        bagpipes_components, drawn = photometric_components(cube, redshifts[mock_index], sfh_model)
        model = updated_photometry_model(model, bagpipes_components)
        magnitudes[mock_index] = ab_magnitudes(model.photometry)
        truths.append(pd.Series(sample_label_row(model, bagpipes_components, drawn, sfh_model)))
    return magnitudes, truths


def photometric_mock_jobs(job_ids, noiseless_magnitudes, redshifts, magnitude_offsets, truths, noise_table, rng,
                          sfh_model):
    """Photometry jobs of mocks at given brightness, noised once as the training data, with their true labels.

    Each mock moves by its magnitude offset (observed minus model, the same in every band), then gets S-PLUS
    errors, flux noise and non-detections (add_splus_noise); its true mass-like labels move with the offset
    (rescale_mass_labels). A mock left with fewer than MIN_DETECTED_BANDS detected bands cannot be fitted and is
    left out.
    Inputs: job_ids: one id per mock (str); noiseless_magnitudes: [mag] (n_mocks, 12) in BANDS order, for
            10^10 M_sun formed; redshifts: one per mock (dimensionless); magnitude_offsets: [mag] one per mock
            (NaN gives a mock with no detection, so it is left out); truths: dict job_id -> Series of the true
            labels for 10^10 M_sun formed; noise_table: output of fit_noise_model; rng: numpy Generator (only
            add_splus_noise draws from it); sfh_model: the model the jobs are fitted with.
    Returns: (jobs, rescaled_truths): jobs is a list of photometry job dicts as fit_and_save takes (job_id, kind,
             sfh_model, magnitudes, errors [mag], redshift), in input order without the left-out mocks;
             rescaled_truths maps each kept job_id to its truth at the mock's brightness.
    """
    offsets = np.asarray(magnitude_offsets, dtype=float)
    scaled_magnitudes = noiseless_magnitudes + offsets[:, np.newaxis]
    noisy_magnitudes, noisy_errors, _ = add_splus_noise(scaled_magnitudes, noise_table, rng)

    jobs = []
    rescaled_truths = {}
    for mock_index, job_id in enumerate(job_ids):
        detected = detected_bands(noisy_magnitudes[mock_index], noisy_errors[mock_index])
        if detected.sum() < MIN_DETECTED_BANDS:
            continue
        jobs.append({"job_id": job_id, "kind": PHOTOMETRY_JOB, "sfh_model": sfh_model,
                     "magnitudes": noisy_magnitudes[mock_index], "errors": noisy_errors[mock_index],
                     "redshift": float(redshifts[mock_index])})
        truth = truths[job_id].to_frame().T
        rescaled_truths[job_id] = rescale_mass_labels(truth, offsets[mock_index]).iloc[0]
    return jobs, rescaled_truths


def draw_photometric_calibration_mocks(n_mocks, sfh_model, donor_redshifts, donor_r, noise_table, seed):
    """Photometric calibration mocks from a model's prior, each at the redshift and brightness of a real galaxy.

    Why a donor galaxy: the mock then has a real redshift and aperture brightness, and so realistic S-PLUS errors.
    The unit-cube points come from draw_photometric_cubes(seed); the donors and the noise from
    numpy.random.default_rng(seed), donors first.
    Inputs: n_mocks: number of mocks (count); sfh_model: "delayed" or "dblplaw"; donor_redshifts: spectroscopic
            redshifts of the candidate donors (dimensionless); donor_r: their aperture r magnitudes [mag] (a donor
            with a non-finite r, or r at or above UNDETECTED_MAGNITUDE, is never chosen); noise_table: output of
            fit_noise_model; seed: integer.
    Returns: (jobs, truths) as photometric_mock_jobs returns them, with ids calib_000, calib_001, ... (a mock with
             too few detected bands is left out, so an id can be missing); truths at the donor's brightness.
    Raises: ValueError when no donor is usable.
    """
    donor_redshifts = np.asarray(donor_redshifts, dtype=float)
    donor_r = np.asarray(donor_r, dtype=float)
    usable = np.isfinite(donor_r) & (donor_r < UNDETECTED_MAGNITUDE) & np.isfinite(donor_redshifts)
    if not usable.any():
        raise ValueError("no usable brightness donor (finite redshift and detected r)")
    usable_redshifts = donor_redshifts[usable]
    usable_r = donor_r[usable]

    # 1. Points of the prior, and one donor galaxy per mock
    cubes = draw_photometric_cubes(n_mocks, sfh_model, seed)
    rng = np.random.default_rng(seed)
    donors = rng.integers(len(usable_r), size=n_mocks)
    redshifts = usable_redshifts[donors]

    # 2. Noiseless photometry for 10^10 M_sun formed, moved to the donor's r magnitude, then noised once
    noiseless_magnitudes, truth_list = noiseless_photometric_mocks(cubes, redshifts, sfh_model)
    magnitude_offsets = usable_r[donors] - noiseless_magnitudes[:, BANDS.index("r")]
    job_ids = []
    truths = {}
    for mock_index in range(n_mocks):
        job_id = f"{CALIBRATION_ID_PREFIX}{mock_index:03d}"
        job_ids.append(job_id)
        truths[job_id] = truth_list[mock_index]
    return photometric_mock_jobs(job_ids, noiseless_magnitudes, redshifts, magnitude_offsets, truths, noise_table,
                                 rng, sfh_model)


def cross_fit_id(mock_id):
    """Job id of the cross-fit of a calibration mock: calib_007 -> cross_calib_007.

    Inputs: mock_id: calibration-mock id (str, CALIBRATION_ID_PREFIX followed by the three-digit index).
    Returns: the cross-fit id (str).
    Raises: ValueError for an id that is not a calibration-mock id.
    """
    if not mock_id.startswith(CALIBRATION_ID_PREFIX):
        raise ValueError(f"{mock_id} is not a calibration-mock id ({CALIBRATION_ID_PREFIX}<k>)")
    return CROSS_CALIBRATION_ID_PREFIX + mock_id[len(CALIBRATION_ID_PREFIX):]


def cross_fit_jobs(mock_jobs, mock_truths, fitted_model):
    """One model's photometric mocks as jobs fitted with the other model, for the evidence test.

    Why: with the truth known, the evidence of the true model should win; if it does not, the colours cannot
    tell the two models apart. The photometry, redshift and truth of each mock are kept unchanged.
    Inputs: mock_jobs: photometry jobs of one model's calibration mocks (as in its mock file; a job without
            sfh_model is delayed-tau); mock_truths: dict mock job_id -> Series of the true labels;
            fitted_model: the model to fit them with, not the mocks' own.
    Returns: (jobs, truths): one job per mock with job_id cross_fit_id(mock id), sfh_model = fitted_model and
             truth_model = the mock's own model; truths maps the new ids to the mocks' truths. The input jobs are
             not changed.
    Raises: ValueError when a mock is already of fitted_model.
    """
    jobs = []
    truths = {}
    for mock_job in mock_jobs:
        truth_model = job_model(mock_job)
        if truth_model == fitted_model:
            raise ValueError(f"{mock_job['job_id']} is a {truth_model} mock: a cross-fit uses the other model")
        job = dict(mock_job)
        job["job_id"] = cross_fit_id(mock_job["job_id"])
        job["sfh_model"] = fitted_model
        job["truth_model"] = truth_model
        jobs.append(job)
        truths[job["job_id"]] = mock_truths[mock_job["job_id"]]
    return jobs, truths


def photometric_summary_labels(sfh_model):
    """Labels of a model's photometric fits that the reference CSV summarises.

    Inputs: sfh_model: "delayed" or "dblplaw".
    Returns: list of label names, each once: the model's calibration labels, then PHOTOMETRIC_EXTRA_LABELS, then
             the SHARED_LABELS not listed yet (the labels the model comparison uses).
    """
    labels = list(MODELS[sfh_model].calibration_labels)
    for label in [*PHOTOMETRIC_EXTRA_LABELS, *SHARED_LABELS]:
        if label not in labels:
            labels.append(label)
    return labels
```

(d) Section 8. After `CALIBRATION_MOCK_FILE = ...`, add:

```python
CROSS_FIT_MOCK_FILE = "cross_fit_mocks.joblib"   # (jobs, truths) of the cross-fits exported into a fits folder
```

After `calibration_mocks_are_current`, add:

```python
def cross_fit_mocks_are_current(mock_file, sfh_model):
    """Whether a stored cross-fit mock file belongs to a fits folder of sfh_model and has current truths.

    Inputs: mock_file: path of a joblib file with (jobs, truths) as cross_fit_jobs returns them;
            sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: bool: True when there is at least one job, every job is fitted with sfh_model, and every truth has
             the calibration labels of its job's truth_model.
    """
    jobs, truths = joblib.load(mock_file)
    if len(jobs) == 0:
        return False
    for job in jobs:
        truth_model = job.get("truth_model")
        if job_model(job) != sfh_model or truth_model not in MODELS or job["job_id"] not in truths:
            return False
        if not set(MODELS[truth_model].calibration_labels) <= set(truths[job["job_id"]].index):
            return False
    return True
```

Then replace `supersede_outdated_fits` with this final version. Its docstring and steps 2–3 add the cross-fit mock file; the rest is as in Task 4.

```python
def supersede_outdated_fits(fits_directory, sfh_model):
    """Move every result of an earlier model or run out of a fits folder, into its SUPERSEDED_SUBDIRECTORY.

    Why: run_spectral_fits.py counts any top-level .h5 as done, and new calibration mocks reuse the calib_<k>
    ids, so an old file left in place would be skipped as finished or scored against a new truth. Each folder is
    judged against its own SFH model (FIT_FOLDER_MODELS): a DPL fit stays in photometric_fits_dblplaw and leaves
    photometric_fits, and the reverse.
    Moved, in this order:
      1. every readable .h5 fit that fails fit_is_current(metadata, sfh_model), with its MultiNest folder (the
         name without .h5); an unreadable .h5 stays, because fit_and_save fits that job again;
      2. CALIBRATION_MOCK_FILE when its truths lack the model's calibration labels (16-bin mocks) or it holds none,
         and CROSS_FIT_MOCK_FILE when cross_fit_mocks_are_current says it does not belong here;
      3. while no mock file is left, the jobs/calib_*.npz files, which then belong to no current mock, and while
         no cross-fit mock file is left, the jobs/cross_calib_*.npz files;
      4. when any fit moved, the run logs RUN_LOG_NAMES (notebook 2 reads failures from them, and
         fetch_fits.sh drops new rows of ids already logged).
    Inputs: fits_directory: folder with the .h5 fits, jobs/ and the mock files (str); a missing folder is fine;
            sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: list of the moved names (str), relative to fits_directory; empty when everything is current.
    Raises: ValueError when the folder is one of FIT_FOLDER_MODELS and sfh_model is not its model.
    """
    folder_name = os.path.basename(os.path.normpath(fits_directory))
    if folder_name in FIT_FOLDER_MODELS and FIT_FOLDER_MODELS[folder_name] != sfh_model:
        raise ValueError(f"{fits_directory} holds {FIT_FOLDER_MODELS[folder_name]} fits, not {sfh_model}")
    if not os.path.isdir(fits_directory):
        return []
    moved = []

    # 1. Fits of an earlier model or run, or of the other model, with their MultiNest folders
    for name in sorted(os.listdir(fits_directory)):
        if not name.endswith(".h5"):
            continue
        try:
            _, metadata = load_fit(os.path.join(fits_directory, name))
        except (OSError, KeyError):
            continue
        if fit_is_current(metadata, sfh_model):
            continue
        for fit_name in (name, name[:-len(".h5")]):
            if move_to_superseded(fits_directory, fit_name) is not None:
                moved.append(fit_name)
    fits_were_moved = len(moved) > 0

    # 2. Mock files of an earlier model, or made for the other folder
    mock_file = os.path.join(fits_directory, CALIBRATION_MOCK_FILE)
    if os.path.exists(mock_file) and not calibration_mocks_are_current(mock_file, sfh_model):
        move_to_superseded(fits_directory, CALIBRATION_MOCK_FILE)
        moved.append(CALIBRATION_MOCK_FILE)
    cross_mock_file = os.path.join(fits_directory, CROSS_FIT_MOCK_FILE)
    if os.path.exists(cross_mock_file) and not cross_fit_mocks_are_current(cross_mock_file, sfh_model):
        move_to_superseded(fits_directory, CROSS_FIT_MOCK_FILE)
        moved.append(CROSS_FIT_MOCK_FILE)

    # 3. Job files of mocks that no stored mock file holds (the notebook stores mocks before exporting them)
    orphan_prefixes = []
    if not os.path.exists(mock_file):
        orphan_prefixes.append(CALIBRATION_ID_PREFIX)
    if not os.path.exists(cross_mock_file):
        orphan_prefixes.append(CROSS_CALIBRATION_ID_PREFIX)
    jobs_directory = os.path.join(fits_directory, JOBS_SUBDIRECTORY)
    if os.path.isdir(jobs_directory):
        for name in sorted(os.listdir(jobs_directory)):
            for prefix in orphan_prefixes:
                if name.startswith(prefix):
                    job_name = os.path.join(JOBS_SUBDIRECTORY, name)
                    move_to_superseded(fits_directory, job_name)
                    moved.append(job_name)

    # 4. Run logs of the moved fits
    if fits_were_moved:
        for name in RUN_LOG_NAMES:
            if move_to_superseded(fits_directory, name) is not None:
                moved.append(name)
    return moved
```

`"cross_calib_..."` does not start with `"calib_"`, so the two job prefixes never catch each other's files. The `glob("calib_*.h5")` guards of notebook 2 also miss the cross-fits.

- [ ] **Step 4: Run the tests and see them pass**

Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py -q -k "photometric_cubes or noiseless_photometric or photometric_mock_jobs or draw_photometric or dblplaw_labels or cross or summary_labels or supersede"`

Expected: PASS.

Then run the whole suite: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q`

Expected: `220 passed`: Task 4's 211, plus the 9 new tests of this task. As in Task 4, the result is `219 passed, 1 skipped` when the short MultiNest run of `test_recomputed_evidence_error_agrees_with_multinest` does not converge.

- [ ] **Step 5: Commit**

```bash
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: photometric calibration mocks from any SFH model's prior, and cross-fit jobs

DPL mocks: a Latin hypercube over the DPL prior, each at the redshift and aperture r of a
validation galaxy, noised once. Cross-fits (cross_calib_<k>) fit each model's mocks with the
other model; each fits folder stores its own cross-fit mocks and supersedes them by its model.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Notebook 2 for the DPL photometric fits and the cross-fits; remote scripts for `photometric_fits_dblplaw`

**Files:**
- Modify: `Synthetic_domain/spectroscopic_reference.ipynb`:
  - cell 0 (markdown, Outputs list);
  - cell 1 (parameter cell: the new `EXPORT_LOCKED_TEST_JOBS` switch after `SEED = 0`);
  - cell 9 (markdown `## 4. Fits (outside the notebook)`, the "Superseding" paragraph);
  - cells 17–20 (section 6b: markdown, mocks/jobs, calibration, summaries), replaced;
  - five new cells after cell 20 (section 6c).
- Modify: `Synthetic_domain/remote/submit_fits.sh` (header comment lines 8–10 only; the `cp` on line 58 already ships `sfh_models.py` since Task 1 and is not touched again), `Synthetic_domain/remote/fit_job.sh` (line 13 comment), `Synthetic_domain/remote/fetch_fits.sh` (line 6 comment), `Synthetic_domain/remote/README.md` (photometric-fits paragraph lines 28–36, end of file), `.gitignore` (repo root).
- Test: notebook validation and compile check; `bash -n`; a `BUNDLE_ONLY` dry run of `submit_fits.sh`.

**Interfaces:**
- **Consumes:**
  - Task 4: `fit_is_current(metadata, sfh_model)`, `FIT_FOLDER_MODELS`, `supersede_outdated_fits(dir, model)`, `fit_evidence`, `PHOTOMETRY_JOB`, `export_job`;
  - Task 5: `noiseless_photometric_mocks`, `photometric_mock_jobs`, `draw_photometric_calibration_mocks`, `cross_fit_jobs`, `cross_fit_id`, `CROSS_FIT_MOCK_FILE`, `CROSS_CALIBRATION_ID_PREFIX`, `PHOTOMETRIC_PREFIXES`, `photometric_summary_labels`;
  - `sm.model_file`, `sm.MODELS`.
- **Produces (files, written by notebook 2 on the PC):**
  - `photometric_fits_dblplaw/jobs/*.npz` (DPL mocks `calib_<k>`, validation galaxies (plus the test half with `EXPORT_LOCKED_TEST_JOBS`), delayed-τ mocks' cross-fits `cross_calib_<k>`);
  - `photometric_fits_dblplaw/calibration_mocks.joblib`;
  - `photometric_fits_dblplaw/cross_fit_mocks.joblib`;
  - `photometric_fits/jobs/cross_calib_<k>.npz` (DPL mocks fitted with delayed-τ) and `photometric_fits/cross_fit_mocks.joblib`;
  - `photometric_cross_fits.csv`.
- **Produces (`spectroscopic_reference.csv` columns):**
  - `phot_<label>_p16/_median/_p84` for `sm.photometric_summary_labels("delayed")` (now also `av_ism`, `log_ssfr`, `f_young`, `f_mid`, `f_old`), plus `phot_converged`, `phot_log_evidence`, `phot_log_evidence_error`, `phot_fitter_calibrated`;
  - `phot_dblplaw_<label>_p16/_median/_p84` for `sm.photometric_summary_labels("dblplaw")`, plus `phot_dblplaw_converged`, `phot_dblplaw_log_evidence`, `phot_dblplaw_log_evidence_error`, `phot_dblplaw_fitter_calibrated`;
  - these columns have values for the fitted validation galaxies. With `EXPORT_LOCKED_TEST_JOBS = True` (once, for the locked test), the fitted test-half galaxies get values too. Otherwise the test half's rows are NaN, although the columns exist.
- **Produces (`photometric_cross_fits.csv`, one row per mock of either model, 66 when complete):** `mock_id` (`calib_<k>`), `truth_model`, `phot_converged`, `phot_log_evidence`, `phot_log_evidence_error`, `phot_dblplaw_converged`, `phot_dblplaw_log_evidence`, `phot_dblplaw_log_evidence_error`. The `phot_*` columns come from the delayed-τ fit of that mock and the `phot_dblplaw_*` columns from the DPL fit, whether that fit is the mock's own (`calib_<k>`) or a cross-fit (`cross_calib_<k>`). The column names are those of the reference CSV, not those of Task 9: `ev.true_model_recovery` reads `log_evidence_delayed` and `log_evidence_dblplaw`. Task 10 reads this file, keeps the mocks whose two fits both converged, and renames `phot_log_evidence` to `log_evidence_delayed` and `phot_dblplaw_log_evidence` to `log_evidence_dblplaw` before it calls `ev.true_model_recovery`. It also reports `ev.log_bayes_factor_error` of the two `_error` columns. A missing fit gives NaN, and `_converged` is then empty.

**Which jobs go to which folder:**
- `photometric_fits/jobs/` keeps the delayed-τ jobs. It also gets the 33 `cross_calib_<k>` jobs of the DPL mocks; these are `sfh_model = "delayed"`, `truth_model = "dblplaw"`, and so are current delayed-τ fits there.
- `photometric_fits_dblplaw/jobs/` gets the DPL mocks, the validation galaxies, and the 33 `cross_calib_<k>` jobs of the delayed-τ mocks (`sfh_model = "dblplaw"`, `truth_model = "delayed"`).
- After the DPL folder's PC run, a `FITS_SUBDIR=photometric_fits` run fits only the 33 cross-fits, because the finished delayed-τ fits travel in its bundle.
- **Locked test (spec amendment 2).** `EXPORT_LOCKED_TEST_JOBS = False` in the parameter cell. Set to True once, together with `RUN_LOCKED_TEST` of `ml_model.ipynb`, it does two things:
  - it exports the test-half galaxy jobs of both models into the same two folders, so both fits runs then include them;
  - the reading cells fill the `phot_*` and `phot_dblplaw_*` columns of the test galaxies.

  The brightness donors of the mocks stay as they are: the delayed-τ mocks exclude the test half and the DPL mocks use the validation half. So the stored mocks do not change.

- [ ] **Step 1: Edit the notebook**

Run this script from `Synthetic_domain/` with `../.venv/bin/python <script>`. Save it in your scratchpad, not in the repository.

```python
"""Task 6: the locked-test switch, section 6b (delayed-tau) rebuilt on shared helpers, section 6c (DPL and
cross-fits) added."""
import json
import uuid

NOTEBOOK = "spectroscopic_reference.ipynb"
with open(NOTEBOOK) as notebook_file:
    notebook = json.load(notebook_file)
cells = notebook["cells"]


def source_of(cell):
    return "".join(cell["source"])


def find_cell(first_line):
    """Index of the one cell whose source starts with first_line."""
    matches = [index for index, cell in enumerate(cells) if source_of(cell).startswith(first_line)]
    assert len(matches) == 1, (first_line, matches)
    return matches[0]


def set_source(cell, text):
    cell["source"] = text.splitlines(keepends=True)
    if cell["cell_type"] == "code":
        cell["outputs"] = []
        cell["execution_count"] = None


def new_cell(cell_type, text):
    cell = {"cell_type": cell_type, "id": uuid.uuid4().hex[:8], "metadata": {}, "source": []}
    if cell_type == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    set_source(cell, text)
    return cell


# ---------------------------------------------------------------------------------------------------------------
PARAMETERS_OLD = "\nSEED = 0\n"
PARAMETERS_NEW = ("\nSEED = 0\n"
                  "# Locked test (DPL spec, amendment 2). False: the photometric fits cover the validation half only.\n"
                  "# Set it True once, together with RUN_LOCKED_TEST of ml_model.ipynb, to also export and read the test\n"
                  "# half's photometric jobs of both SFH models; keep it True afterwards, or the next run empties the\n"
                  "# test half's phot_* and phot_dblplaw_* columns again\n"
                  "EXPORT_LOCKED_TEST_JOBS = False\n")

OUTPUTS_OLD = ("`photometric_fits/jobs/` and `photometric_fits/calibration_mocks.joblib` (photometric fits, "
               "section 6b), ")
OUTPUTS_NEW = ("`photometric_fits/jobs/` and `photometric_fits/calibration_mocks.joblib` (photometric fits, "
               "section 6b), `photometric_fits_dblplaw/jobs/` and `photometric_fits_dblplaw/calibration_mocks.joblib` "
               "(double-power-law photometric fits, section 6c), `cross_fit_mocks.joblib` in both photometric "
               "folders and `photometric_cross_fits.csv` (cross-fits, section 6c), ")

SUPERSEDING = r'''**Superseding:** results of an earlier model or run are first moved to `superseded/` inside every fits folder of `sm.FIT_FOLDER_MODELS` (`spectral_fits/`, `photometric_fits/`, `photometric_fits_dblplaw/`), each judged against its own SFH model (`sm.supersede_outdated_fits`): every fit that fails `sm.fit_is_current` (another `sfh_model` than the folder's, as for the 16-bin fits that have none, or a lower iteration cap than `sm.MAX_ITERATIONS`) with its MultiNest folder; the calibration-mock file when its truths lack the model's calibration labels; then the `jobs/calib_*.npz` files of those mocks (the new mocks reuse the `calib_<k>` ids); likewise a cross-fit mock file (`cross_fit_mocks.joblib`, section 6c) made for the other folder, and then its `jobs/cross_calib_*.npz`; and the run logs. A name already taken in `superseded/` gets the suffix `_16bin`, so nothing is skipped or overwritten. The cell only moves files: run it where the fits are kept, because a notebook job moves its own copies and its moves are not fetched back.'''

SECTION_6B_MARKDOWN = r'''## 6b. Photometric Bagpipes fits (spec amendment A)
The ML is calibrated on synthetic galaxies but disagrees with the spectral fits on real ones. These fits ask whether it agrees with Bagpipes given only what it sees: the 12-band 3″ aperture photometry at the spectroscopic redshift, under the delayed-τ training prior (`sm.MODELS["delayed"].photometric_dims`) and the likelihood of the simulator that made the training set (Gaussian in magnitudes with the catalogue errors, no error floor, detected bands only, brightness offset marginalised analytically: `sm.fit_photometry`). Only the `validation` half is fitted until the locked test. `EXPORT_LOCKED_TEST_JOBS = True` (parameter cell) is set once, together with `RUN_LOCKED_TEST` of `ml_model.ipynb`. It adds the `test` half for both SFH models, and the reading cells then fill the test galaxies' `phot_*` and `phot_dblplaw_*` columns too.

**Calibration:** the spectral calibration mocks' own points and redshifts, as photometry at the aperture brightness of a reference galaxy near that redshift, noised once with `sm.add_splus_noise` (`sm.photometric_mock_jobs`). The fitter passes when at least `sm.MIN_COVERAGE_95` of the truths lie inside the 95% intervals. Jobs go to `photometric_fits/jobs/` and run on the PC (`FITS_SUBDIR=photometric_fits remote/submit_fits.sh`). Finished fits become `phot_<label>_p16/_median/_p84` for `sm.photometric_summary_labels("delayed")` (the calibration labels, `mw_age`, `log_sfr` and the labels shared with the DPL model), `phot_converged`, `phot_log_evidence` and `phot_log_evidence_error` (ln Z and MultiNest's √(H / n_live); fits saved before the DPL model stored no error, so `sm.fit_evidence` recomputes it from their posterior samples, a few minutes in all), and `phot_fitter_calibrated` in the reference CSV, which `ml_model.ipynb` compares with the ML.

The section needs the split, so it waits for the spectral fits: while no galaxy is in `validation`, it draws no mocks (their brightness donors would include the test half, and the stored mock file would keep them) and exports no jobs. Run the notebook again once the spectral fits are back.'''

SECTION_6B_JOBS = r'''NOISE_MODEL_FILE = "splus_noise_model.csv"           # the training set's S-PLUS error model (training_set.ipynb)
PHOTOMETRIC_MOCK_SEED = 1                            # brightness and noise of the delayed-tau photometric mocks
# Halves of the reference split whose photometry is fitted (both models): the test half only for the locked test
PHOTOMETRIC_HALVES = ["validation", "test"] if EXPORT_LOCKED_TEST_JOBS else ["validation"]
os.makedirs(PHOTOMETRIC_FITS_DIRECTORY, exist_ok=True)
noise_table = pd.read_csv(NOISE_MODEL_FILE)


def enough_detected_bands(magnitudes, errors):
    """Whether a photometry job can be fitted: at least sm.MIN_DETECTED_BANDS detected bands.

    Inputs: magnitudes, errors: [mag], 12 values in sm.BANDS order. Returns: bool.
    """
    return sm.detected_bands(np.asarray(magnitudes), np.asarray(errors)).sum() >= sm.MIN_DETECTED_BANDS


def reference_galaxy_jobs(half, sfh_model):
    """Photometric fit jobs of the galaxies of one half of the reference split, under one SFH model.

    Inputs: half: "validation", or "test" (only with EXPORT_LOCKED_TEST_JOBS: the test half stays locked until the
            final report); sfh_model: "delayed" or "dblplaw".
    Returns: list of job dicts as sm.fit_and_save takes: job_id (the specobjid), kind sm.PHOTOMETRY_JOB, sfh_model,
             dereddened 3" magnitudes and errors [mag] and the spectroscopic redshift. Galaxies with fewer than
             sm.MIN_DETECTED_BANDS detected bands are left out, with a printed line.
    """
    jobs = []
    for row_index in reference.index[reference["split"] == half]:
        magnitudes = reference.loc[row_index, aperture_columns].to_numpy(dtype=float)
        errors = reference.loc[row_index, aperture_error_columns].to_numpy(dtype=float)
        specobjid = str(reference.loc[row_index, "specobjid"])
        if not enough_detected_bands(magnitudes, errors):
            print(f"{specobjid}: fewer than {sm.MIN_DETECTED_BANDS} detected bands, not fitted ({sfh_model})")
            continue
        jobs.append({"job_id": specobjid, "kind": sm.PHOTOMETRY_JOB, "sfh_model": sfh_model,
                     "magnitudes": magnitudes, "errors": errors, "redshift": float(reference.loc[row_index, "z"])})
    return jobs


def fitted_galaxy_jobs(sfh_model):
    """Photometric fit jobs of every half in PHOTOMETRIC_HALVES under one SFH model, the validation half first.

    Inputs: sfh_model: "delayed" or "dblplaw".
    Returns: list of job dicts as reference_galaxy_jobs gives them.
    """
    jobs = []
    for half in PHOTOMETRIC_HALVES:
        jobs.extend(reference_galaxy_jobs(half, sfh_model))
    return jobs


# 1. Calibration mocks: drawn once and stored with their true labels, because their fits are cached on disk.
# Section 4 moves mocks of an earlier model away. They wait for the split (see the note above)
if os.path.exists(PHOTOMETRIC_MOCK_FILE):
    photometric_mock_jobs, photometric_truths = joblib.load(PHOTOMETRIC_MOCK_FILE)
    print(f"Reusing {len(photometric_mock_jobs)} stored photometric calibration mocks "
          f"from {PHOTOMETRIC_MOCK_FILE}")
elif not reference["split"].eq("validation").any():
    photometric_mock_jobs = []
    photometric_truths = {}
    print("No validation galaxy yet (the spectral fits make the split): no photometric mocks or jobs this run")
else:
    # Why refuse: new mocks would differ from the fitted ones (the brightness donors grow between passes), and
    # the old calib_*.h5 fits would be scored against the new truths without any error
    stale_fits = glob.glob(os.path.join(PHOTOMETRIC_FITS_DIRECTORY, f"{sm.CALIBRATION_ID_PREFIX}*.h5"))
    assert len(stale_fits) == 0, (f"{len(stale_fits)} calib_*.h5 fits exist but {PHOTOMETRIC_MOCK_FILE} is "
                                  "missing: bundle it (with the fits), or delete the fits to draw new mocks")
    # 1a. The spectral mocks' unit-cube points, redrawn from SEED as sm.draw_calibration_mocks drew them, and
    # their photometric coordinates (veldisp and noise_scale are spectral only)
    spectral_cubes = sm.draw_calibration_cubes(len(calibration_jobs), np.random.default_rng(SEED))
    delayed_photometric_dims = sm.MODELS["delayed"].photometric_dims
    photometric_cubes = np.zeros((len(calibration_jobs), len(delayed_photometric_dims)))
    for mock_index in range(len(calibration_jobs)):
        cube_by_name = dict(zip(sm.SPECTRAL_DIMS, spectral_cubes[mock_index]))
        for dimension_index, dimension in enumerate(delayed_photometric_dims):
            photometric_cubes[mock_index, dimension_index] = cube_by_name[dimension]
    mock_ids = [job["job_id"] for job in calibration_jobs]
    mock_redshifts = np.array([job["redshift"] for job in calibration_jobs])

    # 1b. Noiseless photometry of the same points for 10^10 M_sun formed; each must be its mock's own point
    noiseless_magnitudes, noiseless_truths = sm.noiseless_photometric_mocks(photometric_cubes, mock_redshifts,
                                                                            "delayed")
    for mock_id, noiseless_truth in zip(mock_ids, noiseless_truths):
        for label in ("age", "tau", "metallicity"):
            assert noiseless_truth[label] == calibration_truths[mock_id][label], (
                f"{mock_id}: the redrawn point is not the mock's")

    # 1c. Brightness of a reference galaxy near the mock's redshift (aperture r, as the fits see; the test half
    # stays untouched), then the training simulator's noise; the spectral mocks' truths move with the brightness
    mock_rng = np.random.default_rng(PHOTOMETRIC_MOCK_SEED)
    brightness_donors = reference[reference["split"] != "test"]
    _, brightness_offsets = sm.scale_to_observed_r(noiseless_magnitudes, mock_redshifts,
                                                   brightness_donors["r_aper"].to_numpy(),
                                                   brightness_donors["z"].to_numpy(), mock_rng)
    spectral_truths = {}
    for mock_id in mock_ids:
        spectral_truths[mock_id] = calibration_truths[mock_id]
    photometric_mock_jobs, photometric_truths = sm.photometric_mock_jobs(
        mock_ids, noiseless_magnitudes, mock_redshifts, brightness_offsets, spectral_truths, noise_table, mock_rng,
        "delayed")
    left_out = sorted(set(mock_ids) - set(photometric_truths))
    if left_out:
        print(f"Left out, fewer than {sm.MIN_DETECTED_BANDS} detected bands: {', '.join(left_out)}")
    joblib.dump((photometric_mock_jobs, photometric_truths), PHOTOMETRIC_MOCK_FILE)
    print(f"Made {len(photometric_mock_jobs)} photometric calibration mocks, stored in {PHOTOMETRIC_MOCK_FILE}")

# 2. One job per galaxy of the fitted halves (the test half only for the locked test): dereddened 3" photometry
# and spectroscopic redshift
photometric_galaxy_jobs = fitted_galaxy_jobs("delayed")

# 3. Job files, and how to run them
for job in [*photometric_mock_jobs, *photometric_galaxy_jobs]:
    sm.export_job(job, PHOTOMETRIC_JOBS_DIRECTORY)
print(f"Exported {len(photometric_mock_jobs)} calibration mocks and {len(photometric_galaxy_jobs)} galaxies "
      f"({' and '.join(PHOTOMETRIC_HALVES)} halves) to {PHOTOMETRIC_JOBS_DIRECTORY}/")
print("Run them on the PC:  FITS_SUBDIR=photometric_fits remote/submit_fits.sh")
print("then fetch them:     remote/fetch_fits.sh <job-id>")'''

SECTION_6B_CALIBRATION = r'''def photometric_calibration(sfh_model, fits_directory, truths):
    """Calibration test of one model's photometric fitter on its stored mocks; prints the table and the verdict.

    Inputs: sfh_model: "delayed" or "dblplaw"; fits_directory: the model's photometric fits folder (str);
            truths: dict mock job_id -> Series of the true labels (from the model's mock file).
    Returns: (fitter_calibrated, coverage table from sm.calibration_coverage). fitter_calibrated is True when the
             95% coverage of every calibration label of the model reaches sm.MIN_COVERAGE_95, False otherwise,
             and pd.NA while a mock is unfitted (photometric fits are fast, so there is no skip list) or when
             there are no mocks.
    """
    # 1. Posteriors of the mocks that have finished
    posteriors = {}
    n_converged = 0
    for job_id in truths:
        fit_file = os.path.join(fits_directory, f"{job_id}.h5")
        if os.path.exists(fit_file):
            posterior, metadata = sm.load_fit(fit_file)
            assert metadata.get("kind") == sm.PHOTOMETRY_JOB, f"{fit_file} is not a photometric fit"
            assert sm.fit_is_current(metadata, sfh_model), f"{fit_file} is from an earlier model or run: run section 4"
            posteriors[job_id] = posterior
            n_converged += int(bool(metadata["converged"]))

    # 2. Coverage per label and interval
    table = sm.calibration_coverage(posteriors, truths, sm.MODELS[sfh_model].calibration_labels)
    display(table.pivot(index="label", columns="level", values=["coverage", "n_jobs"]).round(2))

    # 3. Verdict once every mock is fitted; none before the mocks exist
    coverage_95 = table[table["level"] == 0.95].set_index("label")["coverage"]
    verdict_column = f"{sm.PHOTOMETRIC_PREFIXES[sfh_model]}fitter_calibrated"
    if len(truths) == 0 or len(posteriors) < len(truths):
        print(f"INCOMPLETE ({sfh_model}): {len(posteriors)} of {len(truths)} photometric calibration mocks fitted "
              f"({n_converged} converged); {verdict_column} stays empty until all are done")
        return pd.NA, table
    fitter_calibrated = bool((coverage_95 >= sm.MIN_COVERAGE_95).all())
    verdict = "PASS" if fitter_calibrated else "FAIL"
    print(f"{verdict} ({sfh_model}): 95% coverage {coverage_95.round(2).to_dict()} against the required "
          f"{sm.MIN_COVERAGE_95} ({len(posteriors)} mocks, {n_converged} converged)")
    return fitter_calibrated, table


photometric_fitter_calibrated, photometric_calibration_table = photometric_calibration(
    "delayed", PHOTOMETRIC_FITS_DIRECTORY, photometric_truths)'''

SECTION_6B_SUMMARY = r'''def photometric_fit_columns(fit_file, job, sfh_model):
    """Posterior and run columns of one photometric fit, checked to be a current fit of sfh_model.

    Inputs: fit_file: path of the fit (str); job: the photometry job the fit ran on (magnitudes and errors [mag]);
            sfh_model: the model of the fit's folder.
    Returns: (posterior DataFrame, dict with <prefix>converged (bool), <prefix>log_evidence and
             <prefix>log_evidence_error (dimensionless ln Z and its MultiNest uncertainty, sm.fit_evidence)),
             prefix being sm.PHOTOMETRIC_PREFIXES[sfh_model].
    """
    posterior, metadata = sm.load_fit(fit_file)
    assert metadata.get("kind") == sm.PHOTOMETRY_JOB, f"{fit_file} is not a photometric fit"
    assert sm.fit_is_current(metadata, sfh_model), f"{fit_file} is from an earlier model or run: run section 4"
    log_evidence, log_evidence_error = sm.fit_evidence(fit_file, job)
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    run_columns = {f"{prefix}converged": bool(metadata["converged"]), f"{prefix}log_evidence": log_evidence,
                   f"{prefix}log_evidence_error": log_evidence_error}
    return posterior, run_columns


def photometric_summary(sfh_model, fits_directory, galaxy_jobs):
    """One row per exported galaxy with a photometric fit of sfh_model: summaries, convergence and evidence.

    The stored labels are already at the galaxy's brightness: each sample's mass-like labels moved by its own
    drawn offset.
    Inputs: sfh_model: "delayed" or "dblplaw"; fits_directory: the model's photometric fits folder (str);
            galaxy_jobs: the model's galaxy jobs of the fitted halves (fitted_galaxy_jobs).
    Returns: DataFrame with specobjid, then <prefix><label>_p16/_median/_p84 for sm.photometric_summary_labels,
             then <prefix>converged, <prefix>log_evidence and <prefix>log_evidence_error; galaxies without a fit
             have no row.
    """
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    labels = sm.photometric_summary_labels(sfh_model)
    empty_summary = pd.DataFrame(np.zeros((1, len(labels))), columns=labels)
    columns = ["specobjid"]
    for name in sm.summarise_posterior(empty_summary, labels):
        columns.append(f"{prefix}{name}")
    for name in ("converged", "log_evidence", "log_evidence_error"):
        columns.append(f"{prefix}{name}")

    rows = []
    for job in galaxy_jobs:
        fit_file = os.path.join(fits_directory, f"{job['job_id']}.h5")
        if not os.path.exists(fit_file):
            continue
        posterior, run_columns = photometric_fit_columns(fit_file, job, sfh_model)
        row = {"specobjid": job["job_id"], **run_columns}
        for name, value in sm.summarise_posterior(posterior, labels).items():
            row[f"{prefix}{name}"] = value
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def with_photometric_summary(table, summary, sfh_model, fitter_calibrated):
    """The reference table with one model's photometric columns replaced (re-running adds no copies).

    Why the prefix check: "phot_" (delayed-tau) also starts "phot_dblplaw_", whose columns must stay.
    Inputs: table: the reference DataFrame; summary: photometric_summary of sfh_model; sfh_model: its model;
            fitter_calibrated: the model's calibration verdict (bool or pd.NA).
    Returns: new DataFrame: table without its earlier columns of this model, merged with summary on specobjid
             (NaN for galaxies without a fit), with <prefix>fitter_calibrated on every row.
    """
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    other_prefixes = []
    for other_model, other_prefix in sm.PHOTOMETRIC_PREFIXES.items():
        if other_model != sfh_model and other_prefix.startswith(prefix):
            other_prefixes.append(other_prefix)
    earlier_columns = []
    for column in table.columns:
        belongs_to_other_model = any(column.startswith(other_prefix) for other_prefix in other_prefixes)
        if column.startswith(prefix) and not belongs_to_other_model:
            earlier_columns.append(column)
    merged = table.drop(columns=earlier_columns).merge(summary, on="specobjid", how="left")
    merged[f"{prefix}fitter_calibrated"] = fitter_calibrated
    return merged


# Delayed-tau photometric fits of the exported galaxies (validation half; test half too for the locked test),
# into the reference table
delayed_photometric_summary = photometric_summary("delayed", PHOTOMETRIC_FITS_DIRECTORY, photometric_galaxy_jobs)
reference = with_photometric_summary(reference, delayed_photometric_summary, "delayed",
                                     photometric_fitter_calibrated)
n_photometric_converged = int(delayed_photometric_summary["phot_converged"].sum())
print(f"Delayed-tau photometric fits read for {len(delayed_photometric_summary)} of {len(photometric_galaxy_jobs)} "
      f"exported galaxies ({' and '.join(PHOTOMETRIC_HALVES)} halves; {n_photometric_converged} converged); "
      f"phot_fitter_calibrated = {photometric_fitter_calibrated}")'''

SECTION_6C_MARKDOWN = r'''## 6c. Double-power-law photometric fits and cross-fits (DPL spec, part 2)
Bagpipes' double power law (`dblplaw`, SFR ∝ [(t/τ)^α + (t/τ)^−β]^−1, t and τ counted from the Big Bang) is fitted to the same photometry of the same galaxies (the `validation` half, plus the `test` half with `EXPORT_LOCKED_TEST_JOBS`), with the same likelihood, in its own folder `photometric_fits_dblplaw/`. τ is the turnover time of the double power law. It is the time of the SFR maximum only when α = β; in general the maximum lies at t = τ (β/α)^(1/(α+β)). `turnover_lookback` is the lookback time from the observation back to τ. There are no DPL spectral mocks, so its **calibration mocks** are 33 points of a Latin hypercube over the DPL prior, each at the redshift and aperture r brightness of a random validation galaxy, noised once (`sm.draw_photometric_calibration_mocks`). The fitter passes when the 95% coverage is at least `sm.MIN_COVERAGE_95` for `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity` and `log_mstar`; the 68% coverage is reported too. Results become `phot_dblplaw_<label>_p16/_median/_p84` (`sm.photometric_summary_labels("dblplaw")`), `phot_dblplaw_converged`, `phot_dblplaw_log_evidence`, `phot_dblplaw_log_evidence_error` and `phot_dblplaw_fitter_calibrated`.

**Cross-fits:** each model's mocks are also fitted with the other model, as `cross_calib_<k>` jobs. The delayed-τ mocks are fitted in `photometric_fits_dblplaw/`, and the DPL mocks in `photometric_fits/`. Each folder stores its cross-fit jobs with their truths in `cross_fit_mocks.joblib`. Since the truth is known, these fits show whether the colours can tell the two models apart. Their evidences go to `photometric_cross_fits.csv`, with one row per mock and both models' ln Z, for `ml_model.ipynb`.

**Running:** run `FITS_SUBDIR=photometric_fits_dblplaw remote/submit_fits.sh` first (DPL mocks, galaxies and the delayed-τ mocks' cross-fits). Then run `FITS_SUBDIR=photometric_fits remote/submit_fits.sh`. It fits only the DPL mocks' cross-fits (and, with `EXPORT_LOCKED_TEST_JOBS`, the test-half galaxies), because the finished delayed-τ fits travel in its bundle. Fetch each run with `remote/fetch_fits.sh <job-id>`. Like section 6b, this section waits for the split.'''

SECTION_6C_JOBS = r'''DBLPLAW_FITS_DIRECTORY = sm.model_file(PHOTOMETRIC_FITS_DIRECTORY, "dblplaw")    # photometric_fits_dblplaw
DBLPLAW_JOBS_DIRECTORY = os.path.join(DBLPLAW_FITS_DIRECTORY, sm.JOBS_SUBDIRECTORY)
DBLPLAW_MOCK_FILE = os.path.join(DBLPLAW_FITS_DIRECTORY, sm.CALIBRATION_MOCK_FILE)
DBLPLAW_MOCK_SEED = 2            # Latin hypercube, brightness donors and noise of the DPL calibration mocks
CROSS_FITS_FILE = "photometric_cross_fits.csv"     # ln Z of every mock under both models, for ml_model.ipynb
os.makedirs(DBLPLAW_FITS_DIRECTORY, exist_ok=True)

# 1. DPL calibration mocks: drawn once from the DPL prior and stored, as the delayed-tau ones; each takes the
# redshift and aperture r of a random validation galaxy. They wait for the split too
if os.path.exists(DBLPLAW_MOCK_FILE):
    dblplaw_mock_jobs, dblplaw_truths = joblib.load(DBLPLAW_MOCK_FILE)
    print(f"Reusing {len(dblplaw_mock_jobs)} stored DPL calibration mocks from {DBLPLAW_MOCK_FILE}")
elif not reference["split"].eq("validation").any():
    dblplaw_mock_jobs = []
    dblplaw_truths = {}
    print("No validation galaxy yet (the spectral fits make the split): no DPL mocks or jobs this run")
else:
    stale_fits = glob.glob(os.path.join(DBLPLAW_FITS_DIRECTORY, f"{sm.CALIBRATION_ID_PREFIX}*.h5"))
    assert len(stale_fits) == 0, (f"{len(stale_fits)} calib_*.h5 fits exist but {DBLPLAW_MOCK_FILE} is "
                                  "missing: bundle it (with the fits), or delete the fits to draw new mocks")
    validation_galaxies = reference[reference["split"] == "validation"]
    dblplaw_mock_jobs, dblplaw_truths = sm.draw_photometric_calibration_mocks(
        sm.N_CALIBRATION_MOCKS, "dblplaw", validation_galaxies["z"].to_numpy(),
        validation_galaxies["r_aper"].to_numpy(), noise_table, DBLPLAW_MOCK_SEED)
    joblib.dump((dblplaw_mock_jobs, dblplaw_truths), DBLPLAW_MOCK_FILE)
    n_left_out = sm.N_CALIBRATION_MOCKS - len(dblplaw_mock_jobs)
    print(f"Made {len(dblplaw_mock_jobs)} DPL calibration mocks ({n_left_out} left out for fewer than "
          f"{sm.MIN_DETECTED_BANDS} detected bands), stored in {DBLPLAW_MOCK_FILE}")

# 2. The same galaxies as the delayed-tau photometric fits (the halves of PHOTOMETRIC_HALVES)
dblplaw_galaxy_jobs = fitted_galaxy_jobs("dblplaw")

# 3. Job files
for job in [*dblplaw_mock_jobs, *dblplaw_galaxy_jobs]:
    sm.export_job(job, DBLPLAW_JOBS_DIRECTORY)
print(f"Exported {len(dblplaw_mock_jobs)} DPL calibration mocks and {len(dblplaw_galaxy_jobs)} galaxies "
      f"({' and '.join(PHOTOMETRIC_HALVES)} halves) to {DBLPLAW_JOBS_DIRECTORY}/")'''

SECTION_6C_CROSS_JOBS = r'''def stored_or_new_cross_fit_mocks(source_jobs, source_truths, fitted_model, fits_directory):
    """Cross-fit mocks of one fits folder: the stored ones, or new ones made from the other model's mocks and stored.

    Why stored: their fits are cached on disk and must be compared with the truths of the photometry they fitted,
    even if the source mocks are redrawn later.
    Inputs: source_jobs, source_truths: the other model's calibration mocks (jobs and truths, as in its mock file);
            fitted_model: this folder's model; fits_directory: this folder (str).
    Returns: (jobs, truths) as sm.cross_fit_jobs returns them; empty while the source mocks do not exist yet.
    """
    mock_file = os.path.join(fits_directory, sm.CROSS_FIT_MOCK_FILE)
    if os.path.exists(mock_file):
        jobs, truths = joblib.load(mock_file)
        print(f"Reusing {len(jobs)} stored cross-fit mocks from {mock_file}")
        return jobs, truths
    if len(source_jobs) == 0:
        print(f"No source mocks yet: no cross-fits for {fits_directory} this run")
        return [], {}
    stale_fits = glob.glob(os.path.join(fits_directory, f"{sm.CROSS_CALIBRATION_ID_PREFIX}*.h5"))
    assert len(stale_fits) == 0, (f"{len(stale_fits)} cross-fits exist but {mock_file} is missing: bundle it "
                                  "(with the fits), or delete the fits to make new cross-fit mocks")
    jobs, truths = sm.cross_fit_jobs(source_jobs, source_truths, fitted_model)
    joblib.dump((jobs, truths), mock_file)
    print(f"Made {len(jobs)} cross-fit mocks ({fitted_model} fits of the other model's mocks), stored in {mock_file}")
    return jobs, truths


# 1. The delayed-tau mocks fitted with DPL, and the DPL mocks fitted with delayed-tau
dblplaw_cross_jobs, dblplaw_cross_truths = stored_or_new_cross_fit_mocks(
    photometric_mock_jobs, photometric_truths, "dblplaw", DBLPLAW_FITS_DIRECTORY)
delayed_cross_jobs, delayed_cross_truths = stored_or_new_cross_fit_mocks(
    dblplaw_mock_jobs, dblplaw_truths, "delayed", PHOTOMETRIC_FITS_DIRECTORY)

# 2. Job files, each in the folder of the model it is fitted with
for job in dblplaw_cross_jobs:
    sm.export_job(job, DBLPLAW_JOBS_DIRECTORY)
for job in delayed_cross_jobs:
    sm.export_job(job, PHOTOMETRIC_JOBS_DIRECTORY)
print(f"Exported {len(dblplaw_cross_jobs)} cross-fits to {DBLPLAW_JOBS_DIRECTORY}/ and {len(delayed_cross_jobs)} "
      f"to {PHOTOMETRIC_JOBS_DIRECTORY}/")
print("Run them on the PC:  FITS_SUBDIR=photometric_fits_dblplaw remote/submit_fits.sh")
print("              then:  FITS_SUBDIR=photometric_fits remote/submit_fits.sh   (the DPL mocks' cross-fits, and "
      "the test half when EXPORT_LOCKED_TEST_JOBS is True)")
print("and fetch each:      remote/fetch_fits.sh <job-id>")'''

SECTION_6C_READ = r'''# DPL calibration test and the DPL photometric fits of the exported galaxies, into the reference table
dblplaw_fitter_calibrated, dblplaw_calibration_table = photometric_calibration(
    "dblplaw", DBLPLAW_FITS_DIRECTORY, dblplaw_truths)
dblplaw_photometric_summary = photometric_summary("dblplaw", DBLPLAW_FITS_DIRECTORY, dblplaw_galaxy_jobs)
reference = with_photometric_summary(reference, dblplaw_photometric_summary, "dblplaw", dblplaw_fitter_calibrated)
n_dblplaw_converged = int(dblplaw_photometric_summary["phot_dblplaw_converged"].sum())
print(f"DPL photometric fits read for {len(dblplaw_photometric_summary)} of {len(dblplaw_galaxy_jobs)} exported "
      f"galaxies ({' and '.join(PHOTOMETRIC_HALVES)} halves; {n_dblplaw_converged} converged); "
      f"phot_dblplaw_fitter_calibrated = {dblplaw_fitter_calibrated}")'''

SECTION_6C_CROSS_TABLE = r'''# Every mock of either model under both models: its own-model fit (calib_<k> in its own folder) and its
# cross-fit (cross_calib_<k> in the other folder). The columns keep the reference CSV's names; ml_model.ipynb keeps
# the mocks whose two fits converged and renames the two ln Z columns for ev.true_model_recovery. A missing fit
# leaves its columns empty
fits_directory_of_model = {"delayed": PHOTOMETRIC_FITS_DIRECTORY, "dblplaw": DBLPLAW_FITS_DIRECTORY}
cross_fit_columns = ["mock_id", "truth_model"]
for prefix in sm.PHOTOMETRIC_PREFIXES.values():
    for name in ("converged", "log_evidence", "log_evidence_error"):
        cross_fit_columns.append(f"{prefix}{name}")

cross_fit_rows = []
for truth_model, mock_jobs in (("delayed", photometric_mock_jobs), ("dblplaw", dblplaw_mock_jobs)):
    for mock_job in mock_jobs:
        row = {"mock_id": mock_job["job_id"], "truth_model": truth_model}
        for fitted_model, fits_directory in fits_directory_of_model.items():
            if fitted_model == truth_model:
                job_id = mock_job["job_id"]
            else:
                job_id = sm.cross_fit_id(mock_job["job_id"])
            fit_file = os.path.join(fits_directory, f"{job_id}.h5")
            if os.path.exists(fit_file):
                _, run_columns = photometric_fit_columns(fit_file, mock_job, fitted_model)
                row.update(run_columns)
        cross_fit_rows.append(row)
cross_fits = pd.DataFrame(cross_fit_rows, columns=cross_fit_columns)
cross_fits.to_csv(CROSS_FITS_FILE, index=False)

both_fitted = cross_fits["phot_log_evidence"].notna() & cross_fits["phot_dblplaw_log_evidence"].notna()
print(f"{CROSS_FITS_FILE}: {len(cross_fits)} mocks, {int(both_fitted.sum())} fitted with both models")'''

# ---------------------------------------------------------------------------------------------------------------
# Locate every cell before changing any (first lines are those of today's notebook)
title_index = find_cell("# Spectroscopic reference set (part 2)")
parameters_index = find_cell("import glob")
section_4_index = find_cell("## 4. Fits (outside the notebook)")
section_6b_index = find_cell("## 6b. Photometric Bagpipes fits")
jobs_index = find_cell('NOISE_MODEL_FILE = "splus_noise_model.csv"')
calibration_index = find_cell("# 1. Posteriors of the photometric calibration mocks that have finished")
summary_index = find_cell("# 1. Summary columns, built once: phot_<label>_p16")
assert [section_6b_index, jobs_index, calibration_index, summary_index] == list(range(section_6b_index,
                                                                                         section_6b_index + 4))

parameters = source_of(cells[parameters_index])
assert parameters.count(PARAMETERS_OLD) == 1 and "EXPORT_LOCKED_TEST_JOBS" not in parameters
set_source(cells[parameters_index], parameters.replace(PARAMETERS_OLD, PARAMETERS_NEW))

title = source_of(cells[title_index])
assert title.count(OUTPUTS_OLD) == 1
set_source(cells[title_index], title.replace(OUTPUTS_OLD, OUTPUTS_NEW))

section_4 = source_of(cells[section_4_index])
superseding_start = section_4.index("**Superseding:**")
set_source(cells[section_4_index], section_4[:superseding_start] + SUPERSEDING)

set_source(cells[section_6b_index], SECTION_6B_MARKDOWN)
set_source(cells[jobs_index], SECTION_6B_JOBS)
set_source(cells[calibration_index], SECTION_6B_CALIBRATION)
set_source(cells[summary_index], SECTION_6B_SUMMARY)

new_cells = [new_cell("markdown", SECTION_6C_MARKDOWN), new_cell("code", SECTION_6C_JOBS),
             new_cell("code", SECTION_6C_CROSS_JOBS), new_cell("code", SECTION_6C_READ),
             new_cell("code", SECTION_6C_CROSS_TABLE)]
cells[summary_index + 1:summary_index + 1] = new_cells

with open(NOTEBOOK, "w") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print("notebook updated:", len(cells), "cells")
```

Expected output: `notebook updated: 36 cells`.

- [ ] **Step 2: Check the notebook (never execute it on the laptop)**

Run from `Synthetic_domain/`:

```bash
../.venv/bin/python -c "
import nbformat
nb = nbformat.read('spectroscopic_reference.ipynb', as_version=4); nbformat.validate(nb)
for i, c in enumerate(nb.cells):
    if c.cell_type == 'code': compile(c.source, f'cell {i}', 'exec')
print('ok', len(nb.cells))
print([c.source.splitlines()[0] for c in nb.cells if c.cell_type == 'markdown' and c.source.startswith('## ')])"
grep -c 'PHOTOMETRIC_LABELS\|sm.PHOTOMETRIC_DIMS\|fit_is_current(metadata)\|validation_galaxy_jobs' spectroscopic_reference.ipynb   # expect 0
grep -c 'EXPORT_LOCKED_TEST_JOBS = False' spectroscopic_reference.ipynb   # expect 1
```

Expected: `ok 36`. The headings list then has `## 6b. Photometric Bagpipes fits (spec amendment A)` directly followed by `## 6c. Double-power-law photometric fits and cross-fits (DPL spec, part 2)`, and then `## 7. Validation of the spectral fits`. The first `grep` prints `0` and the second `1`.

- [ ] **Step 3: Remote scripts, README and `.gitignore`**

`remote/submit_fits.sh`. Replace these header lines:

```
# FITS_SUBDIR picks the set of fits: spectral_fits (the spectral fits) or photometric_fits (the photometric fits of
# spec amendment A). The two share job ids (specobjids), so their jobs and results never mix:
#     FITS_SUBDIR=photometric_fits Synthetic_domain/remote/submit_fits.sh
```

with:

```
# FITS_SUBDIR picks the set of fits: spectral_fits (the spectral fits), photometric_fits (the delayed-tau photometric
# fits of spec amendment A) or photometric_fits_dblplaw (the double-power-law photometric fits and the cross-fits of
# the delayed-tau mocks). They share job ids (specobjids, calib_<k>), so their jobs and results never mix:
#     FITS_SUBDIR=photometric_fits_dblplaw Synthetic_domain/remote/submit_fits.sh
```

Leave the rest of the script as it is. Task 1 already changed line 58 to

```
cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/sfh_models.py" "$SYNTHETIC_DIR/run_spectral_fits.py" "$bundle_dir/"
```

and Step 4 below only checks that the bundle carries `sfh_models.py`.

`remote/fit_job.sh`. Replace `#   FITS_SUBDIR      folder of the jobs and results: spectral_fits or photometric_fits   [default spectral_fits]` with `#   FITS_SUBDIR      folder of the jobs and results: spectral_fits, photometric_fits or photometric_fits_dblplaw   [default spectral_fits]`.

`remote/fetch_fits.sh`. Replace `# The fits folder is the one the job recorded in out/fits_subdir.txt (spectral_fits or photometric_fits). Jobs` with `# The fits folder is the one the job recorded in out/fits_subdir.txt (spectral_fits, photometric_fits or photometric_fits_dblplaw). Jobs`.

`remote/README.md`. Replace:

```
    FITS_SUBDIR=photometric_fits Synthetic_domain/remote/submit_fits.sh

They share job ids (specobjids, `calib_<k>`) with the spectral fits, so they never share a folder: the
```

with:

```
    FITS_SUBDIR=photometric_fits Synthetic_domain/remote/submit_fits.sh

The double-power-law photometric fits (section 6c of notebook 2: the DPL calibration mocks, the validation galaxies
and the delayed-tau mocks' cross-fits `cross_calib_<k>`) have their own folder:

    FITS_SUBDIR=photometric_fits_dblplaw Synthetic_domain/remote/submit_fits.sh

The DPL mocks' cross-fits (fitted with delayed-tau) are exported to `photometric_fits/jobs/`, so a
`FITS_SUBDIR=photometric_fits` run after it fits only those: the finished delayed-tau fits travel in its bundle.

They share job ids (specobjids, `calib_<k>`) with the spectral fits, so they never share a folder: the
```

At the end of `remote/README.md`, after the line `but the mock file is missing.`, add:

```

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
```

`.gitignore` (repo root). After the line `Synthetic_domain/photometric_fits/`, add `Synthetic_domain/photometric_fits_dblplaw/`.

- [ ] **Step 4: Check the scripts**

Run from the repository root:

```bash
bash -n Synthetic_domain/remote/submit_fits.sh && bash -n Synthetic_domain/remote/fit_job.sh \
  && bash -n Synthetic_domain/remote/fetch_fits.sh && echo syntax-ok
BUNDLE_ONLY="$SCRATCH/bundle_check" FITS_SUBDIR=photometric_fits Synthetic_domain/remote/submit_fits.sh
tar -tzf "$SCRATCH/bundle_check/bundle.tar.gz" | grep -E '^\./(sfh_models|synthetic_model|run_spectral_fits)\.py$'
FITS_SUBDIR=photometric_fits_dblplaw BUNDLE_ONLY="$SCRATCH/bundle_check" Synthetic_domain/remote/submit_fits.sh; echo "exit $?"
git check-ignore Synthetic_domain/photometric_fits_dblplaw/jobs/x.npz
```

Set `SCRATCH` to your scratchpad directory first. `BUNDLE_ONLY` writes the bundle and submits nothing.

Expected:
- `syntax-ok`;
- a `Bundle (photometric_fits): ...` line;
- the three module paths (`sfh_models.py` comes from Task 1's line 58);
- `No job folder .../photometric_fits_dblplaw/jobs (FITS_SUBDIR=photometric_fits_dblplaw)` with `exit 1` (the folder exists only after notebook 2 has run on the PC);
- the ignored path printed.

If `~/Programs/MultiNest` is missing on this machine, skip the two `BUNDLE_ONLY` lines and rely on `bash -n`.

- [ ] **Step 5: Run the module tests (nothing in them reads the notebook, but the whole suite must stay green)**

Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests -q`

Expected: `220 passed`, the same as after Task 5, because this task adds no test. The result is `219 passed, 1 skipped` when the short MultiNest run of `test_recomputed_evidence_error_agrees_with_multinest` does not converge.

- [ ] **Step 6: Commit**

```bash
git add Synthetic_domain/spectroscopic_reference.ipynb Synthetic_domain/remote/submit_fits.sh \
        Synthetic_domain/remote/fit_job.sh Synthetic_domain/remote/fetch_fits.sh Synthetic_domain/remote/README.md \
        .gitignore
git commit -m "feat: notebook 2 fits the DPL photometry, its mocks and the cross-fits; evidences in the reference CSV

Section 6b reads both models through shared helpers and adds phot_log_evidence(_error) and the
shared labels; new section 6c exports photometric_fits_dblplaw/jobs and the cross-fits, reads the
DPL fits into phot_dblplaw_* and writes photometric_cross_fits.csv. EXPORT_LOCKED_TEST_JOBS
(False) adds the test half's photometric jobs of both models for the locked test. The remote
scripts document FITS_SUBDIR=photometric_fits_dblplaw.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 7: `ml_data` encodes and decodes the targets of either SFH model

**Files:**
- Modify: `Synthetic_domain/ml_data.py`:
  - module docstring, imports and section 1 (lines 1–37);
  - section 3, "Targets and decoding" (lines 87–152, rewritten);
  - section 4: `simulate_observations` (lines 170–197) and `validation_observations` (lines 200–229) take `sfh_model`.
- Modify: `Synthetic_domain/train_ml_model.py` lines 95–96, 108 and 114. Each call passes `"delayed"`, so the CLI keeps working until Task 8 adds `--sfh-model`.
- Test: `Synthetic_domain/tests/test_ml_data.py` (full rewrite below).

**Interfaces:**
- **Consumes (Tasks 1–3):**
  - `sfh_models.MODELS[name]` with `.dims`, `.drawn_labels`, `.log_labels`, `.analytic_labels(drawn_columns, log_mstar)`, `.target_names` and `.training_labels`;
  - the constants `sfh_models.Z_RANGE`, `AGE_MIN`, `TAU_RANGE`, `METALLICITY_RANGE`, `TURNOVER_LOOKBACK_MIN`, `ALPHA_RANGE`, `BETA_RANGE`, `COSMOLOGY` and `UNIVERSE_AGE_MARGIN`;
  - the functions `sfh_models.maximum_age(redshift)` and `sfh_models.model_file(name, sfh_model)`;
  - `sm.simulate_models(unit_vectors, sfh_model)`, which returns the `MODELS[sfh_model].training_labels` columns and then `sm.BANDS`.
- **Produces:**
  - `ml_data.EXTRA_LOG_TARGETS = {"log_tau_bc": "tau_bc"}`;
  - `ml_data.target_names(sfh_model) -> list` (a new list each call);
  - `ml_data.encode_targets(table, sfh_model) -> ndarray (n, len(target_names(sfh_model)))`;
  - `ml_data.delayed_sfh_parameters(targets, redshift) -> dict` with keys age and tau;
  - `ml_data.dblplaw_sfh_parameters(targets, redshift) -> dict` with keys turnover_lookback, tau, alpha and beta;
  - `ml_data.SFH_PARAMETER_DECODERS = {"delayed": ..., "dblplaw": ...}`;
  - `ml_data.decode_rows(theta, sfh_model) -> DataFrame`: the decode itself, in this process, without checks;
  - `ml_data.decode_samples(theta, sfh_model, n_jobs=1) -> DataFrame`, with columns `MODELS[sfh_model].training_labels`. It raises `ValueError` when theta has the wrong number of columns or `n_jobs < 1`. With `n_jobs > 1` it splits the rows into `n_jobs` blocks, decodes them in threads (`joblib.Parallel(prefer="threads")`) and concatenates them in order, so the result equals that of `n_jobs=1` to rounding. Task 10 uses this: a DPL decode sums Bagpipes' age grid for every sample (about 0.06 ms each, Task 2 Step 5).
  - `ml_data.simulate_observations(table, catalogue_r, catalogue_z, noise_table, rng, sfh_model)`;
  - `ml_data.validation_observations(held_out, mode, rng, sfh_model)`.
  - **Temporary:** `ml_data.TARGET_NAMES` and `ml_data.N_TARGETS` stay as delayed-τ aliases. `ml_networks.DEFAULT_CONFIG` reads them at import, so they remain until Task 8 removes them.

- [ ] **Step 1: Write the failing tests.** Replace the whole of `Synthetic_domain/tests/test_ml_data.py` with the file below. The existing delayed-τ tests keep every assertion; they change only through the `sfh_model` argument and `ml_data.target_names("delayed")` in place of `TARGET_NAMES`/`N_TARGETS`, and are parametrised over both models where the check is model-independent. Do not loosen them.

```python
"""Tests for ml_data.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest tests/test_ml_data.py -q"""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy.stats import qmc

import ml_data
import sfh_models
import synthetic_model as sm

VALIDATION_FILE = "synthetic_validation.parquet"
ROUND_TRIP_TOLERANCE = 1e-9      # the decode uses the same label functions as the training labels: exact
# A parallel block's grid sum may span fewer grid ages than the serial one (its oldest redshift can differ),
# so the two decodes agree to rounding, not bit for bit
PARALLEL_TOLERANCE = 1e-12
SFH_MODEL_NAMES = ["delayed", "dblplaw"]
N_SIMULATED = 200                # prior draws per model in the simulated fixture
DELAYED_TARGETS = ["log_age", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU",
                   "log_mass_formed", "log_mstar", "redshift"]
DBLPLAW_TARGETS = ["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "av_ism", "log_tau_bc",
                   "logU", "log_mass_formed", "log_mstar", "redshift"]
# Target values inside both priors at z = 0.1, for the targets a test row does not set
TYPICAL_TARGETS = {"log_age": 0.0, "log_tau": 0.0, "log_turnover_lookback": 0.0, "log_alpha": 0.0,
                   "log_beta": 0.0, "log_metallicity": 0.0, "av_ism": 0.5, "log_tau_bc": 0.0, "logU": -2.5,
                   "log_mass_formed": 10.0, "log_mstar": 9.8, "redshift": 0.1}
# A column only that model's sets hold: the set has been generated for the model
SET_COLUMNS = {"delayed": "age", "dblplaw": "turnover_lookback"}


def parquet_has_column(path, column):
    """True when the parquet file exists and its schema holds the column (read without loading the data)."""
    try:
        return column in pq.read_schema(path).names
    except OSError:
        return False


@pytest.fixture(scope="module")
def simulated_sets():
    """N_SIMULATED prior draws (Latin hypercube) of each SFH model, with labels and noiseless magnitudes."""
    sets = {}
    for sfh_model in SFH_MODEL_NAMES:
        n_dims = len(sfh_models.MODELS[sfh_model].dims)
        unit_vectors = qmc.LatinHypercube(d=n_dims, seed=30).random(N_SIMULATED)
        sets[sfh_model] = sm.simulate_models(unit_vectors, sfh_model)
    return sets


def theta_rows(rows, sfh_model):
    """Flow targets (n, n_targets) of sfh_model from dicts keyed by target name; unset targets are typical."""
    names = ml_data.target_names(sfh_model)
    theta = np.zeros((len(rows), len(names)))
    for row_index, row in enumerate(rows):
        values = dict(TYPICAL_TARGETS, **row)
        for target_index, name in enumerate(names):
            theta[row_index, target_index] = values[name]
    return theta


def expected_targets(table, sfh_model):
    """The flow targets of each row, written out by hand: dict target name -> array."""
    expected = {
        "log_metallicity": np.log10(table["metallicity"].to_numpy(dtype=float)),
        "av_ism": table["av_ism"].to_numpy(dtype=float),
        "log_tau_bc": np.log10(table["tau_bc"].to_numpy(dtype=float)),
        "logU": table["logU"].to_numpy(dtype=float),
        "log_mass_formed": table["log_mass_formed"].to_numpy(dtype=float),
        "log_mstar": table["log_mstar"].to_numpy(dtype=float),
        "redshift": table["redshift"].to_numpy(dtype=float),
    }
    if sfh_model == "delayed":
        expected["log_age"] = np.log10(table["age"].to_numpy(dtype=float))
        expected["log_tau"] = np.log10(table["tau"].to_numpy(dtype=float))
    else:
        expected["log_turnover_lookback"] = np.log10(table["turnover_lookback"].to_numpy(dtype=float))
        expected["log_alpha"] = np.log10(table["alpha"].to_numpy(dtype=float))
        expected["log_beta"] = np.log10(table["beta"].to_numpy(dtype=float))
    return expected


def detections(n_rows, detected_bands):
    """Magnitudes 20 mag with 0.05 mag errors in the listed bands, undetected (NaN) elsewhere."""
    magnitudes = np.full((n_rows, len(sm.BANDS)), np.nan)
    errors = np.full((n_rows, len(sm.BANDS)), np.nan)
    for band in detected_bands:
        band_index = sm.BANDS.index(band)
        magnitudes[:, band_index] = 20.0 + 0.1 * band_index
        errors[:, band_index] = 0.05
    return magnitudes, errors


def test_colours_do_not_change_with_brightness():
    magnitudes, errors = detections(3, sm.BANDS)
    z_obs = np.full(3, 0.1)
    sigma_z = np.zeros(3)
    bright = ml_data.colour_inputs(magnitudes, errors, z_obs, sigma_z)
    faint = ml_data.colour_inputs(magnitudes + 1.5, errors, z_obs, sigma_z)
    assert np.allclose(bright["flat"][:, :33], faint["flat"][:, :33])
    r_column = 33
    assert np.allclose(faint["flat"][:, r_column] - bright["flat"][:, r_column], 1.5)


def test_colour_inputs_without_r_band_is_invalid():
    bands_without_r = [band for band in sm.BANDS if band != "r"]
    magnitudes, errors = detections(2, bands_without_r)
    inputs = ml_data.colour_inputs(magnitudes, errors, np.full(2, 0.1), np.zeros(2))
    assert not inputs["valid"].any()


def test_colour_inputs_only_r_detected_has_no_nan():
    magnitudes, errors = detections(2, ["r"])
    inputs = ml_data.colour_inputs(magnitudes, errors, np.full(2, 0.1), np.zeros(2))
    assert inputs["valid"].all()
    assert np.isfinite(inputs["flat"]).all()
    assert np.isfinite(inputs["banded"]).all()
    assert inputs["flat"][:, 22:33].sum() == 0      # every colour flag is 0


def test_target_names_of_each_model_end_with_the_redshift():
    assert ml_data.target_names("delayed") == DELAYED_TARGETS
    assert ml_data.target_names("dblplaw") == DBLPLAW_TARGETS
    names = ml_data.target_names("dblplaw")
    names.append("extra")
    assert ml_data.target_names("dblplaw") == DBLPLAW_TARGETS      # a copy: a caller cannot change the model


def test_target_names_of_an_unknown_model_raise():
    with pytest.raises(KeyError):
        ml_data.target_names("exponential")


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_encode_targets_takes_logs_of_the_log_uniform_parameters(simulated_sets, sfh_model):
    simulated = simulated_sets[sfh_model]
    theta = ml_data.encode_targets(simulated, sfh_model)
    names = ml_data.target_names(sfh_model)
    assert theta.shape == (len(simulated), len(names))
    expected = expected_targets(simulated, sfh_model)
    for target_index, name in enumerate(names):
        assert np.array_equal(theta[:, target_index], expected[name]), name


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_decoded_labels_equal_the_summary_labels_of_simulated_models(simulated_sets, sfh_model):
    """Encode then decode the simulated truths: every training label (sm.summary_labels made them) comes back."""
    simulated = simulated_sets[sfh_model]
    decoded = ml_data.decode_samples(ml_data.encode_targets(simulated, sfh_model), sfh_model)
    training_labels = sfh_models.MODELS[sfh_model].training_labels
    assert list(decoded.columns) == training_labels
    for label in training_labels:
        expected = simulated[label].to_numpy(dtype=float)
        assert np.allclose(decoded[label], expected, rtol=0, atol=ROUND_TRIP_TOLERANCE), label


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_parallel_decode_equals_the_serial_decode(simulated_sets, sfh_model):
    """n_jobs > 1 decodes blocks of rows in threads: the rows, their order and the labels stay."""
    theta = ml_data.encode_targets(simulated_sets[sfh_model], sfh_model)
    serial = ml_data.decode_samples(theta, sfh_model)
    parallel = ml_data.decode_samples(theta, sfh_model, n_jobs=3)
    assert list(parallel.columns) == list(serial.columns)
    assert list(parallel.index) == list(serial.index)
    assert np.allclose(parallel.to_numpy(), serial.to_numpy(), rtol=PARALLEL_TOLERANCE, atol=PARALLEL_TOLERANCE)
    few_rows = ml_data.decode_samples(theta[:2], sfh_model, n_jobs=4)      # fewer rows than threads
    assert np.allclose(few_rows.to_numpy(), serial.iloc[:2].to_numpy(), rtol=PARALLEL_TOLERANCE,
                       atol=PARALLEL_TOLERANCE)


def test_decode_refuses_fewer_than_one_job():
    with pytest.raises(ValueError, match="n_jobs"):
        ml_data.decode_samples(theta_rows([{}], "delayed"), "delayed", n_jobs=0)


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_decode_refuses_theta_of_the_other_model(sfh_model):
    other_model = "dblplaw" if sfh_model == "delayed" else "delayed"
    theta = theta_rows([{}], other_model)
    with pytest.raises(ValueError, match="targets"):
        ml_data.decode_samples(theta, sfh_model)


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_decode_clips_redshift_to_prior(simulated_sets, sfh_model):
    theta = ml_data.encode_targets(simulated_sets[sfh_model].iloc[:2], sfh_model)
    theta[0, -1] = -0.3
    theta[1, -1] = 1.5
    decoded = ml_data.decode_samples(theta, sfh_model)
    assert decoded["redshift"].between(*sfh_models.Z_RANGE).all()


def test_decode_clips_age_to_the_oldest_age_at_the_sample_redshift():
    log_age_too_old = np.log10(10.0)        # older than t_max at z = 0.7 (7.07 Gyr), not at z = 0.1 (12.07)
    theta = theta_rows([
        {"log_age": log_age_too_old, "redshift": 0.7},
        {"log_age": log_age_too_old, "redshift": 0.1},
        {"log_age": np.log10(sfh_models.maximum_age(0.7)), "redshift": 0.7},
        {"log_age": np.log10(0.01), "redshift": 0.1},
        {"log_age": log_age_too_old, "redshift": 1.5},
    ], "delayed")
    decoded = ml_data.decode_samples(theta, "delayed")
    assert decoded["age"][0] <= sfh_models.maximum_age(0.7)
    assert np.isclose(decoded["age"][0], sfh_models.maximum_age(0.7), rtol=1e-12)
    assert np.isclose(decoded["age"][1], 10.0, rtol=1e-12)
    assert np.isclose(decoded["age"][2], sfh_models.maximum_age(0.7), rtol=1e-12)
    assert decoded["age"][3] == sfh_models.AGE_MIN
    assert np.isclose(decoded["age"][4], sfh_models.maximum_age(sfh_models.Z_RANGE[1]), rtol=1e-12)
    assert np.allclose(decoded["log_age"], np.log10(decoded["age"]))


def test_decode_clips_tau_and_metallicity_to_their_ranges():
    theta = theta_rows([
        {"log_tau": -5.0, "log_metallicity": -6.0},
        {"log_tau": 5.0, "log_metallicity": 3.0},
    ], "delayed")
    decoded = ml_data.decode_samples(theta, "delayed")
    assert list(decoded["tau"]) == list(sfh_models.TAU_RANGE)
    assert list(decoded["metallicity"]) == list(sfh_models.METALLICITY_RANGE)
    assert np.allclose(decoded["log_tau"], np.log10(sfh_models.TAU_RANGE))
    assert np.allclose(decoded["log_metallicity"], np.log10(sfh_models.METALLICITY_RANGE))


def test_decode_clips_turnover_lookback_to_the_prior_at_the_sample_redshift():
    log_too_old = np.log10(10.0)            # older than the oldest allowed turnover at z = 0.7 (7.07 Gyr ago), not at z = 0.1
    theta = theta_rows([
        {"log_turnover_lookback": log_too_old, "redshift": 0.7},
        {"log_turnover_lookback": log_too_old, "redshift": 0.1},
        {"log_turnover_lookback": np.log10(0.001), "redshift": 0.1},
        {"log_turnover_lookback": log_too_old, "redshift": 1.5},
    ], "dblplaw")
    decoded = ml_data.decode_samples(theta, "dblplaw")
    assert np.isclose(decoded["turnover_lookback"][0], sfh_models.maximum_age(0.7), rtol=1e-12)
    assert np.isclose(decoded["turnover_lookback"][1], 10.0, rtol=1e-12)
    assert decoded["turnover_lookback"][2] == sfh_models.TURNOVER_LOOKBACK_MIN
    assert np.isclose(decoded["turnover_lookback"][3], sfh_models.maximum_age(sfh_models.Z_RANGE[1]), rtol=1e-12)
    # Bagpipes' tau is the turnover time since the Big Bang, so it follows the clipped turnover lookback
    universe_age = sfh_models.COSMOLOGY.age(decoded["redshift"].to_numpy()).value
    assert np.allclose(decoded["tau"], universe_age - decoded["turnover_lookback"], rtol=0, atol=1e-12)
    assert (decoded["tau"] >= sfh_models.UNIVERSE_AGE_MARGIN - 1e-12).all()
    assert np.allclose(decoded["log_turnover_lookback"], np.log10(decoded["turnover_lookback"]))


def test_decode_clips_alpha_beta_and_metallicity_to_their_ranges():
    theta = theta_rows([
        {"log_alpha": -5.0, "log_beta": 6.0, "log_metallicity": -6.0},
        {"log_alpha": 6.0, "log_beta": -5.0, "log_metallicity": 3.0},
    ], "dblplaw")
    decoded = ml_data.decode_samples(theta, "dblplaw")
    assert list(decoded["alpha"]) == list(sfh_models.ALPHA_RANGE)
    assert list(decoded["beta"]) == [sfh_models.BETA_RANGE[1], sfh_models.BETA_RANGE[0]]
    assert list(decoded["metallicity"]) == list(sfh_models.METALLICITY_RANGE)
    assert np.allclose(decoded["log_alpha"], np.log10(decoded["alpha"]))
    assert np.allclose(decoded["log_beta"], np.log10(decoded["beta"]))


def test_decode_extreme_delayed_samples_stay_finite():
    theta = theta_rows([
        {"log_age": 50.0, "log_tau": -50.0, "redshift": 0.7},
        {"log_age": -50.0, "log_tau": 50.0, "redshift": 0.002},
        {"log_age": 50.0, "log_tau": 50.0, "log_metallicity": 40.0},
        {"log_age": -50.0, "log_tau": -50.0, "log_metallicity": -40.0},
    ], "delayed")
    with np.errstate(divide="raise", over="raise", invalid="raise"):
        decoded = ml_data.decode_samples(theta, "delayed")
    assert np.isfinite(decoded.to_numpy()).all()
    fractions = decoded[["f_young", "f_mid", "f_old"]].to_numpy()
    assert np.all((fractions >= 0) & (fractions <= 1))
    assert np.allclose(fractions.sum(axis=1), 1.0)
    assert np.allclose(decoded["log_ssfr"], decoded["log_sfr"] - decoded["log_mstar"])


def test_decode_extreme_dblplaw_samples_stay_finite():
    """The corners of the DPL prior: sharp (alpha = beta = 1000) and flat (0.01) histories, oldest and latest turnovers."""
    theta = theta_rows([
        {"log_turnover_lookback": 50.0, "log_alpha": 50.0, "log_beta": -50.0, "redshift": 0.7},
        {"log_turnover_lookback": -50.0, "log_alpha": -50.0, "log_beta": 50.0, "redshift": 0.002},
        {"log_turnover_lookback": 50.0, "log_alpha": 50.0, "log_beta": 50.0, "log_metallicity": 40.0},
        {"log_turnover_lookback": -50.0, "log_alpha": 50.0, "log_beta": 50.0, "redshift": 0.7},
        {"log_turnover_lookback": -50.0, "log_alpha": -50.0, "log_beta": -50.0, "log_metallicity": -40.0},
    ], "dblplaw")
    decoded = ml_data.decode_samples(theta, "dblplaw")
    assert np.isfinite(decoded.to_numpy()).all()
    fractions = decoded[["f_young", "f_mid", "f_old"]].to_numpy()
    assert np.all((fractions >= 0) & (fractions <= 1))
    assert np.allclose(fractions.sum(axis=1), 1.0)
    assert np.allclose(decoded["log_ssfr"], decoded["log_sfr"] - decoded["log_mstar"])


@pytest.fixture(scope="module")
def noise_and_catalogue():
    catalogue = pd.read_csv("dr3_corrected.csv")
    catalogue = catalogue[catalogue["PROB_GAL"] > 0.8]
    noise_table = sm.fit_noise_model(catalogue)
    return noise_table, catalogue["r_auto"].to_numpy(), catalogue["zml"].to_numpy()


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_simulate_observations_shifts_mass_and_draws_redshift(simulated_sets, noise_and_catalogue, sfh_model):
    noise_table, catalogue_r, catalogue_z = noise_and_catalogue
    # Repeat the models 10 times: each copy gets its own brightness and noise draw, so enough
    # galaxies end up with a precise r (the check below needs at least 5)
    repeated_models = pd.concat([simulated_sets[sfh_model]] * 10, ignore_index=True)
    observations = ml_data.simulate_observations(repeated_models, catalogue_r, catalogue_z, noise_table,
                                                 np.random.default_rng(31), sfh_model)
    targets = observations["targets"]
    names = ml_data.target_names(sfh_model)
    assert targets.shape == (len(repeated_models), len(names))
    flat = observations["flat"]
    log_mass_formed_index = names.index("log_mass_formed")
    r_magnitude = flat[:, 33]
    model_r = repeated_models["r"].to_numpy()
    r_error = flat[:, 34]
    precise = observations["valid"] & (r_error < 0.05)
    expected_mass = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * (r_magnitude - model_r)
    # The noisy r differs from the scaled r by its noise: compare only galaxies with a precise r
    assert precise.sum() >= 5
    assert np.allclose(targets[precise, log_mass_formed_index], expected_mass[precise], atol=0.1)
    sigma_z = flat[:, 36]
    z_obs = flat[:, 35]
    spectroscopic = sigma_z == 0
    assert np.allclose(z_obs[spectroscopic], repeated_models["redshift"].to_numpy()[spectroscopic])
    photometric = ~spectroscopic
    assert np.allclose(sigma_z[photometric], sm.SIGMA_Z * (1 + repeated_models["redshift"].to_numpy()[photometric]))


@pytest.mark.parametrize("sfh_model", SFH_MODEL_NAMES)
def test_validation_observations_modes(sfh_model):
    path = sfh_models.model_file(VALIDATION_FILE, sfh_model)
    if not parquet_has_column(path, SET_COLUMNS[sfh_model]):
        pytest.skip(f"{path} not yet generated on the PC")
    held_out = pd.read_parquet(path).iloc[:200]
    spectroscopic = ml_data.validation_observations(held_out, "spectroscopic", np.random.default_rng(32), sfh_model)
    photometric = ml_data.validation_observations(held_out, "photometric", np.random.default_rng(32), sfh_model)
    assert np.all(spectroscopic["flat"][:, 36] == 0)
    assert np.all(photometric["flat"][:, 36] > 0)
    assert np.allclose(spectroscopic["flat"][:, 35], held_out["redshift"])
    log_mass_formed_index = ml_data.target_names(sfh_model).index("log_mass_formed")
    expected = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * held_out["magnitude_offset"].to_numpy()
    assert np.allclose(spectroscopic["targets"][:, log_mass_formed_index], expected)


def test_standardiser_round_trip_and_masked_zero():
    values = np.random.default_rng(33).normal(5, 2, (1000, 3))
    mask = np.ones_like(values, dtype=bool)
    mask[::4, 1] = False
    standardiser = ml_data.Standardiser.fit(values, mask)
    transformed = standardiser.transform(values, mask)
    assert np.allclose(transformed[mask[:, 1], 1].mean(), 0, atol=0.05)
    assert np.all(transformed[~mask[:, 1], 1] == 0)
    restored = ml_data.Standardiser.from_dict(standardiser.to_dict())
    assert np.allclose(restored.inverse_transform(restored.transform(values)), values)
```

- [ ] **Step 2: Run the tests and see them fail.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_data.py -q`
  Expected: the target and decoding tests FAIL with `AttributeError: module 'ml_data' has no attribute 'target_names'`. `test_encode_*`, `test_parallel_decode_*`, `test_simulate_observations_*` and `test_decode_clips_redshift_to_prior` fail with `TypeError: encode_targets() takes 1 positional argument but 2 were given` (or the same error for `simulate_observations`). The colour and standardiser tests pass.

- [ ] **Step 3: Implement `ml_data.py`.**

  Replace lines 1–37 (the module docstring through `N_TARGETS = len(TARGET_NAMES)  # 9`) with:

```python
"""
Network inputs and targets for the SFH posterior estimator (Phase 3).

Model tables (training set, held-out sets, real galaxies) become:
- inputs: 11 colours (band - r) with errors and detection flags, plus r, its error, z_obs and sigma_z;
- targets: the parameters a flow models for one SFH model (sfh_models.MODELS: 9 for delayed-tau, 10 for the
  double power law): the SFH shape, metallicity, dust, nebular emission, masses and redshift; every other
  label is derived from them when decoding.
Every physical definition (prior, labels, noise, brightness) comes from sfh_models and synthetic_model.
"""

import joblib
import numpy as np
import pandas as pd

import sfh_models
import synthetic_model as sm


# =============================================================================
# 1. Constants
# =============================================================================
R_INDEX = sm.BANDS.index("r")
COLOUR_BANDS = []
for band_name in sm.BANDS:
    if band_name != "r":
        COLOUR_BANDS.append(band_name)
N_COLOURS = len(COLOUR_BANDS)                    # 11
N_INPUTS = 3 * N_COLOURS + 4                     # colours, colour errors, flags, then 4 scalars = 37
SCALAR_INPUT_SLICE = slice(3 * N_COLOURS, N_INPUTS)   # [m_r, sigma_r, z_obs, sigma_z]
N_BAND_CHANNELS = 4                              # [m_b - m_r, sigma_b, detected, log10 pivot]

P_SPECTROSCOPIC = 0.5                            # share of simulated galaxies with an exact redshift

# Flow targets that are the log of a drawn label but not among the model's derived labels (each model's
# log_labels give the others, e.g. log_age is the log of age)
EXTRA_LOG_TARGETS = {"log_tau_bc": "tau_bc"}
```

  Replace section 3 (lines 87–152, from the `# 3. Targets and decoding` banner through the end of `decode_samples`) with:

```python
# =============================================================================
# 3. Targets and decoding
# =============================================================================
def target_names(sfh_model):
    """The flow targets of an SFH model, in the order the flow models them.

    Why these: parameters drawn log-uniform in the prior enter as their logs. log_sfr is not a target,
    because a fifth of the delayed-tau prior sits on its floor, a point mass a flow cannot learn. The
    redshift is last, because the decoder needs it first: the oldest allowed age or turnover lookback depends on it.
    Inputs: sfh_model: "delayed" or "dblplaw" (a key of sfh_models.MODELS).
    Returns: list of target names, a new list on every call.
    Raises: KeyError for an unknown model.
    """
    return list(sfh_models.MODELS[sfh_model].target_names)


def encode_targets(table, sfh_model):
    """The flow targets of each row.

    Inputs: table with the model's drawn labels (sfh_models.MODELS[sfh_model].drawn_labels: the redshift,
            the SFH parameters [Gyr or dimensionless], metallicity [Z/Z_sun], av_ism [mag], tau_bc, logU and
            log_mass_formed [log M_sun]) and log_mstar [log M_sun]; sfh_model: "delayed" or "dblplaw".
    Returns: array (n, len(target_names(sfh_model))) in target_names order. The log-uniform drawn
             parameters (the model's log_labels, and tau_bc) enter as log10; the other targets enter
             unchanged.
    """
    log_sources = dict(sfh_models.MODELS[sfh_model].log_labels)
    log_sources.update(EXTRA_LOG_TARGETS)
    columns = []
    for name in target_names(sfh_model):
        if name in log_sources:
            columns.append(np.log10(table[log_sources[name]].to_numpy(dtype=float)))
        else:
            columns.append(table[name].to_numpy(dtype=float))
    return np.column_stack(columns)


def delayed_sfh_parameters(targets, redshift):
    """Delayed-tau SFH parameters of flow samples, inside the prior at each sample's redshift.

    Inputs: targets: dict target name -> array (n,), with log_age and log_tau [log Gyr];
            redshift: array (n,), already inside sfh_models.Z_RANGE.
    Returns: dict with age clipped to [AGE_MIN, maximum_age(redshift)] and tau clipped to TAU_RANGE [Gyr].
             Values inside the prior pass unchanged.
    """
    age = np.clip(10 ** targets["log_age"], sfh_models.AGE_MIN, sfh_models.maximum_age(redshift))
    tau = np.clip(10 ** targets["log_tau"], *sfh_models.TAU_RANGE)
    return {"age": age, "tau": tau}


def dblplaw_sfh_parameters(targets, redshift):
    """Double-power-law SFH parameters of flow samples, inside the prior at each sample's redshift.

    Why tau is recomputed: Bagpipes' tau is the turnover time of the double power law, counted from the Big
    Bang (the SFR maximum only when alpha = beta; in general the maximum is at t = tau (beta/alpha)**(1/(alpha+beta))).
    turnover_lookback is the lookback time to tau, so tau = age of the Universe - turnover_lookback (as in
    sfh_models.dblplaw_components). The flow models turnover_lookback, and tau follows from it.
    Inputs: targets: dict target name -> array (n,), with log_turnover_lookback [log Gyr], log_alpha and log_beta;
            redshift: array (n,), already inside sfh_models.Z_RANGE.
    Returns: dict with
             - turnover_lookback [Gyr], clipped to [TURNOVER_LOOKBACK_MIN, maximum_age(redshift)];
             - tau [Gyr since the Big Bang];
             - alpha clipped to ALPHA_RANGE and beta clipped to BETA_RANGE (dimensionless).
             Values inside the prior pass unchanged.
    """
    turnover_lookback = np.clip(10 ** targets["log_turnover_lookback"], sfh_models.TURNOVER_LOOKBACK_MIN,
                            sfh_models.maximum_age(redshift))
    tau = sfh_models.COSMOLOGY.age(redshift).value - turnover_lookback
    alpha = np.clip(10 ** targets["log_alpha"], *sfh_models.ALPHA_RANGE)
    beta = np.clip(10 ** targets["log_beta"], *sfh_models.BETA_RANGE)
    return {"turnover_lookback": turnover_lookback, "tau": tau, "alpha": alpha, "beta": beta}


# The SFH parameters of each model, decoded from its flow targets
SFH_PARAMETER_DECODERS = {"delayed": delayed_sfh_parameters, "dblplaw": dblplaw_sfh_parameters}


def decode_samples(theta, sfh_model, n_jobs=1):
    """Every training label of flow samples, with the same definitions as the training labels.

    Each sample is decoded at its own redshift. The flow's tails can leave the prior, so:
    - the redshift is clipped to Z_RANGE;
    - the SFH parameters are clipped to their prior at that redshift (SFH_PARAMETER_DECODERS);
    - the metallicity is clipped to METALLICITY_RANGE.
    Values inside the prior pass unchanged, so encoded truths decode exactly.
    Why n_jobs: the double-power-law labels are summed over Bagpipes' age grid, about 0.06 ms per sample, so
    tens of millions of posterior samples take an hour or more in one thread. With n_jobs > 1 the rows are split
    into n_jobs blocks, decoded in threads and concatenated in order. Threads, not processes: numpy releases
    the GIL inside the grid sums, and threads share one copy of Bagpipes, which matters in a 7 GB GPU job.
    Each row is decoded on its own, so the result equals that of n_jobs = 1 to rounding.
    Inputs: theta: array (n, len(target_names(sfh_model))) in target_names order (posterior samples or
            encoded truths); sfh_model: "delayed" or "dblplaw"; n_jobs: threads (1 decodes in this thread).
    Returns: DataFrame with the model's training_labels columns, in that order, one row per row of theta:
             - the drawn labels and their logs;
             - log_mstar [log M_sun];
             - from the model's analytic_labels: mw_age [Gyr ago], f_young, f_mid and f_old,
               log_sfr [log M_sun/yr] and log_ssfr [log yr^-1].
    Raises: ValueError when theta does not have one column per target of the model, or n_jobs < 1.
    """
    names = target_names(sfh_model)
    theta = np.asarray(theta, dtype=float)
    if theta.ndim != 2 or theta.shape[1] != len(names):
        raise ValueError(f"theta of shape {theta.shape} does not hold the {len(names)} targets of the "
                         f"{sfh_model!r} model")
    if n_jobs < 1:
        raise ValueError(f"n_jobs = {n_jobs}: decode with 1 or more threads")
    if n_jobs == 1 or len(theta) < 2:
        return decode_rows(theta, sfh_model)

    tasks = []
    for block in np.array_split(theta, min(n_jobs, len(theta))):
        tasks.append(joblib.delayed(decode_rows)(block, sfh_model))
    decoded_blocks = joblib.Parallel(n_jobs=n_jobs, prefer="threads")(tasks)
    return pd.concat(decoded_blocks, ignore_index=True)


def decode_rows(theta, sfh_model):
    """The decode of decode_samples, in this process and without its checks.

    Inputs: theta: array (n, len(target_names(sfh_model))) in target_names order; sfh_model: "delayed" or
            "dblplaw".
    Returns: DataFrame with the model's training_labels columns (see decode_samples).
    """
    model = sfh_models.MODELS[sfh_model]
    names = target_names(sfh_model)
    targets = {}
    for target_index, name in enumerate(names):
        targets[name] = theta[:, target_index]
    redshift = np.clip(targets["redshift"], *sfh_models.Z_RANGE)

    # 1. The drawn labels, inside the prior at each sample's redshift
    drawn = {"redshift": redshift}
    drawn.update(SFH_PARAMETER_DECODERS[sfh_model](targets, redshift))
    drawn["metallicity"] = np.clip(10 ** targets["log_metallicity"], *sfh_models.METALLICITY_RANGE)
    drawn["av_ism"] = targets["av_ism"]
    drawn["tau_bc"] = 10 ** targets["log_tau_bc"]
    drawn["logU"] = targets["logU"]
    drawn["log_mass_formed"] = targets["log_mass_formed"]

    # 2. Labels derived from them, as sm.summary_labels defines them
    labels = dict(drawn)
    for log_label, drawn_label in model.log_labels.items():
        labels[log_label] = np.log10(drawn[drawn_label])
    labels["log_mstar"] = targets["log_mstar"]
    labels.update(model.analytic_labels(drawn, targets["log_mstar"]))
    return pd.DataFrame(labels, columns=model.training_labels)


# delayed-tau aliases for ml_networks.DEFAULT_CONFIG, removed in Task 8
TARGET_NAMES = target_names("delayed")
N_TARGETS = len(TARGET_NAMES)
```

  In section 4, replace `simulate_observations` (lines 170–197) with:

```python
def simulate_observations(table, catalogue_r, catalogue_z, noise_table, rng, sfh_model):
    """One fresh noisy realisation of every noise-free model: brightness, S-PLUS noise and redshift.

    Inputs: table: noise-free models of sfh_model (its training-set columns);
            catalogue_r [mag], catalogue_z: DR3 r_auto and zml (the brightness donors);
            noise_table: splus_noise_model.csv; rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: dict(flat, banded, valid, targets):
             - flat, banded and valid: the inputs as colour_inputs gives them;
             - targets: array (n, len(target_names(sfh_model))), with the mass-like labels rescaled to the
               drawn brightness.
             valid is False where r is undetected or no DR3 galaxy lies at the model's redshift.
    """
    true_redshift = table["redshift"].to_numpy(dtype=float)
    noiseless = table[sm.BANDS].to_numpy(dtype=float)

    # 1. Realistic brightness, and the masses that go with it
    scaled, offsets = sm.scale_to_observed_r(noiseless, true_redshift, catalogue_r, catalogue_z, rng)
    labels = table.copy()
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED     # the noise-free models' value, rescaled below
    labels = sm.rescale_mass_labels(labels, offsets)

    # 2. S-PLUS noise and non-detections
    noisy, errors, _ = sm.add_splus_noise(scaled, noise_table, rng)

    # 3. Observed redshift
    z_obs, sigma_z = draw_observed_redshift(true_redshift, rng)

    observations = colour_inputs(noisy, errors, z_obs, sigma_z)
    observations["valid"] = observations["valid"] & np.isfinite(offsets)
    observations["targets"] = encode_targets(labels, sfh_model)
    return observations
```

  Replace `validation_observations` (lines 200–229) with:

```python
def validation_observations(held_out, mode, rng, sfh_model):
    """Inputs and targets of a held-out set whose noise is already fixed (synthetic_validation/test).

    Inputs: held_out: DataFrame of sfh_model with noisy BANDS magnitudes, e_<band> errors, magnitude_offset
            and labels (mass labels already rescaled);
            mode: "spectroscopic" (sigma_z = 0) or "photometric" (sigma_z = sm.SIGMA_Z * (1 + z), with z_obs
            drawn from rng);
            rng: numpy Generator (seeded by the caller); sfh_model: "delayed" or "dblplaw".
    Returns: dict(flat, banded, valid, targets), with targets (n, len(target_names(sfh_model))).
    """
    true_redshift = held_out["redshift"].to_numpy(dtype=float)
    if mode == "spectroscopic":
        sigma_z = np.zeros(len(true_redshift))
        z_obs = true_redshift.copy()
    elif mode == "photometric":
        sigma_z = sm.SIGMA_Z * (1 + true_redshift)
        z_obs = true_redshift + sigma_z * rng.normal(0, 1, len(true_redshift))
    else:
        raise ValueError(f"unknown mode {mode!r}: use 'spectroscopic' or 'photometric'")

    error_columns = []
    for band in sm.BANDS:
        error_columns.append(f"e_{band}")
    magnitudes = held_out[sm.BANDS].to_numpy(dtype=float)
    errors = held_out[error_columns].to_numpy(dtype=float)

    labels = held_out.copy()
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * held_out["magnitude_offset"]
    observations = colour_inputs(magnitudes, errors, z_obs, sigma_z)
    observations["valid"] = observations["valid"] & np.isfinite(held_out["magnitude_offset"].to_numpy())
    observations["targets"] = encode_targets(labels, sfh_model)
    return observations
```

  In `Synthetic_domain/train_ml_model.py`, make three edits. Task 8 replaces `"delayed"` with the `--sfh-model` argument.
  - Lines 95–96: replace
    ```python
        first_epoch = to_tensors(ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table,
                                                               np.random.default_rng(arguments.seed)))
    ```
    with
    ```python
        first_epoch = to_tensors(ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table,
                                                               np.random.default_rng(arguments.seed), "delayed"))
    ```
  - Line 108: replace `observations = ml_data.validation_observations(validation, mode, np.random.default_rng(seed))` with `observations = ml_data.validation_observations(validation, mode, np.random.default_rng(seed), "delayed")`.
  - Line 114: replace `observations = ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table, rng)` with `observations = ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table, rng, "delayed")`.

- [ ] **Step 4: Run the tests and see them pass.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_data.py -q`
  Expected: all pass. `test_validation_observations_modes[dblplaw]` SKIPS ("synthetic_validation_dblplaw.parquet not yet generated on the PC"). `test_validation_observations_modes[delayed]` passes, because the stored delayed-τ set has the `age` column. The simulated fixture builds 400 Bagpipes models, which takes under a minute.
  `test_parallel_decode_equals_the_serial_decode` decodes in joblib threads; it adds a few seconds.
  If `test_decode_extreme_dblplaw_samples_stay_finite` fails, `sfh_models.dblplaw_labels` (Task 2) produced NaN at a prior corner. Fix it there, not here. Task 2 tests the likeliest corner, α = β = 1000 with τ = 0.1 Gyr: Bagpipes' grid step near a lookback of 13 Gyr is about 0.08 Gyr, so the sharp peak can fall between grid points.

- [ ] **Step 5: Run the affected test files.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_data.py tests/test_ml_networks.py tests/test_train_ml_model.py -q`
  Expected: all pass. `test_training_cli_writes_a_run` runs on the stored delayed-τ sets and passes through the `"delayed"` calls; it takes a few minutes on the CPU.
  `ml_model.ipynb` still calls the old signatures. Task 10 rewrites it, and nothing executes it before then.

- [ ] **Step 6: Commit.**
```bash
git add Synthetic_domain/ml_data.py Synthetic_domain/train_ml_model.py Synthetic_domain/tests/test_ml_data.py
git commit -m "feat: ml_data encodes and decodes the flow targets of either SFH model

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 8: `train_ml_model.py --sfh-model` and the checkpoint guard

**Files:**
- Modify: `Synthetic_domain/ml_networks.py`:
  - section 1, `DEFAULT_CONFIG` (lines 18–42): add `model_config` and `LEGACY_SFH_MODEL`;
  - section 4, `save_run` and `load_run` (lines 214–269).
- Modify: `Synthetic_domain/train_ml_model.py` (full rewrite below).
- Modify: `Synthetic_domain/ml_data.py`: delete the two alias lines added in Task 7, at the end of section 3.
- Modify: `Synthetic_domain/remote/submit_training.sh` (on top of Task 1's version, which already ships `sfh_models.py`: this task adds `SFH_MODEL` and the model's set files only), `Synthetic_domain/remote/train_job.sh` (78 lines: header line 7, `EMBEDDINGS` doc line 10, `EMBEDDINGS=` line 15, training loop lines 58–70) and `Synthetic_domain/remote/README.md` (lines 87–135 after Task 6 (78–126 before it; find the block by its first line `# Training the ML posterior estimator (GPU)` and its last line `        --max-train-rows 3000 --max-validation-rows 500 --device cpu`), the training section).
- Test:
  - `Synthetic_domain/tests/test_ml_networks.py` (full rewrite below);
  - `Synthetic_domain/tests/test_train_ml_model.py` (full rewrite below);
  - `Synthetic_domain/tests/test_ml_data.py` (one test appended).

**Interfaces:**
- **Consumes:** `ml_data.target_names(sfh_model)` (Task 7), `sfh_models.MODELS` and `sfh_models.model_file(name, sfh_model)` (Task 1).
- **Produces:**
  - `ml_networks.DEFAULT_CONFIG`, which no longer holds `n_targets` or `target_names`;
  - `ml_networks.LEGACY_SFH_MODEL = "delayed"`;
  - `ml_networks.model_config(sfh_model) -> dict` (`DEFAULT_CONFIG` + `sfh_model`, `target_names`, `n_targets`);
  - `ml_networks.save_run(directory, estimator, config, history, standardisers)`. Its signature is unchanged. It raises `ValueError` when the config has no `sfh_model` or `target_names`.
  - `ml_networks.load_run(directory, sfh_model, device="cpu", target_names=None) -> (estimator, config, standardisers)`:
    - `config["sfh_model"]` is always set, to `"delayed"` for a run saved without one;
    - it raises `ValueError("run ... was trained for the ... SFH model ...")` on a model mismatch;
    - it raises `ValueError("run was trained for targets ...")` on a target mismatch.
  - `train_ml_model.DEFAULT_SFH_MODEL = "delayed"`;
  - `train_ml_model.data_file(data_dir, name, sfh_model) -> str`;
  - `train_ml_model.run_directory(output_dir, embedding, sfh_model) -> str`;
  - `train_ml_model.parse_arguments(argv=None)` and `train_ml_model.main(argv=None)`;
  - the CLI flag `--sfh-model {dblplaw,delayed}` (default `delayed`);
  - the remote setting `SFH_MODEL` (default `delayed`) in `submit_training.sh` and `train_job.sh`.
- **Removes:** `ml_data.TARGET_NAMES` and `ml_data.N_TARGETS`.

- [ ] **Step 1: Write the failing tests.**

  Replace the whole of `Synthetic_domain/tests/test_ml_networks.py` with:

```python
"""Tests for ml_networks.py (CPU, small). Run from Synthetic_domain/."""

import json
import os

import numpy as np
import pytest
import torch

import ml_data
import ml_networks

NOISE_SIGMA = 0.5        # toy problem: x = theta + N(0, NOISE_SIGMA), theta ~ N(0, 1)
TOY_TARGET_NAMES = ["toy_target_0", "toy_target_1"]
STORED_DELAYED_RUN = "ml_runs/mlp"   # trained before runs recorded their SFH model


def toy_problem(n_rows, generator):
    """Gaussian toy: two targets observed with noise in the first inputs (flat) and band channel 0."""
    theta = torch.randn(n_rows, 2, generator=generator)
    observed = theta + NOISE_SIGMA * torch.randn(n_rows, 2, generator=generator)
    flat = torch.zeros(n_rows, ml_data.N_INPUTS)
    flat[:, :2] = observed
    banded = torch.zeros(n_rows, 12, ml_data.N_BAND_CHANNELS)
    banded[:, :2, 0] = observed
    return {"flat": flat, "banded": banded, "targets": theta}


def small_config(embedding, sfh_model="delayed"):
    """A small, fast configuration for tests: the toy problem's two targets, recorded under sfh_model."""
    config = ml_networks.model_config(sfh_model)
    config.update(embedding=embedding, n_targets=2, target_names=TOY_TARGET_NAMES, summary_size=16,
                  hidden_size=64, flow_transforms=3, flow_hidden=(64, 64), epochs=40, patience=40,
                  batch_size=512, learning_rate=3e-3, seed=0)
    return config


def unit_standardisers(n_targets):
    """Standardisers that leave inputs and n_targets targets unchanged."""
    return {
        "flat": ml_data.Standardiser(np.zeros(ml_data.N_INPUTS), np.ones(ml_data.N_INPUTS)),
        "banded": ml_data.Standardiser(np.zeros(ml_data.N_BAND_CHANNELS), np.ones(ml_data.N_BAND_CHANNELS)),
        "targets": ml_data.Standardiser(np.zeros(n_targets), np.ones(n_targets)),
    }


def save_toy_run(directory, sfh_model="delayed"):
    """Write an untrained toy run (2 targets named TOY_TARGET_NAMES, recorded under sfh_model) with save_run."""
    config = small_config("mlp", sfh_model)
    history = [{"epoch": 0, "train_loss": 1.0, "validation_loss": 1.1}]
    ml_networks.save_run(str(directory), ml_networks.build_estimator(config), config, history,
                         unit_standardisers(2))


@pytest.mark.parametrize("embedding", ["mlp", "cnn"])
def test_flow_recovers_gaussian_posterior(embedding):
    generator = torch.Generator().manual_seed(0)
    config = small_config(embedding)
    estimator = ml_networks.build_estimator(config)
    validation = toy_problem(2000, generator)

    def make_epoch(epoch_index):
        return toy_problem(8192, generator)

    ml_networks.train_estimator(estimator, make_epoch, [validation], config, "cpu", log=lambda message: None)

    observed_value = 1.0
    flat = torch.zeros(1, ml_data.N_INPUTS)
    flat[0, :2] = observed_value
    banded = torch.zeros(1, 12, ml_data.N_BAND_CHANNELS)
    banded[0, :2, 0] = observed_value
    torch.manual_seed(1)
    samples = estimator.sample(flat, banded, 4000)[0].detach().numpy()
    posterior_mean = observed_value / (1 + NOISE_SIGMA**2)
    posterior_std = np.sqrt(NOISE_SIGMA**2 / (1 + NOISE_SIGMA**2))
    assert np.allclose(samples.mean(axis=0), posterior_mean, atol=0.12)
    assert np.allclose(samples.std(axis=0), posterior_std, rtol=0.2)


def test_estimator_samples_finite_for_fully_masked_input():
    estimator = ml_networks.build_estimator(small_config("cnn"))
    flat = torch.zeros(3, ml_data.N_INPUTS)
    banded = torch.zeros(3, 12, ml_data.N_BAND_CHANNELS)
    samples = estimator.sample(flat, banded, 10)
    assert samples.shape == (3, 10, 2)
    assert torch.isfinite(samples).all()


@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_model_config_records_the_model_and_its_targets(sfh_model):
    config = ml_networks.model_config(sfh_model)
    assert config["sfh_model"] == sfh_model
    assert config["target_names"] == ml_data.target_names(sfh_model)
    assert config["n_targets"] == len(ml_data.target_names(sfh_model))
    assert "sfh_model" not in ml_networks.DEFAULT_CONFIG        # the defaults belong to no model
    assert "target_names" not in ml_networks.DEFAULT_CONFIG


def test_dblplaw_estimator_samples_its_ten_targets():
    config = ml_networks.model_config("dblplaw")
    config.update(summary_size=16, hidden_size=64, flow_transforms=2, flow_hidden=(32, 32))
    estimator = ml_networks.build_estimator(config)
    samples = estimator.sample(torch.zeros(2, ml_data.N_INPUTS), torch.zeros(2, 12, ml_data.N_BAND_CHANNELS), 5)
    assert samples.shape == (2, 5, 10)


def test_save_load_round_trip_on_cpu(tmp_path):
    config = small_config("mlp")
    estimator = ml_networks.build_estimator(config)
    standardisers = unit_standardisers(2)
    standardisers["targets"] = ml_data.Standardiser(np.array([1.0, 2.0]), np.array([0.5, 3.0]))
    history = [{"epoch": 0, "train_loss": 1.0, "validation_loss": 1.1}]
    ml_networks.save_run(str(tmp_path), estimator, config, history, standardisers)
    loaded, loaded_config, loaded_standardisers = ml_networks.load_run(str(tmp_path), "delayed", device="cpu",
                                                                       target_names=TOY_TARGET_NAMES)

    flat = np.random.default_rng(0).normal(size=(4, ml_data.N_INPUTS))
    banded = np.random.default_rng(1).normal(size=(4, 12, ml_data.N_BAND_CHANNELS))
    before = ml_networks.posterior_samples(estimator, standardisers, flat, banded, 50, seed=3)
    after = ml_networks.posterior_samples(loaded, loaded_standardisers, flat, banded, 50, seed=3)
    assert loaded_config["embedding"] == "mlp"
    assert loaded_config["sfh_model"] == "delayed"
    assert loaded_config["target_names"] == TOY_TARGET_NAMES
    assert before.shape == (4, 50, 2)
    assert np.allclose(before, after)


def test_load_run_of_a_checkpoint_without_target_names_raises(tmp_path):
    save_toy_run(tmp_path)
    config_path = tmp_path / ml_networks.CONFIG_FILE
    stored = json.loads(config_path.read_text())
    del stored["target_names"]                       # a checkpoint written before target names were recorded
    config_path.write_text(json.dumps(stored))
    with pytest.raises(ValueError, match="run was trained for targets"):
        ml_networks.load_run(str(tmp_path), "delayed", target_names=TOY_TARGET_NAMES)


def test_load_run_of_other_targets_raises(tmp_path):
    save_toy_run(tmp_path)
    with pytest.raises(ValueError, match="run was trained for targets"):
        ml_networks.load_run(str(tmp_path), "delayed")          # expects the delayed-tau targets


def test_load_run_of_another_sfh_model_raises(tmp_path):
    save_toy_run(tmp_path, sfh_model="dblplaw")
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(str(tmp_path), "delayed", target_names=TOY_TARGET_NAMES)


def test_load_run_without_sfh_model_counts_as_delayed(tmp_path):
    save_toy_run(tmp_path)
    config_path = tmp_path / ml_networks.CONFIG_FILE
    stored = json.loads(config_path.read_text())
    del stored["sfh_model"]                          # a run saved before runs recorded their model
    config_path.write_text(json.dumps(stored))
    _, config, _ = ml_networks.load_run(str(tmp_path), "delayed", target_names=TOY_TARGET_NAMES)
    assert config["sfh_model"] == "delayed"
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(str(tmp_path), "dblplaw", target_names=TOY_TARGET_NAMES)


@pytest.mark.skipif(not os.path.exists(os.path.join(STORED_DELAYED_RUN, ml_networks.MODEL_FILE)),
                    reason="no ml_runs/mlp weights here (model.pt is git-ignored)")
def test_the_stored_delayed_run_loads_only_as_delayed():
    _, config, _ = ml_networks.load_run(STORED_DELAYED_RUN, "delayed")
    assert config["sfh_model"] == "delayed"
    assert config["target_names"] == ml_data.target_names("delayed")
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(STORED_DELAYED_RUN, "dblplaw")


def test_save_run_rejects_a_config_without_sfh_model(tmp_path):
    config = small_config("mlp")
    del config["sfh_model"]
    with pytest.raises(ValueError, match="sfh_model"):
        ml_networks.save_run(str(tmp_path), ml_networks.build_estimator(small_config("mlp")), config, [], {})


def test_save_run_rejects_target_names_that_do_not_match_n_targets(tmp_path):
    config = small_config("mlp")
    config["target_names"] = ml_data.target_names("delayed")   # 9 names for a 2-target flow
    with pytest.raises(ValueError, match="target_names"):
        ml_networks.save_run(str(tmp_path), ml_networks.build_estimator(small_config("mlp")), config, [], {})
```

  Replace the whole of `Synthetic_domain/tests/test_train_ml_model.py` with:

```python
"""Tests of the training CLI: paths per SFH model, and a smoke run on a tiny subset (CPU). Run from Synthetic_domain/."""

import json
import os
import subprocess
import sys

import pandas as pd
import pyarrow.parquet as pq
import pytest

import ml_data
import ml_networks
import sfh_models
import train_ml_model

DATA_FILES = ("training_set.parquet", "synthetic_validation.parquet")
SET_COLUMNS = {"delayed": "age", "dblplaw": "turnover_lookback"}   # a column only that model's sets hold
OTHER_MODEL = {"delayed": "dblplaw", "dblplaw": "delayed"}


def sets_are_generated(sfh_model):
    """True when both sets of sfh_model exist and hold its labels (the schema has SET_COLUMNS[sfh_model])."""
    for name in DATA_FILES:
        path = sfh_models.model_file(name, sfh_model)
        try:
            if SET_COLUMNS[sfh_model] not in pq.read_schema(path).names:
                return False
        except OSError:
            return False
    return True


def test_data_files_and_run_directories_follow_the_model_suffix():
    assert train_ml_model.data_file("data", "training_set.parquet", "delayed") == os.path.join(
        "data", "training_set.parquet")
    assert train_ml_model.data_file("data", "training_set.parquet", "dblplaw") == os.path.join(
        "data", "training_set_dblplaw.parquet")
    assert train_ml_model.run_directory("ml_runs", "mlp", "delayed") == os.path.join("ml_runs", "mlp")
    assert train_ml_model.run_directory("ml_runs", "mlp", "dblplaw") == os.path.join("ml_runs", "mlp_dblplaw")


def test_sfh_model_argument_defaults_to_delayed():
    assert train_ml_model.parse_arguments([]).sfh_model == "delayed"
    assert train_ml_model.parse_arguments(["--sfh-model", "dblplaw"]).sfh_model == "dblplaw"
    with pytest.raises(SystemExit):
        train_ml_model.parse_arguments(["--sfh-model", "exponential"])


@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_training_cli_writes_a_run(tmp_path, sfh_model):
    if not sets_are_generated(sfh_model):
        pytest.skip(f"the {sfh_model} sets are not yet generated on the PC")
    command = [sys.executable, "train_ml_model.py", "--sfh-model", sfh_model, "--embedding", "mlp",
               "--output-dir", str(tmp_path), "--epochs", "2", "--max-train-rows", "3000",
               "--max-validation-rows", "500", "--device", "cpu"]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=900)
    assert completed.returncode == 0, completed.stderr[-2000:]
    run_directory = train_ml_model.run_directory(str(tmp_path), "mlp", sfh_model)
    for name in ("model.pt", "config.json", "history.csv"):
        assert os.path.exists(os.path.join(run_directory, name)), name
    history = pd.read_csv(os.path.join(run_directory, "history.csv"))
    assert len(history) == 2
    with open(os.path.join(run_directory, "config.json")) as source:
        config = json.load(source)
    assert config["sfh_model"] == sfh_model
    assert config["target_names"] == ml_data.target_names(sfh_model)
    assert config["n_targets"] == len(ml_data.target_names(sfh_model))
    ml_networks.load_run(run_directory, sfh_model)
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(run_directory, OTHER_MODEL[sfh_model])
```

  Append to `Synthetic_domain/tests/test_ml_data.py`:

```python
def test_no_global_target_list_remains():
    """The targets depend on the SFH model: callers ask target_names(sfh_model)."""
    assert not hasattr(ml_data, "TARGET_NAMES")
    assert not hasattr(ml_data, "N_TARGETS")
```

- [ ] **Step 2: Run the tests and see them fail.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_networks.py tests/test_train_ml_model.py tests/test_ml_data.py -q`
  Expected failures:
  - `AttributeError: module 'ml_networks' has no attribute 'model_config'`;
  - `AttributeError: module 'train_ml_model' has no attribute 'data_file'` and `... 'parse_arguments'`;
  - `test_training_cli_writes_a_run[delayed]`: the script has no `--sfh-model` option ("unrecognized arguments");
  - `test_no_global_target_list_remains`, which fails its assertion.

- [ ] **Step 3: Implement.**

  **`ml_networks.py`, section 1.** Replace lines 21–39, the `DEFAULT_CONFIG` dict, with:

```python
DEFAULT_CONFIG = {
    "embedding": "mlp",              # "mlp" or "cnn"; only the MLP is trained (the CNN lost in 16-bin runs)
    "summary_size": 64,              # length of the embedding output
    "hidden_size": 256,              # width of the dense layers
    "mlp_layers": 3,
    "cnn_channels": 64,
    "cnn_kernel": 3,
    "cnn_layers": 3,
    "head_layers": 2,                # dense layers after the CNN pooling
    "flow_transforms": 8,
    "flow_hidden": (256, 256),
    "learning_rate": 1e-3,
    "batch_size": 4096,
    "epochs": 200,
    "patience": 15,                  # epochs without validation improvement before stopping
    "seed": 0,
}
LEGACY_SFH_MODEL = "delayed"         # SFH model of a run saved before runs recorded one (ml_runs/mlp)
```

  Insert after the `HISTORY_FILE = "history.csv"` line (line 42):

```python


def model_config(sfh_model):
    """DEFAULT_CONFIG for the flow of one SFH model, with its sfh_model, target_names and n_targets.

    Why stored: config.json is all a later reader has, so a run must say which SFH model and targets it
    learnt. load_run checks both.
    Inputs: sfh_model: "delayed" or "dblplaw".
    Returns: a new dict.
    """
    config = dict(DEFAULT_CONFIG)
    names = ml_data.target_names(sfh_model)
    config["sfh_model"] = sfh_model
    config["target_names"] = names
    config["n_targets"] = len(names)
    return config
```

  **`ml_networks.py`, section 4.** Replace `save_run` and `load_run` (lines 214–269) with:

```python
def save_run(directory, estimator, config, history, standardisers):
    """Write model.pt (weights), config.json (config, SFH model, target names, standardisers) and history.csv.

    Inputs: directory: output folder (created if missing); estimator: PosteriorEstimator;
            config: its configuration (model_config output, updated as needed). sfh_model names the SFH
            model, and target_names names the n_targets flow dimensions in order;
            history: list of per-epoch dicts; standardisers: {name: ml_data.Standardiser}.
    Raises: ValueError when config has no sfh_model or target_names, or when target_names does not have
            n_targets entries.
    """
    for key in ("sfh_model", "target_names"):
        if key not in config:
            raise ValueError(f"config has no {key!r}: build it with model_config(sfh_model)")
    target_names = list(config["target_names"])
    if len(target_names) != config["n_targets"]:
        raise ValueError(f"target_names has {len(target_names)} entries but the flow has "
                         f"n_targets = {config['n_targets']}")
    os.makedirs(directory, exist_ok=True)
    torch.save(estimator.state_dict(), os.path.join(directory, MODEL_FILE))
    stored = dict(config)
    stored["target_names"] = target_names
    stored["flow_hidden"] = list(config["flow_hidden"])
    stored["standardisers"] = {}
    for name, standardiser in standardisers.items():
        stored["standardisers"][name] = standardiser.to_dict()
    with open(os.path.join(directory, CONFIG_FILE), "w") as output:
        json.dump(stored, output, indent=1)
    with open(os.path.join(directory, HISTORY_FILE), "w") as output:
        output.write("epoch,train_loss,validation_loss,learning_rate\n")
        for row in history:
            output.write(f"{row['epoch']},{row['train_loss']},{row['validation_loss']},{row.get('learning_rate', '')}\n")


def load_run(directory, sfh_model, device="cpu", target_names=None):
    """Read a run written by save_run, refusing a run of another SFH model or of other targets.

    Why the checks: a run of the other SFH model, or an old checkpoint (e.g. the 19 targets of the 16-bin
    model), would otherwise load. Its samples would then be decoded as the wrong targets without any error.
    Inputs: directory: run folder; sfh_model: the SFH model the caller decodes ("delayed" or "dblplaw");
            device: torch device; target_names: the targets the caller decodes, by default
            ml_data.target_names(sfh_model).
    Returns: (estimator on device, config, standardisers dict). config["sfh_model"] is always set:
             LEGACY_SFH_MODEL for a run saved before models were recorded.
    Raises: ValueError("run ... was trained for the ... SFH model ...") when the stored sfh_model (or
            LEGACY_SFH_MODEL when there is none) differs from sfh_model.
            ValueError("run was trained for targets ...") when the stored target_names (None for a run
            saved before they were recorded) differ from target_names.
    """
    if target_names is None:
        target_names = ml_data.target_names(sfh_model)
    with open(os.path.join(directory, CONFIG_FILE)) as source:
        stored = json.load(source)
    stored_model = stored.get("sfh_model", LEGACY_SFH_MODEL)
    if stored_model != sfh_model:
        raise ValueError(f"run {directory} was trained for the {stored_model!r} SFH model, not {sfh_model!r}: "
                         "load the run of that model (sfh_models.model_file)")
    if stored.get("target_names") != list(target_names):
        raise ValueError(f"run was trained for targets {stored.get('target_names')}, "
                         f"not {list(target_names)}: retrain it, or move {directory} to ml_runs/superseded/")
    stored["sfh_model"] = stored_model
    standardisers = {}
    for name, statistics in stored.pop("standardisers").items():
        standardisers[name] = ml_data.Standardiser.from_dict(statistics)
    estimator = PosteriorEstimator(stored)
    weights = torch.load(os.path.join(directory, MODEL_FILE), map_location=device, weights_only=True)
    estimator.load_state_dict(weights)
    estimator.to(device)
    estimator.eval()
    return estimator, stored, standardisers
```

  **`ml_data.py`.** Delete these three lines at the end of section 3:

```python
# delayed-tau aliases for ml_networks.DEFAULT_CONFIG, removed in Task 8
TARGET_NAMES = target_names("delayed")
N_TARGETS = len(TARGET_NAMES)
```

  **`train_ml_model.py`.** Replace the whole file with:

```python
"""
Train the SFH posterior estimator of one SFH model.

Usage (from Synthetic_domain/):
    python train_ml_model.py                         # the MLP of the delayed-tau model
    python train_ml_model.py --sfh-model dblplaw     # the MLP of the double-power-law model
    python train_ml_model.py --embedding cnn         # kept for comparison; not part of the pipeline
Each run reads its model's sets. sfh_models.model_file names them: training_set.parquet and
synthetic_validation.parquet for delayed-tau, training_set_dblplaw.parquet and
synthetic_validation_dblplaw.parquet for the double power law.
Each run writes ml_runs/<embedding>/ (delayed-tau) or ml_runs/<embedding>_dblplaw/, with model.pt,
config.json and history.csv. config.json records the sfh_model and the ml_data.target_names the flow was
trained for.
The training set is noise-free. Every epoch draws a fresh brightness, S-PLUS noise and observed redshift for
every model (ml_data).
"""

import argparse
import os
import subprocess
import time

import numpy as np
import pandas as pd
import torch

import ml_data
import ml_networks
import sfh_models

TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
NOISE_MODEL_FILE = "splus_noise_model.csv"
CATALOGUE_FILE = "dr3_corrected.csv"
MIN_PROB_GALAXY = 0.8
VALIDATION_SEEDS = {"spectroscopic": 100, "photometric": 101}
DEFAULT_EMBEDDING = "mlp"        # the CNN lost the 16-bin comparison, so only the MLP is trained
DEFAULT_SFH_MODEL = "delayed"    # the model of the existing run, ml_runs/mlp


def to_tensors(observations):
    """Keep valid rows and convert flat, banded, targets to float32 tensors (not yet standardised)."""
    valid = observations["valid"]
    tensors = {}
    for name in ("flat", "banded", "targets"):
        tensors[name] = torch.as_tensor(observations[name][valid], dtype=torch.float32)
    return tensors


def standardise(tensors, standardisers):
    """Standardised copy of a tensor dict (masked colours become 0)."""
    flat = tensors["flat"].numpy()
    banded = tensors["banded"].numpy()
    flat_mask, banded_mask = ml_data.input_masks(flat, banded)
    return {
        "flat": torch.as_tensor(standardisers["flat"].transform(flat, flat_mask), dtype=torch.float32),
        "banded": torch.as_tensor(standardisers["banded"].transform(banded, banded_mask), dtype=torch.float32),
        "targets": torch.as_tensor(standardisers["targets"].transform(tensors["targets"].numpy()),
                                   dtype=torch.float32),
    }


def git_commit():
    """Current git commit hash, or 'unknown' outside a repository (e.g. inside a job)."""
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def data_file(data_dir, name, sfh_model):
    """Path of one of an SFH model's data sets: name with the model's suffix (sfh_models.model_file), in data_dir."""
    return os.path.join(data_dir, sfh_models.model_file(name, sfh_model))


def run_directory(output_dir, embedding, sfh_model):
    """Folder of a run: <output_dir>/<embedding> for delayed-tau (as ml_runs/mlp), with the model's suffix otherwise."""
    return sfh_models.model_file(os.path.join(output_dir, embedding), sfh_model)


def parse_arguments(argv=None):
    """Command-line arguments (argv: list of strings, default sys.argv[1:])."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sfh-model", choices=sorted(sfh_models.MODELS), default=DEFAULT_SFH_MODEL)
    parser.add_argument("--embedding", choices=sorted(ml_networks.EMBEDDINGS), default=DEFAULT_EMBEDDING)
    parser.add_argument("--data-dir", default=".")
    parser.add_argument("--output-dir", default="ml_runs")
    parser.add_argument("--epochs", type=int, default=ml_networks.DEFAULT_CONFIG["epochs"])
    parser.add_argument("--max-train-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--max-validation-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--seed", type=int, default=ml_networks.DEFAULT_CONFIG["seed"])
    return parser.parse_args(argv)


def main(argv=None):
    """Parse arguments, build the data pipeline of one SFH model, train one embedding and save the run."""
    arguments = parse_arguments(argv)
    sfh_model = arguments.sfh_model

    device = arguments.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}, SFH model: {sfh_model}", flush=True)

    # 1. Data
    training = pd.read_parquet(data_file(arguments.data_dir, TRAINING_FILE, sfh_model))
    if arguments.max_train_rows is not None:
        training = training.iloc[:arguments.max_train_rows]
    validation = pd.read_parquet(data_file(arguments.data_dir, VALIDATION_FILE, sfh_model))
    if arguments.max_validation_rows is not None:
        validation = validation.iloc[:arguments.max_validation_rows]
    noise_table = pd.read_csv(os.path.join(arguments.data_dir, NOISE_MODEL_FILE))
    catalogue = pd.read_csv(os.path.join(arguments.data_dir, CATALOGUE_FILE), usecols=["PROB_GAL", "r_auto", "zml"])
    catalogue = catalogue[catalogue["PROB_GAL"] > MIN_PROB_GALAXY]
    catalogue_r = catalogue["r_auto"].to_numpy()
    catalogue_z = catalogue["zml"].to_numpy()

    # 2. Standardisation statistics from one simulated epoch (seeded)
    first_epoch = to_tensors(ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table,
                                                           np.random.default_rng(arguments.seed), sfh_model))
    flat_mask, banded_mask = ml_data.input_masks(first_epoch["flat"].numpy(), first_epoch["banded"].numpy())
    standardisers = {
        "flat": ml_data.Standardiser.fit(first_epoch["flat"].numpy(), flat_mask),
        "banded": ml_data.Standardiser.fit(first_epoch["banded"].numpy().reshape(-1, ml_data.N_BAND_CHANNELS),
                                           banded_mask.reshape(-1, ml_data.N_BAND_CHANNELS)),
        "targets": ml_data.Standardiser.fit(first_epoch["targets"].numpy()),
    }

    # 3. Fixed validation copies: spectroscopic and photometric redshifts
    validation_sets = []
    for mode, seed in VALIDATION_SEEDS.items():
        observations = ml_data.validation_observations(validation, mode, np.random.default_rng(seed), sfh_model)
        validation_sets.append(standardise(to_tensors(observations), standardisers))

    def make_epoch(epoch_index):
        """Fresh noisy realisation for this epoch (seeded by run seed and epoch)."""
        rng = np.random.default_rng([arguments.seed, epoch_index])
        observations = ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table, rng,
                                                     sfh_model)
        return standardise(to_tensors(observations), standardisers)

    # 4. Train and save
    config = ml_networks.model_config(sfh_model)
    config.update(embedding=arguments.embedding, epochs=arguments.epochs, seed=arguments.seed,
                  training_commit=git_commit(), torch_version=torch.__version__, n_training_models=len(training),
                  device=device)
    estimator = ml_networks.build_estimator(config)
    start = time.time()
    history = ml_networks.train_estimator(estimator, make_epoch, validation_sets, config, device,
                                          log=lambda message: print(message, flush=True))
    config["training_minutes"] = (time.time() - start) / 60
    output_directory = run_directory(arguments.output_dir, arguments.embedding, sfh_model)
    ml_networks.save_run(output_directory, estimator, config, history, standardisers)
    print(f"saved {output_directory} after {config['training_minutes']:.1f} min", flush=True)


if __name__ == "__main__":
    main()
```

  **`remote/submit_training.sh`.** Task 1 already put `sfh_models.py` in this file's `cp` (lines 39–41, same line count); the edits below apply on top of that version and do not touch `sfh_models.py` again. The line numbers are those of the file as Task 1 left it; make the edits in the order listed (bottom up), so each number still holds when it is used.
  - Lines 68–72, the submit and its echo: replace
    ```bash
    job_id="$(hljob run --detach --gpu -n "$JOB_NAME" -c "$CORES" -m "$MEMORY_MB" -t "$TIME_LIMIT" \
        -e EMBEDDINGS="$EMBEDDINGS" -e ALLOW_CPU_TRAINING="$ALLOW_CPU_TRAINING" "$@" \
        -f "$scratch_dir/bundle.tar.gz" -f "$REMOTE_DIR/train_job.sh" -- bash inputs/train_job.sh)"

    echo "Submitted job $job_id (GPU, $CORES cores, $MEMORY_MB MB, $TIME_LIMIT, embeddings: $EMBEDDINGS)"
    ```
    with
    ```bash
    job_id="$(hljob run --detach --gpu -n "$JOB_NAME" -c "$CORES" -m "$MEMORY_MB" -t "$TIME_LIMIT" \
        -e EMBEDDINGS="$EMBEDDINGS" -e SFH_MODEL="$SFH_MODEL" -e ALLOW_CPU_TRAINING="$ALLOW_CPU_TRAINING" "$@" \
        -f "$scratch_dir/bundle.tar.gz" -f "$REMOTE_DIR/train_job.sh" -- bash inputs/train_job.sh)"

    echo "Submitted job $job_id (GPU, $CORES cores, $MEMORY_MB MB, $TIME_LIMIT, SFH model: $SFH_MODEL, embeddings: $EMBEDDINGS)"
    ```
  - Lines 39–41, the `cp` of the code and sets as Task 1 wrote it:
    ```bash
    cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/sfh_models.py" "$SYNTHETIC_DIR/ml_data.py" \
       "$SYNTHETIC_DIR/ml_networks.py" "$SYNTHETIC_DIR/train_ml_model.py" "$SYNTHETIC_DIR/training_set.parquet" \
       "$SYNTHETIC_DIR/synthetic_validation.parquet" "$SYNTHETIC_DIR/splus_noise_model.csv" "$bundle_dir/"
    ```
    Only the two set names change, to the model's files:
    ```bash
    cp "$SYNTHETIC_DIR/synthetic_model.py" "$SYNTHETIC_DIR/sfh_models.py" "$SYNTHETIC_DIR/ml_data.py" \
       "$SYNTHETIC_DIR/ml_networks.py" "$SYNTHETIC_DIR/train_ml_model.py" "$SYNTHETIC_DIR/training_set${FILE_SUFFIX}.parquet" \
       "$SYNTHETIC_DIR/synthetic_validation${FILE_SUFFIX}.parquet" "$SYNTHETIC_DIR/splus_noise_model.csv" "$bundle_dir/"
    ```
  - After line 22 (`EMBEDDINGS="${EMBEDDINGS:-mlp}"`), insert:
    ```bash
    SFH_MODEL="${SFH_MODEL:-delayed}"
    # The model's file suffix, as sfh_models.MODELS[SFH_MODEL].file_suffix: its sets are training_set<suffix>.parquet
    case "$SFH_MODEL" in
        delayed) FILE_SUFFIX="" ;;
        dblplaw) FILE_SUFFIX="_dblplaw" ;;
        *) echo "ERROR: unknown SFH_MODEL '$SFH_MODEL' (use delayed or dblplaw)" >&2; exit 1 ;;
    esac
    ```
  - After line 10 (`#   EMBEDDINGS       embeddings to train, one after the other ...`), insert:
    ```bash
    #   SFH_MODEL        SFH model to train: delayed or dblplaw (sfh_models.MODELS)   [default delayed]
    ```
  - Line 4, header: replace `# at the end. The trained runs come back in out/ml_runs/<embedding>/.` with `# at the end. The trained runs come back in out/ml_runs/<embedding>/ (delayed-tau) or out/ml_runs/<embedding>_dblplaw/.`.

  **`remote/train_job.sh`** (78 lines). Make the edits in the order listed (bottom up), with today's line numbers:
  - Lines 58–70, the training loop with the `failed_embeddings=""` line before it: replace them with
    ```bash
    failed_embeddings=""
    for embedding in $EMBEDDINGS; do
        echo "== 5. Training the '$embedding' embedding of the '$SFH_MODEL' SFH model"
        set_stage "train_${embedding}_${SFH_MODEL}" "$OUT_DIR/memory_samples.csv"
        training_exit_code=0
        python train_ml_model.py --sfh-model "$SFH_MODEL" --embedding "$embedding" --output-dir "$OUT_DIR/ml_runs" \
            --device "$DEVICE" \
            >> "$OUT_DIR/train_${embedding}_${SFH_MODEL}.log" 2>&1 || training_exit_code=$?
        echo "'$embedding' finished with exit code $training_exit_code"
        if [[ "$training_exit_code" -ne 0 ]]; then
            failed_embeddings="$failed_embeddings $embedding"
        fi
    done
    ```
  - After line 15 (`EMBEDDINGS="${EMBEDDINGS:-mlp}"`), insert `SFH_MODEL="${SFH_MODEL:-delayed}"`.
  - After line 10 (`#   EMBEDDINGS   embeddings to train, one after the other   [default "mlp"]`), insert `#   SFH_MODEL    SFH model to train: delayed or dblplaw          [default delayed]`.
  - Line 7, header: replace `# out/ml_runs/<embedding>/, so a run that finished before the limit is never lost.` with `# out/ml_runs/<embedding>/ (with the model's suffix, e.g. mlp_dblplaw), so a run that finished before the limit is never lost.`.

  **`remote/README.md`.** Replace lines 87–135 after Task 6 (78–126 before it; find the block by its first line `# Training the ML posterior estimator (GPU)` and its last line `        --max-train-rows 3000 --max-validation-rows 500 --device cpu`), from `# Training the ML posterior estimator (GPU)` through the line `        --max-train-rows 3000 --max-validation-rows 500 --device cpu`, with:

````markdown
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
````

- [ ] **Step 4: Run the tests and see them pass.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_networks.py tests/test_train_ml_model.py tests/test_ml_data.py -q`
  Expected:
  - all pass;
  - `test_training_cli_writes_a_run[dblplaw]` SKIPS until the DPL sets are generated on the PC;
  - `test_training_cli_writes_a_run[delayed]` and `test_the_stored_delayed_run_loads_only_as_delayed` run and pass. The stored `ml_runs/mlp/config.json` has no `sfh_model` and still loads as delayed-τ.

- [ ] **Step 5: Check the scripts and that no global target list is left.**
  Run: `cd Synthetic_domain && bash -n remote/submit_training.sh && bash -n remote/train_job.sh && echo syntax ok`
  Expected: `syntax ok`.
  Run: `cd Synthetic_domain && grep -c 'sfh_models.py' remote/submit_training.sh && grep -n 'FILE_SUFFIX}.parquet' remote/submit_training.sh`
  Expected: `1` (Task 1's single copy of `sfh_models.py`), then the two `cp` lines with `training_set${FILE_SUFFIX}.parquet` and `synthetic_validation${FILE_SUFFIX}.parquet`.
  Run: `cd Synthetic_domain && SFH_MODEL=exponential BUNDLE_ONLY=/tmp/no_bundle bash remote/submit_training.sh; echo "exit $?"`
  Expected: `ERROR: unknown SFH_MODEL 'exponential' (use delayed or dblplaw)`, then `exit 1`. Nothing is packed or submitted: the check runs before the bundle is built.
  Run: `cd Synthetic_domain && grep -rnw -e TARGET_NAMES -e N_TARGETS --include=*.py .`
  Expected: only the `hasattr(ml_data, "TARGET_NAMES")` and `hasattr(ml_data, "N_TARGETS")` lines of `tests/test_ml_data.py`. `ml_model.ipynb` still uses them; Task 10 rewrites it.

- [ ] **Step 6: Commit.**
```bash
git add Synthetic_domain/ml_networks.py Synthetic_domain/train_ml_model.py Synthetic_domain/ml_data.py \
    Synthetic_domain/remote/submit_training.sh Synthetic_domain/remote/train_job.sh Synthetic_domain/remote/README.md \
    Synthetic_domain/tests/test_ml_networks.py Synthetic_domain/tests/test_train_ml_model.py \
    Synthetic_domain/tests/test_ml_data.py
git commit -m "feat: train_ml_model --sfh-model; runs record their SFH model and refuse to load as another

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 9: `ml_evaluation` helpers for the model comparison and the evidence

**Files:**
- Modify: `Synthetic_domain/ml_evaluation.py`:
  - the module docstring (lines 1–5);
  - section 1, new constants (after line 25);
  - two new sections at the end (after line 205): `5. Model comparison` and `6. Bayesian evidence`.
- Test: `Synthetic_domain/tests/test_ml_evaluation.py` (tests appended).

**Interfaces:**
- **Consumes:** the existing `central_interval`, `bias_and_nmad`, `agreement_verdict`, `ALL_LABELS_ROW`, `ONE_SIGMA_LEVEL` and `TWO_SIGMA_LEVEL`. None of them changes.
- **Produces, constants:**
  - `SUMMARY_SUFFIXES = ("p16", "median", "p84")`, `GAUSSIAN_95_WIDTHS = 1.959964`;
  - `VERDICT_COLUMNS = ["n_galaxies", "inside_95", "median_offset", "agrees"]`;
  - `INCONCLUSIVE_LOG_BAYES = 1.0`, `STRONG_LOG_BAYES = 3.0`;
  - `DELAYED_MODEL = "delayed"`, `DBLPLAW_MODEL = "dblplaw"`, `NO_MODEL = "neither"`;
  - `TRUTH_SIGNS`, `ALL_MOCKS_ROW = "all"`.
- **Produces, functions:**
  - `percentile_summary(samples) -> (p16, median, p84)`;
  - `compare_summaries_with_reference(p16, median, p84, reference_p16, reference_median, reference_p84) -> dict(offset_in_widths, inside_68, inside_95)`;
  - `summary_comparison_rows(galaxies, labels, estimate_prefix, reference_prefix="") -> list of dicts`. The dicts have the keys of `reference_comparison_rows`.
  - `model_difference_statistics(first_p16, first_median, first_p84, second_p16, second_median, second_p84) -> dict(n_galaxies, median_difference, scatter, overlap_68)`.
    - The scatter is the NMAD (1.4826 × the median absolute deviation) of the per-galaxy median differences around their median.
    - overlap_68 is the share of galaxies whose two [p16, p84] intervals overlap; touching intervals count as overlapping.
  - `model_comparison_rows(galaxies, labels, first_prefix, second_prefix) -> DataFrame(label, n_galaxies, median_difference, scatter, overlap_68)`;
  - `side_by_side_verdicts(verdicts: dict name -> agreement_verdict table) -> DataFrame`, indexed by label, with MultiIndex columns (name, `VERDICT_COLUMNS`);
  - `log_bayes_factor(log_evidence_delayed, log_evidence_dblplaw) -> ndarray` (positive favours delayed-τ);
  - `nested_sampling_log_evidence_error(information, n_live_points) -> sqrt(H / n_live)`;
  - `log_bayes_factor_error(error_delayed, error_dblplaw) -> sqrt(a**2 + b**2)`;
  - `evidence_strength(log_bayes) -> object array`: "inconclusive" for |ln B| < 1, "positive" for 1 ≤ |ln B| ≤ 3, "strong" for |ln B| > 3, None for NaN;
  - `favoured_model(log_bayes) -> object array` of DELAYED_MODEL, DBLPLAW_MODEL, NO_MODEL or None;
  - `evidence_category_fractions(log_bayes) -> DataFrame(favoured_model, strength, n_galaxies, fraction)`, five rows;
  - `true_model_recovery(cross_fits) -> DataFrame(truth_model, n_mocks, prefers_truth, strong_for_truth, inconclusive, strong_for_other)`, with rows delayed, dblplaw and all.

- [ ] **Step 1: Write the failing tests.** Append to `Synthetic_domain/tests/test_ml_evaluation.py`:

```python
# =============================================================================
# Model comparison and Bayesian evidence (delayed-tau against the double power law)
# =============================================================================
def test_percentile_summary_gives_the_16th_50th_and_84th_percentiles():
    samples = np.tile(np.arange(101.0), (2, 1))
    p16, median, p84 = ev.percentile_summary(samples)
    assert np.allclose(p16, 16.0)
    assert np.allclose(median, 50.0)
    assert np.allclose(p84, 84.0)


def test_compare_summaries_with_reference_extends_each_side_as_a_gaussian():
    """Estimate median 0 with p16 = -1 and p84 = +1 (+2 in the last row); reference sigma 1."""
    reference_median = np.array([0.0, 1.5, 2.5, 3.5])
    result = ev.compare_summaries_with_reference(np.full(4, -1.0), np.zeros(4), np.array([1.0, 1.0, 1.0, 2.0]),
                                                 reference_median - 1, reference_median, reference_median + 1)
    assert list(result["inside_68"]) == [True, False, False, False]
    # 95% upper bound: 1.96 for the first three rows, 3.92 for the wider upper side of the last
    assert list(result["inside_95"]) == [True, True, False, True]
    assert np.isclose(result["offset_in_widths"][0], 0.0)
    assert np.isclose(result["offset_in_widths"][1], -1.5 / np.sqrt(2))


def test_summary_comparison_rows_leave_out_galaxies_without_both_fits():
    galaxies = pd.DataFrame({
        "specobjid": ["a", "b", "c"],
        "phot_dblplaw_mw_age_p16": [-1.0, -1.0, np.nan],
        "phot_dblplaw_mw_age_median": [0.0, 0.0, np.nan],
        "phot_dblplaw_mw_age_p84": [1.0, 1.0, np.nan],
        "mw_age_p16": [-1.0, 4.0, -1.0],
        "mw_age_median": [0.0, 5.0, 0.0],
        "mw_age_p84": [1.0, 6.0, 1.0],
    })
    rows = ev.summary_comparison_rows(galaxies, ["mw_age", "log_ssfr"], "phot_dblplaw_")   # no log_ssfr columns
    assert [row["specobjid"] for row in rows] == ["a", "b"]
    assert bool(rows[0]["inside_95"]) and not bool(rows[1]["inside_95"])
    assert set(rows[0]) == {"label", "specobjid", "offset_in_widths", "inside_68", "inside_95"}


def test_model_difference_statistics_recover_a_known_offset_and_scatter():
    rng = np.random.default_rng(50)
    second_median = rng.normal(10, 1, 20000)
    first_median = second_median + 0.3 + rng.normal(0, 0.1, len(second_median))
    half_width = 0.05
    statistics = ev.model_difference_statistics(first_median - half_width, first_median, first_median + half_width,
                                                second_median - half_width, second_median,
                                                second_median + half_width)
    assert statistics["n_galaxies"] == 20000
    assert np.isclose(statistics["median_difference"], 0.3, atol=0.005)
    assert np.isclose(statistics["scatter"], 0.1, rtol=0.05)
    # The intervals overlap when the difference is at most 0.1, i.e. 2 sigma below its mean: Phi(-2) = 0.023
    assert np.isclose(statistics["overlap_68"], 0.0228, atol=0.005)


def test_model_difference_statistics_on_hand_made_intervals():
    """[0, 1] against [0.5, 2] overlaps, against [1, 2] touches (counts), against [1.1, 2] does not."""
    statistics = ev.model_difference_statistics(
        np.array([0.0, 0.0, 0.0, np.nan]), np.array([0.5, 0.5, 0.5, 0.5]), np.array([1.0, 1.0, 1.0, 1.0]),
        np.array([0.5, 1.0, 1.1, 0.0]), np.array([1.0, 1.5, 1.5, 0.5]), np.array([2.0, 2.0, 2.0, 2.0]))
    assert statistics["n_galaxies"] == 3                       # the galaxy with a NaN is left out
    assert np.isclose(statistics["overlap_68"], 2 / 3)
    assert np.isclose(statistics["median_difference"], -1.0)   # differences -0.5, -1.0, -1.0


def test_model_difference_statistics_without_galaxies_are_nan():
    missing = np.array([np.nan])
    statistics = ev.model_difference_statistics(missing, missing, missing, missing, missing, missing)
    assert statistics["n_galaxies"] == 0
    assert np.isnan(statistics["median_difference"])
    assert np.isnan(statistics["scatter"])
    assert np.isnan(statistics["overlap_68"])


def paired_galaxies(first_prefix, second_prefix, first_medians, second_medians):
    """Galaxies with {prefix}mw_age_p16/_median/_p84 (sigma 1) of two estimates."""
    columns = {}
    for prefix, medians in ((first_prefix, first_medians), (second_prefix, second_medians)):
        medians = np.asarray(medians, dtype=float)
        columns[f"{prefix}mw_age_p16"] = medians - 1
        columns[f"{prefix}mw_age_median"] = medians
        columns[f"{prefix}mw_age_p84"] = medians + 1
    return pd.DataFrame(columns)


def test_model_comparison_rows_read_both_prefixes_and_mark_missing_labels():
    galaxies = paired_galaxies("phot_", "phot_dblplaw_", [1.0, 2.0, 3.0], [0.5, 1.5, 2.5])
    table = ev.model_comparison_rows(galaxies, ["mw_age", "av_ism"], "phot_", "phot_dblplaw_")
    assert list(table.columns) == ["label", "n_galaxies", "median_difference", "scatter", "overlap_68"]
    assert list(table["label"]) == ["mw_age", "av_ism"]
    by_label = table.set_index("label")
    assert by_label.loc["mw_age", "n_galaxies"] == 3
    assert np.isclose(by_label.loc["mw_age", "median_difference"], 0.5)
    assert np.isclose(by_label.loc["mw_age", "scatter"], 0.0)
    assert by_label.loc["mw_age", "overlap_68"] == 1.0
    assert by_label.loc["av_ism", "n_galaxies"] == 0              # no av_ism columns
    assert np.isnan(by_label.loc["av_ism", "median_difference"])


def test_side_by_side_verdicts_put_each_estimator_in_its_own_column_group():
    rows = [{"label": "mw_age", "offset_in_widths": 0.0, "inside_95": True},
            {"label": "mw_age", "offset_in_widths": 0.5, "inside_95": True}]
    delayed_verdict = ev.agreement_verdict(rows, ["log_age", "mw_age"])
    dblplaw_verdict = ev.agreement_verdict(rows, ["log_alpha", "mw_age"])
    table = ev.side_by_side_verdicts({"delayed": delayed_verdict, "dblplaw": dblplaw_verdict})
    assert list(table.index) == ["log_age", "mw_age", "log_alpha", ev.ALL_LABELS_ROW]
    assert list(table.columns.get_level_values(0).unique()) == ["delayed", "dblplaw"]
    assert list(table["delayed"].columns) == ev.VERDICT_COLUMNS
    assert table.loc["mw_age", ("delayed", "inside_95")] == 1.0
    assert np.isclose(table.loc["mw_age", ("dblplaw", "median_offset")], 0.25)
    assert np.isnan(table.loc["log_alpha", ("delayed", "inside_95")])     # a label the delayed table lacks


def test_log_bayes_factor_is_delayed_minus_dblplaw_and_keeps_nan():
    log_bayes = ev.log_bayes_factor([-100.0, -50.0, np.nan], [-102.0, -45.0, -10.0])
    assert np.allclose(log_bayes[:2], [2.0, -5.0])
    assert np.isnan(log_bayes[2])


def test_evidence_strength_and_favoured_model_at_the_category_edges():
    log_bayes = np.array([0.0, 0.99, 1.0, -1.0, 2.5, 3.0, 3.01, -4.0, np.nan])
    assert list(ev.evidence_strength(log_bayes)) == ["inconclusive", "inconclusive", "positive", "positive",
                                                     "positive", "positive", "strong", "strong", None]
    assert list(ev.favoured_model(log_bayes)) == ["neither", "neither", "delayed", "dblplaw", "delayed",
                                                  "delayed", "delayed", "dblplaw", None]


def test_evidence_category_fractions_of_a_hand_made_table():
    log_bayes = np.array([0.5, -0.2, 0.0, 2.0, 5.0, 4.0, -1.5, -3.5, -10.0, np.nan])
    table = ev.evidence_category_fractions(log_bayes)
    assert list(zip(table["favoured_model"], table["strength"])) == [
        ("neither", "inconclusive"), ("delayed", "positive"), ("delayed", "strong"),
        ("dblplaw", "positive"), ("dblplaw", "strong")]
    assert list(table["n_galaxies"]) == [3, 1, 2, 1, 2]
    assert np.allclose(table["fraction"], np.array([3, 1, 2, 1, 2]) / 9)     # the NaN galaxy is left out
    assert np.isclose(table["fraction"].sum(), 1.0)


def test_evidence_category_fractions_without_finite_values_are_nan():
    table = ev.evidence_category_fractions(np.array([np.nan, np.nan]))
    assert (table["n_galaxies"] == 0).all()
    assert table["fraction"].isna().all()


def test_evidence_uncertainties():
    assert np.isclose(ev.nested_sampling_log_evidence_error(12.0, 300), 0.2)
    assert np.isclose(ev.log_bayes_factor_error(0.3, 0.4), 0.5)
    assert np.allclose(ev.log_bayes_factor_error(np.array([0.3, 0.0]), np.array([0.4, 0.1])), [0.5, 0.1])


def cross_fit_table(truth_models, log_bayes):
    """Cross-fitted mocks with the given true models and ln B (ln Z of the DPL fit fixed at -100)."""
    log_bayes = np.asarray(log_bayes, dtype=float)
    mock_ids = []
    for index in range(len(log_bayes)):
        mock_ids.append(f"calib_{index:03d}")
    return pd.DataFrame({"mock_id": mock_ids, "truth_model": truth_models,
                         "log_evidence_delayed": -100.0 + log_bayes,
                         "log_evidence_dblplaw": np.full(len(log_bayes), -100.0)})


def test_true_model_recovery_of_hand_made_mocks():
    table = cross_fit_table(["delayed"] * 4 + ["dblplaw"] * 4, [5.0, 2.0, -0.5, -4.0, -6.0, -0.2, 1.5, np.nan])
    recovery = ev.true_model_recovery(table).set_index("truth_model")
    assert list(recovery.index) == ["delayed", "dblplaw", ev.ALL_MOCKS_ROW]
    assert list(recovery["n_mocks"]) == [4, 3, 7]                            # the NaN mock is left out
    assert np.allclose(recovery["prefers_truth"], [2 / 4, 2 / 3, 4 / 7])
    assert np.allclose(recovery["strong_for_truth"], [1 / 4, 1 / 3, 2 / 7])
    assert np.allclose(recovery["inconclusive"], [1 / 4, 1 / 3, 2 / 7])
    assert np.allclose(recovery["strong_for_other"], [1 / 4, 0.0, 1 / 7])


def test_true_model_recovery_without_mocks_is_nan():
    recovery = ev.true_model_recovery(cross_fit_table([], []))
    assert (recovery["n_mocks"] == 0).all()
    assert recovery["prefers_truth"].isna().all()
```

- [ ] **Step 2: Run the tests and see them fail.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_evaluation.py -q`
  Expected: `15 failed, 15 passed`. The 15 new tests FAIL with `AttributeError: module 'ml_evaluation' has no attribute 'percentile_summary'` (or `compare_summaries_with_reference`, `model_difference_statistics`, ...). The 15 existing tests pass.

- [ ] **Step 3: Implement `ml_evaluation.py`.**

  Replace the module docstring (lines 1–5) with:

```python
"""
Metrics for the SFH posterior estimator: calibration (coverage, ranks), accuracy (bias, NMAD),
information per label (posterior/prior width) and agreement with the Bagpipes reference posteriors
(spectral or photometric fits). Also the comparison of the two SFH models (delayed-tau against the double
power law) on their shared labels, and their Bayesian evidence.
"""
```

  After line 25 (`ALL_LABELS_ROW = "all"  ...`), insert:

```python
SUMMARY_SUFFIXES = ("p16", "median", "p84")  # column suffixes of a posterior stored as percentiles
GAUSSIAN_95_WIDTHS = 1.959964                # [1-sigma widths] half-width of the central 95% of a Gaussian
VERDICT_COLUMNS = ["n_galaxies", "inside_95", "median_offset", "agrees"]   # agreement_verdict per label

# Evidence scale on |ln B| (Jeffreys, as in Trotta 2008): below 1 inconclusive, 1 to 3 positive, above 3 strong
INCONCLUSIVE_LOG_BAYES = 1.0                 # [nats] |ln B| below this is inconclusive
STRONG_LOG_BAYES = 3.0                       # [nats] |ln B| above this is strong
DELAYED_MODEL = "delayed"                    # favoured by ln B >= INCONCLUSIVE_LOG_BAYES
DBLPLAW_MODEL = "dblplaw"                    # favoured by ln B <= -INCONCLUSIVE_LOG_BAYES
NO_MODEL = "neither"                         # an inconclusive ln B favours no model
TRUTH_SIGNS = {DELAYED_MODEL: 1.0, DBLPLAW_MODEL: -1.0}   # sign of ln B that favours each true model
ALL_MOCKS_ROW = "all"                        # true_model_recovery row over every mock
```

  Append at the end of the file (after `agreement_verdict`):

```python
# =============================================================================
# 5. Model comparison: delayed-tau against the double power law
# =============================================================================
def percentile_summary(samples):
    """16th, 50th and 84th percentiles of each galaxy's samples (axis -1), as the reference CSV stores fits.

    Inputs: samples (n, m) of one label. Returns: (p16, median, p84), arrays (n,) in the label's units.
    """
    low, high = central_interval(samples, ONE_SIGMA_LEVEL)
    return low, np.nanmedian(samples, axis=-1), high


def compare_summaries_with_reference(p16, median, p84, reference_p16, reference_median, reference_p84):
    """compare_with_spectral for an estimate known only by its 16th, 50th and 84th percentiles.

    Why the Gaussian extension: a Bagpipes fit in the reference CSV keeps only these percentiles. Its 95%
    interval is therefore taken as median - 1.96 (median - p16) to median + 1.96 (p84 - median): each side is
    scaled from its own 1-sigma width.
    Inputs: the estimate's and the reference's percentiles, arrays (n,) in the label's units.
    Returns: dict with
             - offset_in_widths: estimate minus reference median, in combined 1-sigma widths;
             - inside_68 and inside_95: whether the reference median lies inside the estimate's interval.
    """
    p16 = np.asarray(p16, dtype=float)
    median = np.asarray(median, dtype=float)
    p84 = np.asarray(p84, dtype=float)
    reference_median = np.asarray(reference_median, dtype=float)
    sigma = (p84 - p16) / 2
    reference_sigma = (np.asarray(reference_p84, dtype=float) - np.asarray(reference_p16, dtype=float)) / 2
    combined = np.sqrt(sigma**2 + reference_sigma**2)
    offset = (median - reference_median) / combined

    low_95 = median - GAUSSIAN_95_WIDTHS * (median - p16)
    high_95 = median + GAUSSIAN_95_WIDTHS * (p84 - median)
    inside_68 = (reference_median >= p16) & (reference_median <= p84)
    inside_95 = (reference_median >= low_95) & (reference_median <= high_95)
    return {"offset_in_widths": offset, "inside_68": inside_68, "inside_95": inside_95}


def summary_columns(prefix, label):
    """The {prefix}<label>_p16 / _median / _p84 column names of one stored posterior."""
    columns = []
    for suffix in SUMMARY_SUFFIXES:
        columns.append(f"{prefix}{label}_{suffix}")
    return columns


def summary_comparison_rows(galaxies, labels, estimate_prefix, reference_prefix=""):
    """reference_comparison_rows for an estimate stored as percentiles (e.g. a photometric Bagpipes fit).

    Inputs: galaxies: DataFrame with specobjid, the estimate's {estimate_prefix}<label>_p16/_median/_p84
            columns and the reference's {reference_prefix}<label>_* columns (by default the spectral fit's);
            labels: the labels to compare.
    Returns: list of dicts with label, specobjid, offset_in_widths, inside_68 and inside_95, labels outer and
             galaxies inner. A galaxy without both medians of a label is left out of that label. A label
             whose columns are missing gives no rows.
    """
    rows = []
    for label in labels:
        estimate_columns = summary_columns(estimate_prefix, label)
        reference_columns = summary_columns(reference_prefix, label)
        if not set(estimate_columns + reference_columns) <= set(galaxies.columns):
            continue
        both = galaxies[estimate_columns[1]].notna() & galaxies[reference_columns[1]].notna()
        compared = galaxies[both].reset_index(drop=True)
        comparison = compare_summaries_with_reference(
            compared[estimate_columns[0]].to_numpy(), compared[estimate_columns[1]].to_numpy(),
            compared[estimate_columns[2]].to_numpy(), compared[reference_columns[0]].to_numpy(),
            compared[reference_columns[1]].to_numpy(), compared[reference_columns[2]].to_numpy())
        for galaxy_index in range(len(compared)):
            rows.append({"label": label, "specobjid": compared.loc[galaxy_index, "specobjid"],
                         "offset_in_widths": comparison["offset_in_widths"][galaxy_index],
                         "inside_68": comparison["inside_68"][galaxy_index],
                         "inside_95": comparison["inside_95"][galaxy_index]})
    return rows


def model_difference_statistics(first_p16, first_median, first_p84, second_p16, second_median, second_p84):
    """How two estimates of one label differ over galaxies, first minus second.

    Why NMAD for the scatter: a few catastrophic galaxies should not set it (as in bias_and_nmad).
    Inputs: the 16th, 50th and 84th percentiles of each estimate, arrays (n,) in the label's units, paired by
            galaxy. A galaxy missing any value is left out.
    Returns: dict with
             - n_galaxies (count);
             - median_difference: median of first minus second medians, in label units;
             - scatter: NMAD of that difference around its median, in label units;
             - overlap_68: share of galaxies whose two 16-84% intervals overlap (touching counts).
             The statistics are NaN without galaxies.
    """
    first_p16 = np.asarray(first_p16, dtype=float)
    first_median = np.asarray(first_median, dtype=float)
    first_p84 = np.asarray(first_p84, dtype=float)
    second_p16 = np.asarray(second_p16, dtype=float)
    second_median = np.asarray(second_median, dtype=float)
    second_p84 = np.asarray(second_p84, dtype=float)
    usable = (np.isfinite(first_p16) & np.isfinite(first_median) & np.isfinite(first_p84)
              & np.isfinite(second_p16) & np.isfinite(second_median) & np.isfinite(second_p84))
    n_galaxies = int(usable.sum())
    if n_galaxies == 0:
        return {"n_galaxies": 0, "median_difference": np.nan, "scatter": np.nan, "overlap_68": np.nan}

    median_difference, scatter = bias_and_nmad(first_median[usable], second_median[usable])
    overlap = (first_p16[usable] <= second_p84[usable]) & (second_p16[usable] <= first_p84[usable])
    return {"n_galaxies": n_galaxies, "median_difference": float(median_difference), "scatter": float(scatter),
            "overlap_68": float(np.mean(overlap))}


def model_comparison_rows(galaxies, labels, first_prefix, second_prefix):
    """model_difference_statistics of each label, for two estimates stored side by side in one table.

    Inputs: galaxies: DataFrame with {prefix}<label>_p16/_median/_p84 of both estimates, one row per galaxy;
            labels: the labels to compare; first_prefix, second_prefix: their column prefixes (e.g. "phot_" and
            "phot_dblplaw_"). The difference is first minus second.
    Returns: DataFrame with one row per label: label, n_galaxies, median_difference, scatter and overlap_68.
             A label whose columns are missing for either estimate gets n_galaxies 0 and NaN statistics.
    """
    rows = []
    for label in labels:
        first_columns = summary_columns(first_prefix, label)
        second_columns = summary_columns(second_prefix, label)
        if not set(first_columns + second_columns) <= set(galaxies.columns):
            rows.append({"label": label, "n_galaxies": 0, "median_difference": np.nan, "scatter": np.nan,
                         "overlap_68": np.nan})
            continue
        statistics = model_difference_statistics(
            galaxies[first_columns[0]], galaxies[first_columns[1]], galaxies[first_columns[2]],
            galaxies[second_columns[0]], galaxies[second_columns[1]], galaxies[second_columns[2]])
        rows.append({"label": label, **statistics})
    return pd.DataFrame(rows, columns=["label", "n_galaxies", "median_difference", "scatter", "overlap_68"])


def side_by_side_verdicts(verdicts):
    """agreement_verdict tables of several estimators in one table: one row per label, one column group each.

    Inputs: verdicts: dict name -> agreement_verdict DataFrame, in display order.
    Returns: DataFrame indexed by label (first-seen order, ALL_LABELS_ROW last) with MultiIndex columns
             (name, VERDICT_COLUMNS). A label an estimator's table lacks is NaN there.
    """
    label_order = []
    has_all_row = False
    for verdict in verdicts.values():
        for label in verdict["label"]:
            if label == ALL_LABELS_ROW:
                has_all_row = True
            elif label not in label_order:
                label_order.append(label)
    if has_all_row:
        label_order.append(ALL_LABELS_ROW)

    tables = {}
    for name, verdict in verdicts.items():
        tables[name] = verdict.set_index("label")[VERDICT_COLUMNS].reindex(label_order)
    return pd.concat(tables, axis=1)


# =============================================================================
# 6. Bayesian evidence
# =============================================================================
def log_bayes_factor(log_evidence_delayed, log_evidence_dblplaw):
    """ln B = ln Z(delayed-tau) - ln Z(double power law) of each galaxy; positive favours delayed-tau.

    Inputs: the natural-log evidences of the two models' fits of the same data, arrays (n,) [nats].
    Returns: array (n,) [nats]; NaN where either evidence is missing.
    """
    return np.asarray(log_evidence_delayed, dtype=float) - np.asarray(log_evidence_dblplaw, dtype=float)


def nested_sampling_log_evidence_error(information, n_live_points):
    """Statistical 1-sigma uncertainty of a nested-sampling ln Z: sqrt(H / n_live) (Skilling 2006).

    Inputs: information: H, the Kullback-Leibler divergence from prior to posterior [nats];
            n_live_points: live points of the run (count).
    Returns: uncertainty of ln Z [nats].
    """
    return np.sqrt(np.asarray(information, dtype=float) / n_live_points)


def log_bayes_factor_error(log_evidence_error_delayed, log_evidence_error_dblplaw):
    """1-sigma uncertainty of ln B [nats]: the two fits are independent runs, so their errors add in quadrature."""
    error_delayed = np.asarray(log_evidence_error_delayed, dtype=float)
    error_dblplaw = np.asarray(log_evidence_error_dblplaw, dtype=float)
    return np.sqrt(error_delayed**2 + error_dblplaw**2)


def evidence_strength(log_bayes):
    """Strength of each ln B on the |ln B| scale.

    The categories are "inconclusive" (below INCONCLUSIVE_LOG_BAYES), "positive" (from it to
    STRONG_LOG_BAYES, both inclusive) and "strong" (above STRONG_LOG_BAYES).
    Inputs: log_bayes: array (n,) of ln B [nats].
    Returns: object array (n,) of strengths; None where ln B is NaN.
    """
    log_bayes = np.asarray(log_bayes, dtype=float)
    size = np.abs(log_bayes)
    strengths = np.full(len(log_bayes), None, dtype=object)
    strengths[size < INCONCLUSIVE_LOG_BAYES] = "inconclusive"
    strengths[(size >= INCONCLUSIVE_LOG_BAYES) & (size <= STRONG_LOG_BAYES)] = "positive"
    strengths[size > STRONG_LOG_BAYES] = "strong"
    return strengths


def favoured_model(log_bayes):
    """The model each ln B favours.

    Inputs: log_bayes: array (n,) of ln B [nats].
    Returns: object array (n,):
             - DELAYED_MODEL for ln B >= INCONCLUSIVE_LOG_BAYES;
             - DBLPLAW_MODEL for ln B <= -INCONCLUSIVE_LOG_BAYES;
             - NO_MODEL in between;
             - None where ln B is NaN.
    """
    log_bayes = np.asarray(log_bayes, dtype=float)
    models = np.full(len(log_bayes), None, dtype=object)
    models[np.abs(log_bayes) < INCONCLUSIVE_LOG_BAYES] = NO_MODEL
    models[log_bayes >= INCONCLUSIVE_LOG_BAYES] = DELAYED_MODEL
    models[log_bayes <= -INCONCLUSIVE_LOG_BAYES] = DBLPLAW_MODEL
    return models


def evidence_category_fractions(log_bayes):
    """Share of galaxies in each evidence category: inconclusive, then positive and strong for each model.

    Inputs: log_bayes: array (n,) of ln B (positive favours delayed-tau). NaN entries are left out.
    Returns: DataFrame with 5 rows (favoured_model, strength, n_galaxies, fraction), in this order:
             (NO_MODEL, inconclusive), (DELAYED_MODEL, positive), (DELAYED_MODEL, strong),
             (DBLPLAW_MODEL, positive), (DBLPLAW_MODEL, strong). Fractions are of the finite ln B and sum
             to 1; they are NaN when no ln B is finite.
    """
    log_bayes = np.asarray(log_bayes, dtype=float)
    finite = log_bayes[np.isfinite(log_bayes)]
    strengths = evidence_strength(finite)
    models = favoured_model(finite)
    categories = [(NO_MODEL, "inconclusive"), (DELAYED_MODEL, "positive"), (DELAYED_MODEL, "strong"),
                  (DBLPLAW_MODEL, "positive"), (DBLPLAW_MODEL, "strong")]
    rows = []
    for model, strength in categories:
        n_galaxies = int(np.sum((models == model) & (strengths == strength)))
        fraction = np.nan
        if len(finite) > 0:
            fraction = n_galaxies / len(finite)
        rows.append({"favoured_model": model, "strength": strength, "n_galaxies": n_galaxies,
                     "fraction": fraction})
    return pd.DataFrame(rows)


def true_model_recovery(cross_fits):
    """How often the evidence prefers the true model of mock galaxies fitted with both models.

    Inputs: cross_fits: DataFrame with one row per mock and the columns
            - truth_model (DELAYED_MODEL or DBLPLAW_MODEL);
            - log_evidence_delayed and log_evidence_dblplaw: ln Z of its fit with each model [nats]
              (ml_model.ipynb renames phot_log_evidence and phot_dblplaw_log_evidence of
              photometric_cross_fits.csv to these).
            Rows with a missing evidence are left out.
    Returns: DataFrame with one row per true model (DELAYED_MODEL, DBLPLAW_MODEL), then ALL_MOCKS_ROW.
             Its columns are truth_model and n_mocks, then four fractions of the mocks:
             - prefers_truth: ln B has the true model's sign;
             - strong_for_truth: |ln B| > STRONG_LOG_BAYES in favour of the true model;
             - inconclusive: |ln B| < INCONCLUSIVE_LOG_BAYES;
             - strong_for_other: |ln B| > STRONG_LOG_BAYES in favour of the wrong model.
             The fractions are NaN for a group without mocks.
    """
    log_bayes = log_bayes_factor(cross_fits["log_evidence_delayed"], cross_fits["log_evidence_dblplaw"])
    truth = cross_fits["truth_model"].to_numpy(dtype=object)
    usable = np.isfinite(log_bayes)

    # ln B with the sign that favours the true model: positive = the evidence prefers the truth
    signed = np.full(len(log_bayes), np.nan)
    for model, sign in TRUTH_SIGNS.items():
        is_model = truth == model
        signed[is_model] = sign * log_bayes[is_model]

    groups = {DELAYED_MODEL: usable & (truth == DELAYED_MODEL),
              DBLPLAW_MODEL: usable & (truth == DBLPLAW_MODEL),
              ALL_MOCKS_ROW: usable & np.isin(truth, list(TRUTH_SIGNS))}
    rows = []
    for name, in_group in groups.items():
        n_mocks = int(in_group.sum())
        if n_mocks == 0:
            rows.append({"truth_model": name, "n_mocks": 0, "prefers_truth": np.nan, "strong_for_truth": np.nan,
                         "inconclusive": np.nan, "strong_for_other": np.nan})
            continue
        group_signed = signed[in_group]
        rows.append({"truth_model": name, "n_mocks": n_mocks,
                     "prefers_truth": float(np.mean(group_signed > 0)),
                     "strong_for_truth": float(np.mean(group_signed > STRONG_LOG_BAYES)),
                     "inconclusive": float(np.mean(np.abs(group_signed) < INCONCLUSIVE_LOG_BAYES)),
                     "strong_for_other": float(np.mean(group_signed < -STRONG_LOG_BAYES))})
    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run the tests and see them pass.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_evaluation.py -q`
  Expected: `30 passed` (the 15 existing tests and the 15 new ones).

- [ ] **Step 5: Run the other ML test files.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_evaluation.py tests/test_ml_data.py tests/test_ml_networks.py -q`
  Expected: all pass. Skips are listed in Tasks 7 and 8.

- [ ] **Step 6: Commit.**
```bash
git add Synthetic_domain/ml_evaluation.py Synthetic_domain/tests/test_ml_evaluation.py
git commit -m "feat: ml_evaluation helpers for the delayed-tau vs double-power-law comparison and the Bayesian evidence

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 10: `ml_model.ipynb` evaluates both SFH models and compares them

**Files:**
- Modify: `Synthetic_domain/ml_model.ipynb`:
  - cells 0–13 are rewritten in place, by cell id;
  - 11 new cells are inserted after the last cell (id `2e021085`).
- No test file. The check is that the notebook is valid JSON, every code cell compiles and no name is undefined (Step 3). **Do not execute the notebook on the laptop.** It runs on the PC (GPU) as an `hljob` notebook job with `GPU=1 CORES=16 TIME_LIMIT=12:00:00` (`remote/submit_notebook.sh`). Why: decoding DPL labels sums Bagpipes' age grid for every sample (about 0.06 ms each, Task 2 Step 5), and the notebook decodes about 8e7 DPL samples (3.3e7 for the synthetic validation set, 4.5e7 for DR3), another 5e7 with `RUN_LOCKED_TEST`: about 1.5 h in one thread, 2.5 h with the locked test. `N_DECODE_JOBS = 15` decoding threads (`ml_data.decode_samples(..., n_jobs=N_DECODE_JOBS)`) bring that under half an hour. The 12 h limit, not the 04:00:00 default of `submit_notebook.sh`, leaves room for the GPU sampling, the two DR3 coverage tests and the locked test: a timeout loses the whole GPU run. The notebook prints the decoding rate (samples/s) of each model, so the first run shows the real cost. GPU jobs are capped at about 7.1 GB (`MEMORY_MB=7100`). The decoding threads share the process, so they add only their 32-sample grid blocks (0.42 MB per array each). To stay under the cap, the notebook keeps only one model's synthetic posteriors at a time: after a model's metrics and figures, it deletes that model's posterior arrays (`del` then `gc.collect()`) and keeps only the summaries. `out/memory_peak_bytes.txt` of the job shows the real peak.

**Interfaces:**
- **Consumes:**
  - `ml_data.target_names`, `ml_data.decode_samples(theta, sfh_model, n_jobs=1)` and `ml_data.validation_observations(held_out, mode, rng, sfh_model)` (Task 7);
  - `ml_networks.load_run(directory, sfh_model, device=...)`, whose returned config has `sfh_model` and `n_targets` (Task 8);
  - the `ev.*` helpers of Task 9;
  - `sfh_models.MODELS[...]` with `.training_labels` and `.calibration_labels`, `sfh_models.SHARED_LABELS` and `sfh_models.model_file` (Tasks 1–2);
  - `sm.PHOTOMETRIC_PREFIXES = {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}` (Task 5), the one source of the photometric column prefixes;
  - **the reference-CSV columns of Task 6:**
    - `phot_dblplaw_<label>_p16/_median/_p84`;
    - `phot_dblplaw_converged` and `phot_dblplaw_fitter_calibrated`;
    - `phot_log_evidence` and `phot_dblplaw_log_evidence`;
    - optionally `phot_log_evidence_error` and `phot_dblplaw_log_evidence_error`;
    - for the `test` half, these columns hold values only when notebook 2 ran with `EXPORT_LOCKED_TEST_JOBS = True` and read the test-half fits (Task 6). Section 5 checks for them and otherwise covers the `validation` half only;
  - **`photometric_cross_fits.csv` of Task 6** (notebook 2, section 6c), one row per calibration mock of either model: `mock_id`, `truth_model`, `phot_converged`, `phot_log_evidence`, `phot_log_evidence_error`, `phot_dblplaw_converged`, `phot_dblplaw_log_evidence`, `phot_dblplaw_log_evidence_error`. Notebook 2 has already checked `sm.fit_is_current` for every fit it read, so this notebook reads no `.h5` file.
- **Produces (outputs written when the notebook runs on the PC):**
  - `ml_evaluation_summary.csv` and `ml_evaluation_summary_dblplaw.csv`;
  - `ml_sfh_catalogue.parquet` and `ml_sfh_catalogue_dblplaw.parquet`;
  - `ml_model_comparison.csv`: a new long table with a `section` column (own_photometric, model_pairs, against_spectra, evidence_distribution, evidence_categories, true_model_recovery).

- [ ] **Step 1: Rewrite cells 0–13 in place.** Use NotebookEdit with `edit_mode: replace` and the `cell_id` given; keep each cell's type. Each source below is the complete new cell.

  **Cell 0** (markdown, id `55011031`, "# SFH posterior from S-PLUS colours (part 3)"):

```markdown
# SFH posterior from S-PLUS colours (part 3)

Two SFH models, each with its own normalizing flow (trained by `train_ml_model.py --sfh-model <model>` on the 3090):
- **delayed-τ** (`delayed`): 9 targets, which are log age, log τ, log metallicity, the dust (`av_ism`, `log_tau_bc`), `logU`, the formed and the surviving mass, and the redshift;
- **double power law** (`dblplaw`): 10 targets, which are log turnover lookback time (the lookback time to τ, the turnover of the double power law), log α and log β, then the delayed-τ targets from log metallicity on.

Every other label (`mw_age`, `f_young`, `f_mid`, `f_old`, `log_sfr`, `log_ssfr`, ...) is derived from the samples with the same functions as the training labels (`ml_data.decode_samples`). Inputs are the 11 colours against r, the r magnitude (only for the absolute masses), and z_obs with its uncertainty (0 for spectroscopic redshifts).

Sections 1–4 run once per model in `SFH_MODELS`. Section 5 compares the two models and needs both.

**Inputs:** delayed-τ keeps the plain names; DPL adds `_dblplaw` (`sfh_models.model_file`).
- Per model: `ml_runs/mlp[_dblplaw]/`, `synthetic_validation[_dblplaw].parquet`, `synthetic_test[_dblplaw].parquet`, `training_set[_dblplaw].parquet` (prior draws).
- Shared: `spectroscopic_reference.csv`, `photometric_cross_fits.csv` (ln Z of the cross-fitted calibration mocks, notebook 2 section 6c; section 5.4), `dr3_corrected.csv`.

**Outputs:** `ml_evaluation_summary[_dblplaw].csv`, `ml_sfh_catalogue[_dblplaw].parquet`, `ml_model_comparison.csv`.

**Run:** top to bottom with the *Inference4ETGs (MultiNest)* kernel, on the PC, as a GPU notebook job with 16 cores and a 12 h limit (`GPU=1 CORES=16 TIME_LIMIT=12:00:00 remote/submit_notebook.sh ml_model.ipynb ...`). Decoding a DPL sample sums Bagpipes' age grid (about 0.06 ms), so the decoding runs in `N_DECODE_JOBS` threads. Sections marked **LOCKED** run only with `RUN_LOCKED_TEST = True`, once, for both models, after the user has inspected the validation results. They use the synthetic test sets and the reference set's test half. The photometric Bagpipes fits of the test half exist only when notebook 2 exported and read them (`EXPORT_LOCKED_TEST_JOBS = True`); without them the photometric and evidence comparisons of section 5 cover the validation half.

Posterior samples for the same seed differ between CUDA and CPU (`DEVICE`): numbers are reproducible per device, not across devices. `ml_networks.load_run` refuses a checkpoint trained for another SFH model or other targets (e.g. the 16-bin runs): retrain it.
```

  **Cell 1** (code, id `b6229f82`, "import time"):

```python
import os
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import torch

import ml_data
import ml_evaluation as ev
import ml_networks
import sfh_models
import synthetic_model as sm

RUN_LOCKED_TEST = False            # set True once, for the final report (both models)
SFH_MODELS = ["delayed", "dblplaw"]   # evaluated in sections 1-4; section 5 needs both
EMBEDDING = "mlp"                  # the CNN lost in the 16-bin comparison, so only the MLP is trained
RUN_DIRECTORIES = {}
for sfh_model in SFH_MODELS:
    RUN_DIRECTORIES[sfh_model] = sfh_models.model_file(f"ml_runs/{EMBEDDING}", sfh_model)
# Per-model files: these names are the delayed-tau ones; sfh_models.model_file adds the other models' suffix
VALIDATION_FILE = "synthetic_validation.parquet"
TEST_FILE = "synthetic_test.parquet"
TRAINING_FILE = "training_set.parquet"
SUMMARY_FILE = "ml_evaluation_summary.csv"
CATALOGUE_FILE = "ml_sfh_catalogue.parquet"
# Shared files. REFERENCE_FILE holds each model's photometric Bagpipes fits under sm.PHOTOMETRIC_PREFIXES
REFERENCE_FILE = "spectroscopic_reference.csv"
CROSS_FITS_FILE = "photometric_cross_fits.csv"   # ln Z of every calibration mock under both models (notebook 2, 6c)
SPLUS_FILE = "dr3_corrected.csv"
COMPARISON_FILE = "ml_model_comparison.csv"

N_POSTERIOR_SAMPLES = 1000         # per galaxy for evaluation
N_CATALOGUE_SAMPLES = 500          # per galaxy for the DR3 catalogue
SAMPLING_SEED = 7
VALIDATION_SEEDS = {"spectroscopic": 100, "photometric": 101}   # same copies as training used
TEST_SEEDS = {"spectroscopic": 200, "photometric": 201}
PRIOR_WINDOW_DZ = 0.05             # prior draws within this redshift of a galaxy define its prior width
MIN_PROB_GALAXY = 0.8
PASSIVE_SSFR_LIMIT = -11            # [log10 yr^-1] spectral sSFR below this = passive galaxy
PRIOR_SEED = 11                     # brightness offsets of the prior draws
N_JOBS = 11                         # processes of the DR3 coverage test
N_DECODE_JOBS = 15                  # threads decoding posterior samples: a DPL sample sums Bagpipes' age grid
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"   # flow sampling is autoregressive: slow on CPU


def evaluated_labels(sfh_model):
    """The label lists one SFH model is evaluated on.

    Returns: dict with
      flow_target_labels: its flow targets as decoded labels. decode_samples returns tau_bc but not
        log_tau_bc, so that target is evaluated as tau_bc;
      key_labels: those, then every derived label not among them (mw_age, f_young, f_mid, f_old, log_sfr,
        log_ssfr). The linear SFH parameters are left out: they are monotone maps of the log targets, so
        their ranks are the same and their coverage nearly so;
      judged_labels: the model's calibration labels, which the "works" test judges;
      reported_labels: the shared labels (sfh_models.SHARED_LABELS) not among them, reported but not judged;
      reference_labels: the judged labels, then the reported ones, compared with the Bagpipes posteriors.
    """
    flow_target_labels = []
    for name in ml_data.target_names(sfh_model):
        if name == "log_tau_bc":
            flow_target_labels.append("tau_bc")
        else:
            flow_target_labels.append(name)
    key_labels = [*flow_target_labels, "mw_age", "f_young", "f_mid", "f_old", "log_sfr", "log_ssfr"]
    judged_labels = list(sfh_models.MODELS[sfh_model].calibration_labels)
    reported_labels = []
    for label in sfh_models.SHARED_LABELS:
        if label not in judged_labels:
            reported_labels.append(label)
    return {"flow_target_labels": flow_target_labels, "key_labels": key_labels,
            "judged_labels": judged_labels, "reported_labels": reported_labels,
            "reference_labels": judged_labels + reported_labels}


LABELS = {}
runs = {}
for sfh_model in SFH_MODELS:
    LABELS[sfh_model] = evaluated_labels(sfh_model)
    wanted = set(LABELS[sfh_model]["key_labels"]) | set(LABELS[sfh_model]["reference_labels"])
    undecoded = wanted - set(sfh_models.MODELS[sfh_model].training_labels)
    if undecoded:                  # fail here, before any GPU sampling, not with a KeyError in the evaluation
        raise ValueError(f"{sfh_model}: labels the decoder does not return: {sorted(undecoded)}")
    runs[sfh_model] = ml_networks.load_run(RUN_DIRECTORIES[sfh_model], sfh_model, device=DEVICE)
```

  **Cell 2** (markdown, id `b39e2823`, "## 1. Training and validation"):

```markdown
## 1. Training and validation

The MLP is the only embedding (the CNN lost in the 16-bin comparison). For each SFH model, this section shows the training history and the coverage of the key labels on the model's synthetic validation set. The two models' losses are not comparable: they model different targets.
```

  **Cell 3** (code, id `59be2b0a`, "GALAXY_CHUNK = 500"):

```python
GALAXY_CHUNK = 500   # galaxies decoded at once: all labels of 16k x 1000 samples would need ~2.6 GB


def labelled_posteriors(run, observations, n_samples, seed, labels=None):
    """Decoded posterior labels for every valid galaxy: dict label -> array (n_valid, n_samples).

    run: (estimator, config, standardisers) from ml_networks.load_run. Its config names the SFH model the
    samples are decoded with. labels: names to keep, by default all of the model's training labels.
    Decoding runs in galaxy chunks, so only the kept labels are ever held for the whole set; each chunk is
    decoded in N_DECODE_JOBS threads (a DPL sample sums Bagpipes' age grid, about 0.1 ms).
    Returns: (posteriors, truths); truths are the decoded observations["targets"] of the valid galaxies.
    Prints the sampling and decoding time.
    """
    estimator, config, standardisers = run
    sfh_model = config["sfh_model"]
    valid = observations["valid"]
    if not valid.any():
        raise ValueError("no valid galaxy (every galaxy lacks a detected r band or a brightness donor): "
                         "nothing to predict")
    flat = observations["flat"][valid]
    banded = observations["banded"][valid]
    n_galaxies = len(flat)

    posteriors = {}
    sampling_seconds = 0.0
    decoding_seconds = 0.0
    for chunk_start in range(0, n_galaxies, GALAXY_CHUNK):
        chunk = slice(chunk_start, chunk_start + GALAXY_CHUNK)
        start = time.time()
        theta = ml_networks.posterior_samples(estimator, standardisers, flat[chunk], banded[chunk],
                                              n_samples, seed + chunk_start, DEVICE)
        sampling_seconds += time.time() - start
        n_chunk = theta.shape[0]
        start = time.time()
        decoded = ml_data.decode_samples(theta.reshape(-1, config["n_targets"]), sfh_model, n_jobs=N_DECODE_JOBS)
        decoding_seconds += time.time() - start
        for column in decoded.columns:
            if labels is not None and column not in labels:
                continue
            if column not in posteriors:
                posteriors[column] = np.empty((n_galaxies, n_samples))
            posteriors[column][chunk] = decoded[column].to_numpy().reshape(n_chunk, n_samples)
    truths = ml_data.decode_samples(observations["targets"][valid], sfh_model, n_jobs=N_DECODE_JOBS)
    n_decoded = n_galaxies * n_samples
    print(f"{sfh_model}: {n_galaxies} galaxies x {n_samples} samples, sampled in {sampling_seconds:.0f} s on "
          f"{DEVICE}, decoded in {decoding_seconds:.0f} s ({n_decoded / max(decoding_seconds, 1e-9):.0f} samples/s, "
          f"{N_DECODE_JOBS} threads)")
    return posteriors, truths


figure, axis = plt.subplots(figsize=(6, 4))
validation_tables = {}
choices = {}                       # sfh_model -> coverage of the key labels on the synthetic validation set
for sfh_model, run in runs.items():
    validation_tables[sfh_model] = pd.read_parquet(sfh_models.model_file(VALIDATION_FILE, sfh_model))
    history = pd.read_csv(os.path.join(RUN_DIRECTORIES[sfh_model], ml_networks.HISTORY_FILE))
    axis.plot(history["epoch"], history["validation_loss"], label=f"{sfh_model} validation")
    axis.plot(history["epoch"], history["train_loss"], ls=":", label=f"{sfh_model} train")
    choice_rows = []
    for mode, seed in VALIDATION_SEEDS.items():
        observations = ml_data.validation_observations(validation_tables[sfh_model], mode,
                                                       np.random.default_rng(seed), sfh_model)
        posteriors, truths = labelled_posteriors(run, observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED,
                                                 LABELS[sfh_model]["key_labels"])
        for label in LABELS[sfh_model]["key_labels"]:
            coverage = ev.coverage_curve(posteriors[label], truths[label].to_numpy(), [0.68, 0.95])
            choice_rows.append({"sfh_model": sfh_model, "embedding": EMBEDDING, "mode": mode, "label": label,
                                "coverage_68": coverage[0], "coverage_95": coverage[1],
                                "best_validation_loss": history["validation_loss"].min()})
    choice = pd.DataFrame(choice_rows)
    choice["coverage_error"] = np.abs(choice["coverage_68"] - 0.68) + np.abs(choice["coverage_95"] - 0.95)
    choices[sfh_model] = choice
axis.set(xlabel="epoch", ylabel="mean −log p", title="training history")
axis.legend()
by_model = pd.concat(choices.values(), ignore_index=True).groupby("sfh_model").agg(
    best_validation_loss=("best_validation_loss", "first"), mean_coverage_error=("coverage_error", "mean"))
display(by_model.round(4))
plt.show()
```

  **Cell 4** (markdown, id `66ee9249`, "## 2. Calibration, accuracy and information per parameter (LOCKED: synthetic test set)"):

```markdown
## 2. Calibration, accuracy and information per parameter (LOCKED: synthetic test set)

For each SFH model and redshift mode:
- the left panel is the coverage curve of the key labels (a calibrated posterior lies on the diagonal);
- the right panel is the histogram of the rank of the truth among the posterior samples (flat means calibrated, U-shaped means too narrow).

The table gives coverage, bias, NMAD and posterior/prior width for every key label (`LABELS[model]["key_labels"]`):
- the flow targets (9 for delayed-τ, 10 for DPL). `log_tau_bc` is evaluated as the decoded `tau_bc`, so its bias and NMAD are linear, not in dex like the log rows;
- the derived labels `mw_age` [Gyr], `f_young`, `f_mid`, `f_old`, `log_sfr` and `log_ssfr`.

The linear SFH parameters (`age`, `tau`, `turnover_lookback`, `alpha`, `beta`, `metallicity`) are not listed: they are monotone maps of the log targets, so their ranks are the same and their coverage nearly so.

Several derived labels have point masses in the prior:
- `f_young` and `f_old` are often exactly 0 or 1;
- `log_sfr` and `log_ssfr` sit on their floor when nothing formed in the last 100 Myr.

Ties make the inclusive coverage too high and pile the ranks at one edge, and a zero prior width gives a NaN width ratio. Read those rows as report-only.
```

  **Cell 5** (code, id `19fe68be`, "def mass_rescaled_prior(prior_table, catalogue, seed):"):

```python
def mass_rescaled_prior(prior_table, catalogue, seed):
    """Prior draws sorted by redshift, with the mass-like labels at observed brightness.

    Why: the training sets hold models at a fixed mass. An observed galaxy has the mass of its brightness
    (sm.scale_to_observed_r, as in training), so the prior width of log_mstar is ~1 dex, not 0.08.
    Inputs: prior_table: a model's training set; catalogue: DataFrame with r_auto and zml of DR3 galaxies;
            seed: for the brightness donors. Returns: DataFrame without rows lacking a donor, sorted by redshift.
    """
    offsets = sm.scale_to_observed_r(prior_table[sm.BANDS].to_numpy(dtype=float),
                                     prior_table["redshift"].to_numpy(dtype=float),
                                     catalogue["r_auto"].to_numpy(), catalogue["zml"].to_numpy(),
                                     np.random.default_rng(seed))[1]
    rescaled = sm.rescale_mass_labels(prior_table, offsets)
    return rescaled[np.isfinite(offsets)].sort_values("redshift").reset_index(drop=True)


def windowed_width_ratio(posterior, prior_table, label, true_redshift):
    """Posterior/prior width of a label, the prior taken within PRIOR_WINDOW_DZ of each galaxy's redshift."""
    prior_widths = ev.windowed_prior_widths(prior_table["redshift"].to_numpy(), prior_table[label].to_numpy(),
                                            true_redshift, PRIOR_WINDOW_DZ)
    return ev.width_ratio(posterior, None, prior_widths=prior_widths)


synthetic_test_summaries = {}      # sfh_model -> metrics on its synthetic test set (RUN_LOCKED_TEST only)
test_tables = {}
if RUN_LOCKED_TEST:
    dr3_brightness = pd.read_csv(SPLUS_FILE, usecols=["PROB_GAL", "r_auto", "zml"])
    dr3_brightness = dr3_brightness[dr3_brightness["PROB_GAL"] > MIN_PROB_GALAXY]
    for sfh_model in SFH_MODELS:
        test_tables[sfh_model] = pd.read_parquet(sfh_models.model_file(TEST_FILE, sfh_model))
        prior_table = mass_rescaled_prior(pd.read_parquet(sfh_models.model_file(TRAINING_FILE, sfh_model)),
                                          dr3_brightness, PRIOR_SEED)
        key_labels = LABELS[sfh_model]["key_labels"]
        summary_rows = []
        for mode, seed in TEST_SEEDS.items():
            observations = ml_data.validation_observations(test_tables[sfh_model], mode,
                                                           np.random.default_rng(seed), sfh_model)
            posteriors, truths = labelled_posteriors(runs[sfh_model], observations, N_POSTERIOR_SAMPLES,
                                                     SAMPLING_SEED, key_labels)
            figure, axes = plt.subplots(1, 2, figsize=(12, 5))
            for label in key_labels:
                truth = truths[label].to_numpy()
                coverage = ev.coverage_curve(posteriors[label], truth)
                bias, nmad = ev.bias_and_nmad(np.nanmedian(posteriors[label], axis=1), truth)
                ratio = windowed_width_ratio(posteriors[label], prior_table, label, truths["redshift"].to_numpy())
                summary_rows.append({"sfh_model": sfh_model, "section": "synthetic_test", "mode": mode,
                                     "label": label,
                                     "coverage_68": ev.coverage_curve(posteriors[label], truth, [0.68])[0],
                                     "coverage_95": ev.coverage_curve(posteriors[label], truth, [0.95])[0],
                                     "bias": bias, "nmad": nmad, "width_ratio": ratio})
                axes[0].plot(ev.NOMINAL_LEVELS, coverage, marker="o", label=label)
                axes[1].hist(ev.truth_ranks(posteriors[label], truth), bins=20, histtype="step", label=label)
            axes[0].plot([0, 1], [0, 1], color="black", lw=1)
            axes[0].set(xlabel="nominal coverage", ylabel="empirical coverage",
                        title=f"{sfh_model}, {mode} redshift")
            axes[1].set(xlabel="rank of truth", ylabel="galaxies")
            axes[0].legend(fontsize=6, ncol=2)
            plt.show()
        synthetic_test_summaries[sfh_model] = pd.DataFrame(summary_rows)
        display(synthetic_test_summaries[sfh_model].round(3))
else:
    print("synthetic test evaluation: skipped (locked)")
```

  **Cell 6** (code, id `d51a8bf8`, "if RUN_LOCKED_TEST:" with "Information per parameter"):

```python
if RUN_LOCKED_TEST:
    BAR_WIDTH = 0.4
    for sfh_model in SFH_MODELS:
        # Information per parameter: posterior width / prior width of each flow target
        # (1 = the colours say nothing about that parameter)
        summary = synthetic_test_summaries[sfh_model]
        flow_target_labels = LABELS[sfh_model]["flow_target_labels"]
        figure, axis = plt.subplots(figsize=(9, 4))
        for mode_index, mode in enumerate(TEST_SEEDS):
            mode_rows = summary[(summary["mode"] == mode) & summary["label"].isin(flow_target_labels)]
            width_ratios = mode_rows.set_index("label").loc[flow_target_labels, "width_ratio"]
            axis.bar(np.arange(len(flow_target_labels)) + (mode_index - 0.5) * BAR_WIDTH, width_ratios,
                     width=BAR_WIDTH, label=f"{mode} redshift")
        axis.axhline(1, color="grey", ls=":")
        axis.set_xticks(range(len(flow_target_labels)))
        axis.set_xticklabels(flow_target_labels, rotation=45, ha="right")
        axis.set(ylabel="posterior / prior 68% width",
                 title=f"{sfh_model}: information per parameter (1 = prior only)")
        axis.legend()
        plt.show()

        # Accuracy against brightness, spectroscopic mode
        spectroscopic_rng = np.random.default_rng(TEST_SEEDS["spectroscopic"])
        observations = ml_data.validation_observations(test_tables[sfh_model], "spectroscopic", spectroscopic_rng,
                                                       sfh_model)
        posteriors, truths = labelled_posteriors(runs[sfh_model], observations, N_POSTERIOR_SAMPLES,
                                                 SAMPLING_SEED, ["mw_age", "log_mstar"])
        r_magnitude = observations["flat"][observations["valid"], ml_data.SCALAR_INPUT_SLICE.start]
        r_quartile = pd.qcut(r_magnitude, 4, labels=["brightest", "second", "third", "faintest"])
        accuracy_rows = []
        for label in ("mw_age", "log_mstar"):
            medians = np.nanmedian(posteriors[label], axis=1)
            for quartile in r_quartile.categories:
                in_quartile = np.asarray(r_quartile == quartile)
                bias, nmad = ev.bias_and_nmad(medians[in_quartile], truths[label].to_numpy()[in_quartile])
                accuracy_rows.append({"sfh_model": sfh_model, "label": label, "r_quartile": quartile,
                                      "bias": bias, "nmad": nmad})
        display(pd.DataFrame(accuracy_rows).round(3))
else:
    print("information per parameter and accuracy against brightness: skipped (locked)")
```

  **Cell 7** (markdown, id `fe480605`, "## 3. Agreement with the Bagpipes spectral posteriors"):

```markdown
## 3. Agreement with the Bagpipes spectral posteriors

For each SFH model: ML posteriors from the 3″ aperture colours and spectroscopic redshift, against each reference galaxy's spectral posterior. The spectral fits are delayed-τ, so they hold the shared labels and the delayed-τ parameters; a DPL flow is compared on the labels the spectral fits hold. Development uses the `validation` half; the `test` half is LOCKED.
```

  **Cell 8** (code, id `e81a48b9`, "PHOTOMETRIC_PREFIX = "phot_""):

```python
reference = pd.read_csv(REFERENCE_FILE, dtype={"specobjid": str})
halves = ["validation", "test"] if RUN_LOCKED_TEST else ["validation"]
aperture_columns = []
aperture_error_columns = []
for band in sm.BANDS:
    aperture_columns.append(f"{band}_aper")
    aperture_error_columns.append(f"e_{band}_aper")


def available_labels(table, labels, prefix=""):
    """The labels whose {prefix}<label>_median column the table holds (a reference set may lack some)."""
    kept = []
    for label in labels:
        if f"{prefix}{label}_median" in table.columns:
            kept.append(label)
    return kept


def spectral_galaxy_types(galaxies):
    """'star-forming' where the spectral fit's log_ssfr_median exceeds PASSIVE_SSFR_LIMIT, else 'passive' (array)."""
    is_star_forming = galaxies["log_ssfr_median"].to_numpy(dtype=float) > PASSIVE_SSFR_LIMIT
    return np.where(is_star_forming, "star-forming", "passive")


def converged_fits(galaxies, prefix):
    """True where the galaxy's photometric fit with this column prefix exists and converged (array).

    A table without the {prefix}converged column (e.g. before the DPL fits are read) gives False everywhere.
    """
    column = f"{prefix}converged"
    if column not in galaxies.columns:
        return np.zeros(len(galaxies), dtype=bool)
    return galaxies[column].astype("boolean").fillna(False).to_numpy(dtype=bool)


# Why one sampling pass per model and half: the spectral, photometric and cross-model comparisons must use
# the same ML draws of each galaxy, so each model's posterior is drawn once for every galaxy with a
# spectral fit or a photometric fit of that model
has_spectral_fit = reference["fit_status"] == "fitted"
ml_posteriors_by_half = {}         # (sfh_model, half) -> (galaxies with an ML posterior, dict label -> (n, samples))
agreements = {}                    # sfh_model -> ML against the spectral fits, one row per (galaxy, label)
for sfh_model in SFH_MODELS:
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    has_photometric_fit = pd.Series(False, index=reference.index)
    if f"{prefix}mw_age_median" in reference.columns:
        has_photometric_fit = reference[f"{prefix}mw_age_median"].notna()
    spectral_labels = available_labels(reference, LABELS[sfh_model]["reference_labels"])
    n_targets = runs[sfh_model][1]["n_targets"]
    agreement_rows = []
    for half in halves:
        galaxies = reference[(reference["split"] == half) & (has_spectral_fit | has_photometric_fit)]
        galaxies = galaxies.reset_index(drop=True)
        if len(galaxies) == 0:
            print(f"{sfh_model}, {half}: no fitted galaxies yet")
            continue
        observations = ml_data.colour_inputs(galaxies[aperture_columns].to_numpy(),
                                             galaxies[aperture_error_columns].to_numpy(),
                                             galaxies["z"].to_numpy(), np.zeros(len(galaxies)))
        if not observations["valid"].any():
            print(f"{sfh_model}, {half}: no galaxy with a detected r band in the aperture photometry")
            continue
        observations["targets"] = np.zeros((len(galaxies), n_targets))   # unknown; not used
        posteriors, _ = labelled_posteriors(runs[sfh_model], observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED)
        kept = galaxies[observations["valid"]].reset_index(drop=True)
        ml_posteriors_by_half[(sfh_model, half)] = (kept, posteriors)

        # The spectral comparison: galaxies with a spectral fit only, on the labels the spectral fits hold
        spectral = (kept["fit_status"] == "fitted").to_numpy()
        spectral_galaxies = kept[spectral].reset_index(drop=True)
        spectral_posteriors = {label: posteriors[label][spectral] for label in spectral_labels}
        galaxy_types = dict(zip(spectral_galaxies["specobjid"], spectral_galaxy_types(spectral_galaxies)))
        galaxy_details = spectral_galaxies.set_index("specobjid")
        for row in ev.reference_comparison_rows(spectral_posteriors, spectral_galaxies, spectral_labels):
            details = galaxy_details.loc[row["specobjid"]]
            agreement_rows.append({"sfh_model": sfh_model, "half": half, **row,
                                   "galaxy_type": galaxy_types[row["specobjid"]],
                                   "fibre": details["fiber_diameter"],
                                   "calibrated_fitter": details["fitter_calibrated"]})
    agreement = pd.DataFrame(agreement_rows)
    agreements[sfh_model] = agreement
    if len(agreement) > 0:
        print(f"{sfh_model}: ML against the spectral Bagpipes posteriors")
        display(agreement.groupby(["half", "label"]).agg(median_offset=("offset_in_widths", "median"),
                                                          inside_68=("inside_68", "mean"),
                                                          inside_95=("inside_95", "mean"),
                                                          galaxies=("specobjid", "count")).round(3))
        display(agreement.groupby(["half", "label", "galaxy_type"])["inside_95"].mean().unstack().round(3))
```

  **Cell 9** (markdown, id `e05c91f8`, "### 3b. The "works" test: photometric and spectral Bagpipes"):

```markdown
### 3b. The "works" test: photometric and spectral Bagpipes

For each SFH model, the ML is compared with two Bagpipes references, and each answers a different question:

- **Photometric Bagpipes (emulator check).** These are the model's own photometric fits (notebook 2, section 6b for delayed-τ and 6c for DPL): `phot_<label>` for delayed-τ, `phot_dblplaw_<label>` for DPL (`sm.PHOTOMETRIC_PREFIXES`). They fit the same 3″ photometry the ML sees, under that model's training prior and the simulator's likelihood, so a faithful emulator must agree with them.
- **Spectral Bagpipes (physical check).** The spectral fits use the SDSS spectra; agreement tells whether the colours recover what the spectra say. They are delayed-τ fits: they hold the shared labels and the delayed-τ parameters, not the DPL ones. Section 5.3 compares both models with them on the shared labels.

A label agrees with a reference set when two conditions hold:
- at least `ev.MIN_INSIDE_95` (80%) of the reference medians lie inside the ML 95% interval;
- the median offset is within ±`ev.MAX_ABS_MEDIAN_OFFSET` combined 1σ widths.

The judged labels are the model's calibration labels. For delayed-τ they are `log_age`, `log_tau`, `log_metallicity` and `log_mstar`. For DPL they are `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity` and `log_mstar`.

A reference set judges only the judged labels it holds. The delayed-τ spectral fits have no `log_turnover_lookback`, `log_alpha` or `log_beta`, so against the spectra the DPL flow is judged on `log_metallicity` and `log_mstar` only. Its other judged labels are listed with the status "not in the reference" and stay out of the `all` row. The photometric set judges the full list.

The other shared labels (`sfh_models.SHARED_LABELS`) are reported, not judged. A reported label a reference set does not hold has an empty row. Only the `validation` half is used. The verdicts are printed side by side, and the notebook draws no conclusion beyond them.
```

  **Cell 10** (code, id `74230b0e`, "REFERENCE_SET_NAMES = {"):

```python
REFERENCE_SET_NAMES = {"photometric": "photometric Bagpipes (emulator check)",
                       "spectral": "spectral Bagpipes (physical check)"}


NOT_IN_REFERENCE = "not in the reference"   # status of a judged label the reference set does not hold


def verdict_table(rows, judged_labels, reported_labels, not_in_reference=()):
    """One verdict table with a status column.

    Rows: agreement_verdict over the judged labels (status "judged", with its 'all' row, which covers only
    them); then the reported labels (status "reported", without their 'all' row); then the labels the
    reference set cannot judge (status NOT_IN_REFERENCE, no galaxies, agrees left empty).
    """
    tables = [ev.agreement_verdict(rows, judged_labels).assign(status="judged"),
              ev.agreement_verdict(rows, reported_labels).iloc[:-1].assign(status="reported")]
    if not_in_reference:
        missing_rows = []
        for label in not_in_reference:
            missing_rows.append({"label": label, "n_galaxies": 0, "inside_95": np.nan, "median_offset": np.nan,
                                 "agrees": None, "status": NOT_IN_REFERENCE})
        tables.append(pd.DataFrame(missing_rows))
    return pd.concat(tables, ignore_index=True)


n_validation = int((reference["split"] == "validation").sum())
verdicts = {}                      # (sfh_model, reference set) -> verdict table, validation half
for sfh_model in SFH_MODELS:
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    judged_labels = LABELS[sfh_model]["judged_labels"]
    reported_labels = LABELS[sfh_model]["reported_labels"]
    fitter_calibration_columns = {"photometric": f"{prefix}fitter_calibrated", "spectral": "fitter_calibrated"}

    # 1. The validation galaxies with a converged photometric fit of this model, with the ML draws of cell 8
    # (the test half stays locked)
    photometric_labels = available_labels(reference, LABELS[sfh_model]["reference_labels"], prefix)
    photometric_agreement_rows = []
    n_photometric_fits = 0
    n_photometric_not_converged = 0
    if (sfh_model, "validation") in ml_posteriors_by_half and f"{prefix}mw_age_median" in reference.columns:
        validation_galaxies, validation_posteriors = ml_posteriors_by_half[(sfh_model, "validation")]
        photometric = validation_galaxies[f"{prefix}mw_age_median"].notna().to_numpy()
        n_photometric_fits = int(photometric.sum())
        if n_photometric_fits > 0:
            converged = converged_fits(validation_galaxies, prefix)
            n_photometric_not_converged = int((photometric & ~converged).sum())
            compared = photometric & converged
            photometric_galaxies = validation_galaxies[compared].reset_index(drop=True)
            photometric_posteriors = {label: validation_posteriors[label][compared] for label in photometric_labels}
            photometric_agreement_rows = ev.reference_comparison_rows(photometric_posteriors, photometric_galaxies,
                                                                      photometric_labels, prefix=prefix)
    n_photometric_compared = len({row["specobjid"] for row in photometric_agreement_rows})
    print(f"{sfh_model}: photometric fits of {n_photometric_compared} of {n_validation} validation galaxies "
          f"compared ({n_photometric_fits} with a fit and an ML posterior, {n_photometric_not_converged} not "
          "converged and left out)")
    if n_photometric_fits == 0:
        print(f"  no {sfh_model} photometric Bagpipes fits in {REFERENCE_FILE} yet: run notebook 2 after the fits")
    missing_labels = [label for label in LABELS[sfh_model]["reference_labels"] if label not in photometric_labels]
    if missing_labels:
        print(f"  {REFERENCE_FILE} has no {prefix}<label> columns for {missing_labels}: their rows stay empty")

    # 2. The rows of both reference sets, validation half only
    spectral_validation_rows = agreements[sfh_model]
    if len(spectral_validation_rows) > 0:
        spectral_validation_rows = spectral_validation_rows[spectral_validation_rows["half"] == "validation"]
    comparison_rows = {"photometric": photometric_agreement_rows, "spectral": spectral_validation_rows}

    # 3. The judged labels each reference set holds: the spectral fits are delayed-tau, so they cannot judge
    # the DPL shape labels (log_turnover_lookback, log_alpha, log_beta); the photometric set judges every one
    judged_by_set = {"photometric": judged_labels, "spectral": available_labels(reference, judged_labels)}

    # 4. Both verdicts side by side. Why the fitter flag: a reference whose fitter failed its calibration test
    # (or whose calibration mocks are not all fitted) is a weaker reference, so it is shown next to the verdict
    for reference_set, rows in comparison_rows.items():
        set_judged_labels = judged_by_set[reference_set]
        not_in_reference = [label for label in judged_labels if label not in set_judged_labels]
        verdict = verdict_table(rows, set_judged_labels, reported_labels, not_in_reference)
        verdicts[(sfh_model, reference_set)] = verdict
        fitter_flags = pd.Series(dtype=object)
        if fitter_calibration_columns[reference_set] in reference.columns:
            fitter_flags = reference[fitter_calibration_columns[reference_set]].dropna()
        fitter_calibrated = "unknown (calibration mocks not all fitted)"
        if len(fitter_flags) > 0:
            fitter_calibrated = fitter_flags.iloc[0]
        print(f"{sfh_model}: agreement with {REFERENCE_SET_NAMES[reference_set]} (validation half; "
              f"{', '.join(reported_labels)} reported, not judged); fitter calibrated: {fitter_calibrated}")
        display(verdict.round(3))
        if not_in_reference:
            print(f"  {', '.join(not_in_reference)}: {NOT_IN_REFERENCE}, so not judged here")
        judged_rows = verdict[verdict["status"] == "judged"]
        if judged_rows["n_galaxies"].sum() == 0:
            print(f"  no galaxy for the judged labels {set_judged_labels} in this reference set")
        else:
            all_row = verdict[verdict["label"] == ev.ALL_LABELS_ROW]
            print(f"  the ML agrees on {', '.join(set_judged_labels)}: {bool(all_row['agrees'].iloc[0])}")
```

  **Cell 11** (markdown, id `84911f5f`, "## 4. S-PLUS DR3 catalogue"):

```markdown
## 4. S-PLUS DR3 catalogues

One catalogue per SFH model (`ml_sfh_catalogue.parquet`, `ml_sfh_catalogue_dblplaw.parquet`). Each holds photometric-z posteriors for every DR3 galaxy that passes that model's coverage test (against its own training set); the others are flagged.
```

  **Cell 12** (code, id `46b23dbe`, "splus = pd.read_csv(SPLUS_FILE)"):

```python
splus = pd.read_csv(SPLUS_FILE)
splus = splus[splus["PROB_GAL"] > MIN_PROB_GALAXY].reset_index(drop=True)
magnitudes = splus[[f"{band}_auto" for band in sm.BANDS]].to_numpy()
errors = splus[[f"e_{band}_auto" for band in sm.BANDS]].to_numpy()
detected = sm.detected_bands(magnitudes, errors)
zml = splus["zml"].to_numpy()
sigma_z = sm.SIGMA_Z * (1 + zml)
inputs = ml_data.colour_inputs(magnitudes, errors, zml, sigma_z)
blocks = np.array_split(np.arange(len(splus)), 4 * N_JOBS)
CATALOGUE_CHUNK = 2000             # galaxies per chunk: all samples at once would need ~4 GB
SUMMARY_PERCENTILES = ((16, "p16"), (50, "median"), (84, "p84"))

catalogues = {}                    # sfh_model -> DR3 catalogue (also written to its parquet file)
for sfh_model in SFH_MODELS:
    # 1. Coverage test against this model's training set (colour-based through its free offset), in parallel
    training_models = pd.read_parquet(sfh_models.model_file(TRAINING_FILE, sfh_model),
                                      columns=["redshift", *sm.BANDS])
    tasks = []
    for block in blocks:
        tasks.append(delayed(sm.coverage_test)(magnitudes[block], errors[block], detected[block], zml[block],
                                               sm.redshift_window(zml[block]), training_models[sm.BANDS].to_numpy(),
                                               training_models["redshift"].to_numpy()))
    coverage = pd.concat(Parallel(n_jobs=N_JOBS)(tasks), ignore_index=True)

    # 2. Posterior summaries (16th, 50th, 84th percentile of every decoded label) for covered galaxies with a
    # detected r band
    predictable = inputs["valid"] & coverage["covered"].to_numpy()
    estimator, config, standardisers = runs[sfh_model]
    predictable_indices = np.flatnonzero(predictable)
    summaries = {}                 # column name -> values for every DR3 galaxy
    start = time.time()
    decoding_seconds = 0.0
    for chunk_start in range(0, len(predictable_indices), CATALOGUE_CHUNK):
        chunk = predictable_indices[chunk_start:chunk_start + CATALOGUE_CHUNK]
        theta = ml_networks.posterior_samples(estimator, standardisers, inputs["flat"][chunk],
                                              inputs["banded"][chunk], N_CATALOGUE_SAMPLES,
                                              SAMPLING_SEED + chunk_start, DEVICE)
        decoding_start = time.time()
        decoded = ml_data.decode_samples(theta.reshape(-1, config["n_targets"]), sfh_model, n_jobs=N_DECODE_JOBS)
        decoding_seconds += time.time() - decoding_start
        for column in decoded.columns:
            values = decoded[column].to_numpy().reshape(len(chunk), N_CATALOGUE_SAMPLES)
            for percentile, suffix in SUMMARY_PERCENTILES:
                name = f"{column}_{suffix}"
                if name not in summaries:
                    summaries[name] = np.full(len(splus), np.nan)
                summaries[name][chunk] = np.nanpercentile(values, percentile, axis=1)
    seconds = time.time() - start

    catalogue = splus[["RA", "DEC", "zml"]].copy()
    catalogue["no_r_band"] = ~inputs["valid"]
    catalogue["out_of_domain"] = inputs["valid"] & ~coverage["covered"].to_numpy()
    # Why one concat: adding ~60 columns one by one fragments the frame (a pandas PerformanceWarning)
    catalogue = pd.concat([catalogue, pd.DataFrame(summaries, index=catalogue.index)], axis=1)
    catalogue_file = sfh_models.model_file(CATALOGUE_FILE, sfh_model)
    catalogue.to_parquet(catalogue_file, index=False)
    catalogues[sfh_model] = catalogue
    n_predicted = len(predictable_indices)
    print(f"{sfh_model}: {n_predicted} of {len(splus)} galaxies predicted "
          f"({catalogue['out_of_domain'].mean():.1%} out of domain, {catalogue['no_r_band'].sum()} without r); "
          f"{n_predicted / seconds:.0f} galaxies/s on {DEVICE}; wrote {catalogue_file}")
    n_decoded = n_predicted * N_CATALOGUE_SAMPLES
    print(f"  decoding took {decoding_seconds / 60:.1f} of {seconds / 60:.1f} min "
          f"({n_decoded / max(decoding_seconds, 1e-9):.0f} samples/s, {N_DECODE_JOBS} threads)")
```

  **Cell 13** (code, id `2e021085`, "summary_tables = [summary] if RUN_LOCKED_TEST else []"):

```python
for sfh_model in SFH_MODELS:
    summary_tables = []
    if RUN_LOCKED_TEST:
        summary_tables.append(synthetic_test_summaries[sfh_model])
    summary_tables.append(choices[sfh_model].assign(section="validation_choice"))
    for reference_set in ("photometric", "spectral"):
        verdict = verdicts[(sfh_model, reference_set)]
        summary_tables.append(verdict.assign(section=f"agreement_{reference_set}", sfh_model=sfh_model))
    summary_file = sfh_models.model_file(SUMMARY_FILE, sfh_model)
    pd.concat(summary_tables, ignore_index=True).to_csv(summary_file, index=False)
    print(f"wrote {summary_file}")
```

- [ ] **Step 2: Insert the new section 5.** Use NotebookEdit with `edit_mode: insert`, `cell_type` as given. Insert the first new cell after id `2e021085`, and each following cell after the one inserted before it, in this order.

  **New cell A** (markdown):

```markdown
## 5. Delayed-τ against double power law

The two SFH models side by side, in four parts (spec 2026-10-28, Part 4):
1. each model's ML against its own photometric Bagpipes (the "works" tests of section 3b, side by side);
2. delayed-τ against DPL on the shared labels:
   - the ML of each model, on the validation galaxies and on DR3;
   - the photometric Bagpipes fits of each model, on the validation galaxies;
3. four estimators against the spectral Bagpipes fits, passive and star-forming galaxies apart;
4. the Bayesian evidence of each model, galaxy by galaxy and on the cross-fitted mocks.

The shared labels are model-independent: `log_mstar`, `mw_age`, `log_metallicity`, `av_ism`, `log_ssfr`, `f_young`, `f_mid`, `f_old` (`sfh_models.SHARED_LABELS`).

Every difference is delayed-τ minus DPL. Only the validation half is used, plus the test half with `RUN_LOCKED_TEST`. The photometric Bagpipes and evidence comparisons use the test half only when notebook 2 has read its photometric fits (`EXPORT_LOCKED_TEST_JOBS = True` there); otherwise they cover the validation half, and the cell below says so. The section is report-only: the evidence is reported, not used to choose a model per galaxy.

**Caveat on the DPL shape parameters.** In the DPL, τ is the turnover time of the double power law and `turnover_lookback` the lookback time to it; the SFR peaks there only when α = β (in general the maximum lies at t = τ (β/α)^{1/(α+β)}), so `log_turnover_lookback` is not the age of peak star formation and is not compared with the delayed-τ age.
```

  **New cell B** (markdown):

```markdown
### 5.1 ML against its own photometric Bagpipes

The verdicts of section 3b, label × model. Each model's `all` row covers its own calibration labels. The shared labels are reported, and a label one model does not have is empty in its column.
```

  **New cell C** (code):

```python
COMPARED_MODELS = ["delayed", "dblplaw"]   # ln B and every difference below: delayed-tau minus DPL
SHARED_LABELS = sfh_models.SHARED_LABELS
compare_models = set(COMPARED_MODELS) <= set(SFH_MODELS)
comparison_tables = []             # long tables with a "section" column, written to COMPARISON_FILE at the end
if not compare_models:
    print(f"Section 5 needs both models in SFH_MODELS (now {SFH_MODELS}): skipped")
else:
    own_photometric = {}
    for sfh_model in COMPARED_MODELS:
        own_photometric[sfh_model] = verdicts[(sfh_model, "photometric")]
        comparison_tables.append(own_photometric[sfh_model].assign(section="own_photometric",
                                                                   sfh_model=sfh_model))
    print("Each model's ML against its own photometric Bagpipes (validation half): n_galaxies, inside_95, "
          "median offset [combined 1-sigma widths], agrees")
    display(ev.side_by_side_verdicts(own_photometric).round(3))

    # The halves of the photometric Bagpipes and evidence comparisons (5.2-5.4). The test half has photometric
    # fits only when notebook 2 exported and read them (EXPORT_LOCKED_TEST_JOBS); the validation half always
    # takes part, so its tables exist even before the fits are read
    photometric_halves = ["validation"]
    if RUN_LOCKED_TEST:
        test_galaxies = reference[reference["split"] == "test"]
        test_has_both_fits = (converged_fits(test_galaxies, sm.PHOTOMETRIC_PREFIXES["delayed"])
                              & converged_fits(test_galaxies, sm.PHOTOMETRIC_PREFIXES["dblplaw"]))
        if test_has_both_fits.any():
            photometric_halves.append("test")
            print(f"Photometric Bagpipes and evidence comparisons: validation and test halves "
                  f"({int(test_has_both_fits.sum())} test galaxies with both fits converged)")
        else:
            print("Photometric Bagpipes and evidence comparisons: validation half only. The reference set has no "
                  "test-half photometric fits of both models: set EXPORT_LOCKED_TEST_JOBS = True in notebook 2, "
                  "fit those jobs and read them, then run this notebook again")
    else:
        print("Photometric Bagpipes and evidence comparisons: validation half (the test half is locked)")
```

  **New cell D** (markdown):

```markdown
### 5.2 Delayed-τ against DPL on the shared labels

For each pair of estimates, the cells below compute three statistics per label, over the galaxies both estimates have:
- the median of the per-galaxy difference of medians (delayed-τ minus DPL);
- its scatter: the NMAD around that median, 1.4826 × the median absolute deviation;
- the share of galaxies whose two 16–84% intervals overlap.

Pairs:
- **ML:** validation galaxies (the draws of cell 8), and DR3 galaxies predicted by both catalogues;
- **photometric Bagpipes:** galaxies where both photometric fits converged, in the halves of `photometric_halves`.
```

  **New cell E** (code):

```python
if compare_models:
    model_pairs = []
    validation_ml_pairs = None
    for half in halves:
        # 1. The ML of each model, on the galaxies with both posteriors (cell 8)
        ml_summaries = []
        for sfh_model in COMPARED_MODELS:
            if (sfh_model, half) not in ml_posteriors_by_half:
                continue
            galaxies, posteriors = ml_posteriors_by_half[(sfh_model, half)]
            summary = pd.DataFrame({"specobjid": galaxies["specobjid"].to_numpy()})
            for label in SHARED_LABELS:
                p16, median, p84 = ev.percentile_summary(posteriors[label])
                summary[f"ml_{sfh_model}_{label}_p16"] = p16
                summary[f"ml_{sfh_model}_{label}_median"] = median
                summary[f"ml_{sfh_model}_{label}_p84"] = p84
            ml_summaries.append(summary)
        if len(ml_summaries) == len(COMPARED_MODELS):
            ml_pairs = ml_summaries[0].merge(ml_summaries[1], on="specobjid")
            rows = ev.model_comparison_rows(ml_pairs, SHARED_LABELS, "ml_delayed_", "ml_dblplaw_")
            model_pairs.append(rows.assign(pair="ML", sample=half))
            if half == "validation":
                validation_ml_pairs = ml_pairs

    # 2. The photometric Bagpipes fits of each model, on the galaxies where both converged
    for half in photometric_halves:
        half_galaxies = reference[reference["split"] == half].reset_index(drop=True)
        both_converged = (converged_fits(half_galaxies, sm.PHOTOMETRIC_PREFIXES["delayed"])
                          & converged_fits(half_galaxies, sm.PHOTOMETRIC_PREFIXES["dblplaw"]))
        rows = ev.model_comparison_rows(half_galaxies[both_converged], SHARED_LABELS,
                                        sm.PHOTOMETRIC_PREFIXES["delayed"], sm.PHOTOMETRIC_PREFIXES["dblplaw"])
        model_pairs.append(rows.assign(pair="photometric Bagpipes", sample=half))

    # 3. The ML of each model on DR3: both catalogues have the same rows (cell 12); compared where both predicted
    dr3_pairs = pd.concat([catalogues["delayed"].add_prefix("ml_delayed_"),
                           catalogues["dblplaw"].add_prefix("ml_dblplaw_")], axis=1)
    rows = ev.model_comparison_rows(dr3_pairs, SHARED_LABELS, "ml_delayed_", "ml_dblplaw_")
    model_pairs.append(rows.assign(pair="ML", sample="DR3"))

    model_pairs = pd.concat(model_pairs, ignore_index=True)
    comparison_tables.append(model_pairs.assign(section="model_pairs"))
    print("Delayed-tau minus DPL on the shared labels: median difference, scatter (NMAD around it) and the "
          "share of galaxies whose 16-84% intervals overlap")
    display(model_pairs.set_index(["pair", "sample", "label"]).round(3))

    # 4. ML medians of the two models on the validation galaxies, one panel per shared label
    if validation_ml_pairs is not None:
        figure, axes = plt.subplots(2, 4, figsize=(16, 7))
        for axis, label in zip(axes.flat, SHARED_LABELS):
            delayed_median = validation_ml_pairs[f"ml_delayed_{label}_median"]
            dblplaw_median = validation_ml_pairs[f"ml_dblplaw_{label}_median"]
            axis.scatter(delayed_median, dblplaw_median, s=4)
            low = min(delayed_median.min(), dblplaw_median.min())
            high = max(delayed_median.max(), dblplaw_median.max())
            axis.plot([low, high], [low, high], color="black", lw=1)
            axis.set(xlabel=f"delayed-τ {label}", ylabel=f"DPL {label}")
        figure.suptitle("ML posterior medians, validation galaxies")
        plt.tight_layout()
        plt.show()
```

  **New cell F** (markdown):

```markdown
### 5.3 Against the spectra

The "works" test of four estimators against the spectral Bagpipes fits, on the shared labels: the ML of each model and the photometric Bagpipes fits of each model.
- Passive and star-forming galaxies are judged apart, split at a spectral log sSFR of `PASSIVE_SSFR_LIMIT` (as in section 3).
- The photometric Bagpipes fits keep only their 16th, 50th and 84th percentiles. Their 95% interval is therefore extended from these as a Gaussian, each side by its own width (`ev.compare_summaries_with_reference`).

**Caveat: the spectral reference assumes a delayed-τ SFH.** The spectral fits and the spectral calibration mocks are delayed-τ only, so the reference `mw_age` and mass fractions `f_young`, `f_mid`, `f_old` carry the delayed-τ shape. A DPL estimate that disagrees on them may describe the same galaxy with another history rather than be wrong. `log_mstar`, `log_metallicity`, `av_ism` and `log_ssfr` depend less on the shape.
```

  **New cell G** (code):

```python
if compare_models:
    ML_ESTIMATORS = []
    PHOTOMETRIC_ESTIMATORS = []
    for sfh_model in COMPARED_MODELS:
        ML_ESTIMATORS.append(f"ML {sfh_model}")
        PHOTOMETRIC_ESTIMATORS.append(f"photometric Bagpipes {sfh_model}")

    # 1. One row per (estimator, galaxy, label), with the galaxy type
    spectral_tables = []
    for half in halves:
        for sfh_model in COMPARED_MODELS:
            agreement = agreements[sfh_model]
            if len(agreement) == 0:
                continue
            ml_rows = agreement[(agreement["half"] == half) & agreement["label"].isin(SHARED_LABELS)]
            spectral_tables.append(ml_rows.assign(estimator=f"ML {sfh_model}"))
        if half not in photometric_halves:
            continue
        spectral_galaxies = reference[(reference["split"] == half) & (reference["fit_status"] == "fitted")]
        spectral_galaxies = spectral_galaxies.reset_index(drop=True)
        galaxy_types = dict(zip(spectral_galaxies["specobjid"], spectral_galaxy_types(spectral_galaxies)))
        for sfh_model in COMPARED_MODELS:
            prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
            fitted = spectral_galaxies[converged_fits(spectral_galaxies, prefix)].reset_index(drop=True)
            labels = available_labels(fitted, SHARED_LABELS, prefix)
            photometric_rows = pd.DataFrame(ev.summary_comparison_rows(fitted, labels, prefix))
            if len(photometric_rows) == 0:
                continue
            photometric_rows["half"] = half
            photometric_rows["galaxy_type"] = photometric_rows["specobjid"].map(galaxy_types)
            spectral_tables.append(photometric_rows.assign(estimator=f"photometric Bagpipes {sfh_model}"))
    spectral_rows = pd.DataFrame(columns=["half", "galaxy_type", "estimator", "label", "offset_in_widths",
                                          "inside_95"])
    if spectral_tables:
        spectral_rows = pd.concat(spectral_tables, ignore_index=True)

    # 2. The verdict of each estimator, passive and star-forming apart; all shared labels are judged here. A half
    # without photometric fits (see photometric_halves) judges the ML estimators only
    spectral_verdicts = []
    for half in halves:
        half_estimators = list(ML_ESTIMATORS)
        if half in photometric_halves:
            half_estimators.extend(PHOTOMETRIC_ESTIMATORS)
        for galaxy_type in ("passive", "star-forming"):
            by_estimator = {}
            for estimator in half_estimators:
                rows = spectral_rows[(spectral_rows["half"] == half) & (spectral_rows["galaxy_type"] == galaxy_type)
                                     & (spectral_rows["estimator"] == estimator)]
                by_estimator[estimator] = ev.agreement_verdict(rows, SHARED_LABELS)
                spectral_verdicts.append(by_estimator[estimator].assign(half=half, galaxy_type=galaxy_type,
                                                                        estimator=estimator))
            print(f"Against the spectral Bagpipes fits (delayed-tau, see the caveat): {half} half, "
                  f"{galaxy_type} galaxies")
            display(ev.side_by_side_verdicts(by_estimator).round(3))
    comparison_tables.append(pd.concat(spectral_verdicts, ignore_index=True).assign(section="against_spectra"))
```

  **New cell H** (markdown):

```markdown
### 5.4 Bayesian evidence

ln B = `phot_log_evidence` − `phot_dblplaw_log_evidence` for each galaxy whose two photometric fits converged, in the halves of `photometric_halves`; positive favours delayed-τ. Both fits see the same photometry through the same likelihood, so ln B compares the two SFH priors on these data. The categories on |ln B| are: below 1 inconclusive, 1–3 positive, above 3 strong (`ev.evidence_strength`).

The statistical uncertainty of each MultiNest ln Z is √(H / n_live):
- H is the information, the prior-to-posterior Kullback–Leibler divergence;
- n_live = `sm.N_LIVE_POINTS`.

The two fits are independent, so their errors add in quadrature. The reference CSV carries them as `phot_log_evidence_error` and `phot_dblplaw_log_evidence_error` when the fits stored them; older fits have none.

**Cross-fitted mocks.** The 33 delayed-τ photometric mocks were also fitted with DPL, and the 33 DPL mocks with delayed-τ (`cross_calib_<k>`, notebook 2 section 6c). Notebook 2 writes the ln Z of every mock under both models to `photometric_cross_fits.csv`, which the last code cell reads. The true model is known, so the share of mocks whose ln B prefers it shows whether the colours can separate the two models at all. The median ln B uncertainty of the mocks is printed with it.
```

  **New cell I** (code):

```python
if compare_models:
    # 1. ln B and its uncertainty per galaxy, where both photometric fits converged
    evidence_tables = []
    for half in photometric_halves:
        half_galaxies = reference[reference["split"] == half].reset_index(drop=True)
        both_converged = (converged_fits(half_galaxies, sm.PHOTOMETRIC_PREFIXES["delayed"])
                          & converged_fits(half_galaxies, sm.PHOTOMETRIC_PREFIXES["dblplaw"]))
        galaxies = half_galaxies[both_converged].reset_index(drop=True)
        table = pd.DataFrame({"specobjid": galaxies["specobjid"], "half": half,
                              "galaxy_type": spectral_galaxy_types(galaxies)})
        table["log_bayes"] = ev.log_bayes_factor(galaxies["phot_log_evidence"], galaxies["phot_dblplaw_log_evidence"])
        table["log_bayes_error"] = np.nan
        if {"phot_log_evidence_error", "phot_dblplaw_log_evidence_error"} <= set(galaxies.columns):
            table["log_bayes_error"] = ev.log_bayes_factor_error(galaxies["phot_log_evidence_error"],
                                                                 galaxies["phot_dblplaw_log_evidence_error"])
        evidence_tables.append(table)
    evidence = pd.concat(evidence_tables, ignore_index=True)

    # 2. The distribution of ln B and the share of each category, overall and by galaxy type
    distribution_rows = []
    category_tables = []
    for half in photometric_halves:
        in_half = evidence["half"] == half
        groups = {"all": in_half, "passive": in_half & (evidence["galaxy_type"] == "passive"),
                  "star-forming": in_half & (evidence["galaxy_type"] == "star-forming")}
        for group, in_group in groups.items():
            log_bayes = evidence.loc[in_group, "log_bayes"].to_numpy(dtype=float)
            finite = log_bayes[np.isfinite(log_bayes)]
            percentiles = np.full(3, np.nan)
            if len(finite) > 0:
                percentiles = np.percentile(finite, [16, 50, 84])
            distribution_rows.append({"half": half, "group": group, "n_galaxies": len(finite),
                                      "ln_B_p16": percentiles[0], "ln_B_median": percentiles[1],
                                      "ln_B_p84": percentiles[2]})
            category_tables.append(ev.evidence_category_fractions(log_bayes).assign(half=half, group=group))
    distribution = pd.DataFrame(distribution_rows)
    categories = pd.concat(category_tables, ignore_index=True)
    comparison_tables.append(distribution.assign(section="evidence_distribution"))
    comparison_tables.append(categories.assign(section="evidence_categories"))
    print("ln B = ln Z(delayed-tau) - ln Z(DPL) per galaxy (positive favours delayed-tau)")
    display(distribution.round(2))
    display(categories.set_index(["half", "group", "favoured_model", "strength"])["fraction"]
            .unstack(["favoured_model", "strength"]).round(3))

    figure, axis = plt.subplots(figsize=(7, 4))
    validation_evidence = evidence[evidence["half"] == "validation"]
    for group in ("passive", "star-forming"):
        values = validation_evidence.loc[validation_evidence["galaxy_type"] == group, "log_bayes"].dropna()
        axis.hist(values, bins=30, histtype="step", label=f"{group} ({len(values)})")
    for edge in (-ev.STRONG_LOG_BAYES, -ev.INCONCLUSIVE_LOG_BAYES, ev.INCONCLUSIVE_LOG_BAYES, ev.STRONG_LOG_BAYES):
        axis.axvline(edge, color="grey", ls=":")
    axis.set(xlabel="ln B = ln Z(delayed-τ) − ln Z(DPL)", ylabel="galaxies",
             title="Evidence per validation galaxy (positive favours delayed-τ)")
    axis.legend()
    plt.show()

    # 3. The evidence uncertainty: stored per fit when MultiNest gave it, otherwise its scale sqrt(H / n_live)
    with_error = evidence["log_bayes_error"].notna()
    if with_error.any():
        below_error = (evidence.loc[with_error, "log_bayes"].abs() < evidence.loc[with_error, "log_bayes_error"])
        print(f"ln B uncertainty: median {evidence.loc[with_error, 'log_bayes_error'].median():.2f} over "
              f"{int(with_error.sum())} galaxies; {below_error.mean():.0%} of them have |ln B| below it")
    else:
        print(f"No stored evidence uncertainty in {REFERENCE_FILE}. Its scale for n_live = {sm.N_LIVE_POINTS}:")
        for information in (5.0, 10.0, 20.0):
            log_evidence_error = ev.nested_sampling_log_evidence_error(information, sm.N_LIVE_POINTS)
            log_bayes_error = ev.log_bayes_factor_error(log_evidence_error, log_evidence_error)
            print(f"  H = {information:.0f} nats: sigma(ln Z) = {log_evidence_error:.2f}, "
                  f"sigma(ln B) = {log_bayes_error:.2f}")
```

  **New cell J** (code):

```python
if compare_models:
    if not os.path.exists(CROSS_FITS_FILE):
        print(f"No {CROSS_FITS_FILE} yet: run notebook 2 (section 6c) after the cross-fits. "
              "True-model recovery skipped")
    else:
        # 1. The mocks whose own fit and cross-fit both converged (notebook 2 checked that every fit is current)
        cross_fits = pd.read_csv(CROSS_FITS_FILE)
        both_converged = (converged_fits(cross_fits, sm.PHOTOMETRIC_PREFIXES["delayed"])
                          & converged_fits(cross_fits, sm.PHOTOMETRIC_PREFIXES["dblplaw"]))
        # ev.true_model_recovery reads log_evidence_<model>; the CSV names them <prefix>log_evidence
        evidence_columns = {}
        for sfh_model, prefix in sm.PHOTOMETRIC_PREFIXES.items():
            evidence_columns[f"{prefix}log_evidence"] = f"log_evidence_{sfh_model}"
        converged_cross_fits = cross_fits[both_converged].rename(columns=evidence_columns)
        print(f"Cross-fitted mocks: {len(cross_fits)} in {CROSS_FITS_FILE}, {int(both_converged.sum())} with both "
              f"fits converged; the other {len(cross_fits) - int(both_converged.sum())} are left out")

        # 2. How often ln B prefers the true model
        recovery = ev.true_model_recovery(converged_cross_fits)
        comparison_tables.append(recovery.assign(section="true_model_recovery"))
        print("How often ln B prefers the true model (fractions of the mocks of each true model)")
        display(recovery.round(3))

        # 3. The ln B uncertainty of the mocks, from the stored ln Z errors of both fits
        mock_log_bayes_error = ev.log_bayes_factor_error(
            converged_cross_fits[f"{sm.PHOTOMETRIC_PREFIXES['delayed']}log_evidence_error"],
            converged_cross_fits[f"{sm.PHOTOMETRIC_PREFIXES['dblplaw']}log_evidence_error"])
        if np.isfinite(mock_log_bayes_error).any():
            print(f"ln B uncertainty of the mocks: median {np.nanmedian(mock_log_bayes_error):.2f} nats over "
                  f"{int(np.isfinite(mock_log_bayes_error).sum())} mocks")
        else:
            print("No stored ln Z uncertainty for the mocks (fits made before MultiNest's error was stored)")
```

  **New cell K** (code):

```python
if compare_models:
    pd.concat(comparison_tables, ignore_index=True).to_csv(COMPARISON_FILE, index=False)
    print(f"wrote {COMPARISON_FILE}")
```

- [ ] **Step 3: Clear the stale outputs, validate and compile; never execute.**
  Run:
```bash
cd Synthetic_domain && ../.venv/bin/python - <<'EOF'
import nbformat

notebook = nbformat.read("ml_model.ipynb", as_version=4)
for cell in notebook.cells:
    if cell.cell_type == "code":
        cell.outputs = []
        cell.execution_count = None
nbformat.validate(notebook)
nbformat.write(notebook, "ml_model.ipynb")
print(len(notebook.cells), "cells")
EOF
```
  Expected: `25 cells` (14 rewritten + 11 new).
  Then compile every code cell and look for undefined names across cells:
```bash
cd Synthetic_domain && cells_file="$(mktemp --suffix=.py)" && ../.venv/bin/python - "$cells_file" <<'EOF'
import sys

import nbformat

notebook = nbformat.read("ml_model.ipynb", as_version=4)
nbformat.validate(notebook)
sources = []
for index, cell in enumerate(notebook.cells):
    if cell.cell_type == "code":
        compile(cell.source, f"cell {index}", "exec")
        sources.append(cell.source)
with open(sys.argv[1], "w") as output:
    output.write("\n\n".join(sources))
print("valid; every code cell compiles")
EOF
../.venv/bin/flake8 --select=F821 --builtins=display "$cells_file" && echo "no undefined names"
```
  Expected: `valid; every code cell compiles`, then `no undefined names`.
  Check that no old global, no notebook copy of a `synthetic_model` name and no fit-folder read is left:
  ```bash
  cd Synthetic_domain && grep -n -e "ml_data.N_TARGETS" -e "ml_data.TARGET_NAMES" -e "CHOSEN" -e "PHOTOMETRIC_PREFIX\b" \
      -e "PHOTOMETRIC_PREFIXES = " -e "CROSS_FIT_ID_PREFIX" -e "import glob" -e "sm.load_fit" ml_model.ipynb
  ```
  Expected: no output.
  Check that every decode runs in parallel:
  ```bash
  cd Synthetic_domain && grep -c "decode_samples(" ml_model.ipynb && grep -c "decode_samples(.*n_jobs=N_DECODE_JOBS)" ml_model.ipynb
  ```
  Expected: `3`, then `3` (the chunk and truth decodes of `labelled_posteriors`, and the DR3 decode; each call is on one source line).

- [ ] **Step 4: Run the ML unit tests once more.**
  Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_ml_data.py tests/test_ml_networks.py tests/test_ml_evaluation.py tests/test_train_ml_model.py -q`
  Expected: all pass, with the skips of Tasks 7–8. The notebook itself runs only on the PC, after the DPL sets, the DPL MLP and the DPL photometric fits and cross-fits exist (controller tasks 12+).

- [ ] **Step 5: Commit.**
```bash
git add Synthetic_domain/ml_model.ipynb
git commit -m "feat: ml_model.ipynb evaluates both SFH models and compares them (spec Part 4)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 11: Remove the delayed-τ aliases; documentation

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (the lines marked `# delayed-tau alias, removed in Task 11`)
- Modify: any remaining user found by the grep below (modules, tests, notebooks)
- Modify: `README.md`, `Synthetic_domain/remote/README.md`
- Test: `Synthetic_domain/tests/test_sfh_models.py`

**Interfaces:**
- Consumes: everything from Tasks 1–10.
- Produces: `synthetic_model` without the delayed-only globals `DIMS`, `PHOTOMETRIC_DIMS`, `DRAWN_LABELS`, `DERIVED_LABELS`, `TRAINING_LABELS`, `CALIBRATION_LABELS`, `SFH_MODEL` and `components`. Code reads these from `MODELS[sfh_model]` instead.

- [ ] **Step 1: Write the failing test.** Append to `tests/test_sfh_models.py`:

```python
def test_synthetic_model_has_no_delayed_only_globals():
    """Task 11: every model-dependent name comes from sfh_models.MODELS, never from a delayed-tau alias."""
    import synthetic_model as sm
    removed = ["DIMS", "PHOTOMETRIC_DIMS", "DRAWN_LABELS", "DERIVED_LABELS", "TRAINING_LABELS",
               "CALIBRATION_LABELS", "SFH_MODEL", "components"]
    for name in removed:
        assert not hasattr(sm, name), name
```

- [ ] **Step 2: Run it to check it fails.**
  - Run: `cd Synthetic_domain && ../.venv/bin/python -m pytest tests/test_sfh_models.py::test_synthetic_model_has_no_delayed_only_globals -q`
  - Expected: FAIL, with the first alias name in the assertion message.

- [ ] **Step 3: Find the remaining users.**
  - Run (from `Synthetic_domain`):

    ```bash
    grep -n "delayed-tau alias, removed in Task 11" synthetic_model.py
    grep -rnE "sm\.(DIMS|PHOTOMETRIC_DIMS|DRAWN_LABELS|DERIVED_LABELS|TRAINING_LABELS|CALIBRATION_LABELS|SFH_MODEL|components)\b" \
        --include=*.py --include=*.ipynb . | grep -v "/pipes/"
    grep -rnE "^from synthetic_model import|synthetic_model\.(DIMS|TRAINING_LABELS|components)" --include=*.py . | grep -v "/pipes/"
    ```

  - Each hit in a module or test should read `sm.MODELS[sfh_model].<attribute>` with the model available at that point (`dims`, `photometric_dims`, `drawn_labels`, `derived_labels`, `training_labels`, `calibration_labels`, `name`, `components`).
  - Delayed-only code (the spectral fits) uses `sm.MODELS[sm.SPECTRAL_MODEL]`.
  - In notebooks, edit the cell source through a JSON script, as in Tasks 3, 6 and 10.

- [ ] **Step 4: Delete the alias lines** in `synthetic_model.py`, i.e. every line that the first grep of Step 3 printed.

- [ ] **Step 5: Run the whole suite and check the notebooks compile.**

```bash
cd Synthetic_domain && ../.venv/bin/python -m pytest -q
../.venv/bin/python - <<'EOF'
import json
for name in ["training_set.ipynb", "spectroscopic_reference.ipynb", "ml_model.ipynb"]:
    notebook = json.load(open(name))
    for index, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            lines = []
            for line in source.splitlines():
                lines.append("# " + line if line.lstrip().startswith(("%", "!")) else line)
            compile("\n".join(lines), f"{name}[{index}]", "exec")
    print(name, "compiles")
EOF
```

Expected:
- every test passes, with skips only for the tests that skip themselves today;
- the three notebooks compile;
- both greps of Step 3 print nothing.

- [ ] **Step 6: Documentation.**
  - **Root `README.md`:**
    - module map: add `sfh_models.py`, the shared prior plus one section per SFH model, and `MODELS`;
    - update the sentence calling `synthetic_model.py` the single home of the prior;
    - describe the two SFH models, with the DPL prior table from the spec;
    - ML targets: the 9 delayed-τ targets and the 10 DPL targets, from `MODELS[...].target_names`;
    - new outputs: `training_set_dblplaw.parquet`, `synthetic_validation_dblplaw.parquet`, `synthetic_test_dblplaw.parquet`, `photometric_fits_dblplaw/`, `photometric_cross_fits.csv`, `ml_runs/mlp_dblplaw/`, `ml_sfh_catalogue_dblplaw.parquet`, `ml_evaluation_summary_dblplaw.csv` and `ml_model_comparison.csv`;
    - the test count printed by Step 5;
    - one line in the design-choice table: "DPL added next to delayed-τ: compare SFH assumptions on the same colours".
  - **`Synthetic_domain/remote/README.md`:**
    - `FITS_SUBDIR=photometric_fits_dblplaw` for the DPL photometric fits;
    - `SFH_MODEL=dblplaw` for `submit_training.sh`;
    - the local supersede command covers the three fit folders, each with its model (Task 12, step 1).
  - Keep the README style: short sentences and tables; no results yet, which are filled in after the PC runs.

- [ ] **Step 7: Commit.**

```bash
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_sfh_models.py README.md Synthetic_domain/remote/README.md
git add -u Synthetic_domain   # the modules, tests and notebooks fixed in Step 3
git commit -m "refactor: drop the delayed-tau aliases; docs describe both SFH models

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

### Task 12: PC runs (controller)

Run each step only after the previous job has finished and been fetched; the PC runs one job at a time. Every command runs from `Synthetic_domain`.

- **Waiting:** wait with a background until-loop on `hljob ls`, never with `sleep` in the foreground.
- **Fetching:** fetch notebook jobs with `remote/fetch_notebook.sh <job-id>` and fit jobs with `remote/fetch_fits.sh <job-id>`.
- **Committing:** after each fetch, commit the tracked results with a `data:` message, as on 2026-10-27; fit folders are git-ignored.
- **If the PC is offline or slept,** resubmit the same command. The fit runner skips finished fits, and the training restarts from scratch.

- [ ] **Step 1: Supersede locally before every notebook-2 job** (file moves only, seconds):

```bash
cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -c "
import synthetic_model as sm
for folder, model in sm.FIT_FOLDER_MODELS.items():
    print(folder, len(sm.supersede_outdated_fits(folder, model)))"
```

Expected:
- 0 moves in `spectral_fits` and `photometric_fits`, since the delayed-τ fits are current;
- `photometric_fits_dblplaw` is created empty or holds only DPL fits;
- a repeat moves 0.

- [ ] **Step 2: DPL sets (CPU notebook job).**

```bash
CORES=16 MEMORY_MB=16000 remote/submit_notebook.sh training_set.ipynb dr3_corrected.csv splus_noise_model.csv \
    training_set.parquet synthetic_validation.parquet synthetic_test.parquet
```

- The delayed-τ sets travel only so that the notebook sees them and skips that model.
- Check that the notebook output says the delayed-τ sets were kept.
- Check that `training_set_dblplaw.parquet` (2**19 rows) and the two DPL synthetic sets (2**14 rows each) came back.
- Check that `git status` shows no change to the delayed-τ parquet files.
- If the upload passes about 1.5 GB, add `FLOAT32_FILES="training_set.parquet"`, so hljobd (3 GB limit) is not killed.

- [ ] **Step 3: DPL flow (GPU).**

```bash
SFH_MODEL=dblplaw MEMORY_MB=7100 remote/submit_training.sh
```

- When the job has finished, fetch it with `hljob fetch <job-id> -o /tmp/ml_job`.
- Then copy `/tmp/ml_job/out/ml_runs/mlp_dblplaw` to `ml_runs/mlp_dblplaw`.
- Check that `config.json` has `"sfh_model": "dblplaw"` and the 10 DPL `target_names`.

- [ ] **Step 4: Notebook 2 exports the DPL jobs (CPU notebook job).** Run Step 1 first. The bundle carries only the spectra of fitted galaxies, so it stays under hljobd's 3 GB limit:

```bash
fitted_spectra=()
for fit in spectral_fits/*.h5; do
    id="$(basename "$fit" .h5)"
    [[ -f "spectra/$id.fits" ]] && fitted_spectra+=("spectra/$id.fits")
done
CORES=16 MEMORY_MB=16000 TIME_LIMIT=08:00:00 FLOAT32_FILES=training_set.parquet \
    remote/submit_notebook.sh spectroscopic_reference.ipynb \
    dr3_corrected.csv ../Data/sdss_eboss_firefly-dr16_slim.fits ../dustMaps/csfd training_set.parquet \
    "${fitted_spectra[@]}" spectral_fits photometric_fits photometric_fits_dblplaw skipped_fit_jobs.txt \
    splus_aperture_photometry.csv splus_noise_model.csv spectroscopic_reference.csv literature
```

- Before submitting, check the bundle size with `BUNDLE_ONLY=/tmp/nb2_bundle` and the same command. It must stay below about 1.5 GB.
- After the fetch, check these exports:
  - `photometric_fits_dblplaw/jobs/` holds one job per validation galaxy, 33 `calib_*` DPL mocks and 33 `cross_calib_*` DPL fits of the delayed-τ mocks;
  - `photometric_fits/jobs/` gained 33 `cross_calib_*` delayed-τ fits of the DPL mocks.
- The notebook also recomputes the evidence errors of the old delayed-τ fits (`sm.fit_evidence`), which takes about 4 minutes.

- [ ] **Step 5: Photometric fits, DPL folder, then delayed-τ cross-fits (CPU).**

```bash
WORKERS=15 CORES=16 FITS_SUBDIR=photometric_fits_dblplaw remote/submit_fits.sh
# after it has been fetched (remote/fetch_fits.sh <job-id>):
WORKERS=15 CORES=16 FITS_SUBDIR=photometric_fits remote/submit_fits.sh
```

- The second job fits only the 33 new `cross_calib_*` jobs; the finished delayed-τ fits travel in the bundle and are skipped.
- Check with `run_spectral_fits.py --status` (printed in the job log) that every job has a converged fit, or record the ones that did not converge.

- [ ] **Step 6: Notebook 2 reads every fit (CPU notebook job).** Run Step 1, then the Step 4 command again.
- Check that `spectroscopic_reference.csv` gained these columns:
  - `phot_dblplaw_*`, `phot_dblplaw_fitter_calibrated`;
  - `phot_log_evidence(_error)` and `phot_dblplaw_log_evidence(_error)`.
- Check that `photometric_cross_fits.csv` exists with 66 rows.
- Note the DPL calibration result (95% coverage per label, which passes at ≥ 0.85) for the report.

- [ ] **Step 7: ML evaluation and comparison (GPU notebook job).**

```bash
GPU=1 CORES=16 MEMORY_MB=7100 TIME_LIMIT=12:00:00 FLOAT32_FILES="training_set.parquet training_set_dblplaw.parquet" \
    remote/submit_notebook.sh ml_model.ipynb ml_runs/mlp ml_runs/mlp_dblplaw \
    synthetic_validation.parquet synthetic_test.parquet training_set.parquet \
    synthetic_validation_dblplaw.parquet synthetic_test_dblplaw.parquet training_set_dblplaw.parquet \
    spectroscopic_reference.csv photometric_cross_fits.csv dr3_corrected.csv splus_noise_model.csv
```

- `RUN_LOCKED_TEST` stays `False`.
- Check the outputs: `ml_evaluation_summary.csv` and `ml_evaluation_summary_dblplaw.csv`, `ml_sfh_catalogue_dblplaw.parquet`, and `ml_model_comparison.csv`.
- Check that the delayed-τ results match the 2026-10-27 run (commit `9499b37`): same verdicts, and photometric `inside_95` within ±0.01.
- Commit and report to the user, covering:
  - the DPL calibration;
  - each model's works test against its own Bagpipes;
  - the model-against-model statistics;
  - all four estimators against the spectra, by galaxy type;
  - the evidence categories and the cross-fit recovery rate.

- [ ] **Step 8: Locked test (only when the user asks for it, once, for both models).**
  1. Set `EXPORT_LOCKED_TEST_JOBS = True` in notebook 2, then run Steps 4–6 for the test-half photometric jobs.
  2. Set `RUN_LOCKED_TEST = True` in `ml_model.ipynb`, then run Step 7.
  3. Commit and report.

The remaining delayed-τ spectral fits (604 galaxies, about 12 h) are independent of this plan. They run only when the user asks (`WORKERS=15 CORES=16 FITS_SUBDIR=spectral_fits remote/submit_fits.sh`). Afterwards, notebook 2 grows the validation half, and both models' photometric jobs grow with it (Steps 4–6).

---

## Part 2 of the plan: DPL spectral fits and the locked test (spec Amendment 2, 2026-10-29)

Tasks 1–12 are done. These tasks follow the contract in the ledger workspace (`contract2.md`). The delayed-τ spectral path must keep producing identical fits, because a PC job is fitting the remaining delayed-τ spectra from a bundle made at `9b4676e`.

**Controller rulings on the plan review:**
- **Spectral likelihood guard:** the try/except that returns LOG_ZERO where Bagpipes raises also covers delayed-τ. This is accepted: every value that computes is unchanged.
- **Stub check:** Task 15 Step 8 runs notebook cells on fake posteriors for about 30 s on the CPU. This is accepted as a quick check; it runs no fits and no flow.
- **Notebook-2 bundle:** when all spectra are fitted, the notebook-2 bundle sends only as many spectra as keep it under about 1 GB, and the notebook downloads the rest, as on 2026-10-27.


### Task 13: Spectral fits of either SFH model (`synthetic_model.py`), the DPL spectral jobs and mocks, and the remote scripts for `spectral_fits_dblplaw`

**Why:** Amendment 2 lifts decision 1. Every galaxy with a prepared spectrum is also fitted with the double power law, with the machinery of the delayed-τ spectral fits (calibration polynomial, likelihood, MultiNest settings). Today the spectral path is hard-wired to delayed-τ (`SPECTRAL_MODEL`, `SPECTRAL_DIMS`) and `fit_spectrum_job` refuses a DPL job. After this task every spectral function takes the model by name, as the photometric ones already do.

**The delayed-τ spectral path must not change.** A PC job is fitting the remaining delayed-τ spectral fits from a bundle made at commit 9b4676e. Its fits are read by notebook 2 together with the 788 already fetched. So for delayed-τ:
- the dims are the same 8 names in the same order (`spectral_dims("delayed")`);
- the components, drawn labels, posterior labels and mock truths are the same numbers. Two new regression tests relabel a stored fit and the stored mocks and compare them with the files to 1e-10 (they skip on a machine without `spectral_fits/`);
- the likelihood gives the same value at every point. The one change is that a point where Bagpipes raises now returns `LOG_ZERO`, as in the photometric fit, instead of killing the fit through PyMultiNest's `sys.exit`;
- the delayed-τ mocks are drawn from the same uniforms (`rng.random((n, 8))`). New mock jobs also carry `kind` and `sfh_model` keys; `export_job` wrote both before too, so the job files are the same.

**Existing tests (spec amendment 4).** They change only through the new `sfh_model` argument. Two tests stated the old decision and are replaced:
- `test_spectral_dims_stay_the_delayed_photometric_dims_plus_the_nuisances` asserted `sm.SPECTRAL_MODEL == "delayed"`. That constant is removed. The new `test_spectral_dims_are_each_model_s_photometric_dims_plus_the_nuisances` keeps its delayed-τ assertions (the same dims, now also as the literal list of 8 names) and its `len(... photometric_dims) == 7`;
- `test_fit_and_save_fails_a_spectral_job_of_the_dblplaw_model` asserted that a DPL spectral job fails. It becomes `test_fit_and_save_fits_each_spectral_job_with_its_own_model`.

`sm.SPECTRAL_DIMS` in a test becomes `N_DELAYED_SPECTRAL_DIMS` (`= len(sm.spectral_dims("delayed"))`) in `tests/test_synthetic_model.py` and a module constant `SPECTRAL_DIMS = sm.spectral_dims("delayed")` in `tests/test_spectral_validation.py`. Both have the old value, and the short names keep the lines under 110 characters.

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py`:
  - section 5: `SPECTRAL_MODEL` and `SPECTRAL_DIMS` (today lines 607–613) become `SPECTRAL_NUISANCE_DIMS` and `spectral_dims()`; `spectral_components` (636–656);
  - section 6: `fit_spectrum` (736–769), `posterior_labels` (823–845), `make_mock_spectrum` (1130–1161), the `save_fit` docstring (1172), `fit_spectrum_job` (1245–1258);
  - section 7: `SPECTRAL_PREFIXES` after `PHOTOMETRIC_PREFIXES` (1377), `draw_calibration_cubes` and `draw_calibration_mocks` (1439–1468), new `dblplaw_spectral_jobs` after them;
  - section 8: `FIT_FOLDER_MODELS` (1681–1682) and the `supersede_outdated_fits` docstring (1754–1755).
- Modify: `Synthetic_domain/spectral_validation.py`: new constant `RULER_SFH_MODEL = "delayed"` and the three `sm.spectral_components` calls (lines 193, 211, 343) with their docstrings (204–205, 333). The rulers stay on the delayed-τ fits (amendment 2, item 3).
- Modify: `Synthetic_domain/run_spectral_fits.py` (module docstring, lines 7–10).
- Modify: `Synthetic_domain/spectroscopic_reference.ipynb`, names only: in cells 11, 14, 16 and 18 the changed functions get `"delayed"`, and `sm.SPECTRAL_MODEL` / `sm.SPECTRAL_DIMS` give way to their values. Without this step `tests/test_notebooks.py` fails, and cell 11 would call `draw_calibration_mocks` without its model. Task 14 then builds on this version.
- Modify: `Synthetic_domain/remote/submit_fits.sh` (header lines 8–12, `SKIP_LIST` help lines 24–25, the `SKIP_LIST` branch lines 44–48), `Synthetic_domain/remote/fit_job.sh` (lines 13–14), `Synthetic_domain/remote/fetch_fits.sh` (lines 6–8), `Synthetic_domain/remote/README.md` (after line 40, and lines 45–46), `.gitignore` (repo root, after line 21).
- Create: `Synthetic_domain/skipped_fit_jobs_dblplaw.txt` (comments only).
- Test: `Synthetic_domain/tests/test_synthetic_model.py`, `Synthetic_domain/tests/test_spectral_validation.py`; `tests/test_notebooks.py` must stay green.

**Interfaces:**
- **Consumes:** `MODELS[...]` (`photometric_dims`, `components`, `training_labels`, `drawn_labels`, `calibration_labels`), `training_unit_vector`, `summary_labels` through `sample_label_row`, `run_multinest`, `job_kind`, `job_model`, `SPECTRUM_JOB`, `CALIBRATION_ID_PREFIX`, `calibration_mocks_are_current`, `supersede_outdated_fits`, `model_file`.
- **Produces (`synthetic_model`):**
  - `SPECTRAL_NUISANCE_DIMS = ["veldisp", "noise_scale"]`;
  - `spectral_dims(sfh_model) -> list`: `MODELS[sfh_model].photometric_dims + SPECTRAL_NUISANCE_DIMS` (8 names for delayed-τ, 9 for DPL);
  - `spectral_components(cube, redshift, sfh_model) -> (bagpipes_components, drawn, noise_scale)`;
  - `fit_spectrum(spectrum, redshift, sfh_model, output_basename, max_iterations, n_live_points=N_LIVE_POINTS) -> dict` (the model comes before `output_basename`, as in `fit_photometry`);
  - `posterior_labels(unit_cube_samples, redshift, sfh_model) -> DataFrame`;
  - `make_mock_spectrum(cube, redshift, template_spectrum, rng, sfh_model) -> (spectrum, truth)`;
  - `draw_calibration_cubes(n_mocks, rng, sfh_model) -> ndarray (n_mocks, len(spectral_dims(sfh_model)))`;
  - `draw_calibration_mocks(n_mocks, templates, rng, sfh_model) -> (jobs, truths)`; the jobs now carry `kind = SPECTRUM_JOB` and `sfh_model`;
  - `dblplaw_spectral_jobs(delayed_jobs) -> list`: one new job per galaxy, with the same `job_id`, `spectrum` and `redshift`, `kind = SPECTRUM_JOB` and `sfh_model = "dblplaw"`. Calibration mocks (`calib_*`) are left out, and anything but a delayed-τ spectrum job raises `ValueError`;
  - `fit_spectrum_job(job, ...)` fits with `job_model(job)`; a DPL spectrum job is no longer refused;
  - `SPECTRAL_PREFIXES = {"delayed": "", "dblplaw": "spec_dblplaw_"}`: the reference-CSV prefix of each model's spectral columns. Task 14 uses it, and Task 15 must use it too, not define it again;
  - `FIT_FOLDER_MODELS["spectral_fits_dblplaw"] = "dblplaw"`; `fit_is_current` and `supersede_outdated_fits` work for that folder unchanged;
  - removed: `SPECTRAL_MODEL` and `SPECTRAL_DIMS` (the notebooks, modules and tests no longer use them).
- **Produces (`spectral_validation`):** `RULER_SFH_MODEL = "delayed"`; every ruler reads the delayed-τ spectral fits.
- **Produces (remote):** `FITS_SUBDIR=spectral_fits_dblplaw remote/submit_fits.sh` bundles `spectral_fits_dblplaw/jobs/` and the finished `spectral_fits_dblplaw/*.h5`. Its `SKIP_LIST` default is the new `skipped_fit_jobs_dblplaw.txt` (comments only, so no job is left out). Each model gets its own list because the DPL mocks reuse the `calib_<k>` ids of the delayed-τ mocks: the delayed-τ list would drop DPL mocks, and the reverse. The photometric folders keep having no list. `fit_job.sh` records the folder in `out/fits_subdir.txt`, and `fetch_fits.sh` copies into the recorded folder (and stops if a `FITS_SUBDIR` you set differs from it). Neither needs a logic change, because neither has a branch on the folder name; only their comments name the new folder.

**Running the scripts.** Steps 1, 3, 4, 5 and 7 are Python scripts. Save each one in your scratchpad, not in the repository, under the name the step gives. Each replacement asserts that its old text occurs exactly as often as stated, so a script stops without writing the file in hand when that file differs from what this plan expects. Files it patched before that one stay patched: restore them with `git checkout -- <file>` before running it again. Set `SCRATCH` to your scratchpad directory first. Every command runs from `Synthetic_domain/` unless the step says otherwise.

- [ ] **Step 1: Write the failing tests**

Save as `$SCRATCH/task13_tests.py` and run `../.venv/bin/python "$SCRATCH/task13_tests.py" .`:

```python
"""Task 13, step 1: the tests of the spectral fits of either SFH model. Argument: the Synthetic_domain folder."""
import sys

root = sys.argv[1]


def patch(path, pairs):
    text = open(path).read()
    for old, new, count in pairs:
        assert text.count(old) == count, (path, text.count(old), old[:80])
        text = text.replace(old, new)
    open(path, "w").write(text)
    print(path, "patched")


# =============================================================================================================
# tests/test_spectral_validation.py: only the new sfh_model argument (the rulers read delayed-tau fits)
# =============================================================================================================
patch(root + "/tests/test_spectral_validation.py", [
    ("sm.spectral_components(cube, redshift)", "sm.spectral_components(cube, redshift, \"delayed\")", 2),
    ("sm.spectral_components(cube, 0.1)", "sm.spectral_components(cube, 0.1, \"delayed\")", 1),
    ("sm.SPECTRAL_DIMS", "SPECTRAL_DIMS", 6),
    ("""GRID = np.arange(3700.0, 5600.0, 1.0)      # rest-frame air Angstrom, 1 A pixels
""", """GRID = np.arange(3700.0, 5600.0, 1.0)      # rest-frame air Angstrom, 1 A pixels
SPECTRAL_DIMS = sm.spectral_dims("delayed")  # unit-cube dims of the delayed-tau fits the rulers read
""", 1),
])

# =============================================================================================================
# tests/test_synthetic_model.py
# =============================================================================================================
NEW_SPECTRAL_COMPONENT_TESTS = '''

def test_dblplaw_spectral_components_are_its_photometric_components_plus_the_nuisances():
    """A DPL spectral point is the DPL photometric point at fixed redshift, plus veldisp and noise_scale."""
    dims = sm.spectral_dims("dblplaw")
    cube = np.random.default_rng(50).random(len(dims))
    photometric_part = cube[:len(sm.MODELS["dblplaw"].photometric_dims)]
    expected_components, expected_drawn = sm.photometric_components(photometric_part, 0.08, "dblplaw")

    bagpipes_components, drawn, noise_scale = sm.spectral_components(cube, 0.08, "dblplaw")

    veldisp = sm.uniform(cube[dims.index("veldisp")], *sm.VELDISP_RANGE)
    assert bagpipes_components["veldisp"] == veldisp
    assert drawn["veldisp"] == veldisp
    assert noise_scale == 10 ** sm.uniform(cube[dims.index("noise_scale")], *sm.LOG_NOISE_SCALE_RANGE)
    del bagpipes_components["veldisp"]
    del drawn["veldisp"]
    assert bagpipes_components == expected_components
    assert drawn == expected_drawn
    assert "dblplaw" in bagpipes_components and "delayed" not in bagpipes_components
'''

NEW_FIT_FILE_TESTS = '''

def test_spectral_dims_are_each_model_s_photometric_dims_plus_the_nuisances():
    """Each SFH model has its own spectral fits (DPL spec, amendment 2); delayed-tau keeps its 8 dims."""
    assert sm.spectral_dims("delayed") == [*sm.MODELS["delayed"].photometric_dims, "veldisp", "noise_scale"]
    assert sm.spectral_dims("delayed") == ["log_age_fraction", "log_tau", "log_metallicity", "av_ism",
                                           "log_tau_bc", "logU", "veldisp", "noise_scale"]
    assert sm.spectral_dims("dblplaw") == [*sm.MODELS["dblplaw"].photometric_dims, "veldisp", "noise_scale"]
    assert len(sm.MODELS["dblplaw"].photometric_dims) == 7
    assert not hasattr(sm, "SPECTRAL_DIMS") and not hasattr(sm, "SPECTRAL_MODEL")


STORED_SPECTRAL_FITS = "spectral_fits"          # the delayed-tau spectral fits (git-ignored; on this machine)
N_REGRESSION_SAMPLES = 3                        # samples relabelled from the stored fit (~0.1 s each)
STORED_LABEL_RTOL = 1e-10                       # "exactly", up to floating-point rounding


def stored_spectral_fit_files():
    """Fit files of galaxies (not mocks) in STORED_SPECTRAL_FITS, sorted; empty when the folder is missing."""
    if not os.path.isdir(STORED_SPECTRAL_FITS):
        return []
    names = sorted(name for name in os.listdir(STORED_SPECTRAL_FITS)
                   if name.endswith(".h5") and not name.startswith(sm.CALIBRATION_ID_PREFIX))
    return [os.path.join(STORED_SPECTRAL_FITS, name) for name in names]


@pytest.mark.skipif(len(stored_spectral_fit_files()) == 0, reason="no stored delayed-tau spectral fit here")
def test_delayed_spectral_labels_reproduce_a_stored_spectral_fit():
    """The refactored delayed-tau spectral path labels a stored fit's samples exactly as its fit did."""
    fit_file = stored_spectral_fit_files()[0]
    stored_posterior, metadata = sm.load_fit(fit_file)
    assert metadata["sfh_model"] == "delayed"
    # The stored rows are the labels of the thinned samples, in order
    kept_samples = sm.thinned_samples(sm.load_unit_cube_samples(fit_file))
    assert kept_samples.shape[1] == len(sm.spectral_dims("delayed"))

    relabelled = sm.posterior_labels(kept_samples[:N_REGRESSION_SAMPLES], metadata["redshift"], "delayed")

    assert set(relabelled.columns) == set(stored_posterior.columns)
    for column in stored_posterior.columns:
        assert np.allclose(relabelled[column], stored_posterior[column].iloc[:N_REGRESSION_SAMPLES],
                           rtol=STORED_LABEL_RTOL, atol=0, equal_nan=True), column


STORED_MOCK_FILE = os.path.join(STORED_SPECTRAL_FITS, sm.CALIBRATION_MOCK_FILE)
NOTEBOOK_SEED = 0                               # SEED of spectroscopic_reference.ipynb, which drew the mocks


@pytest.mark.skipif(not os.path.exists(STORED_MOCK_FILE), reason="no stored delayed-tau spectral mocks here")
def test_delayed_calibration_mock_truths_reproduce_the_stored_mocks():
    """The delayed-tau mock points, redrawn from the notebook's seed, give the stored truths exactly."""
    mock_jobs, truths = joblib.load(STORED_MOCK_FILE)
    cubes = sm.draw_calibration_cubes(len(mock_jobs), np.random.default_rng(NOTEBOOK_SEED), "delayed")
    noise_index = sm.spectral_dims("delayed").index("noise_scale")
    for mock_index in (0, len(mock_jobs) - 1):
        cube = cubes[mock_index].copy()
        cube[noise_index] = sm.UNIT_NOISE_CUBE_VALUE       # as make_mock_spectrum sets it
        mock_job = mock_jobs[mock_index]
        truth = sm.posterior_labels(cube[np.newaxis, :], mock_job["redshift"], "delayed").iloc[0]
        stored_truth = truths[mock_job["job_id"]]
        assert list(truth.index) == list(stored_truth.index)
        assert np.allclose(truth.to_numpy(dtype=float), stored_truth.to_numpy(dtype=float),
                           rtol=STORED_LABEL_RTOL, atol=0)


def test_dblplaw_posterior_labels_give_every_dblplaw_training_label():
    """Three DPL spectral samples: one row each, every DPL training label finite and as the prior maps it."""
    cubes = np.random.default_rng(51).random((3, len(sm.spectral_dims("dblplaw"))))
    posterior = sm.posterior_labels(cubes, 0.1, "dblplaw")

    training_labels = sm.MODELS["dblplaw"].training_labels
    assert len(posterior) == 3
    assert list(posterior.columns[:len(training_labels)]) == training_labels
    for column in (*training_labels, "veldisp", "noise_scale", "model_r"):
        assert np.isfinite(posterior[column]).all(), column
    for row_index, cube in enumerate(cubes):
        _, drawn, noise_scale = sm.spectral_components(cube, 0.1, "dblplaw")
        for label in sm.MODELS["dblplaw"].drawn_labels:
            assert posterior[label].iloc[row_index] == drawn[label], label
        assert posterior["noise_scale"].iloc[row_index] == noise_scale
    fractions = posterior[["f_young", "f_mid", "f_old"]].sum(axis=1)
    assert np.allclose(fractions, 1.0)


def test_make_mock_spectrum_of_the_dblplaw_model_returns_its_truth_labels():
    """A DPL mock spectrum is finite and its truth carries every DPL training and calibration label."""
    job, truth = mock_job("mock_dblplaw", np.random.default_rng(52), "dblplaw")
    assert job["sfh_model"] == "dblplaw"
    assert np.isfinite(job["spectrum"]["flux"]).all()
    for column in (*sm.MODELS["dblplaw"].training_labels, "veldisp", "noise_scale", "model_r"):
        assert np.isfinite(truth[column]), column
    assert set(sm.MODELS["dblplaw"].calibration_labels) <= set(truth.index)
    assert truth["noise_scale"] == 1.0    # the mock noise is the template's own
'''

NEW_JOB_TESTS = '''

def test_calibration_cubes_of_each_model_fill_its_spectral_cube():
    """Delayed-tau cubes are the same uniforms as before; DPL cubes have the 9 DPL spectral dimensions."""
    delayed_cubes = sm.draw_calibration_cubes(33, np.random.default_rng(0), "delayed")
    assert np.array_equal(delayed_cubes, np.random.default_rng(0).random((33, 8)))
    dblplaw_cubes = sm.draw_calibration_cubes(33, np.random.default_rng(0), "dblplaw")
    assert dblplaw_cubes.shape == (33, len(sm.spectral_dims("dblplaw"))) == (33, 9)
    assert np.all((dblplaw_cubes >= 0) & (dblplaw_cubes < 1))


def test_draw_dblplaw_calibration_mocks_are_dblplaw_spectrum_jobs_with_current_truths(tmp_path):
    """DPL mocks are DPL spectrum jobs, and a mock file of them is current in a DPL folder only."""
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    template = sm.prepare_spectrum(wavelength_vacuum, np.ones_like(wavelength_vacuum),
                                   np.full_like(wavelength_vacuum, 100.0), ebv=0.0)

    jobs, truths = sm.draw_calibration_mocks(2, [(template, 0.1)], np.random.default_rng(53), "dblplaw")

    assert [job["job_id"] for job in jobs] == ["calib_000", "calib_001"]
    for job in jobs:
        assert job["kind"] == sm.SPECTRUM_JOB
        assert sm.job_model(job) == "dblplaw"
        assert set(sm.MODELS["dblplaw"].calibration_labels) <= set(truths[job["job_id"]].index)
    mock_file = str(tmp_path / sm.CALIBRATION_MOCK_FILE)
    joblib.dump((jobs, truths), mock_file)
    assert sm.calibration_mocks_are_current(mock_file, "dblplaw")
    assert not sm.calibration_mocks_are_current(mock_file, "delayed")


def test_dblplaw_spectral_jobs_copy_the_galaxies_and_leave_out_the_mocks():
    """The galaxies' delayed-tau spectral jobs become DPL jobs with the same id, spectrum and redshift."""
    spectrum = {"wavelength": np.arange(4000.0, 4010.0), "flux": np.ones(10), "ivar": np.ones(10)}
    delayed_jobs = [{"job_id": "calib_000", "spectrum": spectrum, "redshift": 0.2},
                    {"job_id": "1105771463723476992", "spectrum": spectrum, "redshift": 0.1},
                    {"job_id": "10311134961714089984", "kind": "spectrum", "sfh_model": "delayed",
                     "spectrum": spectrum, "redshift": 0.3}]

    jobs = sm.dblplaw_spectral_jobs(delayed_jobs)

    assert [job["job_id"] for job in jobs] == ["1105771463723476992", "10311134961714089984"]
    for job, delayed_job in zip(jobs, delayed_jobs[1:]):
        assert job["kind"] == sm.SPECTRUM_JOB
        assert job["sfh_model"] == "dblplaw"
        assert job["redshift"] == delayed_job["redshift"]
        assert job["spectrum"] is delayed_job["spectrum"]
    assert "sfh_model" not in delayed_jobs[1]          # the delayed-tau jobs are left alone
    assert delayed_jobs[2]["sfh_model"] == "delayed"


def test_dblplaw_spectral_jobs_refuse_a_job_that_is_not_a_delayed_spectrum_job():
    spectrum = {"wavelength": np.arange(4000.0, 4010.0), "flux": np.ones(10), "ivar": np.ones(10)}
    photometry = {"job_id": "123", "kind": "photometry", "magnitudes": np.full(12, 20.0),
                  "errors": np.full(12, 0.05), "redshift": 0.1}
    dblplaw_spectrum = {"job_id": "456", "sfh_model": "dblplaw", "spectrum": spectrum, "redshift": 0.1}

    for job in (photometry, dblplaw_spectrum):
        with pytest.raises(ValueError, match="not a delayed-tau spectrum job"):
            sm.dblplaw_spectral_jobs([job])
'''

NEW_FIT_ROUTING_TESTS = '''def test_fit_and_save_fits_each_spectral_job_with_its_own_model(tmp_path, monkeypatch):
    """A DPL spectral job is fitted and labelled as DPL (DPL spec, amendment 2), a job without sfh_model as
    delayed-tau; both store their model."""
    fitted_models = []
    labelled_models = []

    def capturing_fit_spectrum(spectrum, redshift, sfh_model, output_basename, max_iterations, n_live_points):
        fitted_models.append(sfh_model)
        return fake_fit(len(sm.spectral_dims(sfh_model)))

    def capturing_labels(unit_cube_samples, redshift, sfh_model):
        labelled_models.append(sfh_model)
        return pd.DataFrame({"mw_age": [1.0, 2.0]})

    monkeypatch.setattr(sm, "fit_spectrum", capturing_fit_spectrum)
    monkeypatch.setattr(sm, "posterior_labels", capturing_labels)
    spectrum = {"wavelength": np.arange(4000.0, 4010.0), "flux": np.ones(10), "ivar": np.ones(10)}
    dblplaw_job = {"job_id": "spectral_dpl", "kind": "spectrum", "sfh_model": "dblplaw", "redshift": 0.1,
                   "spectrum": spectrum}
    old_job = {"job_id": "spectral_old", "redshift": 0.1, "spectrum": spectrum}   # as every delayed-tau job

    for job in (dblplaw_job, old_job):
        summary = sm.fit_and_save(job, str(tmp_path), max_iterations=10, n_live_points=25)
        assert summary["status"] == "fitted", summary["error"]

    assert fitted_models == ["dblplaw", "delayed"]
    assert labelled_models == ["dblplaw", "delayed"]
    _, dblplaw_metadata = sm.load_fit(str(tmp_path / "spectral_dpl.h5"))
    assert dblplaw_metadata["sfh_model"] == "dblplaw"
    assert dblplaw_metadata["kind"] == "spectrum"
    assert dblplaw_metadata["log_evidence_error"] == 0.1
    assert sm.load_unit_cube_samples(str(tmp_path / "spectral_dpl.h5")).shape[1] == 9
    _, old_metadata = sm.load_fit(str(tmp_path / "spectral_old.h5"))
    assert old_metadata["sfh_model"] == "delayed"


class ZeroNoise:
    """Stands in for the numpy Generator of make_mock_spectrum: no noise, so the mock is the model itself."""

    def normal(self, mean, sigma, size):
        """Zeros of the requested length (dimensionless), whatever mean and sigma are."""
        return np.zeros(size)


@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_fit_and_save_dblplaw_spectral_smoke_run_on_a_noiseless_mock(tmp_path):
    """A short MultiNest run of a DPL spectrum job (noiseless mock) is saved with DPL labels and model."""
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    template = sm.prepare_spectrum(wavelength_vacuum, np.ones_like(wavelength_vacuum),
                                   np.full_like(wavelength_vacuum, 100.0), ebv=0.0)
    cube = np.random.default_rng(54).random(len(sm.spectral_dims("dblplaw")))
    spectrum, _ = sm.make_mock_spectrum(cube, 0.1, template, ZeroNoise(), "dblplaw")
    job = {"job_id": "dpl_spectral_smoke", "kind": "spectrum", "sfh_model": "dblplaw", "spectrum": spectrum,
           "redshift": 0.1}

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=100, n_live_points=40)

    assert summary["status"] == "fitted", summary["error"]
    posterior, metadata = sm.load_fit(str(tmp_path / "dpl_spectral_smoke.h5"))
    assert metadata["sfh_model"] == "dblplaw"
    assert metadata["kind"] == "spectrum"
    assert "log_evidence_error" in metadata
    assert sm.load_unit_cube_samples(str(tmp_path / "dpl_spectral_smoke.h5")).shape[1] == 9
    for column in (*sm.MODELS["dblplaw"].training_labels, "veldisp", "noise_scale", "model_r"):
        assert column in posterior.columns, column
        assert np.isfinite(posterior[column]).all(), column
'''

NEW_SUPERSEDE_TESTS = '''

def test_supersede_keeps_dblplaw_spectral_fits_only_in_their_folder(tmp_path):
    """spectral_fits_dblplaw is a DPL folder: a delayed-tau fit leaves it; spectral_fits does the reverse."""
    dblplaw_folder = tmp_path / "spectral_fits_dblplaw"
    delayed_folder = tmp_path / "spectral_fits"
    for folder in (dblplaw_folder, delayed_folder):
        folder.mkdir()
        write_fit(folder / "dpl.h5", sfh_model="dblplaw")
        write_fit(folder / "delayed.h5", sfh_model="delayed")

    assert sm.supersede_outdated_fits(str(dblplaw_folder), "dblplaw") == ["delayed.h5"]
    assert sm.supersede_outdated_fits(str(delayed_folder), "delayed") == ["dpl.h5"]
    assert (dblplaw_folder / "dpl.h5").exists() and (delayed_folder / "delayed.h5").exists()
    with pytest.raises(ValueError, match="dblplaw"):
        sm.supersede_outdated_fits(str(dblplaw_folder), "delayed")
'''

patch(root + "/tests/test_synthetic_model.py", [
    # --- Spectra: the new sfh_model argument, and the DPL spectral components
    ('''def test_spectral_components_fix_redshift_and_add_veldisp():
    cube = np.full(len(sm.SPECTRAL_DIMS), 0.5)
    bagpipes_components, drawn, noise_scale = sm.spectral_components(cube, redshift=0.08)

    assert "z" not in sm.SPECTRAL_DIMS
''', '''def test_spectral_components_fix_redshift_and_add_veldisp():
    cube = np.full(len(sm.spectral_dims("delayed")), 0.5)
    bagpipes_components, drawn, noise_scale = sm.spectral_components(cube, redshift=0.08, sfh_model="delayed")

    assert "z" not in sm.spectral_dims("delayed")
''', 1),
    ('''    assert bagpipes_components["delayed"]["age"] == drawn["age"] <= sm.maximum_age(0.08)
''', '''    assert bagpipes_components["delayed"]["age"] == drawn["age"] <= sm.maximum_age(0.08)
''' + NEW_SPECTRAL_COMPONENT_TESTS, 1),
    # --- mock_job of either model
    ('''def mock_job(job_id, rng):
    """A small mock spectrum job at z = 0.1.

    Inputs: job_id: identifier (str); rng: numpy Generator.
    Returns: (job dict with job_id, spectrum and redshift, truth labels as a Series).
    """''', '''def mock_job(job_id, rng, sfh_model="delayed"):
    """A small mock spectrum job of sfh_model at z = 0.1.

    Inputs: job_id: identifier (str); rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: (job dict with job_id, spectrum and redshift, and sfh_model for a DPL job (a delayed-tau job has
             no sfh_model, as every job exported before the DPL spectral fits), truth labels as a Series).
    """''', 1),
    ('''    cube = rng.random(len(sm.SPECTRAL_DIMS))
    spectrum, truth = sm.make_mock_spectrum(cube, 0.1, template, rng)
    return {"job_id": job_id, "spectrum": spectrum, "redshift": 0.1}, truth''',
     '''    cube = rng.random(len(sm.spectral_dims(sfh_model)))
    spectrum, truth = sm.make_mock_spectrum(cube, 0.1, template, rng, sfh_model)
    job = {"job_id": job_id, "spectrum": spectrum, "redshift": 0.1}
    if sfh_model != "delayed":
        job["sfh_model"] = sfh_model
    return job, truth''', 1),
    ("sm.MODELS[sm.SPECTRAL_MODEL]", 'sm.MODELS["delayed"]', 5),
    ("len(sm.SPECTRAL_DIMS)", "N_DELAYED_SPECTRAL_DIMS", 13),
    ('''LOW_REDSHIFT = 0.002             # Bagpipes needs z > 0 to build a star formation history
''', '''LOW_REDSHIFT = 0.002             # Bagpipes needs z > 0 to build a star formation history
N_DELAYED_SPECTRAL_DIMS = len(sm.spectral_dims("delayed"))   # dims of a delayed-tau spectral fit (8)
''', 1),
    # --- the fit folders
    ('''    assert sm.FIT_FOLDER_MODELS == {"spectral_fits": "delayed", "photometric_fits": "delayed",
                                    "photometric_fits_dblplaw": "dblplaw"}''',
     '''    assert sm.FIT_FOLDER_MODELS == {"spectral_fits": "delayed", "photometric_fits": "delayed",
                                    "photometric_fits_dblplaw": "dblplaw", "spectral_fits_dblplaw": "dblplaw"}''',
     1),
    # --- spectral dims, stored-fit regression and DPL labels and mocks
    ('''def test_spectral_dims_stay_the_delayed_photometric_dims_plus_the_nuisances():
    """Spectral fits are delayed-tau only (DPL spec, decision 1); the DPL photometric fit has 7 dimensions."""
    assert sm.SPECTRAL_MODEL == "delayed"
    assert sm.SPECTRAL_DIMS == [*sm.MODELS["delayed"].photometric_dims, "veldisp", "noise_scale"]
    assert len(sm.MODELS["dblplaw"].photometric_dims) == 7
''', NEW_FIT_FILE_TESTS.lstrip("\n"), 1),
    ('''    posterior = sm.posterior_labels(cubes, 0.1)

    training_labels = sm.MODELS["delayed"].training_labels''',
     '''    posterior = sm.posterior_labels(cubes, 0.1, "delayed")

    training_labels = sm.MODELS["delayed"].training_labels''', 1),
    ('''    _, drawn, _ = sm.spectral_components(cubes[0], 0.1)''',
     '''    _, drawn, _ = sm.spectral_components(cubes[0], 0.1, "delayed")''', 1),
    # --- the fake fitters take the model
    ('''    def fake_fit_spectrum(spectrum, redshift, output_basename, max_iterations, n_live_points):
        called.append("spectrum")
        return fake_fit(N_DELAYED_SPECTRAL_DIMS)''',
     '''    def fake_fit_spectrum(spectrum, redshift, sfh_model, output_basename, max_iterations, n_live_points):
        called.append("spectrum")
        return fake_fit(len(sm.spectral_dims(sfh_model)))''', 1),
    ('''    def fake_spectral_labels(unit_cube_samples, redshift):''',
     '''    def fake_spectral_labels(unit_cube_samples, redshift, sfh_model):''', 1),
    # --- a DPL spectral job is now fitted (amendment 2 lifts decision 1)
    ('''def test_fit_and_save_fails_a_spectral_job_of_the_dblplaw_model(tmp_path, fake_fitters):
    """Spectral fits are delayed-tau only (DPL spec, decision 1)."""
    job = {"job_id": "spectral_dpl", "kind": "spectrum", "sfh_model": "dblplaw", "redshift": 0.1,
           "spectrum": {"wavelength": np.arange(4000.0, 4010.0), "flux": np.ones(10), "ivar": np.ones(10)}}

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=10, n_live_points=25)

    assert summary["status"] == "failed"
    assert "delayed only" in summary["error"]
    assert fake_fitters == []
''', NEW_FIT_ROUTING_TESTS, 1),
    # --- calibration mocks of each model, and the DPL spectral jobs
    ('''    first = sm.draw_calibration_cubes(60, np.random.default_rng(0))
    second = sm.draw_calibration_cubes(60, np.random.default_rng(0))''',
     '''    first = sm.draw_calibration_cubes(60, np.random.default_rng(0), "delayed")
    second = sm.draw_calibration_cubes(60, np.random.default_rng(0), "delayed")''', 1),
    ('''    posterior = sm.posterior_labels(cubes, 0.1)
    truth = posterior.iloc[1]''', '''    posterior = sm.posterior_labels(cubes, 0.1, "delayed")
    truth = posterior.iloc[1]''', 1),
    ('''    jobs, truths = sm.draw_calibration_mocks(3, templates, np.random.default_rng(4))
    repeat_jobs, _ = sm.draw_calibration_mocks(3, templates, np.random.default_rng(4))''',
     '''    jobs, truths = sm.draw_calibration_mocks(3, templates, np.random.default_rng(4), "delayed")
    repeat_jobs, _ = sm.draw_calibration_mocks(3, templates, np.random.default_rng(4), "delayed")''', 1),
    ('''        assert np.array_equal(job["spectrum"]["flux"], repeat["spectrum"]["flux"])
        assert job["redshift"] == repeat["redshift"]
''', '''        assert np.array_equal(job["spectrum"]["flux"], repeat["spectrum"]["flux"])
        assert job["redshift"] == repeat["redshift"]
''' + NEW_JOB_TESTS, 1),
    # --- the spectral columns of the reference CSV
    ('''    assert sm.PHOTOMETRIC_PREFIXES == {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}
''', '''    assert sm.PHOTOMETRIC_PREFIXES == {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}


def test_spectral_prefixes_keep_the_delayed_names_plain():
    """The delayed-tau spectral columns keep the names of the first fits; the DPL ones are spec_dblplaw_*."""
    assert sm.SPECTRAL_PREFIXES == {"delayed": "", "dblplaw": "spec_dblplaw_"}
    assert set(sm.SPECTRAL_PREFIXES) == set(sm.MODELS)
''', 1),
    # --- superseding: the DPL spectral folder
    ('''    with pytest.raises(ValueError, match="dblplaw"):
        sm.supersede_outdated_fits(str(folder), "delayed")
    assert (folder / "dpl.h5").exists()
''', '''    with pytest.raises(ValueError, match="dblplaw"):
        sm.supersede_outdated_fits(str(folder), "delayed")
    assert (folder / "dpl.h5").exists()
''' + NEW_SUPERSEDE_TESTS, 1),
])
```

Expected output: `./tests/test_spectral_validation.py patched` and `./tests/test_synthetic_model.py patched`.

- [ ] **Step 2: Run the tests and watch them fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py tests/test_spectral_validation.py -q`

Expected: two collection errors, `AttributeError: module 'synthetic_model' has no attribute 'spectral_dims'` (both files call it at import).

- [ ] **Step 3: `synthetic_model.py`: spectral dims, components, fit, labels and mocks per SFH model; the DPL spectral jobs; the new folder**

Save as `$SCRATCH/task13_synthetic_model.py` and run `../.venv/bin/python "$SCRATCH/task13_synthetic_model.py" synthetic_model.py`:

```python
"""Task 13, step 3: synthetic_model.py fits, labels and mocks spectra of either SFH model. Argument: the file."""
import sys

PATH = sys.argv[1]
text = open(PATH).read()


def replace(old, new):
    global text
    assert text.count(old) == 1, (text.count(old), old[:80])
    text = text.replace(old, new)


# 1. spectral_dims instead of SPECTRAL_MODEL / SPECTRAL_DIMS
replace('''# Spectral fits are delayed-tau only (DPL spec, decision 1): the DPL model is fitted to the photometry alone
SPECTRAL_MODEL = "delayed"

# Redshift is fixed to the spectroscopic value, so it leaves the list; veldisp and noise_scale join it
SPECTRAL_DIMS = list(MODELS[SPECTRAL_MODEL].photometric_dims)
SPECTRAL_DIMS.append("veldisp")
SPECTRAL_DIMS.append("noise_scale")
''', '''# The spectral fit's own nuisance parameters, after the SFH model's dims: velocity dispersion, noise scale
SPECTRAL_NUISANCE_DIMS = ["veldisp", "noise_scale"]


def spectral_dims(sfh_model):
    """Unit-cube dimensions of a spectral fit of one SFH model, in order.

    Why: the redshift is fixed to the spectroscopic value, so it leaves the model's dims, and the spectral
    nuisance parameters join them. Each SFH model has its own spectral fits (DPL spec, amendment 2).
    Inputs: sfh_model: "delayed" or "dblplaw".
    Returns: new list of names: MODELS[sfh_model].photometric_dims, then SPECTRAL_NUISANCE_DIMS (8 names for
             delayed-tau, 9 for the double power law).
    """
    return [*MODELS[sfh_model].photometric_dims, *SPECTRAL_NUISANCE_DIMS]
''')

replace('''def spectral_components(cube, redshift):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit unit-cube point.

    The SFH, metallicity, dust and nebular part is the delayed-tau training prior (MODELS[SPECTRAL_MODEL]) at
    the fixed spectroscopic redshift, which also sets the oldest allowed age; veldisp and noise_scale are the
    spectral nuisance parameters.
    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as SPECTRAL_DIMS;
            redshift: spectroscopic redshift (dimensionless).
    Returns: (bagpipes_components, drawn, noise_scale). The components are the delayed-tau model's plus
             veldisp [km/s]; drawn holds the model's drawn labels plus veldisp [km/s]; noise_scale is the
             error multiplier (dimensionless, 1 to 5).
    """
    cube_by_name = dict(zip(SPECTRAL_DIMS, cube))
    unit_vector = training_unit_vector(cube_by_name, SPECTRAL_MODEL)
    bagpipes_components, drawn = MODELS[SPECTRAL_MODEL].components(unit_vector, fixed_redshift=redshift)
''', '''def spectral_components(cube, redshift, sfh_model):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit point of an SFH model.

    The SFH, metallicity, dust and nebular part is the model's training prior (MODELS[sfh_model]) at the fixed
    spectroscopic redshift, which also sets the oldest allowed age (delayed-tau) or turnover (double power
    law); veldisp and noise_scale are the spectral nuisance parameters, the same for every model.
    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as spectral_dims(sfh_model);
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw".
    Returns: (bagpipes_components, drawn, noise_scale). The components are the model's plus veldisp [km/s];
             drawn holds the model's drawn labels plus veldisp [km/s]; noise_scale is the error multiplier
             (dimensionless, 1 to 5).
    """
    cube_by_name = dict(zip(spectral_dims(sfh_model), cube))
    unit_vector = training_unit_vector(cube_by_name, sfh_model)
    bagpipes_components, drawn = MODELS[sfh_model].components(unit_vector, fixed_redshift=redshift)
''')

# 2. fit_spectrum
replace('''def fit_spectrum(spectrum, redshift, output_basename, max_iterations, n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one prepared spectrum under the training prior at a fixed redshift.

    The sampler works directly on the unit cube (SPECTRAL_DIMS): spectral_components maps each point
    to physical values, exactly as the training set does. output_basename must be short (MultiNest
    limit of 100 characters), e.g. "spectral_fits/<id>/".
    Inputs: spectrum: dict with wavelength [Angstrom, observed], flux and ivar, as from prepare_spectrum;
            redshift: spectroscopic redshift (dimensionless); output_basename: directory prefix for the
            MultiNest files (str); max_iterations: MultiNest iteration limit (count);
            n_live_points: live points (count).
    Returns: dict with samples (unit-cube points, shape (n_samples, len(SPECTRAL_DIMS))), log_evidence
             (dimensionless ln Z), log_evidence_error (MultiNest's uncertainty of ln Z, sqrt(H / n_live)),
             n_iterations (count, from multinest_iterations), converged (bool: the run certainly stopped
             before max_iterations, judged by the upper bound of the count) and run_time [s].
    """
    # 1. Bagpipes model built once; each likelihood call updates it in place, which is much faster
    midpoint_cube = np.full(len(SPECTRAL_DIMS), MIDPOINT_DRAW)
    starting_components, _, _ = spectral_components(midpoint_cube, redshift)
    model = pipes.model_galaxy(starting_components, spec_wavs=spectrum["wavelength"])

    # 2. Likelihood of one unit-cube point: model spectrum vs data, with the S/N-scaling noise multiplier
    def log_likelihood(cube):
        """Log-likelihood of the spectrum for one unit-cube point.

        Inputs: cube: unit-cube point (dimensionless, SPECTRAL_DIMS order).
        Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
        """
        bagpipes_components, _, noise_scale = spectral_components(cube, redshift)
        model.update(bagpipes_components)
        model_flux = model.spectrum[:, 1]
        return safe_spectral_log_likelihood(model_flux, spectrum["flux"], spectrum["ivar"], noise_scale)

    # 3. MultiNest over the unit cube
    return run_multinest(log_likelihood, len(SPECTRAL_DIMS), output_basename, max_iterations, n_live_points)
''', '''def fit_spectrum(spectrum, redshift, sfh_model, output_basename, max_iterations, n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one prepared spectrum under an SFH model's training prior at a fixed redshift.

    The sampler works directly on the unit cube (spectral_dims(sfh_model)): spectral_components maps each
    point to physical values, exactly as the training set does. Every SFH model gets the same likelihood
    (calibration polynomial, noise scale) and MultiNest settings, so their evidences can be compared.
    output_basename must be short (MultiNest limit of 100 characters), e.g. "spectral_fits/<id>/".
    Inputs: spectrum: dict with wavelength [Angstrom, observed], flux and ivar, as from prepare_spectrum;
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw";
            output_basename: directory prefix for the MultiNest files (str); max_iterations: MultiNest
            iteration limit (count); n_live_points: live points (count).
    Returns: dict with samples (unit-cube points, shape (n_samples, len(spectral_dims(sfh_model)))),
             log_evidence (dimensionless ln Z), log_evidence_error (MultiNest's uncertainty of ln Z,
             sqrt(H / n_live)), n_iterations (count, from multinest_iterations), converged (bool: the run
             certainly stopped before max_iterations, judged by the upper bound of the count) and
             run_time [s].
    """
    n_dims = len(spectral_dims(sfh_model))

    # 1. Bagpipes model built once; each likelihood call updates it in place, which is much faster
    midpoint_cube = np.full(n_dims, MIDPOINT_DRAW)
    starting_components, _, _ = spectral_components(midpoint_cube, redshift, sfh_model)
    model = pipes.model_galaxy(starting_components, spec_wavs=spectrum["wavelength"])

    # 2. Likelihood of one unit-cube point: model spectrum vs data, with the S/N-scaling noise multiplier
    def log_likelihood(cube):
        """Log-likelihood of the spectrum for one unit-cube point.

        Inputs: cube: unit-cube point (dimensionless, spectral_dims(sfh_model) order).
        Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
        """
        bagpipes_components, _, noise_scale = spectral_components(cube, redshift, sfh_model)
        # PyMultiNest calls sys.exit if the likelihood raises: a model Bagpipes cannot compute is impossible
        # (as in photometric_point_log_likelihood; a point that computes gets the same value as before)
        try:
            model.update(bagpipes_components)
            model_flux = model.spectrum[:, 1]
        except Exception:
            return LOG_ZERO
        return safe_spectral_log_likelihood(model_flux, spectrum["flux"], spectrum["ivar"], noise_scale)

    # 3. MultiNest over the unit cube
    return run_multinest(log_likelihood, n_dims, output_basename, max_iterations, n_live_points)
''')

# 3. posterior_labels
replace('''def posterior_labels(unit_cube_samples, redshift):
    """Labels of each spectral-fit posterior sample, from a Bagpipes model that also gives S-PLUS photometry.

    Keeps at most MAX_POSTERIOR_SAMPLES evenly spaced samples, because building labels is the slow part.
    model_r is the model r magnitude for 10^10 M_sun formed, used later to rescale mass-like labels.
    Inputs: unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, len(SPECTRAL_DIMS)));
            redshift: spectroscopic redshift (dimensionless).
    Returns: DataFrame with one row per kept sample: the delayed-tau training labels (mass-like ones for
             10^10 M_sun formed), then veldisp [km/s], noise_scale (dimensionless) and model_r [mag].
    """''', '''def posterior_labels(unit_cube_samples, redshift, sfh_model):
    """Labels of each spectral-fit posterior sample, from a Bagpipes model that also gives S-PLUS photometry.

    Keeps at most MAX_POSTERIOR_SAMPLES evenly spaced samples, because building labels is the slow part.
    model_r is the model r magnitude for 10^10 M_sun formed, used later to rescale mass-like labels.
    Inputs: unit_cube_samples: unit-cube points (dimensionless, shape (n_samples,
            len(spectral_dims(sfh_model)))); redshift: spectroscopic redshift (dimensionless);
            sfh_model: the fit's model, "delayed" or "dblplaw".
    Returns: DataFrame with one row per kept sample: the model's training labels (mass-like ones for
             10^10 M_sun formed), then veldisp [km/s], noise_scale (dimensionless) and model_r [mag].
    """''')
replace('''        bagpipes_components, drawn, noise_scale = spectral_components(cube, redshift)
        model = updated_photometry_model(model, bagpipes_components)
        drawn_with_noise = dict(drawn)
        drawn_with_noise["noise_scale"] = noise_scale
        rows.append(sample_label_row(model, bagpipes_components, drawn_with_noise, SPECTRAL_MODEL))''',
        '''        bagpipes_components, drawn, noise_scale = spectral_components(cube, redshift, sfh_model)
        model = updated_photometry_model(model, bagpipes_components)
        drawn_with_noise = dict(drawn)
        drawn_with_noise["noise_scale"] = noise_scale
        rows.append(sample_label_row(model, bagpipes_components, drawn_with_noise, sfh_model))''')

# 4. make_mock_spectrum
replace('''def make_mock_spectrum(cube, redshift, template_spectrum, rng):
    """Mock spectrum from a known unit-cube point, with the wavelengths, mask and errors of a real one.

    The model is scaled to the template's median flux so the mock has the template's S/N. The
    noise_scale coordinate is set to UNIT_NOISE_CUBE_VALUE (error multiplier 1) because the mock
    noise is the template's own.
    Inputs: cube: unit-cube point (dimensionless, SPECTRAL_DIMS order); redshift: dimensionless;
            template_spectrum: dict from prepare_spectrum (wavelength [Angstrom], flux, ivar);
            rng: numpy Generator.
    Returns: (spectrum dict with wavelength, flux and ivar, truth labels as a pandas Series with the
             same labels as posterior_labels).
    """
    # 1. Known model on the template's wavelength grid (copied, so the caller's point is not changed)
    cube = np.array(cube, dtype=float)
    cube[SPECTRAL_DIMS.index("noise_scale")] = UNIT_NOISE_CUBE_VALUE
    bagpipes_components, _, _ = spectral_components(cube, redshift)''', '''def make_mock_spectrum(cube, redshift, template_spectrum, rng, sfh_model):
    """Mock spectrum of a known point of an SFH model, with the wavelengths, mask and errors of a real one.

    The model is scaled to the template's median flux so the mock has the template's S/N. The
    noise_scale coordinate is set to UNIT_NOISE_CUBE_VALUE (error multiplier 1) because the mock
    noise is the template's own.
    Inputs: cube: unit-cube point (dimensionless, spectral_dims(sfh_model) order); redshift: dimensionless;
            template_spectrum: dict from prepare_spectrum (wavelength [Angstrom], flux, ivar);
            rng: numpy Generator (only its normal() is used, once); sfh_model: "delayed" or "dblplaw".
    Returns: (spectrum dict with wavelength, flux and ivar, truth labels as a pandas Series with the
             same labels as posterior_labels of sfh_model).
    """
    # 1. Known model on the template's wavelength grid (copied, so the caller's point is not changed)
    cube = np.array(cube, dtype=float)
    cube[spectral_dims(sfh_model).index("noise_scale")] = UNIT_NOISE_CUBE_VALUE
    bagpipes_components, _, _ = spectral_components(cube, redshift, sfh_model)''')
replace('''    truth = posterior_labels(cube[np.newaxis, :], redshift).iloc[0]''',
        '''    truth = posterior_labels(cube[np.newaxis, :], redshift, sfh_model).iloc[0]''')

# 5. save_fit docstring
replace('''            fit: len(SPECTRAL_DIMS) or the model's photometric_dims);''',
        '''            fit: len(spectral_dims(sfh_model)) or the model's photometric_dims);''')

# 6. fit_spectrum_job
replace('''def fit_spectrum_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a spectrum job and label its posterior samples.

    Inputs: job: dict with spectrum (as from prepare_spectrum) and redshift (dimensionless); multinest_basename:
            directory prefix for the MultiNest files (str); max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_spectrum, posterior DataFrame from posterior_labels).
    Raises: ValueError for a job of a model other than SPECTRAL_MODEL (spectral fits are delayed-tau only).
    """
    if job_model(job) != SPECTRAL_MODEL:
        raise ValueError(f"spectral fits are {SPECTRAL_MODEL} only, not {job_model(job)} "
                         "(DPL spec, decision 1)")
    fit = fit_spectrum(job["spectrum"], job["redshift"], multinest_basename, max_iterations, n_live_points)
    posterior = posterior_labels(fit["samples"], job["redshift"])
    return fit, posterior''', '''def fit_spectrum_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a spectrum job with its SFH model and label its posterior samples.

    Inputs: job: dict with spectrum (as from prepare_spectrum), redshift (dimensionless) and sfh_model
            (job_model: "delayed" for a job without one, as every job exported before the DPL spectral fits);
            multinest_basename: directory prefix for the MultiNest files (str); max_iterations, n_live_points:
            counts.
    Returns: (fit dict from fit_spectrum, posterior DataFrame from posterior_labels).
    """
    sfh_model = job_model(job)
    fit = fit_spectrum(job["spectrum"], job["redshift"], sfh_model, multinest_basename, max_iterations,
                       n_live_points)
    posterior = posterior_labels(fit["samples"], job["redshift"], sfh_model)
    return fit, posterior''')

# 7. Reference-CSV prefixes of the spectral fits
replace('''PHOTOMETRIC_PREFIXES = {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}
''', '''PHOTOMETRIC_PREFIXES = {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}
# ... and of each model's spectral fits: <prefix><label>_p16, <prefix>converged, ...; delayed-tau keeps the
# plain names of the first spectral fits (only its ln Z columns are named spec_log_evidence(_error))
SPECTRAL_PREFIXES = {"delayed": "", "dblplaw": "spec_dblplaw_"}
''')

# 8. Calibration mocks per model, and the DPL spectral jobs
replace('''def draw_calibration_cubes(n_mocks, rng):
    """Unit-cube points of the calibration mocks: uniform over every SPECTRAL_DIMS coordinate.

    Uniform on the unit cube is the training prior at the template's redshift (log-uniform age, tau and Z).
    draw_calibration_mocks calls this first on its rng, so the cubes can be redrawn from the same seed.
    Inputs: n_mocks: number of points (count); rng: numpy Generator.
    Returns: array of shape (n_mocks, len(SPECTRAL_DIMS)), values in [0, 1) (dimensionless).
    """
    return rng.random((n_mocks, len(SPECTRAL_DIMS)))


def draw_calibration_mocks(n_mocks, templates, rng):
    """Mock spectra with known labels, each on the wavelengths, mask and errors of a random real galaxy.

    Inputs: n_mocks: number of mocks (count); templates: list of (spectrum dict from prepare_spectrum,
            redshift (dimensionless)); rng: numpy Generator.
    Returns: (jobs, truths): jobs is a list of dicts as fit_and_save takes (job_id calib_000, calib_001, ...);
             truths maps each job_id to its true labels (Series, as posterior_labels columns: the
             spectral model's training_labels, veldisp, noise_scale and model_r).
    """
    cubes = draw_calibration_cubes(n_mocks, rng)
    jobs = []
    truths = {}
    for mock_index, cube in enumerate(cubes):
        job_id = f"{CALIBRATION_ID_PREFIX}{mock_index:03d}"
        template_spectrum, redshift = templates[rng.integers(len(templates))]
        spectrum, truth = make_mock_spectrum(cube, redshift, template_spectrum, rng)
        jobs.append({"job_id": job_id, "spectrum": spectrum, "redshift": redshift})
        truths[job_id] = truth
    return jobs, truths
''', '''def draw_calibration_cubes(n_mocks, rng, sfh_model):
    """Unit-cube points of an SFH model's spectral calibration mocks: uniform over its spectral_dims.

    Uniform on the unit cube is the model's training prior at the template's redshift (log-uniform age, tau
    and Z for delayed-tau; log-uniform turnover lookback, alpha, beta and Z for the double power law).
    draw_calibration_mocks calls this first on its rng, so the cubes can be redrawn from the same seed.
    Inputs: n_mocks: number of points (count); rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: array of shape (n_mocks, len(spectral_dims(sfh_model))), values in [0, 1) (dimensionless).
    """
    return rng.random((n_mocks, len(spectral_dims(sfh_model))))


def draw_calibration_mocks(n_mocks, templates, rng, sfh_model):
    """Mock spectra of an SFH model with known labels, each with the wavelengths, mask and errors of a galaxy.

    Inputs: n_mocks: number of mocks (count); templates: list of (spectrum dict from prepare_spectrum,
            redshift (dimensionless)); rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: (jobs, truths): jobs is a list of spectrum job dicts as fit_and_save takes (job_id calib_000,
             calib_001, ..., kind SPECTRUM_JOB, sfh_model, spectrum, redshift); truths maps each job_id to its
             true labels (Series, as posterior_labels columns: the model's training_labels, veldisp,
             noise_scale and model_r).
    """
    cubes = draw_calibration_cubes(n_mocks, rng, sfh_model)
    jobs = []
    truths = {}
    for mock_index, cube in enumerate(cubes):
        job_id = f"{CALIBRATION_ID_PREFIX}{mock_index:03d}"
        template_spectrum, redshift = templates[rng.integers(len(templates))]
        spectrum, truth = make_mock_spectrum(cube, redshift, template_spectrum, rng, sfh_model)
        jobs.append({"job_id": job_id, "kind": SPECTRUM_JOB, "sfh_model": sfh_model, "spectrum": spectrum,
                     "redshift": redshift})
        truths[job_id] = truth
    return jobs, truths


def dblplaw_spectral_jobs(delayed_jobs):
    """The galaxies' prepared spectra as spectral jobs of the double power law (DPL spec, amendment 2).

    Why the calibration mocks are left out: they were drawn from the delayed-tau prior, so their truths have
    no DPL labels; the DPL fitter gets its own mocks (draw_calibration_mocks with "dblplaw").
    Inputs: delayed_jobs: delayed-tau spectral jobs, as notebook 2 exports them to spectral_fits/jobs/ (a job
            without kind is a spectrum job, one without sfh_model is delayed-tau): galaxies, and possibly
            calibration mocks (job ids starting with CALIBRATION_ID_PREFIX).
    Returns: list of new job dicts, one per galaxy, in input order: the same job_id, spectrum and redshift,
             with kind SPECTRUM_JOB and sfh_model "dblplaw". The input jobs are not changed; their spectrum
             dicts are shared, not copied (nothing writes to them).
    Raises: ValueError for a job that is not a delayed-tau spectrum job.
    """
    jobs = []
    for delayed_job in delayed_jobs:
        if job_kind(delayed_job) != SPECTRUM_JOB or job_model(delayed_job) != "delayed":
            raise ValueError(f"{delayed_job['job_id']} is not a delayed-tau spectrum job "
                             f"({job_kind(delayed_job)}, {job_model(delayed_job)})")
        if str(delayed_job["job_id"]).startswith(CALIBRATION_ID_PREFIX):
            continue
        job = dict(delayed_job)
        job["kind"] = SPECTRUM_JOB
        job["sfh_model"] = "dblplaw"
        jobs.append(job)
    return jobs
''')

# 9. Fit folders
replace('''FIT_FOLDER_MODELS = {"spectral_fits": "delayed", "photometric_fits": "delayed",
                     "photometric_fits_dblplaw": "dblplaw"}''', '''FIT_FOLDER_MODELS = {"spectral_fits": "delayed", "photometric_fits": "delayed",
                     "photometric_fits_dblplaw": "dblplaw", "spectral_fits_dblplaw": "dblplaw"}''')
replace('''    is judged against its own SFH model (FIT_FOLDER_MODELS): a DPL fit stays in photometric_fits_dblplaw and
    leaves photometric_fits, and the reverse.''', '''    is judged against its own SFH model (FIT_FOLDER_MODELS): a DPL fit stays in photometric_fits_dblplaw or
    spectral_fits_dblplaw and leaves photometric_fits or spectral_fits, and the reverse.''')

open(PATH, "w").write(text)
print("synthetic_model.py patched")
```

Expected output: `synthetic_model.py patched`.

- [ ] **Step 4: `spectral_validation.py` names its model; `run_spectral_fits.py` names the new folder**

Save as `$SCRATCH/task13_validation.py` and run `../.venv/bin/python "$SCRATCH/task13_validation.py" .`:

```python
"""Task 13, step 4: spectral_validation.py names its model, run_spectral_fits.py names the new folder."""
import sys

root = sys.argv[1]


def patch(path, pairs):
    text = open(path).read()
    for old, new in pairs:
        assert text.count(old) == 1, (path, text.count(old), old[:80])
        text = text.replace(old, new)
    open(path, "w").write(text)
    print(path, "patched")


patch(root + "/spectral_validation.py", [
    ('''N_PREDICTIVE_SAMPLES = 100     # posterior samples per galaxy for the model indices
MAX_Z_SCORE = 2.0              # an index agrees when |z| < MAX_Z_SCORE
''', '''N_PREDICTIVE_SAMPLES = 100     # posterior samples per galaxy for the model indices
MAX_Z_SCORE = 2.0              # an index agrees when |z| < MAX_Z_SCORE
# The rulers read the delayed-tau spectral fits: the validation stays on them (DPL spec, amendment 2)
RULER_SFH_MODEL = "delayed"
'''),
    ('''    components, _, _ = sm.spectral_components(cube, redshift)
    model.update(components)''', '''    components, _, _ = sm.spectral_components(cube, redshift, RULER_SFH_MODEL)
    model.update(components)'''),
    ('''    Inputs: unit-cube posterior samples (n, SPECTRAL_DIMS); redshift; the galaxy's native observed air
            wavelength grid [Angstrom] and its flux/ivar on that grid (restricted to sm.OBSERVED_RANGE by the caller).''',
     '''    Inputs: unit-cube posterior samples of a delayed-tau spectral fit (n, sm.spectral_dims(RULER_SFH_MODEL));
            redshift; the galaxy's native observed air wavelength grid [Angstrom] and its flux/ivar on that
            grid (restricted to sm.OBSERVED_RANGE by the caller).'''),
    ('''    first_components, _, _ = sm.spectral_components(kept[0], redshift)''',
     '''    first_components, _, _ = sm.spectral_components(kept[0], redshift, RULER_SFH_MODEL)'''),
    ('''    Inputs: unit_cube_samples: unit-cube posterior samples (dimensionless, shape (n, len(SPECTRAL_DIMS)));''',
     '''    Inputs: unit_cube_samples: unit-cube posterior samples of a delayed-tau spectral fit (dimensionless,
            shape (n, len(sm.spectral_dims(RULER_SFH_MODEL))));'''),
    ('''        components, drawn, _ = sm.spectral_components(cube, redshift)''',
     '''        components, drawn, _ = sm.spectral_components(cube, redshift, RULER_SFH_MODEL)'''),
])

patch(root + "/run_spectral_fits.py", [
    ('''its true model (truth_model). The two kinds share job ids (the specobjid), so they live in separate folders:
spectral_fits/, photometric_fits/ (delayed-tau) and photometric_fits_dblplaw/ (double power law); pass e.g.
--jobs-dir photometric_fits_dblplaw/jobs --output-dir photometric_fits_dblplaw.''',
     '''its true model (truth_model). The kinds and models share job ids (the specobjid, calib_<k>), so each has its
own folder: spectral_fits/ and photometric_fits/ (delayed-tau), spectral_fits_dblplaw/ and
photometric_fits_dblplaw/ (double power law); pass e.g. --jobs-dir spectral_fits_dblplaw/jobs --output-dir
spectral_fits_dblplaw.'''),
])
```

Expected output: `./spectral_validation.py patched` and `./run_spectral_fits.py patched`.

- [ ] **Step 5: Notebook 2 passes `"delayed"` to the changed functions (names only)**

Save as `$SCRATCH/task13_notebook.py` and run `../.venv/bin/python "$SCRATCH/task13_notebook.py"`:

```python
"""Task 13, step 5: notebook 2 calls the changed spectral functions with the delayed-tau model (names only)."""
import json

NOTEBOOK = "spectroscopic_reference.ipynb"
with open(NOTEBOOK) as notebook_file:
    notebook = json.load(notebook_file)

# (old text, new text, expected number of occurrences in the notebook's code cells)
REPLACEMENTS = [
    ("sm.N_CALIBRATION_MOCKS, templates, np.random.default_rng(SEED))",
     "sm.N_CALIBRATION_MOCKS, templates, np.random.default_rng(SEED), \"delayed\")", 1),
    ("sm.MODELS[sm.SPECTRAL_MODEL]", "sm.MODELS[\"delayed\"]", 2),
    ("sm.draw_calibration_cubes(len(calibration_jobs), np.random.default_rng(SEED))",
     "sm.draw_calibration_cubes(len(calibration_jobs), np.random.default_rng(SEED), \"delayed\")", 1),
    ("dict(zip(sm.SPECTRAL_DIMS, spectral_cubes[mock_index]))",
     "dict(zip(sm.spectral_dims(\"delayed\"), spectral_cubes[mock_index]))", 1),
]

for old, new, expected_count in REPLACEMENTS:
    found = 0
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if old in source:
            found += source.count(old)
            cell["source"] = source.replace(old, new).splitlines(keepends=True)
    assert found == expected_count, (old, found)

with open(NOTEBOOK, "w") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print("notebook 2 updated")
```

Expected output: `notebook 2 updated`. The cells keep their outputs: only names change, and the PC run of Task 14 replaces them. Check it (never execute the notebook on the laptop):

```bash
../.venv/bin/python -c "
import nbformat
nb = nbformat.read('spectroscopic_reference.ipynb', as_version=4); nbformat.validate(nb)
for i, c in enumerate(nb.cells):
    if c.cell_type == 'code': compile(c.source, f'cell {i}', 'exec')
print('ok', len(nb.cells))"
grep -rn "SPECTRAL_DIMS\|SPECTRAL_MODEL" --include=*.py --include=*.ipynb . | grep -v "^./pipes/" | grep -v "tests/"
```

Expected: `ok 36`, and the `grep` prints nothing. The tests still use the names: `hasattr(sm, "SPECTRAL_DIMS")` in the new dims test, and the test constants `SPECTRAL_DIMS` / `N_DELAYED_SPECTRAL_DIMS`.

- [ ] **Step 6: Run the tests and watch them pass**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py tests/test_spectral_validation.py tests/test_notebooks.py -q`

Expected: `168 passed` (156 before this task: 12 new tests, and 2 replaced). The run takes about 2.5 minutes. The DPL spectral smoke run takes about 11 s, and the two regression tests against the stored delayed-τ fit and mocks about 1 s. On a machine without `spectral_fits/`, those two skip (`166 passed, 2 skipped`). The DPL smoke run skips when pymultinest cannot import.

Then run the whole suite: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -q`. Expected: `282 passed` (270 + 12), with the Bagpipes overflow warnings as before. It is `281 passed, 1 skipped` when the short MultiNest run of `test_recomputed_evidence_error_agrees_with_multinest` does not converge (it skips itself then).

- [ ] **Step 7: Remote scripts, README, the DPL skip list and `.gitignore`**

Save as `$SCRATCH/task13_remote.py` and run it from the repository root: `.venv/bin/python "$SCRATCH/task13_remote.py" .`

```python
"""Task 13, step 7: remote scripts, README, skip list and .gitignore. Argument: the repository root."""
import sys

root = sys.argv[1]


def patch(path, pairs):
    text = open(path).read()
    for old, new in pairs:
        assert text.count(old) == 1, (path, text.count(old), old[:80])
        text = text.replace(old, new)
    open(path, "w").write(text)
    print(path, "patched")


patch(root + "/Synthetic_domain/remote/submit_fits.sh", [
    ('''# FITS_SUBDIR picks the set of fits: spectral_fits (the spectral fits), photometric_fits (the delayed-tau
# photometric fits of spec amendment A) or photometric_fits_dblplaw (the double-power-law photometric fits and
# the cross-fits of the delayed-tau mocks). They share job ids (specobjids, calib_<k>), so their jobs and
# results never mix:
#     FITS_SUBDIR=photometric_fits_dblplaw Synthetic_domain/remote/submit_fits.sh''',
     '''# FITS_SUBDIR picks the set of fits: spectral_fits (the delayed-tau spectral fits), spectral_fits_dblplaw (the
# double-power-law spectral fits), photometric_fits (the delayed-tau photometric fits of spec amendment A) or
# photometric_fits_dblplaw (the double-power-law photometric fits and the cross-fits of the delayed-tau mocks).
# They share job ids (specobjids, calib_<k>), so their jobs and results never mix:
#     FITS_SUBDIR=spectral_fits_dblplaw Synthetic_domain/remote/submit_fits.sh'''),
    ('''#   SKIP_LIST        file of job ids left out of the bundle   [default skipped_fit_jobs.txt for spectral_fits, none
#                    otherwise: the photometric calibration mocks reuse the ids calib_<k> of the spectral ones]''',
     '''#   SKIP_LIST        file of job ids left out of the bundle   [default skipped_fit_jobs.txt for spectral_fits,
#                    skipped_fit_jobs_dblplaw.txt for spectral_fits_dblplaw, none for the photometric folders:
#                    every folder's calibration mocks reuse the ids calib_<k>, so each needs its own list]'''),
    ('''if [[ "$FITS_SUBDIR" == "spectral_fits" ]]; then
    SKIP_LIST="${SKIP_LIST-$SYNTHETIC_DIR/skipped_fit_jobs.txt}"
else''', '''if [[ "$FITS_SUBDIR" == "spectral_fits" ]]; then
    SKIP_LIST="${SKIP_LIST-$SYNTHETIC_DIR/skipped_fit_jobs.txt}"
elif [[ "$FITS_SUBDIR" == "spectral_fits_dblplaw" ]]; then
    SKIP_LIST="${SKIP_LIST-$SYNTHETIC_DIR/skipped_fit_jobs_dblplaw.txt}"
else'''),
])

patch(root + "/Synthetic_domain/remote/fit_job.sh", [
    ('''#   FITS_SUBDIR      folder of the jobs and results: spectral_fits, photometric_fits or
#                    photometric_fits_dblplaw   [default spectral_fits]''',
     '''#   FITS_SUBDIR      folder of the jobs and results: spectral_fits, spectral_fits_dblplaw, photometric_fits
#                    or photometric_fits_dblplaw   [default spectral_fits]'''),
])

patch(root + "/Synthetic_domain/remote/fetch_fits.sh", [
    ('''# The fits folder is the one the job recorded in out/fits_subdir.txt (spectral_fits, photometric_fits or
# photometric_fits_dblplaw). Jobs made before it was recorded are spectral fits; for them FITS_SUBDIR may name
# the folder [default spectral_fits].''',
     '''# The fits folder is the one the job recorded in out/fits_subdir.txt (spectral_fits, spectral_fits_dblplaw,
# photometric_fits or photometric_fits_dblplaw). Jobs made before it was recorded are delayed-tau spectral fits;
# for them FITS_SUBDIR may name the folder [default spectral_fits].'''),
])

patch(root + "/Synthetic_domain/remote/README.md", [
    ('''The DPL mocks' cross-fits (fitted with delayed-tau) are exported to `photometric_fits/jobs/`, so a
`FITS_SUBDIR=photometric_fits` run after it fits only those: the finished delayed-tau fits travel in its bundle.
''', '''The DPL mocks' cross-fits (fitted with delayed-tau) are exported to `photometric_fits/jobs/`, so a
`FITS_SUBDIR=photometric_fits` run after it fits only those: the finished delayed-tau fits travel in its bundle.

The double-power-law spectral fits (section 6d of notebook 2: every galaxy with a prepared spectrum, refitted
with the DPL prior, and 33 DPL calibration mocks, about 1400 fits) have their own folder too:

    FITS_SUBDIR=spectral_fits_dblplaw Synthetic_domain/remote/submit_fits.sh

They take about 30 h, so two 24 h jobs: fetch the first (`fetch_fits.sh ID`) and submit again; the second fits
only the rest, because the finished fits travel in its bundle.
'''),
    ('''`skipped_fit_jobs.txt` applies to the spectral fits only (`SKIP_LIST` overrides it). It is empty for the delayed-tau
model: the ids it listed named 16-bin mocks, and the new mocks reuse the `calib_<k>` ids. `BUNDLE_ONLY=<folder>`''',
     '''`skipped_fit_jobs.txt` applies to the delayed-tau spectral fits and `skipped_fit_jobs_dblplaw.txt` to the DPL
spectral fits (`SKIP_LIST` overrides either; the photometric folders have none). Each folder needs its own list,
because every folder's mocks reuse the `calib_<k>` ids. Both are empty: the ids `skipped_fit_jobs.txt` listed
named 16-bin mocks. `BUNDLE_ONLY=<folder>`'''),
])

patch(root + "/.gitignore", [
    ('''Synthetic_domain/photometric_fits_dblplaw/
''', '''Synthetic_domain/photometric_fits_dblplaw/
Synthetic_domain/spectral_fits_dblplaw/
'''),
])

with open(root + "/Synthetic_domain/skipped_fit_jobs_dblplaw.txt", "w") as skip_file:
    skip_file.write(
        "# DPL spectral fit jobs deliberately left unfitted, one job_id per line ('#' starts a comment).\n"
        "# remote/submit_fits.sh does not send them (FITS_SUBDIR=spectral_fits_dblplaw), and the DPL spectral\n"
        "# calibration test of spectroscopic_reference.ipynb (section 6d) gives its verdict on the other mocks,\n"
        "# listing these as skipped. The DPL mocks reuse the calib_<k> ids of the delayed-tau mocks, so the two\n"
        "# models keep separate lists (skipped_fit_jobs.txt is the delayed-tau one).\n")
print("skipped_fit_jobs_dblplaw.txt written")
```

Expected output: five `patched` lines and `skipped_fit_jobs_dblplaw.txt written`.

The resulting `SKIP_LIST` branch of `submit_fits.sh` (lines 45–51) is:

```bash
if [[ "$FITS_SUBDIR" == "spectral_fits" ]]; then
    SKIP_LIST="${SKIP_LIST-$SYNTHETIC_DIR/skipped_fit_jobs.txt}"
elif [[ "$FITS_SUBDIR" == "spectral_fits_dblplaw" ]]; then
    SKIP_LIST="${SKIP_LIST-$SYNTHETIC_DIR/skipped_fit_jobs_dblplaw.txt}"
else
    SKIP_LIST="${SKIP_LIST-}"
fi
```

- [ ] **Step 8: Check the scripts on a mirror, never in the repository's fit folders**

`submit_fits.sh` needs `$FITS_SUBDIR/jobs/` next to itself, and `spectral_fits_dblplaw/` must not be created in the repository by hand: notebook 2 makes it on the PC. So the check runs on a copy of the scripts and modules in the scratchpad. Run from the repository root:

```bash
bash -n Synthetic_domain/remote/submit_fits.sh && bash -n Synthetic_domain/remote/fit_job.sh \
  && bash -n Synthetic_domain/remote/fetch_fits.sh && echo syntax-ok
MIRROR="$SCRATCH/mirror/Synthetic_domain"
rm -rf "$SCRATCH/mirror" && mkdir -p "$MIRROR/spectral_fits_dblplaw/jobs"
cp -r Synthetic_domain/remote Synthetic_domain/*.py Synthetic_domain/skipped_fit_jobs*.txt "$MIRROR/"
cp -rL Synthetic_domain/splus-filters "$MIRROR/"
touch "$MIRROR/spectral_fits_dblplaw/jobs/calib_000.npz" "$MIRROR/spectral_fits_dblplaw/jobs/1105771463723476992.npz"
BUNDLE_ONLY="$SCRATCH/bundle_check" FITS_SUBDIR=spectral_fits_dblplaw "$MIRROR/remote/submit_fits.sh" | grep Bundle
tar -tzf "$SCRATCH/bundle_check/bundle.tar.gz" | grep -E 'spectral_fits_dblplaw/jobs/.|^\./(sfh_models|synthetic_model|run_spectral_fits)\.py$'
echo "calib_000" >> "$MIRROR/skipped_fit_jobs_dblplaw.txt"
BUNDLE_ONLY="$SCRATCH/bundle_check" FITS_SUBDIR=spectral_fits_dblplaw "$MIRROR/remote/submit_fits.sh" | grep Bundle
BUNDLE_ONLY="$SCRATCH/bundle_check" FITS_SUBDIR=spectral_fits_dblplaw SKIP_LIST= "$MIRROR/remote/submit_fits.sh" | grep Bundle
git check-ignore Synthetic_domain/spectral_fits_dblplaw/jobs/x.npz
git status --short Synthetic_domain/skipped_fit_jobs_dblplaw.txt
```

Expected:
- `syntax-ok`;
- `Bundle (spectral_fits_dblplaw): ..., with 0 finished fit(s), 0 skipped job(s) left out`;
- the three module paths and the two job files;
- after `calib_000` is added to the mirror's list, `1 skipped job(s) left out`: the default list is the DPL one;
- with `SKIP_LIST=` (empty), `0 skipped job(s)`;
- the ignored path printed, and `?? Synthetic_domain/skipped_fit_jobs_dblplaw.txt` (a new file, tracked like `skipped_fit_jobs.txt`).

The mirror only exists in the scratchpad. If `~/Programs/MultiNest` is missing on this machine, skip the `BUNDLE_ONLY` lines and rely on `bash -n`.

- [ ] **Step 9: Commit**

```bash
git add Synthetic_domain/synthetic_model.py Synthetic_domain/spectral_validation.py Synthetic_domain/run_spectral_fits.py \
        Synthetic_domain/spectroscopic_reference.ipynb Synthetic_domain/tests/test_synthetic_model.py \
        Synthetic_domain/tests/test_spectral_validation.py Synthetic_domain/remote/submit_fits.sh \
        Synthetic_domain/remote/fit_job.sh Synthetic_domain/remote/fetch_fits.sh Synthetic_domain/remote/README.md \
        Synthetic_domain/skipped_fit_jobs_dblplaw.txt .gitignore
git commit -m "feat: spectral fits of either SFH model; DPL spectral jobs and mocks; spectral_fits_dblplaw

spectral_dims, spectral_components, fit_spectrum, posterior_labels, make_mock_spectrum and the
calibration mocks take the SFH model; fit_spectrum_job fits a job with its own model, so a DPL
spectrum job is fitted (amendment 2 lifts decision 1). dblplaw_spectral_jobs copies the galaxy
jobs; SPECTRAL_PREFIXES names the reference columns; FIT_FOLDER_MODELS gains spectral_fits_dblplaw.
The delayed-tau spectral path is unchanged: regression tests relabel a stored fit and the stored
mocks exactly. The spectral likelihood returns LOG_ZERO where Bagpipes raises, as the photometric
one does; values elsewhere unchanged. submit_fits.sh skips the jobs of skipped_fit_jobs_dblplaw.txt
for the new folder.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 14: Notebook 2 exports and reads the DPL spectral fits (section 6d), with the calibration verdict and the spectral evidence of both models

**Why:** Amendment 2, items 2 and 3. The jobs of Task 13 must reach the PC, and their fits must reach `spectroscopic_reference.csv`, where `ml_model.ipynb` (Task 15) reads them. The new section 6d follows the delayed-τ path step by step:
- section 4 exports the jobs and draws the mocks;
- section 5 gives the calibration verdict;
- section 6 summarises the labels at the aperture brightness.

**What changes in the notebook:**
- **Section 5** (cell 14): the delayed-τ calibration test becomes the function `spectral_calibration(sfh_model, fits_directory, truths, skipped_jobs_file)`, and section 6d calls it for DPL. The delayed-τ verdict, table and printout stay the same. Two things differ: the lines name the model, and a missing skip file now counts as an empty list, with a printed `WARNING` that names the file to add to the notebook-2 inputs, instead of raising (a listed mock then counts as unfitted, so the verdict stays empty, never a false PASS). `converged_from_multinest` gains `fits_directory=FITS_DIRECTORY`, so the delayed-τ calls do not change.
- **Section 6** (cell 16): the delayed-τ rows gain `spec_log_evidence` and `spec_log_evidence_error`. The error is NaN for the 788 fits made before it was stored; the 604 fits of the running PC job store it. It is not recomputed: the photometric `sm.fit_evidence` recomputes errors with the photometric likelihood only, and the difference between two spectral ln Z values is expected to dwarf errors of about 0.1–0.3.
- **Section 6d** (new, four cells after the cross-fit table of section 6c):
  - export: `spectral_fits_dblplaw/calibration_mocks.joblib` (33 DPL mocks, drawn once with `DBLPLAW_SPECTRAL_MOCK_SEED = 3` on the templates of the selected galaxies, then reused) and `spectral_fits_dblplaw/jobs/*.npz` (the mocks and `sm.dblplaw_spectral_jobs(real_jobs)`). This happens on every run, without waiting for the split;
  - the calibration verdict (`spectral_calibration("dblplaw", ...)` with `skipped_fit_jobs_dblplaw.txt`);
  - the reading of the fits.
- **Markdown:** the title cell, sections 4, 5, 6 and 6c, and "Next pass" mention the new section. Two stale names in the markdown are fixed: `sm.CALIBRATION_LABELS` and `sm.TRAINING_LABELS`.
- **Superseding** (cell 10) loops over `sm.FIT_FOLDER_MODELS`, so it covers `spectral_fits_dblplaw/` through Task 13 without a code change. Only the section-4 markdown names the new folder.

**Interfaces:**
- **Consumes (Task 13):** `sm.spectral_dims`, `sm.draw_calibration_mocks(..., "dblplaw")`, `sm.dblplaw_spectral_jobs`, `sm.SPECTRAL_PREFIXES`, `sm.FIT_FOLDER_MODELS`, `sm.model_file`, `sm.fit_is_current`, `sm.SPECTRUM_JOB`, `skipped_fit_jobs_dblplaw.txt`. From the notebook itself: `real_jobs` (cell 11), `reference` with `r_aper` and `converged` (cell 16), `FITS_DIRECTORY` and `SKIPPED_JOBS_FILE` (cell 1).
- **Produces (notebook functions and variables):**
  - `converged_from_multinest(job_id, metadata, fits_directory=FITS_DIRECTORY)`, `skipped_job_ids(skipped_jobs_file) -> set`;
  - `spectral_calibration(sfh_model, fits_directory, truths, skipped_jobs_file) -> (verdict, table)`;
  - `dblplaw_spectral_summary(galaxies) -> DataFrame`;
  - `DBLPLAW_SPECTRAL_FITS_DIRECTORY = "spectral_fits_dblplaw"`, `DBLPLAW_SKIPPED_JOBS_FILE = "skipped_fit_jobs_dblplaw.txt"`, `DBLPLAW_SPECTRAL_LABELS`, `dblplaw_spectral_fitter_calibrated`.
- **Produces (files, written by notebook 2 on the PC):** `spectral_fits_dblplaw/jobs/*.npz` (33 `calib_<k>` DPL mocks, and one job per prepared spectrum, named by its specobjid as in `spectral_fits/jobs/`) and `spectral_fits_dblplaw/calibration_mocks.joblib`.
- **Produces (`spectroscopic_reference.csv` columns):**
  - delayed-τ: `spec_log_evidence`, `spec_log_evidence_error`; the existing plain columns (`<label>_p16/_median/_p84`, `converged`, `fitter_calibrated`, `split`) do not change;
  - DPL: `spec_dblplaw_<label>_p16/_median/_p84` for `[*sm.MODELS["dblplaw"].training_labels, "veldisp", "noise_scale"]`. That is the label set the delayed-τ spectral columns have (`LABEL_COLUMNS`), with mass-like labels at the aperture r brightness;
  - DPL: `spec_dblplaw_converged`, `spec_dblplaw_log_evidence`, `spec_dblplaw_log_evidence_error`, and `spec_dblplaw_fitter_calibrated` on every row;
  - galaxies without a DPL fit have NaN. Like the delayed-τ spectral columns, the DPL ones cover both halves of the split. The split stays the delayed-τ one, so `ml_model.ipynb` must use the test half only under `RUN_LOCKED_TEST`.

- [ ] **Step 1: Edit the notebook**

Save as `$SCRATCH/task14_notebook.py` and run `../.venv/bin/python "$SCRATCH/task14_notebook.py"` from `Synthetic_domain/`:

```python
"""Task 14: notebook 2 exports and reads the DPL spectral fits (section 6d), gives the DPL spectral calibration
verdict and the spectral evidence columns of both models."""
import json
import uuid

NOTEBOOK = "spectroscopic_reference.ipynb"
with open(NOTEBOOK) as notebook_file:
    notebook = json.load(notebook_file)
cells = notebook["cells"]


def source_of(cell):
    return "".join(cell["source"])


def find_cell(first_line):
    """Index of the one cell whose source starts with first_line."""
    matches = [index for index, cell in enumerate(cells) if source_of(cell).startswith(first_line)]
    assert len(matches) == 1, (first_line, matches)
    return matches[0]


def set_source(cell, text):
    cell["source"] = text.splitlines(keepends=True)
    if cell["cell_type"] == "code":
        cell["outputs"] = []
        cell["execution_count"] = None


def replace_in_cell(index, old, new):
    """Replace the one occurrence of old in cell index (markdown keeps its rendering; code loses its outputs)."""
    text = source_of(cells[index])
    assert text.count(old) == 1, (index, old[:60], text.count(old))
    set_source(cells[index], text.replace(old, new))


def new_cell(cell_type, text):
    cell = {"cell_type": cell_type, "id": uuid.uuid4().hex[:8], "metadata": {}, "source": []}
    if cell_type == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    set_source(cell, text)
    return cell


# ---------------------------------------------------------------------------------------------------------------
# Section 5: the calibration test as a function of the model, so section 6d runs it for the DPL fits
CALIBRATION_CELL = r'''def converged_from_multinest(job_id, metadata, fits_directory=FITS_DIRECTORY):
    """Whether a fit certainly stopped before its iteration limit, judged from its MultiNest folder.

    Uses the upper bound of the outer-loop iteration count (sm.multinest_iteration_bounds), so a capped
    run is never called converged. The limit is the one stored in the fit's metadata. A fit saved before
    the limit was stored (or without its MultiNest folder) keeps the flag stored when it ran: that flag
    came from the ev.dat row count, which overcounts iterations, so a stored True is already safe and a
    stored False is never upgraded against a limit that may differ from the one the fit used.
    Inputs: job_id: fit identifier, the name of its .h5 file and MultiNest folder (str); metadata: dict from
            sm.load_fit; fits_directory: the fit's folder (str; by default the delayed-tau spectral fits).
    Returns: True if the run converged, else False (bool).
    """
    iteration_limit = metadata.get("max_iterations")
    multinest_basename = os.path.join(fits_directory, job_id) + "/"
    if iteration_limit is None or not os.path.exists(multinest_basename + "ev.dat"):
        return bool(metadata["converged"])
    return sm.multinest_iteration_bounds(multinest_basename)[1] < iteration_limit


def skipped_job_ids(skipped_jobs_file):
    """Job ids listed in a skip file: one per line, '#' starts a comment.

    Why a missing file lists none: a skipped mock then counts as unfitted, so the verdict stays empty instead
    of passing on fewer mocks.
    Inputs: skipped_jobs_file: path of the list (str).
    Returns: set of job ids (str); empty, with a printed warning, when the file does not exist.
    """
    if not os.path.exists(skipped_jobs_file):
        print(f"WARNING: {skipped_jobs_file} not found (add it to the notebook-2 inputs): no mock counts "
              "as skipped")
        return set()
    with open(skipped_jobs_file) as skipped_file:
        return {line.split("#")[0].strip() for line in skipped_file} - {""}


def spectral_calibration(sfh_model, fits_directory, truths, skipped_jobs_file):
    """Calibration test of one model's spectral fitter on its stored mocks; prints the table and the verdict.

    Inputs: sfh_model: "delayed" or "dblplaw"; fits_directory: the model's spectral fits folder (str);
            truths: dict mock job_id -> Series of the true labels (from the model's mock file);
            skipped_jobs_file: the model's list of mocks deliberately left unfitted (str).
    Returns: (fitter_calibrated, coverage table from sm.calibration_coverage). fitter_calibrated is True when
             the 95% coverage of every calibration label of the model reaches sm.MIN_COVERAGE_95, False
             otherwise, and pd.NA while a mock is neither fitted nor skipped, or when there are no mocks.
    """
    # 1. Posteriors of the calibration mocks that have finished
    posteriors = {}
    n_converged = 0
    for job_id in truths:
        fit_file = os.path.join(fits_directory, f"{job_id}.h5")
        if os.path.exists(fit_file):
            posteriors[job_id], metadata = sm.load_fit(fit_file)
            assert sm.fit_is_current(metadata, sfh_model), (
                f"{fit_file} is from an earlier model or run: run section 4")
            n_converged += int(converged_from_multinest(job_id, metadata, fits_directory))

    # 2. Coverage per label and interval
    n_finished = len(posteriors)
    table = sm.calibration_coverage(posteriors, truths, sm.MODELS[sfh_model].calibration_labels)
    display(table.pivot(index="label", columns="level", values=["coverage", "n_jobs"]).round(2))

    # 3. Verdict: every mock must be fitted or deliberately skipped (listed in skipped_jobs_file)
    skipped_mocks = sorted(skipped_job_ids(skipped_jobs_file) & set(truths) - set(posteriors))
    coverage_95 = table[table["level"] == 0.95].set_index("label")["coverage"]
    if skipped_mocks:
        # Worst case: every skipped mock would have fallen outside its 95% interval
        worst_case_95 = coverage_95 * n_finished / (n_finished + len(skipped_mocks))
        print(f"Skipped mocks ({sfh_model}): {', '.join(skipped_mocks)}; 95% coverage if all of them had "
              f"missed: {worst_case_95.round(2).to_dict()}")
    if len(truths) == 0 or n_finished + len(skipped_mocks) < len(truths):
        print(f"INCOMPLETE ({sfh_model}): {n_finished} of {len(truths)} calibration mocks fitted "
              f"({n_converged} converged); the verdict stays empty until all are done")
        return pd.NA, table
    fitter_calibrated = bool((coverage_95 >= sm.MIN_COVERAGE_95).all())
    verdict = "PASS" if fitter_calibrated else "FAIL"
    print(f"{verdict} ({sfh_model}): 95% coverage {coverage_95.round(2).to_dict()} against the required "
          f"{sm.MIN_COVERAGE_95} ({n_finished} mocks, {n_converged} converged, {len(skipped_mocks)} skipped)")
    return fitter_calibrated, table


# The delayed-tau spectral fitter, on the mocks of section 4
fitter_calibrated, calibration_table = spectral_calibration("delayed", FITS_DIRECTORY, calibration_truths,
                                                            SKIPPED_JOBS_FILE)
'''

# ---------------------------------------------------------------------------------------------------------------
# Section 6d: the DPL spectral fits
SECTION_6D_MARKDOWN = r'''## 6d. Double-power-law spectral fits (DPL spec, amendment 2)
Every galaxy with a prepared spectrum is fitted again under the double-power-law prior. The spectrum, redshift, calibration polynomial, likelihood and MultiNest settings are those of section 4; only the SFH model changes, so the unit cube is `sm.spectral_dims("dblplaw")` (the DPL photometric dimensions plus `veldisp` and `noise_scale`). The jobs (`sm.dblplaw_spectral_jobs`: the galaxy jobs of section 4, same ids) go to their own folder, `spectral_fits_dblplaw/jobs/`. This section does not wait for the split: it fits every prepared spectrum.

**Calibration:** 33 mocks drawn from the DPL prior (`sm.draw_calibration_mocks(..., "dblplaw")`), each with the wavelengths, mask and errors of a random selected galaxy, as for the delayed-τ mocks. They are stored in `spectral_fits_dblplaw/calibration_mocks.joblib` and reused. The fitter passes when at least `sm.MIN_COVERAGE_95` of the truths lie inside the 95% intervals of `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity` and `log_mstar` (`sm.MODELS["dblplaw"].calibration_labels`); mocks listed in `skipped_fit_jobs_dblplaw.txt` are handled as in section 5.

**Results:** `spec_dblplaw_<label>_p16/_median/_p84` for the DPL training labels plus `veldisp` and `noise_scale` (mass-like labels at the aperture r brightness, as in section 6), `spec_dblplaw_converged`, `spec_dblplaw_log_evidence`, `spec_dblplaw_log_evidence_error` and `spec_dblplaw_fitter_calibrated`. With `spec_log_evidence` of section 6, they give the spectral Bayes factor of the two models in `ml_model.ipynb`. The delayed-τ spectral columns keep their plain names; the split, the photometric fits and the validation rulers of section 7 stay on the delayed-τ fits.

**Running:** `FITS_SUBDIR=spectral_fits_dblplaw remote/submit_fits.sh` (about 1400 fits, about 30 h: two 24 h jobs, the second resuming from the fetched fits), and `remote/fetch_fits.sh <job-id>` after each. From the first export on, every notebook-2 run must bundle `spectral_fits_dblplaw` and `skipped_fit_jobs_dblplaw.txt` (see `remote/README.md`).'''

SECTION_6D_JOBS = r'''DBLPLAW_SPECTRAL_FITS_DIRECTORY = sm.model_file(FITS_DIRECTORY, "dblplaw")    # spectral_fits_dblplaw
DBLPLAW_SPECTRAL_JOBS_DIRECTORY = os.path.join(DBLPLAW_SPECTRAL_FITS_DIRECTORY, sm.JOBS_SUBDIRECTORY)
DBLPLAW_SPECTRAL_MOCK_FILE = os.path.join(DBLPLAW_SPECTRAL_FITS_DIRECTORY, sm.CALIBRATION_MOCK_FILE)
DBLPLAW_SKIPPED_JOBS_FILE = sm.model_file(SKIPPED_JOBS_FILE, "dblplaw")      # skipped_fit_jobs_dblplaw.txt
DBLPLAW_SPECTRAL_MOCK_SEED = 3   # templates and noise of the DPL spectral mocks (SEED: delayed-tau)
os.makedirs(DBLPLAW_SPECTRAL_FITS_DIRECTORY, exist_ok=True)

# 1. DPL spectral calibration mocks: drawn once and stored with their true labels, as the delayed-tau mocks of
# section 4, because their fits are cached on disk and must be compared with the mocks they were made from
if os.path.exists(DBLPLAW_SPECTRAL_MOCK_FILE):
    dblplaw_spectral_mock_jobs, dblplaw_spectral_truths = joblib.load(DBLPLAW_SPECTRAL_MOCK_FILE)
    print(f"Reusing {len(dblplaw_spectral_mock_jobs)} stored DPL spectral calibration mocks from "
          f"{DBLPLAW_SPECTRAL_MOCK_FILE}")
else:
    # Why refuse: old calib_*.h5 fits would be scored against the new truths without any error. calib_* also
    # matches the fits' MultiNest folders
    stale_fits = glob.glob(os.path.join(DBLPLAW_SPECTRAL_FITS_DIRECTORY, f"{sm.CALIBRATION_ID_PREFIX}*"))
    assert len(stale_fits) == 0, (f"{len(stale_fits)} calib_* fits or MultiNest folders exist but "
                                  f"{DBLPLAW_SPECTRAL_MOCK_FILE} is missing: bundle it (with the fits), or "
                                  "move the fits to superseded/ to draw new mocks")
    templates = [(job["spectrum"], job["redshift"]) for job in real_jobs]
    dblplaw_spectral_mock_jobs, dblplaw_spectral_truths = sm.draw_calibration_mocks(
        sm.N_CALIBRATION_MOCKS, templates, np.random.default_rng(DBLPLAW_SPECTRAL_MOCK_SEED), "dblplaw")
    joblib.dump((dblplaw_spectral_mock_jobs, dblplaw_spectral_truths), DBLPLAW_SPECTRAL_MOCK_FILE)
    print(f"Drew {len(dblplaw_spectral_mock_jobs)} DPL spectral calibration mocks from the DPL prior, stored "
          f"in {DBLPLAW_SPECTRAL_MOCK_FILE}")

# 2. The galaxies of section 4, each fitted again with the DPL prior (same job id, spectrum and redshift)
dblplaw_spectral_galaxy_jobs = sm.dblplaw_spectral_jobs(real_jobs)

# 3. Job files, and how to run them
for job in [*dblplaw_spectral_mock_jobs, *dblplaw_spectral_galaxy_jobs]:
    sm.export_job(job, DBLPLAW_SPECTRAL_JOBS_DIRECTORY)
print(f"Exported {len(dblplaw_spectral_mock_jobs)} DPL calibration mocks and "
      f"{len(dblplaw_spectral_galaxy_jobs)} galaxies to {DBLPLAW_SPECTRAL_JOBS_DIRECTORY}/")
print("Run them on the PC:  FITS_SUBDIR=spectral_fits_dblplaw remote/submit_fits.sh   (about 30 h: two 24 h "
      "jobs, the second resuming from the fetched fits)")
print("and fetch each:      remote/fetch_fits.sh <job-id>")
'''

SECTION_6D_CALIBRATION = r'''# The DPL spectral fitter, on its own mocks: section 5's test with the DPL calibration labels
dblplaw_spectral_fitter_calibrated, dblplaw_spectral_calibration_table = spectral_calibration(
    "dblplaw", DBLPLAW_SPECTRAL_FITS_DIRECTORY, dblplaw_spectral_truths, DBLPLAW_SKIPPED_JOBS_FILE)
'''

SECTION_6D_READ = r'''# Summaries of every DPL training label, then the spectral fit's own parameters (as LABEL_COLUMNS, section 6)
DBLPLAW_SPECTRAL_LABELS = [*sm.MODELS["dblplaw"].training_labels, "veldisp", "noise_scale"]
DBLPLAW_SPECTRAL_PREFIX = sm.SPECTRAL_PREFIXES["dblplaw"]       # spec_dblplaw_


def dblplaw_spectral_summary(galaxies):
    """One row per galaxy with a DPL spectral fit: label summaries, convergence and evidence.

    The mass-like labels are moved to the aperture brightness first, as for the delayed-tau fits.
    Inputs: galaxies: the reference table (specobjid, and r_aper [mag]: the dereddened 3" r magnitude).
    Returns: DataFrame with specobjid, then spec_dblplaw_<label>_p16/_median/_p84 for DBLPLAW_SPECTRAL_LABELS,
             spec_dblplaw_converged (bool), spec_dblplaw_log_evidence and spec_dblplaw_log_evidence_error
             (dimensionless ln Z and MultiNest's uncertainty of it); galaxies without a fit have no row.
    """
    empty_summary = pd.DataFrame(np.zeros((1, len(DBLPLAW_SPECTRAL_LABELS))), columns=DBLPLAW_SPECTRAL_LABELS)
    columns = ["specobjid"]
    for name in sm.summarise_posterior(empty_summary, DBLPLAW_SPECTRAL_LABELS):
        columns.append(f"{DBLPLAW_SPECTRAL_PREFIX}{name}")
    for name in ("converged", "log_evidence", "log_evidence_error"):
        columns.append(f"{DBLPLAW_SPECTRAL_PREFIX}{name}")

    rows = []
    for galaxy in galaxies.itertuples():
        fit_file = os.path.join(DBLPLAW_SPECTRAL_FITS_DIRECTORY, f"{galaxy.specobjid}.h5")
        if not os.path.exists(fit_file):
            continue
        posterior, metadata = sm.load_fit(fit_file)
        assert metadata.get("kind", sm.SPECTRUM_JOB) == sm.SPECTRUM_JOB, f"{fit_file} is not a spectral fit"
        assert sm.fit_is_current(metadata, "dblplaw"), (
            f"{fit_file} is from an earlier model or run: run section 4")

        # Mass-like labels follow the aperture r magnitude, as for the delayed-tau fits (section 6)
        posterior = sm.rescale_mass_labels(posterior, galaxy.r_aper - posterior["model_r"])

        converged = converged_from_multinest(galaxy.specobjid, metadata, DBLPLAW_SPECTRAL_FITS_DIRECTORY)
        row = {"specobjid": galaxy.specobjid, f"{DBLPLAW_SPECTRAL_PREFIX}converged": converged,
               f"{DBLPLAW_SPECTRAL_PREFIX}log_evidence": metadata["log_evidence"],
               f"{DBLPLAW_SPECTRAL_PREFIX}log_evidence_error": metadata.get("log_evidence_error", np.nan)}
        for name, value in sm.summarise_posterior(posterior, DBLPLAW_SPECTRAL_LABELS).items():
            row[f"{DBLPLAW_SPECTRAL_PREFIX}{name}"] = value
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


# The DPL spectral columns replace any earlier ones (re-running the cell adds no copies); the verdict goes on
# every row
dblplaw_spectral = dblplaw_spectral_summary(reference)
earlier_columns = [column for column in reference.columns if column.startswith(DBLPLAW_SPECTRAL_PREFIX)]
reference = reference.drop(columns=earlier_columns).merge(dblplaw_spectral, on="specobjid", how="left")
reference[f"{DBLPLAW_SPECTRAL_PREFIX}fitter_calibrated"] = dblplaw_spectral_fitter_calibrated
n_dblplaw_spectral_converged = int(dblplaw_spectral[f"{DBLPLAW_SPECTRAL_PREFIX}converged"].sum())
both_converged = reference["converged"].eq(True) & reference[f"{DBLPLAW_SPECTRAL_PREFIX}converged"].eq(True)
print(f"DPL spectral fits read for {len(dblplaw_spectral)} of {len(dblplaw_spectral_galaxy_jobs)} exported "
      f"galaxies ({n_dblplaw_spectral_converged} converged; {int(both_converged.sum())} galaxies with both "
      f"spectral fits converged); spec_dblplaw_fitter_calibrated = {dblplaw_spectral_fitter_calibrated}")
'''

# ---------------------------------------------------------------------------------------------------------------
title_index = find_cell("# Spectroscopic reference set (part 2)")
section_4_index = find_cell("## 4. Fits (outside the notebook)")
section_5_index = find_cell("## 5. Calibration test")
calibration_index = find_cell("def converged_from_multinest(job_id, metadata):")
section_6_index = find_cell("## 6. Labels, coverage and split")
labels_index = find_cell("# Summaries of every training-set label, then the spectral fit's own parameters")
section_6c_index = find_cell("## 6c. Double-power-law photometric fits and cross-fits (DPL spec, part 2)")
cross_table_index = find_cell("# Every mock of either model under both models")
next_pass_index = find_cell("## Next pass")

# Title cell: the DPL spectral fits and their outputs, and the label list by its registry name
replace_in_cell(title_index, "(`sm.TRAINING_LABELS`)", "(`sm.MODELS[\"delayed\"].training_labels`)")
replace_in_cell(title_index, "between 0.1 and 1 Gyr and older.",
                "between 0.1 and 1 Gyr and older. Section 6d fits the same spectra again with Bagpipes' double "
                "power law (`dblplaw`, DPL spec, amendment 2), into the `spec_dblplaw_*` columns; the split and "
                "the validation stay on the delayed-τ fits.")
replace_in_cell(title_index,
                "(double-power-law photometric fits, section 6c), ",
                "(double-power-law photometric fits, section 6c), `spectral_fits_dblplaw/jobs/` and "
                "`spectral_fits_dblplaw/calibration_mocks.joblib` (double-power-law spectral fits, section 6d), ")

# Section 4: the superseding covers the DPL spectral folder (it loops over sm.FIT_FOLDER_MODELS)
replace_in_cell(section_4_index, "(`spectral_fits/`, `photometric_fits/`, `photometric_fits_dblplaw/`)",
                "(`spectral_fits/`, `photometric_fits/`, `photometric_fits_dblplaw/`, `spectral_fits_dblplaw/`)")

# Section 5: the label list by its registry name, and the DPL test of section 6d
replace_in_cell(section_5_index, "(`sm.CALIBRATION_LABELS`; age, τ and Z in log",
                "(`sm.MODELS[\"delayed\"].calibration_labels`; age, τ and Z in log")
replace_in_cell(section_5_index, "A failure is a finding about the prior or the model, not a crash.",
                "A failure is a finding about the prior or the model, not a crash. Section 6d runs the same test "
                "(`spectral_calibration`) on the DPL spectral fits, with their own mocks and "
                "`skipped_fit_jobs_dblplaw.txt`.")
set_source(cells[calibration_index], CALIBRATION_CELL)

# Section 6: the delayed-tau spectral evidence
replace_in_cell(section_6_index,
                "Each posterior sample carries the training set's labels (`sm.TRAINING_LABELS`) plus the fit's "
                "velocity dispersion and noise scale.",
                "Each posterior sample carries the delayed-τ training labels (`sm.MODELS[\"delayed\"]"
                ".training_labels`) plus the fit's velocity dispersion and noise scale. The fit's ln Z and "
                "MultiNest's error of it become `spec_log_evidence` and `spec_log_evidence_error` (the error is "
                "empty for the fits saved before it was stored, and is not recomputed).")
replace_in_cell(labels_index, '''    row["fit_hours"] = metadata["run_time"] / 3600
''', '''    row["fit_hours"] = metadata["run_time"] / 3600
    # ln Z of the fit and MultiNest's error of it, stored by the fits made since the DPL model; earlier fits
    # have no error, and it is not recomputed (the ln Z difference to the DPL fits is expected to dwarf it)
    row["spec_log_evidence"] = metadata["log_evidence"]
    row["spec_log_evidence_error"] = metadata.get("log_evidence_error", np.nan)
''')
replace_in_cell(labels_index,
                'expected_columns = ["converged", "fit_hours", *sm.summarise_posterior(empty_summary, LABEL_COLUMNS)]',
                'expected_columns = ["converged", "fit_hours", "spec_log_evidence", "spec_log_evidence_error",\n'
                '                    *sm.summarise_posterior(empty_summary, LABEL_COLUMNS)]')

# Section 6c: its mocks predate the DPL spectral fits
replace_in_cell(section_6c_index,
                "There are no DPL spectral mocks, so its **calibration mocks** are 33 points",
                "Its **calibration mocks** (drawn before the DPL spectral fits of section 6d existed, and kept) "
                "are 33 points")

# Next pass: the DPL spectral jobs follow
replace_in_cell(next_pass_index, "Add DR4 Stripe 82 galaxies",
                "Section 6d exports the DPL spectral jobs of the same prepared spectra, so a new pass adds their DPL "
                "jobs too, and its stored mocks are reused likewise. Add DR4 Stripe 82 galaxies")

# Section 6d, after the cross-fit table and before section 7
new_cells = [new_cell("markdown", SECTION_6D_MARKDOWN), new_cell("code", SECTION_6D_JOBS),
             new_cell("code", SECTION_6D_CALIBRATION), new_cell("code", SECTION_6D_READ)]
cells[cross_table_index + 1:cross_table_index + 1] = new_cells

with open(NOTEBOOK, "w") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print("notebook updated:", len(cells), "cells")
```

Expected output: `notebook updated: 40 cells`.

- [ ] **Step 2: Check the notebook (never execute it on the laptop)**

```bash
../.venv/bin/python -c "
import nbformat
nb = nbformat.read('spectroscopic_reference.ipynb', as_version=4); nbformat.validate(nb)
for i, c in enumerate(nb.cells):
    if c.cell_type == 'code': compile(c.source, f'cell {i}', 'exec')
print('ok', len(nb.cells))
print([c.source.splitlines()[0] for c in nb.cells if c.cell_type == 'markdown' and c.source.startswith('## ')])"
grep -c 'SPECTRAL_DIMS\|SPECTRAL_MODEL\|sm.CALIBRATION_LABELS\|sm.TRAINING_LABELS' spectroscopic_reference.ipynb
grep -c 'spectral_calibration(' spectroscopic_reference.ipynb
LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_notebooks.py -q
```

Expected:
- `ok 40`;
- the headings list has `## 6c. Double-power-law photometric fits and cross-fits (DPL spec, part 2)`, then `## 6d. Double-power-law spectral fits (DPL spec, amendment 2)`, then `## 7. Validation of the spectral fits`;
- the first `grep` prints `0`;
- the second prints `3` (the definition, the delayed-τ call and the DPL call: they are in three different cells, and `grep -c` counts lines);
- `3 passed` (every `sm.` name the notebooks use exists).

- [ ] **Step 3: `remote/README.md`: the notebook-2 runs carry the DPL spectral folder**

Save as `$SCRATCH/task14_readme.py` and run `../.venv/bin/python "$SCRATCH/task14_readme.py"` from `Synthetic_domain/`:

```python
"""Task 14, step 3: remote/README.md says how the notebook-2 runs carry the DPL spectral folder."""

README = "remote/README.md"
text = open(README).read()
# The new text is a raw string: its shell command ends lines with a backslash
pairs = [
    ('''`phot_*` and `phot_dblplaw_*` columns. Leave it True afterwards: with False, the next reading run leaves those
columns empty for the test half again.
''', r'''`phot_*` and `phot_dblplaw_*` columns. Leave it True afterwards: with False, the next reading run leaves those
columns empty for the test half again.

**DPL spectral runs (DPL spec, amendment 2).** Section 6d of notebook 2 exports `spectral_fits_dblplaw/jobs/`
and `spectral_fits_dblplaw/calibration_mocks.joblib` in every run; the first run that has it draws the mocks.
From then on, every notebook-2 run must bundle `spectral_fits_dblplaw` (the mock file and the fits) and
`skipped_fit_jobs_dblplaw.txt`. A run without them does not fail: section 6d draws throw-away mocks
(`fetch_notebook.sh` keeps the local ones) and leaves the `spec_dblplaw_*` columns empty. Between the export and
the reading run, fit `spectral_fits_dblplaw` (above) and fetch it. The reading run writes the `spec_dblplaw_*`
columns, and section 6 writes `spec_log_evidence` and `spec_log_evidence_error` of the delayed-tau spectral fits.
Before the first export, make the empty folder (`mkdir -p spectral_fits_dblplaw`): the bundle stops on a missing
input. From then on every notebook-2 run uses this command, with all four fit folders and both skip lists:

    cd Synthetic_domain
    EXCLUDE_NAMES="superseded pilot jobs" FLOAT32_FILES=training_set.parquet remote/submit_notebook.sh \
        spectroscopic_reference.ipynb dr3_corrected.csv ../Data/sdss_eboss_firefly-dr16_slim.fits \
        ../dustMaps/csfd training_set.parquet spectra spectral_fits spectral_fits_dblplaw photometric_fits \
        photometric_fits_dblplaw skipped_fit_jobs.txt skipped_fit_jobs_dblplaw.txt splus_aperture_photometry.csv \
        splus_noise_model.csv spectroscopic_reference.csv literature
    remote/fetch_notebook.sh ID

Notebook 2 reads no `jobs/` folder (it exports every job again), so `jobs` in `EXCLUDE_NAMES` keeps all the
`jobs/` folders out of the bundle (about 80 MB: 39 MB per spectral jobs folder); `fetch_notebook.sh` brings back
the jobs the run exported and still keeps the local `calib_*.npz`. The bundle is still large: 1.5 GB on
2026-10-29 with 788 delayed-tau spectral fits and no DPL spectral folder (`spectra/` 1.4 GB, the CSFD map
289 MB), and about 2 GB once both spectral folders are finished (about 300 MB of fits each). That is close to
hljobd's 3 GB limit, so check it first with `BUNDLE_ONLY=<folder>` in front of the same command.
'''),
]
for old, new in pairs:
    assert text.count(old) == 1, old[:80]
    text = text.replace(old, new)
open(README, "w").write(text)
print(README, "updated")
```

Expected output: `remote/README.md updated`.

- [ ] **Step 4: Run the whole suite**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -q`

Expected: `282 passed` (or `281 passed, 1 skipped`, as in Task 13), the same as after Task 13, because this task adds no test. On a checkout that lacks git-ignored data, more tests skip themselves (a review copy without the data gave `275 passed, 7 skipped`); compare passed + skipped, which is 282.

- [ ] **Step 5: Commit**

```bash
git add Synthetic_domain/spectroscopic_reference.ipynb Synthetic_domain/remote/README.md
git commit -m "feat: notebook 2 exports and reads the DPL spectral fits, with their calibration and evidence

New section 6d exports spectral_fits_dblplaw/jobs (every prepared spectrum and 33 DPL mocks),
runs the spectral calibration test for DPL (section 5's test, now a function of the model) and
reads the fits into spec_dblplaw_* columns, with the DPL ln Z and its error. Section 6 adds
spec_log_evidence(_error) of the delayed-tau spectral fits; the split stays on them.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

**PC order after Task 14 (for Task 16, the controller).** This is Amendment 2, item 6, step by step; Task 15's handoff continues it with the locked run of `ml_model.ipynb`. Every command runs from `Synthetic_domain/`. Rules for every notebook-2 run:
- supersede locally first (`remote/README.md`, "Notebook 2");
- use the command of `remote/README.md`, "DPL spectral runs" (all four fit folders, both skip lists, `EXCLUDE_NAMES="superseded pilot jobs"`), and check its size with `BUNDLE_ONLY=<folder>` first;
- fetch it with `remote/fetch_notebook.sh <id>`.

**(a) The remaining delayed-τ spectral fits.** Wait for the running delayed-τ job (bundle of 9b4676e) to finish, fetch it with `remote/fetch_fits.sh <id>`, then supersede locally. Its fits grow both halves of the split (no galaxy changes half). If it stopped at its time limit before every job was fitted, run `remote/submit_fits.sh` again and fetch that job too: it resumes from the fetched fits.

**(b) Unlock the test half's photometric jobs.** Set `EXPORT_LOCKED_TEST_JOBS = True` in notebook 2's parameter cell, keeping the JSON layout, and commit the change:

```bash
LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python - <<'PYTHON'
import json

with open("spectroscopic_reference.ipynb", encoding="utf-8") as notebook_file:
    notebook = json.load(notebook_file)
cell = next(cell for cell in notebook["cells"] if cell["id"] == "2f7e52ed")
old_line = "EXPORT_LOCKED_TEST_JOBS = False\n"
new_line = "EXPORT_LOCKED_TEST_JOBS = True\n"
assert cell["source"].count(old_line) == 1, "the EXPORT_LOCKED_TEST_JOBS line is not the expected one"
cell["source"][cell["source"].index(old_line)] = new_line
with open("spectroscopic_reference.ipynb", "w", encoding="utf-8") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print("EXPORT_LOCKED_TEST_JOBS = True")
PYTHON
git diff --stat
git add spectroscopic_reference.ipynb
git commit -m "run: notebook 2 exports and reads the test half's photometric jobs (EXPORT_LOCKED_TEST_JOBS = True)

Set once for the locked test of both SFH models (DPL spec, amendment 2, item 5); it stays True.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Expected: `EXPORT_LOCKED_TEST_JOBS = True`, and `git diff --stat` shows one file with 1 insertion and 1 deletion.

**(c) The export run.** First `mkdir -p spectral_fits_dblplaw` (the bundle stops on a missing input), then the notebook-2 run and its fetch. With the grown split and `EXPORT_LOCKED_TEST_JOBS = True`, this one run exports:
- the DPL spectral jobs of every prepared spectrum, and draws and stores the 33 DPL spectral mocks (section 6d);
- the photometric jobs of both models for the validation half (its new galaxies are the ones not yet fitted) and the whole test half (sections 6b and 6c);
- the cross-fits (section 6c; they are fitted already, so the fit runs skip them).

Check in the fetched notebook: "Exported 33 DPL calibration mocks and <n> galaxies to spectral_fits_dblplaw/jobs/" (about 1360 galaxies), and "(validation and test halves)" in the export lines of `photometric_fits/jobs/` and `photometric_fits_dblplaw/jobs/`. A line with "(validation halves)" means the flag of (b) did not reach the run: stop there.

**(d) The fits, each run followed by `remote/fetch_fits.sh <id>`:**
1. `FITS_SUBDIR=spectral_fits_dblplaw remote/submit_fits.sh` (about 1400 fits, about 30 h), fetch, then the same command again for the rest, fetch. The second job resumes from the fetched fits.
2. `FITS_SUBDIR=photometric_fits_dblplaw remote/submit_fits.sh`, fetch.
3. `FITS_SUBDIR=photometric_fits remote/submit_fits.sh`, fetch.

Any run that stops at its 24 h limit is submitted again the same way after its fetch.

**(e) The reading run:** the same notebook-2 command, after the fits of (d) are fetched. Check in the fetched notebook:
- the DPL spectral verdict is `PASS (dblplaw)` or `FAIL (dblplaw)`, not `INCOMPLETE (dblplaw)`;
- "DPL spectral fits read for <n> of <n> exported galaxies" has nearly equal numbers;
- the photometric readings of both models cover the validation and the test half.

Then go on with the pre-flight check of Task 15's handoff, and the locked run of `ml_model.ipynb`.

**Optional overlap.** The DPL spectral jobs do not depend on the split, so an extra export run during (a) could start step (d) 1 before the delayed-τ fits are done. That run keeps `EXPORT_LOCKED_TEST_JOBS = False` and bundles the 788 fetched delayed-τ fits; its photometric exports and its CSV are replaced by (c) and (e). It saves time only if the job server runs two 16-core jobs side by side; otherwise the two fit jobs queue anyway. The user asked for the delayed-τ fits first, so the default is the order above.

### Task 15: `ml_model.ipynb` compares both models with both spectral references, adds the spectral evidence, and runs the locked test for both models

**Files:**
- Modify: `Synthetic_domain/ml_evaluation.py`:
  - a new constant `FAVOURED_MODEL_ORDER` (section 1);
  - two new functions `log_bayes_sign_agreement` and `favoured_model_crosstab` (section 6, before `true_model_recovery`);
  - two docstrings name the spectral prefixes.
- Test: `Synthetic_domain/tests/test_ml_evaluation.py` (3 tests appended; 31 → 34 tests in the file).
- Modify: `Synthetic_domain/ml_model.ipynb`:
  - 13 cells are rewritten in place, by cell id: `55011031`, `fe480605`, `e81a48b9`, `e05c91f8`, `74230b0e`, `2e021085`, `3e17feaa`, `2887477a`, `b051ebb0`, `89339c2d`, `79e11a83`, `10a86e1d`, `0267e995`;
  - 2 new cells (section 5.5: markdown `a7c3e5d1`, code `b8d4f6e2`) are inserted after the cross-fit cell `fc196542`, before the last cell `3dfd5bb5` (which writes `COMPARISON_FILE`);
  - 25 → 27 cells. Cells not listed keep their source and outputs; the edited code cells lose their outputs (they are stale).
- Temporary, deleted in Step 10: `Synthetic_domain/task15_cells/` (the new cell sources and the script that puts them into the notebook).
- Unchanged: every cell not listed, among them cell 1 (`b6229f82`, the parameters), sections 1, 2 and 4, the photometric evidence (`b1c2a26c`), the cross-fits (`fc196542`) and the last cell (`3dfd5bb5`).

**Do not execute the notebook on the laptop.** It runs on the PC as a GPU notebook job (handoff at the end of this task). Step 8 runs the changed cells' table code with stubbed posteriors: no notebook is executed with nbconvert, and there is no Bagpipes or flow work (no sampling, no fit, about 30 s on the CPU). **Controller approval needed:** Step 8 `exec()`s code cells of `ml_model.ipynb` on the laptop, which the rule "never execute notebooks on the laptop" could be read to forbid. It is the only check that the delayed-τ tables stay unchanged, so this plan keeps it, with the controller's explicit approval. Without that approval, skip Step 8: Step 7 (nbformat, compile, flake8 F821, `tests/test_notebooks.py`) still runs, and the equality of the delayed-τ tables is then checked on the PC output against the previous run's CSVs instead.

**Interfaces:**
- **Consumes:**
  - `sm.SPECTRAL_PREFIXES = {"delayed": "", "dblplaw": "spec_dblplaw_"}` (Task 13), next to `sm.PHOTOMETRIC_PREFIXES`. With prefix `""` the delayed-τ spectral columns keep their plain names: `<label>_p16/_median/_p84`, `converged`, `fitter_calibrated`;
  - **the reference-CSV columns of Task 14** (`spectroscopic_reference.csv`):
    - `spec_dblplaw_<label>_p16/_median/_p84` for `sfh_models.MODELS["dblplaw"].training_labels`;
    - `spec_dblplaw_converged`, `spec_dblplaw_fitter_calibrated`;
    - `spec_log_evidence`, `spec_log_evidence_error` (delayed-τ spectral fits; the error is NaN for the fits made before it was stored) and `spec_dblplaw_log_evidence`, `spec_dblplaw_log_evidence_error`;
    - the `split` column stays the delayed-τ spectral split: every galaxy of a half has a converged delayed-τ spectral fit;
  - the existing `ev.*` helpers (`reference_comparison_rows`, `agreement_verdict`, `side_by_side_verdicts`, `summary_comparison_rows`, `model_comparison_rows`, `log_bayes_factor`, `log_bayes_factor_error`, `evidence_category_fractions`, `favoured_model`) and the notebook's own helpers (`available_labels`, `converged_fits`, `spectral_galaxy_types`, `labelled_posteriors`).
- **Produces:**
  - `ev.FAVOURED_MODEL_ORDER = [ev.DELAYED_MODEL, ev.NO_MODEL, ev.DBLPLAW_MODEL]`;
  - `ev.log_bayes_sign_agreement(first_log_bayes, second_log_bayes) -> dict(n_galaxies, same_sign, same_favoured_model)`;
  - `ev.favoured_model_crosstab(first_log_bayes, second_log_bayes) -> DataFrame` (3 × 3 counts, index name `first`, columns name `second`, both in `FAVOURED_MODEL_ORDER`);
  - notebook globals:
    - `SPECTRAL_REFERENCES = list(sm.SPECTRAL_PREFIXES)`;
    - `spectral_fits(galaxies, spectral_model) -> bool array`: fit_status `"fitted"` and the model's spectral `converged` flag;
    - `agreements[(ml_model, spectral_model)]`: ML rows against one spectral reference, with a new `spectral_model` column (it was `agreements[sfh_model]`, against the delayed-τ spectra only);
    - `verdicts[(sfh_model, reference_set, half)]`: the works-test verdicts (it was `verdicts[(sfh_model, reference_set)]`, validation only). `reference_set` "spectral" is now the model's **own** spectral fits;
    - `spectral_evidence`: one row per galaxy with both spectral fits, `specobjid`, `half`, `galaxy_type`, `spectral_log_bayes`, `spectral_log_bayes_error`;
  - outputs (written when the notebook runs on the PC):
    - `ml_evaluation_summary[_dblplaw].csv`: the agreement rows gain a `half` column (test-half rows with `RUN_LOCKED_TEST`);
    - `ml_model_comparison.csv`:
      - `own_photometric` rows gain `half`; new section `own_spectral`;
      - `model_pairs` gains the pair `spectral Bagpipes`;
      - `against_spectra` gains `spectral_model` (rows against both references);
      - new sections `spectral_evidence_distribution`, `spectral_evidence_categories`, `evidence_agreement`, `evidence_favoured_models`.

**What stays identical.** For the same ML draws and the same reference CSV, every delayed-τ table the notebook printed before is unchanged: the photometric verdicts, the delayed-τ ML against the delayed-τ spectra, the DPL ML against the delayed-τ spectra, the model pairs, `against_spectra` against the delayed-τ reference, the photometric evidence and the cross-fit recovery. Step 8 checks this. The PC numbers still change, because the remaining delayed-τ spectral fits grow both halves (Amendment 2, item 1).

**Design choices (why):**
- **No extra sampling.** The ML draws of cell 8 (one pass per model and half) serve every reference, so the four (ML model, spectral model) comparisons cost only table work, and the memory rule (one model's synthetic posteriors at a time, GPU job cap 7.1 GB) is untouched.
- **One galaxy type for both references.** The passive/star-forming split stays the delayed-τ spectral sSFR (as the validation rulers of notebook 2), so a galaxy is in the same group against both references.
- **The `all` row without judged galaxies is empty.** Before notebook 2 has read the DPL spectral fits, the DPL own-spectral verdict has no judged label; `ev.agreement_verdict`'s `all` row would then say `True` (an `all()` over nothing). `verdict_table` sets it to `None` ("no verdict") whenever the judged labels hold no galaxy. Delayed-τ verdicts always have judged galaxies, so they do not change.
- **Locked test.** With `RUN_LOCKED_TEST = True`, section 3b judges the `test` half too: against the model's own spectral fits always, against its own photometric fits only when notebook 2 has read the test-half photometric fits (`EXPORT_LOCKED_TEST_JOBS`); otherwise that verdict is skipped with a message. The spectral parts of section 5 (5.1 spectral, the 5.2 spectral pair, 5.3, 5.5) use `halves`; the photometric parts keep `photometric_halves`.
- **The spectral pair in 5.2** (delayed-τ spectral against DPL spectral fits, shared labels) is not listed in Amendment 2. It costs one `ev.model_comparison_rows` call and shows how far apart the two spectral references are, which the 5.3 caveat needs.

- [ ] **Step 1: Write the failing tests.** Append to `Synthetic_domain/tests/test_ml_evaluation.py`:

```python


def test_log_bayes_sign_agreement_of_a_hand_made_table():
    photometric = np.array([2.0, -5.0, 0.5, 4.0, np.nan, -0.2])
    spectral = np.array([10.0, -50.0, -30.0, -2.0, 1.0, 0.3])
    agreement = ev.log_bayes_sign_agreement(photometric, spectral)
    assert agreement["n_galaxies"] == 5                                     # the NaN galaxy is left out
    assert np.isclose(agreement["same_sign"], 2 / 5)                        # galaxies 0 and 1
    assert np.isclose(agreement["same_favoured_model"], 3 / 5)              # 0, 1, and 5 (both inconclusive)


def test_favoured_model_crosstab_of_a_hand_made_table():
    photometric = np.array([2.0, -5.0, 0.5, 4.0, np.nan, -0.2])
    spectral = np.array([10.0, -50.0, -30.0, -2.0, 1.0, 0.3])
    crosstab = ev.favoured_model_crosstab(photometric, spectral)
    assert list(crosstab.index) == ["delayed", "neither", "dblplaw"]
    assert list(crosstab.columns) == ["delayed", "neither", "dblplaw"]
    assert (crosstab.index.name, crosstab.columns.name) == ("first", "second")
    expected = np.array([[1, 0, 1],      # first favours delayed: second delayed (0), second DPL (3)
                         [0, 1, 1],      # first inconclusive: second inconclusive (5), second DPL (2)
                         [0, 0, 1]])     # first favours DPL: second DPL (1)
    assert np.array_equal(crosstab.to_numpy(), expected)


def test_evidence_agreement_without_galaxies_is_nan_and_zero():
    agreement = ev.log_bayes_sign_agreement([np.nan, 1.0], [2.0, np.nan])
    assert agreement["n_galaxies"] == 0
    assert np.isnan(agreement["same_sign"]) and np.isnan(agreement["same_favoured_model"])
    assert (ev.favoured_model_crosstab([np.nan], [np.nan]).to_numpy() == 0).all()
```

- [ ] **Step 2: Run them to check they fail.**
  - Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_evaluation.py -q`
  - Expected: `3 failed, 31 passed`, each failure an `AttributeError: module 'ml_evaluation' has no attribute 'log_bayes_sign_agreement'` (or `'favoured_model_crosstab'`).

- [ ] **Step 3: Add the constant and the two functions to `ml_evaluation.py`.**
  1. In section 1 (Constants), after the line `ALL_MOCKS_ROW = "all"                        # true_model_recovery row over every mock`, add:

```python
FAVOURED_MODEL_ORDER = [DELAYED_MODEL, NO_MODEL, DBLPLAW_MODEL]   # order of favoured_model_crosstab
```

  2. In section 6 (Bayesian evidence), between the end of `evidence_category_fractions` (its last line is `    return pd.DataFrame(rows)`) and `def true_model_recovery(cross_fits):`, insert (two blank lines before and after, as between the other functions):

```python
def log_bayes_sign_agreement(first_log_bayes, second_log_bayes):
    """How often two ln B of the same galaxies point the same way (e.g. photometric against spectral).

    Inputs: first_log_bayes, second_log_bayes: arrays (n,) of ln B [nats], paired by galaxy (positive favours
            delayed-tau). A galaxy with either ln B NaN is left out.
    Returns: dict with
             - n_galaxies: galaxies with both ln B finite (count);
             - same_sign: share of them whose two ln B have the same sign (a ln B of exactly 0 matches only
               another 0);
             - same_favoured_model: share of them whose two ln B favour the same model (favoured_model; two
               inconclusive ln B both favour NO_MODEL, so they count as the same).
             The shares are NaN without galaxies.
    """
    first = np.asarray(first_log_bayes, dtype=float)
    second = np.asarray(second_log_bayes, dtype=float)
    usable = np.isfinite(first) & np.isfinite(second)
    n_galaxies = int(usable.sum())
    if n_galaxies == 0:
        return {"n_galaxies": 0, "same_sign": np.nan, "same_favoured_model": np.nan}
    same_sign = np.sign(first[usable]) == np.sign(second[usable])
    same_favoured_model = favoured_model(first[usable]) == favoured_model(second[usable])
    return {"n_galaxies": n_galaxies, "same_sign": float(np.mean(same_sign)),
            "same_favoured_model": float(np.mean(same_favoured_model))}


def favoured_model_crosstab(first_log_bayes, second_log_bayes):
    """Counts of galaxies by the model each of two ln B favours (e.g. photometric against spectral).

    Inputs: first_log_bayes, second_log_bayes: arrays (n,) of ln B [nats], paired by galaxy. A galaxy with
            either ln B NaN is left out.
    Returns: DataFrame of counts, 3 x 3: rows (index name "first") are the model the first ln B favours,
             columns (name "second") the model the second favours, both in FAVOURED_MODEL_ORDER
             (favoured_model: DELAYED_MODEL, NO_MODEL for an inconclusive ln B, DBLPLAW_MODEL). An empty
             cell is 0.
    """
    first = np.asarray(first_log_bayes, dtype=float)
    second = np.asarray(second_log_bayes, dtype=float)
    usable = np.isfinite(first) & np.isfinite(second)
    first_models = favoured_model(first[usable])
    second_models = favoured_model(second[usable])
    counts = pd.DataFrame(0, index=pd.Index(FAVOURED_MODEL_ORDER, name="first"),
                          columns=pd.Index(FAVOURED_MODEL_ORDER, name="second"))
    for first_model, second_model in zip(first_models, second_models):
        counts.loc[first_model, second_model] += 1
    return counts
```

  3. Two docstrings name the spectral prefixes. In `reference_comparison_rows`, replace

```text
            labels: the labels to compare; prefix: "" (spectral fits) or "phot_" (photometric fits).
```

  with

```text
            labels: the labels to compare; prefix: a spectral prefix (sm.SPECTRAL_PREFIXES: "" for the
            delayed-tau spectral fits, "spec_dblplaw_" for the DPL ones) or a photometric one
            (sm.PHOTOMETRIC_PREFIXES: "phot_", "phot_dblplaw_").
```

  In `summary_comparison_rows`, replace

```text
            columns and the reference's {reference_prefix}<label>_* columns (by default the spectral fit's);
```

  with

```text
            columns and the reference's {reference_prefix}<label>_* columns (by default the delayed-tau
            spectral fit's; "spec_dblplaw_" for the DPL spectral fit's);
```

- [ ] **Step 4: Run the tests; commit.**
  - Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_evaluation.py -q`
  - Expected: `34 passed`.

```bash
git add Synthetic_domain/ml_evaluation.py Synthetic_domain/tests/test_ml_evaluation.py
git commit -m "feat: ml_evaluation compares two ln B per galaxy (sign agreement, favoured-model counts)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 5: Write the new cell sources.** Create the folder `Synthetic_domain/task15_cells/` and write each block below, verbatim, to the file named in its heading (`.md` for markdown cells, `.py` for code cells). Each file is the complete new source of that cell. The files end with one newline; the script of Step 6 strips it.

  **`task15_cells/55011031.md`** (cell 0, markdown, the notebook header):

```markdown
# SFH posterior from S-PLUS colours (part 3)

Two SFH models, each with its own normalizing flow (trained by `train_ml_model.py --sfh-model <model>` on the 3090):
- **delayed-τ** (`delayed`): 9 targets, which are log age, log τ, log metallicity, the dust (`av_ism`, `log_tau_bc`), `logU`, the formed and the surviving mass, and the redshift;
- **double power law** (`dblplaw`): 10 targets, which are log turnover lookback time (the lookback time to τ, the turnover of the double power law), log α and log β, then the delayed-τ targets from log metallicity on.

Every other label (`mw_age`, `f_young`, `f_mid`, `f_old`, `log_sfr`, `log_ssfr`, ...) is derived from the samples with the same functions as the training labels (`ml_data.decode_samples`). Inputs are the 11 colours against r, the r magnitude (only for the absolute masses), and z_obs with its uncertainty (0 for spectroscopic redshifts).

Sections 1–4 run once per model in `SFH_MODELS`. Section 5 compares the two models and needs both.

**Bagpipes references.** Each model has its own photometric fits (columns `sm.PHOTOMETRIC_PREFIXES[model]`) and its own spectral fits (columns `sm.SPECTRAL_PREFIXES[model]`: plain names for delayed-τ, `spec_dblplaw_` for DPL) in `spectroscopic_reference.csv`. The split into `validation` and `test` halves is the delayed-τ spectral split, for both models.

**Inputs:** delayed-τ keeps the plain names; DPL adds `_dblplaw` (`sfh_models.model_file`).
- Per model: `ml_runs/mlp[_dblplaw]/`, `synthetic_validation[_dblplaw].parquet`, `synthetic_test[_dblplaw].parquet`, `training_set[_dblplaw].parquet` (prior draws).
- Shared: `spectroscopic_reference.csv`, `photometric_cross_fits.csv` (ln Z of the cross-fitted calibration mocks, notebook 2 section 6c; section 5.4), `dr3_corrected.csv`.

**Outputs:** `ml_evaluation_summary[_dblplaw].csv`, `ml_sfh_catalogue[_dblplaw].parquet`, `ml_model_comparison.csv`.

**Run:** top to bottom with the *Inference4ETGs (MultiNest)* kernel, on the PC, as a GPU notebook job with 16 cores and a 12 h limit (`GPU=1 CORES=16 MEMORY_MB=7100 TIME_LIMIT=12:00:00 HLJOB_OPTIONS=--pause-llm remote/submit_notebook.sh ml_model.ipynb <inputs>`, with the inputs listed above). Decoding a DPL sample sums Bagpipes' age grid (about 0.06 ms), so the decoding runs in `N_DECODE_JOBS` threads.

**Locked test.** Sections marked **LOCKED** run only with `RUN_LOCKED_TEST = True`, once, for both models, after the user has inspected the validation results. They use the synthetic test sets and the reference set's test half:
- the spectral comparisons (sections 3, 3b, 5.1, 5.3, 5.5) cover the test half whenever it is unlocked, since both models' spectral fits cover every galaxy with a spectrum;
- the photometric Bagpipes fits of the test half exist only when notebook 2 exported and read them (`EXPORT_LOCKED_TEST_JOBS = True`); without them the photometric and photometric-evidence comparisons cover the validation half, and the notebook says so.

Posterior samples for the same seed differ between CUDA and CPU (`DEVICE`): numbers are reproducible per device, not across devices. `ml_networks.load_run` refuses a checkpoint trained for another SFH model or other targets (e.g. the 16-bin runs): retrain it.
```

  **`task15_cells/fe480605.md`** (section 3 header, markdown):

```markdown
## 3. Agreement with the Bagpipes spectral posteriors

For each SFH model: ML posteriors from the 3″ aperture colours and spectroscopic redshift, against each reference galaxy's spectral posterior. There are two spectral references, one per SFH model (`sm.SPECTRAL_PREFIXES`):
- the **delayed-τ spectral fits** (plain column names) hold the shared labels and the delayed-τ parameters;
- the **DPL spectral fits** (`spec_dblplaw_<label>`) hold the shared labels and the DPL parameters.

Every ML model is compared with both references, on the labels each reference holds. A galaxy enters a reference's comparison when its spectral fit of that model converged. The galaxy type (passive or star-forming) comes from the delayed-τ spectral sSFR for both references, so a galaxy has the same type in every table. Development uses the `validation` half; the `test` half is LOCKED.
```

  **`task15_cells/e81a48b9.py`** (section 3, code: the ML draws and the ML against both spectral references). Changes against the current cell: `SPECTRAL_REFERENCES`, the new helper `spectral_fits`, `agreements` keyed by (ML model, spectral model) with a `spectral_model` column, the fitter flag read from the reference's own `fitter_calibrated` column; the sampling is unchanged.

```python
reference = pd.read_csv(REFERENCE_FILE, dtype={"specobjid": str})
halves = ["validation", "test"] if RUN_LOCKED_TEST else ["validation"]
SPECTRAL_REFERENCES = list(sm.SPECTRAL_PREFIXES)   # one spectral Bagpipes reference per SFH model
aperture_columns = []
aperture_error_columns = []
for band in sm.BANDS:
    aperture_columns.append(f"{band}_aper")
    aperture_error_columns.append(f"e_{band}_aper")


def available_labels(table, labels, prefix=""):
    """The labels whose {prefix}<label>_median column the table holds (a reference set may lack some)."""
    kept = []
    for label in labels:
        if f"{prefix}{label}_median" in table.columns:
            kept.append(label)
    return kept


def spectral_galaxy_types(galaxies):
    """'star-forming' where the delayed-tau spectral fit's log_ssfr_median exceeds PASSIVE_SSFR_LIMIT.

    The other galaxies are 'passive'. Both spectral references use this one split, so a galaxy has the same
    type in every comparison.
    Returns: array of str, one per galaxy.
    """
    is_star_forming = galaxies["log_ssfr_median"].to_numpy(dtype=float) > PASSIVE_SSFR_LIMIT
    return np.where(is_star_forming, "star-forming", "passive")


def converged_fits(galaxies, prefix):
    """True where the galaxy's Bagpipes fit with this column prefix exists and converged (array).

    prefix: a photometric prefix (sm.PHOTOMETRIC_PREFIXES) or a spectral one (sm.SPECTRAL_PREFIXES; "" reads
    the delayed-tau spectral column "converged"). A table without the {prefix}converged column (e.g. before
    the DPL fits are read) gives False everywhere.
    """
    column = f"{prefix}converged"
    if column not in galaxies.columns:
        return np.zeros(len(galaxies), dtype=bool)
    return galaxies[column].astype("boolean").fillna(False).to_numpy(dtype=bool)


def spectral_fits(galaxies, spectral_model):
    """True where the galaxy has a converged spectral fit of this SFH model (array).

    The galaxy needs its delayed-tau spectral fit (fit_status "fitted"), which defines the split, and the
    converged flag of the model's own spectral fit. Every galaxy of a half has a converged delayed-tau fit
    (notebook 2 splits only those), so for delayed-tau this is every galaxy of the half.
    """
    has_delayed_fit = (galaxies["fit_status"] == "fitted").to_numpy()
    return has_delayed_fit & converged_fits(galaxies, sm.SPECTRAL_PREFIXES[spectral_model])


# Why one sampling pass per model and half: the spectral, photometric and cross-model comparisons must use
# the same ML draws of each galaxy, so each model's posterior is drawn once for every galaxy with a
# spectral fit or a photometric fit of that model
has_spectral_fit = reference["fit_status"] == "fitted"
# (sfh_model, half) -> (galaxies with an ML posterior, dict label -> (n, samples))
ml_posteriors_by_half = {}
# (ML model, spectral model) -> the ML of one model against the spectral fits of one model, one row per
# (half, galaxy, label)
agreements = {}
for sfh_model in SFH_MODELS:
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    has_photometric_fit = pd.Series(False, index=reference.index)
    if f"{prefix}mw_age_median" in reference.columns:
        has_photometric_fit = reference[f"{prefix}mw_age_median"].notna()
    n_targets = runs[sfh_model][1]["n_targets"]
    agreement_rows = {}
    for spectral_model in SPECTRAL_REFERENCES:
        agreement_rows[spectral_model] = []
    for half in halves:
        galaxies = reference[(reference["split"] == half) & (has_spectral_fit | has_photometric_fit)]
        galaxies = galaxies.reset_index(drop=True)
        if len(galaxies) == 0:
            print(f"{sfh_model}, {half}: no fitted galaxies yet")
            continue
        observations = ml_data.colour_inputs(galaxies[aperture_columns].to_numpy(),
                                             galaxies[aperture_error_columns].to_numpy(),
                                             galaxies["z"].to_numpy(), np.zeros(len(galaxies)))
        if not observations["valid"].any():
            print(f"{sfh_model}, {half}: no galaxy with a detected r band in the aperture photometry")
            continue
        observations["targets"] = np.zeros((len(galaxies), n_targets))   # unknown; not used
        posteriors, _ = labelled_posteriors(runs[sfh_model], observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED)
        kept = galaxies[observations["valid"]].reset_index(drop=True)
        ml_posteriors_by_half[(sfh_model, half)] = (kept, posteriors)

        # The spectral comparisons: for each spectral reference, the galaxies whose fit of that model
        # converged, on the labels that reference holds
        for spectral_model in SPECTRAL_REFERENCES:
            spectral_prefix = sm.SPECTRAL_PREFIXES[spectral_model]
            spectral_labels = available_labels(reference, LABELS[sfh_model]["reference_labels"],
                                               spectral_prefix)
            spectral = spectral_fits(kept, spectral_model)
            if not spectral.any() or not spectral_labels:
                continue
            spectral_galaxies = kept[spectral].reset_index(drop=True)
            spectral_posteriors = {label: posteriors[label][spectral] for label in spectral_labels}
            galaxy_types = dict(zip(spectral_galaxies["specobjid"], spectral_galaxy_types(spectral_galaxies)))
            galaxy_details = spectral_galaxies.set_index("specobjid")
            fitter_column = f"{spectral_prefix}fitter_calibrated"
            for row in ev.reference_comparison_rows(spectral_posteriors, spectral_galaxies, spectral_labels,
                                                    prefix=spectral_prefix):
                details = galaxy_details.loc[row["specobjid"]]
                agreement_rows[spectral_model].append({
                    "sfh_model": sfh_model, "spectral_model": spectral_model, "half": half, **row,
                    "galaxy_type": galaxy_types[row["specobjid"]], "fibre": details["fiber_diameter"],
                    "calibrated_fitter": details.get(fitter_column, np.nan)})
    for spectral_model in SPECTRAL_REFERENCES:
        agreement = pd.DataFrame(agreement_rows[spectral_model])
        agreements[(sfh_model, spectral_model)] = agreement
        if len(agreement) == 0:
            print(f"{sfh_model} ML: no converged {spectral_model} spectral fits to compare with in "
                  f"{REFERENCE_FILE} yet")
            continue
        print(f"{sfh_model} ML against the {spectral_model} spectral Bagpipes posteriors")
        display(agreement.groupby(["half", "label"]).agg(median_offset=("offset_in_widths", "median"),
                                                          inside_68=("inside_68", "mean"),
                                                          inside_95=("inside_95", "mean"),
                                                          galaxies=("specobjid", "count")).round(3))
        display(agreement.groupby(["half", "label", "galaxy_type"])["inside_95"].mean().unstack().round(3))
```

  **`task15_cells/e05c91f8.md`** (section 3b header, markdown):

```markdown
### 3b. The "works" test: the model's own photometric and spectral Bagpipes

For each SFH model, the ML is compared with two Bagpipes references of the same model, and each answers a different question:

- **Photometric Bagpipes (emulator check).** These are the model's own photometric fits (notebook 2, section 6b for delayed-τ and 6c for DPL): `phot_<label>` for delayed-τ, `phot_dblplaw_<label>` for DPL (`sm.PHOTOMETRIC_PREFIXES`). They fit the same 3″ photometry the ML sees, under that model's training prior and the simulator's likelihood, so a faithful emulator must agree with them.
- **Spectral Bagpipes (physical check).** These are the model's own spectral fits to the SDSS spectra: plain `<label>` columns for delayed-τ, `spec_dblplaw_<label>` for DPL (`sm.SPECTRAL_PREFIXES`). Agreement tells whether the colours recover what the spectra say under the same SFH model. Section 5.3 compares every estimator with both spectral references on the shared labels.

A label agrees with a reference set when two conditions hold:
- at least `ev.MIN_INSIDE_95` (80%) of the reference medians lie inside the ML 95% interval;
- the median offset is within ±`ev.MAX_ABS_MEDIAN_OFFSET` combined 1σ widths.

The judged labels are the model's calibration labels. For delayed-τ they are `log_age`, `log_tau`, `log_metallicity` and `log_mstar`. For DPL they are `log_turnover_lookback`, `log_alpha`, `log_beta`, `log_metallicity` and `log_mstar`. Both of a model's references hold all of them. A judged label a reference set lacks (e.g. before notebook 2 has read the DPL spectral fits) is listed with the status "not in the reference" and stays out of the `all` row; without any judged label holding galaxies, the `all` row gives no verdict (empty).

The other shared labels (`sfh_models.SHARED_LABELS`) are reported, not judged. A reported label a reference set does not hold has an empty row. The `validation` half is always judged. With `RUN_LOCKED_TEST`, the `test` half is judged too (LOCKED): against the spectral fits always, against the photometric fits only when notebook 2 has read the test-half photometric fits. The verdicts are printed side by side, and the notebook draws no conclusion beyond them.

From Amendment 2 on, `agreement_spectral` in `ml_evaluation_summary_dblplaw.csv` is the DPL ML against the DPL spectral fits; the DPL ML against the delayed-τ spectra is in `ml_model_comparison.csv`, section `against_spectra`, `spectral_model` = `delayed`.
```

  **`task15_cells/74230b0e.py`** (section 3b, code: the works test against the model's own photometric and own spectral fits, per half). Changes: the loop over `halves`, `verdicts[(sfh_model, reference_set, half)]`, the own spectral reference (`agreements[(sfh_model, sfh_model)]`, `sm.SPECTRAL_PREFIXES`), the empty `all` row in `verdict_table`, and two helpers (`photometric_comparison_rows`, `fitter_calibration_flag`) that hold the old inline code.

```python
REFERENCE_SET_NAMES = {"photometric": "its own photometric Bagpipes (emulator check)",
                       "spectral": "its own spectral Bagpipes (physical check)"}


NOT_IN_REFERENCE = "not in the reference"   # status of a judged label the reference set does not hold
# What gives a reference set fits in the test half (section 3b judges the test half only with RUN_LOCKED_TEST)
TEST_HALF_HINTS = {
    "photometric": "set EXPORT_LOCKED_TEST_JOBS = True in notebook 2, fit those jobs and read them",
    "spectral": "run notebook 2 after the spectral fits of the test half"}


def verdict_table(rows, judged_labels, reported_labels, not_in_reference=()):
    """One verdict table with a status column.

    Rows: agreement_verdict over the judged labels (status "judged", with its 'all' row, which covers only
    them); then the reported labels (status "reported", without their 'all' row); then the labels the
    reference set cannot judge (status NOT_IN_REFERENCE, no galaxies, agrees left empty).
    Why the empty 'all' row: when no judged label holds a galaxy (no judged label in the reference, or no
    fit read yet), agreement_verdict's 'all' row is no verdict, so its agrees is None rather than True or
    False.
    """
    judged = ev.agreement_verdict(rows, judged_labels).assign(status="judged")
    judged_label_rows = judged[judged["label"] != ev.ALL_LABELS_ROW]
    if judged_label_rows["n_galaxies"].sum() == 0:
        judged["agrees"] = judged["agrees"].astype(object)
        judged.loc[judged["label"] == ev.ALL_LABELS_ROW, "agrees"] = None
    tables = [judged, ev.agreement_verdict(rows, reported_labels).iloc[:-1].assign(status="reported")]
    if not_in_reference:
        missing_rows = []
        for label in not_in_reference:
            missing_rows.append({"label": label, "n_galaxies": 0, "inside_95": np.nan,
                                 "median_offset": np.nan, "agrees": None, "status": NOT_IN_REFERENCE})
        tables.append(pd.DataFrame(missing_rows))
    return pd.concat(tables, ignore_index=True)


def photometric_comparison_rows(sfh_model, half, labels):
    """The ML of one model against its own photometric Bagpipes fits in one half, with the draws of cell 8.

    Returns: (rows, n_fits, n_not_converged): reference_comparison_rows over the galaxies whose photometric
    fit converged; the galaxies with a fit and an ML posterior; those of them whose fit did not converge
    (left out).
    """
    prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
    if (sfh_model, half) not in ml_posteriors_by_half or f"{prefix}mw_age_median" not in reference.columns:
        return [], 0, 0
    half_galaxies, half_posteriors = ml_posteriors_by_half[(sfh_model, half)]
    photometric = half_galaxies[f"{prefix}mw_age_median"].notna().to_numpy()
    n_fits = int(photometric.sum())
    if n_fits == 0:
        return [], 0, 0
    converged = converged_fits(half_galaxies, prefix)
    n_not_converged = int((photometric & ~converged).sum())
    compared = photometric & converged
    photometric_galaxies = half_galaxies[compared].reset_index(drop=True)
    photometric_posteriors = {label: half_posteriors[label][compared] for label in labels}
    rows = ev.reference_comparison_rows(photometric_posteriors, photometric_galaxies, labels, prefix=prefix)
    return rows, n_fits, n_not_converged


def fitter_calibration_flag(column):
    """The fitter_calibrated flag a reference set's column holds, or a note when it is empty or missing."""
    if column in reference.columns:
        flags = reference[column].dropna()
        if len(flags) > 0:
            return flags.iloc[0]
    return "unknown (calibration mocks not all fitted, or fits not read)"


verdicts = {}                      # (sfh_model, reference set, half) -> verdict table
for sfh_model in SFH_MODELS:
    prefixes = {"photometric": sm.PHOTOMETRIC_PREFIXES[sfh_model],
                "spectral": sm.SPECTRAL_PREFIXES[sfh_model]}
    judged_labels = LABELS[sfh_model]["judged_labels"]
    reported_labels = LABELS[sfh_model]["reported_labels"]
    reference_labels = LABELS[sfh_model]["reference_labels"]
    labels_by_set = {}
    for reference_set, prefix in prefixes.items():
        labels_by_set[reference_set] = available_labels(reference, reference_labels, prefix)
        missing_labels = [label for label in reference_labels if label not in labels_by_set[reference_set]]
        if missing_labels:
            print(f"{sfh_model}: {REFERENCE_FILE} has no {reference_set} {prefix}<label> columns for "
                  f"{missing_labels}: their rows stay empty")
    # The judged labels each reference set holds. The photometric set judges the full list: a missing
    # column there means the fits are not read yet, so its labels stay judged and empty
    judged_by_set = {"photometric": judged_labels,
                     "spectral": available_labels(reference, judged_labels, prefixes["spectral"])}

    for half in halves:
        n_half = int((reference["split"] == half).sum())

        # 1. The galaxies of this half with a converged photometric fit of this model
        photometric_rows, n_photometric_fits, n_not_converged = photometric_comparison_rows(
            sfh_model, half, labels_by_set["photometric"])
        n_photometric_compared = len({row["specobjid"] for row in photometric_rows})
        print(f"{sfh_model}, {half} half: photometric fits of {n_photometric_compared} of {n_half} galaxies "
              f"compared ({n_photometric_fits} with a fit and an ML posterior, {n_not_converged} not "
              "converged and left out)")

        # 2. The galaxies of this half with a converged spectral fit of this model (cell 8)
        spectral_rows = agreements[(sfh_model, sfh_model)]
        if len(spectral_rows) > 0:
            spectral_rows = spectral_rows[spectral_rows["half"] == half]
        n_spectral_compared = 0
        if len(spectral_rows) > 0:
            n_spectral_compared = spectral_rows["specobjid"].nunique()
        print(f"{sfh_model}, {half} half: spectral fits of {n_spectral_compared} of {n_half} galaxies "
              "compared")
        comparison_rows = {"photometric": photometric_rows, "spectral": spectral_rows}

        # 3. Both verdicts side by side. Why the fitter flag: a reference whose fitter failed its calibration
        # test (or whose calibration mocks are not all fitted) is a weaker reference, so it is shown next to
        # the verdict
        for reference_set, rows in comparison_rows.items():
            if len(rows) == 0 and half != "validation":
                # The test half has photometric fits only with EXPORT_LOCKED_TEST_JOBS in notebook 2
                print(f"  {half} half: no {sfh_model} {reference_set} Bagpipes fits to compare with, so no "
                      f"verdict ({TEST_HALF_HINTS[reference_set]})")
                continue
            if len(rows) == 0:
                print(f"  no {sfh_model} {reference_set} Bagpipes fits in {REFERENCE_FILE} yet: "
                      "run notebook 2 after the fits")
            set_judged_labels = judged_by_set[reference_set]
            not_in_reference = [label for label in judged_labels if label not in set_judged_labels]
            verdict = verdict_table(rows, set_judged_labels, reported_labels, not_in_reference)
            verdicts[(sfh_model, reference_set, half)] = verdict
            fitter_calibrated = fitter_calibration_flag(f"{prefixes[reference_set]}fitter_calibrated")
            print(f"{sfh_model}: agreement with {REFERENCE_SET_NAMES[reference_set]} ({half} half; "
                  f"{', '.join(reported_labels)} reported, not judged); "
                  f"fitter calibrated: {fitter_calibrated}")
            display(verdict.round(3))
            if not_in_reference:
                print(f"  {', '.join(not_in_reference)}: {NOT_IN_REFERENCE}, so not judged here")
            judged_rows = verdict[verdict["status"] == "judged"]
            if judged_rows["n_galaxies"].sum() == 0:
                print(f"  no galaxy for the judged labels {set_judged_labels} in this reference set")
            else:
                all_row = verdict[verdict["label"] == ev.ALL_LABELS_ROW]
                print(f"  the ML agrees on {', '.join(set_judged_labels)}: {bool(all_row['agrees'].iloc[0])}")
```

  **`task15_cells/2e021085.py`** (end of section 4, code: the per-model summary files, now with `half`):

```python
for sfh_model in SFH_MODELS:
    summary_tables = []
    if RUN_LOCKED_TEST:
        summary_tables.append(synthetic_test_summaries[sfh_model])
    summary_tables.append(choices[sfh_model].assign(section="validation_choice"))
    # The works-test verdicts of section 3b, in their order there: per half, photometric then spectral
    for (verdict_model, reference_set, half), verdict in verdicts.items():
        if verdict_model != sfh_model:
            continue
        summary_tables.append(verdict.assign(section=f"agreement_{reference_set}", half=half,
                                             sfh_model=sfh_model))
    summary_file = sfh_models.model_file(SUMMARY_FILE, sfh_model)
    pd.concat(summary_tables, ignore_index=True).to_csv(summary_file, index=False)
    print(f"wrote {summary_file}")
```

  **`task15_cells/3e17feaa.md`** (section 5 header, markdown):

```markdown
## 5. Delayed-τ against double power law

The two SFH models side by side, in five parts (spec 2026-10-28, Part 4 and Amendment 2):
1. each model's ML against its own photometric and its own spectral Bagpipes (the "works" tests of section 3b, side by side);
2. delayed-τ against DPL on the shared labels:
   - the ML of each model, on the validation galaxies and on DR3;
   - the photometric Bagpipes fits of each model, on the validation galaxies;
   - the spectral Bagpipes fits of each model, on the validation galaxies;
3. four estimators against both spectral references (delayed-τ and DPL spectral Bagpipes), passive and star-forming galaxies apart;
4. the photometric Bayesian evidence of each model, galaxy by galaxy and on the cross-fitted mocks;
5. the spectral Bayesian evidence of each model, galaxy by galaxy, and how it relates to the photometric one.

The shared labels are model-independent: `log_mstar`, `mw_age`, `log_metallicity`, `av_ism`, `log_ssfr`, `f_young`, `f_mid`, `f_old` (`sfh_models.SHARED_LABELS`).

Every difference is delayed-τ minus DPL. Only the validation half is used, plus the test half with `RUN_LOCKED_TEST` (LOCKED):
- the spectral comparisons (5.1 spectral, 5.2 spectral pair, 5.3, 5.5) take the test half whenever it is unlocked;
- the photometric Bagpipes and photometric evidence comparisons (5.1 photometric, 5.2 photometric pair, the photometric estimators of 5.3, 5.4, and the photometric side of 5.5) use the test half only when notebook 2 has read its photometric fits (`EXPORT_LOCKED_TEST_JOBS = True` there); otherwise they cover the validation half, and the cell below says so.

The section is report-only: the evidence is reported, not used to choose a model per galaxy.

**Caveat on the DPL shape parameters.** In the DPL, τ is the turnover time of the double power law and `turnover_lookback` the lookback time to it; the SFR peaks there only when α = β (in general the maximum lies at t = τ (β/α)^{1/(α+β)}), so `log_turnover_lookback` is not the age of peak star formation and is not compared with the delayed-τ age.
```

  **`task15_cells/2887477a.md`** (section 5.1 header, markdown):

```markdown
### 5.1 ML against its own photometric and spectral Bagpipes

The verdicts of section 3b, label × model, one table per reference kind and half: first each model against its own photometric fits, then each model against its own spectral fits (the DPL ML against the DPL spectral fits, beside the delayed-τ ML against the delayed-τ spectral fits). Each model's `all` row covers its own calibration labels. The shared labels are reported, and a label one model does not have is empty in its column.
```

  **`task15_cells/b051ebb0.py`** (section 5.1, code: own photometric and own spectral verdicts side by side, per half; `photometric_halves` unchanged):

```python
COMPARED_MODELS = ["delayed", "dblplaw"]   # ln B and every difference below: delayed-tau minus DPL
SHARED_LABELS = sfh_models.SHARED_LABELS
compare_models = set(COMPARED_MODELS) <= set(SFH_MODELS)
# Long tables with a "section" column, written to COMPARISON_FILE at the end
comparison_tables = []
if not compare_models:
    print(f"Section 5 needs both models in SFH_MODELS (now {SFH_MODELS}): skipped")
else:
    # Each model's works-test verdicts of section 3b, side by side: own photometric, then own spectral fits
    for reference_set in ("photometric", "spectral"):
        for half in halves:
            own_verdicts = {}
            for sfh_model in COMPARED_MODELS:
                key = (sfh_model, reference_set, half)
                if key not in verdicts:   # e.g. no test-half photometric fits (section 3b said so)
                    continue
                own_verdicts[sfh_model] = verdicts[key]
                comparison_tables.append(verdicts[key].assign(section=f"own_{reference_set}", half=half,
                                                              sfh_model=sfh_model))
            if not own_verdicts:
                continue
            print(f"Each model's ML against its own {reference_set} Bagpipes ({half} half): n_galaxies, "
                  "inside_95, median offset [combined 1-sigma widths], agrees")
            display(ev.side_by_side_verdicts(own_verdicts).round(3))

    # The halves of the photometric Bagpipes and evidence comparisons (5.2-5.5). The test half has photometric
    # fits only when notebook 2 exported and read them (EXPORT_LOCKED_TEST_JOBS); the validation half always
    # takes part, so its tables exist even before the fits are read
    photometric_halves = ["validation"]
    if RUN_LOCKED_TEST:
        test_galaxies = reference[reference["split"] == "test"]
        test_has_both_fits = (converged_fits(test_galaxies, sm.PHOTOMETRIC_PREFIXES["delayed"])
                              & converged_fits(test_galaxies, sm.PHOTOMETRIC_PREFIXES["dblplaw"]))
        if test_has_both_fits.any():
            photometric_halves.append("test")
            print(f"Photometric Bagpipes and evidence comparisons: validation and test halves "
                  f"({int(test_has_both_fits.sum())} test galaxies with both fits converged)")
        else:
            print("Photometric Bagpipes and evidence comparisons: validation half only. The reference set "
                  "has no test-half photometric fits of both models: set EXPORT_LOCKED_TEST_JOBS = True in "
                  "notebook 2, fit those jobs and read them, then run this notebook again")
    else:
        print("Photometric Bagpipes and evidence comparisons: validation half (the test half is locked)")
```

  **`task15_cells/89339c2d.md`** (section 5.2 header, markdown):

```markdown
### 5.2 Delayed-τ against DPL on the shared labels

For each pair of estimates, the cells below compute three statistics per label, over the galaxies both estimates have:
- the median of the per-galaxy difference of medians (delayed-τ minus DPL);
- its scatter: the NMAD around that median, 1.4826 × the median absolute deviation;
- the share of galaxies whose two 16–84% intervals overlap.

Pairs:
- **ML:** validation galaxies (the draws of cell 8; test galaxies too with `RUN_LOCKED_TEST`), and DR3 galaxies predicted by both catalogues;
- **photometric Bagpipes:** galaxies where both photometric fits converged, in the halves of `photometric_halves`;
- **spectral Bagpipes:** galaxies where both spectral fits converged, in the halves of `halves`. This pair shows how much the assumed SFH shape moves the shared labels when the data are the spectra.
```

  **`task15_cells/79e11a83.py`** (section 5.2, code: the new spectral pair as step 3; the other pairs unchanged):

```python
if compare_models:
    model_pairs = []
    validation_ml_pairs = None
    for half in halves:
        # 1. The ML of each model, on the galaxies with both posteriors (cell 8)
        ml_summaries = []
        for sfh_model in COMPARED_MODELS:
            if (sfh_model, half) not in ml_posteriors_by_half:
                continue
            galaxies, posteriors = ml_posteriors_by_half[(sfh_model, half)]
            summary = pd.DataFrame({"specobjid": galaxies["specobjid"].to_numpy()})
            for label in SHARED_LABELS:
                p16, median, p84 = ev.percentile_summary(posteriors[label])
                summary[f"ml_{sfh_model}_{label}_p16"] = p16
                summary[f"ml_{sfh_model}_{label}_median"] = median
                summary[f"ml_{sfh_model}_{label}_p84"] = p84
            ml_summaries.append(summary)
        if len(ml_summaries) == len(COMPARED_MODELS):
            ml_pairs = ml_summaries[0].merge(ml_summaries[1], on="specobjid")
            rows = ev.model_comparison_rows(ml_pairs, SHARED_LABELS, "ml_delayed_", "ml_dblplaw_")
            model_pairs.append(rows.assign(pair="ML", sample=half))
            if half == "validation":
                validation_ml_pairs = ml_pairs

    # 2. The photometric Bagpipes fits of each model, on the galaxies where both converged
    for half in photometric_halves:
        half_galaxies = reference[reference["split"] == half].reset_index(drop=True)
        both_converged = (converged_fits(half_galaxies, sm.PHOTOMETRIC_PREFIXES["delayed"])
                          & converged_fits(half_galaxies, sm.PHOTOMETRIC_PREFIXES["dblplaw"]))
        rows = ev.model_comparison_rows(half_galaxies[both_converged], SHARED_LABELS,
                                        sm.PHOTOMETRIC_PREFIXES["delayed"],
                                        sm.PHOTOMETRIC_PREFIXES["dblplaw"])
        model_pairs.append(rows.assign(pair="photometric Bagpipes", sample=half))

    # 3. The spectral Bagpipes fits of each model, on the galaxies where both converged. Every half of
    # halves: the spectral fits of both models cover every galaxy with a spectrum. Before notebook 2 has read
    # the DPL spectral fits, the rows have no galaxies
    for half in halves:
        half_galaxies = reference[reference["split"] == half].reset_index(drop=True)
        both_converged = spectral_fits(half_galaxies, "delayed") & spectral_fits(half_galaxies, "dblplaw")
        rows = ev.model_comparison_rows(half_galaxies[both_converged], SHARED_LABELS,
                                        sm.SPECTRAL_PREFIXES["delayed"], sm.SPECTRAL_PREFIXES["dblplaw"])
        model_pairs.append(rows.assign(pair="spectral Bagpipes", sample=half))

    # 4. The ML of each model on DR3: both catalogues have the same rows (cell 12); compared where
    # both predicted
    dr3_pairs = pd.concat([catalogues["delayed"].add_prefix("ml_delayed_"),
                           catalogues["dblplaw"].add_prefix("ml_dblplaw_")], axis=1)
    rows = ev.model_comparison_rows(dr3_pairs, SHARED_LABELS, "ml_delayed_", "ml_dblplaw_")
    model_pairs.append(rows.assign(pair="ML", sample="DR3"))

    model_pairs = pd.concat(model_pairs, ignore_index=True)
    comparison_tables.append(model_pairs.assign(section="model_pairs"))
    print("Delayed-tau minus DPL on the shared labels: median difference, scatter (NMAD around it) and the "
          "share of galaxies whose 16-84% intervals overlap")
    display(model_pairs.set_index(["pair", "sample", "label"]).round(3))

    # 5. ML medians of the two models on the validation galaxies, one panel per shared label
    if validation_ml_pairs is not None:
        figure, axes = plt.subplots(2, 4, figsize=(16, 7))
        for axis, label in zip(axes.flat, SHARED_LABELS):
            delayed_median = validation_ml_pairs[f"ml_delayed_{label}_median"]
            dblplaw_median = validation_ml_pairs[f"ml_dblplaw_{label}_median"]
            axis.scatter(delayed_median, dblplaw_median, s=4)
            low = min(delayed_median.min(), dblplaw_median.min())
            high = max(delayed_median.max(), dblplaw_median.max())
            axis.plot([low, high], [low, high], color="black", lw=1)
            axis.set(xlabel=f"delayed-τ {label}", ylabel=f"DPL {label}")
        figure.suptitle("ML posterior medians, validation galaxies")
        plt.tight_layout()
        plt.show()
```

  **`task15_cells/10a86e1d.md`** (section 5.3 header, markdown):

```markdown
### 5.3 Against the spectra

The "works" test of four estimators against each of the two spectral references (the delayed-τ and the DPL spectral Bagpipes fits), on the shared labels: the ML of each model and the photometric Bagpipes fits of each model.
- Passive and star-forming galaxies are judged apart, split at the delayed-τ spectral log sSFR of `PASSIVE_SSFR_LIMIT` (as in section 3), for both references.
- A reference judges the galaxies whose spectral fit of its model converged, on the shared labels it holds.
- The photometric Bagpipes fits keep only their 16th, 50th and 84th percentiles. Their 95% interval is therefore extended from these as a Gaussian, each side by its own width (`ev.compare_summaries_with_reference`).

**Caveat: each spectral reference assumes its own SFH shape.** The reference `mw_age` and mass fractions `f_young`, `f_mid`, `f_old` carry the delayed-τ shape in the delayed-τ fits and the DPL shape in the DPL fits. An estimate that agrees with only one reference on them may describe the galaxy with the other history rather than be wrong, and the spectral pair of section 5.2 shows how far apart the two references are. `log_mstar`, `log_metallicity`, `av_ism` and `log_ssfr` depend less on the shape.
```

  **`task15_cells/0267e995.py`** (section 5.3, code: the four estimators against both spectral references):

```python
if compare_models:
    ML_ESTIMATORS = []
    PHOTOMETRIC_ESTIMATORS = []
    for sfh_model in COMPARED_MODELS:
        ML_ESTIMATORS.append(f"ML {sfh_model}")
        PHOTOMETRIC_ESTIMATORS.append(f"photometric Bagpipes {sfh_model}")

    # 1. One row per (spectral reference, estimator, galaxy, label), with the galaxy type
    spectral_tables = []
    for spectral_model in SPECTRAL_REFERENCES:
        spectral_prefix = sm.SPECTRAL_PREFIXES[spectral_model]
        for half in halves:
            # The ML rows of cell 8 (they carry spectral_model, half and galaxy_type)
            for sfh_model in COMPARED_MODELS:
                agreement = agreements[(sfh_model, spectral_model)]
                if len(agreement) == 0:
                    continue
                ml_rows = agreement[(agreement["half"] == half) & agreement["label"].isin(SHARED_LABELS)]
                spectral_tables.append(ml_rows.assign(estimator=f"ML {sfh_model}"))
            if half not in photometric_halves:
                continue
            # The photometric Bagpipes rows, on the galaxies with this reference's converged spectral fit
            half_galaxies = reference[reference["split"] == half].reset_index(drop=True)
            spectral_galaxies = half_galaxies[spectral_fits(half_galaxies, spectral_model)]
            spectral_galaxies = spectral_galaxies.reset_index(drop=True)
            galaxy_types = dict(zip(spectral_galaxies["specobjid"], spectral_galaxy_types(spectral_galaxies)))
            for sfh_model in COMPARED_MODELS:
                prefix = sm.PHOTOMETRIC_PREFIXES[sfh_model]
                fitted = spectral_galaxies[converged_fits(spectral_galaxies, prefix)].reset_index(drop=True)
                # The shared labels both the photometric fit and this spectral reference hold
                labels = available_labels(fitted, available_labels(fitted, SHARED_LABELS, prefix),
                                          spectral_prefix)
                photometric_rows = pd.DataFrame(ev.summary_comparison_rows(fitted, labels, prefix,
                                                                           spectral_prefix))
                if len(photometric_rows) == 0:
                    continue
                photometric_rows["spectral_model"] = spectral_model
                photometric_rows["half"] = half
                photometric_rows["galaxy_type"] = photometric_rows["specobjid"].map(galaxy_types)
                spectral_tables.append(photometric_rows.assign(estimator=f"photometric Bagpipes {sfh_model}"))
    spectral_rows = pd.DataFrame(columns=["spectral_model", "half", "galaxy_type", "estimator", "label",
                                          "offset_in_widths", "inside_95"])
    if spectral_tables:
        spectral_rows = pd.concat(spectral_tables, ignore_index=True)

    # 2. The verdict of each estimator against each reference, passive and star-forming apart, on the shared
    # labels the reference holds. A half without photometric fits (see photometric_halves) judges the ML
    # estimators only
    spectral_verdicts = []
    for spectral_model in SPECTRAL_REFERENCES:
        spectral_shared_labels = available_labels(reference, SHARED_LABELS,
                                                  sm.SPECTRAL_PREFIXES[spectral_model])
        if not spectral_shared_labels:
            print(f"No {spectral_model} spectral Bagpipes fits in {REFERENCE_FILE} yet: run notebook 2 after "
                  "the fits. Comparison against them skipped")
            continue
        for half in halves:
            half_estimators = list(ML_ESTIMATORS)
            if half in photometric_halves:
                half_estimators.extend(PHOTOMETRIC_ESTIMATORS)
            for galaxy_type in ("passive", "star-forming"):
                by_estimator = {}
                for estimator in half_estimators:
                    rows = spectral_rows[(spectral_rows["spectral_model"] == spectral_model)
                                         & (spectral_rows["half"] == half)
                                         & (spectral_rows["galaxy_type"] == galaxy_type)
                                         & (spectral_rows["estimator"] == estimator)]
                    by_estimator[estimator] = ev.agreement_verdict(rows, spectral_shared_labels)
                    spectral_verdicts.append(by_estimator[estimator].assign(
                        spectral_model=spectral_model, half=half, galaxy_type=galaxy_type,
                        estimator=estimator))
                print(f"Against the {spectral_model} spectral Bagpipes fits (see the caveat): {half} half, "
                      f"{galaxy_type} galaxies")
                display(ev.side_by_side_verdicts(by_estimator).round(3))
    if spectral_verdicts:
        comparison_tables.append(pd.concat(spectral_verdicts, ignore_index=True)
                                 .assign(section="against_spectra"))
```

  **`task15_cells/a7c3e5d1.md`** (new section 5.5 header, markdown):

```markdown
### 5.5 Spectral Bayesian evidence

The spectral ln B = `spec_log_evidence` − `spec_dblplaw_log_evidence` for each galaxy whose two spectral fits converged, in the halves of `halves`; positive favours delayed-τ. Both fits see the same spectrum through the same likelihood (calibration polynomial, noise scale and MultiNest settings), so ln B compares the two SFH priors on these data. The spectra carry far more information than 12 colours, so this tests whether they prefer one SFH shape. The categories are those of section 5.4 (`ev.evidence_category_fractions`: below 1 inconclusive, 1–3 positive, above 3 strong), overall and by galaxy type.

**Uncertainty.** The DPL spectral fits store MultiNest's ln Z error (`spec_dblplaw_log_evidence_error`). The delayed-τ spectral fits made before the error was stored have none (`spec_log_evidence_error` is empty), so their ln B error is unknown. ln Z differences between spectral models are expected to dwarf the usual 0.1–0.3 nats.

**Against the photometric ln B.** For the galaxies with both ln B (the halves of `photometric_halves`):
- the share whose two ln B have the same sign (`ev.log_bayes_sign_agreement`), and the share whose two ln B favour the same model, an inconclusive one counting as favouring neither;
- the counts by the model each favours (`ev.favoured_model_crosstab`, rows photometric, columns spectral);
- a scatter plot of the two.

The spectral ln B axis is drawn as sign(ln B) log10(1 + |ln B|), since spectral ln B may reach hundreds of nats.
```

  **`task15_cells/b8d4f6e2.py`** (new section 5.5, code: the spectral ln B, its categories by galaxy type, and its relation to the photometric ln B of cell 22):

```python
SPECTRAL_EVIDENCE_COLUMNS = {"delayed": "spec_log_evidence", "dblplaw": "spec_dblplaw_log_evidence"}
EVIDENCE_ERROR_SUFFIX = "_error"   # spec_log_evidence_error, spec_dblplaw_log_evidence_error


def signed_log(values):
    """sign(x) log10(1 + |x|): a monotone axis for ln B from 0.1 to thousands of nats, either sign."""
    values = np.asarray(values, dtype=float)
    return np.sign(values) * np.log10(1 + np.abs(values))


def log_bayes_summary(table, column, halves_used):
    """Percentiles and evidence categories of one ln B column, per half, overall and by galaxy type.

    Inputs: table: DataFrame with half, galaxy_type and the ln B column [nats]; halves_used: halves, in order.
    Returns: (distribution, categories): distribution has one row per (half, group) with n_galaxies and
    ln_B_p16 / ln_B_median / ln_B_p84 of the finite ln B [nats]; categories stacks
    ev.evidence_category_fractions per (half, group).
    """
    distribution_rows = []
    category_tables = []
    for half in halves_used:
        in_half = table["half"] == half
        groups = {"all": in_half, "passive": in_half & (table["galaxy_type"] == "passive"),
                  "star-forming": in_half & (table["galaxy_type"] == "star-forming")}
        for group, in_group in groups.items():
            log_bayes = table.loc[in_group, column].to_numpy(dtype=float)
            finite = log_bayes[np.isfinite(log_bayes)]
            percentiles = np.full(3, np.nan)
            if len(finite) > 0:
                percentiles = np.percentile(finite, [16, 50, 84])
            distribution_rows.append({"half": half, "group": group, "n_galaxies": len(finite),
                                      "ln_B_p16": percentiles[0], "ln_B_median": percentiles[1],
                                      "ln_B_p84": percentiles[2]})
            category_tables.append(ev.evidence_category_fractions(log_bayes).assign(half=half, group=group))
    return pd.DataFrame(distribution_rows), pd.concat(category_tables, ignore_index=True)


if compare_models:
    # 1. Spectral ln B and its uncertainty per galaxy, where both spectral fits converged
    spectral_evidence_columns = []   # ln Z of each model, then their errors (columns of REFERENCE_FILE)
    for suffix in ("", EVIDENCE_ERROR_SUFFIX):
        for sfh_model in COMPARED_MODELS:
            spectral_evidence_columns.append(f"{SPECTRAL_EVIDENCE_COLUMNS[sfh_model]}{suffix}")
    spectral_evidence_tables = []
    for half in halves:
        half_galaxies = reference[reference["split"] == half].reset_index(drop=True)
        both_converged = spectral_fits(half_galaxies, "delayed") & spectral_fits(half_galaxies, "dblplaw")
        galaxies = half_galaxies[both_converged].reset_index(drop=True)
        table = pd.DataFrame({"specobjid": galaxies["specobjid"], "half": half,
                              "galaxy_type": spectral_galaxy_types(galaxies)})
        # Why reindex: a column the reference set lacks (no DPL spectral fits read yet, or delayed-tau fits
        # without a stored ln Z error) becomes NaN instead of a KeyError
        log_evidence_delayed, log_evidence_dblplaw, error_delayed, error_dblplaw = (
            galaxies.reindex(columns=spectral_evidence_columns).to_numpy(dtype=float).T)
        table["spectral_log_bayes"] = ev.log_bayes_factor(log_evidence_delayed, log_evidence_dblplaw)
        table["spectral_log_bayes_error"] = ev.log_bayes_factor_error(error_delayed, error_dblplaw)
        spectral_evidence_tables.append(table)
    spectral_evidence = pd.concat(spectral_evidence_tables, ignore_index=True)

    if not np.isfinite(spectral_evidence["spectral_log_bayes"].to_numpy(dtype=float)).any():
        print(f"No galaxy with both spectral ln Z in {REFERENCE_FILE} yet: run notebook 2 after the DPL "
              "spectral fits. Spectral evidence skipped")
    else:
        # 2. The distribution of the spectral ln B and the share of each category, overall and by galaxy type
        spectral_distribution, spectral_categories = log_bayes_summary(spectral_evidence,
                                                                       "spectral_log_bayes", halves)
        comparison_tables.append(spectral_distribution.assign(section="spectral_evidence_distribution"))
        comparison_tables.append(spectral_categories.assign(section="spectral_evidence_categories"))
        print("Spectral ln B = ln Z(delayed-tau spectral) - ln Z(DPL spectral) per galaxy "
              "(positive favours delayed-tau)")
        display(spectral_distribution.round(2))
        display(spectral_categories.set_index(["half", "group", "favoured_model", "strength"])["fraction"]
                .unstack(["favoured_model", "strength"]).round(3))

        figure, axis = plt.subplots(figsize=(7, 4))
        validation_spectral = spectral_evidence[spectral_evidence["half"] == "validation"]
        for group in ("passive", "star-forming"):
            values = validation_spectral.loc[validation_spectral["galaxy_type"] == group,
                                             "spectral_log_bayes"].dropna()
            axis.hist(signed_log(values), bins=30, histtype="step", label=f"{group} ({len(values)})")
        for edge in (-ev.STRONG_LOG_BAYES, -ev.INCONCLUSIVE_LOG_BAYES, ev.INCONCLUSIVE_LOG_BAYES,
                     ev.STRONG_LOG_BAYES):
            axis.axvline(signed_log(edge), color="grey", ls=":")
        axis.set(xlabel="sign(ln B) log10(1 + |ln B|), spectral ln B = ln Z(delayed-τ) − ln Z(DPL)",
                 ylabel="galaxies",
                 title="Spectral evidence per validation galaxy (positive favours delayed-τ)")
        axis.legend()
        plt.show()

        # 3. The spectral ln B uncertainty, where both fits stored their MultiNest ln Z error
        with_error = spectral_evidence["spectral_log_bayes_error"].notna()
        if with_error.any():
            errors = spectral_evidence.loc[with_error, "spectral_log_bayes_error"]
            below_error = spectral_evidence.loc[with_error, "spectral_log_bayes"].abs() < errors
            print(f"Spectral ln B uncertainty: median {errors.median():.2f} over {int(with_error.sum())} "
                  f"galaxies; {below_error.mean():.0%} of them have |ln B| below it")
        else:
            print("No spectral ln B uncertainty: the delayed-tau spectral fits store no ln Z error "
                  f"(spec_log_evidence_error is empty); the median |spectral ln B| is "
                  f"{spectral_evidence['spectral_log_bayes'].abs().median():.1f} nats")

        # 4. Against the photometric ln B of section 5.4, on the galaxies with both
        both_evidence = evidence.merge(spectral_evidence[["specobjid", "half", "spectral_log_bayes"]],
                                       on=["specobjid", "half"])
        agreement_rows = []
        crosstab_tables = []
        for half in photometric_halves:
            in_half = both_evidence["half"] == half
            groups = {"all": in_half, "passive": in_half & (both_evidence["galaxy_type"] == "passive"),
                      "star-forming": in_half & (both_evidence["galaxy_type"] == "star-forming")}
            for group, in_group in groups.items():
                photometric_log_bayes = both_evidence.loc[in_group, "log_bayes"]
                spectral_log_bayes = both_evidence.loc[in_group, "spectral_log_bayes"]
                sign_agreement = ev.log_bayes_sign_agreement(photometric_log_bayes, spectral_log_bayes)
                agreement_rows.append({"half": half, "group": group, **sign_agreement})
                crosstab = ev.favoured_model_crosstab(photometric_log_bayes, spectral_log_bayes)
                crosstab_rows = crosstab.stack().rename("n_galaxies").reset_index()
                crosstab_rows = crosstab_rows.rename(columns={"first": "photometric_favours",
                                                              "second": "spectral_favours"})
                crosstab_tables.append(crosstab_rows.assign(half=half, group=group))
                if group == "all":
                    print(f"Galaxies by the model each ln B favours ({half} half; rows photometric, "
                          "columns spectral)")
                    display(crosstab.rename_axis(index="photometric", columns="spectral"))
        evidence_agreement = pd.DataFrame(agreement_rows)
        comparison_tables.append(evidence_agreement.assign(section="evidence_agreement"))
        comparison_tables.append(pd.concat(crosstab_tables, ignore_index=True)
                                 .assign(section="evidence_favoured_models"))
        print("Photometric against spectral ln B: share with the same sign and share favouring the "
              "same model")
        display(evidence_agreement.round(3))

        validation_both = both_evidence[both_evidence["half"] == "validation"]
        if len(validation_both) > 0:
            figure, axis = plt.subplots(figsize=(6, 5))
            for group in ("passive", "star-forming"):
                in_group = validation_both["galaxy_type"] == group
                axis.scatter(validation_both.loc[in_group, "log_bayes"],
                             signed_log(validation_both.loc[in_group, "spectral_log_bayes"]), s=6,
                             label=f"{group} ({int(in_group.sum())})")
            axis.axhline(0, color="grey", lw=1)
            axis.axvline(0, color="grey", lw=1)
            axis.set(xlabel="photometric ln B [nats]", ylabel="sign(ln B) log10(1 + |ln B|), spectral",
                     title="Evidence per validation galaxy (positive favours delayed-τ)")
            axis.legend()
            plt.show()
```

- [ ] **Step 6: Put the cells into the notebook.** Write this script to `Synthetic_domain/task15_cells/edit_notebook.py`:

```python
"""Task 15: put the new cell sources of task15_cells/ into ml_model.ipynb (run from Synthetic_domain/).

Each file task15_cells/<cell id>.py (code) or .md (markdown) is the complete new source of that cell.
REPLACED cells keep their id and type; NEW cells are inserted, in order, after the cell INSERT_AFTER.
Edited code cells lose their stale outputs. The file is written as Jupyter writes it (indent 1, UTF-8,
trailing newline), so the diff shows only the changed cells.
"""

import json
import os
import sys

NOTEBOOK = sys.argv[1] if len(sys.argv) > 1 else "ml_model.ipynb"
CELLS_DIRECTORY = sys.argv[2] if len(sys.argv) > 2 else "task15_cells"
REPLACED = ["55011031", "fe480605", "e81a48b9", "e05c91f8", "74230b0e", "2e021085", "3e17feaa",
            "2887477a", "b051ebb0", "89339c2d", "79e11a83", "10a86e1d", "0267e995"]
INSERT_AFTER = "fc196542"           # the cross-fit cell of section 5.4
NEW = [("a7c3e5d1", "markdown"), ("b8d4f6e2", "code")]   # section 5.5, before the COMPARISON_FILE cell
EXTENSIONS = {"code": ".py", "markdown": ".md"}


def cell_source(cell_id, cell_type):
    """The new source of one cell as notebook lines (each ends with a newline except the last)."""
    path = os.path.join(CELLS_DIRECTORY, cell_id + EXTENSIONS[cell_type])
    with open(path, encoding="utf-8") as source_file:
        text = source_file.read().rstrip("\n")
    return text.splitlines(keepends=True)


with open(NOTEBOOK, encoding="utf-8") as notebook_file:
    notebook = json.load(notebook_file)
cells_by_id = {cell["id"]: cell for cell in notebook["cells"]}
for cell_id in REPLACED:
    cell = cells_by_id[cell_id]
    cell["source"] = cell_source(cell_id, cell["cell_type"])
    if cell["cell_type"] == "code":
        cell["outputs"] = []
        cell["execution_count"] = None

assert not set(cell_id for cell_id, _ in NEW) & set(cells_by_id), "a new cell id is already in the notebook"
insert_at = [cell["id"] for cell in notebook["cells"]].index(INSERT_AFTER) + 1
for offset, (cell_id, cell_type) in enumerate(NEW):
    if cell_type == "code":
        new_cell = {"cell_type": "code", "execution_count": None, "id": cell_id, "metadata": {},
                    "outputs": [], "source": cell_source(cell_id, cell_type)}
    else:
        new_cell = {"cell_type": "markdown", "id": cell_id, "metadata": {},
                    "source": cell_source(cell_id, cell_type)}
    notebook["cells"].insert(insert_at + offset, new_cell)

with open(NOTEBOOK, "w", encoding="utf-8") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print(f"{NOTEBOOK}: {len(REPLACED)} cells replaced, {len(NEW)} inserted, {len(notebook['cells'])} cells")
```

  - Run: `cd Synthetic_domain && ../.venv/bin/python task15_cells/edit_notebook.py`
  - Expected: `ml_model.ipynb: 13 cells replaced, 2 inserted, 27 cells`. A `KeyError` names a cell id that is not in the notebook: stop and report (the notebook is not the one this task was written against).

- [ ] **Step 7: Validate and compile; never execute.**
  1. JSON, nbformat, every code cell compiles, no undefined name across cells:

```bash
cd Synthetic_domain && cells_file="$(mktemp --suffix=.py)" && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python - "$cells_file" <<'EOF'
import sys

import nbformat

notebook = nbformat.read("ml_model.ipynb", as_version=4)
nbformat.validate(notebook)
sources = []
for index, cell in enumerate(notebook.cells):
    if cell.cell_type == "code":
        compile(cell.source, f"cell {index}", "exec")
        sources.append(cell.source)
with open(sys.argv[1], "w") as output:
    output.write("\n\n".join(sources))
print(len(notebook.cells), "cells; valid; every code cell compiles")
EOF
../.venv/bin/flake8 --select=F821 --builtins=display "$cells_file" && echo "no undefined names"
awk 'length > 110 {print FILENAME ": " FNR ": " length " characters"}' task15_cells/*.py
```

  - Expected: `27 cells; valid; every code cell compiles`, then `no undefined names`, and no line from `awk`.

  2. Every module name the notebooks use exists (this needs Task 13's `sm.SPECTRAL_PREFIXES` and Step 3's functions):
     - Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_notebooks.py -q`
     - Expected: `3 passed`.

  3. Only the intended places changed:

```bash
cd Synthetic_domain && git diff --stat -- ml_model.ipynb && grep -c '"id": "' ml_model.ipynb \
    && grep -n -e 'agreements\[sfh_model\]' -e 'verdicts\[(sfh_model, "photometric")\]' \
       -e 'verdicts\[(sfh_model, reference_set)\]' ml_model.ipynb
```

  - Expected: one changed file; `27`; and no line from the last `grep`, which therefore exits with status 1 (no old key of `agreements` or `verdicts` is left; the current notebook has 5 such lines).

- [ ] **Step 8: Stub check of the changed tables (CPU, about 30 s, outside the repository; needs the controller's approval, see the top of this task).** It runs the table code of the changed cells with stubbed posteriors: no notebook is executed with nbconvert, and there is no Bagpipes or flow work. It runs the code cells of sections 3, 3b and 5 of the old and the new notebook with the same fake ML posteriors (`labelled_posteriors` and `ml_networks.load_run` are stubs; sections 1, 2 and 4 are skipped, the DR3 catalogues are small random tables). It checks that every delayed-τ table is unchanged, and that the new tables appear with fake DPL spectral columns and stay empty without them, unlocked and locked. It writes only into its own folder.
  Every command of this step runs from the repository root, in one shell (it keeps `$check_dir`). The check reads the laptop's `Synthetic_domain/spectroscopic_reference.csv` and `photometric_cross_fits.csv`, which exist since Task 12.
  1. Make a scratch folder outside the repository and save the old notebook there (HEAD still holds it: Step 4 committed only `ml_evaluation.py` and its tests):

```bash
check_dir="$(mktemp -d)" && echo "$check_dir"
git show HEAD:Synthetic_domain/ml_model.ipynb > "$check_dir/ml_model_before_task15.ipynb"
```

  2. Write this script to `"$check_dir"/task15_stub_check.py`:

```python
"""Task 15 stub check: run ml_model.ipynb sections 3, 3b and 5 on fake ML posteriors, old cells against new.

No flow, no GPU, no DR3: labelled_posteriors and ml_networks.load_run are stubs, the DR3 catalogues are
small random tables, and sections 1, 2 and 4 are skipped. The same stub draws feed the old and the new
notebook, so every delayed-tau table (and the DPL tables against the delayed-tau spectra) must be equal.
Usage: python task15_stub_check.py <Synthetic_domain directory> <old ml_model.ipynb> <new ml_model.ipynb>
Writes its CSV files (reference copies, the notebook's outputs) into run/ next to this script.
"""

import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

REPO, OLD_NOTEBOOK, NEW_NOTEBOOK = (os.path.abspath(path) for path in sys.argv[1:4])
RUN_DIRECTORY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run")
sys.path.insert(0, REPO)
os.chdir(REPO)                                   # synthetic_model reads its grids from here at import
import ml_data  # noqa: E402
import ml_networks  # noqa: E402
import sfh_models  # noqa: E402

os.makedirs(RUN_DIRECTORY, exist_ok=True)
shutil.copy(os.path.join(REPO, "photometric_cross_fits.csv"), RUN_DIRECTORY)
os.chdir(RUN_DIRECTORY)                          # the notebook writes its outputs here, not in the repo
STUB_SAMPLES = 60
SKIPPED_CELLS = {"b6229f82", "59be2b0a", "19fe68be", "d51a8bf8", "46b23dbe"}   # run by hand, or not at all


def stub_load_run(directory, sfh_model, device="cpu"):
    """ml_networks.load_run without a checkpoint: only the config fields the notebook reads."""
    return None, {"sfh_model": sfh_model, "n_targets": len(ml_data.target_names(sfh_model))}, None


def stub_labelled_posteriors(run, observations, n_samples, seed, labels=None):
    """Random posteriors of every training label; the same draws for the same model and galaxy count."""
    sfh_model = run[1]["sfh_model"]
    n_valid = int(observations["valid"].sum())
    rng = np.random.default_rng([n_valid, len(sfh_model)])
    posteriors = {}
    for label in sfh_models.MODELS[sfh_model].training_labels:
        centre = rng.normal(0, 1, size=(n_valid, 1))
        posteriors[label] = centre + rng.normal(0, 1, size=(n_valid, n_samples))
    return posteriors, None


def run_notebook(notebook_path, reference_path, locked):
    """Execute the notebook's code cells (except SKIPPED_CELLS) with the stubs; returns its namespace."""
    ml_networks.load_run = stub_load_run
    notebook = json.load(open(notebook_path))
    sources = {cell["id"]: "".join(cell["source"]) for cell in notebook["cells"]}
    code_ids = [cell["id"] for cell in notebook["cells"] if cell["cell_type"] == "code"]
    namespace = {"display": lambda value: None, "__name__": "__main__"}
    exec(sources["b6229f82"].replace("RUN_LOCKED_TEST = False", f"RUN_LOCKED_TEST = {locked}"), namespace)
    namespace.update(REFERENCE_FILE=reference_path, N_POSTERIOR_SAMPLES=STUB_SAMPLES,
                     labelled_posteriors=stub_labelled_posteriors)
    namespace["choices"] = {model: pd.DataFrame({"label": ["stub"]}) for model in namespace["SFH_MODELS"]}
    namespace["synthetic_test_summaries"] = dict(namespace["choices"])
    rng = np.random.default_rng(3)
    namespace["catalogues"] = {}
    for sfh_model in namespace["SFH_MODELS"]:
        catalogue = {}
        for label in sfh_models.SHARED_LABELS:
            median = rng.normal(size=50)
            catalogue.update({f"{label}_p16": median - 1, f"{label}_median": median,
                              f"{label}_p84": median + 1})
        namespace["catalogues"][sfh_model] = pd.DataFrame(catalogue)
    plt = namespace["plt"]
    plt.show = lambda: plt.close("all")
    for cell_id in code_ids:
        if cell_id not in SKIPPED_CELLS:
            exec(sources[cell_id], namespace)
    return namespace


def same(first, second):
    """Assert two tables hold the same values (compared as text, so None, NaN and dtypes match exactly)."""
    pd.testing.assert_frame_equal(first.reset_index(drop=True).astype(str),
                                  second.reset_index(drop=True).astype(str))


def tables_by_section(namespace):
    """The notebook's comparison_tables, grouped by their section column."""
    tables = {}
    for table in namespace["comparison_tables"]:
        tables.setdefault(table["section"].iloc[0], []).append(table)
    return tables


# 1. Two reference sets: without any DPL spectral column (the old CSV), and with fake DPL spectral fits
reference = pd.read_csv(os.path.join(REPO, "spectroscopic_reference.csv"), dtype={"specobjid": str})
plain = reference.drop(columns=[column for column in reference.columns if column.startswith("spec_")])
plain.to_csv("reference_plain.csv", index=False)
rng = np.random.default_rng(5)
fitted = (plain["fit_status"] == "fitted").to_numpy()
fake = {}
for label in sfh_models.MODELS["dblplaw"].training_labels:
    median = np.where(fitted, rng.normal(size=len(plain)), np.nan)
    fake.update({f"spec_dblplaw_{label}_p16": median - 1, f"spec_dblplaw_{label}_median": median,
                 f"spec_dblplaw_{label}_p84": median + 1})
converged = rng.random(len(plain)) < 0.9
fake["spec_dblplaw_converged"] = pd.array([bool(c) if f else None for c, f in zip(converged, fitted)],
                                          dtype=object)
fake["spec_dblplaw_fitter_calibrated"] = np.where(fitted, True, None)
fake["spec_log_evidence"] = np.where(fitted, rng.normal(-500, 50, len(plain)), np.nan)
fake["spec_log_evidence_error"] = np.nan
fake["spec_dblplaw_log_evidence"] = np.where(fitted, rng.normal(-500, 50, len(plain)), np.nan)
fake["spec_dblplaw_log_evidence_error"] = np.where(fitted, 0.2, np.nan)
with_dpl = pd.concat([plain, pd.DataFrame(fake, index=plain.index)], axis=1)
with_dpl.to_csv("reference_with_dpl.csv", index=False)

# 2. Old cells against new cells, unlocked and locked
for locked in (False, True):
    old = run_notebook(OLD_NOTEBOOK, "reference_plain.csv", locked)
    old_tables = {section: tables[0] for section, tables in tables_by_section(old).items()}
    for reference_file in ("reference_plain.csv", "reference_with_dpl.csv"):
        new = run_notebook(NEW_NOTEBOOK, reference_file, locked)
        new_tables = tables_by_section(new)
        for sfh_model in ("delayed", "dblplaw"):
            same(old["verdicts"][(sfh_model, "photometric")],
                 new["verdicts"][(sfh_model, "photometric", "validation")])
            same(old["agreements"][sfh_model],
                 new["agreements"][(sfh_model, "delayed")].drop(columns="spectral_model"))
        same(old["verdicts"][("delayed", "spectral")], new["verdicts"][("delayed", "spectral", "validation")])
        same(old["model_pairs"], new["model_pairs"][new["model_pairs"]["pair"] != "spectral Bagpipes"])
        new_against_spectra = new_tables["against_spectra"][0]
        new_against_delayed = new_against_spectra[new_against_spectra["spectral_model"] == "delayed"]
        same(old_tables["against_spectra"], new_against_delayed.drop(columns="spectral_model"))
        for section in ("evidence_distribution", "evidence_categories", "true_model_recovery"):
            same(old_tables[section], new_tables[section][0])
        dpl_own_spectral = new["verdicts"][("dblplaw", "spectral", "validation")]
        all_row = dpl_own_spectral[dpl_own_spectral["label"] == "all"]["agrees"].iloc[0]
        if reference_file == "reference_plain.csv":
            assert all_row is None, all_row                  # no DPL spectral fit: no verdict
            assert "spectral_evidence_distribution" not in new_tables
        else:
            assert all_row in (True, False), all_row
            assert dpl_own_spectral["n_galaxies"].iloc[0] > 0
            for section in ("own_spectral", "spectral_evidence_distribution", "spectral_evidence_categories",
                            "evidence_agreement", "evidence_favoured_models"):
                assert section in new_tables, section
            assert set(new_tables["against_spectra"][0]["spectral_model"]) == {"delayed", "dblplaw"}
            assert (new["model_pairs"]["pair"] == "spectral Bagpipes").any()
        expected_halves = {"validation", "test"} if locked else {"validation"}
        assert {key[2] for key in new["verdicts"] if key[1] == "spectral"} == expected_halves
        print(f"locked={locked}, {reference_file}: delayed-tau tables reproduce; new sections present")
print("stub check passed")
```

  3. Run it:

```bash
MPLBACKEND=Agg LD_LIBRARY_PATH=.venv/lib .venv/bin/python -W ignore "$check_dir/task15_stub_check.py" \
    "$PWD/Synthetic_domain" "$check_dir/ml_model_before_task15.ipynb" "$PWD/Synthetic_domain/ml_model.ipynb" \
    | grep -e "reproduce" -e "stub check"
git status --short Synthetic_domain
```

  - Expected:

```text
locked=False, reference_plain.csv: delayed-tau tables reproduce; new sections present
locked=False, reference_with_dpl.csv: delayed-tau tables reproduce; new sections present
locked=True, reference_plain.csv: delayed-tau tables reproduce; new sections present
locked=True, reference_with_dpl.csv: delayed-tau tables reproduce; new sections present
stub check passed
```

  - `git status` lists only `ml_model.ipynb` (modified), `task15_cells/` and the untracked files that were there before (`__pycache__/`, `pipes/`, the two training logs): the check writes nothing into the repository.
  - An `AssertionError` in `pd.testing.assert_frame_equal` names the table and column that changed: fix the cell, rerun Step 6 on a fresh copy (`git checkout Synthetic_domain/ml_model.ipynb`, then Step 6), Step 7 and this step.
  - Remove the scratch folder: `rm -r "$check_dir"`.

- [ ] **Step 9: Run the whole suite.**
  - Run: `cd Synthetic_domain && LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -q`
  - Expected: all pass: the total after Task 14 plus the 3 tests of Step 1, with the 12 Bagpipes overflow warnings and the self-skipping MultiNest smoke test. Passed + skipped is 285 (282 after Task 14, plus 3); tests that need git-ignored data skip themselves where it is missing.

- [ ] **Step 10: Remove the temporary folder; commit.**

```bash
rm -r Synthetic_domain/task15_cells
git add Synthetic_domain/ml_model.ipynb
git commit -m "feat: ml_model.ipynb judges each model against its own spectral fits, compares all estimators
with both spectral references, reports the spectral evidence against the photometric one, and judges the
test half in the locked test

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

**Handoff to the controller (PC runs, not part of this task).** The order is that of Amendment 2, item 6. First the "PC order after Task 14" at the end of Task 14:
- (a) the remaining delayed-τ spectral fits;
- (b) `EXPORT_LOCKED_TEST_JOBS = True` in notebook 2, committed;
- (c) the notebook-2 export run;
- (d) the DPL spectral fits, then the photometric fits of both models for the new validation and the test galaxies;
- (e) the notebook-2 reading run.

Then this notebook with the locked test, from `Synthetic_domain/`. The locked test runs once, so its inputs are checked before it runs:

0. **Pre-flight check of the fetched `spectroscopic_reference.csv`.** A reading run without one of the fit folders does not fail; it leaves that model's columns empty, and the locked run would then only print messages such as "dblplaw ML: no converged dblplaw spectral fits to compare with", "No dblplaw spectral Bagpipes fits in ..." or "validation half only" instead of failing. Save as `$SCRATCH/locked_test_preflight.py` and run `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python "$SCRATCH/locked_test_preflight.py"`:

```python
"""Pre-flight check of the locked test: spectroscopic_reference.csv must hold every input of both models."""
import sys

import numpy as np
import pandas as pd

MIN_FIT_FRACTION = 0.9   # converged fits per half and model, as a fraction of the half's galaxies
HALVES = ["validation", "test"]
VERDICT_COLUMNS = ["fitter_calibrated", "spec_dblplaw_fitter_calibrated", "phot_fitter_calibrated",
                   "phot_dblplaw_fitter_calibrated"]
CONVERGED_COLUMNS = ["spec_dblplaw_converged", "phot_converged", "phot_dblplaw_converged"]
EVIDENCE_COLUMNS = ["spec_log_evidence", "spec_dblplaw_log_evidence"]

reference = pd.read_csv("spectroscopic_reference.csv")

# 1. Every column the locked run reads; a missing one means notebook 2 did not read that fit folder
missing_columns = [column for column in ["split", *VERDICT_COLUMNS, *CONVERGED_COLUMNS, *EVIDENCE_COLUMNS]
                   if column not in reference.columns]
if missing_columns:
    sys.exit(f"STOP: spectroscopic_reference.csv lacks {missing_columns}; do not run the locked test")

failures = []


def check(passed, message):
    """Print one check and remember a failure."""
    print(("ok    " if passed else "FAIL  ") + message)
    if not passed:
        failures.append(message)


# 2. Fitter verdicts: True or False is a result; empty means a calibration test was incomplete
for column in VERDICT_COLUMNS:
    check(bool(reference[column].notna().all()), f"{column} filled ({reference[column].iloc[0]})")

# 3. Per half: spectral fits of both models with their ln Z, and photometric fits of both models
for half in HALVES:
    galaxies = reference[reference["split"] == half]
    n_galaxies = len(galaxies)
    check(n_galaxies > 0, f"{half}: {n_galaxies} galaxies")
    check(bool(np.isfinite(galaxies["spec_log_evidence"]).all()),
          f"{half}: spec_log_evidence finite for every galaxy")
    for column in CONVERGED_COLUMNS:
        n_converged = int(galaxies[column].eq(True).sum())
        check(n_converged >= MIN_FIT_FRACTION * n_galaxies,
              f"{half}: {column} True for {n_converged} of {n_galaxies}")
    dblplaw_converged = galaxies["spec_dblplaw_converged"].eq(True)
    check(bool(np.isfinite(galaxies.loc[dblplaw_converged, "spec_dblplaw_log_evidence"]).all()),
          f"{half}: spec_dblplaw_log_evidence finite for every converged DPL spectral fit")

if failures:
    sys.exit(f"STOP: {len(failures)} check(s) failed; do not run the locked test")
print("pre-flight passed: run the locked test")
```

  - Expected: every line `ok` and `pre-flight passed: run the locked test`. A fitter verdict of `False` is a result, not a failure; an empty one means that calibration test was incomplete.
  - On `STOP`, do not flip `RUN_LOCKED_TEST`: redo the step of the Task 14 PC order that the failing line names (a missing or empty `spec_dblplaw_*` column: (e) without `spectral_fits_dblplaw`, or (d) 1 unfinished; too few test-half `phot_*`/`phot_dblplaw_*` fits: (b)-(d)), then (e), then this check again. If the DPL spectral fits themselves converge for fewer than 90% of a half, report it to the user before going on.

1. Set `RUN_LOCKED_TEST = True` in cell 1, keeping the JSON layout (from `Synthetic_domain`):

```bash
LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python - <<'EOF'
import json

with open("ml_model.ipynb", encoding="utf-8") as notebook_file:
    notebook = json.load(notebook_file)
cell = next(cell for cell in notebook["cells"] if cell["id"] == "b6229f82")
old_line = "RUN_LOCKED_TEST = False            # set True once, for the final report (both models)\n"
new_line = "RUN_LOCKED_TEST = True             # set True once, for the final report (both models)\n"
assert cell["source"].count(old_line) == 1, "the RUN_LOCKED_TEST line is not the expected one"
cell["source"][cell["source"].index(old_line)] = new_line
with open("ml_model.ipynb", "w", encoding="utf-8") as notebook_file:
    json.dump(notebook, notebook_file, indent=1, ensure_ascii=False)
    notebook_file.write("\n")
print("RUN_LOCKED_TEST = True")
EOF
```

2. Run the GPU notebook job with the same inputs as Task 12 Step 7 (`spectroscopic_reference.csv` now carries the `spec_dblplaw_*` and `spec_*log_evidence*` columns of Task 14):

```bash
GPU=1 CORES=16 MEMORY_MB=7100 TIME_LIMIT=12:00:00 FLOAT32_FILES="training_set.parquet training_set_dblplaw.parquet" \
    HLJOB_OPTIONS=--pause-llm remote/submit_notebook.sh ml_model.ipynb ml_runs/mlp ml_runs/mlp_dblplaw \
    synthetic_validation.parquet synthetic_test.parquet training_set.parquet \
    synthetic_validation_dblplaw.parquet synthetic_test_dblplaw.parquet training_set_dblplaw.parquet \
    spectroscopic_reference.csv photometric_cross_fits.csv dr3_corrected.csv splus_noise_model.csv
```

3. Check in the outputs: "Photometric Bagpipes and evidence comparisons: validation and test halves" (otherwise the test-half photometric fits were not read); the `own_spectral` verdicts of both models for both halves; `against_spectra` rows for both `spectral_model` values; the spectral evidence tables and the `evidence_agreement` rows; `out/memory_peak_bytes.txt` below the 7.1 GB cap. The locked test adds the synthetic test sets of both models (sections 2) to the run time; the 12 h limit covers it.

4. **README notes for Task 16 (repository `README.md`):**
   - Replace the bullet "**The spectral reference is delayed-τ only.** Its `mw_age` and mass fractions carry the delayed-τ SFH shape, so part of a DPL disagreement with the spectra is a model difference, not an ML error." (line 537 at 9b4676e) with the caveat of the new `10a86e1d` markdown: each spectral reference carries its own SFH shape (`mw_age`, `f_young`, `f_mid`, `f_old` follow delayed-τ in the delayed-τ fits and DPL in the DPL fits), and the "spectral Bagpipes" pair of section 5.2 shows how far apart the two references are.
   - Say that, from Amendment 2 on, `agreement_spectral` in `ml_evaluation_summary_dblplaw.csv` is the DPL ML against the DPL spectral fits; the DPL ML against the delayed-τ spectra is in `ml_model_comparison.csv`, section `against_spectra`, `spectral_model` = `delayed`. Before Amendment 2 those rows were the DPL ML against the delayed-τ spectral fits, judged on `log_metallicity` and `log_mstar`.

### Task 16: README and PC runs of Part 2 (controller)

- [ ] **Step 1:** After Tasks 13–15, update `README.md`:
  - module map and §4: DPL spectral fits;
  - §7: the run order of Part 2;
  - §8: the DPL spectral calibration, the comparison against both spectral references, the spectral evidence, and the locked-test results;
  - §9: replace "The spectral reference is delayed-τ only".
- [ ] **Step 2:** On the PC, in order, with fetches and commits between the steps:
  1. The remaining delayed-τ spectral fits (job `351292`).
  2. Notebook 2 exports the DPL spectral jobs and mocks.
  3. `FITS_SUBDIR=spectral_fits_dblplaw remote/submit_fits.sh`, twice if the 24 h limit is reached.
  4. Notebook 2 with `EXPORT_LOCKED_TEST_JOBS = True`: reads everything, makes the grown split, and exports the photometric jobs of both models for the new validation galaxies and the test half.
  5. The photometric fits of both folders.
  6. Notebook 2 reads them; then the Task 15 pre-flight check.
  7. `ml_model.ipynb` with `RUN_LOCKED_TEST = True`, once.
