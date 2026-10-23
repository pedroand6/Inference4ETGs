# Training Set v2 and Spectroscopic Reference Sets Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the data foundation for an ML model that predicts 16-bin SFHs from S-PLUS photometry: a Bagpipes training set with 16-bin labels and a realistic noise model, synthetic validation/test sets, a χ² coverage test, and a real-galaxy reference set whose labels come from Bagpipes fits to SDSS spectra under the training prior.

**Architecture:** One shared module, `Synthetic_domain/synthetic_model.py`, holds the prior and every piece of logic both notebooks use, so the training set and the spectral fits cannot drift apart. `training_set.ipynb` (part 1) and `spectroscopic_reference.ipynb` (part 2) only orchestrate, display and save. Spectral fits run MultiNest with a custom likelihood whose prior is the training unit cube.

**Tech Stack:** Python 3.10 venv at `.venv`, Bagpipes 1.3.5 (`bc03_miles` grids), PyMultiNest 2.12 (library in `.venv/lib`), numpy/scipy/pandas/astropy, joblib, h5py, pyarrow, astroquery, splusdata, dustmaps (CSFD), extinction, pytest, nbconvert.

**Spec:** `docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md`

## Global Constraints

- SPS: Bagpipes `bc03_miles` (installed default). Cosmology: `FlatLambdaCDM(H0=70, Om0=0.3)`.
- `Z_RANGE = (0.002, 0.7)`; `N_MODELS = 2**19`; synthetic validation and test sets `2**14` models each.
- Spectral fits: MultiNest through PyMultiNest; kernel must have `LD_LIBRARY_PATH=/home/pedroand/Projects/Inference4ETGs/.venv/lib`.
- Spectral fitting budget: at most 3 hours wall time on 11 workers.
- SDSS spectra are corrected for Milky Way extinction with CSFD E(B−V), CCM89, R_V = 3.1 (same as `extinctions.py`).
- Spectral fits use the spectrum only; S-PLUS photometry never enters the spectral likelihood.
- `spectra/`, `spectral_fits/`, `*.parquet` are git-ignored; the old tracked `Synthetic_domain/training_set.csv` is removed from the repository.
- **Code style** (applies to every file):
  1. English names that say what they hold; no single letters except loop indices and standard maths (`z`).
  2. One step per block, numbered section banners (`# 1. Load the S-PLUS catalogue`).
  3. Comments explain *why* or the physics, never restate the code.
  4. No compact code: no nested comprehensions, no chained one-liners, no lambdas, no clever indexing without a named intermediate. Prefer explicit loops when they read better.
  5. Every function has a docstring: purpose, inputs with units, return value.
  6. UPPER_CASE constants at the top with units in a comment; no inline magic numbers.
  7. Each notebook opens with a markdown cell (goal, inputs, outputs, run order); each section starts with a short markdown cell.
  8. Notebook cells orchestrate and display; shared or unit-tested logic lives in the module.
  9. Parallel cases follow one pattern (e.g. the three SFH branches share one signature).
  10. Be concise.

**Working directory for every command below:** `/home/pedroand/Projects/Inference4ETGs/Synthetic_domain` unless stated. Python is `../.venv/bin/python`.

## Review Focus

1. **A spectrum with most pixels masked** (sky residuals, `ivar = 0` stretches). Expected: `prepare_spectrum` raises a clear `ValueError`, the notebook logs it to `failures.csv` and the batch continues. Test: Task 6, `test_prepare_spectrum_rejects_mostly_masked_spectrum`.
2. **A galaxy near the top of `Z_RANGE`**, where the oldest age bins do not exist. Expected: labels are padded (`frac = 0`, `logsfr = NaN`, `zmet = NaN`) at the old end and the fractions still sum to 1. Test: Task 3, `test_bin_labels_pad_missing_old_bins`.
3. **A galaxy with no models in its redshift window** (zml outside the prior). Expected: reported uncovered with `chi2_min = NaN`, not a crash. Test: Task 5, `test_coverage_galaxy_without_models_in_window`.
4. **Firefly missing values (−9999).** Expected: `firefly_consistent` is NA, not False, so the galaxy is not counted as a failed check. Test: Task 6, `test_compare_with_firefly_missing_values_are_na`.
5. **Re-running part 2 with more galaxies.** Expected: finished fits are skipped and no existing galaxy changes split. Tests: Task 6, `test_split_by_hash_is_stable`; Task 7, `test_fit_and_save_skips_existing_output`.

---

## File Structure

| File | Responsibility |
|---|---|
| `Synthetic_domain/synthetic_model.py` (create) | Prior, SFH branches, enrichment, labels, magnitudes, model simulation, detection rule, noise model, coverage test, spectrum preparation, spectral likelihood, MultiNest fit, posterior labels, Firefly comparison, split. |
| `Synthetic_domain/pytest.ini` (create) | Lets tests import the module from `Synthetic_domain/`. |
| `Synthetic_domain/tests/make_legacy_fixture.py` (create) | One-off: freezes the original notebook prior as a JSON fixture. |
| `Synthetic_domain/tests/legacy_components.json` (generated) | Regression fixture. |
| `Synthetic_domain/tests/test_synthetic_model.py` (create) | All module tests. |
| `Synthetic_domain/training_set.ipynb` (rewrite) | Part 1. |
| `Synthetic_domain/spectroscopic_reference.ipynb` (create) | Part 2. |
| `.gitignore`, `requirements.txt` (modify) | Ignore generated data; list new dependencies. |

---

### Task 1: Environment

**Files:**
- Modify: `requirements.txt`, `.gitignore`
- Create: `Synthetic_domain/pytest.ini`

**Interfaces:**
- Produces: Jupyter kernel `inference4etgs` with `LD_LIBRARY_PATH` set; `pytest` runnable from `Synthetic_domain/`.

- [ ] **Step 1: Install test and notebook-execution tools**

Run (from repo root): `.venv/bin/pip install pytest nbconvert`
Expected: both install without errors.

- [ ] **Step 2: Register the MultiNest-enabled kernel**

Run (from repo root):
```bash
.venv/bin/python -m ipykernel install --user --name inference4etgs \
  --display-name "Inference4ETGs (MultiNest)" \
  --env LD_LIBRARY_PATH /home/pedroand/Projects/Inference4ETGs/.venv/lib
```
Expected: `Installed kernelspec inference4etgs in ...`.

- [ ] **Step 3: Verify PyMultiNest loads inside that kernel**

Create a scratch notebook check (from repo root):
```bash
cat > /tmp/claude-1000/kernel_check.py <<'EOF'
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
notebook = nbformat.v4.new_notebook()
notebook.cells.append(nbformat.v4.new_code_cell("import pymultinest; print('pymultinest OK')"))
ExecutePreprocessor(timeout=120, kernel_name="inference4etgs").preprocess(notebook, {})
print(notebook.cells[0].outputs[0]["text"])
EOF
.venv/bin/python /tmp/claude-1000/kernel_check.py
```
Expected: `pymultinest OK`.

- [ ] **Step 4: Update `requirements.txt`**

Append these lines:
```
pymultinest
h5py
pyarrow
splusdata
pytest
nbconvert
ipykernel
```

- [ ] **Step 5: Update `.gitignore`**

Append:
```
Synthetic_domain/spectra/
Synthetic_domain/spectral_fits/
Synthetic_domain/*.parquet
```

- [ ] **Step 6: Create `Synthetic_domain/pytest.ini`**

```ini
[pytest]
pythonpath = .
testpaths = tests
```

- [ ] **Step 7: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add requirements.txt .gitignore Synthetic_domain/pytest.ini
git commit -m "chore: test tooling, MultiNest kernel docs and ignores for generated data"
```

---

### Task 2: Shared prior in `synthetic_model.py`

**Files:**
- Create: `Synthetic_domain/tests/make_legacy_fixture.py`, `Synthetic_domain/tests/legacy_components.json` (generated), `Synthetic_domain/synthetic_model.py`, `Synthetic_domain/tests/test_synthetic_model.py`

**Interfaces:**
- Produces:
  - Constants `COSMOLOGY`, `AGE_BIN_CENTRES`, `N_BINS`, `Z_RANGE`, `DIMS`, `LOG_MASS_FORMED`, `MIN_FRACTION` and the other prior constants below.
  - `uniform(draw, low, high) -> float`
  - `positive_normal(draw, mean, sigma) -> float`
  - `bin_edges(redshift: float) -> np.ndarray` (Gyr ago, youngest first)
  - `dirichlet_fractions(params: dict, edges: np.ndarray) -> np.ndarray`
  - `skew_normal_fractions(params: dict, edges: np.ndarray) -> np.ndarray`
  - `skew_normal_burst_fractions(params: dict, edges: np.ndarray) -> np.ndarray`
  - `SFH_BRANCHES: dict[int, function]`
  - `choose_sfh_branch(sfh_type_draw: float) -> int` (0 Dirichlet, 1 skew-normal, 2 skew-normal + burst)
  - `bin_metallicities(fractions: np.ndarray, z_final: float) -> np.ndarray` (Z/Z☉)
  - `components(u: np.ndarray, fixed_redshift: float | None = None) -> tuple[dict, dict]` (Bagpipes components; drawn parameters `z, sfh_type, Z_fin, Av_ism, Av_bc, logU`)
  - `bins_from_components(bagpipes_components: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]` (fractions, metallicities, edges)

- [ ] **Step 1: Write the legacy fixture script**

Create `Synthetic_domain/tests/make_legacy_fixture.py`:
```python
"""Freeze the prior of the original training_set.ipynb (commit c04268e) as a test fixture.

Run once from Synthetic_domain/:  ../.venv/bin/python tests/make_legacy_fixture.py
The refactored module must reproduce these components exactly.
"""

import json
import subprocess

import numpy as np

LEGACY_COMMIT = "c04268e"
NOTEBOOK_PATH = "Synthetic_domain/training_set.ipynb"
CONSTANTS_CELL = 1      # imports and prior constants
PRIOR_CELL = 3          # bin_edges, mass_fractions, components
NEW_Z_RANGE = (0.002, 0.7)   # the only prior change in the new design
N_VECTORS = 30
SEED = 123
OUTPUT_FILE = "tests/legacy_components.json"


# 1. Read the original notebook from git
notebook_text = subprocess.run(
    ["git", "show", f"{LEGACY_COMMIT}:{NOTEBOOK_PATH}"],
    capture_output=True, text=True, check=True,
).stdout
notebook = json.loads(notebook_text)

# 2. Run its constants and prior cells in a private namespace
namespace = {}
for cell_index in (CONSTANTS_CELL, PRIOR_CELL):
    source = "".join(notebook["cells"][cell_index]["source"])
    exec(source, namespace)
namespace["Z_RANGE"] = NEW_Z_RANGE

# 3. Evaluate the legacy prior on fixed unit-cube vectors
rng = np.random.default_rng(SEED)
unit_vectors = rng.random((N_VECTORS, len(namespace["DIMS"])))
cases = []
for u in unit_vectors:
    legacy_components, legacy_params = namespace["components"](u)
    cases.append({"u": u.tolist(), "components": legacy_components, "params": legacy_params})

# 4. Save (numpy floats are Python floats; anything else is cast to float)
with open(OUTPUT_FILE, "w") as output:
    json.dump(cases, output, default=float, indent=1)
print(f"{len(cases)} legacy cases written to {OUTPUT_FILE}")
```

- [ ] **Step 2: Generate the fixture**

Run: `../.venv/bin/python tests/make_legacy_fixture.py`
Expected: `30 legacy cases written to tests/legacy_components.json` (Bagpipes banner lines may print first).

- [ ] **Step 3: Write the failing tests**

Create `Synthetic_domain/tests/test_synthetic_model.py`:
```python
"""Tests for synthetic_model.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest -v"""

import json

import numpy as np
import pandas as pd
import pytest

import synthetic_model as sm

LEGACY_FIXTURE = "tests/legacy_components.json"


# =============================================================================
# Prior
# =============================================================================
def test_components_match_legacy_notebook():
    with open(LEGACY_FIXTURE) as fixture:
        cases = json.load(fixture)

    for case in cases:
        new_components, new_params = sm.components(np.array(case["u"]))
        legacy_components = case["components"]

        assert set(new_components) == set(legacy_components)
        assert np.isclose(new_components["redshift"], legacy_components["redshift"], rtol=1e-12)
        assert np.isclose(new_components["t_bc"], legacy_components["t_bc"], rtol=1e-12)
        for key in ("Av", "eta", "n"):
            assert np.isclose(new_components["dust"][key], legacy_components["dust"][key], rtol=1e-12)
        assert np.isclose(new_components["nebular"]["logU"], legacy_components["nebular"]["logU"], rtol=1e-12)

        for name in legacy_components:
            if not name.startswith("constant"):
                continue
            for key in ("age_min", "age_max", "massformed", "metallicity"):
                assert np.isclose(new_components[name][key], legacy_components[name][key], rtol=1e-12), (name, key)

        for key, legacy_value in case["params"].items():
            assert np.isclose(new_params[key], legacy_value, rtol=1e-12), key


@pytest.mark.parametrize("branch, sfh_type_draw", [(0, 0.25), (1, 0.75), (2, 0.99)])
def test_each_branch_gives_fractions_that_sum_to_one(branch, sfh_type_draw):
    rng = np.random.default_rng(1)
    u = rng.random(len(sm.DIMS))
    u[sm.DIMS.index("sfh_type")] = sfh_type_draw

    assert sm.choose_sfh_branch(sfh_type_draw) == branch
    bagpipes_components, drawn = sm.components(u)
    fractions, metallicities, edges = sm.bins_from_components(bagpipes_components)

    assert drawn["sfh_type"] == branch
    assert np.isclose(fractions.sum(), 1.0, atol=1e-4)    # MIN_FRACTION floor adds at most 16e-6
    assert len(fractions) == len(edges) - 1


def test_old_bins_are_dropped_at_high_redshift():
    low_redshift_edges = sm.bin_edges(0.01)
    high_redshift_edges = sm.bin_edges(0.69)

    # The 14.2 Gyr bin centre is older than the Universe (13.37 Gyr at most here), so the
    # prior never has more than N_BINS - 1 bins; its labels are always padding.
    assert len(low_redshift_edges) - 1 == sm.N_BINS - 1
    assert len(high_redshift_edges) - 1 < sm.N_BINS - 1
    universe_age = sm.COSMOLOGY.age(0.69).value
    assert high_redshift_edges[-1] < universe_age


def test_metallicity_rises_from_oldest_to_youngest_bin():
    fractions = np.full(sm.N_BINS, 1 / sm.N_BINS)
    metallicities = sm.bin_metallicities(fractions, z_final=1.0)

    # Bins run youngest first, so metallicity must fall along the array
    assert np.all(np.diff(metallicities) < 0)
    assert metallicities[0] < 1.0
    assert metallicities[-1] > sm.ZINI_FRACTION * 1.0


def test_fixed_redshift_overrides_the_drawn_redshift():
    u = np.full(len(sm.DIMS), 0.5)
    bagpipes_components, drawn = sm.components(u, fixed_redshift=0.123)

    assert bagpipes_components["redshift"] == 0.123
    assert drawn["z"] == 0.123
```

- [ ] **Step 4: Run tests to verify they fail**

Run: `../.venv/bin/python -m pytest -v`
Expected: collection error `ModuleNotFoundError: No module named 'synthetic_model'`.

- [ ] **Step 5: Create the module with the prior**

Create `Synthetic_domain/synthetic_model.py`:
```python
"""
Shared Bagpipes model of the S-PLUS synthetic domain.

training_set.ipynb and spectroscopic_reference.ipynb both import this module,
so the training set and the spectral fits use exactly the same prior.
"""

import hashlib
import os
import time
import warnings

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simpson
from scipy.stats import chi2 as chi_squared
from scipy.stats import norm, skewnorm

import bagpipes as pipes

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1. Prior (fixed before looking at the data)
# =============================================================================
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3)
UNIVERSE_AGE_MARGIN = 0.1        # [Gyr] oldest bin ends this long after the Big Bang

AGE_BIN_CENTRES = np.array([0.007, 0.021, 0.0425, 0.076, 0.1145, 0.1615, 0.23, 0.345,
                            0.545, 0.865, 1.37, 2.17, 3.44, 5.455, 8.645, 14.2])   # [Gyr ago]
N_BINS = len(AGE_BIN_CENTRES)
LABEL_EDGES = [0.1, 1.0]         # [Gyr] f_young < 0.1 < f_mid < 1.0 < f_old

Z_RANGE = (0.002, 0.7)                               # redshift
LOG_ZFIN_RANGE = (-1.0, 0.4)                         # log(Z_fin / Z_sun): 0.1 to 2.5
ZINI_FRACTION = 0.4                                  # Z_ini = ZINI_FRACTION * Z_fin
AV_ISM_RANGE = (0.01, 2.0)                           # [mag]
LOG_TAU_BC_RANGE = (np.log10(0.5), np.log10(4.0))    # birth-cloud optical depth
LOGU_RANGE = (-4.0, -1.0)                            # nebular ionisation parameter
DUST_SLOPE = 0.7                                     # CF00 attenuation-law exponent
BIRTH_CLOUD_AGE = 0.01                               # [Gyr]
LOG_MASS_FORMED = 10.0                               # [log M_sun]
TAU_TO_AV = 1.086                                    # A_V = 1.086 * tau_V

P_SKEW = 0.5                     # share of skew-normal SFHs (with or without burst)
P_BURST = 6 / 165                # share of skew-normal + burst SFHs
MIN_FRACTION = 1e-6              # keeps every bin as a Bagpipes component
SFR_FLOOR = 1e-14                # [yr^-1] sSFR floor for SFHs with no recent star formation

SKEW_PEAK_MEAN, SKEW_PEAK_SIGMA = 12.0, 3.0          # [Gyr ago]
SKEW_WIDTH_MEAN, SKEW_WIDTH_SIGMA = 1.0, 0.4         # [Gyr]
SKEW_SHAPE_RANGE = (-0.25, 0.75)                     # skew-normal shape parameter
LOG_BURST_AGE_RANGE = (np.log10(0.01), np.log10(0.5))   # [log Gyr ago]
BURST_WIDTH_RANGE = (0.040, 0.050)                   # [Gyr] full width
BURST_FRACTION_RANGE = (0.01, 0.10)                  # share of the mass in the burst

DIRICHLET_DIMS = []
for bin_index in range(N_BINS):
    DIRICHLET_DIMS.append(f"dir{bin_index}")
DIMS = ["z", "sfh_type", *DIRICHLET_DIMS, "peak", "width", "skew",
        "burst_age", "burst_width", "burst_frac", "log_zfin", "av_ism", "log_tau_bc", "logU"]


def uniform(draw, low, high):
    """Map a unit-interval draw onto [low, high)."""
    return low + (high - low) * draw


def positive_normal(draw, mean, sigma):
    """Map a unit-interval draw onto N(mean, sigma) truncated to positive values (inverse CDF)."""
    probability_below_zero = norm.cdf(0, mean, sigma)
    return norm.ppf(probability_below_zero + draw * (1 - probability_below_zero), mean, sigma)


def bin_edges(redshift):
    """Age-bin edges [Gyr ago], youngest first.

    Bins centred beyond the age of the Universe (minus the margin) are dropped,
    so high-redshift galaxies have fewer bins.
    """
    maximum_age = COSMOLOGY.age(redshift).value - UNIVERSE_AGE_MARGIN
    centres = AGE_BIN_CENTRES[AGE_BIN_CENTRES <= maximum_age]
    midpoints = (centres[:-1] + centres[1:]) / 2
    return np.concatenate([[0.0], midpoints, [maximum_age]])


# -----------------------------------------------------------------------------
# SFH branches: same signature, each returns mass fractions per bin (sum 1)
# -----------------------------------------------------------------------------
def dirichlet_fractions(params, edges):
    """Dirichlet(1) fractions: normalised exponential draws, one per bin."""
    n_bins = len(edges) - 1
    exponential_draws = np.zeros(n_bins)
    for bin_index in range(n_bins):
        exponential_draws[bin_index] = -np.log1p(-params[f"dir{bin_index}"])
    return exponential_draws / exponential_draws.sum()


def skew_normal_fractions(params, edges):
    """Fractions from one skew-normal episode in lookback time (typically old)."""
    peak = norm.ppf(params["peak"], SKEW_PEAK_MEAN, SKEW_PEAK_SIGMA)
    width = positive_normal(params["width"], SKEW_WIDTH_MEAN, SKEW_WIDTH_SIGMA)
    shape = uniform(params["skew"], *SKEW_SHAPE_RANGE)

    cumulative_mass = skewnorm.cdf(edges, shape, peak, width)
    fractions = np.diff(cumulative_mass)
    if fractions.sum() > 0:
        return fractions / fractions.sum()

    # The whole episode falls outside the bins: fall back to a flat SFH
    n_bins = len(edges) - 1
    return np.full(n_bins, 1 / n_bins)


def skew_normal_burst_fractions(params, edges):
    """Skew-normal fractions plus a short recent burst holding 1-10% of the mass."""
    base_fractions = skew_normal_fractions(params, edges)
    burst_age = 10 ** uniform(params["burst_age"], *LOG_BURST_AGE_RANGE)
    burst_half_width = uniform(params["burst_width"], *BURST_WIDTH_RANGE) / 2
    burst_fraction = uniform(params["burst_frac"], *BURST_FRACTION_RANGE)

    # Time the burst overlaps each bin, so its mass is shared in proportion
    bin_starts = edges[:-1]
    bin_ends = edges[1:]
    overlap_start = np.maximum(bin_starts, burst_age - burst_half_width)
    overlap_end = np.minimum(bin_ends, burst_age + burst_half_width)
    overlap = np.clip(overlap_end - overlap_start, 0, None)

    fractions = (1 - burst_fraction) * base_fractions + burst_fraction * overlap / overlap.sum()
    return fractions / fractions.sum()


SFH_BRANCHES = {
    0: dirichlet_fractions,
    1: skew_normal_fractions,
    2: skew_normal_burst_fractions,
}


def choose_sfh_branch(sfh_type_draw):
    """SFH branch for a unit draw: 0 Dirichlet, 1 skew-normal, 2 skew-normal + burst."""
    if sfh_type_draw < 1 - P_SKEW:
        return 0
    if sfh_type_draw < 1 - P_BURST:
        return 1
    return 2


def bin_metallicities(fractions, z_final):
    """Z/Z_sun of each bin: rises from Z_ini to Z_final with the mass already formed."""
    z_initial = ZINI_FRACTION * z_final
    # Bins run youngest first, so the mass formed before a bin is the sum of the older ones
    mass_before = np.cumsum(fractions[::-1])[::-1] - fractions
    enrichment = (mass_before + fractions / 2) / fractions.sum()
    return z_initial + (z_final - z_initial) * enrichment


def components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one unit-cube vector u (order of DIMS).

    fixed_redshift: if given, replaces the redshift drawn from u (spectral fits).
    Returns (bagpipes_components, drawn) where drawn holds z, sfh_type, Z_fin, Av_ism, Av_bc, logU.
    """
    params = dict(zip(DIMS, u))

    # 1. Redshift and the age bins that exist at that redshift
    if fixed_redshift is None:
        redshift = uniform(params["z"], *Z_RANGE)
    else:
        redshift = fixed_redshift
    edges = bin_edges(redshift)

    # 2. Star formation history
    sfh_type = choose_sfh_branch(params["sfh_type"])
    fractions = SFH_BRANCHES[sfh_type](params, edges)
    fractions = np.fmax(fractions, MIN_FRACTION)

    # 3. Chemical enrichment
    z_final = 10 ** uniform(params["log_zfin"], *LOG_ZFIN_RANGE)
    metallicities = bin_metallicities(fractions, z_final)

    # 4. Dust and nebular emission
    av_ism = uniform(params["av_ism"], *AV_ISM_RANGE)
    av_birth_cloud = TAU_TO_AV * 10 ** uniform(params["log_tau_bc"], *LOG_TAU_BC_RANGE)
    log_u = uniform(params["logU"], *LOGU_RANGE)

    bagpipes_components = {
        "redshift": redshift,
        "t_bc": BIRTH_CLOUD_AGE,
        "dust": {"type": "CF00", "Av": av_ism, "eta": 1 + av_birth_cloud / av_ism, "n": DUST_SLOPE},
        "nebular": {"logU": log_u},
    }
    # One "constant" component per bin; Bagpipes reads the type from the name minus its last letter
    for bin_index in range(len(fractions)):
        name = "constant" + chr(ord("A") + bin_index)
        bagpipes_components[name] = {
            "age_min": edges[bin_index],
            "age_max": edges[bin_index + 1],
            "massformed": LOG_MASS_FORMED + np.log10(fractions[bin_index]),
            "metallicity": metallicities[bin_index],
        }

    drawn = {"z": redshift, "sfh_type": sfh_type, "Z_fin": z_final,
             "Av_ism": av_ism, "Av_bc": av_birth_cloud, "logU": log_u}
    return bagpipes_components, drawn


def bins_from_components(bagpipes_components):
    """Mass fractions, metallicities and edges [Gyr ago] read back from the constant components."""
    names = []
    for name in bagpipes_components:
        if name.startswith("constant"):
            names.append(name)
    names.sort()

    n_bins = len(names)
    fractions = np.zeros(n_bins)
    metallicities = np.zeros(n_bins)
    edges = np.zeros(n_bins + 1)
    for bin_index, name in enumerate(names):
        component = bagpipes_components[name]
        fractions[bin_index] = 10 ** (component["massformed"] - LOG_MASS_FORMED)
        metallicities[bin_index] = component["metallicity"]
        edges[bin_index] = component["age_min"]
        edges[bin_index + 1] = component["age_max"]
    return fractions, metallicities, edges
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `../.venv/bin/python -m pytest -v`
Expected: 7 passed.

- [ ] **Step 7: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests
git commit -m "feat: shared prior module reproducing the training_set notebook prior"
```

---

### Task 3: Labels, magnitudes and model simulation

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (append section 2)
- Modify: `Synthetic_domain/tests/test_synthetic_model.py` (append)

**Interfaces:**
- Consumes: `components`, `bins_from_components`, `bin_edges`, `bin_metallicities` (Task 2).
- Produces:
  - `BANDS = ["u", "J0378", "J0395", "J0410", "J0430", "g", "J0515", "r", "J0660", "i", "J0861", "z"]`
  - `FILTER_FILES` (absolute paths, same order as `BANDS`), `PIVOT_WAVELENGTHS` (Å)
  - `BIN_LABEL_COLUMNS` (48 names), `SUMMARY_LABEL_COLUMNS`, `MASS_DEPENDENT_LABELS`
  - `bin_labels(fractions, metallicities, edges, log_mass_formed=LOG_MASS_FORMED) -> dict`
  - `summary_labels(bagpipes_sfh) -> dict`
  - `mass_weighted_metallicity(fractions, metallicities) -> float`
  - `ab_magnitudes(flux_per_angstrom: np.ndarray) -> np.ndarray` (NaN if any flux ≤ 0)
  - `simulate_models(unit_vectors: np.ndarray) -> pd.DataFrame`: columns drawn params, summary labels, bin labels, one magnitude column per name in `BANDS`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_synthetic_model.py`:
```python
# =============================================================================
# Labels and magnitudes
# =============================================================================
def test_bin_labels_pad_missing_old_bins():
    u = np.full(len(sm.DIMS), 0.5)
    bagpipes_components, _ = sm.components(u, fixed_redshift=0.69)
    fractions, metallicities, edges = sm.bins_from_components(bagpipes_components)
    labels = sm.bin_labels(fractions, metallicities, edges)

    n_existing = len(fractions)
    assert n_existing < sm.N_BINS
    total = 0.0
    for bin_index in range(sm.N_BINS):
        suffix = f"{bin_index:02d}"
        total += labels[f"frac_bin_{suffix}"]
        if bin_index < n_existing:
            assert np.isfinite(labels[f"logsfr_bin_{suffix}"])
            assert np.isfinite(labels[f"zmet_bin_{suffix}"])
        else:
            assert labels[f"frac_bin_{suffix}"] == 0.0
            assert np.isnan(labels[f"logsfr_bin_{suffix}"])
            assert np.isnan(labels[f"zmet_bin_{suffix}"])
    assert np.isclose(total, 1.0)
    assert list(labels) == sm.BIN_LABEL_COLUMNS


def test_mass_weighted_metallicity_of_equal_bins_is_the_mean():
    fractions = np.array([0.5, 0.5])
    metallicities = np.array([1.0, 2.0])
    assert np.isclose(sm.mass_weighted_metallicity(fractions, metallicities), 1.5)


def test_ab_magnitudes_are_nan_when_any_flux_is_not_positive():
    flux = np.full(len(sm.BANDS), 1e-17)
    assert np.all(np.isfinite(sm.ab_magnitudes(flux)))
    flux[3] = 0.0
    assert np.all(np.isnan(sm.ab_magnitudes(flux)))


def test_simulate_models_returns_labels_and_magnitudes():
    unit_vectors = np.random.default_rng(2).random((4, len(sm.DIMS)))
    models = sm.simulate_models(unit_vectors)

    assert len(models) == 4
    for column in ["z", "sfh_type", *sm.SUMMARY_LABEL_COLUMNS, *sm.BIN_LABEL_COLUMNS, *sm.BANDS]:
        assert column in models.columns, column
    fraction_columns = []
    for bin_index in range(sm.N_BINS):
        fraction_columns.append(f"frac_bin_{bin_index:02d}")
    assert np.allclose(models[fraction_columns].sum(axis=1), 1.0)
    assert np.isfinite(models[sm.BANDS].to_numpy()).all()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `../.venv/bin/python -m pytest -v -k "labels or magnitudes or simulate or mass_weighted"`
Expected: FAIL with `AttributeError: module 'synthetic_model' has no attribute 'bin_labels'` (and similar).

- [ ] **Step 3: Append section 2 to the module**

Append to `synthetic_model.py`:
```python
# =============================================================================
# 2. Labels, S-PLUS magnitudes and model simulation
# =============================================================================
MODULE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
FILTER_LIST_FILE = os.path.join(MODULE_DIRECTORY, "splus-filters", "splus_filters.txt")
SPEED_OF_LIGHT = 2.99792458e18   # [Angstrom / s]
AB_ZERO_POINT = 48.60

# S-PLUS bands, in the order of splus_filters.txt
BANDS = ["u", "J0378", "J0395", "J0410", "J0430", "g", "J0515", "r", "J0660", "i", "J0861", "z"]

FILTER_FILES = []
for relative_path in np.loadtxt(FILTER_LIST_FILE, dtype="str"):
    FILTER_FILES.append(os.path.join(MODULE_DIRECTORY, relative_path))


def pivot_wavelength(filter_file):
    """Pivot wavelength [Angstrom] of a filter curve file (wavelength, transmission)."""
    wavelength, transmission = np.loadtxt(filter_file, unpack=True)
    numerator = simpson(transmission * wavelength, x=wavelength)
    denominator = simpson(transmission / wavelength, x=wavelength)
    return np.sqrt(numerator / denominator)


PIVOT_WAVELENGTHS = np.zeros(len(FILTER_FILES))
for band_index, filter_file in enumerate(FILTER_FILES):
    PIVOT_WAVELENGTHS[band_index] = pivot_wavelength(filter_file)

SUMMARY_LABEL_COLUMNS = ["log_mstar", "log_sfr", "mw_age", "f_young", "f_mid", "f_old"]

BIN_LABEL_COLUMNS = []
for label_kind in ("frac_bin", "logsfr_bin", "zmet_bin"):
    for bin_index in range(N_BINS):
        BIN_LABEL_COLUMNS.append(f"{label_kind}_{bin_index:02d}")

# Labels that scale with brightness: log value - 0.4 * (m_obs - m_model)
MASS_DEPENDENT_LABELS = ["log_mstar", "log_sfr"]
for bin_index in range(N_BINS):
    MASS_DEPENDENT_LABELS.append(f"logsfr_bin_{bin_index:02d}")


def bin_labels(fractions, metallicities, edges, log_mass_formed=LOG_MASS_FORMED):
    """16-bin labels, youngest bin first.

    fractions: mass formed per existing bin (any normalisation); metallicities: Z/Z_sun;
    edges: [Gyr ago]. Bins that do not exist at this redshift get frac 0, logsfr and zmet NaN.
    Returns a dict ordered as BIN_LABEL_COLUMNS.
    """
    n_existing = len(fractions)
    normalised_fractions = fractions / fractions.sum()
    widths_in_years = np.diff(edges) * 1e9

    labels = {}
    for bin_index in range(N_BINS):
        if bin_index < n_existing:
            labels[f"frac_bin_{bin_index:02d}"] = normalised_fractions[bin_index]
        else:
            labels[f"frac_bin_{bin_index:02d}"] = 0.0

    for bin_index in range(N_BINS):
        if bin_index < n_existing:
            mass_in_bin = 10**log_mass_formed * normalised_fractions[bin_index]
            labels[f"logsfr_bin_{bin_index:02d}"] = np.log10(mass_in_bin / widths_in_years[bin_index])
        else:
            labels[f"logsfr_bin_{bin_index:02d}"] = np.nan

    for bin_index in range(N_BINS):
        if bin_index < n_existing:
            labels[f"zmet_bin_{bin_index:02d}"] = metallicities[bin_index]
        else:
            labels[f"zmet_bin_{bin_index:02d}"] = np.nan
    return labels


def summary_labels(bagpipes_sfh):
    """log_mstar (living), log_sfr (last 100 Myr), mw_age [Gyr] and f_young/f_mid/f_old."""
    mass_per_age_step = bagpipes_sfh.sfh * bagpipes_sfh.age_widths
    total_mass = mass_per_age_step.sum()
    fraction_younger_than = []
    for edge in LABEL_EDGES:
        younger = bagpipes_sfh.ages < edge * 1e9
        fraction_younger_than.append(mass_per_age_step[younger].sum() / total_mass)

    sfr_floor = SFR_FLOOR * 10**bagpipes_sfh.stellar_mass
    log_sfr = np.log10(max(bagpipes_sfh.sfr, sfr_floor))
    return {
        "log_mstar": bagpipes_sfh.stellar_mass,
        "log_sfr": log_sfr,
        "mw_age": bagpipes_sfh.mass_weighted_age,
        "f_young": fraction_younger_than[0],
        "f_mid": fraction_younger_than[1] - fraction_younger_than[0],
        "f_old": 1 - fraction_younger_than[1],
    }


def mass_weighted_metallicity(fractions, metallicities):
    """Mean Z/Z_sun weighted by the mass formed in each bin."""
    return np.sum(fractions * metallicities) / np.sum(fractions)


def ab_magnitudes(flux_per_angstrom):
    """AB magnitudes from Bagpipes photometry [erg/s/cm2/A]; all NaN if any flux is not positive."""
    if not np.all(flux_per_angstrom > 0):
        return np.full(len(PIVOT_WAVELENGTHS), np.nan)
    flux_per_hertz = flux_per_angstrom * PIVOT_WAVELENGTHS**2 / SPEED_OF_LIGHT
    return -2.5 * np.log10(flux_per_hertz) - AB_ZERO_POINT


def simulate_models(unit_vectors):
    """Drawn parameters, labels and AB magnitudes (10^10 M_sun formed) for a block of unit-cube vectors.

    One Bagpipes model is kept per number of bins: update() is ~70x faster than building a new one.
    """
    models_by_bin_count = {}
    rows = []
    for u in unit_vectors:
        bagpipes_components, drawn = components(u)
        fractions, metallicities, edges = bins_from_components(bagpipes_components)

        n_bins = len(fractions)
        if n_bins in models_by_bin_count:
            model = models_by_bin_count[n_bins]
            model.update(bagpipes_components)
        else:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
            models_by_bin_count[n_bins] = model

        row = dict(drawn)
        row.update(summary_labels(model.sfh))
        row.update(bin_labels(fractions, metallicities, edges))
        magnitudes = ab_magnitudes(model.photometry)
        for band, magnitude in zip(BANDS, magnitudes):
            row[band] = magnitude
        rows.append(row)
    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `../.venv/bin/python -m pytest -v`
Expected: 11 passed.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: 16-bin SFR and metallicity labels and model simulation"
```

---

### Task 4: Detection rule and S-PLUS noise model

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (append section 3)
- Modify: `Synthetic_domain/tests/test_synthetic_model.py` (append)

**Interfaces:**
- Consumes: `BANDS` (Task 3).
- Produces:
  - Constants `UNDETECTED_MAGNITUDE = 50`, `MAX_DETECTED_ERROR = 1.0`, `ERROR_FLOOR = 0.05`, `MIN_DETECTED_BANDS = 6`, `Z_SIGMA`, `DZ_MIN`, `DZ_MAX`, `SIGMA_Z`
  - `detected_bands(magnitudes, errors) -> np.ndarray[bool]` (n, 12)
  - `redshift_window(photometric_z) -> np.ndarray` (half-width dz)
  - `fit_noise_model(catalogue: pd.DataFrame) -> pd.DataFrame` (columns `band, magnitude_low, magnitude_high, n_galaxies, log_error_q00 … log_error_q20`; catalogue has `<band>_auto`, `e_<band>_auto`)
  - `add_splus_noise(magnitudes, noise_table, rng) -> (noisy_magnitudes, errors, detected)`; undetected entries are NaN in the first two
  - `scale_to_observed_r(model_magnitudes, model_z, catalogue_r, catalogue_z, rng, max_dz=0.05) -> (scaled_magnitudes, offsets)`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_synthetic_model.py`:
```python
# =============================================================================
# Detection rule and noise model
# =============================================================================
DR3_FILE = "dr3_corrected.csv"


@pytest.fixture(scope="module")
def dr3_galaxies():
    catalogue = pd.read_csv(DR3_FILE)
    return catalogue[catalogue["PROB_GAL"] > 0.8].reset_index(drop=True)


def test_detected_bands_rule():
    magnitudes = np.array([[20.0, 99.0, 20.0, np.nan]])
    errors = np.array([[0.1, 0.1, 1.5, 0.1]])
    assert sm.detected_bands(magnitudes, errors).tolist() == [[True, False, False, False]]


def test_noise_model_reproduces_dr3_detection_rate_and_median_error(dr3_galaxies):
    noise_table = sm.fit_noise_model(dr3_galaxies)
    rng = np.random.default_rng(3)

    # Use DR3's own measured magnitudes as "true" magnitudes
    true_magnitudes = np.full((len(dr3_galaxies), len(sm.BANDS)), np.nan)
    observed_errors = np.full((len(dr3_galaxies), len(sm.BANDS)), np.nan)
    for band_index, band in enumerate(sm.BANDS):
        magnitudes = dr3_galaxies[f"{band}_auto"].to_numpy()
        errors = dr3_galaxies[f"e_{band}_auto"].to_numpy()
        measured = (magnitudes > 0) & (magnitudes < sm.UNDETECTED_MAGNITUDE) & (errors > 0)
        true_magnitudes[measured, band_index] = magnitudes[measured]
        observed_errors[measured, band_index] = errors[measured]

    _, simulated_errors, simulated_detected = sm.add_splus_noise(true_magnitudes, noise_table, rng)

    for band_index, band in enumerate(sm.BANDS):
        measured = np.isfinite(true_magnitudes[:, band_index])
        observed_detected = observed_errors[measured, band_index] < sm.MAX_DETECTED_ERROR
        observed_rate = observed_detected.mean()
        simulated_rate = simulated_detected[measured, band_index].mean()
        assert abs(simulated_rate - observed_rate) < 0.05, band

        observed_median = np.median(observed_errors[measured, band_index][observed_detected])
        simulated_median = np.nanmedian(simulated_errors[measured & simulated_detected[:, band_index], band_index])
        assert abs(simulated_median / observed_median - 1) < 0.05, band


def test_noise_model_marks_very_faint_models_undetected(dr3_galaxies):
    noise_table = sm.fit_noise_model(dr3_galaxies)
    very_faint = np.full((10, len(sm.BANDS)), 35.0)
    noisy, errors, detected = sm.add_splus_noise(very_faint, noise_table, np.random.default_rng(4))
    assert not detected.any()
    assert np.isnan(noisy).all()


def test_scale_to_observed_r_takes_r_from_galaxies_at_similar_redshift():
    model_magnitudes = np.full((5, len(sm.BANDS)), 25.0)
    model_z = np.full(5, 0.3)
    catalogue_r = np.array([18.0, 19.0, 30.0])
    catalogue_z = np.array([0.31, 0.29, 0.9])     # the third galaxy is outside the window
    scaled, offsets = sm.scale_to_observed_r(model_magnitudes, model_z, catalogue_r, catalogue_z,
                                             np.random.default_rng(5))
    r_index = sm.BANDS.index("r")
    assert set(scaled[:, r_index]).issubset({18.0, 19.0})
    assert np.allclose(scaled[:, 0] - model_magnitudes[:, 0], offsets)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `../.venv/bin/python -m pytest -v -k "detected or noise or scale"`
Expected: FAIL with `AttributeError` for `detected_bands`, `fit_noise_model`, `scale_to_observed_r`.

- [ ] **Step 3: Append section 3 to the module**

Append to `synthetic_model.py`:
```python
# =============================================================================
# 3. Detection rule and S-PLUS noise model
# =============================================================================
UNDETECTED_MAGNITUDE = 50.0      # S-PLUS writes 99 for non-detections
MAX_DETECTED_ERROR = 1.0         # [mag] error >= 1 mag (S/N < 1) counts as a non-detection
ERROR_FLOOR = 0.05               # [mag] added in quadrature; stands in for model error
MIN_DETECTED_BANDS = 6

Z_SIGMA = 2.0                    # redshift window: |z_model - zml| < Z_SIGMA * zml_error,
DZ_MIN, DZ_MAX = 0.02, 0.10      # limited to [DZ_MIN, DZ_MAX]
SIGMA_Z = 0.025                  # assumed zml_error = SIGMA_Z * (1 + zml)

NOISE_MAGNITUDE_RANGE = (10.0, 30.0)   # [mag] range of the noise-model bins
NOISE_BIN_WIDTH = 0.25                 # [mag]
NOISE_MIN_GALAXIES = 20                # bins with fewer galaxies are not used
NOISE_N_QUANTILES = 21                 # quantiles of log10(error) stored per bin
MAGNITUDE_TO_RELATIVE_FLUX = 2.5 / np.log(10)   # sigma_mag = 1.0857 * sigma_flux / flux
SCALING_MAX_DZ = 0.05                  # redshift window for borrowing a real r magnitude


def detected_bands(magnitudes, errors):
    """True where a band is detected: magnitude measured and 0 < error < MAX_DETECTED_ERROR (NaN -> False)."""
    measured = magnitudes < UNDETECTED_MAGNITUDE
    good_error = (errors > 0) & (errors < MAX_DETECTED_ERROR)
    return measured & good_error


def redshift_window(photometric_z):
    """Half-width of the model redshift window around each photometric redshift."""
    assumed_error = SIGMA_Z * (1 + photometric_z)
    return np.clip(Z_SIGMA * assumed_error, DZ_MIN, DZ_MAX)


def noise_bin_index(magnitudes):
    """Index of the noise-model bin of each magnitude (may fall outside the table)."""
    return np.floor((magnitudes - NOISE_MAGNITUDE_RANGE[0]) / NOISE_BIN_WIDTH).astype(int)


def quantile_column_names():
    """Names of the stored log-error quantile columns: log_error_q00 ... log_error_q20."""
    names = []
    for quantile_index in range(NOISE_N_QUANTILES):
        names.append(f"log_error_q{quantile_index:02d}")
    return names


def fit_noise_model(catalogue):
    """Empirical error distribution per band and 0.25 mag bin of measured magnitude.

    catalogue: S-PLUS table with <band>_auto and e_<band>_auto columns.
    Returns one row per populated (band, bin) with 21 quantiles of log10(error).
    A band is later detected when its drawn error is below MAX_DETECTED_ERROR, so each bin's
    detection probability is the share of its errors below 1 mag.
    """
    quantile_levels = np.linspace(0, 1, NOISE_N_QUANTILES)
    n_bins = int(round((NOISE_MAGNITUDE_RANGE[1] - NOISE_MAGNITUDE_RANGE[0]) / NOISE_BIN_WIDTH))
    rows = []
    for band in BANDS:
        magnitudes = catalogue[f"{band}_auto"].to_numpy()
        errors = catalogue[f"e_{band}_auto"].to_numpy()
        measured = (magnitudes > 0) & (magnitudes < UNDETECTED_MAGNITUDE) & (errors > 0)
        bin_of_galaxy = noise_bin_index(magnitudes[measured])
        log_errors = np.log10(errors[measured])

        for bin_index in range(n_bins):
            in_bin = bin_of_galaxy == bin_index
            if in_bin.sum() < NOISE_MIN_GALAXIES:
                continue
            row = {
                "band": band,
                "magnitude_low": NOISE_MAGNITUDE_RANGE[0] + bin_index * NOISE_BIN_WIDTH,
                "magnitude_high": NOISE_MAGNITUDE_RANGE[0] + (bin_index + 1) * NOISE_BIN_WIDTH,
                "n_galaxies": int(in_bin.sum()),
            }
            quantiles = np.quantile(log_errors[in_bin], quantile_levels)
            for name, value in zip(quantile_column_names(), quantiles):
                row[name] = value
            rows.append(row)
    return pd.DataFrame(rows)


def draw_log_errors(magnitudes, band_table, rng):
    """Draw log10(error) for each magnitude from its bin's quantiles (inverse-CDF sampling).

    Magnitudes brighter than the first populated bin use that bin; fainter than the last one,
    or NaN, get NaN (undetected).
    """
    table_bins = noise_bin_index(band_table["magnitude_low"].to_numpy() + NOISE_BIN_WIDTH / 2)
    quantile_table = band_table[quantile_column_names()].to_numpy()

    log_errors = np.full(len(magnitudes), np.nan)
    finite = np.isfinite(magnitudes)
    model_bins = np.full(len(magnitudes), -1)
    model_bins[finite] = noise_bin_index(magnitudes[finite])

    for model_index in np.flatnonzero(finite):
        if model_bins[model_index] > table_bins[-1]:
            continue    # fainter than any S-PLUS galaxy with a measurement
        # Nearest populated bin at or below this one; the brightest bin covers anything brighter
        row_index = np.searchsorted(table_bins, model_bins[model_index], side="right") - 1
        row_index = max(row_index, 0)
        quantile_position = rng.random() * (NOISE_N_QUANTILES - 1)
        lower = int(np.floor(quantile_position))
        upper = min(lower + 1, NOISE_N_QUANTILES - 1)
        weight = quantile_position - lower
        row = quantile_table[row_index]
        log_errors[model_index] = row[lower] + weight * (row[upper] - row[lower])
    return log_errors


def add_splus_noise(magnitudes, noise_table, rng):
    """S-PLUS-like errors, flux noise and non-detections for model magnitudes (n_models, 12).

    Returns (noisy_magnitudes, errors, detected); undetected entries are NaN in the first two.
    """
    n_models = len(magnitudes)
    noisy_magnitudes = np.full(magnitudes.shape, np.nan)
    errors = np.full(magnitudes.shape, np.nan)
    detected = np.zeros(magnitudes.shape, dtype=bool)

    for band_index, band in enumerate(BANDS):
        band_table = noise_table[noise_table["band"] == band]
        true_magnitudes = magnitudes[:, band_index]
        band_errors = 10 ** draw_log_errors(true_magnitudes, band_table, rng)

        # Gaussian noise in flux, converted back to magnitudes
        relative_flux_error = band_errors / MAGNITUDE_TO_RELATIVE_FLUX
        flux_ratio = 1 + rng.normal(0, 1, n_models) * relative_flux_error
        positive_flux = flux_ratio > 0
        band_detected = (band_errors < MAX_DETECTED_ERROR) & positive_flux

        noisy_band = np.full(n_models, np.nan)
        noisy_band[band_detected] = true_magnitudes[band_detected] - 2.5 * np.log10(flux_ratio[band_detected])
        noisy_magnitudes[:, band_index] = noisy_band
        errors[band_detected, band_index] = band_errors[band_detected]
        detected[:, band_index] = band_detected
    return noisy_magnitudes, errors, detected


def scale_to_observed_r(model_magnitudes, model_z, catalogue_r, catalogue_z, rng, max_dz=SCALING_MAX_DZ):
    """Shift each model to the r magnitude of a random real galaxy with |zml - z| < max_dz.

    Gives synthetic galaxies a realistic brightness (and so realistic errors) at each redshift.
    Returns (scaled_magnitudes, offsets) with offsets = m_scaled - m_model; NaN if no galaxy is close.
    """
    r_index = BANDS.index("r")
    usable = np.isfinite(catalogue_r) & (catalogue_r < UNDETECTED_MAGNITUDE) & np.isfinite(catalogue_z)
    order = np.argsort(catalogue_z[usable])
    sorted_z = catalogue_z[usable][order]
    sorted_r = catalogue_r[usable][order]

    offsets = np.full(len(model_magnitudes), np.nan)
    for model_index in range(len(model_magnitudes)):
        first = np.searchsorted(sorted_z, model_z[model_index] - max_dz)
        last = np.searchsorted(sorted_z, model_z[model_index] + max_dz)
        if last == first:
            continue
        chosen = rng.integers(first, last)
        offsets[model_index] = sorted_r[chosen] - model_magnitudes[model_index, r_index]

    scaled_magnitudes = model_magnitudes + offsets[:, np.newaxis]
    return scaled_magnitudes, offsets
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `../.venv/bin/python -m pytest -v`
Expected: 15 passed. If `test_noise_model_reproduces_dr3...` fails for one band by a small margin, print that band's observed vs simulated values and check `draw_log_errors` row lookup before changing any tolerance.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: S-PLUS detection rule and empirical noise model"
```

---

### Task 5: Coverage test v2

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (append section 4)
- Modify: `Synthetic_domain/tests/test_synthetic_model.py` (append)

**Interfaces:**
- Consumes: `ERROR_FLOOR` (Task 4).
- Produces: `P_MIN = 0.01`; `coverage_test(observed_magnitudes, observed_errors, detected, observed_z, observed_dz, model_magnitudes, model_z) -> pd.DataFrame` with columns `chi2_min, n_dof, p_value, covered, best_model` (row per galaxy; `best_model` indexes `model_magnitudes`, −1 if none).

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_synthetic_model.py`:
```python
# =============================================================================
# Coverage test
# =============================================================================
@pytest.fixture(scope="module")
def random_model_grid():
    rng = np.random.default_rng(6)
    model_magnitudes = rng.uniform(20, 23, size=(2000, len(sm.BANDS)))
    model_z = rng.uniform(0.0, 0.5, size=2000)
    return model_magnitudes, model_z


def galaxy_from_model(model_magnitudes, model_z, model_index, rng):
    """A noisy copy of one model, 2 mag fainter, with 0.05 mag errors in all bands."""
    errors = np.full((1, len(sm.BANDS)), 0.05)
    magnitudes = model_magnitudes[[model_index]] + 2.0 + rng.normal(0, 0.05, size=(1, len(sm.BANDS)))
    detected = np.ones((1, len(sm.BANDS)), dtype=bool)
    redshift = model_z[[model_index]]
    return magnitudes, errors, detected, redshift


def test_coverage_noisy_model_is_covered(random_model_grid):
    model_magnitudes, model_z = random_model_grid
    magnitudes, errors, detected, redshift = galaxy_from_model(model_magnitudes, model_z, 7, np.random.default_rng(7))

    result = sm.coverage_test(magnitudes, errors, detected, redshift, np.array([0.05]), model_magnitudes, model_z)

    assert bool(result.loc[0, "covered"])
    assert result.loc[0, "best_model"] == 7
    assert result.loc[0, "n_dof"] == len(sm.BANDS) - 1


def test_coverage_model_with_shifted_colour_is_not_covered(random_model_grid):
    model_magnitudes, model_z = random_model_grid
    magnitudes, errors, detected, redshift = galaxy_from_model(model_magnitudes, model_z, 7, np.random.default_rng(7))
    magnitudes[0, sm.BANDS.index("g")] += 0.5

    result = sm.coverage_test(magnitudes, errors, detected, redshift, np.array([0.05]), model_magnitudes, model_z)

    assert not bool(result.loc[0, "covered"])


def test_coverage_galaxy_without_models_in_window(random_model_grid):
    model_magnitudes, model_z = random_model_grid
    magnitudes, errors, detected, _ = galaxy_from_model(model_magnitudes, model_z, 7, np.random.default_rng(7))

    result = sm.coverage_test(magnitudes, errors, detected, np.array([0.9]), np.array([0.05]),
                              model_magnitudes, model_z)

    assert not bool(result.loc[0, "covered"])
    assert np.isnan(result.loc[0, "chi2_min"])
    assert result.loc[0, "best_model"] == -1
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `../.venv/bin/python -m pytest -v -k coverage`
Expected: FAIL with `AttributeError: module 'synthetic_model' has no attribute 'coverage_test'`.

- [ ] **Step 3: Append section 4 to the module**

Append to `synthetic_model.py`:
```python
# =============================================================================
# 4. Coverage test v2
# =============================================================================
P_MIN = 0.01    # a galaxy is covered if its best chi^2 has p-value above this


def coverage_test(observed_magnitudes, observed_errors, detected, observed_z, observed_dz,
                  model_magnitudes, model_z):
    """Is each galaxy consistent with some model at its redshift?

    For every model with |z_model - z| < dz: chi^2 over the detected bands, with the magnitude
    offset (the model's mass normalisation) solved analytically. Covered if the minimum chi^2
    has p-value > P_MIN for n_detected - 1 degrees of freedom.
    Arrays: magnitudes and errors (n, 12) [mag]; detected (n, 12) bool; z and dz (n,).
    Returns a DataFrame with chi2_min, n_dof, p_value, covered, best_model (index into
    model_magnitudes, -1 when no model is in the window).
    """
    # 1. Models with finite magnitudes, sorted by redshift for fast windows
    finite_models = np.flatnonzero(np.isfinite(model_magnitudes).all(axis=1))
    order = finite_models[np.argsort(model_z[finite_models])]
    sorted_z = model_z[order]
    sorted_magnitudes = model_magnitudes[order]

    rows = []
    for galaxy_index in range(len(observed_magnitudes)):
        bands = np.flatnonzero(detected[galaxy_index])
        n_dof = len(bands) - 1
        first = np.searchsorted(sorted_z, observed_z[galaxy_index] - observed_dz[galaxy_index])
        last = np.searchsorted(sorted_z, observed_z[galaxy_index] + observed_dz[galaxy_index])
        if last == first or n_dof < 1:
            rows.append({"chi2_min": np.nan, "n_dof": n_dof, "p_value": np.nan,
                         "covered": False, "best_model": -1})
            continue

        # 2. chi^2 with free offset: sum w d^2 - (sum w d)^2 / sum w, where d = m_obs - m_model
        weights = 1 / (observed_errors[galaxy_index, bands] ** 2 + ERROR_FLOOR**2)
        differences = observed_magnitudes[galaxy_index, bands] - sorted_magnitudes[first:last][:, bands]
        weighted_sum = differences @ weights
        best_offsets = weighted_sum / weights.sum()
        chi2_values = (differences**2) @ weights - weighted_sum * best_offsets

        # 3. Best model and its p-value
        best_in_window = np.argmin(chi2_values)
        chi2_min = chi2_values[best_in_window]
        p_value = chi_squared.sf(chi2_min, n_dof)
        rows.append({"chi2_min": chi2_min, "n_dof": n_dof, "p_value": p_value,
                     "covered": p_value > P_MIN, "best_model": int(order[first + best_in_window])})
    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `../.venv/bin/python -m pytest -v`
Expected: 18 passed.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: chi-square coverage test with redshift window"
```

---

### Task 6: Rewrite `training_set.ipynb` (part 1) and run it

**Files:**
- Rewrite: `Synthetic_domain/training_set.ipynb`
- Remove from git: `Synthetic_domain/training_set.csv`
- Generated (ignored): `training_set.parquet`, `synthetic_validation.parquet`, `synthetic_test.parquet`
- Generated (committed): `splus_noise_model.csv`, `color_color.png`, `coverage_criterion.png`

**Interfaces:**
- Consumes: everything from Tasks 2–5.
- Produces: `training_set.parquet` (columns: drawn params, `SUMMARY_LABEL_COLUMNS`, `BIN_LABEL_COLUMNS`, `BANDS`), `synthetic_validation.parquet` and `synthetic_test.parquet` (same plus `e_<band>` and `magnitude_offset`; undetected bands NaN), `splus_noise_model.csv`. Task 9 reads `training_set.parquet`.

Replace all cells of the notebook with the cells below, in order (use the NotebookEdit tool or `nbformat`). Set the notebook kernel to `inference4etgs`.

- [ ] **Step 1: Cell 0 (markdown)**

```markdown
# Synthetic S-PLUS training set (part 1)

Bagpipes (BC03-MILES) models in the 12 S-PLUS bands, drawn by Latin Hypercube from wide priors fixed **before** looking at the data. All prior code lives in `synthetic_model.py`, shared with `spectroscopic_reference.ipynb`.

**Inputs:** `splus-filters/`, `dr3_corrected.csv` (S-PLUS DR3 corrected for Galactic extinction).
**Outputs:** `training_set.parquet` (noise-free), `synthetic_validation.parquet` and `synthetic_test.parquet` (S-PLUS noise and non-detections), `splus_noise_model.csv`, `color_color.png`, `coverage_criterion.png`.
**Run:** top to bottom with the *Inference4ETGs (MultiNest)* kernel, ~35 min on 11 workers.

**SFH:** 16 constant age bins up to the age of the Universe at $z$: ~50% Dirichlet(1) fractions; ~46% skew-normal (peak $\sim N(12, 3)$ Gyr ago, width $\sim N(1, 0.4)$, shape $U(-0.25, 0.75)$); ~4% skew-normal + burst (10–500 Myr ago, 40–50 Myr wide, 1–10% of the mass). Metallicity rises from $Z_{\rm ini} = 0.4 Z_{\rm fin}$ (oldest stars) to $Z_{\rm fin}$ with the mass formed, so the 16 `zmet_bin_*` labels carry a single free number, `Z_fin`.

| parameter | prior |
|---|---|
| $z$ | $U(0.002, 0.7)$ |
| $Z_{\rm fin}$ | log-uniform, $0.1$–$2.5\,Z_\odot$ |
| ISM $A_V$ | $U(0.01, 2)$ mag |
| extra birth-cloud $A_V$ (< 10 Myr) | $1.086\,\tau_{\rm BC}$, $\tau_{\rm BC}$ log-uniform 0.5–4 |
| $\log U$ | $U(-4, -1)$ |

**Labels:** `frac_bin_00…15` (mass fraction per bin, youngest first), `logsfr_bin_00…15` (log mean SFR per bin for $10^{10}\,M_\odot$ formed), `zmet_bin_00…15` ($Z/Z_\odot$), plus `log_mstar`, `log_sfr`, `mw_age`, `f_young`, `f_mid`, `f_old`. Bins older than the Universe have `frac = 0` and NaN SFR and metallicity. For an observed galaxy, $\Delta m = m_{\rm obs} - m_{\rm model}$ rescales every mass-like label by $-0.4\,\Delta m$ dex.
```

- [ ] **Step 2: Cell 1 (code): imports and constants**

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from scipy.stats import qmc

import synthetic_model as sm

# Sizes and seeds
N_MODELS = 2**19                 # training set
N_HELD_OUT = 2**14               # each synthetic validation / test set
TRAINING_SEED = 0
VALIDATION_SEED = 1
TEST_SEED = 2
N_JOBS = 11
BLOCKS_PER_JOB = 20              # small blocks keep the workers busy until the end

# Files
SPLUS_FILE = "dr3_corrected.csv"
TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
TEST_FILE = "synthetic_test.parquet"
NOISE_MODEL_FILE = "splus_noise_model.csv"

MIN_PROB_GALAXY = 0.8
OLD_GEOMETRIC_COVERAGE = 0.991   # geometric test of commit c04268e (K = 2), kept for one comparison
```

- [ ] **Step 3: Cell 2 (markdown) and Cell 3 (code): training set**

Cell 2:
```markdown
## 1. Training set
Each model is a Latin Hypercube vector $u \in [0, 1)^{28}$; `sm.components` maps it to Bagpipes components by inverse CDFs, so the SFH shares and priors are exactly those of the table.
```

Cell 3:
```python
def simulate_in_parallel(n_models, seed):
    """Latin Hypercube draw of n_models unit-cube vectors, simulated on N_JOBS workers."""
    sampler = qmc.LatinHypercube(len(sm.DIMS), seed=seed)
    unit_vectors = sampler.random(n_models)
    blocks = np.array_split(unit_vectors, BLOCKS_PER_JOB * N_JOBS)
    tables = Parallel(n_jobs=N_JOBS, verbose=5)(delayed(sm.simulate_models)(block) for block in blocks)
    return pd.concat(tables, ignore_index=True)


training_models = simulate_in_parallel(N_MODELS, TRAINING_SEED)
training_models.to_parquet(TRAINING_FILE, index=False)

# Check: Bagpipes returned positive fluxes in every band for every model
finite_magnitudes = np.isfinite(training_models[sm.BANDS].to_numpy()).all(axis=1)
print(f"{len(training_models)} models saved to {TRAINING_FILE}")
assert finite_magnitudes.all(), f"{(~finite_magnitudes).sum()} models have non-finite magnitudes"
display(training_models[["z", "sfh_type", "Z_fin", "log_mstar", "log_sfr", "mw_age", "f_old"]].describe().T.round(3))
```

- [ ] **Step 4: Cell 4 (markdown) and Cell 5 (code): S-PLUS catalogue**

Cell 4:
```markdown
## 2. S-PLUS catalogue
A band is detected when its magnitude is measured and its error is below 1 mag. Galaxies need at least 6 detected bands and a photometric redshift.
```

Cell 5:
```python
galaxy_catalogue = pd.read_csv(SPLUS_FILE)
galaxy_catalogue = galaxy_catalogue[galaxy_catalogue["PROB_GAL"] > MIN_PROB_GALAXY].reset_index(drop=True)

magnitude_columns = []
error_columns = []
for band in sm.BANDS:
    magnitude_columns.append(f"{band}_auto")
    error_columns.append(f"e_{band}_auto")

all_magnitudes = galaxy_catalogue[magnitude_columns].to_numpy()
all_errors = galaxy_catalogue[error_columns].to_numpy()
all_detected = sm.detected_bands(all_magnitudes, all_errors)

# Galaxies for the coverage test
enough_bands = all_detected.sum(axis=1) >= sm.MIN_DETECTED_BANDS
has_redshift = np.isfinite(galaxy_catalogue["zml"].to_numpy())
keep = enough_bands & has_redshift

splus = galaxy_catalogue[keep].reset_index(drop=True)
observed_magnitudes = all_magnitudes[keep]
observed_errors = all_errors[keep]
detected = all_detected[keep]
observed_z = splus["zml"].to_numpy()
observed_dz = sm.redshift_window(observed_z)

outside_prior = (observed_z < sm.Z_RANGE[0]) | (observed_z > sm.Z_RANGE[1])
print(f"{len(splus)} galaxies, {detected.sum(axis=1).mean():.1f} detected bands on average; "
      f"{outside_prior.sum()} with zml outside Z_RANGE")
```

- [ ] **Step 5: Cell 6 (markdown) and Cell 7 (code): noise model**

Cell 6:
```markdown
## 3. Noise model
Per band and 0.25 mag bin of measured magnitude, the distribution of DR3 errors (21 quantiles of $\log_{10}$ error). A synthetic band is detected when its drawn error is below 1 mag, so faint bins reproduce DR3's non-detection rate.
```

Cell 7:
```python
noise_table = sm.fit_noise_model(galaxy_catalogue)
noise_table.to_csv(NOISE_MODEL_FILE, index=False)
print(f"{len(noise_table)} (band, magnitude bin) rows saved to {NOISE_MODEL_FILE}")

# Median and 10-90% range of the errors in each band
figure, axes = plt.subplots(3, 4, figsize=(14, 9), sharex=True, sharey=True)
for axis, band in zip(axes.flat, sm.BANDS):
    band_rows = noise_table[noise_table["band"] == band]
    bin_centres = (band_rows["magnitude_low"] + band_rows["magnitude_high"]) / 2
    axis.fill_between(bin_centres, 10 ** band_rows["log_error_q02"], 10 ** band_rows["log_error_q18"],
                      color="grey", alpha=0.4, lw=0)
    axis.plot(bin_centres, 10 ** band_rows["log_error_q10"], color="black")
    axis.axhline(sm.MAX_DETECTED_ERROR, color="red", ls=":")
    axis.set(yscale="log", title=band)
for axis in axes[-1]:
    axis.set_xlabel("magnitude")
for axis in axes[:, 0]:
    axis.set_ylabel("error (mag)")
figure.tight_layout()
plt.show()
```

- [ ] **Step 6: Cell 8 (markdown) and Cell 9 (code): synthetic validation and test sets**

Cell 8:
```markdown
## 4. Synthetic validation and test sets
Independent draws with their own seeds. Each model takes the $r$ magnitude of a random DR3 galaxy at similar redshift, then gets S-PLUS errors, flux noise and non-detections. Mass-like labels are rescaled by the same magnitude offset.
```

Cell 9:
```python
catalogue_r = galaxy_catalogue["r_auto"].to_numpy()
catalogue_z = galaxy_catalogue["zml"].to_numpy()
held_out_seeds = {VALIDATION_FILE: VALIDATION_SEED, TEST_FILE: TEST_SEED}

for output_file, seed in held_out_seeds.items():
    rng = np.random.default_rng(seed)
    models = simulate_in_parallel(N_HELD_OUT, seed)
    model_magnitudes = models[sm.BANDS].to_numpy()

    # 1. Realistic brightness at each redshift
    scaled_magnitudes, offsets = sm.scale_to_observed_r(model_magnitudes, models["z"].to_numpy(),
                                                        catalogue_r, catalogue_z, rng)
    for label in sm.MASS_DEPENDENT_LABELS:
        models[label] = models[label] - 0.4 * offsets
    models["magnitude_offset"] = offsets

    # 2. S-PLUS noise and non-detections
    noisy_magnitudes, errors, detected_mask = sm.add_splus_noise(scaled_magnitudes, noise_table, rng)
    for band_index, band in enumerate(sm.BANDS):
        models[band] = noisy_magnitudes[:, band_index]
        models[f"e_{band}"] = errors[:, band_index]
    models.to_parquet(output_file, index=False)

    print(f"{output_file}: {len(models)} models, {detected_mask.sum(axis=1).mean():.1f} detected bands "
          f"on average (DR3: {all_detected.sum(axis=1).mean():.1f})")
```

- [ ] **Step 7: Cell 10 (markdown) and Cell 11 (code): coverage test v2**

Cell 10:
```markdown
## 5. Coverage test v2
For each galaxy, only models with $|z_{\rm model} - z_{\rm ml}| < \Delta z$ are compared. Over the detected bands, $\chi^2$ is minimised analytically over the magnitude offset (the model's mass). The galaxy is covered if the best $\chi^2$ has $p > 0.01$ for $n_{\rm detected} - 1$ degrees of freedom; the 0.05 mag error floor stands in for model error. Noisy synthetic galaxies should be covered ~99% of the time; that is the test's calibration check.
```

Cell 11:
```python
training_magnitudes = training_models[sm.BANDS].to_numpy()
training_z = training_models["z"].to_numpy()


def coverage_in_parallel(magnitudes, errors, detected_mask, z, dz):
    """sm.coverage_test on blocks of galaxies, on N_JOBS workers."""
    blocks = np.array_split(np.arange(len(z)), 4 * N_JOBS)
    tables = Parallel(n_jobs=N_JOBS)(
        delayed(sm.coverage_test)(magnitudes[block], errors[block], detected_mask[block], z[block], dz[block],
                                  training_magnitudes, training_z)
        for block in blocks
    )
    return pd.concat(tables, ignore_index=True)


# 1. S-PLUS galaxies
coverage = coverage_in_parallel(observed_magnitudes, observed_errors, detected, observed_z, observed_dz)
for column in coverage.columns:
    splus[column] = coverage[column].to_numpy()
print(f"Coverage v2 (p > {sm.P_MIN}): {splus['covered'].mean():.1%} of {len(splus)} galaxies")
print(f"Old geometric test (K = 2): {OLD_GEOMETRIC_COVERAGE:.1%}")

# 2. Calibration check on noisy synthetic galaxies
validation = pd.read_parquet(VALIDATION_FILE)
validation_magnitudes = validation[sm.BANDS].to_numpy()
validation_error_columns = []
for band in sm.BANDS:
    validation_error_columns.append(f"e_{band}")
validation_errors = validation[validation_error_columns].to_numpy()
validation_detected = sm.detected_bands(validation_magnitudes, validation_errors)
usable = validation_detected.sum(axis=1) >= sm.MIN_DETECTED_BANDS
validation_z = validation["z"].to_numpy()[usable]
validation_coverage = coverage_in_parallel(validation_magnitudes[usable], validation_errors[usable],
                                           validation_detected[usable], validation_z,
                                           sm.redshift_window(validation_z))
print(f"Synthetic validation galaxies covered: {validation_coverage['covered'].mean():.1%} "
      f"(expected about {1 - sm.P_MIN:.0%})")

# 3. Coverage by redshift, magnitude and number of detected bands
groupings = {
    "zml": pd.cut(splus["zml"], [-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 10]),
    "r_auto": pd.cut(splus["r_auto"], [0, 17, 18, 19, 20, 20.5, 21, 21.5, 30]),
    "n_detected": splus["n_dof"] + 1,
}
for name, grouping in groupings.items():
    grouped = splus.groupby(grouping, observed=True)
    table = pd.DataFrame({
        "covered": grouped["covered"].mean(),
        "median p": grouped["p_value"].median(),
        "galaxies": grouped.size(),
    })
    display(table.round(3))
```

- [ ] **Step 8: Cell 12 (markdown) and Cell 13 (code): figures**

Cell 12:
```markdown
## 6. Figures
Colour-colour diagrams of the synthetic domain with covered and uncovered galaxies, and the p-value distributions of S-PLUS and synthetic galaxies (uniform if models and errors are right).
```

Cell 13:
```python
from matplotlib.colors import LinearSegmentedColormap

# Poster palette
POSTER = {"panel": "#1b1b3a", "grid": "#2a2a52", "edge": "#3d3d70", "text": "#d6d6ec",
          "muted": "#8a8ab8", "cyan": "#5fd0f5", "red": "#e0524f", "red_edge": "#8c2424",
          "covered": "#b8b8d0"}
DOMAIN_COLOURMAP = LinearSegmentedColormap.from_list("poster_blue", ["#22305f", "#2f5fa0", "#4aa6e0", POSTER["cyan"]])
POSTER_STYLE = {"figure.facecolor": POSTER["panel"], "axes.facecolor": POSTER["panel"],
                "axes.edgecolor": POSTER["edge"], "axes.labelcolor": POSTER["text"], "axes.labelsize": 20,
                "xtick.color": POSTER["muted"], "ytick.color": POSTER["muted"], "text.color": POSTER["text"],
                "font.size": 16, "legend.frameon": False}
PANELS = [("g", "i", "r", "z"), ("u", "g", "g", "r"), ("g", "r", "r", "i"), ("r", "i", "i", "z")]
N_GALAXIES_SHOWN = 1500


def colour(magnitudes, first_band, second_band):
    """first_band - second_band from a (n, 12) magnitude array."""
    return magnitudes[:, sm.BANDS.index(first_band)] - magnitudes[:, sm.BANDS.index(second_band)]


covered = splus["covered"].to_numpy()
shown = np.random.default_rng(TRAINING_SEED).random(len(splus)) < N_GALAXIES_SHOWN / len(splus)

# Figure 1: colour-colour diagrams
with plt.rc_context(POSTER_STYLE):
    figure, axes = plt.subplots(2, 2, figsize=(12, 10))
    for axis, (band_a, band_b, band_c, band_d) in zip(axes.flat, PANELS):
        panel_detected = shown.copy()
        for band in (band_a, band_b, band_c, band_d):
            panel_detected = panel_detected & detected[:, sm.BANDS.index(band)]
        galaxy_x = colour(observed_magnitudes, band_a, band_b)[panel_detected]
        galaxy_y = colour(observed_magnitudes, band_c, band_d)[panel_detected]
        panel_covered = covered[panel_detected]

        x_limits = np.percentile(galaxy_x, [0.5, 99.5]) + [-0.8, 1.7]
        y_limits = np.percentile(galaxy_y, [0.5, 99.5]) + [-0.8, 0.95]
        axis.hexbin(colour(training_magnitudes, band_a, band_b), colour(training_magnitudes, band_c, band_d),
                    extent=(*x_limits, *y_limits), bins="log", mincnt=1, cmap=DOMAIN_COLOURMAP, zorder=1)
        axis.scatter(galaxy_x[panel_covered], galaxy_y[panel_covered], s=8, c=POSTER["covered"],
                     alpha=0.45, lw=0, zorder=2)
        axis.scatter(galaxy_x[~panel_covered], galaxy_y[~panel_covered], s=10, c=POSTER["red"],
                     edgecolors=POSTER["red_edge"], lw=0.5, zorder=3)
        axis.set(xlim=x_limits, ylim=y_limits, xlabel=f"${band_a} - {band_b}$", ylabel=f"${band_c} - {band_d}$")

    handles = [
        plt.Line2D([], [], ls="", marker="h", ms=10, color=POSTER["cyan"], label="Synthetic Domain"),
        plt.Line2D([], [], ls="", marker="o", ms=10, color=POSTER["covered"], label="Covered Galaxies"),
        plt.Line2D([], [], ls="", marker="o", ms=10, color=POSTER["red"], mec=POSTER["red_edge"], label="Uncovered Galaxies"),
    ]
    figure.legend(handles=handles, loc="upper center", ncol=3, labelcolor=POSTER["text"])
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig("color_color.png", dpi=300, facecolor=figure.get_facecolor())
    plt.show()

# Figure 2: p-value distributions
log_p_bins = np.linspace(-8, 0, 41)
with plt.rc_context(POSTER_STYLE):
    figure, axis = plt.subplots(figsize=(10, 7))
    splus_log_p = np.log10(np.clip(splus["p_value"].dropna(), 1e-8, 1))
    validation_log_p = np.log10(np.clip(validation_coverage["p_value"].dropna(), 1e-8, 1))
    axis.hist(splus_log_p, bins=log_p_bins, density=True, histtype="step", lw=2, color=POSTER["covered"], label="S-PLUS DR3")
    axis.hist(validation_log_p, bins=log_p_bins, density=True, histtype="step", lw=2, color=POSTER["cyan"], label="Synthetic validation")
    axis.axvline(np.log10(sm.P_MIN), color=POSTER["red"], ls="--", lw=2, label=f"$p = {sm.P_MIN}$")
    axis.set(xlabel=r"$\log_{10} p$", ylabel="density")
    axis.grid(True, color=POSTER["grid"])
    axis.legend(labelcolor=POSTER["text"])
    figure.tight_layout()
    figure.savefig("coverage_criterion.png", dpi=300, facecolor=figure.get_facecolor())
    plt.show()
```

- [ ] **Step 9: Run the notebook end to end (~35 min, run in background)**

Run:
```bash
../.venv/bin/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=inference4etgs training_set.ipynb
```
Expected: exits 0. In the saved outputs: `524288 models saved`, `Synthetic validation galaxies covered:` between 90% and 100%, coverage tables displayed, both PNGs written. If the synthetic coverage is far below 99%, report it (it measures grid density plus noise-model realism; it is a result, not a bug, unless the code is wrong).

- [ ] **Step 10: Run the full test suite**

Run: `../.venv/bin/python -m pytest -v`
Expected: 18 passed.

- [ ] **Step 11: Commit (and stop tracking the old CSV)**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git rm Synthetic_domain/training_set.csv
git add Synthetic_domain/training_set.ipynb Synthetic_domain/splus_noise_model.csv \
        Synthetic_domain/color_color.png Synthetic_domain/coverage_criterion.png
git commit -m "feat: training set v2 with 16-bin labels, noise model, held-out sets and chi-square coverage"
```

---

### Task 7: Spectrum preparation, likelihood, Firefly comparison and split

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (append section 5)
- Modify: `Synthetic_domain/tests/test_synthetic_model.py` (append)

**Interfaces:**
- Consumes: `components`, `uniform`, `DIMS` (Task 2).
- Produces:
  - `EXTINCTION_WAVELENGTHS` (12 values from `extinctions.py`), `R_V = 3.1`
  - `vacuum_to_air(wavelength) -> np.ndarray`
  - `deredden_spectrum(wavelength, flux, ivar, ebv) -> (flux, ivar)`
  - `deredden_magnitudes(magnitudes (n,12), ebv (n,)) -> np.ndarray`
  - `prepare_spectrum(wavelength_vacuum, flux, ivar, ebv) -> dict(wavelength, flux, ivar)`; raises `ValueError` if fewer than `MIN_GOOD_BINS` bins have data
  - `calibration_polynomial(model_flux, data_flux, data_ivar, order=CALIBRATION_ORDER) -> np.ndarray`
  - `spectral_log_likelihood(model_flux, data_flux, data_ivar, noise_scale) -> float`
  - `SPECTRAL_DIMS` (DIMS without `z`, plus `veldisp`, `noise_scale`)
  - `spectral_components(cube, redshift) -> (bagpipes_components, drawn, noise_scale)`; `drawn` also holds `veldisp`
  - `compare_with_firefly(bagpipes_age_gyr, bagpipes_metallicity, firefly_age_years, firefly_metallicity) -> pd.DataFrame` (`delta_log_age`, `delta_log_metallicity`, `firefly_consistent` as pandas `boolean` with NA)
  - `split_by_hash(specobjid) -> str` (`"validation"` or `"test"`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_synthetic_model.py`:
```python
# =============================================================================
# Spectra, likelihood, Firefly and split
# =============================================================================
import extinction


def test_vacuum_to_air_halpha():
    assert abs(sm.vacuum_to_air(np.array([6564.61]))[0] - 6562.80) < 0.02


def test_deredden_spectrum_undoes_ccm89_extinction():
    wavelength = np.linspace(3800, 9200, 500)
    intrinsic_flux = np.ones_like(wavelength)
    ebv = 0.1
    extinction_magnitudes = extinction.ccm89(wavelength, sm.R_V * ebv, sm.R_V)
    observed_flux = intrinsic_flux * 10 ** (-0.4 * extinction_magnitudes)
    observed_ivar = np.ones_like(wavelength)

    flux, ivar = sm.deredden_spectrum(wavelength, observed_flux, observed_ivar, ebv)

    assert np.allclose(flux, intrinsic_flux)
    assert np.all(ivar < observed_ivar)


def test_prepare_spectrum_rebins_and_masks():
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    flux = np.ones_like(wavelength_vacuum)
    ivar = np.ones_like(wavelength_vacuum)

    spectrum = sm.prepare_spectrum(wavelength_vacuum, flux, ivar, ebv=0.0)

    assert np.allclose(np.diff(spectrum["wavelength"]), sm.REBIN_WIDTH)
    assert spectrum["wavelength"][0] >= sm.OBSERVED_RANGE[0]
    assert spectrum["wavelength"][-1] <= sm.OBSERVED_RANGE[1]
    good = spectrum["ivar"] > 0
    assert np.allclose(spectrum["flux"][good], 1.0)
    sky_bin = np.argmin(np.abs(spectrum["wavelength"] - 5577.0))
    assert spectrum["ivar"][sky_bin] == 0


def test_prepare_spectrum_rejects_mostly_masked_spectrum():
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    flux = np.ones_like(wavelength_vacuum)
    ivar = np.zeros_like(wavelength_vacuum)
    ivar[:200] = 1.0
    with pytest.raises(ValueError, match="good bins"):
        sm.prepare_spectrum(wavelength_vacuum, flux, ivar, ebv=0.0)


def test_calibration_polynomial_recovers_known_polynomial():
    n_bins = 1000
    model_flux = 1 + 0.5 * np.sin(np.linspace(0, 20, n_bins))
    x = np.linspace(-1, 1, n_bins)
    true_polynomial = 2.0 + 0.3 * x - 0.1 * x**2
    data_flux = model_flux * true_polynomial
    data_ivar = np.ones(n_bins)

    polynomial = sm.calibration_polynomial(model_flux, data_flux, data_ivar)

    assert np.allclose(polynomial, true_polynomial)


def test_log_likelihood_is_highest_for_the_true_model():
    rng = np.random.default_rng(8)
    model_flux = 1 + 0.5 * np.sin(np.linspace(0, 20, 1000))
    data_flux = model_flux + rng.normal(0, 0.01, 1000)
    data_ivar = np.full(1000, 1 / 0.01**2)
    wrong_model = 1 + 0.5 * np.sin(np.linspace(0, 21, 1000))

    true_value = sm.spectral_log_likelihood(model_flux, data_flux, data_ivar, noise_scale=1.0)
    wrong_value = sm.spectral_log_likelihood(wrong_model, data_flux, data_ivar, noise_scale=1.0)

    assert true_value > wrong_value


def test_spectral_components_fix_redshift_and_add_veldisp():
    cube = np.full(len(sm.SPECTRAL_DIMS), 0.5)
    bagpipes_components, drawn, noise_scale = sm.spectral_components(cube, redshift=0.08)

    assert "z" not in sm.SPECTRAL_DIMS
    assert bagpipes_components["redshift"] == 0.08
    assert sm.VELDISP_RANGE[0] <= bagpipes_components["veldisp"] <= sm.VELDISP_RANGE[1]
    assert drawn["veldisp"] == bagpipes_components["veldisp"]
    assert 1.0 <= noise_scale <= 5.0


def test_compare_with_firefly_flags_consistency():
    result = sm.compare_with_firefly(
        bagpipes_age_gyr=np.array([5.0, 5.0]),
        bagpipes_metallicity=np.array([1.0, 1.0]),
        firefly_age_years=np.array([5.5e9, 20e9]),
        firefly_metallicity=np.array([0.9, 1.0]),
    )
    assert bool(result.loc[0, "firefly_consistent"])
    assert not bool(result.loc[1, "firefly_consistent"])


def test_compare_with_firefly_missing_values_are_na():
    result = sm.compare_with_firefly(
        bagpipes_age_gyr=np.array([5.0]),
        bagpipes_metallicity=np.array([1.0]),
        firefly_age_years=np.array([-9999.0]),
        firefly_metallicity=np.array([1.0]),
    )
    assert pd.isna(result.loc[0, "firefly_consistent"])
    assert np.isnan(result.loc[0, "delta_log_age"])


def test_split_by_hash_is_stable():
    identifiers = [299489677444933632, 299490502078654464, 1237645879578460255, 42]
    first = []
    for identifier in identifiers:
        first.append(sm.split_by_hash(identifier))
    # Same answer when asked again, in another order, with other galaxies around
    for identifier, expected in zip(reversed(identifiers), reversed(first)):
        assert sm.split_by_hash(identifier) == expected
    assert set(first).issubset({"validation", "test"})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `../.venv/bin/python -m pytest -v -k "vacuum or deredden or prepare or calibration or likelihood or spectral_components or firefly or split"`
Expected: FAIL with `AttributeError` for the new functions.

- [ ] **Step 3: Append section 5 to the module**

Append to `synthetic_model.py` (add `import extinction` to the imports at the top of the file):
```python
# =============================================================================
# 5. Spectra, spectral likelihood, Firefly comparison and split
# =============================================================================
R_V = 3.1
# Same S-PLUS wavelengths [Angstrom] and CCM89 law as extinctions.py
EXTINCTION_WAVELENGTHS = np.array([3536, 3770, 3940, 4094, 4292, 4751, 5133, 6258, 6614, 7690, 8611, 8831],
                                  dtype=float)

OBSERVED_RANGE = (3800.0, 9200.0)       # [Angstrom, observed] kept for the fit
SKY_LINES = (5577.0, 6300.0)            # [Angstrom] strong sky residuals
SKY_LINE_HALF_WIDTH = 10.0              # [Angstrom]
REBIN_WIDTH = 5.0                       # [Angstrom]
MIN_GOOD_BINS = 300                     # out of ~1080 bins
CALIBRATION_ORDER = 2

VELDISP_RANGE = (40.0, 400.0)                  # [km/s]
LOG_NOISE_SCALE_RANGE = (0.0, np.log10(5.0))   # error multiplier 1 to 5

FIREFLY_MISSING = -999.0                # Firefly writes -9999 for missing values
MAX_DELTA_LOG_AGE = 0.3                 # [dex]
MAX_DELTA_LOG_METALLICITY = 0.3         # [dex]


def vacuum_to_air(wavelength):
    """Air wavelengths [Angstrom] from vacuum ones (IAU standard, Morton 2000).

    SDSS spectra are in vacuum; the BC03-MILES models are in air.
    """
    wavenumber_squared = (1e4 / wavelength) ** 2
    refractive_index = (1 + 0.0000834254 + 0.02406147 / (130 - wavenumber_squared)
                        + 0.00015998 / (38.9 - wavenumber_squared))
    return wavelength / refractive_index


def deredden_spectrum(wavelength, flux, ivar, ebv):
    """Remove Milky Way extinction (CCM89, R_V = 3.1) from a spectrum. Returns (flux, ivar)."""
    extinction_magnitudes = extinction.ccm89(wavelength.astype(float), R_V * ebv, R_V)
    correction = 10 ** (0.4 * extinction_magnitudes)
    return flux * correction, ivar / correction**2


def deredden_magnitudes(magnitudes, ebv):
    """Remove Milky Way extinction from S-PLUS magnitudes (n, 12), as extinctions.py does."""
    corrected = np.array(magnitudes, dtype=float)
    for galaxy_index in range(len(corrected)):
        band_extinction = extinction.ccm89(EXTINCTION_WAVELENGTHS, R_V * ebv[galaxy_index], R_V)
        corrected[galaxy_index] = corrected[galaxy_index] - band_extinction
    return corrected


def prepare_spectrum(wavelength_vacuum, flux, ivar, ebv):
    """SDSS spectrum ready for the fit: air wavelengths, dereddened, masked, rebinned to 5 Angstrom.

    Each output bin is the inverse-variance-weighted mean of its unmasked pixels; bins without
    any get ivar = 0. Raises ValueError when fewer than MIN_GOOD_BINS bins have data.
    Returns dict(wavelength [Angstrom, air, observed], flux, ivar).
    """
    # 1. Air wavelengths and Milky Way extinction
    wavelength = vacuum_to_air(wavelength_vacuum)
    flux, ivar = deredden_spectrum(wavelength, flux, ivar, ebv)

    # 2. Pixels to ignore
    good = (ivar > 0) & np.isfinite(flux)
    for sky_line in SKY_LINES:
        near_sky_line = np.abs(wavelength - sky_line) < SKY_LINE_HALF_WIDTH
        good = good & ~near_sky_line

    # 3. Rebin onto a regular grid
    bin_centres = np.arange(OBSERVED_RANGE[0] + REBIN_WIDTH / 2, OBSERVED_RANGE[1], REBIN_WIDTH)
    rebinned_flux = np.zeros(len(bin_centres))
    rebinned_ivar = np.zeros(len(bin_centres))
    pixel_bin = np.floor((wavelength - OBSERVED_RANGE[0]) / REBIN_WIDTH).astype(int)
    for bin_index in range(len(bin_centres)):
        in_bin = good & (pixel_bin == bin_index)
        if not in_bin.any():
            continue
        weights = ivar[in_bin]
        rebinned_ivar[bin_index] = weights.sum()
        rebinned_flux[bin_index] = np.sum(weights * flux[in_bin]) / weights.sum()

    n_good_bins = int((rebinned_ivar > 0).sum())
    if n_good_bins < MIN_GOOD_BINS:
        raise ValueError(f"only {n_good_bins} good bins (need {MIN_GOOD_BINS})")
    return {"wavelength": bin_centres, "flux": rebinned_flux, "ivar": rebinned_ivar}


def calibration_polynomial(model_flux, data_flux, data_ivar, order=CALIBRATION_ORDER):
    """Multiplicative polynomial (in normalised wavelength) that best maps the model onto the data.

    Weighted linear least squares; the constant term absorbs the model normalisation, so the
    spectrum constrains the shape of the SFH, not the mass. Returns the polynomial per bin.
    """
    n_bins = len(model_flux)
    x = np.linspace(-1, 1, n_bins)
    design = np.zeros((n_bins, order + 1))
    for power in range(order + 1):
        design[:, power] = model_flux * x**power

    weighted_design = design * data_ivar[:, np.newaxis]
    normal_matrix = design.T @ weighted_design
    right_hand_side = weighted_design.T @ data_flux
    coefficients = np.linalg.solve(normal_matrix, right_hand_side)

    polynomial = np.zeros(n_bins)
    for power in range(order + 1):
        polynomial = polynomial + coefficients[power] * x**power
    return polynomial


def spectral_log_likelihood(model_flux, data_flux, data_ivar, noise_scale):
    """Gaussian log-likelihood of the data given a calibrated model, with errors x noise_scale.

    The calibration polynomial does not depend on noise_scale (it scales every weight equally).
    """
    polynomial = calibration_polynomial(model_flux, data_flux, data_ivar)
    residuals = data_flux - model_flux * polynomial
    good = data_ivar > 0
    chi2_value = np.sum(residuals[good] ** 2 * data_ivar[good]) / noise_scale**2
    return -0.5 * chi2_value - good.sum() * np.log(noise_scale)


SPECTRAL_DIMS = []
for dimension in DIMS:
    if dimension != "z":            # redshift is fixed to the spectroscopic value
        SPECTRAL_DIMS.append(dimension)
SPECTRAL_DIMS.append("veldisp")
SPECTRAL_DIMS.append("noise_scale")


def spectral_components(cube, redshift):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit unit-cube point.

    The SFH, metallicity, dust and nebular part is the training prior (sm.components) at the
    fixed spectroscopic redshift; veldisp and noise_scale are the spectral nuisance parameters.
    """
    cube_by_name = dict(zip(SPECTRAL_DIMS, cube))
    u = np.zeros(len(DIMS))
    for dimension_index, dimension in enumerate(DIMS):
        if dimension == "z":
            u[dimension_index] = 0.5    # ignored: fixed_redshift replaces it
        else:
            u[dimension_index] = cube_by_name[dimension]

    bagpipes_components, drawn = components(u, fixed_redshift=redshift)
    veldisp = uniform(cube_by_name["veldisp"], *VELDISP_RANGE)
    bagpipes_components["veldisp"] = veldisp
    drawn["veldisp"] = veldisp
    noise_scale = 10 ** uniform(cube_by_name["noise_scale"], *LOG_NOISE_SCALE_RANGE)
    return bagpipes_components, drawn, noise_scale


def compare_with_firefly(bagpipes_age_gyr, bagpipes_metallicity, firefly_age_years, firefly_metallicity):
    """Bagpipes vs Firefly mass-weighted age and metallicity (Z/Z_sun, both linear).

    Firefly missing values (-9999) become NaN, and the consistency flag is then NA, not False.
    Returns a DataFrame with delta_log_age, delta_log_metallicity [dex] and firefly_consistent.
    """
    firefly_age_years = np.where(firefly_age_years > FIREFLY_MISSING, firefly_age_years, np.nan)
    firefly_metallicity = np.where(firefly_metallicity > FIREFLY_MISSING, firefly_metallicity, np.nan)

    delta_log_age = np.log10(bagpipes_age_gyr * 1e9) - np.log10(firefly_age_years)
    delta_log_metallicity = np.log10(bagpipes_metallicity) - np.log10(firefly_metallicity)

    consistent = []
    for age_offset, metallicity_offset in zip(delta_log_age, delta_log_metallicity):
        if np.isnan(age_offset) or np.isnan(metallicity_offset):
            consistent.append(pd.NA)
        else:
            consistent.append(bool(abs(age_offset) < MAX_DELTA_LOG_AGE
                                   and abs(metallicity_offset) < MAX_DELTA_LOG_METALLICITY))

    return pd.DataFrame({
        "delta_log_age": delta_log_age,
        "delta_log_metallicity": delta_log_metallicity,
        "firefly_consistent": pd.array(consistent, dtype="boolean"),
    })


def split_by_hash(specobjid):
    """'validation' or 'test' from a hash of the SDSS specobjid.

    Depends only on the galaxy itself, so adding galaxies later never moves an existing one.
    """
    digest = hashlib.md5(str(int(specobjid)).encode()).hexdigest()
    if int(digest[:8], 16) % 2 == 0:
        return "validation"
    return "test"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `../.venv/bin/python -m pytest -v`
Expected: 28 passed.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: spectrum preparation, spectral likelihood, Firefly comparison and hash split"
```

---

### Task 8: MultiNest fit, posterior labels, mocks and fit files

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (append section 6)
- Modify: `Synthetic_domain/tests/test_synthetic_model.py` (append)

**Interfaces:**
- Consumes: `spectral_components`, `spectral_log_likelihood`, `bins_from_components`, `bin_labels`, `summary_labels`, `mass_weighted_metallicity`, `ab_magnitudes`, `FILTER_FILES`, `BANDS` (Tasks 2–7).
- Produces:
  - Constants `N_LIVE_POINTS = 300`, `SAMPLING_EFFICIENCY = 0.8`, `MAX_POSTERIOR_SAMPLES = 500`
  - `fit_spectrum(spectrum, redshift, output_basename, max_iterations, n_live_points=N_LIVE_POINTS) -> dict` (`samples`, `log_evidence`, `n_iterations`, `converged`, `run_time`)
  - `posterior_labels(unit_cube_samples, redshift) -> pd.DataFrame` (drawn params incl. `veldisp`, `noise_scale`, summary labels, bin labels, `mw_metallicity`, `model_r`)
  - `make_mock_spectrum(cube, redshift, template_spectrum, rng) -> (spectrum dict, truth pd.Series)`
  - `save_fit(path, posterior, unit_cube_samples, metadata)` / `load_fit(path) -> (pd.DataFrame, dict)`
  - `fit_and_save(job, output_directory, max_iterations, n_live_points=N_LIVE_POINTS) -> dict` with `job_id`, `status` (`"fitted"`, `"skipped"`, `"failed"`), `converged`, `run_time`, `n_iterations`, `error`; `job` is a dict with `job_id` (str), `spectrum` (dict), `redshift` (float)
  - `summarise_posterior(posterior, columns) -> dict` (`<column>_p16`, `<column>_median`, `<column>_p84`)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_synthetic_model.py`:
```python
# =============================================================================
# MultiNest fit and fit files (run with LD_LIBRARY_PATH=../.venv/lib)
# =============================================================================
def mock_job(job_id, rng):
    """A small mock spectrum job at z = 0.1."""
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    template = sm.prepare_spectrum(wavelength_vacuum, np.ones_like(wavelength_vacuum),
                                   np.full_like(wavelength_vacuum, 100.0), ebv=0.0)
    cube = rng.random(len(sm.SPECTRAL_DIMS))
    spectrum, truth = sm.make_mock_spectrum(cube, 0.1, template, rng)
    return {"job_id": job_id, "spectrum": spectrum, "redshift": 0.1}, truth


def test_make_mock_spectrum_returns_truth_labels():
    job, truth = mock_job("mock_test", np.random.default_rng(9))
    assert np.isfinite(job["spectrum"]["flux"]).all()
    for column in ("mw_age", "Z_fin", "f_old", "model_r"):
        assert np.isfinite(truth[column]), column


def test_fit_and_save_smoke_run(tmp_path):
    job, _ = mock_job("smoke", np.random.default_rng(10))
    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=100, n_live_points=40)

    assert summary["status"] == "fitted", summary["error"]
    posterior, metadata = sm.load_fit(str(tmp_path / "smoke.h5"))
    assert len(posterior) > 0
    for column in ("mw_age", "Z_fin", "f_old", "frac_bin_00", "model_r", "veldisp"):
        assert column in posterior.columns
    assert metadata["converged"] in (True, False)


def test_fit_and_save_skips_existing_output(tmp_path):
    job, _ = mock_job("skip", np.random.default_rng(11))
    posterior = pd.DataFrame({"mw_age": [1.0, 2.0]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "log_evidence": -1.0}
    sm.save_fit(str(tmp_path / "skip.h5"), posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), metadata)

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=100, n_live_points=40)

    assert summary["status"] == "skipped"
    assert summary["run_time"] == 5.0


def test_summarise_posterior_percentiles():
    posterior = pd.DataFrame({"mw_age": np.arange(101.0)})
    summary = sm.summarise_posterior(posterior, ["mw_age"])
    assert summary["mw_age_median"] == 50.0
    assert summary["mw_age_p16"] == 16.0
    assert summary["mw_age_p84"] == 84.0
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -v -k "mock or fit_and_save or summarise"`
Expected: FAIL with `AttributeError` for `make_mock_spectrum`, `fit_and_save`, `save_fit`, `summarise_posterior`.

- [ ] **Step 3: Append section 6 to the module**

Append to `synthetic_model.py` (add `import h5py` to the imports at the top):
```python
# =============================================================================
# 6. MultiNest fit, posterior labels, mocks and fit files
# =============================================================================
N_LIVE_POINTS = 300
SAMPLING_EFFICIENCY = 0.8        # MultiNest parameter-estimation mode
MAX_POSTERIOR_SAMPLES = 500      # posterior samples kept for labels


def fit_spectrum(spectrum, redshift, output_basename, max_iterations, n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one prepared spectrum under the training prior at a fixed redshift.

    The sampler works directly on the unit cube (SPECTRAL_DIMS): spectral_components maps it to
    physical values, exactly as the training set does. output_basename must be short (MultiNest
    limit of 100 characters), e.g. "spectral_fits/<id>/".
    Returns dict(samples, log_evidence, n_iterations, converged, run_time [s]).
    """
    import pymultinest   # needs LD_LIBRARY_PATH=.venv/lib (the inference4etgs kernel sets it)

    os.makedirs(output_basename, exist_ok=True)
    starting_components, _, _ = spectral_components(np.full(len(SPECTRAL_DIMS), 0.5), redshift)
    model = pipes.model_galaxy(starting_components, spec_wavs=spectrum["wavelength"])

    def log_likelihood(cube):
        bagpipes_components, _, noise_scale = spectral_components(cube, redshift)
        model.update(bagpipes_components)
        model_flux = model.spectrum[:, 1]
        return spectral_log_likelihood(model_flux, spectrum["flux"], spectrum["ivar"], noise_scale)

    def unit_cube_prior(cube):
        return cube

    start = time.time()
    result = pymultinest.solve(
        LogLikelihood=log_likelihood, Prior=unit_cube_prior, n_dims=len(SPECTRAL_DIMS),
        outputfiles_basename=output_basename, n_live_points=n_live_points,
        sampling_efficiency=SAMPLING_EFFICIENCY, multimodal=True, max_iter=max_iterations,
        resume=False, verbose=False, use_MPI=False,
    )
    run_time = time.time() - start

    # Each line of ev.dat is one MultiNest iteration (one dead point)
    with open(output_basename + "ev.dat") as dead_points:
        n_iterations = sum(1 for _ in dead_points)
    return {
        "samples": result["samples"],
        "log_evidence": result["logZ"],
        "n_iterations": n_iterations,
        "converged": n_iterations < max_iterations,
        "run_time": run_time,
    }


def posterior_labels(unit_cube_samples, redshift):
    """Labels of each posterior sample, from a Bagpipes model that also computes S-PLUS photometry.

    Keeps at most MAX_POSTERIOR_SAMPLES evenly spaced samples. model_r is the model r magnitude
    for 10^10 M_sun formed, used later to rescale mass-like labels.
    """
    step = max(1, int(np.ceil(len(unit_cube_samples) / MAX_POSTERIOR_SAMPLES)))
    kept_samples = unit_cube_samples[::step]

    model = None
    rows = []
    for cube in kept_samples:
        bagpipes_components, drawn, noise_scale = spectral_components(cube, redshift)
        if model is None:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
        else:
            model.update(bagpipes_components)
        fractions, metallicities, edges = bins_from_components(bagpipes_components)

        row = dict(drawn)
        row["noise_scale"] = noise_scale
        row.update(summary_labels(model.sfh))
        row.update(bin_labels(fractions, metallicities, edges))
        row["mw_metallicity"] = mass_weighted_metallicity(fractions, metallicities)
        row["model_r"] = ab_magnitudes(model.photometry)[BANDS.index("r")]
        rows.append(row)
    return pd.DataFrame(rows)


def make_mock_spectrum(cube, redshift, template_spectrum, rng):
    """Mock spectrum from a known unit-cube point, with the wavelengths, mask and errors of a real one.

    The model is scaled to the template's median flux so the mock has the template's S/N.
    The noise_scale coordinate is set to 0 (scale 1) because the mock noise is the template's.
    Returns (spectrum dict, truth labels as a Series).
    """
    cube = np.array(cube, dtype=float)
    cube[SPECTRAL_DIMS.index("noise_scale")] = 0.0
    bagpipes_components, _, _ = spectral_components(cube, redshift)
    model = pipes.model_galaxy(bagpipes_components, spec_wavs=template_spectrum["wavelength"])
    model_flux = model.spectrum[:, 1]

    good = template_spectrum["ivar"] > 0
    scale = np.median(template_spectrum["flux"][good]) / np.median(model_flux[good])
    errors = np.zeros(len(model_flux))
    errors[good] = 1 / np.sqrt(template_spectrum["ivar"][good])
    mock_flux = scale * model_flux + rng.normal(0, 1, len(model_flux)) * errors
    mock_flux[~good] = 0.0

    spectrum = {"wavelength": template_spectrum["wavelength"], "flux": mock_flux,
                "ivar": template_spectrum["ivar"].copy()}
    truth = posterior_labels(cube[np.newaxis, :], redshift).iloc[0]
    return spectrum, truth


def save_fit(path, posterior, unit_cube_samples, metadata):
    """Write posterior labels, unit-cube samples and fit metadata to one HDF5 file."""
    with h5py.File(path, "w") as output:
        output.create_dataset("unit_cube_samples", data=unit_cube_samples)
        labels_group = output.create_group("posterior")
        for column in posterior.columns:
            labels_group.create_dataset(column, data=posterior[column].to_numpy(dtype=float))
        for key, value in metadata.items():
            output.attrs[key] = value


def load_fit(path):
    """Read a file written by save_fit. Returns (posterior DataFrame, metadata dict)."""
    with h5py.File(path, "r") as source:
        columns = {}
        for column in source["posterior"]:
            columns[column] = source["posterior"][column][()]
        metadata = {}
        for key, value in source.attrs.items():
            metadata[key] = value.item() if hasattr(value, "item") else value
    return pd.DataFrame(columns), metadata


def fit_and_save(job, output_directory, max_iterations, n_live_points=N_LIVE_POINTS):
    """Fit one job (dict with job_id, spectrum, redshift) and save <output_directory>/<job_id>.h5.

    Skips jobs whose file already exists, so a batch can be re-run with more galaxies.
    Never raises: failures come back with status "failed" and the error text.
    """
    output_path = os.path.join(output_directory, f"{job['job_id']}.h5")
    summary = {"job_id": job["job_id"], "status": None, "converged": None,
               "run_time": None, "n_iterations": None, "error": None}

    if os.path.exists(output_path):
        _, metadata = load_fit(output_path)
        summary.update(status="skipped", converged=metadata["converged"],
                       run_time=metadata["run_time"], n_iterations=metadata["n_iterations"])
        return summary

    try:
        multinest_basename = os.path.join(output_directory, str(job["job_id"])) + "/"
        fit = fit_spectrum(job["spectrum"], job["redshift"], multinest_basename, max_iterations, n_live_points)
        posterior = posterior_labels(fit["samples"], job["redshift"])
        metadata = {"converged": fit["converged"], "run_time": fit["run_time"],
                    "n_iterations": fit["n_iterations"], "log_evidence": fit["log_evidence"],
                    "redshift": job["redshift"]}
        save_fit(output_path, posterior, fit["samples"], metadata)
        summary.update(status="fitted", converged=fit["converged"],
                       run_time=fit["run_time"], n_iterations=fit["n_iterations"])
    except Exception as error:
        summary.update(status="failed", error=repr(error))
    return summary


def summarise_posterior(posterior, columns):
    """16th, 50th and 84th percentiles of each column: <column>_p16, _median, _p84."""
    summary = {}
    for column in columns:
        values = posterior[column].to_numpy(dtype=float)
        summary[f"{column}_p16"] = np.nanpercentile(values, 16)
        summary[f"{column}_median"] = np.nanpercentile(values, 50)
        summary[f"{column}_p84"] = np.nanpercentile(values, 84)
    return summary
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -v`
Expected: 32 passed (the smoke fit takes about a minute).

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py
git commit -m "feat: MultiNest spectral fit under the training prior, posterior labels and fit files"
```

---

### Task 9: `spectroscopic_reference.ipynb`: sample, spectra and aperture photometry

**Files:**
- Create: `Synthetic_domain/spectroscopic_reference.ipynb` (cells 0–9)

**Interfaces:**
- Consumes: `sm.prepare_spectrum`, `sm.deredden_magnitudes`, `sm.Z_RANGE`, `sm.BANDS` (Tasks 2–7); `../Data/sdss_eboss_firefly-dr16.fits`; `dr3_corrected.csv`; `../dustMaps/csfd/`.
- Produces (notebook variables used by Task 10's cells): `selected` (DataFrame, one row per galaxy: `specobjid, ra, dec, z, plate, mjd, fiberid, fiber_diameter, sn_median, firefly_age, firefly_metallicity, firefly_mass, <band>_aper, e_<band>_aper`), `prepared_spectra` (dict specobjid → spectrum dict), `failures` (list of dicts).

The S-PLUS Cloud needs credentials: set `SPLUS_USERNAME` and `SPLUS_PASSWORD` in the environment before starting Jupyter. If they are not set, ask the human partner for them; do not hard-code them.

- [ ] **Step 1: Discover the S-PLUS aperture table and columns**

Run (repo root, with the credentials in the environment):
```bash
.venv/bin/python - <<'EOF'
import os
import splusdata
connection = splusdata.Core(os.environ["SPLUS_USERNAME"], os.environ["SPLUS_PASSWORD"])
tables = connection.query("SELECT table_schema, table_name FROM TAP_SCHEMA.tables")
print(type(tables)); print(tables)
columns = connection.query("SELECT table_name, column_name FROM TAP_SCHEMA.columns "
                           "WHERE column_name LIKE '%aper_3%' OR column_name LIKE '%APER_3%'")
print(columns)
EOF
```
Expected: a list of tables and the 3″-aperture magnitude/error columns. Record: (a) the table(s) holding them in DR4 (prefer DR4 over DR3 if both exist), (b) the exact magnitude and error column names per band, (c) the RA/Dec column names, (d) the Python type `query` returns. Fill these into the `APERTURE_SOURCES`, `SPLUS_RA_COLUMN` and `SPLUS_DEC_COLUMN` constants of cell 1 (the values shown there are the expected pattern; replace any that differ). If the return type is not an astropy Table, adapt `query_to_dataframe` in cell 9.

- [ ] **Step 2: Cell 0 (markdown)**

```markdown
# Spectroscopic reference set (part 2)

Real S-PLUS galaxies whose SFH labels come from **Bagpipes fits to their SDSS spectra under exactly the training prior** (`synthetic_model.py`), so they test the ML model on real data with labels it can match. Firefly (SDSS DR16) is an independent check of those fits, never a label.

**Inputs:** `dr3_corrected.csv`, `../Data/sdss_eboss_firefly-dr16.fits`, `../dustMaps/csfd/`, `training_set.parquet` (part 1), SDSS spectra (downloaded), S-PLUS Cloud aperture photometry (needs `SPLUS_USERNAME` / `SPLUS_PASSWORD`).
**Outputs:** `spectra/`, `spectral_fits/` (one `.h5` per fit), `spectral_fits/failures.csv`, `spectroscopic_reference.csv`.
**Run:** top to bottom with the *Inference4ETGs (MultiNest)* kernel, after `training_set.ipynb`. About 3 h; re-running skips finished fits, so raising `N_REAL_GALAXIES` later only fits the new galaxies.

Spectra fit 3″ (SDSS) or 2″ (BOSS) fibres, so photometry is the 3″ S-PLUS aperture, and BOSS galaxies are reported separately. The **test** half of the split is used once, for the final ML evaluation; with ~30/30 galaxies this is a pilot-size set.
```

- [ ] **Step 3: Cell 1 (code): imports and constants**

```python
import math
import os
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import fits
from astroquery.sdss import SDSS
from dustmaps.csfd import CSFDQuery
from joblib import Parallel, delayed

import synthetic_model as sm

# Files and folders
SPLUS_FILE = "dr3_corrected.csv"
FIREFLY_FILE = "../Data/sdss_eboss_firefly-dr16.fits"
DUST_MAP_FILE = "../dustMaps/csfd/csfd_ebv.fits"
DUST_MASK_FILE = "../dustMaps/csfd/mask.fits"
TRAINING_FILE = "training_set.parquet"
SPECTRA_DIRECTORY = "spectra"
FITS_DIRECTORY = "spectral_fits"
PILOT_DIRECTORY = "spectral_fits/pilot"
FAILURES_FILE = "spectral_fits/failures.csv"
OUTPUT_FILE = "spectroscopic_reference.csv"

# Sample
MIN_PROB_GALAXY = 0.8
MATCH_RADIUS = 1.5               # [arcsec]
N_REDSHIFT_BINS = 4
N_AGE_BINS = 3
GALAXIES_PER_CELL = 5
N_REAL_GALAXIES = N_REDSHIFT_BINS * N_AGE_BINS * GALAXIES_PER_CELL   # 60
FIBER_DIAMETER = {"26": 3.0, "v5_13_0": 2.0}   # [arcsec] SDSS legacy and BOSS
SQL_CHUNK = 200                  # specobjids per SDSS SQL query

# Fits
N_WORKERS = 11
N_MOCKS = 6                      # 2 per SFH branch
MOCK_SFH_DRAWS = [0.25, 0.25, 0.75, 0.75, 0.99, 0.99]   # sfh_type draws: Dirichlet, skew-normal, + burst
MIN_MOCKS_RECOVERED = 5
MOCK_CHECK_LABELS = ["mw_age", "Z_fin", "f_old"]
FIRST_WAVE_REAL_GALAXIES = N_WORKERS - N_MOCKS            # 5 real galaxies share the first wave
BUDGET_HOURS = 3.0
TARGET_MINUTES_PER_FIT = 25.0
ITERATION_SAFETY = 0.8           # MultiNest slows down as it converges
PILOT_ITERATIONS = 1000
SEED = 0

# S-PLUS aperture photometry (from Step 1's discovery; replace names that differ)
APERTURE_TABLE = "dr4_dual.dr4_dual"
SPLUS_RA_COLUMN = "RA"
SPLUS_DEC_COLUMN = "DEC"
APERTURE_SOURCES = {}
for band in sm.BANDS:
    APERTURE_SOURCES[band] = (APERTURE_TABLE, f"{band}_aper_3", f"e_{band}_aper_3")

SSFR_PASSIVE = -11.0             # [log yr^-1] star-forming above, passive below

os.makedirs(SPECTRA_DIRECTORY, exist_ok=True)
os.makedirs(PILOT_DIRECTORY, exist_ok=True)
failures = []
```

- [ ] **Step 4: Cell 2 (markdown) and Cell 3 (code): crossmatch with Firefly**

Cell 2:
```markdown
## 1. Sample
S-PLUS DR3 galaxies matched to Firefly DR16 spectra within 1.5″, with spectroscopic redshift inside the prior. One spectrum per S-PLUS galaxy (the closest match).
```

Cell 3:
```python
# 1. S-PLUS galaxies
splus = pd.read_csv(SPLUS_FILE)
splus = splus[splus["PROB_GAL"] > MIN_PROB_GALAXY].reset_index(drop=True)

# 2. Firefly galaxies with a redshift inside the prior, in the S-PLUS declination range
FIREFLY_COLUMNS = ["SPECOBJID", "RUN2D", "PLATE", "MJD", "FIBERID", "RA", "DEC", "Z_NOQSO", "CLASS_NOQSO",
                   "Chabrier_MILES_age_massW", "Chabrier_MILES_metallicity_massW", "Chabrier_MILES_total_mass"]
firefly_data = fits.getdata(FIREFLY_FILE, 1)
firefly = {}
for column in FIREFLY_COLUMNS:
    firefly[column] = np.asarray(firefly_data[column])
firefly = pd.DataFrame(firefly)
firefly["CLASS_NOQSO"] = firefly["CLASS_NOQSO"].str.strip()
firefly["RUN2D"] = firefly["RUN2D"].str.strip()
is_galaxy = firefly["CLASS_NOQSO"] == "GALAXY"
in_prior = (firefly["Z_NOQSO"] > sm.Z_RANGE[0]) & (firefly["Z_NOQSO"] < sm.Z_RANGE[1])
in_footprint = (firefly["DEC"] > splus["DEC"].min() - 1) & (firefly["DEC"] < splus["DEC"].max() + 1)
firefly = firefly[is_galaxy & in_prior & in_footprint].reset_index(drop=True)
del firefly_data

# 3. Nearest S-PLUS galaxy of each spectrum
spectrum_coordinates = SkyCoord(firefly["RA"].to_numpy() * u.deg, firefly["DEC"].to_numpy() * u.deg)
splus_coordinates = SkyCoord(splus["RA"].to_numpy() * u.deg, splus["DEC"].to_numpy() * u.deg)
nearest_splus, separation, _ = match_coordinates_sky(spectrum_coordinates, splus_coordinates)
firefly["splus_index"] = nearest_splus
firefly["separation"] = separation.arcsec
matched = firefly[firefly["separation"] < MATCH_RADIUS]

# 4. One spectrum per S-PLUS galaxy: the closest
matched = matched.sort_values("separation").drop_duplicates("splus_index").reset_index(drop=True)
# Kept as text: eBOSS specobjids exceed the int64 range
matched["specobjid"] = matched["SPECOBJID"].str.strip()
matched["fiber_diameter"] = matched["RUN2D"].map(FIBER_DIAMETER)
for column in ("Chabrier_MILES_age_massW", "Chabrier_MILES_metallicity_massW", "Chabrier_MILES_total_mass"):
    missing = matched[column] <= sm.FIREFLY_MISSING
    matched.loc[missing, column] = np.nan
print(f"{len(matched)} S-PLUS galaxies with an SDSS/BOSS spectrum; "
      f"{matched['fiber_diameter'].eq(3.0).sum()} SDSS, {matched['fiber_diameter'].eq(2.0).sum()} BOSS")
```

- [ ] **Step 5: Cell 4 (code): S/N and stratified selection**

```python
# 1. Median S/N of each spectrum from the SDSS database
sn_by_id = {}
identifiers = matched["specobjid"].tolist()
for start in range(0, len(identifiers), SQL_CHUNK):
    chunk = identifiers[start:start + SQL_CHUNK]
    id_list = ", ".join(chunk)
    query = ("SELECT CAST(specObjID AS VARCHAR(20)) AS specobjid_text, snMedian "
             f"FROM SpecObjAll WHERE specObjID IN ({id_list})")
    result = SDSS.query_sql(query, data_release=16)
    for row in result:
        sn_by_id[str(row["specobjid_text"]).strip()] = float(row["snMedian"])
matched["sn_median"] = matched["specobjid"].map(sn_by_id)

# 2. Grid of redshift x Firefly age quantile bins; highest S/N first in each cell
candidates = matched.dropna(subset=["Chabrier_MILES_age_massW", "sn_median"]).copy()
candidates["redshift_bin"] = pd.qcut(candidates["Z_NOQSO"], N_REDSHIFT_BINS, labels=False)
candidates["age_bin"] = pd.qcut(candidates["Chabrier_MILES_age_massW"], N_AGE_BINS, labels=False)

chosen = []
for (redshift_bin, age_bin), cell in candidates.groupby(["redshift_bin", "age_bin"]):
    best = cell.sort_values("sn_median", ascending=False).head(GALAXIES_PER_CELL).copy()
    best["rank_in_cell"] = np.arange(len(best))
    chosen.append(best)
# Round-robin order (best of every cell first), so a budget cut never empties a whole cell
selected = pd.concat(chosen, ignore_index=True)
selected = selected.sort_values(["rank_in_cell", "redshift_bin", "age_bin"]).reset_index(drop=True)
selected = selected.rename(columns={"RA": "ra", "DEC": "dec", "Z_NOQSO": "z", "PLATE": "plate", "MJD": "mjd",
                                    "FIBERID": "fiberid", "Chabrier_MILES_age_massW": "firefly_age",
                                    "Chabrier_MILES_metallicity_massW": "firefly_metallicity",
                                    "Chabrier_MILES_total_mass": "firefly_mass"})
print(f"{len(selected)} galaxies selected (target {N_REAL_GALAXIES})")
display(pd.crosstab(selected["redshift_bin"], selected["age_bin"]))
```

- [ ] **Step 6: Cell 5 (markdown) and Cell 6 (code): download and prepare spectra**

Cell 5:
```markdown
## 2. Spectra
SDSS spectra are **not** corrected for Galactic extinction, so each is dereddened with the same map and law as the S-PLUS photometry (CSFD, CCM89, $R_V = 3.1$), converted to air wavelengths (the models' frame), masked (bad pixels, 5577/6300 Å sky lines) and rebinned to 5 Å.
```

Cell 6:
```python
dust_map = CSFDQuery(map_fname=DUST_MAP_FILE, mask_fname=DUST_MASK_FILE)
galaxy_coordinates = SkyCoord(selected["ra"].to_numpy() * u.deg, selected["dec"].to_numpy() * u.deg)
selected["ebv"] = dust_map(galaxy_coordinates)

prepared_spectra = {}
for row in selected.itertuples():
    spectrum_file = os.path.join(SPECTRA_DIRECTORY, f"{row.specobjid}.fits")

    # 1. Download once
    if not os.path.exists(spectrum_file):
        try:
            downloaded = SDSS.get_spectra(plate=int(row.plate), mjd=int(row.mjd), fiberID=int(row.fiberid),
                                          data_release=16)
            downloaded[0].writeto(spectrum_file)
        except Exception as error:
            failures.append({"job_id": str(row.specobjid), "stage": "download", "error": repr(error)})
            continue

    # 2. Prepare for the fit
    with fits.open(spectrum_file) as spectrum_hdus:
        table = spectrum_hdus[1].data
        wavelength_vacuum = 10 ** table["loglam"].astype(float)
        flux = table["flux"].astype(float)
        ivar = table["ivar"].astype(float)
    try:
        prepared_spectra[row.specobjid] = sm.prepare_spectrum(wavelength_vacuum, flux, ivar, row.ebv)
    except ValueError as error:
        failures.append({"job_id": str(row.specobjid), "stage": "prepare", "error": repr(error)})

print(f"{len(prepared_spectra)} spectra ready; {len(failures)} failures so far")
```

- [ ] **Step 7: Cell 7 (markdown) and Cell 8 (code): aperture photometry**

Cell 7:
```markdown
## 3. Aperture photometry
3″-diameter S-PLUS magnitudes, matching the SDSS fibre, dereddened like `dr3_corrected.csv`. They set the mass normalisation of the labels and feed the coverage test; they never enter the spectral fit.
```

Cell 8:
```python
import splusdata


def query_to_dataframe(result):
    """S-PLUS Cloud query result as a DataFrame."""
    if hasattr(result, "to_pandas"):
        return result.to_pandas()
    return pd.DataFrame(result)


def cone_query(table, columns, ras, decs, radius_arcsec):
    """ADQL selecting columns of every object within radius of any of the given positions."""
    radius_degrees = radius_arcsec / 3600
    circles = []
    for ra, dec in zip(ras, decs):
        circles.append(f"1 = CONTAINS(POINT('ICRS', {SPLUS_RA_COLUMN}, {SPLUS_DEC_COLUMN}), "
                       f"CIRCLE('ICRS', {ra}, {dec}, {radius_degrees}))")
    column_list = ", ".join([SPLUS_RA_COLUMN, SPLUS_DEC_COLUMN, *columns])
    return f"SELECT {column_list} FROM {table} WHERE " + " OR ".join(circles)


connection = splusdata.Core(os.environ["SPLUS_USERNAME"], os.environ["SPLUS_PASSWORD"])

# 1. One query per table, then match back to the selected galaxies
columns_by_table = {}
for band, (table, magnitude_column, error_column) in APERTURE_SOURCES.items():
    columns_by_table.setdefault(table, []).extend([magnitude_column, error_column])

for table, columns in columns_by_table.items():
    query = cone_query(table, columns, selected["ra"], selected["dec"], MATCH_RADIUS)
    photometry = query_to_dataframe(connection.query(query))
    photometry_coordinates = SkyCoord(photometry[SPLUS_RA_COLUMN].to_numpy() * u.deg,
                                      photometry[SPLUS_DEC_COLUMN].to_numpy() * u.deg)
    nearest, separation, _ = match_coordinates_sky(galaxy_coordinates, photometry_coordinates)
    found = separation.arcsec < MATCH_RADIUS
    for band, (band_table, magnitude_column, error_column) in APERTURE_SOURCES.items():
        if band_table != table:
            continue
        magnitudes = np.full(len(selected), np.nan)
        errors = np.full(len(selected), np.nan)
        magnitudes[found] = photometry[magnitude_column].to_numpy()[nearest[found]]
        errors[found] = photometry[error_column].to_numpy()[nearest[found]]
        selected[f"{band}_aper"] = magnitudes
        selected[f"e_{band}_aper"] = errors

# 2. Galactic extinction, as for dr3_corrected.csv (non-detections keep their 99)
aperture_columns = []
aperture_error_columns = []
for band in sm.BANDS:
    aperture_columns.append(f"{band}_aper")
    aperture_error_columns.append(f"e_{band}_aper")
raw_magnitudes = selected[aperture_columns].to_numpy()
dereddened = sm.deredden_magnitudes(raw_magnitudes, selected["ebv"].to_numpy())
measured = raw_magnitudes < sm.UNDETECTED_MAGNITUDE
selected[aperture_columns] = np.where(measured, dereddened, raw_magnitudes)

missing_photometry = selected["r_aper"].isna()
print(f"Aperture photometry found for {(~missing_photometry).sum()} of {len(selected)} galaxies")
assert (~missing_photometry).sum() >= 0.9 * len(selected), "Aperture photometry missing for >10% of galaxies"
```

- [ ] **Step 8: Run cells 0–8 and check**

Run the notebook up to cell 8 in VS Code with the *Inference4ETGs (MultiNest)* kernel (or with nbconvert after temporarily keeping only cells 0–8 in a scratch copy).
Expected: ~1,400 matched galaxies printed; ~60 selected with a 4×3 crosstab of mostly 5s; ≥ 55 spectra ready; aperture photometry found for ≥ 90%.

- [ ] **Step 9: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/spectroscopic_reference.ipynb
git commit -m "feat: spectroscopic reference sample, dereddened SDSS spectra and S-PLUS aperture photometry"
```

---

### Task 10: `spectroscopic_reference.ipynb`: fits, Firefly check, coverage and split

**Files:**
- Modify: `Synthetic_domain/spectroscopic_reference.ipynb` (append cells 9–16)
- Generated (committed): `Synthetic_domain/spectroscopic_reference.csv`, `Synthetic_domain/firefly_comparison.png`

**Interfaces:**
- Consumes: `selected`, `prepared_spectra`, `failures` (Task 9); `sm.fit_and_save`, `sm.make_mock_spectrum`, `sm.load_fit`, `sm.summarise_posterior`, `sm.compare_with_firefly`, `sm.coverage_test`, `sm.split_by_hash` (Tasks 5, 7, 8); `training_set.parquet` (Task 6).
- Produces: `spectroscopic_reference.csv` (one row per selected galaxy: identifiers, Firefly values, aperture photometry, fit status, label percentiles, Firefly deltas and flag, coverage columns, `split`). The ML spec reads it.

- [ ] **Step 1: Cell 9 (markdown) and Cell 10 (code): pilot timing**

Cell 9:
```markdown
## 4. Fits
**Pilot:** one short fit measures the time per MultiNest iteration, which sets the iteration cap for ~25 min fits. **Inject-and-recover:** 6 mocks (2 per SFH branch) with real wavelengths, masks and errors run in the first wave with 5 real galaxies; at least 5 must recover the true `mw_age`, `Z_fin` and `f_old` inside their 95% intervals before the remaining fits start. Fits that hit the cap are flagged `converged = False`.
```

Cell 10:
```python
rng = np.random.default_rng(SEED)
pilot_galaxy = next(iter(prepared_spectra))
pilot_job = {"job_id": f"pilot_{pilot_galaxy}", "spectrum": prepared_spectra[pilot_galaxy],
             "redshift": float(selected.loc[selected["specobjid"] == pilot_galaxy, "z"].iloc[0])}
pilot = sm.fit_and_save(pilot_job, PILOT_DIRECTORY, PILOT_ITERATIONS)
assert pilot["status"] in ("fitted", "skipped"), pilot["error"]

seconds_per_iteration = pilot["run_time"] / pilot["n_iterations"]
max_iterations = int(ITERATION_SAFETY * TARGET_MINUTES_PER_FIT * 60 / seconds_per_iteration)
n_fits = N_MOCKS + len(prepared_spectra)
n_waves = math.ceil(n_fits / N_WORKERS)
print(f"Pilot: {seconds_per_iteration:.2f} s per iteration -> cap of {max_iterations} iterations per fit")
print(f"Projected: {n_fits} fits in {n_waves} waves of ~{TARGET_MINUTES_PER_FIT:.0f} min "
      f"= {n_waves * TARGET_MINUTES_PER_FIT / 60:.1f} h (budget {BUDGET_HOURS} h)")
```

- [ ] **Step 2: Cell 11 (code): first wave (mocks + 5 real galaxies) and the inject-and-recover check**

```python
def run_wave(jobs):
    """Fit a list of jobs on N_WORKERS processes; returns their summaries."""
    return Parallel(n_jobs=N_WORKERS, verbose=5)(
        delayed(sm.fit_and_save)(job, FITS_DIRECTORY, max_iterations) for job in jobs
    )


def real_galaxy_job(specobjid):
    """Fit job of one selected galaxy."""
    redshift = float(selected.loc[selected["specobjid"] == specobjid, "z"].iloc[0])
    return {"job_id": str(specobjid), "spectrum": prepared_spectra[specobjid], "redshift": redshift}


# 1. Mock jobs: real wavelengths, masks and errors; known unit-cube points
real_identifiers = list(prepared_spectra)
mock_jobs = []
mock_truths = {}
sfh_type_index = sm.SPECTRAL_DIMS.index("sfh_type")
for mock_index, sfh_type_draw in enumerate(MOCK_SFH_DRAWS):
    template_identifier = real_identifiers[rng.integers(len(real_identifiers))]
    template_job = real_galaxy_job(template_identifier)
    cube = rng.random(len(sm.SPECTRAL_DIMS))
    cube[sfh_type_index] = sfh_type_draw
    mock_spectrum, truth = sm.make_mock_spectrum(cube, template_job["redshift"], template_job["spectrum"], rng)
    job_id = f"mock_{mock_index}"
    mock_jobs.append({"job_id": job_id, "spectrum": mock_spectrum, "redshift": template_job["redshift"]})
    mock_truths[job_id] = truth

# 2. First wave
fitting_start = time.time()
first_wave_jobs = list(mock_jobs)
for specobjid in real_identifiers[:FIRST_WAVE_REAL_GALAXIES]:
    first_wave_jobs.append(real_galaxy_job(specobjid))
fit_summaries = run_wave(first_wave_jobs)
first_wave_hours = (time.time() - fitting_start) / 3600

# 3. Inject-and-recover: truth inside the 95% interval for every checked label
recovered = 0
for job in mock_jobs:
    posterior, _ = sm.load_fit(os.path.join(FITS_DIRECTORY, f"{job['job_id']}.h5"))
    inside_all = True
    for label in MOCK_CHECK_LABELS:
        low, high = np.nanpercentile(posterior[label], [2.5, 97.5])
        truth_value = mock_truths[job["job_id"]][label]
        inside = low <= truth_value <= high
        inside_all = inside_all and inside
        print(f"{job['job_id']} {label}: truth {truth_value:.3f}, 95% [{low:.3f}, {high:.3f}] {'ok' if inside else 'MISS'}")
    recovered += int(inside_all)
print(f"{recovered} of {len(mock_jobs)} mocks recovered; first wave took {first_wave_hours:.2f} h")
assert recovered >= MIN_MOCKS_RECOVERED, "Inject-and-recover failed: do not trust the real fits; investigate first"
```

- [ ] **Step 3: Cell 12 (code): remaining waves within the budget**

```python
# Fit as many of the remaining galaxies as the budget allows, at the measured wave time
remaining = real_identifiers[FIRST_WAVE_REAL_GALAXIES:]
hours_left = BUDGET_HOURS - (time.time() - fitting_start) / 3600
# On a re-run the first wave is skipped and takes no time; fall back to the target fit time
if first_wave_hours < 0.05:
    wave_hours = TARGET_MINUTES_PER_FIT / 60
else:
    wave_hours = first_wave_hours
waves_left = int(hours_left // wave_hours)
n_allowed = waves_left * N_WORKERS
if n_allowed < len(remaining):
    print(f"Budget allows {n_allowed} of the {len(remaining)} remaining galaxies; "
          f"fitting the first {n_allowed}. Raise BUDGET_HOURS and re-run to fit the rest.")
    remaining = remaining[:n_allowed]

remaining_jobs = []
for specobjid in remaining:
    remaining_jobs.append(real_galaxy_job(specobjid))
fit_summaries = fit_summaries + run_wave(remaining_jobs)

summary_table = pd.DataFrame(fit_summaries)
for record in fit_summaries:
    if record["status"] == "failed":
        failures.append({"job_id": record["job_id"], "stage": "fit", "error": record["error"]})
pd.DataFrame(failures).to_csv(FAILURES_FILE, index=False)
print(f"Fitting took {(time.time() - fitting_start) / 3600:.2f} h")
display(summary_table["status"].value_counts())
display(summary_table.groupby("status")[["run_time", "n_iterations"]].median())
```

- [ ] **Step 4: Cell 13 (markdown) and Cell 14 (code): labels, Firefly check, coverage and split**

Cell 13:
```markdown
## 5. Labels, Firefly check, coverage and split
Each posterior sample's mass-like labels are rescaled by its aperture $r$ offset ($-0.4\,\Delta m$ dex), then summarised as 16/50/84th percentiles. A fit is **Firefly-consistent** if mass-weighted age and metallicity agree within 0.3 dex; failures are flagged, not deleted. Coverage uses the aperture photometry at the spectroscopic redshift. Galaxies that converged, are covered and are Firefly-consistent go to `validation` or `test` by a hash of `specobjid`. Read the output with `dtype={"specobjid": str}` (eBOSS ids exceed int64).
```

Cell 14:
```python
LABEL_COLUMNS = ["log_ssfr", "mw_metallicity", "veldisp", "noise_scale", "Z_fin", "Av_ism",
                 *sm.SUMMARY_LABEL_COLUMNS, *sm.BIN_LABEL_COLUMNS]

# 1. Posterior summaries per galaxy
rows = []
for galaxy in selected.itertuples():
    fit_file = os.path.join(FITS_DIRECTORY, f"{galaxy.specobjid}.h5")
    row = {"specobjid": galaxy.specobjid}
    if not os.path.exists(fit_file):
        row["fit_status"] = "not_fitted"
        rows.append(row)
        continue
    posterior, metadata = sm.load_fit(fit_file)

    # Mass-like labels follow the aperture r magnitude: log X = log X_model - 0.4 (r_obs - r_model)
    offsets = galaxy.r_aper - posterior["model_r"]
    for label in sm.MASS_DEPENDENT_LABELS:
        posterior[label] = posterior[label] - 0.4 * offsets
    posterior["log_ssfr"] = posterior["log_sfr"] - posterior["log_mstar"]   # normalisation-free

    row["fit_status"] = "fitted"
    row["converged"] = bool(metadata["converged"])
    row["fit_hours"] = metadata["run_time"] / 3600
    row.update(sm.summarise_posterior(posterior, LABEL_COLUMNS))
    rows.append(row)
reference = selected.merge(pd.DataFrame(rows), on="specobjid", how="left")

# 2. Firefly check
firefly_check = sm.compare_with_firefly(
    reference["mw_age_median"].to_numpy(), reference["mw_metallicity_median"].to_numpy(),
    reference["firefly_age"].to_numpy(), reference["firefly_metallicity"].to_numpy(),
)
for column in firefly_check.columns:
    reference[column] = firefly_check[column].to_numpy()
reference["delta_log_mass"] = reference["log_mstar_median"] - np.log10(reference["firefly_mass"])

# 3. Coverage of the aperture photometry at the spectroscopic redshift
training = pd.read_parquet(TRAINING_FILE, columns=["z", *sm.BANDS])
aperture_magnitudes = reference[aperture_columns].to_numpy()
aperture_errors = reference[aperture_error_columns].to_numpy()
aperture_detected = sm.detected_bands(aperture_magnitudes, aperture_errors)
coverage = sm.coverage_test(aperture_magnitudes, aperture_errors, aperture_detected,
                            reference["z"].to_numpy(), np.full(len(reference), sm.DZ_MIN),
                            training[sm.BANDS].to_numpy(), training["z"].to_numpy())
for column in coverage.columns:
    reference[column] = coverage[column].to_numpy()

# 4. Locked split
usable = (reference["converged"].fillna(False).astype(bool)
          & reference["covered"].astype(bool)
          & reference["firefly_consistent"].fillna(False).astype(bool))
reference["split"] = "excluded"
for row_index in np.flatnonzero(usable.to_numpy()):
    reference.loc[row_index, "split"] = sm.split_by_hash(reference.loc[row_index, "specobjid"])

reference.to_csv(OUTPUT_FILE, index=False)
fitted = reference["fit_status"] == "fitted"
assert reference.loc[fitted, "mw_age_median"].notna().all(), "A fitted galaxy has no labels in the output"
print(f"{OUTPUT_FILE}: {len(reference)} galaxies, {fitted.sum()} fitted, "
      f"{reference['converged'].fillna(False).sum()} converged, {reference['covered'].sum()} covered, "
      f"{reference['firefly_consistent'].fillna(False).sum()} Firefly-consistent")
display(reference["split"].value_counts())
```

- [ ] **Step 5: Cell 15 (code): Firefly offsets by fibre and galaxy type, and figure**

```python
reference["fibre"] = reference["fiber_diameter"].map({3.0: "SDSS 3\"", 2.0: "BOSS 2\""})
reference["galaxy_type"] = np.where(reference["log_ssfr_median"] > SSFR_PASSIVE, "star-forming", "passive")

for grouping in ("fibre", "galaxy_type"):
    grouped = reference[fitted].groupby(grouping)
    table = pd.DataFrame({
        "median Δlog age": grouped["delta_log_age"].median(),
        "scatter Δlog age": grouped["delta_log_age"].std(),
        "median Δlog Z": grouped["delta_log_metallicity"].median(),
        "scatter Δlog Z": grouped["delta_log_metallicity"].std(),
        "median Δlog M": grouped["delta_log_mass"].median(),
        "consistent": grouped["firefly_consistent"].mean(),
        "galaxies": grouped.size(),
    })
    display(table.round(3))

figure, axes = plt.subplots(1, 2, figsize=(12, 5.5))
comparisons = [
    (np.log10(reference["firefly_age"]), np.log10(reference["mw_age_median"] * 1e9), "log age [yr]"),
    (np.log10(reference["firefly_metallicity"]), np.log10(reference["mw_metallicity_median"]), r"log $Z/Z_\odot$"),
]
for axis, (firefly_values, bagpipes_values, label) in zip(axes, comparisons):
    for galaxy_type, marker in (("star-forming", "o"), ("passive", "s")):
        chosen = fitted & (reference["galaxy_type"] == galaxy_type)
        axis.scatter(firefly_values[chosen], bagpipes_values[chosen], marker=marker, label=galaxy_type)
    limits = [np.nanmin([firefly_values, bagpipes_values]), np.nanmax([firefly_values, bagpipes_values])]
    axis.plot(limits, limits, color="black", lw=1)
    axis.set(xlabel=f"Firefly {label}", ylabel=f"Bagpipes {label}")
    axis.legend()
figure.tight_layout()
figure.savefig("firefly_comparison.png", dpi=200)
plt.show()
```

- [ ] **Step 6: Cell 16 (markdown): what the next pass changes**

```markdown
## Next pass
Raise `GALAXIES_PER_CELL` (and `BUDGET_HOURS`) and re-run: finished fits are skipped and the hash split never moves an existing galaxy. Add DR4 Stripe 82 galaxies by pointing `SPLUS_FILE` at a DR4 catalogue with the same columns.
```

- [ ] **Step 7: Run the whole notebook (~3.5 h, in the background)**

Run (credentials in the environment):
```bash
../.venv/bin/jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=inference4etgs spectroscopic_reference.ipynb
```
Expected: exits 0; `≥ 5 of 6 mocks recovered`; fitting time ≤ 3 h; `spectroscopic_reference.csv` written with a `split` value count of roughly equal `validation`/`test` plus `excluded`. If the inject-and-recover assertion fails, stop and report the printed per-label lines; do not loosen `MIN_MOCKS_RECOVERED`.

- [ ] **Step 8: Run the full test suite**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -v`
Expected: 32 passed.

- [ ] **Step 9: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/spectroscopic_reference.ipynb Synthetic_domain/spectroscopic_reference.csv \
        Synthetic_domain/firefly_comparison.png
git commit -m "feat: Bagpipes spectral fits of SDSS x S-PLUS galaxies, Firefly check and locked split"
```
