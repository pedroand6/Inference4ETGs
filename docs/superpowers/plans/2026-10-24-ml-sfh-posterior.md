# Phase 3 — ML Posterior for 16-bin SFHs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Train and evaluate a neural posterior estimator (conditional normalizing flow) that maps S-PLUS colours, r magnitude and redshift to a full posterior over each galaxy's SFH (15 populated age bins, metallicity, dust, masses), and apply it to S-PLUS DR3.

**Architecture:**
- `ml_data.py` turns model tables into network inputs and targets, simulating observations fresh each epoch with the shared `synthetic_model.py` noise model.
- `ml_networks.py` holds two swappable embeddings (MLP, 1-D CNN), a `zuko` neural spline flow and the training loop.
- `ml_evaluation.py` holds the metrics.
- `train_ml_model.py` trains on the user's RTX 3090 through `hljob`.
- `ml_model.ipynb` loads the trained runs, chooses the embedding on validation data, then evaluates once on the locked test data and writes the DR3 catalogue.

**Tech Stack:** Python 3.10 venv (`.venv`), numpy 1.26, pandas, pyarrow, torch 2.5.1, zuko, Bagpipes 1.3.5 (only through `synthetic_model.py`), pytest, `hljob` for GPU jobs.

**Spec:** `docs/superpowers/specs/2026-10-23-ml-sfh-posterior-design.md` (including its 2026-10-24 amendment: the redshift is the 19th target).

## Global Constraints

- **Code style:** the ten rules of the parts 1–2 spec apply verbatim:
  - English names; numbered banners; why-comments; no compact code;
  - docstrings with purpose, inputs with units and returns; named constants;
  - an opening markdown cell per notebook; logic in modules, orchestration in notebooks;
  - parallel cases in one pattern; concise.
- **Libraries:** `torch` (CUDA build on the GPU host, CPU build locally) and `zuko`. No other ML framework.
- **No re-implementation:** prior, noise model and label definitions come only from `synthetic_model.py` (imported as `sm`). Phase 3 never re-implements them; vectorised versions are added to `synthetic_model.py` itself, with tests against the originals.
- **Locked data:** `synthetic_test.parquet` and the reference set's `test` half are used once, in the final evaluation, behind `RUN_LOCKED_TEST = False` in `ml_model.ipynb`. No architecture, epoch or threshold choice may look at them.
- **Git:** weights (`ml_runs/*/model.pt`) are git-ignored. Configs, training histories and the catalogue summary are committed. The catalogue parquet itself is git-ignored.
- **Bins:** with the current prior, 15 age bins are populated (bin index 15 is always padding), and bin 14 exists only below z ≈ 0.48.
- **Inputs (37 numbers):**
  - 11 colours `m_b − m_r` (bands in `sm.BANDS` order, r excluded);
  - 11 colour errors `sqrt(σ_b² + σ_r²)`;
  - 11 detection flags;
  - `m_r`, `σ_r`, `z_obs`, `σ_z`.
  Galaxies without a detected r band are `no_r_band` and are never predicted.
- **Targets (19 numbers):**
  - `log10(f_i / f_13)` for bins i = 0..14, i ≠ 13 (14 numbers), with a non-existent bin set to `sm.MIN_FRACTION` before the transform;
  - `log10 Z_fin`, `Av_ism`, `log10 M_formed`, `log10 M_living` (`log_mstar`), and the true `redshift`.
- **Simulated observations per epoch:**
  1. brightness via `sm.scale_to_observed_r`;
  2. noise via `sm.add_splus_noise`;
  3. redshift: σ_z = 0 with probability 0.5, otherwise σ_z = `sm.SIGMA_Z`·(1 + z) and z_obs = z + N(0, σ_z).
- **Network:**
  - MLP embedding: 37 → 3 × 256 SiLU → 64.
  - CNN embedding: 12 bands × 4 channels [m_b − m_r, σ_b, detected, log10 pivot] → 3 × Conv1d(64, kernel 3, SiLU) → mean over bands, then concatenated with [m_r, σ_r, z_obs, σ_z] → 2 × 256 SiLU → 64.
  - Flow: `zuko.flows.NSF(features=19, context=64, transforms=8, hidden_features=(256, 256))`.
- **Training:** AdamW with lr 1e-3, cosine decay, batch 4096, ≤ 200 epochs, early stopping with patience 15 on the mean validation loss over the spectroscopic and photometric copies of `synthetic_validation.parquet`, keeping the best weights. Fixed seeds.
- **Model choice:** lower validation loss, with validation 68% / 95% coverage closest to nominal as the tie-breaker.
- **Evaluation is report-only:** no pass/fail gates.
- **Working directory for commands:** `/home/pedroand/Projects/Inference4ETGs/Synthetic_domain`. Python is `../.venv/bin/python`. Tests run with `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest`.

## Review Focus

1. **A galaxy with r undetected** (`no_r_band`): `colour_inputs` marks it invalid, it is never fed to the network, and the catalogue carries the flag. Test: Task 2, `test_colour_inputs_without_r_band_is_invalid`.
2. **A posterior redshift sample outside the prior** (photo-z near 0 or above 0.7): decoding clips it to `sm.Z_RANGE` and never crashes in `bin_edges_array`. Test: Task 2, `test_decode_clips_redshift_to_prior`.
3. **Extreme log-ratio samples** from the flow's tails: decoded fractions stay finite and sum to 1 (log-ratios clipped). Test: Task 2, `test_decode_extreme_log_ratios_stay_finite`.
4. **A galaxy detected only in r** (every colour masked): inputs contain no NaN and the network still returns finite samples. Tests: Task 2, `test_colour_inputs_only_r_detected_has_no_nan`; Task 3, `test_estimator_samples_finite_for_fully_masked_input`.
5. **A run trained on the GPU loaded on a CPU-only machine:** `load_run(..., device="cpu")` works and gives the same samples as before saving for a fixed seed. Test: Task 3, `test_save_load_round_trip_on_cpu`.

---

## File Structure

| File | Responsibility |
|---|---|
| `Synthetic_domain/synthetic_model.py` (modify) | Vectorised `draw_log_errors` and `scale_to_observed_r` (same API); new `bin_edges_array`, `bin_metallicities_array`. |
| `Synthetic_domain/ml_data.py` (create) | Inputs, targets, decoding, simulated observations, validation copies, `Standardiser`. |
| `Synthetic_domain/ml_networks.py` (create) | Embeddings, estimator, training loop, checkpoints, `posterior_samples`. |
| `Synthetic_domain/ml_evaluation.py` (create) | Coverage, ranks, bias / NMAD, width ratios, reference comparison. |
| `Synthetic_domain/train_ml_model.py` (create) | Training CLI. |
| `Synthetic_domain/remote/train_job.sh`, `remote/submit_training.sh` (create); `remote/build_env.sh` (modify) | GPU training jobs via `hljob`. |
| `Synthetic_domain/ml_model.ipynb` (create) | Part 3 story notebook. |
| `Synthetic_domain/tests/test_ml_data.py`, `test_ml_networks.py`, `test_ml_evaluation.py`, `test_train_ml_model.py` (create); `tests/test_synthetic_model.py` (modify) | Tests. |
| `.gitignore`, `requirements.txt` (modify) | Ignore weights / catalogue; add torch, zuko. |

---

### Task 1: Vectorised helpers in `synthetic_model.py` and the ML environment

**Files:**
- Modify: `Synthetic_domain/synthetic_model.py` (sections 1 and 3), `requirements.txt`, `.gitignore`
- Test: `Synthetic_domain/tests/test_synthetic_model.py`

**Interfaces:**
- Produces:
  - `sm.bin_edges_array(redshifts: np.ndarray) -> np.ndarray` of shape (n, N_BINS + 1), with edges [Gyr ago] youngest first and NaN beyond the last existing edge;
  - `sm.bin_metallicities_array(fractions: np.ndarray (n, N_BINS), n_existing: np.ndarray (n,), z_final: np.ndarray (n,)) -> np.ndarray (n, N_BINS)` in Z/Z☉, NaN for non-existent bins;
  - `sm.draw_log_errors` and `sm.scale_to_observed_r` keep their signatures and become vectorised.

- [ ] **Step 1: Install torch and zuko locally (CPU build)**

Run (from repo root):
```bash
.venv/bin/pip install "torch==2.5.1" --index-url https://download.pytorch.org/whl/cpu
.venv/bin/pip install "zuko==1.2.0"
.venv/bin/python -c "import torch, zuko; print(torch.__version__, zuko.__version__)"
```
Expected: `2.5.1+cpu 1.2.0`. If `zuko==1.2.0` is unavailable, install the newest zuko that supports torch 2.5.1, and write the version used in the report.

Append to `requirements.txt`:
```
torch==2.5.1
zuko==1.2.0
```
(Use the zuko version actually installed.)

Append to `.gitignore`:
```
Synthetic_domain/ml_runs/*/model.pt
Synthetic_domain/ml_sfh_catalogue.parquet
```

- [ ] **Step 2: Write the failing tests**

Append to `tests/test_synthetic_model.py`:
```python
# =============================================================================
# Vectorised helpers for the ML phase
# =============================================================================
def legacy_draw_log_errors(magnitudes, band_table, rng):
    """The original per-model loop of sm.draw_log_errors, kept as the reference."""
    table_bins = sm.noise_bin_index(band_table["magnitude_low"].to_numpy() + sm.NOISE_BIN_WIDTH / 2)
    quantile_table = band_table[sm.quantile_column_names()].to_numpy()
    log_errors = np.full(len(magnitudes), np.nan)
    finite = np.isfinite(magnitudes)
    model_bins = np.full(len(magnitudes), -1)
    model_bins[finite] = sm.noise_bin_index(magnitudes[finite])
    for model_index in np.flatnonzero(finite):
        if model_bins[model_index] > table_bins[-1]:
            continue
        row_index = max(np.searchsorted(table_bins, model_bins[model_index], side="right") - 1, 0)
        quantile_position = rng.random() * (sm.NOISE_N_QUANTILES - 1)
        lower = int(np.floor(quantile_position))
        upper = min(lower + 1, sm.NOISE_N_QUANTILES - 1)
        weight = quantile_position - lower
        row = quantile_table[row_index]
        log_errors[model_index] = row[lower] + weight * (row[upper] - row[lower])
    return log_errors


def test_vectorised_draw_log_errors_matches_the_legacy_loop(dr3_galaxies):
    noise_table = sm.fit_noise_model(dr3_galaxies)
    band_table = noise_table[noise_table["band"] == "g"]
    magnitudes = np.random.default_rng(20).uniform(8, 32, 5000)
    magnitudes[::17] = np.nan
    expected = legacy_draw_log_errors(magnitudes, band_table, np.random.default_rng(21))
    result = sm.draw_log_errors(magnitudes, band_table, np.random.default_rng(21))
    assert np.array_equal(np.isnan(expected), np.isnan(result))
    assert np.allclose(expected[np.isfinite(expected)], result[np.isfinite(result)])


def test_vectorised_scale_to_observed_r_draws_from_the_window():
    rng = np.random.default_rng(22)
    catalogue_z = rng.uniform(0, 1, 3000)
    catalogue_r = 15 + 5 * catalogue_z + rng.normal(0, 0.1, 3000)
    model_magnitudes = np.full((400, len(sm.BANDS)), 25.0)
    model_z = rng.uniform(0, 0.7, 400)
    scaled, offsets = sm.scale_to_observed_r(model_magnitudes, model_z, catalogue_r, catalogue_z, rng)
    r_index = sm.BANDS.index("r")
    for model_index in range(len(model_z)):
        in_window = np.abs(catalogue_z - model_z[model_index]) < sm.SCALING_MAX_DZ
        assert np.any(np.isclose(catalogue_r[in_window], scaled[model_index, r_index]))
    assert np.allclose(scaled[:, 0] - model_magnitudes[:, 0], offsets)


def test_bin_edges_array_matches_bin_edges():
    redshifts = np.random.default_rng(23).uniform(*sm.Z_RANGE, 300)
    edges_array = sm.bin_edges_array(redshifts)
    assert edges_array.shape == (300, sm.N_BINS + 1)
    for row_index, redshift in enumerate(redshifts):
        expected = sm.bin_edges(redshift)
        assert np.allclose(edges_array[row_index, :len(expected)], expected)
        assert np.all(np.isnan(edges_array[row_index, len(expected):]))


def test_bin_metallicities_array_matches_bin_metallicities():
    rng = np.random.default_rng(24)
    redshifts = rng.uniform(*sm.Z_RANGE, 200)
    n_existing = np.isfinite(sm.bin_edges_array(redshifts)).sum(axis=1) - 1
    fractions = rng.dirichlet(np.ones(sm.N_BINS), 200)
    z_final = rng.uniform(0.1, 2.5, 200)
    result = sm.bin_metallicities_array(fractions, n_existing, z_final)
    for row_index in range(200):
        n_bins = n_existing[row_index]
        expected = sm.bin_metallicities(fractions[row_index, :n_bins], z_final[row_index])
        assert np.allclose(result[row_index, :n_bins], expected)
        assert np.all(np.isnan(result[row_index, n_bins:]))
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_synthetic_model.py -v -k "vectorised or bin_edges_array or bin_metallicities_array"`
Expected: FAIL with `AttributeError: module 'synthetic_model' has no attribute 'bin_edges_array'`. The two vectorised tests may already pass (the API is unchanged); that is fine.

- [ ] **Step 4: Implement**

In `synthetic_model.py`, replace the body of `draw_log_errors` (keep the docstring and add one line saying it is vectorised and draws in the same order as the original loop):
```python
    table_bins = noise_bin_index(band_table["magnitude_low"].to_numpy() + NOISE_BIN_WIDTH / 2)
    quantile_table = band_table[quantile_column_names()].to_numpy()

    log_errors = np.full(len(magnitudes), np.nan)
    finite = np.isfinite(magnitudes)
    model_bins = np.full(len(magnitudes), -1)
    model_bins[finite] = noise_bin_index(magnitudes[finite])

    # Models fainter than the faintest populated bin get no error (undetected)
    drawable = finite & (model_bins <= table_bins[-1])
    drawable_indices = np.flatnonzero(drawable)

    # Nearest populated bin at or below each model's bin; the brightest bin covers anything brighter
    row_indices = np.searchsorted(table_bins, model_bins[drawable_indices], side="right") - 1
    row_indices = np.maximum(row_indices, 0)

    # One uniform draw per drawable model, in the same order as the original loop
    quantile_positions = rng.random(len(drawable_indices)) * (NOISE_N_QUANTILES - 1)
    lower = np.floor(quantile_positions).astype(int)
    upper = np.minimum(lower + 1, NOISE_N_QUANTILES - 1)
    weights = quantile_positions - lower
    lower_values = quantile_table[row_indices, lower]
    upper_values = quantile_table[row_indices, upper]
    log_errors[drawable_indices] = lower_values + weights * (upper_values - lower_values)
    return log_errors
```

Replace the per-model loop of `scale_to_observed_r` (everything from `offsets = np.full(...)` to the end) with:
```python
    # Window of real galaxies around each model's redshift
    first = np.searchsorted(sorted_z, model_z - max_dz)
    last = np.searchsorted(sorted_z, model_z + max_dz)
    has_neighbours = last > first

    offsets = np.full(len(model_magnitudes), np.nan)
    chosen = rng.integers(first[has_neighbours], last[has_neighbours])
    model_r = model_magnitudes[has_neighbours, r_index]
    offsets[has_neighbours] = sorted_r[chosen] - model_r

    scaled_magnitudes = model_magnitudes + offsets[:, np.newaxis]
    return scaled_magnitudes, offsets
```

Append to section 1 (after `bin_edges`) of `synthetic_model.py`:
```python
def bin_edges_array(redshifts):
    """Vectorised bin_edges: one row of age-bin edges per redshift.

    Inputs: redshifts: array of dimensionless redshifts (n,).
    Returns: array (n, N_BINS + 1) of edges [Gyr ago], youngest first, as bin_edges gives them;
             entries beyond a row's last edge (bins that do not exist at that redshift) are NaN.
    """
    redshifts = np.asarray(redshifts, dtype=float)
    maximum_ages = COSMOLOGY.age(redshifts).value - UNIVERSE_AGE_MARGIN
    all_midpoints = (AGE_BIN_CENTRES[:-1] + AGE_BIN_CENTRES[1:]) / 2

    edges = np.full((len(redshifts), N_BINS + 1), np.nan)
    for row_index, maximum_age in enumerate(maximum_ages):
        n_existing = int(np.sum(AGE_BIN_CENTRES <= maximum_age))
        edges[row_index, 0] = 0.0
        edges[row_index, 1:n_existing] = all_midpoints[:n_existing - 1]
        edges[row_index, n_existing] = maximum_age
    return edges


def bin_metallicities_array(fractions, n_existing, z_final):
    """Vectorised bin_metallicities for many galaxies whose number of existing bins differs.

    Inputs: fractions: mass fraction per bin (n, N_BINS), youngest first (entries beyond n_existing ignored);
            n_existing: number of existing bins per row (n,); z_final: Z_fin per row in Z/Z_sun (n,).
    Returns: Z/Z_sun per bin (n, N_BINS); NaN for bins that do not exist.
    """
    fractions = np.asarray(fractions, dtype=float)
    bin_index = np.arange(fractions.shape[1])
    exists = bin_index[np.newaxis, :] < np.asarray(n_existing)[:, np.newaxis]
    existing_fractions = np.where(exists, fractions, 0.0)

    # Mass formed before each bin = sum of the older existing bins (bins run youngest first)
    older_mass = np.cumsum(existing_fractions[:, ::-1], axis=1)[:, ::-1] - existing_fractions
    total = existing_fractions.sum(axis=1, keepdims=True)
    enrichment = (older_mass + existing_fractions / 2) / total

    z_final = np.asarray(z_final, dtype=float)[:, np.newaxis]
    z_initial = ZINI_FRACTION * z_final
    metallicities = z_initial + (z_final - z_initial) * enrichment
    return np.where(exists, metallicities, np.nan)
```

- [ ] **Step 5: Run the full suite**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -v`
Expected: all previous tests plus the 4 new ones pass (the noise-model statistics test still passes with the vectorised draw).

- [ ] **Step 6: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py requirements.txt .gitignore
git commit -m "feat: vectorised noise draws, brightness scaling, bin edges and metallicities for the ML phase"
```

---

### Task 2: `ml_data.py` — inputs, targets, decoding, simulated observations

**Files:**
- Create: `Synthetic_domain/ml_data.py`, `Synthetic_domain/tests/test_ml_data.py`

**Interfaces:**
- Consumes: `sm.BANDS`, `sm.PIVOT_WAVELENGTHS`, `sm.detected_bands`, `sm.scale_to_observed_r`, `sm.add_splus_noise`, `sm.rescale_mass_labels`, `sm.bin_edges_array`, `sm.bin_metallicities_array`, `sm.MIN_FRACTION`, `sm.LOG_MASS_FORMED`, `sm.SIGMA_Z`, `sm.Z_RANGE`, `sm.LABEL_EDGES`, `sm.SFR_FLOOR`, `sm.GYR_TO_YEARS`, `sm.DEX_PER_MAGNITUDE`, `sm.N_BINS`.
- Produces:
  - Constants: `COLOUR_BANDS` (11 names), `N_INPUTS = 37`, `N_BAND_CHANNELS = 4`, `N_FRACTION_BINS = 15`, `REFERENCE_BIN = 13`, `TARGET_NAMES` (19 names), `N_TARGETS = 19`, `SCALAR_INPUT_SLICE = slice(33, 37)`, `LOG_RATIO_LIMIT = 12.0`.
  - `colour_inputs(magnitudes, errors, z_obs, sigma_z) -> dict(flat=(n,37), banded=(n,12,4), valid=(n,) bool)`.
  - `encode_targets(table) -> np.ndarray (n, 19)` (the table needs `frac_bin_00..14`, `Z_fin`, `Av_ism`, `log_mass_formed`, `log_mstar`, `redshift`).
  - `decode_samples(theta) -> pd.DataFrame` with columns `redshift, Z_fin, Av_ism, log_mass_formed, log_mstar, log_sfr, mw_age, f_young, f_mid, f_old` + `sm.BIN_LABEL_COLUMNS`.
  - `simulate_observations(table, catalogue_r, catalogue_z, noise_table, rng) -> dict(flat, banded, valid, targets)`.
  - `validation_observations(held_out, mode, rng) -> dict(flat, banded, valid, targets)` with `mode` in `{"spectroscopic", "photometric"}`.
  - `class Standardiser` with `fit(values, mask=None)` (classmethod), `transform(values, mask=None)`, `inverse_transform(values)`, `to_dict()`, `from_dict(d)` (classmethod).

- [ ] **Step 1: Write the failing tests**

Create `Synthetic_domain/tests/test_ml_data.py`:
```python
"""Tests for ml_data.py. Run from Synthetic_domain/: LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -v"""

import numpy as np
import pandas as pd
import pytest

import ml_data
import synthetic_model as sm


@pytest.fixture(scope="module")
def simulated_models():
    """Twenty prior draws with labels and noiseless magnitudes."""
    unit_vectors = np.random.default_rng(30).random((20, len(sm.DIMS)))
    models = sm.simulate_models(unit_vectors)
    models["log_mass_formed"] = sm.LOG_MASS_FORMED
    return models


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


def test_encode_decode_round_trip(simulated_models):
    theta = ml_data.encode_targets(simulated_models)
    assert theta.shape == (len(simulated_models), ml_data.N_TARGETS)
    decoded = ml_data.decode_samples(theta)
    for column in ["redshift", "Z_fin", "Av_ism", "log_mstar"]:
        assert np.allclose(decoded[column], simulated_models[column])
    for bin_index in range(sm.N_BINS):
        for kind in ("frac_bin", "logsfr_bin", "zmet_bin"):
            column = f"{kind}_{bin_index:02d}"
            expected = simulated_models[column].to_numpy(dtype=float)
            result = decoded[column].to_numpy(dtype=float)
            assert np.array_equal(np.isnan(expected), np.isnan(result)), column
            finite = np.isfinite(expected)
            assert np.allclose(result[finite], expected[finite], rtol=1e-6, atol=1e-9), column


def test_decoded_summary_labels_match_bagpipes(simulated_models):
    decoded = ml_data.decode_samples(ml_data.encode_targets(simulated_models))
    assert np.allclose(decoded["mw_age"], simulated_models["mw_age"], rtol=0.02)
    for column in ("f_young", "f_mid", "f_old"):
        assert np.allclose(decoded[column], simulated_models[column], atol=0.01), column
    star_forming = simulated_models["log_sfr"] > np.log10(sm.SFR_FLOOR) + simulated_models["log_mstar"] + 0.5
    assert np.allclose(decoded["log_sfr"][star_forming], simulated_models["log_sfr"][star_forming], atol=0.05)


def test_decode_clips_redshift_to_prior(simulated_models):
    theta = ml_data.encode_targets(simulated_models.iloc[:2])
    theta[0, -1] = -0.3
    theta[1, -1] = 1.5
    decoded = ml_data.decode_samples(theta)
    assert decoded["redshift"].between(*sm.Z_RANGE).all()


def test_decode_extreme_log_ratios_stay_finite(simulated_models):
    theta = ml_data.encode_targets(simulated_models.iloc[:3])
    theta[:, :14] = np.array([[500.0], [-500.0], [0.0]])
    decoded = ml_data.decode_samples(theta)
    fraction_columns = [f"frac_bin_{bin_index:02d}" for bin_index in range(sm.N_BINS)]
    sums = decoded[fraction_columns].sum(axis=1)
    assert np.allclose(sums, 1.0)
    assert np.isfinite(decoded[fraction_columns].to_numpy()).all()


@pytest.fixture(scope="module")
def noise_and_catalogue():
    catalogue = pd.read_csv("dr3_corrected.csv")
    catalogue = catalogue[catalogue["PROB_GAL"] > 0.8]
    noise_table = sm.fit_noise_model(catalogue)
    return noise_table, catalogue["r_auto"].to_numpy(), catalogue["zml"].to_numpy()


def test_simulate_observations_shifts_mass_and_draws_redshift(simulated_models, noise_and_catalogue):
    noise_table, catalogue_r, catalogue_z = noise_and_catalogue
    observations = ml_data.simulate_observations(simulated_models, catalogue_r, catalogue_z, noise_table,
                                                 np.random.default_rng(31))
    targets = observations["targets"]
    flat = observations["flat"]
    log_mass_formed_index = ml_data.TARGET_NAMES.index("log_mass_formed")
    r_magnitude = flat[:, 33]
    model_r = simulated_models["r"].to_numpy()
    r_error = flat[:, 34]
    precise = observations["valid"] & (r_error < 0.05)
    expected_mass = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * (r_magnitude - model_r)
    # The noisy r differs from the scaled r by its noise: compare only galaxies with a precise r
    assert precise.sum() >= 5
    assert np.allclose(targets[precise, log_mass_formed_index], expected_mass[precise], atol=0.1)
    sigma_z = flat[:, 36]
    z_obs = flat[:, 35]
    spectroscopic = sigma_z == 0
    assert np.allclose(z_obs[spectroscopic], simulated_models["redshift"].to_numpy()[spectroscopic])
    photometric = ~spectroscopic
    assert np.allclose(sigma_z[photometric], sm.SIGMA_Z * (1 + simulated_models["redshift"].to_numpy()[photometric]))


def test_validation_observations_modes():
    held_out = pd.read_parquet("synthetic_validation.parquet").iloc[:200]
    spectroscopic = ml_data.validation_observations(held_out, "spectroscopic", np.random.default_rng(32))
    photometric = ml_data.validation_observations(held_out, "photometric", np.random.default_rng(32))
    assert np.all(spectroscopic["flat"][:, 36] == 0)
    assert np.all(photometric["flat"][:, 36] > 0)
    assert np.allclose(spectroscopic["flat"][:, 35], held_out["redshift"])
    log_mass_formed_index = ml_data.TARGET_NAMES.index("log_mass_formed")
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

- [ ] **Step 2: Run the tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_data.py -v`
Expected: collection error `ModuleNotFoundError: No module named 'ml_data'`.

- [ ] **Step 3: Implement `ml_data.py`**

Create `Synthetic_domain/ml_data.py`:
```python
"""
Network inputs and targets for the SFH posterior estimator (Phase 3).

Model tables (training set, held-out sets, real galaxies) become:
- inputs: 11 colours (band - r) with errors and detection flags, plus r, its error, z_obs and sigma_z;
- targets: 19 numbers (bin log-ratios, metallicity, dust, masses, redshift) that a flow can model.
Every physical definition (noise, brightness, bins, enrichment) comes from synthetic_model.
"""

import numpy as np
import pandas as pd

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

N_FRACTION_BINS = sm.N_BINS - 1                  # 15 populated bins (index 15 never exists)
REFERENCE_BIN = 13                               # oldest bin that exists at every redshift of the prior
LOG_RATIO_LIMIT = 12.0                           # [dex] clip of sampled log-ratios (flow tails)
P_SPECTROSCOPIC = 0.5                            # share of simulated galaxies with an exact redshift

RATIO_BINS = []
for bin_index in range(N_FRACTION_BINS):
    if bin_index != REFERENCE_BIN:
        RATIO_BINS.append(bin_index)
TARGET_NAMES = []
for bin_index in RATIO_BINS:
    TARGET_NAMES.append(f"log_ratio_bin_{bin_index:02d}")
TARGET_NAMES += ["log_Z_fin", "Av_ism", "log_mass_formed", "log_mstar", "redshift"]
N_TARGETS = len(TARGET_NAMES)                    # 19


# =============================================================================
# 2. Inputs
# =============================================================================
def colour_inputs(magnitudes, errors, z_obs, sigma_z):
    """Network inputs from magnitudes: colours against r, their errors, flags and the scalars.

    Inputs: magnitudes, errors: [mag] arrays (n, 12) in sm.BANDS order (NaN or 99 when undetected);
            z_obs: observed redshift (n,); sigma_z: its uncertainty (n,), 0 for spectroscopic redshifts.
    Returns: dict with
             flat (n, 37): [11 colours, 11 colour errors, 11 flags, m_r, sigma_r, z_obs, sigma_z] for the MLP;
             banded (n, 12, 4): per band [m_b - m_r, sigma_b, detected, log10 pivot] for the CNN;
             valid (n,): True where r is detected (only those galaxies may be predicted).
    Unusable entries are 0, never NaN.
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    n_galaxies = len(magnitudes)
    detected = sm.detected_bands(magnitudes, errors)
    valid = detected[:, R_INDEX]

    r_magnitude = np.where(valid, magnitudes[:, R_INDEX], 0.0)
    r_error = np.where(valid, errors[:, R_INDEX], 0.0)

    flat = np.zeros((n_galaxies, N_INPUTS))
    banded = np.zeros((n_galaxies, len(sm.BANDS), N_BAND_CHANNELS))
    log_pivots = np.log10(sm.PIVOT_WAVELENGTHS)

    for colour_index, band in enumerate(COLOUR_BANDS):
        band_index = sm.BANDS.index(band)
        usable = detected[:, band_index] & valid
        colour = np.where(usable, magnitudes[:, band_index] - r_magnitude, 0.0)
        colour_error = np.where(usable, np.sqrt(errors[:, band_index] ** 2 + r_error**2), 0.0)
        flat[:, colour_index] = colour
        flat[:, N_COLOURS + colour_index] = colour_error
        flat[:, 2 * N_COLOURS + colour_index] = usable.astype(float)

    flat[:, SCALAR_INPUT_SLICE] = np.column_stack([r_magnitude, r_error, z_obs, sigma_z])

    for band_index in range(len(sm.BANDS)):
        usable = detected[:, band_index] & valid
        banded[:, band_index, 0] = np.where(usable, magnitudes[:, band_index] - r_magnitude, 0.0)
        banded[:, band_index, 1] = np.where(usable, errors[:, band_index], 0.0)
        banded[:, band_index, 2] = usable.astype(float)
        banded[:, band_index, 3] = log_pivots[band_index]
    return {"flat": flat, "banded": banded, "valid": valid}


# =============================================================================
# 3. Targets and decoding
# =============================================================================
def fraction_columns(n_bins):
    """Names frac_bin_00 ... of the first n_bins bins."""
    names = []
    for bin_index in range(n_bins):
        names.append(f"frac_bin_{bin_index:02d}")
    return names


def encode_targets(table):
    """The 19 flow targets of each row.

    Inputs: table with frac_bin_00..14, Z_fin [Z_sun], Av_ism [mag], log_mass_formed and log_mstar [log M_sun],
            redshift.
    Returns: array (n, 19) in TARGET_NAMES order. Bins that do not exist (fraction 0) take sm.MIN_FRACTION
             before the log-ratio, so every row has the same dimension.
    """
    fractions = table[fraction_columns(N_FRACTION_BINS)].to_numpy(dtype=float)
    fractions = np.where(fractions > 0, fractions, sm.MIN_FRACTION)
    reference = fractions[:, REFERENCE_BIN]

    columns = []
    for bin_index in RATIO_BINS:
        columns.append(np.log10(fractions[:, bin_index] / reference))
    columns.append(np.log10(table["Z_fin"].to_numpy(dtype=float)))
    columns.append(table["Av_ism"].to_numpy(dtype=float))
    columns.append(table["log_mass_formed"].to_numpy(dtype=float))
    columns.append(table["log_mstar"].to_numpy(dtype=float))
    columns.append(table["redshift"].to_numpy(dtype=float))
    return np.column_stack(columns)


def mass_formed_younger_than(fractions, edges, age_limit):
    """Share of the formed mass younger than age_limit, with constant SFR inside each bin.

    Inputs: fractions (n, N_BINS) summing to 1; edges (n, N_BINS + 1) [Gyr ago], NaN beyond the last edge;
            age_limit [Gyr].
    Returns: (n,) fraction of the mass formed more recently than age_limit.
    """
    starts = edges[:, :-1]
    ends = edges[:, 1:]
    widths = ends - starts
    overlap = np.clip(np.minimum(ends, age_limit) - starts, 0, None)
    share_of_bin = np.where(widths > 0, overlap / widths, 0.0)
    share_of_bin = np.nan_to_num(share_of_bin)
    return np.sum(fractions * share_of_bin, axis=1)


def decode_samples(theta):
    """Physical labels of flow samples, with the same definitions as the training labels.

    Inputs: theta: array (n, 19) in TARGET_NAMES order (posterior samples or encoded truths).
    Returns: DataFrame with redshift, Z_fin, Av_ism, log_mass_formed, log_mstar, log_sfr, mw_age, f_young,
             f_mid, f_old and the 48 sm.BIN_LABEL_COLUMNS. The redshift is clipped to the prior; bins that do
             not exist at the sampled redshift get frac 0, logsfr and zmet NaN, as in sm.bin_labels.
    """
    theta = np.asarray(theta, dtype=float)
    n_samples = len(theta)
    redshift = np.clip(theta[:, TARGET_NAMES.index("redshift")], *sm.Z_RANGE)
    z_final = 10 ** theta[:, TARGET_NAMES.index("log_Z_fin")]
    log_mass_formed = theta[:, TARGET_NAMES.index("log_mass_formed")]
    log_mstar = theta[:, TARGET_NAMES.index("log_mstar")]

    # 1. Bins that exist at each sampled redshift
    edges = sm.bin_edges_array(redshift)
    n_existing = np.isfinite(edges).sum(axis=1) - 1
    exists = np.arange(sm.N_BINS)[np.newaxis, :] < n_existing[:, np.newaxis]

    # 2. Fractions from the log-ratios (clipped: the flow's tails can be extreme)
    log_ratios = np.clip(theta[:, :len(RATIO_BINS)], -LOG_RATIO_LIMIT, LOG_RATIO_LIMIT)
    relative_mass = np.zeros((n_samples, sm.N_BINS))
    relative_mass[:, REFERENCE_BIN] = 1.0
    for ratio_index, bin_index in enumerate(RATIO_BINS):
        relative_mass[:, bin_index] = 10 ** log_ratios[:, ratio_index]
    relative_mass = np.where(exists, relative_mass, 0.0)
    fractions = relative_mass / relative_mass.sum(axis=1, keepdims=True)

    # 3. Per-bin labels, as sm.bin_labels defines them
    widths_in_years = np.diff(edges, axis=1) * sm.GYR_TO_YEARS
    with np.errstate(divide="ignore", invalid="ignore"):
        log_sfr_bins = log_mass_formed[:, np.newaxis] + np.log10(fractions) - np.log10(widths_in_years)
    log_sfr_bins = np.where(exists, log_sfr_bins, np.nan)
    metallicities = sm.bin_metallicities_array(fractions, n_existing, z_final)

    # 4. Summary labels from the bins (constant SFR inside each bin)
    midpoints = np.nan_to_num((edges[:, :-1] + edges[:, 1:]) / 2)
    mw_age = np.sum(fractions * midpoints, axis=1)
    younger_than_first = mass_formed_younger_than(fractions, edges, sm.LABEL_EDGES[0])
    younger_than_second = mass_formed_younger_than(fractions, edges, sm.LABEL_EDGES[1])
    recent_mass = 10**log_mass_formed * younger_than_first
    recent_sfr = recent_mass / (sm.LABEL_EDGES[0] * sm.GYR_TO_YEARS)
    sfr_floor = sm.SFR_FLOOR * 10**log_mstar
    log_sfr = np.log10(np.maximum(recent_sfr, sfr_floor))

    labels = {
        "redshift": redshift,
        "Z_fin": z_final,
        "Av_ism": theta[:, TARGET_NAMES.index("Av_ism")],
        "log_mass_formed": log_mass_formed,
        "log_mstar": log_mstar,
        "log_sfr": log_sfr,
        "mw_age": mw_age,
        "f_young": younger_than_first,
        "f_mid": younger_than_second - younger_than_first,
        "f_old": 1 - younger_than_second,
    }
    for bin_index in range(sm.N_BINS):
        labels[f"frac_bin_{bin_index:02d}"] = fractions[:, bin_index]
    for bin_index in range(sm.N_BINS):
        labels[f"logsfr_bin_{bin_index:02d}"] = log_sfr_bins[:, bin_index]
    for bin_index in range(sm.N_BINS):
        labels[f"zmet_bin_{bin_index:02d}"] = metallicities[:, bin_index]
    return pd.DataFrame(labels)


# =============================================================================
# 4. Simulated observations
# =============================================================================
def draw_observed_redshift(true_redshift, rng):
    """Observed redshift and its uncertainty: exact for a P_SPECTROSCOPIC share, photometric otherwise.

    Inputs: true_redshift (n,); rng: numpy Generator.
    Returns: (z_obs, sigma_z), each (n,); sigma_z = sm.SIGMA_Z * (1 + z) for photometric redshifts.
    """
    spectroscopic = rng.random(len(true_redshift)) < P_SPECTROSCOPIC
    sigma_z = np.where(spectroscopic, 0.0, sm.SIGMA_Z * (1 + true_redshift))
    z_obs = true_redshift + sigma_z * rng.normal(0, 1, len(true_redshift))
    return z_obs, sigma_z


def simulate_observations(table, catalogue_r, catalogue_z, noise_table, rng):
    """One fresh noisy realisation of every noise-free model: brightness, S-PLUS noise and redshift.

    Inputs: table: noise-free models (training_set.parquet columns); catalogue_r [mag], catalogue_z: DR3 r_auto
            and zml (brightness donors); noise_table: splus_noise_model.csv; rng: numpy Generator.
    Returns: dict(flat, banded, valid, targets) — inputs as colour_inputs gives them and the (n, 19) targets.
             valid is False where r is undetected or no DR3 galaxy lies at the model's redshift.
    """
    true_redshift = table["redshift"].to_numpy(dtype=float)
    noiseless = table[sm.BANDS].to_numpy(dtype=float)

    # 1. Realistic brightness, and the masses that go with it
    scaled, offsets = sm.scale_to_observed_r(noiseless, true_redshift, catalogue_r, catalogue_z, rng)
    labels = table.copy()
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED
    labels = sm.rescale_mass_labels(labels, offsets)
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * offsets

    # 2. S-PLUS noise and non-detections
    noisy, errors, _ = sm.add_splus_noise(scaled, noise_table, rng)

    # 3. Observed redshift
    z_obs, sigma_z = draw_observed_redshift(true_redshift, rng)

    observations = colour_inputs(noisy, errors, z_obs, sigma_z)
    observations["valid"] = observations["valid"] & np.isfinite(offsets)
    observations["targets"] = encode_targets(labels)
    return observations


def validation_observations(held_out, mode, rng):
    """Inputs and targets of a held-out set whose noise is already fixed (synthetic_validation/test).

    Inputs: held_out: DataFrame with noisy BANDS magnitudes, e_<band> errors, magnitude_offset and labels
            (mass labels already rescaled); mode: "spectroscopic" (sigma_z = 0) or "photometric"
            (sigma_z = sm.SIGMA_Z * (1 + z), z_obs drawn with rng); rng: numpy Generator (seeded by the caller).
    Returns: dict(flat, banded, valid, targets).
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
    observations["targets"] = encode_targets(labels)
    return observations


# =============================================================================
# 5. Standardisation
# =============================================================================
class Standardiser:
    """Column-wise (x - mean) / std, fitted on training data and stored with the model.

    Entries whose mask is False are excluded from the fit and set to 0 after transforming
    (0 = 'no information' for a masked colour).
    """

    def __init__(self, mean, std):
        """Inputs: mean, std: arrays with one value per column (std > 0)."""
        self.mean = np.asarray(mean, dtype=float)
        self.std = np.asarray(std, dtype=float)

    @classmethod
    def fit(cls, values, mask=None):
        """Fit on values (n, ...) using only entries where mask is True; returns a Standardiser."""
        values = np.asarray(values, dtype=float)
        if mask is None:
            mask = np.ones(values.shape, dtype=bool)
        masked = np.where(mask, values, np.nan)
        mean = np.nanmean(masked, axis=0)
        std = np.nanstd(masked, axis=0)
        std = np.where(std > 0, std, 1.0)        # constant columns (e.g. log pivots) stay unscaled
        mean = np.nan_to_num(mean)
        return cls(mean, std)

    def transform(self, values, mask=None):
        """Standardised copy of values; masked entries become 0."""
        transformed = (np.asarray(values, dtype=float) - self.mean) / self.std
        if mask is not None:
            transformed = np.where(mask, transformed, 0.0)
        return transformed

    def inverse_transform(self, values):
        """Undo transform (for unmasked entries)."""
        return np.asarray(values, dtype=float) * self.std + self.mean

    def to_dict(self):
        """JSON-friendly dict of the statistics."""
        return {"mean": self.mean.tolist(), "std": self.std.tolist()}

    @classmethod
    def from_dict(cls, statistics):
        """Standardiser from to_dict output."""
        return cls(statistics["mean"], statistics["std"])


def input_masks(flat, banded):
    """Which entries of flat and banded carry information (for Standardiser fit/transform).

    Inputs: flat (n, 37), banded (n, 12, 4) as colour_inputs returns them.
    Returns: (flat_mask, banded_mask) boolean arrays of the same shapes: colours and colour errors follow
             their flags, flags and scalars are always informative; band channels 0-1 follow channel 2.
    """
    flat_mask = np.ones(flat.shape, dtype=bool)
    flags = flat[:, 2 * N_COLOURS:3 * N_COLOURS] > 0
    flat_mask[:, :N_COLOURS] = flags
    flat_mask[:, N_COLOURS:2 * N_COLOURS] = flags
    banded_mask = np.ones(banded.shape, dtype=bool)
    band_detected = banded[:, :, 2] > 0
    banded_mask[:, :, 0] = band_detected
    banded_mask[:, :, 1] = band_detected
    return flat_mask, banded_mask
```

- [ ] **Step 4: Run the tests**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_data.py -v`
Expected: 10 passed. If `test_decoded_summary_labels_match_bagpipes` fails only for `mw_age` by a small margin, print the per-row relative differences in the report before touching the tolerance: Bagpipes integrates on its own age grid, and the bin-midpoint formula must agree within 2% for constant-SFR bins.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/ml_data.py Synthetic_domain/tests/test_ml_data.py
git commit -m "feat: ML inputs from colours, 19-number targets, decoding and simulated observations"
```

---

### Task 3: `ml_networks.py` — embeddings, flow, training, checkpoints, sampling

**Files:**
- Create: `Synthetic_domain/ml_networks.py`, `Synthetic_domain/tests/test_ml_networks.py`

**Interfaces:**
- Consumes: `ml_data.N_INPUTS`, `ml_data.N_BAND_CHANNELS`, `ml_data.SCALAR_INPUT_SLICE`, `ml_data.N_TARGETS`, `ml_data.Standardiser`, `ml_data.input_masks`.
- Produces:
  - `DEFAULT_CONFIG` (dict);
  - `build_estimator(config) -> PosteriorEstimator`;
  - `PosteriorEstimator.log_prob(theta, flat, banded) -> Tensor (batch,)`;
  - `PosteriorEstimator.sample(flat, banded, n_samples) -> Tensor (batch, n_samples, n_targets)`;
  - `train_estimator(estimator, make_epoch, validation_sets, config, device, log=print) -> list[dict]` (history; restores the best weights);
  - `save_run(directory, estimator, config, history, standardisers)`;
  - `load_run(directory, device="cpu") -> (estimator, config, standardisers)`;
  - `posterior_samples(estimator, standardisers, flat, banded, n_samples, seed, device="cpu", batch_size=512) -> np.ndarray (n, n_samples, n_targets)` in physical target units.
  - Epoch and validation sets are dicts of float32 tensors `{"flat", "banded", "targets"}` that are already standardised. `standardisers` is a dict `{"flat", "banded", "targets"}` of `ml_data.Standardiser`.

- [ ] **Step 1: Write the failing tests**

Create `Synthetic_domain/tests/test_ml_networks.py`:
```python
"""Tests for ml_networks.py (CPU, small). Run from Synthetic_domain/."""

import numpy as np
import pytest
import torch

import ml_data
import ml_networks

NOISE_SIGMA = 0.5        # toy problem: x = theta + N(0, NOISE_SIGMA), theta ~ N(0, 1)


def toy_problem(n_rows, generator):
    """Gaussian toy: two targets observed with noise in the first inputs (flat) and band channel 0."""
    theta = torch.randn(n_rows, 2, generator=generator)
    observed = theta + NOISE_SIGMA * torch.randn(n_rows, 2, generator=generator)
    flat = torch.zeros(n_rows, ml_data.N_INPUTS)
    flat[:, :2] = observed
    banded = torch.zeros(n_rows, 12, ml_data.N_BAND_CHANNELS)
    banded[:, :2, 0] = observed
    return {"flat": flat, "banded": banded, "targets": theta}


def small_config(embedding):
    """A small, fast configuration for tests."""
    config = dict(ml_networks.DEFAULT_CONFIG)
    config.update(embedding=embedding, n_targets=2, summary_size=16, hidden_size=64,
                  flow_transforms=3, flow_hidden=(64, 64), epochs=40, patience=40, batch_size=512,
                  learning_rate=3e-3, seed=0)
    return config


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


def test_save_load_round_trip_on_cpu(tmp_path):
    config = small_config("mlp")
    estimator = ml_networks.build_estimator(config)
    standardisers = {
        "flat": ml_data.Standardiser(np.zeros(ml_data.N_INPUTS), np.ones(ml_data.N_INPUTS)),
        "banded": ml_data.Standardiser(np.zeros(ml_data.N_BAND_CHANNELS), np.ones(ml_data.N_BAND_CHANNELS)),
        "targets": ml_data.Standardiser(np.array([1.0, 2.0]), np.array([0.5, 3.0])),
    }
    history = [{"epoch": 0, "train_loss": 1.0, "validation_loss": 1.1}]
    ml_networks.save_run(str(tmp_path), estimator, config, history, standardisers)
    loaded, loaded_config, loaded_standardisers = ml_networks.load_run(str(tmp_path), device="cpu")

    flat = np.random.default_rng(0).normal(size=(4, ml_data.N_INPUTS))
    banded = np.random.default_rng(1).normal(size=(4, 12, ml_data.N_BAND_CHANNELS))
    before = ml_networks.posterior_samples(estimator, standardisers, flat, banded, 50, seed=3)
    after = ml_networks.posterior_samples(loaded, loaded_standardisers, flat, banded, 50, seed=3)
    assert loaded_config["embedding"] == "mlp"
    assert before.shape == (4, 50, 2)
    assert np.allclose(before, after)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_networks.py -v`
Expected: `ModuleNotFoundError: No module named 'ml_networks'`.

- [ ] **Step 3: Implement `ml_networks.py`**

Create `Synthetic_domain/ml_networks.py`:
```python
"""
Posterior estimator for the SFH targets: an embedding network (MLP or 1-D CNN) summarises the inputs,
and a zuko neural spline flow models p(targets | summary). Training maximises log p(true targets).
"""

import json
import math
import os

import numpy as np
import torch
from torch import nn
import zuko

import ml_data


# =============================================================================
# 1. Configuration
# =============================================================================
DEFAULT_CONFIG = {
    "embedding": "mlp",              # "mlp" or "cnn"
    "n_targets": ml_data.N_TARGETS,
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
MODEL_FILE = "model.pt"
CONFIG_FILE = "config.json"
HISTORY_FILE = "history.csv"


# =============================================================================
# 2. Networks
# =============================================================================
def dense_stack(input_size, hidden_size, n_layers, output_size):
    """n_layers of Linear + SiLU, then a Linear to output_size."""
    layers = []
    current_size = input_size
    for _ in range(n_layers):
        layers.append(nn.Linear(current_size, hidden_size))
        layers.append(nn.SiLU())
        current_size = hidden_size
    layers.append(nn.Linear(current_size, output_size))
    return nn.Sequential(*layers)


class MLPEmbedding(nn.Module):
    """Dense network on the 37 flat inputs (colours, errors, flags, scalars)."""

    def __init__(self, config):
        """Inputs: config with hidden_size, mlp_layers, summary_size."""
        super().__init__()
        self.network = dense_stack(ml_data.N_INPUTS, config["hidden_size"], config["mlp_layers"],
                                   config["summary_size"])

    def forward(self, flat, banded):
        """Inputs: flat (batch, 37); banded is ignored. Returns: summary (batch, summary_size)."""
        return self.network(flat)


class CNNEmbedding(nn.Module):
    """1-D convolutions over the 12 bands ordered by wavelength, then a dense head with the scalars."""

    def __init__(self, config):
        """Inputs: config with cnn_channels, cnn_kernel, cnn_layers, hidden_size, head_layers, summary_size."""
        super().__init__()
        layers = []
        in_channels = ml_data.N_BAND_CHANNELS
        for _ in range(config["cnn_layers"]):
            layers.append(nn.Conv1d(in_channels, config["cnn_channels"], config["cnn_kernel"],
                                    padding=config["cnn_kernel"] // 2))
            layers.append(nn.SiLU())
            in_channels = config["cnn_channels"]
        self.convolutions = nn.Sequential(*layers)
        n_scalars = ml_data.SCALAR_INPUT_SLICE.stop - ml_data.SCALAR_INPUT_SLICE.start
        self.head = dense_stack(config["cnn_channels"] + n_scalars, config["hidden_size"],
                                config["head_layers"], config["summary_size"])

    def forward(self, flat, banded):
        """Inputs: flat (batch, 37) for the scalars; banded (batch, 12, 4). Returns: summary (batch, summary_size)."""
        channels_first = banded.transpose(1, 2)              # (batch, 4, 12): Conv1d wants channels first
        features = self.convolutions(channels_first).mean(dim=2)
        scalars = flat[:, ml_data.SCALAR_INPUT_SLICE]
        return self.head(torch.cat([features, scalars], dim=1))


EMBEDDINGS = {"mlp": MLPEmbedding, "cnn": CNNEmbedding}


class PosteriorEstimator(nn.Module):
    """Embedding + conditional neural spline flow over the standardised targets."""

    def __init__(self, config):
        """Inputs: config (see DEFAULT_CONFIG)."""
        super().__init__()
        self.embedding = EMBEDDINGS[config["embedding"]](config)
        self.flow = zuko.flows.NSF(features=config["n_targets"], context=config["summary_size"],
                                   transforms=config["flow_transforms"],
                                   hidden_features=tuple(config["flow_hidden"]))

    def log_prob(self, theta, flat, banded):
        """log p(theta | inputs) per row. Inputs: tensors (batch, n_targets), (batch, 37), (batch, 12, 4)."""
        summary = self.embedding(flat, banded)
        return self.flow(summary).log_prob(theta)

    def sample(self, flat, banded, n_samples):
        """n_samples posterior draws per row. Returns: tensor (batch, n_samples, n_targets)."""
        summary = self.embedding(flat, banded)
        draws = self.flow(summary).sample((n_samples,))       # (n_samples, batch, n_targets)
        return draws.transpose(0, 1)


def build_estimator(config):
    """A fresh PosteriorEstimator with weights seeded by config['seed']."""
    torch.manual_seed(config["seed"])
    return PosteriorEstimator(config)


# =============================================================================
# 3. Training
# =============================================================================
def mean_loss(estimator, data_set, batch_size, device):
    """Mean negative log-probability over a data set (no gradients)."""
    estimator.eval()
    total = 0.0
    n_rows = len(data_set["targets"])
    with torch.no_grad():
        for start in range(0, n_rows, batch_size):
            stop = start + batch_size
            theta = data_set["targets"][start:stop].to(device)
            flat = data_set["flat"][start:stop].to(device)
            banded = data_set["banded"][start:stop].to(device)
            total += -estimator.log_prob(theta, flat, banded).sum().item()
    return total / n_rows


def train_estimator(estimator, make_epoch, validation_sets, config, device, log=print):
    """Train with AdamW + cosine decay and early stopping; keeps the best weights.

    Inputs: estimator; make_epoch(epoch_index) -> dict of standardised tensors {flat, banded, targets}
            (a fresh noisy realisation each epoch); validation_sets: list of such dicts (fixed);
            config: batch_size, learning_rate, epochs, patience, seed; device: "cpu" or "cuda";
            log: function called with one progress line per epoch.
    Returns: history, a list of dicts {epoch, train_loss, validation_loss, learning_rate}.
    """
    estimator.to(device)
    optimiser = torch.optim.AdamW(estimator.parameters(), lr=config["learning_rate"])
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimiser, T_max=config["epochs"])
    shuffle_generator = torch.Generator().manual_seed(config["seed"])

    best_loss = math.inf
    best_state = None
    epochs_without_improvement = 0
    history = []
    for epoch_index in range(config["epochs"]):
        # 1. One pass over a fresh noisy realisation of the training set
        epoch = make_epoch(epoch_index)
        n_rows = len(epoch["targets"])
        order = torch.randperm(n_rows, generator=shuffle_generator)
        estimator.train()
        running = 0.0
        for start in range(0, n_rows, config["batch_size"]):
            batch = order[start:start + config["batch_size"]]
            theta = epoch["targets"][batch].to(device)
            flat = epoch["flat"][batch].to(device)
            banded = epoch["banded"][batch].to(device)
            loss = -estimator.log_prob(theta, flat, banded).mean()
            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            running += loss.item() * len(batch)
        scheduler.step()

        # 2. Validation loss: mean over the validation sets (e.g. spectroscopic and photometric z)
        validation_losses = []
        for validation_set in validation_sets:
            validation_losses.append(mean_loss(estimator, validation_set, config["batch_size"], device))
        validation_loss = float(np.mean(validation_losses))
        history.append({"epoch": epoch_index, "train_loss": running / n_rows,
                        "validation_loss": validation_loss, "learning_rate": scheduler.get_last_lr()[0]})
        log(f"epoch {epoch_index}: train {running / n_rows:.4f}  validation {validation_loss:.4f}")

        # 3. Early stopping on the validation loss
        if validation_loss < best_loss:
            best_loss = validation_loss
            best_state = {name: value.detach().cpu().clone() for name, value in estimator.state_dict().items()}
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= config["patience"]:
                log(f"early stop after epoch {epoch_index}: best validation {best_loss:.4f}")
                break

    estimator.load_state_dict(best_state)
    return history


# =============================================================================
# 4. Checkpoints and sampling
# =============================================================================
def save_run(directory, estimator, config, history, standardisers):
    """Write model.pt (weights), config.json (config + standardisation statistics) and history.csv."""
    os.makedirs(directory, exist_ok=True)
    torch.save(estimator.state_dict(), os.path.join(directory, MODEL_FILE))
    stored = dict(config)
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


def load_run(directory, device="cpu"):
    """Read a run written by save_run. Returns: (estimator on device, config, standardisers dict)."""
    with open(os.path.join(directory, CONFIG_FILE)) as source:
        stored = json.load(source)
    standardisers = {}
    for name, statistics in stored.pop("standardisers").items():
        standardisers[name] = ml_data.Standardiser.from_dict(statistics)
    estimator = PosteriorEstimator(stored)
    weights = torch.load(os.path.join(directory, MODEL_FILE), map_location=device, weights_only=True)
    estimator.load_state_dict(weights)
    estimator.to(device)
    estimator.eval()
    return estimator, stored, standardisers


def posterior_samples(estimator, standardisers, flat, banded, n_samples, seed, device="cpu", batch_size=512):
    """Posterior draws in physical target units for raw (unstandardised) inputs.

    Inputs: estimator; standardisers {flat, banded, targets}; flat (n, 37) and banded (n, 12, 4) as
            ml_data.colour_inputs returns them; n_samples per galaxy; seed (torch); device; batch_size (galaxies).
    Returns: array (n, n_samples, n_targets).
    """
    flat_mask, banded_mask = ml_data.input_masks(flat, banded)
    flat_standard = standardisers["flat"].transform(flat, flat_mask)
    banded_standard = standardisers["banded"].transform(banded, banded_mask)

    torch.manual_seed(seed)
    estimator.to(device)
    estimator.eval()
    pieces = []
    with torch.no_grad():
        for start in range(0, len(flat), batch_size):
            flat_batch = torch.as_tensor(flat_standard[start:start + batch_size], dtype=torch.float32, device=device)
            banded_batch = torch.as_tensor(banded_standard[start:start + batch_size], dtype=torch.float32,
                                           device=device)
            draws = estimator.sample(flat_batch, banded_batch, n_samples).cpu().numpy()
            pieces.append(standardisers["targets"].inverse_transform(draws))
    return np.concatenate(pieces, axis=0)
```

- [ ] **Step 4: Run the tests**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_networks.py -v`
Expected: 4 passed, in under ~3 min on CPU. If the Gaussian recovery misses by a small margin, raise `epochs` in `small_config` (not the tolerances) and report it.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/ml_networks.py Synthetic_domain/tests/test_ml_networks.py
git commit -m "feat: MLP and CNN embeddings with a zuko spline flow, training loop and checkpoints"
```

---

### Task 4: `ml_evaluation.py` — metrics

**Files:**
- Create: `Synthetic_domain/ml_evaluation.py`, `Synthetic_domain/tests/test_ml_evaluation.py`

**Interfaces:**
- Produces:
  - `coverage_curve(samples, truths, levels) -> np.ndarray` (empirical central-interval coverage per level; samples (n, m), truths (n,));
  - `truth_ranks(samples, truths) -> np.ndarray` (n,) in [0, 1];
  - `bias_and_nmad(estimates, truths) -> (bias, nmad)`;
  - `width_ratio(posterior_samples, prior_samples, level=0.68) -> float` (median over galaxies of posterior width / prior width);
  - `compare_with_spectral(ml_samples, spectral_p16, spectral_median, spectral_p84) -> dict(offset_in_widths (n,), inside_68 (n,) bool, inside_95 (n,) bool)`;
  - `NOMINAL_LEVELS = np.linspace(0.1, 0.9, 9)`; `NMAD_SCALE = 1.4826`.

- [ ] **Step 1: Write the failing tests**

Create `Synthetic_domain/tests/test_ml_evaluation.py`:
```python
"""Tests for ml_evaluation.py. Run from Synthetic_domain/."""

import numpy as np

import ml_evaluation as ev


def test_calibrated_posteriors_give_a_diagonal_coverage_curve():
    rng = np.random.default_rng(40)
    n_galaxies, n_samples = 4000, 400
    centres = rng.normal(0, 1, n_galaxies)
    truths = centres + rng.normal(0, 1, n_galaxies)         # truth drawn from the posterior itself
    samples = centres[:, None] + rng.normal(0, 1, (n_galaxies, n_samples))
    coverage = ev.coverage_curve(samples, truths, ev.NOMINAL_LEVELS)
    assert np.allclose(coverage, ev.NOMINAL_LEVELS, atol=0.03)


def test_truth_ranks_are_uniform_for_calibrated_posteriors():
    rng = np.random.default_rng(41)
    samples = rng.normal(0, 1, (5000, 200))
    truths = rng.normal(0, 1, 5000)
    ranks = ev.truth_ranks(samples, truths)
    histogram, _ = np.histogram(ranks, bins=10, range=(0, 1))
    assert np.all(np.abs(histogram / 5000 - 0.1) < 0.02)


def test_bias_and_nmad():
    truths = np.linspace(0, 1, 101)
    bias, nmad = ev.bias_and_nmad(truths + 0.2, truths)
    assert np.isclose(bias, 0.2)
    assert np.isclose(nmad, 0.0)


def test_width_ratio_is_one_when_posterior_equals_prior():
    rng = np.random.default_rng(42)
    prior = rng.normal(0, 1, 20000)
    posterior = rng.normal(0, 1, (300, 2000))
    assert np.isclose(ev.width_ratio(posterior, prior), 1.0, atol=0.05)
    assert ev.width_ratio(0.1 * posterior, prior) < 0.15


def test_compare_with_spectral():
    rng = np.random.default_rng(43)
    ml_samples = rng.normal(0, 1, (2, 5000))
    result = ev.compare_with_spectral(ml_samples, np.array([-0.1, 2.9]), np.array([0.0, 3.0]), np.array([0.1, 3.1]))
    assert bool(result["inside_68"][0]) and bool(result["inside_95"][0])
    assert not bool(result["inside_95"][1])
    assert abs(result["offset_in_widths"][0]) < 0.05
    assert result["offset_in_widths"][1] < -2.5
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_evaluation.py -v`
Expected: `ModuleNotFoundError: No module named 'ml_evaluation'`.

- [ ] **Step 3: Implement**

Create `Synthetic_domain/ml_evaluation.py`:
```python
"""
Metrics for the SFH posterior estimator: calibration (coverage, ranks), accuracy (bias, NMAD),
information per bin (posterior/prior width) and agreement with the Bagpipes spectral posteriors.
"""

import numpy as np

NOMINAL_LEVELS = np.linspace(0.1, 0.9, 9)   # central-interval probabilities of the coverage curve
NMAD_SCALE = 1.4826                          # makes the median absolute deviation a Gaussian sigma
ONE_SIGMA_LEVEL = 0.68
TWO_SIGMA_LEVEL = 0.95
PERCENT = 100


def central_interval(samples, level):
    """Lower and upper bounds of the central interval holding `level` of the samples (axis -1)."""
    tail = (1 - level) / 2 * PERCENT
    low = np.nanpercentile(samples, tail, axis=-1)
    high = np.nanpercentile(samples, PERCENT - tail, axis=-1)
    return low, high


def coverage_curve(samples, truths, levels=NOMINAL_LEVELS):
    """Empirical fraction of truths inside the central interval, for each nominal level.

    Inputs: samples (n, m) posterior draws of one label; truths (n,); levels: nominal probabilities.
    Returns: array of empirical coverages, one per level (equal to levels for a calibrated estimator).
    """
    coverages = []
    for level in levels:
        low, high = central_interval(samples, level)
        inside = (truths >= low) & (truths <= high)
        coverages.append(np.mean(inside))
    return np.array(coverages)


def truth_ranks(samples, truths):
    """Fraction of each galaxy's samples below its truth (uniform on [0, 1] when calibrated)."""
    return np.mean(samples < truths[:, np.newaxis], axis=1)


def bias_and_nmad(estimates, truths):
    """Median offset and robust scatter (NMAD) of estimates around truths. Returns: (bias, nmad)."""
    residuals = np.asarray(estimates) - np.asarray(truths)
    bias = np.nanmedian(residuals)
    nmad = NMAD_SCALE * np.nanmedian(np.abs(residuals - bias))
    return bias, nmad


def width_ratio(posterior_samples, prior_samples, level=ONE_SIGMA_LEVEL):
    """How much the data shrink a label: median over galaxies of posterior width / prior width.

    Inputs: posterior_samples (n, m); prior_samples (k,) draws of the same label from the prior;
            level: central-interval probability.
    Returns: ratio near 1 = prior-dominated, near 0 = well constrained.
    """
    prior_low, prior_high = central_interval(prior_samples, level)
    posterior_low, posterior_high = central_interval(posterior_samples, level)
    return float(np.nanmedian((posterior_high - posterior_low) / (prior_high - prior_low)))


def compare_with_spectral(ml_samples, spectral_p16, spectral_median, spectral_p84):
    """ML posterior vs the Bagpipes spectral posterior of each reference galaxy.

    Inputs: ml_samples (n, m) of one label; spectral 16th, 50th and 84th percentiles (n,).
    Returns: dict with offset_in_widths (median difference / combined 1-sigma width, ML minus spectral),
             inside_68 and inside_95 (spectral median inside the ML central interval).
    """
    ml_low, ml_high = central_interval(ml_samples, ONE_SIGMA_LEVEL)
    ml_median = np.nanmedian(ml_samples, axis=1)
    ml_sigma = (ml_high - ml_low) / 2
    spectral_sigma = (np.asarray(spectral_p84) - np.asarray(spectral_p16)) / 2
    combined = np.sqrt(ml_sigma**2 + spectral_sigma**2)
    offset = (ml_median - np.asarray(spectral_median)) / combined

    low_95, high_95 = central_interval(ml_samples, TWO_SIGMA_LEVEL)
    inside_68 = (spectral_median >= ml_low) & (spectral_median <= ml_high)
    inside_95 = (spectral_median >= low_95) & (spectral_median <= high_95)
    return {"offset_in_widths": offset, "inside_68": inside_68, "inside_95": inside_95}
```

- [ ] **Step 4: Run the tests**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_ml_evaluation.py -v`
Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/ml_evaluation.py Synthetic_domain/tests/test_ml_evaluation.py
git commit -m "feat: calibration, accuracy, information and spectral-agreement metrics"
```

---

### Task 5: Training CLI and GPU jobs

**Files:**
- Create: `Synthetic_domain/train_ml_model.py`, `Synthetic_domain/tests/test_train_ml_model.py`, `Synthetic_domain/remote/train_job.sh`, `Synthetic_domain/remote/submit_training.sh`
- Modify: `Synthetic_domain/remote/build_env.sh` (optional torch install), `Synthetic_domain/remote/README.md`

**Interfaces:**
- Consumes: `ml_data.simulate_observations`, `ml_data.validation_observations`, `ml_data.Standardiser`, `ml_data.input_masks`, `ml_networks.*`.
- Produces: `python train_ml_model.py --embedding mlp|cnn [--data-dir .] [--output-dir ml_runs] [--epochs N] [--max-train-rows N] [--device auto|cpu|cuda] [--seed S]`, which writes `<output-dir>/<embedding>/{model.pt, config.json, history.csv}`. `remote/submit_training.sh` submits one GPU job that trains both embeddings and returns `ml_runs/` in `out/`.

- [ ] **Step 1: Write the failing smoke test**

Create `Synthetic_domain/tests/test_train_ml_model.py`:
```python
"""Smoke test of the training CLI on a tiny subset (CPU). Run from Synthetic_domain/."""

import os
import subprocess
import sys

import pandas as pd


def test_training_cli_writes_a_run(tmp_path):
    command = [sys.executable, "train_ml_model.py", "--embedding", "mlp", "--output-dir", str(tmp_path),
               "--epochs", "2", "--max-train-rows", "3000", "--max-validation-rows", "500", "--device", "cpu"]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=900)
    assert completed.returncode == 0, completed.stderr[-2000:]
    run_directory = tmp_path / "mlp"
    for name in ("model.pt", "config.json", "history.csv"):
        assert (run_directory / name).exists(), name
    history = pd.read_csv(run_directory / "history.csv")
    assert len(history) == 2
```

- [ ] **Step 2: Run it to verify it fails**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_train_ml_model.py -v`
Expected: FAIL (non-zero return code: `train_ml_model.py` does not exist).

- [ ] **Step 3: Implement `train_ml_model.py`**

Create `Synthetic_domain/train_ml_model.py`:
```python
"""
Train the SFH posterior estimator.

Usage (from Synthetic_domain/):
    python train_ml_model.py --embedding mlp
    python train_ml_model.py --embedding cnn
Each run writes ml_runs/<embedding>/{model.pt, config.json, history.csv}. The training set is noise-free;
every epoch draws a fresh brightness, S-PLUS noise and observed redshift for every model (ml_data).
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
import synthetic_model as sm

TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
NOISE_MODEL_FILE = "splus_noise_model.csv"
CATALOGUE_FILE = "dr3_corrected.csv"
MIN_PROB_GALAXY = 0.8
VALIDATION_SEEDS = {"spectroscopic": 100, "photometric": 101}


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


def main():
    """Parse arguments, build the data pipeline, train one embedding and save the run."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--embedding", choices=sorted(ml_networks.EMBEDDINGS), required=True)
    parser.add_argument("--data-dir", default=".")
    parser.add_argument("--output-dir", default="ml_runs")
    parser.add_argument("--epochs", type=int, default=ml_networks.DEFAULT_CONFIG["epochs"])
    parser.add_argument("--max-train-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--max-validation-rows", type=int, default=None, help="subset for smoke tests")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--seed", type=int, default=ml_networks.DEFAULT_CONFIG["seed"])
    arguments = parser.parse_args()

    device = arguments.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {device}", flush=True)

    # 1. Data
    training = pd.read_parquet(os.path.join(arguments.data_dir, TRAINING_FILE))
    if arguments.max_train_rows is not None:
        training = training.iloc[:arguments.max_train_rows]
    validation = pd.read_parquet(os.path.join(arguments.data_dir, VALIDATION_FILE))
    if arguments.max_validation_rows is not None:
        validation = validation.iloc[:arguments.max_validation_rows]
    noise_table = pd.read_csv(os.path.join(arguments.data_dir, NOISE_MODEL_FILE))
    catalogue = pd.read_csv(os.path.join(arguments.data_dir, CATALOGUE_FILE), usecols=["PROB_GAL", "r_auto", "zml"])
    catalogue = catalogue[catalogue["PROB_GAL"] > MIN_PROB_GALAXY]
    catalogue_r = catalogue["r_auto"].to_numpy()
    catalogue_z = catalogue["zml"].to_numpy()

    # 2. Standardisation statistics from one simulated epoch (seeded)
    first_epoch = to_tensors(ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table,
                                                           np.random.default_rng(arguments.seed)))
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
        observations = ml_data.validation_observations(validation, mode, np.random.default_rng(seed))
        validation_sets.append(standardise(to_tensors(observations), standardisers))

    def make_epoch(epoch_index):
        """Fresh noisy realisation for this epoch (seeded by run seed and epoch)."""
        rng = np.random.default_rng([arguments.seed, epoch_index])
        observations = ml_data.simulate_observations(training, catalogue_r, catalogue_z, noise_table, rng)
        return standardise(to_tensors(observations), standardisers)

    # 4. Train and save
    config = dict(ml_networks.DEFAULT_CONFIG)
    config.update(embedding=arguments.embedding, epochs=arguments.epochs, seed=arguments.seed,
                  training_commit=git_commit(), torch_version=torch.__version__,
                  n_training_models=len(training), device=device)
    estimator = ml_networks.build_estimator(config)
    start = time.time()
    history = ml_networks.train_estimator(estimator, make_epoch, validation_sets, config, device,
                                          log=lambda message: print(message, flush=True))
    config["training_minutes"] = (time.time() - start) / 60
    output_directory = os.path.join(arguments.output_dir, arguments.embedding)
    ml_networks.save_run(output_directory, estimator, config, history, standardisers)
    print(f"saved {output_directory} after {config['training_minutes']:.1f} min", flush=True)


if __name__ == "__main__":
    main()
```

Note for the implementer: the `banded` standardiser is fitted per channel over all bands (reshape to `(-1, 4)`). `Standardiser.transform` on a `(n, 12, 4)` array then broadcasts its 4-element mean and std over the last axis, which is what the CNN needs.

- [ ] **Step 4: Run the smoke test and the whole suite**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_train_ml_model.py -v`, then the full suite.
Expected: PASS, and all tests pass.

- [ ] **Step 5: GPU job scripts**

Modify `remote/build_env.sh`: when the environment variable `INSTALL_TORCH=1` is set, after the pinned packages run
```bash
"$VENV_PYTHON" -m pip install "torch==2.5.1" --index-url https://download.pytorch.org/whl/cu121
"$VENV_PYTHON" -m pip install "zuko==<the version pinned in requirements.txt>"
"$VENV_PYTHON" -c "import torch; print('cuda available:', torch.cuda.is_available())"
```
Use the variable and helper names that `build_env.sh` already defines for the venv Python and logging.

Create `remote/train_job.sh` (runs inside the job; follow the structure and the shared helpers of `fit_job.sh`/`job_common.sh`):
1. unpack the bundle;
2. `INSTALL_TORCH=1` run `build_env.sh`;
3. for each embedding in `${EMBEDDINGS:-mlp cnn}`, run `python train_ml_model.py --embedding "$embedding" --output-dir out/ml_runs` from the unpacked `Synthetic_domain` folder, appending to `out/train_<embedding>.log`;
4. copy `out/ml_runs` at the end (it is written directly to `out/`).

Create `remote/submit_training.sh` (runs locally, modelled on `submit_fits.sh`):
1. bundle `synthetic_model.py`, `ml_data.py`, `ml_networks.py`, `train_ml_model.py`, `splus-filters/`, `training_set.parquet`, `synthetic_validation.parquet`, `splus_noise_model.csv`, `dr3_corrected.csv` (only the `PROB_GAL`, `r_auto`, `zml` columns, written as a small CSV under the same name), `build_env.sh`, `job_common.sh`, `grid_md5.txt` and the MultiNest source (`build_env.sh` builds it);
2. submit with `hljob run --detach --gpu -c 8 -m 7000 -t 08:00:00 -e EMBEDDINGS="mlp cnn" -f bundle.tar.gz -f train_job.sh -- bash inputs/train_job.sh`;
3. print the job id and the logs / wait / fetch commands.

Add a "Training" section to `remote/README.md` (how to submit, follow and fetch, and where the runs land locally: `Synthetic_domain/ml_runs/<embedding>/`).

- [ ] **Step 6: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/train_ml_model.py Synthetic_domain/tests/test_train_ml_model.py Synthetic_domain/remote/
git commit -m "feat: training CLI and GPU training jobs for both embeddings"
```

- [ ] **Step 7 (controller): run the GPU training**

The controller runs `remote/submit_training.sh` once the job server is free (spectral fits have priority), fetches the job, and copies `out/ml_runs/*` into `Synthetic_domain/ml_runs/`. Then it commits `ml_runs/*/config.json` and `ml_runs/*/history.csv` (weights stay ignored).

---

### Task 6: `ml_model.ipynb` — model choice, evaluation and DR3 catalogue

**Files:**
- Create: `Synthetic_domain/ml_model.ipynb`

**Interfaces:**
- Consumes: `ml_runs/{mlp,cnn}/` (Task 5), `ml_data`, `ml_networks`, `ml_evaluation`, `synthetic_model`, `synthetic_validation.parquet`, `synthetic_test.parquet`, `training_set.parquet`, `spectroscopic_reference.csv` (read with `dtype={"specobjid": str}`), `dr3_corrected.csv`.
- Produces: figures inline, `ml_evaluation_summary.csv` (committed), `ml_sfh_catalogue.parquet` (ignored).

Build the notebook cell by cell (NotebookEdit or nbformat). Set the kernel to `inference4etgs`. Every cell follows the code-style rules.

- [ ] **Step 1: Cell 0 (markdown) and Cell 1 (constants)**

Cell 0:
```markdown
# SFH posterior from S-PLUS colours (part 3)

A normalizing flow (trained by `train_ml_model.py` on the 3090) gives each galaxy a posterior over its SFH: 15 bin mass fractions, metallicity (Z_fin, hence every bin), dust, masses and redshift. Inputs are the 11 colours against r, the r magnitude (only for the absolute masses), and z_obs with its uncertainty (0 for spectroscopic redshifts).

**Inputs:** `ml_runs/mlp/`, `ml_runs/cnn/`, `synthetic_validation.parquet`, `synthetic_test.parquet`, `training_set.parquet` (prior draws), `spectroscopic_reference.csv`, `dr3_corrected.csv`.
**Outputs:** `ml_evaluation_summary.csv`, `ml_sfh_catalogue.parquet`.
**Run:** top to bottom with the *Inference4ETGs (MultiNest)* kernel. Sections marked **LOCKED** run only with `RUN_LOCKED_TEST = True`, once, after the embedding is chosen: they use the synthetic test set and the reference set's test half.
```

Cell 1:
```python
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

import ml_data
import ml_evaluation as ev
import ml_networks
import synthetic_model as sm

RUN_LOCKED_TEST = False            # set True once, for the final report
RUN_DIRECTORIES = {"mlp": "ml_runs/mlp", "cnn": "ml_runs/cnn"}
VALIDATION_FILE = "synthetic_validation.parquet"
TEST_FILE = "synthetic_test.parquet"
TRAINING_FILE = "training_set.parquet"
REFERENCE_FILE = "spectroscopic_reference.csv"
SPLUS_FILE = "dr3_corrected.csv"
SUMMARY_FILE = "ml_evaluation_summary.csv"
CATALOGUE_FILE = "ml_sfh_catalogue.parquet"

N_POSTERIOR_SAMPLES = 1000         # per galaxy for evaluation
N_CATALOGUE_SAMPLES = 500          # per galaxy for the DR3 catalogue
SAMPLING_SEED = 7
VALIDATION_SEEDS = {"spectroscopic": 100, "photometric": 101}   # same copies as training used
TEST_SEEDS = {"spectroscopic": 200, "photometric": 201}
KEY_LABELS = ["mw_age", "Z_fin", "f_young", "f_mid", "f_old", "log_mstar", "log_sfr"]
REFERENCE_LABELS = ["mw_age", "Z_fin", "f_old", "log_mstar"]
PRIOR_WINDOW_DZ = 0.05             # prior draws within this redshift of a galaxy define its prior width
MIN_PROB_GALAXY = 0.8
N_JOBS = 11

runs = {}
for name, directory in RUN_DIRECTORIES.items():
    runs[name] = ml_networks.load_run(directory, device="cpu")
```

- [ ] **Step 2: Cells 2–3: training histories and model choice (validation only)**

Cell 2 (markdown): `## 1. Training and model choice` — "Both embeddings saw the same data, budget and seeds. The choice uses the validation set only: lower loss, then coverage closer to nominal."

Cell 3:
```python
def labelled_posteriors(run, observations, n_samples, seed):
    """Decoded posterior labels for every valid galaxy: dict label -> array (n_valid, n_samples)."""
    estimator, config, standardisers = run
    valid = observations["valid"]
    theta = ml_networks.posterior_samples(estimator, standardisers, observations["flat"][valid],
                                          observations["banded"][valid], n_samples, seed)
    n_galaxies = theta.shape[0]
    decoded = ml_data.decode_samples(theta.reshape(-1, ml_data.N_TARGETS))
    posteriors = {}
    for column in decoded.columns:
        posteriors[column] = decoded[column].to_numpy().reshape(n_galaxies, n_samples)
    truths = ml_data.decode_samples(observations["targets"][valid])
    return posteriors, truths


validation_table = pd.read_parquet(VALIDATION_FILE)
figure, axes = plt.subplots(1, 2, figsize=(12, 4))
choice_rows = []
for name, run in runs.items():
    history = pd.read_csv(f"{RUN_DIRECTORIES[name]}/{ml_networks.HISTORY_FILE}")
    axes[0].plot(history["epoch"], history["validation_loss"], label=f"{name} validation")
    axes[0].plot(history["epoch"], history["train_loss"], ls=":", label=f"{name} train")
    for mode, seed in VALIDATION_SEEDS.items():
        observations = ml_data.validation_observations(validation_table, mode, np.random.default_rng(seed))
        posteriors, truths = labelled_posteriors(run, observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED)
        for label in KEY_LABELS:
            coverage = ev.coverage_curve(posteriors[label], truths[label].to_numpy(), [0.68, 0.95])
            choice_rows.append({"embedding": name, "mode": mode, "label": label,
                                "coverage_68": coverage[0], "coverage_95": coverage[1],
                                "best_validation_loss": history["validation_loss"].min()})
axes[0].set(xlabel="epoch", ylabel="mean −log p", title="training histories")
axes[0].legend()
choice = pd.DataFrame(choice_rows)
choice["coverage_error"] = np.abs(choice["coverage_68"] - 0.68) + np.abs(choice["coverage_95"] - 0.95)
by_embedding = choice.groupby("embedding").agg(best_validation_loss=("best_validation_loss", "first"),
                                               mean_coverage_error=("coverage_error", "mean"))
display(by_embedding.round(4))
CHOSEN = by_embedding.sort_values(["best_validation_loss", "mean_coverage_error"]).index[0]
print(f"Chosen embedding: {CHOSEN}")
axes[1].axis("off")
plt.show()
```

- [ ] **Step 3: Cells 4–5: calibration, accuracy and information (LOCKED test, both redshift modes)**

Cell 4 (markdown): `## 2. Calibration, accuracy and information per bin (LOCKED: synthetic test set)` — one or two sentences on what each plot shows.

Cell 5:
```python
assert RUN_LOCKED_TEST, "Locked section: set RUN_LOCKED_TEST = True only for the final report"
test_table = pd.read_parquet(TEST_FILE)
prior_table = pd.read_parquet(TRAINING_FILE)
summary_rows = []
bin_labels = []
for bin_index in range(ml_data.N_FRACTION_BINS):
    bin_labels.append(f"frac_bin_{bin_index:02d}")
for mode, seed in TEST_SEEDS.items():
    observations = ml_data.validation_observations(test_table, mode, np.random.default_rng(seed))
    posteriors, truths = labelled_posteriors(runs[CHOSEN], observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED)
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    for label in KEY_LABELS + bin_labels:
        truth = truths[label].to_numpy()
        coverage = ev.coverage_curve(posteriors[label], truth)
        bias, nmad = ev.bias_and_nmad(np.nanmedian(posteriors[label], axis=1), truth)
        prior_draws = prior_table[label].to_numpy()
        ratio = ev.width_ratio(posteriors[label], prior_draws)
        summary_rows.append({"section": "synthetic_test", "mode": mode, "label": label,
                             "coverage_68": ev.coverage_curve(posteriors[label], truth, [0.68])[0],
                             "coverage_95": ev.coverage_curve(posteriors[label], truth, [0.95])[0],
                             "bias": bias, "nmad": nmad, "width_ratio": ratio})
        if label in KEY_LABELS:
            axes[0].plot(ev.NOMINAL_LEVELS, coverage, marker="o", label=label)
            axes[1].hist(ev.truth_ranks(posteriors[label], truth), bins=20, histtype="step", label=label)
    axes[0].plot([0, 1], [0, 1], color="black", lw=1)
    axes[0].set(xlabel="nominal coverage", ylabel="empirical coverage", title=f"{mode} redshift")
    axes[1].set(xlabel="rank of truth", ylabel="galaxies")
    axes[0].legend(fontsize=8)
    plt.show()
summary = pd.DataFrame(summary_rows)
display(summary.round(3))
```
Also add, in the same cell, a plot of `width_ratio` against bin index for `frac_bin_*` (spectroscopic mode) titled "information per age bin (1 = prior only)". Add bias / NMAD against r magnitude in 4 r-magnitude quartile bins for `mw_age` and `log_mstar`, as a small table.

- [ ] **Step 4: Cells 6–7: agreement with the Bagpipes spectral posteriors and the Firefly cross-check**

Cell 6 (markdown): `## 3. Agreement with the Bagpipes spectral posteriors` — "ML posteriors from the 3″ aperture colours and spectroscopic redshift, against each reference galaxy's spectral posterior. Development uses the `validation` half; the `test` half is LOCKED."

Cell 7:
```python
reference = pd.read_csv(REFERENCE_FILE, dtype={"specobjid": str})
halves = ["validation", "test"] if RUN_LOCKED_TEST else ["validation"]
aperture_columns = []
aperture_error_columns = []
for band in sm.BANDS:
    aperture_columns.append(f"{band}_aper")
    aperture_error_columns.append(f"e_{band}_aper")

agreement_rows = []
for half in halves:
    galaxies = reference[(reference["split"] == half) & (reference["fit_status"] == "fitted")].reset_index(drop=True)
    if len(galaxies) == 0:
        print(f"{half}: no fitted galaxies yet")
        continue
    observations = ml_data.colour_inputs(galaxies[aperture_columns].to_numpy(),
                                         galaxies[aperture_error_columns].to_numpy(),
                                         galaxies["z"].to_numpy(), np.zeros(len(galaxies)))
    observations["targets"] = np.zeros((len(galaxies), ml_data.N_TARGETS))   # unknown; not used
    posteriors, _ = labelled_posteriors(runs[CHOSEN], observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED)
    kept = galaxies[observations["valid"]].reset_index(drop=True)
    for label in REFERENCE_LABELS:
        comparison = ev.compare_with_spectral(posteriors[label], kept[f"{label}_p16"], kept[f"{label}_median"],
                                              kept[f"{label}_p84"])
        for galaxy_index in range(len(kept)):
            agreement_rows.append({"half": half, "label": label, "specobjid": kept.loc[galaxy_index, "specobjid"],
                                   "galaxy_type": "star-forming" if kept.loc[galaxy_index, "log_ssfr_median"] > -11 else "passive",
                                   "fibre": kept.loc[galaxy_index, "fiber_diameter"],
                                   "calibrated_fitter": kept.loc[galaxy_index, "fitter_calibrated"],
                                   "offset_in_widths": comparison["offset_in_widths"][galaxy_index],
                                   "inside_68": comparison["inside_68"][galaxy_index],
                                   "inside_95": comparison["inside_95"][galaxy_index]})
agreement = pd.DataFrame(agreement_rows)
if len(agreement) > 0:
    display(agreement.groupby(["half", "label"]).agg(median_offset=("offset_in_widths", "median"),
                                                      inside_68=("inside_68", "mean"),
                                                      inside_95=("inside_95", "mean"),
                                                      galaxies=("specobjid", "count")).round(3))
    display(agreement.groupby(["half", "label", "galaxy_type"])["inside_95"].mean().unstack().round(3))
```
Then add Cell 7b (Firefly cross-check of the ML posteriors):
```python
firefly_rows = []
for half in halves:
    galaxies = reference[(reference["split"] == half) & (reference["fit_status"] == "fitted")].reset_index(drop=True)
    if len(galaxies) == 0:
        continue
    observations = ml_data.colour_inputs(galaxies[aperture_columns].to_numpy(),
                                         galaxies[aperture_error_columns].to_numpy(),
                                         galaxies["z"].to_numpy(), np.zeros(len(galaxies)))
    observations["targets"] = np.zeros((len(galaxies), ml_data.N_TARGETS))
    posteriors, _ = labelled_posteriors(runs[CHOSEN], observations, N_POSTERIOR_SAMPLES, SAMPLING_SEED)
    kept = galaxies[observations["valid"]].reset_index(drop=True)

    # Mass-weighted metallicity per sample: sum(f * Z) / sum(f) over the existing bins
    weighted_sum = np.zeros_like(posteriors["frac_bin_00"])
    fraction_sum = np.zeros_like(posteriors["frac_bin_00"])
    for bin_index in range(ml_data.N_FRACTION_BINS):
        fractions = posteriors[f"frac_bin_{bin_index:02d}"]
        metallicities = np.nan_to_num(posteriors[f"zmet_bin_{bin_index:02d}"])
        weighted_sum += fractions * metallicities
        fraction_sum += fractions
    mass_weighted_metallicity = np.nanmedian(weighted_sum / fraction_sum, axis=1)
    mass_weighted_age = np.nanmedian(posteriors["mw_age"], axis=1)

    comparison = sm.compare_with_firefly(mass_weighted_age, mass_weighted_metallicity,
                                         kept["firefly_age"].to_numpy(), kept["firefly_metallicity"].to_numpy())
    firefly_rows.append({"half": half, "galaxies": len(kept),
                         "median_delta_log_age": np.nanmedian(comparison["delta_log_age"]),
                         "median_delta_log_metallicity": np.nanmedian(comparison["delta_log_metallicity"]),
                         "consistent_fraction": comparison["firefly_consistent"].mean()})
display(pd.DataFrame(firefly_rows).round(3))
```

- [ ] **Step 5: Cells 8–9: DR3 catalogue and speed**

Cell 8 (markdown): `## 4. S-PLUS DR3 catalogue` — "Photometric-z posteriors for every DR3 galaxy that passes the coverage test; the others are flagged."

Cell 9:
```python
splus = pd.read_csv(SPLUS_FILE)
splus = splus[splus["PROB_GAL"] > MIN_PROB_GALAXY].reset_index(drop=True)
magnitudes = splus[[f"{band}_auto" for band in sm.BANDS]].to_numpy()
errors = splus[[f"e_{band}_auto" for band in sm.BANDS]].to_numpy()
detected = sm.detected_bands(magnitudes, errors)
zml = splus["zml"].to_numpy()
sigma_z = sm.SIGMA_Z * (1 + zml)

# 1. Coverage test (colour-based through its free offset), in parallel blocks
training_models = pd.read_parquet(TRAINING_FILE, columns=["redshift", *sm.BANDS])
blocks = np.array_split(np.arange(len(splus)), 4 * N_JOBS)
tasks = []
for block in blocks:
    tasks.append(delayed(sm.coverage_test)(magnitudes[block], errors[block], detected[block], zml[block],
                                           sm.redshift_window(zml[block]), training_models[sm.BANDS].to_numpy(),
                                           training_models["redshift"].to_numpy()))
coverage = pd.concat(Parallel(n_jobs=N_JOBS)(tasks), ignore_index=True)

# 2. Posterior summaries for covered galaxies with a detected r band
inputs = ml_data.colour_inputs(magnitudes, errors, zml, sigma_z)
predictable = inputs["valid"] & coverage["covered"].to_numpy()
estimator, config, standardisers = runs[CHOSEN]
CATALOGUE_CHUNK = 2000             # galaxies per chunk: all samples at once would need ~4 GB
SUMMARY_PERCENTILES = ((16, "p16"), (50, "median"), (84, "p84"))
predictable_indices = np.flatnonzero(predictable)
summaries = {}
start = time.time()
for chunk_start in range(0, len(predictable_indices), CATALOGUE_CHUNK):
    chunk = predictable_indices[chunk_start:chunk_start + CATALOGUE_CHUNK]
    theta = ml_networks.posterior_samples(estimator, standardisers, inputs["flat"][chunk], inputs["banded"][chunk],
                                          N_CATALOGUE_SAMPLES, SAMPLING_SEED + chunk_start)
    decoded = ml_data.decode_samples(theta.reshape(-1, ml_data.N_TARGETS))
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
for name, values in summaries.items():
    catalogue[name] = values
n_predicted = len(predictable_indices)
catalogue.to_parquet(CATALOGUE_FILE, index=False)
print(f"{n_predicted} of {len(splus)} galaxies predicted ({catalogue['out_of_domain'].mean():.1%} out of domain, "
      f"{catalogue['no_r_band'].sum()} without r); {n_predicted / seconds:.0f} galaxies/s on CPU")
```

- [ ] **Step 6: Cell 10: save the summary**

```python
summary_tables = [summary] if RUN_LOCKED_TEST else []
summary_tables.append(choice.assign(section="validation_choice"))
pd.concat(summary_tables, ignore_index=True).to_csv(SUMMARY_FILE, index=False)
print(f"wrote {SUMMARY_FILE}")
```

- [ ] **Step 7: Validate on scratch runs, then commit**

Before real runs exist, validate the notebook with two tiny CPU runs:
- `python train_ml_model.py --embedding mlp --epochs 2 --max-train-rows 3000 --max-validation-rows 500 --device cpu --output-dir _scratch_runs`, and the same for `cnn`;
- execute a scratch copy of the notebook with `RUN_DIRECTORIES` pointed at `_scratch_runs`, `RUN_LOCKED_TEST = True`, and the DR3 section limited to its first 2,000 galaxies (scratch copy only);
- confirm every cell runs and report the printed lines;
- delete the scratch copy, `_scratch_runs` and any scratch outputs.

Commit the notebook with outputs cleared:
```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/ml_model.ipynb
git commit -m "feat: part 3 notebook — model choice, calibration, spectral agreement and DR3 catalogue"
```

- [ ] **Step 8 (controller): final runs**

After Task 5 Step 7 provides real runs:
1. execute the notebook with `RUN_LOCKED_TEST = False` (choice on validation) and check the choice;
2. execute it once with `RUN_LOCKED_TEST = True` and commit `ml_evaluation_summary.csv` plus the executed notebook.

The spectroscopic reference comparison uses whatever fits exist. Run the final report after the spectral-fit job and notebook 2's analysis have produced `spectroscopic_reference.csv` with the fitted galaxies.
