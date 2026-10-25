# Bagpipes Validation Rulers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Validate the Bagpipes spectral fits of the reference galaxies against rulers independent of full-spectrum SED fitting — spectral indices (posterior predictive), Hα SFR and Gallazzi et al. (2005) ages/metallicities, and literature scaling relations — and remove Firefly as a validation.

**Architecture:** A new focused module `Synthetic_domain/spectral_validation.py` holds the index definitions and measurement, the posterior-predictive model indices, light-weighted ages/metallicities, the Speagle main sequence, relation bands and the MPA-JHU SQL builder (all unit-tested). Notebook 2 (`spectroscopic_reference.ipynb`) fetches the literature data, computes the three layers for every fitted galaxy and writes per-galaxy columns, `validation_summary.csv` and figures; its Firefly comparison is removed (Firefly stays only as the selection catalogue).

**Tech Stack:** Python 3.10 venv (`.venv`), numpy, pandas, astropy, Bagpipes 1.3.5 (`bc03_miles`), astroquery (SDSS SkyServer SQL), pytest.

**Spec:** `docs/superpowers/specs/2026-10-23-training-and-reference-sets-design.md` — section "Amendment (2026-10-25): validating the Bagpipes spectral fits against independent rulers; Firefly dropped as a validation".

## Global Constraints

- Code style: the ten rules of the parts 1–2 spec (English names; numbered banners; why-comments; no compact code; docstrings with purpose, inputs with units, returns; named UPPER_CASE constants with units; notebook opening markdown; logic in modules, orchestration in notebooks; parallel cases in one pattern; concise).
- Index passbands (air wavelengths, rest frame, Å), verbatim from the spec: D4000ₙ blue 3850–3950, red 4000–4100 (ratio of mean F_ν, red/blue); Hδ_A feature 4083.500–4122.250, blue 4041.600–4079.750, red 4128.500–4161.000; Mgb 5160.125–5192.625 / 5142.625–5161.375 / 5191.375–5206.375; Fe5270 5245.650–5285.650 / 5233.150–5248.150 / 5285.650–5318.150; Fe5335 5312.125–5352.125 / 5304.625–5315.875 / 5353.375–5363.375; [MgFe]′ = √(Mgb · (0.72 Fe5270 + 0.28 Fe5335)); equivalent widths against a straight pseudo-continuum through the blue and red band means.
- Observed indices: native SDSS resolution, air wavelengths, Milky-Way dereddened, rest frame; errors by 200 Monte-Carlo noise draws from ivar; any ivar = 0 pixel inside an index's bands → that index is NaN for that galaxy.
- Model indices: 100 equal-weight posterior samples per galaxy, Bagpipes model spectrum on the galaxy's native observed grid (same components incl. `veldisp` and nebular emission) times the calibration polynomial fitted to that galaxy's data; same index code.
- Statistics: z = (observed − model median) / √(σ_obs² + σ_model²); two-sided posterior-predictive p; report fraction |z| < 2, median z and z scatter by star-forming/passive and SDSS/BOSS.
- Layer 2: MPA-JHU `galSpecIndx` (code cross-check), `galSpecExtra` `sfr_fib_p50` for BPT star-forming galaxies only vs Bagpipes `logsfr_bin_00` (aperture-normalised), no IMF offset (all sources Kroupa; ruling Q6); Gallazzi et al. (2005) r-band light-weighted log age and log Z/Z☉ vs Bagpipes light-weighted values from per-bin r-band luminosities. BOSS galaxies: "not available".
- Layer 3: stellar mass–metallicity (Gallazzi catalogue + its masses), Speagle et al. (2014) log SFR = (0.84 − 0.026 t) log M* − (6.51 − 0.11 t), ±0.2 dex, t = age of the Universe [Gyr]; Gallazzi age vs SDSS `velDisp` for passive galaxies; fraction of our galaxies inside each 16–84th band.
- Report, never filter: the validation/test split stays "converged and covered". Firefly is removed from validation (columns, figure, flags, check functions) and stays only as the selection catalogue.
- If the Gallazzi catalogue cannot be obtained, its comparisons are reported as unavailable — never replaced by an SED-fitting catalogue.
- Working directory for commands: `/home/pedroand/Projects/Inference4ETGs/Synthetic_domain`; Python `../.venv/bin/python`; tests `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest`.
- Do not touch `spectral_fits/` (a remote job writes fits there), `ml_runs/`, or files another task owns.

## Review Focus

1. **An index band falling outside the observed spectrum or on masked pixels** (high-z galaxy whose Mgb region is redshifted beyond 9200 Å, or a sky-line gap) → that index is NaN for that galaxy (observed and model alike) and is left out of the fractions, never a crash or a zero. Test: Task 1, `test_index_outside_coverage_is_nan`.
2. **A galaxy whose model or observed [MgFe]′ argument is negative** (weak, noisy Fe/Mg lines) → [MgFe]′ is NaN, not a complex number or an exception. Test: Task 1, `test_mgfe_prime_negative_argument_is_nan`.
3. **A galaxy with no MPA-JHU row** (BOSS spectra, or SDSS spectra missing from galSpec) → its layer-2 columns are NaN and labelled "not available"; the summary counts only galaxies with data. Test: Task 4, notebook check cell asserting BOSS rows have NaN `mpa_sfr_fib` and are excluded from the Hα summary count.
4. **Model spectrum with all-zero flux in an index band** (e.g. a degenerate posterior sample) → the index for that sample is NaN and the sample is skipped in the median/std, never inf. Test: Task 2, `test_predictive_indices_skip_degenerate_samples`.
5. **A galaxy in a relation's x-range gap** (mass outside the binned relation) → `inside_band` returns NA for it and it is excluded from the "inside" fraction denominator. Test: Task 2, `test_inside_band_outside_range_is_na`.

---

## File Structure

| File | Responsibility |
|---|---|
| `Synthetic_domain/spectral_validation.py` (create) | Index definitions/measurement, rest-frame spectrum prep, posterior-predictive indices and z-scores, light-weighted age/metallicity, Speagle MS, relation bands, MPA-JHU SQL. |
| `Synthetic_domain/tests/test_spectral_validation.py` (create) | Tests of the module. |
| `Synthetic_domain/literature/` (create) | `mpa_jhu_reference.csv` (per reference galaxy, committed), `relations/*.csv` (binned relation bands, committed), `raw/` (downloaded catalogues, git-ignored), `SOURCES.md` (where each file came from). |
| `Synthetic_domain/spectroscopic_reference.ipynb` (modify) | Remove Firefly validation; add literature fetch and the three validation layers. |
| `Synthetic_domain/ml_model.ipynb` (modify, Task 5) | Remove the Firefly cross-check cell. |
| `Synthetic_domain/synthetic_model.py`, `tests/test_synthetic_model.py` (modify, Task 5) | Delete `compare_with_firefly`, `firefly_scatter_consistency` and their tests (keep `FIREFLY_MISSING`, used by the selection). |
| `.gitignore` (modify) | `Synthetic_domain/literature/raw/`. |

---

### Task 1: Spectral indices — definitions and measurement

**Files:**
- Create: `Synthetic_domain/spectral_validation.py` (sections 1–2), `Synthetic_domain/tests/test_spectral_validation.py`

**Interfaces:**
- Consumes: `sm.vacuum_to_air`, `sm.deredden_spectrum`.
- Produces:
  - `INDEX_NAMES = ["D4000_n", "Hdelta_A", "Mgb", "Fe5270", "Fe5335", "MgFe_prime"]`, `EQUIVALENT_WIDTH_INDICES` (dict name → {"blue","feature","red"} tuples), `D4000_BLUE_BAND`, `D4000_RED_BAND`, `N_NOISE_DRAWS = 200`.
  - `index_values(wavelength, flux) -> dict[str, float]` (rest-frame air Å, any flux units; NaN for an index whose bands hold no pixel).
  - `index_validity(wavelength, ivar) -> dict[str, bool]` (True when every band of the index has pixels and all of them have ivar > 0).
  - `measure_indices(wavelength, flux, ivar, rng, n_draws=N_NOISE_DRAWS) -> dict` with keys `<name>` and `<name>_err` for every name (NaN where invalid).
  - `rest_frame_spectrum(sdss_fits_path, redshift, ebv) -> (wavelength_rest_air, flux, ivar)`.

- [ ] **Step 1: Write the failing tests**

Create `Synthetic_domain/tests/test_spectral_validation.py`:
```python
"""Tests for spectral_validation.py. Run from Synthetic_domain/."""

import numpy as np
import pytest

import spectral_validation as sv

GRID = np.arange(3700.0, 5600.0, 1.0)      # rest-frame air Angstrom, 1 A pixels


def flat_spectrum(level=1.0):
    """Constant F_lambda on GRID."""
    return np.full(GRID.shape, level)


def test_flat_spectrum_has_zero_equivalent_widths():
    values = sv.index_values(GRID, flat_spectrum())
    for name in ("Hdelta_A", "Mgb", "Fe5270", "Fe5335"):
        assert abs(values[name]) < 1e-9, name


def test_d4000_of_flat_f_lambda_is_the_lambda_squared_ratio():
    values = sv.index_values(GRID, flat_spectrum())
    red = (GRID >= 4000) & (GRID <= 4100)
    blue = (GRID >= 3850) & (GRID <= 3950)
    expected = np.mean(GRID[red] ** 2) / np.mean(GRID[blue] ** 2)
    assert np.isclose(values["D4000_n"], expected)


def test_gaussian_absorption_line_equivalent_width():
    depth, sigma, centre = 0.3, 3.0, 5176.0               # inside the Mgb feature band
    flux = 1 - depth * np.exp(-0.5 * ((GRID - centre) / sigma) ** 2)
    expected = depth * sigma * np.sqrt(2 * np.pi)
    assert np.isclose(sv.index_values(GRID, flux)["Mgb"], expected, rtol=0.01)


def test_mgfe_prime_negative_argument_is_nan():
    assert np.isnan(sv.mgfe_prime(-1.0, 2.0, 2.0))
    assert np.isclose(sv.mgfe_prime(4.0, 2.0, 2.0), np.sqrt(4.0 * 2.0))


def test_index_outside_coverage_is_nan():
    short_grid = np.arange(3700.0, 4500.0, 1.0)         # stops before the Lick Mg/Fe region
    values = sv.index_values(short_grid, np.ones_like(short_grid))
    assert np.isnan(values["Mgb"]) and np.isnan(values["MgFe_prime"])
    assert np.isfinite(values["D4000_n"])


def test_masked_pixel_in_band_invalidates_only_that_index():
    ivar = np.ones_like(GRID)
    ivar[(GRID > 5170) & (GRID < 5172)] = 0              # inside the Mgb feature
    measured = sv.measure_indices(GRID, flat_spectrum(), ivar, np.random.default_rng(0), n_draws=20)
    assert np.isnan(measured["Mgb"]) and np.isnan(measured["Mgb_err"]) and np.isnan(measured["MgFe_prime"])
    assert np.isfinite(measured["Hdelta_A"]) and np.isfinite(measured["D4000_n"])


def test_measurement_error_shrinks_with_signal_to_noise():
    rng = np.random.default_rng(1)
    noisy = sv.measure_indices(GRID, flat_spectrum(), np.full(GRID.shape, 100.0), rng, n_draws=200)
    precise = sv.measure_indices(GRID, flat_spectrum(), np.full(GRID.shape, 10000.0), rng, n_draws=200)
    assert precise["Hdelta_A_err"] < noisy["Hdelta_A_err"] / 5
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_spectral_validation.py -v`
Expected: `ModuleNotFoundError: No module named 'spectral_validation'`.

- [ ] **Step 3: Implement**

Create `Synthetic_domain/spectral_validation.py`:
```python
"""
Rulers for the Bagpipes spectral fits that do not rely on full-spectrum SED fitting:
model-independent spectral indices (compared with the fits' posterior predictions), light-weighted ages and
metallicities (for Gallazzi et al. 2005), the Speagle et al. (2014) main sequence and percentile bands of
literature scaling relations.
"""

import numpy as np
import pandas as pd
from astropy.io import fits

import synthetic_model as sm


# =============================================================================
# 1. Index definitions (air wavelengths, rest frame, Angstrom)
# =============================================================================
D4000_BLUE_BAND = (3850.0, 3950.0)      # Balogh et al. (1999), narrow 4000 A break
D4000_RED_BAND = (4000.0, 4100.0)
EQUIVALENT_WIDTH_INDICES = {
    "Hdelta_A": {"blue": (4041.600, 4079.750), "feature": (4083.500, 4122.250), "red": (4128.500, 4161.000)},
    "Mgb": {"blue": (5142.625, 5161.375), "feature": (5160.125, 5192.625), "red": (5191.375, 5206.375)},
    "Fe5270": {"blue": (5233.150, 5248.150), "feature": (5245.650, 5285.650), "red": (5285.650, 5318.150)},
    "Fe5335": {"blue": (5304.625, 5315.875), "feature": (5312.125, 5352.125), "red": (5353.375, 5363.375)},
}
MGFE_FE5270_WEIGHT = 0.72               # Thomas et al. (2003) [MgFe]' weights
MGFE_FE5335_WEIGHT = 0.28
INDEX_NAMES = ["D4000_n", "Hdelta_A", "Mgb", "Fe5270", "Fe5335", "MgFe_prime"]
MGFE_COMPONENTS = ["Mgb", "Fe5270", "Fe5335"]
N_NOISE_DRAWS = 200                     # Monte-Carlo draws for observed-index errors


def in_band(wavelength, band):
    """True for pixels whose centre lies inside band = (low, high) [Angstrom]."""
    return (wavelength >= band[0]) & (wavelength <= band[1])


def band_mean(wavelength, values, band):
    """Mean of values over the pixels inside band; NaN when the band holds no pixel."""
    inside = in_band(wavelength, band)
    if not inside.any():
        return np.nan
    return np.mean(values[inside])


def d4000_n(wavelength, flux):
    """Narrow 4000 A break: mean F_nu in 4000-4100 A over mean F_nu in 3850-3950 A (F_nu ∝ F_lambda lambda^2).

    Inputs: rest-frame air wavelength [Angstrom]; flux F_lambda (any units). Returns: the ratio (NaN if a band is empty).
    """
    flux_per_frequency = flux * wavelength**2
    red = band_mean(wavelength, flux_per_frequency, D4000_RED_BAND)
    blue = band_mean(wavelength, flux_per_frequency, D4000_BLUE_BAND)
    return red / blue


def equivalent_width(wavelength, flux, definition):
    """Equivalent width [Angstrom] against a straight pseudo-continuum through the blue and red band means.

    Inputs: rest-frame air wavelength [Angstrom]; flux F_lambda; definition: dict with blue, feature, red bands.
    Returns: sum over feature pixels of (1 - F / continuum) * pixel width; NaN if a band holds no pixel.
    """
    blue_flux = band_mean(wavelength, flux, definition["blue"])
    red_flux = band_mean(wavelength, flux, definition["red"])
    feature = in_band(wavelength, definition["feature"])
    if np.isnan(blue_flux) or np.isnan(red_flux) or not feature.any():
        return np.nan
    blue_centre = np.mean(definition["blue"])
    red_centre = np.mean(definition["red"])
    slope = (red_flux - blue_flux) / (red_centre - blue_centre)
    continuum = blue_flux + slope * (wavelength[feature] - blue_centre)
    pixel_widths = np.gradient(wavelength)[feature]
    return np.sum((1 - flux[feature] / continuum) * pixel_widths)


def mgfe_prime(mgb, fe5270, fe5335):
    """[MgFe]' = sqrt(Mgb (0.72 Fe5270 + 0.28 Fe5335)) [Angstrom]; NaN when the argument is negative or NaN."""
    argument = mgb * (MGFE_FE5270_WEIGHT * fe5270 + MGFE_FE5335_WEIGHT * fe5335)
    if not np.isfinite(argument) or argument < 0:
        return np.nan
    return np.sqrt(argument)


def index_values(wavelength, flux):
    """All indices of one spectrum. Inputs: rest-frame air wavelength [Angstrom], flux F_lambda.

    Returns: dict name -> value for INDEX_NAMES (D4000_n dimensionless, the others in Angstrom).
    """
    values = {"D4000_n": d4000_n(wavelength, flux)}
    for name, definition in EQUIVALENT_WIDTH_INDICES.items():
        values[name] = equivalent_width(wavelength, flux, definition)
    values["MgFe_prime"] = mgfe_prime(values["Mgb"], values["Fe5270"], values["Fe5335"])
    return values


def index_bands(name):
    """The wavelength bands an index uses."""
    if name == "D4000_n":
        return [D4000_BLUE_BAND, D4000_RED_BAND]
    return list(EQUIVALENT_WIDTH_INDICES[name].values())


def index_validity(wavelength, ivar):
    """Whether each index can be measured: every band has pixels and none of them is masked (ivar = 0).

    Returns: dict name -> bool for INDEX_NAMES ([MgFe]' needs its three components).
    """
    validity = {}
    for name in INDEX_NAMES:
        if name == "MgFe_prime":
            continue
        valid = True
        for band in index_bands(name):
            inside = in_band(wavelength, band)
            if not inside.any() or np.any(ivar[inside] <= 0):
                valid = False
        validity[name] = valid
    components_valid = []
    for component in MGFE_COMPONENTS:
        components_valid.append(validity[component])
    validity["MgFe_prime"] = all(components_valid)
    return validity


def measure_indices(wavelength, flux, ivar, rng, n_draws=N_NOISE_DRAWS):
    """Observed indices with Monte-Carlo errors from the inverse variance.

    Inputs: rest-frame air wavelength [Angstrom]; flux and ivar (same units as the spectrum); rng; n_draws.
    Returns: dict with <name> and <name>_err for every name in INDEX_NAMES; NaN where index_validity is False.
    """
    validity = index_validity(wavelength, ivar)
    values = index_values(wavelength, flux)
    sigma = np.zeros_like(flux)
    good = ivar > 0
    sigma[good] = 1 / np.sqrt(ivar[good])

    draws = {name: [] for name in INDEX_NAMES}
    for _ in range(n_draws):
        noisy_flux = flux + rng.normal(0, 1, len(flux)) * sigma
        draw_values = index_values(wavelength, noisy_flux)
        for name in INDEX_NAMES:
            draws[name].append(draw_values[name])

    measured = {}
    for name in INDEX_NAMES:
        if validity[name]:
            measured[name] = values[name]
            measured[f"{name}_err"] = np.nanstd(draws[name])
        else:
            measured[name] = np.nan
            measured[f"{name}_err"] = np.nan
    return measured


# =============================================================================
# 2. Observed spectra for the indices
# =============================================================================
def rest_frame_spectrum(sdss_fits_path, redshift, ebv):
    """Native-resolution SDSS spectrum ready for index measurement.

    Inputs: path of an SDSS spec FITS file (HDU 1 with loglam, flux, ivar); redshift; Milky-Way E(B-V) [mag].
    Returns: (rest-frame air wavelength [Angstrom], dereddened flux, ivar), the flux in SDSS units.
    """
    with fits.open(sdss_fits_path) as spectrum_file:
        table = spectrum_file[1].data
        wavelength_vacuum = 10 ** table["loglam"].astype(float)
        flux = table["flux"].astype(float)
        ivar = table["ivar"].astype(float)
    wavelength_air = sm.vacuum_to_air(wavelength_vacuum)
    flux, ivar = sm.deredden_spectrum(wavelength_air, flux, ivar, ebv)
    return wavelength_air / (1 + redshift), flux, ivar
```

- [ ] **Step 4: Run the tests**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_spectral_validation.py -v`
Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/spectral_validation.py Synthetic_domain/tests/test_spectral_validation.py
git commit -m "feat: model-independent spectral indices (D4000n, Hdelta_A, Lick Mg/Fe, [MgFe]') with Monte-Carlo errors"
```

---

### Task 2: Posterior predictive, light-weighted quantities, relations and MPA-JHU SQL

**Files:**
- Modify: `Synthetic_domain/spectral_validation.py` (append sections 3–6), `Synthetic_domain/tests/test_spectral_validation.py` (append)

**Interfaces:**
- Consumes: Task 1 functions; `sm.spectral_components`, `sm.calibration_polynomial`, `sm.bins_from_components`, `sm.OBSERVED_RANGE`, `sm.FILTER_FILES`, `sm.BANDS`, `sm.COSMOLOGY`, `sm.GYR_TO_YEARS`, `sm.LOG_MASS_FORMED`; `bagpipes`.
- Produces:
  - `N_PREDICTIVE_SAMPLES = 100`, `N_LIGHT_WEIGHT_SAMPLES = 50`, `MAX_Z_SCORE = 2.0`, `SPEAGLE_SCATTER = 0.2` (no IMF constant: ruling Q6).
  - `predictive_indices(unit_cube_samples, redshift, observed_wavelength_air, data_flux, data_ivar, n_samples=N_PREDICTIVE_SAMPLES) -> pd.DataFrame` (rows = samples, columns = INDEX_NAMES).
  - `index_comparison(measured, predictive) -> dict` with `<name>_model_median`, `<name>_model_std`, `<name>_z`, `<name>_p` for every name.
  - `light_weighted_age_metallicity(unit_cube_samples, redshift, n_samples=N_LIGHT_WEIGHT_SAMPLES) -> pd.DataFrame` with `log_age_light` [log yr] and `log_metallicity_light` [log Z/Z☉].
  - `speagle_log_sfr(log_mstar, age_gyr) -> np.ndarray`; `speagle_main_sequence(log_mstar, redshift) -> np.ndarray`.
  - `relation_band(x, y, bin_edges) -> pd.DataFrame` (`x_low, x_high, n, y_p16, y_median, y_p84`); `inside_band(x, y, band) -> pd.array("boolean")`.
  - `mpa_jhu_query(specobjids) -> str`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_spectral_validation.py`:
```python
import synthetic_model as sm


def model_galaxy_spectrum(cube, redshift):
    """Noise-free model spectrum of one cube point on a 1 A observed air grid, and that grid."""
    import bagpipes as pipes
    observed_grid = np.arange(sm.OBSERVED_RANGE[0], sm.OBSERVED_RANGE[1], 1.0)
    components, _, _ = sm.spectral_components(cube, redshift)
    model = pipes.model_galaxy(components, spec_wavs=observed_grid)
    return observed_grid, model.spectrum[:, 1]


def test_predictive_indices_reproduce_a_noise_free_model():
    cube = np.full(len(sm.SPECTRAL_DIMS), 0.5)
    redshift = 0.05
    observed_grid, flux = model_galaxy_spectrum(cube, redshift)
    ivar = np.full(flux.shape, 1.0 / (0.01 * np.median(flux)) ** 2)
    predictive = sv.predictive_indices(cube[np.newaxis, :], redshift, observed_grid, flux, ivar, n_samples=1)
    expected = sv.index_values(observed_grid / (1 + redshift), flux)
    for name in sv.INDEX_NAMES:
        assert np.isclose(predictive[name].iloc[0], expected[name], rtol=1e-3, atol=1e-3), name


def test_predictive_indices_skip_degenerate_samples(monkeypatch):
    cube = np.full((2, len(sm.SPECTRAL_DIMS)), 0.5)
    redshift = 0.05
    observed_grid, flux = model_galaxy_spectrum(cube[0], redshift)
    ivar = np.ones_like(flux)
    original = sv.model_flux_for_sample
    calls = {"count": 0}

    def zero_second(model, sample, redshift_value, data_flux, data_ivar):
        calls["count"] += 1
        result = original(model, sample, redshift_value, data_flux, data_ivar)
        return result * 0 if calls["count"] == 2 else result

    monkeypatch.setattr(sv, "model_flux_for_sample", zero_second)
    predictive = sv.predictive_indices(cube, redshift, observed_grid, flux, ivar, n_samples=2)
    assert np.isfinite(predictive["D4000_n"].iloc[0])
    assert np.isnan(predictive["D4000_n"].iloc[1])


def test_index_comparison_z_and_p():
    measured = {name: 1.0 for name in sv.INDEX_NAMES}
    for name in sv.INDEX_NAMES:
        measured[f"{name}_err"] = 0.1
    rng = np.random.default_rng(2)
    predictive = pd.DataFrame({name: rng.normal(1.0, 0.1, 2000) for name in sv.INDEX_NAMES})
    comparison = sv.index_comparison(measured, predictive)
    assert abs(comparison["D4000_n_z"]) < 0.1
    assert comparison["D4000_n_p"] > 0.8


def test_light_weighted_age_of_a_single_populated_bin_lies_in_that_bin():
    cube = np.full(len(sm.SPECTRAL_DIMS), 0.5)
    cube[sm.SPECTRAL_DIMS.index("sfh_type")] = 0.25            # Dirichlet branch
    for dimension_index, dimension in enumerate(sm.SPECTRAL_DIMS):
        if dimension.startswith("dir"):
            cube[dimension_index] = 1e-9                     # almost no mass in every bin ...
    cube[sm.SPECTRAL_DIMS.index("dir10")] = 1 - 1e-9          # ... except bin 10
    result = sv.light_weighted_age_and_metallicity(cube[np.newaxis, :], 0.05, n_samples=1)
    edges = sm.bin_edges(0.05)
    log_age = result["log_age_light"].iloc[0]
    assert np.log10(edges[10] * sm.GYR_TO_YEARS) - 0.1 <= log_age <= np.log10(edges[11] * sm.GYR_TO_YEARS) + 0.1


def test_speagle_main_sequence_value():
    assert np.isclose(sv.speagle_log_sfr(10.0, 13.7), (0.84 - 0.026 * 13.7) * 10 - (6.51 - 0.11 * 13.7))


def test_relation_band_and_inside_band():
    rng = np.random.default_rng(3)
    x = rng.uniform(9, 12, 20000)
    y = 0.5 * x + rng.normal(0, 0.1, 20000)
    band = sv.relation_band(x, y, np.arange(9, 12.01, 0.5))
    bin_centres = (band["x_low"] + band["x_high"]) / 2
    assert np.allclose(band["y_median"], 0.5 * bin_centres, atol=0.03)
    inside = sv.inside_band(np.array([10.25, 10.25]), np.array([5.125, 6.0]), band)
    assert bool(inside[0]) and not bool(inside[1])


def test_inside_band_outside_range_is_na():
    band = pd.DataFrame({"x_low": [9.0], "x_high": [10.0], "n": [100], "y_p16": [0.0], "y_median": [0.5], "y_p84": [1.0]})
    inside = sv.inside_band(np.array([12.0]), np.array([0.5]), band)
    assert pd.isna(inside[0])


def test_mpa_jhu_query_contains_ids_and_columns():
    query = sv.mpa_jhu_query(["123", "10310136606464692224"])
    assert "10310136606464692224" in query and "galSpecIndx" in query and "galSpecExtra" in query
    assert "sfr_fib_p50" in query and "lick_hd_a" in query and "velDisp" in query
```
(Also add `import pandas as pd` at the top of the test file.)

- [ ] **Step 2: Run the tests to verify they fail**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_spectral_validation.py -v`
Expected: the new tests fail with `AttributeError` (`predictive_indices`, …).

- [ ] **Step 3: Implement** — append to `spectral_validation.py`:
```python
import bagpipes as pipes


# =============================================================================
# 3. Posterior-predictive indices
# =============================================================================
N_PREDICTIVE_SAMPLES = 100     # posterior samples per galaxy for the model indices
MAX_Z_SCORE = 2.0              # an index agrees when |z| < MAX_Z_SCORE


def evenly_spaced_rows(samples, n_samples):
    """At most n_samples rows of samples, evenly spaced (posterior samples are equal-weight)."""
    step = max(1, int(np.ceil(len(samples) / n_samples)))
    return samples[::step][:n_samples]


def model_flux_for_sample(model, cube, redshift, data_flux, data_ivar):
    """Model spectrum of one posterior sample on the model's grid, times the calibration polynomial fitted to the data."""
    components, _, _ = sm.spectral_components(cube, redshift)
    model.update(components)
    flux = model.spectrum[:, 1]
    polynomial = sm.calibration_polynomial(flux, data_flux, data_ivar)
    return flux * polynomial


def predictive_indices(unit_cube_samples, redshift, observed_wavelength_air, data_flux, data_ivar,
                       n_samples=N_PREDICTIVE_SAMPLES):
    """Indices of the Bagpipes model spectra of posterior samples, measured exactly like the observed ones.

    Inputs: unit-cube posterior samples (n, SPECTRAL_DIMS); redshift; the galaxy's native observed air
            wavelength grid [Angstrom] and its flux/ivar on that grid (restricted to sm.OBSERVED_RANGE by the caller).
    Returns: DataFrame (samples x INDEX_NAMES); a sample whose model flux is not positive in a band gets NaN there.
    """
    kept = evenly_spaced_rows(unit_cube_samples, n_samples)
    first_components, _, _ = sm.spectral_components(kept[0], redshift)
    model = pipes.model_galaxy(first_components, spec_wavs=observed_wavelength_air)
    rest_wavelength = observed_wavelength_air / (1 + redshift)
    rows = []
    for cube in kept:
        flux = model_flux_for_sample(model, cube, redshift, data_flux, data_ivar)
        positive = np.where(flux > 0, 1.0, 0.0)
        values = index_values(rest_wavelength, flux)
        validity = index_validity(rest_wavelength, positive)    # zero or negative model flux invalidates a band
        for name in INDEX_NAMES:
            if not validity[name] or not np.isfinite(values[name]):
                values[name] = np.nan
        rows.append(values)
    return pd.DataFrame(rows, columns=INDEX_NAMES)


def index_comparison(measured, predictive):
    """Observed vs posterior-predictive indices of one galaxy.

    Inputs: measured (measure_indices output); predictive (predictive_indices output).
    Returns: dict with <name>_model_median, <name>_model_std, <name>_z = (obs - median)/sqrt(err^2 + std^2) and
             <name>_p (two-sided posterior-predictive p-value) for every name; NaN where either side is missing.
    """
    comparison = {}
    for name in INDEX_NAMES:
        model_values = predictive[name].to_numpy(dtype=float)
        model_values = model_values[np.isfinite(model_values)]
        observed = measured[name]
        if len(model_values) == 0 or not np.isfinite(observed):
            for suffix in ("model_median", "model_std", "z", "p"):
                comparison[f"{name}_{suffix}"] = np.nan
            continue
        median = np.median(model_values)
        spread = np.std(model_values)
        comparison[f"{name}_model_median"] = median
        comparison[f"{name}_model_std"] = spread
        comparison[f"{name}_z"] = (observed - median) / np.sqrt(measured[f"{name}_err"] ** 2 + spread**2)
        below = np.mean(model_values < observed)
        comparison[f"{name}_p"] = 2 * min(below, 1 - below)
    return comparison


# =============================================================================
# 4. Light-weighted age and metallicity (for Gallazzi et al. 2005)
# =============================================================================
N_LIGHT_WEIGHT_SAMPLES = 50
LIGHT_WEIGHT_REDSHIFT = 0.001   # near rest frame: the r band then samples rest-frame r-band light


def light_weighted_age_and_metallicity(unit_cube_samples, redshift, n_samples=N_LIGHT_WEIGHT_SAMPLES):
    """r-band light-weighted log age and log metallicity of posterior samples.

    Each existing age bin is modelled alone (dust-free, no nebular emission, near rest frame) to get its r-band
    luminosity L_i; the light-weighted values are sum(L_i * log10 x_i) / sum(L_i) with x_i the bin's mid age [yr]
    or metallicity [Z/Z_sun]. Check this convention against Gallazzi et al. (2005) §2 (they quote r-band
    light-weighted log age and log Z) and keep both definitions identical.
    Inputs: unit-cube samples; the galaxy redshift (sets the bins). Returns: DataFrame with log_age_light [log yr]
            and log_metallicity_light [log Z/Z_sun], one row per kept sample.
    """
    r_filter = [sm.FILTER_FILES[sm.BANDS.index("r")]]
    kept = evenly_spaced_rows(unit_cube_samples, n_samples)
    model = None
    rows = []
    for cube in kept:
        components, _, _ = sm.spectral_components(cube, redshift)
        fractions, metallicities, edges = sm.bins_from_components(components)
        luminosities = np.zeros(len(fractions))
        for bin_index in range(len(fractions)):
            single_bin = {
                "redshift": LIGHT_WEIGHT_REDSHIFT,
                "constant": {"age_min": edges[bin_index], "age_max": edges[bin_index + 1],
                             "massformed": sm.LOG_MASS_FORMED + np.log10(fractions[bin_index]),
                             "metallicity": metallicities[bin_index]},
            }
            if model is None:
                model = pipes.model_galaxy(single_bin, filt_list=r_filter)
            else:
                model.update(single_bin)
            luminosities[bin_index] = model.photometry[0]
        mid_ages = (edges[:-1] + edges[1:]) / 2 * sm.GYR_TO_YEARS
        weights = luminosities / luminosities.sum()
        rows.append({"log_age_light": np.sum(weights * np.log10(mid_ages)),
                     "log_metallicity_light": np.sum(weights * np.log10(metallicities))})
    return pd.DataFrame(rows)


# =============================================================================
# 5. Scaling relations
# =============================================================================
SPEAGLE_SCATTER = 0.2           # [dex] intrinsic scatter of the main sequence (Speagle et al. 2014)
PERCENTILES = (16, 50, 84)


def speagle_log_sfr(log_mstar, age_gyr):
    """Speagle et al. (2014) main sequence: log SFR [M_sun/yr] = (0.84 - 0.026 t) log M* - (6.51 - 0.11 t)."""
    return (0.84 - 0.026 * age_gyr) * np.asarray(log_mstar) - (6.51 - 0.11 * age_gyr)


def speagle_main_sequence(log_mstar, redshift):
    """Main-sequence log SFR at each galaxy's redshift (t = age of the Universe [Gyr] in the module cosmology)."""
    age_gyr = sm.COSMOLOGY.age(np.asarray(redshift)).value
    return speagle_log_sfr(log_mstar, age_gyr)


def relation_band(x, y, bin_edges):
    """Median and 16-84th percentiles of y in bins of x (a literature relation from a large sample).

    Returns: DataFrame with x_low, x_high, n, y_p16, y_median, y_p84 (bins with < 10 points are dropped).
    """
    rows = []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        in_bin = (x >= low) & (x < high) & np.isfinite(y)
        if in_bin.sum() < 10:
            continue
        p16, p50, p84 = np.percentile(y[in_bin], PERCENTILES)
        rows.append({"x_low": low, "x_high": high, "n": int(in_bin.sum()),
                     "y_p16": p16, "y_median": p50, "y_p84": p84})
    return pd.DataFrame(rows)


def inside_band(x, y, band):
    """Whether each (x, y) lies inside the 16-84th band of its x bin; NA outside the binned x range or for NaN."""
    result = []
    for x_value, y_value in zip(x, y):
        matching = band[(band["x_low"] <= x_value) & (band["x_high"] > x_value)]
        if len(matching) == 0 or not np.isfinite(y_value):
            result.append(pd.NA)
        else:
            row = matching.iloc[0]
            result.append(bool(row["y_p16"] <= y_value <= row["y_p84"]))
    return pd.array(result, dtype="boolean")


# =============================================================================
# 6. MPA-JHU (SDSS galSpec) query
# =============================================================================
MPA_JHU_COLUMNS = {
    "galSpecIndx": ["d4000_n", "d4000_n_err", "lick_hd_a", "lick_hd_a_err", "lick_mgb", "lick_mgb_err",
                    "lick_fe5270", "lick_fe5270_err", "lick_fe5335", "lick_fe5335_err"],
    "galSpecExtra": ["bptclass", "sfr_fib_p16", "sfr_fib_p50", "sfr_fib_p84", "lgm_tot_p50", "lgm_fib_p50"],
}


def mpa_jhu_query(specobjids):
    """SkyServer SQL selecting MPA-JHU indices and SFRs plus SDSS velocity dispersions for the given spectra.

    Inputs: SDSS specObjIDs as text (eBOSS ids exceed int64). Returns: the SQL string (one LEFT JOIN per table;
            spectra without galSpec rows come back with NULLs).
    """
    selected = ["CAST(s.specObjID AS VARCHAR(20)) AS specobjid_text", "s.velDisp", "s.velDispErr"]
    for column in MPA_JHU_COLUMNS["galSpecIndx"]:
        selected.append(f"i.{column}")
    for column in MPA_JHU_COLUMNS["galSpecExtra"]:
        selected.append(f"e.{column}")
    id_list = ", ".join(specobjids)
    return ("SELECT " + ", ".join(selected) + " FROM SpecObjAll AS s "
            "LEFT JOIN galSpecIndx AS i ON i.specObjID = s.specObjID "
            "LEFT JOIN galSpecExtra AS e ON e.specObjID = s.specObjID "
            f"WHERE s.specObjID IN ({id_list})")
```
Move the `import bagpipes as pipes` line to the module's import block at the top.

- [ ] **Step 4: Run the tests**

Run: `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest tests/test_spectral_validation.py -v`
Expected: 16 passed. If `test_light_weighted_age_of_a_single_populated_bin...` fails because the other bins are not negligible (MIN_FRACTION floor 1e-6 times young-star luminosity can matter), report the measured value and the per-bin weights before changing anything; the ruling would be to compare against the luminosity-weighted expectation computed from the same weights.

- [ ] **Step 5: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/spectral_validation.py Synthetic_domain/tests/test_spectral_validation.py
git commit -m "feat: posterior-predictive indices, light-weighted age/metallicity, Speagle MS, relation bands, MPA-JHU query"
```

---

### Task 3: Literature data

**Files:**
- Create: `Synthetic_domain/literature/SOURCES.md`, `Synthetic_domain/literature/mpa_jhu_reference.csv`, `Synthetic_domain/literature/relations/mass_metallicity.csv`, `Synthetic_domain/literature/relations/age_sigma_passive.csv`, `Synthetic_domain/literature/gallazzi_reference.csv`; `Synthetic_domain/fetch_literature.py`
- Modify: `.gitignore` (add `Synthetic_domain/literature/raw/`)

**Interfaces:**
- Consumes: `sv.mpa_jhu_query`, `sv.relation_band`; `spectroscopic_reference.csv` (columns `specobjid`, `plate`, `mjd`, `fiberid`, `RUN2D`).
- Produces (read by Task 4):
  - `literature/mpa_jhu_reference.csv`: `specobjid` (text), `velDisp`, `velDispErr`, the `MPA_JHU_COLUMNS` columns.
  - `literature/gallazzi_reference.csv`: `specobjid`, `gallazzi_log_age_p16/p50/p84` [log yr], `gallazzi_log_z_p16/p50/p84` [log Z/Z☉], `gallazzi_log_mstar` [log M☉] (empty table with these columns if unavailable).
  - `literature/relations/mass_metallicity.csv` and `age_sigma_passive.csv`: `relation_band` outputs (x = log M* or log σ [km/s], y = log Z/Z☉ or log age [yr]).
  - `literature/SOURCES.md`: URL, retrieval date, file names, column meanings, conventions (light-weighted definitions, IMF).

`fetch_literature.py` is a small, legible script (run once from `Synthetic_domain/`; re-runnable; it never overwrites `raw/` downloads that exist).

- [ ] **Step 1: Locate the Gallazzi et al. (2005) catalogue.** Search the MPA-JHU DR4 release pages (`https://www.mpa-garching.mpg.de/SDSS/DR4/`) and VizieR (`J/MNRAS/362/41`) for the public table of r-band light-weighted ages and stellar metallicities (percentiles) with plate/mjd/fiber and stellar masses. Download it into `literature/raw/`. Record in `SOURCES.md`: URL, date, file name, the column names used, and how the paper defines the light-weighted quantities (mean of log, or log of mean) — Task 4 uses this to keep the Bagpipes side identical (adjust `sv.light_weighted_age_and_metallicity`'s docstring/formula only if the definition differs, with a test update). If no public table exists, write that in `SOURCES.md` and produce empty `gallazzi_reference.csv` / `mass_metallicity.csv` / `age_sigma_passive.csv` with the documented columns.

- [ ] **Step 2: Write `fetch_literature.py`** with numbered steps:
  1. read `spectroscopic_reference.csv` (dtype specobjid str), take every `specobjid`;
  2. query MPA-JHU in chunks of 200 ids with `astroquery.sdss.SDSS.query_sql(sv.mpa_jhu_query(chunk), data_release=16)`, concatenate, rename `specobjid_text` → `specobjid` (stripped text), write `literature/mpa_jhu_reference.csv`;
  3. if the Gallazzi table exists in `raw/`: map its columns (the names recorded in step 1), match our galaxies by plate/mjd/fiber, write `gallazzi_reference.csv`; build the relation bands from the whole catalogue: mass–metallicity with `sv.relation_band(log M*, log Z p50, np.arange(8.5, 12.01, 0.25))`; for the passive age–σ relation, draw a random subsample of up to 5,000 Gallazzi galaxies with D4000ₙ > 1.6 from MPA-JHU (query their `velDisp` and `d4000_n` by specObjID or plate/mjd/fiber in chunks of 200), then `sv.relation_band(log10 velDisp, log age p50, np.arange(1.8, 2.61, 0.1))`; write both CSVs;
  4. print counts (galaxies with MPA-JHU rows, with Gallazzi rows; relation bins).
  The D4000ₙ > 1.6 passive cut and the bin edges are named constants at the top of the script.

- [ ] **Step 3: Run it** (`../.venv/bin/python fetch_literature.py`) and check: MPA-JHU rows exist for (almost) all SDSS-legacy galaxies (`RUN2D == "26"`) and are NULL/absent for BOSS ones; relation CSVs have ≥ 5 bins each (or are documented as unavailable).

- [ ] **Step 4: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add .gitignore Synthetic_domain/fetch_literature.py Synthetic_domain/literature/SOURCES.md Synthetic_domain/literature/*.csv Synthetic_domain/literature/relations/
git commit -m "data: MPA-JHU and Gallazzi (2005) literature tables and scaling-relation bands for validating the spectral fits"
```

---

### Task 4: Notebook 2 — three validation layers, Firefly removed

**Files:**
- Modify: `Synthetic_domain/spectroscopic_reference.ipynb`

**Interfaces:**
- Consumes: Tasks 1–3 (module functions, literature CSVs); notebook variables `reference`, `fitted`, `FITS_DIRECTORY`, `SPECTRA_DIRECTORY`, `SSFR_PASSIVE`.
- Produces: new columns in `spectroscopic_reference.csv` (`<index>_obs`, `<index>_obs_err`, `<index>_model_median`, `<index>_model_std`, `<index>_z`, `<index>_p` for every index; `mpa_<column>` for the MPA-JHU columns; `bagpipes_log_sfr_young`, `delta_log_sfr_halpha`; `log_age_light_median/_p16/_p84`, `log_metallicity_light_median/_p16/_p84`; `delta_log_age_gallazzi`, `delta_log_z_gallazzi`; `log_mstar_total`; `in_mass_metallicity_band`, `in_main_sequence_band`, `in_age_sigma_band`), `validation_summary.csv`, figures `validation_indices.png`, `validation_halpha_gallazzi.png`, `validation_relations.png`.

- [ ] **Step 1: Remove Firefly as a validation.** In cell 16 delete the Firefly check block (step 2 "Firefly check", the `firefly_scatter_consistency` call and their CSV columns, the two Firefly headline prints); delete cell 17's Firefly tables and figure (keep its `fibre`/`galaxy_type` columns, moving them into cell 16 before the CSV write); remove `FIGURE_FILE` (Firefly figure) from cell 1 and its mention in cell 0; update cell 15's markdown (no Firefly sentence). Firefly stays in cells 3–4 (selection) untouched. Delete `firefly_comparison.png` from the repo (`git rm`).

- [ ] **Step 2: Add section "## 7. Validation of the spectral fits" (markdown + code cells)** after the reference table is built and before the CSV is written. Markdown (≤ 6 lines): the three layers, report-only, which galaxies each covers. Code cell A (constants + layer 1):
```python
import h5py

import spectral_validation as sv

LITERATURE_DIRECTORY = "literature"
VALIDATION_SUMMARY_FILE = "validation_summary.csv"
INDEX_FIGURE_FILE = "validation_indices.png"
TRACER_FIGURE_FILE = "validation_halpha_gallazzi.png"
RELATION_FIGURE_FILE = "validation_relations.png"
VALIDATION_SEED = 11


def load_unit_cube_samples(path):
    """Equal-weight unit-cube posterior samples stored by sm.save_fit (dataset 'unit_cube_samples')."""
    with h5py.File(path, "r") as source:
        return source["unit_cube_samples"][()]


# Layer 1: observed indices vs posterior-predictive model indices, galaxy by galaxy
rng = np.random.default_rng(VALIDATION_SEED)
index_rows = []
for galaxy in reference[fitted].itertuples():
    rest_wavelength, flux, ivar = sv.rest_frame_spectrum(os.path.join(SPECTRA_DIRECTORY, f"{galaxy.specobjid}.fits"),
                                                         galaxy.z, galaxy.ebv)
    measured = sv.measure_indices(rest_wavelength, flux, ivar, rng)
    observed_air = rest_wavelength * (1 + galaxy.z)
    in_fit_range = (observed_air >= sm.OBSERVED_RANGE[0]) & (observed_air <= sm.OBSERVED_RANGE[1])
    posterior_cubes = load_unit_cube_samples(os.path.join(FITS_DIRECTORY, f"{galaxy.specobjid}.h5"))
    predictive = sv.predictive_indices(posterior_cubes, galaxy.z, observed_air[in_fit_range],
                                       flux[in_fit_range], ivar[in_fit_range])
    row = {"specobjid": galaxy.specobjid}
    for name in sv.INDEX_NAMES:
        row[f"{name}_obs"] = measured[name]
        row[f"{name}_obs_err"] = measured[f"{name}_err"]
    row.update(sv.index_comparison(measured, predictive))
    index_rows.append(row)
reference = reference.merge(pd.DataFrame(index_rows), on="specobjid", how="left")
```
Code cell B (layer 2):
```python
# Layer 2a: MPA-JHU indices (cross-check of our index code), Halpha SFR, velocity dispersion
mpa = pd.read_csv(os.path.join(LITERATURE_DIRECTORY, "mpa_jhu_reference.csv"), dtype={"specobjid": str})
mpa = mpa.add_prefix("mpa_").rename(columns={"mpa_specobjid": "specobjid"})
reference = reference.merge(mpa, on="specobjid", how="left")
STAR_FORMING_BPT_CLASSES = (1, 2)   # MPA-JHU BPT classes: 1 star-forming, 2 low-S/N star-forming (ruling Q7)
FIBRE_YOUNG_BIN = "logsfr_bin_00_median"   # Bagpipes SFR over 0-14 Myr, the time scale Halpha traces
reference["bagpipes_log_sfr_young"] = reference[FIBRE_YOUNG_BIN]
reference["delta_log_sfr_halpha"] = (reference["bagpipes_log_sfr_young"]
                                     - (reference["mpa_sfr_fib_p50"]))   # both Kroupa: no IMF offset (ruling Q6)
halpha_usable = (reference["mpa_bptclass"].isin(STAR_FORMING_BPT_CLASSES)
                 & reference["delta_log_sfr_halpha"].notna())

# Layer 2b: Gallazzi et al. (2005) light-weighted ages and metallicities
gallazzi = pd.read_csv(os.path.join(LITERATURE_DIRECTORY, "gallazzi_reference.csv"), dtype={"specobjid": str})
light_rows = []
for galaxy in reference[fitted].itertuples():
    samples = load_unit_cube_samples(os.path.join(FITS_DIRECTORY, f"{galaxy.specobjid}.h5"))
    light = sv.light_weighted_age_and_metallicity(samples, galaxy.z)
    row = {"specobjid": galaxy.specobjid}
    for column in ("log_age_light", "log_metallicity_light"):
        low, middle, high = np.percentile(light[column], [16, 50, 84])
        row[f"{column}_p16"], row[f"{column}_median"], row[f"{column}_p84"] = low, middle, high
    light_rows.append(row)
reference = reference.merge(pd.DataFrame(light_rows), on="specobjid", how="left")
reference = reference.merge(gallazzi, on="specobjid", how="left")
reference["delta_log_age_gallazzi"] = reference["log_age_light_median"] - reference["gallazzi_log_age_p50"]
reference["delta_log_z_gallazzi"] = reference["log_metallicity_light_median"] - reference["gallazzi_log_z_p50"]
```
(If the Gallazzi columns are absent because the table is empty, the merge yields NaN columns; the summary reports "not available".)

Code cell C (layer 3):
```python
# Layer 3: our galaxies on literature scaling relations. Relations use total masses: the 3" aperture mass is
# scaled to the total S-PLUS r magnitude assuming a constant mass-to-light ratio across the galaxy
# splus_index points into the in-memory splus table of cell 3 (filtered on PROB_GAL and re-indexed), not into
# the file's rows, so r_auto is taken from that table; the total magnitude must be the brighter one
reference["r_auto"] = splus["r_auto"].to_numpy()[reference["splus_index"].astype(int)]
total_minus_aperture = (reference["r_auto"] - reference["r_aper"]).dropna()
assert (total_minus_aperture < 0).mean() >= 0.95, "r_auto fainter than r_aper: wrong galaxies (check splus_index)"
aperture_to_total = -sm.DEX_PER_MAGNITUDE * (reference["r_auto"] - reference["r_aper"])
reference["log_mstar_total"] = reference["log_mstar_median"] + aperture_to_total
reference["log_sfr_total"] = reference["log_sfr_median"] + aperture_to_total
star_forming = reference["galaxy_type"].eq("star-forming")
passive = reference["galaxy_type"].eq("passive")

relations_directory = os.path.join(LITERATURE_DIRECTORY, "relations")
mass_metallicity_band = pd.read_csv(os.path.join(relations_directory, "mass_metallicity.csv"))
age_sigma_band = pd.read_csv(os.path.join(relations_directory, "age_sigma_passive.csv"))

reference["in_mass_metallicity_band"] = sv.inside_band(reference["log_mstar_total"].to_numpy(),
                                                       reference["log_metallicity_light_median"].to_numpy(),
                                                       mass_metallicity_band) if len(mass_metallicity_band) else pd.NA
main_sequence = sv.speagle_main_sequence(reference["log_mstar_total"].to_numpy(), reference["z"].to_numpy())
offset_from_sequence = reference["log_sfr_total"].to_numpy() - main_sequence
on_sequence = np.abs(offset_from_sequence) <= sv.SPEAGLE_SCATTER
reference["in_main_sequence_band"] = pd.array(np.where(star_forming, on_sequence, pd.NA), dtype="boolean")
log_sigma = np.log10(reference["mpa_velDisp"].to_numpy(dtype=float))
age_sigma_flags = sv.inside_band(log_sigma, reference["log_age_light_median"].to_numpy(), age_sigma_band) \
    if len(age_sigma_band) else pd.array([pd.NA] * len(reference), dtype="boolean")
reference["in_age_sigma_band"] = pd.array(np.where(passive, age_sigma_flags, pd.NA), dtype="boolean")
```

Code cell D (summary + figures):
```python
GROUPS = {
    "all": reference["fit_status"].eq("fitted"),
    "star-forming": star_forming,
    "passive": passive,
    "SDSS 3\"": reference["fiber_diameter"].eq(3.0),
    "BOSS 2\"": reference["fiber_diameter"].eq(2.0),
}
summary_rows = []
for group_name, in_group in GROUPS.items():
    group = reference[in_group]
    for name in sv.INDEX_NAMES:
        z_values = group[f"{name}_z"].dropna()
        summary_rows.append({"layer": "1 indices", "quantity": name, "group": group_name, "n": len(z_values),
                             "available": len(z_values) > 0,
                             "fraction_abs_z_below_2": np.mean(np.abs(z_values) < sv.MAX_Z_SCORE) if len(z_values) else np.nan,
                             "median_z": z_values.median() if len(z_values) else np.nan,
                             "z_scatter": 1.4826 * (z_values - z_values.median()).abs().median() if len(z_values) else np.nan})
    tracer_columns = {"Halpha SFR": group.loc[halpha_usable[in_group], "delta_log_sfr_halpha"],
                      "Gallazzi age": group["delta_log_age_gallazzi"], "Gallazzi metallicity": group["delta_log_z_gallazzi"]}
    for quantity, offsets in tracer_columns.items():
        offsets = offsets.dropna()
        summary_rows.append({"layer": "2 tracers", "quantity": quantity, "group": group_name, "n": len(offsets),
                             "available": len(offsets) > 0,
                             "median_offset": offsets.median() if len(offsets) else np.nan,
                             "scatter": 1.4826 * (offsets - offsets.median()).abs().median() if len(offsets) else np.nan})
    for quantity, column in (("mass-metallicity", "in_mass_metallicity_band"),
                             ("main sequence", "in_main_sequence_band"), ("age-sigma", "in_age_sigma_band")):
        flags = group[column].dropna()
        summary_rows.append({"layer": "3 relations", "quantity": quantity, "group": group_name, "n": len(flags),
                             "available": len(flags) > 0,
                             "fraction_inside_band": flags.astype(bool).mean() if len(flags) else np.nan})
validation_summary = pd.DataFrame(summary_rows)
validation_summary.to_csv(VALIDATION_SUMMARY_FILE, index=False)
display(validation_summary[validation_summary["group"] == "all"].round(3))

# Figure 1: indices, observed against posterior-predictive model
figure, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].errorbar(reference["D4000_n_obs"], reference["Hdelta_A_obs"], xerr=reference["D4000_n_obs_err"],
                 yerr=reference["Hdelta_A_obs_err"], fmt="o", ms=4, label="observed")
axes[0].scatter(reference["D4000_n_model_median"], reference["Hdelta_A_model_median"], marker="x", label="Bagpipes model")
for galaxy in reference.itertuples():
    axes[0].plot([galaxy.D4000_n_obs, galaxy.D4000_n_model_median], [galaxy.Hdelta_A_obs, galaxy.Hdelta_A_model_median],
                 color="grey", lw=0.5)
axes[0].set(xlabel="D4000$_n$", ylabel=r"H$\delta_A$ [Å]")
axes[0].legend()
axes[1].errorbar(reference["MgFe_prime_model_median"], reference["MgFe_prime_obs"], yerr=reference["MgFe_prime_obs_err"],
                 xerr=reference["MgFe_prime_model_std"], fmt="o", ms=4)
limits = [np.nanmin(reference["MgFe_prime_obs"]), np.nanmax(reference["MgFe_prime_obs"])]
axes[1].plot(limits, limits, color="black", lw=1)
axes[1].set(xlabel="[MgFe]' Bagpipes model [Å]", ylabel="[MgFe]' observed [Å]")
figure.tight_layout()
figure.savefig(INDEX_FIGURE_FILE, dpi=200)
plt.show()

# Figure 2: Halpha SFR and Gallazzi light-weighted age and metallicity
figure, axes = plt.subplots(1, 3, figsize=(16, 5))
comparisons = [
    (reference.loc[halpha_usable, "mpa_sfr_fib_p50"],
     reference.loc[halpha_usable, "bagpipes_log_sfr_young"], "log SFR (Halpha, 3\" fibre)", "log SFR Bagpipes 0-14 Myr"),
    (reference["gallazzi_log_age_p50"], reference["log_age_light_median"], "Gallazzi log age_r [yr]", "Bagpipes log age_r [yr]"),
    (reference["gallazzi_log_z_p50"], reference["log_metallicity_light_median"], "Gallazzi log Z/Z_sun", "Bagpipes log Z/Z_sun"),
]
for axis, (literature_values, bagpipes_values, x_label, y_label) in zip(axes, comparisons):
    axis.scatter(literature_values, bagpipes_values, s=12)
    finite = np.isfinite(literature_values) & np.isfinite(bagpipes_values)
    if finite.any():
        limits = [min(literature_values[finite].min(), bagpipes_values[finite].min()),
                  max(literature_values[finite].max(), bagpipes_values[finite].max())]
        axis.plot(limits, limits, color="black", lw=1)
    else:
        axis.text(0.5, 0.5, "not available", transform=axis.transAxes, ha="center")
    axis.set(xlabel=x_label, ylabel=y_label)
figure.tight_layout()
figure.savefig(TRACER_FIGURE_FILE, dpi=200)
plt.show()

# Figure 3: scaling relations with their 16-84th bands
figure, axes = plt.subplots(1, 3, figsize=(16, 5))
panels = [(mass_metallicity_band, reference["log_mstar_total"], reference["log_metallicity_light_median"],
           "log M* total [M_sun]", "log Z/Z_sun (light-weighted)"),
          (None, reference.loc[star_forming, "log_mstar_total"], reference.loc[star_forming, "log_sfr_total"],
           "log M* total [M_sun]", "log SFR [M_sun/yr]"),
          (age_sigma_band, pd.Series(log_sigma)[passive.to_numpy()], reference.loc[passive, "log_age_light_median"],
           "log sigma [km/s]", "log age_r [yr]")]
for axis, (band, x_values, y_values, x_label, y_label) in zip(axes, panels):
    if band is not None and len(band):
        centres = (band["x_low"] + band["x_high"]) / 2
        axis.fill_between(centres, band["y_p16"], band["y_p84"], color="grey", alpha=0.4, label="literature 16-84%")
        axis.plot(centres, band["y_median"], color="black")
    axis.scatter(x_values, y_values, s=14, label="this work")
    axis.set(xlabel=x_label, ylabel=y_label)
mass_grid = np.linspace(9, 11.5, 50)
median_redshift = reference.loc[star_forming, "z"].median()
axes[1].fill_between(mass_grid, sv.speagle_main_sequence(mass_grid, median_redshift) - sv.SPEAGLE_SCATTER,
                     sv.speagle_main_sequence(mass_grid, median_redshift) + sv.SPEAGLE_SCATTER, color="grey", alpha=0.4,
                     label=f"Speagle+14 at z={median_redshift:.2f}")
for axis in axes:
    axis.legend(fontsize=8)
figure.tight_layout()
figure.savefig(RELATION_FIGURE_FILE, dpi=200)
plt.show()
```

Code cell E (checks):
```python
boss = reference["fiber_diameter"].eq(2.0)
assert reference.loc[boss, "mpa_sfr_fib_p50"].isna().all(), "BOSS spectra have no MPA-JHU row"
mpa_index_offset = (reference["D4000_n_obs"] - reference["mpa_d4000_n"]).median()
print(f"Index code check: median D4000_n difference to MPA-JHU = {mpa_index_offset:.3f}")
```
plus the same median difference for Hδ_A, Mgb, Fe5270, Fe5335 against their MPA-JHU columns, printed.

- [ ] **Step 3: Run notebook 2** (all fits are present; nothing refits; no S-PLUS queries):
`../.venv/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=inference4etgs spectroscopic_reference.ipynb` (background + bounded waits; layer 1 and the light-weighted pass take ~30–60 min on CPU). Check: no error outputs; the index-code check medians are within the typical MPA-JHU errors (report them); summary printed.

- [ ] **Step 4: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/spectroscopic_reference.ipynb Synthetic_domain/spectroscopic_reference.csv Synthetic_domain/validation_summary.csv Synthetic_domain/validation_*.png
git rm Synthetic_domain/firefly_comparison.png
git commit -m "feat: validate the spectral fits with spectral indices, Halpha SFR, Gallazzi (2005) and scaling relations; drop Firefly validation"
```

---

### Task 5: Remove the remaining Firefly validation code

**Files:**
- Modify: `Synthetic_domain/ml_model.ipynb` (delete the Firefly cross-check cell, "7b", and its mention in markdown), `Synthetic_domain/synthetic_model.py` (delete `compare_with_firefly`, `firefly_scatter_consistency`, `MAX_DELTA_LOG_AGE`, `MAX_DELTA_LOG_METALLICITY`; keep `FIREFLY_MISSING`, used by the selection), `Synthetic_domain/tests/test_synthetic_model.py` (delete their tests).

**Interfaces:**
- Consumes: nothing new. **Ordering:** run only after the Phase-3 notebook fix workflow on `ml_model.ipynb` has finished (it edits that notebook).

- [ ] **Step 1:** `grep -rn "compare_with_firefly\|firefly_scatter_consistency\|MAX_DELTA_LOG" Synthetic_domain --include=*.py --include=*.ipynb` — list every use; after Task 4 only `ml_model.ipynb` and the module/tests should remain.
- [ ] **Step 2:** delete the uses and definitions listed above (and the `ml_model.ipynb` cell); keep the notebook's other cells unchanged.
- [ ] **Step 3:** run `LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -q` (all pass) and validate `ml_model.ipynb` with nbformat.
- [ ] **Step 4: Commit**

```bash
cd /home/pedroand/Projects/Inference4ETGs
git add Synthetic_domain/synthetic_model.py Synthetic_domain/tests/test_synthetic_model.py Synthetic_domain/ml_model.ipynb
git commit -m "refactor: remove the Firefly comparison code (Firefly no longer validates anything)"
```
