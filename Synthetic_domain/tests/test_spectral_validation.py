"""Tests for spectral_validation.py. Run from Synthetic_domain/."""

import numpy as np
import pandas as pd
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
