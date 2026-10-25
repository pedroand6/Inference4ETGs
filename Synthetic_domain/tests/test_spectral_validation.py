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
