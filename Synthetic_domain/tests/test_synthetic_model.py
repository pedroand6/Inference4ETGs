"""Tests for synthetic_model.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest -v"""

import json

import extinction
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


# =============================================================================
# Spectra, likelihood, Firefly and split
# =============================================================================
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
