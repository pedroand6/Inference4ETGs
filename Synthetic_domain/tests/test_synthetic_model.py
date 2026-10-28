"""Tests for synthetic_model.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest -v"""

import os
import subprocess
import sys

import extinction
import h5py
import joblib
import numpy as np
import pandas as pd
import pytest
from scipy.stats import qmc

import bagpipes as pipes
import sfh_models
import synthetic_model as sm

MIN_PLAUSIBLE_MAGNITUDE = 5.0    # [mag] AB, for 10^10 M_sun formed
MAX_PLAUSIBLE_MAGNITUDE = 40.0   # [mag] AB
LOW_REDSHIFT = 0.002             # Bagpipes needs z > 0 to build a star formation history


# =============================================================================
# Prior
# =============================================================================
def test_prior_draws_stay_inside_bounds_and_age_below_universe_age():
    rng = np.random.default_rng(0)
    for u in rng.random((500, len(sm.MODELS["delayed"].dims))):
        _, drawn = sm.MODELS["delayed"].components(u)
        assert sm.Z_RANGE[0] <= drawn["redshift"] <= sm.Z_RANGE[1]
        assert sm.AGE_MIN <= drawn["age"] <= sm.maximum_age(drawn["redshift"]) + 1e-12
        assert sm.TAU_RANGE[0] <= drawn["tau"] <= sm.TAU_RANGE[1]
        assert sm.METALLICITY_RANGE[0] <= drawn["metallicity"] <= sm.METALLICITY_RANGE[1]
        assert sm.AV_ISM_RANGE[0] <= drawn["av_ism"] <= sm.AV_ISM_RANGE[1]
        assert 10 ** sm.LOG_TAU_BC_RANGE[0] <= drawn["tau_bc"] <= 10 ** sm.LOG_TAU_BC_RANGE[1]
        assert sm.LOGU_RANGE[0] <= drawn["logU"] <= sm.LOGU_RANGE[1]
        assert drawn["log_mass_formed"] == sm.LOG_MASS_FORMED


def test_age_at_the_top_of_the_cube_never_exceeds_the_universe_age():
    """Edge of the prior: z = 0.7 and the largest age draw, drawn or fixed redshift."""
    cube = np.full(len(sm.MODELS["delayed"].dims), 1.0)
    for fixed_redshift in (None, sm.Z_RANGE[1]):
        _, drawn = sm.MODELS["delayed"].components(cube, fixed_redshift=fixed_redshift)
        assert drawn["age"] <= sm.maximum_age(sm.Z_RANGE[1])
    _, drawn = sm.MODELS["delayed"].components(cube, fixed_redshift=0.3)
    assert drawn["age"] <= sm.maximum_age(0.3)


def test_maximum_age_is_the_universe_age_minus_the_margin_for_scalars_and_arrays():
    redshifts = np.array([sm.Z_RANGE[0], 0.3, sm.Z_RANGE[1]])
    expected = sm.COSMOLOGY.age(redshifts).value - sm.UNIVERSE_AGE_MARGIN
    assert np.allclose(sm.maximum_age(redshifts), expected)
    assert np.isclose(sm.maximum_age(0.3), expected[1])
    assert np.isclose(sm.maximum_age(sm.Z_RANGE[1]), 7.065, atol=1e-3)


def test_uniform_cube_gives_log_uniform_tau_and_metallicity():
    cube = np.full(len(sm.MODELS["delayed"].dims), 0.5)
    _, drawn = sm.MODELS["delayed"].components(cube)
    assert np.isclose(np.log10(drawn["tau"]), np.mean(np.log10(sm.TAU_RANGE)))
    assert np.isclose(np.log10(drawn["metallicity"]), np.mean(np.log10(sm.METALLICITY_RANGE)))
    log_age_range = [np.log10(sm.AGE_MIN), np.log10(sm.maximum_age(drawn["redshift"]))]
    assert np.isclose(np.log10(drawn["age"]), np.mean(log_age_range))


def test_latin_hypercube_is_stratified_in_log_tau():
    """Each of n equal slices of log tau holds exactly one of n Latin hypercube draws."""
    n_draws = 50
    cube = qmc.LatinHypercube(len(sm.MODELS["delayed"].dims), seed=4).random(n_draws)
    log_taus = np.zeros(n_draws)
    for row_index, u in enumerate(cube):
        _, drawn = sm.MODELS["delayed"].components(u)
        log_taus[row_index] = np.log10(drawn["tau"])
    log_tau_range = np.log10(sm.TAU_RANGE)
    slices = np.floor((log_taus - log_tau_range[0]) / np.diff(log_tau_range) * n_draws).astype(int)
    assert np.array_equal(np.sort(slices), np.arange(n_draws))


def test_components_hold_one_delayed_component_and_the_dust_and_nebular_entries():
    u = np.random.default_rng(5).random(len(sm.MODELS["delayed"].dims))
    bagpipes_components, drawn = sm.MODELS["delayed"].components(u)

    assert list(drawn) == sm.MODELS["delayed"].drawn_labels
    delayed = bagpipes_components["delayed"]
    assert delayed == {"age": drawn["age"], "tau": drawn["tau"], "massformed": sm.LOG_MASS_FORMED,
                       "metallicity": drawn["metallicity"]}
    dust = bagpipes_components["dust"]
    assert dust["type"] == "CF00" and dust["n"] == sm.DUST_SLOPE
    assert dust["Av"] == drawn["av_ism"]
    assert np.isclose(dust["eta"], 1 + sm.TAU_TO_AV * drawn["tau_bc"] / drawn["av_ism"])
    assert bagpipes_components["nebular"] == {"logU": drawn["logU"]}
    assert bagpipes_components["t_bc"] == sm.BIRTH_CLOUD_AGE


def test_fixed_redshift_overrides_the_drawn_redshift():
    u = np.full(len(sm.MODELS["delayed"].dims), 0.5)
    bagpipes_components, drawn = sm.MODELS["delayed"].components(u, fixed_redshift=0.123)

    assert bagpipes_components["redshift"] == 0.123
    assert drawn["redshift"] == 0.123


# =============================================================================
# Labels and magnitudes
# =============================================================================
def bagpipes_delayed_sfh(age, tau):
    """Bagpipes SFH of one delayed-tau component (10^10 M_sun formed, Z_sun) at LOW_REDSHIFT.

    Inputs: age, tau [Gyr]. Returns: bagpipes star_formation_history object.
    """
    component = {"age": age, "tau": tau, "massformed": 10.0, "metallicity": 1.0}
    return pipes.models.star_formation_history({"redshift": LOW_REDSHIFT, "delayed": component})


# tau << age, tau ~ age and tau >> age; ages away from LABEL_EDGES, where Bagpipes' grid step is coarse
@pytest.mark.parametrize("age, tau", [(13.0, 0.3), (5.0, 1.0), (2.0, 20.0), (8.0, 3.0), (0.5, 0.1),
                                      (0.05, 0.01)])
def test_analytic_labels_match_bagpipes_sfh(age, tau):
    sfh = bagpipes_delayed_sfh(age, tau)
    mass = sfh.sfh * sfh.age_widths
    ages = sfh.ages / 1e9
    assert np.isclose(sm.delayed_mass_weighted_age(age, tau), sfh.mass_weighted_age, rtol=0.01)
    for limit in (0.1, 1.0):
        younger = mass[ages < limit].sum() / mass.sum()
        assert abs(sm.delayed_mass_fraction_younger_than(age, tau, limit) - younger) < 0.01


@pytest.mark.parametrize("age, tau", [(5.0, 1.0), (2.0, 20.0)])
def test_analytic_log_sfr_matches_bagpipes_sfr(age, tau):
    sfh = bagpipes_delayed_sfh(age, tau)
    labels = sm.delayed_tau_labels(np.array([age]), np.array([tau]), np.array([10.0]),
                                   np.array([sfh.stellar_mass]))
    assert abs(labels["log_sfr"][0] - np.log10(sfh.sfr)) < 0.02
    assert np.isclose(labels["log_ssfr"][0], labels["log_sfr"][0] - sfh.stellar_mass)


def test_analytic_labels_are_finite_in_the_prior_corners():
    t_max = sm.maximum_age(LOW_REDSHIFT)
    for age, tau in [(t_max, 0.01), (0.05, 20.0), (t_max, 20.0), (0.05, 0.01)]:
        labels = sm.delayed_tau_labels(np.array([age]), np.array([tau]), np.array([10.0]), np.array([9.8]))
        assert all(np.isfinite(value).all() for value in labels.values())
        for fraction in ("f_young", "f_mid", "f_old"):
            assert ((labels[fraction] >= 0) & (labels[fraction] <= 1)).all(), fraction
        fractions = labels["f_young"] + labels["f_mid"] + labels["f_old"]
        assert np.allclose(fractions, 1.0)
        assert 0 < labels["mw_age"][0] <= age


def test_analytic_labels_reach_their_limits():
    """tau << age gives mw_age = age - 2 tau; tau >> age gives age / 3; a young galaxy is all f_young."""
    assert np.isclose(sm.delayed_mass_weighted_age(13.0, 0.01), 13.0 - 0.02)
    assert np.isclose(sm.delayed_mass_weighted_age(0.05, 20.0), 0.05 / 3, rtol=1e-3)
    assert sm.delayed_mass_fraction_younger_than(0.05, 1.0, 0.1) == 1.0


def test_analytic_labels_are_vectorised():
    ages = np.array([13.0, 5.0, 0.05])
    taus = np.array([0.3, 1.0, 20.0])
    labels = sm.delayed_tau_labels(ages, taus, np.full(3, 10.0), np.full(3, 9.8))
    for label in ("mw_age", "f_young", "f_mid", "f_old", "log_sfr", "log_ssfr"):
        assert labels[label].shape == (3,), label
        for row_index in range(3):
            single = sm.delayed_tau_labels(ages[row_index:row_index + 1], taus[row_index:row_index + 1],
                                           np.array([10.0]), np.array([9.8]))
            assert np.isclose(labels[label][row_index], single[label][0]), label


def test_no_recent_star_formation_sits_on_the_sfr_floor():
    labels = sm.delayed_tau_labels(np.array([13.0]), np.array([0.01]), np.array([10.0]), np.array([9.7]))
    assert np.isclose(labels["log_sfr"][0], np.log10(sm.SFR_FLOOR) + 9.7)
    assert np.isfinite(labels["log_ssfr"][0])


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


def test_ab_magnitudes_are_nan_when_any_flux_is_not_positive():
    flux = np.full(len(sm.BANDS), 1e-17)
    assert np.all(np.isfinite(sm.ab_magnitudes(flux)))
    flux[3] = 0.0
    assert np.all(np.isnan(sm.ab_magnitudes(flux)))


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
# Spectra, likelihood and split
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


def test_spectral_log_likelihood_guard_returns_log_zero():
    model_flux = 1 + 0.5 * np.sin(np.linspace(0, 20, 1000))
    data_flux = model_flux
    data_ivar = np.ones(1000)
    nan_model = model_flux.copy()
    nan_model[10] = np.nan

    # An all-zero model makes the normal matrix singular, so np.linalg.solve raises LinAlgError
    assert sm.safe_spectral_log_likelihood(np.zeros(1000), data_flux, data_ivar, 1.0) == sm.LOG_ZERO
    assert sm.safe_spectral_log_likelihood(nan_model, data_flux, data_ivar, 1.0) == sm.LOG_ZERO
    expected = sm.spectral_log_likelihood(model_flux, data_flux, data_ivar, 1.0)
    assert sm.safe_spectral_log_likelihood(model_flux, data_flux, data_ivar, 1.0) == expected


def test_spectral_components_fix_redshift_and_add_veldisp():
    cube = np.full(len(sm.SPECTRAL_DIMS), 0.5)
    bagpipes_components, drawn, noise_scale = sm.spectral_components(cube, redshift=0.08)

    assert "z" not in sm.SPECTRAL_DIMS
    assert bagpipes_components["redshift"] == 0.08
    assert sm.VELDISP_RANGE[0] <= bagpipes_components["veldisp"] <= sm.VELDISP_RANGE[1]
    assert drawn["veldisp"] == bagpipes_components["veldisp"]
    assert 1.0 <= noise_scale <= 5.0
    assert bagpipes_components["delayed"]["age"] == drawn["age"] <= sm.maximum_age(0.08)


def test_split_by_hash_is_stable():
    """Known string ids (as the notebook passes them, one above the int64 range) keep their label forever."""
    expected_labels = {
        "299489677444933632": "validation",
        "299490502078654464": "validation",
        "1237645879578460255": "test",
        "42": "validation",
        "10310136606464692224": "test",    # eBOSS id above the int64 maximum (9.22e18)
    }
    for identifier, expected in expected_labels.items():
        assert sm.split_by_hash(identifier) == expected, identifier
    # Same answer in another order and as an integer, with other galaxies around
    for identifier in reversed(list(expected_labels)):
        assert sm.split_by_hash(int(identifier)) == expected_labels[identifier]


def test_split_by_hash_is_balanced():
    """About half of 1,000+ ids go to validation, as the notebook's split assumes."""
    n_ids = 2000
    labels = [sm.split_by_hash(str(299489677444933632 + 1000003 * index)) for index in range(n_ids)]
    validation_fraction = labels.count("validation") / n_ids
    assert set(labels) == {"validation", "test"}
    assert 0.4 < validation_fraction < 0.6


def write_ev_file(basename, modes):
    """Write an ev.dat in MultiNest's format: parameters, log-likelihood, log-weight, then the mode number.

    Inputs: basename: directory prefix (str); modes: mode number of each row (list of int), in file order.
    Returns: None.
    """
    n_parameters = 3
    rows = []
    for row_index, mode in enumerate(modes):
        rows.append([0.5] * n_parameters + [-100.0 + row_index, -5.0, mode])
    np.savetxt(basename + "ev.dat", np.array(rows), fmt="%28.18E")


def test_multinest_iterations_of_a_single_mode_run_is_the_row_count(tmp_path):
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1] * 7)
    assert sm.multinest_iterations(basename) == 7


def test_multinest_iterations_counts_one_iteration_for_all_modes_of_a_step(tmp_path):
    """Two modes alive for 3 iterations, then mode 1 stops and mode 2 runs 2 more: 5 iterations, 8 rows."""
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1, 2, 1, 2, 1, 2, 2, 2])
    assert sm.multinest_iterations(basename) == 5


def test_multinest_iterations_handles_a_mode_that_appears_later(tmp_path):
    """A mode that splits off in iteration 3 adds rows, not iterations: 1,1 then 1,2 then 1,2."""
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1, 1, 1, 2, 1, 2])
    assert sm.multinest_iterations(basename) == 4


def test_multinest_iteration_bounds_when_every_mode_is_replaced(tmp_path):
    """Mode 1 splits into 2 and 3 and is emptied: rows 1 | 1 | 2 3 | 2 3. Only 3 -> 2 shows a boundary
    in the mode numbers, so the lower bound is 3 and the upper bound is 4 (1 -> 2 may be a
    boundary since mode 1 never returns; 2 -> 3 is not, mode 2 returns), which is the true count."""
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1, 1, 2, 3, 2, 3])
    assert sm.multinest_iteration_bounds(basename) == (3, 4)
    assert sm.multinest_iterations(basename) == 3


def test_multinest_iteration_bounds_are_equal_when_modes_keep_their_points(tmp_path):
    """A new mode next to a surviving parent is not a possible boundary: 1 1 | 1 2 | 1 2."""
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1, 1, 1, 2, 1, 2])
    assert sm.multinest_iteration_bounds(basename) == (4, 4)


def test_multinest_iteration_bounds_of_a_one_row_file(tmp_path):
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1])
    assert sm.multinest_iteration_bounds(basename) == (1, 1)


def test_multinest_iterations_of_a_one_row_file(tmp_path):
    basename = str(tmp_path) + "/"
    write_ev_file(basename, [1])
    assert sm.multinest_iterations(basename) == 1


# =============================================================================
# Mass rescaling
# =============================================================================
def label_table(n_rows):
    """Table with every MASS_DEPENDENT_LABELS column at 10.0 and three columns that must not move."""
    table = pd.DataFrame({"mw_age": np.full(n_rows, 5.0), "metallicity": np.full(n_rows, 1.0),
                          "log_ssfr": np.full(n_rows, -10.0)})
    for label in sm.MASS_DEPENDENT_LABELS:
        table[label] = 10.0
    return table


def test_rescale_mass_labels_fainter_galaxy_has_lower_mass():
    """An observed galaxy 2.5 mag fainter than its model has 10**-1 times the mass: -1 dex."""
    table = label_table(3)
    rescaled = sm.rescale_mass_labels(table, 2.5)
    for label in sm.MASS_DEPENDENT_LABELS:
        assert np.allclose(rescaled[label], 9.0), label


def test_rescale_mass_labels_leaves_other_columns_and_the_input_alone():
    table = label_table(3)
    rescaled = sm.rescale_mass_labels(table, 2.5)
    assert (rescaled["mw_age"] == 5.0).all()
    assert (rescaled["metallicity"] == 1.0).all()
    assert (rescaled["log_ssfr"] == -10.0).all()
    for label in sm.MASS_DEPENDENT_LABELS:
        assert (table[label] == 10.0).all(), label    # the input is not modified


def test_rescale_mass_labels_takes_one_offset_per_row_and_nan_gives_nan():
    table = label_table(3)
    rescaled = sm.rescale_mass_labels(table, np.array([-2.5, 0.0, np.nan]))
    assert np.allclose(rescaled["log_mstar"].iloc[:2], [11.0, 10.0])
    assert np.isnan(rescaled["log_mstar"].iloc[2])
    assert np.isnan(rescaled["log_sfr"].iloc[2])
    assert np.isnan(rescaled["log_mass_formed"].iloc[2])
    assert rescaled["mw_age"].iloc[2] == 5.0


# =============================================================================
# MultiNest fit and fit files (run with LD_LIBRARY_PATH=../.venv/lib)
# =============================================================================
def mock_job(job_id, rng):
    """A small mock spectrum job at z = 0.1.

    Inputs: job_id: identifier (str); rng: numpy Generator.
    Returns: (job dict with job_id, spectrum and redshift, truth labels as a Series).
    """
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    template = sm.prepare_spectrum(wavelength_vacuum, np.ones_like(wavelength_vacuum),
                                   np.full_like(wavelength_vacuum, 100.0), ebv=0.0)
    cube = rng.random(len(sm.SPECTRAL_DIMS))
    spectrum, truth = sm.make_mock_spectrum(cube, 0.1, template, rng)
    return {"job_id": job_id, "spectrum": spectrum, "redshift": 0.1}, truth


def test_make_mock_spectrum_returns_truth_labels():
    """The mock spectrum is finite and its truth carries every training label, finite."""
    job, truth = mock_job("mock_test", np.random.default_rng(9))
    assert np.isfinite(job["spectrum"]["flux"]).all()
    for column in (*sm.MODELS[sm.SPECTRAL_MODEL].training_labels, "veldisp", "noise_scale", "model_r"):
        assert np.isfinite(truth[column]), column
    assert truth["noise_scale"] == 1.0    # the mock noise is the template's own


def test_fit_and_save_smoke_run(tmp_path):
    """A short MultiNest run on a mock spectrum is saved with posterior labels and metadata."""
    job, _ = mock_job("smoke", np.random.default_rng(10))
    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=100, n_live_points=40)

    assert summary["status"] == "fitted", summary["error"]
    posterior, metadata = sm.load_fit(str(tmp_path / "smoke.h5"))
    assert len(posterior) > 0
    for column in (*sm.MODELS[sm.SPECTRAL_MODEL].training_labels, "model_r", "veldisp", "noise_scale"):
        assert column in posterior.columns, column
    assert metadata["converged"] in (True, False)
    assert metadata["max_iterations"] == 100
    assert metadata["n_iterations"] <= 100
    assert metadata["sfh_model"] == "delayed"


def test_fit_is_current_needs_the_delayed_model_and_the_production_iteration_limit():
    """Fits without sfh_model (16-bin era), of another model or with a smaller limit are not current."""
    assert not sm.fit_is_current({"converged": True, "max_iterations": sm.MAX_ITERATIONS}, "delayed")
    assert sm.fit_is_current({"sfh_model": "delayed", "max_iterations": sm.MAX_ITERATIONS}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "nonparametric", "max_iterations": sm.MAX_ITERATIONS}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "delayed", "max_iterations": sm.MAX_ITERATIONS - 1}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "delayed"}, "delayed")
    assert sm.MODELS["delayed"].name == "delayed"


def test_fit_is_current_compares_the_fit_model_with_the_folder_model():
    """A fit is current only in a folder of its own model; FIT_FOLDER_MODELS names each folder's model."""
    full_run = sm.MAX_ITERATIONS
    assert not sm.fit_is_current({"sfh_model": "delayed", "max_iterations": full_run}, "dblplaw")
    assert sm.fit_is_current({"sfh_model": "dblplaw", "max_iterations": full_run}, "dblplaw")
    assert not sm.fit_is_current({"sfh_model": "dblplaw", "max_iterations": full_run}, "delayed")
    assert not sm.fit_is_current({"sfh_model": "dblplaw", "max_iterations": full_run - 1}, "dblplaw")
    assert sm.FIT_FOLDER_MODELS == {"spectral_fits": "delayed", "photometric_fits": "delayed",
                                    "photometric_fits_dblplaw": "dblplaw"}


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


def test_a_job_without_sfh_model_is_delayed():
    """Every job exported before the DPL model has no sfh_model: it is a delayed-tau job."""
    assert sm.job_model({"job_id": "123"}) == "delayed"
    assert sm.job_model({"job_id": "cross_calib_000", "sfh_model": "dblplaw"}) == "dblplaw"


def test_spectral_dims_stay_the_delayed_photometric_dims_plus_the_nuisances():
    """Spectral fits are delayed-tau only (DPL spec, decision 1); the DPL photometric fit has 7 dimensions."""
    assert sm.SPECTRAL_MODEL == "delayed"
    assert sm.SPECTRAL_DIMS == [*sm.MODELS["delayed"].photometric_dims, "veldisp", "noise_scale"]
    assert len(sm.MODELS["dblplaw"].photometric_dims) == 7


def test_posterior_labels_of_midpoint_cubes_give_every_training_label():
    """Three midpoint samples: one row each, every training-label column finite and as the prior maps it."""
    cubes = np.full((3, len(sm.SPECTRAL_DIMS)), 0.5)
    posterior = sm.posterior_labels(cubes, 0.1)

    training_labels = sm.MODELS[sm.SPECTRAL_MODEL].training_labels
    assert len(posterior) == 3
    assert list(posterior.columns[:len(training_labels)]) == training_labels
    for column in (*training_labels, "veldisp", "noise_scale", "model_r"):
        assert np.isfinite(posterior[column]).all(), column
    _, drawn, _ = sm.spectral_components(cubes[0], 0.1)
    for label in ("redshift", "age", "tau", "metallicity", "av_ism", "tau_bc", "logU", "log_mass_formed"):
        assert posterior[label].iloc[0] == drawn[label], label
    fractions = posterior[["f_young", "f_mid", "f_old"]].sum(axis=1)
    assert np.allclose(fractions, 1.0)


def test_fit_and_save_skips_existing_output(tmp_path):
    """A job whose .h5 file already exists is skipped and its run time is read from the file."""
    job, _ = mock_job("skip", np.random.default_rng(11))
    posterior = pd.DataFrame({"mw_age": [1.0, 2.0]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "log_evidence": -1.0}
    sm.save_fit(str(tmp_path / "skip.h5"), posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), metadata)

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=100, n_live_points=40)

    assert summary["status"] == "skipped"
    assert summary["run_time"] == 5.0


def test_save_fit_failure_leaves_no_file_at_the_final_path(tmp_path):
    """A write that fails part-way leaves no file at the final path, so the job is not taken as done."""
    posterior = pd.DataFrame({"mw_age": [1.0, 2.0]})
    path = str(tmp_path / "failed_write.h5")
    with pytest.raises(TypeError):
        sm.save_fit(path, posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), {"converged": object()})
    assert not os.path.exists(path)


def test_fit_and_save_refits_a_truncated_output_file(tmp_path):
    """An unreadable .h5 left by an interrupted run is not counted as done: the job is fitted again."""
    job, _ = mock_job("truncated", np.random.default_rng(12))
    posterior = pd.DataFrame({"mw_age": [1.0, 2.0]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "log_evidence": -1.0}
    path = str(tmp_path / "truncated.h5")
    sm.save_fit(path, posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), metadata)
    with open(path, "r+b") as damaged:
        damaged.truncate(os.path.getsize(path) // 2)

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=100, n_live_points=40)

    assert summary["status"] == "fitted", summary["error"]
    posterior, _ = sm.load_fit(path)
    assert len(posterior) > 0


def test_summarise_posterior_percentiles():
    """The 16th, 50th and 84th percentiles of 0..100 are 16, 50 and 84."""
    posterior = pd.DataFrame({"mw_age": np.arange(101.0)})
    summary = sm.summarise_posterior(posterior, ["mw_age"])
    assert summary["mw_age_median"] == 50.0
    assert summary["mw_age_p16"] == 16.0
    assert summary["mw_age_p84"] == 84.0


# =============================================================================
# Photometric fit (spec amendment A; run with LD_LIBRARY_PATH=../.venv/lib)
# =============================================================================
TOY_MODEL_MAGNITUDES = np.array([20.0, 19.5, 19.2])     # [mag]
TOY_OBSERVED_MAGNITUDES = np.array([21.1, 20.4, 20.35])  # [mag] about 1 mag fainter, with colour noise
TOY_ERRORS = np.array([0.05, 0.08, 0.1])                 # [mag]
OFFSET_INTEGRAL_HALF_WIDTH = 20                          # [offset sigmas] the integrand is ~exp(-200) beyond


def multinest_available():
    """Whether pymultinest imports (it needs libmultinest on LD_LIBRARY_PATH). Returns: bool."""
    try:
        import pymultinest    # noqa: F401
    except (ImportError, OSError, SystemExit):    # pymultinest calls sys.exit without the library
        return False
    return True


def noiseless_photometry(cube, redshift, sfh_model="delayed"):
    """Model AB magnitudes [mag] of one photometric-fit point of sfh_model, and the unscaled stellar mass [log M_sun]."""
    import bagpipes as pipes

    bagpipes_components, _ = sm.photometric_components(cube, redshift, sfh_model)
    model = pipes.model_galaxy(bagpipes_components, filt_list=sm.FILTER_FILES)
    return sm.ab_magnitudes(model.photometry), model.sfh.stellar_mass


def test_offset_marginal_equals_a_numerical_integral_over_the_offset():
    """The analytic marginal over a flat-prior offset equals a brute-force integral (3-band toy)."""
    from scipy.integrate import quad

    detected = np.ones(3, dtype=bool)
    log_likelihood, best_offset, offset_sigma = sm.offset_marginal_log_likelihood(
        TOY_MODEL_MAGNITUDES, TOY_OBSERVED_MAGNITUDES, TOY_ERRORS, detected)

    weights = 1 / TOY_ERRORS**2

    def likelihood_at_offset(offset):
        residuals = TOY_OBSERVED_MAGNITUDES - TOY_MODEL_MAGNITUDES - offset
        return np.exp(-0.5 * np.sum(weights * residuals**2))

    half_width = OFFSET_INTEGRAL_HALF_WIDTH * offset_sigma
    integral, _ = quad(likelihood_at_offset, best_offset - half_width, best_offset + half_width,
                       points=[best_offset], epsabs=0, epsrel=1e-12)
    assert np.isclose(log_likelihood, np.log(integral), rtol=0, atol=1e-6)
    differences = TOY_OBSERVED_MAGNITUDES - TOY_MODEL_MAGNITUDES
    assert np.isclose(best_offset, np.sum(weights * differences) / np.sum(weights))
    assert np.isclose(offset_sigma, 1 / np.sqrt(np.sum(weights)))


def test_offset_marginal_ignores_undetected_bands():
    """A fourth, undetected band (99 mag, NaN model, zero error) leaves the marginal unchanged."""
    three_bands = sm.offset_marginal_log_likelihood(TOY_MODEL_MAGNITUDES, TOY_OBSERVED_MAGNITUDES, TOY_ERRORS,
                                                    np.ones(3, dtype=bool))
    model_magnitudes = np.append(TOY_MODEL_MAGNITUDES, np.nan)
    observed_magnitudes = np.append(TOY_OBSERVED_MAGNITUDES, 99.0)
    errors = np.append(TOY_ERRORS, 0.0)
    detected = np.array([True, True, True, False])
    four_bands = sm.offset_marginal_log_likelihood(model_magnitudes, observed_magnitudes, errors, detected)
    assert np.allclose(four_bands, three_bands)


def test_photometric_log_likelihood_guard_returns_log_zero():
    """A model with a NaN magnitude in a detected band is an impossible point; one in an undetected band is not."""
    detected = np.ones(3, dtype=bool)
    nan_model = TOY_MODEL_MAGNITUDES.copy()
    nan_model[1] = np.nan
    assert sm.safe_photometric_log_likelihood(nan_model, TOY_OBSERVED_MAGNITUDES, TOY_ERRORS, detected) == sm.LOG_ZERO

    expected = sm.offset_marginal_log_likelihood(TOY_MODEL_MAGNITUDES, TOY_OBSERVED_MAGNITUDES, TOY_ERRORS,
                                                 detected)[0]
    safe = sm.safe_photometric_log_likelihood(TOY_MODEL_MAGNITUDES, TOY_OBSERVED_MAGNITUDES, TOY_ERRORS, detected)
    assert safe == expected
    two_detected = np.array([True, False, True])
    undetected_nan = sm.safe_photometric_log_likelihood(nan_model, TOY_OBSERVED_MAGNITUDES, TOY_ERRORS, two_detected)
    assert undetected_nan > sm.LOG_ZERO


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


def test_fit_evidence_refuses_a_fit_that_is_not_photometric(tmp_path):
    """A spectral fit, or one stored without a kind (so spectral), has no photometric likelihood."""
    posterior = pd.DataFrame({"log_age": [0.5]})
    metadata = {"converged": True, "run_time": 1.0, "n_iterations": 5, "max_iterations": sm.MAX_ITERATIONS,
                "log_evidence": -5.0, "log_evidence_error": 0.25, "redshift": 0.1}
    samples = np.zeros((1, len(sm.SPECTRAL_DIMS)))
    sm.save_fit(str(tmp_path / "spectrum.h5"), posterior, samples, dict(metadata, kind="spectrum"))
    sm.save_fit(str(tmp_path / "no_kind.h5"), posterior, samples, metadata)
    job = {"magnitudes": np.full(12, 20.0), "errors": np.full(12, 0.05), "redshift": 0.1}

    for name in ("spectrum.h5", "no_kind.h5"):
        with pytest.raises(ValueError, match="not a photometric fit"):
            sm.fit_evidence(str(tmp_path / name), job)


# =============================================================================
# Job files and the calibration test
# =============================================================================
def test_export_and_load_job_round_trip(tmp_path):
    job, _ = mock_job("round_trip", np.random.default_rng(13))
    path = sm.export_job(job, str(tmp_path / "jobs"))    # the folder does not exist yet

    assert path == str(tmp_path / "jobs" / "round_trip.npz")
    loaded = sm.load_job(path)
    assert loaded["job_id"] == "round_trip"
    assert loaded["redshift"] == job["redshift"]
    for key in ("wavelength", "flux", "ivar"):
        assert np.array_equal(loaded["spectrum"][key], job["spectrum"][key]), key
    assert loaded["kind"] == "spectrum"


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


def test_export_and_load_photometry_job_round_trip(tmp_path):
    job = photometry_job("10311134961714089984", np.random.default_rng(30))
    path = sm.export_job(job, str(tmp_path / "jobs"))

    loaded = sm.load_job(path)
    assert loaded["kind"] == "photometry"
    assert loaded["job_id"] == "10311134961714089984"
    assert loaded["redshift"] == 0.1
    assert np.array_equal(loaded["magnitudes"], job["magnitudes"])
    assert np.array_equal(loaded["errors"], job["errors"])
    assert "spectrum" not in loaded


def test_job_file_written_before_kinds_loads_as_a_spectrum_job(tmp_path):
    """Job files exported before the photometric fits existed have no kind: they are spectral jobs."""
    job, _ = mock_job("old_job", np.random.default_rng(31))
    path = str(tmp_path / "old_job.npz")
    np.savez(path, wavelength=job["spectrum"]["wavelength"], flux=job["spectrum"]["flux"],
             ivar=job["spectrum"]["ivar"], redshift=job["redshift"], job_id="old_job")

    loaded = sm.load_job(path)
    assert loaded["kind"] == "spectrum"
    assert np.array_equal(loaded["spectrum"]["flux"], job["spectrum"]["flux"])


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


def test_fit_and_save_routes_each_job_to_its_fitter_and_stores_the_kind(tmp_path, fake_fitters):
    spectral_job, _ = mock_job("spectral", np.random.default_rng(32))    # no kind: a spectral job by default
    photometric_job = {"job_id": "photometric", "kind": "photometry", "magnitudes": np.full(12, 20.0),
                       "errors": np.full(12, 0.05), "redshift": 0.1}

    spectral_summary = sm.fit_and_save(spectral_job, str(tmp_path), max_iterations=10, n_live_points=25)
    photometric_summary = sm.fit_and_save(photometric_job, str(tmp_path), max_iterations=10, n_live_points=25)

    assert fake_fitters == ["spectrum", "photometry"]
    assert spectral_summary["status"] == "fitted", spectral_summary["error"]
    assert photometric_summary["status"] == "fitted", photometric_summary["error"]
    spectral_posterior, spectral_metadata = sm.load_fit(str(tmp_path / "spectral.h5"))
    photometric_posterior, photometric_metadata = sm.load_fit(str(tmp_path / "photometric.h5"))
    assert spectral_metadata["kind"] == "spectrum"
    assert photometric_metadata["kind"] == "photometry"
    assert list(photometric_posterior["mw_age"]) == [3.0, 4.0]


def test_fit_and_save_never_replaces_a_fit_of_the_other_kind(tmp_path, fake_fitters):
    """Job ids are specobjids, shared by the spectral and photometric fits of one galaxy."""
    posterior = pd.DataFrame({"mw_age": [1.0, 2.0]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "log_evidence": -1.0}
    path = str(tmp_path / "1105771463723476992.h5")
    sm.save_fit(path, posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), metadata)   # no kind: a spectral fit
    photometric_job = {"job_id": "1105771463723476992", "kind": "photometry", "magnitudes": np.full(12, 20.0),
                       "errors": np.full(12, 0.05), "redshift": 0.1}

    summary = sm.fit_and_save(photometric_job, str(tmp_path), max_iterations=10, n_live_points=25)

    assert summary["status"] == "failed"
    assert "spectrum fit" in summary["error"]
    assert fake_fitters == []
    _, stored_metadata = sm.load_fit(path)
    assert "kind" not in stored_metadata    # the spectral file is untouched


@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_fit_and_save_photometry_smoke_run(tmp_path):
    """A short MultiNest run on a photometry job is saved with rescaled labels, the offsets and the kind."""
    job = photometry_job("phot_smoke", np.random.default_rng(33))
    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=200, n_live_points=25)

    assert summary["status"] == "fitted", summary["error"]
    posterior, metadata = sm.load_fit(str(tmp_path / "phot_smoke.h5"))
    assert metadata["kind"] == "photometry"
    for column in (*sm.MODELS["delayed"].training_labels, "magnitude_offset"):
        assert column in posterior.columns, column
        assert np.isfinite(posterior[column]).all(), column
    assert metadata["sfh_model"] == "delayed"


def test_fit_and_save_never_reads_or_replaces_a_fit_of_an_older_sfh_model(tmp_path, fake_fitters):
    """A same-kind file without sfh_model (a 16-bin fit) is neither skipped as done nor overwritten."""
    path = str(tmp_path / "old_model.h5")
    posterior = pd.DataFrame({"mw_age": [1.0, 2.0]})
    metadata = {"converged": True, "run_time": 5.0, "n_iterations": 10, "log_evidence": -1.0}
    sm.save_fit(path, posterior, np.zeros((2, len(sm.SPECTRAL_DIMS))), metadata)
    with h5py.File(path, "a") as stored:
        del stored.attrs["sfh_model"]    # what every fit written before the delayed-tau model looks like
    job, _ = mock_job("old_model", np.random.default_rng(35))

    summary = sm.fit_and_save(job, str(tmp_path), max_iterations=10, n_live_points=25)

    assert summary["status"] == "failed"
    assert "superseded" in summary["error"]
    assert fake_fitters == []
    _, stored_metadata = sm.load_fit(path)
    assert "sfh_model" not in stored_metadata    # the old file is untouched


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


def test_calibration_cubes_are_deterministic_and_fill_the_spectral_cube():
    first = sm.draw_calibration_cubes(60, np.random.default_rng(0))
    second = sm.draw_calibration_cubes(60, np.random.default_rng(0))
    assert np.array_equal(first, second)
    assert first.shape == (60, len(sm.SPECTRAL_DIMS))
    assert np.all((first >= 0) & (first < 1))


def test_calibration_labels_are_the_log_parameters_and_the_stellar_mass():
    calibration_labels = sm.MODELS[sm.SPECTRAL_MODEL].calibration_labels
    assert calibration_labels == ["log_age", "log_tau", "log_metallicity", "log_mstar"]
    assert len(sm.SPECTRAL_DIMS) == 8 and len(sm.MODELS["delayed"].photometric_dims) == 6


def test_calibration_coverage_runs_on_the_calibration_labels():
    """Truth = one posterior sample: every calibration label has one usable mock and a finite coverage."""
    cubes = np.random.default_rng(36).random((3, len(sm.SPECTRAL_DIMS)))
    posterior = sm.posterior_labels(cubes, 0.1)
    truth = posterior.iloc[1]

    calibration_labels = sm.MODELS[sm.SPECTRAL_MODEL].calibration_labels
    table = sm.calibration_coverage({"calib_000": posterior}, {"calib_000": truth}, calibration_labels)

    assert len(table) == len(calibration_labels) * len(sm.CALIBRATION_LEVELS)
    assert (table["n_jobs"] == 1).all()
    assert table["coverage"].notna().all()


def test_draw_calibration_mocks_ids_count_and_seed():
    wavelength_vacuum = 10 ** np.arange(np.log10(3700), np.log10(9300), 1e-4)
    template = sm.prepare_spectrum(wavelength_vacuum, np.ones_like(wavelength_vacuum),
                                   np.full_like(wavelength_vacuum, 100.0), ebv=0.0)
    templates = [(template, 0.1), (template, 0.2)]

    jobs, truths = sm.draw_calibration_mocks(3, templates, np.random.default_rng(4))
    repeat_jobs, _ = sm.draw_calibration_mocks(3, templates, np.random.default_rng(4))

    assert [job["job_id"] for job in jobs] == ["calib_000", "calib_001", "calib_002"]
    assert set(truths) == {job["job_id"] for job in jobs}
    assert {job["redshift"] for job in jobs} <= {0.1, 0.2}
    for job, repeat in zip(jobs, repeat_jobs):
        assert np.array_equal(job["spectrum"]["flux"], repeat["spectrum"]["flux"])
        assert job["redshift"] == repeat["redshift"]


def test_calibration_coverage_of_a_calibrated_and_a_biased_label():
    """Truths drawn from the posterior itself are covered at about the level; shifted truths never are."""
    rng = np.random.default_rng(5)
    posteriors, truths = {}, {}
    for job_index in range(400):
        posteriors[f"job_{job_index}"] = pd.DataFrame({"good": rng.normal(0, 1, 2000),
                                                       "biased": rng.normal(0, 1, 2000)})
        truths[f"job_{job_index}"] = pd.Series({"good": rng.normal(0, 1), "biased": 10.0})

    table = sm.calibration_coverage(posteriors, truths, ["good", "biased"], levels=(0.68, 0.95))

    assert list(table.columns) == ["label", "level", "coverage", "n_jobs"]
    assert len(table) == 4
    good = table[table["label"] == "good"].set_index("level")
    assert abs(good.loc[0.68, "coverage"] - 0.68) < 0.07
    assert abs(good.loc[0.95, "coverage"] - 0.95) < 0.04
    assert (table[table["label"] == "biased"]["coverage"] == 0.0).all()
    assert (table["n_jobs"] == 400).all()


def test_calibration_coverage_known_small_case_and_unusable_jobs():
    samples = pd.DataFrame({"mw_age": np.arange(101.0)})    # 68% interval: 16 to 84; 95%: 2.5 to 97.5
    posteriors = {"a": samples, "b": samples, "c": samples, "d": samples}
    truths = {"a": pd.Series({"mw_age": 50.0}), "b": pd.Series({"mw_age": 99.0}),
              "c": pd.Series({"mw_age": np.nan})}    # c has no finite truth and d has no truth at all

    table = sm.calibration_coverage(posteriors, truths, ["mw_age"]).set_index("level")

    assert table.loc[0.68, "coverage"] == 0.5 and table.loc[0.68, "n_jobs"] == 2
    assert table.loc[0.95, "coverage"] == 0.5 and table.loc[0.95, "n_jobs"] == 2


def test_calibration_coverage_without_jobs_is_nan():
    table = sm.calibration_coverage({}, {}, ["mw_age"])
    assert table["coverage"].isna().all()
    assert (table["n_jobs"] == 0).all()


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


def test_calibration_mock_donors_are_independent_of_the_latin_hypercube(monkeypatch):
    """With one shared seed, mock 1's donor tracked mock 0's first-coordinate jitter (correlation 0.99999)."""
    n_mocks, n_donors = 33, 200
    captured = {}

    def recorded_noiseless_mocks(cubes, redshifts, sfh_model):
        captured["cubes"], captured["redshifts"] = cubes, redshifts
        truth = {**dict.fromkeys(sm.MASS_DEPENDENT_LABELS, 10.0), "model_r": 18.0}
        truths = [pd.Series(truth) for _ in cubes]
        return np.full((len(cubes), len(sm.BANDS)), 18.0), truths

    monkeypatch.setattr(sm, "noiseless_photometric_mocks", recorded_noiseless_mocks)
    noise_table = pd.read_csv(NOISE_MODEL_FILE)
    donor_redshifts = np.arange(n_donors) / 1000    # the redshift gives back the donor index
    donor_r = np.full(n_donors, 17.0)
    donor_fractions, jitters = [], []
    for seed in range(200):
        sm.draw_photometric_calibration_mocks(n_mocks, "dblplaw", donor_redshifts, donor_r, noise_table, seed)
        donor_fractions.append(np.rint(captured["redshifts"][1] * 1000) / n_donors)
        jitters.append(1 - np.modf(n_mocks * captured["cubes"][0, 0])[0])

    assert abs(np.corrcoef(donor_fractions, jitters)[0, 1]) < 0.3


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


# =============================================================================
# Superseding the results of an earlier model (spectroscopic_reference.ipynb, section 4)
# =============================================================================
def write_fit(path, max_iterations=sm.MAX_ITERATIONS, sfh_model="delayed"):
    """A small fit file; sfh_model=None writes it as a 16-bin fit did, without the key."""
    metadata = {"converged": True, "run_time": 1.0, "n_iterations": 5, "max_iterations": max_iterations}
    sm.save_fit(str(path), pd.DataFrame({"log_age": [0.5]}), np.zeros((1, len(sm.SPECTRAL_DIMS))), metadata)
    with h5py.File(path, "a") as fit_file:
        if sfh_model is None:
            del fit_file.attrs["sfh_model"]
        else:
            fit_file.attrs["sfh_model"] = sfh_model


def write_mock_file(path, truth_labels):
    """A calibration-mock file as the notebook stores it: (jobs, truths), each truth a Series."""
    truths = {"calib_000": pd.Series(np.zeros(len(truth_labels)), index=truth_labels)}
    joblib.dump(([{"job_id": "calib_000"}], truths), path)


def test_supersede_moves_old_fits_with_their_folders_and_keeps_current_ones(tmp_path):
    """A fit without sfh_model = "delayed" (or from a pilot run) is superseded; a current one stays."""
    write_fit(tmp_path / "123.h5", sfh_model=None)
    (tmp_path / "123").mkdir()
    (tmp_path / "123" / "ev.dat").write_text("multinest output")
    write_fit(tmp_path / "456.h5", sfh_model="nonparametric")
    write_fit(tmp_path / "pilot_789.h5", max_iterations=100)
    write_fit(tmp_path / "current.h5")
    (tmp_path / "broken.h5").write_bytes(b"not an hdf5 file")

    moved = sm.supersede_outdated_fits(str(tmp_path), "delayed")

    superseded = tmp_path / sm.SUPERSEDED_SUBDIRECTORY
    assert sorted(moved) == ["123", "123.h5", "456.h5", "pilot_789.h5"]
    assert (superseded / "123.h5").exists() and (superseded / "123" / "ev.dat").exists()
    assert (superseded / "456.h5").exists() and (superseded / "pilot_789.h5").exists()
    assert (tmp_path / "current.h5").exists()
    assert (tmp_path / "broken.h5").exists()    # unreadable: fit_and_save fits that job again
    assert not (tmp_path / "123.h5").exists() and not (tmp_path / "123").exists()


def test_supersede_renames_a_name_already_taken_in_superseded(tmp_path):
    """Earlier superseded fits of the same job id get the _16bin suffix: nothing is skipped or overwritten."""
    superseded = tmp_path / sm.SUPERSEDED_SUBDIRECTORY
    superseded.mkdir()
    (superseded / "123.h5").write_text("first-run fit")
    (superseded / "123").mkdir()
    write_fit(tmp_path / "123.h5", sfh_model=None)
    (tmp_path / "123").mkdir()

    moved = sm.supersede_outdated_fits(str(tmp_path), "delayed")

    assert sorted(moved) == ["123", "123.h5"]
    assert (superseded / "123.h5").read_text() == "first-run fit"
    assert (superseded / f"123{sm.SUPERSEDED_NAME_SUFFIX}.h5").exists()
    assert (superseded / f"123{sm.SUPERSEDED_NAME_SUFFIX}").is_dir()
    assert not (tmp_path / "123.h5").exists()
    assert sm.SUPERSEDED_NAME_SUFFIX == "_16bin"


def test_move_to_superseded_counts_up_when_the_suffixed_name_is_taken_too(tmp_path):
    """A third copy of one name gets _16bin_2, a fourth _16bin_3: no earlier superseded copy is overwritten."""
    superseded = tmp_path / sm.SUPERSEDED_SUBDIRECTORY
    destinations = []
    for copy in ("first", "second", "third", "fourth"):
        (tmp_path / "123.h5").write_text(copy)
        destinations.append(sm.move_to_superseded(str(tmp_path), "123.h5"))

    suffix = sm.SUPERSEDED_NAME_SUFFIX
    names = ["123.h5", f"123{suffix}.h5", f"123{suffix}_2.h5", f"123{suffix}_3.h5"]
    assert destinations == [str(superseded / name) for name in names]
    assert [(superseded / name).read_text() for name in names] == ["first", "second", "third", "fourth"]
    assert sm.move_to_superseded(str(tmp_path), "123.h5") is None    # nothing left to move


def test_calibration_mocks_without_truths_are_not_current(tmp_path):
    """A mock file holding no mock (e.g. every photometric mock lacked detected bands) is redrawn, not reused."""
    joblib.dump(([], {}), tmp_path / sm.CALIBRATION_MOCK_FILE)
    assert not sm.calibration_mocks_are_current(str(tmp_path / sm.CALIBRATION_MOCK_FILE), "delayed")
    write_mock_file(tmp_path / "current.joblib", sm.MODELS["delayed"].training_labels)
    assert sm.calibration_mocks_are_current(str(tmp_path / "current.joblib"), "delayed")


def test_calibration_mocks_are_current_only_with_their_model_labels(tmp_path):
    """A mock file is current only in a folder of the model whose calibration labels its truths carry."""
    write_mock_file(tmp_path / "delayed.joblib", sm.MODELS["delayed"].training_labels)
    write_mock_file(tmp_path / "dblplaw.joblib", sm.MODELS["dblplaw"].training_labels)
    assert not sm.calibration_mocks_are_current(str(tmp_path / "delayed.joblib"), "dblplaw")
    assert sm.calibration_mocks_are_current(str(tmp_path / "dblplaw.joblib"), "dblplaw")
    assert not sm.calibration_mocks_are_current(str(tmp_path / "dblplaw.joblib"), "delayed")


def test_supersede_moves_old_mocks_their_jobs_and_the_run_logs(tmp_path):
    """16-bin mock truths lack log_age: the mock file, its calib_* jobs and the logs of the old fits move."""
    write_mock_file(tmp_path / sm.CALIBRATION_MOCK_FILE, ["mw_age", "Z_fin", "f_old", "log_mstar"])
    jobs = tmp_path / "jobs"
    jobs.mkdir()
    (jobs / "calib_000.npz").write_bytes(b"old mock job")
    (jobs / "123.npz").write_bytes(b"galaxy job")
    write_fit(tmp_path / "calib_000.h5", sfh_model=None)
    (tmp_path / "fit_log.csv").write_text("job_id,status\ncalib_000,fitted\n")
    (tmp_path / "fit_logs").mkdir()

    moved = sm.supersede_outdated_fits(str(tmp_path), "delayed")

    superseded = tmp_path / sm.SUPERSEDED_SUBDIRECTORY
    old_mock_job = os.path.join("jobs", "calib_000.npz")
    expected = ["calib_000.h5", sm.CALIBRATION_MOCK_FILE, old_mock_job, "fit_log.csv", "fit_logs"]
    assert sorted(moved) == sorted(expected)
    assert (superseded / sm.CALIBRATION_MOCK_FILE).exists()
    assert (superseded / "jobs" / "calib_000.npz").read_bytes() == b"old mock job"
    assert (jobs / "123.npz").exists()    # galaxy jobs do not depend on the model
    assert (superseded / "fit_log.csv").exists() and (superseded / "fit_logs").is_dir()


def test_supersede_keeps_current_mocks_their_jobs_and_logs(tmp_path):
    """A second run after the new mocks and fits exist moves nothing."""
    write_mock_file(tmp_path / sm.CALIBRATION_MOCK_FILE, sm.MODELS["delayed"].training_labels)
    jobs = tmp_path / "jobs"
    jobs.mkdir()
    (jobs / "calib_000.npz").write_bytes(b"new mock job")
    write_fit(tmp_path / "calib_000.h5")
    (tmp_path / "fit_log.csv").write_text("job_id,status\ncalib_000,fitted\n")

    moved = sm.supersede_outdated_fits(str(tmp_path), "delayed")

    assert moved == []
    assert (jobs / "calib_000.npz").exists() and (tmp_path / "fit_log.csv").exists()
    assert not (tmp_path / sm.SUPERSEDED_SUBDIRECTORY).exists()


def test_supersede_moves_jobs_of_mocks_whose_file_is_missing(tmp_path):
    """calib_* jobs without any mock file belong to no current mocks; a missing folder is no error."""
    jobs = tmp_path / "jobs"
    jobs.mkdir()
    (jobs / "calib_004.npz").write_bytes(b"orphan job")

    moved = sm.supersede_outdated_fits(str(tmp_path), "delayed")

    assert moved == [os.path.join("jobs", "calib_004.npz")]
    assert sm.supersede_outdated_fits(str(tmp_path / "missing"), "delayed") == []


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


# =============================================================================
# run_spectral_fits.py (run with LD_LIBRARY_PATH=../.venv/lib)
# =============================================================================
def run_script(*arguments):
    """Run run_spectral_fits.py in a subprocess, as a user would.

    Inputs: arguments: command-line options (str).
    Returns: the CompletedProcess, with stdout and stderr as text.
    """
    return subprocess.run([sys.executable, "run_spectral_fits.py", *arguments],
                          capture_output=True, text=True, timeout=600)


def test_run_script_fits_an_exported_job_and_logs_it(tmp_path):
    job, _ = mock_job("calib_smoke", np.random.default_rng(14))
    sm.export_job(job, str(tmp_path / "jobs"))

    result = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", str(tmp_path),
                        "--max-iterations", "100", "--n-live", "40", "--workers", "1")

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "calib_smoke.h5").exists()
    log = pd.read_csv(tmp_path / "fit_log.csv")
    assert list(log["job_id"]) == ["calib_smoke"]
    assert list(log["status"]) == ["fitted"]
    assert log["n_iterations"].iloc[0] <= 100


@pytest.mark.skipif(not multinest_available(), reason="pymultinest cannot import (LD_LIBRARY_PATH?)")
def test_run_script_fits_an_exported_photometry_job(tmp_path):
    sm.export_job(photometry_job("calib_000", np.random.default_rng(34)), str(tmp_path / "jobs"))

    result = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", str(tmp_path),
                        "--max-iterations", "200", "--n-live", "25", "--workers", "1")

    assert result.returncode == 0, result.stderr
    _, metadata = sm.load_fit(str(tmp_path / "calib_000.h5"))
    assert metadata["kind"] == "photometry"


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


def test_run_script_status_counts_done_and_pending_jobs(tmp_path):
    for job_id in ("calib_000", "calib_001", "galaxy_1"):
        sm.export_job(mock_job(job_id, np.random.default_rng(15))[0], str(tmp_path / "jobs"))
    (tmp_path / "calib_000.h5").write_bytes(b"")    # presence is what counts as done

    result = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", str(tmp_path), "--status")

    assert result.returncode == 0, result.stderr
    assert "3 jobs: 1 done, 2 pending" in result.stdout


def test_run_script_interleaves_calibration_mocks_and_galaxies():
    import run_spectral_fits
    order = run_spectral_fits.interleave_calibration_and_galaxies(
        ["calib_000", "calib_001", "calib_002", "galaxy_a", "galaxy_b"])
    assert order == ["calib_000", "galaxy_a", "calib_001", "galaxy_b", "calib_002"]


def test_run_script_rejects_too_long_output_path_and_bad_task_index(tmp_path):
    sm.export_job(mock_job("calib_000", np.random.default_rng(16))[0], str(tmp_path / "jobs"))

    long_output = str(tmp_path / ("x" * 100))    # MultiNest accepts under 100 characters in all
    too_long = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", long_output)
    bad_index = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", str(tmp_path), "--task-index", "5")

    assert too_long.returncode != 0 and "shorter relative --output-dir" in too_long.stderr
    assert bad_index.returncode != 0 and "--task-index 5 is outside 0..0" in bad_index.stderr


def cpuinfo_text(n_sockets, cores_per_socket, threads_per_core):
    """Build the text of a Linux /proc/cpuinfo for a machine with this topology (processors numbered from 0)."""
    blocks = []
    processor = 0
    for thread in range(threads_per_core):
        for socket_id in range(n_sockets):
            for core_id in range(cores_per_socket):
                blocks.append(f"processor\t: {processor}\nmodel name\t: test cpu\n"
                              f"physical id\t: {socket_id}\ncore id\t\t: {core_id}\n")
                processor += 1
    return "\n".join(blocks)


def test_physical_core_count_ignores_hardware_threads():
    import run_spectral_fits
    four_cores_two_threads = cpuinfo_text(1, 4, 2)                 # 8 processors
    two_sockets = cpuinfo_text(2, 3, 2)                            # 12 processors, 6 cores
    assert run_spectral_fits.count_physical_cores(four_cores_two_threads, set(range(8))) == 4
    assert run_spectral_fits.count_physical_cores(two_sockets, set(range(12))) == 6


def test_physical_core_count_follows_the_allowed_processors():
    import run_spectral_fits
    text = cpuinfo_text(1, 4, 2)    # processors 0-3 are the first threads, 4-7 their siblings
    assert run_spectral_fits.count_physical_cores(text, {0, 1, 4}) == 2     # cores 0 and 1 only
    assert run_spectral_fits.count_physical_cores(text, {0, 1, 2, 3}) == 4


def test_physical_core_count_is_none_without_topology_fields():
    import run_spectral_fits
    assert run_spectral_fits.count_physical_cores("processor\t: 0\nmodel name\t: arm\n", {0}) is None
    assert run_spectral_fits.count_physical_cores("", {0}) is None


def test_slurm_tasks_write_one_log_file_each(tmp_path):
    for job_id in ("calib_000", "calib_001"):
        sm.export_job(mock_job(job_id, np.random.default_rng(17))[0], str(tmp_path / "jobs"))

    for task_index in ("0", "1"):
        result = run_script("--jobs-dir", str(tmp_path / "jobs"), "--output-dir", str(tmp_path),
                            "--max-iterations", "100", "--n-live", "40", "--task-index", task_index)
        assert result.returncode == 0, result.stderr

    log_files = sorted(path.name for path in (tmp_path / "fit_logs").iterdir())
    assert log_files == ["calib_000.csv", "calib_001.csv"]
    assert not (tmp_path / "fit_log.csv").exists()
    for job_id in ("calib_000", "calib_001"):
        assert list(pd.read_csv(tmp_path / "fit_logs" / f"{job_id}.csv")["job_id"]) == [job_id]


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
