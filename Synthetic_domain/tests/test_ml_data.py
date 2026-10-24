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
    # Repeat the 20 models 10 times: each copy gets its own brightness and noise draw, so enough
    # galaxies end up with a precise r (the check below needs at least 5)
    repeated_models = pd.concat([simulated_models] * 10, ignore_index=True)
    observations = ml_data.simulate_observations(repeated_models, catalogue_r, catalogue_z, noise_table,
                                                 np.random.default_rng(31))
    targets = observations["targets"]
    flat = observations["flat"]
    log_mass_formed_index = ml_data.TARGET_NAMES.index("log_mass_formed")
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
