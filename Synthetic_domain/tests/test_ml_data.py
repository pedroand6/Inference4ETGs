"""Tests for ml_data.py. Run from Synthetic_domain/: LD_LIBRARY_PATH=../.venv/lib ../.venv/bin/python -m pytest -v"""

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy.stats import qmc

import ml_data
import synthetic_model as sm

VALIDATION_FILE = "synthetic_validation.parquet"
ROUND_TRIP_TOLERANCE = 1e-9      # the decode uses the same analytic functions as the labels: exact


def parquet_has_column(path, column):
    """True when the parquet file exists and its schema holds the column (read without loading the data)."""
    try:
        return column in pq.read_schema(path).names
    except OSError:
        return False


@pytest.fixture(scope="module")
def simulated_models():
    """Two hundred prior draws (Latin hypercube) with labels and noiseless magnitudes."""
    unit_vectors = qmc.LatinHypercube(d=len(sm.DIMS), seed=30).random(200)
    return sm.simulate_models(unit_vectors)


def theta_rows(rows):
    """Flow targets (n, N_TARGETS) from dicts keyed by TARGET_NAMES; missing targets take typical values."""
    typical = {"log_age": 0.0, "log_tau": 0.0, "log_metallicity": 0.0, "av_ism": 0.5, "log_tau_bc": 0.0,
               "logU": -2.5, "log_mass_formed": 10.0, "log_mstar": 9.8, "redshift": 0.1}
    theta = np.zeros((len(rows), ml_data.N_TARGETS))
    for row_index, row in enumerate(rows):
        values = dict(typical, **row)
        for target_index, name in enumerate(ml_data.TARGET_NAMES):
            theta[row_index, target_index] = values[name]
    return theta


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


def test_targets_are_the_delayed_tau_parameters_with_redshift_last():
    assert ml_data.TARGET_NAMES == ["log_age", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU",
                                    "log_mass_formed", "log_mstar", "redshift"]
    assert ml_data.N_TARGETS == 9


def test_encode_targets_takes_logs_of_age_tau_metallicity_and_tau_bc(simulated_models):
    theta = ml_data.encode_targets(simulated_models)
    assert theta.shape == (len(simulated_models), ml_data.N_TARGETS)
    expected = {
        "log_age": np.log10(simulated_models["age"]),
        "log_tau": np.log10(simulated_models["tau"]),
        "log_metallicity": np.log10(simulated_models["metallicity"]),
        "av_ism": simulated_models["av_ism"],
        "log_tau_bc": np.log10(simulated_models["tau_bc"]),
        "logU": simulated_models["logU"],
        "log_mass_formed": simulated_models["log_mass_formed"],
        "log_mstar": simulated_models["log_mstar"],
        "redshift": simulated_models["redshift"],
    }
    for target_index, name in enumerate(ml_data.TARGET_NAMES):
        assert np.array_equal(theta[:, target_index], expected[name].to_numpy(dtype=float)), name


def test_encode_decode_round_trip_is_exact_for_every_training_label(simulated_models):
    decoded = ml_data.decode_samples(ml_data.encode_targets(simulated_models))
    assert list(decoded.columns) == sm.TRAINING_LABELS
    for label in sm.TRAINING_LABELS:
        expected = simulated_models[label].to_numpy(dtype=float)
        assert np.allclose(decoded[label], expected, rtol=0, atol=ROUND_TRIP_TOLERANCE), label


def test_decode_clips_redshift_to_prior(simulated_models):
    theta = ml_data.encode_targets(simulated_models.iloc[:2])
    theta[0, -1] = -0.3
    theta[1, -1] = 1.5
    decoded = ml_data.decode_samples(theta)
    assert decoded["redshift"].between(*sm.Z_RANGE).all()


def test_decode_clips_age_to_the_oldest_age_at_the_sample_redshift():
    log_age_too_old = np.log10(10.0)        # older than t_max at z = 0.7 (7.07 Gyr), not at z = 0.1 (12.07)
    theta = theta_rows([
        {"log_age": log_age_too_old, "redshift": 0.7},
        {"log_age": log_age_too_old, "redshift": 0.1},
        {"log_age": np.log10(sm.maximum_age(0.7)), "redshift": 0.7},
        {"log_age": np.log10(0.01), "redshift": 0.1},
        {"log_age": log_age_too_old, "redshift": 1.5},
    ])
    decoded = ml_data.decode_samples(theta)
    assert decoded["age"][0] <= sm.maximum_age(0.7)
    assert np.isclose(decoded["age"][0], sm.maximum_age(0.7), rtol=1e-12)
    assert np.isclose(decoded["age"][1], 10.0, rtol=1e-12)
    assert np.isclose(decoded["age"][2], sm.maximum_age(0.7), rtol=1e-12)
    assert decoded["age"][3] == sm.AGE_MIN
    assert np.isclose(decoded["age"][4], sm.maximum_age(sm.Z_RANGE[1]), rtol=1e-12)
    assert np.allclose(decoded["log_age"], np.log10(decoded["age"]))


def test_decode_clips_tau_and_metallicity_to_their_ranges():
    theta = theta_rows([
        {"log_tau": -5.0, "log_metallicity": -6.0},
        {"log_tau": 5.0, "log_metallicity": 3.0},
    ])
    decoded = ml_data.decode_samples(theta)
    assert list(decoded["tau"]) == list(sm.TAU_RANGE)
    assert list(decoded["metallicity"]) == list(sm.METALLICITY_RANGE)
    assert np.allclose(decoded["log_tau"], np.log10(sm.TAU_RANGE))
    assert np.allclose(decoded["log_metallicity"], np.log10(sm.METALLICITY_RANGE))


def test_decode_extreme_samples_stay_finite():
    theta = theta_rows([
        {"log_age": 50.0, "log_tau": -50.0, "redshift": 0.7},
        {"log_age": -50.0, "log_tau": 50.0, "redshift": 0.002},
        {"log_age": 50.0, "log_tau": 50.0, "log_metallicity": 40.0},
        {"log_age": -50.0, "log_tau": -50.0, "log_metallicity": -40.0},
    ])
    with np.errstate(divide="raise", over="raise", invalid="raise"):
        decoded = ml_data.decode_samples(theta)
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


def test_simulate_observations_shifts_mass_and_draws_redshift(simulated_models, noise_and_catalogue):
    noise_table, catalogue_r, catalogue_z = noise_and_catalogue
    # Repeat the models 10 times: each copy gets its own brightness and noise draw, so enough
    # galaxies end up with a precise r (the check below needs at least 5)
    repeated_models = pd.concat([simulated_models] * 10, ignore_index=True)
    observations = ml_data.simulate_observations(repeated_models, catalogue_r, catalogue_z, noise_table,
                                                 np.random.default_rng(31))
    targets = observations["targets"]
    assert targets.shape == (len(repeated_models), ml_data.N_TARGETS)
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


@pytest.mark.skipif(not parquet_has_column(VALIDATION_FILE, "age"),
                    reason="synthetic sets not yet regenerated on the PC")
def test_validation_observations_modes():
    held_out = pd.read_parquet(VALIDATION_FILE).iloc[:200]
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
