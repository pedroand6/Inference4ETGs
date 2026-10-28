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


def test_no_global_target_list_remains():
    """The targets depend on the SFH model: callers ask target_names(sfh_model)."""
    assert not hasattr(ml_data, "TARGET_NAMES")
    assert not hasattr(ml_data, "N_TARGETS")
