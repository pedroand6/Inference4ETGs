"""Tests for sfh_models.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest tests/test_sfh_models.py -q"""

import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy.stats import qmc

import bagpipes as pipes
import sfh_models as sfh
import synthetic_model as sm

# Sizes and seeds of the stored delayed-tau sets, as training_set.ipynb draws them
TRAINING_FILE = "training_set.parquet"
VALIDATION_FILE = "synthetic_validation.parquet"
N_TRAINING_MODELS = 2**19
N_HELD_OUT_MODELS = 2**14
TRAINING_SEED = 0
VALIDATION_SEED = 1
N_REGRESSION_ROWS = 4            # first rows of each stored set that are simulated again (Bagpipes: ~0.3 s each)
REGRESSION_RTOL = 1e-10          # "exactly", up to floating-point rounding
REGRESSION_ATOL = 1e-12


def parquet_has_column(path, column):
    """True when the parquet file exists and its schema holds the column (read without loading the data)."""
    try:
        return column in pq.read_schema(path).names
    except OSError:
        return False


def first_unit_vectors(seed, n_models):
    """First N_REGRESSION_ROWS Latin hypercube vectors of a delayed-tau set, drawn as training_set.ipynb does.

    Why the full draw: the Latin hypercube permutes all n_models strata, so its first rows depend on n_models.
    Inputs: seed: Latin hypercube seed; n_models: size of the stored set (count).
    Returns: array (N_REGRESSION_ROWS, len(DELAYED_DIMS)) of unit-cube values (dimensionless).
    """
    sampler = qmc.LatinHypercube(len(sfh.DELAYED_DIMS), seed=seed)
    return sampler.random(n_models)[:N_REGRESSION_ROWS]


# =============================================================================
# Delayed-tau unchanged by the move to sfh_models.py
# =============================================================================
@pytest.mark.skipif(not parquet_has_column(TRAINING_FILE, "age"), reason="no delayed-tau training set here")
def test_delayed_path_regenerates_the_first_rows_of_the_stored_training_set():
    simulated = sm.simulate_models(first_unit_vectors(TRAINING_SEED, N_TRAINING_MODELS), "delayed")
    stored = pd.read_parquet(TRAINING_FILE).head(N_REGRESSION_ROWS)
    columns = sfh.DELAYED.training_labels + sm.BANDS

    assert list(stored.columns) == columns
    np.testing.assert_allclose(simulated[columns].to_numpy(), stored[columns].to_numpy(),
                               rtol=REGRESSION_RTOL, atol=REGRESSION_ATOL)


@pytest.mark.skipif(not parquet_has_column(VALIDATION_FILE, "age"), reason="no delayed-tau validation set here")
def test_delayed_path_regenerates_the_labels_of_the_stored_validation_set():
    """The stored labels are rescaled to an observed brightness: undo that with the stored magnitude offset."""
    simulated = sm.simulate_models(first_unit_vectors(VALIDATION_SEED, N_HELD_OUT_MODELS), "delayed")
    stored = pd.read_parquet(VALIDATION_FILE).head(N_REGRESSION_ROWS)
    offsets = stored["magnitude_offset"].to_numpy()
    assert np.isfinite(offsets).all()
    unscaled = sm.rescale_mass_labels(stored, -offsets)
    columns = sfh.DELAYED.training_labels

    np.testing.assert_allclose(simulated[columns].to_numpy(), unscaled[columns].to_numpy(),
                               rtol=REGRESSION_RTOL, atol=REGRESSION_ATOL)


def test_synthetic_model_reexports_the_shared_prior_and_the_delayed_aliases():
    assert sm.Z_RANGE == sfh.Z_RANGE and sm.COSMOLOGY is sfh.COSMOLOGY
    assert sm.maximum_age is sfh.maximum_age and sm.delayed_tau_labels is sfh.delayed_tau_labels
    assert sm.GYR_TO_YEARS == sfh.GYR_TO_YEARS == 1e9
    delayed = sm.MODELS["delayed"]
    assert sm.DIMS == delayed.dims
    assert sm.PHOTOMETRIC_DIMS == delayed.photometric_dims
    assert sm.DRAWN_LABELS == delayed.drawn_labels
    assert sm.DERIVED_LABELS == delayed.derived_labels
    assert sm.TRAINING_LABELS == delayed.training_labels
    assert sm.CALIBRATION_LABELS == delayed.calibration_labels
    assert sm.SFH_MODEL == "delayed"
    assert sm.components is sfh.delayed_components


def test_sfh_models_does_not_import_synthetic_model():
    """synthetic_model imports sfh_models; the reverse would be a circular import."""
    command = "import sys, sfh_models; assert 'synthetic_model' not in sys.modules"
    subprocess.run([sys.executable, "-c", command], check=True)


# =============================================================================
# Shared prior and the registry
# =============================================================================
def test_dust_and_nebular_follow_the_shared_prior():
    params = {"av_ism": 0.5, "log_tau_bc": 0.5, "logU": 0.5}
    bagpipes_parts, drawn = sfh.dust_and_nebular(params)

    assert list(bagpipes_parts) == ["t_bc", "dust", "nebular"]
    assert np.isclose(drawn["av_ism"], np.mean(sfh.AV_ISM_RANGE))
    assert np.isclose(np.log10(drawn["tau_bc"]), np.mean(sfh.LOG_TAU_BC_RANGE))
    assert np.isclose(drawn["logU"], np.mean(sfh.LOGU_RANGE))
    dust = bagpipes_parts["dust"]
    assert dust["type"] == "CF00" and dust["n"] == sfh.DUST_SLOPE and dust["Av"] == drawn["av_ism"]
    assert np.isclose(dust["eta"], 1 + sfh.TAU_TO_AV * drawn["tau_bc"] / drawn["av_ism"])
    assert bagpipes_parts["nebular"] == {"logU": drawn["logU"]}
    assert bagpipes_parts["t_bc"] == sfh.BIRTH_CLOUD_AGE


@pytest.mark.parametrize("name", list(sfh.MODELS))
def test_registry_entries_are_consistent(name):
    model = sfh.MODELS[name]
    u = np.full(len(model.dims), 0.5)
    _, drawn = model.components(u)

    assert model.name == name
    assert model.dims[0] == "z"
    assert model.photometric_dims == model.dims[1:]
    assert list(drawn) == model.drawn_labels
    assert model.training_labels == model.drawn_labels + model.derived_labels
    for log_label, drawn_label in model.log_labels.items():
        assert log_label in model.derived_labels and drawn_label in model.drawn_labels
    assert set(sfh.SHARED_LABELS) <= set(model.training_labels)
    assert set(sfh.HISTORY_LABELS) <= set(model.derived_labels)
    assert model.target_names[-1] == "redshift"
    assert set(model.calibration_labels) <= set(model.training_labels)


def test_model_file_keeps_the_delayed_tau_names():
    for name in ("training_set.parquet", "photometric_fits", "ml_runs/mlp", "ml_sfh_catalogue.parquet"):
        assert sfh.model_file(name, "delayed") == name


# =============================================================================
# Double power law: prior
# =============================================================================
LOW_REDSHIFT = 0.002             # Bagpipes needs z > 0 to build a star formation history


def test_dblplaw_draws_stay_inside_bounds_with_the_turnover_after_the_big_bang():
    rng = np.random.default_rng(0)
    for u in rng.random((500, len(sfh.DBLPLAW_DIMS))):
        bagpipes_components, drawn = sfh.dblplaw_components(u)
        redshift = drawn["redshift"]
        assert sfh.Z_RANGE[0] <= redshift <= sfh.Z_RANGE[1]
        assert sfh.TURNOVER_LOOKBACK_MIN <= drawn["turnover_lookback"] <= sfh.maximum_age(redshift)
        assert np.isclose(drawn["tau"], sfh.COSMOLOGY.age(redshift).value - drawn["turnover_lookback"])
        assert drawn["tau"] > 0
        # Bagpipes flags tau > 1.5 x the age of the Universe as unphysical
        assert drawn["tau"] < 1.5 * sfh.bagpipes_universe_age(redshift)
        assert sfh.ALPHA_RANGE[0] <= drawn["alpha"] <= sfh.ALPHA_RANGE[1]
        assert sfh.BETA_RANGE[0] <= drawn["beta"] <= sfh.BETA_RANGE[1]
        assert sfh.METALLICITY_RANGE[0] <= drawn["metallicity"] <= sfh.METALLICITY_RANGE[1]
        assert sfh.AV_ISM_RANGE[0] <= drawn["av_ism"] <= sfh.AV_ISM_RANGE[1]
        assert 10 ** sfh.LOG_TAU_BC_RANGE[0] <= drawn["tau_bc"] <= 10 ** sfh.LOG_TAU_BC_RANGE[1]
        assert sfh.LOGU_RANGE[0] <= drawn["logU"] <= sfh.LOGU_RANGE[1]
        assert drawn["log_mass_formed"] == sfh.LOG_MASS_FORMED


def test_dblplaw_turnover_at_the_top_of_the_cube_never_precedes_the_big_bang():
    cube = np.full(len(sfh.DBLPLAW_DIMS), 1.0)
    for fixed_redshift in (None, sfh.Z_RANGE[1], 0.3):
        _, drawn = sfh.dblplaw_components(cube, fixed_redshift=fixed_redshift)
        assert drawn["turnover_lookback"] <= sfh.maximum_age(drawn["redshift"])
        assert drawn["tau"] >= sfh.UNIVERSE_AGE_MARGIN - 1e-9


def test_uniform_cube_gives_log_uniform_turnover_lookback_alpha_beta_and_metallicity():
    cube = np.full(len(sfh.DBLPLAW_DIMS), 0.5)
    _, drawn = sfh.dblplaw_components(cube)
    log_lookback_range = [np.log10(sfh.TURNOVER_LOOKBACK_MIN), np.log10(sfh.maximum_age(drawn["redshift"]))]
    assert np.isclose(np.log10(drawn["turnover_lookback"]), np.mean(log_lookback_range))
    assert np.isclose(np.log10(drawn["alpha"]), np.mean(np.log10(sfh.ALPHA_RANGE)))
    assert np.isclose(np.log10(drawn["beta"]), np.mean(np.log10(sfh.BETA_RANGE)))
    assert np.isclose(np.log10(drawn["metallicity"]), np.mean(np.log10(sfh.METALLICITY_RANGE)))


@pytest.mark.parametrize("label, value_range", [("alpha", sfh.ALPHA_RANGE), ("beta", sfh.BETA_RANGE)])
def test_latin_hypercube_is_stratified_in_log_alpha_and_log_beta(label, value_range):
    """Each of n equal slices of the log range holds exactly one of n Latin hypercube draws."""
    n_draws = 50
    cube = qmc.LatinHypercube(len(sfh.DBLPLAW_DIMS), seed=4).random(n_draws)
    log_values = np.zeros(n_draws)
    for row_index, u in enumerate(cube):
        _, drawn = sfh.dblplaw_components(u)
        log_values[row_index] = np.log10(drawn[label])
    log_range = np.log10(value_range)
    slices = np.floor((log_values - log_range[0]) / np.diff(log_range) * n_draws).astype(int)
    assert np.array_equal(np.sort(slices), np.arange(n_draws))


def test_dblplaw_components_hold_one_dblplaw_component_and_the_shared_dust():
    u = np.random.default_rng(5).random(len(sfh.DBLPLAW_DIMS))
    bagpipes_components, drawn = sfh.dblplaw_components(u)

    assert list(drawn) == sfh.DBLPLAW.drawn_labels
    assert bagpipes_components["dblplaw"] == {"tau": drawn["tau"], "alpha": drawn["alpha"], "beta": drawn["beta"],
                                              "massformed": sfh.LOG_MASS_FORMED,
                                              "metallicity": drawn["metallicity"]}
    params = dict(zip(sfh.DBLPLAW_DIMS, u))
    shared_parts, _ = sfh.dust_and_nebular(params)
    for key, value in shared_parts.items():
        assert bagpipes_components[key] == value
    assert "delayed" not in bagpipes_components


def test_dblplaw_fixed_redshift_overrides_the_drawn_redshift():
    bagpipes_components, drawn = sfh.dblplaw_components(np.full(len(sfh.DBLPLAW_DIMS), 0.5), fixed_redshift=0.123)
    assert bagpipes_components["redshift"] == 0.123 and drawn["redshift"] == 0.123
    assert np.isclose(drawn["tau"], sfh.COSMOLOGY.age(0.123).value - drawn["turnover_lookback"])


@pytest.mark.parametrize("cube_value", [0.0, 0.5, 0.999999])
def test_dblplaw_prior_corners_are_physical_for_bagpipes(cube_value):
    bagpipes_components, _ = sfh.dblplaw_components(np.full(len(sfh.DBLPLAW_DIMS), cube_value))
    history = pipes.models.star_formation_history(bagpipes_components)
    assert not history.unphysical
    assert np.isfinite(history.stellar_mass)


# =============================================================================
# Double power law: labels summed on Bagpipes' grid
# =============================================================================
def bagpipes_dblplaw_sfh(tau, alpha, beta, redshift):
    """Bagpipes SFH of one dblplaw component (10^10 M_sun formed, Z_sun).

    Inputs: tau [Gyr since the Big Bang]; alpha, beta, redshift (dimensionless).
    Returns: bagpipes star_formation_history object.
    """
    component = {"tau": tau, "alpha": alpha, "beta": beta, "massformed": 10.0, "metallicity": 1.0}
    return pipes.models.star_formation_history({"redshift": redshift, "dblplaw": component})


def test_bagpipes_grid_is_rebuilt_exactly():
    history = bagpipes_dblplaw_sfh(5.0, 2.0, 3.0, 0.1)
    np.testing.assert_allclose(sfh.BAGPIPES_AGES * 1e9, history.ages, rtol=1e-12, atol=0)
    np.testing.assert_allclose(sfh.BAGPIPES_AGE_WIDTHS * 1e9, history.age_widths, rtol=1e-12, atol=0)
    assert np.isclose(sfh.bagpipes_universe_age(0.1) * 1e9, history.age_of_universe, rtol=1e-12)


# (alpha, beta, redshift, where the turnover tau lies): a sharp peak (alpha = beta = 1000, where the SFR maximum is
# at tau), a flat SFH, early and recent turnovers, asymmetric slopes. "oldest" puts tau at UNIVERSE_AGE_MARGIN.
DBLPLAW_CASES = [
    (1000.0, 1000.0, LOW_REDSHIFT, "early"),
    (1000.0, 1000.0, LOW_REDSHIFT, "recent"),
    (1000.0, 1000.0, LOW_REDSHIFT, "oldest"),
    (1000.0, 1000.0, 0.7, "oldest"),
    (0.01, 0.01, 0.3, "middle"),
    (0.5, 10.0, LOW_REDSHIFT, "early"),
    (10.0, 0.5, 0.7, "recent"),
    (0.01, 1000.0, LOW_REDSHIFT, "recent"),
]
TURNOVER_LOOKBACKS = {"recent": 0.06, "middle": 3.0}     # [Gyr]; "early" is 1 Gyr after the Big Bang
FRACTION_TOLERANCE = 0.01        # absolute, on fractions
DEX_PER_PERCENT = 0.0043         # log10(1.01): a 1% difference in SFR [dex]


@pytest.mark.parametrize("alpha, beta, redshift, where", DBLPLAW_CASES)
def test_dblplaw_labels_match_bagpipes_sfh(alpha, beta, redshift, where):
    universe_age = sfh.COSMOLOGY.age(redshift).value
    if where == "early":
        turnover_lookback = universe_age - 1.0
    elif where == "oldest":
        turnover_lookback = sfh.maximum_age(redshift)
    else:
        turnover_lookback = TURNOVER_LOOKBACKS[where]
    tau = universe_age - turnover_lookback
    history = bagpipes_dblplaw_sfh(tau, alpha, beta, redshift)
    mass = history.sfh * history.age_widths
    ages = history.ages / 1e9
    labels = sfh.dblplaw_labels(tau, alpha, beta, redshift, 10.0, history.stellar_mass)

    assert np.isclose(labels["mw_age"], history.mass_weighted_age, rtol=0.01)
    assert abs(labels["f_young"] - mass[ages < 0.1].sum() / mass.sum()) < FRACTION_TOLERANCE
    assert abs(labels["f_mid"] - mass[(ages >= 0.1) & (ages < 1.0)].sum() / mass.sum()) < FRACTION_TOLERANCE
    assert abs(labels["f_old"] - mass[ages >= 1.0].sum() / mass.sum()) < FRACTION_TOLERANCE
    assert np.isclose(labels["f_young"] + labels["f_mid"] + labels["f_old"], 1.0)

    recent_sfr = mass[ages < sfh.SFR_TIMESCALE].sum() / (sfh.SFR_TIMESCALE * 1e9)
    sfr_floor = sfh.SFR_FLOOR * 10**history.stellar_mass
    assert abs(labels["log_sfr"] - np.log10(max(recent_sfr, sfr_floor))) < DEX_PER_PERCENT
    if history.sfr > sfr_floor:
        assert abs(labels["log_sfr"] - np.log10(history.sfr)) < DEX_PER_PERCENT
    assert np.isclose(labels["log_ssfr"], labels["log_sfr"] - history.stellar_mass)


def test_dblplaw_labels_are_finite_in_the_prior_corners():
    corners = []
    for redshift in (LOW_REDSHIFT, sfh.Z_RANGE[1]):
        for turnover_lookback in (sfh.TURNOVER_LOOKBACK_MIN, sfh.maximum_age(redshift)):
            for alpha in sfh.ALPHA_RANGE:
                for beta in sfh.BETA_RANGE:
                    tau = sfh.COSMOLOGY.age(redshift).value - turnover_lookback
                    corners.append((tau, alpha, beta, redshift))
    tau, alpha, beta, redshift = np.array(corners).T
    labels = sfh.dblplaw_labels(tau, alpha, beta, redshift, np.full(len(tau), 10.0), np.full(len(tau), 9.8))

    assert all(np.isfinite(value).all() for value in labels.values())
    for fraction in ("f_young", "f_mid", "f_old"):
        assert ((labels[fraction] >= 0) & (labels[fraction] <= 1)).all(), fraction
    assert np.allclose(labels["f_young"] + labels["f_mid"] + labels["f_old"], 1.0)
    assert ((labels["mw_age"] > 0) & (labels["mw_age"] <= sfh.bagpipes_universe_age(redshift))).all()


SHARP_TURNOVER_AGE_TOLERANCE = 0.08  # [Gyr] one step of Bagpipes' grid near a lookback time of 13 Gyr


@pytest.mark.parametrize("redshift", [LOW_REDSHIFT, sfh.Z_RANGE[1]])
def test_sharp_turnover_at_the_oldest_lookback_gives_finite_labels(redshift):
    """alpha = beta = 1000 with tau = 0.1 Gyr, the corner the ML decoding reaches when it clips to the prior.

    Why: the SFR is a spike narrower than the grid step there, so every grid SFR is tiny (about e^-300).
    Normalising by each sample's largest grid SFR keeps the sums finite: no 0 / 0.
    """
    turnover_lookback = sfh.maximum_age(redshift)
    tau = sfh.COSMOLOGY.age(redshift).value - turnover_lookback
    assert np.isclose(tau, sfh.UNIVERSE_AGE_MARGIN)
    labels = sfh.dblplaw_labels(tau, 1000.0, 1000.0, redshift, 10.0, 9.8)

    assert all(np.isfinite(value) for value in labels.values())
    assert np.isclose(labels["f_young"] + labels["f_mid"] + labels["f_old"], 1.0)
    # All the mass forms within a grid step of the turnover, more than 1 Gyr ago
    assert np.isclose(labels["f_old"], 1.0)
    assert abs(labels["mw_age"] - turnover_lookback) < SHARP_TURNOVER_AGE_TOLERANCE
    # No star formation in the last 100 Myr, so log_ssfr sits on the floor
    assert np.isclose(labels["log_ssfr"], np.log10(sfh.SFR_FLOOR))


def test_dblplaw_labels_do_not_depend_on_the_chunk_size(monkeypatch):
    """Blocks of 3 samples stop at their own oldest Universe, so they use fewer grid ages than one block of 10.

    This also checks that the grid ages each block skips hold no mass.
    """
    rng = np.random.default_rng(7)
    n_samples = 10
    redshift = rng.uniform(*sfh.Z_RANGE, n_samples)
    tau = sfh.COSMOLOGY.age(redshift).value - rng.uniform(0.1, 5.0, n_samples)
    alpha = 10 ** rng.uniform(-2, 3, n_samples)
    beta = 10 ** rng.uniform(-2, 3, n_samples)
    one_block = sfh.dblplaw_labels(tau, alpha, beta, redshift, np.full(n_samples, 10.0), np.full(n_samples, 9.8))
    monkeypatch.setattr(sfh, "DBLPLAW_CHUNK", 3)
    small_blocks = sfh.dblplaw_labels(tau, alpha, beta, redshift, np.full(n_samples, 10.0),
                                      np.full(n_samples, 9.8))
    for label in sfh.HISTORY_LABELS:
        np.testing.assert_allclose(small_blocks[label], one_block[label], rtol=1e-12, atol=1e-15)


def test_dblplaw_labels_keep_the_input_shape():
    scalar = sfh.dblplaw_labels(5.0, 2.0, 3.0, 0.1, 10.0, 9.8)
    assert all(np.shape(value) == () for value in scalar.values())
    grid = sfh.dblplaw_labels(np.full((2, 3), 5.0), 2.0, 3.0, 0.1, 10.0, 9.8)
    assert all(value.shape == (2, 3) for value in grid.values())


def test_dblplaw_analytic_labels_read_the_drawn_columns():
    _, drawn = sfh.dblplaw_components(np.full(len(sfh.DBLPLAW_DIMS), 0.5))
    drawn_columns = {}
    for label, value in drawn.items():
        drawn_columns[label] = np.array([value])
    from_registry = sfh.DBLPLAW.analytic_labels(drawn_columns, np.array([9.8]))
    direct = sfh.dblplaw_labels(drawn["tau"], drawn["alpha"], drawn["beta"], drawn["redshift"],
                                drawn["log_mass_formed"], 9.8)
    assert list(from_registry) == sfh.HISTORY_LABELS
    for label in sfh.HISTORY_LABELS:
        assert np.isclose(from_registry[label][0], direct[label])


def test_model_file_adds_the_dblplaw_suffix():
    assert sfh.model_file("training_set.parquet", "dblplaw") == "training_set_dblplaw.parquet"
    assert sfh.model_file("photometric_fits", "dblplaw") == "photometric_fits_dblplaw"
    assert sfh.model_file("ml_runs/mlp", "dblplaw") == "ml_runs/mlp_dblplaw"
    assert sfh.model_file("ml_sfh_catalogue.parquet", "dblplaw") == "ml_sfh_catalogue_dblplaw.parquet"
