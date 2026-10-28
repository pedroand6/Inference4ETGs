"""Tests for sfh_models.py. Run from Synthetic_domain/: ../.venv/bin/python -m pytest tests/test_sfh_models.py -q"""

import subprocess
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest
from scipy.stats import qmc

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
    simulated = sm.simulate_models(first_unit_vectors(TRAINING_SEED, N_TRAINING_MODELS))
    stored = pd.read_parquet(TRAINING_FILE).head(N_REGRESSION_ROWS)
    columns = sfh.DELAYED.training_labels + sm.BANDS

    assert list(stored.columns) == columns
    np.testing.assert_allclose(simulated[columns].to_numpy(), stored[columns].to_numpy(),
                               rtol=REGRESSION_RTOL, atol=REGRESSION_ATOL)


@pytest.mark.skipif(not parquet_has_column(VALIDATION_FILE, "age"), reason="no delayed-tau validation set here")
def test_delayed_path_regenerates_the_labels_of_the_stored_validation_set():
    """The stored labels are rescaled to an observed brightness: undo that with the stored magnitude offset."""
    simulated = sm.simulate_models(first_unit_vectors(VALIDATION_SEED, N_HELD_OUT_MODELS))
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
