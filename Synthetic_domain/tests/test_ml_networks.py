"""Tests for ml_networks.py (CPU, small). Run from Synthetic_domain/."""

import json
import os

import numpy as np
import pytest
import torch

import ml_data
import ml_networks

NOISE_SIGMA = 0.5        # toy problem: x = theta + N(0, NOISE_SIGMA), theta ~ N(0, 1)
TOY_TARGET_NAMES = ["toy_target_0", "toy_target_1"]
STORED_DELAYED_RUN = "ml_runs/mlp"   # trained before runs recorded their SFH model


def toy_problem(n_rows, generator):
    """Gaussian toy: two targets observed with noise in the first inputs (flat) and band channel 0."""
    theta = torch.randn(n_rows, 2, generator=generator)
    observed = theta + NOISE_SIGMA * torch.randn(n_rows, 2, generator=generator)
    flat = torch.zeros(n_rows, ml_data.N_INPUTS)
    flat[:, :2] = observed
    banded = torch.zeros(n_rows, 12, ml_data.N_BAND_CHANNELS)
    banded[:, :2, 0] = observed
    return {"flat": flat, "banded": banded, "targets": theta}


def small_config(embedding, sfh_model="delayed"):
    """A small, fast configuration for tests: the toy problem's two targets, recorded under sfh_model."""
    config = ml_networks.model_config(sfh_model)
    config.update(embedding=embedding, n_targets=2, target_names=TOY_TARGET_NAMES, summary_size=16,
                  hidden_size=64, flow_transforms=3, flow_hidden=(64, 64), epochs=40, patience=40,
                  batch_size=512, learning_rate=3e-3, seed=0)
    return config


def unit_standardisers(n_targets):
    """Standardisers that leave inputs and n_targets targets unchanged."""
    return {
        "flat": ml_data.Standardiser(np.zeros(ml_data.N_INPUTS), np.ones(ml_data.N_INPUTS)),
        "banded": ml_data.Standardiser(np.zeros(ml_data.N_BAND_CHANNELS), np.ones(ml_data.N_BAND_CHANNELS)),
        "targets": ml_data.Standardiser(np.zeros(n_targets), np.ones(n_targets)),
    }


def save_toy_run(directory, sfh_model="delayed"):
    """Write an untrained toy run (2 targets, TOY_TARGET_NAMES, recorded under sfh_model) with save_run."""
    config = small_config("mlp", sfh_model)
    history = [{"epoch": 0, "train_loss": 1.0, "validation_loss": 1.1}]
    ml_networks.save_run(str(directory), ml_networks.build_estimator(config), config, history,
                         unit_standardisers(2))


@pytest.mark.parametrize("embedding", ["mlp", "cnn"])
def test_flow_recovers_gaussian_posterior(embedding):
    generator = torch.Generator().manual_seed(0)
    config = small_config(embedding)
    estimator = ml_networks.build_estimator(config)
    validation = toy_problem(2000, generator)

    def make_epoch(epoch_index):
        return toy_problem(8192, generator)

    ml_networks.train_estimator(estimator, make_epoch, [validation], config, "cpu", log=lambda message: None)

    observed_value = 1.0
    flat = torch.zeros(1, ml_data.N_INPUTS)
    flat[0, :2] = observed_value
    banded = torch.zeros(1, 12, ml_data.N_BAND_CHANNELS)
    banded[0, :2, 0] = observed_value
    torch.manual_seed(1)
    samples = estimator.sample(flat, banded, 4000)[0].detach().numpy()
    posterior_mean = observed_value / (1 + NOISE_SIGMA**2)
    posterior_std = np.sqrt(NOISE_SIGMA**2 / (1 + NOISE_SIGMA**2))
    assert np.allclose(samples.mean(axis=0), posterior_mean, atol=0.12)
    assert np.allclose(samples.std(axis=0), posterior_std, rtol=0.2)


def test_estimator_samples_finite_for_fully_masked_input():
    estimator = ml_networks.build_estimator(small_config("cnn"))
    flat = torch.zeros(3, ml_data.N_INPUTS)
    banded = torch.zeros(3, 12, ml_data.N_BAND_CHANNELS)
    samples = estimator.sample(flat, banded, 10)
    assert samples.shape == (3, 10, 2)
    assert torch.isfinite(samples).all()


@pytest.mark.parametrize("sfh_model", ["delayed", "dblplaw"])
def test_model_config_records_the_model_and_its_targets(sfh_model):
    config = ml_networks.model_config(sfh_model)
    assert config["sfh_model"] == sfh_model
    assert config["target_names"] == ml_data.target_names(sfh_model)
    assert config["n_targets"] == len(ml_data.target_names(sfh_model))
    assert "sfh_model" not in ml_networks.DEFAULT_CONFIG        # the defaults belong to no model
    assert "target_names" not in ml_networks.DEFAULT_CONFIG


def test_dblplaw_estimator_samples_its_ten_targets():
    config = ml_networks.model_config("dblplaw")
    config.update(summary_size=16, hidden_size=64, flow_transforms=2, flow_hidden=(32, 32))
    estimator = ml_networks.build_estimator(config)
    samples = estimator.sample(torch.zeros(2, ml_data.N_INPUTS), torch.zeros(2, 12, ml_data.N_BAND_CHANNELS),
                               5)
    assert samples.shape == (2, 5, 10)


def test_save_load_round_trip_on_cpu(tmp_path):
    config = small_config("mlp")
    estimator = ml_networks.build_estimator(config)
    standardisers = unit_standardisers(2)
    standardisers["targets"] = ml_data.Standardiser(np.array([1.0, 2.0]), np.array([0.5, 3.0]))
    history = [{"epoch": 0, "train_loss": 1.0, "validation_loss": 1.1}]
    ml_networks.save_run(str(tmp_path), estimator, config, history, standardisers)
    loaded, loaded_config, loaded_standardisers = ml_networks.load_run(str(tmp_path), "delayed", device="cpu",
                                                                       target_names=TOY_TARGET_NAMES)

    flat = np.random.default_rng(0).normal(size=(4, ml_data.N_INPUTS))
    banded = np.random.default_rng(1).normal(size=(4, 12, ml_data.N_BAND_CHANNELS))
    before = ml_networks.posterior_samples(estimator, standardisers, flat, banded, 50, seed=3)
    after = ml_networks.posterior_samples(loaded, loaded_standardisers, flat, banded, 50, seed=3)
    assert loaded_config["embedding"] == "mlp"
    assert loaded_config["sfh_model"] == "delayed"
    assert loaded_config["target_names"] == TOY_TARGET_NAMES
    assert before.shape == (4, 50, 2)
    assert np.allclose(before, after)


def test_load_run_of_a_checkpoint_without_target_names_raises(tmp_path):
    save_toy_run(tmp_path)
    config_path = tmp_path / ml_networks.CONFIG_FILE
    stored = json.loads(config_path.read_text())
    del stored["target_names"]                       # a checkpoint written before target names were recorded
    config_path.write_text(json.dumps(stored))
    with pytest.raises(ValueError, match="run was trained for targets"):
        ml_networks.load_run(str(tmp_path), "delayed", target_names=TOY_TARGET_NAMES)


def test_load_run_of_other_targets_raises(tmp_path):
    save_toy_run(tmp_path)
    with pytest.raises(ValueError, match="run was trained for targets"):
        ml_networks.load_run(str(tmp_path), "delayed")          # expects the delayed-tau targets


def test_load_run_of_another_sfh_model_raises(tmp_path):
    save_toy_run(tmp_path, sfh_model="dblplaw")
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(str(tmp_path), "delayed", target_names=TOY_TARGET_NAMES)


def test_load_run_without_sfh_model_counts_as_delayed(tmp_path):
    save_toy_run(tmp_path)
    config_path = tmp_path / ml_networks.CONFIG_FILE
    stored = json.loads(config_path.read_text())
    del stored["sfh_model"]                          # a run saved before runs recorded their model
    config_path.write_text(json.dumps(stored))
    _, config, _ = ml_networks.load_run(str(tmp_path), "delayed", target_names=TOY_TARGET_NAMES)
    assert config["sfh_model"] == "delayed"
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(str(tmp_path), "dblplaw", target_names=TOY_TARGET_NAMES)


@pytest.mark.skipif(not os.path.exists(os.path.join(STORED_DELAYED_RUN, ml_networks.MODEL_FILE)),
                    reason="no ml_runs/mlp weights here (model.pt is git-ignored)")
def test_the_stored_delayed_run_loads_only_as_delayed():
    _, config, _ = ml_networks.load_run(STORED_DELAYED_RUN, "delayed")
    assert config["sfh_model"] == "delayed"
    assert config["target_names"] == ml_data.target_names("delayed")
    with pytest.raises(ValueError, match="SFH model"):
        ml_networks.load_run(STORED_DELAYED_RUN, "dblplaw")


def test_save_run_rejects_a_config_without_sfh_model(tmp_path):
    config = small_config("mlp")
    del config["sfh_model"]
    with pytest.raises(ValueError, match="sfh_model"):
        ml_networks.save_run(str(tmp_path), ml_networks.build_estimator(small_config("mlp")), config, [], {})


def test_save_run_rejects_target_names_that_do_not_match_n_targets(tmp_path):
    config = small_config("mlp")
    config["target_names"] = ml_data.target_names("delayed")   # 9 names for a 2-target flow
    with pytest.raises(ValueError, match="target_names"):
        ml_networks.save_run(str(tmp_path), ml_networks.build_estimator(small_config("mlp")), config, [], {})
