"""Tests for ml_networks.py (CPU, small). Run from Synthetic_domain/."""

import numpy as np
import pytest
import torch

import ml_data
import ml_networks

NOISE_SIGMA = 0.5        # toy problem: x = theta + N(0, NOISE_SIGMA), theta ~ N(0, 1)


def toy_problem(n_rows, generator):
    """Gaussian toy: two targets observed with noise in the first inputs (flat) and band channel 0."""
    theta = torch.randn(n_rows, 2, generator=generator)
    observed = theta + NOISE_SIGMA * torch.randn(n_rows, 2, generator=generator)
    flat = torch.zeros(n_rows, ml_data.N_INPUTS)
    flat[:, :2] = observed
    banded = torch.zeros(n_rows, 12, ml_data.N_BAND_CHANNELS)
    banded[:, :2, 0] = observed
    return {"flat": flat, "banded": banded, "targets": theta}


def small_config(embedding):
    """A small, fast configuration for tests."""
    config = dict(ml_networks.DEFAULT_CONFIG)
    config.update(embedding=embedding, n_targets=2, summary_size=16, hidden_size=64,
                  flow_transforms=3, flow_hidden=(64, 64), epochs=40, patience=40, batch_size=512,
                  learning_rate=3e-3, seed=0)
    return config


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


def test_save_load_round_trip_on_cpu(tmp_path):
    config = small_config("mlp")
    estimator = ml_networks.build_estimator(config)
    standardisers = {
        "flat": ml_data.Standardiser(np.zeros(ml_data.N_INPUTS), np.ones(ml_data.N_INPUTS)),
        "banded": ml_data.Standardiser(np.zeros(ml_data.N_BAND_CHANNELS), np.ones(ml_data.N_BAND_CHANNELS)),
        "targets": ml_data.Standardiser(np.array([1.0, 2.0]), np.array([0.5, 3.0])),
    }
    history = [{"epoch": 0, "train_loss": 1.0, "validation_loss": 1.1}]
    ml_networks.save_run(str(tmp_path), estimator, config, history, standardisers)
    loaded, loaded_config, loaded_standardisers = ml_networks.load_run(str(tmp_path), device="cpu")

    flat = np.random.default_rng(0).normal(size=(4, ml_data.N_INPUTS))
    banded = np.random.default_rng(1).normal(size=(4, 12, ml_data.N_BAND_CHANNELS))
    before = ml_networks.posterior_samples(estimator, standardisers, flat, banded, 50, seed=3)
    after = ml_networks.posterior_samples(loaded, loaded_standardisers, flat, banded, 50, seed=3)
    assert loaded_config["embedding"] == "mlp"
    assert before.shape == (4, 50, 2)
    assert np.allclose(before, after)
