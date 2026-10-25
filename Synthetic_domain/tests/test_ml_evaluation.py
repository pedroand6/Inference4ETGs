"""Tests for ml_evaluation.py. Run from Synthetic_domain/."""

import numpy as np

import ml_evaluation as ev


def test_calibrated_posteriors_give_a_diagonal_coverage_curve():
    rng = np.random.default_rng(40)
    n_galaxies, n_samples = 4000, 400
    centres = rng.normal(0, 1, n_galaxies)
    truths = centres + rng.normal(0, 1, n_galaxies)         # truth drawn from the posterior itself
    samples = centres[:, None] + rng.normal(0, 1, (n_galaxies, n_samples))
    coverage = ev.coverage_curve(samples, truths, ev.NOMINAL_LEVELS)
    assert np.allclose(coverage, ev.NOMINAL_LEVELS, atol=0.03)


def test_truth_ranks_are_uniform_for_calibrated_posteriors():
    rng = np.random.default_rng(41)
    samples = rng.normal(0, 1, (5000, 200))
    truths = rng.normal(0, 1, 5000)
    ranks = ev.truth_ranks(samples, truths)
    histogram, _ = np.histogram(ranks, bins=10, range=(0, 1))
    assert np.all(np.abs(histogram / 5000 - 0.1) < 0.02)


def test_bias_and_nmad():
    truths = np.linspace(0, 1, 101)
    bias, nmad = ev.bias_and_nmad(truths + 0.2, truths)
    assert np.isclose(bias, 0.2)
    assert np.isclose(nmad, 0.0)


def test_width_ratio_is_one_when_posterior_equals_prior():
    rng = np.random.default_rng(42)
    prior = rng.normal(0, 1, 20000)
    posterior = rng.normal(0, 1, (300, 2000))
    assert np.isclose(ev.width_ratio(posterior, prior), 1.0, atol=0.05)
    assert ev.width_ratio(0.1 * posterior, prior) < 0.15


def test_compare_with_spectral():
    rng = np.random.default_rng(43)
    ml_samples = rng.normal(0, 1, (2, 5000))
    result = ev.compare_with_spectral(ml_samples, np.array([-0.1, 2.9]), np.array([0.0, 3.0]), np.array([0.1, 3.1]))
    assert bool(result["inside_68"][0]) and bool(result["inside_95"][0])
    assert not bool(result["inside_95"][1])
    assert abs(result["offset_in_widths"][0]) < 0.05
    assert result["offset_in_widths"][1] < -2.5


def test_windowed_prior_widths_follow_the_local_prior():
    rng = np.random.default_rng(44)
    redshift = np.sort(rng.uniform(0, 1, 60000))
    samples = rng.normal(0, 1, len(redshift)) * (1 + 2 * redshift)      # prior widens with redshift
    widths = ev.windowed_prior_widths(redshift, samples, np.array([0.1, 0.9, 5.0]), window=0.05)
    expected = 2 * 0.994 * np.array([1.2, 2.8])                          # 68% central interval of a normal
    assert np.allclose(widths[:2], expected, rtol=0.1)
    assert np.isnan(widths[2])                                           # no prior draws near z = 5


def test_width_ratio_uses_per_galaxy_prior_widths():
    rng = np.random.default_rng(45)
    posterior = rng.normal(0, 1, (300, 2000))
    wide = np.full(300, 2 * 0.994 * 2)                                   # prior twice as wide as the posterior
    assert np.isclose(ev.width_ratio(posterior, None, prior_widths=wide), 0.5, atol=0.05)


def test_zero_prior_width_is_left_out_not_infinite():
    rng = np.random.default_rng(46)
    posterior = rng.normal(0, 1e-3, (300, 2000))                         # a bin absent at high redshift
    widths = np.full(300, 2 * 0.994 * 2e-3)
    widths[:200] = 0.0                                                   # most priors collapsed: inf would win
    ratio = ev.width_ratio(posterior, None, prior_widths=widths)
    assert np.isfinite(ratio)
    assert np.isclose(ratio, 0.5, atol=0.05)
