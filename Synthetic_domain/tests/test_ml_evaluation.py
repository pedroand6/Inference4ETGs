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


def comparison_rows(labels, reference_median, n_galaxies=60, seed=47):
    """Per-galaxy rows of ML vs reference comparisons: ML posteriors N(0, 1), reference median given, sigma 1.

    Inputs: labels: label names; reference_median: the reference median of every galaxy (label units);
            n_galaxies: count; seed: rng seed.
    Returns: list of dicts with label, offset_in_widths and inside_95, as ml_model.ipynb builds them.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for label in labels:
        ml_samples = rng.normal(0, 1, (n_galaxies, 4000))
        medians = np.full(n_galaxies, reference_median)
        comparison = ev.compare_with_spectral(ml_samples, medians - 1, medians, medians + 1)
        for galaxy_index in range(n_galaxies):
            rows.append({"label": label, "offset_in_widths": comparison["offset_in_widths"][galaxy_index],
                         "inside_95": comparison["inside_95"][galaxy_index]})
    return rows


def test_agreement_verdict_of_an_agreeing_set_is_true():
    verdict = ev.agreement_verdict(comparison_rows(["mw_age", "Z_fin"], 0.0), ["mw_age", "Z_fin"])
    by_label = verdict.set_index("label")
    assert list(verdict["label"]) == ["mw_age", "Z_fin", ev.ALL_LABELS_ROW]
    assert by_label.loc["mw_age", "inside_95"] == 1.0
    assert abs(by_label.loc["mw_age", "median_offset"]) < 0.05
    assert by_label.loc["mw_age", "n_galaxies"] == 60
    assert bool(by_label.loc["Z_fin", "agrees"])
    assert bool(by_label.loc[ev.ALL_LABELS_ROW, "agrees"])


def test_agreement_verdict_of_a_set_shifted_by_three_sigma_is_false():
    """Reference medians 3 sigma away: outside the ML 95% interval, offset -3 / sqrt(2) combined widths."""
    verdict = ev.agreement_verdict(comparison_rows(["mw_age"], 3.0), ["mw_age"])
    by_label = verdict.set_index("label")
    assert by_label.loc["mw_age", "inside_95"] < 0.80
    assert by_label.loc["mw_age", "median_offset"] < -1.0
    assert not bool(by_label.loc["mw_age", "agrees"])
    assert not bool(by_label.loc[ev.ALL_LABELS_ROW, "agrees"])


def test_agreement_verdict_thresholds_are_inclusive_and_settable():
    """A label passes when exactly at the thresholds; stricter thresholds fail it."""
    rows = []
    for offset, inside in ((0.0, True), (1.0, True), (1.0, True), (2.0, False), (-1.0, True)):
        rows.append({"label": "f_old", "offset_in_widths": offset, "inside_95": inside})
    verdict = ev.agreement_verdict(rows, ["f_old"]).set_index("label")
    assert verdict.loc["f_old", "inside_95"] == 0.80 and verdict.loc["f_old", "median_offset"] == 1.0
    assert bool(verdict.loc["f_old", "agrees"])
    strict = ev.agreement_verdict(rows, ["f_old"], min_inside_95=0.9).set_index("label")
    assert not bool(strict.loc["f_old", "agrees"])


def test_agreement_verdict_of_a_label_without_galaxies_is_nan_and_false():
    verdict = ev.agreement_verdict(comparison_rows(["mw_age"], 0.0), ["mw_age", "log_mstar"])
    by_label = verdict.set_index("label")
    assert np.isnan(by_label.loc["log_mstar", "inside_95"])
    assert np.isnan(by_label.loc["log_mstar", "median_offset"])
    assert by_label.loc["log_mstar", "n_galaxies"] == 0
    assert not bool(by_label.loc["log_mstar", "agrees"])
    assert not bool(by_label.loc[ev.ALL_LABELS_ROW, "agrees"])
    empty = ev.agreement_verdict([], ["mw_age"]).set_index("label")    # no comparison at all
    assert not bool(empty.loc[ev.ALL_LABELS_ROW, "agrees"])
