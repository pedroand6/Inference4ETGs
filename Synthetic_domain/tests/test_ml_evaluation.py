"""Tests for ml_evaluation.py. Run from Synthetic_domain/."""

import numpy as np
import pandas as pd
import pytest

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
    posterior = rng.normal(0, 1e-3, (300, 2000))                         # a label the prior pins to one value
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
    labels = ["mw_age", "log_metallicity"]
    verdict = ev.agreement_verdict(comparison_rows(labels, 0.0), labels)
    by_label = verdict.set_index("label")
    assert list(verdict["label"]) == ["mw_age", "log_metallicity", ev.ALL_LABELS_ROW]
    assert by_label.loc["mw_age", "inside_95"] == 1.0
    assert abs(by_label.loc["mw_age", "median_offset"]) < 0.05
    assert by_label.loc["mw_age", "n_galaxies"] == 60
    assert bool(by_label.loc["log_metallicity", "agrees"])
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


def test_agreement_verdict_refuses_rows_without_a_label_column():
    rows = [{"offset_in_widths": 0.0, "inside_95": True}]
    with pytest.raises(KeyError, match="label"):
        ev.agreement_verdict(rows, ["mw_age"])


def reference_galaxies(prefix, medians):
    """Reference rows with {prefix}mw_age_p16/_median/_p84 (sigma 1) for the given medians."""
    medians = np.asarray(medians, dtype=float)
    return pd.DataFrame({"specobjid": [f"id{index}" for index in range(len(medians))],
                         f"{prefix}mw_age_p16": medians - 1, f"{prefix}mw_age_median": medians,
                         f"{prefix}mw_age_p84": medians + 1})


def test_reference_comparison_rows_reads_the_prefixed_columns():
    rng = np.random.default_rng(48)
    posteriors = {"mw_age": rng.normal(0, 1, (2, 4000))}
    galaxies = reference_galaxies("phot_", [0.0, 5.0])
    rows = ev.reference_comparison_rows(posteriors, galaxies, ["mw_age"], prefix="phot_")

    assert [row["specobjid"] for row in rows] == ["id0", "id1"]
    assert all(row["label"] == "mw_age" for row in rows)
    assert bool(rows[0]["inside_95"]) and not bool(rows[1]["inside_95"])
    assert abs(rows[0]["offset_in_widths"]) < 0.1
    assert set(rows[0]) == {"label", "specobjid", "offset_in_widths", "inside_68", "inside_95"}


def test_reference_comparison_rows_with_the_default_prefix_matches_compare_with_spectral():
    rng = np.random.default_rng(49)
    posteriors = {"mw_age": rng.normal(0, 1, (3, 1000))}
    galaxies = reference_galaxies("", [0.0, 1.0, -2.0])
    rows = ev.reference_comparison_rows(posteriors, galaxies, ["mw_age"])
    direct = ev.compare_with_spectral(posteriors["mw_age"], galaxies["mw_age_p16"], galaxies["mw_age_median"],
                                      galaxies["mw_age_p84"])
    assert np.allclose([row["offset_in_widths"] for row in rows], direct["offset_in_widths"])


# =============================================================================
# Model comparison and Bayesian evidence (delayed-tau against the double power law)
# =============================================================================
def test_percentile_summary_gives_the_16th_50th_and_84th_percentiles():
    samples = np.tile(np.arange(101.0), (2, 1))
    p16, median, p84 = ev.percentile_summary(samples)
    assert np.allclose(p16, 16.0)
    assert np.allclose(median, 50.0)
    assert np.allclose(p84, 84.0)


def test_compare_summaries_with_reference_extends_each_side_as_a_gaussian():
    """Estimate median 0 with p16 = -1 and p84 = +1 (+2 in the last row); reference sigma 1."""
    reference_median = np.array([0.0, 1.5, 2.5, 3.5])
    result = ev.compare_summaries_with_reference(np.full(4, -1.0), np.zeros(4), np.array([1.0, 1.0, 1.0, 2.0]),
                                                 reference_median - 1, reference_median, reference_median + 1)
    assert list(result["inside_68"]) == [True, False, False, False]
    # 95% upper bound: 1.96 for the first three rows, 3.92 for the wider upper side of the last
    assert list(result["inside_95"]) == [True, True, False, True]
    assert np.isclose(result["offset_in_widths"][0], 0.0)
    assert np.isclose(result["offset_in_widths"][1], -1.5 / np.sqrt(2))


def test_summary_comparison_rows_leave_out_galaxies_without_both_fits():
    galaxies = pd.DataFrame({
        "specobjid": ["a", "b", "c"],
        "phot_dblplaw_mw_age_p16": [-1.0, -1.0, np.nan],
        "phot_dblplaw_mw_age_median": [0.0, 0.0, np.nan],
        "phot_dblplaw_mw_age_p84": [1.0, 1.0, np.nan],
        "mw_age_p16": [-1.0, 4.0, -1.0],
        "mw_age_median": [0.0, 5.0, 0.0],
        "mw_age_p84": [1.0, 6.0, 1.0],
    })
    rows = ev.summary_comparison_rows(galaxies, ["mw_age", "log_ssfr"], "phot_dblplaw_")   # no log_ssfr columns
    assert [row["specobjid"] for row in rows] == ["a", "b"]
    assert bool(rows[0]["inside_95"]) and not bool(rows[1]["inside_95"])
    assert set(rows[0]) == {"label", "specobjid", "offset_in_widths", "inside_68", "inside_95"}


def test_model_difference_statistics_recover_a_known_offset_and_scatter():
    rng = np.random.default_rng(50)
    second_median = rng.normal(10, 1, 20000)
    first_median = second_median + 0.3 + rng.normal(0, 0.1, len(second_median))
    half_width = 0.05
    statistics = ev.model_difference_statistics(first_median - half_width, first_median, first_median + half_width,
                                                second_median - half_width, second_median,
                                                second_median + half_width)
    assert statistics["n_galaxies"] == 20000
    assert np.isclose(statistics["median_difference"], 0.3, atol=0.005)
    assert np.isclose(statistics["scatter"], 0.1, rtol=0.05)
    # The intervals overlap when the difference is at most 0.1, i.e. 2 sigma below its mean: Phi(-2) = 0.023
    assert np.isclose(statistics["overlap_68"], 0.0228, atol=0.005)


def test_model_difference_statistics_on_hand_made_intervals():
    """[0, 1] against [0.5, 2] overlaps, against [1, 2] touches (counts), against [1.1, 2] does not."""
    statistics = ev.model_difference_statistics(
        np.array([0.0, 0.0, 0.0, np.nan]), np.array([0.5, 0.5, 0.5, 0.5]), np.array([1.0, 1.0, 1.0, 1.0]),
        np.array([0.5, 1.0, 1.1, 0.0]), np.array([1.0, 1.5, 1.5, 0.5]), np.array([2.0, 2.0, 2.0, 2.0]))
    assert statistics["n_galaxies"] == 3                       # the galaxy with a NaN is left out
    assert np.isclose(statistics["overlap_68"], 2 / 3)
    assert np.isclose(statistics["median_difference"], -1.0)   # differences -0.5, -1.0, -1.0


def test_model_difference_statistics_without_galaxies_are_nan():
    missing = np.array([np.nan])
    statistics = ev.model_difference_statistics(missing, missing, missing, missing, missing, missing)
    assert statistics["n_galaxies"] == 0
    assert np.isnan(statistics["median_difference"])
    assert np.isnan(statistics["scatter"])
    assert np.isnan(statistics["overlap_68"])


def paired_galaxies(first_prefix, second_prefix, first_medians, second_medians):
    """Galaxies with {prefix}mw_age_p16/_median/_p84 (sigma 1) of two estimates."""
    columns = {}
    for prefix, medians in ((first_prefix, first_medians), (second_prefix, second_medians)):
        medians = np.asarray(medians, dtype=float)
        columns[f"{prefix}mw_age_p16"] = medians - 1
        columns[f"{prefix}mw_age_median"] = medians
        columns[f"{prefix}mw_age_p84"] = medians + 1
    return pd.DataFrame(columns)


def test_model_comparison_rows_read_both_prefixes_and_mark_missing_labels():
    galaxies = paired_galaxies("phot_", "phot_dblplaw_", [1.0, 2.0, 3.0], [0.5, 1.5, 2.5])
    table = ev.model_comparison_rows(galaxies, ["mw_age", "av_ism"], "phot_", "phot_dblplaw_")
    assert list(table.columns) == ["label", "n_galaxies", "median_difference", "scatter", "overlap_68"]
    assert list(table["label"]) == ["mw_age", "av_ism"]
    by_label = table.set_index("label")
    assert by_label.loc["mw_age", "n_galaxies"] == 3
    assert np.isclose(by_label.loc["mw_age", "median_difference"], 0.5)
    assert np.isclose(by_label.loc["mw_age", "scatter"], 0.0)
    assert by_label.loc["mw_age", "overlap_68"] == 1.0
    assert by_label.loc["av_ism", "n_galaxies"] == 0              # no av_ism columns
    assert np.isnan(by_label.loc["av_ism", "median_difference"])


def test_side_by_side_verdicts_put_each_estimator_in_its_own_column_group():
    rows = [{"label": "mw_age", "offset_in_widths": 0.0, "inside_95": True},
            {"label": "mw_age", "offset_in_widths": 0.5, "inside_95": True}]
    delayed_verdict = ev.agreement_verdict(rows, ["log_age", "mw_age"])
    dblplaw_verdict = ev.agreement_verdict(rows, ["log_alpha", "mw_age"])
    table = ev.side_by_side_verdicts({"delayed": delayed_verdict, "dblplaw": dblplaw_verdict})
    assert list(table.index) == ["log_age", "mw_age", "log_alpha", ev.ALL_LABELS_ROW]
    assert list(table.columns.get_level_values(0).unique()) == ["delayed", "dblplaw"]
    assert list(table["delayed"].columns) == ev.VERDICT_COLUMNS
    assert table.loc["mw_age", ("delayed", "inside_95")] == 1.0
    assert np.isclose(table.loc["mw_age", ("dblplaw", "median_offset")], 0.25)
    assert np.isnan(table.loc["log_alpha", ("delayed", "inside_95")])     # a label the delayed table lacks


def test_log_bayes_factor_is_delayed_minus_dblplaw_and_keeps_nan():
    log_bayes = ev.log_bayes_factor([-100.0, -50.0, np.nan], [-102.0, -45.0, -10.0])
    assert np.allclose(log_bayes[:2], [2.0, -5.0])
    assert np.isnan(log_bayes[2])


def test_evidence_strength_and_favoured_model_at_the_category_edges():
    log_bayes = np.array([0.0, 0.99, 1.0, -1.0, 2.5, 3.0, 3.01, -4.0, np.nan])
    assert list(ev.evidence_strength(log_bayes)) == ["inconclusive", "inconclusive", "positive", "positive",
                                                     "positive", "positive", "strong", "strong", None]
    assert list(ev.favoured_model(log_bayes)) == ["neither", "neither", "delayed", "dblplaw", "delayed",
                                                  "delayed", "delayed", "dblplaw", None]


def test_evidence_category_fractions_of_a_hand_made_table():
    log_bayes = np.array([0.5, -0.2, 0.0, 2.0, 5.0, 4.0, -1.5, -3.5, -10.0, np.nan])
    table = ev.evidence_category_fractions(log_bayes)
    assert list(zip(table["favoured_model"], table["strength"])) == [
        ("neither", "inconclusive"), ("delayed", "positive"), ("delayed", "strong"),
        ("dblplaw", "positive"), ("dblplaw", "strong")]
    assert list(table["n_galaxies"]) == [3, 1, 2, 1, 2]
    assert np.allclose(table["fraction"], np.array([3, 1, 2, 1, 2]) / 9)     # the NaN galaxy is left out
    assert np.isclose(table["fraction"].sum(), 1.0)


def test_evidence_category_fractions_without_finite_values_are_nan():
    table = ev.evidence_category_fractions(np.array([np.nan, np.nan]))
    assert (table["n_galaxies"] == 0).all()
    assert table["fraction"].isna().all()


def test_evidence_uncertainties():
    assert np.isclose(ev.nested_sampling_log_evidence_error(12.0, 300), 0.2)
    assert np.isclose(ev.log_bayes_factor_error(0.3, 0.4), 0.5)
    assert np.allclose(ev.log_bayes_factor_error(np.array([0.3, 0.0]), np.array([0.4, 0.1])), [0.5, 0.1])


def cross_fit_table(truth_models, log_bayes):
    """Cross-fitted mocks with the given true models and ln B (ln Z of the DPL fit fixed at -100)."""
    log_bayes = np.asarray(log_bayes, dtype=float)
    mock_ids = []
    for index in range(len(log_bayes)):
        mock_ids.append(f"calib_{index:03d}")
    return pd.DataFrame({"mock_id": mock_ids, "truth_model": truth_models,
                         "log_evidence_delayed": -100.0 + log_bayes,
                         "log_evidence_dblplaw": np.full(len(log_bayes), -100.0)})


def test_true_model_recovery_of_hand_made_mocks():
    table = cross_fit_table(["delayed"] * 4 + ["dblplaw"] * 4, [5.0, 2.0, -0.5, -4.0, -6.0, -0.2, 1.5, np.nan])
    recovery = ev.true_model_recovery(table).set_index("truth_model")
    assert list(recovery.index) == ["delayed", "dblplaw", ev.ALL_MOCKS_ROW]
    assert list(recovery["n_mocks"]) == [4, 3, 7]                            # the NaN mock is left out
    assert np.allclose(recovery["prefers_truth"], [2 / 4, 2 / 3, 4 / 7])
    assert np.allclose(recovery["strong_for_truth"], [1 / 4, 1 / 3, 2 / 7])
    assert np.allclose(recovery["inconclusive"], [1 / 4, 1 / 3, 2 / 7])
    assert np.allclose(recovery["strong_for_other"], [1 / 4, 0.0, 1 / 7])


def test_true_model_recovery_without_mocks_is_nan():
    recovery = ev.true_model_recovery(cross_fit_table([], []))
    assert (recovery["n_mocks"] == 0).all()
    assert recovery["prefers_truth"].isna().all()
