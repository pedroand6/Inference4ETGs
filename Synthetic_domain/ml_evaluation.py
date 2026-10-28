"""
Metrics for the SFH posterior estimator: calibration (coverage, ranks), accuracy (bias, NMAD),
information per label (posterior/prior width) and agreement with the Bagpipes reference posteriors
(spectral or photometric fits). Also the comparison of the two SFH models (delayed-tau against the double
power law) on their shared labels, and their Bayesian evidence.
"""

import numpy as np
import pandas as pd

# =============================================================================
# 1. Constants
# =============================================================================
NOMINAL_LEVELS = np.linspace(0.1, 0.9, 9)   # central-interval probabilities of the coverage curve
NMAD_SCALE = 1.4826                          # makes the median absolute deviation a Gaussian sigma
ONE_SIGMA_LEVEL = 0.68                       # central probability of a "1 sigma" interval
TWO_SIGMA_LEVEL = 0.95                       # central probability of a "2 sigma" interval
PERCENT = 100
WINDOW_STEPS = 5                             # cell width = window / 5, so a galaxy's prior window is centred
                                             # within window / 10 of its redshift
# The "works" test of the delayed-tau spec (2026-10-26): a label agrees with a reference set when at least
# MIN_INSIDE_95 of the reference medians lie inside the ML 95% interval and the median offset is at most
# MAX_ABS_MEDIAN_OFFSET combined 1-sigma widths
MIN_INSIDE_95 = 0.80                         # fraction of galaxies
MAX_ABS_MEDIAN_OFFSET = 1.0                  # [combined 1-sigma widths]
ALL_LABELS_ROW = "all"                       # label of the verdict row that holds for every label at once
SUMMARY_SUFFIXES = ("p16", "median", "p84")  # column suffixes of a posterior stored as percentiles
GAUSSIAN_95_WIDTHS = 1.959964                # [1-sigma widths] half-width of the central 95% of a Gaussian
VERDICT_COLUMNS = ["n_galaxies", "inside_95", "median_offset", "agrees"]   # agreement_verdict per label

# Evidence scale on |ln B| (Jeffreys, as in Trotta 2008): below 1 inconclusive, 1 to 3 positive, above 3 strong
INCONCLUSIVE_LOG_BAYES = 1.0                 # [nats] |ln B| below this is inconclusive
STRONG_LOG_BAYES = 3.0                       # [nats] |ln B| above this is strong
DELAYED_MODEL = "delayed"                    # favoured by ln B >= INCONCLUSIVE_LOG_BAYES
DBLPLAW_MODEL = "dblplaw"                    # favoured by ln B <= -INCONCLUSIVE_LOG_BAYES
NO_MODEL = "neither"                         # an inconclusive ln B favours no model
TRUTH_SIGNS = {DELAYED_MODEL: 1.0, DBLPLAW_MODEL: -1.0}   # sign of ln B that favours each true model
ALL_MOCKS_ROW = "all"                        # true_model_recovery row over every mock


# =============================================================================
# 2. Intervals and calibration
# =============================================================================
def central_interval(samples, level):
    """Lower and upper bounds of the central interval holding `level` of the samples (axis -1).

    Why nanpercentile: a few galaxies may have NaN draws (failed decodes); they should not
    poison the whole interval.
    """
    tail = (1 - level) / 2 * PERCENT
    low = np.nanpercentile(samples, tail, axis=-1)
    high = np.nanpercentile(samples, PERCENT - tail, axis=-1)
    return low, high


def coverage_curve(samples, truths, levels=NOMINAL_LEVELS):
    """Empirical fraction of truths inside the central interval, for each nominal level.

    Inputs: samples (n, m) posterior draws of one label; truths (n,); levels: nominal probabilities.
    Returns: array of empirical coverages, one per level (equal to levels for a calibrated estimator).
    """
    coverages = []
    for level in levels:
        low, high = central_interval(samples, level)
        inside = (truths >= low) & (truths <= high)
        coverages.append(np.mean(inside))
    return np.array(coverages)


def truth_ranks(samples, truths):
    """Fraction of each galaxy's samples below its truth (uniform on [0, 1] when calibrated)."""
    return np.mean(samples < truths[:, np.newaxis], axis=1)


# =============================================================================
# 3. Accuracy and information
# =============================================================================
def bias_and_nmad(estimates, truths):
    """Median offset and robust scatter (NMAD) of estimates around truths. Returns: (bias, nmad).

    Why NMAD rather than the standard deviation: one catastrophic galaxy should not dominate the
    scatter; NMAD measures the typical spread, scaled to match a Gaussian sigma.
    """
    residuals = np.asarray(estimates) - np.asarray(truths)
    bias = np.nanmedian(residuals)
    nmad = NMAD_SCALE * np.nanmedian(np.abs(residuals - bias))
    return bias, nmad


def windowed_prior_widths(prior_redshift, prior_samples, galaxy_redshift, window, level=ONE_SIGMA_LEVEL):
    """Prior central-interval width of one label around each galaxy's redshift.

    Inputs: prior_redshift (k,) sorted ascending; prior_samples (k,) draws of the label, in the same order;
            galaxy_redshift (n,); window: half-width in redshift (dimensionless); level: central probability.
    Returns: widths (n,) in the label's units; NaN where no prior draw lies within the window.

    Why a grid: galaxies are grouped in cells of window / WINDOW_STEPS and every cell uses the draws within
    `window` of its centre, so the cost scales with the number of occupied cells, not galaxies.
    """
    cell_width = window / WINDOW_STEPS
    cells = np.round(np.asarray(galaxy_redshift) / cell_width).astype(int)
    widths = np.full(len(cells), np.nan)
    for cell in np.unique(cells):
        first = np.searchsorted(prior_redshift, cell * cell_width - window)
        last = np.searchsorted(prior_redshift, cell * cell_width + window)
        if last > first:
            low, high = central_interval(prior_samples[first:last], level)
            widths[cells == cell] = high - low
    return widths


def width_ratio(posterior_samples, prior_samples, level=ONE_SIGMA_LEVEL, prior_widths=None):
    """How much the data shrink a label: median over galaxies of posterior width / prior width.

    Inputs: posterior_samples (n, m); prior_samples (k,) draws of the same label from the prior;
            level: central-interval probability; prior_widths: optional (n,) prior width of each galaxy
            (e.g. windowed_prior_widths), used instead of the width of all prior_samples.
    Returns: ratio near 1 = prior-dominated, near 0 = well constrained.

    A zero prior width (a label the prior pins to one value in the galaxy's window, e.g. f_young of SFHs
    with no recent star formation) carries no information and would give inf, so it becomes NaN and is
    left out of the median.
    """
    if prior_widths is None:
        prior_low, prior_high = central_interval(prior_samples, level)
        prior_widths = prior_high - prior_low
    prior_widths = np.where(prior_widths > 0, prior_widths, np.nan)
    posterior_low, posterior_high = central_interval(posterior_samples, level)
    return float(np.nanmedian((posterior_high - posterior_low) / prior_widths))


# =============================================================================
# 4. Agreement with the spectral posteriors
# =============================================================================
def compare_with_spectral(ml_samples, spectral_p16, spectral_median, spectral_p84):
    """ML posterior vs the Bagpipes spectral posterior of each reference galaxy.

    Inputs: ml_samples (n, m) of one label; spectral 16th, 50th and 84th percentiles (n,).
    Returns: dict with offset_in_widths (median difference / combined 1-sigma width, ML minus spectral),
             inside_68 and inside_95 (spectral median inside the ML central interval).
    """
    ml_low, ml_high = central_interval(ml_samples, ONE_SIGMA_LEVEL)
    ml_median = np.nanmedian(ml_samples, axis=1)
    ml_sigma = (ml_high - ml_low) / 2
    spectral_sigma = (np.asarray(spectral_p84) - np.asarray(spectral_p16)) / 2
    # Why add the two sigmas in quadrature: the offset is judged against the total uncertainty of both methods.
    combined = np.sqrt(ml_sigma**2 + spectral_sigma**2)
    offset = (ml_median - np.asarray(spectral_median)) / combined

    low_95, high_95 = central_interval(ml_samples, TWO_SIGMA_LEVEL)
    inside_68 = (spectral_median >= ml_low) & (spectral_median <= ml_high)
    inside_95 = (spectral_median >= low_95) & (spectral_median <= high_95)
    return {"offset_in_widths": offset, "inside_68": inside_68, "inside_95": inside_95}


def reference_comparison_rows(ml_posteriors, galaxies, labels, prefix=""):
    """ML posterior against a Bagpipes reference posterior of each galaxy, one row per (galaxy, label).

    One helper for both reference sets, so the spectral and photometric verdicts use the same ML draws.
    Inputs: ml_posteriors: dict label -> ML samples (n, m), rows in the order of galaxies;
            galaxies: DataFrame (n rows) with specobjid and {prefix}<label>_p16/_median/_p84;
            labels: the labels to compare; prefix: "" (spectral fits) or "phot_" (photometric fits).
    Returns: list of dicts with label, specobjid, offset_in_widths, inside_68 and inside_95
             (as compare_with_spectral), labels outer and galaxies inner.
    """
    galaxies = galaxies.reset_index(drop=True)
    rows = []
    for label in labels:
        comparison = compare_with_spectral(ml_posteriors[label], galaxies[f"{prefix}{label}_p16"].to_numpy(),
                                           galaxies[f"{prefix}{label}_median"].to_numpy(),
                                           galaxies[f"{prefix}{label}_p84"].to_numpy())
        for galaxy_index in range(len(galaxies)):
            rows.append({"label": label, "specobjid": galaxies.loc[galaxy_index, "specobjid"],
                         "offset_in_widths": comparison["offset_in_widths"][galaxy_index],
                         "inside_68": comparison["inside_68"][galaxy_index],
                         "inside_95": comparison["inside_95"][galaxy_index]})
    return rows


def agreement_verdict(comparison_rows, labels, min_inside_95=MIN_INSIDE_95,
                      max_abs_median_offset=MAX_ABS_MEDIAN_OFFSET):
    """Whether the ML agrees with a reference label set, per label and for all labels at once.

    A label agrees when at least min_inside_95 of the reference medians lie inside the ML 95% interval
    and the median of offset_in_widths is within +-max_abs_median_offset (both limits inclusive).
    Inputs: comparison_rows: one row per (galaxy, label) with label, offset_in_widths (ML minus reference
            median, in combined 1-sigma widths) and inside_95 (bool), as built from compare_with_spectral
            (DataFrame or list of dicts); labels: the labels to judge, in output order;
            min_inside_95: fraction; max_abs_median_offset: combined 1-sigma widths.
    Returns: DataFrame with one row per label and a last row ALL_LABELS_ROW: label, n_galaxies (count),
             inside_95 (fraction), median_offset (combined widths) and agrees (bool). A label without
             galaxies has NaN fractions and does not agree; the last row agrees only if every label does.
    """
    rows = pd.DataFrame(comparison_rows)
    if not rows.empty and "label" not in rows.columns:
        raise KeyError("comparison_rows have no 'label' column")
    verdict_rows = []
    for label in labels:
        # Why the empty check: an empty comparison has no columns at all
        if rows.empty:
            label_rows = rows
        else:
            label_rows = rows[rows["label"] == label]
        n_galaxies = len(label_rows)
        if n_galaxies == 0:
            verdict_rows.append({"label": label, "n_galaxies": 0, "inside_95": np.nan,
                                 "median_offset": np.nan, "agrees": False})
            continue
        inside_95 = float(np.mean(label_rows["inside_95"].to_numpy(dtype=bool)))
        median_offset = float(np.nanmedian(label_rows["offset_in_widths"].to_numpy(dtype=float)))
        agrees = inside_95 >= min_inside_95 and abs(median_offset) <= max_abs_median_offset
        verdict_rows.append({"label": label, "n_galaxies": n_galaxies, "inside_95": inside_95,
                             "median_offset": median_offset, "agrees": agrees})

    every_label_agrees = all(row["agrees"] for row in verdict_rows)
    verdict_rows.append({"label": ALL_LABELS_ROW, "n_galaxies": np.nan, "inside_95": np.nan,
                         "median_offset": np.nan, "agrees": every_label_agrees})
    return pd.DataFrame(verdict_rows)


# =============================================================================
# 5. Model comparison: delayed-tau against the double power law
# =============================================================================
def percentile_summary(samples):
    """16th, 50th and 84th percentiles of each galaxy's samples (axis -1), as the reference CSV stores fits.

    Inputs: samples (n, m) of one label. Returns: (p16, median, p84), arrays (n,) in the label's units.
    """
    low, high = central_interval(samples, ONE_SIGMA_LEVEL)
    return low, np.nanmedian(samples, axis=-1), high


def compare_summaries_with_reference(p16, median, p84, reference_p16, reference_median, reference_p84):
    """compare_with_spectral for an estimate known only by its 16th, 50th and 84th percentiles.

    Why the Gaussian extension: a Bagpipes fit in the reference CSV keeps only these percentiles. Its 95%
    interval is therefore taken as median - 1.96 (median - p16) to median + 1.96 (p84 - median): each side is
    scaled from its own 1-sigma width.
    Inputs: the estimate's and the reference's percentiles, arrays (n,) in the label's units.
    Returns: dict with
             - offset_in_widths: estimate minus reference median, in combined 1-sigma widths;
             - inside_68 and inside_95: whether the reference median lies inside the estimate's interval.
    """
    p16 = np.asarray(p16, dtype=float)
    median = np.asarray(median, dtype=float)
    p84 = np.asarray(p84, dtype=float)
    reference_median = np.asarray(reference_median, dtype=float)
    sigma = (p84 - p16) / 2
    reference_sigma = (np.asarray(reference_p84, dtype=float) - np.asarray(reference_p16, dtype=float)) / 2
    combined = np.sqrt(sigma**2 + reference_sigma**2)
    offset = (median - reference_median) / combined

    low_95 = median - GAUSSIAN_95_WIDTHS * (median - p16)
    high_95 = median + GAUSSIAN_95_WIDTHS * (p84 - median)
    inside_68 = (reference_median >= p16) & (reference_median <= p84)
    inside_95 = (reference_median >= low_95) & (reference_median <= high_95)
    return {"offset_in_widths": offset, "inside_68": inside_68, "inside_95": inside_95}


def summary_columns(prefix, label):
    """The {prefix}<label>_p16 / _median / _p84 column names of one stored posterior."""
    columns = []
    for suffix in SUMMARY_SUFFIXES:
        columns.append(f"{prefix}{label}_{suffix}")
    return columns


def summary_comparison_rows(galaxies, labels, estimate_prefix, reference_prefix=""):
    """reference_comparison_rows for an estimate stored as percentiles (e.g. a photometric Bagpipes fit).

    Inputs: galaxies: DataFrame with specobjid, the estimate's {estimate_prefix}<label>_p16/_median/_p84
            columns and the reference's {reference_prefix}<label>_* columns (by default the spectral fit's);
            labels: the labels to compare.
    Returns: list of dicts with label, specobjid, offset_in_widths, inside_68 and inside_95, labels outer and
             galaxies inner. A galaxy missing any of the six percentiles of a label (the estimate's or the
             reference's p16, median or p84) is left out of that label. A label whose columns are missing
             gives no rows.
    """
    rows = []
    for label in labels:
        estimate_columns = summary_columns(estimate_prefix, label)
        reference_columns = summary_columns(reference_prefix, label)
        if not set(estimate_columns + reference_columns) <= set(galaxies.columns):
            continue
        complete = galaxies[estimate_columns + reference_columns].notna().all(axis=1)
        compared = galaxies[complete].reset_index(drop=True)
        comparison = compare_summaries_with_reference(
            compared[estimate_columns[0]].to_numpy(), compared[estimate_columns[1]].to_numpy(),
            compared[estimate_columns[2]].to_numpy(), compared[reference_columns[0]].to_numpy(),
            compared[reference_columns[1]].to_numpy(), compared[reference_columns[2]].to_numpy())
        for galaxy_index in range(len(compared)):
            rows.append({"label": label, "specobjid": compared.loc[galaxy_index, "specobjid"],
                         "offset_in_widths": comparison["offset_in_widths"][galaxy_index],
                         "inside_68": comparison["inside_68"][galaxy_index],
                         "inside_95": comparison["inside_95"][galaxy_index]})
    return rows


def model_difference_statistics(first_p16, first_median, first_p84, second_p16, second_median, second_p84):
    """How two estimates of one label differ over galaxies, first minus second.

    Why NMAD for the scatter: a few catastrophic galaxies should not set it (as in bias_and_nmad).
    Inputs: the 16th, 50th and 84th percentiles of each estimate, arrays (n,) in the label's units, paired by
            galaxy. A galaxy missing any value is left out.
    Returns: dict with
             - n_galaxies (count);
             - median_difference: median of first minus second medians, in label units;
             - scatter: NMAD of that difference around its median, in label units;
             - overlap_68: share of galaxies whose two 16-84% intervals overlap (touching counts).
             The statistics are NaN without galaxies.
    """
    first_p16 = np.asarray(first_p16, dtype=float)
    first_median = np.asarray(first_median, dtype=float)
    first_p84 = np.asarray(first_p84, dtype=float)
    second_p16 = np.asarray(second_p16, dtype=float)
    second_median = np.asarray(second_median, dtype=float)
    second_p84 = np.asarray(second_p84, dtype=float)
    usable = (np.isfinite(first_p16) & np.isfinite(first_median) & np.isfinite(first_p84)
              & np.isfinite(second_p16) & np.isfinite(second_median) & np.isfinite(second_p84))
    n_galaxies = int(usable.sum())
    if n_galaxies == 0:
        return {"n_galaxies": 0, "median_difference": np.nan, "scatter": np.nan, "overlap_68": np.nan}

    median_difference, scatter = bias_and_nmad(first_median[usable], second_median[usable])
    overlap = (first_p16[usable] <= second_p84[usable]) & (second_p16[usable] <= first_p84[usable])
    return {"n_galaxies": n_galaxies, "median_difference": float(median_difference), "scatter": float(scatter),
            "overlap_68": float(np.mean(overlap))}


def model_comparison_rows(galaxies, labels, first_prefix, second_prefix):
    """model_difference_statistics of each label, for two estimates stored side by side in one table.

    Inputs: galaxies: DataFrame with {prefix}<label>_p16/_median/_p84 of both estimates, one row per galaxy;
            labels: the labels to compare; first_prefix, second_prefix: their column prefixes (e.g. "phot_" and
            "phot_dblplaw_"). The difference is first minus second.
    Returns: DataFrame with one row per label: label, n_galaxies, median_difference, scatter and overlap_68.
             A label whose columns are missing for either estimate gets n_galaxies 0 and NaN statistics.
    """
    rows = []
    for label in labels:
        first_columns = summary_columns(first_prefix, label)
        second_columns = summary_columns(second_prefix, label)
        if not set(first_columns + second_columns) <= set(galaxies.columns):
            rows.append({"label": label, "n_galaxies": 0, "median_difference": np.nan, "scatter": np.nan,
                         "overlap_68": np.nan})
            continue
        statistics = model_difference_statistics(
            galaxies[first_columns[0]], galaxies[first_columns[1]], galaxies[first_columns[2]],
            galaxies[second_columns[0]], galaxies[second_columns[1]], galaxies[second_columns[2]])
        rows.append({"label": label, **statistics})
    return pd.DataFrame(rows, columns=["label", "n_galaxies", "median_difference", "scatter", "overlap_68"])


def side_by_side_verdicts(verdicts):
    """agreement_verdict tables of several estimators in one table: one row per label, one column group each.

    Inputs: verdicts: dict name -> agreement_verdict DataFrame, in display order.
    Returns: DataFrame indexed by label (first-seen order, ALL_LABELS_ROW last) with MultiIndex columns
             (name, VERDICT_COLUMNS). A label an estimator's table lacks is NaN there.
    """
    label_order = []
    has_all_row = False
    for verdict in verdicts.values():
        for label in verdict["label"]:
            if label == ALL_LABELS_ROW:
                has_all_row = True
            elif label not in label_order:
                label_order.append(label)
    if has_all_row:
        label_order.append(ALL_LABELS_ROW)

    tables = {}
    for name, verdict in verdicts.items():
        tables[name] = verdict.set_index("label")[VERDICT_COLUMNS].reindex(label_order)
    return pd.concat(tables, axis=1)


# =============================================================================
# 6. Bayesian evidence
# =============================================================================
def log_bayes_factor(log_evidence_delayed, log_evidence_dblplaw):
    """ln B = ln Z(delayed-tau) - ln Z(double power law) of each galaxy; positive favours delayed-tau.

    Inputs: the natural-log evidences of the two models' fits of the same data, arrays (n,) [nats].
    Returns: array (n,) [nats]; NaN where either evidence is missing.
    """
    return np.asarray(log_evidence_delayed, dtype=float) - np.asarray(log_evidence_dblplaw, dtype=float)


def nested_sampling_log_evidence_error(information, n_live_points):
    """Statistical 1-sigma uncertainty of a nested-sampling ln Z: sqrt(H / n_live) (Skilling 2006).

    Inputs: information: H, the Kullback-Leibler divergence from prior to posterior [nats];
            n_live_points: live points of the run (count).
    Returns: uncertainty of ln Z [nats].
    """
    return np.sqrt(np.asarray(information, dtype=float) / n_live_points)


def log_bayes_factor_error(log_evidence_error_delayed, log_evidence_error_dblplaw):
    """1-sigma uncertainty of ln B [nats]: the two fits are independent runs, so their errors add in quadrature."""
    error_delayed = np.asarray(log_evidence_error_delayed, dtype=float)
    error_dblplaw = np.asarray(log_evidence_error_dblplaw, dtype=float)
    return np.sqrt(error_delayed**2 + error_dblplaw**2)


def evidence_strength(log_bayes):
    """Strength of each ln B on the |ln B| scale.

    The categories are "inconclusive" (below INCONCLUSIVE_LOG_BAYES), "positive" (from it to
    STRONG_LOG_BAYES, both inclusive) and "strong" (above STRONG_LOG_BAYES).
    Inputs: log_bayes: array (n,) of ln B [nats].
    Returns: object array (n,) of strengths; None where ln B is NaN.
    """
    log_bayes = np.asarray(log_bayes, dtype=float)
    size = np.abs(log_bayes)
    strengths = np.full(len(log_bayes), None, dtype=object)
    strengths[size < INCONCLUSIVE_LOG_BAYES] = "inconclusive"
    strengths[(size >= INCONCLUSIVE_LOG_BAYES) & (size <= STRONG_LOG_BAYES)] = "positive"
    strengths[size > STRONG_LOG_BAYES] = "strong"
    return strengths


def favoured_model(log_bayes):
    """The model each ln B favours.

    Inputs: log_bayes: array (n,) of ln B [nats].
    Returns: object array (n,):
             - DELAYED_MODEL for ln B >= INCONCLUSIVE_LOG_BAYES;
             - DBLPLAW_MODEL for ln B <= -INCONCLUSIVE_LOG_BAYES;
             - NO_MODEL in between;
             - None where ln B is NaN.
    """
    log_bayes = np.asarray(log_bayes, dtype=float)
    models = np.full(len(log_bayes), None, dtype=object)
    models[np.abs(log_bayes) < INCONCLUSIVE_LOG_BAYES] = NO_MODEL
    models[log_bayes >= INCONCLUSIVE_LOG_BAYES] = DELAYED_MODEL
    models[log_bayes <= -INCONCLUSIVE_LOG_BAYES] = DBLPLAW_MODEL
    return models


def evidence_category_fractions(log_bayes):
    """Share of galaxies in each evidence category: inconclusive, then positive and strong for each model.

    Inputs: log_bayes: array (n,) of ln B (positive favours delayed-tau). NaN entries are left out.
    Returns: DataFrame with 5 rows (favoured_model, strength, n_galaxies, fraction), in this order:
             (NO_MODEL, inconclusive), (DELAYED_MODEL, positive), (DELAYED_MODEL, strong),
             (DBLPLAW_MODEL, positive), (DBLPLAW_MODEL, strong). Fractions are of the finite ln B and sum
             to 1; they are NaN when no ln B is finite.
    """
    log_bayes = np.asarray(log_bayes, dtype=float)
    finite = log_bayes[np.isfinite(log_bayes)]
    strengths = evidence_strength(finite)
    models = favoured_model(finite)
    categories = [(NO_MODEL, "inconclusive"), (DELAYED_MODEL, "positive"), (DELAYED_MODEL, "strong"),
                  (DBLPLAW_MODEL, "positive"), (DBLPLAW_MODEL, "strong")]
    rows = []
    for model, strength in categories:
        n_galaxies = int(np.sum((models == model) & (strengths == strength)))
        fraction = np.nan
        if len(finite) > 0:
            fraction = n_galaxies / len(finite)
        rows.append({"favoured_model": model, "strength": strength, "n_galaxies": n_galaxies,
                     "fraction": fraction})
    return pd.DataFrame(rows)


def true_model_recovery(cross_fits):
    """How often the evidence prefers the true model of mock galaxies fitted with both models.

    Inputs: cross_fits: DataFrame with one row per mock and the columns
            - truth_model (DELAYED_MODEL or DBLPLAW_MODEL);
            - log_evidence_delayed and log_evidence_dblplaw: ln Z of its fit with each model [nats]
              (ml_model.ipynb renames phot_log_evidence and phot_dblplaw_log_evidence of
              photometric_cross_fits.csv to these).
            Rows with a missing evidence are left out.
    Returns: DataFrame with one row per true model (DELAYED_MODEL, DBLPLAW_MODEL), then ALL_MOCKS_ROW.
             Its columns are truth_model and n_mocks, then four fractions of the mocks:
             - prefers_truth: ln B has the true model's sign;
             - strong_for_truth: |ln B| > STRONG_LOG_BAYES in favour of the true model;
             - inconclusive: |ln B| < INCONCLUSIVE_LOG_BAYES;
             - strong_for_other: |ln B| > STRONG_LOG_BAYES in favour of the wrong model.
             The fractions are NaN for a group without mocks.
    """
    log_bayes = log_bayes_factor(cross_fits["log_evidence_delayed"], cross_fits["log_evidence_dblplaw"])
    truth = cross_fits["truth_model"].to_numpy(dtype=object)
    usable = np.isfinite(log_bayes)

    # ln B with the sign that favours the true model: positive = the evidence prefers the truth
    signed = np.full(len(log_bayes), np.nan)
    for model, sign in TRUTH_SIGNS.items():
        is_model = truth == model
        signed[is_model] = sign * log_bayes[is_model]

    groups = {DELAYED_MODEL: usable & (truth == DELAYED_MODEL),
              DBLPLAW_MODEL: usable & (truth == DBLPLAW_MODEL),
              ALL_MOCKS_ROW: usable & np.isin(truth, list(TRUTH_SIGNS))}
    rows = []
    for name, in_group in groups.items():
        n_mocks = int(in_group.sum())
        if n_mocks == 0:
            rows.append({"truth_model": name, "n_mocks": 0, "prefers_truth": np.nan, "strong_for_truth": np.nan,
                         "inconclusive": np.nan, "strong_for_other": np.nan})
            continue
        group_signed = signed[in_group]
        rows.append({"truth_model": name, "n_mocks": n_mocks,
                     "prefers_truth": float(np.mean(group_signed > 0)),
                     "strong_for_truth": float(np.mean(group_signed > STRONG_LOG_BAYES)),
                     "inconclusive": float(np.mean(np.abs(group_signed) < INCONCLUSIVE_LOG_BAYES)),
                     "strong_for_other": float(np.mean(group_signed < -STRONG_LOG_BAYES))})
    return pd.DataFrame(rows)
