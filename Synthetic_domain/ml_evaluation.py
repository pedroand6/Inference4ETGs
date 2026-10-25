"""
Metrics for the SFH posterior estimator: calibration (coverage, ranks), accuracy (bias, NMAD),
information per bin (posterior/prior width) and agreement with the Bagpipes spectral posteriors.
"""

import numpy as np

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

    A zero prior width (a bin that cannot exist at the galaxy's redshift, e.g. frac_bin_14 above z ~ 0.48)
    carries no information and would give inf, so it becomes NaN and is left out of the median.
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
