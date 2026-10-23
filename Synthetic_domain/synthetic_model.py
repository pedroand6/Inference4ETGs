"""
Shared Bagpipes model of the S-PLUS synthetic domain.

training_set.ipynb and spectroscopic_reference.ipynb both import this module,
so the training set and the spectral fits use exactly the same prior.
"""

import hashlib
import os
import time
import warnings

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simpson
from scipy.stats import chi2 as chi_squared
from scipy.stats import norm, skewnorm

import bagpipes as pipes

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1. Prior (fixed before looking at the data)
# =============================================================================
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3)
UNIVERSE_AGE_MARGIN = 0.1        # [Gyr] oldest bin ends this long after the Big Bang

AGE_BIN_CENTRES = np.array([0.007, 0.021, 0.0425, 0.076, 0.1145, 0.1615, 0.23, 0.345,
                            0.545, 0.865, 1.37, 2.17, 3.44, 5.455, 8.645, 14.2])   # [Gyr ago]
N_BINS = len(AGE_BIN_CENTRES)
LABEL_EDGES = [0.1, 1.0]         # [Gyr] f_young < 0.1 < f_mid < 1.0 < f_old

Z_RANGE = (0.002, 0.7)                               # redshift
LOG_ZFIN_RANGE = (-1.0, 0.4)                         # log(Z_fin / Z_sun): 0.1 to 2.5
ZINI_FRACTION = 0.4                                  # Z_ini = ZINI_FRACTION * Z_fin
AV_ISM_RANGE = (0.01, 2.0)                           # [mag]
LOG_TAU_BC_RANGE = (np.log10(0.5), np.log10(4.0))    # birth-cloud optical depth
LOGU_RANGE = (-4.0, -1.0)                            # nebular ionisation parameter
DUST_SLOPE = 0.7                                     # CF00 attenuation-law exponent
BIRTH_CLOUD_AGE = 0.01                               # [Gyr]
LOG_MASS_FORMED = 10.0                               # [log M_sun]
TAU_TO_AV = 1.086                                    # A_V = 1.086 * tau_V

P_SKEW = 0.5                     # share of skew-normal SFHs (with or without burst)
P_BURST = 6 / 165                # share of skew-normal + burst SFHs
MIN_FRACTION = 1e-6              # keeps every bin as a Bagpipes component
SFR_FLOOR = 1e-14                # [yr^-1] sSFR floor for SFHs with no recent star formation

SKEW_PEAK_MEAN, SKEW_PEAK_SIGMA = 12.0, 3.0          # [Gyr ago]
SKEW_WIDTH_MEAN, SKEW_WIDTH_SIGMA = 1.0, 0.4         # [Gyr]
SKEW_SHAPE_RANGE = (-0.25, 0.75)                     # skew-normal shape parameter
LOG_BURST_AGE_RANGE = (np.log10(0.01), np.log10(0.5))   # [log Gyr ago]
BURST_WIDTH_RANGE = (0.040, 0.050)                   # [Gyr] full width
BURST_FRACTION_RANGE = (0.01, 0.10)                  # share of the mass in the burst

DIRICHLET_DIMS = []
for bin_index in range(N_BINS):
    DIRICHLET_DIMS.append(f"dir{bin_index}")
DIMS = ["z", "sfh_type", *DIRICHLET_DIMS, "peak", "width", "skew",
        "burst_age", "burst_width", "burst_frac", "log_zfin", "av_ism", "log_tau_bc", "logU"]


def uniform(draw, low, high):
    """Map a unit-interval draw onto [low, high)."""
    return low + (high - low) * draw


def positive_normal(draw, mean, sigma):
    """Map a unit-interval draw onto N(mean, sigma) truncated to positive values (inverse CDF)."""
    probability_below_zero = norm.cdf(0, mean, sigma)
    return norm.ppf(probability_below_zero + draw * (1 - probability_below_zero), mean, sigma)


def bin_edges(redshift):
    """Age-bin edges [Gyr ago], youngest first.

    Bins centred beyond the age of the Universe (minus the margin) are dropped,
    so high-redshift galaxies have fewer bins.
    """
    maximum_age = COSMOLOGY.age(redshift).value - UNIVERSE_AGE_MARGIN
    centres = AGE_BIN_CENTRES[AGE_BIN_CENTRES <= maximum_age]
    midpoints = (centres[:-1] + centres[1:]) / 2
    return np.concatenate([[0.0], midpoints, [maximum_age]])


# -----------------------------------------------------------------------------
# SFH branches: same signature, each returns mass fractions per bin (sum 1)
# -----------------------------------------------------------------------------
def dirichlet_fractions(params, edges):
    """Dirichlet(1) fractions: normalised exponential draws, one per bin."""
    n_bins = len(edges) - 1
    exponential_draws = np.zeros(n_bins)
    for bin_index in range(n_bins):
        exponential_draws[bin_index] = -np.log1p(-params[f"dir{bin_index}"])
    return exponential_draws / exponential_draws.sum()


def skew_normal_fractions(params, edges):
    """Fractions from one skew-normal episode in lookback time (typically old)."""
    peak = norm.ppf(params["peak"], SKEW_PEAK_MEAN, SKEW_PEAK_SIGMA)
    width = positive_normal(params["width"], SKEW_WIDTH_MEAN, SKEW_WIDTH_SIGMA)
    shape = uniform(params["skew"], *SKEW_SHAPE_RANGE)

    cumulative_mass = skewnorm.cdf(edges, shape, peak, width)
    fractions = np.diff(cumulative_mass)
    if fractions.sum() > 0:
        return fractions / fractions.sum()

    # The whole episode falls outside the bins: fall back to a flat SFH
    n_bins = len(edges) - 1
    return np.full(n_bins, 1 / n_bins)


def skew_normal_burst_fractions(params, edges):
    """Skew-normal fractions plus a short recent burst holding 1-10% of the mass."""
    base_fractions = skew_normal_fractions(params, edges)
    burst_age = 10 ** uniform(params["burst_age"], *LOG_BURST_AGE_RANGE)
    burst_half_width = uniform(params["burst_width"], *BURST_WIDTH_RANGE) / 2
    burst_fraction = uniform(params["burst_frac"], *BURST_FRACTION_RANGE)

    # Time the burst overlaps each bin, so its mass is shared in proportion
    bin_starts = edges[:-1]
    bin_ends = edges[1:]
    overlap_start = np.maximum(bin_starts, burst_age - burst_half_width)
    overlap_end = np.minimum(bin_ends, burst_age + burst_half_width)
    overlap = np.clip(overlap_end - overlap_start, 0, None)

    fractions = (1 - burst_fraction) * base_fractions + burst_fraction * overlap / overlap.sum()
    return fractions / fractions.sum()


SFH_BRANCHES = {
    0: dirichlet_fractions,
    1: skew_normal_fractions,
    2: skew_normal_burst_fractions,
}


def choose_sfh_branch(sfh_type_draw):
    """SFH branch for a unit draw: 0 Dirichlet, 1 skew-normal, 2 skew-normal + burst."""
    if sfh_type_draw < 1 - P_SKEW:
        return 0
    if sfh_type_draw < 1 - P_BURST:
        return 1
    return 2


def bin_metallicities(fractions, z_final):
    """Z/Z_sun of each bin: rises from Z_ini to Z_final with the mass already formed."""
    z_initial = ZINI_FRACTION * z_final
    # Bins run youngest first, so the mass formed before a bin is the sum of the older ones
    mass_before = np.cumsum(fractions[::-1])[::-1] - fractions
    enrichment = (mass_before + fractions / 2) / fractions.sum()
    return z_initial + (z_final - z_initial) * enrichment


def components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one unit-cube vector u (order of DIMS).

    fixed_redshift: if given, replaces the redshift drawn from u (spectral fits).
    Returns (bagpipes_components, drawn) where drawn holds z, sfh_type, Z_fin, Av_ism, Av_bc, logU.
    """
    params = dict(zip(DIMS, u))

    # 1. Redshift and the age bins that exist at that redshift
    if fixed_redshift is None:
        redshift = uniform(params["z"], *Z_RANGE)
    else:
        redshift = fixed_redshift
    edges = bin_edges(redshift)

    # 2. Star formation history
    sfh_type = choose_sfh_branch(params["sfh_type"])
    fractions = SFH_BRANCHES[sfh_type](params, edges)
    fractions = np.fmax(fractions, MIN_FRACTION)

    # 3. Chemical enrichment
    z_final = 10 ** uniform(params["log_zfin"], *LOG_ZFIN_RANGE)
    metallicities = bin_metallicities(fractions, z_final)

    # 4. Dust and nebular emission
    av_ism = uniform(params["av_ism"], *AV_ISM_RANGE)
    av_birth_cloud = TAU_TO_AV * 10 ** uniform(params["log_tau_bc"], *LOG_TAU_BC_RANGE)
    log_u = uniform(params["logU"], *LOGU_RANGE)

    bagpipes_components = {
        "redshift": redshift,
        "t_bc": BIRTH_CLOUD_AGE,
        "dust": {"type": "CF00", "Av": av_ism, "eta": 1 + av_birth_cloud / av_ism, "n": DUST_SLOPE},
        "nebular": {"logU": log_u},
    }
    # One "constant" component per bin; Bagpipes reads the type from the name minus its last letter
    for bin_index in range(len(fractions)):
        name = "constant" + chr(ord("A") + bin_index)
        bagpipes_components[name] = {
            "age_min": edges[bin_index],
            "age_max": edges[bin_index + 1],
            "massformed": LOG_MASS_FORMED + np.log10(fractions[bin_index]),
            "metallicity": metallicities[bin_index],
        }

    drawn = {"z": redshift, "sfh_type": sfh_type, "Z_fin": z_final,
             "Av_ism": av_ism, "Av_bc": av_birth_cloud, "logU": log_u}
    return bagpipes_components, drawn


def bins_from_components(bagpipes_components):
    """Mass fractions, metallicities and edges [Gyr ago] read back from the constant components."""
    names = []
    for name in bagpipes_components:
        if name.startswith("constant"):
            names.append(name)
    names.sort()

    n_bins = len(names)
    fractions = np.zeros(n_bins)
    metallicities = np.zeros(n_bins)
    edges = np.zeros(n_bins + 1)
    for bin_index, name in enumerate(names):
        component = bagpipes_components[name]
        fractions[bin_index] = 10 ** (component["massformed"] - LOG_MASS_FORMED)
        metallicities[bin_index] = component["metallicity"]
        edges[bin_index] = component["age_min"]
        edges[bin_index + 1] = component["age_max"]
    return fractions, metallicities, edges


# =============================================================================
# 2. Labels, S-PLUS magnitudes and model simulation
# =============================================================================
MODULE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
FILTER_LIST_FILE = os.path.join(MODULE_DIRECTORY, "splus-filters", "splus_filters.txt")
SPEED_OF_LIGHT = 2.99792458e18   # [Angstrom / s]
AB_ZERO_POINT = 48.60

# S-PLUS bands, in the order of splus_filters.txt
BANDS = ["u", "J0378", "J0395", "J0410", "J0430", "g", "J0515", "r", "J0660", "i", "J0861", "z"]

FILTER_FILES = []
for relative_path in np.loadtxt(FILTER_LIST_FILE, dtype="str"):
    FILTER_FILES.append(os.path.join(MODULE_DIRECTORY, relative_path))


def pivot_wavelength(filter_file):
    """Pivot wavelength [Angstrom] of a filter curve file (wavelength, transmission)."""
    wavelength, transmission = np.loadtxt(filter_file, unpack=True)
    numerator = simpson(transmission * wavelength, x=wavelength)
    denominator = simpson(transmission / wavelength, x=wavelength)
    return np.sqrt(numerator / denominator)


PIVOT_WAVELENGTHS = np.zeros(len(FILTER_FILES))
for band_index, filter_file in enumerate(FILTER_FILES):
    PIVOT_WAVELENGTHS[band_index] = pivot_wavelength(filter_file)

SUMMARY_LABEL_COLUMNS = ["log_mstar", "log_sfr", "mw_age", "f_young", "f_mid", "f_old"]

BIN_LABEL_COLUMNS = []
for label_kind in ("frac_bin", "logsfr_bin", "zmet_bin"):
    for bin_index in range(N_BINS):
        BIN_LABEL_COLUMNS.append(f"{label_kind}_{bin_index:02d}")

# Labels that scale with brightness: log value - 0.4 * (m_obs - m_model)
MASS_DEPENDENT_LABELS = ["log_mstar", "log_sfr"]
for bin_index in range(N_BINS):
    MASS_DEPENDENT_LABELS.append(f"logsfr_bin_{bin_index:02d}")


def bin_labels(fractions, metallicities, edges, log_mass_formed=LOG_MASS_FORMED):
    """16-bin labels, youngest bin first.

    fractions: mass formed per existing bin (any normalisation); metallicities: Z/Z_sun;
    edges: [Gyr ago]. Bins that do not exist at this redshift get frac 0, logsfr and zmet NaN.
    Returns a dict ordered as BIN_LABEL_COLUMNS.
    """
    n_existing = len(fractions)
    normalised_fractions = fractions / fractions.sum()
    widths_in_years = np.diff(edges) * 1e9

    labels = {}
    for bin_index in range(N_BINS):
        if bin_index < n_existing:
            labels[f"frac_bin_{bin_index:02d}"] = normalised_fractions[bin_index]
        else:
            labels[f"frac_bin_{bin_index:02d}"] = 0.0

    for bin_index in range(N_BINS):
        if bin_index < n_existing:
            mass_in_bin = 10**log_mass_formed * normalised_fractions[bin_index]
            labels[f"logsfr_bin_{bin_index:02d}"] = np.log10(mass_in_bin / widths_in_years[bin_index])
        else:
            labels[f"logsfr_bin_{bin_index:02d}"] = np.nan

    for bin_index in range(N_BINS):
        if bin_index < n_existing:
            labels[f"zmet_bin_{bin_index:02d}"] = metallicities[bin_index]
        else:
            labels[f"zmet_bin_{bin_index:02d}"] = np.nan
    return labels


def summary_labels(bagpipes_sfh):
    """log_mstar (living), log_sfr (last 100 Myr), mw_age [Gyr] and f_young/f_mid/f_old."""
    mass_per_age_step = bagpipes_sfh.sfh * bagpipes_sfh.age_widths
    total_mass = mass_per_age_step.sum()
    fraction_younger_than = []
    for edge in LABEL_EDGES:
        younger = bagpipes_sfh.ages < edge * 1e9
        fraction_younger_than.append(mass_per_age_step[younger].sum() / total_mass)

    sfr_floor = SFR_FLOOR * 10**bagpipes_sfh.stellar_mass
    log_sfr = np.log10(max(bagpipes_sfh.sfr, sfr_floor))
    return {
        "log_mstar": bagpipes_sfh.stellar_mass,
        "log_sfr": log_sfr,
        "mw_age": bagpipes_sfh.mass_weighted_age,
        "f_young": fraction_younger_than[0],
        "f_mid": fraction_younger_than[1] - fraction_younger_than[0],
        "f_old": 1 - fraction_younger_than[1],
    }


def mass_weighted_metallicity(fractions, metallicities):
    """Mean Z/Z_sun weighted by the mass formed in each bin."""
    return np.sum(fractions * metallicities) / np.sum(fractions)


def ab_magnitudes(flux_per_angstrom):
    """AB magnitudes from Bagpipes photometry [erg/s/cm2/A]; all NaN if any flux is not positive."""
    if not np.all(flux_per_angstrom > 0):
        return np.full(len(PIVOT_WAVELENGTHS), np.nan)
    flux_per_hertz = flux_per_angstrom * PIVOT_WAVELENGTHS**2 / SPEED_OF_LIGHT
    return -2.5 * np.log10(flux_per_hertz) - AB_ZERO_POINT


def simulate_models(unit_vectors):
    """Drawn parameters, labels and AB magnitudes (10^10 M_sun formed) for a block of unit-cube vectors.

    One Bagpipes model is kept per number of bins: update() is ~70x faster than building a new one.
    """
    models_by_bin_count = {}
    rows = []
    for u in unit_vectors:
        bagpipes_components, drawn = components(u)
        fractions, metallicities, edges = bins_from_components(bagpipes_components)

        n_bins = len(fractions)
        if n_bins in models_by_bin_count:
            model = models_by_bin_count[n_bins]
            model.update(bagpipes_components)
        else:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
            models_by_bin_count[n_bins] = model

        row = dict(drawn)
        row.update(summary_labels(model.sfh))
        row.update(bin_labels(fractions, metallicities, edges))
        magnitudes = ab_magnitudes(model.photometry)
        for band, magnitude in zip(BANDS, magnitudes):
            row[band] = magnitude
        rows.append(row)
    return pd.DataFrame(rows)
