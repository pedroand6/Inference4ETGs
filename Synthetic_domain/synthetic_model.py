"""
Shared Bagpipes model of the S-PLUS synthetic domain.

training_set.ipynb and spectroscopic_reference.ipynb both import this module,
so the training set and the spectral fits use exactly the same prior.
"""

import hashlib
import os
import time
import warnings

import extinction
import h5py
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
N_BINS = len(AGE_BIN_CENTRES)    # 16 labelled bins, but only 15 are ever populated (see below)
# The last centre (14.2 Gyr) is older than the Universe for every z >= 0.002 (13.34 Gyr minus the margin), so
# *_bin_15 is always padding. The oldest populated bin runs from 7.05 Gyr to the age of the Universe, which
# caps the mass-weighted age (mw_age) at about 10.2 Gyr at z = 0.002. Kept as in the legacy notebook.
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
    """Map a unit-interval draw onto [low, high).

    Inputs: draw in [0, 1) (dimensionless); low, high: bounds in the units of the parameter.
    Returns: parameter value in [low, high), in the same units as low and high.
    """
    return low + (high - low) * draw


def positive_normal(draw, mean, sigma):
    """Map a unit-interval draw onto N(mean, sigma) truncated to positive values (inverse CDF).

    Inputs: draw in [0, 1) (dimensionless); mean, sigma: of the untruncated normal, in the parameter's units.
    Returns: positive value in the same units as mean and sigma.
    """
    probability_below_zero = norm.cdf(0, mean, sigma)
    return norm.ppf(probability_below_zero + draw * (1 - probability_below_zero), mean, sigma)


def bin_edges(redshift):
    """Age-bin edges [Gyr ago], youngest first.

    Bins centred beyond the age of the Universe (minus the margin) are dropped,
    so high-redshift galaxies have fewer bins.

    Inputs: redshift: dimensionless redshift z.
    Returns: edges [Gyr ago], youngest first: 0.0 (today) up to the maximum age (Universe age minus margin).
    """
    maximum_age = COSMOLOGY.age(redshift).value - UNIVERSE_AGE_MARGIN
    centres = AGE_BIN_CENTRES[AGE_BIN_CENTRES <= maximum_age]
    midpoints = (centres[:-1] + centres[1:]) / 2
    return np.concatenate([[0.0], midpoints, [maximum_age]])


# -----------------------------------------------------------------------------
# SFH branches: same signature, each returns mass fractions per bin (sum 1)
# -----------------------------------------------------------------------------
def dirichlet_fractions(params, edges):
    """Dirichlet(1) fractions: normalised exponential draws, one per bin.

    Inputs: params: dict of unit-cube draws keyed by DIMS names (dimensionless);
            edges: bin edges [Gyr ago] from bin_edges.
    Returns: mass fraction per bin (dimensionless, sums to 1), youngest bin first.
    """
    n_bins = len(edges) - 1
    exponential_draws = np.zeros(n_bins)
    for bin_index in range(n_bins):
        exponential_draws[bin_index] = -np.log1p(-params[f"dir{bin_index}"])
    return exponential_draws / exponential_draws.sum()


def skew_normal_fractions(params, edges):
    """Fractions from one skew-normal episode in lookback time (typically old).

    Inputs: params: dict of unit-cube draws (uses peak, width, skew; dimensionless);
            edges: bin edges [Gyr ago] from bin_edges.
    Returns: mass fraction per bin (dimensionless, sums to 1); flat if the episode lies outside all bins.
    """
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
    """Skew-normal fractions plus a short recent burst holding 1-10% of the mass.

    Inputs: params: dict of unit-cube draws (the skew-normal keys plus burst_age, burst_width,
            burst_frac; dimensionless); edges: bin edges [Gyr ago] from bin_edges.
    Returns: mass fraction per bin (dimensionless, sums to 1), youngest bin first.
    """
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
    """SFH branch for a unit draw: 0 Dirichlet, 1 skew-normal, 2 skew-normal + burst.

    Inputs: sfh_type_draw: unit-cube draw in [0, 1) (dimensionless).
    Returns: integer branch index 0, 1 or 2, a key of SFH_BRANCHES.
    """
    if sfh_type_draw < 1 - P_SKEW:
        return 0
    if sfh_type_draw < 1 - P_BURST:
        return 1
    return 2


def bin_metallicities(fractions, z_final):
    """Z/Z_sun of each bin: rises from Z_ini to Z_final with the mass already formed.

    Inputs: fractions: mass fraction per bin, youngest first (any normalisation);
            z_final: final metallicity Z_fin in Z/Z_sun.
    Returns: Z/Z_sun per bin (same length as fractions), youngest first.
    """
    z_initial = ZINI_FRACTION * z_final
    # Bins run youngest first, so the mass formed before a bin is the sum of the older ones
    mass_before = np.cumsum(fractions[::-1])[::-1] - fractions
    enrichment = (mass_before + fractions / 2) / fractions.sum()
    return z_initial + (z_final - z_initial) * enrichment


def components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one unit-cube vector u (order of DIMS).

    Inputs: u: unit-cube vector (dimensionless, in [0, 1)) ordered as DIMS;
            fixed_redshift: if given, replaces the redshift drawn from u (spectral fits); dimensionless z.
    Returns: (bagpipes_components, drawn). drawn holds redshift (dimensionless), sfh_type (0, 1 or 2),
             Z_fin (Z/Z_sun), Av_ism and Av_bc (mag), logU (dimensionless).
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

    drawn = {"redshift": redshift, "sfh_type": sfh_type, "Z_fin": z_final,
             "Av_ism": av_ism, "Av_bc": av_birth_cloud, "logU": log_u}
    return bagpipes_components, drawn


def bins_from_components(bagpipes_components):
    """Mass fractions, metallicities and edges [Gyr ago] read back from the constant components.

    Inputs: bagpipes_components: dict from components() with constant* entries (massformed in log M_sun).
    Returns: (fractions, metallicities, edges): mass fraction per bin (dimensionless), Z/Z_sun per bin,
             and n_bins + 1 bin edges [Gyr ago], youngest first.
    """
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
MAGNITUDES_PER_DEX = 2.5         # [mag/dex] magnitude change per decade of flux
DEX_PER_MAGNITUDE = 1 / MAGNITUDES_PER_DEX   # [dex/mag] = 0.4: flux ratio is 10**(-0.4 * magnitude)
AB_ZERO_POINT = 48.60            # [mag] AB offset for flux density f_nu in erg/s/cm2/Hz

# S-PLUS bands, in the order of splus_filters.txt
BANDS = ["u", "J0378", "J0395", "J0410", "J0430", "g", "J0515", "r", "J0660", "i", "J0861", "z"]

FILTER_FILES = []
for relative_path in np.loadtxt(FILTER_LIST_FILE, dtype="str"):
    FILTER_FILES.append(os.path.join(MODULE_DIRECTORY, relative_path))


def pivot_wavelength(filter_file):
    """Pivot wavelength [Angstrom] of a filter curve file (wavelength, transmission).

    Inputs: filter_file: path to a text file with wavelength [Angstrom] and transmission (dimensionless) columns.
    Returns: pivot wavelength [Angstrom] as a float.
    """
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


def rescale_mass_labels(table, magnitude_offset):
    """Copy of a label table with the mass-like labels moved to the brightness of an observed galaxy.

    A galaxy fainter than its model by magnitude_offset has 10**(-0.4 * offset) times the mass,
    so every MASS_DEPENDENT_LABELS column changes by -0.4 * offset dex.
    Inputs: table: DataFrame with the MASS_DEPENDENT_LABELS columns [log M_sun or log M_sun/yr];
            magnitude_offset: m_observed - m_model [mag], one number for every row or one per row.
    Returns: new DataFrame, same columns and index; the mass-like columns shifted, the others
             unchanged. A NaN offset gives NaN labels in that row.
    """
    shifted = table.copy()
    offset = np.asarray(magnitude_offset, dtype=float)
    for label in MASS_DEPENDENT_LABELS:
        shifted[label] = table[label].to_numpy(dtype=float) - DEX_PER_MAGNITUDE * offset
    return shifted


def bin_labels(fractions, metallicities, edges, log_mass_formed=LOG_MASS_FORMED):
    """Age-bin labels (N_BINS columns of each kind), youngest bin first.

    Inputs: fractions: mass formed per existing bin (any normalisation, dimensionless);
            metallicities: Z/Z_sun per existing bin; edges: bin edges [Gyr ago];
            log_mass_formed: total mass formed [log M_sun].
    Returns: dict ordered as BIN_LABEL_COLUMNS: frac (dimensionless, sums to 1), logsfr [log M_sun/yr]
             and zmet (Z/Z_sun). Bins that do not exist at this redshift get frac 0 and NaN for logsfr and zmet.
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
    """Summary labels of a Bagpipes star formation history.

    Inputs: bagpipes_sfh: fitted model's sfh object (ages in yr, sfh in M_sun/yr, stellar_mass in log M_sun).
    Returns: dict ordered as SUMMARY_LABEL_COLUMNS: log_mstar [log M_sun] (living stars), log_sfr [log M_sun/yr]
             (last 100 Myr), mw_age [Gyr], and dimensionless f_young (< 0.1 Gyr), f_mid (0.1-1 Gyr), f_old (> 1 Gyr).
    """
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
    """Mean Z/Z_sun weighted by the mass formed in each bin.

    Inputs: fractions: mass formed per bin (any normalisation); metallicities: Z/Z_sun per bin.
    Returns: mass-weighted Z/Z_sun (dimensionless).
    """
    return np.sum(fractions * metallicities) / np.sum(fractions)


def ab_magnitudes(flux_per_angstrom):
    """AB magnitudes from Bagpipes photometry [erg/s/cm2/A]; all NaN if any flux is not positive.

    Inputs: flux_per_angstrom: F_lambda [erg/s/cm2/Angstrom], one per filter in PIVOT_WAVELENGTHS.
    Returns: AB magnitude per filter [mag], in the order of PIVOT_WAVELENGTHS.
    """
    if not np.all(flux_per_angstrom > 0):
        return np.full(len(PIVOT_WAVELENGTHS), np.nan)
    flux_per_hertz = flux_per_angstrom * PIVOT_WAVELENGTHS**2 / SPEED_OF_LIGHT
    return -MAGNITUDES_PER_DEX * np.log10(flux_per_hertz) - AB_ZERO_POINT


def simulate_models(unit_vectors):
    """Drawn parameters, labels and AB magnitudes (10^10 M_sun formed) for a block of unit-cube vectors.

    Inputs: unit_vectors: iterable of unit-cube vectors (dimensionless, in [0, 1)), each ordered as DIMS.
    Returns: DataFrame with one row per model: the drawn parameters (redshift, sfh_type, Z_fin, ...),
             SUMMARY_LABEL_COLUMNS, BIN_LABEL_COLUMNS and the AB magnitudes [mag] in BANDS order
             (NaN where a flux is not positive). The "z" column is the z-band magnitude, not the redshift.

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


# =============================================================================
# 3. Detection rule and S-PLUS noise model
# =============================================================================
UNDETECTED_MAGNITUDE = 50.0      # S-PLUS writes 99 for non-detections
MAX_DETECTED_ERROR = 1.0         # [mag] error >= 1 mag (S/N < 1) counts as a non-detection
ERROR_FLOOR = 0.05               # [mag] added in quadrature; stands in for model error
MIN_DETECTED_BANDS = 6

Z_SIGMA = 2.0                    # redshift window: |z_model - zml| < Z_SIGMA * zml_error,
DZ_MIN, DZ_MAX = 0.02, 0.10      # limited to [DZ_MIN, DZ_MAX]
SIGMA_Z = 0.025                  # assumed zml_error = SIGMA_Z * (1 + zml)

NOISE_MAGNITUDE_RANGE = (10.0, 30.0)   # [mag] range of the noise-model bins
NOISE_BIN_WIDTH = 0.25                 # [mag]
NOISE_MIN_GALAXIES = 20                # bins with fewer galaxies are not used
NOISE_N_QUANTILES = 21                 # quantiles of log10(error) stored per bin
MAGNITUDE_TO_RELATIVE_FLUX = MAGNITUDES_PER_DEX / np.log(10)   # sigma_mag = 1.0857 * sigma_flux / flux
SCALING_MAX_DZ = 0.05                  # redshift window for borrowing a real r magnitude


def detected_bands(magnitudes, errors):
    """True where a band is detected: magnitude measured and 0 < error < MAX_DETECTED_ERROR (NaN -> False).

    Inputs: magnitudes: [mag] array (n, 12), non-detections at or above UNDETECTED_MAGNITUDE;
            errors: [mag] array of the same shape.
    Returns: boolean array of the same shape, True where the band is detected.
    """
    measured = magnitudes < UNDETECTED_MAGNITUDE
    good_error = (errors > 0) & (errors < MAX_DETECTED_ERROR)
    return measured & good_error


def redshift_window(photometric_z):
    """Half-width of the model redshift window around each photometric redshift.

    Inputs: photometric_z: array of photometric redshifts (dimensionless).
    Returns: half-width dz (dimensionless) per galaxy, clipped to [DZ_MIN, DZ_MAX].
    """
    assumed_error = SIGMA_Z * (1 + photometric_z)
    return np.clip(Z_SIGMA * assumed_error, DZ_MIN, DZ_MAX)


def noise_bin_index(magnitudes):
    """Index of the noise-model bin of each magnitude (may fall outside the table).

    Inputs: magnitudes: finite magnitudes [mag] (array).
    Returns: integer bin index per magnitude, NOISE_BIN_WIDTH wide from NOISE_MAGNITUDE_RANGE[0];
             can be negative or beyond the last bin.
    """
    return np.floor((magnitudes - NOISE_MAGNITUDE_RANGE[0]) / NOISE_BIN_WIDTH).astype(int)


def quantile_column_names():
    """Names of the stored log-error quantile columns: log_error_q00 ... log_error_q20.

    Inputs: none.
    Returns: list of NOISE_N_QUANTILES column names (str).
    """
    names = []
    for quantile_index in range(NOISE_N_QUANTILES):
        names.append(f"log_error_q{quantile_index:02d}")
    return names


def fit_noise_model(catalogue):
    """Empirical error distribution per band and 0.25 mag bin of measured magnitude.

    Inputs: catalogue: S-PLUS DataFrame with <band>_auto magnitudes [mag] and e_<band>_auto errors [mag].
    Returns: DataFrame with one row per populated (band, bin): band, magnitude_low and magnitude_high [mag],
             n_galaxies, and 21 quantiles of log10(error [mag]) (log_error_q00 ... log_error_q20).

    A band is later detected when its drawn error is below MAX_DETECTED_ERROR, so each bin's
    detection probability is the share of its errors below 1 mag.
    """
    quantile_levels = np.linspace(0, 1, NOISE_N_QUANTILES)
    n_bins = int(round((NOISE_MAGNITUDE_RANGE[1] - NOISE_MAGNITUDE_RANGE[0]) / NOISE_BIN_WIDTH))
    rows = []
    for band in BANDS:
        magnitudes = catalogue[f"{band}_auto"].to_numpy()
        errors = catalogue[f"e_{band}_auto"].to_numpy()
        measured = (magnitudes > 0) & (magnitudes < UNDETECTED_MAGNITUDE) & (errors > 0)
        bin_of_galaxy = noise_bin_index(magnitudes[measured])
        log_errors = np.log10(errors[measured])

        for bin_index in range(n_bins):
            in_bin = bin_of_galaxy == bin_index
            if in_bin.sum() < NOISE_MIN_GALAXIES:
                continue
            row = {
                "band": band,
                "magnitude_low": NOISE_MAGNITUDE_RANGE[0] + bin_index * NOISE_BIN_WIDTH,
                "magnitude_high": NOISE_MAGNITUDE_RANGE[0] + (bin_index + 1) * NOISE_BIN_WIDTH,
                "n_galaxies": int(in_bin.sum()),
            }
            quantiles = np.quantile(log_errors[in_bin], quantile_levels)
            for name, value in zip(quantile_column_names(), quantiles):
                row[name] = value
            rows.append(row)
    return pd.DataFrame(rows)


def draw_log_errors(magnitudes, band_table, rng):
    """Draw log10(error) for each magnitude from its bin's quantiles (inverse-CDF sampling).

    Inputs: magnitudes: model magnitudes [mag] for one band (array, NaN allowed);
            band_table: rows of the noise table for that band; rng: numpy Generator.
    Returns: log10(error [mag]) per magnitude; NaN where the magnitude is NaN or fainter than the last populated bin.

    Magnitudes brighter than the first populated bin use that bin.
    """
    table_bins = noise_bin_index(band_table["magnitude_low"].to_numpy() + NOISE_BIN_WIDTH / 2)
    quantile_table = band_table[quantile_column_names()].to_numpy()

    log_errors = np.full(len(magnitudes), np.nan)
    finite = np.isfinite(magnitudes)
    model_bins = np.full(len(magnitudes), -1)
    model_bins[finite] = noise_bin_index(magnitudes[finite])

    for model_index in np.flatnonzero(finite):
        if model_bins[model_index] > table_bins[-1]:
            continue    # fainter than any S-PLUS galaxy with a measurement
        # Nearest populated bin at or below this one; the brightest bin covers anything brighter
        row_index = np.searchsorted(table_bins, model_bins[model_index], side="right") - 1
        row_index = max(row_index, 0)
        quantile_position = rng.random() * (NOISE_N_QUANTILES - 1)
        lower = int(np.floor(quantile_position))
        upper = min(lower + 1, NOISE_N_QUANTILES - 1)
        weight = quantile_position - lower
        row = quantile_table[row_index]
        log_errors[model_index] = row[lower] + weight * (row[upper] - row[lower])
    return log_errors


def add_splus_noise(magnitudes, noise_table, rng):
    """S-PLUS-like errors, flux noise and non-detections for model magnitudes (n_models, 12).

    Inputs: magnitudes: noiseless model magnitudes [mag], shape (n_models, 12) in BANDS order;
            noise_table: output of fit_noise_model; rng: numpy Generator.
    Returns: (noisy_magnitudes [mag], errors [mag], detected), each (n_models, 12); detected is a bool array.
             Undetected entries are NaN in the first two.
    """
    n_models = len(magnitudes)
    noisy_magnitudes = np.full(magnitudes.shape, np.nan)
    errors = np.full(magnitudes.shape, np.nan)
    detected = np.zeros(magnitudes.shape, dtype=bool)

    for band_index, band in enumerate(BANDS):
        band_table = noise_table[noise_table["band"] == band]
        true_magnitudes = magnitudes[:, band_index]
        band_errors = 10 ** draw_log_errors(true_magnitudes, band_table, rng)

        # Gaussian noise in flux, converted back to magnitudes
        relative_flux_error = band_errors / MAGNITUDE_TO_RELATIVE_FLUX
        flux_ratio = 1 + rng.normal(0, 1, n_models) * relative_flux_error
        positive_flux = flux_ratio > 0
        band_detected = (band_errors < MAX_DETECTED_ERROR) & positive_flux

        noisy_band = np.full(n_models, np.nan)
        noisy_band[band_detected] = true_magnitudes[band_detected] - MAGNITUDES_PER_DEX * np.log10(flux_ratio[band_detected])
        noisy_magnitudes[:, band_index] = noisy_band
        errors[band_detected, band_index] = band_errors[band_detected]
        detected[:, band_index] = band_detected
    return noisy_magnitudes, errors, detected


def scale_to_observed_r(model_magnitudes, model_z, catalogue_r, catalogue_z, rng, max_dz=SCALING_MAX_DZ):
    """Shift each model to the r magnitude of a random real galaxy with |zml - z| < max_dz.

    Inputs: model_magnitudes: [mag] array (n_models, 12) in BANDS order; model_z: model redshifts (dimensionless);
            catalogue_r: real r magnitudes [mag]; catalogue_z: photometric redshifts of those galaxies (dimensionless);
            rng: numpy Generator; max_dz: redshift half-window (dimensionless).
    Returns: (scaled_magnitudes [mag], offsets [mag]); offsets is the shift applied to every band (real r minus
             model r) and both are NaN where no real galaxy lies within max_dz.

    Gives synthetic galaxies a realistic brightness (and so realistic errors) at each redshift.
    """
    r_index = BANDS.index("r")
    usable = np.isfinite(catalogue_r) & (catalogue_r < UNDETECTED_MAGNITUDE) & np.isfinite(catalogue_z)
    order = np.argsort(catalogue_z[usable])
    sorted_z = catalogue_z[usable][order]
    sorted_r = catalogue_r[usable][order]

    offsets = np.full(len(model_magnitudes), np.nan)
    for model_index in range(len(model_magnitudes)):
        first = np.searchsorted(sorted_z, model_z[model_index] - max_dz)
        last = np.searchsorted(sorted_z, model_z[model_index] + max_dz)
        if last == first:
            continue
        chosen = rng.integers(first, last)
        offsets[model_index] = sorted_r[chosen] - model_magnitudes[model_index, r_index]

    scaled_magnitudes = model_magnitudes + offsets[:, np.newaxis]
    return scaled_magnitudes, offsets


# =============================================================================
# 4. Coverage test v2
# =============================================================================
P_MIN = 0.01                 # probability (dimensionless): covered if the best p-value exceeds this
NO_MODEL = -1                # best_model value when no model lies in the redshift window
FREE_OFFSETS = 1             # the magnitude normalisation is fitted, so it uses one degree of freedom


def coverage_test(observed_magnitudes, observed_errors, detected, observed_z, observed_dz,
                  model_magnitudes, model_z):
    """Is each galaxy consistent with some model at its redshift?

    Purpose: for every model with |z_model - z| < dz, chi^2 over the detected bands, with the magnitude
    offset (the model's mass normalisation) solved analytically. A galaxy is covered if its minimum chi^2
    has p-value > P_MIN for n_detected - FREE_OFFSETS degrees of freedom.
    Inputs: observed_magnitudes and observed_errors [mag], arrays (n_galaxies, 12) in BANDS order;
            detected: bool array (n_galaxies, 12); observed_z and observed_dz: dimensionless, shape (n_galaxies,);
            model_magnitudes [mag], array (n_models, 12); model_z: dimensionless, shape (n_models,).
    Returns: DataFrame with one row per galaxy and columns chi2_min, n_dof, p_value, covered (bool) and
             best_model (index into model_magnitudes; NO_MODEL when no model is in the window).
    """
    # 1. Models with finite magnitudes, sorted by redshift so each galaxy's window is a contiguous slice
    finite_models = np.flatnonzero(np.isfinite(model_magnitudes).all(axis=1))
    order = finite_models[np.argsort(model_z[finite_models])]
    sorted_z = model_z[order]
    sorted_magnitudes = model_magnitudes[order]

    rows = []
    for galaxy_index in range(len(observed_magnitudes)):
        bands = np.flatnonzero(detected[galaxy_index])
        n_dof = len(bands) - FREE_OFFSETS
        first = np.searchsorted(sorted_z, observed_z[galaxy_index] - observed_dz[galaxy_index])
        last = np.searchsorted(sorted_z, observed_z[galaxy_index] + observed_dz[galaxy_index])
        if last == first or n_dof < 1:    # a chi^2 needs at least one degree of freedom
            rows.append({"chi2_min": np.nan, "n_dof": n_dof, "p_value": np.nan,
                         "covered": False, "best_model": NO_MODEL})
            continue

        # 2. chi^2 with free offset: sum w d^2 - (sum w d)^2 / sum w, where d = m_obs - m_model.
        # Solving for the offset analytically avoids a grid over normalisations.
        # The variance is the error squared plus ERROR_FLOOR squared (model error in quadrature).
        weights = 1 / (observed_errors[galaxy_index, bands] ** 2 + ERROR_FLOOR**2)
        window_magnitudes = sorted_magnitudes[first:last]
        differences = observed_magnitudes[galaxy_index, bands] - window_magnitudes[:, bands]
        weighted_sum = differences @ weights
        best_offsets = weighted_sum / weights.sum()
        weighted_square_sum = (differences**2) @ weights
        chi2_values = weighted_square_sum - weighted_sum * best_offsets

        # 3. Best model in the window and its p-value
        best_in_window = np.argmin(chi2_values)
        chi2_min = chi2_values[best_in_window]
        p_value = chi_squared.sf(chi2_min, n_dof)
        rows.append({"chi2_min": chi2_min, "n_dof": n_dof, "p_value": p_value,
                     "covered": p_value > P_MIN, "best_model": int(order[first + best_in_window])})
    return pd.DataFrame(rows)


# =============================================================================
# 5. Spectra, spectral likelihood, Firefly comparison and split
# =============================================================================
R_V = 3.1                                # [dimensionless] CCM89 R_V = A_V / E(B-V), as in extinctions.py
GYR_TO_YEARS = 1e9                       # [yr/Gyr]

# S-PLUS wavelengths [Angstrom] at which the CCM89 law is evaluated, as in extinctions.py
EXTINCTION_WAVELENGTHS = np.array([3536, 3770, 3940, 4094, 4292, 4751, 5133, 6258, 6614, 7690, 8611, 8831],
                                  dtype=float)

# Refractive index of air (Morton 2000, IAU standard); s is the wavenumber in 1/micron
ANGSTROM_PER_MICRON = 1e4                # [Angstrom/micron]
AIR_INDEX_OFFSET = 8.34254e-5            # dimensionless constant part of (n - 1)
AIR_INDEX_NUMERATOR_1 = 0.02406147       # dimensionless
AIR_INDEX_POLE_1 = 130.0                 # [1/micron^2] pole of the first term in s^2
AIR_INDEX_NUMERATOR_2 = 0.00015998       # dimensionless
AIR_INDEX_POLE_2 = 38.9                  # [1/micron^2] pole of the second term in s^2

OBSERVED_RANGE = (3800.0, 9200.0)        # [Angstrom, observed] kept for the fit
SKY_LINES = (5577.0, 6300.0)             # [Angstrom] strong sky residuals
SKY_LINE_HALF_WIDTH = 10.0               # [Angstrom]
REBIN_WIDTH = 5.0                        # [Angstrom]
MIN_GOOD_BINS = 300                      # bins with data, out of about 1080
CALIBRATION_ORDER = 2                    # degree of the multiplicative calibration polynomial

VELDISP_RANGE = (40.0, 400.0)                  # [km/s]
LOG_NOISE_SCALE_RANGE = (0.0, np.log10(5.0))   # [dex] log10 of the error multiplier, 1 to 5
IGNORED_REDSHIFT_DRAW = 0.5              # unit-cube value for the redshift slot, which fixed_redshift replaces

FIREFLY_MISSING = -999.0                 # [same units as Firefly] Firefly writes -9999 for missing values
MAX_DELTA_LOG_AGE = 0.3                  # [dex]
MAX_DELTA_LOG_METALLICITY = 0.3          # [dex]

LOG_ZERO = -1e100                        # log-likelihood of an impossible point (MultiNest convention)
HASH_HEX_DIGITS = 8                      # digits of the md5 digest used for the split (32 bits)
SPLIT_MODULUS = 2                        # even hash goes to validation, odd to test: half each


def vacuum_to_air(wavelength):
    """Air wavelengths [Angstrom] from vacuum ones (IAU standard, Morton 2000).

    SDSS spectra are in vacuum; the BC03-MILES models are in air.
    Inputs: wavelength: vacuum wavelengths [Angstrom] (array).
    Returns: air wavelengths [Angstrom], same shape as wavelength.
    """
    # The index of air departs from 1 by an amount that grows as the wavelength shortens
    wavenumber_squared = (ANGSTROM_PER_MICRON / wavelength) ** 2
    first_term = AIR_INDEX_NUMERATOR_1 / (AIR_INDEX_POLE_1 - wavenumber_squared)
    second_term = AIR_INDEX_NUMERATOR_2 / (AIR_INDEX_POLE_2 - wavenumber_squared)
    refractive_index = 1 + AIR_INDEX_OFFSET + first_term + second_term
    return wavelength / refractive_index


def deredden_spectrum(wavelength, flux, ivar, ebv):
    """Remove Milky Way extinction (CCM89, R_V = 3.1) from a spectrum.

    Inputs: wavelength [Angstrom], flux and ivar (inverse variance, 1/flux^2) as arrays of equal length;
            ebv: Milky Way E(B-V) [mag] at the galaxy.
    Returns: (flux, ivar), both dereddened, in the same units as the inputs.
    """
    extinction_magnitudes = extinction.ccm89(wavelength.astype(float), R_V * ebv, R_V)
    # Observed flux = intrinsic flux * 10**(-0.4 A), so the correction multiplies by the inverse
    correction = 10 ** (DEX_PER_MAGNITUDE * extinction_magnitudes)
    # Errors scale with flux, so the inverse variance divides by the square of the correction
    return flux * correction, ivar / correction**2


def deredden_magnitudes(magnitudes, ebv):
    """Remove Milky Way extinction from S-PLUS magnitudes, as extinctions.py does.

    Inputs: magnitudes: observed [mag], array (n_galaxies, 12) in EXTINCTION_WAVELENGTHS order;
            ebv: Milky Way E(B-V) [mag] per galaxy, shape (n_galaxies,).
    Returns: dereddened magnitudes [mag], same shape as magnitudes.
    """
    corrected = np.array(magnitudes, dtype=float)
    for galaxy_index in range(len(corrected)):
        # Magnitudes are brighter than intrinsic by A_lambda, so the extinction is subtracted
        band_extinction = extinction.ccm89(EXTINCTION_WAVELENGTHS, R_V * ebv[galaxy_index], R_V)
        corrected[galaxy_index] = corrected[galaxy_index] - band_extinction
    return corrected


def prepare_spectrum(wavelength_vacuum, flux, ivar, ebv):
    """SDSS spectrum ready for the fit: air wavelengths, dereddened, masked, rebinned to REBIN_WIDTH.

    Each output bin is the inverse-variance-weighted mean of its unmasked pixels; bins without any
    get ivar = 0. Raises ValueError when fewer than MIN_GOOD_BINS bins have data.
    Inputs: wavelength_vacuum [Angstrom, vacuum]; flux [erg/s/cm2/Angstrom, any scale] and ivar
            (inverse variance) per pixel, arrays of equal length; ebv: Milky Way E(B-V) [mag].
    Returns: dict(wavelength [Angstrom, air, observed bin centres], flux, ivar), one value per bin.
    """
    # 1. Air wavelengths and Milky Way extinction
    wavelength = vacuum_to_air(wavelength_vacuum)
    flux, ivar = deredden_spectrum(wavelength, flux, ivar, ebv)

    # 2. Pixels to ignore: no data, non-finite flux, and strong sky residuals
    good = (ivar > 0) & np.isfinite(flux)
    for sky_line in SKY_LINES:
        near_sky_line = np.abs(wavelength - sky_line) < SKY_LINE_HALF_WIDTH
        good = good & ~near_sky_line

    # 3. Rebin onto a regular grid, weighting each pixel by its inverse variance
    bin_centres = np.arange(OBSERVED_RANGE[0] + REBIN_WIDTH / 2, OBSERVED_RANGE[1], REBIN_WIDTH)
    rebinned_flux = np.zeros(len(bin_centres))
    rebinned_ivar = np.zeros(len(bin_centres))
    pixel_bin = np.floor((wavelength - OBSERVED_RANGE[0]) / REBIN_WIDTH).astype(int)
    for bin_index in range(len(bin_centres)):
        in_bin = good & (pixel_bin == bin_index)
        if not in_bin.any():
            continue
        weights = ivar[in_bin]
        rebinned_ivar[bin_index] = weights.sum()
        rebinned_flux[bin_index] = np.sum(weights * flux[in_bin]) / weights.sum()

    # 4. Refuse spectra with too little usable wavelength range
    n_good_bins = int((rebinned_ivar > 0).sum())
    if n_good_bins < MIN_GOOD_BINS:
        raise ValueError(f"only {n_good_bins} good bins (need {MIN_GOOD_BINS})")
    return {"wavelength": bin_centres, "flux": rebinned_flux, "ivar": rebinned_ivar}


def calibration_polynomial(model_flux, data_flux, data_ivar, order=CALIBRATION_ORDER):
    """Multiplicative polynomial (in normalised wavelength) that best maps the model onto the data.

    Weighted linear least squares; the constant term absorbs the model normalisation, so the
    spectrum constrains the shape of the SFH, not the mass.
    Inputs: model_flux, data_flux: flux per bin (arrays of equal length, any common scale);
            data_ivar: inverse variance per bin (1/flux^2); order: polynomial degree (dimensionless).
    Returns: the polynomial evaluated on each bin (dimensionless), same length as model_flux.
    """
    n_bins = len(model_flux)
    # Wavelength axis rescaled to [-1, 1] keeps the polynomial terms well conditioned
    x = np.linspace(-1, 1, n_bins)
    design = np.zeros((n_bins, order + 1))
    for power in range(order + 1):
        design[:, power] = model_flux * x**power

    weighted_design = design * data_ivar[:, np.newaxis]
    normal_matrix = design.T @ weighted_design
    right_hand_side = weighted_design.T @ data_flux
    coefficients = np.linalg.solve(normal_matrix, right_hand_side)

    polynomial = np.zeros(n_bins)
    for power in range(order + 1):
        polynomial = polynomial + coefficients[power] * x**power
    return polynomial


def spectral_log_likelihood(model_flux, data_flux, data_ivar, noise_scale):
    """Gaussian log-likelihood of the data given the calibrated model, with errors x noise_scale.

    The calibration polynomial is refitted for each model. It does not depend on noise_scale,
    because scaling every weight by the same factor leaves the best fit unchanged. Terms that
    do not depend on the model are dropped.
    Inputs: model_flux, data_flux, data_ivar: as in calibration_polynomial;
            noise_scale: error multiplier (dimensionless, positive).
    Returns: log-likelihood (dimensionless) summed over bins with data_ivar > 0.
    """
    polynomial = calibration_polynomial(model_flux, data_flux, data_ivar)
    residuals = data_flux - model_flux * polynomial
    good = data_ivar > 0
    chi2_value = np.sum(residuals[good] ** 2 * data_ivar[good]) / noise_scale**2
    return -0.5 * chi2_value - good.sum() * np.log(noise_scale)


def safe_spectral_log_likelihood(model_flux, data_flux, data_ivar, noise_scale):
    """spectral_log_likelihood that returns LOG_ZERO instead of raising for an impossible model.

    PyMultiNest calls sys.exit if the likelihood raises, which would kill a whole fitting batch.
    A model with non-finite flux, or one whose calibration fit is singular (for example all zeros),
    is an impossible point, not an error.
    Inputs: as in spectral_log_likelihood.
    Returns: log-likelihood (dimensionless), or LOG_ZERO for an impossible model.
    """
    if not np.all(np.isfinite(model_flux)):
        return LOG_ZERO
    try:
        return spectral_log_likelihood(model_flux, data_flux, data_ivar, noise_scale)
    except np.linalg.LinAlgError:
        return LOG_ZERO


# Redshift is fixed to the spectroscopic value, so it leaves the list; veldisp and noise_scale join it
SPECTRAL_DIMS = []
for dimension in DIMS:
    if dimension != "z":
        SPECTRAL_DIMS.append(dimension)
SPECTRAL_DIMS.append("veldisp")
SPECTRAL_DIMS.append("noise_scale")


def spectral_components(cube, redshift):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit unit-cube point.

    The SFH, metallicity, dust and nebular part is the training prior (components) at the fixed
    spectroscopic redshift; veldisp and noise_scale are the spectral nuisance parameters.
    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as SPECTRAL_DIMS;
            redshift: spectroscopic redshift (dimensionless).
    Returns: (bagpipes_components, drawn, noise_scale). The components are as in components() plus
             veldisp [km/s]; drawn is as in components() plus veldisp [km/s]; noise_scale is the error
             multiplier (dimensionless, 1 to 5).
    """
    cube_by_name = dict(zip(SPECTRAL_DIMS, cube))
    unit_vector = np.zeros(len(DIMS))
    for dimension_index, dimension in enumerate(DIMS):
        if dimension == "z":
            # components() ignores this value, because fixed_redshift replaces the drawn redshift
            unit_vector[dimension_index] = IGNORED_REDSHIFT_DRAW
        else:
            unit_vector[dimension_index] = cube_by_name[dimension]

    bagpipes_components, drawn = components(unit_vector, fixed_redshift=redshift)
    veldisp = uniform(cube_by_name["veldisp"], *VELDISP_RANGE)
    bagpipes_components["veldisp"] = veldisp
    drawn["veldisp"] = veldisp
    # Uniform in log10 so that small and large error multipliers are equally likely per decade
    noise_scale = 10 ** uniform(cube_by_name["noise_scale"], *LOG_NOISE_SCALE_RANGE)
    return bagpipes_components, drawn, noise_scale


def compare_with_firefly(bagpipes_age_gyr, bagpipes_metallicity, firefly_age_years, firefly_metallicity):
    """Bagpipes vs Firefly mass-weighted age and metallicity (Z/Z_sun, both linear).

    Firefly missing values (-9999) become NaN, and the consistency flag is then NA, not False.
    Inputs: bagpipes_age_gyr [Gyr]; bagpipes_metallicity [Z/Z_sun]; firefly_age_years [yr];
            firefly_metallicity [Z/Z_sun]; arrays of equal length.
    Returns: DataFrame with delta_log_age and delta_log_metallicity [dex] (Bagpipes minus Firefly)
             and firefly_consistent (pandas boolean, NA where Firefly is missing).
    """
    firefly_age_years = np.where(firefly_age_years > FIREFLY_MISSING, firefly_age_years, np.nan)
    firefly_metallicity = np.where(firefly_metallicity > FIREFLY_MISSING, firefly_metallicity, np.nan)

    delta_log_age = np.log10(bagpipes_age_gyr * GYR_TO_YEARS) - np.log10(firefly_age_years)
    delta_log_metallicity = np.log10(bagpipes_metallicity) - np.log10(firefly_metallicity)

    consistent = []
    for age_offset, metallicity_offset in zip(delta_log_age, delta_log_metallicity):
        if np.isnan(age_offset) or np.isnan(metallicity_offset):
            consistent.append(pd.NA)
        else:
            consistent.append(bool(abs(age_offset) < MAX_DELTA_LOG_AGE
                                   and abs(metallicity_offset) < MAX_DELTA_LOG_METALLICITY))

    return pd.DataFrame({
        "delta_log_age": delta_log_age,
        "delta_log_metallicity": delta_log_metallicity,
        "firefly_consistent": pd.array(consistent, dtype="boolean"),
    })


def split_by_hash(specobjid):
    """'validation' or 'test' from a hash of the SDSS specobjid.

    Depends only on the galaxy itself, so adding galaxies later never moves an existing one.
    Inputs: specobjid: SDSS spectroscopic object identifier (integer).
    Returns: the string "validation" or "test", about half of the galaxies each.
    """
    digest = hashlib.md5(str(int(specobjid)).encode()).hexdigest()
    hash_value = int(digest[:HASH_HEX_DIGITS], 16)     # base 16: the digest is hexadecimal
    if hash_value % SPLIT_MODULUS == 0:
        return "validation"
    return "test"


# =============================================================================
# 6. MultiNest fit, posterior labels, mocks and fit files
# =============================================================================
N_LIVE_POINTS = 300              # live points: more gives finer posteriors but a slower fit (count)
SAMPLING_EFFICIENCY = 0.8        # MultiNest parameter-estimation mode (dimensionless)
MAX_POSTERIOR_SAMPLES = 500      # posterior samples kept for labels (count)
MIDPOINT_DRAW = 0.5              # unit-cube value of every parameter, for the starting model (dimensionless)
UNIT_NOISE_CUBE_VALUE = 0.0      # unit-cube value of noise_scale that gives an error multiplier of 1
LOW_PERCENTILE = 16              # lower posterior summary (percent), 1 sigma
MEDIAN_PERCENTILE = 50           # central posterior summary (percent)
HIGH_PERCENTILE = 84             # upper posterior summary (percent), 1 sigma
MAX_ITERATIONS = 50000           # MultiNest iteration limit of a production fit (count)


def multinest_iteration_bounds(output_basename):
    """Lower and upper bound on the outer-loop iterations of a finished MultiNest run.

    MultiNest writes one row to ev.dat for every dead point, and with multimodal=True it removes one
    dead point per active mode in each iteration, so counting rows overstates the iterations as soon
    as the posterior splits into several modes. Inside one iteration the rows follow the mode order
    (the last column of ev.dat is the mode number, increasing). A row whose mode number is not larger
    than the previous row's certainly starts a new iteration, which gives the lower bound. The count
    can be missed when every mode of one iteration is replaced in the next (nested.F90, isolateModes2:
    new modes get the next free number and an emptied parent is dropped), e.g. rows 1 | 2 3: there the
    mode numbers only rise. Such a boundary is an ascending step k -> j into a mode number j that never
    appeared before, with k never occurring again (k is then dead); the upper bound counts each of
    those steps. The exact value is
    in between and is the loop counter that the max_iter limit is compared with in nested.F90 (no
    MultiNest file stores it).
    Inputs: output_basename: directory prefix of the MultiNest files (str), as given to fit_spectrum.
    Returns: (lower, upper) iteration counts (counts); (0, 0) if ev.dat has no rows.
    """
    mode_of_each_row = np.loadtxt(output_basename + "ev.dat", usecols=-1, ndmin=1)
    if len(mode_of_each_row) == 0:
        return 0, 0
    certain_starts = int((mode_of_each_row[1:] <= mode_of_each_row[:-1]).sum())
    possible_starts = 0
    modes_seen = {mode_of_each_row[0]}
    for row_index in range(1, len(mode_of_each_row)):
        previous_mode, mode = mode_of_each_row[row_index - 1], mode_of_each_row[row_index]
        previous_mode_is_dead = not (mode_of_each_row[row_index:] == previous_mode).any()
        if mode > previous_mode and mode not in modes_seen and previous_mode_is_dead:
            possible_starts += 1
        modes_seen.add(mode)
    return 1 + certain_starts, 1 + certain_starts + possible_starts


def multinest_iterations(output_basename):
    """Best (lower-bound) count of outer-loop iterations of a finished MultiNest run.

    Exact unless a whole set of modes is replaced at once (see multinest_iteration_bounds); use the
    upper bound from multinest_iteration_bounds to decide whether the run hit max_iter.
    Inputs: output_basename: directory prefix of the MultiNest files (str), as given to fit_spectrum.
    Returns: iterations (count); 0 if ev.dat has no rows.
    """
    return multinest_iteration_bounds(output_basename)[0]


def fit_spectrum(spectrum, redshift, output_basename, max_iterations, n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one prepared spectrum under the training prior at a fixed redshift.

    The sampler works directly on the unit cube (SPECTRAL_DIMS): spectral_components maps each point
    to physical values, exactly as the training set does. output_basename must be short (MultiNest
    limit of 100 characters), e.g. "spectral_fits/<id>/".
    Inputs: spectrum: dict with wavelength [Angstrom, observed], flux and ivar, as from prepare_spectrum;
            redshift: spectroscopic redshift (dimensionless); output_basename: directory prefix for the
            MultiNest files (str); max_iterations: MultiNest iteration limit (count);
            n_live_points: live points (count).
    Returns: dict with samples (unit-cube points, shape (n_samples, len(SPECTRAL_DIMS))), log_evidence
             (dimensionless ln Z), n_iterations (count, from multinest_iterations), converged (bool: the
             run certainly stopped before max_iterations, judged by the upper bound of the count) and run_time [s].
    """
    # pymultinest loads the MultiNest C library on import, so it needs LD_LIBRARY_PATH (set by the kernel)
    import pymultinest

    os.makedirs(output_basename, exist_ok=True)

    # 1. Bagpipes model built once; each likelihood call updates it in place, which is much faster
    midpoint_cube = np.full(len(SPECTRAL_DIMS), MIDPOINT_DRAW)
    starting_components, _, _ = spectral_components(midpoint_cube, redshift)
    model = pipes.model_galaxy(starting_components, spec_wavs=spectrum["wavelength"])

    # 2. Likelihood of one unit-cube point: model spectrum vs data, with the S/N-scaling noise multiplier
    def log_likelihood(cube):
        """Log-likelihood of the spectrum for one unit-cube point.

        Inputs: cube: unit-cube point (dimensionless, SPECTRAL_DIMS order).
        Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
        """
        bagpipes_components, _, noise_scale = spectral_components(cube, redshift)
        model.update(bagpipes_components)
        model_flux = model.spectrum[:, 1]
        return safe_spectral_log_likelihood(model_flux, spectrum["flux"], spectrum["ivar"], noise_scale)

    def unit_cube_prior(cube):
        """Identity map: the training prior is already uniform on the unit cube.

        Inputs: cube: unit-cube point (dimensionless). Returns: the same point (dimensionless).
        """
        return cube

    # 3. Run MultiNest, timing only the sampling itself
    start = time.time()
    result = pymultinest.solve(
        LogLikelihood=log_likelihood, Prior=unit_cube_prior, n_dims=len(SPECTRAL_DIMS),
        outputfiles_basename=output_basename, n_live_points=n_live_points,
        sampling_efficiency=SAMPLING_EFFICIENCY, multimodal=True, max_iter=max_iterations,
        # Standard nested-sampling evidence: the importance-sampled one underflows in short runs and
        # MultiNest then writes it in a form pymultinest cannot parse
        importance_nested_sampling=False,
        resume=False, verbose=False, use_MPI=False,
    )
    run_time = time.time() - start

    # 4. Hitting the iteration limit means the run stopped before it converged; the upper bound of the
    #    count is used so that a capped run is never reported as converged
    n_iterations, max_possible_iterations = multinest_iteration_bounds(output_basename)
    return {
        "samples": result["samples"],
        "log_evidence": result["logZ"],
        "n_iterations": n_iterations,
        "converged": max_possible_iterations < max_iterations,
        "run_time": run_time,
    }


def posterior_labels(unit_cube_samples, redshift):
    """Labels of each posterior sample, from a Bagpipes model that also computes S-PLUS photometry.

    Keeps at most MAX_POSTERIOR_SAMPLES evenly spaced samples, because building labels is the slow part.
    model_r is the model r magnitude for 10^10 M_sun formed, used later to rescale mass-like labels.
    Inputs: unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, len(SPECTRAL_DIMS)));
            redshift: spectroscopic redshift (dimensionless).
    Returns: DataFrame with one row per kept sample: the drawn parameters (veldisp [km/s], Z_fin [Z/Z_sun], ...),
             noise_scale (dimensionless), SUMMARY_LABEL_COLUMNS, BIN_LABEL_COLUMNS, mw_metallicity (Z/Z_sun)
             and model_r [mag].
    """
    # 1. Thin the samples evenly: neighbouring MultiNest samples are nearly identical
    step = max(1, int(np.ceil(len(unit_cube_samples) / MAX_POSTERIOR_SAMPLES)))
    kept_samples = unit_cube_samples[::step]

    # 2. Labels per kept sample, with one photometry model updated in place
    model = None
    rows = []
    for cube in kept_samples:
        bagpipes_components, drawn, noise_scale = spectral_components(cube, redshift)
        if model is None:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
        else:
            model.update(bagpipes_components)
        fractions, metallicities, edges = bins_from_components(bagpipes_components)

        row = dict(drawn)
        row["noise_scale"] = noise_scale
        row.update(summary_labels(model.sfh))
        row.update(bin_labels(fractions, metallicities, edges))
        row["mw_metallicity"] = mass_weighted_metallicity(fractions, metallicities)
        row["model_r"] = ab_magnitudes(model.photometry)[BANDS.index("r")]
        rows.append(row)
    return pd.DataFrame(rows)


def make_mock_spectrum(cube, redshift, template_spectrum, rng):
    """Mock spectrum from a known unit-cube point, with the wavelengths, mask and errors of a real one.

    The model is scaled to the template's median flux so the mock has the template's S/N. The
    noise_scale coordinate is set to UNIT_NOISE_CUBE_VALUE (error multiplier 1) because the mock
    noise is the template's own.
    Inputs: cube: unit-cube point (dimensionless, SPECTRAL_DIMS order); redshift: dimensionless;
            template_spectrum: dict from prepare_spectrum (wavelength [Angstrom], flux, ivar);
            rng: numpy Generator.
    Returns: (spectrum dict with wavelength, flux and ivar, truth labels as a pandas Series with the
             same labels as posterior_labels).
    """
    # 1. Known model on the template's wavelength grid (copied, so the caller's point is not changed)
    cube = np.array(cube, dtype=float)
    cube[SPECTRAL_DIMS.index("noise_scale")] = UNIT_NOISE_CUBE_VALUE
    bagpipes_components, _, _ = spectral_components(cube, redshift)
    model = pipes.model_galaxy(bagpipes_components, spec_wavs=template_spectrum["wavelength"])
    model_flux = model.spectrum[:, 1]

    # 2. Match the template's flux scale, then add Gaussian noise with the template's errors
    good = template_spectrum["ivar"] > 0
    scale = np.median(template_spectrum["flux"][good]) / np.median(model_flux[good])
    errors = np.zeros(len(model_flux))
    errors[good] = 1 / np.sqrt(template_spectrum["ivar"][good])
    mock_flux = scale * model_flux + rng.normal(0, 1, len(model_flux)) * errors
    mock_flux[~good] = 0.0    # bins without data stay empty, as in the template

    # 3. Spectrum and truth labels of the same point, so a fit can be checked against them
    spectrum = {"wavelength": template_spectrum["wavelength"], "flux": mock_flux,
                "ivar": template_spectrum["ivar"].copy()}
    truth = posterior_labels(cube[np.newaxis, :], redshift).iloc[0]
    return spectrum, truth


def save_fit(path, posterior, unit_cube_samples, metadata):
    """Write posterior labels, unit-cube samples and fit metadata to one HDF5 file, all or nothing.

    The file is written as path + ".tmp" and renamed to path only once complete, so an interrupted
    or failed write never leaves a partial file at path (which would count as a finished job).
    Inputs: path: output file path (str); posterior: DataFrame of labels, one column per label;
            unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, len(SPECTRAL_DIMS)));
            metadata: dict of scalars (converged, run_time [s], n_iterations, max_iterations, log_evidence, redshift).
    Returns: None. An existing file at path is replaced.
    """
    temporary_path = path + ".tmp"
    with h5py.File(temporary_path, "w") as output:
        output.create_dataset("unit_cube_samples", data=unit_cube_samples)
        labels_group = output.create_group("posterior")
        for column in posterior.columns:
            labels_group.create_dataset(column, data=posterior[column].to_numpy(dtype=float))
        for key, value in metadata.items():
            output.attrs[key] = value
    os.replace(temporary_path, path)


def load_fit(path):
    """Read a file written by save_fit.

    Inputs: path: file path (str).
    Returns: (posterior DataFrame with one column per label, metadata dict with the attributes save_fit wrote).
    """
    with h5py.File(path, "r") as source:
        columns = {}
        for column in source["posterior"]:
            columns[column] = source["posterior"][column][()]
        metadata = {}
        for key, value in source.attrs.items():
            # HDF5 returns numpy scalars; .item() gives back the plain Python value
            metadata[key] = value.item() if hasattr(value, "item") else value
    return pd.DataFrame(columns), metadata


def fit_and_save(job, output_directory, max_iterations, n_live_points=N_LIVE_POINTS):
    """Fit one job and save <output_directory>/<job_id>.h5. Never raises.

    Skips jobs whose readable file already exists, so a batch can be re-run with more galaxies. An
    unreadable file does not count as finished: the job is fitted again and the file replaced. A failed
    fit comes back with status "failed" and the error text instead of an exception.
    Inputs: job: dict with job_id (str), spectrum (dict as from prepare_spectrum) and redshift (dimensionless);
            output_directory: folder for the .h5 file and the MultiNest files (str); max_iterations: MultiNest
            iteration limit (count); n_live_points: live points (count).
    Returns: dict with job_id, status ("fitted", "skipped" or "failed"), converged (bool or None),
             run_time [s] (or None), n_iterations (count or None) and error (str or None).
    """
    # 1. Output path and a summary that the steps below fill in
    output_path = os.path.join(output_directory, f"{job['job_id']}.h5")
    summary = {"job_id": job["job_id"], "status": None, "converged": None,
               "run_time": None, "n_iterations": None, "error": None}

    # 2. Skip finished jobs; their run details come from the file. Any failure to read it (damaged,
    # truncated or missing fields) means the job is not finished, so it is fitted again below
    if os.path.exists(output_path):
        try:
            _, metadata = load_fit(output_path)
            summary.update(status="skipped", converged=metadata["converged"],
                           run_time=metadata["run_time"], n_iterations=metadata["n_iterations"])
            return summary
        except Exception:
            pass

    # 3. Fit, label the posterior and save; any error is recorded rather than raised
    try:
        multinest_basename = os.path.join(output_directory, str(job["job_id"])) + "/"
        fit = fit_spectrum(job["spectrum"], job["redshift"], multinest_basename, max_iterations, n_live_points)
        posterior = posterior_labels(fit["samples"], job["redshift"])
        metadata = {"converged": fit["converged"], "run_time": fit["run_time"],
                    "n_iterations": fit["n_iterations"], "max_iterations": max_iterations,
                    "log_evidence": fit["log_evidence"], "redshift": job["redshift"]}
        save_fit(output_path, posterior, fit["samples"], metadata)
        summary.update(status="fitted", converged=fit["converged"],
                       run_time=fit["run_time"], n_iterations=fit["n_iterations"])
    except Exception as error:
        summary.update(status="failed", error=repr(error))
    return summary


def summarise_posterior(posterior, columns):
    """16th, 50th and 84th percentiles of each listed posterior column.

    Inputs: posterior: DataFrame of labels, one row per posterior sample; columns: names of the columns to summarise.
    Returns: dict with <column>_p16, <column>_median and <column>_p84 for each column, in the column's units.
    """
    summary = {}
    for column in columns:
        values = posterior[column].to_numpy(dtype=float)
        summary[f"{column}_p{LOW_PERCENTILE}"] = np.nanpercentile(values, LOW_PERCENTILE)
        summary[f"{column}_median"] = np.nanpercentile(values, MEDIAN_PERCENTILE)
        summary[f"{column}_p{HIGH_PERCENTILE}"] = np.nanpercentile(values, HIGH_PERCENTILE)
    return summary


# =============================================================================
# 7. Job files and the calibration test
# =============================================================================
N_CALIBRATION_MOCKS = 33         # mocks drawn from the prior to test the fitter (count)
MIN_COVERAGE_95 = 0.85           # smallest acceptable share of truths inside the 95% interval (fraction)
CALIBRATION_LABELS = ["mw_age", "Z_fin", "f_old"]   # labels the calibration test checks
CALIBRATION_LEVELS = (0.68, 0.95)                   # central-interval probabilities reported (fractions)
CALIBRATION_ID_PREFIX = "calib_"                    # job ids are calib_000, calib_001, ...
PERCENT = 100                    # percent per unit fraction


def export_job(job, jobs_directory):
    """Write one fit job to <jobs_directory>/<job_id>.npz, a portable file for run_spectral_fits.py.

    Inputs: job: dict with job_id (str), spectrum (dict with wavelength [Angstrom, observed], flux and
            ivar) and redshift (dimensionless); jobs_directory: folder for the file (str), created if missing.
    Returns: path of the written file (str).
    """
    os.makedirs(jobs_directory, exist_ok=True)
    path = os.path.join(jobs_directory, f"{job['job_id']}.npz")
    np.savez(path, wavelength=job["spectrum"]["wavelength"], flux=job["spectrum"]["flux"],
             ivar=job["spectrum"]["ivar"], redshift=job["redshift"], job_id=str(job["job_id"]))
    return path


def load_job(path):
    """Read a file written by export_job.

    Inputs: path: file path (str).
    Returns: job dict in the format fit_and_save takes: job_id (str), spectrum (wavelength, flux, ivar)
             and redshift (dimensionless float).
    """
    with np.load(path) as source:
        spectrum = {"wavelength": source["wavelength"], "flux": source["flux"], "ivar": source["ivar"]}
        return {"job_id": str(source["job_id"]), "spectrum": spectrum, "redshift": float(source["redshift"])}


def draw_calibration_cubes(n_mocks, rng):
    """Unit-cube points of the calibration mocks: uniform over every SPECTRAL_DIMS coordinate.

    sfh_type is drawn too, so the mocks follow the training prior and use all three SFH branches.
    Inputs: n_mocks: number of points (count); rng: numpy Generator.
    Returns: array of shape (n_mocks, len(SPECTRAL_DIMS)), values in [0, 1) (dimensionless).
    """
    return rng.random((n_mocks, len(SPECTRAL_DIMS)))


def draw_calibration_mocks(n_mocks, templates, rng):
    """Mock spectra with known labels, each on the wavelengths, mask and errors of a random real galaxy.

    Inputs: n_mocks: number of mocks (count); templates: list of (spectrum dict from prepare_spectrum,
            redshift (dimensionless)); rng: numpy Generator.
    Returns: (jobs, truths): jobs is a list of dicts as fit_and_save takes (job_id calib_000, calib_001, ...);
             truths maps each job_id to its true labels (Series, as posterior_labels columns).
    """
    cubes = draw_calibration_cubes(n_mocks, rng)
    jobs = []
    truths = {}
    for mock_index, cube in enumerate(cubes):
        job_id = f"{CALIBRATION_ID_PREFIX}{mock_index:03d}"
        template_spectrum, redshift = templates[rng.integers(len(templates))]
        spectrum, truth = make_mock_spectrum(cube, redshift, template_spectrum, rng)
        jobs.append({"job_id": job_id, "spectrum": spectrum, "redshift": redshift})
        truths[job_id] = truth
    return jobs, truths


def calibration_coverage(posteriors, truths, labels, levels=CALIBRATION_LEVELS):
    """Share of mocks whose true label lies inside the central posterior interval, per label and level.

    A calibrated fitter gives a share close to the level (0.95 of truths inside the 95% interval).
    A mock counts for a label only if its truth and its posterior samples for that label are finite.
    Inputs: posteriors: dict job_id -> posterior DataFrame (one row per sample); truths: dict job_id -> Series
            of true labels; labels: label names to check; levels: interval probabilities (fractions).
    Returns: DataFrame with one row per (label, level): label, level, coverage (fraction, NaN if no mock
             was usable) and n_jobs (count of mocks used).
    """
    rows = []
    for label in labels:
        for level in levels:
            half_width = PERCENT * level / 2
            lower_percentile = MEDIAN_PERCENTILE - half_width
            upper_percentile = MEDIAN_PERCENTILE + half_width
            inside = []
            for job_id, posterior in posteriors.items():
                samples = posterior[label].to_numpy(dtype=float)
                samples = samples[np.isfinite(samples)]
                if job_id not in truths or len(samples) == 0 or not np.isfinite(truths[job_id][label]):
                    continue
                lower, upper = np.percentile(samples, [lower_percentile, upper_percentile])
                inside.append(lower <= truths[job_id][label] <= upper)
            coverage = np.mean(inside) if inside else np.nan
            rows.append({"label": label, "level": level, "coverage": coverage, "n_jobs": len(inside)})
    return pd.DataFrame(rows)
