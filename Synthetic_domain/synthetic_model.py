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
    Returns: (bagpipes_components, drawn). drawn holds z (dimensionless), sfh_type (0, 1 or 2),
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

    drawn = {"z": redshift, "sfh_type": sfh_type, "Z_fin": z_final,
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


def bin_labels(fractions, metallicities, edges, log_mass_formed=LOG_MASS_FORMED):
    """16-bin labels, youngest bin first.

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
    Returns: DataFrame with one row per model: the drawn parameters, SUMMARY_LABEL_COLUMNS, BIN_LABEL_COLUMNS
             and the AB magnitudes [mag] in BANDS order (NaN where a flux is not positive).

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
