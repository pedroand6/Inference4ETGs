"""
Shared Bagpipes model of the S-PLUS synthetic domain.

training_set.ipynb and spectroscopic_reference.ipynb both import this module,
so the training set and the spectral fits use exactly the same prior.
"""

import hashlib
import os
import shutil
import time
import warnings

import extinction
import h5py
import joblib
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import simpson
from scipy.special import gammainc, gammaincc
from scipy.stats import chi2 as chi_squared

import bagpipes as pipes

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1. Prior (fixed before looking at the data)
# =============================================================================
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3)
UNIVERSE_AGE_MARGIN = 0.1        # [Gyr] age upper bound = Universe age at z minus this
LABEL_EDGES = [0.1, 1.0]         # [Gyr] f_young < 0.1 < f_mid < 1.0 < f_old
SFR_TIMESCALE = 0.1              # [Gyr] log_sfr is the mean SFR over this lookback time, as Bagpipes' sfr

# SFH: one Bagpipes "delayed" component, SFR(t) ~ t exp(-t / tau) since star formation began,
# with a single metallicity. Ages and tau are drawn log-uniform, because they span several decades.
Z_RANGE = (0.002, 0.7)                               # redshift
AGE_MIN = 0.05                                       # [Gyr] youngest age (time since star formation began)
TAU_RANGE = (0.01, 20.0)                             # [Gyr] e-folding time
METALLICITY_RANGE = (0.005, 2.5)                     # [Z/Z_sun] the bc03_miles grid range
AV_ISM_RANGE = (0.01, 2.0)                           # [mag]
LOG_TAU_BC_RANGE = (np.log10(0.5), np.log10(4.0))    # birth-cloud optical depth
LOGU_RANGE = (-4.0, -1.0)                            # nebular ionisation parameter
DUST_SLOPE = 0.7                                     # CF00 attenuation-law exponent
BIRTH_CLOUD_AGE = 0.01                               # [Gyr]
LOG_MASS_FORMED = 10.0                               # [log M_sun]
TAU_TO_AV = 1.086                                    # A_V = 1.086 * tau_V
SFR_FLOOR = 1e-14                # [yr^-1] sSFR floor for SFHs with no recent star formation

# "log_age_fraction" places the age inside [log AGE_MIN, log t_max(z)], so one cube works at every z
DIMS = ["z", "log_age_fraction", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU"]


def uniform(draw, low, high):
    """Map a unit-interval draw onto [low, high).

    Inputs: draw in [0, 1) (dimensionless); low, high: bounds in the units of the parameter.
    Returns: parameter value in [low, high), in the same units as low and high.
    """
    return low + (high - low) * draw


def maximum_age(redshift):
    """Oldest allowed age [Gyr]: the age of the Universe at the redshift minus UNIVERSE_AGE_MARGIN.

    Why the margin: Bagpipes flags a model that forms stars before the Big Bang as unphysical (NaN magnitudes).
    Inputs: redshift: dimensionless z, scalar or array.
    Returns: maximum age [Gyr], a float for a scalar redshift or an array of the same shape.
    """
    return COSMOLOGY.age(redshift).value - UNIVERSE_AGE_MARGIN


def components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one unit-cube vector u (order of DIMS).

    Inputs: u: unit-cube vector (dimensionless, in [0, 1)) ordered as DIMS;
            fixed_redshift: if given, replaces the redshift drawn from u (spectral fits); dimensionless z.
    Returns: (bagpipes_components, drawn). drawn holds the DRAWN_LABELS: redshift (dimensionless), age and
             tau [Gyr], metallicity [Z/Z_sun], av_ism [mag], tau_bc and logU (dimensionless),
             log_mass_formed [log M_sun].
    """
    params = dict(zip(DIMS, u))

    # 1. Redshift, which sets the oldest allowed age
    if fixed_redshift is None:
        redshift = uniform(params["z"], *Z_RANGE)
    else:
        redshift = fixed_redshift
    oldest_age = maximum_age(redshift)

    # 2. Star formation history and metallicity
    age = 10 ** uniform(params["log_age_fraction"], np.log10(AGE_MIN), np.log10(oldest_age))
    # Rounding in 10**log10 can push the top draw a hair above the bound
    age = min(age, oldest_age)
    tau = 10 ** uniform(params["log_tau"], *np.log10(TAU_RANGE))
    metallicity = 10 ** uniform(params["log_metallicity"], *np.log10(METALLICITY_RANGE))

    # 3. Dust and nebular emission
    av_ism = uniform(params["av_ism"], *AV_ISM_RANGE)
    tau_birth_cloud = 10 ** uniform(params["log_tau_bc"], *LOG_TAU_BC_RANGE)
    av_birth_cloud = TAU_TO_AV * tau_birth_cloud
    log_u = uniform(params["logU"], *LOGU_RANGE)

    bagpipes_components = {
        "redshift": redshift,
        "t_bc": BIRTH_CLOUD_AGE,
        "dust": {"type": "CF00", "Av": av_ism, "eta": 1 + av_birth_cloud / av_ism, "n": DUST_SLOPE},
        "nebular": {"logU": log_u},
        "delayed": {"age": age, "tau": tau, "massformed": LOG_MASS_FORMED, "metallicity": metallicity},
    }
    drawn = {"redshift": redshift, "age": age, "tau": tau, "metallicity": metallicity, "av_ism": av_ism,
             "tau_bc": tau_birth_cloud, "logU": log_u, "log_mass_formed": LOG_MASS_FORMED}
    return bagpipes_components, drawn


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

# Label columns of the parquet sets and the fits: the drawn parameters, then the derived labels
DRAWN_LABELS = ["redshift", "age", "tau", "metallicity", "av_ism", "tau_bc", "logU", "log_mass_formed"]
DERIVED_LABELS = ["log_age", "log_tau", "log_metallicity", "log_mstar", "log_sfr", "log_ssfr",
                  "mw_age", "f_young", "f_mid", "f_old"]
TRAINING_LABELS = DRAWN_LABELS + DERIVED_LABELS

# Labels that scale with brightness: log value - 0.4 * (m_obs - m_model)
MASS_DEPENDENT_LABELS = ["log_mass_formed", "log_mstar", "log_sfr"]


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


def delayed_mass_fraction_younger_than(age, tau, lookback):
    """Share of the formed mass with lookback time below `lookback` for SFR(t) ~ t exp(-t / tau).

    t runs from the start of star formation, `age` ago. With P and Q the regularised lower and upper
    incomplete gamma functions, the mass formed up to t is proportional to P(2, t / tau).
    Why the Q form: Q(2, (A - L) / tau) - Q(2, A / tau) avoids the cancellation of 1 - P / P when tau << A.
    Inputs: age A, tau, lookback L [Gyr], scalars or arrays that broadcast together.
    Returns: fraction in [0, 1] (dimensionless), an array of the broadcast shape; 1 where L >= A.
    """
    age, tau, lookback = np.broadcast_arrays(np.asarray(age, dtype=float), np.asarray(tau, dtype=float),
                                             np.asarray(lookback, dtype=float))
    formed_before_now = gammainc(2, age / tau)
    older_start = np.clip(age - lookback, 0, None) / tau
    fraction = (gammaincc(2, older_start) - gammaincc(2, age / tau)) / formed_before_now
    return np.where(lookback >= age, 1.0, np.clip(fraction, 0.0, 1.0))


def delayed_mass_weighted_age(age, tau):
    """Formed-mass-weighted lookback age of SFR(t) ~ t exp(-t / tau): A - 2 tau P(3, A / tau) / P(2, A / tau).

    Limits: A - 2 tau when tau << A, and A / 3 when tau >> A.
    Inputs: age A and tau [Gyr], scalars or arrays that broadcast together.
    Returns: mass-weighted age [Gyr ago], an array of the broadcast shape.
    """
    age = np.asarray(age, dtype=float)
    tau = np.asarray(tau, dtype=float)
    return age - 2 * tau * gammainc(3, age / tau) / gammainc(2, age / tau)


def delayed_tau_labels(age, tau, log_mass_formed, log_mstar):
    """Analytic labels of a delayed-tau SFH, shared by the training sets and the ML decoding.

    log_sfr is the mean SFR over the last SFR_TIMESCALE, floored at SFR_FLOOR * M_star: about a fifth
    of the prior has no star formation that recent, and log(0) would be -inf.
    Inputs: age, tau [Gyr]; log_mass_formed, log_mstar [log M_sun]; arrays of one shape (or scalars).
    Returns: dict of arrays: mw_age [Gyr ago], f_young, f_mid, f_old (dimensionless, sum to 1; split at
             LABEL_EDGES), log_sfr [log M_sun/yr] and log_ssfr [log yr^-1].
    """
    log_mass_formed = np.asarray(log_mass_formed, dtype=float)
    log_mstar = np.asarray(log_mstar, dtype=float)
    younger_than_young_edge = delayed_mass_fraction_younger_than(age, tau, LABEL_EDGES[0])
    younger_than_old_edge = delayed_mass_fraction_younger_than(age, tau, LABEL_EDGES[1])
    recent_fraction = delayed_mass_fraction_younger_than(age, tau, SFR_TIMESCALE)

    recent_sfr = 10**log_mass_formed * recent_fraction / (SFR_TIMESCALE * GYR_TO_YEARS)
    sfr_floor = SFR_FLOOR * 10**log_mstar
    log_sfr = np.log10(np.maximum(recent_sfr, sfr_floor))
    return {
        "mw_age": delayed_mass_weighted_age(age, tau),
        "f_young": younger_than_young_edge,
        "f_mid": younger_than_old_edge - younger_than_young_edge,
        "f_old": 1 - younger_than_old_edge,
        "log_sfr": log_sfr,
        "log_ssfr": log_sfr - log_mstar,
    }


def summary_labels(bagpipes_sfh, drawn):
    """Every TRAINING_LABELS value of one model.

    Inputs: bagpipes_sfh: the model's sfh object (stellar_mass in log M_sun, living stars);
            drawn: dict from components() with the DRAWN_LABELS (extra keys are ignored).
    Returns: dict ordered as TRAINING_LABELS: the drawn values; log_age and log_tau [log Gyr] and
             log_metallicity [log Z/Z_sun]; log_mstar from Bagpipes; the delayed_tau_labels values as floats.
    """
    labels = {}
    for label in DRAWN_LABELS:
        labels[label] = drawn[label]
    labels["log_age"] = np.log10(drawn["age"])
    labels["log_tau"] = np.log10(drawn["tau"])
    labels["log_metallicity"] = np.log10(drawn["metallicity"])
    labels["log_mstar"] = bagpipes_sfh.stellar_mass

    analytic_labels = delayed_tau_labels(drawn["age"], drawn["tau"], drawn["log_mass_formed"],
                                         bagpipes_sfh.stellar_mass)
    for label in ("log_sfr", "log_ssfr", "mw_age", "f_young", "f_mid", "f_old"):
        labels[label] = float(analytic_labels[label])
    return labels


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
    """Labels and AB magnitudes (10^10 M_sun formed) for a block of unit-cube vectors.

    Inputs: unit_vectors: iterable of unit-cube vectors (dimensionless, in [0, 1)), each ordered as DIMS.
    Returns: DataFrame with one row per model: the TRAINING_LABELS, then the AB magnitudes [mag] in BANDS
             order (NaN where a flux is not positive). The "z" column is the z-band magnitude, not the redshift.

    One Bagpipes model is built and then updated: update() is ~70x faster than building a new one, and safe
    because every model has the same component set (Bagpipes fixes the components when it builds a model).
    """
    model = None
    rows = []
    for u in unit_vectors:
        bagpipes_components, drawn = components(u)
        if model is None:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
        else:
            model.update(bagpipes_components)

        row = summary_labels(model.sfh, drawn)
        magnitudes = ab_magnitudes(model.photometry)
        for band, magnitude in zip(BANDS, magnitudes):
            row[band] = magnitude
        rows.append(row)
    return pd.DataFrame(rows, columns=TRAINING_LABELS + BANDS)


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
    Vectorised: it draws in the same order as the original per-model loop, so the same rng gives the same errors.
    """
    table_bins = noise_bin_index(band_table["magnitude_low"].to_numpy() + NOISE_BIN_WIDTH / 2)
    quantile_table = band_table[quantile_column_names()].to_numpy()

    log_errors = np.full(len(magnitudes), np.nan)
    finite = np.isfinite(magnitudes)
    model_bins = np.full(len(magnitudes), -1)
    model_bins[finite] = noise_bin_index(magnitudes[finite])

    # Models fainter than the faintest populated bin get no error (undetected)
    drawable = finite & (model_bins <= table_bins[-1])
    drawable_indices = np.flatnonzero(drawable)

    # Nearest populated bin at or below each model's bin; the brightest bin covers anything brighter
    row_indices = np.searchsorted(table_bins, model_bins[drawable_indices], side="right") - 1
    row_indices = np.maximum(row_indices, 0)

    # One uniform draw per drawable model, in the same order as the original loop
    quantile_positions = rng.random(len(drawable_indices)) * (NOISE_N_QUANTILES - 1)
    lower = np.floor(quantile_positions).astype(int)
    upper = np.minimum(lower + 1, NOISE_N_QUANTILES - 1)
    weights = quantile_positions - lower
    lower_values = quantile_table[row_indices, lower]
    upper_values = quantile_table[row_indices, upper]
    log_errors[drawable_indices] = lower_values + weights * (upper_values - lower_values)
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

    # Window of real galaxies around each model's redshift
    first = np.searchsorted(sorted_z, model_z - max_dz)
    last = np.searchsorted(sorted_z, model_z + max_dz)
    has_neighbours = last > first

    offsets = np.full(len(model_magnitudes), np.nan)
    chosen = rng.integers(first[has_neighbours], last[has_neighbours])
    model_r = model_magnitudes[has_neighbours, r_index]
    offsets[has_neighbours] = sorted_r[chosen] - model_r

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
# 5. Spectra, spectral likelihood and split
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


def training_unit_vector(cube_by_name):
    """Unit-cube vector in DIMS order for a fit at fixed redshift, from the fit's own named coordinates.

    Inputs: cube_by_name: dict of unit-cube values (dimensionless, in [0, 1)) holding every DIMS name but "z";
            other names (the fit's nuisance parameters) are ignored.
    Returns: array of len(DIMS) values (dimensionless); the "z" slot holds IGNORED_REDSHIFT_DRAW.
    """
    unit_vector = np.zeros(len(DIMS))
    for dimension_index, dimension in enumerate(DIMS):
        if dimension == "z":
            # components() ignores this value, because fixed_redshift replaces the drawn redshift
            unit_vector[dimension_index] = IGNORED_REDSHIFT_DRAW
        else:
            unit_vector[dimension_index] = cube_by_name[dimension]
    return unit_vector


def spectral_components(cube, redshift):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit unit-cube point.

    The delayed-tau SFH, the single metallicity, dust and nebular part is the training prior (components) at
    the fixed spectroscopic redshift, which also sets the oldest allowed age; veldisp and noise_scale are the
    spectral nuisance parameters.
    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as SPECTRAL_DIMS;
            redshift: spectroscopic redshift (dimensionless).
    Returns: (bagpipes_components, drawn, noise_scale). The components are as in components() plus
             veldisp [km/s]; drawn is as in components() plus veldisp [km/s]; noise_scale is the error
             multiplier (dimensionless, 1 to 5).
    """
    cube_by_name = dict(zip(SPECTRAL_DIMS, cube))
    unit_vector = training_unit_vector(cube_by_name)
    bagpipes_components, drawn = components(unit_vector, fixed_redshift=redshift)
    veldisp = uniform(cube_by_name["veldisp"], *VELDISP_RANGE)
    bagpipes_components["veldisp"] = veldisp
    drawn["veldisp"] = veldisp
    # Uniform in log10 so that small and large error multipliers are equally likely per decade
    noise_scale = 10 ** uniform(cube_by_name["noise_scale"], *LOG_NOISE_SCALE_RANGE)
    return bagpipes_components, drawn, noise_scale


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
OFFSET_DRAW_SEED = 0             # seed of the per-sample magnitude offsets of a photometric fit's labels
SPECTRUM_JOB = "spectrum"        # job kind of a spectral fit; also every job or fit file stored without a kind
PHOTOMETRY_JOB = "photometry"    # job kind of a photometric fit (spec amendment A)
# SFH model of every fit this code writes, stored in each fit file: fits made under an earlier model
# (the 16-bin non-parametric SFH, whose files have no sfh_model) must never be read as current ones
SFH_MODEL = "delayed"

# The photometric fit (spec amendment A) samples the training prior at the spectroscopic redshift and nothing
# else: its brightness offset is marginalised analytically, and the simulator has no error multiplier
PHOTOMETRIC_DIMS = []
for dimension in DIMS:
    if dimension != "z":
        PHOTOMETRIC_DIMS.append(dimension)


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

    # 3. MultiNest over the unit cube
    return run_multinest(log_likelihood, len(SPECTRAL_DIMS), output_basename, max_iterations, n_live_points)


def run_multinest(log_likelihood, n_dims, output_basename, max_iterations, n_live_points):
    """One MultiNest run over the unit cube, with the settings shared by the spectral and photometric fits.

    The prior is the identity map: each fit's components function turns a unit-cube point into physical
    values, exactly as the training set does.
    Inputs: log_likelihood: function of one unit-cube point (dimensionless) giving a log-likelihood
            (dimensionless); it must never raise; n_dims: number of unit-cube dimensions (count);
            output_basename: directory prefix for the MultiNest files (str, under 100 characters);
            max_iterations: MultiNest iteration limit (count); n_live_points: live points (count).
    Returns: dict as fit_spectrum returns: samples (n_samples, n_dims), log_evidence, n_iterations, converged
             and run_time [s].
    """
    # pymultinest loads the MultiNest C library on import, so it needs LD_LIBRARY_PATH (set by the kernel)
    import pymultinest

    os.makedirs(output_basename, exist_ok=True)

    def unit_cube_prior(cube):
        """Identity map: the training prior is already uniform on the unit cube.

        Inputs: cube: unit-cube point (dimensionless). Returns: the same point (dimensionless).
        """
        return cube

    # 1. Run MultiNest, timing only the sampling itself
    start = time.time()
    result = pymultinest.solve(
        LogLikelihood=log_likelihood, Prior=unit_cube_prior, n_dims=n_dims,
        outputfiles_basename=output_basename, n_live_points=n_live_points,
        sampling_efficiency=SAMPLING_EFFICIENCY, multimodal=True, max_iter=max_iterations,
        # Standard nested-sampling evidence: the importance-sampled one underflows in short runs and
        # MultiNest then writes it in a form pymultinest cannot parse
        importance_nested_sampling=False,
        resume=False, verbose=False, use_MPI=False,
    )
    run_time = time.time() - start

    # 2. Hitting the iteration limit means the run stopped before it converged; the upper bound of the
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
    Returns: DataFrame with one row per kept sample: the TRAINING_LABELS (mass-like ones for 10^10 M_sun
             formed), then veldisp [km/s], noise_scale (dimensionless) and model_r [mag].
    """
    # 1. Thin the samples evenly: neighbouring MultiNest samples are nearly identical
    kept_samples = thinned_samples(unit_cube_samples)

    # 2. Labels per kept sample, with one photometry model updated in place
    model = None
    rows = []
    for cube in kept_samples:
        bagpipes_components, drawn, noise_scale = spectral_components(cube, redshift)
        model = updated_photometry_model(model, bagpipes_components)
        drawn_with_noise = dict(drawn)
        drawn_with_noise["noise_scale"] = noise_scale
        rows.append(sample_label_row(model, bagpipes_components, drawn_with_noise))
    return pd.DataFrame(rows)


def thinned_samples(unit_cube_samples):
    """At most MAX_POSTERIOR_SAMPLES evenly spaced posterior samples, because building labels is the slow part.

    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, n_dims).
    Returns: every step-th row (dimensionless), step chosen so that no more than MAX_POSTERIOR_SAMPLES remain.
    """
    step = max(1, int(np.ceil(len(unit_cube_samples) / MAX_POSTERIOR_SAMPLES)))
    return unit_cube_samples[::step]


def updated_photometry_model(model, bagpipes_components):
    """A Bagpipes model with the S-PLUS filters at new components: built the first time, then updated in place.

    Inputs: model: an earlier photometry model_galaxy, or None; bagpipes_components: dict as from components().
    Returns: the model_galaxy (filt_list=FILTER_FILES), now at bagpipes_components.
    """
    if model is None:
        return pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
    model.update(bagpipes_components)
    return model


def sample_label_row(model, bagpipes_components, drawn):
    """Labels of one posterior sample, read from its photometry model as the training set labels a model.

    Inputs: model: photometry model_galaxy already at bagpipes_components; bagpipes_components: dict as from
            components(), the point the model is at (not read: the labels come from drawn and the model);
            drawn: dict from components() with the DRAWN_LABELS, plus any extra parameters of the fit
            (veldisp [km/s], noise_scale).
    Returns: dict: the TRAINING_LABELS in order (summary_labels), then the extra drawn parameters, then model_r
             [mag], the model r magnitude for 10^10 M_sun formed. Mass-like labels are for that mass.
    """
    row = summary_labels(model.sfh, drawn)
    for key, value in drawn.items():
        if key not in row:
            row[key] = value
    row["model_r"] = ab_magnitudes(model.photometry)[BANDS.index("r")]
    return row


# -----------------------------------------------------------------------------
# Photometric fit (spec amendment A): the ML's own input, fitted by Bagpipes with the simulator's likelihood
# -----------------------------------------------------------------------------
def photometric_components(cube, redshift):
    """Bagpipes components and drawn parameters for one photometric-fit unit-cube point.

    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as PHOTOMETRIC_DIMS;
            redshift: spectroscopic redshift (dimensionless).
    Returns: (bagpipes_components, drawn), exactly as components(..., fixed_redshift=redshift) gives them.
    """
    cube_by_name = dict(zip(PHOTOMETRIC_DIMS, cube))
    unit_vector = training_unit_vector(cube_by_name)
    return components(unit_vector, fixed_redshift=redshift)


def offset_marginal_log_likelihood(model_magnitudes, observed_magnitudes, errors, detected):
    """Gaussian log-likelihood of the photometry in magnitudes, integrated over a flat-prior magnitude offset.

    The offset (observed minus model, the same in every band) is the galaxy's brightness against a model of
    10^10 M_sun formed. With weights w = 1 / error^2 and d = m_obs - m_model over the detected bands, the
    likelihood is Gaussian in the offset, so the integral is exact:
        ln L = -chi2_min / 2 + ln(2 pi / sum w) / 2,   chi2_min = sum w d^2 - (sum w d)^2 / sum w.
    No ERROR_FLOOR is added, because add_splus_noise, which made the training data, adds none. The Gaussian
    normalisation of each band does not depend on the model, so it is dropped.
    Inputs: model_magnitudes, observed_magnitudes, errors: [mag], 1-D over the bands; detected: bool array of
            the same length with at least one True (detected_bands); undetected entries may hold anything, NaN too.
    Returns: (log_likelihood (dimensionless), best_offset [mag], offset_sigma [mag]); given the model, the
             offset is Gaussian with mean best_offset and standard deviation offset_sigma.
    """
    model_magnitudes = np.asarray(model_magnitudes, dtype=float)
    observed_magnitudes = np.asarray(observed_magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)

    differences = observed_magnitudes[detected] - model_magnitudes[detected]
    weights = 1 / errors[detected] ** 2
    weight_sum = np.sum(weights)
    weighted_difference_sum = np.sum(weights * differences)

    best_offset = weighted_difference_sum / weight_sum
    chi2_min = np.sum(weights * differences**2) - weighted_difference_sum * best_offset
    log_likelihood = -0.5 * chi2_min + 0.5 * np.log(2 * np.pi / weight_sum)
    offset_sigma = 1 / np.sqrt(weight_sum)
    return log_likelihood, best_offset, offset_sigma


def safe_photometric_log_likelihood(model_magnitudes, observed_magnitudes, errors, detected):
    """offset_marginal_log_likelihood that returns LOG_ZERO for a model with no magnitude in a detected band.

    ab_magnitudes gives NaN when a model flux is not positive; such a model cannot have made the data.
    Inputs: as in offset_marginal_log_likelihood.
    Returns: log-likelihood (dimensionless), or LOG_ZERO for an impossible model.
    """
    model_magnitudes = np.asarray(model_magnitudes, dtype=float)
    if not np.all(np.isfinite(model_magnitudes[detected])):
        return LOG_ZERO
    return offset_marginal_log_likelihood(model_magnitudes, observed_magnitudes, errors, detected)[0]


def fit_photometry(magnitudes, errors, redshift, output_basename, max_iterations, n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one galaxy's 12-band photometry under the training prior at a fixed redshift.

    It fits the ML's own target: the training prior (PHOTOMETRIC_DIMS) and the likelihood of the simulator
    that made the training data (offset_marginal_log_likelihood), over the detected bands only.
    Inputs: magnitudes, errors: [mag], 12 values in BANDS order; a band is left out when detected_bands says
            it is not detected (99 or NaN magnitude, or error >= MAX_DETECTED_ERROR);
            redshift: spectroscopic redshift (dimensionless); output_basename, max_iterations and n_live_points:
            as in fit_spectrum.
    Returns: dict as fit_spectrum returns, with samples of shape (n_samples, len(PHOTOMETRIC_DIMS)).
    Raises: ValueError when fewer than MIN_DETECTED_BANDS bands are detected.
    """
    # 1. Bands that enter the likelihood
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)
    if detected.sum() < MIN_DETECTED_BANDS:
        raise ValueError(f"{detected.sum()} detected bands, fewer than MIN_DETECTED_BANDS = {MIN_DETECTED_BANDS}")

    # 2. Bagpipes model built once with the S-PLUS filters; every point has the same component set
    midpoint_cube = np.full(len(PHOTOMETRIC_DIMS), MIDPOINT_DRAW)
    starting_components, _ = photometric_components(midpoint_cube, redshift)
    model = updated_photometry_model(None, starting_components)

    # 3. Likelihood of one unit-cube point
    def log_likelihood(cube):
        """Log-likelihood of the photometry for one unit-cube point.

        Inputs: cube: unit-cube point (dimensionless, PHOTOMETRIC_DIMS order).
        Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
        """
        bagpipes_components, _ = photometric_components(cube, redshift)
        # PyMultiNest calls sys.exit if the likelihood raises: a model Bagpipes cannot compute is impossible
        try:
            model.update(bagpipes_components)
            model_magnitudes = ab_magnitudes(model.photometry)
        except Exception:
            return LOG_ZERO
        return safe_photometric_log_likelihood(model_magnitudes, magnitudes, errors, detected)

    # 4. MultiNest over the unit cube
    return run_multinest(log_likelihood, len(PHOTOMETRIC_DIMS), output_basename, max_iterations, n_live_points)


def photometric_posterior_labels(unit_cube_samples, redshift, magnitudes, errors, rng):
    """Labels of each photometric-fit posterior sample, at the brightness of the observed galaxy.

    Each kept sample gets its own magnitude offset, drawn from the offset's Gaussian given that sample's
    model (offset_marginal_log_likelihood), and its mass-like labels move by it (rescale_mass_labels), as the
    training labels move with the brightness offset of each model.
    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, len(PHOTOMETRIC_DIMS));
            redshift: spectroscopic redshift (dimensionless); magnitudes, errors: [mag] the fitted photometry,
            as given to fit_photometry; rng: numpy Generator.
    Returns: DataFrame with one row per kept sample (at most MAX_POSTERIOR_SAMPLES): the TRAINING_LABELS
             (MASS_DEPENDENT_LABELS already rescaled), model_r [mag, for 10^10 M_sun formed] and
             magnitude_offset [mag, observed minus model].
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)

    # 1. Labels per kept sample, each with an offset drawn given its own model
    kept_samples = thinned_samples(unit_cube_samples)
    model = None
    rows = []
    for cube in kept_samples:
        bagpipes_components, drawn = photometric_components(cube, redshift)
        model = updated_photometry_model(model, bagpipes_components)
        row = sample_label_row(model, bagpipes_components, drawn)
        model_magnitudes = ab_magnitudes(model.photometry)
        _, best_offset, offset_sigma = offset_marginal_log_likelihood(model_magnitudes, magnitudes, errors, detected)
        row["magnitude_offset"] = rng.normal(best_offset, offset_sigma)
        rows.append(row)

    # 2. Mass-like labels at the observed brightness
    labels = pd.DataFrame(rows)
    return rescale_mass_labels(labels, labels["magnitude_offset"].to_numpy())


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
    or failed write never leaves a partial file at path (which would count as a finished job). Every file
    also records sfh_model = SFH_MODEL, so fits of an earlier SFH model can be told apart (fit_is_current).
    Inputs: path: output file path (str); posterior: DataFrame of labels, one column per label;
            unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, n_dims), n_dims that of the
            fit's kind: len(SPECTRAL_DIMS) or len(PHOTOMETRIC_DIMS));
            metadata: dict of scalars (converged, run_time [s], n_iterations, max_iterations, log_evidence, redshift,
            kind: SPECTRUM_JOB or PHOTOMETRY_JOB; files written before kinds existed have none and are spectral).
            The caller's dict is not changed.
    Returns: None. An existing file at path is replaced.
    """
    stored_metadata = dict(metadata)
    stored_metadata["sfh_model"] = SFH_MODEL
    temporary_path = path + ".tmp"
    with h5py.File(temporary_path, "w") as output:
        output.create_dataset("unit_cube_samples", data=unit_cube_samples)
        labels_group = output.create_group("posterior")
        for column in posterior.columns:
            labels_group.create_dataset(column, data=posterior[column].to_numpy(dtype=float))
        for key, value in stored_metadata.items():
            output.attrs[key] = value
    os.replace(temporary_path, path)


def load_fit(path):
    """Read a file written by save_fit.

    Inputs: path: file path (str).
    Returns: (posterior DataFrame with one column per label (in alphabetical order, as HDF5 stores them),
             metadata dict with the attributes save_fit wrote, sfh_model included; files written before the
             delayed-tau model have none).
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


def fit_is_current(metadata):
    """Whether a stored fit belongs to the current production run: SFH_MODEL and the full iteration limit.

    Fits that fail this (16-bin fits without sfh_model, pilot runs with a smaller limit) are moved to
    superseded/ by supersede_outdated_fits (spectroscopic_reference.ipynb, section 4) and never read as results.
    Inputs: metadata: dict as load_fit returns it.
    Returns: bool: True only when sfh_model == SFH_MODEL and max_iterations >= MAX_ITERATIONS.
    """
    same_model = metadata.get("sfh_model") == SFH_MODEL
    full_run = metadata.get("max_iterations", 0) >= MAX_ITERATIONS
    return same_model and full_run


def job_kind(job):
    """Kind of a fit job: SPECTRUM_JOB or PHOTOMETRY_JOB.

    Inputs: job: job dict (as fit_and_save takes).
    Returns: job["kind"], or SPECTRUM_JOB for a job without one (jobs made before the photometric fits existed).
    """
    return job.get("kind", SPECTRUM_JOB)


def fit_spectrum_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a spectrum job and label its posterior samples.

    Inputs: job: dict with spectrum (as from prepare_spectrum) and redshift (dimensionless); multinest_basename:
            directory prefix for the MultiNest files (str); max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_spectrum, posterior DataFrame from posterior_labels).
    """
    fit = fit_spectrum(job["spectrum"], job["redshift"], multinest_basename, max_iterations, n_live_points)
    posterior = posterior_labels(fit["samples"], job["redshift"])
    return fit, posterior


def fit_photometry_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a photometry job and label its posterior samples at the galaxy's brightness.

    Inputs: job: dict with magnitudes and errors ([mag], 12 in BANDS order) and redshift (dimensionless);
            multinest_basename: directory prefix for the MultiNest files (str); max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_photometry, posterior DataFrame from photometric_posterior_labels).
    """
    fit = fit_photometry(job["magnitudes"], job["errors"], job["redshift"], multinest_basename,
                         max_iterations, n_live_points)
    # A fixed seed, so labelling the same samples again gives the same offsets
    rng = np.random.default_rng(OFFSET_DRAW_SEED)
    posterior = photometric_posterior_labels(fit["samples"], job["redshift"], job["magnitudes"], job["errors"], rng)
    return fit, posterior


# The fitter of each job kind: same signature, (job, multinest_basename, max_iterations, n_live_points)
JOB_FITTERS = {SPECTRUM_JOB: fit_spectrum_job, PHOTOMETRY_JOB: fit_photometry_job}


def fit_and_save(job, output_directory, max_iterations, n_live_points=N_LIVE_POINTS):
    """Fit one job and save <output_directory>/<job_id>.h5. Never raises.

    Skips jobs whose readable file of the same kind and SFH_MODEL already exists, so a batch can be re-run with
    more galaxies. An unreadable file does not count as finished: the job is fitted again and the file replaced.
    A readable file of the other kind is never replaced (job ids are specobjids, shared by the spectral and
    photometric fits of one galaxy), and neither is a fit of an earlier SFH model (the only copy of that result,
    to be moved to superseded/): the job fails instead. A failed fit comes back with status "failed" and the
    error text instead of an exception.
    Inputs: job: dict with job_id (str), redshift (dimensionless) and kind (SPECTRUM_JOB when absent); a
            spectrum job holds spectrum (dict as from prepare_spectrum), a photometry job magnitudes and
            errors ([mag], 12 in BANDS order); output_directory: folder for the .h5 file and the MultiNest
            files (str); max_iterations: MultiNest iteration limit (count); n_live_points: live points (count).
    Returns: dict with job_id, status ("fitted", "skipped" or "failed"), converged (bool or None),
             run_time [s] (or None), n_iterations (count or None) and error (str or None).
    """
    # 1. Output path and a summary that the steps below fill in
    output_path = os.path.join(output_directory, f"{job['job_id']}.h5")
    kind = job_kind(job)
    summary = {"job_id": job["job_id"], "status": None, "converged": None,
               "run_time": None, "n_iterations": None, "error": None}

    # 2. Skip finished jobs; their run details come from the file. Any failure to read it (damaged,
    # truncated or missing fields) means the job is not finished, so it is fitted again below. A readable fit of
    # an earlier SFH model is neither skipped (it is not a result of this model) nor overwritten (its only copy)
    stored_kind = None
    stored_model = None
    if os.path.exists(output_path):
        try:
            _, metadata = load_fit(output_path)
            stored_details = {"converged": metadata["converged"], "run_time": metadata["run_time"],
                              "n_iterations": metadata["n_iterations"]}
            stored_kind = metadata.get("kind", SPECTRUM_JOB)
            stored_model = metadata.get("sfh_model")
        except Exception:
            pass
    if stored_kind == kind and stored_model != SFH_MODEL:
        summary.update(status="failed", error=f"{output_path} holds a fit of an earlier SFH model "
                                              f"({stored_model}, not {SFH_MODEL}): move it to superseded/")
        return summary
    if stored_kind == kind:
        summary.update(status="skipped", **stored_details)
        return summary
    if stored_kind is not None:
        summary.update(status="failed", error=f"{output_path} holds a {stored_kind} fit and this job is a "
                                              f"{kind} fit: use a separate output folder for each kind")
        return summary

    # 3. Fit, label the posterior and save; any error is recorded rather than raised
    try:
        multinest_basename = os.path.join(output_directory, str(job["job_id"])) + "/"
        fit, posterior = JOB_FITTERS[kind](job, multinest_basename, max_iterations, n_live_points)
        metadata = {"converged": fit["converged"], "run_time": fit["run_time"],
                    "n_iterations": fit["n_iterations"], "max_iterations": max_iterations,
                    "log_evidence": fit["log_evidence"], "redshift": job["redshift"], "kind": kind}
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
# Labels the calibration test checks; age, tau and Z in log, because the prior spans several decades of each
CALIBRATION_LABELS = ["log_age", "log_tau", "log_metallicity", "log_mstar"]
CALIBRATION_LEVELS = (0.68, 0.95)                   # central-interval probabilities reported (fractions)
CALIBRATION_ID_PREFIX = "calib_"                    # job ids are calib_000, calib_001, ...
PERCENT = 100                    # percent per unit fraction


def export_job(job, jobs_directory):
    """Write one fit job to <jobs_directory>/<job_id>.npz, a portable file for run_spectral_fits.py.

    Inputs: job: dict with job_id (str), redshift (dimensionless) and kind (SPECTRUM_JOB when absent); a
            spectrum job holds spectrum (dict with wavelength [Angstrom, observed], flux and ivar), a photometry
            job magnitudes and errors ([mag], 12 in BANDS order); jobs_directory: folder for the file (str),
            created if missing.
    Returns: path of the written file (str).
    """
    os.makedirs(jobs_directory, exist_ok=True)
    path = os.path.join(jobs_directory, f"{job['job_id']}.npz")
    kind = job_kind(job)
    arrays = {"job_id": str(job["job_id"]), "redshift": job["redshift"], "kind": kind}
    if kind == PHOTOMETRY_JOB:
        arrays["magnitudes"] = job["magnitudes"]
        arrays["errors"] = job["errors"]
    else:
        arrays["wavelength"] = job["spectrum"]["wavelength"]
        arrays["flux"] = job["spectrum"]["flux"]
        arrays["ivar"] = job["spectrum"]["ivar"]
    np.savez(path, **arrays)
    return path


def load_job(path):
    """Read a file written by export_job.

    Inputs: path: file path (str).
    Returns: job dict in the format fit_and_save takes: job_id (str), kind, redshift (dimensionless float) and
             either spectrum (wavelength, flux, ivar) or magnitudes and errors [mag]. A file without a kind
             (exported before the photometric fits existed) is a spectrum job.
    """
    with np.load(path) as source:
        kind = SPECTRUM_JOB
        if "kind" in source.files:
            kind = str(source["kind"])
        job = {"job_id": str(source["job_id"]), "kind": kind, "redshift": float(source["redshift"])}
        if kind == PHOTOMETRY_JOB:
            job["magnitudes"] = source["magnitudes"]
            job["errors"] = source["errors"]
        else:
            job["spectrum"] = {"wavelength": source["wavelength"], "flux": source["flux"], "ivar": source["ivar"]}
        return job


def draw_calibration_cubes(n_mocks, rng):
    """Unit-cube points of the calibration mocks: uniform over every SPECTRAL_DIMS coordinate.

    Uniform on the unit cube is the training prior at the template's redshift (log-uniform age, tau and Z).
    draw_calibration_mocks calls this first on its rng, so the cubes can be redrawn from the same seed.
    Inputs: n_mocks: number of points (count); rng: numpy Generator.
    Returns: array of shape (n_mocks, len(SPECTRAL_DIMS)), values in [0, 1) (dimensionless).
    """
    return rng.random((n_mocks, len(SPECTRAL_DIMS)))


def draw_calibration_mocks(n_mocks, templates, rng):
    """Mock spectra with known labels, each on the wavelengths, mask and errors of a random real galaxy.

    Inputs: n_mocks: number of mocks (count); templates: list of (spectrum dict from prepare_spectrum,
            redshift (dimensionless)); rng: numpy Generator.
    Returns: (jobs, truths): jobs is a list of dicts as fit_and_save takes (job_id calib_000, calib_001, ...);
             truths maps each job_id to its true labels (Series, as posterior_labels columns: the
             TRAINING_LABELS, veldisp, noise_scale and model_r).
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


# =============================================================================
# 8. Superseding the results of an earlier model or run (spectroscopic_reference.ipynb, section 4)
# =============================================================================
SUPERSEDED_SUBDIRECTORY = "superseded"     # inside each fits folder: results of earlier, replaced runs
SUPERSEDED_NAME_SUFFIX = "_16bin"          # added to a name already taken there (first-run fits of the same id)
CALIBRATION_MOCK_FILE = "calibration_mocks.joblib"   # (jobs, truths) of the calibration mocks of a fits folder
JOBS_SUBDIRECTORY = "jobs"                 # exported fit jobs, one .npz each
RUN_LOG_NAMES = ("fit_log.csv", "fit_logs")   # logs of run_spectral_fits.py; their rows name the fits that ran


def move_to_superseded(fits_directory, name):
    """Move a file or folder of fits_directory into its SUPERSEDED_SUBDIRECTORY, keeping its relative path.

    A name already taken there gets SUPERSEDED_NAME_SUFFIX before its extension (then _2, _3, ... if that is
    taken too), so an earlier superseded copy is never overwritten and nothing is left in place.
    Inputs: fits_directory: folder of the fits (str); name: path relative to it (str), e.g. "123.h5",
            "123" (a MultiNest folder) or "jobs/calib_000.npz".
    Returns: the destination path (str), or None when the source does not exist.
    """
    source = os.path.join(fits_directory, name)
    if not os.path.exists(source):
        return None
    superseded_directory = os.path.join(fits_directory, SUPERSEDED_SUBDIRECTORY)
    stem, extension = os.path.splitext(name)
    destination = os.path.join(superseded_directory, name)
    copy_number = 1
    while os.path.exists(destination):
        counter = "" if copy_number == 1 else f"_{copy_number}"
        destination = os.path.join(superseded_directory, f"{stem}{SUPERSEDED_NAME_SUFFIX}{counter}{extension}")
        copy_number += 1
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    shutil.move(source, destination)
    return destination


def calibration_mocks_are_current(mock_file):
    """Whether a stored calibration-mock file holds truths with the current calibration labels.

    Inputs: mock_file: path of a joblib file with (jobs, truths), truths mapping job_id -> Series of true labels.
    Returns: bool: True when every truth has every CALIBRATION_LABELS value (16-bin truths lack log_age).
    """
    _, truths = joblib.load(mock_file)
    for truth in truths.values():
        if not set(CALIBRATION_LABELS) <= set(truth.index):
            return False
    return True


def supersede_outdated_fits(fits_directory):
    """Move every result of an earlier model or run out of a fits folder, into its SUPERSEDED_SUBDIRECTORY.

    Why: run_spectral_fits.py counts any top-level .h5 as done, and the new calibration mocks reuse the
    calib_<k> ids, so an old file left in place would be skipped as finished or scored against a new truth.
    Moved, in this order:
      1. every readable .h5 fit that fails fit_is_current, with its MultiNest folder (the name without .h5);
         an unreadable .h5 stays, because fit_and_save fits that job again;
      2. CALIBRATION_MOCK_FILE when its truths lack the CALIBRATION_LABELS (mocks of the 16-bin model);
      3. while no mock file is left, the jobs/calib_*.npz files, which then belong to no current mock;
      4. when any fit moved, the run logs RUN_LOG_NAMES (notebook 2 reads failures from them, and
         fetch_fits.sh drops new rows of ids already logged).
    Inputs: fits_directory: folder with the .h5 fits, jobs/ and the mock file (str); a missing folder is fine.
    Returns: list of the moved names (str), relative to fits_directory; empty when everything is current.
    """
    if not os.path.isdir(fits_directory):
        return []
    moved = []

    # 1. Fits of an earlier model or run, with their MultiNest folders
    for name in sorted(os.listdir(fits_directory)):
        if not name.endswith(".h5"):
            continue
        try:
            _, metadata = load_fit(os.path.join(fits_directory, name))
        except (OSError, KeyError):
            continue
        if fit_is_current(metadata):
            continue
        for fit_name in (name, name[:-len(".h5")]):
            if move_to_superseded(fits_directory, fit_name) is not None:
                moved.append(fit_name)
    fits_were_moved = len(moved) > 0

    # 2. Calibration mocks of an earlier model
    mock_file = os.path.join(fits_directory, CALIBRATION_MOCK_FILE)
    if os.path.exists(mock_file) and not calibration_mocks_are_current(mock_file):
        move_to_superseded(fits_directory, CALIBRATION_MOCK_FILE)
        moved.append(CALIBRATION_MOCK_FILE)

    # 3. Job files of mocks that no stored mock file holds (the notebook stores the mocks before exporting them)
    jobs_directory = os.path.join(fits_directory, JOBS_SUBDIRECTORY)
    if not os.path.exists(mock_file) and os.path.isdir(jobs_directory):
        for name in sorted(os.listdir(jobs_directory)):
            if name.startswith(CALIBRATION_ID_PREFIX):
                job_name = os.path.join(JOBS_SUBDIRECTORY, name)
                move_to_superseded(fits_directory, job_name)
                moved.append(job_name)

    # 4. Run logs of the moved fits
    if fits_were_moved:
        for name in RUN_LOG_NAMES:
            if move_to_superseded(fits_directory, name) is not None:
                moved.append(name)
    return moved
