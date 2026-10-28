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
from scipy.integrate import simpson
from scipy.stats import chi2 as chi_squared
from scipy.stats import qmc

import bagpipes as pipes
# The shared prior and the SFH models live in sfh_models.py. They are re-exported here, so sm.Z_RANGE,
# sm.maximum_age, sm.delayed_tau_labels and the rest keep working
from sfh_models import (AGE_MIN, AV_ISM_RANGE, BIRTH_CLOUD_AGE, COSMOLOGY, DUST_SLOPE, GYR_TO_YEARS,
                        HISTORY_LABELS, LABEL_EDGES, LOG_MASS_FORMED, LOG_TAU_BC_RANGE, LOGU_RANGE,
                        METALLICITY_RANGE, MODELS, SFR_FLOOR, SFR_TIMESCALE, SHARED_LABELS, TAU_RANGE,
                        TAU_TO_AV, UNIVERSE_AGE_MARGIN, Z_RANGE, delayed_mass_fraction_younger_than,
                        delayed_mass_weighted_age, delayed_tau_labels, maximum_age, model_file, uniform)

warnings.filterwarnings("ignore", category=RuntimeWarning)


# =============================================================================
# 1. Prior (fixed before looking at the data): see sfh_models.py
# =============================================================================


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


def summary_labels(bagpipes_sfh, drawn, sfh_model):
    """Every training label of one model of the given SFH model.

    Inputs: bagpipes_sfh: the model's sfh object (stellar_mass in log M_sun, living stars);
            drawn: dict from the model's components() with its drawn labels (extra keys are ignored);
            sfh_model: a name in MODELS ("delayed" or "dblplaw").
    Returns: dict ordered as MODELS[sfh_model].training_labels: the drawn values; the log labels (base-10 logs
             of drawn values); log_mstar from Bagpipes; the HISTORY_LABELS from the analytic_labels, as floats.
    """
    model_definition = MODELS[sfh_model]
    labels = {}
    for label in model_definition.drawn_labels:
        labels[label] = drawn[label]
    for log_label, drawn_label in model_definition.log_labels.items():
        labels[log_label] = np.log10(drawn[drawn_label])
    labels["log_mstar"] = bagpipes_sfh.stellar_mass

    # analytic_labels is vectorised: give it one-element columns
    drawn_columns = {}
    for label in model_definition.drawn_labels:
        drawn_columns[label] = np.array([drawn[label]])
    history_labels = model_definition.analytic_labels(drawn_columns, np.array([bagpipes_sfh.stellar_mass]))
    for label in HISTORY_LABELS:
        labels[label] = float(history_labels[label][0])

    ordered_labels = {}
    for label in model_definition.training_labels:
        ordered_labels[label] = labels[label]
    return ordered_labels


def ab_magnitudes(flux_per_angstrom):
    """AB magnitudes from Bagpipes photometry [erg/s/cm2/A]; all NaN if any flux is not positive.

    Inputs: flux_per_angstrom: F_lambda [erg/s/cm2/Angstrom], one per filter in PIVOT_WAVELENGTHS.
    Returns: AB magnitude per filter [mag], in the order of PIVOT_WAVELENGTHS.
    """
    if not np.all(flux_per_angstrom > 0):
        return np.full(len(PIVOT_WAVELENGTHS), np.nan)
    flux_per_hertz = flux_per_angstrom * PIVOT_WAVELENGTHS**2 / SPEED_OF_LIGHT
    return -MAGNITUDES_PER_DEX * np.log10(flux_per_hertz) - AB_ZERO_POINT


def simulate_models(unit_vectors, sfh_model):
    """Labels and AB magnitudes (10^10 M_sun formed) for a block of unit-cube vectors of one SFH model.

    Inputs: unit_vectors: iterable of unit-cube vectors (dimensionless, in [0, 1)), each ordered as
            MODELS[sfh_model].dims; sfh_model: a name in MODELS ("delayed" or "dblplaw").
    Returns: DataFrame with one row per model: the model's training_labels, then the AB magnitudes [mag] in
             BANDS order (NaN where a flux is not positive). The "z" column is the z-band magnitude, not the
             redshift.

    One Bagpipes model is built and then updated: update() is ~70x faster than building a new one, and safe
    because every model has the same component set (Bagpipes fixes the components when it builds a model).
    """
    model_definition = MODELS[sfh_model]
    model = None
    rows = []
    for u in unit_vectors:
        bagpipes_components, drawn = model_definition.components(u)
        if model is None:
            model = pipes.model_galaxy(bagpipes_components, filt_list=FILTER_FILES)
        else:
            model.update(bagpipes_components)

        row = summary_labels(model.sfh, drawn, sfh_model)
        magnitudes = ab_magnitudes(model.photometry)
        for band, magnitude in zip(BANDS, magnitudes):
            row[band] = magnitude
        rows.append(row)
    return pd.DataFrame(rows, columns=model_definition.training_labels + BANDS)


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


# Spectral fits are delayed-tau only (DPL spec, decision 1): the DPL model is fitted to the photometry alone
SPECTRAL_MODEL = "delayed"

# Redshift is fixed to the spectroscopic value, so it leaves the list; veldisp and noise_scale join it
SPECTRAL_DIMS = list(MODELS[SPECTRAL_MODEL].photometric_dims)
SPECTRAL_DIMS.append("veldisp")
SPECTRAL_DIMS.append("noise_scale")


def training_unit_vector(cube_by_name, sfh_model):
    """Unit-cube vector in a model's dims order for a fit at fixed redshift, from the fit's own named coordinates.

    Inputs: cube_by_name: dict of unit-cube values (dimensionless, in [0, 1)) holding every dimension of the model
            but "z"; other names (the fit's nuisance parameters) are ignored; sfh_model: "delayed" or "dblplaw".
    Returns: array of len(MODELS[sfh_model].dims) values (dimensionless); the "z" slot holds IGNORED_REDSHIFT_DRAW.
    """
    dims = MODELS[sfh_model].dims
    unit_vector = np.zeros(len(dims))
    for dimension_index, dimension in enumerate(dims):
        if dimension == "z":
            # components() ignores this value, because fixed_redshift replaces the drawn redshift
            unit_vector[dimension_index] = IGNORED_REDSHIFT_DRAW
        else:
            unit_vector[dimension_index] = cube_by_name[dimension]
    return unit_vector


def spectral_components(cube, redshift):
    """Bagpipes components, drawn parameters and noise scale for one spectral-fit unit-cube point.

    The SFH, metallicity, dust and nebular part is the delayed-tau training prior (MODELS[SPECTRAL_MODEL]) at the
    fixed spectroscopic redshift, which also sets the oldest allowed age; veldisp and noise_scale are the
    spectral nuisance parameters.
    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as SPECTRAL_DIMS;
            redshift: spectroscopic redshift (dimensionless).
    Returns: (bagpipes_components, drawn, noise_scale). The components are the delayed-tau model's plus
             veldisp [km/s]; drawn holds the model's drawn labels plus veldisp [km/s]; noise_scale is the error
             multiplier (dimensionless, 1 to 5).
    """
    cube_by_name = dict(zip(SPECTRAL_DIMS, cube))
    unit_vector = training_unit_vector(cube_by_name, SPECTRAL_MODEL)
    bagpipes_components, drawn = MODELS[SPECTRAL_MODEL].components(unit_vector, fixed_redshift=redshift)
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
JOB_DEFAULT_SFH_MODEL = "delayed"   # SFH model of a job or save_fit call naming none (every job before the DPL model)
EVIDENCE_ERROR_SAMPLES = 200        # posterior samples that recompute the ln Z error of a fit that stored none (count)


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
             (dimensionless ln Z), log_evidence_error (MultiNest's uncertainty of ln Z, sqrt(H / n_live)),
             n_iterations (count, from multinest_iterations), converged (bool: the run certainly stopped before
             max_iterations, judged by the upper bound of the count) and run_time [s].
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
    Returns: dict as fit_spectrum returns: samples (n_samples, n_dims), log_evidence, log_evidence_error,
             n_iterations, converged and run_time [s].
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
        # MultiNest's "global log-evidence error", sqrt(H / n_live), read by pymultinest from stats.dat
        "log_evidence_error": result["logZerr"],
        "n_iterations": n_iterations,
        "converged": max_possible_iterations < max_iterations,
        "run_time": run_time,
    }


def posterior_labels(unit_cube_samples, redshift):
    """Labels of each spectral-fit posterior sample, from a Bagpipes model that also computes S-PLUS photometry.

    Keeps at most MAX_POSTERIOR_SAMPLES evenly spaced samples, because building labels is the slow part.
    model_r is the model r magnitude for 10^10 M_sun formed, used later to rescale mass-like labels.
    Inputs: unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, len(SPECTRAL_DIMS)));
            redshift: spectroscopic redshift (dimensionless).
    Returns: DataFrame with one row per kept sample: the delayed-tau training labels (mass-like ones for
             10^10 M_sun formed), then veldisp [km/s], noise_scale (dimensionless) and model_r [mag].
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
        rows.append(sample_label_row(model, bagpipes_components, drawn_with_noise, SPECTRAL_MODEL))
    return pd.DataFrame(rows)


def thinned_samples(unit_cube_samples, max_samples=MAX_POSTERIOR_SAMPLES):
    """At most max_samples evenly spaced posterior samples, because computing anything per sample is the slow part.

    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, n_dims);
            max_samples: largest number of samples kept (count).
    Returns: every step-th row (dimensionless), step chosen so that no more than max_samples remain.
    """
    step = max(1, int(np.ceil(len(unit_cube_samples) / max_samples)))
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


def sample_label_row(model, bagpipes_components, drawn, sfh_model):
    """Labels of one posterior sample, computed from its photometry model as for a training-set model.

    Inputs: model: photometry model_galaxy already at bagpipes_components; bagpipes_components: dict as from the
            model's components(), the point the model is at (not read: the labels come from drawn and the model);
            drawn: dict from components() with the model's drawn labels, plus any extra parameters of the fit
            (veldisp [km/s], noise_scale); sfh_model: "delayed" or "dblplaw".
    Returns: dict: the model's training labels in order (summary_labels), then the extra drawn parameters, then
             model_r [mag], the model r magnitude for 10^10 M_sun formed. Mass-like labels are for that mass.
    """
    row = summary_labels(model.sfh, drawn, sfh_model)
    for key, value in drawn.items():
        if key not in row:
            row[key] = value
    row["model_r"] = ab_magnitudes(model.photometry)[BANDS.index("r")]
    return row


# -----------------------------------------------------------------------------
# Photometric fit (spec amendment A): the ML's own input, fitted by Bagpipes with the simulator's likelihood
# -----------------------------------------------------------------------------
def photometric_components(cube, redshift, sfh_model):
    """Bagpipes components and drawn parameters for one photometric-fit unit-cube point of one SFH model.

    Inputs: cube: unit-cube values (dimensionless, in [0, 1)) ordered as MODELS[sfh_model].photometric_dims;
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw".
    Returns: (bagpipes_components, drawn), exactly as MODELS[sfh_model].components(..., fixed_redshift=redshift)
             gives them.
    """
    cube_by_name = dict(zip(MODELS[sfh_model].photometric_dims, cube))
    unit_vector = training_unit_vector(cube_by_name, sfh_model)
    return MODELS[sfh_model].components(unit_vector, fixed_redshift=redshift)


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


def starting_photometry_model(redshift, sfh_model):
    """The photometry model a photometric fit starts from: built once, at the midpoint of the unit cube.

    Inputs: redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw".
    Returns: a model_galaxy with the S-PLUS filters. Every point of a fit has the same component set, so the
             model is then updated in place (photometric_point_log_likelihood).
    """
    midpoint_cube = np.full(len(MODELS[sfh_model].photometric_dims), MIDPOINT_DRAW)
    starting_components, _ = photometric_components(midpoint_cube, redshift, sfh_model)
    return updated_photometry_model(None, starting_components)


def photometric_point_log_likelihood(model, cube, redshift, sfh_model, magnitudes, errors, detected):
    """Log-likelihood of the photometry at one unit-cube point, as the photometric fit's sampler computes it.

    Why one function: fit_photometry's likelihood and photometric_log_likelihoods (the evidence error) must
    give the same number at the same point.
    Inputs: model: photometry model_galaxy from starting_photometry_model, updated in place to the point;
            cube: unit-cube point (dimensionless, in MODELS[sfh_model].photometric_dims order);
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw";
            magnitudes, errors: [mag] the fitted photometry; detected: boolean mask from detected_bands.
    Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
    """
    bagpipes_components, _ = photometric_components(cube, redshift, sfh_model)
    # PyMultiNest calls sys.exit if the likelihood raises: a model Bagpipes cannot compute is impossible
    try:
        model.update(bagpipes_components)
        model_magnitudes = ab_magnitudes(model.photometry)
    except Exception:
        return LOG_ZERO
    return safe_photometric_log_likelihood(model_magnitudes, magnitudes, errors, detected)


def fit_photometry(magnitudes, errors, redshift, sfh_model, output_basename, max_iterations,
                   n_live_points=N_LIVE_POINTS):
    """MultiNest fit of one galaxy's 12-band photometry under one SFH model's training prior at a fixed redshift.

    It fits the ML's own target: the model's training prior (its photometric_dims) and the likelihood of the
    simulator that made the training data (offset_marginal_log_likelihood), over the detected bands only.
    Inputs: magnitudes, errors: [mag], 12 values in BANDS order; a band is left out when detected_bands says
            it is not detected (99 or NaN magnitude, or error >= MAX_DETECTED_ERROR);
            redshift: spectroscopic redshift (dimensionless); sfh_model: "delayed" or "dblplaw";
            output_basename, max_iterations and n_live_points: as in fit_spectrum.
    Returns: dict as fit_spectrum returns, with samples of shape (n_samples, len(MODELS[sfh_model].photometric_dims)).
    Raises: ValueError when fewer than MIN_DETECTED_BANDS bands are detected.
    """
    # 1. Bands that enter the likelihood
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)
    if detected.sum() < MIN_DETECTED_BANDS:
        raise ValueError(f"{detected.sum()} detected bands, fewer than MIN_DETECTED_BANDS = {MIN_DETECTED_BANDS}")

    # 2. Bagpipes model built once with the S-PLUS filters; every point has the same component set
    model = starting_photometry_model(redshift, sfh_model)

    # 3. Likelihood of one unit-cube point
    def log_likelihood(cube):
        """Log-likelihood of the photometry for one unit-cube point (photometric_point_log_likelihood).

        Inputs: cube: unit-cube point (dimensionless, in the model's photometric_dims order).
        Returns: log-likelihood (dimensionless), LOG_ZERO for an impossible model.
        """
        return photometric_point_log_likelihood(model, cube, redshift, sfh_model, magnitudes, errors,
                                                detected)

    # 4. MultiNest over the unit cube
    n_dims = len(MODELS[sfh_model].photometric_dims)
    return run_multinest(log_likelihood, n_dims, output_basename, max_iterations, n_live_points)


def photometric_posterior_labels(unit_cube_samples, redshift, magnitudes, errors, rng, sfh_model):
    """Labels of each photometric-fit posterior sample, at the brightness of the observed galaxy.

    Each kept sample gets its own magnitude offset, drawn from the offset's Gaussian given that sample's
    model (offset_marginal_log_likelihood), and its mass-like labels move by it (rescale_mass_labels), as the
    training labels move with the brightness offset of each model.
    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, len(MODELS[sfh_model].photometric_dims));
            redshift: spectroscopic redshift (dimensionless); magnitudes, errors: [mag] the fitted photometry,
            as given to fit_photometry; rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: DataFrame with one row per kept sample (at most MAX_POSTERIOR_SAMPLES): the model's training labels
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
        bagpipes_components, drawn = photometric_components(cube, redshift, sfh_model)
        model = updated_photometry_model(model, bagpipes_components)
        row = sample_label_row(model, bagpipes_components, drawn, sfh_model)
        model_magnitudes = ab_magnitudes(model.photometry)
        _, best_offset, offset_sigma = offset_marginal_log_likelihood(model_magnitudes, magnitudes, errors, detected)
        row["magnitude_offset"] = rng.normal(best_offset, offset_sigma)
        rows.append(row)

    # 2. Mass-like labels at the observed brightness
    labels = pd.DataFrame(rows)
    return rescale_mass_labels(labels, labels["magnitude_offset"].to_numpy())


# -----------------------------------------------------------------------------
# Evidence uncertainty (DPL spec, part 4): stored by every new fit, recomputed for fits saved without it
# -----------------------------------------------------------------------------
def evidence_error(log_likelihoods, log_evidence, n_live_points):
    """MultiNest's uncertainty of ln Z, sqrt(H / n_live), with the information H taken from posterior samples.

    H = E_posterior[ln L] - ln Z is the information gained from prior to posterior [nats]; MultiNest reports
    sqrt(H / n_live) as the error of its ln Z.
    Inputs: log_likelihoods: ln L of equal-weight posterior samples (dimensionless), with the fit's own likelihood;
            log_evidence: the fit's ln Z (dimensionless); n_live_points: live points of the fit (count).
    Returns: the ln Z uncertainty (dimensionless float); 0 when the samples give H <= 0 (no information).
    """
    information = np.mean(np.asarray(log_likelihoods, dtype=float)) - log_evidence
    return float(np.sqrt(max(information, 0.0) / n_live_points))


def photometric_log_likelihoods(unit_cube_samples, redshift, magnitudes, errors, sfh_model,
                                max_samples=EVIDENCE_ERROR_SAMPLES):
    """The photometric fit's log-likelihood at its posterior samples, as fit_photometry computed it.

    Inputs: unit_cube_samples: unit-cube points (dimensionless), shape (n_samples, len(MODELS[sfh_model].photometric_dims));
            redshift: spectroscopic redshift (dimensionless); magnitudes, errors: [mag] the fitted photometry;
            sfh_model: the fit's model; max_samples: at most this many evenly spaced samples are used (count).
    Returns: array of log-likelihoods (dimensionless), one per used sample.
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    detected = detected_bands(magnitudes, errors)
    model = starting_photometry_model(redshift, sfh_model)
    log_likelihoods = []
    for cube in thinned_samples(unit_cube_samples, max_samples):
        log_likelihoods.append(photometric_point_log_likelihood(model, cube, redshift, sfh_model, magnitudes,
                                                                errors, detected))
    return np.array(log_likelihoods)


def load_unit_cube_samples(path):
    """Unit-cube posterior samples of a file written by save_fit.

    Inputs: path: file path (str). Returns: array (n_samples, n_dims), dimensionless.
    """
    with h5py.File(path, "r") as source:
        return source["unit_cube_samples"][()]


def fit_evidence(fit_file, job):
    """ln Z of a photometric fit and its uncertainty: the stored one, or sqrt(H / n_live) recomputed.

    Why recompute: fits saved before the DPL model stored no evidence error, and their MultiNest folders stayed
    in the PC job; the recomputed value is MultiNest's own estimate, from the stored posterior samples. A fit
    without n_live_points ran with N_LIVE_POINTS (submit_fits.sh never sets --n-live).
    Inputs: fit_file: path of a photometric fit (str); job: its photometry job (magnitudes and errors [mag],
            the photometry the fit used).
    Returns: (log_evidence, log_evidence_error), dimensionless floats.
    Raises: ValueError when the fit's kind is not PHOTOMETRY_JOB (a fit stored without a kind is spectral):
            its error could not be recomputed with the photometric likelihood.
    """
    _, metadata = load_fit(fit_file)
    stored_kind = metadata.get("kind", SPECTRUM_JOB)
    if stored_kind != PHOTOMETRY_JOB:
        raise ValueError(f"{fit_file} is not a photometric fit (kind {stored_kind!r})")
    if "log_evidence_error" in metadata:
        return float(metadata["log_evidence"]), float(metadata["log_evidence_error"])
    samples = load_unit_cube_samples(fit_file)
    log_likelihoods = photometric_log_likelihoods(samples, metadata["redshift"], job["magnitudes"], job["errors"],
                                                  metadata["sfh_model"])
    n_live_points = metadata.get("n_live_points", N_LIVE_POINTS)
    return float(metadata["log_evidence"]), evidence_error(log_likelihoods, metadata["log_evidence"], n_live_points)


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
    records its sfh_model, so fits of an earlier or of the other SFH model can be told apart (fit_is_current).
    Inputs: path: output file path (str); posterior: DataFrame of labels, one column per label;
            unit_cube_samples: unit-cube points (dimensionless, shape (n_samples, n_dims), n_dims that of the
            fit: len(SPECTRAL_DIMS) or the model's photometric_dims);
            metadata: dict of scalars (converged, run_time [s], n_iterations, max_iterations, n_live_points,
            log_evidence, log_evidence_error, redshift, kind, sfh_model and, for a cross-fit, truth_model; files
            written before kinds existed have no kind and are spectral). Without sfh_model, the fit is stored as
            JOB_DEFAULT_SFH_MODEL. The caller's dict is not changed.
    Returns: None. An existing file at path is replaced.
    """
    stored_metadata = dict(metadata)
    stored_metadata["sfh_model"] = metadata.get("sfh_model", JOB_DEFAULT_SFH_MODEL)
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
             delayed-tau model have no sfh_model, and files written before the DPL model have no
             log_evidence_error and no n_live_points).
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


def fit_is_current(metadata, sfh_model):
    """Whether a stored fit belongs to the current production run of a fits folder's SFH model.

    Fits that fail this (16-bin fits without sfh_model, fits of the other model, pilot runs with a smaller limit)
    are moved to superseded/ by supersede_outdated_fits (spectroscopic_reference.ipynb, section 4), never read
    as results. The folder's model comes from FIT_FOLDER_MODELS.
    Inputs: metadata: dict as load_fit returns it; sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: bool: True only when metadata["sfh_model"] == sfh_model and max_iterations >= MAX_ITERATIONS.
    """
    same_model = metadata.get("sfh_model") == sfh_model
    full_run = metadata.get("max_iterations", 0) >= MAX_ITERATIONS
    return same_model and full_run


def job_kind(job):
    """Kind of a fit job: SPECTRUM_JOB or PHOTOMETRY_JOB.

    Inputs: job: job dict (as fit_and_save takes).
    Returns: job["kind"], or SPECTRUM_JOB for a job without one (jobs made before the photometric fits existed).
    """
    return job.get("kind", SPECTRUM_JOB)


def job_model(job):
    """SFH model a fit job is fitted with.

    Inputs: job: job dict (as fit_and_save takes).
    Returns: job["sfh_model"], or JOB_DEFAULT_SFH_MODEL ("delayed") for a job without one (every job exported
             before the DPL model).
    """
    return job.get("sfh_model", JOB_DEFAULT_SFH_MODEL)


def fit_spectrum_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a spectrum job and label its posterior samples.

    Inputs: job: dict with spectrum (as from prepare_spectrum) and redshift (dimensionless); multinest_basename:
            directory prefix for the MultiNest files (str); max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_spectrum, posterior DataFrame from posterior_labels).
    Raises: ValueError for a job of a model other than SPECTRAL_MODEL (spectral fits are delayed-tau only).
    """
    if job_model(job) != SPECTRAL_MODEL:
        raise ValueError(f"spectral fits are {SPECTRAL_MODEL} only, not {job_model(job)} (DPL spec, decision 1)")
    fit = fit_spectrum(job["spectrum"], job["redshift"], multinest_basename, max_iterations, n_live_points)
    posterior = posterior_labels(fit["samples"], job["redshift"])
    return fit, posterior


def fit_photometry_job(job, multinest_basename, max_iterations, n_live_points):
    """Fit a photometry job with its SFH model and label its posterior samples at the galaxy's brightness.

    Inputs: job: dict with magnitudes and errors ([mag], 12 in BANDS order), redshift (dimensionless) and
            sfh_model (job_model); multinest_basename: directory prefix for the MultiNest files (str);
            max_iterations, n_live_points: counts.
    Returns: (fit dict from fit_photometry, posterior DataFrame from photometric_posterior_labels).
    """
    sfh_model = job_model(job)
    fit = fit_photometry(job["magnitudes"], job["errors"], job["redshift"], sfh_model, multinest_basename,
                         max_iterations, n_live_points)
    # A fixed seed, so labelling the same samples again gives the same offsets
    rng = np.random.default_rng(OFFSET_DRAW_SEED)
    posterior = photometric_posterior_labels(fit["samples"], job["redshift"], job["magnitudes"], job["errors"],
                                             rng, sfh_model)
    return fit, posterior


# The fitter of each job kind: same signature, (job, multinest_basename, max_iterations, n_live_points)
JOB_FITTERS = {SPECTRUM_JOB: fit_spectrum_job, PHOTOMETRY_JOB: fit_photometry_job}


def fit_and_save(job, output_directory, max_iterations, n_live_points=N_LIVE_POINTS):
    """Fit one job and save <output_directory>/<job_id>.h5. Never raises.

    Skips jobs whose readable file of the same kind and SFH model already exists, so a batch can be re-run
    with more galaxies. An unreadable file does not count as finished: the job is fitted again and the file
    replaced. A readable file of the other kind is never replaced (job ids are specobjids, shared by the
    spectral and photometric fits of one galaxy), and neither is a fit of another SFH model (an earlier model's
    result, the only copy, to be moved to superseded/): the job fails instead. A failed fit comes back with
    status "failed" and the error text instead of an exception.
    Inputs: job: dict with job_id (str), redshift (dimensionless), kind (SPECTRUM_JOB when absent), sfh_model
            (JOB_DEFAULT_SFH_MODEL when absent) and, for a cross-fit, truth_model; a spectrum job holds spectrum
            (dict as from prepare_spectrum), a photometry job magnitudes and errors ([mag], 12 in BANDS order);
            output_directory: folder for the .h5 file and the MultiNest files (str); max_iterations: MultiNest
            iteration limit (count); n_live_points: live points (count).
    Returns: dict with job_id, status ("fitted", "skipped" or "failed"), converged (bool or None),
             run_time [s] (or None), n_iterations (count or None) and error (str or None).
    """
    # 1. Output path and a summary that the steps below fill in
    output_path = os.path.join(output_directory, f"{job['job_id']}.h5")
    kind = job_kind(job)
    sfh_model = job_model(job)
    summary = {"job_id": job["job_id"], "status": None, "converged": None,
               "run_time": None, "n_iterations": None, "error": None}

    # 2. Skip finished jobs; their run details come from the file. Any failure to read it (damaged,
    # truncated or missing fields) means the job is not finished, so it is fitted again below. A readable
    # fit of another SFH model is neither skipped (not a result of this model) nor overwritten (only copy)
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
    if stored_kind == kind and stored_model != sfh_model:
        summary.update(status="failed", error=f"{output_path} holds a fit of another SFH model "
                                              f"({stored_model}, not {sfh_model}): move it to superseded/")
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
                    "n_live_points": n_live_points, "log_evidence": fit["log_evidence"],
                    "log_evidence_error": fit["log_evidence_error"], "redshift": job["redshift"], "kind": kind,
                    "sfh_model": sfh_model}
        if "truth_model" in job:
            metadata["truth_model"] = job["truth_model"]
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
CALIBRATION_LEVELS = (0.68, 0.95)                   # central-interval probabilities reported (fractions)
CALIBRATION_ID_PREFIX = "calib_"                    # job ids are calib_000, calib_001, ...
CROSS_CALIBRATION_ID_PREFIX = "cross_calib_"   # cross-fit ids: cross_calib_007 is mock calib_007 fitted with the other model
# Reference-CSV columns of each model's photometric fits: <prefix><label>_p16, <prefix>log_evidence, ...
PHOTOMETRIC_PREFIXES = {"delayed": "phot_", "dblplaw": "phot_dblplaw_"}
PHOTOMETRIC_EXTRA_LABELS = ["mw_age", "log_sfr"]   # summarised besides the calibration and the shared labels
PERCENT = 100                    # percent per unit fraction


def export_job(job, jobs_directory):
    """Write one fit job to <jobs_directory>/<job_id>.npz, a portable file for run_spectral_fits.py.

    Inputs: job: dict with job_id (str), redshift (dimensionless), kind (SPECTRUM_JOB when absent), sfh_model
            (JOB_DEFAULT_SFH_MODEL when absent) and, for a cross-fit, truth_model; a spectrum job holds spectrum
            (dict with wavelength [Angstrom, observed], flux and ivar), a photometry job magnitudes and errors
            ([mag], 12 in BANDS order); jobs_directory: folder for the file (str), created if missing.
    Returns: path of the written file (str).
    """
    os.makedirs(jobs_directory, exist_ok=True)
    path = os.path.join(jobs_directory, f"{job['job_id']}.npz")
    kind = job_kind(job)
    arrays = {"job_id": str(job["job_id"]), "redshift": job["redshift"], "kind": kind, "sfh_model": job_model(job)}
    if "truth_model" in job:
        arrays["truth_model"] = job["truth_model"]
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
    Returns: job dict in the format fit_and_save takes: job_id (str), kind, sfh_model, redshift (dimensionless
             float), truth_model when the file has one, and either spectrum (wavelength, flux, ivar) or magnitudes
             and errors [mag]. A file without a kind (exported before the photometric fits existed) is a spectrum
             job; a file without sfh_model (exported before the DPL model) is a JOB_DEFAULT_SFH_MODEL job.
    """
    with np.load(path) as source:
        kind = SPECTRUM_JOB
        if "kind" in source.files:
            kind = str(source["kind"])
        sfh_model = JOB_DEFAULT_SFH_MODEL
        if "sfh_model" in source.files:
            sfh_model = str(source["sfh_model"])
        job = {"job_id": str(source["job_id"]), "kind": kind, "sfh_model": sfh_model,
               "redshift": float(source["redshift"])}
        if "truth_model" in source.files:
            job["truth_model"] = str(source["truth_model"])
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
             spectral model's training_labels, veldisp, noise_scale and model_r).
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


# -----------------------------------------------------------------------------
# Photometric calibration mocks of any SFH model, and the cross-fits of the evidence test (DPL spec, part 2)
# -----------------------------------------------------------------------------
def draw_photometric_cubes(n_mocks, sfh_model, seed):
    """Unit-cube points of photometric calibration mocks: a Latin hypercube over the model's photometric dims.

    Why a Latin hypercube: as for the training sets, so that 33 mocks still cover every parameter range evenly.
    Inputs: n_mocks: number of points (count); sfh_model: "delayed" or "dblplaw";
            seed: Latin hypercube seed (integer or numpy Generator).
    Returns: array (n_mocks, len(MODELS[sfh_model].photometric_dims)), values in [0, 1) (dimensionless).
    """
    sampler = qmc.LatinHypercube(len(MODELS[sfh_model].photometric_dims), seed=seed)
    return sampler.random(n_mocks)


def noiseless_photometric_mocks(cubes, redshifts, sfh_model):
    """Noiseless S-PLUS photometry and true labels of known photometric-fit points, for 10^10 M_sun formed.

    Inputs: cubes: unit-cube points (dimensionless), shape (n_mocks, len(MODELS[sfh_model].photometric_dims));
            redshifts: one redshift per point (dimensionless); sfh_model: "delayed" or "dblplaw".
    Returns: (magnitudes [mag], shape (n_mocks, 12) in BANDS order, NaN where a flux is not positive;
             truths: list of Series, one per point, as sample_label_row gives them: the model's training labels,
             mass-like ones for 10^10 M_sun formed, and model_r [mag]).
    """
    magnitudes = np.zeros((len(cubes), len(BANDS)))
    truths = []
    model = None
    for mock_index, cube in enumerate(cubes):
        bagpipes_components, drawn = photometric_components(cube, redshifts[mock_index], sfh_model)
        model = updated_photometry_model(model, bagpipes_components)
        magnitudes[mock_index] = ab_magnitudes(model.photometry)
        truths.append(pd.Series(sample_label_row(model, bagpipes_components, drawn, sfh_model)))
    return magnitudes, truths


def photometric_mock_jobs(job_ids, noiseless_magnitudes, redshifts, magnitude_offsets, truths, noise_table, rng,
                          sfh_model):
    """Photometry jobs of mocks at given brightness, noised once as the training data, with their true labels.

    Each mock moves by its magnitude offset (observed minus model, the same in every band), then gets S-PLUS
    errors, flux noise and non-detections (add_splus_noise); its true mass-like labels move with the offset
    (rescale_mass_labels). A mock left with fewer than MIN_DETECTED_BANDS detected bands cannot be fitted and is
    left out.
    Inputs: job_ids: one id per mock (str); noiseless_magnitudes: [mag] (n_mocks, 12) in BANDS order, for
            10^10 M_sun formed; redshifts: one per mock (dimensionless); magnitude_offsets: [mag] one per mock
            (NaN gives a mock with no detection, so it is left out); truths: dict job_id -> Series of the true
            labels for 10^10 M_sun formed; noise_table: output of fit_noise_model; rng: numpy Generator (only
            add_splus_noise draws from it); sfh_model: the model the jobs are fitted with.
    Returns: (jobs, rescaled_truths): jobs is a list of photometry job dicts as fit_and_save takes (job_id, kind,
             sfh_model, magnitudes, errors [mag], redshift), in input order without the left-out mocks;
             rescaled_truths maps each kept job_id to its truth at the mock's brightness.
    """
    offsets = np.asarray(magnitude_offsets, dtype=float)
    scaled_magnitudes = noiseless_magnitudes + offsets[:, np.newaxis]
    noisy_magnitudes, noisy_errors, _ = add_splus_noise(scaled_magnitudes, noise_table, rng)

    jobs = []
    rescaled_truths = {}
    for mock_index, job_id in enumerate(job_ids):
        detected = detected_bands(noisy_magnitudes[mock_index], noisy_errors[mock_index])
        if detected.sum() < MIN_DETECTED_BANDS:
            continue
        jobs.append({"job_id": job_id, "kind": PHOTOMETRY_JOB, "sfh_model": sfh_model,
                     "magnitudes": noisy_magnitudes[mock_index], "errors": noisy_errors[mock_index],
                     "redshift": float(redshifts[mock_index])})
        truth = truths[job_id].to_frame().T
        rescaled_truths[job_id] = rescale_mass_labels(truth, offsets[mock_index]).iloc[0]
    return jobs, rescaled_truths


def draw_photometric_calibration_mocks(n_mocks, sfh_model, donor_redshifts, donor_r, noise_table, seed):
    """Photometric calibration mocks from a model's prior, each at the redshift and brightness of a real galaxy.

    Why a donor galaxy: the mock then has a real redshift and aperture brightness, and so realistic S-PLUS errors.
    The unit-cube points (draw_photometric_cubes) and the donors+noise come from two independent streams
    spawned from numpy.random.SeedSequence(seed), donors first in the second. Why two streams: with one seed
    for both, the donors and noise reuse the Latin hypercube's own uniforms, so the mocks are not independent.
    Inputs: n_mocks: number of mocks (count); sfh_model: "delayed" or "dblplaw"; donor_redshifts: spectroscopic
            redshifts of the candidate donors (dimensionless); donor_r: their aperture r magnitudes [mag] (a donor
            with a non-finite r, or r at or above UNDETECTED_MAGNITUDE, is never chosen); noise_table: output of
            fit_noise_model; seed: integer.
    Returns: (jobs, truths) as photometric_mock_jobs returns them, with ids calib_000, calib_001, ... (a mock with
             too few detected bands is left out, so an id can be missing); truths at the donor's brightness.
    Raises: ValueError when no donor is usable.
    """
    donor_redshifts = np.asarray(donor_redshifts, dtype=float)
    donor_r = np.asarray(donor_r, dtype=float)
    usable = np.isfinite(donor_r) & (donor_r < UNDETECTED_MAGNITUDE) & np.isfinite(donor_redshifts)
    if not usable.any():
        raise ValueError("no usable brightness donor (finite redshift and detected r)")
    usable_redshifts = donor_redshifts[usable]
    usable_r = donor_r[usable]

    # 1. Points of the prior, and one donor galaxy per mock
    cube_seed, donor_seed = np.random.SeedSequence(seed).spawn(2)
    cubes = draw_photometric_cubes(n_mocks, sfh_model, np.random.default_rng(cube_seed))
    rng = np.random.default_rng(donor_seed)
    donors = rng.integers(len(usable_r), size=n_mocks)
    redshifts = usable_redshifts[donors]

    # 2. Noiseless photometry for 10^10 M_sun formed, moved to the donor's r magnitude, then noised once
    noiseless_magnitudes, truth_list = noiseless_photometric_mocks(cubes, redshifts, sfh_model)
    magnitude_offsets = usable_r[donors] - noiseless_magnitudes[:, BANDS.index("r")]
    job_ids = []
    truths = {}
    for mock_index in range(n_mocks):
        job_id = f"{CALIBRATION_ID_PREFIX}{mock_index:03d}"
        job_ids.append(job_id)
        truths[job_id] = truth_list[mock_index]
    return photometric_mock_jobs(job_ids, noiseless_magnitudes, redshifts, magnitude_offsets, truths, noise_table,
                                 rng, sfh_model)


def cross_fit_id(mock_id):
    """Job id of the cross-fit of a calibration mock: calib_007 -> cross_calib_007.

    Inputs: mock_id: calibration-mock id (str, CALIBRATION_ID_PREFIX followed by the three-digit index).
    Returns: the cross-fit id (str).
    Raises: ValueError for an id that is not a calibration-mock id.
    """
    if not mock_id.startswith(CALIBRATION_ID_PREFIX):
        raise ValueError(f"{mock_id} is not a calibration-mock id ({CALIBRATION_ID_PREFIX}<k>)")
    return CROSS_CALIBRATION_ID_PREFIX + mock_id[len(CALIBRATION_ID_PREFIX):]


def cross_fit_jobs(mock_jobs, mock_truths, fitted_model):
    """One model's photometric mocks as jobs fitted with the other model, for the evidence test.

    Why: with the truth known, the evidence of the true model should win; if it does not, the colours cannot
    tell the two models apart. The photometry, redshift and truth of each mock are kept unchanged.
    Inputs: mock_jobs: photometry jobs of one model's calibration mocks (as in its mock file; a job without
            sfh_model is delayed-tau); mock_truths: dict mock job_id -> Series of the true labels;
            fitted_model: the model to fit them with, not the mocks' own.
    Returns: (jobs, truths): one job per mock with job_id cross_fit_id(mock id), sfh_model = fitted_model and
             truth_model = the mock's own model; truths maps the new ids to the mocks' truths. The input jobs are
             not changed.
    Raises: ValueError when a mock is already of fitted_model.
    """
    jobs = []
    truths = {}
    for mock_job in mock_jobs:
        truth_model = job_model(mock_job)
        if truth_model == fitted_model:
            raise ValueError(f"{mock_job['job_id']} is a {truth_model} mock: a cross-fit uses the other model")
        job = dict(mock_job)
        job["job_id"] = cross_fit_id(mock_job["job_id"])
        job["sfh_model"] = fitted_model
        job["truth_model"] = truth_model
        jobs.append(job)
        truths[job["job_id"]] = mock_truths[mock_job["job_id"]]
    return jobs, truths


def photometric_summary_labels(sfh_model):
    """Labels of a model's photometric fits that the reference CSV summarises.

    Inputs: sfh_model: "delayed" or "dblplaw".
    Returns: list of label names, each once: the model's calibration labels, then PHOTOMETRIC_EXTRA_LABELS, then
             the SHARED_LABELS not listed yet (the labels the model comparison uses).
    """
    labels = list(MODELS[sfh_model].calibration_labels)
    for label in [*PHOTOMETRIC_EXTRA_LABELS, *SHARED_LABELS]:
        if label not in labels:
            labels.append(label)
    return labels


# =============================================================================
# 8. Superseding the results of an earlier model or run (spectroscopic_reference.ipynb, section 4)
# =============================================================================
SUPERSEDED_SUBDIRECTORY = "superseded"     # inside each fits folder: results of earlier, replaced runs
SUPERSEDED_NAME_SUFFIX = "_16bin"          # added to a name already taken there (first-run fits, same id)
CALIBRATION_MOCK_FILE = "calibration_mocks.joblib"   # (jobs, truths) of a fits folder's calibration mocks
CROSS_FIT_MOCK_FILE = "cross_fit_mocks.joblib"   # (jobs, truths) of the cross-fits exported into a fits folder
JOBS_SUBDIRECTORY = "jobs"                 # exported fit jobs, one .npz each
RUN_LOG_NAMES = ("fit_log.csv", "fit_logs")   # logs of run_spectral_fits.py; their rows name the fits run
# The SFH model of each fits folder: a fit of another model in it is superseded, never read
FIT_FOLDER_MODELS = {"spectral_fits": "delayed", "photometric_fits": "delayed", "photometric_fits_dblplaw": "dblplaw"}


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
        new_name = f"{stem}{SUPERSEDED_NAME_SUFFIX}{counter}{extension}"
        destination = os.path.join(superseded_directory, new_name)
        copy_number += 1
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    shutil.move(source, destination)
    return destination


def calibration_mocks_are_current(mock_file, sfh_model):
    """Whether a stored calibration-mock file holds truths with the calibration labels of a folder's model.

    Inputs: mock_file: path of a joblib file with (jobs, truths), truths mapping job_id -> Series of the
            true labels; sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: bool: True when there is at least one truth and every truth has every
             MODELS[sfh_model].calibration_labels value (16-bin truths lack log_age, delayed-tau truths lack
             log_turnover_lookback). A file without mocks is not current, so it is redrawn, not reused.
    """
    _, truths = joblib.load(mock_file)
    if len(truths) == 0:
        return False
    for truth in truths.values():
        if not set(MODELS[sfh_model].calibration_labels) <= set(truth.index):
            return False
    return True


def cross_fit_mocks_are_current(mock_file, sfh_model):
    """Whether a stored cross-fit mock file belongs to a fits folder of sfh_model and has current truths.

    Inputs: mock_file: path of a joblib file with (jobs, truths) as cross_fit_jobs returns them;
            sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: bool: True when there is at least one job, every job is fitted with sfh_model, and every truth has
             the calibration labels of its job's truth_model.
    """
    jobs, truths = joblib.load(mock_file)
    if len(jobs) == 0:
        return False
    for job in jobs:
        truth_model = job.get("truth_model")
        if job_model(job) != sfh_model or truth_model not in MODELS or job["job_id"] not in truths:
            return False
        if not set(MODELS[truth_model].calibration_labels) <= set(truths[job["job_id"]].index):
            return False
    return True


def supersede_outdated_fits(fits_directory, sfh_model):
    """Move every result of an earlier model or run out of a fits folder, into its SUPERSEDED_SUBDIRECTORY.

    Why: run_spectral_fits.py counts any top-level .h5 as done, and new calibration mocks reuse the calib_<k>
    ids, so an old file left in place would be skipped as finished or scored against a new truth. Each folder is
    judged against its own SFH model (FIT_FOLDER_MODELS): a DPL fit stays in photometric_fits_dblplaw and leaves
    photometric_fits, and the reverse.
    Moved, in this order:
      1. every readable .h5 fit that fails fit_is_current(metadata, sfh_model), with its MultiNest folder (the
         name without .h5); an unreadable .h5 stays, because fit_and_save fits that job again;
      2. CALIBRATION_MOCK_FILE when its truths lack the model's calibration labels (16-bin mocks) or it holds none,
         and CROSS_FIT_MOCK_FILE when cross_fit_mocks_are_current says it does not belong here;
      3. while no mock file is left, the jobs/calib_*.npz files, which then belong to no current mock, and while
         no cross-fit mock file is left, the jobs/cross_calib_*.npz files;
      4. when any fit moved, the run logs RUN_LOG_NAMES (notebook 2 reads failures from them, and
         fetch_fits.sh drops new rows of ids already logged).
    Inputs: fits_directory: folder with the .h5 fits, jobs/ and the mock files (str); a missing folder is fine;
            sfh_model: the folder's model, "delayed" or "dblplaw".
    Returns: list of the moved names (str), relative to fits_directory; empty when everything is current.
    Raises: ValueError when the folder is one of FIT_FOLDER_MODELS and sfh_model is not its model.
    """
    folder_name = os.path.basename(os.path.normpath(fits_directory))
    if folder_name in FIT_FOLDER_MODELS and FIT_FOLDER_MODELS[folder_name] != sfh_model:
        raise ValueError(f"{fits_directory} holds {FIT_FOLDER_MODELS[folder_name]} fits, not {sfh_model}")
    if not os.path.isdir(fits_directory):
        return []
    moved = []

    # 1. Fits of an earlier model or run, or of the other model, with their MultiNest folders
    for name in sorted(os.listdir(fits_directory)):
        if not name.endswith(".h5"):
            continue
        try:
            _, metadata = load_fit(os.path.join(fits_directory, name))
        except (OSError, KeyError):
            continue
        if fit_is_current(metadata, sfh_model):
            continue
        for fit_name in (name, name[:-len(".h5")]):
            if move_to_superseded(fits_directory, fit_name) is not None:
                moved.append(fit_name)
    fits_were_moved = len(moved) > 0

    # 2. Mock files of an earlier model, or made for the other folder
    mock_file = os.path.join(fits_directory, CALIBRATION_MOCK_FILE)
    if os.path.exists(mock_file) and not calibration_mocks_are_current(mock_file, sfh_model):
        move_to_superseded(fits_directory, CALIBRATION_MOCK_FILE)
        moved.append(CALIBRATION_MOCK_FILE)
    cross_mock_file = os.path.join(fits_directory, CROSS_FIT_MOCK_FILE)
    if os.path.exists(cross_mock_file) and not cross_fit_mocks_are_current(cross_mock_file, sfh_model):
        move_to_superseded(fits_directory, CROSS_FIT_MOCK_FILE)
        moved.append(CROSS_FIT_MOCK_FILE)

    # 3. Job files of mocks that no stored mock file holds (the notebook stores mocks before exporting them)
    orphan_prefixes = []
    if not os.path.exists(mock_file):
        orphan_prefixes.append(CALIBRATION_ID_PREFIX)
    if not os.path.exists(cross_mock_file):
        orphan_prefixes.append(CROSS_CALIBRATION_ID_PREFIX)
    jobs_directory = os.path.join(fits_directory, JOBS_SUBDIRECTORY)
    if os.path.isdir(jobs_directory):
        for name in sorted(os.listdir(jobs_directory)):
            for prefix in orphan_prefixes:
                if name.startswith(prefix):
                    job_name = os.path.join(JOBS_SUBDIRECTORY, name)
                    move_to_superseded(fits_directory, job_name)
                    moved.append(job_name)

    # 4. Run logs of the moved fits
    if fits_were_moved:
        for name in RUN_LOG_NAMES:
            if move_to_superseded(fits_directory, name) is not None:
                moved.append(name)
    return moved
