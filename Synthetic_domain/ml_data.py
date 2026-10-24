"""
Network inputs and targets for the SFH posterior estimator (Phase 3).

Model tables (training set, held-out sets, real galaxies) become:
- inputs: 11 colours (band - r) with errors and detection flags, plus r, its error, z_obs and sigma_z;
- targets: 19 numbers (bin log-ratios, metallicity, dust, masses, redshift) that a flow can model.
Every physical definition (noise, brightness, bins, enrichment) comes from synthetic_model.
"""

import numpy as np
import pandas as pd

import synthetic_model as sm


# =============================================================================
# 1. Constants
# =============================================================================
R_INDEX = sm.BANDS.index("r")
COLOUR_BANDS = []
for band_name in sm.BANDS:
    if band_name != "r":
        COLOUR_BANDS.append(band_name)
N_COLOURS = len(COLOUR_BANDS)                    # 11
N_INPUTS = 3 * N_COLOURS + 4                     # colours, colour errors, flags, then 4 scalars = 37
SCALAR_INPUT_SLICE = slice(3 * N_COLOURS, N_INPUTS)   # [m_r, sigma_r, z_obs, sigma_z]
N_BAND_CHANNELS = 4                              # [m_b - m_r, sigma_b, detected, log10 pivot]

N_FRACTION_BINS = sm.N_BINS - 1                  # 15 populated bins (index 15 never exists)
REFERENCE_BIN = 13                               # oldest bin that exists at every redshift of the prior
LOG_RATIO_LIMIT = 12.0                           # [dex] clip of sampled log-ratios (flow tails)
P_SPECTROSCOPIC = 0.5                            # share of simulated galaxies with an exact redshift

RATIO_BINS = []
for bin_index in range(N_FRACTION_BINS):
    if bin_index != REFERENCE_BIN:
        RATIO_BINS.append(bin_index)
TARGET_NAMES = []
for bin_index in RATIO_BINS:
    TARGET_NAMES.append(f"log_ratio_bin_{bin_index:02d}")
TARGET_NAMES += ["log_Z_fin", "Av_ism", "log_mass_formed", "log_mstar", "redshift"]
N_TARGETS = len(TARGET_NAMES)                    # 19


# =============================================================================
# 2. Inputs
# =============================================================================
def colour_inputs(magnitudes, errors, z_obs, sigma_z):
    """Network inputs from magnitudes: colours against r, their errors, flags and the scalars.

    Inputs: magnitudes, errors: [mag] arrays (n, 12) in sm.BANDS order (NaN or 99 when undetected);
            z_obs: observed redshift (n,); sigma_z: its uncertainty (n,), 0 for spectroscopic redshifts.
    Returns: dict with
             flat (n, 37): [11 colours, 11 colour errors, 11 flags, m_r, sigma_r, z_obs, sigma_z] for the MLP;
             banded (n, 12, 4): per band [m_b - m_r, sigma_b, detected, log10 pivot] for the CNN;
             valid (n,): True where r is detected (only those galaxies may be predicted).
    Unusable entries are 0, never NaN.
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    n_galaxies = len(magnitudes)
    detected = sm.detected_bands(magnitudes, errors)
    valid = detected[:, R_INDEX]

    r_magnitude = np.where(valid, magnitudes[:, R_INDEX], 0.0)
    r_error = np.where(valid, errors[:, R_INDEX], 0.0)

    flat = np.zeros((n_galaxies, N_INPUTS))
    banded = np.zeros((n_galaxies, len(sm.BANDS), N_BAND_CHANNELS))
    log_pivots = np.log10(sm.PIVOT_WAVELENGTHS)

    for colour_index, band in enumerate(COLOUR_BANDS):
        band_index = sm.BANDS.index(band)
        usable = detected[:, band_index] & valid
        colour = np.where(usable, magnitudes[:, band_index] - r_magnitude, 0.0)
        colour_error = np.where(usable, np.sqrt(errors[:, band_index] ** 2 + r_error**2), 0.0)
        flat[:, colour_index] = colour
        flat[:, N_COLOURS + colour_index] = colour_error
        flat[:, 2 * N_COLOURS + colour_index] = usable.astype(float)

    flat[:, SCALAR_INPUT_SLICE] = np.column_stack([r_magnitude, r_error, z_obs, sigma_z])

    for band_index in range(len(sm.BANDS)):
        usable = detected[:, band_index] & valid
        banded[:, band_index, 0] = np.where(usable, magnitudes[:, band_index] - r_magnitude, 0.0)
        banded[:, band_index, 1] = np.where(usable, errors[:, band_index], 0.0)
        banded[:, band_index, 2] = usable.astype(float)
        banded[:, band_index, 3] = log_pivots[band_index]
    return {"flat": flat, "banded": banded, "valid": valid}


# =============================================================================
# 3. Targets and decoding
# =============================================================================
def fraction_columns(n_bins):
    """Names frac_bin_00 ... of the first n_bins bins."""
    names = []
    for bin_index in range(n_bins):
        names.append(f"frac_bin_{bin_index:02d}")
    return names


def encode_targets(table):
    """The 19 flow targets of each row.

    Inputs: table with frac_bin_00..14, Z_fin [Z_sun], Av_ism [mag], log_mass_formed and log_mstar [log M_sun],
            redshift.
    Returns: array (n, 19) in TARGET_NAMES order. Bins that do not exist (fraction 0) take sm.MIN_FRACTION
             before the log-ratio, so every row has the same dimension.
    """
    fractions = table[fraction_columns(N_FRACTION_BINS)].to_numpy(dtype=float)
    fractions = np.where(fractions > 0, fractions, sm.MIN_FRACTION)
    reference = fractions[:, REFERENCE_BIN]

    columns = []
    for bin_index in RATIO_BINS:
        columns.append(np.log10(fractions[:, bin_index] / reference))
    columns.append(np.log10(table["Z_fin"].to_numpy(dtype=float)))
    columns.append(table["Av_ism"].to_numpy(dtype=float))
    columns.append(table["log_mass_formed"].to_numpy(dtype=float))
    columns.append(table["log_mstar"].to_numpy(dtype=float))
    columns.append(table["redshift"].to_numpy(dtype=float))
    return np.column_stack(columns)


def mass_formed_younger_than(fractions, edges, age_limit):
    """Share of the formed mass younger than age_limit, with constant SFR inside each bin.

    Inputs: fractions (n, N_BINS) summing to 1; edges (n, N_BINS + 1) [Gyr ago], NaN beyond the last edge;
            age_limit [Gyr].
    Returns: (n,) fraction of the mass formed more recently than age_limit.
    """
    starts = edges[:, :-1]
    ends = edges[:, 1:]
    widths = ends - starts
    overlap = np.clip(np.minimum(ends, age_limit) - starts, 0, None)
    share_of_bin = np.where(widths > 0, overlap / widths, 0.0)
    share_of_bin = np.nan_to_num(share_of_bin)
    return np.sum(fractions * share_of_bin, axis=1)


def decode_samples(theta):
    """Physical labels of flow samples, with the same definitions as the training labels.

    Inputs: theta: array (n, 19) in TARGET_NAMES order (posterior samples or encoded truths).
    Returns: DataFrame with redshift, Z_fin, Av_ism, log_mass_formed, log_mstar, log_sfr, mw_age, f_young,
             f_mid, f_old and the 48 sm.BIN_LABEL_COLUMNS. The redshift is clipped to the prior; bins that do
             not exist at the sampled redshift get frac 0, logsfr and zmet NaN, as in sm.bin_labels.
    """
    theta = np.asarray(theta, dtype=float)
    n_samples = len(theta)
    redshift = np.clip(theta[:, TARGET_NAMES.index("redshift")], *sm.Z_RANGE)
    z_final = 10 ** theta[:, TARGET_NAMES.index("log_Z_fin")]
    log_mass_formed = theta[:, TARGET_NAMES.index("log_mass_formed")]
    log_mstar = theta[:, TARGET_NAMES.index("log_mstar")]

    # 1. Bins that exist at each sampled redshift
    edges = sm.bin_edges_array(redshift)
    n_existing = np.isfinite(edges).sum(axis=1) - 1
    exists = np.arange(sm.N_BINS)[np.newaxis, :] < n_existing[:, np.newaxis]

    # 2. Fractions from the log-ratios (clipped: the flow's tails can be extreme)
    log_ratios = np.clip(theta[:, :len(RATIO_BINS)], -LOG_RATIO_LIMIT, LOG_RATIO_LIMIT)
    relative_mass = np.zeros((n_samples, sm.N_BINS))
    relative_mass[:, REFERENCE_BIN] = 1.0
    for ratio_index, bin_index in enumerate(RATIO_BINS):
        relative_mass[:, bin_index] = 10 ** log_ratios[:, ratio_index]
    relative_mass = np.where(exists, relative_mass, 0.0)
    fractions = relative_mass / relative_mass.sum(axis=1, keepdims=True)

    # 3. Per-bin labels, as sm.bin_labels defines them
    widths_in_years = np.diff(edges, axis=1) * sm.GYR_TO_YEARS
    with np.errstate(divide="ignore", invalid="ignore"):
        log_sfr_bins = log_mass_formed[:, np.newaxis] + np.log10(fractions) - np.log10(widths_in_years)
    log_sfr_bins = np.where(exists, log_sfr_bins, np.nan)
    metallicities = sm.bin_metallicities_array(fractions, n_existing, z_final)

    # 4. Summary labels from the bins (constant SFR inside each bin)
    midpoints = np.nan_to_num((edges[:, :-1] + edges[:, 1:]) / 2)
    mw_age = np.sum(fractions * midpoints, axis=1)
    younger_than_first = mass_formed_younger_than(fractions, edges, sm.LABEL_EDGES[0])
    younger_than_second = mass_formed_younger_than(fractions, edges, sm.LABEL_EDGES[1])
    recent_mass = 10**log_mass_formed * younger_than_first
    recent_sfr = recent_mass / (sm.LABEL_EDGES[0] * sm.GYR_TO_YEARS)
    sfr_floor = sm.SFR_FLOOR * 10**log_mstar
    log_sfr = np.log10(np.maximum(recent_sfr, sfr_floor))

    labels = {
        "redshift": redshift,
        "Z_fin": z_final,
        "Av_ism": theta[:, TARGET_NAMES.index("Av_ism")],
        "log_mass_formed": log_mass_formed,
        "log_mstar": log_mstar,
        "log_sfr": log_sfr,
        "mw_age": mw_age,
        "f_young": younger_than_first,
        "f_mid": younger_than_second - younger_than_first,
        "f_old": 1 - younger_than_second,
    }
    for bin_index in range(sm.N_BINS):
        labels[f"frac_bin_{bin_index:02d}"] = fractions[:, bin_index]
    for bin_index in range(sm.N_BINS):
        labels[f"logsfr_bin_{bin_index:02d}"] = log_sfr_bins[:, bin_index]
    for bin_index in range(sm.N_BINS):
        labels[f"zmet_bin_{bin_index:02d}"] = metallicities[:, bin_index]
    return pd.DataFrame(labels)


# =============================================================================
# 4. Simulated observations
# =============================================================================
def draw_observed_redshift(true_redshift, rng):
    """Observed redshift and its uncertainty: exact for a P_SPECTROSCOPIC share, photometric otherwise.

    Inputs: true_redshift (n,); rng: numpy Generator.
    Returns: (z_obs, sigma_z), each (n,); sigma_z = sm.SIGMA_Z * (1 + z) for photometric redshifts.
    """
    spectroscopic = rng.random(len(true_redshift)) < P_SPECTROSCOPIC
    sigma_z = np.where(spectroscopic, 0.0, sm.SIGMA_Z * (1 + true_redshift))
    z_obs = true_redshift + sigma_z * rng.normal(0, 1, len(true_redshift))
    return z_obs, sigma_z


def simulate_observations(table, catalogue_r, catalogue_z, noise_table, rng):
    """One fresh noisy realisation of every noise-free model: brightness, S-PLUS noise and redshift.

    Inputs: table: noise-free models (training_set.parquet columns); catalogue_r [mag], catalogue_z: DR3 r_auto
            and zml (brightness donors); noise_table: splus_noise_model.csv; rng: numpy Generator.
    Returns: dict(flat, banded, valid, targets) — inputs as colour_inputs gives them and the (n, 19) targets.
             valid is False where r is undetected or no DR3 galaxy lies at the model's redshift.
    """
    true_redshift = table["redshift"].to_numpy(dtype=float)
    noiseless = table[sm.BANDS].to_numpy(dtype=float)

    # 1. Realistic brightness, and the masses that go with it
    scaled, offsets = sm.scale_to_observed_r(noiseless, true_redshift, catalogue_r, catalogue_z, rng)
    labels = table.copy()
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED
    labels = sm.rescale_mass_labels(labels, offsets)
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * offsets

    # 2. S-PLUS noise and non-detections
    noisy, errors, _ = sm.add_splus_noise(scaled, noise_table, rng)

    # 3. Observed redshift
    z_obs, sigma_z = draw_observed_redshift(true_redshift, rng)

    observations = colour_inputs(noisy, errors, z_obs, sigma_z)
    observations["valid"] = observations["valid"] & np.isfinite(offsets)
    observations["targets"] = encode_targets(labels)
    return observations


def validation_observations(held_out, mode, rng):
    """Inputs and targets of a held-out set whose noise is already fixed (synthetic_validation/test).

    Inputs: held_out: DataFrame with noisy BANDS magnitudes, e_<band> errors, magnitude_offset and labels
            (mass labels already rescaled); mode: "spectroscopic" (sigma_z = 0) or "photometric"
            (sigma_z = sm.SIGMA_Z * (1 + z), z_obs drawn with rng); rng: numpy Generator (seeded by the caller).
    Returns: dict(flat, banded, valid, targets).
    """
    true_redshift = held_out["redshift"].to_numpy(dtype=float)
    if mode == "spectroscopic":
        sigma_z = np.zeros(len(true_redshift))
        z_obs = true_redshift.copy()
    elif mode == "photometric":
        sigma_z = sm.SIGMA_Z * (1 + true_redshift)
        z_obs = true_redshift + sigma_z * rng.normal(0, 1, len(true_redshift))
    else:
        raise ValueError(f"unknown mode {mode!r}: use 'spectroscopic' or 'photometric'")

    error_columns = []
    for band in sm.BANDS:
        error_columns.append(f"e_{band}")
    magnitudes = held_out[sm.BANDS].to_numpy(dtype=float)
    errors = held_out[error_columns].to_numpy(dtype=float)

    labels = held_out.copy()
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED - sm.DEX_PER_MAGNITUDE * held_out["magnitude_offset"]
    observations = colour_inputs(magnitudes, errors, z_obs, sigma_z)
    observations["valid"] = observations["valid"] & np.isfinite(held_out["magnitude_offset"].to_numpy())
    observations["targets"] = encode_targets(labels)
    return observations


# =============================================================================
# 5. Standardisation
# =============================================================================
class Standardiser:
    """Column-wise (x - mean) / std, fitted on training data and stored with the model.

    Entries whose mask is False are excluded from the fit and set to 0 after transforming
    (0 = 'no information' for a masked colour).
    """

    def __init__(self, mean, std):
        """Inputs: mean, std: arrays with one value per column (std > 0)."""
        self.mean = np.asarray(mean, dtype=float)
        self.std = np.asarray(std, dtype=float)

    @classmethod
    def fit(cls, values, mask=None):
        """Fit on values (n, ...) using only entries where mask is True; returns a Standardiser."""
        values = np.asarray(values, dtype=float)
        if mask is None:
            mask = np.ones(values.shape, dtype=bool)
        masked = np.where(mask, values, np.nan)
        mean = np.nanmean(masked, axis=0)
        std = np.nanstd(masked, axis=0)
        std = np.where(std > 0, std, 1.0)        # constant columns (e.g. log pivots) stay unscaled
        mean = np.nan_to_num(mean)
        return cls(mean, std)

    def transform(self, values, mask=None):
        """Standardised copy of values; masked entries become 0."""
        transformed = (np.asarray(values, dtype=float) - self.mean) / self.std
        if mask is not None:
            transformed = np.where(mask, transformed, 0.0)
        return transformed

    def inverse_transform(self, values):
        """Undo transform (for unmasked entries)."""
        return np.asarray(values, dtype=float) * self.std + self.mean

    def to_dict(self):
        """JSON-friendly dict of the statistics."""
        return {"mean": self.mean.tolist(), "std": self.std.tolist()}

    @classmethod
    def from_dict(cls, statistics):
        """Standardiser from to_dict output."""
        return cls(statistics["mean"], statistics["std"])


def input_masks(flat, banded):
    """Which entries of flat and banded carry information (for Standardiser fit/transform).

    Inputs: flat (n, 37), banded (n, 12, 4) as colour_inputs returns them.
    Returns: (flat_mask, banded_mask) boolean arrays of the same shapes: colours and colour errors follow
             their flags, flags and scalars are always informative; band channels 0-1 follow channel 2.
    """
    flat_mask = np.ones(flat.shape, dtype=bool)
    flags = flat[:, 2 * N_COLOURS:3 * N_COLOURS] > 0
    flat_mask[:, :N_COLOURS] = flags
    flat_mask[:, N_COLOURS:2 * N_COLOURS] = flags
    banded_mask = np.ones(banded.shape, dtype=bool)
    band_detected = banded[:, :, 2] > 0
    banded_mask[:, :, 0] = band_detected
    banded_mask[:, :, 1] = band_detected
    return flat_mask, banded_mask
