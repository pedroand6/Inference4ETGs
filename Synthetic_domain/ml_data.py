"""
Network inputs and targets for the SFH posterior estimator (Phase 3).

Model tables (training set, held-out sets, real galaxies) become:
- inputs: 11 colours (band - r) with errors and detection flags, plus r, its error, z_obs and sigma_z;
- targets: the 9 delayed-tau parameters a flow models (log age, log tau, log metallicity, dust, nebular,
  masses, redshift); every other label is derived from them when decoding.
Every physical definition (prior, analytic labels, noise, brightness) comes from synthetic_model.
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

P_SPECTROSCOPIC = 0.5                            # share of simulated galaxies with an exact redshift

# Flow targets. Age, tau, metallicity and tau_bc are log-uniform in the prior, so the flow models their logs.
# log_sfr is not a target: a fifth of the prior sits on its floor, a point mass a flow cannot learn.
# The redshift stays last (decoded first, since the oldest allowed age depends on it).
TARGET_NAMES = ["log_age", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU",
                "log_mass_formed", "log_mstar", "redshift"]
N_TARGETS = len(TARGET_NAMES)                    # 9


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
def encode_targets(table):
    """The N_TARGETS flow targets of each row.

    Inputs: table with redshift, age and tau [Gyr], metallicity [Z/Z_sun], av_ism [mag], tau_bc, logU,
            log_mass_formed and log_mstar [log M_sun] (the sm.TRAINING_LABELS names).
    Returns: array (n, N_TARGETS) in TARGET_NAMES order: log10 of age, tau, metallicity and tau_bc;
             the other targets as they are.
    """
    targets = {
        "log_age": np.log10(table["age"].to_numpy(dtype=float)),
        "log_tau": np.log10(table["tau"].to_numpy(dtype=float)),
        "log_metallicity": np.log10(table["metallicity"].to_numpy(dtype=float)),
        "av_ism": table["av_ism"].to_numpy(dtype=float),
        "log_tau_bc": np.log10(table["tau_bc"].to_numpy(dtype=float)),
        "logU": table["logU"].to_numpy(dtype=float),
        "log_mass_formed": table["log_mass_formed"].to_numpy(dtype=float),
        "log_mstar": table["log_mstar"].to_numpy(dtype=float),
        "redshift": table["redshift"].to_numpy(dtype=float),
    }
    columns = []
    for name in TARGET_NAMES:
        columns.append(targets[name])
    return np.column_stack(columns)


def decode_samples(theta):
    """Every training label of flow samples, with the same analytic definitions as the training labels.

    Each sample is decoded at its own redshift. The flow's tails can leave the prior, so the redshift is
    clipped to sm.Z_RANGE, then age to [sm.AGE_MIN, sm.maximum_age(redshift)] and tau and metallicity to
    their prior ranges; values inside the prior pass unchanged, so encoded truths decode exactly.
    Inputs: theta: array (n, N_TARGETS) in TARGET_NAMES order (posterior samples or encoded truths).
    Returns: DataFrame with the sm.TRAINING_LABELS columns, in that order: redshift, age and tau [Gyr],
             metallicity [Z/Z_sun], av_ism [mag], tau_bc, logU, log_mass_formed, their logs, log_mstar
             [log M_sun], and log_sfr, log_ssfr, mw_age and f_young/f_mid/f_old from sm.delayed_tau_labels.
    """
    theta = np.asarray(theta, dtype=float)
    redshift = np.clip(theta[:, TARGET_NAMES.index("redshift")], *sm.Z_RANGE)

    # 1. The SFH parameters, inside the prior at each sample's redshift
    age = np.clip(10 ** theta[:, TARGET_NAMES.index("log_age")], sm.AGE_MIN, sm.maximum_age(redshift))
    tau = np.clip(10 ** theta[:, TARGET_NAMES.index("log_tau")], *sm.TAU_RANGE)
    metallicity = np.clip(10 ** theta[:, TARGET_NAMES.index("log_metallicity")], *sm.METALLICITY_RANGE)
    log_mass_formed = theta[:, TARGET_NAMES.index("log_mass_formed")]
    log_mstar = theta[:, TARGET_NAMES.index("log_mstar")]

    # 2. Labels derived from them, as sm.summary_labels defines them
    labels = {
        "redshift": redshift,
        "age": age,
        "tau": tau,
        "metallicity": metallicity,
        "av_ism": theta[:, TARGET_NAMES.index("av_ism")],
        "tau_bc": 10 ** theta[:, TARGET_NAMES.index("log_tau_bc")],
        "logU": theta[:, TARGET_NAMES.index("logU")],
        "log_mass_formed": log_mass_formed,
        "log_age": np.log10(age),
        "log_tau": np.log10(tau),
        "log_metallicity": np.log10(metallicity),
        "log_mstar": log_mstar,
    }
    labels.update(sm.delayed_tau_labels(age, tau, log_mass_formed, log_mstar))
    return pd.DataFrame(labels, columns=sm.TRAINING_LABELS)


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
    Returns: dict(flat, banded, valid, targets) — inputs as colour_inputs gives them and the targets
             (n, N_TARGETS), with the mass-like labels rescaled to the drawn brightness.
             valid is False where r is undetected or no DR3 galaxy lies at the model's redshift.
    """
    true_redshift = table["redshift"].to_numpy(dtype=float)
    noiseless = table[sm.BANDS].to_numpy(dtype=float)

    # 1. Realistic brightness, and the masses that go with it
    scaled, offsets = sm.scale_to_observed_r(noiseless, true_redshift, catalogue_r, catalogue_z, rng)
    labels = table.copy()
    labels["log_mass_formed"] = sm.LOG_MASS_FORMED     # the noise-free models' value, rescaled below
    labels = sm.rescale_mass_labels(labels, offsets)

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
