"""
Network inputs and targets for the SFH posterior estimator (Phase 3).

Model tables (training set, held-out sets, real galaxies) become:
- inputs: 11 colours (band - r) with errors and detection flags, plus r, its error, z_obs and sigma_z;
- targets: the parameters a flow models for one SFH model (sfh_models.MODELS: 9 for delayed-tau, 10 for the
  double power law): the SFH shape, metallicity, dust, nebular emission, masses and redshift; every other
  label is derived from them when decoding.
Every physical definition (prior, labels, noise, brightness) comes from sfh_models and synthetic_model.
"""

import joblib
import numpy as np
import pandas as pd

import sfh_models
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

# Flow targets that are the log of a drawn label but not among the model's derived labels (each model's
# log_labels give the others, e.g. log_age is the log of age)
EXTRA_LOG_TARGETS = {"log_tau_bc": "tau_bc"}


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
def target_names(sfh_model):
    """The flow targets of an SFH model, in the order the flow models them.

    Why these: parameters drawn log-uniform in the prior enter as their logs. log_sfr is not a target,
    because a fifth of the delayed-tau prior sits on its floor, a point mass a flow cannot learn. The
    redshift is last, because the decoder needs it first: the oldest allowed age or turnover lookback
    depends on it.
    Inputs: sfh_model: "delayed" or "dblplaw" (a key of sfh_models.MODELS).
    Returns: list of target names, a new list on every call.
    Raises: KeyError for an unknown model.
    """
    return list(sfh_models.MODELS[sfh_model].target_names)


def encode_targets(table, sfh_model):
    """The flow targets of each row.

    Inputs: table with the model's drawn labels (sfh_models.MODELS[sfh_model].drawn_labels: the redshift,
            the SFH parameters [Gyr or dimensionless], metallicity [Z/Z_sun], av_ism [mag], tau_bc, logU and
            log_mass_formed [log M_sun]) and log_mstar [log M_sun]; sfh_model: "delayed" or "dblplaw".
    Returns: array (n, len(target_names(sfh_model))) in target_names order. The log-uniform drawn
             parameters (the model's log_labels, and tau_bc) enter as log10; the other targets enter
             unchanged.
    """
    log_sources = dict(sfh_models.MODELS[sfh_model].log_labels)
    log_sources.update(EXTRA_LOG_TARGETS)
    columns = []
    for name in target_names(sfh_model):
        if name in log_sources:
            columns.append(np.log10(table[log_sources[name]].to_numpy(dtype=float)))
        else:
            columns.append(table[name].to_numpy(dtype=float))
    return np.column_stack(columns)


def delayed_sfh_parameters(targets, redshift):
    """Delayed-tau SFH parameters of flow samples, inside the prior at each sample's redshift.

    Inputs: targets: dict target name -> array (n,), with log_age and log_tau [log Gyr];
            redshift: array (n,), already inside sfh_models.Z_RANGE.
    Returns: dict with age clipped to [AGE_MIN, maximum_age(redshift)] and tau clipped to TAU_RANGE [Gyr].
             Values inside the prior pass unchanged.
    """
    age = np.clip(10 ** targets["log_age"], sfh_models.AGE_MIN, sfh_models.maximum_age(redshift))
    tau = np.clip(10 ** targets["log_tau"], *sfh_models.TAU_RANGE)
    return {"age": age, "tau": tau}


def dblplaw_sfh_parameters(targets, redshift):
    """Double-power-law SFH parameters of flow samples, inside the prior at each sample's redshift.

    Why tau is recomputed: Bagpipes' tau is the turnover time of the double power law, counted from the Big
    Bang (the SFR maximum only when alpha = beta; in general the maximum is at
    t = tau (beta/alpha)**(1/(alpha+beta))). turnover_lookback is the lookback time to tau, so
    tau = age of the Universe - turnover_lookback (as in sfh_models.dblplaw_components). The flow models
    turnover_lookback, and tau follows from it.
    Inputs: targets: dict target name -> array (n,), with log_turnover_lookback [log Gyr], log_alpha and
            log_beta;
            redshift: array (n,), already inside sfh_models.Z_RANGE.
    Returns: dict with
             - turnover_lookback [Gyr], clipped to [TURNOVER_LOOKBACK_MIN, maximum_age(redshift)];
             - tau [Gyr since the Big Bang];
             - alpha clipped to ALPHA_RANGE and beta clipped to BETA_RANGE (dimensionless).
             Values inside the prior pass unchanged.
    """
    turnover_lookback = np.clip(10 ** targets["log_turnover_lookback"], sfh_models.TURNOVER_LOOKBACK_MIN,
                                sfh_models.maximum_age(redshift))
    tau = sfh_models.COSMOLOGY.age(redshift).value - turnover_lookback
    alpha = np.clip(10 ** targets["log_alpha"], *sfh_models.ALPHA_RANGE)
    beta = np.clip(10 ** targets["log_beta"], *sfh_models.BETA_RANGE)
    return {"turnover_lookback": turnover_lookback, "tau": tau, "alpha": alpha, "beta": beta}


# The SFH parameters of each model, decoded from its flow targets
SFH_PARAMETER_DECODERS = {"delayed": delayed_sfh_parameters, "dblplaw": dblplaw_sfh_parameters}


def decode_samples(theta, sfh_model, n_jobs=1):
    """Every training label of flow samples, with the same definitions as the training labels.

    Each sample is decoded at its own redshift. The flow's tails can leave the prior, so:
    - the redshift is clipped to Z_RANGE;
    - the SFH parameters are clipped to their prior at that redshift (SFH_PARAMETER_DECODERS);
    - the metallicity is clipped to METALLICITY_RANGE.
    Values inside the prior pass unchanged, so encoded truths decode exactly.
    Why n_jobs: the double-power-law labels are summed over Bagpipes' age grid, about 0.06 ms per sample,
    so tens of millions of posterior samples take an hour or more in one thread. With n_jobs > 1 the rows
    are split into n_jobs blocks, decoded in threads and concatenated in order. Threads, not processes:
    numpy releases the GIL inside the grid sums, and threads share one copy of Bagpipes, which matters in a
    7 GB GPU job.
    Each row is decoded on its own, so the result equals that of n_jobs = 1 to rounding.
    Inputs: theta: array (n, len(target_names(sfh_model))) in target_names order (posterior samples or
            encoded truths); sfh_model: "delayed" or "dblplaw"; n_jobs: threads (1 decodes in this thread).
    Returns: DataFrame with the model's training_labels columns, in that order, one row per row of theta:
             - the drawn labels and their logs;
             - log_mstar [log M_sun];
             - from the model's analytic_labels: mw_age [Gyr ago], f_young, f_mid and f_old,
               log_sfr [log M_sun/yr] and log_ssfr [log yr^-1].
    Raises: ValueError when theta does not have one column per target of the model, or n_jobs < 1.
    """
    names = target_names(sfh_model)
    theta = np.asarray(theta, dtype=float)
    if theta.ndim != 2 or theta.shape[1] != len(names):
        raise ValueError(f"theta of shape {theta.shape} does not hold the {len(names)} targets of the "
                         f"{sfh_model!r} model")
    if n_jobs < 1:
        raise ValueError(f"n_jobs = {n_jobs}: decode with 1 or more threads")
    if n_jobs == 1 or len(theta) < 2:
        return decode_rows(theta, sfh_model)

    tasks = []
    for block in np.array_split(theta, min(n_jobs, len(theta))):
        tasks.append(joblib.delayed(decode_rows)(block, sfh_model))
    decoded_blocks = joblib.Parallel(n_jobs=n_jobs, prefer="threads")(tasks)
    return pd.concat(decoded_blocks, ignore_index=True)


def decode_rows(theta, sfh_model):
    """The decode of decode_samples, in this process and without its checks.

    Inputs: theta: array (n, len(target_names(sfh_model))) in target_names order; sfh_model: "delayed" or
            "dblplaw".
    Returns: DataFrame with the model's training_labels columns (see decode_samples).
    """
    model = sfh_models.MODELS[sfh_model]
    names = target_names(sfh_model)
    targets = {}
    for target_index, name in enumerate(names):
        targets[name] = theta[:, target_index]
    redshift = np.clip(targets["redshift"], *sfh_models.Z_RANGE)

    # 1. The drawn labels, inside the prior at each sample's redshift
    drawn = {"redshift": redshift}
    drawn.update(SFH_PARAMETER_DECODERS[sfh_model](targets, redshift))
    drawn["metallicity"] = np.clip(10 ** targets["log_metallicity"], *sfh_models.METALLICITY_RANGE)
    drawn["av_ism"] = targets["av_ism"]
    drawn["tau_bc"] = 10 ** targets["log_tau_bc"]
    drawn["logU"] = targets["logU"]
    drawn["log_mass_formed"] = targets["log_mass_formed"]

    # 2. Labels derived from them, as sm.summary_labels defines them
    labels = dict(drawn)
    for log_label, drawn_label in model.log_labels.items():
        labels[log_label] = np.log10(drawn[drawn_label])
    labels["log_mstar"] = targets["log_mstar"]
    labels.update(model.analytic_labels(drawn, targets["log_mstar"]))
    return pd.DataFrame(labels, columns=model.training_labels)


# delayed-tau aliases for ml_networks.DEFAULT_CONFIG, removed in Task 8
TARGET_NAMES = target_names("delayed")
N_TARGETS = len(TARGET_NAMES)


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


def simulate_observations(table, catalogue_r, catalogue_z, noise_table, rng, sfh_model):
    """One fresh noisy realisation of every noise-free model: brightness, S-PLUS noise and redshift.

    Inputs: table: noise-free models of sfh_model (its training-set columns);
            catalogue_r [mag], catalogue_z: DR3 r_auto and zml (the brightness donors);
            noise_table: splus_noise_model.csv; rng: numpy Generator; sfh_model: "delayed" or "dblplaw".
    Returns: dict(flat, banded, valid, targets):
             - flat, banded and valid: the inputs as colour_inputs gives them;
             - targets: array (n, len(target_names(sfh_model))), with the mass-like labels rescaled to the
               drawn brightness.
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
    observations["targets"] = encode_targets(labels, sfh_model)
    return observations

def validation_observations(held_out, mode, rng, sfh_model):
    """Inputs and targets of a held-out set whose noise is already fixed (synthetic_validation/test).

    Inputs: held_out: DataFrame of sfh_model with noisy BANDS magnitudes, e_<band> errors, magnitude_offset
            and labels (mass labels already rescaled);
            mode: "spectroscopic" (sigma_z = 0) or "photometric" (sigma_z = sm.SIGMA_Z * (1 + z), with z_obs
            drawn from rng);
            rng: numpy Generator (seeded by the caller); sfh_model: "delayed" or "dblplaw".
    Returns: dict(flat, banded, valid, targets), with targets (n, len(target_names(sfh_model))).
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
    observations["targets"] = encode_targets(labels, sfh_model)
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
