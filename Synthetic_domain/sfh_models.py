"""
SFH models of the S-PLUS synthetic domain: the shared prior, one section per Bagpipes SFH, and the registry.

synthetic_model.py imports this module (never the reverse), so the training sets, the fits and the ML use
exactly the same prior and labels for each SFH model.
"""

import os
from dataclasses import dataclass
from typing import Callable

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.special import gammainc, gammaincc

import bagpipes as pipes


# =============================================================================
# 1. Shared prior (fixed before looking at the data)
# =============================================================================
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3)
UNIVERSE_AGE_MARGIN = 0.1        # [Gyr] age upper bound = Universe age at z minus this
LABEL_EDGES = [0.1, 1.0]         # [Gyr] f_young < 0.1 < f_mid < 1.0 < f_old
SFR_TIMESCALE = 0.1              # [Gyr] log_sfr is the mean SFR over this lookback time, as Bagpipes' sfr
GYR_TO_YEARS = 1e9               # [yr/Gyr]

# Parameters every SFH model shares: redshift, one metallicity, dust, nebular emission and the formed mass
Z_RANGE = (0.002, 0.7)                               # redshift
METALLICITY_RANGE = (0.005, 2.5)                     # [Z/Z_sun] the bc03_miles grid range
AV_ISM_RANGE = (0.01, 2.0)                           # [mag]
LOG_TAU_BC_RANGE = (np.log10(0.5), np.log10(4.0))    # birth-cloud optical depth
LOGU_RANGE = (-4.0, -1.0)                            # nebular ionisation parameter
DUST_SLOPE = 0.7                                     # CF00 attenuation-law exponent
BIRTH_CLOUD_AGE = 0.01                               # [Gyr]
LOG_MASS_FORMED = 10.0                               # [log M_sun]
TAU_TO_AV = 1.086                                    # A_V = 1.086 * tau_V
SFR_FLOOR = 1e-14                # [yr^-1] sSFR floor for SFHs with no recent star formation

# Labels every model's analytic_labels returns, computed from the shape of its SFH
HISTORY_LABELS = ["mw_age", "f_young", "f_mid", "f_old", "log_sfr", "log_ssfr"]


def uniform(draw, low, high):
    """Map a unit-interval draw onto [low, high).

    Inputs: draw in [0, 1) (dimensionless); low, high: bounds in the units of the parameter.
    Returns: parameter value in [low, high), in the same units as low and high.
    """
    return low + (high - low) * draw


def maximum_age(redshift):
    """Oldest allowed age [Gyr]: the age of the Universe at the redshift minus UNIVERSE_AGE_MARGIN.

    Why the margin: Bagpipes flags a model forming stars before the Big Bang as unphysical (NaN magnitudes).
    Inputs: redshift: dimensionless z, scalar or array.
    Returns: maximum age [Gyr], a float for a scalar redshift or an array of the same shape.
    """
    return COSMOLOGY.age(redshift).value - UNIVERSE_AGE_MARGIN


def dust_and_nebular(params):
    """Bagpipes dust, birth-cloud and nebular entries of one model, the same for every SFH model.

    Inputs: params: dict of unit-cube values (dimensionless, in [0, 1)) holding "av_ism", "log_tau_bc" and "logU".
    Returns: (bagpipes_parts, drawn). bagpipes_parts holds "t_bc" [Gyr], "dust" (CF00) and "nebular";
             drawn holds av_ism [mag], tau_bc and logU (dimensionless).
    """
    av_ism = uniform(params["av_ism"], *AV_ISM_RANGE)
    tau_birth_cloud = 10 ** uniform(params["log_tau_bc"], *LOG_TAU_BC_RANGE)
    av_birth_cloud = TAU_TO_AV * tau_birth_cloud
    log_u = uniform(params["logU"], *LOGU_RANGE)

    bagpipes_parts = {
        "t_bc": BIRTH_CLOUD_AGE,
        "dust": {"type": "CF00", "Av": av_ism, "eta": 1 + av_birth_cloud / av_ism, "n": DUST_SLOPE},
        "nebular": {"logU": log_u},
    }
    drawn = {"av_ism": av_ism, "tau_bc": tau_birth_cloud, "logU": log_u}
    return bagpipes_parts, drawn


# =============================================================================
# 2. Delayed-tau SFH: SFR(t) ~ t exp(-t / tau) since star formation began
# =============================================================================
# One Bagpipes "delayed" component with a single metallicity. Ages and tau are drawn log-uniform,
# because they span several decades.
AGE_MIN = 0.05                   # [Gyr] youngest age (time since star formation began)
TAU_RANGE = (0.01, 20.0)         # [Gyr] e-folding time

# "log_age_fraction" places the age inside [log AGE_MIN, log t_max(z)], so one cube works at every z
DELAYED_DIMS = ["z", "log_age_fraction", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU"]


def delayed_components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one delayed-tau unit-cube vector u (order of DELAYED_DIMS).

    Inputs: u: unit-cube vector (dimensionless, in [0, 1)) ordered as DELAYED_DIMS;
            fixed_redshift: if given, replaces the redshift drawn from u (fits); dimensionless z.
    Returns: (bagpipes_components, drawn). drawn holds the delayed-tau drawn labels: redshift (dimensionless),
             age and tau [Gyr], metallicity [Z/Z_sun], av_ism [mag], tau_bc and logU (dimensionless),
             log_mass_formed [log M_sun].
    """
    params = dict(zip(DELAYED_DIMS, u))

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
    dust_parts, dust_drawn = dust_and_nebular(params)

    bagpipes_components = {"redshift": redshift}
    bagpipes_components.update(dust_parts)
    bagpipes_components["delayed"] = {"age": age, "tau": tau, "massformed": LOG_MASS_FORMED,
                                      "metallicity": metallicity}
    drawn = {"redshift": redshift, "age": age, "tau": tau, "metallicity": metallicity,
             "av_ism": dust_drawn["av_ism"], "tau_bc": dust_drawn["tau_bc"], "logU": dust_drawn["logU"],
             "log_mass_formed": LOG_MASS_FORMED}
    return bagpipes_components, drawn


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


def delayed_analytic_labels(drawn_columns, log_mstar):
    """The HISTORY_LABELS of delayed-tau models from their drawn labels (the registry's analytic_labels).

    Inputs: drawn_columns: dict of arrays holding at least "age", "tau" [Gyr] and "log_mass_formed"
            [log M_sun]; log_mstar: array [log M_sun] of the same shape.
    Returns: dict of arrays, as delayed_tau_labels.
    """
    return delayed_tau_labels(drawn_columns["age"], drawn_columns["tau"], drawn_columns["log_mass_formed"],
                              log_mstar)


# =============================================================================
# 3. Double power law: SFR(t) ~ [(t / tau)**alpha + (t / tau)**-beta]**-1, t since the Big Bang
# =============================================================================
# One Bagpipes "dblplaw" component with a single metallicity. Star formation starts at the Big Bang, so there is
# no age. tau is the turnover time since the Big Bang, alpha the falling slope after it and beta the rising slope
# before it. The SFR maximum lies at tau only when alpha = beta; in general it lies at t = tau (beta / alpha)**(1 /
# (alpha + beta)). turnover_lookback, the lookback time to tau, is drawn as the delayed-tau age, log-uniform in
# lookback time; alpha and beta are log-uniform (Carnall et al. 2018).
TURNOVER_LOOKBACK_MIN = 0.05     # [Gyr] most recent turnover, as AGE_MIN
ALPHA_RANGE = (0.01, 1000.0)     # falling slope after the turnover (dimensionless)
BETA_RANGE = (0.01, 1000.0)      # rising slope before the turnover (dimensionless)

# "log_turnover_lookback_fraction" places the turnover inside [log TURNOVER_LOOKBACK_MIN, log t_max(z)], as the
# delayed-tau age
DBLPLAW_DIMS = ["z", "log_turnover_lookback_fraction", "log_alpha", "log_beta", "log_metallicity", "av_ism",
                "log_tau_bc", "logU"]

# Bagpipes' own SFH age grid (star_formation_history.__init__, default log_sampling)
BAGPIPES_LOG_AGE_START = 6.0     # [log yr] youngest grid age
BAGPIPES_LOG_AGE_STEP = 0.0025   # [dex] grid step
LOG_GYR_TO_YEARS = 9.0           # [dex] log10 of GYR_TO_YEARS, which Bagpipes adds to the log Hubble time
DBLPLAW_CHUNK = 32               # samples per block of the grid sum: 32 x 1654 grid ages, 0.42 MB per array


def dblplaw_components(u, fixed_redshift=None):
    """Bagpipes components and drawn parameters for one DPL unit-cube vector u (order of DBLPLAW_DIMS).

    Why tau >= UNIVERSE_AGE_MARGIN: turnover_lookback is at most maximum_age(z), so the turnover always lies after
    the Big Bang, and tau stays below Bagpipes' "unphysical" limit of 1.5 times the age of the Universe.
    Inputs: u: unit-cube vector (dimensionless, in [0, 1)) ordered as DBLPLAW_DIMS;
            fixed_redshift: if given, replaces the redshift drawn from u (fits); dimensionless z.
    Returns: (bagpipes_components, drawn). drawn holds the DPL drawn labels: redshift (dimensionless),
             turnover_lookback [Gyr ago], tau [Gyr since the Big Bang], alpha and beta (dimensionless),
             metallicity [Z/Z_sun], av_ism [mag], tau_bc and logU (dimensionless), log_mass_formed [log M_sun].
    """
    params = dict(zip(DBLPLAW_DIMS, u))

    # 1. Redshift, which sets the oldest allowed turnover
    if fixed_redshift is None:
        redshift = uniform(params["z"], *Z_RANGE)
    else:
        redshift = fixed_redshift
    oldest_turnover = maximum_age(redshift)

    # 2. Star formation history and metallicity
    turnover_lookback = 10 ** uniform(params["log_turnover_lookback_fraction"], np.log10(TURNOVER_LOOKBACK_MIN),
                                      np.log10(oldest_turnover))
    # Rounding in 10**log10 can push the top draw a hair above the bound
    turnover_lookback = min(turnover_lookback, oldest_turnover)
    tau = COSMOLOGY.age(redshift).value - turnover_lookback
    alpha = 10 ** uniform(params["log_alpha"], *np.log10(ALPHA_RANGE))
    beta = 10 ** uniform(params["log_beta"], *np.log10(BETA_RANGE))
    metallicity = 10 ** uniform(params["log_metallicity"], *np.log10(METALLICITY_RANGE))

    # 3. Dust and nebular emission
    dust_parts, dust_drawn = dust_and_nebular(params)

    bagpipes_components = {"redshift": redshift}
    bagpipes_components.update(dust_parts)
    bagpipes_components["dblplaw"] = {"tau": tau, "alpha": alpha, "beta": beta, "massformed": LOG_MASS_FORMED,
                                      "metallicity": metallicity}
    drawn = {"redshift": redshift, "turnover_lookback": turnover_lookback, "tau": tau, "alpha": alpha,
             "beta": beta, "metallicity": metallicity, "av_ism": dust_drawn["av_ism"],
             "tau_bc": dust_drawn["tau_bc"], "logU": dust_drawn["logU"], "log_mass_formed": LOG_MASS_FORMED}
    return bagpipes_components, drawn


def bagpipes_age_grid():
    """Bagpipes' SFH age grid, rebuilt exactly as star_formation_history.__init__ builds it.

    Why: the DPL labels have no closed form. Summed on Bagpipes' own grid they match its sfh arrays, including
    its discretisation of sharp peaks (alpha or beta near 1000).
    Returns: (ages, widths) [Gyr], two arrays of the grid length: the lookback age of each grid point and the
             width of its bin. The first bin starts at 0 and the last ends at the Hubble time.
    """
    hubble_time = pipes.utils.age_at_z[pipes.utils.z_array == 0.0][0]     # [Gyr]
    log_age_max = np.log10(hubble_time) + LOG_GYR_TO_YEARS + 2 * BAGPIPES_LOG_AGE_STEP
    log_ages = np.arange(BAGPIPES_LOG_AGE_START, log_age_max, BAGPIPES_LOG_AGE_STEP)    # [log yr]
    log_edges = pipes.utils.make_bins(log_ages, make_rhs=True)[0]                       # [log yr]
    edges = 10**log_edges                                                                # [yr]
    edges[0] = 0.0
    edges[-1] = GYR_TO_YEARS * hubble_time
    ages = 10**log_ages / GYR_TO_YEARS
    widths = (edges[1:] - edges[:-1]) / GYR_TO_YEARS
    return ages, widths


BAGPIPES_AGES, BAGPIPES_AGE_WIDTHS = bagpipes_age_grid()      # [Gyr]

# The grid ages increase, so each label is a sum over one slice of the grid. Each index is the first grid age at
# or above the edge, which matches Bagpipes' "ages < edge"
YOUNG_EDGE_INDEX = int(np.searchsorted(BAGPIPES_AGES, LABEL_EDGES[0]))     # f_young: grid ages below 0.1 Gyr
OLD_EDGE_INDEX = int(np.searchsorted(BAGPIPES_AGES, LABEL_EDGES[1]))       # f_old: grid ages from 1 Gyr on
# log_sfr: grid ages below SFR_TIMESCALE
RECENT_EDGE_INDEX = int(np.searchsorted(BAGPIPES_AGES, SFR_TIMESCALE))


def bagpipes_universe_age(redshift):
    """Age of the Universe [Gyr] at the redshift, interpolated in Bagpipes' redshift table as Bagpipes does.

    Why not COSMOLOGY.age: the cosmology is the same, but Bagpipes interpolates; its value decides which grid
    ages lie after the Big Bang, exactly as in Bagpipes.
    Inputs: redshift (dimensionless), scalar or array.
    Returns: age [Gyr], same shape.
    """
    return np.interp(redshift, pipes.utils.z_array, pipes.utils.age_at_z)


def dblplaw_grid_sums(tau, alpha, beta, redshift):
    """Mass-weighted age and formed-mass fractions of double-power-law SFHs, summed on Bagpipes' age grid.

    The SFR is evaluated at every grid age younger than the Universe at the redshift (t = age of the Universe
    - lookback age), as Bagpipes' dblplaw does. Samples go DBLPLAW_CHUNK at a time, vectorised over the block,
    so the arrays stay in the CPU cache. A block uses only the grid ages below the age of the Universe at its
    lowest redshift, because older ones hold no stars.
    Why in logs: alpha and beta reach 1000, so (t / tau)**alpha overflows. With x = ln(t / tau), the exponents
    are alpha x and -beta x, the larger of the two is L >= 0, and
    SFR = 1 / (e^(alpha x) + e^(-beta x)) = e^-L / (1 + e^-((alpha + beta) |x|)).
    e^-L is divided by its largest value in each sample's row, which leaves the normalised SFH unchanged and
    keeps the sum finite when the SFR is a spike between two grid ages.
    Why no matrix product: a threaded BLAS product of these small blocks is slower than numpy's sums.
    Inputs: tau [Gyr since the Big Bang], alpha, beta, redshift (dimensionless): 1-D arrays of one length.
    Returns: dict of 1-D arrays: "mw_age" [Gyr ago]; "f_young", "f_mid", "f_old" (formed-mass fractions split
             at LABEL_EDGES, sum to 1); "recent" (formed-mass fraction younger than SFR_TIMESCALE).
    """
    n_samples = len(tau)
    sums = {}
    for key in ("mw_age", "f_young", "f_mid", "f_old", "recent"):
        sums[key] = np.zeros(n_samples)

    for start in range(0, n_samples, DBLPLAW_CHUNK):
        block = slice(start, start + DBLPLAW_CHUNK)
        block_tau = tau[block, np.newaxis]             # [Gyr], (block, 1)
        block_alpha = alpha[block, np.newaxis]
        block_beta = beta[block, np.newaxis]
        universe_age = bagpipes_universe_age(redshift[block])      # [Gyr], (block,)

        # 1. Grid ages that the oldest Universe of the block reaches
        n_ages = int(np.searchsorted(BAGPIPES_AGES, universe_age.max()))
        ages = BAGPIPES_AGES[:n_ages]
        widths = BAGPIPES_AGE_WIDTHS[:n_ages]

        # 2. x = ln(t / tau), with t the time since the Big Bang of each grid age, (block, n_ages). Grid ages
        #    before the Big Bang get t = tau (x = 0), so the log stays finite; step 4 gives them no mass
        log_time_ratio = np.subtract.outer(universe_age, ages)     # t [Gyr]
        before_big_bang = log_time_ratio <= 0
        np.copyto(log_time_ratio, block_tau, where=before_big_bang)
        log_time_ratio /= block_tau
        np.log(log_time_ratio, out=log_time_ratio)

        # 3. The two exponents of SFR = 1 / (e^(alpha x) + e^(-beta x))
        falling_exponent = block_alpha * log_time_ratio
        rising_exponent = -block_beta * log_time_ratio
        larger_exponent = np.maximum(falling_exponent, rising_exponent)          # L >= 0
        exponent_gap = np.abs(falling_exponent - rising_exponent)                # (alpha + beta) |x|

        # 4. Formed mass per grid bin, relative to the largest grid SFR of each sample
        larger_exponent[before_big_bang] = np.inf
        smallest_exponent = larger_exponent.min(axis=1, keepdims=True)
        formed_mass = np.exp(smallest_exponent - larger_exponent) / (1 + np.exp(-exponent_gap)) * widths

        # 5. The sums, as Bagpipes' mass_weighted_age and the label edges
        total_mass = formed_mass.sum(axis=1)
        sums["mw_age"][block] = (formed_mass * ages).sum(axis=1) / total_mass
        sums["f_young"][block] = formed_mass[:, :YOUNG_EDGE_INDEX].sum(axis=1) / total_mass
        sums["f_mid"][block] = formed_mass[:, YOUNG_EDGE_INDEX:OLD_EDGE_INDEX].sum(axis=1) / total_mass
        sums["f_old"][block] = formed_mass[:, OLD_EDGE_INDEX:].sum(axis=1) / total_mass
        sums["recent"][block] = formed_mass[:, :RECENT_EDGE_INDEX].sum(axis=1) / total_mass
    return sums


def dblplaw_labels(tau, alpha, beta, redshift, log_mass_formed, log_mstar):
    """Labels of a double-power-law SFH, shared by the training sets and the ML decoding.

    mw_age and the mass fractions have no closed form: they are summed on Bagpipes' age grid (dblplaw_grid_sums).
    log_sfr is the formed mass younger than SFR_TIMESCALE over SFR_TIMESCALE, floored at SFR_FLOOR * M_star,
    as for delayed-tau.
    Inputs: tau [Gyr since the Big Bang], alpha, beta, redshift (dimensionless), log_mass_formed and log_mstar
            [log M_sun]; arrays that broadcast together, or scalars.
    Returns: dict of arrays of the broadcast shape, the keys of delayed_tau_labels: mw_age [Gyr ago], f_young,
             f_mid, f_old (dimensionless, sum to 1), log_sfr [log M_sun/yr] and log_ssfr [log yr^-1].
    """
    tau, alpha, beta, redshift, log_mass_formed, log_mstar = np.broadcast_arrays(
        np.asarray(tau, dtype=float), np.asarray(alpha, dtype=float), np.asarray(beta, dtype=float),
        np.asarray(redshift, dtype=float), np.asarray(log_mass_formed, dtype=float),
        np.asarray(log_mstar, dtype=float))
    shape = tau.shape
    sums = dblplaw_grid_sums(tau.ravel(), alpha.ravel(), beta.ravel(), redshift.ravel())

    recent_sfr = 10**log_mass_formed.ravel() * sums["recent"] / (SFR_TIMESCALE * GYR_TO_YEARS)
    sfr_floor = SFR_FLOOR * 10**log_mstar.ravel()
    log_sfr = np.log10(np.maximum(recent_sfr, sfr_floor))
    labels = {
        "mw_age": sums["mw_age"],
        "f_young": sums["f_young"],
        "f_mid": sums["f_mid"],
        "f_old": sums["f_old"],
        "log_sfr": log_sfr,
        "log_ssfr": log_sfr - log_mstar.ravel(),
    }
    for label in labels:
        labels[label] = labels[label].reshape(shape)
    return labels


def dblplaw_analytic_labels(drawn_columns, log_mstar):
    """The HISTORY_LABELS of double-power-law models from their drawn labels (the registry's analytic_labels).

    Inputs: drawn_columns: dict of arrays holding at least "tau" [Gyr since the Big Bang], "alpha", "beta",
            "redshift" (dimensionless) and "log_mass_formed" [log M_sun]; log_mstar: array [log M_sun].
    Returns: dict of arrays, as dblplaw_labels.
    """
    return dblplaw_labels(drawn_columns["tau"], drawn_columns["alpha"], drawn_columns["beta"],
                          drawn_columns["redshift"], drawn_columns["log_mass_formed"], log_mstar)


# =============================================================================
# 4. The model registry: everything that differs between SFH models
# =============================================================================
# Labels with the same meaning under every SFH model, compared between the models
SHARED_LABELS = ["log_mstar", "mw_age", "log_metallicity", "av_ism", "log_ssfr", "f_young", "f_mid", "f_old"]


@dataclass(frozen=True)
class SfhModel:
    """One SFH model: its unit cube, Bagpipes components, labels, flow targets and file names."""

    name: str                    # "delayed" or "dblplaw", stored in every job, fit and checkpoint
    dims: list                   # unit-cube dimensions, "z" first
    components: Callable         # components(u, fixed_redshift=None) -> (bagpipes_components, drawn)
    drawn_labels: list           # keys of drawn, in this order
    derived_labels: list         # labels computed from drawn and from Bagpipes' stellar mass
    log_labels: dict             # derived log label -> drawn label, e.g. {"log_age": "age"}
    analytic_labels: Callable    # analytic_labels(drawn_columns, log_mstar) -> dict with the HISTORY_LABELS
    target_names: list           # flow targets, redshift last
    calibration_labels: list     # labels of the calibration (coverage) test of the photometric fits
    file_suffix: str             # "" for delayed-tau (today's file names), "_dblplaw" for the double power law

    @property
    def training_labels(self):
        """Label columns of the sets and the fits: the drawn labels, then the derived ones. Returns: list."""
        return self.drawn_labels + self.derived_labels

    @property
    def photometric_dims(self):
        """Unit-cube dimensions of a fit at fixed redshift: dims without "z". Returns: list."""
        photometric = []
        for dimension in self.dims:
            if dimension != "z":
                photometric.append(dimension)
        return photometric


DELAYED = SfhModel(
    name="delayed",
    dims=DELAYED_DIMS,
    components=delayed_components,
    drawn_labels=["redshift", "age", "tau", "metallicity", "av_ism", "tau_bc", "logU", "log_mass_formed"],
    derived_labels=["log_age", "log_tau", "log_metallicity", "log_mstar", "log_sfr", "log_ssfr",
                    "mw_age", "f_young", "f_mid", "f_old"],
    log_labels={"log_age": "age", "log_tau": "tau", "log_metallicity": "metallicity"},
    analytic_labels=delayed_analytic_labels,
    target_names=["log_age", "log_tau", "log_metallicity", "av_ism", "log_tau_bc", "logU", "log_mass_formed",
                  "log_mstar", "redshift"],
    calibration_labels=["log_age", "log_tau", "log_metallicity", "log_mstar"],
    file_suffix="",
)

DBLPLAW = SfhModel(
    name="dblplaw",
    dims=DBLPLAW_DIMS,
    components=dblplaw_components,
    drawn_labels=["redshift", "turnover_lookback", "tau", "alpha", "beta", "metallicity", "av_ism", "tau_bc", "logU",
                  "log_mass_formed"],
    derived_labels=["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "log_mstar", "log_sfr",
                    "log_ssfr", "mw_age", "f_young", "f_mid", "f_old"],
    log_labels={"log_turnover_lookback": "turnover_lookback", "log_alpha": "alpha", "log_beta": "beta",
                "log_metallicity": "metallicity"},
    analytic_labels=dblplaw_analytic_labels,
    target_names=["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "av_ism", "log_tau_bc", "logU",
                  "log_mass_formed", "log_mstar", "redshift"],
    calibration_labels=["log_turnover_lookback", "log_alpha", "log_beta", "log_metallicity", "log_mstar"],
    file_suffix="_dblplaw",
)

MODELS = {"delayed": DELAYED, "dblplaw": DBLPLAW}


def model_file(name, sfh_model):
    """File or folder name of one SFH model's product: the delayed-tau name with the model's suffix.

    Inputs: name: the delayed-tau file or folder name ("training_set.parquet", "photometric_fits", "ml_runs/mlp");
            sfh_model: a name in MODELS.
    Returns: the name with the suffix before a file's extension, or after a folder name
             ("training_set_dblplaw.parquet", "photometric_fits_dblplaw", "ml_runs/mlp_dblplaw");
             unchanged for delayed-tau.
    """
    suffix = MODELS[sfh_model].file_suffix
    stem, extension = os.path.splitext(name)
    return stem + suffix + extension
