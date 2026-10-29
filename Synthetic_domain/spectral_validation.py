"""
Rulers for the Bagpipes spectral fits that do not rely on full-spectrum SED fitting:
model-independent spectral indices (compared with the fits' posterior predictions), light-weighted ages and
metallicities (for Gallazzi et al. 2005), the Speagle et al. (2014) main sequence and percentile bands of
literature scaling relations.
"""

import functools

import bagpipes as pipes
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.stats import norm

import synthetic_model as sm


# =============================================================================
# 1. Index definitions (air wavelengths, rest frame, Angstrom)
# =============================================================================
D4000_BLUE_BAND = (3850.0, 3950.0)      # Balogh et al. (1999), narrow 4000 A break
D4000_RED_BAND = (4000.0, 4100.0)
EQUIVALENT_WIDTH_INDICES = {
    "Hdelta_A": {"blue": (4041.600, 4079.750), "feature": (4083.500, 4122.250), "red": (4128.500, 4161.000)},
    "Mgb": {"blue": (5142.625, 5161.375), "feature": (5160.125, 5192.625), "red": (5191.375, 5206.375)},
    "Fe5270": {"blue": (5233.150, 5248.150), "feature": (5245.650, 5285.650), "red": (5285.650, 5318.150)},
    "Fe5335": {"blue": (5304.625, 5315.875), "feature": (5312.125, 5352.125), "red": (5353.375, 5363.375)},
}
MGFE_FE5270_WEIGHT = 0.72               # Thomas et al. (2003) [MgFe]' weights
MGFE_FE5335_WEIGHT = 0.28
INDEX_NAMES = ["D4000_n", "Hdelta_A", "Mgb", "Fe5270", "Fe5335", "MgFe_prime"]
MGFE_COMPONENTS = ["Mgb", "Fe5270", "Fe5335"]
N_NOISE_DRAWS = 200                     # Monte-Carlo draws for observed-index errors


def in_band(wavelength, band):
    """True for pixels whose centre lies inside band = (low, high) [Angstrom]."""
    return (wavelength >= band[0]) & (wavelength <= band[1])


def band_mean(wavelength, values, band):
    """Mean of values over the pixels inside band; NaN when the band holds no pixel."""
    inside = in_band(wavelength, band)
    if not inside.any():
        return np.nan
    return np.mean(values[inside])


def d4000_n(wavelength, flux):
    """Narrow 4000 A break: mean F_nu in 4000-4100 A over mean F_nu in 3850-3950 A (F_nu ∝ F_lambda lambda^2).

    Inputs: rest-frame air wavelength [Angstrom]; flux F_lambda (any units). Returns: the ratio (NaN if a band is empty).
    """
    flux_per_frequency = flux * wavelength**2
    red = band_mean(wavelength, flux_per_frequency, D4000_RED_BAND)
    blue = band_mean(wavelength, flux_per_frequency, D4000_BLUE_BAND)
    return red / blue


def equivalent_width(wavelength, flux, definition):
    """Equivalent width [Angstrom] against a straight pseudo-continuum through the blue and red band means.

    Inputs: rest-frame air wavelength [Angstrom]; flux F_lambda; definition: dict with blue, feature, red bands.
    Returns: sum over feature pixels of (1 - F / continuum) * pixel width; NaN if a band holds no pixel.
    """
    blue_flux = band_mean(wavelength, flux, definition["blue"])
    red_flux = band_mean(wavelength, flux, definition["red"])
    feature = in_band(wavelength, definition["feature"])
    if np.isnan(blue_flux) or np.isnan(red_flux) or not feature.any():
        return np.nan
    blue_centre = np.mean(definition["blue"])
    red_centre = np.mean(definition["red"])
    slope = (red_flux - blue_flux) / (red_centre - blue_centre)
    continuum = blue_flux + slope * (wavelength[feature] - blue_centre)
    pixel_widths = np.gradient(wavelength)[feature]
    return np.sum((1 - flux[feature] / continuum) * pixel_widths)


def mgfe_prime(mgb, fe5270, fe5335):
    """[MgFe]' = sqrt(Mgb (0.72 Fe5270 + 0.28 Fe5335)) [Angstrom]; NaN when the argument is negative or NaN."""
    argument = mgb * (MGFE_FE5270_WEIGHT * fe5270 + MGFE_FE5335_WEIGHT * fe5335)
    if not np.isfinite(argument) or argument < 0:
        return np.nan
    return np.sqrt(argument)


def index_values(wavelength, flux):
    """All indices of one spectrum. Inputs: rest-frame air wavelength [Angstrom], flux F_lambda.

    Returns: dict name -> value for INDEX_NAMES (D4000_n dimensionless, the others in Angstrom).
    """
    values = {"D4000_n": d4000_n(wavelength, flux)}
    for name, definition in EQUIVALENT_WIDTH_INDICES.items():
        values[name] = equivalent_width(wavelength, flux, definition)
    values["MgFe_prime"] = mgfe_prime(values["Mgb"], values["Fe5270"], values["Fe5335"])
    return values


def index_bands(name):
    """The wavelength bands an index uses."""
    if name == "D4000_n":
        return [D4000_BLUE_BAND, D4000_RED_BAND]
    return list(EQUIVALENT_WIDTH_INDICES[name].values())


def index_validity(wavelength, ivar):
    """Whether each index can be measured: every band has pixels and none of them is masked (ivar = 0).

    Returns: dict name -> bool for INDEX_NAMES ([MgFe]' needs its three components).
    """
    validity = {}
    for name in INDEX_NAMES:
        if name == "MgFe_prime":
            continue
        valid = True
        for band in index_bands(name):
            inside = in_band(wavelength, band)
            if not inside.any() or np.any(ivar[inside] <= 0):
                valid = False
        validity[name] = valid
    components_valid = []
    for component in MGFE_COMPONENTS:
        components_valid.append(validity[component])
    validity["MgFe_prime"] = all(components_valid)
    return validity


def measure_indices(wavelength, flux, ivar, rng, n_draws=N_NOISE_DRAWS):
    """Observed indices with Monte-Carlo errors from the inverse variance.

    Inputs: rest-frame air wavelength [Angstrom]; flux and ivar (same units as the spectrum); rng; n_draws.
    Returns: dict with <name> and <name>_err for every name in INDEX_NAMES; NaN where index_validity is False.
    """
    validity = index_validity(wavelength, ivar)
    values = index_values(wavelength, flux)
    sigma = np.zeros_like(flux)
    good = ivar > 0
    sigma[good] = 1 / np.sqrt(ivar[good])

    draws = {name: [] for name in INDEX_NAMES}
    for _ in range(n_draws):
        noisy_flux = flux + rng.normal(0, 1, len(flux)) * sigma
        draw_values = index_values(wavelength, noisy_flux)
        for name in INDEX_NAMES:
            draws[name].append(draw_values[name])

    measured = {}
    for name in INDEX_NAMES:
        if validity[name]:
            measured[name] = values[name]
            measured[f"{name}_err"] = np.nanstd(draws[name])
        else:
            measured[name] = np.nan
            measured[f"{name}_err"] = np.nan
    return measured


# =============================================================================
# 2. Observed spectra for the indices
# =============================================================================
def rest_frame_spectrum(sdss_fits_path, redshift, ebv):
    """Native-resolution SDSS spectrum ready for index measurement.

    Inputs: path of an SDSS spec FITS file (HDU 1 with loglam, flux, ivar); redshift; Milky-Way E(B-V) [mag].
    Returns: (rest-frame air wavelength [Angstrom], dereddened flux, ivar), the flux in SDSS units.
    """
    with fits.open(sdss_fits_path) as spectrum_file:
        table = spectrum_file[1].data
        wavelength_vacuum = 10 ** table["loglam"].astype(float)
        flux = table["flux"].astype(float)
        ivar = table["ivar"].astype(float)
    wavelength_air = sm.vacuum_to_air(wavelength_vacuum)
    flux, ivar = sm.deredden_spectrum(wavelength_air, flux, ivar, ebv)
    return wavelength_air / (1 + redshift), flux, ivar


# =============================================================================
# 3. Posterior-predictive indices
# =============================================================================
N_PREDICTIVE_SAMPLES = 100     # posterior samples per galaxy for the model indices
MAX_Z_SCORE = 2.0              # an index agrees when |z| < MAX_Z_SCORE
# The rulers read the delayed-tau spectral fits: the validation stays on them (DPL spec, amendment 2)
RULER_SFH_MODEL = "delayed"


def evenly_spaced_rows(samples, n_samples):
    """At most n_samples rows of samples, evenly spaced (posterior samples are equal-weight)."""
    step = max(1, int(np.ceil(len(samples) / n_samples)))
    return samples[::step][:n_samples]


def model_flux_for_sample(model, cube, redshift, data_flux, data_ivar):
    """Model spectrum of one posterior sample on the model's grid, times the calibration polynomial fitted to the data."""
    components, _, _ = sm.spectral_components(cube, redshift, RULER_SFH_MODEL)
    model.update(components)
    flux = model.spectrum[:, 1]
    polynomial = sm.calibration_polynomial(flux, data_flux, data_ivar)
    return flux * polynomial


def predictive_indices(unit_cube_samples, redshift, observed_wavelength_air, data_flux, data_ivar,
                       n_samples=N_PREDICTIVE_SAMPLES):
    """Indices of the Bagpipes model spectra of posterior samples, measured exactly like the observed ones.

    Inputs: unit-cube posterior samples of a delayed-tau spectral fit (n, sm.spectral_dims(RULER_SFH_MODEL));
            redshift; the galaxy's native observed air wavelength grid [Angstrom] and its flux/ivar on that
            grid (restricted to sm.OBSERVED_RANGE by the caller).
    The calibration polynomial basis is uniform in pixel index, which differs slightly from the fit's 5 A bins
    (negligible for the indices).
    Returns: DataFrame (samples x INDEX_NAMES); a sample whose model flux is not positive in a band gets NaN there.
    """
    kept = evenly_spaced_rows(unit_cube_samples, n_samples)
    first_components, _, _ = sm.spectral_components(kept[0], redshift, RULER_SFH_MODEL)
    model = pipes.model_galaxy(first_components, spec_wavs=observed_wavelength_air)
    rest_wavelength = observed_wavelength_air / (1 + redshift)
    rows = []
    for cube in kept:
        flux = model_flux_for_sample(model, cube, redshift, data_flux, data_ivar)
        positive = np.where(flux > 0, 1.0, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):    # a zero-flux sample gives NaN, handled below
            values = index_values(rest_wavelength, flux)
        validity = index_validity(rest_wavelength, positive)    # zero or negative model flux invalidates a band
        for name in INDEX_NAMES:
            if not validity[name] or not np.isfinite(values[name]):
                values[name] = np.nan
        rows.append(values)
    return pd.DataFrame(rows, columns=INDEX_NAMES)


def index_comparison(measured, predictive):
    """Observed vs posterior-predictive indices of one galaxy.

    Inputs: measured (measure_indices output); predictive (predictive_indices output).
    Returns: dict with <name>_model_median, <name>_model_std, <name>_z = (obs - median)/sqrt(err^2 + std^2) and
             <name>_p (two-sided posterior-predictive p-value) for every name; NaN where either side is missing.
    The p-value is checked against replicated data = model + observational noise: for each posterior sample the
    replicate is N(model, err), so P(replicate < observed) = Phi((observed - model) / err), averaged over samples.
    Comparing with the noise-free model alone would flag every well-measured index as a misfit.
    """
    comparison = {}
    for name in INDEX_NAMES:
        model_values = predictive[name].to_numpy(dtype=float)
        model_values = model_values[np.isfinite(model_values)]
        observed = measured[name]
        if len(model_values) == 0 or not np.isfinite(observed):
            for suffix in ("model_median", "model_std", "z", "p"):
                comparison[f"{name}_{suffix}"] = np.nan
            continue
        median = np.median(model_values)
        spread = np.std(model_values)
        comparison[f"{name}_model_median"] = median
        comparison[f"{name}_model_std"] = spread
        comparison[f"{name}_z"] = (observed - median) / np.sqrt(measured[f"{name}_err"] ** 2 + spread**2)
        below = np.mean(norm.cdf((observed - model_values) / measured[f"{name}_err"]))
        comparison[f"{name}_p"] = 2 * min(below, 1 - below)
    return comparison


# =============================================================================
# 4. Light-weighted age and metallicity (for Gallazzi et al. 2005)
# =============================================================================
N_LIGHT_WEIGHT_SAMPLES = 50
LIGHT_WEIGHT_REDSHIFT = 0.001   # near rest frame: the r band then samples rest-frame r-band light


def light_weighted_log_mean(values, luminosities):
    """log10 of the luminosity-weighted linear mean, log10(sum(L v) / sum(L)): the log of the mean, not the mean of the log."""
    values = np.asarray(values, dtype=float)
    luminosities = np.asarray(luminosities, dtype=float)
    return np.log10(np.sum(luminosities * values) / np.sum(luminosities))


@functools.lru_cache(maxsize=1)
def ssp_r_band_luminosities():
    """r-band luminosity per unit formed mass of every Bagpipes SSP: grid metallicity x Bagpipes age bin.

    Built once from a dust-free, nebular-free model near rest frame (LIGHT_WEIGHT_REDSHIFT): its stellar grid
    holds each SSP spectrum per M_sun formed on Bagpipes' 0.1-dex age bins (config.age_sampling), and its r
    filter measures them exactly as it measures the model's own spectrum.
    Inputs: none.
    Returns: array of shape (len(pipes.config.metallicities), len(pipes.config.age_sampling)): filter-averaged
             F_lambda per M_sun formed, in units common to every entry (only ratios are used).
    """
    r_filter = [sm.FILTER_FILES[sm.BANDS.index("r")]]
    any_sfh = {"age": 1.0, "tau": 1.0, "massformed": sm.LOG_MASS_FORMED, "metallicity": 1.0}
    any_population = {"redshift": LIGHT_WEIGHT_REDSHIFT, "delayed": any_sfh}
    model = pipes.model_galaxy(any_population, filt_list=r_filter)
    n_metallicities, n_ages = model.stellar.grid.shape[1:]
    luminosities = np.zeros((n_metallicities, n_ages))
    for metallicity_index in range(n_metallicities):
        for age_index in range(n_ages):
            ssp_spectrum = model.stellar.grid[:, metallicity_index, age_index]
            photometry = model.filter_set.get_photometry(ssp_spectrum, LIGHT_WEIGHT_REDSHIFT)
            # Bagpipes returns a 0-d array for a one-filter list, so promote it before indexing
            luminosities[metallicity_index, age_index] = np.atleast_1d(photometry)[0]
    return luminosities


def r_band_light_per_age_step(star_formation, metallicity):
    """r-band light of each age step of a Bagpipes SFH: formed mass x SSP r luminosity per unit mass.

    It follows Bagpipes' own model: the SSP of a step is that of the 0.1-dex age bin it falls in (Bagpipes
    histograms the SFH into config.age_bins), interpolated linearly in Z between the grid metallicities and
    clamped at the ends (chemical_enrichment_history.delta). So the steps add up to the dust-free model's
    r flux.
    Inputs: star_formation: a Bagpipes star_formation_history (ages [yr], age_widths [yr], sfh [M_sun/yr]);
            metallicity: the single metallicity [Z/Z_sun].
    Returns: (ages [yr], r-band light per step [ssp_r_band_luminosities units x M_sun]), arrays of
             len(star_formation.ages).
    """
    formed_mass = star_formation.sfh * star_formation.age_widths

    # Linear interpolation is linear in the tabulated values: interpolating each unit vector gives the weights
    grid_metallicities = pipes.config.metallicities
    metallicity_weights = np.zeros(len(grid_metallicities))
    for grid_index, unit_vector in enumerate(np.eye(len(grid_metallicities))):
        metallicity_weights[grid_index] = np.interp(metallicity, grid_metallicities, unit_vector)
    luminosity_per_age_bin = metallicity_weights @ ssp_r_band_luminosities()

    # Steps past the last bin edge (the age of the Universe today) hold no mass: every age is below it
    age_bin = np.searchsorted(pipes.config.age_bins, star_formation.ages, side="right") - 1
    age_bin = np.clip(age_bin, 0, len(luminosity_per_age_bin) - 1)
    return star_formation.ages, formed_mass * luminosity_per_age_bin[age_bin]


def light_weighted_age_and_metallicity(unit_cube_samples, redshift, n_samples=N_LIGHT_WEIGHT_SAMPLES):
    """r-band light-weighted log age and log metallicity of posterior samples (Gallazzi et al. 2005 definition).

    Gallazzi et al. (2005, section 3.1) define the r-band light-weighted age as the linear mean
    t_r = sum(L_i t_i) / sum(L_i), quoted as log10(t_r). Here i runs over Bagpipes' SFH age steps
    (0.0025 dex) and L_i is the step's formed mass times the SSP r-band luminosity per unit mass at its age
    and the sample's metallicity (r_band_light_per_age_step; dust-free, without nebular emission). With a
    single metallicity the light-weighted Z is that metallicity. The S-PLUS r filter stands in for SDSS r
    (negligible difference).
    Inputs: unit_cube_samples: unit-cube posterior samples of a delayed-tau spectral fit (dimensionless,
            shape (n, len(sm.spectral_dims(RULER_SFH_MODEL))));
            redshift: the galaxy redshift (dimensionless; it sets the oldest allowed age);
            n_samples: at most this many evenly spaced samples are used (count).
    Returns: DataFrame with log_age_light [log yr] and log_metallicity_light [log Z/Z_sun], one row per kept
             sample.
    """
    kept = evenly_spaced_rows(unit_cube_samples, n_samples)
    star_formation = None
    rows = []
    for cube in kept:
        components, drawn, _ = sm.spectral_components(cube, redshift, RULER_SFH_MODEL)
        # One SFH object updated in place; only the SFH is needed, not a full model
        if star_formation is None:
            star_formation = pipes.models.star_formation_history(components)
        else:
            star_formation.update(components)
        ages, light = r_band_light_per_age_step(star_formation, drawn["metallicity"])
        rows.append({"log_age_light": light_weighted_log_mean(ages, light),
                     "log_metallicity_light": np.log10(drawn["metallicity"])})
    return pd.DataFrame(rows, columns=["log_age_light", "log_metallicity_light"])


# =============================================================================
# 5. Scaling relations
# =============================================================================
SPEAGLE_SCATTER = 0.2           # [dex] intrinsic scatter of the main sequence (Speagle et al. 2014)
MIN_POINTS_PER_BIN = 10         # relation bins with fewer points are dropped (percentiles too noisy)
PERCENTILES = (16, 50, 84)


def speagle_log_sfr(log_mstar, age_gyr):
    """Speagle et al. (2014) main sequence: log SFR [M_sun/yr] = (0.84 - 0.026 t) log M* - (6.51 - 0.11 t)."""
    return (0.84 - 0.026 * age_gyr) * np.asarray(log_mstar) - (6.51 - 0.11 * age_gyr)


def speagle_main_sequence(log_mstar, redshift):
    """Main-sequence log SFR at each galaxy's redshift (t = age of the Universe [Gyr] in the module cosmology)."""
    age_gyr = sm.COSMOLOGY.age(np.asarray(redshift)).value
    return speagle_log_sfr(log_mstar, age_gyr)


def relation_band(x, y, bin_edges):
    """Median and 16-84th percentiles of y in bins of x (a literature relation from a large sample).

    Returns: DataFrame with x_low, x_high, n, y_p16, y_median, y_p84 (bins with < MIN_POINTS_PER_BIN points are dropped).
    """
    rows = []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        in_bin = (x >= low) & (x < high) & np.isfinite(y)
        if in_bin.sum() < MIN_POINTS_PER_BIN:
            continue
        p16, p50, p84 = np.percentile(y[in_bin], PERCENTILES)
        rows.append({"x_low": low, "x_high": high, "n": int(in_bin.sum()),
                     "y_p16": p16, "y_median": p50, "y_p84": p84})
    return pd.DataFrame(rows)


def inside_band(x, y, band):
    """Whether each (x, y) lies inside the 16-84th band of its x bin; NA outside the binned x range or for NaN."""
    x = np.asarray(x, dtype=float)    # pandas nullable inputs become float with NaN
    y = np.asarray(y, dtype=float)
    result = []
    for x_value, y_value in zip(x, y):
        matching = band[(band["x_low"] <= x_value) & (band["x_high"] > x_value)]
        if len(matching) == 0 or not np.isfinite(y_value):
            result.append(pd.NA)
        else:
            row = matching.iloc[0]
            result.append(bool(row["y_p16"] <= y_value <= row["y_p84"]))
    return pd.array(result, dtype="boolean")


# =============================================================================
# 6. MPA-JHU (SDSS galSpec) query
# =============================================================================
MPA_JHU_COLUMNS = {
    "galSpecIndx": ["d4000_n", "d4000_n_err", "lick_hd_a", "lick_hd_a_err", "lick_mgb", "lick_mgb_err",
                    "lick_fe5270", "lick_fe5270_err", "lick_fe5335", "lick_fe5335_err"],
    "galSpecExtra": ["bptclass", "sfr_fib_p16", "sfr_fib_p50", "sfr_fib_p84", "lgm_tot_p50", "lgm_fib_p50"],
}


def mpa_jhu_query(specobjids):
    """SkyServer SQL selecting MPA-JHU indices and SFRs plus SDSS velocity dispersions for the given spectra.

    Inputs: SDSS specObjIDs as text (eBOSS ids exceed int64). Returns: the SQL string (one LEFT JOIN per table;
            spectra without galSpec rows come back with NULLs). Callers should chunk long lists (<= 200 ids).
    """
    for specobjid in specobjids:
        assert str(specobjid).isdigit(), f"specObjID must be a digit string, got {specobjid!r}"
    selected = ["CAST(s.specObjID AS VARCHAR(20)) AS specobjid_text", "s.velDisp", "s.velDispErr"]
    for column in MPA_JHU_COLUMNS["galSpecIndx"]:
        selected.append(f"i.{column}")
    for column in MPA_JHU_COLUMNS["galSpecExtra"]:
        selected.append(f"e.{column}")
    id_list = ", ".join(specobjids)
    return ("SELECT " + ", ".join(selected) + " FROM SpecObjAll AS s "
            "LEFT JOIN galSpecIndx AS i ON i.specObjID = s.specObjID "
            "LEFT JOIN galSpecExtra AS e ON e.specObjID = s.specObjID "
            f"WHERE s.specObjID IN ({id_list})")
