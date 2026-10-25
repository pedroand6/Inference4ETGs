"""
Rulers for the Bagpipes spectral fits that do not rely on full-spectrum SED fitting:
model-independent spectral indices (compared with the fits' posterior predictions), light-weighted ages and
metallicities (for Gallazzi et al. 2005), the Speagle et al. (2014) main sequence and percentile bands of
literature scaling relations.
"""

import numpy as np
import pandas as pd
from astropy.io import fits

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
