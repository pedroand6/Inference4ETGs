"""
Convert BPASS v2.2.1 single-star (sin) dat.gz files into a BAGPIPES-compatible
stellar grids FITS file.

BAGPIPES FITS structure (reverse-engineered from bc03_miles_stellar_grids.fits):
  HDU 0  : PrimaryHDU (empty)
  HDU 1  : ImageHDU  ZMET_X.XXXZSOL  shape (n_ages, n_wavelengths)  flux per unit mass [Lsun/AA/Msun]
  ...      one HDU per metallicity, in ascending order
  HDU -3 : ImageHDU  LIV_MSTAR_FRAC  shape (n_ages, n_metallicities+1)
  HDU -2 : ImageHDU  STELLAR_AGE_YR  shape (n_ages,)   ages in YEARS
  HDU -1 : ImageHDU  WAVELENGTHS_AA  shape (n_wavelengths,)  wavelengths in Angstroms

BPASS v2.2.1 file naming convention (single-star, Chabrier IMF, Mup=300):
  Spectra : spectra-sin-imf_chab300.z{ZZZ}.dat.gz
  Starmass: starmass-sin-imf_chab300.z{ZZZ}.dat.gz  (columns: log_age, total_mass, alive_mass, ...)

Metallicity codes vs Z/Zsun (Zsun=0.020):
  z001 = 0.001  → 0.05  Zsun
  z002 = 0.002  → 0.10  Zsun
  z003 = 0.003  → 0.15  Zsun
  z004 = 0.004  → 0.20  Zsun
  z006 = 0.006  → 0.30  Zsun
  z008 = 0.008  → 0.40  Zsun
  z010 = 0.010  → 0.50  Zsun
  z014 = 0.014  → 0.70  Zsun  (≈ solar, Asplund 2009)
  z020 = 0.020  → 1.00  Zsun  (Anders & Grevesse solar)
  z030 = 0.030  → 1.50  Zsun
  z040 = 0.040  → 2.00  Zsun

Usage:
    python bpass_to_bagpipes.py --bpass_dir /path/to/bpass_v2.2.1_sin \
                                 --output    /path/to/bagpipes/models/grids/bpass_sin_chab300_stellar_grids.fits

Then in your script, BEFORE importing bagpipes or creating any model:
    import bagpipes as pipes
    pipes.config.stellar_file = "bpass_sin_chab300_stellar_grids.fits"
    pipes.config.metallicities = np.array([0.05, 0.10, 0.20, 0.30, 0.40,
                                            0.50, 0.70, 1.00, 1.50, 2.00])
"""

import os
import gzip
import argparse
import numpy as np
from astropy.io import fits

# ---------------------------------------------------------------------------
# BPASS metallicity table
# ---------------------------------------------------------------------------
BPASS_ZMETS = {
    "z001": (0.001, 0.05),
    "z002": (0.002, 0.10),
    "z003": (0.003, 0.15),
    "z004": (0.004, 0.20),
    "z006": (0.006, 0.30),
    "z008": (0.008, 0.40),
    "z010": (0.010, 0.50),
    "z014": (0.014, 0.70),
    "z020": (0.020, 1.00),
    "z030": (0.030, 1.50),
    "z040": (0.040, 2.00),
}

# BPASS log10(age/yr) grid — fixed across all versions
BPASS_LOG_AGES = np.arange(6.0, 11.1, 0.1)          # 6.0 … 11.0  → 51 steps
BPASS_AGES_YR  = 10.0 ** BPASS_LOG_AGES              # years


def read_spectra(path):
    """
    Read a BPASS spectra-*.dat.gz file.

    Format: first column = wavelength (AA), subsequent columns = flux at each
    age step in units of Lsun / AA / 1e6 Msun (i.e. per 10^6 Msun formed).
    We convert to Lsun / AA / Msun by dividing by 1e6.

    Returns
    -------
    wavelengths : (n_wav,)
    fluxes      : (n_ages, n_wav)   Lsun/AA/Msun
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        data = np.loadtxt(fh)

    wavelengths = data[:, 0]                  # Angstroms
    # columns 1..51 are the 51 age steps; transpose → (n_ages, n_wav)
    fluxes = data[:, 1:].T / 1.0e6           # Lsun/AA/Msun
    return wavelengths, fluxes


def read_starmass(path, n_ages=51):
    """
    Read a BPASS starmass-*.dat.gz file.

    Columns (BPASS 2.2.1):
        0  log10(age/yr)
        1  total_mass        (mass still in stars + remnants, per 1e6 Msun)
        2  alive_mass        (mass in living stars only, per 1e6 Msun)
        3  remnant_mass
        4  ...

    Returns
    -------
    live_frac : (n_ages,)  fraction of initial mass in living stars
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        data = np.loadtxt(fh)

    # alive_mass / 1e6 gives fraction of 1 Msun formed that is alive
    alive = data[:n_ages, 2] / 1.0e6
    alive = np.clip(alive, 0.0, 1.0)
    return alive


def build_fits(bpass_dir, output_path, imf_tag="imf_chab300", kind="sin",
               zmet_keys=None):
    """
    Build the BAGPIPES-format FITS from BPASS dat.gz files.

    Parameters
    ----------
    bpass_dir : str   directory containing BPASS dat.gz files
    output_path : str output FITS path
    imf_tag : str     IMF string as it appears in BPASS filenames
    kind : str        'sin' or 'bin'
    zmet_keys : list  subset of BPASS_ZMETS keys to include (None = all)
    """
    if zmet_keys is None:
        zmet_keys = sorted(BPASS_ZMETS.keys())

    n_ages = len(BPASS_AGES_YR)

    # ------------------------------------------------------------------ #
    # 1. Read all metallicities                                            #
    # ------------------------------------------------------------------ #
    grids      = []   # list of (n_ages, n_wav) arrays, one per metallicity
    live_fracs = []   # list of (n_ages,) arrays
    wavelengths = None
    valid_keys  = []
    zsol_values = []

    for zkey in zmet_keys:
        z_abs, z_sol = BPASS_ZMETS[zkey]

        spec_name  = f"spectra-{kind}-{imf_tag}.{zkey}.dat.gz"
        mass_name  = f"starmass-{kind}-{imf_tag}.{zkey}.dat.gz"

        spec_path  = os.path.join(bpass_dir, spec_name)
        mass_path  = os.path.join(bpass_dir, mass_name)

        # Try without .gz as fallback
        if not os.path.exists(spec_path):
            spec_path = spec_path[:-3]
        if not os.path.exists(mass_path):
            mass_path = mass_path[:-3]

        if not os.path.exists(spec_path):
            print(f"  WARNING: spectra file not found for {zkey}, skipping.")
            continue
        if not os.path.exists(mass_path):
            print(f"  WARNING: starmass file not found for {zkey}, skipping.")
            continue

        print(f"  Reading {zkey}  (Z/Zsun = {z_sol:.2f}) ...", end=" ", flush=True)
        wavs, flux = read_spectra(spec_path)
        live       = read_starmass(mass_path, n_ages=n_ages)

        if wavelengths is None:
            wavelengths = wavs
        else:
            assert np.allclose(wavelengths, wavs), \
                f"Wavelength mismatch at {zkey}"

        # Verify age dimension
        assert flux.shape[0] == n_ages, \
            f"Expected {n_ages} ages, got {flux.shape[0]} in {spec_name}"

        grids.append(flux)
        live_fracs.append(live)
        valid_keys.append(zkey)
        zsol_values.append(z_sol)
        print("OK")

    n_met = len(valid_keys)
    n_wav = len(wavelengths)
    print(f"\nIncluded metallicities ({n_met}): {zsol_values}")
    print(f"Wavelength grid: {n_wav} points, {wavelengths[0]:.1f}–{wavelengths[-1]:.1f} AA")
    print(f"Age grid: {n_ages} points, {BPASS_AGES_YR[0]:.2e}–{BPASS_AGES_YR[-1]:.2e} yr")

    # ------------------------------------------------------------------ #
    # 2. Build LIV_MSTAR_FRAC array                                        #
    #    BAGPIPES shape: (n_ages, n_met+1)                                 #
    #    Column 0 is unused padding (set to 1); columns 1..n_met are the   #
    #    live fractions per metallicity (matching bc03 convention).         #
    # ------------------------------------------------------------------ #
    live_arr = np.ones((n_ages, n_met + 1), dtype=np.float64)
    for j, lf in enumerate(live_fracs):
        live_arr[:, j + 1] = lf

    # ------------------------------------------------------------------ #
    # 3. Assemble HDUList                                                  #
    # ------------------------------------------------------------------ #
    hdulist = fits.HDUList()
    hdulist.append(fits.PrimaryHDU())   # HDU 0 — empty primary

    for j, (zkey, zsol, grid) in enumerate(zip(valid_keys, zsol_values, grids)):
        name = f"ZMET_{zsol:.3f}ZSOL"
        hdu  = fits.ImageHDU(data=grid.astype(np.float64), name=name)
        hdu.header["BUNIT"]   = "Lsun/AA/Msun"
        hdu.header["ZSOL"]    = zsol
        hdu.header["ZABS"]    = BPASS_ZMETS[zkey][0]
        hdu.header["COMMENT"] = f"BPASS {kind} {imf_tag} - {zkey}"
        hdulist.append(hdu)

    hdulist.append(fits.ImageHDU(data=live_arr.astype(np.float64),
                                  name="LIV_MSTAR_FRAC"))
    hdulist.append(fits.ImageHDU(data=BPASS_AGES_YR.astype(np.float64),
                                  name="STELLAR_AGE_YR"))
    hdulist.append(fits.ImageHDU(data=wavelengths.astype(np.float64),
                                  name="WAVELENGTHS_AA"))

    # ------------------------------------------------------------------ #
    # 4. Write                                                             #
    # ------------------------------------------------------------------ #
    hdulist.writeto(output_path, overwrite=True)
    print(f"\nWritten: {output_path}")
    print(f"Total HDUs: {len(hdulist)}  "
          f"(1 primary + {n_met} metallicity grids + 3 auxiliary)")

    # ------------------------------------------------------------------ #
    # 5. Print config snippet                                              #
    # ------------------------------------------------------------------ #
    fname = os.path.basename(output_path)
    print("\n" + "="*60)
    print("Add this to your script BEFORE importing bagpipes models:")
    print("="*60)
    print("import bagpipes as pipes")
    print("import numpy as np")
    print(f'pipes.config.stellar_file = "{fname}"')
    met_str = ", ".join(f"{z:.3f}" for z in zsol_values)
    print(f"pipes.config.metallicities = np.array([{met_str}])")
    print("# Also update metallicity_bins:")
    print("from bagpipes.utils import make_bins")
    print("pipes.config.metallicity_bins = make_bins(pipes.config.metallicities, make_rhs=True)[0]")
    print("pipes.config.metallicity_bins[0]  = 0.")
    print("pipes.config.metallicity_bins[-1] = 10.")
    print("="*60)

    return zsol_values


# ---------------------------------------------------------------------------
# Validation helper
# ---------------------------------------------------------------------------
def validate_fits(output_path, metallicities):
    """Quick sanity check on the written FITS file."""
    print("\nValidating output FITS...")
    f = fits.open(output_path)
    n_met = len(metallicities)

    assert len(f) == n_met + 4, f"Expected {n_met+4} HDUs, got {len(f)}"

    ages = f[-2].data
    wavs = f[-1].data
    live = f[-3].data

    assert ages.ndim == 1, "Ages must be 1D"
    assert wavs.ndim == 1, "Wavelengths must be 1D"
    assert live.shape == (len(ages), n_met + 1), \
        f"live_frac shape mismatch: {live.shape}"

    for i in range(1, n_met + 1):
        grid = f[i].data
        assert grid.shape == (len(ages), len(wavs)), \
            f"HDU {i} shape {grid.shape} != ({len(ages)}, {len(wavs)})"
        assert np.all(np.isfinite(grid)), f"NaN/Inf in HDU {i}"
        assert np.all(grid >= 0), f"Negative flux in HDU {i}"

    print(f"  HDUs         : {len(f)}")
    print(f"  Metallicities: {n_met}")
    print(f"  Ages         : {len(ages)}  ({ages[0]:.2e} – {ages[-1]:.2e} yr)")
    print(f"  Wavelengths  : {len(wavs)}  ({wavs[0]:.1f} – {wavs[-1]:.1f} AA)")
    print(f"  Grid shape   : {f[1].data.shape}")
    print(f"  live_frac    : min={live[:,1:].min():.4f}  max={live[:,1:].max():.4f}")
    print("  All checks PASSED")
    f.close()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert BPASS v2.2.1 dat.gz files to BAGPIPES FITS grid")
    parser.add_argument("--bpass_dir", required=True,
                        help="Directory containing BPASS dat.gz files")
    parser.add_argument("--output", required=True,
                        help="Output FITS path (e.g. .../grids/bpass_sin_chab300_stellar_grids.fits)")
    parser.add_argument("--kind", default="sin", choices=["sin", "bin"],
                        help="Single (sin) or binary (bin) stellar evolution")
    parser.add_argument("--imf", default="imf_chab300",
                        help="IMF tag as it appears in BPASS filenames")
    parser.add_argument("--zmets", nargs="+", default=None,
                        choices=list(BPASS_ZMETS.keys()),
                        help="Metallicity codes to include (default: all available)")
    parser.add_argument("--validate", action="store_true",
                        help="Run validation after writing")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    zsol_vals = build_fits(
        bpass_dir   = args.bpass_dir,
        output_path = args.output,
        imf_tag     = args.imf,
        kind        = args.kind,
        zmet_keys   = args.zmets,
    )

    if args.validate:
        validate_fits(args.output, zsol_vals)