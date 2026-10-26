"""
Writes a slim copy of the Firefly DR16 catalogue for spectroscopic_reference.ipynb: only the columns its cell 3
reads and only the spectra that can match an S-PLUS galaxy. The full file is 2.4 GB, too large for the job
server, which buffers whole uploads in a small memory limit; the slim one is a few MB.

The slim rows are a superset of what cell 3 keeps (galaxies inside sm.Z_RANGE within SLIM_MATCH_RADIUS of an
S-PLUS galaxy with PROB_GAL above MIN_PROB_GALAXY), and every spectrum is matched on its own, so the notebook
gives the same result with either file. Remake it whenever dr3_corrected.csv changes.

Usage (from Synthetic_domain/):  python make_slim_firefly.py
Then submit notebook 2 with ../Data/sdss_eboss_firefly-dr16_slim.fits instead of the full file
(remote/README.md); cell 1 reads the slim file whenever it exists.
"""

import argparse

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.io import fits

import synthetic_model as sm

FULL_FILE = "../Data/sdss_eboss_firefly-dr16.fits"
SLIM_FILE = "../Data/sdss_eboss_firefly-dr16_slim.fits"
SPLUS_FILE = "dr3_corrected.csv"
# The columns cell 3 of spectroscopic_reference.ipynb reads (its FIREFLY_COLUMNS)
FIREFLY_COLUMNS = ["SPECOBJID", "RUN2D", "PLATE", "MJD", "FIBERID", "RA", "DEC", "Z_NOQSO", "CLASS_NOQSO",
                   "Chabrier_MILES_age_massW"]
MIN_PROB_GALAXY = 0.8           # as the notebook's MIN_PROB_GALAXY
SLIM_MATCH_RADIUS = 5.0         # [arcsec] wider than the notebook's MATCH_RADIUS (1.5), so the rows are a superset


def write_slim_firefly(full_file, splus, slim_file, min_prob_galaxy=MIN_PROB_GALAXY,
                       match_radius=SLIM_MATCH_RADIUS):
    """Write the Firefly rows and columns spectroscopic_reference.ipynb can use to slim_file.

    Inputs: full_file, slim_file: FITS paths; splus: DataFrame with RA, DEC [deg] and PROB_GAL;
            min_prob_galaxy: PROB_GAL limit; match_radius: [arcsec].
    Returns: the number of rows written.
    """
    with fits.open(full_file, memmap=True) as hdus:
        data = hdus[1].data
        kind = np.char.strip(np.asarray(data["CLASS_NOQSO"]).astype(str))
        redshift = np.asarray(data["Z_NOQSO"], dtype=float)
        keep = (kind == "GALAXY") & (redshift > sm.Z_RANGE[0]) & (redshift < sm.Z_RANGE[1])

        galaxies = splus[splus["PROB_GAL"] > min_prob_galaxy]
        candidates = np.flatnonzero(keep)
        spectrum_coordinates = SkyCoord(np.asarray(data["RA"])[candidates] * u.deg,
                                        np.asarray(data["DEC"])[candidates] * u.deg)
        splus_coordinates = SkyCoord(galaxies["RA"].to_numpy() * u.deg, galaxies["DEC"].to_numpy() * u.deg)
        _, separation, _ = match_coordinates_sky(spectrum_coordinates, splus_coordinates)
        rows = candidates[separation.arcsec < match_radius]

        columns = []
        for column in hdus[1].columns:
            if column.name in FIREFLY_COLUMNS:
                columns.append(fits.Column(name=column.name, format=column.format,
                                           array=np.asarray(data[column.name])[rows]))
    columns.sort(key=lambda column: FIREFLY_COLUMNS.index(column.name))
    fits.BinTableHDU.from_columns(columns).writeto(slim_file, overwrite=True)
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description="Write the slim Firefly catalogue of notebook 2.")
    parser.add_argument("--full", default=FULL_FILE)
    parser.add_argument("--splus", default=SPLUS_FILE)
    parser.add_argument("--output", default=SLIM_FILE)
    arguments = parser.parse_args()
    n_rows = write_slim_firefly(arguments.full, pd.read_csv(arguments.splus), arguments.output)
    print(f"Wrote {n_rows} Firefly rows to {arguments.output}")


if __name__ == "__main__":
    main()
