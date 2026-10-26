"""Tests for make_slim_firefly.py. Run from Synthetic_domain/."""

import numpy as np
import pandas as pd
from astropy.io import fits

import make_slim_firefly as slim


def write_firefly(path, rows):
    """A tiny Firefly-like FITS table with the slim columns plus one extra column."""
    columns = []
    for name in slim.FIREFLY_COLUMNS:
        values = np.array([row[name] for row in rows])
        fmt = "20A" if values.dtype.kind in "US" else "D"
        columns.append(fits.Column(name=name, format=fmt, array=values))
    columns.append(fits.Column(name="UNUSED", format="D", array=np.zeros(len(rows))))
    fits.BinTableHDU.from_columns(columns).writeto(path)


def firefly_row(specobjid, ra, dec, z, kind="GALAXY"):
    return {"SPECOBJID": specobjid, "RUN2D": "26", "PLATE": 1, "MJD": 2, "FIBERID": 3, "RA": ra, "DEC": dec,
            "Z_NOQSO": z, "CLASS_NOQSO": kind, "Chabrier_MILES_age_massW": 1.0}


def test_slim_keeps_only_matchable_galaxies_and_the_read_columns(tmp_path):
    splus = pd.DataFrame({"RA": [10.0, 20.0], "DEC": [0.0, 0.0], "PROB_GAL": [0.9, 0.1]})
    rows = [firefly_row("kept", 10.0, 0.0, 0.1),
            firefly_row("star", 10.0, 0.0, 0.1, kind="STAR"),
            firefly_row("out_of_prior", 10.0, 0.0, 5.0),
            firefly_row("far", 11.0, 0.0, 0.1),
            firefly_row("low_prob_gal", 20.0, 0.0, 0.1)]
    write_firefly(tmp_path / "full.fits", rows)

    n_kept = slim.write_slim_firefly(str(tmp_path / "full.fits"), splus, str(tmp_path / "slim.fits"))

    data = fits.getdata(tmp_path / "slim.fits", 1)
    assert n_kept == 1
    assert [value.strip() for value in data["SPECOBJID"]] == ["kept"]
    assert list(data.columns.names) == slim.FIREFLY_COLUMNS
