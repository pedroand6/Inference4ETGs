"""Fetch the literature tables used to validate the Bagpipes fits (run once from Synthetic_domain/; re-runnable).

Produces, under literature/:
  mpa_jhu_reference.csv           MPA-JHU / SDSS rows for our reference galaxies (SDSS-I/II spectra only)
  gallazzi_reference.csv          Gallazzi et al. (2005) age, metallicity and mass for our reference galaxies
  relations/mass_metallicity.csv  stellar mass-metallicity band from the whole Gallazzi catalogue
  relations/age_sigma_passive.csv age versus velocity dispersion band for passive Gallazzi galaxies
Downloads in literature/raw/ are never overwritten. Sources and conventions: literature/SOURCES.md.
"""
import io
import os
import time
import urllib.request

import numpy as np
import pandas as pd
from astroquery.sdss import SDSS

import spectral_validation as sv

# =============================================================================
# 0. Constants
# =============================================================================
LITERATURE_DIR = "literature"
RAW_DIR = os.path.join(LITERATURE_DIR, "raw")
RELATIONS_DIR = os.path.join(LITERATURE_DIR, "relations")
REFERENCE_CSV = "spectroscopic_reference.csv"

GALLAZZI_BASE_URL = "https://wwwmpa.mpa-garching.mpg.de/SDSS/DR4/Data/Gallazzi/"
GALLAZZI_FILES = {"age": "all_stat_age.dat.gz", "log_z": "all_stat_z_log.dat.gz", "log_mstar": "all_stat_mstar.dat.gz"}
# Columns of each Gallazzi file: plate, mjd, fiber, then the percentiles below, then the mode and the row index
GALLAZZI_COLUMN_NAMES = ["plate", "mjd", "fiberid", "p2_5", "p16", "p50", "p84", "p97_5", "mode", "index"]
GALLAZZI_NO_ESTIMATE = -99.0    # sentinel the catalogue uses where a galaxy has no estimate
LOG_SOLAR_Z_GALLAZZI = np.log10(0.02)   # catalogue log Z is not in solar units; subtract this (Z_sun = 0.02)

SDSS_DATA_RELEASE = 16
IDS_PER_QUERY = 100             # polite size; 200 ids make the GET URL too long for SkyServer (HTML error page)
RETRY_WAIT_SECONDS = 30         # wait before the single retry of a failed chunk
LEGACY_RUN2D = 26               # SDSS-I/II spectra reprocessed in DR16 (the ones MPA-JHU covers)

MISSING_VALUE_SENTINELS = [-9999.0, -99.0]   # MPA-JHU placeholders for "not measured"
# Index columns whose value is a placeholder (0 with a non-positive error) when the index was not measured
INDEX_COLUMNS = ["d4000_n", "lick_hd_a", "lick_mgb", "lick_fe5270", "lick_fe5335"]
DROPPED_IDS_CSV = os.path.join(LITERATURE_DIR, "dropped_ids.csv")

PASSIVE_D4000N_MIN = 1.6        # passive cut on the 4000 A break
MAX_PASSIVE_SAMPLE = 5000       # random Gallazzi galaxies queried for the age-sigma relation
RANDOM_SEED = 20261025
MASS_BIN_EDGES = np.round(np.arange(8.5, 12.01, 0.25), 3)    # log M* [Msun]
SIGMA_BIN_EDGES = np.round(np.arange(1.8, 2.61, 0.1), 3)     # log10 velocity dispersion [km/s]


def run_sql(sql):
    """One SkyServer SQL query as a DataFrame with the specobjid_text column kept as text.

    astropy's CSV reader (behind SDSS.query_sql) overflows on 19-digit ids, so the raw response is parsed here.
    Raises on an HTTP error or a response that is not a CSV table.
    """
    response = SDSS.query_sql_async(sql, data_release=SDSS_DATA_RELEASE)
    response.raise_for_status()
    lines = response.text.splitlines()
    assert lines and lines[0].startswith("#Table1"), f"unexpected SkyServer response: {response.text[:200]!r}"
    return pd.read_csv(io.StringIO("\n".join(lines[1:])), dtype={"specobjid_text": str})


def query_with_splitting(chunk, dropped):
    """run_sql for one chunk; a chunk that still fails is split in halves until the failing ids stand alone.

    SkyServer returns an error page (not a table) when a row holds a non-finite float it cannot write as CSV, so one
    such spectrum poisons its whole chunk. Failing single ids are appended to dropped. Returns: list of DataFrames.
    """
    try:
        return [run_sql(sv.mpa_jhu_query(chunk))]
    except Exception as error:
        if len(chunk) == 1:
            dropped.append(chunk[0])
            return []
        print(f"query of {len(chunk)} ids failed ({str(error)[:40]}); splitting")
    middle = len(chunk) // 2
    return query_with_splitting(chunk[:middle], dropped) + query_with_splitting(chunk[middle:], dropped)


def clean_mpa_values(rows):
    """Turn MPA-JHU missing-value placeholders into NaN.

    -9999 and -99 become NaN in every numeric column except bptclass, whose -1 is a genuine "unclassified" code.
    For each index (D4000n, Lick) a non-positive error marks a placeholder, so the value and its error become NaN.
    Returns: the cleaned DataFrame.
    """
    rows = rows.copy()
    numeric_columns = [c for c in rows.columns if c not in ["specobjid", "bptclass"]]
    rows[numeric_columns] = rows[numeric_columns].replace(MISSING_VALUE_SENTINELS, np.nan)
    for index in INDEX_COLUMNS:
        not_measured = rows[f"{index}_err"] <= 0
        rows.loc[not_measured, [index, f"{index}_err"]] = np.nan
    return rows


def query_in_chunks(specobjids):
    """MPA-JHU / SDSS rows for the given specObjIDs (digit strings), IDS_PER_QUERY per SkyServer query.

    A failed chunk is retried once after RETRY_WAIT_SECONDS, then split (see query_with_splitting); ids that fail
    alone are reported and left out. Returns: DataFrame with specobjid (stripped text) and the columns of
    sv.mpa_jhu_query with placeholders as NaN; spectra without a row are absent.
            Also the list of ids that could not be queried.
    """
    tables = []
    dropped = []
    for start in range(0, len(specobjids), IDS_PER_QUERY):
        chunk = list(specobjids[start:start + IDS_PER_QUERY])
        try:
            tables.append(run_sql(sv.mpa_jhu_query(chunk)))
        except Exception as error:
            print(f"chunk at {start} failed ({str(error)[:40]}); retrying once")
            time.sleep(RETRY_WAIT_SECONDS)
            try:
                tables.append(run_sql(sv.mpa_jhu_query(chunk)))
            except Exception:
                tables.extend(query_with_splitting(chunk, dropped))
    if dropped:
        print(f"{len(dropped)} ids could not be queried (SkyServer error on their rows) and are left out")
    rows = pd.concat(tables, ignore_index=True)
    rows = rows.rename(columns={"specobjid_text": "specobjid"})
    rows["specobjid"] = rows["specobjid"].astype(str).str.strip()
    return clean_mpa_values(rows), dropped


def legacy_specobjid(plate, mjd, fiberid):
    """SDSS specObjID of a legacy (RUN2D 26) spectrum, as text. Layout: plate<<50 | fiber<<38 | (mjd-50000)<<24 | run2d<<10."""
    value = (int(plate) << 50) + (int(fiberid) << 38) + ((int(mjd) - 50000) << 24) + (LEGACY_RUN2D << 10)
    return str(value)


def read_gallazzi_catalogue():
    """All three Gallazzi files merged on row position into one DataFrame (plate, mjd, fiberid, age/log_z/log_mstar p16/p50/p84).

    Ages are log yr, log_z is converted to log Z/Z_sun, log_mstar is log Msun; entries without an estimate become NaN.
    """
    catalogue = None
    for name, file_name in GALLAZZI_FILES.items():
        table = pd.read_csv(os.path.join(RAW_DIR, file_name), sep=r"\s+", header=None, names=GALLAZZI_COLUMN_NAMES)
        if catalogue is None:
            catalogue = table[["plate", "mjd", "fiberid"]].copy()
        same_order = (table[["plate", "mjd", "fiberid"]].values == catalogue[["plate", "mjd", "fiberid"]].values).all()
        assert same_order, f"{file_name} is not row-aligned with the other Gallazzi files"
        for percentile in ["p16", "p50", "p84"]:
            values = table[percentile].where(table[percentile] != GALLAZZI_NO_ESTIMATE)
            if name == "log_z":
                values = values - LOG_SOLAR_Z_GALLAZZI
            catalogue[f"gallazzi_{name}_{percentile}"] = values
    return catalogue


# =============================================================================
# 1. Reference galaxies
# =============================================================================
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(RELATIONS_DIR, exist_ok=True)
reference = pd.read_csv(REFERENCE_CSV, dtype={"specobjid": str})
reference["specobjid"] = reference["specobjid"].str.strip()
is_legacy = reference["RUN2D"] == str(LEGACY_RUN2D)

# =============================================================================
# 2. MPA-JHU rows (BOSS spectra come back with NULL galSpec columns, which documents their absence)
# =============================================================================
mpa_rows, dropped_reference = query_in_chunks(list(reference["specobjid"]))
mpa_rows.to_csv(os.path.join(LITERATURE_DIR, "mpa_jhu_reference.csv"), index=False)
# SpecObjAll rows exist for every spectrum, so count galaxies with real MPA-JHU measurements (galSpecExtra join)
mass_ids = set(mpa_rows.loc[mpa_rows["lgm_tot_p50"].notna(), "specobjid"])
sfr_ids = set(mpa_rows.loc[mpa_rows["sfr_fib_p50"].notna(), "specobjid"])

# =============================================================================
# 3. Gallazzi et al. (2005) catalogue: download (once), match to our galaxies, build the relation bands
# =============================================================================
for file_name in GALLAZZI_FILES.values():
    destination = os.path.join(RAW_DIR, file_name)
    if not os.path.exists(destination):
        urllib.request.urlretrieve(GALLAZZI_BASE_URL + file_name, destination)
gallazzi = read_gallazzi_catalogue()

ours = reference[["specobjid", "plate", "mjd", "fiberid"]].merge(gallazzi, on=["plate", "mjd", "fiberid"], how="inner")
gallazzi_columns = ["specobjid"] + [c for c in gallazzi.columns if c.startswith("gallazzi_")]
gallazzi_out = ours[gallazzi_columns].rename(columns={
    "gallazzi_age_p16": "gallazzi_log_age_p16", "gallazzi_age_p50": "gallazzi_log_age_p50",
    "gallazzi_age_p84": "gallazzi_log_age_p84", "gallazzi_log_mstar_p50": "gallazzi_log_mstar"})
gallazzi_out = gallazzi_out.drop(columns=["gallazzi_log_mstar_p16", "gallazzi_log_mstar_p84"])
gallazzi_out.to_csv(os.path.join(LITERATURE_DIR, "gallazzi_reference.csv"), index=False)

# 3a. Mass-metallicity band from the whole catalogue
mass_metallicity = sv.relation_band(gallazzi["gallazzi_log_mstar_p50"].values, gallazzi["gallazzi_log_z_p50"].values,
                                    MASS_BIN_EDGES)
mass_metallicity.to_csv(os.path.join(RELATIONS_DIR, "mass_metallicity.csv"), index=False)

# 3b. Age-sigma band for passive galaxies. A random sample is queried (not the whole catalogue) to keep the
#     SkyServer load modest, then cut on D4000n; the age column is the Gallazzi r-band light-weighted age.
has_age = gallazzi.dropna(subset=["gallazzi_age_p50"])
sample = has_age.sample(n=min(MAX_PASSIVE_SAMPLE, len(has_age)), random_state=RANDOM_SEED)
sample_ids = [legacy_specobjid(p, m, f) for p, m, f in zip(sample["plate"], sample["mjd"], sample["fiberid"])]
sample = sample.assign(specobjid=sample_ids)
sdss_rows, dropped_sample = query_in_chunks(sample_ids)
sdss_rows = sdss_rows[["specobjid", "velDisp", "d4000_n"]]
dropped = [(i, "reference") for i in dropped_reference] + [(i, "age_sigma_sample") for i in dropped_sample]
pd.DataFrame(dropped, columns=["id", "query"]).to_csv(DROPPED_IDS_CSV, index=False)
passive = sample.merge(sdss_rows, on="specobjid", how="inner")
passive = passive[(passive["d4000_n"] > PASSIVE_D4000N_MIN) & (passive["velDisp"] > 0)]
age_sigma = sv.relation_band(np.log10(passive["velDisp"].values), passive["gallazzi_age_p50"].values, SIGMA_BIN_EDGES)
age_sigma.to_csv(os.path.join(RELATIONS_DIR, "age_sigma_passive.csv"), index=False)

# =============================================================================
# 4. Counts
# =============================================================================
print(f"reference galaxies: {len(reference)} ({is_legacy.sum()} legacy, {(~is_legacy).sum()} BOSS)")
for label, ids in [("finite lgm_tot_p50", mass_ids), ("finite sfr_fib_p50", sfr_ids)]:
    print(f"with {label}: {reference['specobjid'].isin(ids).sum()} "
          f"(legacy {reference.loc[is_legacy, 'specobjid'].isin(ids).sum()}, "
          f"BOSS {reference.loc[~is_legacy, 'specobjid'].isin(ids).sum()})")
print(f"ids not queryable: {len(dropped_reference)} reference, {len(dropped_sample)} age-sigma sample")
print(f"with Gallazzi rows: {len(gallazzi_out)} ({gallazzi_out['gallazzi_log_age_p50'].notna().sum()} with an age)")
print(f"mass-metallicity bins: {len(mass_metallicity)}; passive galaxies in age-sigma: {len(passive)}, "
      f"bins: {len(age_sigma)}")
