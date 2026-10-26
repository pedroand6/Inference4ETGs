# Literature data used to validate the Bagpipes fits

Retrieved on 2026-10-25 by `fetch_literature.py` (run from `Synthetic_domain/`). Raw downloads live in `literature/raw/`
(git-ignored, never overwritten by the script).

Refreshed on 2026-10-26 for all 168 reference galaxies (153 SDSS legacy, 15 BOSS), up from the 60 of the first run:
- `mpa_jhu_reference.csv`: 168 rows; 142 with a finite `lgm_tot_p50` and 118 with a finite `sfr_fib_p50` (all legacy;
  the BOSS rows have empty galSpec columns);
- `gallazzi_reference.csv`: 144 rows (legacy spectra on DR4 plates), 124 of them with an age;
- the 60 and 51 rows of the first run came back with identical values;
- `relations/*.csv` and `dropped_ids.csv` came back identical to the 2026-10-25 files.

## Gallazzi et al. (2005) stellar ages, metallicities and masses

- Paper: Gallazzi, Charlot, Brinchmann, White, Tremonti 2005, MNRAS 362, 41 (arXiv:astro-ph/0506539).
- Public tables (DR4 version, 567 486 spectra), MPA-JHU page
  https://wwwmpa.mpa-garching.mpg.de/SDSS/DR4/Data/stellarmet.html (VizieR was not needed):
  - `https://wwwmpa.mpa-garching.mpg.de/SDSS/DR4/Data/Gallazzi/all_stat_age.dat.gz` r-band light-weighted log age [log yr]
  - `.../Gallazzi/all_stat_z_log.dat.gz` stellar metallicity, log Z (NOT in solar units; the script subtracts
    log10(0.02) so that the values are log Z/Z_sun, as the page instructs)
  - `.../Gallazzi/all_stat_mstar.dat.gz` log stellar mass [log Msun] from z-band model magnitudes (total mass).
    The Petrosian-magnitude masses (`all_stat_mstar.petromag.dat.gz`, about 0.1 dex lower) are not used.
- File columns (whitespace separated, no header; `-99` = no estimate, converted to NaN): plate, mjd, fiberid,
  2.5th, 16th, 50th, 84th, 97.5th percentile of the likelihood distribution, mode, row index into the MPA-JHU raw files.
  We match to our galaxies by plate, mjd and fiberid.
- Light-weighting definition (paper, section 3.1, items (i)-(ii), and its discussion of ages): the age is
  `t_r = int dtau psi(t - tau) f_r(tau) tau / int dtau psi(t - tau) f_r(tau)`, the r-band-flux-weighted LINEAR mean of
  stellar age (f_r = r-band flux of stars of age tau), quoted as log10. That is the log of the mean, not the mean of the
  log. Each model has a single metallicity for all its stars, interpreted as the "(optical) light-weighted" metallicity.
- `sv.light_weighted_age_and_metallicity` follows this definition: the age is `log10(sum(L_i t_i) / sum(L_i))`
  (`sv.light_weighted_log_mean`), and the metallicity is the model's single Z.

## MPA-JHU (SDSS galSpec) values

- Queried from SkyServer DR16 (https://skyserver.sdss.org/dr16, tables `SpecObjAll`, `galSpecIndx`, `galSpecExtra`)
  with `sv.mpa_jhu_query`, 100 specObjIDs per query. The queries run against DR16 (`SDSS_DATA_RELEASE = 16`); the
  documentation URL https://www.sdss4.org/dr17/spectro/galaxy_mpajhu/ is the DR17 page (the galSpec tables are the same).
- `mpa_jhu_reference.csv`: `specobjid` (text), `velDisp`, `velDispErr` [km/s] (SDSS pipeline), galSpecIndx `d4000_n`,
  `lick_hd_a`, `lick_mgb`, `lick_fe5270`, `lick_fe5335` (+ `_err`; Angstrom except D4000n), galSpecExtra `bptclass`
  (1 = star forming), `sfr_fib_p16/p50/p84` [log Msun/yr, fibre], `lgm_tot_p50`, `lgm_fib_p50` [log Msun].
  SDSS-I/II spectra (RUN2D 26) have these rows; BOSS spectra (RUN2D v5_13_0) come back with NULL (empty) galSpec columns
  but still carry `velDisp`. A legacy spectrum with no row (a low-S/N one) also has empty galSpec columns.

Missing values: the script converts the placeholders -9999 and -99 to NaN in every numeric MPA-JHU column, and for D4000n
and each Lick index sets value and error to NaN where the error is <= 0 (placeholder 0.0 with error -1). `bptclass` keeps
-1, which is a genuine "unclassified" code. Spectra the SkyServer could not return (error pages on rows with non-finite
floats) are listed in `dropped_ids.csv` (columns id, query), so one can check whether a reference galaxy was ever lost.

## Passive age-sigma relation

Reproducibility: `age_sigma_passive.csv` is reproducible only to a few galaxies, because SkyServer chunk failures vary from run to
run (a re-run differed by about 2 galaxies, about 0.001 dex).


`relations/age_sigma_passive.csv`: 5 000 random Gallazzi galaxies (seed in the script) were looked up in SkyServer by their
specObjID (rebuilt from plate, mjd, fiber with RUN2D 26); kept those with `d4000_n` > 1.6 and `velDisp` > 0.
x = log10 velDisp [km/s], y = Gallazzi log age p50 [log yr]; bins of 0.1 dex (bins with < 10 galaxies dropped).
`relations/mass_metallicity.csv`: whole Gallazzi catalogue, x = log M* (model-magnitude mass) [log Msun], y = log Z/Z_sun
p50, bins of 0.25 dex. Both files hold x_low, x_high, n, y_p16, y_median, y_p84 (`sv.relation_band`).

## Initial mass functions

| Quantity | IMF | Where it comes from |
|---|---|---|
| MPA-JHU stellar masses (`lgm_*`) | Kroupa (2001) | stated in the SDSS DR17 MPA-JHU page |
| MPA-JHU SFRs (`sfr_fib_*`) | Kroupa (2001), by inheritance | Brinchmann et al. (2004) method; the SDSS page does not name an IMF for SFRs (not verified there) |
| Gallazzi (2005) ages and metallicities | IMF-independent in practice | single-metallicity BC03 models; no IMF dependence stated |
| Gallazzi (2005) stellar masses | Kroupa (2001), INFERRED | The paper text and the DR4 pages do not state the IMF. The model library follows Kauffmann et al. (2003), which uses a Kroupa (2001) IMF (the same convention as the MPA-JHU masses above). Not confirmed in the Gallazzi paper itself. |
| Bagpipes `bc03_miles`, Speagle et al. (2014) | Kroupa | project convention |

Task 4 applies the Kroupa to Chabrier offset only where our masses meet a Chabrier-based axis. The Gallazzi mass-metallicity axis
is Kroupa by the inference above, so none is needed unless the inference is overturned.
