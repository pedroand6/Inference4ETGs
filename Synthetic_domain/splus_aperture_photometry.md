# S-PLUS 3" aperture photometry for the SDSS/eBOSS Firefly reference galaxies
- Data release: S-PLUS DR4 (dual-mode catalogue), schema `dr4_dual`. Values are RAW (no Galactic-extinction correction); non-detections kept as returned (~99).
- Position match table: `dr4_dual.dr4_dual_detection` (columns `ID`, `RA`, `DEC`). Output `splus_ra`/`splus_dec` are its `RA`/`DEC`.
- Photometry tables (one per band, joined on `ID`): `dr4_dual.dr4_dual_{u,j0378,j0395,j0410,j0430,g,j0515,j0660,j0861,r,i,z}`.
- Columns per band: `<band>_aper_3`, `e_<band>_aper_3` with band in u, J0378, J0395, J0410, J0430, g, J0515, r, J0660, i, J0861, z (same names as the output CSV).
- Method: uploaded-table crossmatch (`upload.<name>`, 1=CONTAINS(POINT(RA,DEC), CIRCLE(qra,qdec,1.5/3600))), 10 chunks of 150 positions on the detection table, then ID-join upload per band. Closest S-PLUS object kept (match_separation_arcsec).
- Input list: dr3_corrected.csv (PROB_GAL>0.8) x Data/sdss_eboss_firefly-dr16.fits (CLASS_NOQSO=GALAXY, 0.002<Z_NOQSO<0.7), nearest within 1.5", one spectrum per S-PLUS galaxy; `specobjid` is a string.
- Counts: 1362 galaxies queried, 1360 with photometry (2 have no DR4 object within 1.5"; those rows are NaN), 0 chunks failed.
- Run time: ~17 min for 4 parallel detection chunks x 3 waves (~50 min total, server ~5 s/row) plus ~15 min for the 12 band queries; wall clock about 1h15 including a failed first attempt.
- Note: DR4 aperture mags differ from DR3 (`dr3.all_dr3`, same column names) by roughly 0.05-0.1 mag (recalibration); dr3_corrected.csv is DR3-based. Many 99s in the blue bands: u 195, J0378 215, J0395 246, J0410 202, J0430 155, J0515 58, g 16; ~0 in r, J0660, i, J0861, z.
