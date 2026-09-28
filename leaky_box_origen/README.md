# Leaky Box (ORIGEN)

Leaky-box utilities built around SCALE/ORIGEN. This is where leakage tests and
F71-based leakage simulations live.

**Primary entrypoints**
- `LeakyBox.py` is the maintained implementation. Box A leaks into box B, and box B
  leaks into box C, which collects the release.
- `AnalysisRerun.py` rebuilds the spreadsheet of a single-isotope test from its saved
  JSON outputs:
  `python leaky_box_origen/AnalysisRerun.py --run-dir <dir> --prefix xe135|xe136 [--volume 1.0]`.
  It reads `boxA_<prefix>.json5`, `boxB_<prefix>.json5` and `boxC_<prefix>.json5`
  (falling back to the unprefixed names), adds the analytic solutions, and writes
  `leaky_boxes_rerun_<prefix>.xlsx` into the run directory.
- `Case1F71DecaySeries.py` runs a per-timestep decay scan from an F71 file
  (for each case-1 step: 6.5 L sample decayed for 2 days by default), and writes
  `case_decay_activity.csv` + `case_decay_activity.png`. The CSV has the total
  activity and the activity of one tracked nuclide (`tracked_nuclide`,
  `tracked_activity_bq`). The tracked nuclide is F-18 by default (`--track-nuclide`),
  the fluorine activation product in fluoride salt. F-19 is stable.
- `extract_fgr11_dcf.py` regenerates the two FGR-11 DCF CSVs in `data/`.
- `extract_icrp72_dcf.py` regenerates the four ICRP-72 DCF CSVs in `data/`.

**Outputs**
- Runs write all artifacts under `leaky_box_origen/run_YYYY-MM-DD/`.
- ORIGEN working directories (`_box_*`), JSON, Excel, CSV, and plots are all written to that run directory.
- A run with an output prefix gets its own case directories (`_box_A_<prefix>_leak`,
  `_box_B_<prefix>_NNNN_leak`, `_box_C_<prefix>_NNNN`), so the xe135 and xe136 tests
  can share a run directory. Runs without a prefix keep the names `_box_A_leak`,
  `_box_B_NNNN_leak` and `_box_C_NNNN`. The activity readers take the prefix from the
  box JSON name (`boxB_xe135.json5` -> `xe135`) and fall back to the unprefixed names.

**Model conventions**
- Leaked elements: each box removes H, He, O, N, I, Br, Ar, Kr and Xe
  (`LEAKED_ELEMENTS`, one ORIGEN `removal` entry per element). All other elements,
  including decay daughters such as Rb, Sr, Cs and Ba, stay in the box.
  `compute_dose_from_box` converts inventories to release rates for those elements
  only (`released_elements` overrides the list).
- Normalization: boxes B and C are 1 cm^3 volumes fed by the leak from 1 cm^3 of box A
  (default `apply_volume_scaling=False`). Their inventories, activities and the doses
  computed from them are therefore per cm^3 of box A. Pass
  `activity_scale=<box-A volume in cm^3>` to `compute_dose_from_box` for the dose of
  the whole inventory. The factor is written to the `activity_scale [-]` column.
- Missing coefficients: the dose output has the column
  `activity_fraction_without_dcf [-]`, the fraction of the release rate at each time
  step carried by nuclides with neither an inhalation nor an immersion DCF. A warning
  names the largest of them. `max_missing_dcf_fraction` turns the warning into an error
  above a chosen fraction.
- DCF loaders keep values between 1e-20 and 1e-2 (Sv/Bq, or Sv/day per Bq/m^3) and
  drop anything outside that window as an extraction error.

**Dose coefficient tables (`data/`)**

| File | Source | Quantity | Selection per nuclide | Rows |
|---|---|---|---|---|
| `dcf_fgr11_inhalation_worker.csv` | FGR-11 Table 2.1 | h_E,50 [Sv/Bq] (`dcf_sv_bq`) | largest value over the listed lung clearance classes (D/W/Y), vapour and chemical forms | 733 |
| `dcf_fgr11_submersion_worker.csv` | FGR-11 Table 2.3 | h_E,ext [Sv/day per Bq/m^3] (`dcf_sv_per_bq_m3_day`) | tabulated value x 24 (FGR-11 prints Sv/hr per Bq/m^3) | 27 |
| `dcf_icrp72_inhalation_adult.csv` | ICRP-72 Table A.2, adult | e(50) [Sv/Bq] | Type F where listed, otherwise the fastest listed type (M, then S) | 738 |
| `dcf_icrp72_anib26_1_adult.csv` | ICRP-72 Table A.2, adult | e(50) [Sv/Bq] | Type M where listed, otherwise the fastest listed type | 738 |
| `dcf_icrp72_inhalation_adult_max.csv` | ICRP-72 Table A.2, adult | e(50) [Sv/Bq] | largest value over the listed types (F/M/S) | 738 |
| `dcf_icrp72_immersion_adult.csv` | ICRP-72 Table A.4 | effective dose rate [Sv/day per Bq/m^3] | tabulated value | 26 |

- Sources: EPA Federal Guidance Report No. 11, EPA-520/1-88-020 (1988), and ICRP
  Publication 72, Annals of the ICRP 26(1) (1996). ANIB_26_1 is the file name of that
  Annals issue.
- FGR-11 values are for workers (ICRP-30 lung model). ICRP-72 values are for adult
  members of the public (ICRP-66 lung model). Mercury in ICRP-72 lists organic and
  inorganic forms; the larger value of each type is kept.
- The largest-over-types ICRP table is the conservative choice when the chemical form
  of the release is unknown. The Type F table has lower values for insoluble forms.
  For Cs-137 the three ICRP tables give 4.6e-9, 9.7e-9 and 3.9e-8 Sv/Bq.
- H-3: FGR-11 inhalation is tritiated water vapour (class V, 1.73e-11 Sv/Bq), and
  FGR-11 submersion is elemental tritium (2.86e-14 Sv/day per Bq/m^3). The three ICRP
  inhalation tables use tritiated water vapour from Table A.3 (1.8e-11 Sv/Bq). Table
  A.2 "tritium compounds" (particulates) are not used. ICRP-72 Table A.4 has no tritium.
- ICRP-72 Table A.3 lists other gas and vapour forms (elemental iodine, methyl iodide,
  carbon compounds and others). Apart from HTO they are not in these tables.
- FGR-11 separates some isomers by half-life instead of an `m` suffix. The CSVs name
  them after the ground or metastable state: Nb-89 (122 min) is `nb-89` and Nb-89
  (66 min) is `nb-89m`; the same holds for In-110, Sb-120, Sb-128, Eu-150 and
  Tb-156m (`tb-156m`, 24.4 h, and `tb-156m2`, 5.0 h).
- FGR-11 and ICRP-72 inhalation values agree within a factor of 30 for all 733 shared
  nuclides when the largest class and the largest type are compared. Submersion and
  immersion values agree within 0.5--1.25 for 23 of the 26 noble gases. Ar-37, Ar-39
  and Kr-83m differ by factors of 0.12--2.7. These three emit little penetrating
  radiation, so the skin, lens and lung dose conventions set their coefficients, and
  FGR-11 h_E excludes skin and lens.

**Regenerating the tables**
- The source PDFs are not distributed with the repository. Place them under `PDF/`
  (gitignored) or pass their paths.
- FGR-11: `python -m leaky_box_origen.extract_fgr11_dcf [--pdf <FGR-11 PDF> ...] [--audit <dir>]`.
  It needs `pdftoppm`, `pdftotext`, `pdfinfo` and `tesseract`, plus numpy, scipy and
  Pillow. It reads the OSTI copy (`6294233.pdf`) as the primary scan and the EPA copy
  (`EPA 1988_FGR11_0.pdf`) as a second reading. Both PDF text layers are OCR output
  whose superscript exponents are unreadable, so the extractor locates each
  `m.mm 10^-e` value on the page image and reads the exponent digits separately. The
  text-layer mantissas are used as additional votes. Every row must satisfy the
  definition of the effective value, h_E = sum_T w_T h_T with the ICRP-26 weights, to
  within 1.2%. Rows the automatic pass cannot settle are read by eye from the page
  images and listed with their page number in `MANUAL_CORRECTIONS` (17 of 1351 Table
  2.1 rows).
- ICRP-72: `python -m leaky_box_origen.extract_icrp72_dcf [--pdf PDF/ANIB_26_1.pdf] [--audit <dir>]`.
  It needs `pdftoppm` and `tesseract`, plus numpy, scipy and Pillow. Each Table A.2 page
  is read five times (three full-page readings and two readings of the value columns).
  The majority value is kept, with the 15-year column as a plausibility check. Values
  read by eye are listed in `VALUE_FIXES` (26 of 1683 rows), and OCR misreads of mass
  numbers identified by the printed half-life in `NAME_FIXES`.
- `--audit <dir>` writes the per-row readings, the weighted-sum check and the source of
  every value. Review the output before using the DCFs for licensing decisions.
- `test/test_dcf_tables.py` checks anchor nuclides of every table against the published
  values within a factor of 1.5, and rejects impossible nuclide names.

**Notes**
- Defaults (F71 path, case/time selection, leak elements) are set in `LeakyBox.py`.
- Override F71 defaults with environment variables:
  - `LEAKYBOX_F71_PATH`
  - `LEAKYBOX_F71_CASE`

Example:

```bash
python leaky_box_origen/Case1F71DecaySeries.py --f71 /path/to/msrr.f71 --case 1
```
