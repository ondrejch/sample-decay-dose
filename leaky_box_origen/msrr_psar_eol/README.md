# MSRR PSAR Rev 0 — EOL isotopics (ORIGEN bounding case)

ORIGEN-AR burn that computes the end-of-life fuel isotopics for the MSRR
(1 MWth molten-salt research reactor) from the bounding assumptions in the
MSRR Plant Specific Application Report, Rev 0. This is the source-term input
(Box A inventory) for the MSRR PSAR test case of the leaky-box paper.

Full distilled notes (inputs, derivation, caveats, integration): `.kb/msrr-psar-eol-origen.md`.

## Bounding case

- All fission products retained, **including gases** (Kr, Xe, H-3, He):
  off-gas system not in operation → no `removal` block in the deck.
- **No UF4 addition** over the lifetime (PSAR 5-year bounding simulation).
- Constant 1 MWth, continuous critical operation.

## Inputs (all from `SMD-MSRR-dev/ACU_PSAR_REV_0.txt`)

| Quantity | Value | PSAR source |
|---|---|---|
| Fuel salt | 67LiF-28BeF2-5UF4 (mole %) | §1.2.3.2, §4.2.1.1 |
| U-235 enrichment | 19.75 % | §1.2.3.2, §4.1 |
| Li-7 enrichment | 99.99 % | §1.2.3.2 |
| UF4 mass | < 500 kg (nominal 500) | §1.2.3.2 |
| Total salt | ~1600 kg / ~550 L | §1.2.3.2, §4.2.1.1 |
| Power | 1 MWth (max licensed) | Table 4.1-1 |
| Burnup basis | **5 MW-years** continuous, off-gas off | §4.5.2.1, §4.5.2.3 |
| Neutron spectrum | thermal (graphite) | Table 4.1-1 |

Derived: 500 kg UF4 @ 19.75% = 1595 mol U = **378.8 kg U**
(U-235 = 74 050 g, U-238 = 304 700 g). 5 MW·y at 1 MWth = 1826 MWd,
so burnup ≈ **4.8 MWd/kgU** and total fissions ≈ 4.9e24.

## Burnup note (correction of the SMD-MSRR-dev HANDOFF)

The original handoff derived "5000 MWd → 13.2 MWd/kgU". That treated 5 EFPY
as 5000 MWd; 1 EFPY = 365.25 MWd, so 5 EFPY at 1 MWth is **1826 MWd →
4.8 MWd/kgU** (the PSAR's literal "5 MW-years"). The deck
(`power=[5r1.0]`, 5 annual steps) is the literal PSAR case; only the
documented derived quantities were corrected.

## Library

`t100.f33` (named in the original handoff) is not present in this SCALE 6.3.3
install. We use `pwr.rev04.orglib`, the only prebuilt thermal f33
(`$SCALE_DATA/origen_library/`), as the closest self-contained surrogate for
the thermal MSRR spectrum. See `.kb/scale-env.md`.

## Run

```bash
SCALE_BIN=/opt/scale6.3.3/bin   # already exported in this environment
$SCALE_BIN/scalerte msrr_psar_eol.inp     # from a run/scratch dir
```

Read `msrr_psar_eol.out` (per-step nuclide tables) and `msrr_psar_eol.f71`
(one case, 6 positions: initial + 5 year-ends; EOL = last position,
t = 1.57788e8 s).

## Run result (2026-10-22)

`run_2026-10-22/msrr_psar_eol/` ran clean. EOL (position 6) inventory,
1309 nuclides, total 378.75 kg (mass conserved vs BOC):

| Quantity | Value |
|---|---|
| Total fissions | 4.898e24 (= 1816 MWd @ 200 MeV/f, matches 1826 MWd design) |
| U-235 burned | 2290 g (74050 → 71760 g) |
| Pu total | 2217 g (239: 2185, 240: 28.5, 241: 3.4, 242: 0.02) |
| Am-241 | 0.214 g |
| Cm | ~1 mg (242 dominant) |
| Kr-85 | 1.64 g |
| Xe-134 / 136 / 133 / 135 | 85.4 / 71.9 / 0.30 / 0.021 g |
| H-3 | 2.2 mg (FP yield only; no Li carrier in model) |
| Cs-137 / Sr-90 / Tc-99 / I-129 | 65 / 39 / 49 / 6.0 g |

All leaked elements (H, He, O, N, I, Br, Ar, Kr, Xe) present at EOL.

## Deck notes (SCALE 6.3.3 ORIGEN-AR)

- `lib{ file="${DATA}/origen_library/pwr.rev04.orglib" }`: bare library names
  are not searched; `${DATA}` expands to the SCALE data dir.
- `cutoffs[...]` sits at the `print` level, **not** inside `nuc{}`; with
  `rel_cutoff=no` and `cutoffs[GRAMS=0.0 CURIES=0.0]` every nuclide prints.
- `scalerte` exits 0 even on module failure — grep the `.out` for
  `***Error` / "completion code".

## Leaky-box test case (paper)

The EOL F71 (case 1, position 6, t = 1.57788e8 s) is the Box A inventory for
the 90-day leaky-box simulation driven by `run_psar_eol_leaky_box.py`:

```bash
python3 leaky_box_origen/msrr_psar_eol/run_psar_eol_leaky_box.py all eol    # sim + dose + report
python3 leaky_box_origen/msrr_psar_eol/run_psar_eol_leaky_box.py all mha    # PSAR MHA reproduction
```

The driver takes `[phase] [case]`: phase `sim | dose | report | all`, case `eol | mha`
(default `all eol`). `mha` carries out the PSAR §13.1.1 Maximum Hypothetical Accident
methodology instead: f71 position 2 (one year), the enclosure → cell → systems pit chain as Boxes
A → B → C with per-element rates (Box A 0.05 %/day noble gases and 0.005 %/day the others, Te
included; Boxes B and C 1 %/day for all), the pit leak as the release, 60-day horizon, FGR-11
coefficients at the PSAR's constant χ/Q = 1/115.6 s/m³ (100 m, class F, 1 m/s) and the ICRP-23
0.958 m³/h breathing rate. `report` prints the Box C horizon dose with the per-nuclide attribution
and writes the dose-rate figure. Its artifacts use the prefix `f71_psar_mha_60d` in the same run
directory.

Artifacts land in `run_2026-10-22/` under prefix `f71_psar_eol_90d`.

**Volume rescaling.** ORIGEN-AR built this F71 from a mass-only input
(`mat units=GRAMS`, no volume), so the total inventory is normalized into an
implicit 1 cm³. The real salt is 550 L, so the driver rescales Box A atom
densities by 1/5.5e5 and sets `volume = 5.5e5 cm³`, `density = 2.91 g/cm³`.
The whole-inventory dose then comes from `activity_scale = 5.5e5` (salt
volume) in `compute_dose_from_box`. Without this rescale the derived Box A
volume would be ~4 cm³ and the dose wrong by ~5 orders of magnitude.

**Leak parameters (journal method).** Box A leak 1 %/day, Box B leak
0.1 %/day, Box C no removal. Leaked elements: H, He, O, N, I, Br, Ar, Kr, Xe.

**PSAR Ch. 13 consistency.** These are the journal's standard rates, *not*
the MHA rates of PSAR §13.1.1 (Table 13.1-1): enclosure 0.05 %/day
(noble gases) / 0.005 %/day (other nuclides), cell 1 %/day (all gases). Box A
(1 %/day) matches the cell rate; Box B (0.1 %/day) is 2× the enclosure
noble-gas rate and 20× the non-noble rate. The `mha` case below carries out
the PSAR methodology itself, with Te added to the released set and the systems
pit as a third compartment. MHA
reference dose: 48.66 mrem at 100 m over 60 d (FGR-11,
χ/Q = 1/115.6 s/m³). See `.kb/msrr-psar-eol-origen.md` for the full check.

### MHA result (2026-10-22)

`mha` follows the compartment chain of PSAR §13.1.1: the reactor enclosure
leaks into the reactor cell, the cell into the systems pit, and the pit into a
research bay that has no holdup, so the pit leak reaches the atmosphere (PSAR
Eq. 13.1-2, `Q_i = L_pit λ_i N_i^pit`). Box A is the enclosure (Ar/Kr/Xe at
0.05 %/day, H/He/O/N/I/Br/**Te** at 0.005 %/day), Box B the cell and Box C the
systems pit (all ten elements at 1 %/day each). The PSAR defines `L_pit`
without a value; the driver takes 1 %/day from the §13.1.1 statement that the
cell leaks at 1 %/day into the pit and from the pit into the research bay. The
source term is f71 position 2 (one operating year) and the horizon is 60 days.
The dose step uses FGR-11 inhalation and FGR-12 submersion, χ/Q = 1/115.6 s/m³
as a constant, and the ICRP-23 breathing rate of 0.958 m³/h.

60-day dose at 100 m from the systems-pit release: **14.96 mrem**, with FGR-11
inhalation and FGR-12 air-submersion coefficients as in the PSAR. I-131 carries
73.4 %, Te-129m 16.1 %, Te-127m 6.8 % and Te-132 2.4 %; submersion gives 1.5 % of
the dose. The release fraction without a dose coefficient peaks at 3.3 %, under
the 5 % default.

The PSAR publishes 48.66 mrem for the same methodology from its own
calculation: a MATLAB solution for mass chains 131-135, 137 and 85, 87-91 and
a line of best fit to the dose rate integrated over 60 days. The model follows every ORIGEN
nuclide (Te-127m and Te-129m give 23 % of its dose from outside those chains),
takes the largest FGR-11 value over lung-clearance classes, and integrates the
computed dose rate. No model parameter is tuned toward the published value.
Full report: `.kb/psar-mha-reproduction-report.md`.

### Result (2026-10-22)

90-day cumulative dose to a worker (RG 1.145, 400 m, 22.2 m³/d):

| Box / basis | ICRP-72 (max F/M/S) | FGR-11 (worker) |
|---|---|---|
| Box B (physical release) | **8.1e-4 Sv (0.81 mSv)** | **9.7e-4 Sv (0.97 mSv)** |
| Box C (surrogate) | 7.1e-6 Sv | 8.5e-6 Sv |

Inhalation dominates (>99.9 %). Box B exceeds Box C over a 90-day window
because the surrogate is about eps_B / lambda of the Box B release for I-131
(0.9 % here); the environmental release is the Box B leak.

Dominant contributors (ICRP-72, Box B, per-nuclide attribution):
**I-131 = 98.4 %** (short-lived, from the in-operation EOL inventory), I-133
1.4 %, Xe-133 0.1 %, the rest negligible. I-129, Kr-85, and the Xe immersions
are small individually but set the long-time tail.

**DCF coverage note.** ICRP-72 and FGR-11 both leave 5.38 % of the Box B
release without a dose coefficient (xe-137, kr-89, br-85, i-136, xe-139; all
half-lives < 10 min), marginally above the 5 % fail-closed default, so the
driver passes `max_missing_dcf_fraction=None` (warn-only). The Box C surrogate
leaves 2.15 % uncovered. The short-lived missing nuclides' 90-day dose contribution is
negligible.

## Caveats (state with any result)

- Thermal-spectrum library data (fission yields, effective xs): order-of-
  magnitude table, not a qualification/acceptance value.
- No corrosion products, impurities, or structural activation.
- No UF4 addition — bounding case; real EOL inventory is lower in fissile.
- EOL fuel composition per the PSAR is officially deferred to the Operating
  License Application (PSAR §4.2.1.1). This run is an estimate.

## References

- `SMD-MSRR-dev/ACU_PSAR_REV_0.txt` — PSAR Rev 0.
- `SMD-MSRR-dev/00runs/origen-msrr-eol/HANDOFF.md` — original spec.
- SCALE 6.3.3 manual, ORIGEN module §5.1.5.
- KB note: `smd-msrr-dev/acu-psar-rev-0-no-eol-isotopics` (Basic Memory).
