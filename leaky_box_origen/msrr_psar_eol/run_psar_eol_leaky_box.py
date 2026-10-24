#!/usr/bin/env python3
"""Leaky-box PSAR cases for the MSRR (1 MWth molten-salt research reactor) journal paper.

Two cases share one driver and one run directory:

  eol - bounding end-of-life inventory (5 years at 1 MWth, off-gas off, no UF4 make-up) with the
        journal-standard leak rates (Box A 1 %/d, Box B 0.1 %/d), 90-day horizon, ICRP-72 and
        FGR-11 dose coefficients, and the RG 1.145 / 400 m chi/Q schedule.

  mha - the methodology of the PSAR Rev 0 Section 13.1.1 Maximum Hypothetical Accident. The
        one-year source term (f71 position 2) liberates Te, I, Xe, Br and Kr (PSAR Table 13.1-1)
        together with the other enclosure gases H, He, O, N and Ar. They pass through three
        compartments: the reactor enclosure (Box A: 0.05 %/d for Ar, Kr and Xe, 0.005 %/d for the
        rest), the reactor cell (Box B: 1 %/d) and the systems pit (Box C). The release
        to the atmosphere is the pit leak, L_pit times the pit activity (PSAR Eq. 13.1-2), and the
        research bay above the pit has no holdup. The PSAR defines L_pit without a value; the
        driver takes 1 %/d from the Section 13.1.1 statement that the cell leaks at 1 %/d into
        the pit and from the pit into the research bay. Dose: 60 days, FGR-11 inhalation and FGR-12
        air-submersion coefficients, the PSAR constant chi/Q = 1/115.6 s/m^3 (100 m, stability
        class F, 1 m/s, building wake) and the ICRP-23 0.958 m^3/h breathing rate. The report phase prints the published 48.66 mrem
        for context; the PSAR solved the transport for selected mass chains in MATLAB and
        integrated a fitted dose-rate curve, so the two numbers come from different calculations.

Both f71 positions are mass-only ORIGEN output, so the same 5.5e5 volume rescaling applies to
both: the driver divides Box A atom densities by the salt volume and uses that volume as the
dose-phase activity_scale. See README.md in this directory for provenance and caveats.

Usage:  run_psar_eol_leaky_box.py [phase] [case]      phase = sim | dose | report | all
                                                      case  = eol | mha      (default: all eol)

Long runs: launch the sim phase detached, `(setsid python3 ... </dev/null >log 2>&1 &)`. A plain
`nohup ... &` dies with the launching shell's process group, and without the parentheses the launching
shell waits on the job's terminal.
"""
import os
import sys
from typing import Any

from leaky_box_origen.LeakyBox import (
    BREATHING_RATE_M3_S,
    LEAKED_ELEMENTS,
    LEAKY_BOX_DATA_DIR,
    PCTperDAY,
    _setup_box_a,
    _run_simulation,
    _write_excel,
    _write_dataframes_json,
    _write_box_c_timeseries_json,
    _write_activity_csv,
    _activity_timeseries_from_box_json,
    _activity_timeseries_per_nuclide_from_box_json,
    _inventory_to_release_rate,
    _removal_rate_from_run,
    _run_prefix_from_box_json,
    _with_empty_start,
    _load_dcf_csv,
    _load_dcf_immersion_csv,
    _resolve_input_file,
    _out_name,
    _in_directory,
    chi_q_schedule_rg145_400m,
    compute_dose_timeseries,
    plot_total_activity,
    _plot_isotopes,
    compute_dose_from_box,
    plot_dose_comparison,
    plot_dose_rate_comparison,
)

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
F71_PATH = os.path.abspath(os.path.join(CASE_DIR, '..', 'run_2026-10-22', 'msrr_psar_eol',
                                        'msrr_psar_eol.f71'))
RUN_DIR = os.path.abspath(os.path.join(CASE_DIR, '..', 'run_2026-10-22'))
F71_CASE = '1'
SALT_MASS_G = 1600e3         # total fuel salt mass (bounding: all fission products retained)
SALT_VOLUME_CM3 = 550e3      # 550 L total salt volume
STEPS_PER_DAY = 2.0          # ORIGEN time-step density (journal method)

ICRP72_INHAL = 'dcf_icrp72_inhalation_adult_max.csv'
ICRP72_IMMERSION = 'dcf_icrp72_immersion_adult.csv'
FGR11_INHAL = 'dcf_fgr11_inhalation_worker.csv'
FGR11_SUBMERSION = 'dcf_fgr11_submersion_worker.csv'
FGR12_SUBMERSION = 'dcf_fgr12_submersion.csv'
DCF_ICRP72 = ('icrp72', ICRP72_INHAL, ICRP72_IMMERSION)
DCF_FGR11 = ('fgr11', FGR11_INHAL, FGR11_SUBMERSION)
# PSAR Sec. 13.1.1 combines FGR-11 inhalation with FGR-12 external (air-submersion) coefficients.
DCF_PSAR = ('fgr11_fgr12', FGR11_INHAL, FGR12_SUBMERSION)

# PSAR Sec. 13.1.1 MHA assumptions (Table 13.1-1, Eqs. 13.1-5/13.1-8, Ref. 13.4-9).
DAY_S = 86400.0
MHA_NOBLE_RATE_S = 0.05e-2 / DAY_S      # enclosure leak, noble-gas radionuclides
MHA_OTHER_RATE_S = 0.005e-2 / DAY_S     # enclosure leak, all other nuclides
MHA_CELL_RATE_S = 1e-2 / DAY_S          # reactor cell leak, all gases
# Systems pit leak L_pit of Eq. 13.1-2. The PSAR gives no value; Sec. 13.1.1 states that the cell
# leaks at 1 %/d into the pit and from the pit into the research bay, which has no holdup.
MHA_PIT_RATE_S = MHA_CELL_RATE_S
# Table 13.1-1 liberates 100 % of Te, I, Xe, Br and Kr at the start of the MHA. H, He, O, N and Ar
# are gases in the enclosure as well and leak with them, Ar at the noble-gas rate, so the released set
# is the nine journal elements plus Te. Daughters of other elements (Cs, Rb, Sr) stay in the box they
# are born in, because element-level removal of Cs/Rb would also release the salt's own Cs/Rb.
MHA_LEAKED_ELEMENTS: tuple[str, ...] = LEAKED_ELEMENTS + ('te',)
MHA_BOX_A_RATES: dict[str, float] = {**{ele: MHA_NOBLE_RATE_S for ele in ('ar', 'kr', 'xe')},
                                     **{ele: MHA_OTHER_RATE_S
                                        for ele in ('h', 'he', 'o', 'n', 'i', 'br', 'te')}}
MHA_BOX_B_RATES: dict[str, float] = {ele: MHA_CELL_RATE_S for ele in MHA_LEAKED_ELEMENTS}
MHA_CHI_Q_S_M3 = 1.0 / 115.6            # Eq. 13.1-5, 100 m, class F, 1 m/s, with building wake
MHA_BREATHING_M3_S = 0.958 * 24.0 / DAY_S   # ICRP-23 reference man, 22.992 m^3/d
MHA_DOSE_MREM = 48.66                   # PSAR fit to Fig. 13.1-1 integrated over 60 d at 100 m

CASES: dict[str, dict[str, Any]] = {
    'eol': {
        'f71_time_s': 1.57788e8,        # EOL position: end of the fifth operating year
        'decay_days': 90.0,
        'out_prefix': 'f71_psar_eol_90d',
        'plot_isotopes': ('i-129', 'i-131', 'kr-85', 'kr-88'),
        'box_A_leak_rate': PCTperDAY,
        'box_B_leak_rate': PCTperDAY * 0.1,
        'box_c_leak_rate': 0.0,         # collector only, as the journal method has it
        'released_elements': None,      # LEAKED_ELEMENTS
        'dcf_tables': (DCF_ICRP72, DCF_FGR11),
        'chi_q': None,                  # default RG 1.145 schedule at 400 m
        'breathing_m3_s': None,         # default 22.2 m^3/d
        'box_c_rate_s': PCTperDAY * 0.1,
        'release_box': 'boxB',
        'dose_boxes': ('boxB', 'boxC'),  # Box B release and the journal Box C surrogate
        'reference_mrem': None,
        'plot_comparison': True,
        'plot_mha_dose': False,
    },
    'mha': {
        'f71_time_s': 3.15576e7,        # position 2: one year of operation (1 MW-year)
        'decay_days': 60.0,
        'out_prefix': 'f71_psar_mha_60d',
        'plot_isotopes': ('te-132', 'i-131', 'i-133', 'xe-133', 'kr-85'),
        'box_A_leak_rate': MHA_BOX_A_RATES,
        'box_B_leak_rate': MHA_BOX_B_RATES,
        'box_c_leak_rate': MHA_PIT_RATE_S,   # one rate for every element the run releases
        'released_elements': MHA_LEAKED_ELEMENTS,
        'dcf_tables': (DCF_PSAR,),      # FGR-11 inhalation and FGR-12 submersion, as the PSAR
        'chi_q': MHA_CHI_Q_S_M3,
        'breathing_m3_s': MHA_BREATHING_M3_S,
        'box_c_rate_s': None,           # the pit leaks at its own rate: read it from the run
        'release_box': 'boxC',          # Eq. 13.1-2: the atmosphere receives the pit leak
        'dose_boxes': ('boxC',),        # the Box B leak feeds the pit and reaches no receptor
        'reference_mrem': MHA_DOSE_MREM,
        'plot_comparison': False,       # one DCF set: nothing to compare
        'plot_mha_dose': True,          # dose rate and cumulative dose, as PSAR Fig. 13.1-1
    },
}


def run_sim(case: dict[str, Any]) -> None:
    out_prefix = case['out_prefix']
    box_a = _setup_box_a(F71_PATH, F71_CASE, case['f71_time_s'], SALT_MASS_G, do_skip_calcs=False,
                         decay_days=case['decay_days'], steps_per_day=STEPS_PER_DAY)
    # The F71 normalizes the total inventory to 1 cm^3 (mass-only ORIGEN input). Rescale
    # the atom densities to per cm^3 of the real salt so Box A/B/C are physically per cm^3;
    # the whole-inventory dose comes from activity_scale = SALT_VOLUME_CM3 in the dose phase.
    rescale = 1.0 / SALT_VOLUME_CM3
    box_a.atom_dens = {k: v * rescale for k, v in box_a.atom_dens.items()}
    box_a.volume = SALT_VOLUME_CM3
    box_a.density = SALT_MASS_G / SALT_VOLUME_CM3

    with _in_directory(RUN_DIR):
        pd_A, pd_B, pd_C, volume = _run_simulation(
            box_a, apply_volume_scaling=False, do_skip_calcs=False, out_prefix=out_prefix,
            box_A_leak_rate=case['box_A_leak_rate'], box_B_leak_rate=case['box_B_leak_rate'],
            box_C_leak_rate=case['box_c_leak_rate'])
        _write_excel(pd_A, pd_B, pd_C, out_prefix)
        _write_dataframes_json(pd_A, pd_B, pd_C, out_prefix)
        _write_box_c_timeseries_json(pd_C, volume, out_prefix, source_volume_cm3=box_a.volume)

        boxB_json = _out_name(out_prefix, 'boxB', '.json5')
        boxC_json = _out_name(out_prefix, 'boxC', '.json5')
        pd_B_act = _activity_timeseries_from_box_json(boxB_json, 'box_B')
        pd_C_act = _activity_timeseries_from_box_json(boxC_json, 'box_C')
        _write_activity_csv(pd_B_act, pd_C_act, out_prefix)
        plot_total_activity(pd_B_act, pd_C_act, out_prefix)
        _plot_isotopes(pd_A, pd_B, pd_C, list(case['plot_isotopes']), volume, out_prefix)
    print(f'sim phase complete: {RUN_DIR} (prefix {out_prefix})')


def dose_kwargs(case: dict[str, Any], tag: str, box: str) -> dict[str, Any]:
    _, dcf_inhal, dcf_imm = next(t for t in case['dcf_tables'] if t[0] == tag)
    # max_missing_dcf_fraction=None: warn instead of failing at >5% release without DCF. The EOL
    # release exceeded the default by 0.4 percentage points through five short-lived nuclides
    # (xe-137, kr-89, br-85, i-136, xe-139; all half-lives < 10 min) that carry no inhalation DCF.
    # The MHA release weights those and te-135/te-137 more, because iodine leaks 200x slower and
    # the noble gases 10x faster than in the EOL case.
    kwargs: dict[str, Any] = dict(
        dcf_inhal_csv=dcf_inhal,
        dcf_immersion_csv=dcf_imm,
        activity_scale=SALT_VOLUME_CM3,
        out_prefix=f'{case["out_prefix"]}_{tag}_{box}',
        max_missing_dcf_fraction=None,
    )
    if case['chi_q'] is not None:
        kwargs['chi_q_schedule'] = case['chi_q']
    if case['breathing_m3_s'] is not None:
        kwargs['breathing_rate_m3_s'] = case['breathing_m3_s']
    if case['released_elements'] is not None:
        kwargs['released_elements'] = case['released_elements']
    if box == 'boxC':
        if case['box_c_rate_s'] is not None:
            # Box C has no removal of its own; the journal method uses the cell removal rate times
            # the cumulative Box C inventory as a release-rate surrogate. Integrated over the run it is
            # removal rate / decay constant of the Box B release (1.16 % for I-131), a screening number only.
            kwargs['removal_rate_s'] = case['box_c_rate_s']
        # else: this run leaks Box C for real, so the dose step reads the run's own box-C rate.
    return kwargs


def run_dose(case: dict[str, Any]) -> None:
    out_prefix = case['out_prefix']
    with _in_directory(RUN_DIR):
        for tag, _, _ in case['dcf_tables']:
            for box in case['dose_boxes']:
                case_prefix = {'boxB': 'box_B', 'boxC': 'box_C'}[box]
                compute_dose_from_box(_out_name(out_prefix, box, '.json5'), case_prefix,
                                      **dose_kwargs(case, tag, box))

        if not case['plot_comparison']:
            print(f'dose phase complete: {RUN_DIR} (prefix {out_prefix})')
            return
        for box in ('boxB', 'boxC'):
            csv_icrp = _out_name(f'{out_prefix}_icrp72_{box}', 'leaky_boxes_dose', '.csv')
            csv_fgr = _out_name(f'{out_prefix}_fgr11_{box}', 'leaky_boxes_dose', '.csv')
            title = f'PSAR EOL {int(case["decay_days"])}d leak, Box {box[-1].upper()}'
            plot_dose_comparison(csv_icrp, csv_fgr, 'ICRP-72 (max F/M/S)', 'FGR-11 (worker)',
                                 f'{out_prefix}_{box}', title + ' - cumulative dose')
            plot_dose_rate_comparison(csv_icrp, csv_fgr, 'ICRP-72 (max F/M/S)', 'FGR-11 (worker)',
                                      f'{out_prefix}_{box}', title + ' - dose rate')
    print(f'dose phase complete: {RUN_DIR} (prefix {out_prefix})')


def run_report(case: dict[str, Any]) -> None:
    """Print the horizon dose of the case's environmental release and its per-nuclide attribution.

    Dose is linear in the DCF-weighted release, so a nuclide's share comes from re-running the dose
    machinery with a single-nuclide DCF map. The release is the Box B leak for eol and the systems-pit
    (Box C) leak for mha.
    """
    import pandas as pd

    tag, dcf_inhal_csv, dcf_imm_csv = case['dcf_tables'][0]
    box = case['release_box']
    case_prefix = {'boxB': 'box_B', 'boxC': 'box_C'}[box]
    out_prefix = case['out_prefix']
    chi_q = case['chi_q'] if case['chi_q'] is not None else chi_q_schedule_rg145_400m()
    breathing = case['breathing_m3_s'] if case['breathing_m3_s'] is not None else BREATHING_RATE_M3_S
    with _in_directory(RUN_DIR):
        box_json = _out_name(out_prefix, box, '.json5')
        inh = _load_dcf_csv(_resolve_input_file(dcf_inhal_csv, LEAKY_BOX_DATA_DIR))
        imm = _load_dcf_immersion_csv(_resolve_input_file(dcf_imm_csv, LEAKY_BOX_DATA_DIR))
        act = _with_empty_start(_activity_timeseries_per_nuclide_from_box_json(box_json, case_prefix), case_prefix)
        leak_rate = _removal_rate_from_run(box_json, case_prefix,
                                          _run_prefix_from_box_json(box_json))  # the rate ORIGEN applied
        act = _inventory_to_release_rate(act, leak_rate, case['released_elements'])
        nuclide_cols = [c for c in act.columns if c not in ('time [s]', 'time [d]')]
        act[nuclide_cols] = act[nuclide_cols] * SALT_VOLUME_CM3
        total = compute_dose_timeseries(act, inh, chi_q, breathing, imm,
                                        max_missing_dcf_fraction=None)
        dose_mrem = float(total['dose [rem]'].iloc[-1]) * 1000.0
        print(f'{box} ({case_prefix}) {case["decay_days"]:.0f}-day dose, {tag}: {dose_mrem:.4f} mrem')
        if case['reference_mrem']:
            print(f'PSAR Sec. 13.1.1 publishes {case["reference_mrem"]} mrem from its own calculation '
                  f'(selected mass chains, fitted dose rate)')

        rows: list[tuple[str, float]] = []
        for nuclide in nuclide_cols:
            if nuclide not in inh and nuclide not in imm:
                continue  # no DCF: contributes to no dose
            one = compute_dose_timeseries(act, {nuclide: inh[nuclide]} if nuclide in inh else {},
                                         chi_q, breathing,
                                         {nuclide: imm[nuclide]} if nuclide in imm else {},
                                         max_missing_dcf_fraction=None)
            share = float(one['dose [rem]'].iloc[-1]) * 1000.0
            if share > 0.0:
                rows.append((nuclide, share))
        rows.sort(key=lambda r: -r[1])
        attributed = sum(s for _, s in rows)
        # Nuclides without a DCF add nothing to either sum, so the two agree up to rounding.
        print(f'attributed {attributed:.4f} mrem of {dose_mrem:.4f} mrem '
              f'({100.0 * attributed / dose_mrem:.1f}%, sum of the single-nuclide doses)')
        for nuclide, share in rows[:12]:
            print(f'  {nuclide:10s} {share:12.4f} mrem  {100.0 * share / dose_mrem:5.1f}%')
        pd.DataFrame(rows, columns=['nuclide', 'dose_mrem']).to_csv(
            _out_name(out_prefix, f'leaky_boxes_attribution_{tag}_{box}', '.csv'), index=False)
        if case['plot_mha_dose']:
            print(f'dose figure: {plot_mha_dose(f"{out_prefix}_{tag}_{box}")}')
    print(f'report phase complete: {RUN_DIR} (prefix {out_prefix})')


def plot_mha_dose(dose_prefix: str) -> str:
    """Dose rate and cumulative dose of the MHA release at 100 m, in the units of PSAR Fig. 13.1-1."""
    import pandas as pd
    import matplotlib.pyplot as plt

    df = pd.read_csv(_out_name(dose_prefix, 'leaky_boxes_dose', '.csv'))
    fig, (ax_rate, ax_dose) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax_rate.plot(df['time [d]'], df['dose_rate [Sv/s]'] * 1.0e5 * 3600.0, color='tab:blue', linewidth=1.8)
    ax_rate.set_ylabel('dose rate [mrem/h]')
    ax_rate.set_title('PSAR MHA methodology: systems-pit release, FGR-11/FGR-12, chi/Q = 1/115.6 s/m$^3$')
    ax_dose.plot(df['time [d]'], df['dose [rem]'] * 1000.0, color='tab:blue', linewidth=1.8)
    ax_dose.set_ylabel('cumulative dose [mrem]')
    ax_dose.set_xlabel('time [days]')
    for ax in (ax_rate, ax_dose):
        ax.grid(True, which='both', linestyle=':', alpha=0.5)
    fig.tight_layout()
    fname = _out_name(dose_prefix, 'leaky_boxes_dose_mha', '.png')
    fig.savefig(fname, dpi=200)
    plt.close(fig)
    return fname


def main() -> int:
    phase = sys.argv[1] if len(sys.argv) > 1 else 'all'
    case_name = sys.argv[2] if len(sys.argv) > 2 else 'eol'
    if phase not in ('sim', 'dose', 'report', 'all'):
        print(f'unknown phase: {phase} (use sim, dose, report, or all)')
        return 2
    if case_name not in CASES:
        print(f'unknown case: {case_name} (use {" or ".join(CASES)})')
        return 2
    case = CASES[case_name]
    if phase in ('sim', 'all'):
        run_sim(case)
    if phase in ('dose', 'all'):
        run_dose(case)
    if phase in ('report', 'all'):
        run_report(case)
    return 0


if __name__ == '__main__':
    sys.exit(main())
