#!/usr/bin/env python3
"""Numbers of the MSRR applications (Part II) of the leaky-box paper, from the existing PSAR runs.

Reads `run_2026-10-22/` (the EOL and MHA runs of `msrr_psar_eol/run_psar_eol_leaky_box.py` and the
`msrr_psar_eol.f71` inventory) with obiwan; it runs no ORIGEN case. Sections, each printing what the paper quotes:

  totals          60-day MHA and 90-day EOL doses of the runs, integrated from t = 0
  attribution     per-nuclide 60-day MHA doses of the cell exit (Box B) and the systems-pit release (Box C), with
                  FGR-11 + FGR-12 and with FGR-11 alone, and their restriction to the PSAR mass chains
  two-hour        2-hour doses of the 90-day EOL case, integrated over the Box B and Box C sub-run positions
  positions       dose-driver activities at F71 positions 2 and 6 and the share-weighted dose change
  pit-rate        systems-pit leak rate that reproduces the PSAR 48.66 mrem, from parent-free compartment chains
  daughters       bound on the Cs, Rb and Sr daughters that the PSAR chains transport and this model retains
  hand-check      analytic three-compartment I-131 dose, with I-131 bred from Te-131m and the A = 131 precursors
  iodine-vapour   EOL and MHA doses with the ICRP-72 elemental iodine vapour coefficients

The numbers also go to `paper_numbers_psar.json` and the attribution to `paper_mha_attribution.csv` in the run
directory.

Usage:  python -m leaky_box_origen.paper.psar_numbers [section ...]      (default: all sections)
"""
import json
import math
import os
import re
import sys
import warnings

import json5
import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.optimize import brentq

from leaky_box_origen.LeakyBox import (
    BREATHING_RATE_M3_S, ICRP72_IODINE_VAPOUR_CSV, LEAKY_BOX_DATA_DIR,
    _activity_timeseries_per_nuclide_from_box_json, _element_of, _find_case_f71, _inventory_to_release_rate,
    _load_dcf_csv, _load_dcf_immersion_csv, _out_name, _removal_rate_from_run, _resolve_input_file,
    _step_activity_from_f71, _step_case_dirs, _with_empty_start, chi_q_schedule_rg145_400m, compute_dose_timeseries,
)
from leaky_box_origen.msrr_psar_eol import run_psar_eol_leaky_box as drv
from sample_decay_dose.utils import get_burned_nuclide_atom_dens, get_burned_nuclide_data, get_f71_positions_index

DAY_S = 86400.0
MREM_PER_SV = 1.0e5
RUN_DIR = drv.RUN_DIR
MHA = drv.CASES['mha']
EOL = drv.CASES['eol']
# PSAR Sec. 13.1.1 follows the Te-I-Xe-Cs chains A = 131-135, 137 and the Br-Kr-Rb-Sr chains A = 85, 87-91.
PSAR_CHAINS = frozenset((131, 132, 133, 134, 135, 137, 85, 87, 88, 89, 90, 91))
PSAR_MHA_MREM = drv.MHA_DOSE_MREM
# Short-lived daughters carry the holdup of the parent that regenerates them in each box.
PARENT_OF = {'i-132': 'te-132', 'te-127': 'te-127m', 'te-129': 'te-129m', 'te-131': 'te-131m'}
NOBLE = ('ar', 'kr', 'xe')
# Atomic numbers along the beta-decay chains of the liberated elements.
Z_OF = {'br': 35, 'kr': 36, 'rb': 37, 'sr': 38, 'y': 39, 'sb': 51, 'te': 52, 'i': 53, 'xe': 54, 'cs': 55, 'ba': 56}
SYMBOL_OF = {z: s for s, z in Z_OF.items()}

results: dict = {}
_release_cache: dict = {}


def _dcf(name: str) -> dict[str, float]:
    path = _resolve_input_file(name, LEAKY_BOX_DATA_DIR)
    return _load_dcf_immersion_csv(path) if 'submersion' in name or 'immersion' in name else _load_dcf_csv(path)


INH_FGR11 = _dcf(drv.FGR11_INHAL)
SUB_FGR11 = _dcf(drv.FGR11_SUBMERSION)
SUB_FGR12 = _dcf(drv.FGR12_SUBMERSION)
INH_ICRP72 = _dcf(drv.ICRP72_INHAL)
IMM_ICRP72 = _dcf(drv.ICRP72_IMMERSION)
IODINE_VAPOUR = _dcf(ICRP72_IODINE_VAPOUR_CSV)
DCF_SETS = {'fgr11_fgr12': (INH_FGR11, SUB_FGR12), 'fgr11': (INH_FGR11, SUB_FGR11),
            'icrp72': (INH_ICRP72, IMM_ICRP72)}


def _mass_number(nuclide: str) -> int:
    return int(re.sub(r'\D', '', nuclide.split('-')[1]))


def _nuclide_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c not in ('time [s]', 'time [d]')]


def _chi_q(case: dict):
    return case['chi_q'] if case['chi_q'] is not None else chi_q_schedule_rg145_400m()


def _breathing(case: dict) -> float:
    return case['breathing_m3_s'] if case['breathing_m3_s'] is not None else BREATHING_RATE_M3_S


def release(case: dict, box: str, rate: float | dict | None = None) -> pd.DataFrame:
    """ Whole-inventory release rate [Bq/s] of a box: its own removal rates from the run unless `rate` is given """
    key = (case['out_prefix'], box, rate)
    if key in _release_cache:
        return _release_cache[key].copy()
    case_prefix = {'boxB': 'box_B', 'boxC': 'box_C'}[box]
    box_json = os.path.join(RUN_DIR, _out_name(case['out_prefix'], box, '.json5'))
    act = _with_empty_start(_activity_timeseries_per_nuclide_from_box_json(
        box_json, case_prefix, run_prefix=case['out_prefix']), case_prefix)
    if rate is None:
        rate = _removal_rate_from_run(box_json, case_prefix, case['out_prefix'])
    rel = _inventory_to_release_rate(act, rate, case['released_elements'])
    cols = _nuclide_cols(rel)
    rel[cols] = rel[cols] * drv.SALT_VOLUME_CM3
    _release_cache[key] = rel
    return rel.copy()


def dose_mrem(rel: pd.DataFrame, case: dict, inh: dict, sub: dict, until_s: float | None = None) -> float:
    """ Cumulative dose [mrem] of a release frame at its last row, or at `until_s` (release linear in between) """
    if until_s is not None and until_s not in set(rel['time [s]']):
        t = rel['time [s]'].to_numpy(dtype=float)
        row = {c: float(np.interp(until_s, t, rel[c].to_numpy(dtype=float))) for c in rel.columns}
        rel = pd.concat([rel, pd.DataFrame([row])], ignore_index=True).sort_values('time [s]')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out = compute_dose_timeseries(rel, inh, _chi_q(case), _breathing(case), sub, max_missing_dcf_fraction=None)
    if until_s is not None:
        out = out[out['time [s]'] <= until_s]
    return float(out['dose [rem]'].iloc[-1]) * 1000.0


def per_nuclide_mrem(rel: pd.DataFrame, case: dict, inh: dict, sub: dict) -> pd.Series:
    """ Single-nuclide doses [mrem]; the dose is linear in each nuclide's release, so they add to the total """
    doses = {}
    for n in _nuclide_cols(rel):
        if n in inh or n in sub:
            one = rel[['time [s]', n]]
            doses[n] = dose_mrem(one, case, {n: inh[n]} if n in inh else {}, {n: sub[n]} if n in sub else {})
    return pd.Series(doses).sort_values(ascending=False)


def f71_activity_and_lambda(position: int) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    """ Whole-inventory atoms, activities [Bq] and decay constants [1/s] at an F71 position (1 cm^3 normalization) """
    atoms = {n: v * 1.0e24 for n, v in get_burned_nuclide_atom_dens(drv.F71_PATH, position).items()}
    becq = get_burned_nuclide_data(drv.F71_PATH, position, f71units='becq')
    lam = {n: becq[n] / atoms[n] for n in becq if atoms.get(n, 0.0) > 0.0}
    return atoms, becq, lam


def chain_release_fraction(lam: float, eps_a: float, eps_b: float, eps_c: float, horizon_s: float) -> float:
    """ Activity released from Box C over the horizon per initial Box A atom [Bq s], three-box parent-free chain """
    m = np.zeros((4, 4))
    m[0, 0] = -(lam + eps_a)
    m[1, 0], m[1, 1] = eps_a, -(lam + eps_b)
    m[2, 1], m[2, 2] = eps_b, -(lam + eps_c)
    m[3, 2] = eps_c * lam
    return float((expm(m * horizon_s) @ np.array([1.0, 0.0, 0.0, 0.0]))[3])


def mha_rates(element: str) -> tuple[float, float, float]:
    eps_a = drv.MHA_NOBLE_RATE_S if element in NOBLE else drv.MHA_OTHER_RATE_S
    return eps_a, drv.MHA_CELL_RATE_S, drv.MHA_PIT_RATE_S


# ------------------------------------------------------------------------------------------------------------------

def section_totals() -> None:
    for name, case, box, tag in (('mha_pit', MHA, 'boxC', 'fgr11_fgr12'), ('mha_pit', MHA, 'boxC', 'fgr11'),
                                 ('mha_cell', MHA, 'boxB', 'fgr11_fgr12'), ('eol_boxB', EOL, 'boxB', 'icrp72'),
                                 ('eol_boxB', EOL, 'boxB', 'fgr11')):
        rel = release(case, box)
        d = dose_mrem(rel, case, *DCF_SETS[tag])
        results[f'{name}_{tag}_mrem'] = d
        print(f'{name:9s} {tag:12s} {d:12.4f} mrem')
    for tag in ('icrp72', 'fgr11'):
        rel = release(EOL, 'boxC', rate=EOL['box_c_rate_s'])
        d = dose_mrem(rel, EOL, *DCF_SETS[tag])
        results[f'eol_boxC_surrogate_{tag}_mrem'] = d
        print(f'eol_boxC  {tag:12s} {d:12.4f} mrem (surrogate eps_B * A_C), '
              f'{100.0 * d / results[f"eol_boxB_{tag}_mrem"]:.3f} % of Box B')


def section_attribution() -> pd.DataFrame:
    table = {}
    for box in ('boxB', 'boxC'):
        rel = release(MHA, box)
        for tag in ('fgr11_fgr12', 'fgr11'):
            table[f'{box}_{tag}'] = per_nuclide_mrem(rel, MHA, *DCF_SETS[tag])
    df = pd.DataFrame(table).fillna(0.0).sort_values('boxC_fgr11_fgr12', ascending=False)
    df['psar_chain'] = [_mass_number(n) in PSAR_CHAINS for n in df.index]
    df['attenuation_B_over_C'] = df['boxB_fgr11_fgr12'] / df['boxC_fgr11_fgr12']
    df.index.name = 'nuclide'
    df.to_csv(os.path.join(RUN_DIR, 'paper_mha_attribution.csv'))
    print(df.head(12).to_string(float_format=lambda v: f'{v:.4f}'))
    cols = ['boxC_fgr11_fgr12', 'boxC_fgr11', 'boxB_fgr11_fgr12', 'boxB_fgr11']
    top = df.index[:10]
    rows = {'all others': df.loc[~df.index.isin(top), cols].sum(), 'total': df[cols].sum(),
            'PSAR chains': df.loc[df['psar_chain'], cols].sum(), 'outside chains': df.loc[~df['psar_chain'], cols].sum()}
    summary = pd.DataFrame(rows).T
    summary['attenuation_B_over_C'] = summary['boxB_fgr11_fgr12'] / summary['boxC_fgr11_fgr12']
    print(summary.to_string(float_format=lambda v: f'{v:.4f}'))
    te = df.loc[[n for n in df.index if n.startswith('te-')], 'boxC_fgr11_fgr12'].sum()
    total = summary.loc['total', 'boxC_fgr11_fgr12']
    print(f'tellurium in the pit release: {te:.4f} mrem, {100.0 * te / total:.2f} %')
    results.update({'mha_attribution_totals': summary.to_dict(orient='index'), 'mha_tellurium_mrem': te})
    return df


def section_two_hour() -> None:
    """ D_2h of the 90-day EOL case from the sub-run positions of Box A steps 2 (0 to 0.01 d) and 3 (to 0.51 d) """
    for box, case_prefix in (('boxB', 'box_B'), ('boxC', 'box_C')):
        with open(os.path.join(RUN_DIR, _out_name(EOL['out_prefix'], box, '.json5'))) as f:
            step_end = {int(k): float(v['time']) for k, v in json5.load(f).items()}
        rows: dict[float, dict] = {}
        for step in (2, 3):
            f71 = _find_case_f71(_step_case_dirs(case_prefix, step, EOL['out_prefix']), search_roots=[RUN_DIR])
            index = get_f71_positions_index(f71)
            # F71 times count from the start of each sub-run; the box json holds the absolute end time.
            last = max(float(info['time']) for info in index.values())
            for pos, info in index.items():
                rows[step_end[step] - last + float(info['time'])] = _step_activity_from_f71(f71, pos)
        act = pd.DataFrame([{'time [s]': t, 'time [d]': t / DAY_S, **a} for t, a in sorted(rows.items())]).fillna(0.0)
        rate = EOL['box_B_leak_rate']  # eps_B: the Box B leak, and the journal surrogate rate for Box C
        rel = _inventory_to_release_rate(act, rate, EOL['released_elements'])
        cols = _nuclide_cols(rel)
        rel[cols] = rel[cols] * drv.SALT_VOLUME_CM3
        spacing = np.diff(rel['time [s]'].to_numpy())
        for tag in ('icrp72', 'fgr11'):
            d2h = dose_mrem(rel, EOL, *DCF_SETS[tag], until_s=7200.0) / 1000.0
            results[f'eol_{box}_{tag}_D2h_rem'] = d2h
            print(f'{box} {tag:6s} D_2h = {d2h:.3e} rem (sub-run spacing {spacing.min():.0f} to {spacing.max():.0f} s)')


def section_positions(attribution: pd.DataFrame | None) -> None:
    _, a2, _ = f71_activity_and_lambda(2)
    _, a6, _ = f71_activity_and_lambda(6)
    for n in ('i-131', 'i-133', 'xe-133', 'te-132', 'te-129m', 'te-127m', 'kr-85', 'i-129', 'cs-137', 'h-3'):
        print(f'{n:8s} position 6 / position 2 = {a6[n] / a2[n]:.4f}')
    if attribution is None:
        attribution = pd.read_csv(os.path.join(RUN_DIR, 'paper_mha_attribution.csv'), index_col='nuclide')
    shares = attribution['boxC_fgr11_fgr12']

    def weighted(doses: pd.Series, num: dict, den: dict) -> float:
        # First-order change: each nuclide's dose scales with its own source activity (daughters with the parent).
        ratio = {n: num[PARENT_OF.get(n, n)] / den[PARENT_OF.get(n, n)] for n in doses.index
                 if PARENT_OF.get(n, n) in num and PARENT_OF.get(n, n) in den}
        return 100.0 * sum(doses[n] * (r - 1.0) for n, r in ratio.items()) / doses.sum()

    mha = weighted(shares, a6, a2)
    rel = release(EOL, 'boxB')
    eol = weighted(per_nuclide_mrem(rel, EOL, *DCF_SETS['icrp72']), a2, a6)
    results.update({'mha_dose_change_pos6_over_pos2_pct': mha, 'eol_dose_change_pos2_over_pos6_pct': eol})
    print(f'MHA dose with position 6 instead of 2: {mha:+.2f} %; EOL dose with position 2 instead of 6: {eol:+.2f} %')


def section_pit_rate(attribution: pd.DataFrame | None) -> None:
    if attribution is None:
        attribution = pd.read_csv(os.path.join(RUN_DIR, 'paper_mha_attribution.csv'), index_col='nuclide')
    _, _, lam = f71_activity_and_lambda(2)
    doses = attribution['boxC_fgr11_fgr12']
    doses = doses[doses > 1e-6 * doses.sum()]
    horizon = MHA['decay_days'] * DAY_S

    def nuclide_lambda(n: str) -> float:
        return lam.get(PARENT_OF.get(n, n), lam.get(n, 0.0))

    def estimate(pit_pct_per_day: float, chains_only: bool = False) -> float:
        eps_c = pit_pct_per_day * 1e-2 / DAY_S
        total = 0.0
        for n, d in doses.items():
            if chains_only and _mass_number(n) not in PSAR_CHAINS:
                continue
            src = PARENT_OF.get(n, n)
            eps_a, eps_b, eps_1 = mha_rates(_element_of(src))
            lam_n = nuclide_lambda(n)
            total += d * (chain_release_fraction(lam_n, eps_a, eps_b, eps_c, horizon)
                          / chain_release_fraction(lam_n, eps_a, eps_b, eps_1, horizon))
        return total

    open_pit = estimate(1.0e4)
    results['pit_rate_estimate_open_pit_mrem'] = open_pit
    print(f'estimate at 1 %/d: {estimate(1.0):.3f} mrem; at a 1e4 %/d pit: {open_pit:.2f} mrem against the cell-exit '
          f'dose of the run, {attribution["boxB_fgr11_fgr12"].sum():.2f} mrem')
    for pct in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0):
        print(f'  L_pit = {pct:4.1f} %/d: {estimate(pct):6.2f} mrem (full), {estimate(pct, True):6.2f} mrem (chains)')
    full = brentq(lambda p: estimate(p) - PSAR_MHA_MREM, 1.0, 1.0e3)
    chains = brentq(lambda p: estimate(p, True) - PSAR_MHA_MREM, 1.0, 1.0e3)
    results.update({'pit_rate_for_psar_dose_full_pct_per_day': full, 'pit_rate_for_psar_dose_chains_pct_per_day': chains})
    print(f'L_pit reproducing {PSAR_MHA_MREM} mrem: {full:.2f} %/d (full inventory), {chains:.2f} %/d (PSAR chains)')


def section_daughters() -> None:
    """ Upper estimate of the Cs/Rb/Sr daughters of the liberated inventory that PSAR Eq. 13.1-1 would transport """
    atoms, _, lam = f71_activity_and_lambda(2)
    liberated = ('te', 'i', 'xe', 'br', 'kr')
    daughters: dict[str, float] = {}
    for n, count in atoms.items():
        ele = _element_of(n)
        a = _mass_number(n)
        if ele not in liberated or a not in PSAR_CHAINS:
            continue
        z = Z_OF[ele]
        current = f'{ele}-{a}'  # an isomer enters its chain through the ground state of its element
        while current in lam and lam[current] > 0.0 and z + 1 in SYMBOL_OF:
            z += 1
            current = f'{SYMBOL_OF[z]}-{a}'
            if SYMBOL_OF[z] in ('cs', 'rb', 'sr') and current in lam:
                daughters[current] = daughters.get(current, 0.0) + count
    chi_q, v = _chi_q(MHA), _breathing(MHA)
    horizon = MHA['decay_days'] * DAY_S
    total = 0.0
    for d, n0 in sorted(daughters.items(), key=lambda kv: -kv[1]):
        # Born at t = 0 in the enclosure and leaked at the PSAR rate of "all other nuclides".
        released_bq = n0 * chain_release_fraction(lam[d], drv.MHA_OTHER_RATE_S, drv.MHA_CELL_RATE_S, drv.MHA_PIT_RATE_S,
                                           horizon)
        mrem = chi_q * (v * INH_FGR11.get(d, 0.0) + SUB_FGR12.get(d, 0.0) / DAY_S) * released_bq * MREM_PER_SV
        total += mrem
        if mrem > 1e-6:
            print(f'  {d:7s} {n0:.3e} atoms at T0 -> {mrem:.2e} mrem')
    results['daughter_bound_mrem'] = total
    print(f'Cs/Rb/Sr daughters of the liberated inventory: at most {total:.3f} mrem')


def section_hand_check() -> None:
    """ I-131 through enclosure, cell and pit: initial I-131 plus I-131 bred from Te-131m and the A = 131 chain """
    atoms, becq, lam = f71_activity_and_lambda(2)
    eps = (drv.MHA_OTHER_RATE_S, drv.MHA_CELL_RATE_S, drv.MHA_PIT_RATE_S)  # Te and I share the "other" rates
    lam_i, lam_te = lam['i-131'], lam['te-131m']
    # Sb-131, Sn-131 and Te-131 decay to I-131 within an hour; they enter as I-131 at t = 0.
    fast = sum(atoms.get(n, 0.0) for n in ('in-131', 'sn-131', 'sb-131', 'te-131'))
    horizon = MHA['decay_days'] * DAY_S
    # States: Te-131m in A, B, C; I-131 in A, B, C; released I-131 activity [Bq s].
    m = np.zeros((7, 7))
    for l_, off in ((lam_te, 0), (lam_i, 3)):
        m[off, off] = -(l_ + eps[0])
        m[off + 1, off], m[off + 1, off + 1] = eps[0], -(l_ + eps[1])
        m[off + 2, off + 1], m[off + 2, off + 2] = eps[1], -(l_ + eps[2])
    for box in range(3):
        m[3 + box, box] += lam_te  # Te-131m -> I-131 in the same box (beta and IT through Te-131)
    m[6, 5] = eps[2] * lam_i

    def released(n_te: float, n_i: float) -> float:
        return float((expm(m * horizon) @ np.array([n_te, 0, 0, n_i, 0, 0, 0.0]))[6])

    factor = _chi_q(MHA) * _breathing(MHA) * INH_FGR11['i-131'] * MREM_PER_SV
    direct = released(0.0, atoms['i-131'])
    bred = released(atoms['te-131m'], fast)
    frac = direct / becq['i-131']
    results.update({'i131_release_fraction': frac, 'i131_direct_mrem': factor * direct,
                    'i131_bred_mrem': factor * bred, 'te131m_half_life_h': math.log(2) / lam_te / 3600.0})
    print(f'I-131 A0 = {becq["i-131"]:.4e} Bq, released fraction {frac:.4e} -> {factor * direct:.3f} mrem; bred from '
          f'Te-131m ({math.log(2) / lam_te / 3600.0:.2f} h) and the fast A = 131 precursors: {factor * bred:.3f} mrem; '
          f'sum {factor * (direct + bred):.3f} mrem')


def section_iodine_vapour() -> None:
    for name, case, box, tag in (('eol_boxB', EOL, 'boxB', 'icrp72'), ('mha_pit', MHA, 'boxC', 'fgr11_fgr12')):
        inh, sub = DCF_SETS[tag]
        vapour = {**inh, **IODINE_VAPOUR}
        rel = release(case, box)
        base, vap = dose_mrem(rel, case, inh, sub), dose_mrem(rel, case, vapour, sub)
        results.update({f'{name}_{tag}_iodine_vapour_mrem': vap})
        print(f'{name} {tag}: {base:.3f} mrem particulate iodine, {vap:.3f} mrem elemental iodine vapour '
              f'(x {vap / base:.3f})')


SECTIONS = ('totals', 'attribution', 'two-hour', 'positions', 'pit-rate', 'daughters', 'hand-check', 'iodine-vapour')


def main(argv: list[str]) -> int:
    chosen = argv or list(SECTIONS)
    unknown = [s for s in chosen if s not in SECTIONS]
    if unknown:
        print(f'unknown section(s) {unknown}; choose from {SECTIONS}')
        return 2
    attribution = None
    for s in chosen:
        print(f'=== {s}', flush=True)
        if s == 'totals':
            section_totals()
        elif s == 'attribution':
            attribution = section_attribution()
        elif s == 'two-hour':
            section_two_hour()
        elif s == 'positions':
            section_positions(attribution)
        elif s == 'pit-rate':
            section_pit_rate(attribution)
        elif s == 'daughters':
            section_daughters()
        elif s == 'hand-check':
            section_hand_check()
        elif s == 'iodine-vapour':
            section_iodine_vapour()
    out = os.path.join(RUN_DIR, 'paper_numbers_psar.json')
    previous = json.load(open(out)) if os.path.isfile(out) else {}
    with open(out, 'w') as f:
        json.dump({**previous, **results}, f, indent=2, default=float)
    print(f'numbers: {out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
