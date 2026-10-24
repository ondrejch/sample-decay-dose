#!/usr/bin/env python3
"""Xenon verification of the leaky-box paper: Xe-135 feed-step sweep and Xe-136 atom conservation.

Phases:
  run      - ORIGEN runs (SCALE): Xe-135 at 0.5, 1, 2, 4 and 8 steps/day and Xe-136 at 2 steps/day, 12 days each.
             Box A takes its density from the PSAR EOL F71 only to size the volume; the test inventory is
             1 atom/b-cm of the xenon isotope, so the result does not depend on that F71.
  analyse  - (default) reads the leaky_boxes_<prefix>.json snapshots that `run` leaves and writes
             convergence.csv (Table tab:convergence), conservation.csv and the convergence figure.

convergence.csv holds, per step density, the 12-day mean absolute relative difference (Table tab:analytic) and
the largest error normalized to the initial inventory, which stays well conditioned as the analytic solution
decays. conservation.csv holds (Box A + Box B + Box C) / N0 for stable Xe-136; the analysis prints the
trapezoidal excess that endpoint-averaging the Box A leak predicts on the actual output grid.

Usage:  python -m leaky_box_origen.paper.xe_convergence [run|analyse|all] [--run-dir DIR] [--fig-dir DIR]
"""
import argparse
import os

import json5
import numpy as np
import pandas as pd

from leaky_box_origen.LeakyBox import LEAKY_BOX_ORIGEN_DIR, PCTperDAY, run_test

DEFAULT_RUN_DIR = LEAKY_BOX_ORIGEN_DIR / 'run_2026-10-23'
DEFAULT_FIG_DIR = LEAKY_BOX_ORIGEN_DIR / 'journal_paper' / 'figs'
F71 = LEAKY_BOX_ORIGEN_DIR / 'run_2026-10-22' / 'msrr_psar_eol' / 'msrr_psar_eol.f71'
XE135_LAMBDA = 2.106574217602e-5  # [1/s], T1/2 = 9.14 h
DAYS = 12.0
STEPS_PER_DAY = (0.5, 1.0, 2.0, 4.0, 8.0)


def _prefix(isotope: str, spd: float) -> str:
    return f"{isotope.replace('-', '')}_spd{spd:g}".replace('.', 'p')


def run(run_dir: str) -> None:
    """ ORIGEN runs of the sweep and of the conservation test (SCALE, about an hour in total) """
    for spd in STEPS_PER_DAY:
        run_test('xe-135', XE135_LAMBDA, False, False, out_prefix=_prefix('xe-135', spd), f71_path=str(F71),
                 f71_case='1', sample_mass_g=1200e3, decay_days=DAYS, steps_per_day=spd, run_dir=run_dir)
    run_test('xe-136', None, False, False, out_prefix=_prefix('xe-136', 2.0), f71_path=str(F71), f71_case='1',
             sample_mass_g=1200e3, decay_days=DAYS, steps_per_day=2.0, run_dir=run_dir)


def _frames(run_dir: str, prefix: str) -> dict[str, pd.DataFrame]:
    with open(os.path.join(run_dir, f'leaky_boxes_{prefix}.json')) as f:
        return {box: pd.DataFrame(cols) for box, cols in json5.load(f).items()}


def _max_norm_err_pct(frame: pd.DataFrame, n0: float) -> float:
    return float(np.nanmax(np.abs(frame['total'] - frame['analytic'])) / n0 * 100.0)


def _mean_abs_rel_pct(frame: pd.DataFrame) -> float:
    return float(np.nanmean(np.abs(frame['diff [%]'])))


def analyse(run_dir: str, fig_dir: str | None) -> None:
    rows = []
    for spd in STEPS_PER_DAY:
        f = _frames(run_dir, _prefix('xe-135', spd))
        n0 = float(f['box_A']['xe-135'].iloc[0])
        rows.append((spd, _mean_abs_rel_pct(f['box_B']), _mean_abs_rel_pct(f['box_C']),
                     _max_norm_err_pct(f['box_B'], n0), _max_norm_err_pct(f['box_C'], n0)))
    conv = pd.DataFrame(rows, columns=['steps_per_day', 'boxB_mean_abs_rel_pct_all', 'boxC_mean_abs_rel_pct_all',
                                       'boxB_max_abs_norm_pct', 'boxC_max_abs_norm_pct'])
    conv.to_csv(os.path.join(run_dir, 'convergence.csv'), index=False)
    print(conv.to_string(index=False))
    ratios = conv['boxB_max_abs_norm_pct'].to_numpy()
    print('Box B error ratio per step halving:', ' '.join(f'{a / b:.2f}' for a, b in zip(ratios, ratios[1:])))
    for spd in STEPS_PER_DAY:
        x = XE135_LAMBDA * 86400.0 / spd
        print(f'  {24.0 / spd:4g} h: lambda dt = {x:.3f}, sinh(x)/x - 1 = {100.0 * (np.sinh(x) / x - 1.0):.3f} %')

    # Box B and Box C start one step after Box A, so join on time rather than adding the frames by position.
    f = _frames(run_dir, _prefix('xe-136', 2.0))
    n0 = float(f['box_A']['xe-136'].iloc[0])
    joined = (f['box_A'][['time [s]', 'xe-136']]
              .merge(f['box_B'][['time [s]', 'total']], on='time [s]')
              .merge(f['box_C'][['time [s]', 'total']], on='time [s]'))
    conserved = pd.DataFrame({'time [d]': joined['time [s]'] / 86400.0,
                              'fraction_of_initial': (joined['xe-136'] + joined['total_x'] + joined['total_y']) / n0})
    conserved.to_csv(os.path.join(run_dir, 'conservation.csv'), index=False)
    # Endpoint averaging of the convex leak eps N0 exp(-eps t) over a step h supplies (eps h)^2 / 12 of the atoms
    # that leave Box A in that step in excess; summed over the actual grid of Box A output times.
    t = f['box_A']['time [s]'].to_numpy(dtype=float)
    eps = PCTperDAY
    excess = sum((eps * (b - a)) ** 2 / 12.0 * (np.exp(-eps * a) - np.exp(-eps * b)) for a, b in zip(t, t[1:]))
    h = 0.5 * 86400.0
    print(f'Xe-136 conservation: max deviation {conserved["fraction_of_initial"].max() - 1.0:.4e}; '
          f'trapezoid excess on the output grid {excess:.4e} '
          f'(closed form with dt = 0.5 d: {(eps * h) ** 2 / 12.0 * (1.0 - np.exp(-eps * t[-1])):.4e})')

    if fig_dir:
        print(f'figure: {plot(conv, fig_dir)}')


def plot(conv: pd.DataFrame, fig_dir: str) -> str:
    """ Largest normalized error against the feed interval, with a second-order reference line """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    conv = conv.sort_values('steps_per_day')
    step_h = 24.0 / conv['steps_per_day']
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(step_h, conv['boxB_max_abs_norm_pct'], 'o-', color='tab:blue', label='Box B (holdup)')
    ax.loglog(step_h, conv['boxC_max_abs_norm_pct'], 's-', color='tab:orange', label='Box C (collector)')
    ref = conv['boxB_max_abs_norm_pct'].iloc[-1] * (step_h / step_h.iloc[-1]) ** 2
    ax.loglog(step_h, ref, ':', color='0.4', label=r'second order ($\propto \Delta t^2$)')
    ax.set_xlabel('feed interval [h]')
    ax.set_ylabel(r'max $|N_\mathrm{ORIGEN} - N_\mathrm{analytic}| / N_0$ [%]')
    ax.set_title('Xe-135, 12-day horizon: error against the feed interval')
    ax.set_xticks([3, 6, 12, 24, 48])
    ax.set_xticklabels(['3', '6', '12', '24', '48'])
    ax.set_xlim(2.2, 60)
    ax.grid(True, which='both', linestyle=':', alpha=0.5)
    ax.legend()
    fig.tight_layout()
    os.makedirs(fig_dir, exist_ok=True)
    out = os.path.join(fig_dir, 'leaky_boxes_xe135_convergence.png')
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('phase', nargs='?', default='analyse', choices=('run', 'analyse', 'all'))
    parser.add_argument('--run-dir', default=str(DEFAULT_RUN_DIR))
    parser.add_argument('--fig-dir', default=str(DEFAULT_FIG_DIR) if DEFAULT_FIG_DIR.parent.is_dir() else None,
                        help='where the convergence figure goes (default: the paper figs/ when present)')
    args = parser.parse_args()
    if args.phase in ('run', 'all'):
        run(args.run_dir)
    if args.phase in ('analyse', 'all'):
        analyse(args.run_dir, args.fig_dir)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
