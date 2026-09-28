#!/bin/env python3
"""
Leaky box using SCALE/Origen
Ondrej Chvala <ochvala@utexas.edu>

Load analysis results from JSONS and re-generates the spreadsheet.

Reads boxA_<prefix>.json5, boxB_<prefix>.json5 and boxC_<prefix>.json5 written by
LeakyBox.run_test (falling back to unprefixed boxA.json5 etc.), adds the analytic
single-isotope solutions, and writes leaky_boxes_rerun_<prefix>.xlsx into the run directory.

The analytic solutions are those of LeakyBox._add_analytic_columns:
N_i^A(t) & = N_i^A e^{ - (\\lambda_i + \\epsilon_i^A) t}
N_i^B(t) & = N_i^A \\epsilon_i^A \\frac{e^{-(\\epsilon_i^A + \\lambda_i) t} - e^{-(\\epsilon_i^B + \\lambda_i) t}}
                                    {\\epsilon_i^B - \\epsilon_i^A}
N_i^C(t) & = N_i^A e^{ - \\lambda_i t}
             \\frac{ -\\epsilon_i^B e^{-\\epsilon_i^A t} + \\epsilon_i^A (e^{-\\epsilon_i^B t} - 1) + \\epsilon_i^B}
                  {\\epsilon_i^B - \\epsilon_i^A}
"""
import argparse
from datetime import date
from pathlib import Path

import pandas as pd
import json5
from leaky_box_origen.LeakyBox import get_dataframe, PCTperDAY, _add_analytic_columns, _out_name

# Leak rates used by LeakyBox._run_simulation [1/s].
box_A_leak_rate: float = PCTperDAY
box_B_leak_rate: float = PCTperDAY * 0.1

# Single-isotope test cases of LeakyBox.run_all_tests: prefix -> (isotope, decay constant [1/s]).
TEST_CASES: dict[str, tuple[str, float | None]] = {
    'xe136': ('xe-136', None),
    'xe135': ('xe-135', 2.106574217602 * 1e-5),
}

# LeakyBox._run_simulation models boxes B and C as 1 cm^3 volumes fed by the leak from
# 1 cm^3 of box A, so the B/C totals equal the B/C atom densities.
DEFAULT_VOLUME: float = 1.0
DEFAULT_N0_DENSITY: float = 1.0  # run_test sets box A to 1 atom/b-cm of the test isotope


def _box_json_names(prefix: str | None) -> list[tuple[str, str, str]]:
    # Candidate (A, B, C) json names: prefixed first, then legacy unprefixed.
    names = []
    if prefix:
        names.append(tuple(_out_name(prefix, f'box{b}', '.json5') for b in 'ABC'))
    names.append(tuple(_out_name(None, f'box{b}', '.json5') for b in 'ABC'))
    return names


def _find_box_jsons(base: Path, prefix: str | None) -> tuple[Path, Path, Path] | None:
    for names in _box_json_names(prefix):
        paths = tuple(base / n for n in names)
        if all(p.is_file() for p in paths):
            return paths
    return None


def _resolve_run_dir(run_dir: str | None, prefix: str | None = None) -> Path:
    if run_dir:
        p = Path(run_dir).expanduser()
        if not p.is_absolute():
            p = (Path.cwd() / p).resolve()
        else:
            p = p.resolve()
        return p

    cwd = Path.cwd()
    if _find_box_jsons(cwd, prefix):
        return cwd

    base = Path(__file__).resolve().parent
    today = base / f'run_{date.today().isoformat()}'
    if _find_box_jsons(today, prefix):
        return today

    # Newest run first, but only a run that actually holds this prefix's box JSONs.
    # Returning the newest run blindly can silently rerun the wrong isotope case.
    for run in sorted(base.glob('run_*'), reverse=True):
        if _find_box_jsons(run, prefix):
            return run

    raise FileNotFoundError(f"No boxA/boxB/boxC json5 files for prefix '{prefix}' under {base}")


def main(run_dir: str | None = None, prefix: str = 'xe135', volume: float = DEFAULT_VOLUME,
         n0_density: float = DEFAULT_N0_DENSITY, write_excel: bool = True
         ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if prefix not in TEST_CASES:
        raise ValueError(f"Unknown prefix '{prefix}', expected one of {sorted(TEST_CASES)}")
    isotope, lambda_decay = TEST_CASES[prefix]
    base = _resolve_run_dir(run_dir, prefix)
    paths = _find_box_jsons(base, prefix)
    if paths is None:
        raise FileNotFoundError(f"No boxA/boxB/boxC json5 files for prefix '{prefix}' in {base}")

    box_adens = []
    for path in paths:
        with open(path, 'r') as f:
            box_adens.append(json5.load(f))
    pd_A: pd.DataFrame = get_dataframe(box_adens[0])
    pd_B: pd.DataFrame = get_dataframe(box_adens[1])
    pd_C: pd.DataFrame = get_dataframe(box_adens[2])

    pd_B['total'] = pd_B[isotope] * volume
    pd_C['total'] = pd_C[isotope] * volume
    _add_analytic_columns(pd_A, pd_B, pd_C, isotope, n0_density, n0_density * volume,
                          box_A_leak_rate, box_B_leak_rate, lambda_decay)

    if write_excel:
        writer = pd.ExcelWriter(base / _out_name(prefix, 'leaky_boxes_rerun', '.xlsx'))
        pd_A.to_excel(writer, sheet_name='box A')
        pd_B.to_excel(writer, sheet_name='box B')
        pd_C.to_excel(writer, sheet_name='box C')
        writer.close()
    return pd_A, pd_B, pd_C


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Regenerate leaky-box spreadsheet from saved JSON files.")
    parser.add_argument("--run-dir", default=None,
                        help="Run directory containing boxA_<prefix>.json5/boxB_<prefix>.json5/boxC_<prefix>.json5")
    parser.add_argument("--prefix", default='xe135', choices=sorted(TEST_CASES),
                        help="Test-case output prefix (default: xe135)")
    parser.add_argument("--volume", type=float, default=DEFAULT_VOLUME,
                        help="Box B/C volume [cm^3] used to convert densities to totals (default: 1.0)")
    args = parser.parse_args()
    main(run_dir=args.run_dir, prefix=args.prefix, volume=args.volume)
