#!/bin/env python3
"""Run per-step ORIGEN decays for case-1 records in an input F71 file.

For each time step (position) belonging to case 1 in the source F71:
- read atom densities from that position
- create a dedicated run directory
- run ORIGEN decay for a fixed material volume and decay time
- read final activities from the generated ORIGEN F71

Outputs:
- CSV with total activity and F-19 activity vs source time
- PNG plot with total and F-19 activity curves
"""

from __future__ import annotations

import argparse
import os
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from datetime import date
from pathlib import Path

import pandas as pd

from sample_decay_dose.SampleDose import NOW
from sample_decay_dose.utils import (
    atom_dens_for_origen,
    get_burned_nuclide_atom_dens,
    get_burned_nuclide_data,
    get_f71_positions_index,
    run_scale,
)


DEFAULT_F71_PATH: str = os.getenv(
    "LEAKYBOX_F71_PATH",
    os.path.expanduser("~/0.03/20-burn-MHA/mha-4.5-a4/msrr.f71"),
)
DEFAULT_CASE: int = 1
MODULE_DIR: Path = Path(__file__).resolve().parent


@dataclass
class StepActivity:
    source_position: int
    source_time_s: float
    source_time_d: float
    run_dir: str
    total_activity_bq: float
    f19_activity_bq: float


@contextmanager
def _in_directory(path: Path):
    prev = Path.cwd()
    path.mkdir(parents=True, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def _resolve_output_root(out_root: str | None) -> Path:
    if out_root:
        return Path(out_root).expanduser().resolve()
    return (MODULE_DIR / f"run_{date.today().isoformat()}" / "case1_f71_decay_2d").resolve()


def _case_positions(f71_path: str, case: int) -> list[tuple[int, float]]:
    f71_idx = get_f71_positions_index(f71_path)
    rows: list[tuple[int, float]] = []
    for pos, rec in f71_idx.items():
        if int(rec["case"]) != int(case):
            continue
        rows.append((int(pos), float(rec["time"])))
    return sorted(rows, key=lambda item: item[1])


def _origen_decay_deck(
    *,
    atom_file: str,
    out_f71: str,
    volume_cm3: float,
    decay_days: float,
    decay_steps: int,
) -> str:
    interp_steps = max(int(decay_steps) - 3, 1)
    return f"""=shell
cp -r ${{INPDIR}}/{atom_file} .
end

=origen
' {NOW}
options{{
    digits=6
}}
bounds {{
    neutron="scale.rev13.xn200g47v7.1"
    gamma="scale.rev13.xn200g47v7.1"
    beta=[100L 1.0e7 1.0e-3]
}}
case {{
    gamma=yes
    neutron=yes
    beta=yes
    lib {{
        file="end7dec"
    }}
    mat {{
        iso [
<{atom_file}
        ]
        units=ATOMS-PER-BARN-CM
        volume={volume_cm3}
    }}
    time {{
        units=DAYS
        t=[{interp_steps}I 0.01 {decay_days}]
        start=0
    }}
    save {{
        file="{out_f71}"
    }}
}}
end
"""


def _run_single_step(
    *,
    f71_path: str,
    source_position: int,
    source_time_s: float,
    out_root: Path,
    case: int,
    volume_cm3: float,
    decay_days: float,
    decay_steps: int,
    skip_scale: bool,
) -> StepActivity:
    run_dir_name = f"case{case}_pos{source_position:04d}_t{source_time_s:012.1f}s"
    run_dir = out_root / run_dir_name
    atom_file = "sample_atom_dens.inp"
    deck_file = "origen.inp"
    out_f71 = "origen.f71"

    atom_dens = get_burned_nuclide_atom_dens(f71_path, source_position)
    if not atom_dens:
        raise RuntimeError(f"No atom-density data at position {source_position} in {f71_path}")

    with _in_directory(run_dir):
        with open(atom_file, "w", encoding="utf-8") as fh:
            fh.write(atom_dens_for_origen(atom_dens))
        with open(deck_file, "w", encoding="utf-8") as fh:
            fh.write(
                _origen_decay_deck(
                    atom_file=atom_file,
                    out_f71=out_f71,
                    volume_cm3=volume_cm3,
                    decay_days=decay_days,
                    decay_steps=decay_steps,
                )
            )

        if not skip_scale:
            ok = run_scale(deck_file)
            if not ok:
                raise RuntimeError(f"SCALE run failed for {run_dir_name}")

        if not Path(out_f71).is_file():
            raise RuntimeError(f"Missing ORIGEN output {run_dir / out_f71}")

        activity_bq = get_burned_nuclide_data(out_f71, -1, f71units="becq")

    total_bq = float(sum(activity_bq.values()))
    f19_bq = float(activity_bq.get("f-19", 0.0))
    return StepActivity(
        source_position=source_position,
        source_time_s=source_time_s,
        source_time_d=source_time_s / float(24 * 60 * 60),
        run_dir=str(run_dir),
        total_activity_bq=total_bq,
        f19_activity_bq=f19_bq,
    )


def _plot_activity(df: pd.DataFrame, png_path: Path, *, logy: bool = True) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    ax.plot(df["source_time_d"], df["total_activity_bq"], label="Total activity", color="tab:blue", marker="o")
    ax.plot(df["source_time_d"], df["f19_activity_bq"], label="F-19 activity", color="tab:orange", marker="s")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("Source case time [days]")
    ax.set_ylabel("Activity [Bq]")
    ax.set_title("2-day decay of 6.5 L sample from each case-1 timestep")
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    plt.close(fig)


def run_case1_decay_series(
    *,
    f71_path: str,
    case: int = 1,
    volume_liters: float = 6.5,
    decay_days: float = 2.0,
    decay_steps: int = 30,
    out_root: str | None = None,
    skip_scale: bool = False,
    logy: bool = True,
) -> tuple[pd.DataFrame, Path, Path]:
    out_dir = _resolve_output_root(out_root)
    out_dir.mkdir(parents=True, exist_ok=True)

    positions = _case_positions(f71_path, case)
    if not positions:
        raise ValueError(f"No positions found for case {case} in {f71_path}")

    volume_cm3 = float(volume_liters) * 1000.0
    rows: list[StepActivity] = []
    for pos, time_s in positions:
        print(f"Running case {case}, position {pos}, t={time_s:.3f} s")
        row = _run_single_step(
            f71_path=f71_path,
            source_position=pos,
            source_time_s=time_s,
            out_root=out_dir,
            case=case,
            volume_cm3=volume_cm3,
            decay_days=decay_days,
            decay_steps=decay_steps,
            skip_scale=skip_scale,
        )
        rows.append(row)

    df = pd.DataFrame([asdict(row) for row in rows]).sort_values("source_time_s")
    csv_path = out_dir / "case_decay_activity.csv"
    png_path = out_dir / "case_decay_activity.png"
    df.to_csv(csv_path, index=False)
    _plot_activity(df, png_path, logy=logy)
    return df, csv_path, png_path


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "For each timestep in F71 case 1, run a 2-day ORIGEN decay of a 6.5-liter sample "
            "and plot total/F-19 activity."
        )
    )
    parser.add_argument("--f71", default=DEFAULT_F71_PATH, help="Input F71 path.")
    parser.add_argument("--case", type=int, default=DEFAULT_CASE, help="F71 case to process (default: 1).")
    parser.add_argument("--volume-liters", type=float, default=6.5, help="Sample volume in liters.")
    parser.add_argument("--decay-days", type=float, default=2.0, help="ORIGEN decay time in days.")
    parser.add_argument("--decay-steps", type=int, default=30, help="ORIGEN time steps.")
    parser.add_argument("--out-root", default=None, help="Output root directory.")
    parser.add_argument(
        "--skip-scale",
        action="store_true",
        help="Skip running SCALE and only read existing per-step origen.f71 files.",
    )
    parser.add_argument(
        "--linear-y",
        action="store_true",
        help="Use linear y-axis for activity plot (default is log scale).",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    df, csv_path, png_path = run_case1_decay_series(
        f71_path=args.f71,
        case=args.case,
        volume_liters=args.volume_liters,
        decay_days=args.decay_days,
        decay_steps=args.decay_steps,
        out_root=args.out_root,
        skip_scale=args.skip_scale,
        logy=not args.linear_y,
    )
    print(f"Wrote {len(df)} rows")
    print(f"CSV: {csv_path}")
    print(f"Plot: {png_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
