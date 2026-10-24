#!/usr/bin/env python3
"""
Extract EPA Federal Guidance Report No. 12 (FGR-12, EPA-402-R-93-081, 1993) air-submersion dose
coefficients from the EPA PDF.

Output (in leaky_box_origen/data/ unless --out-dir is given):
  - dcf_fgr12_submersion.csv   nuclide,dcf_sv_per_bq_m3_day
        Table III.1 effective dose equivalent coefficient h_E [Sv per Bq s m^-3] for submersion in a
        semi-infinite cloud (air density 1.2 kg/m^3), stored as Sv/day per Bq/m^3 (x 86400). h_E excludes
        skin. Rows with h_E = 0.0 (no penetrating emission, e.g. Ca-41, Ru-106) are left out. FGR-12 states
        that Table III.1 replaces Table 2.3 of FGR-11, and its H-3 row is elemental tritium.

The EPA PDF (402-r-93-081.pdf) was produced with Acrobat Distiller, so its text layer holds the printed
values exactly and pdftotext -layout reads the table without OCR. Every row is checked against
h_E = sum_T w_T h_T with the ICRP-26 weights (FGR-12 Table II.1). Nuclide labels carry page-local
footnote markers (H-3 is printed 'H-31,2', Ar-37 'Ar-372') and ICRP-38 isomer letters (Nb-89a/b, Sb-124n).
Table A.1 lists every label without markers, which settles the footnotes, and ISOMER_NAMES maps the
lettered isomers to the ground/m/m2 names of ORIGEN and the other CSVs by their Table A.1 half-lives.

Requirements: pdftotext (poppler). Review the output before using the DCFs for licensing decisions.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from leaky_box_origen.extract_fgr11_dcf import (DEFAULT_OUT_DIR, DEFAULT_PDF_DIR, is_plausible_nuclide,
                                                row_is_consistent, write_csv)

DEFAULT_PDF: Path = DEFAULT_PDF_DIR / "402-r-93-081.pdf"
SUBMERSION_FILE: str = "dcf_fgr12_submersion.csv"
SUBMERSION_HEADER: str = "nuclide,dcf_sv_per_bq_m3_day"
SECONDS_PER_DAY: float = 86400.0  # Table III.1 is Sv per Bq s m^-3, i.e. Sv/s per Bq/m^3

# ICRP-38 labels that FGR-12 uses when neither isomer is designated metastable ('a'/'b') or for a second
# metastable state ('n'), mapped by their Table A.1 half-life to the names used by ORIGEN and the
# ICRP-72/FGR-11 CSVs. Ta-178 follows the ORIGEN decay library, whose ground state has T1/2 = 9.31 min.
ISOMER_NAMES: dict[str, str] = {
    'nb-89b': 'nb-89',      # 122 min, ground state
    'nb-89a': 'nb-89m',     # 66 min
    'in-110b': 'in-110',    # 4.9 h, ground state
    'in-110a': 'in-110m',   # 69.1 min
    'sb-120a': 'sb-120',    # 15.89 min, ground state
    'sb-120b': 'sb-120m',   # 5.76 d
    'sb-124n': 'sb-124m2',  # 20.2 min; Sb-124m is the 93 s state
    'sb-128b': 'sb-128',    # 9.01 h, ground state
    'sb-128a': 'sb-128m',   # 10.4 min
    'eu-150b': 'eu-150',    # 34.2 y, ground state
    'eu-150a': 'eu-150m',   # 12.62 h
    'tb-156n': 'tb-156m2',  # 5.0 h; Tb-156m is the 24.4 h state
    'ta-178a': 'ta-178',    # 9.31 min, ground state
    'ta-178b': 'ta-178m',   # 2.2 h
    're-182b': 're-182',    # 64.0 h, ground state
    're-182a': 're-182m',   # 12.7 h
    'ir-186a': 'ir-186',    # 15.8 h, ground state
    'ir-186b': 'ir-186m',   # 1.75 h
    'ir-190n': 'ir-190m2',  # 3.1 h; Ir-190m is the 1.2 h state
    'np-236a': 'np-236',    # 1.15e5 y, ground state
    'np-236b': 'np-236m',   # 22.5 h
}

_VALUE = r'(?:0\.0|\d\.\d\dE[+-]\d\d)'
# Nuclide label followed by the 9 columns Gonad .. Remainder, Effective, Skin.
_ROW_RE = re.compile(r'^\s*(\S+)' + rf'\s+({_VALUE})' * 9 + r'\s*$')
# Table A.1 row: label, then the half-life ('12.35y', '122m', '115E3y'), which can run into the decay mode
# ('35msA', '164.3µsA').
_A1_RE = re.compile(r'^\s*([A-Z][a-z]?-\d+[a-z]?)\s+[\d.]+(?:E\d+)?(?:us|µs|ms|s|m|h|d|y)(?=\s|[A-Z])')
_LABEL_RE = re.compile(r'([A-Z][a-z]?)-(\d+)([a-z]?)')


@dataclass
class Row:
    label: str            # label as printed, footnote markers included
    nuclide: str | None   # resolved CSV name, or None when the label matches no Table A.1 entry
    values: list[float]   # Gonad, Breast, Lung, R Marrow, B Surface, Thyroid, Remainder, Effective, Skin

    @property
    def effective(self) -> float:
        return self.values[7]


def pdf_text(pdf: Path) -> list[str]:
    """ pdftotext -layout output of the whole report, one entry per line """
    if shutil.which("pdftotext") is None:
        raise RuntimeError("pdftotext (poppler) is required")
    out = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, check=True)
    return out.stdout.decode("utf-8", errors="replace").split("\n")


def table_a1_labels(lines: list[str]) -> set[str]:
    """ Nuclide labels of Table A.1, which carry no footnote markers """
    start = next(i for i, line in enumerate(lines) if 'Table A.1.' in line and 'Summary' in line)
    return {m.group(1) for line in lines[start:] if (m := _A1_RE.match(line))}


def resolve_label(label: str, valid: set[str]) -> str | None:
    """ CSV name of a Table III.1 label: strip footnote markers (digits and commas after the label)
    until it matches a Table A.1 label, then lower-case it and map ICRP-38 isomer letters. """
    base = label
    while base not in valid:
        stripped = re.sub(r'[\d,]$', '', base)
        if stripped == base or not stripped:
            return None
        base = stripped
    name = ISOMER_NAMES.get(base.lower(), base.lower())
    return name if is_plausible_nuclide(name) else None


def parse_table_iii1(lines: list[str], valid: set[str]) -> list[Row]:
    """ Rows of Table III.1 (air submersion), from its first page heading to the first page of Table III.2 """
    start = next(i for i, line in enumerate(lines) if 'Table III.1.' in line and 'Air Submersion' in line)
    end = next(i for i, line in enumerate(lines) if i > start and 'Table III.2.' in line)
    rows: list[Row] = []
    for line in lines[start:end]:
        m = _ROW_RE.match(line)
        if m:
            rows.append(Row(m.group(1), resolve_label(m.group(1), valid), [float(v) for v in m.groups()[1:]]))
    return rows


def submersion_table(rows: list[Row]) -> tuple[dict[str, float], list[Row]]:
    """ {nuclide: h_E in Sv/day per Bq/m^3} for the rows with h_E > 0, and the rows that fail a check """
    data: dict[str, float] = {}
    bad: list[Row] = []
    for r in rows:
        if r.effective == 0.0:
            continue
        if r.nuclide is None or r.nuclide in data or not row_is_consistent(r.values[:8]):
            bad.append(r)
            continue
        data[r.nuclide] = r.effective * SECONDS_PER_DAY
    return data, bad


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Extract FGR-12 Table III.1 air-submersion DCFs.")
    parser.add_argument("--pdf", default=str(DEFAULT_PDF), help=f"FGR-12 PDF (default: {DEFAULT_PDF})")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for the CSV")
    args = parser.parse_args(argv)

    pdf = Path(args.pdf)
    if not pdf.is_file():
        print(f"FGR-12 PDF not found: {pdf}. EPA distributes it as 402-r-93-081.pdf.", file=sys.stderr)
        return 1
    lines = pdf_text(pdf)
    rows = parse_table_iii1(lines, table_a1_labels(lines))
    data, bad = submersion_table(rows)
    for r in bad:
        print(f"  check: {r.label} -> {r.nuclide} values={r.values}")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / SUBMERSION_FILE, data, SUBMERSION_HEADER)
    zeros = sum(1 for r in rows if r.effective == 0.0)
    print(f"Table III.1: {len(rows)} rows, {zeros} with h_E = 0, {len(bad)} unresolved, "
          f"{len(data)} nuclides -> {out_dir / SUBMERSION_FILE}")
    return 0 if not bad else 2


if __name__ == "__main__":
    raise SystemExit(main())
