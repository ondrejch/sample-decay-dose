#!/usr/bin/env python3
"""
Extract EPA Federal Guidance Report No. 11 (FGR-11, EPA-520/1-88-020, 1988) dose conversion
factors from scanned copies of the report.

Outputs (in leaky_box_origen/data/ unless --out-dir is given):
  - dcf_fgr11_inhalation_worker.csv   nuclide,dcf_sv_bq
        Table 2.1 effective dose equivalent per unit intake h_E,50 [Sv/Bq]. FGR-11 lists one row per
        lung clearance class (D/W/Y), vapour (V) or chemical form. The CSV keeps the LARGEST
        h_E,50 over the listed rows of each nuclide, which is the conservative choice when the
        chemical form of the release is unknown.
  - dcf_fgr11_submersion_worker.csv   nuclide,dcf_sv_per_bq_m3_day
        Table 2.3 effective dose equivalent rate per unit air concentration h_E,ext. FGR-11 tabulates
        it in Sv/hr per Bq/m^3 (see the table heading), and the CSV stores Sv/day per Bq/m^3 (x 24).
        h_E,ext excludes skin and lens; FGR-11 lists those organ values only where they are limiting.
The ICRP-72 DCF CSVs in data/ are produced by extract_icrp72_dcf.py.

Why not the PDF text layers alone: both known PDF editions (EPA "EPA 1988_FGR11_0.pdf" and the OSTI
copy 6294233.pdf) carry OCR text layers in which the mantissas are mostly right but the superscript
exponents ("1.73 10^-11" printed with a raised "-11") are garbled ("10"", "IO-''", "lo-''"). The
previous extractor decoded those suffixes heuristically and produced errors of up to six orders of
magnitude. This extractor therefore:
  1. renders each table page at 300 dpi (pdftoppm) and finds every "m.mm 10^-e" value cell
     geometrically: glyph components on a text line, the "10", and the raised exponent digits;
  2. reads mantissas with tesseract, and reads exponent digits with a nearest-neighbour classifier
     trained on the mantissa digits of the same scan (plus a tesseract reading as a second vote);
  3. adds the mantissas of both PDF text layers (pdftotext -layout) as further votes;
  4. repeats 1-2 on the second scan when it is available, and aligns the two scans row by row;
  5. checks every Table 2.1 and Table 2.3 row against the definition of the effective value,
     h_E = sum_T w_T h_T with the ICRP-26 weights (gonad 0.25, breast 0.15, lung 0.12, red marrow 0.12,
     bone surface 0.03, thyroid 0.03, remainder 0.30). When the most-voted reading of a row fails the
     check, the lowest-cost combination of alternative readings that satisfies it is used;
  6. applies MANUAL_CORRECTIONS, which hold values read by eye from the page images for rows the
     automatic pass could not settle. Each entry names the page it was read from. A Table 2.1 row
     with an undetected cell counts as unsettled whenever the solver would have to alter another
     reading, because all eight values are printed in every Table 2.1 row.
Rows that still fail the check are listed on stdout and in the audit CSV (--audit).

Requirements: pdftoppm, pdftotext, pdfinfo (poppler), tesseract (eng), numpy, scipy, Pillow.
Review the output before using the DCFs for licensing decisions.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import itertools
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

REPO_DIR: Path = Path(__file__).resolve().parent.parent
DEFAULT_PDF_DIR: Path = REPO_DIR / "PDF"
DEFAULT_OUT_DIR: Path = Path(__file__).resolve().parent / "data"

# ICRP-26 tissue weighting factors in FGR-11 column order (Gonad, Breast, Lung, R Marrow, B Surface,
# Thyroid, Remainder). The Effective column is the weighted sum (FGR-11 p. 6).
ORGAN_WEIGHTS: tuple[float, ...] = (0.25, 0.15, 0.12, 0.12, 0.03, 0.03, 0.30)
N_VALUE_COLUMNS: int = 8  # 7 organs + effective
# Organ values and h_E are printed to 3 significant figures, so the weighted sum reproduces h_E to
# about 1%. Rows outside this tolerance hold a misread value.
CONSISTENCY_TOL: float = 0.012

# Plausibility window for FGR-11 coefficients (Sv/Bq and Sv/hr per Bq/m^3). Values outside it are
# extraction errors: the smallest printed inhalation factor is ~1e-14 Sv/Bq, and no factor reaches
# 1e-2.
MIN_DCF: float = 1e-20
MAX_DCF: float = 1e-2
EXPONENT_RANGE: tuple[int, int] = (2, 20)  # printed exponents run from 10^-2 (bone surface) to 10^-20
# Exponent readings of 2 are almost always a dropped leading '1' (10^-12); only the free-exponent
# pass of build_rows may use 10^-2.
VOTE_EXPONENT_RANGE: tuple[int, int] = (3, 20)
SUBMERSION_HR_TO_DAY: float = 24.0  # Table 2.3 is Sv/hr per Bq/m^3; the CSV stores Sv/day per Bq/m^3

INHALATION_HEADER: str = "nuclide,dcf_sv_bq"
SUBMERSION_HEADER: str = "nuclide,dcf_sv_per_bq_m3_day"


@dataclass(frozen=True)
class Edition:
    """ Page layout of one PDF edition of FGR-11 """
    name: str
    n_pages: int
    inhalation_pages: tuple[int, ...]  # Table 2.1
    submersion_pages: tuple[int, ...]  # Table 2.3
    default_file: str


EDITIONS: dict[str, Edition] = {
    # OSTI scan: clean 300 dpi page images. Used as the primary reading.
    'osti': Edition('osti', 234, tuple(range(131, 163)), (191,), '6294233.pdf'),
    # EPA scan: blurrier, used as the second reading.
    'epa': Edition('epa', 224, tuple(range(127, 159)), (185,), 'EPA 1988_FGR11_0.pdf'),
}

# Chemical elements Z = 1..103: (symbol, name[, alternative spellings]).
ELEMENTS: tuple[tuple[str, ...], ...] = (
    ('H', 'Hydrogen'), ('He', 'Helium'), ('Li', 'Lithium'), ('Be', 'Beryllium'), ('B', 'Boron'),
    ('C', 'Carbon'), ('N', 'Nitrogen'), ('O', 'Oxygen'), ('F', 'Fluorine'), ('Ne', 'Neon'),
    ('Na', 'Sodium'), ('Mg', 'Magnesium'), ('Al', 'Aluminum', 'Aluminium'), ('Si', 'Silicon'),
    ('P', 'Phosphorus'), ('S', 'Sulfur', 'Sulphur'), ('Cl', 'Chlorine'), ('Ar', 'Argon'),
    ('K', 'Potassium'), ('Ca', 'Calcium'), ('Sc', 'Scandium'), ('Ti', 'Titanium'), ('V', 'Vanadium'),
    ('Cr', 'Chromium'), ('Mn', 'Manganese'), ('Fe', 'Iron'), ('Co', 'Cobalt'), ('Ni', 'Nickel'),
    ('Cu', 'Copper'), ('Zn', 'Zinc'), ('Ga', 'Gallium'), ('Ge', 'Germanium'), ('As', 'Arsenic'),
    ('Se', 'Selenium'), ('Br', 'Bromine'), ('Kr', 'Krypton'), ('Rb', 'Rubidium'), ('Sr', 'Strontium'),
    ('Y', 'Yttrium'), ('Zr', 'Zirconium'), ('Nb', 'Niobium'), ('Mo', 'Molybdenum'), ('Tc', 'Technetium'),
    ('Ru', 'Ruthenium'), ('Rh', 'Rhodium'), ('Pd', 'Palladium'), ('Ag', 'Silver'), ('Cd', 'Cadmium'),
    ('In', 'Indium'), ('Sn', 'Tin'), ('Sb', 'Antimony'), ('Te', 'Tellurium'), ('I', 'Iodine'),
    ('Xe', 'Xenon'), ('Cs', 'Cesium', 'Caesium'), ('Ba', 'Barium'), ('La', 'Lanthanum'),
    ('Ce', 'Cerium'), ('Pr', 'Praseodymium'), ('Nd', 'Neodymium'), ('Pm', 'Promethium'),
    ('Sm', 'Samarium'), ('Eu', 'Europium'), ('Gd', 'Gadolinium'), ('Tb', 'Terbium'),
    ('Dy', 'Dysprosium'), ('Ho', 'Holmium'), ('Er', 'Erbium'), ('Tm', 'Thulium'), ('Yb', 'Ytterbium'),
    ('Lu', 'Lutetium'), ('Hf', 'Hafnium'), ('Ta', 'Tantalum'), ('W', 'Tungsten', 'Wolfram'),
    ('Re', 'Rhenium'), ('Os', 'Osmium'), ('Ir', 'Iridium'), ('Pt', 'Platinum'), ('Au', 'Gold'),
    ('Hg', 'Mercury'), ('Tl', 'Thallium'), ('Pb', 'Lead'), ('Bi', 'Bismuth'), ('Po', 'Polonium'),
    ('At', 'Astatine'), ('Rn', 'Radon'), ('Fr', 'Francium'), ('Ra', 'Radium'), ('Ac', 'Actinium'),
    ('Th', 'Thorium'), ('Pa', 'Protactinium'), ('U', 'Uranium'), ('Np', 'Neptunium'),
    ('Pu', 'Plutonium'), ('Am', 'Americium'), ('Cm', 'Curium'), ('Bk', 'Berkelium'),
    ('Cf', 'Californium'), ('Es', 'Einsteinium'), ('Fm', 'Fermium'), ('Md', 'Mendelevium'),
    ('No', 'Nobelium'), ('Lr', 'Lawrencium'),
)
ELEMENT_Z: dict[str, int] = {e[0].lower(): z for z, e in enumerate(ELEMENTS, start=1)}
_ELEMENT_BY_NAME: dict[str, str] = {name.lower(): e[0] for e in ELEMENTS for name in e[1:]}

# Values read by eye from the OSTI page images (6294233.pdf) for rows whose automatic reading
# failed the h_E consistency check. Key: (table, nuclide, row index within the nuclide) with table
# '2.1' or '2.3'; value: the 8 printed values (Gonad ... Effective, None where blank) and the page.
MANUAL_CORRECTIONS: dict[tuple[str, str, int], tuple[tuple[float | None, ...], int]] = {
    # Exponent or effective cells that touch the page frame or a neighbouring glyph in the scan.
    ('2.1', 'te-123', 1): ((3.31e-12, 3.28e-12, 5.19e-10, 2.57e-9, 3.12e-8, 2.22e-12, 2.03e-11, 1.31e-9), 144),
    ('2.1', 'te-123m', 1): ((1.88e-10, 2.04e-10, 1.27e-8, 2.41e-9, 2.40e-8, 1.46e-10, 8.06e-10, 2.86e-9), 144),
    ('2.1', 'te-125m', 1): ((7.93e-11, 7.08e-11, 1.04e-8, 1.15e-9, 1.18e-8, 3.87e-11, 6.75e-10, 1.97e-9), 144),
    ('2.1', 'sm-146', 0): ((0.0, 0.0, 8.40e-6, 3.03e-5, 3.79e-4, 0.0, 2.08e-5, 2.23e-5), 148),
    ('2.1', 'eu-154', 0): ((1.17e-8, 1.55e-8, 7.92e-8, 1.06e-7, 5.23e-7, 7.14e-9, 1.13e-7, 7.73e-8), 148),
    ('2.1', 'dy-166', 0): ((2.86e-11, 8.09e-12, 9.10e-9, 4.37e-10, 1.44e-9, 2.68e-12, 2.75e-9, 2.02e-9), 149),
    ('2.1', 'ac-227', 0): ((3.96e-4, 6.66e-8, 1.23e-7, 2.57e-3, 3.21e-2, 3.59e-8, 1.47e-3, 1.81e-3), 158),
    ('2.1', 'th-229', 0): ((2.76e-6, 2.76e-6, 7.95e-5, 1.15e-3, 1.43e-2, 2.76e-6, 7.05e-6, 5.80e-4), 159),
    ('2.1', 'th-232', 0): ((7.62e-7, 7.72e-7, 1.44e-5, 8.93e-4, 1.11e-2, 7.44e-7, 1.87e-6, 4.43e-4), 159),
    ('2.1', 'u-234', 2): ((2.65e-9, 2.68e-9, 2.98e-4, 7.22e-8, 1.13e-6, 2.65e-9, 1.06e-7, 3.58e-5), 159),
    # Rows with an undetected (bold) organ cell, where the solver had to alter other cells.
    ('2.1', 'mo-99', 1): ((9.51e-11, 2.75e-11, 4.29e-9, 5.24e-11, 4.13e-11, 1.52e-11, 1.74e-9, 1.07e-9), 138),
    ('2.1', 'te-121m', 0): ((1.18e-9, 1.23e-9, 1.41e-9, 9.42e-9, 6.94e-8, 1.12e-9, 1.38e-9, 4.31e-9), 144),
    ('2.1', 'te-123', 0): ((7.21e-12, 6.92e-12, 1.61e-11, 5.86e-9, 7.13e-8, 5.03e-12, 1.15e-11, 2.85e-9), 144),
    ('2.1', 'te-123m', 0): ((2.77e-10, 2.80e-10, 6.05e-10, 5.79e-9, 6.09e-8, 2.40e-10, 4.75e-10, 2.86e-9), 144),
    ('2.1', 'te-125m', 0): ((1.24e-10, 1.07e-10, 4.66e-10, 3.01e-9, 3.21e-8, 9.93e-11, 3.14e-10, 1.52e-9), 144),
    ('2.1', 'pu-239', 1): ((1.20e-5, 3.99e-10, 3.23e-4, 6.57e-5, 8.21e-4, 3.75e-10, 3.02e-5, 8.33e-5), 160),
    ('2.1', 'pu-240', 1): ((1.20e-5, 4.33e-10, 3.23e-4, 6.57e-5, 8.21e-4, 3.76e-10, 3.02e-5, 8.33e-5), 160),
}


# --------------------------------------------------------------------------------------------------
# Pure helpers (unit tested without PDFs)
# --------------------------------------------------------------------------------------------------

def element_symbol(text: str) -> str | None:
    """ Element symbol for an element name ('Hydrogen', 'Sulphur') or symbol ('Cs', 'cs'), else None """
    t = re.sub(r'[^A-Za-z]', '', text or '')
    if not t:
        return None
    if t.lower() in _ELEMENT_BY_NAME:
        return _ELEMENT_BY_NAME[t.lower()]
    if t.lower() in ELEMENT_Z:
        return ELEMENTS[ELEMENT_Z[t.lower()] - 1][0]
    return None


def _edit_distance(a: str, b: str) -> int:
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        cur = [i]
        for j, cb in enumerate(b, start=1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


def element_from_heading(word: str) -> str | None:
    """ Element symbol for an OCR'd element heading ('Iodine', 'lodine', 'Sulphur', 'Hydrogen*').

    Headings are whole words starting with a capital letter (OCR may turn 'I' into 'l'). One
    wrong letter is tolerated for names of five to seven letters, two for longer names.
    """
    w = re.sub(r'[^A-Za-z]', '', word or '')
    if len(w) < 3 or not (w[0].isupper() or w[0] == 'l'):
        return None
    w = w.lower()
    if w in _ELEMENT_BY_NAME:
        return _ELEMENT_BY_NAME[w]
    if len(w) >= 5:
        tol = 2 if len(w) >= 8 else 1
        close = [(_edit_distance(w, name), sym) for name, sym in _ELEMENT_BY_NAME.items()
                 if abs(len(name) - len(w)) <= tol]
        close = [c for c in close if c[0] <= tol]
        if len(close) == 1:
            return close[0][1]
    return None


def is_plausible_nuclide(nuclide: str) -> bool:
    """ True for 'el-A', 'el-Am' or 'el-Am2' with a real element and a mass number in the band of
    known nuclides: A >= Z for Z < 10, A >= 1.8 Z above (e.g. iodine starts near A = 108), and
    A <= 3 Z + 30. This rejects OCR names such as 'i-78m', 'i-94' or 'to-104'.
    """
    m = re.fullmatch(r'([a-z]{1,2})-(\d{1,3})(m\d?)?', nuclide or '')
    if not m or m.group(1) not in ELEMENT_Z:
        return False
    z, a = ELEMENT_Z[m.group(1)], int(m.group(2))
    a_min = z if z < 10 else 1.8 * z
    return a_min <= a <= 3 * z + 30 and (z > 1 or a <= 3)


_DIGITISH = r"0-9IlOoSBZi|!\]\}"
_LABEL_RE = re.compile(rf"^\s*(\S{{1,3}}?)\s*-\s*([{_DIGITISH}](?:[{_DIGITISH}\s]*[{_DIGITISH}])?)\s*(m\s*\d?|rn)?"
                       rf"\s*[^0-9A-Za-z]*$")
_OCR_DIGIT_FIXES = str.maketrans({'I': '1', 'l': '1', 'i': '1', '|': '1', '!': '1', ']': '1', '}': '1',
                                  'O': '0', 'o': '0', 'S': '5', 'B': '8', 'Z': '2'})
_HALF_LIFE_RE = re.compile(r"^\s*(\d+(?:[.,]\d+)?)\s*([mhdy])\b")

# FGR-11 separates some isomer pairs by half-life on a second label line instead of an 'm' suffix.
# (nuclide printed, half-life label) -> name used in the CSVs (ground state / metastable state).
HALF_LIFE_ISOMERS: dict[tuple[str, str], str] = {
    ('nb-89', '122 m'): 'nb-89',     # Nb-89 ground state, T1/2 = 2.03 h
    ('nb-89', '66 m'): 'nb-89m',     # Nb-89m, T1/2 = 66 min
    ('in-110', '4.9 h'): 'in-110',   # In-110 ground state, T1/2 = 4.9 h
    ('in-110', '69.1 m'): 'in-110m',  # In-110m, T1/2 = 69.1 min
    ('sb-120', '15.89 m'): 'sb-120',  # Sb-120 ground state, T1/2 = 15.89 min
    ('sb-120', '5.76 d'): 'sb-120m',  # Sb-120m, T1/2 = 5.76 d
    ('sb-128', '9.01 h'): 'sb-128',   # Sb-128 ground state, T1/2 = 9.01 h
    ('sb-128', '10.4 m'): 'sb-128m',  # Sb-128m, T1/2 = 10.4 min
    ('eu-150', '34.2 y'): 'eu-150',   # Eu-150 ground state (FGR-11 half-life 34.2 y)
    ('eu-150', '12.62 h'): 'eu-150m',  # Eu-150m, T1/2 = 12.6 h
    ('tb-156m', '24.4 h'): 'tb-156m',  # Tb-156m1, T1/2 = 24.4 h
    ('tb-156m', '5.0 h'): 'tb-156m2',  # Tb-156m2, T1/2 = 5.0 h
}


def parse_nuclide_label(label: str, element: str | None = None) -> str | None:
    """ Parse an OCR'd FGR-11 nuclide label ('Be-7', 'Kr-8 1', 'Xe- 129m', 'Sm-14im') into 'be-7'.

    element (symbol from the table's element heading) replaces the label's symbol when given, because
    headings are long words that OCR reliably while one- or two-letter symbols are often misread
    ('C1-38', 'T1-204', '[r-185'). Returns None when no mass number is found or the result is not a
    plausible nuclide.
    """
    if not label:
        return None
    s = label.strip().replace('\u2014', '-').replace('\u2013', '-').replace('_', '-')
    m = _LABEL_RE.match(s)
    if not m:
        return None
    sym = element if element else element_symbol(m.group(1))
    if not sym:
        return None
    mass = re.sub(r'\s+', '', m.group(2)).translate(_OCR_DIGIT_FIXES)
    if not mass.isdigit():
        return None
    meta = ''
    if m.group(3):
        meta = 'm' + re.sub(r'\D', '', m.group(3).replace('rn', 'm'))
    nuc = f'{sym.lower()}-{int(mass)}{meta}'
    return nuc if is_plausible_nuclide(nuc) else None


def half_life_label(label: str) -> str | None:
    """ '122 m', '5.76 d', '49h' -> normalised '122 m' / '5.76 d' / '4.9 h' form, else None """
    m = _HALF_LIFE_RE.match(label or '')
    if not m:
        return None
    num = m.group(1).replace(',', '.')
    if m.group(2) == 'h' and num == '49':
        num = '4.9'  # OCR drops the decimal point of '4.9 h'
    return f'{num} {m.group(2)}'


def text_layer_mantissas(line: str) -> list[str]:
    """ Mantissas 'd.dd' found in one pdftotext line, in order (commas read as decimal points) """
    return [f'{a}.{b}' for a, b in re.findall(r'(?<![\d.,])(\d)[.,](\d\d)(?![\d])', line)]


def effective_from_organs(organs: list[float | None]) -> float:
    """ ICRP-26 weighted sum of the 7 organ values; blank organs (None) count as zero """
    return float(sum(w * (v or 0.0) for w, v in zip(ORGAN_WEIGHTS, organs)))


def row_is_consistent(values: list[float | None], tol: float = CONSISTENCY_TOL) -> bool:
    """ True when values[7] (effective) equals the weighted organ sum of values[0:7] within tol """
    if len(values) != N_VALUE_COLUMNS or values[7] is None or values[7] <= 0.0:
        return False
    return abs(effective_from_organs(list(values[:7])) - values[7]) <= tol * values[7]


def solve_row(options: list[list[tuple[float, float]]], tol: float = CONSISTENCY_TOL,
              max_changes: int = 3) -> tuple[list[float | None], float, bool]:
    """ Pick one value per column so that the row satisfies the h_E consistency check.

    options[i] lists (value, cost) candidates for column i, best first (cost 0 = unanimous reading);
    an empty list means a blank cell. The first-choice combination is tried first, then combinations
    that change up to max_changes columns, cheapest first. Returns (values, cost, consistent). When
    no combination is consistent, the first choices are returned with consistent=False.
    """
    first = [opt[0][0] if opt else None for opt in options]
    base_cost = sum(opt[0][1] for opt in options if opt)
    if row_is_consistent(first, tol):
        return first, base_cost, True
    best: tuple[float, list[float | None]] | None = None
    cols = [i for i, opt in enumerate(options) if len(opt) > 1]
    for n_change in range(1, max_changes + 1):
        for subset in itertools.combinations(cols, n_change):
            alt_lists = [options[i][1:] for i in subset]
            for alts in itertools.product(*alt_lists):
                vals = list(first)
                cost = base_cost
                for i, (v, c) in zip(subset, alts):
                    vals[i] = v
                    cost += c - options[i][0][1]
                if best is not None and cost >= best[0]:
                    continue
                if row_is_consistent(vals, tol):
                    best = (cost, vals)
        if best is not None:
            return best[1], best[0], True
    return first, base_cost, False


def cell_options(mantissa_votes: dict[str, float], exponent_votes: dict[int, float],
                 max_mantissas: int = 3, max_exponents: int = 3,
                 widen_exponents: bool = True) -> list[tuple[float, float]]:
    """ Candidate values for one table cell from weighted mantissa and exponent votes.

    Cost of a candidate = (total vote weight - weight of its mantissa) + (total - weight of its
    exponent), so the unanimous reading costs 0. With widen_exponents, the neighbours +/-1 of the
    best exponent are added at a high cost; they fix single-digit exponent misreads that every
    reader shares. Candidates are sorted by cost.
    """
    if not mantissa_votes:
        return []
    m_tot = sum(mantissa_votes.values())
    mants = sorted(mantissa_votes.items(), key=lambda kv: -kv[1])[:max_mantissas]
    e_votes = dict(exponent_votes)
    e_tot = sum(e_votes.values()) or 1.0
    exps = sorted(e_votes.items(), key=lambda kv: -kv[1])[:max_exponents]
    if widen_exponents:
        e0 = exps[0][0] if exps else None
        extra = [e0 - 1, e0 + 1] if e0 is not None else list(range(EXPONENT_RANGE[0], EXPONENT_RANGE[1] + 1))
        for e in extra:
            if EXPONENT_RANGE[0] <= e <= EXPONENT_RANGE[1] and e not in dict(exps):
                exps.append((e, -2.0 * e_tot))  # penalised: cost = 3 * e_tot
    out: list[tuple[float, float]] = []
    for (m, mw), (e, ew) in itertools.product(mants, exps):
        try:
            v = float(m) * 10.0 ** (-e)
        except ValueError:
            continue
        out.append((v, (m_tot - mw) + (e_tot - ew)))
    out.sort(key=lambda t: t[1])
    return out


def max_by_nuclide(rows: list[tuple[str, float]]) -> dict[str, float]:
    """ Largest value per nuclide over (nuclide, value) rows (e.g. over lung clearance classes) """
    out: dict[str, float] = {}
    for nuc, val in rows:
        if val is None or not math.isfinite(val) or not (MIN_DCF <= val <= MAX_DCF):
            continue
        if val > out.get(nuc, 0.0):
            out[nuc] = val
    return out


def write_csv(path: Path, data: dict[str, float], header: str) -> None:
    """ Write nuclide,value rows sorted by nuclide (%.6g keeps the printed 3 significant figures) """
    with Path(path).open("w") as f:
        f.write(header + "\n")
        for k in sorted(data):
            f.write(f"{k},{data[k]:.6g}\n")


# --------------------------------------------------------------------------------------------------
# Page images: glyph components, text lines and value cells
# --------------------------------------------------------------------------------------------------

@dataclass
class Cell:
    """ One printed 'm.mm 10^-e' value """
    mant: list  # glyph boxes (x0, y0, x1, y1) of the mantissa
    sup: list  # glyph boxes of the exponent digits (minus sign removed)
    x: float  # horizontal centre [px]
    col: int = -1  # 0..7 value column, -1 unassigned (e.g. the f1 column)
    reads: list = field(default_factory=list)  # (reader, mantissa str | None, exponent int | None)


@dataclass
class Line:
    page: int
    baseline: float
    cells: list[Cell]
    label_boxes: list  # glyph boxes of the nuclide field (empty for continuation rows)
    rest_boxes: list  # glyph boxes of the class/f1 field
    text_boxes: list  # all glyph boxes left of the value columns (element headings, footnotes)
    label: str = ''
    rest: str = ''
    text: str = ''


def _check_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found in PATH: {tool}")


def render_page(pdf: Path, page: int, dpi: int = 300) -> np.ndarray:
    """ Render one page to a deskewed boolean ink mask """
    from PIL import Image
    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / 'p'
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", str(dpi), "-gray", "-png",
                        "-singlefile", str(pdf), str(prefix)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        im = Image.open(str(prefix) + '.png').convert('L')
        im.load()
    small = im.resize((im.width // 4, im.height // 4))
    best_angle, best_score = 0.0, -1.0
    for a20 in range(-30, 31):  # -1.5 .. 1.5 degrees
        a = a20 / 20.0
        prof = (np.asarray(small.rotate(a, resample=Image.BILINEAR, fillcolor=255)) < 128).sum(axis=1)
        score = float(np.var(prof.astype(float)))
        if score > best_score:
            best_angle, best_score = a, score
    if best_angle != 0.0:
        im = im.rotate(best_angle, resample=Image.BICUBIC, fillcolor=255)
    return np.asarray(im) < 128


def glyph_boxes(ink: np.ndarray, min_area: int = 9) -> np.ndarray:
    """ Bounding boxes (x0, y0, x1, y1) of 8-connected ink components, specks removed """
    from scipy import ndimage
    lab, n = ndimage.label(ink, structure=np.ones((3, 3)))
    if n == 0:
        return np.zeros((0, 4), dtype=int)
    areas = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1))
    out = [(sl[1].start, sl[0].start, sl[1].stop, sl[0].stop)
           for i, sl in enumerate(ndimage.find_objects(lab)) if areas[i] >= min_area]
    return np.array(out, dtype=int).reshape(-1, 4)


def _text_baselines(boxes: np.ndarray) -> tuple[list[float], float]:
    """ Baselines of text lines (median bottom of full-height glyphs) and the glyph height """
    h = boxes[:, 3] - boxes[:, 1]
    w = boxes[:, 2] - boxes[:, 0]
    # Full-height glyph size from the upper quartile of text-sized components. Exponent digits are
    # ~0.7 of it, often the single most common height, and must not form a line of their own.
    cand = h[(h >= 14) & (h <= 44) & (w >= 4) & (w < 40)]
    if len(cand) == 0:
        return [], 24.0
    ref_h = float(np.percentile(cand, 75))
    main = np.where((h >= 0.82 * ref_h) & (h <= 1.45 * ref_h) & (w < 60))[0]
    glyph_h = float(np.median(h[main]))
    order = main[np.argsort(boxes[main, 3])]
    lines, cur = [], [order[0]]
    for i in order[1:]:
        if boxes[i, 3] - np.median(boxes[cur, 3]) <= 0.35 * glyph_h:
            cur.append(i)
        else:
            lines.append(cur)
            cur = [i]
    lines.append(cur)
    bases = [float(np.median(boxes[L, 3])) for L in lines if len(L) >= 3]
    # Words set mostly in lower case (element headings such as 'Samarium') have fewer than three
    # full-height glyphs. Add rows of x-height glyphs whose bottoms are not the raised exponents of
    # an existing line (those sit 0.4-0.8 glyph heights above its baseline).
    small = np.where((h >= 0.5 * ref_h) & (h < 0.82 * ref_h) & (w < 90))[0]  # bold letters merge
    if len(small):
        order = small[np.argsort(boxes[small, 3])]
        groups, cur = [], [order[0]]
        for i in order[1:]:
            if boxes[i, 3] - np.median(boxes[cur, 3]) <= 0.25 * glyph_h:
                cur.append(i)
            else:
                groups.append(cur)
                cur = [i]
        groups.append(cur)
        for g in groups:
            yb = float(np.median(boxes[g, 3]))
            n_tall = int(np.sum(np.abs(boxes[main, 3] - yb) <= 0.25 * glyph_h))
            if len(g) + n_tall < 3:
                continue
            if all(not (b - 1.0 * glyph_h <= yb <= b + 0.6 * glyph_h) for b in bases):
                bases.append(yb)
    return sorted(bases), glyph_h


def _group_words(items: list, gap: float) -> list[list]:
    items = sorted(items, key=lambda t: t[0])
    out, cur, right = [], [], None
    for t in items:
        if cur and t[0] - right > gap:
            out.append(cur)
            cur = []
        cur.append(t)
        right = t[2] if right is None or len(cur) == 1 else max(right, t[2])
    if cur:
        out.append(cur)
    return out


def _trim_minus(ink: np.ndarray, box: tuple) -> tuple:
    """ Cut a touching exponent minus sign off the left of the first exponent digit.

    The minus is a short bar in the middle band of the glyph box; columns whose ink lies only there
    are dropped. Without this, '-8' merged into one component is classified as '4'.
    """
    x0, y0, x1, y1 = box[:4]
    h = y1 - y0
    sub = ink[y0:y1, x0:x1]
    cut = 0
    for col in range(sub.shape[1]):
        rows = np.flatnonzero(sub[:, col])
        if len(rows) == 0:
            if cut:
                cut = col + 1
                continue
            break
        if rows.min() >= 0.2 * h and rows.max() <= 0.8 * h and len(rows) <= 0.3 * h:
            cut = col + 1
        else:
            break
    if (x1 - x0) <= 0.8 * h or cut < 3 or (x1 - x0 - cut) < 0.25 * h:
        return box  # a lone digit is narrower than 0.8 h; only merged '-d' glyphs are wider
    return (x0 + cut, y0, x1, y1, *box[4:])


def segment_page(ink: np.ndarray, page: int, n_columns: int = N_VALUE_COLUMNS) -> list[Line]:
    """ Find text lines, value cells, and the label fields of one table page """
    boxes = glyph_boxes(ink)
    bases, gh = _text_baselines(boxes)
    if not bases:
        return []
    B = np.array(bases)
    per_line: dict[int, list] = {k: [] for k in range(len(bases))}
    for x0, y0, x1, y1 in boxes:
        d = B - y1
        ok = np.where((d >= -0.5 * gh) & (d <= 1.1 * gh))[0]
        if len(ok) == 0:
            continue
        k = int(ok[np.argmin(np.abs(d[ok] - 0.2 * gh))])
        if y0 < B[k] - 1.9 * gh:
            continue
        rel_bot = y1 - B[k]
        hh = y1 - y0
        if rel_bot >= -0.15 * gh:
            kind = 'm'  # on the baseline (digits, letters, decimal points)
        elif rel_bot >= -0.4 * gh and hh <= 0.3 * gh:
            kind = 'dash'  # hyphen of 'Be-7'
        else:
            kind = 's'  # raised: exponent digits and minus, footnote marks
        per_line[k].append((int(x0), int(y0), int(x1), int(y1), kind))

    lines: list[Line] = []
    for k, base in enumerate(bases):
        words = _group_words(per_line[k], gap=0.42 * gh)
        cells: list[Cell] = []
        used: set[int] = set()
        for j, wd in enumerate(words):
            sup = [t for t in wd if t[4] == 's']
            mains = [t for t in wd if t[4] == 'm']
            if not sup or len(mains) < 2:
                continue
            first_sup_x = min(t[0] for t in sup)
            before = [t for t in mains if t[2] <= first_sup_x + 2]
            if len(before) < 2 or any(t[0] > first_sup_x for t in mains):
                continue
            one, zero = before[-2], before[-1]
            # '1' and '0' of the '10'; the '0' may carry a touching exponent minus (wider box).
            if one[3] - one[1] >= 0.8 * gh and zero[3] - zero[1] >= 0.8 * gh \
                    and (one[2] - one[0]) < 0.55 * gh and 0.4 * gh <= (zero[2] - zero[0]) <= 1.3 * gh:
                mant = before[:-2]
            elif zero[3] - zero[1] >= 0.8 * gh and 0.9 * gh <= (zero[2] - zero[0]) <= 1.9 * gh:
                mant = before[:-1]  # '1' and '0' merged into one component
            else:
                continue
            if not mant:
                if j == 0 or (j - 1) in used:
                    continue
                mant = [t for t in words[j - 1] if t[4] == 'm']
                used.add(j - 1)
            if not mant:
                continue
            digits = sorted((t for t in sup if (t[3] - t[1]) > 0.35 * gh), key=lambda t: t[0])
            if digits:
                digits[0] = _trim_minus(ink, digits[0])
            used.add(j)
            cells.append(Cell(mant=[t[:4] for t in mant], sup=[t[:4] for t in digits],
                              x=(min(t[0] for t in mant) + max(t[2] for t in wd)) / 2.0))
        text_items = [t for j, wd in enumerate(words) if j not in used for t in wd]
        lines.append(Line(page=page, baseline=base, cells=cells, label_boxes=[], rest_boxes=[],
                          text_boxes=[t[:4] for t in text_items]))
    _assign_columns(lines, n_columns, gh)
    return lines


def _assign_columns(lines: list[Line], n_columns: int, gh: float) -> None:
    """ Assign cells to the n_columns rightmost value columns and split label fields """
    xs = sorted(c.x for ln in lines for c in ln.cells)
    if not xs:
        return
    clusters: list[list[float]] = [[xs[0]]]
    for x in xs[1:]:
        if x - clusters[-1][-1] > 2.0 * gh:
            clusters.append([x])
        else:
            clusters[-1].append(x)
    # Keep clusters that hold a reasonable share of the cells (drops stray detections).
    n_max = max(len(c) for c in clusters)
    centers = [float(np.median(c)) for c in clusters if len(c) >= max(2, 0.15 * n_max)]
    value_centers = centers[-n_columns:]
    if len(value_centers) < 2:
        return
    spacing = float(np.median(np.diff(value_centers)))
    left_edge = value_centers[0] - 0.5 * spacing
    for ln in lines:
        for c in ln.cells:
            d = [abs(c.x - vc) for vc in value_centers]
            i = int(np.argmin(d))
            c.col = i if d[i] < 0.4 * spacing else -1
        # Label fields: glyphs left of the first value column, split at the widest gap
        items = [t for t in ln.text_boxes if t[2] < left_edge]
        items += [b for c in ln.cells if c.col < 0 for b in c.mant + c.sup]
        fields = _group_words([(*t, '') for t in items], gap=2.0 * gh)
        ln.text_boxes = [t[:4] for t in items]
        if fields:
            ln.label_boxes = [t[:4] for t in fields[0]]
            ln.rest_boxes = [t[:4] for f in fields[1:] for t in f]
    # The nuclide column starts at the leftmost label field on the page; a first field that starts
    # well to the right of it is the class/f1 field of a continuation row.
    firsts = [min(t[0] for t in ln.label_boxes) for ln in lines if ln.label_boxes and ln.cells]
    if firsts:
        col0 = float(np.percentile(firsts, 5))
        for ln in lines:
            if ln.label_boxes and min(t[0] for t in ln.label_boxes) > col0 + 3.0 * gh:
                ln.rest_boxes = ln.label_boxes + ln.rest_boxes
                ln.label_boxes = []


# --------------------------------------------------------------------------------------------------
# OCR of stacked crops
# --------------------------------------------------------------------------------------------------

ROW_PITCH: int = 64


def _paste(canvas: np.ndarray, patch: np.ndarray, x: int, y_bottom: int) -> int:
    h, w = patch.shape
    w = min(w, canvas.shape[1] - x)
    y0 = max(0, y_bottom - h)
    canvas[y0:y_bottom, x:x + w] |= patch[h - (y_bottom - y0):, :w]
    return x + w


def _crop(ink: np.ndarray, boxes: list, pad: int = 1) -> np.ndarray:
    xa = max(min(b[0] for b in boxes) - pad, 0)
    xb = max(b[2] for b in boxes) + pad
    ya = max(min(b[1] for b in boxes) - pad, 0)
    yb = max(b[3] for b in boxes) + pad
    return ink[ya:yb, xa:xb]


def _scaled(patch: np.ndarray, scale: float) -> np.ndarray:
    from PIL import Image
    if scale == 1.0:
        return patch
    im = Image.fromarray((patch * 255).astype(np.uint8))
    im = im.resize((max(1, int(round(patch.shape[1] * scale))), max(1, int(round(patch.shape[0] * scale)))),
                   Image.LANCZOS)
    return np.asarray(im) > 110


def ocr_stack(ink: np.ndarray, groups: list[list], whitelist: str | None, spaced: bool = False,
              scale: float = 1.0, width: int = 1400, chunk: int = 150) -> list[str]:
    """ OCR many glyph groups at once: one group per canvas row, rows mapped back by position.

    spaced=True places each glyph separately with a gap, which lets tesseract read tight pairs
    such as the exponent '11'.
    """
    from PIL import Image
    out: list[str] = []
    for start in range(0, len(groups), chunk):
        part = groups[start:start + chunk]
        pitch = int(ROW_PITCH * max(1.0, scale))
        canvas = np.zeros((pitch * len(part) + 40, width), dtype=bool)
        for r, g in enumerate(part):
            if not g:
                continue
            y_bottom = 20 + r * pitch + int(0.75 * pitch)
            if spaced:
                x = 30
                ybot = max(b[3] for b in g)
                for b in sorted(g, key=lambda b: b[0]):
                    patch = _scaled(_crop(ink, [b], pad=0), scale)
                    lift = int(round((ybot - b[3]) * scale))
                    x = _paste(canvas, patch, x, y_bottom - lift) + int(10 * scale)
            else:
                _paste(canvas, _scaled(_crop(ink, g), scale), 30, y_bottom)
        img = Image.fromarray(np.where(canvas, 0, 255).astype(np.uint8))
        img = img.resize((img.width * 2, img.height * 2), Image.BICUBIC)
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, 'c.png')
            img.save(p)
            cmd = ['tesseract', p, '-', '--psm', '6']
            if whitelist is not None:
                cmd += ['-c', f'tessedit_char_whitelist={whitelist}']
            # One OpenMP thread per tesseract: pages already run in parallel processes.
            env = dict(os.environ, OMP_THREAD_LIMIT='1')
            tsv = subprocess.run(cmd + ['tsv'], capture_output=True, text=True, env=env).stdout
        rows: list[list[tuple[int, str]]] = [[] for _ in part]
        for rec in tsv.splitlines()[1:]:
            f = rec.split('\t')
            if len(f) < 12 or f[0] != '5' or not f[11].strip():
                continue
            cy = (int(f[7]) + int(f[9]) / 2.0) / 2.0
            r = int((cy - 20) // pitch)
            if 0 <= r < len(part):
                rows[r].append((int(f[6]), f[11].strip()))
        out.extend(' '.join(w for _, w in sorted(rw)) for rw in rows)
    return out


def _glyph_vector(ink: np.ndarray, b, H: int = 20, W: int = 18) -> np.ndarray:
    patch = _crop(ink, [b], pad=0).astype(np.float32)
    h, w = patch.shape
    nw = max(1, min(W, int(round(w * H / h))))
    from PIL import Image
    im = Image.fromarray((patch * 255).astype(np.uint8)).resize((nw, H), Image.BILINEAR)
    out = np.zeros((H, W), dtype=np.float32)
    x = (W - nw) // 2
    out[:, x:x + nw] = np.asarray(im, dtype=np.float32) / 255.0
    return np.concatenate([out.ravel(), [2.0 * w / h]])


class DigitClassifier:
    """ k-nearest-neighbour digit classifier trained on glyphs of one scan """

    def __init__(self, vectors: np.ndarray, labels: np.ndarray, k: int = 5, per_class: int = 600):
        rng = np.random.default_rng(0)
        keep = []
        for d in np.unique(labels):
            idx = np.where(labels == d)[0]
            keep.extend(rng.choice(idx, size=min(per_class, len(idx)), replace=False))
        keep = np.array(sorted(keep))
        self.X = vectors[keep]
        self.y = labels[keep]
        self.k = k
        self.x2 = (self.X ** 2).sum(axis=1)

    def predict(self, vectors: np.ndarray) -> list[str]:
        if len(vectors) == 0:
            return []
        d = self.x2[None, :] - 2.0 * vectors @ self.X.T + (vectors ** 2).sum(axis=1)[:, None]
        nn = np.argsort(d, axis=1)[:, :self.k]
        out = []
        for row in nn:
            vals, cnt = np.unique(self.y[row], return_counts=True)
            out.append(str(vals[np.argmax(cnt)]))
        return out


# --------------------------------------------------------------------------------------------------
# Reading one PDF edition
# --------------------------------------------------------------------------------------------------

def _read_page(args: tuple) -> tuple[list[Line], list[tuple[np.ndarray, str]], list]:
    """ Worker: segment and OCR one page; returns lines, digit training samples, raw glyph data """
    pdf, page, n_columns = args
    ink = render_page(Path(pdf), page)
    lines = segment_page(ink, page, n_columns)
    cells = [c for ln in lines for c in ln.cells if c.col >= 0]
    mant_txt = ocr_stack(ink, [c.mant for c in cells], '0123456789.')
    exp_txt = ocr_stack(ink, [c.sup for c in cells], '0123456789', spaced=True, scale=2.0)
    train: list[tuple[np.ndarray, str]] = []
    sup_vectors = []
    for c, mt, et in zip(cells, mant_txt, exp_txt):
        mt = mt.replace(' ', '')
        et = et.replace(' ', '')
        ms = mt if re.fullmatch(r'\d\.\d\d', mt) else None
        es = int(et) if re.fullmatch(r'\d{1,2}', et) else None
        c.reads.append(('tess', ms, es))
        digits = sorted([b for b in c.mant if (b[3] - b[1]) > 12], key=lambda b: b[0])
        if ms and len(digits) == 3:
            for b, ch in zip(digits, ms.replace('.', '')):
                train.append((_glyph_vector(ink, b), ch))
        sup_vectors.append(([_glyph_vector(ink, b) for b in sorted(c.sup, key=lambda b: b[0])],
                            [_glyph_vector(ink, b) for b in digits]))
    # Labels, class fields and element headings
    lab_groups, lab_ref = [], []
    for ln in lines:
        for attr, boxes in (('label', ln.label_boxes), ('rest', ln.rest_boxes)):
            if boxes:
                lab_groups.append(boxes)
                lab_ref.append((ln, attr))
        if not ln.cells and ln.text_boxes:
            lab_groups.append(ln.text_boxes)
            lab_ref.append((ln, 'text'))
    for (ln, attr), txt in zip(lab_ref, ocr_stack(ink, lab_groups, None, width=2600)):
        setattr(ln, attr, txt)
    for ln in lines:  # drop glyph boxes before pickling back to the parent
        ln.label_boxes, ln.rest_boxes, ln.text_boxes = [], [], []
    return lines, train, sup_vectors


def read_edition(pdf: Path, edition: Edition, table: str, workers: int) -> list[Line]:
    """ Read all pages of one table ('2.1' or '2.3') from one PDF edition """
    pages = edition.inhalation_pages if table == '2.1' else edition.submersion_pages
    jobs = [(str(pdf), p, N_VALUE_COLUMNS) for p in pages]
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as ex:
        results = list(ex.map(_read_page, jobs))
    all_lines: list[Line] = []
    train_X, train_y = [], []
    page_vectors = []
    for lines, train, vecs in results:
        all_lines.extend(lines)
        for v, ch in train:
            train_X.append(v)
            train_y.append(ch)
        page_vectors.append((lines, vecs))
    if not train_X:
        return all_lines
    clf = DigitClassifier(np.array(train_X), np.array(train_y))
    for lines, vecs in page_vectors:
        cells = [c for ln in lines for c in ln.cells if c.col >= 0]
        for c, (sup_v, mant_v) in zip(cells, vecs):
            e = ''.join(clf.predict(np.array(sup_v))) if sup_v else ''
            m = ''.join(clf.predict(np.array(mant_v))) if len(mant_v) == 3 else ''
            c.reads.append(('knn', f'{m[0]}.{m[1:]}' if m else None,
                            int(e) if re.fullmatch(r'\d{1,2}', e) else None))
    return all_lines


def text_layer_pages(pdf: Path, pages: tuple[int, ...]) -> dict[int, list[list[str]]]:
    """ Mantissa sequences of the pdftotext -layout lines of each page """
    out: dict[int, list[list[str]]] = {}
    for p in pages:
        txt = subprocess.run(['pdftotext', '-layout', '-f', str(p), '-l', str(p), str(pdf), '-'],
                             capture_output=True, text=True).stdout
        out[p] = [m for m in (text_layer_mantissas(ln) for ln in txt.splitlines()) if len(m) >= 2]
    return out


def _align(seq_a: list[list[str]], seq_b: list[list[str]]) -> list[tuple[int, int]]:
    """ Monotonic alignment of two line lists by the number of shared mantissas in order """
    def score(a: list[str], b: list[str]) -> int:
        n, m = len(a), len(b)
        dp = [[0] * (m + 1) for _ in range(n + 1)]
        for i in range(n):
            for j in range(m):
                dp[i + 1][j + 1] = dp[i][j] + 1 if a[i] == b[j] else max(dp[i][j + 1], dp[i + 1][j])
        return dp[n][m]
    n, m = len(seq_a), len(seq_b)
    S = [[0] * (m + 1) for _ in range(n + 1)]
    back = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = score(seq_a[i - 1], seq_b[j - 1])
            opts = [(S[i - 1][j], 'a'), (S[i][j - 1], 'b')]
            if s >= max(2, len(seq_a[i - 1]) // 2):
                opts.append((S[i - 1][j - 1] + s, 'ab'))
            S[i][j], back[i][j] = max(opts)
    pairs = []
    i, j = n, m
    while i > 0 and j > 0:
        mv = back[i][j]
        if mv == 'ab':
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif mv == 'a':
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def _line_mantissas(ln: Line) -> list[str]:
    out = []
    for c in sorted((c for c in ln.cells if c.col >= 0), key=lambda c: c.col):
        ms = [r[1] for r in c.reads if r[1]]
        out.append(max(set(ms), key=ms.count) if ms else '?')
    return out


def merge_second_reading(primary: list[Line], secondary: list[Line], page_offset: int, tag: str) -> None:
    """ Attach the cell readings of a second scan to the aligned rows of the primary scan """
    by_page_p: dict[int, list[Line]] = {}
    by_page_s: dict[int, list[Line]] = {}
    for ln in primary:
        if ln.cells:
            by_page_p.setdefault(ln.page, []).append(ln)
    for ln in secondary:
        if ln.cells:
            by_page_s.setdefault(ln.page + page_offset, []).append(ln)
    for page, lines_p in by_page_p.items():
        lines_s = by_page_s.get(page, [])
        pairs = _align([_line_mantissas(ln) for ln in lines_p], [_line_mantissas(ln) for ln in lines_s])
        for i, j in pairs:
            cols_s = {c.col: c for c in lines_s[j].cells if c.col >= 0}
            cols_p = {c.col for c in lines_p[i].cells if c.col >= 0}
            for c in lines_p[i].cells:
                if c.col in cols_s:
                    c.reads.extend((f'{tag}-{r[0]}', r[1], r[2]) for r in cols_s[c.col].reads)
            for col, cs in cols_s.items():  # cells the primary scan missed
                if col not in cols_p:
                    lines_p[i].cells.append(Cell(mant=[], sup=[], x=cs.x, col=col,
                                                 reads=[(f'{tag}-{r[0]}', r[1], r[2]) for r in cs.reads]))


def merge_text_layer(primary: list[Line], text_pages: dict[int, list[list[str]]], page_offset: int,
                     tag: str) -> None:
    """ Attach text-layer mantissas to primary rows whose aligned text line has one per column """
    by_page: dict[int, list[Line]] = {}
    for ln in primary:
        if ln.cells:
            by_page.setdefault(ln.page, []).append(ln)
    for page, lines_p in by_page.items():
        seqs = text_pages.get(page - page_offset, [])
        pairs = _align([_line_mantissas(ln) for ln in lines_p], seqs)
        for i, j in pairs:
            cells = sorted((c for c in lines_p[i].cells if c.col >= 0), key=lambda c: c.col)
            if len(seqs[j]) == len(cells):
                for c, m in zip(cells, seqs[j]):
                    c.reads.append((f'{tag}-text', m, None))


# Vote weights by reader. The OSTI scan is the sharper one.
READER_WEIGHTS: dict[str, float] = {
    'tess': 1.0, 'knn': 1.0, 'epa-tess': 0.6, 'epa-knn': 0.6, 'osti-text': 0.5, 'epa-text': 0.4,
}


def _votes(cell: Cell) -> tuple[dict[str, float], dict[int, float]]:
    mv: dict[str, float] = {}
    ev: dict[int, float] = {}
    for reader, m, e in cell.reads:
        w = READER_WEIGHTS.get(reader, 0.5)
        if m:
            mv[m] = mv.get(m, 0.0) + w
        if e is not None and VOTE_EXPONENT_RANGE[0] <= e <= VOTE_EXPONENT_RANGE[1]:
            ev[e] = ev.get(e, 0.0) + w
    return mv, ev


@dataclass
class TableRow:
    page: int
    nuclide: str | None
    form: str
    values: list
    consistent: bool
    cost: float
    source: str = 'auto'
    label: str = ''
    changed: tuple = ()  # columns whose value differs from the most-voted reading


def build_rows(lines: list[Line], table: str) -> list[TableRow]:
    """ Solve every value row and attach nuclide names from the labels and element headings """
    rows: list[TableRow] = []
    element: str | None = None
    nuclide: str | None = None
    label_row = 0
    for ln in lines:
        if not ln.cells:
            words = ln.text.split()
            hl = half_life_label(ln.text)
            if hl and nuclide and rows:
                # Isomer half-life printed on its own line below the rows it names ('Eu-150' / '12.62 h')
                iso = HALF_LIFE_ISOMERS.get((nuclide, hl))
                nuclide = iso if iso else f'?{nuclide} ({hl})'
                for r in rows[label_row:]:
                    r.nuclide = nuclide
            elif words and len(words) <= 2 and re.fullmatch(r'[A-Zl][a-z]{2,}\*?', words[0]):
                sym = element_from_heading(words[0])  # element heading ('Hydrogen', 'Hydrogen*')
                if sym:
                    element = sym
            continue
        cols = {c.col: c for c in ln.cells if c.col >= 0}
        if not cols:
            continue
        if 7 not in cols and set(cols) <= {6} and not ln.label:
            continue  # limiting-organ annotation rows ('ST wall', 'Skin', 'Lens'): one value, no h_E
        hl = half_life_label(ln.label) if ln.label else None
        if hl:
            # Second label line of an isomer pair printed as 'Nb-89' / '122 m': it names the rows
            # since the last nuclide label, including the current one.
            iso = HALF_LIFE_ISOMERS.get((nuclide, hl)) if nuclide else None
            nuclide = iso if iso else f'?{nuclide} ({hl})'
            for r in rows[label_row:]:
                r.nuclide = nuclide
        elif ln.label:
            label_row = len(rows)
            label_sym = None
            m = re.match(r'\s*([A-Za-z]{1,2})\s*-', ln.label)
            if m:
                label_sym = element_symbol(m.group(1))
            if element is None:
                element = label_sym
            elif label_sym and label_sym != element and _label_confirms(ln.label, label_sym) \
                    and ELEMENT_Z[label_sym.lower()] > ELEMENT_Z[element.lower()] \
                    and parse_nuclide_label(ln.label, element) is None:
                # Element heading missed: the label reads cleanly as a later element and its mass
                # number is implausible for the current one.
                element = label_sym
            nuc = parse_nuclide_label(ln.label, element)
            if nuc:
                nuclide = nuc
            else:
                nuclide = f'?{ln.label.strip()}'
        options = []
        for i in range(N_VALUE_COLUMNS):
            if i in cols:
                mv, ev = _votes(cols[i])
                options.append(cell_options(mv, ev))
            else:
                options.append([])
        first = [opt[0][0] if opt else None for opt in options]
        values, cost, ok = solve_row(options)
        if not ok:
            # Second pass: any exponent for up to two cells (exponent digits shared by all readers
            # can be misread, e.g. a bold '8' as '3').
            wide = []
            for i in range(N_VALUE_COLUMNS):
                if i in cols:
                    mv, ev = _votes(cols[i])
                    opts = cell_options(mv, ev)
                    tot = sum(ev.values()) or 1.0
                    have = {round(math.log10(v), 6) for v, _ in opts if v > 0.0}
                    for m in sorted(mv, key=lambda k: -mv[k])[:1]:
                        for e in range(EXPONENT_RANGE[0], EXPONENT_RANGE[1] + 1):
                            v = float(m) * 10.0 ** (-e)
                            if v > 0.0 and round(math.log10(v), 6) not in have:
                                opts.append((v, 4.0 * tot))
                    wide.append(opts)
                else:
                    wide.append([])
            values2, cost2, ok2 = solve_row(wide, max_changes=2)
            if ok2:
                values, cost, ok = values2, cost2, ok2
        changed = tuple(i for i in range(N_VALUE_COLUMNS) if values[i] != first[i])
        if table == '2.1' and changed and len(cols) < N_VALUE_COLUMNS:
            # Every Table 2.1 row prints all 8 values. With a cell missing, the check is too weak to
            # justify altering the other readings, so the row needs a manual reading.
            ok = False
        rows.append(TableRow(page=ln.page, nuclide=nuclide, form=ln.rest.strip(), values=values,
                             consistent=ok, cost=cost, label=ln.label.strip(), changed=changed))
    return rows


def _label_confirms(label: str, sym: str) -> bool:
    # A label symbol overrides the heading only when it reads cleanly as that symbol.
    return bool(re.match(rf'\s*{sym}\s*-\s*\d', label))


def apply_manual_corrections(rows: list[TableRow], table: str) -> None:
    counters: dict[str, int] = {}
    for r in rows:
        if r.nuclide is None:
            continue
        idx = counters.get(r.nuclide, 0)
        counters[r.nuclide] = idx + 1
        key = (table, r.nuclide, idx)
        if key in MANUAL_CORRECTIONS:
            vals, page = MANUAL_CORRECTIONS[key]
            r.values = list(vals)
            r.consistent = row_is_consistent(r.values)
            r.source = f'manual (p. {page})'


def write_audit(path: Path, rows: list[TableRow]) -> None:
    with Path(path).open('w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['page', 'nuclide', 'label', 'form', 'gonad', 'breast', 'lung', 'r_marrow', 'b_surface',
                     'thyroid', 'remainder', 'effective', 'weighted_sum', 'consistent', 'cost', 'changed',
                     'source'])
        for r in rows:
            vals = [f'{v:.3g}' if v else '' for v in r.values]
            wsum = effective_from_organs(r.values[:7])
            wr.writerow([r.page, r.nuclide, r.label, r.form, *vals, f'{wsum:.4g}', r.consistent,
                         f'{r.cost:.2f}', ' '.join(str(c) for c in r.changed), r.source])


def extract_table(pdfs: dict[str, Path], table: str, workers: int) -> list[TableRow]:
    """ Read one table from all available editions and return the solved rows """
    if 'osti' in pdfs:
        primary_name, secondary_name = 'osti', 'epa'
    else:
        primary_name, secondary_name = 'epa', None
    ed_p = EDITIONS[primary_name]
    primary = read_edition(pdfs[primary_name], ed_p, table, workers)
    pages_p = ed_p.inhalation_pages if table == '2.1' else ed_p.submersion_pages
    if secondary_name and secondary_name in pdfs:
        ed_s = EDITIONS[secondary_name]
        secondary = read_edition(pdfs[secondary_name], ed_s, table, workers)
        pages_s = ed_s.inhalation_pages if table == '2.1' else ed_s.submersion_pages
        offset = pages_p[0] - pages_s[0]
        merge_second_reading(primary, secondary, offset, secondary_name)
        merge_text_layer(primary, text_layer_pages(pdfs[secondary_name], pages_s), offset, secondary_name)
    merge_text_layer(primary, text_layer_pages(pdfs[primary_name], pages_p), 0, primary_name)
    rows = build_rows(primary, table)
    apply_manual_corrections(rows, table)
    return rows


def detect_editions(paths: list[Path]) -> dict[str, Path]:
    """ Map edition name -> PDF path by page count """
    out: dict[str, Path] = {}
    for p in paths:
        info = subprocess.run(['pdfinfo', str(p)], capture_output=True, text=True).stdout
        m = re.search(r'Pages:\s+(\d+)', info)
        n = int(m.group(1)) if m else -1
        for ed in EDITIONS.values():
            if ed.n_pages == n and ed.name not in out:
                out[ed.name] = p
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Extract FGR-11 inhalation and submersion DCFs.")
    parser.add_argument("--pdf", action="append", default=None,
                        help="FGR-11 PDF (repeat for a second edition). Default: known files under PDF/")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for the two CSVs")
    parser.add_argument("--audit", default=None,
                        help="Directory for per-row audit CSVs (fgr11_table21_rows.csv, fgr11_table23_rows.csv)")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args(argv)

    for tool in ("pdftoppm", "pdftotext", "pdfinfo", "tesseract"):
        _check_tool(tool)
    paths = [Path(p) for p in args.pdf] if args.pdf else \
        [DEFAULT_PDF_DIR / ed.default_file for ed in EDITIONS.values() if (DEFAULT_PDF_DIR / ed.default_file).is_file()]
    pdfs = detect_editions([p for p in paths if p.is_file()])
    if not pdfs:
        print(f"No FGR-11 PDF found (looked for: {', '.join(str(p) for p in paths)})", file=sys.stderr)
        return 1
    print(f"Editions: {', '.join(f'{k}={v}' for k, v in pdfs.items())}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    for table, fname, header, scale in (('2.1', 'dcf_fgr11_inhalation_worker.csv', INHALATION_HEADER, 1.0),
                                        ('2.3', 'dcf_fgr11_submersion_worker.csv', SUBMERSION_HEADER,
                                         SUBMERSION_HR_TO_DAY)):
        rows = extract_table(pdfs, table, args.workers)
        if args.audit:
            Path(args.audit).mkdir(parents=True, exist_ok=True)
            write_audit(Path(args.audit) / f"fgr11_table{table.replace('.', '')}_rows.csv", rows)
        bad = [r for r in rows if not r.consistent or not r.nuclide or r.nuclide.startswith('?')]
        for r in bad:
            print(f"  check: Table {table} p.{r.page} {r.nuclide} {r.form!r} values={r.values}")
        good = [(r.nuclide, r.values[7]) for r in rows
                if r.consistent and r.nuclide and not r.nuclide.startswith('?') and r.values[7]]
        data = {k: v * scale for k, v in max_by_nuclide(good).items()}
        write_csv(out_dir / fname, data, header)
        results[table] = (len(rows), len(bad), len(data))
        print(f"Table {table}: {len(rows)} rows, {len(bad)} unresolved, {len(data)} nuclides -> {out_dir / fname}")
    return 0 if all(v[1] == 0 for v in results.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
