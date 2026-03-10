#!/usr/bin/env python3
"""
Extract EPA FGR-11 dose conversion factors (DCF) from the PDF using OCR.

Outputs:
  - leaky_box_origen/data/dcf_fgr11_inhalation_worker.csv
  - leaky_box_origen/data/dcf_fgr11_submersion_worker.csv

This script uses pdftoppm + tesseract. It applies conservative parsing:
  - Inhalation: max DCF across lung clearance classes (D/W/Y) per nuclide.
  - Submersion: effective dose rate per unit concentration (Sv/hr per Bq/m^3),
    converted to Sv/day per Bq/m^3.

Notes:
  - FGR-11 tables are multi-page and OCR is imperfect. We use heuristics
    to decode exponent suffixes and filter non-physical outliers.
  - Review the output if you plan to use the DCFs for licensing decisions.
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


INHAL_PAGES = list(range(126, 159))  # Table 2.1 inhalation pages
SUB_PAGES = list(range(184, 186))    # Table 2.3 submersion pages

# Heuristic exponent maps for OCR suffixes.
INHAL_EXP_MAP = {'°': 9, '%': 8, '"': 10, "'": 7, '!': 8, '?': 9}
SUB_EXP_MAP = {'°"?': 12, '°"': 11, '?)': 10, '"': 12, '°': 11, '%': 11}

# Sanity caps to drop OCR artifacts.
MAX_DCF_SV_BQ = 1e-2
MAX_DCF_SV_DAY = 1e-2


def _check_tool(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Required tool not found in PATH: {tool}")


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _ocr_pages(pdf: Path, pages: list[int], dpi: int, work: Path) -> list[Path]:
    txt_paths: list[Path] = []
    work.mkdir(parents=True, exist_ok=True)
    for page in pages:
        img_prefix = work / f"p{page:03d}"
        img = img_prefix.with_suffix(".png")
        if not img.exists():
            _run(["pdftoppm", "-f", str(page), "-l", str(page),
                  "-r", str(dpi), "-png", "-singlefile", str(pdf), str(img_prefix)])
        txt_path = img_prefix.with_suffix(".txt")
        if not txt_path.exists():
            _run(["tesseract", str(img), str(img_prefix), "-l", "eng", "--psm", "6"])
        txt_paths.append(txt_path)
    return txt_paths


def _parse_mantissa(token: str) -> float | None:
    token = token.strip().replace(",", ".")
    if "." in token:
        try:
            return float(token)
        except ValueError:
            return None
    if token.isdigit():
        if len(token) == 3:
            return float(token[0] + "." + token[1:])
        if len(token) == 2:
            return float(token[0] + "." + token[1])
    try:
        return float(token)
    except ValueError:
        return None


def _exp_from_suffix(suffix: str, mapping: dict[str, int]) -> int | None:
    suffix = suffix.strip()
    digits = "".join(ch for ch in suffix if ch.isdigit())
    if digits:
        try:
            return int(digits)
        except ValueError:
            pass
    for key in sorted(mapping.keys(), key=lambda x: -len(x)):
        if suffix == key:
            return mapping[key]
    for ch in suffix:
        if ch in mapping:
            return mapping[ch]
    return None


def _parse_last_value(line: str, mapping: dict[str, int]) -> float | None:
    # Find last occurrence of mantissa + 10 + suffix (OCR forms like "8.89 10°").
    matches = re.findall(r"(\d+(?:\.\d+)?)\s*10(\S+)", line)
    if not matches:
        return None
    mantissa, suffix = matches[-1]
    mant = _parse_mantissa(mantissa)
    if mant is None:
        return None
    exp = _exp_from_suffix(suffix, mapping)
    if exp is None:
        return None
    return mant * (10 ** (-exp))


def extract_inhalation(txt_paths: list[Path]) -> dict[str, float]:
    nuclide_re = re.compile(r"\b([A-Z][a-z]?|1)\s*[-–]?\s*(\d{2,3}m?)\b")
    dcf: dict[str, float] = {}
    last_nuclide: str | None = None
    for txt_path in txt_paths:
        for raw_line in txt_path.read_text(errors="ignore").splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.lower().startswith("table") or "Committed" in line or "Nuclide" in line:
                continue
            line = line.replace("—", "-").replace("–", "-")
            m = nuclide_re.search(line)
            if m:
                elem, mass = m.group(1), m.group(2)
                if elem == "1":
                    elem = "I"
                nuclide = f"{elem}-{mass}".lower()
                last_nuclide = nuclide
            else:
                if last_nuclide and re.match(r"^[DWMYV]\b", line):
                    nuclide = last_nuclide
                else:
                    continue
            val = _parse_last_value(line, INHAL_EXP_MAP)
            if val is None or not math.isfinite(val) or val <= 0.0 or val > MAX_DCF_SV_BQ:
                continue
            prev = dcf.get(nuclide)
            if prev is None or val > prev:
                dcf[nuclide] = val
    return dcf


def extract_submersion(txt_paths: list[Path]) -> dict[str, float]:
    # Uses effective dose rate per unit concentration (Sv/hr per Bq/m^3),
    # then converts to Sv/day per Bq/m^3.
    nuclide_re = re.compile(r"\b([A-Z][a-z]?|1)\s*[-–]?\s*(\d{2,3}m?)\b")
    dcf: dict[str, float] = {}
    for txt_path in txt_paths:
        for raw_line in txt_path.read_text(errors="ignore").splitlines():
            line = raw_line.strip()
            if not line:
                continue
            if line.lower().startswith("table") or "Dose Equivalent" in line:
                continue
            line = line.replace("—", "-").replace("–", "-")
            m = nuclide_re.search(line)
            if not m:
                continue
            elem, mass = m.group(1), m.group(2)
            if elem == "1":
                elem = "I"
            nuclide = f"{elem}-{mass}".lower()
            val = _parse_last_value(line, SUB_EXP_MAP)
            if val is None:
                continue
            val *= 24.0  # Sv/hr -> Sv/day
            if not math.isfinite(val) or val <= 0.0 or val > MAX_DCF_SV_DAY:
                continue
            dcf[nuclide] = val
    return dcf


def write_csv(path: Path, data: dict[str, float], header: str) -> None:
    with path.open("w") as f:
        f.write(header + "\n")
        for k in sorted(data):
            f.write(f"{k},{data[k]:.6e}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract FGR-11 DCFs via OCR.")
    default_pdf = Path(__file__).resolve().parent.parent / "PDF" / "EPA 1988_FGR11_0.pdf"
    parser.add_argument("--pdf", default=str(default_pdf),
                        help="Path to FGR-11 PDF")
    parser.add_argument("--dpi", type=int, default=300, help="OCR resolution (dpi)")
    parser.add_argument("--keep-ocr", action="store_true",
                        help="Keep OCR intermediates under ./_fgr11_ocr")
    args = parser.parse_args()

    pdf = Path(args.pdf)
    if not pdf.exists():
        print(f"PDF not found: {pdf}", file=sys.stderr)
        return 1

    _check_tool("pdftoppm")
    _check_tool("tesseract")

    out_dir = Path(__file__).parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.keep_ocr:
        work = Path("_fgr11_ocr")
        work.mkdir(parents=True, exist_ok=True)
        txt_inh = _ocr_pages(pdf, INHAL_PAGES, args.dpi, work / "inhalation")
        txt_sub = _ocr_pages(pdf, SUB_PAGES, args.dpi, work / "submersion")
    else:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            txt_inh = _ocr_pages(pdf, INHAL_PAGES, args.dpi, work / "inhalation")
            txt_sub = _ocr_pages(pdf, SUB_PAGES, args.dpi, work / "submersion")

            inhal = extract_inhalation(txt_inh)
            sub = extract_submersion(txt_sub)
            write_csv(out_dir / "dcf_fgr11_inhalation_worker.csv", inhal, "nuclide,dcf_sv_bq")
            write_csv(out_dir / "dcf_fgr11_submersion_worker.csv", sub, "nuclide,dcf_sv_per_bq_m3_day")

            print(f"Inhalation nuclides: {len(inhal)}")
            print(f"Submersion nuclides: {len(sub)}")
            return 0

    # keep-ocr path
    inhal = extract_inhalation(txt_inh)
    sub = extract_submersion(txt_sub)
    write_csv(out_dir / "dcf_fgr11_inhalation_worker.csv", inhal, "nuclide,dcf_sv_bq")
    write_csv(out_dir / "dcf_fgr11_submersion_worker.csv", sub, "nuclide,dcf_sv_per_bq_m3_day")
    print(f"Inhalation nuclides: {len(inhal)}")
    print(f"Submersion nuclides: {len(sub)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
