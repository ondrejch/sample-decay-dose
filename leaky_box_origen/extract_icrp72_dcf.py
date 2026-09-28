#!/usr/bin/env python3
"""
Extract ICRP Publication 72 (Annals of the ICRP 26(1), 1996) adult dose coefficients from the scanned
compilation ANIB_26_1.pdf.

Outputs (in leaky_box_origen/data/ unless --out-dir is given), all from Table A.2 (adult committed
effective dose per unit intake e(50) [Sv/Bq]) except the immersion table:
  - dcf_icrp72_inhalation_adult.csv       nuclide,dcf_sv_bq
        Type F where Table A.2 lists it, otherwise the fastest listed type (M, then S).
  - dcf_icrp72_anib26_1_adult.csv         nuclide,dcf_sv_bq
        Type M where Table A.2 lists it, otherwise the fastest listed type.
  - dcf_icrp72_inhalation_adult_max.csv   nuclide,dcf_sv_bq
        The largest value over the listed types (F, M, S). This is the conservative choice when the
        chemical form of the release is unknown.
  - dcf_icrp72_immersion_adult.csv        nuclide,dcf_sv_per_bq_m3_day
        Table A.4 effective dose rate for adults immersed in inert gases [Sv/day per Bq/m^3].
Mercury rows list organic and inorganic forms separately; the larger value of a type is kept.
The three inhalation tables also carry H-3 as tritiated water (HTO) vapour from Table A.3 (Class SR-2),
1.8E-11 Sv/Bq for adults. Table A.2 lists 'tritium compounds' (particulate forms), which are not used.
Table A.4 has no tritium entry.

The table pages are typewritten (monospace, 'd.dE-dd'), so tesseract reads them well. Each page is
read three times (300, 450 and 600 dpi) and the readings are merged row by row. A value is accepted
when the majority of readings agree, or when a single candidate lies within a factor of 3 of the
15-year value printed next to it. NAME_FIXES and VALUE_FIXES hold readings checked by eye on the page images.

Requirements: pdftoppm, tesseract (eng), Pillow.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from leaky_box_origen.extract_fgr11_dcf import (ELEMENT_Z, _group_words, element_from_heading,  # noqa: E402
                                                element_symbol, glyph_boxes, is_plausible_nuclide,
                                                ocr_stack, render_page, write_csv)

REPO_DIR: Path = Path(__file__).resolve().parent.parent
DEFAULT_PDF: Path = REPO_DIR / "PDF" / "ANIB_26_1.pdf"
DEFAULT_OUT_DIR: Path = Path(__file__).resolve().parent / "data"
TABLE_A2_PAGES: tuple[int, ...] = tuple(range(58, 98))  # PDF pages (printed pages 44-83)
TABLE_A4_PAGE: int = 103  # printed page 89
N_TABLE_A4: int = 26  # Ar-37 ... Xe-138

# Table A.3 (soluble or reactive gases and vapours, PDF p. 98 = printed p. 84), adult column.
HTO_ADULT_SV_BQ: float = 1.8e-11  # tritiated water vapour

# Nuclides that ICRP-72 separates by half-life instead of an 'm' suffix, and OCR misreads of the
# mass number that the printed half-life identifies. (name as parsed, half-life as printed) -> name.
NAME_FIXES: dict[tuple[str, str], str] = {
    ('nb-89', '2.03 h'): 'nb-89', ('nb-89', '1.10 h'): 'nb-89m',
    ('in-110', '4.90 h'): 'in-110', ('in-110', '1.15 h'): 'in-110m',
    ('sb-120', '5.76 d'): 'sb-120m', ('sb-120', '0.265 h'): 'sb-120',
    ('sb-128', '9.01 h'): 'sb-128', ('sb-128', '0.173 h'): 'sb-128m',
    ('eu-150', '34.2 y'): 'eu-150', ('eu-150', '12.6 h'): 'eu-150m',
    ('tb-156m', '1.02 d'): 'tb-156m', ('tb-156m', '5.00 h'): 'tb-156m2',
    ('re-182', '2.67 d'): 're-182', ('re-182', '12.7 h'): 're-182m',
    ('ir-186', '15.6 h'): 'ir-186', ('ir-186', '15.8 h'): 'ir-186', ('ir-186', '1.75 h'): 'ir-186m',
    ('ir-190m', '3.10 h'): 'ir-190m2', ('ir-190m', '1.20 h'): 'ir-190m',
    ('np-236', '1.15E+05 y'): 'np-236', ('np-236', '22.5 h'): 'np-236m',
    # Mass numbers misread by OCR; the printed half-life identifies the nuclide.
    ('rb-61', '4.58 h'): 'rb-81', ('zr-68', '83.4 d'): 'zr-88', ('hf-161', '42.4 d'): 'hf-181',
    ('ta-162m', '0.264 h'): 'ta-182m', ('re-164', '38.0 d'): 're-184', ('re-164', '386.0 d'): 're-184',
    ('os-165', '94.0 d'): 'os-185', ('re-186m', '0.310 h'): 're-188m', ('ir-186', '1.73 d'): 'ir-188',
    ('tl-202', '3.04 d'): 'tl-201', ('zn-711m', '3.92 h'): 'zn-71m', ('nb-33m', '13.6 y'): 'nb-93m',
    ('ge-7', ''): 'ge-77',  # between Ge-75 and Ge-78 on p. 63 (printed 49)
    ('tb-156m', ''): 'tb-156m2',  # second Tb-156m entry (5.00 h) on p. 81; its half-life is not read
}

# Adult e(50) [Sv/Bq] read by eye on the page images for rows whose OCR readings disagree, did not
# parse, or give an adult/15-year ratio outside RATIO_RANGE. Key: (PDF page, nuclide, type,
# occurrence of that nuclide/type on the page; mercury lists organic and inorganic forms).
VALUE_FIXES: dict[tuple[int, str, str, int], float] = {
    (59, 'ca-47', 'F', 0): 5.5e-10, (69, 'rh-105', 'F', 0): 8.2e-11, (70, 'ag-108m', 'S', 0): 3.7e-8,
    (70, 'ag-110m', 'S', 0): 1.2e-8, (71, 'cd-117', 'F', 0): 6.7e-11, (73, 'sb-124', 'M', 0): 6.4e-9,
    (74, 'te-129', 'F', 0): 1.6e-11, (76, 'cs-130', 'S', 0): 1.4e-11, (76, 'ba-133', 'F', 0): 1.5e-9,
    (76, 'ba-133', 'S', 0): 1.0e-8, (78, 'nd-139', 'S', 0): 1.0e-11, (78, 'nd-141', 'M', 0): 4.8e-12,
    (80, 'gd-152', 'M', 0): 8.0e-6, (83, 'hf-175', 'F', 0): 7.2e-10, (87, 'au-198', 'S', 0): 8.6e-10,
    (88, 'au-201', 'F', 0): 8.7e-12, (88, 'hg-197m', 'F', 1): 1.1e-10, (91, 'ra-226', 'F', 0): 3.6e-7,
    (91, 'ra-228', 'F', 0): 9.0e-7, (93, 'np-232', 'M', 0): 5.0e-11, (66, 'zr-93', 'F', 0): 2.5e-8,
    (66, 'zr-93', 'M', 0): 1.0e-8, (84, 'w-181', 'F', 0): 2.7e-11, (94, 'pu-238', 'S', 0): 1.6e-5,
    (95, 'am-241', 'F', 0): 9.6e-5, (96, 'cm-242', 'S', 0): 5.9e-6,
}
# Adult e(50) is normally 0.3-1.0 times the 15-year value; rows outside this range are checked.
RATIO_RANGE: tuple[float, float] = (0.2, 1.5)

_FIX = str.maketrans({'$': '5', 'S': '5', 's': '5', '§': '5', 'O': '0', 'o': '0', 'L': '1', 'l': '1',
                      'I': '1', 'J': '1', '@': '9'})


def parse_coefficient(tok: str) -> float | None:
    """ Parse a typewritten 'd.dE-dd' coefficient with common OCR confusions, else None.

    Every e(t) in Tables A.2 and A.4 is below 1, so 'E+09' is read as 'E-09' (the Ag-110m entry
    of the previous table was 5.5E+09 for this reason).
    """
    t = tok.strip().strip('()._-°|;,:\'"').replace('\u2014', '-').replace('\u2013', '-').replace('+', '-')
    t = t.replace(' ', '')
    m = re.fullmatch(r'([0-9SsOoLlIJ§@$])[.,]+([0-9SsOoLlIJ§@$])(?:[E£Be]|[68](?=[-~]))[E£B8e~\-]*(\d\d)', t)
    if not m:
        return None
    return float(f'{m.group(1).translate(_FIX)}.{m.group(2).translate(_FIX)}e-{m.group(3)}')


def _lines_and_words(ink, word_gap: int = 14) -> list[list[list[tuple]]]:
    """ Text lines (top to bottom) of glyph words (left to right) on a typewritten page """
    import numpy as np
    boxes = glyph_boxes(ink, min_area=6)
    if len(boxes) == 0:
        return []
    yc = (boxes[:, 1] + boxes[:, 3]) / 2.0
    order = np.argsort(yc)
    groups, cur = [], [order[0]]
    for i in order[1:]:
        if yc[i] - np.median(yc[cur]) <= 8:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    out = []
    for g in groups:
        items = [(*map(int, boxes[i]), '') for i in g]
        words = _group_words(items, gap=word_gap)
        out.append([[w[:4] for w in wd] for wd in words])
    return out


def _page_tsv_words(ink, scale: float, psm: str = '6') -> list[tuple[float, float, str]]:
    """ Full-page tesseract reading: (x centre, y centre, text) per word in 300-dpi coordinates """
    import numpy as np
    from PIL import Image
    img = Image.fromarray(np.where(ink, 0, 255).astype(np.uint8))
    if scale != 1.0:
        img = img.resize((int(img.width * scale), int(img.height * scale)), Image.LANCZOS)
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, 'p.png')
        img.save(path)
        env = dict(os.environ, OMP_THREAD_LIMIT='1')
        tsv = subprocess.run(['tesseract', path, '-', '--psm', psm, 'tsv'], capture_output=True, text=True,
                             env=env).stdout
    out = []
    for rec in tsv.splitlines()[1:]:
        f = rec.split('\t')
        if len(f) < 12 or f[0] != '5' or not f[11].strip():
            continue
        x, y, w, h = (int(v) / scale for v in f[6:10])
        out.append((x + w / 2.0, y + h / 2.0, f[11].strip()))
    return out


def read_a2_page(args: tuple[str, int]) -> list[dict]:
    """ Worker: segment one Table A.2 page and read its fields several times.

    Columns are located by position: the Type letter is the one- or two-glyph word at the page's
    Type column, the adult value is the rightmost word and the 15-year value the word before it,
    and the name and half-life are the words left of the Type column. Each geometric word collects
    the text of every full-page tesseract reading (three scales) whose word centre falls inside it,
    plus whitelisted crop readings of the two value columns (three scales).
    """
    import numpy as np
    pdf, page = args
    ink = render_page(Path(pdf), page)
    lines = _lines_and_words(ink)
    wide = [ln for ln in lines if len(ln) >= 8]
    single_x = [min(b[0] for b in w) for ln in wide for w in ln[:-6] if len(w) <= 2]
    right_x = [max(b[2] for b in ln[-1]) for ln in wide]
    if not single_x or not right_x:
        return []
    type_x = float(np.median(single_x))
    adult_x1 = float(np.median(right_x))
    readings = [_page_tsv_words(ink, sc) for sc in (1.0, 1.5, 2.0)]

    def texts_in(boxes, pad=4):
        x0 = min(b[0] for b in boxes) - pad
        x1 = max(b[2] for b in boxes) + pad
        y0 = min(b[1] for b in boxes) - pad
        y1 = max(b[3] for b in boxes) + pad
        return [' '.join(t for cx, cy, t in rd if x0 <= cx <= x1 and y0 <= cy <= y1) for rd in readings]

    recs = []
    for words in lines:
        x0s = [min(b[0] for b in w) for w in words]
        tw = [i for i, w in enumerate(words) if len(w) <= 2 and abs(x0s[i] - type_x) < 30]
        last_x1 = max(b[2] for b in words[-1])
        if tw and len(words) >= 6 and abs(last_x1 - adult_x1) < 60:
            i = tw[0]
            left = [b for w in words[:i] for b in w]
            recs.append({'kind': 'type', 'type_w': words[i], 'e15_w': words[-2], 'adult_w': words[-1],
                         'name_w': left, 'y': (min(b[1] for w in words for b in w),
                                               max(b[3] for w in words for b in w))})
        elif len(words) <= 3 and x0s[0] < type_x - 30:
            recs.append({'kind': 'text', 'text_w': [b for w in words for b in w]})

    def crop_ocr(key, whitelist, scale):
        groups = [r[key] for r in recs if r.get(key)]
        texts = ocr_stack(ink, groups, whitelist, scale=scale, width=1200)
        it = iter(texts)
        return [next(it) if r.get(key) else '' for r in recs]

    crops = {k: [crop_ocr(k, '0123456789.E-+', sc) for sc in (1.0, 1.5)] for k in ('adult_w', 'e15_w')}
    out = []
    for i, r in enumerate(recs):
        if r['kind'] == 'text':
            out.append({'page': page, 'kind': 'text', 'text': texts_in(r['text_w'])[0]})
            continue
        vals = {}
        for k in ('adult_w', 'e15_w'):
            txt = texts_in(r[k], pad=2) + [c[i] for c in crops[k]]
            vals[k] = [parse_coefficient(t) for t in txt]
        types = [t.strip()[:1] for t in texts_in(r['type_w'], pad=2)]
        out.append({'page': page, 'kind': 'type', 'types': types, 'y': r['y'],
                    'names': texts_in(r['name_w'], pad=2) if r['name_w'] else [],
                    'adult': vals['adult_w'], 'e15': vals['e15_w']})
    return out


def _vote(values: list[float | None], e15: float | None) -> tuple[float | None, str]:
    """ Majority of the parsed readings; ties broken by closeness to the 15-year value """
    from collections import Counter
    c = Counter(v for v in values if v is not None).most_common()
    if not c:
        return None, 'check'
    if len(c) == 1 or c[0][1] > c[1][1]:
        return c[0][0], 'agree' if c[0][1] >= 2 else 'single'
    near = [v for v, n in c if n == c[0][1] and e15 and e15 / 5.0 <= v <= 5.0 * e15]
    return (near[0], 'resolved') if len(near) == 1 else (c[0][0], 'check')


def _name_vote(names: list[str]) -> str:
    """ Most common parsed 'Xx-NNN[m] half-life' reading of a name field, else the first reading """
    from collections import Counter
    keyed = []
    for n in names:
        m = re.match(r'^(\S{0,3}?)\s*[-~+]+\s*([0-9ilIO]{1,3})\s*(m\d?)?\s*(.*)$', n.strip())
        if m:
            keyed.append((m.group(2) + (m.group(3) or ''), n.strip()))
    if not keyed:
        return names[0].strip() if names else ''
    best = Counter(k for k, _ in keyed).most_common(1)[0][0]
    return next(n for k, n in keyed if k == best)


def _type_vote(types: list[str]) -> str:
    from collections import Counter
    fixed = [{'s': 'S', '5': 'S', '$': 'S', '8': 'S', '§': 'S', '3': 'S', 'P': 'F', 'E': 'F', '¥': 'F',
              'f': 'F', 'm': 'M'}.get(t, t) for t in types]
    c = Counter(t for t in fixed if t in ('F', 'M', 'S')).most_common(1)
    return c[0][0] if c else '?'


def assemble_a2(page_recs: list[dict]) -> list[dict]:
    """ Attach element headings and nuclide names to the type rows and vote on the values """
    rows = []
    element, nuclide, half_life, last_mass = None, None, '', 0
    for r in page_recs:
        if r['kind'] == 'text':
            hm = re.match(r'^([A-Z][a-z]{2,})', r['text'].strip())
            sym = element_from_heading(hm.group(1)) if hm else None
            if sym:
                element, last_mass = sym, 0
            continue
        r['name'] = _name_vote(r.get('names', []))
        r['type'] = _type_vote(r.get('types', []))
        if r['name'] and not r['name'].startswith(('(', '{', 'Nuclide', 'compounds')):
            if r['name'].startswith('Tritium'):
                nuclide, half_life = 'h-3 compounds', ''
            else:
                nm = re.match(r'^(\S{0,3}?)\s*[-~+=]+\s*([0-9ilIOB]{1,3})\s*(m\d?)?\s*(.*)$', r['name'])
                if nm:
                    mass = int(nm.group(2).translate(str.maketrans('ilIOB', '11108')))
                    lab = element_symbol(nm.group(1)) if re.fullmatch(r'[A-Za-z][a-z]?', nm.group(1)) else None
                    if lab and (element is None or (ELEMENT_Z[lab.lower()] > ELEMENT_Z[element.lower()]
                                                    and mass < last_mass)):
                        element = lab  # heading missed: a later element whose masses restart lower
                    nuclide = f'{(element or lab or "?").lower()}-{mass}{nm.group(3) or ""}'
                    last_mass = mass
                    hl = re.match(r'([\d.,E+]+)\s*([a-z])', nm.group(4))
                    half_life = f'{hl.group(1)} {hl.group(2)}' if hl else ''
                else:
                    nuclide, half_life = '?' + r['name'], ''
        if nuclide is None:
            continue
        if all(v is None for v in r['adult'] + r['e15']):
            continue  # column headings
        typ = r['type']
        e15, _ = _vote(r['e15'], None)
        adult, status = _vote(r['adult'], e15)
        unanimous = adult is not None and all(v == adult for v in r['adult'])
        if status in ('single', 'agree') and e15 and not (RATIO_RANGE[0] <= adult / e15 <= RATIO_RANGE[1]):
            status = 'check'
        elif not unanimous and status == 'single' and not e15:
            status = 'check'  # one reading and nothing to compare it with
        rows.append({'page': r['page'], 'y': r.get('y'), 'nuclide': nuclide, 'half_life': half_life, 'type': typ,
                     'adult': adult, 'e15': e15, 'status': status, 'readings': r['adult']})
    return rows


def resolve_names(rows: list[dict]) -> None:
    """ Apply NAME_FIXES and mark implausible names with a leading '?' """
    for r in rows:
        hl = re.sub(r'([\d.E+]+)\s*([a-z])[a-z]*', r'\1 \2', r['half_life']).strip()
        key = (r['nuclide'], hl)
        if key in NAME_FIXES:
            r['nuclide'] = NAME_FIXES[key]
        elif r['nuclide'] != 'h-3 compounds' and not is_plausible_nuclide(r['nuclide']):
            r['nuclide'] = '?' + r['nuclide']
    seen: dict[tuple, int] = {}
    for r in rows:
        k = (r['page'], r['nuclide'], r['type'])
        seen[k] = seen.get(k, 0) + 1
        fix = VALUE_FIXES.get(k + (seen[k] - 1,))
        if fix is not None:
            r['adult'], r['status'] = fix, 'manual'


def _ocr_page(args: tuple[str, int]) -> tuple[int, str]:
    pdf, page = args
    with tempfile.TemporaryDirectory() as td:
        prefix = os.path.join(td, 'p')
        subprocess.run(['pdftoppm', '-f', str(page), '-l', str(page), '-r', '300', '-gray', '-png',
                        '-singlefile', pdf, prefix], check=True)
        env = dict(os.environ, OMP_THREAD_LIMIT='1')
        txt = subprocess.run(['tesseract', prefix + '.png', '-', '--psm', '6'], capture_output=True,
                             text=True, env=env).stdout
    return page, txt


def parse_table_a4(txt: str) -> dict[str, float]:
    """ Parse Table A.4 (inert gases): 'Kr-85 10.7 y 2.2E-11' rows -> {nuclide: Sv/day per Bq/m^3} """
    out: dict[str, float] = {}
    for line in txt.splitlines():
        line = re.sub(r'(\d[.,])\s+(\d)', r'\1\2', line)  # '5. 9E-10'
        m = re.match(r'^\s*([A-Z][a-z]?)\s*[-~]+\s*(\d{1,3})\s*(m)?\s.*?(\S+)\s*$', line)
        if not m:
            continue
        nuc = f'{m.group(1).lower()}-{int(m.group(2))}{m.group(3) or ""}'
        val = parse_coefficient(m.group(4))
        if val is not None and is_plausible_nuclide(nuc):
            out[nuc] = val
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Extract ICRP-72 adult inhalation and immersion DCFs.")
    parser.add_argument("--pdf", default=str(DEFAULT_PDF), help="ICRP-72 compilation (ANIB_26_1.pdf)")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--audit", default=None, help="Write per-row audit CSV to this directory")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args(argv)
    for tool in ('pdftoppm', 'tesseract'):
        if shutil.which(tool) is None:
            raise RuntimeError(f"Required tool not found in PATH: {tool}")
    pdf = str(Path(args.pdf))
    if not Path(pdf).is_file():
        print(f"PDF not found: {pdf}", file=sys.stderr)
        return 1
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as ex:
        page_recs = list(ex.map(read_a2_page, [(pdf, p) for p in TABLE_A2_PAGES]))
        a4_text = dict(ex.map(_ocr_page, [(pdf, TABLE_A4_PAGE)]))[TABLE_A4_PAGE]
    merged = assemble_a2([r for recs in page_recs for r in recs])
    resolve_names(merged)
    a4 = parse_table_a4(a4_text)

    if args.audit:
        import csv
        Path(args.audit).mkdir(parents=True, exist_ok=True)
        with open(Path(args.audit) / 'icrp72_table_a2_rows.csv', 'w', newline='') as f:
            wr = csv.writer(f)
            wr.writerow(['page', 'nuclide', 'half_life', 'type', 'adult', 'e15', 'status', 'readings'])
            for r in merged:
                wr.writerow([r['page'], r['nuclide'], r['half_life'], r['type'], r['adult'], r['e15'],
                             r['status'], ' '.join(str(v) for v in r['readings'])])

    bad = [r for r in merged if r['status'] == 'check' or r['adult'] is None or r['type'] not in ('F', 'M', 'S')
           or r['nuclide'].startswith('?')]
    for r in bad:
        print(f"  check: p.{r['page']} {r['nuclide']} ({r['half_life']}) type {r['type']}: "
              f"{r['readings']} (15 y: {r['e15']})")
    good = [r for r in merged if r not in bad and r['nuclide'] != 'h-3 compounds']
    per_type: dict[str, dict[str, float]] = {}
    for r in good:
        t = per_type.setdefault(r['nuclide'], {})
        t[r['type']] = max(t.get(r['type'], 0.0), r['adult'])  # Hg: larger of the organic/inorganic rows
    fastest = {n: next(t[k] for k in ('F', 'M', 'S') if k in t) for n, t in per_type.items()}
    type_m = {n: t['M'] if 'M' in t else fastest[n] for n, t in per_type.items()}
    max_type = {n: max(t.values()) for n, t in per_type.items()}
    for table in (fastest, type_m, max_type):
        table['h-3'] = HTO_ADULT_SV_BQ
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / 'dcf_icrp72_inhalation_adult.csv', fastest, 'nuclide,dcf_sv_bq')
    write_csv(out_dir / 'dcf_icrp72_anib26_1_adult.csv', type_m, 'nuclide,dcf_sv_bq')
    write_csv(out_dir / 'dcf_icrp72_inhalation_adult_max.csv', max_type, 'nuclide,dcf_sv_bq')
    write_csv(out_dir / 'dcf_icrp72_immersion_adult.csv', a4, 'nuclide,dcf_sv_per_bq_m3_day')
    if len(a4) != N_TABLE_A4:
        print(f"  check: Table A.4 gave {len(a4)} nuclides, expected {N_TABLE_A4}")
        bad.append({'table': 'A.4'})
    print(f"Table A.2: {len(merged)} rows, {len(bad)} unresolved, {len(max_type)} nuclides; "
          f"Table A.4: {len(a4)} nuclides")
    return 0 if not bad else 2


if __name__ == "__main__":
    raise SystemExit(main())
