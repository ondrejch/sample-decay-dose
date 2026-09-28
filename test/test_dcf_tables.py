"""Regression checks for the dose-coefficient tables in leaky_box_origen/data/ and their extractors.

Anchor values are the published numbers read from the source pages:
  - FGR-11 (EPA-520/1-88-020): Table 2.1 h_E,50 [Sv/Bq] and Table 2.3 h_E,ext [Sv/hr per Bq/m^3]
    (OSTI copy 6294233.pdf, PDF pp. 131-162 and 191).
  - ICRP-72 (Annals of the ICRP 26(1)): Table A.2 adult e(50) [Sv/Bq], Table A.3 HTO, and
    Table A.4 immersion [Sv/day per Bq/m^3] (ANIB_26_1.pdf, printed pp. 44-89).
A bad extraction (exponent misread, row shift, wrong absorption type) moves an anchor by at least
a factor of 1.5, so each anchor must match within that factor.
"""
import csv
import math
from pathlib import Path

import numpy as np
import pytest

from leaky_box_origen import extract_fgr11_dcf as fgr
from leaky_box_origen import extract_icrp72_dcf as icrp
from leaky_box_origen.LeakyBox import LEAKY_BOX_DATA_DIR, MAX_DCF, MIN_DCF, _load_dcf_csv, _load_dcf_immersion_csv

DATA = Path(LEAKY_BOX_DATA_DIR)
HR_TO_DAY = 24.0

# table file -> (header, {nuclide: published value in the CSV units})
ANCHORS = {
    'dcf_fgr11_inhalation_worker.csv': ('nuclide,dcf_sv_bq', {
        'h-3': 1.73e-11,     # water vapour (V)
        'co-60': 5.91e-8,    # class Y (largest listed class)
        'sr-90': 3.51e-7,    # class Y
        'cs-137': 8.63e-9,   # class D
        'i-131': 8.89e-9,    # class D
        'pu-238': 1.06e-4,   # class W
    }),
    'dcf_fgr11_submersion_worker.csv': ('nuclide,dcf_sv_per_bq_m3_day', {
        'h-3': 1.19e-15 * HR_TO_DAY,   # elemental
        'ar-41': 2.17e-10 * HR_TO_DAY,
        'kr-85': 4.70e-13 * HR_TO_DAY,
        'kr-87': 1.42e-10 * HR_TO_DAY,
        'xe-133': 6.07e-12 * HR_TO_DAY,
        'xe-138': 1.92e-10 * HR_TO_DAY,
    }),
    'dcf_icrp72_inhalation_adult.csv': ('nuclide,dcf_sv_bq', {   # Type F where listed
        'h-3': 1.8e-11, 'co-60': 5.2e-9, 'sr-90': 2.4e-8, 'cs-137': 4.6e-9, 'i-131': 7.4e-9,
        'pu-238': 1.1e-4, 'ag-110m': 5.5e-9,
    }),
    'dcf_icrp72_anib26_1_adult.csv': ('nuclide,dcf_sv_bq', {     # Type M where listed
        'h-3': 1.8e-11, 'co-60': 1.0e-8, 'sr-90': 3.6e-8, 'cs-137': 9.7e-9, 'i-131': 2.4e-9,
        'pu-238': 4.6e-5, 'ag-110m': 7.6e-9,
    }),
    'dcf_icrp72_inhalation_adult_max.csv': ('nuclide,dcf_sv_bq', {  # largest over F/M/S
        'h-3': 1.8e-11, 'co-60': 3.1e-8, 'sr-90': 1.6e-7, 'cs-137': 3.9e-8, 'i-131': 7.4e-9,
        'pu-238': 1.1e-4, 'ag-110m': 1.2e-8,
    }),
    'dcf_icrp72_immersion_adult.csv': ('nuclide,dcf_sv_per_bq_m3_day', {
        'ar-41': 5.3e-9, 'kr-85': 2.2e-11, 'kr-87': 3.4e-9, 'xe-133': 1.2e-10, 'xe-138': 4.7e-9,
    }),
}


def _raw_rows(name):
    with open(DATA / name, newline='') as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def _load(name):
    if 'dcf_sv_per_bq_m3_day' in ANCHORS[name][0]:
        return _load_dcf_immersion_csv(str(DATA / name))
    return _load_dcf_csv(str(DATA / name))


@pytest.mark.parametrize('name', sorted(ANCHORS))
def test_header_names_and_value_range(name):
    header, rows = _raw_rows(name)
    assert ','.join(header) == ANCHORS[name][0]
    names = [r[0] for r in rows]
    assert len(names) == len(set(names)), 'duplicate nuclide rows'
    bad = [n for n in names if not fgr.is_plausible_nuclide(n)]
    assert not bad, f'impossible nuclide names: {bad}'
    values = [float(r[1]) for r in rows]
    assert all(MIN_DCF <= v <= MAX_DCF for v in values)
    # The loader keeps every row, so no nuclide silently loses its coefficient.
    assert len(_load(name)) == len(rows)


@pytest.mark.parametrize('name', sorted(ANCHORS))
def test_anchor_values_within_factor_1p5(name):
    table = _load(name)
    for nuclide, published in ANCHORS[name][1].items():
        assert nuclide in table, f'{nuclide} missing from {name}'
        ratio = table[nuclide] / published
        assert 1 / 1.5 <= ratio <= 1.5, f'{name} {nuclide}: {table[nuclide]:.3g} vs published {published:.3g}'


def test_icrp72_type_tables_are_ordered_per_nuclide():
    # Type F / Type M selections can never exceed the maximum over types.
    fastest = _load('dcf_icrp72_inhalation_adult.csv')
    type_m = _load('dcf_icrp72_anib26_1_adult.csv')
    largest = _load('dcf_icrp72_inhalation_adult_max.csv')
    assert set(fastest) == set(type_m) == set(largest)
    for n in largest:
        assert fastest[n] <= largest[n] * (1 + 1e-9)
        assert type_m[n] <= largest[n] * (1 + 1e-9)


def test_fgr11_and_icrp72_agree_within_lung_model_spread():
    # FGR-11 (ICRP-30 lung model, largest class) and ICRP-72 (ICRP-66 model, largest type) differ by
    # up to ~20x for genuine reasons. A transcription error (exponent or row shift) exceeds 30x.
    f11 = _load('dcf_fgr11_inhalation_worker.csv')
    i72 = _load('dcf_icrp72_inhalation_adult_max.csv')
    shared = [n for n in f11 if n in i72]
    assert len(shared) > 700
    outliers = [n for n in shared if max(f11[n] / i72[n], i72[n] / f11[n]) > 30.0]
    assert not outliers, outliers
    sub = _load('dcf_fgr11_submersion_worker.csv')
    imm = _load('dcf_icrp72_immersion_adult.csv')
    assert set(imm) <= set(sub)
    # Skin dose is outside FGR-11 h_E and inside the ICRP-72 immersion value: Ar-39 differs ~8x.
    assert all(0.1 <= sub[n] / imm[n] <= 3.0 for n in imm)


def test_h3_rows_present():
    for name in ('dcf_fgr11_inhalation_worker.csv', 'dcf_fgr11_submersion_worker.csv',
                 'dcf_icrp72_inhalation_adult.csv', 'dcf_icrp72_anib26_1_adult.csv',
                 'dcf_icrp72_inhalation_adult_max.csv'):
        assert 'h-3' in _load(name)


# ------------------------------------------------------------------------------------------------
# FGR-11 extractor helpers
# ------------------------------------------------------------------------------------------------

def test_element_lookup_and_headings():
    assert fgr.element_symbol('Hydrogen') == 'H'
    assert fgr.element_symbol('cs') == 'Cs'
    assert fgr.element_symbol('Xx') is None
    assert fgr.element_from_heading('lodine') == 'I'       # OCR 'l' for 'I'
    assert fgr.element_from_heading('Sulphur') == 'S'
    assert fgr.element_from_heading('Samarlum') == 'Sm'    # one wrong letter
    assert fgr.element_from_heading('Anericiun') == 'Am'   # two wrong letters, long name
    assert fgr.element_from_heading('Labelled') is None
    assert fgr.element_from_heading('Table') is None


@pytest.mark.parametrize('label, element, expected', [
    ('Be-7', None, 'be-7'),
    ('H-3', 'H', 'h-3'),                # one-digit mass numbers
    ('Kr-8 1', 'Kr', 'kr-81'),
    ('C1-38', 'Cl', 'cl-38'),
    ('Sm-14im', 'Sm', 'sm-141m'),
    ('La-14]', 'La', 'la-141'),
    ('Xe- 129m', 'Xe', 'xe-129m'),
    ('1-131', 'I', 'i-131'),
    ('Be-IO', 'Be', 'be-10'),
    ('Nb-89', 'Nb', 'nb-89'),
    ('122 m', 'Nb', None),
    ('Te-310', 'Te', None),             # implausible mass number
])
def test_parse_nuclide_label(label, element, expected):
    assert fgr.parse_nuclide_label(label, element) == expected


def test_plausible_nuclide_rejects_old_junk_names():
    for junk in ('i-00', 'i-78m', 'i-94', 'to-104', 'to-94m', 'ec-161', 'e-10', 'y-310', 'd-110', 'i-01'):
        assert not fgr.is_plausible_nuclide(junk), junk
    for ok in ('h-3', 'tc-104', 'tc-94m', 'ag-110m', 'tb-156m2', 'pu-238'):
        assert fgr.is_plausible_nuclide(ok), ok


def test_half_life_labels():
    assert fgr.half_life_label('122 m') == '122 m'
    assert fgr.half_life_label('49h') == '4.9 h'
    assert fgr.half_life_label('5,76 d') == '5.76 d'
    assert fgr.half_life_label('Nb-89') is None
    assert fgr.HALF_LIFE_ISOMERS[('sb-120', '5.76 d')] == 'sb-120m'


def test_text_layer_mantissas():
    line = 'Be-7   W 5 IO\'   3.72 IQ"   3.12 lu"  2.15 IU\'O  6.37tc"'
    assert fgr.text_layer_mantissas(line) == ['3.72', '3.12', '2.15', '6.37']
    assert fgr.text_layer_mantissas('C-11  1.0t  3,41 lo-"') == ['3.41']


BE7_W = [3.72e-11, 3.12e-11, 2.15e-10, 4.58e-11, 4.09e-11, 2.60e-11, 5.46e-11, 6.37e-11]  # p. 131


def test_effective_is_weighted_organ_sum():
    assert fgr.effective_from_organs(BE7_W[:7]) == pytest.approx(6.366e-11, rel=1e-3)
    assert fgr.row_is_consistent(BE7_W)
    wrong = list(BE7_W)
    wrong[7] = 6.37e-10  # exponent misread in the effective column
    assert not fgr.row_is_consistent(wrong)
    # Table 2.3 rows leave organ cells blank (H-3: lung only, h_E = 0.12 h_lung)
    assert fgr.row_is_consistent([None, None, 9.90e-15, None, None, None, None, 1.19e-15])


def test_solve_row_repairs_single_exponent_misread():
    options = []
    for i, v in enumerate(BE7_W):
        m = f'{v:.2e}'.split('e')[0]
        e = -int(math.floor(math.log10(v)))
        votes_e = {e: 2.0}
        if i == 2:
            votes_e = {e + 4: 1.2, e: 1.0}  # lung exponent misread by most readers ('-10' as '-14')
        options.append(fgr.cell_options({m: 2.0}, votes_e))
    first = [opt[0][0] for opt in options]
    assert not fgr.row_is_consistent(first)
    values, cost, ok = fgr.solve_row(options)
    assert ok and cost > 0
    assert values[2] == pytest.approx(2.15e-10)
    assert values[7] == pytest.approx(6.37e-11)


def test_cell_options_costs():
    opts = fgr.cell_options({'1.73': 2.0, '1.78': 0.5}, {11: 2.5})
    assert opts[0] == (pytest.approx(1.73e-11), 0.5)
    assert any(v == pytest.approx(1.73e-10) for v, _ in opts)  # widened neighbour exponent
    assert fgr.cell_options({}, {11: 1.0}) == []


def test_max_by_nuclide_and_bounds():
    rows = [('co-60', 8.94e-9), ('co-60', 5.91e-8), ('zn-69', 1.06e-97), ('ar-39', 0.013)]
    assert fgr.max_by_nuclide(rows) == {'co-60': 5.91e-8}


def test_write_csv_roundtrip(tmp_path):
    out = tmp_path / 'x.csv'
    fgr.write_csv(out, {'kr-85': 4.70e-13 * 24, 'h-3': 1.73e-11}, fgr.INHALATION_HEADER)
    lines = out.read_text().splitlines()
    assert lines[0] == 'nuclide,dcf_sv_bq'
    assert lines[1:] == ['h-3,1.73e-11', 'kr-85,1.128e-11']


def test_manual_corrections_are_self_consistent():
    for key, (values, page) in fgr.MANUAL_CORRECTIONS.items():
        assert key[0] in ('2.1', '2.3') and fgr.is_plausible_nuclide(key[1])
        assert fgr.EDITIONS['osti'].inhalation_pages[0] <= page <= fgr.EDITIONS['osti'].inhalation_pages[-1]
        assert fgr.row_is_consistent(list(values)), key


def test_trim_minus_cuts_touching_exponent_minus():
    ink = np.zeros((17, 20), dtype=bool)
    ink[7:9, 0:7] = True     # minus bar in the middle band
    ink[0:17, 7:9] = True    # left stroke of the digit
    ink[0:2, 7:18] = True
    ink[15:17, 7:18] = True
    ink[0:17, 16:18] = True  # a '0'-like box touching the bar
    trimmed = fgr._trim_minus(ink, (0, 0, 20, 17))
    assert trimmed[0] == 7
    lone = np.zeros((17, 10), dtype=bool)
    lone[:, 4:6] = True      # a lone '1' is narrower than 0.8 h and stays unchanged
    assert fgr._trim_minus(lone, (0, 0, 10, 17)) == (0, 0, 10, 17)


def test_submersion_conversion_is_per_hour_to_per_day():
    # FGR-11 Table 2.3 is tabulated in Sv/hr per Bq/m^3.
    assert fgr.SUBMERSION_HR_TO_DAY == 24.0


# ------------------------------------------------------------------------------------------------
# ICRP-72 extractor helpers
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('token, expected', [
    ('2.6E-11', 2.6e-11),
    ('5.5E+09', 5.5e-9),     # sign misread, the old Ag-110m entry
    ('1.6£-09', 1.6e-9),
    ('3.68-08', 3.6e-8),     # 'E' read as '8'
    ('4,1E~08', 4.1e-8),
    ('§.6E-11', 5.6e-11),
    ('1.000', None),         # f1 value
    ('Adult', None),
])
def test_parse_coefficient(token, expected):
    got = icrp.parse_coefficient(token)
    if expected is None:
        assert got is None
    else:
        assert got == pytest.approx(expected)


def test_parse_table_a4_lines():
    txt = 'Krypton\nKr~76 14.8 h 1.6£-09\nKr-83m 1.83 h 2.1£-13\nKr-85m 4.48 h 5. 9E~10\nXe-135m 15.3 min 1. 6E-09.\n'
    assert icrp.parse_table_a4(txt) == {'kr-76': pytest.approx(1.6e-9), 'kr-83m': pytest.approx(2.1e-13),
                                        'kr-85m': pytest.approx(5.9e-10), 'xe-135m': pytest.approx(1.6e-9)}


def test_vote_majority_and_tiebreak():
    assert icrp._vote([2.1e-10, 2.1e-10, 2.7e-10], None) == (2.1e-10, 'agree')
    assert icrp._vote([8.2e-11, 6.2e-11], 9.6e-11)[1] == 'check'   # both within 5x of the 15 y value
    assert icrp._vote([1.0e-8, 1.0e-6], 1.1e-8) == (1.0e-8, 'resolved')
    assert icrp._vote([None, None], 1e-9) == (None, 'check')


def test_icrp_fix_tables_are_plausible():
    for (page, nuclide, typ, occ), value in icrp.VALUE_FIXES.items():
        assert page in icrp.TABLE_A2_PAGES and typ in ('F', 'M', 'S') and occ >= 0
        assert fgr.is_plausible_nuclide(nuclide)
        assert 1e-15 <= value <= 1e-2
    for (_, _), name in icrp.NAME_FIXES.items():
        assert fgr.is_plausible_nuclide(name)
    assert icrp.HTO_ADULT_SV_BQ == pytest.approx(1.8e-11)
