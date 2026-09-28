""" Smoke tests: every estimator class must generate a structurally valid MAVRIC deck """
import os
import re
import tempfile
import unittest
from unittest.mock import patch, mock_open

import sample_decay_dose.SampleDose as sd
import sample_decay_dose.HotCell as hc
import sample_decay_dose.Radiator as rad

BODY_KEYWORDS = ('cylinder', 'ycylinder', 'xcylinder', 'cuboid', 'sphere')


class DummyOrigen:
    """ Minimal stand-in for a run Origen object """

    def __init__(self):
        self.debug = 0
        self.cwd = '.'
        self.case_dir = 'run_DUMMY'
        self.sample_weight = 10.0
        self.sample_density = 7.8
        self.sample_volume = 10.0 / 7.8
        self.SAMPLE_F71_file_name = 'decay.f71'
        self.SAMPLE_F71_position = 12
        self.SAMPLE_DECAY_days = 30.0
        self.decayed_atom_dens = {'fe-56': 0.08}
        self.ORIGEN_input_file_name = 'origen.inp'

    def get_beta_to_gamma(self):
        return 0.5

    def get_neutron_integral(self):
        return 1.0e5


def active_lines(text: str) -> list[str]:
    """ Deck lines without SCALE comments (lines starting with an apostrophe) """
    return [ln for ln in text.splitlines() if not ln.strip().startswith("'")]


def block_ids(deck: str, keyword: str) -> list[str]:
    """ IDs of 'keyword N' block openings, e.g. pointDetector 1, adjointSource 5, src 2, distribution 1 """
    return [m.group(1) for ln in active_lines(deck) for m in [re.match(rf'\s*{keyword}\s+(\d+)\s*$', ln)] if m]


def plane_values(deck: str, keyword: str) -> list[float]:
    """ Numbers of all 'keyword ... end' plane lists in the gridGeometry block """
    grid = deck[deck.index('gridGeometry 1'):deck.index('end gridGeometry')]
    values = []
    for m in re.finditer(rf'^\s*{keyword}\s+(.*?)\bend\b', grid, re.M | re.S):
        values += [float(x) for x in m.group(1).split()]
    return values


class DeckTestMixin:
    def check_common(self, deck):
        self.assertIn('=mavric parm=(   )', deck)
        self.assertIn('read parameters', deck)
        self.assertIn('read geometry', deck)
        self.assertIn('read definitions', deck)
        self.assertIn('read importanceMap', deck)
        self.assertIn('read tallies', deck)
        self.assertTrue(deck.rstrip().endswith('end'))
        detectors = block_ids(deck, 'pointDetector')
        self.assertIn('1', detectors)
        self.assertIn('2', detectors)
        adjoint_sources = block_ids(deck, 'adjointSource')
        if 'adjointFluxes=' in '\n'.join(active_lines(deck)):
            self.assertEqual(adjoint_sources, [])  # MAVRIC rejects adjoint sources with an adjoint flux file
        else:
            self.assertEqual(sorted(adjoint_sources), sorted(detectors))  # one adjoint source per detector
        distributions = set(block_ids(deck, 'distribution'))
        for src_id in block_ids(deck, 'src'):
            self.assertIn(src_id, distributions)  # src N uses eDistributionID=N
        self.check_bodies(deck)

    def check_bodies(self, deck):
        """ Every body referenced by media and boundary records is defined in the same unit """
        geometry = deck[deck.index('read geometry') + len('read geometry'):deck.index('end geometry')]
        units = re.split(r'^\s*(?:global\s+)?unit\s+\d+\s*$', geometry, flags=re.M)[1:]
        self.assertTrue(units)
        for unit in units:
            defined, referenced = set(), set()
            for ln in active_lines(unit):
                tokens = ln.split()
                if not tokens:
                    continue
                if tokens[0] in BODY_KEYWORDS:
                    defined.add(int(tokens[1]))
                elif tokens[0] == 'media':
                    referenced |= {abs(int(t)) for t in tokens[3:]}
                elif tokens[0] == 'boundary':
                    referenced.add(int(tokens[1]))
            self.assertTrue(referenced <= defined, f'undefined bodies {referenced - defined} in unit:\n{unit}')


class TestSampleDoseDecks(unittest.TestCase, DeckTestMixin):

    def setUp(self):
        self.o = DummyOrigen()

    def test_dose_estimator(self):
        deck = sd.DoseEstimator(self.o).mavric_deck()
        self.check_common(deck)
        self.assertIn('Sample dose', deck)
        self.assertNotIn('helium', deck)
        self.assertNotIn('adjointSource 5', deck)

    def test_square_tank(self):
        deck = sd.DoseEstimatorSquareTank(self.o).mavric_deck()
        self.check_common(deck)
        self.assertIn('DoseEstimatorSquareTank', deck)
        self.assertIn('media 10', deck)
        self.assertNotIn('helium 2 end', deck)

    def test_storage_tank(self):
        deck = sd.DoseEstimatorStorageTank(self.o).mavric_deck()
        self.check_common(deck)
        self.assertIn('DoseEstimatorStorageTank', deck)
        self.assertIn('helium 2 end', deck)
        self.assertIn('media 10 1 -1 -2 3', deck)  # gas plenum

    def test_generic_tank(self):
        est = sd.DoseEstimatorGenericTank(self.o)
        est.cyl_r = 5.0
        deck = est.mavric_deck()
        self.check_common(deck)
        self.assertIn("DoseEstimatorGenericTank", deck)
        self.assertIn("' helium 2 end", deck)
        self.assertNotIn('adjointSource 5', deck)

    def test_handling_contact_tank(self):
        est = sd.HandlingContactDoseEstimatorGenericTank(self.o)
        est.cyl_r = 5.0
        deck = est.mavric_deck()
        self.check_common(deck)
        self.assertIn('location 2', deck)
        self.assertIn('pointDetector 5', deck)
        self.assertIn('pointDetector 6', deck)
        self.assertIn('adjointSource 6', deck)

    def test_mha_tank(self):
        est = sd.MHATank(self.o)
        est.source_multiplier = 2.0
        deck = est.mavric_deck()
        self.check_common(deck)
        self.assertIn('multiplier=2.0', deck)
        self.assertIn('THIS IS BASICALLY MEANINGLESS', deck)

    def test_no_neutron_source_omitted(self):
        o = DummyOrigen()
        o.get_neutron_integral = lambda: -1.0
        deck = sd.DoseEstimator(o).mavric_deck()
        self.check_common(deck)
        self.assertNotIn('src 1', deck)
        self.assertIn('src 2', deck)
        self.assertNotIn('distribution 1', deck)
        self.assertIn('distribution 2', deck)

    def test_neutron_source_included_when_intensity_unknown(self):
        o = DummyOrigen()
        o.decayed_atom_dens = {}  # ORIGEN not run in this process and no OPUS files on disk

        def missing():
            raise FileNotFoundError('no OPUS file')
        o.get_neutron_integral = missing
        o.get_beta_to_gamma = missing
        est = sd.DoseEstimatorSquareTank(o)
        est.debug = 0
        self.assertIsNone(est.neutron_intensity)
        deck = est.mavric_deck()
        self.check_common(deck)
        self.assertIn('src 1', deck)
        self.assertIn('parameters 12 1 end', deck)

    def test_zero_day_sample_uses_position_1(self):
        o = DummyOrigen()
        o.SAMPLE_DECAY_days = 0.0
        deck = sd.DoseEstimator(o).mavric_deck()
        self.assertIn('parameters 1 5 end', deck)
        self.assertIn('parameters 1 1 end', deck)

    def test_base_detector_outside_sample(self):
        o = DummyOrigen()
        o.sample_volume = 5.0e5  # radius 43 cm, the default detector at 30 cm from the centre is inside
        est = sd.DoseEstimator(o)
        with self.assertRaises(ValueError):
            est.mavric_deck()
        est.det_x = 60.0
        deck = est.mavric_deck()
        self.check_common(deck)
        self.assertIn('position 60.0 0 0', deck)

    def test_mavric_deck_twice_is_identical(self):
        for cls in (sd.DoseEstimatorSquareTank, sd.DoseEstimatorStorageTank, sd.DoseEstimatorGenericTank,
                    sd.HandlingContactDoseEstimatorGenericTank, sd.MHATank, hc.HotCellDoses):
            est = cls(self.o)
            est.debug = 0
            if cls in (sd.DoseEstimatorGenericTank, sd.HandlingContactDoseEstimatorGenericTank):
                est.cyl_r = 5.0
            if cls is hc.HotCellDoses:
                est.layers_thicknesses = [35.0, 7.5]
                est.layers_mats = [sd.ADENS_DRYAIR_COLD, sd.ADENS_LEAD_COLD]
                est.layers_temperature_K = [300.0, 300.0]
            box_a = est.box_a
            self.assertEqual(est.mavric_deck(), est.mavric_deck(), cls.__name__)
            self.assertEqual(est.box_a, box_a, cls.__name__)

    @patch('sample_decay_dose.SampleDose.run_scale', return_value=True)
    @patch('sample_decay_dose.SampleDose.shutil.copy2')
    @patch('builtins.open', new_callable=mock_open)
    @patch('os.chdir')
    @patch('os.mkdir')
    @patch('os.path.exists', return_value=True)
    @patch('os.path.isfile', return_value=True)
    def test_run_mavric_twice_keeps_case_dir(self, *_mocks):
        est = sd.DoseEstimatorSquareTank(self.o)
        est.debug = 0
        est.run_mavric()
        first = est.case_dir
        self.assertEqual(first, 'run_DUMMY_MAVRIC_2.540_5.080_7.620')
        est.run_mavric()
        self.assertEqual(est.case_dir, first)
        est.layers_thicknesses = [1.0, 2.0, 3.0]  # a new layer set replaces the suffix
        est.run_mavric()
        self.assertEqual(est.case_dir, 'run_DUMMY_MAVRIC_1.000_2.000_3.000')

    def test_empty_layers_give_bare_sample(self):
        for cls in (sd.DoseEstimatorSquareTank, sd.DoseEstimatorStorageTank, sd.DoseEstimatorGenericTank,
                    sd.HandlingContactDoseEstimatorGenericTank, sd.MHATank, hc.HotCellDoses):
            est = cls(self.o)
            est.debug = 0
            est.layers_mats, est.layers_thicknesses, est.layers_temperature_K = [], [], []
            if cls in (sd.DoseEstimatorGenericTank, sd.HandlingContactDoseEstimatorGenericTank):
                est.cyl_r = 0.5
            deck = est.mavric_deck()
            self.check_common(deck)  # the outer void references only defined bodies
            self.assertTrue(est.beta_applies, cls.__name__)
            self.assertGreater(est.det_x, est.cyl_r, cls.__name__)  # detector outside the sample

    def test_layer_list_lengths_are_checked_before_the_deck(self):
        est = sd.DoseEstimatorSquareTank(self.o)
        est.layers_thicknesses = [1.0, 2.0, 3.0, 4.0]  # four thicknesses, three materials
        with self.assertRaises(ValueError):
            est.mavric_deck()
        est.layers_thicknesses = [1.0, 2.0, 3.0]
        est.layers_temperature_K = [300.0]
        with self.assertRaises(ValueError):
            est.mavric_deck()

    def _radiator(self, index_record: dict | None = None):
        class DummyOrigenTriton(DummyOrigen):
            def __init__(self):
                super().__init__()
                self.BURNED_MATERIAL_F71_file_name = 'data/burned.f71'
                self.BURNED_MATERIAL_F71_position = 16
                self.BURNED_MATERIAL_F71_index = {
                    15: {'case': '1', 'time': '100.0', 'DCGNAB': 'DC----'},
                    16: index_record or {'case': '2', 'time': '864000.0', 'DCGNAB': 'DCGN--'}}
                self.burned_atom_dens = {'fe-56': 0.08}

        box = rad.RadiatorBox(DummyOrigenTriton())
        box.debug = 0
        box.layers_mats = [sd.ADENS_SS316H_HOT, sd.ADENS_DRYAIR_COLD]
        box.layers_temperature_K = [600.0, 300.0]
        return box

    @patch('sample_decay_dose.Radiator.get_f71_volume', return_value=1.0e6)
    def test_radiator_box(self, mock_volume):
        box = self._radiator()
        deck = box.mavric_deck()
        self.check_common(deck)
        self.assertIn('array', deck)
        # The F71 position of the OrigenFromTriton object feeds both distributions, also without run_mavric()
        self.assertIn('parameters 16 1 end', deck)
        self.assertIn('parameters 16 5 end', deck)
        self.assertIn('src 1', deck)  # the neutron source is present
        self.assertIn('filename="burned.f71"', deck)  # the copied basename, not the source path
        self.assertIn('cp -r ${INPDIR}/burned.f71 .', deck)
        mock_volume.assert_called_with(os.path.join(box.cwd, 'data/burned.f71'), 16)
        self.assertEqual(deck.count(f'multiplier={box.all_pins_volume / 1.0e6}'), 2)
        self.assertIn('RadiatorBox, F71 position 16 (t = 10 d)', deck)
        self.assertNotIn('HotCellDoses', deck)
        self.assertNotIn('30.0 days', deck)
        # The importance-map grid covers both detectors and the problem boundary
        x_planes = plane_values(deck, 'xPlanes')
        self.assertGreaterEqual(max(x_planes), box.handling_det_x + box.planes_xy_around_det)
        self.assertLessEqual(min(x_planes), -(box.radiator_geometry['n_X'] * box.half_pitch + box.box_a))

    @patch('sample_decay_dose.Radiator.get_f71_volume', return_value=1.0e6)
    def test_radiator_needs_f71_spectra(self, _mock_volume):
        box = self._radiator({'case': '1', 'time': '100.0', 'DCGNAB': 'DC----'})  # TRITON-like position
        with self.assertRaises(ValueError):
            box.mavric_deck()
        box = self._radiator({'case': '2', 'time': '100.0', 'DCGNAB': 'DCG---'})  # gamma only
        with self.assertRaises(ValueError):
            box.mavric_deck()
        box.neutron_intensity = 0.0  # explicit photon-only run
        deck = box.mavric_deck()
        self.assertNotIn('src 1', deck)
        self.assertIn('src 2', deck)

    def test_radiator_source_multiplier_override_and_layers(self):
        box = self._radiator()
        box.source_multiplier = 1.0
        deck = box.mavric_deck()  # no obiwan call
        self.assertEqual(deck.count('multiplier=1.0\n'), 2)
        box.layers_mats = [sd.ADENS_SS316H_HOT]
        with self.assertRaises(ValueError):
            box.mavric_deck()

    def test_hot_cell(self):
        class DummyOrigenTriton(DummyOrigen):
            def __init__(self):
                super().__init__()
                self.BURNED_MATERIAL_F71_file_name = 'burned.f71'
                self.BURNED_MATERIAL_F71_position = 16
                self.burned_atom_dens = {'fe-56': 0.08}

        cell = hc.HotCellDoses(DummyOrigenTriton())
        deck = cell.mavric_deck()
        self.check_common(deck)
        self.assertIn('HotCellDoses', deck)
        self.assertIn('cuboid 2', deck)  # square-cuboid shielding layers
        self.assertIn("' helium 2 end", deck)
        self.assertIn('pointDetector 5', deck)
        self.assertIn('pointDetector 6', deck)
        self.assertFalse(cell.beta_applies)

    def test_hot_cell_first_layer_is_half_width(self):
        cell = hc.HotCellDoses(self.o)
        cell.layers_mats = [sd.ADENS_DRYAIR_COLD, sd.ADENS_LEAD_COLD]
        cell.layers_thicknesses = [35.0, 7.5]  # 70 cm cell interior, 7.5 cm of lead
        cell.layers_temperature_K = [300.0, 300.0]
        deck = cell.mavric_deck()
        self.assertIn('cuboid 2 4p 35.0 2p 35.0', deck)
        self.assertIn('cuboid 3 4p 42.5 2p 42.5', deck)
        self.assertAlmostEqual(cell.det_x, 42.6)
        cell.layers_thicknesses = [0.1, 7.5]  # half-width below the sample radius would cut the sample
        with self.assertRaises(ValueError):
            cell.mavric_deck()

    def test_hot_cell_bare_sample(self):
        cell = hc.HotCellDoses(self.o)
        cell.layers_mats, cell.layers_thicknesses, cell.layers_temperature_K = [], [], []
        deck = cell.mavric_deck()
        self.check_common(deck)
        self.assertIn('media 0 1 99999 -1', deck)
        self.assertAlmostEqual(cell.det_x, cell.cyl_r + cell.det_standoff_distance)

    def test_manual_det_x_warns_and_recomputes(self):
        est = sd.DoseEstimatorSquareTank(self.o)
        est.debug = 0
        est.det_x = 999.0  # direct assignment cannot survive the deck build
        with self.assertWarns(UserWarning) as cm:
            deck = est.mavric_deck()
        self.assertIn('det_standoff_distance', str(cm.warning))
        self.assertNotIn('position 999.0 0 0', deck)
        self.assertAlmostEqual(est.det_x,
                               est.cyl_r + sum(est.layers_thicknesses) + est.det_standoff_distance)
        # Without a manual assignment there is no warning.
        est2 = sd.DoseEstimatorSquareTank(self.o)
        est2.debug = 0
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('error', UserWarning)
            est2.mavric_deck()

    def test_manual_handling_det_x_warns(self):
        est = sd.HandlingContactDoseEstimatorGenericTank(self.o)
        est.debug = 0
        est.cyl_r = 5.0
        est.handling_det_x = 999.0
        with self.assertWarns(UserWarning) as cm:
            est.mavric_deck()
        self.assertIn('handling_det_standoff_distance', str(cm.warning))

    def test_stale_opus_spectra_treated_as_unknown(self):
        import time
        import warnings
        with tempfile.TemporaryDirectory() as tmp:
            case = os.path.join(tmp, 'run_DUMMY')
            os.makedirs(case)
            f71 = os.path.join(case, 'decay.f71')
            base = 'origen'
            plt_names = [base + '.000000000000000000.plt', base + '.000000000000000001.plt',
                         base + '.000000000000000002.plt']
            for name in plt_names:
                open(os.path.join(case, name), 'w').close()
            open(f71, 'w').close()
            now = time.time()
            # Fresh spectra: plt files newer than the F71.
            os.utime(f71, (now - 100.0, now - 100.0))
            for name in plt_names:
                os.utime(os.path.join(case, name), (now, now))

            class DiskOrigen(DummyOrigen):
                def __init__(self):
                    super().__init__()
                    self.cwd = tmp
                    self.decayed_atom_dens = {}  # not run in this process: read OPUS from disk

            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                est = sd.DoseEstimator(DiskOrigen())
            self.assertEqual(est.beta_over_gamma, 0.5)
            self.assertTrue(any('on disk' in str(x.message) for x in w))
            # Stale spectra: F71 regenerated after OPUS -> both treated as unknown.
            os.utime(f71, (now + 100.0, now + 100.0))
            with warnings.catch_warnings(record=True) as w2:
                warnings.simplefilter('always')
                est2 = sd.DoseEstimator(DiskOrigen())
            self.assertIsNone(est2.beta_over_gamma)
            self.assertIsNone(est2.neutron_intensity)
            self.assertTrue(any('predate' in str(x.message) for x in w2))

    def test_hot_cell_reuse_adjoint_flux(self):
        with tempfile.TemporaryDirectory() as tmp:
            o = DummyOrigen()
            o.cwd = tmp
            cell = hc.HotCellDoses(o)
            cell.reuse_adjoint_flux = True
            with self.assertRaises(FileNotFoundError):  # nothing to reuse on a fresh run
                cell.mavric_deck()
            case_dir = os.path.join(tmp, cell._layer_case_dir())  # where run_mavric() runs MAVRIC
            os.makedirs(case_dir)
            dff = os.path.join(case_dir, 'my_dose.adjoint.dff')
            open(dff, 'w').close()
            deck = cell.mavric_deck()
            self.check_common(deck)
            self.assertIn(f'adjointFluxes="{dff}"', deck)
            self.assertEqual(block_ids(deck, 'adjointSource'), [])
            cell.adjoint_flux_file = os.path.join(tmp, 'other.adjoint.dff')  # explicit file from another run
            with self.assertRaises(FileNotFoundError):
                cell.mavric_deck()
