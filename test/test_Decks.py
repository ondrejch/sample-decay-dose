""" Smoke tests: every estimator class must generate a structurally valid MAVRIC deck """
import unittest
from unittest.mock import patch

import sample_decay_dose.SampleDose as sd
import sample_decay_dose.HotCell as hc
import sample_decay_dose.Radiator as rad


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


class DeckTestMixin:
    def check_common(self, deck):
        self.assertIn('=mavric parm=(   )', deck)
        self.assertIn('read parameters', deck)
        self.assertIn('read geometry', deck)
        self.assertIn('read definitions', deck)
        self.assertIn('read importanceMap', deck)
        self.assertIn('read tallies', deck)
        self.assertTrue(deck.rstrip().endswith('end'))
        self.assertEqual(deck.count('pointDetector'), 2 * deck.count('adjointSource') / 2)


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

    @patch('sample_decay_dose.Radiator.get_f71_positions_index',
           return_value={1: {'case': '1', 'time': '100.0'}})
    def test_radiator_box(self, _mock_idx):
        class DummyOrigenTriton(DummyOrigen):
            def __init__(self):
                super().__init__()
                self.BURNED_MATERIAL_F71_file_name = 'burned.f71'
                self.BURNED_MATERIAL_F71_position = 16
                self.burned_atom_dens = {'fe-56': 0.08}

        box = rad.RadiatorBox(DummyOrigenTriton())
        deck = box.mavric_deck()
        self.check_common(deck)
        self.assertIn('array', deck)

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
