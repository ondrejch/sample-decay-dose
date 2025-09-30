import unittest
from unittest.mock import patch, MagicMock, mock_open
import numpy as np
import sys
import types

import sample_decay_dose.SampleDose as sd


def ensure_dummy_isotopes():
    mod_name = 'sample_decay_dose.isotopes'
    if mod_name not in sys.modules:
        m = types.ModuleType(mod_name)
        m.rel_iso_mass = {'h-1': 1.007, 'o-16': 15.995}
        m.m_Da = 1.660539e-24
        sys.modules[mod_name] = m


def ensure_dummy_read_opus():
    mod_name = 'sample_decay_dose.read_opus'
    if mod_name not in sys.modules:
        m = types.ModuleType(mod_name)
        m.integrate_opus = lambda filename: 1.0
        sys.modules[mod_name] = m


class TestSampleDoseFunctions(unittest.TestCase):

    def test_scale_adens(self):
        adens = {'h-1': 1.0, 'o-16': 0.5}
        self.assertEqual(sd.scale_adens(adens, 2.0), {'h-1': 2.0, 'o-16': 1.0})
        self.assertEqual(sd.scale_adens(adens), adens)

    def test_get_rho_from_atom_density(self):
        ensure_dummy_isotopes()
        adens = {'h-1': 0.04, 'o-16': 0.02}
        expected = (1.007 * 0.04 + 15.995 * 0.02) * 1.660539
        self.assertAlmostEqual(sd.get_rho_from_atom_density(adens), expected, places=5)

    def test_get_cyl_r(self):
        V = 10.0
        self.assertAlmostEqual(sd.get_cyl_r(V), (V / (2.0 * np.pi)) ** (1 / 3))

    def test_get_cyl_r_4_1(self):
        V = 10.0
        self.assertAlmostEqual(sd.get_cyl_r_4_1(V), (V / (4.0 * np.pi)) ** (1 / 3))

    def test_get_fill_height_4_1(self):
        cyl_V = 100.0
        fill_V = 50.0
        r = sd.get_cyl_r_4_1(cyl_V)
        expected = fill_V / (np.pi * r ** 2)
        self.assertAlmostEqual(sd.get_fill_height_4_1(fill_V, cyl_V), expected)
        with self.assertRaises(ValueError):
            sd.get_fill_height_4_1(200.0, cyl_V)

    def test_get_cyl_h(self):
        V = 100.0
        r = 5.0
        expected = V / (np.pi * r ** 2)
        self.assertAlmostEqual(sd.get_cyl_h(V, r), expected)
        with self.assertRaises(ValueError):
            sd.get_cyl_h(10.0, 0.0)
        with self.assertRaises(ValueError):
            sd.get_cyl_h(10.0, -1.0)

    @patch('sample_decay_dose.SampleDose.subprocess.run')
    def test_run_scale_success(self, mock_run):
        proc = MagicMock()
        proc.stdout.decode.return_value = "All good\n"
        mock_run.return_value = proc
        self.assertTrue(sd.run_scale('deck.inp'))

    @patch('sample_decay_dose.SampleDose.subprocess.run')
    def test_run_scale_failure(self, mock_run):
        proc = MagicMock()
        proc.stdout.decode.return_value = "Error: something failed\n"
        mock_run.return_value = proc
        self.assertFalse(sd.run_scale('deck.inp'))

    def test_atom_dens_for_origen(self):
        adens = {'h-1': 1.0, 'o-16': 0.5}
        expected = 'h-1 = 1.0 \no-16 = 0.5 \n'
        self.assertEqual(sd.atom_dens_for_origen(adens), expected)

    def test_atom_dens_for_mavric(self):
        adens = {'h-1': 1.0, 'o-16': 0.5, 'u-235m': 0.1}
        expected = (
            'h-1 1 0 1.0 873.0 end\n'
            'o-16 1 0 0.5 873.0 end\n'
            'u-235 1 0 0.1 873.0 end\n'
        )
        self.assertEqual(sd.atom_dens_for_mavric(adens), expected)

    @patch('sample_decay_dose.SampleDose.subprocess.run')
    def test_get_f71_positions_index(self, mock_run):
        # Tokens: id time power flux fluence energy initialhm libpos case step DCGNAB
        mock_run.return_value.stdout = (
            b"1 0.0 0.0 0.0 0.0 0.0 0.0 1 1 1 1\n"
            b"2 1.0 1.0 1.0 1.0 1.0 1.0 2 1 2 1\n"
            b"state definition present\n"
        )
        idx = sd.get_f71_positions_index('x.f71')
        self.assertIn(1, idx)
        self.assertIn(2, idx)
        # case is token 8 (zero-based) -> '1' for the second line
        self.assertEqual(idx[2]['case'], '1')
        # optional: also check libpos for clarity
        self.assertEqual(idx[2]['libpos'], '2')

    @patch('sample_decay_dose.SampleDose.subprocess.run')
    def test_get_burned_nuclide_atom_dens(self, mock_run):
        mock_run.return_value.stdout = (
            b"case,1,2,3,4,5\n"
            b"U235,0,0,0,0,1.20e-3\n"
            b"Pu239,0,0,0,0,5.00e-5\n"
        )
        dens = sd.get_burned_nuclide_atom_dens('x.f71', 5)
        self.assertAlmostEqual(dens['u-235'], 1.20e-3)
        self.assertAlmostEqual(dens['pu-239'], 5.00e-5)


class TestOrigenFromTritonMHA(unittest.TestCase):

    @patch('sample_decay_dose.SampleDose.get_last_position_for_case', return_value=20)
    @patch('sample_decay_dose.SampleDose.get_f71_positions_index',
           return_value={i: {'case': '20', 'time': str(i)} for i in range(1, 21)})
    def setUp(self, mock_idx, mock_last):
        self.o = sd.OrigenFromTritonMHA(_f71='core.f71', _MTiHM=0.5, _f71_case=20)
        self.o.debug = 0

    def test_init(self):
        self.assertEqual(self.o.BURNED_MATERIAL_F71_position, 19)
        self.assertEqual(self.o.MTiHM, 0.5)

    def test_origen_deck(self):
        deck = self.o.origen_deck()
        self.assertIn('OrigenFromTritonMHA', deck)
        self.assertIn(f'pos={self.o.BURNED_MATERIAL_F71_position}', deck)
        self.assertIn(f'retained=[Te={self.o.MTiHM}', deck)

    @patch('os.path.exists', return_value=False)
    @patch('os.mkdir')
    @patch('os.chdir')
    @patch('builtins.open', new_callable=mock_open)
    @patch('sample_decay_dose.SampleDose.run_scale', return_value=True)
    @patch('sample_decay_dose.SampleDose.get_burned_nuclide_atom_dens', return_value={'u-235': 1e-3})
    def test_run_decay_sample(self, mock_get, mock_run, mock_file, mock_chdir, mock_mkdir, mock_exists):
        self.o.burned_atom_dens = {'u-238': 0.1}
        self.o.run_decay_sample()
        mock_mkdir.assert_called_with(self.o.case_dir)
        mock_run.assert_called_with(self.o.ORIGEN_input_file_name)
        self.assertIn('u-235', self.o.decayed_atom_dens)


class TestDoseEstimator(unittest.TestCase):

    def setUp(self):
        ensure_dummy_read_opus()
        mock_o = MagicMock(spec=sd.Origen)
        mock_o.sample_weight = 10.0
        mock_o.sample_density = 7.8
        mock_o.sample_volume = 10.0 / 7.8
        mock_o.SAMPLE_F71_file_name = 'decay.f71'
        mock_o.SAMPLE_F71_position = 12
        mock_o.SAMPLE_DECAY_days = 30.0
        mock_o.decayed_atom_dens = {'fe-56': 0.08}
        mock_o.get_beta_to_gamma.return_value = 0.5
        mock_o.get_neutron_integral.return_value = 1.0e5
        mock_o.case_dir = 'run_10.0_g-30_days'
        mock_o.cwd = '/tmp'
        self.d = sd.DoseEstimator(mock_o)
        self.d.debug = 0

    def test_init(self):
        self.assertEqual(self.d.sample_weight, 10.0)
        self.assertEqual(self.d.case_dir, self.d.ORIGEN_dir + '_MAVRIC')
        self.assertEqual(self.d.beta_over_gamma, 0.5)

    @patch('os.path.isfile', return_value=True)
    @patch('os.path.exists', return_value=False)
    @patch('os.mkdir')
    @patch('os.chdir')
    @patch('shutil.copy2')
    @patch('builtins.open', new_callable=mock_open)
    @patch('sample_decay_dose.SampleDose.run_scale', return_value=True)
    def test_run_mavric(self, mock_run, mock_file, mock_copy, mock_chdir, mock_mkdir, mock_exists, mock_isfile):
        self.d.run_mavric()
        mock_mkdir.assert_called_with(self.d.case_dir)
        mock_run.assert_called_with(self.d.MAVRIC_input_file_name)

    @patch('os.path.isfile', return_value=True)
    @patch('os.chdir')
    def test_get_responses(self, mock_chdir, mock_isfile):
        mock_text = (
            " Header stuff\n"
            " Final Tally Results Summary\n"
            " response 1 1.0e-3 1.0e-4\n"
            " response 2 2.0e-3 2.0e-4\n"
        )
        m = mock_open(read_data=mock_text)
        with patch('builtins.open', m):
            self.d.get_responses()
        self.assertAlmostEqual(self.d.responses['1']['value'], 1.0e-3)
        self.assertAlmostEqual(self.d.responses['2']['value'], 2.0e-3)
        self.assertAlmostEqual(self.d.responses['3']['value'], 0.5 * 2.0e-3)

    def test_total_dose(self):
        self.d.responses = {
            '1': {'value': 1.0, 'stdev': 0.1},
            '2': {'value': 2.0, 'stdev': 0.2},
            '3': {'value': 3.0, 'stdev': 0.3},
        }
        total = self.d.total_dose
        self.assertAlmostEqual(total['value'], 6.0)
        expected_sd = 6.0 * np.sqrt((0.1 / 1.0) ** 2 + (0.2 / 2.0) ** 2 + (0.3 / 3.0) ** 2)
        self.assertAlmostEqual(total['stdev'], expected_sd)


if __name__ == '__main__':
    ensure_dummy_isotopes()
    ensure_dummy_read_opus()
    unittest.main()
