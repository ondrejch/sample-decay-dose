import unittest
from unittest.mock import patch, MagicMock, mock_open
import numpy as np
import os
import re
import sys
import types

import sample_decay_dose.SampleDose as sd
from sample_decay_dose import utils
from sample_decay_dose.read_opus import integrate_opus


def dummy_isotopes() -> types.ModuleType:
    """ Stand-in isotopes module with known masses; install it with patch.dict(sys.modules, ...) """
    m = types.ModuleType('sample_decay_dose.isotopes')
    m.rel_iso_mass = {'h-1': 1.007, 'o-16': 15.995}
    m.m_Da = 1.660539e-24
    return m


def ensure_dummy_read_opus():
    mod_name = 'sample_decay_dose.read_opus'
    if mod_name not in sys.modules:
        m = types.ModuleType(mod_name)
        m.integrate_opus = lambda filename: 1.0
        sys.modules[mod_name] = m


def mavric_tally_summary(detectors: dict) -> str:
    """ MAVRIC output tail in the format of a real 'Final Tally Results Summary' block.
    detectors: {detector_id: (particle header, title word, response id, value, stdev)} """
    text = (" pointDetector 1\n        title=\"neutron detector\"\n"  # input echo, not a result
            "\n\n Final Tally Results Summary\n ============================\n\n"
            "     Final Statistical Checks (fits are over the last half of the simulation)\n\n"
            "     quantity                 check                     goal         passing\n"
            "     1 mean                   rel slope of linear fit   0.00  |slope| < 0.10\n\n")
    for det, (header, word, rid, value, stdev) in detectors.items():
        text += (f"\n {header} Point Detector {det}.  {word} detector\n"
                 "                         average      standard     relat      FOM    stat checks\n"
                 "    tally/quantity        value       deviation    uncert   (/min)   1 2 3 4 5 6\n"
                 "    ------------------  -----------  -----------  -------  --------  -----------\n"
                 "    uncollided flux     1.82626E+02  1.39784E-01  0.00077\n"
                 "    total flux          5.34887E+02  2.77958E+01  0.05197  1.44E+02  - - - - - -\n"
                 f"    response {rid}          {value:.5E}  {stdev:.5E}  0.00330  3.57E+04  X - X X - -\n"
                 "    ------------------  -----------  -----------  -------  --------  -----------\n")
    text += "\n Total Monaco cpu time for this problem was  2.65 minutes\n"
    return text


class TestSampleDoseFunctions(unittest.TestCase):

    def test_scale_adens(self):
        adens = {'h-1': 1.0, 'o-16': 0.5}
        self.assertEqual(utils.scale_adens(adens, 2.0), {'h-1': 2.0, 'o-16': 1.0})
        self.assertEqual(utils.scale_adens(adens), adens)

    def test_get_rho_from_atom_density(self):
        adens = {'h-1': 0.04, 'o-16': 0.02}
        expected = (1.007 * 0.04 + 15.995 * 0.02) * 1.660539
        with patch.dict(sys.modules, {'sample_decay_dose.isotopes': dummy_isotopes()}):
            self.assertAlmostEqual(utils.get_rho_from_atom_density(adens), expected, places=5)

    def test_get_cyl_r(self):
        V = 10.0
        self.assertAlmostEqual(utils.get_cyl_r(V), (V / (2.0 * np.pi)) ** (1 / 3))

    def test_get_cyl_r_4_1(self):
        V = 10.0
        self.assertAlmostEqual(utils.get_cyl_r_4_1(V), (V / (4.0 * np.pi)) ** (1 / 3))

    def test_get_fill_height_4_1(self):
        cyl_V = 100.0
        fill_V = 50.0
        r = utils.get_cyl_r_4_1(cyl_V)
        expected = fill_V / (np.pi * r ** 2)
        self.assertAlmostEqual(utils.get_fill_height_4_1(fill_V, cyl_V), expected)
        with self.assertRaises(ValueError):
            utils.get_fill_height_4_1(200.0, cyl_V)

    def test_get_cyl_h(self):
        V = 100.0
        r = 5.0
        expected = V / (np.pi * r ** 2)
        self.assertAlmostEqual(utils.get_cyl_h(V, r), expected)
        with self.assertRaises(ValueError):
            utils.get_cyl_h(10.0, 0.0)
        with self.assertRaises(ValueError):
            utils.get_cyl_h(10.0, -1.0)

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_run_scale_success(self, mock_run):
        proc = MagicMock()
        proc.returncode = 0
        proc.stdout.decode.return_value = "All good\n"
        proc.stderr.decode.return_value = ""
        mock_run.return_value = proc
        self.assertTrue(utils.run_scale('deck.inp'))

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_run_scale_failure(self, mock_run):
        # scalerte exits with code 0 after an ORIGEN abort; the "***Error:" marker is the failure signal
        proc = MagicMock()
        proc.returncode = 0
        proc.stdout.decode.return_value = "  Now executing origen\n***Error: Time step dt=-nan must be positive 2/11\n"
        proc.stderr.decode.return_value = ""
        mock_run.return_value = proc
        self.assertFalse(utils.run_scale('deck.inp'))

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_run_scale_ignores_error_words_in_normal_output(self, mock_run):
        proc = MagicMock()
        proc.returncode = 0
        proc.stdout.decode.return_value = ("relative error of the tally is 0.01\nBacktrace for this error: none\n"
                                           "Scale job deck.inp is finished.\n")
        proc.stderr.decode.return_value = ""
        mock_run.return_value = proc
        self.assertTrue(utils.run_scale('deck.inp'))

    def test_run_scale_reads_failure_from_output_file(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            deck = os.path.join(tmp, 'deck.inp')

            def fake_run(*_args, **_kwargs):
                with open(os.path.join(tmp, 'deck.out'), 'w') as f:
                    f.write("-------------------------- Summary --------------------------\n"
                            "    terminated due to errors. completion code 134. used 1.4 seconds\n"
                            "origen failed. used 1.4474 seconds.\n")
                proc = MagicMock()
                proc.returncode = 0
                proc.stdout = b"Scale job deck.inp is finished.\n"
                proc.stderr = b""
                return proc

            with patch('sample_decay_dose.utils.subprocess.run', side_effect=fake_run):
                self.assertFalse(utils.run_scale(deck))

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_run_scale_nonzero_return_code(self, mock_run):
        proc = MagicMock()
        proc.returncode = 1
        proc.stdout.decode.return_value = "All good\n"
        proc.stderr.decode.return_value = ""
        mock_run.return_value = proc
        self.assertFalse(utils.run_scale('deck.inp'))

    def test_atom_dens_for_origen(self):
        adens = {'h-1': 1.0, 'o-16': 0.5}
        expected = 'h-1 = 1.0 \no-16 = 0.5 \n'
        self.assertEqual(utils.atom_dens_for_origen(adens), expected)

    def test_atom_dens_for_mavric(self):
        adens = {'h-1': 1.0, 'o-16': 0.5, 'u-235m': 0.1}
        expected = (
            'h-1 1 0 1.0 873.0 end\n'
            'o-16 1 0 0.5 873.0 end\n'
            'u-235 1 0 0.1 873.0 end\n'
        )
        self.assertEqual(utils.atom_dens_for_mavric(adens), expected)

    def test_atom_dens_for_mavric_keeps_natural_sm_tm(self):
        out = utils.atom_dens_for_mavric({'sm': 1.0, 'tm': 2.0, 'am242m': 3.0}, 5, 300.0)
        self.assertEqual(out, 'sm 5 0 1.0 300.0 end\ntm 5 0 2.0 300.0 end\nam242 5 0 3.0 300.0 end\n')

    def test_get_rho_handles_bound_hydrogen_and_metastables(self):
        fake = types.ModuleType('sample_decay_dose.isotopes')
        fake.rel_iso_mass = {'h-1': 1.007825, 'c-12': 12.0, 'c-13': 13.003355, 'sm': 150.36, 'am-242': 242.0595}
        fake.m_Da = 1.660539e-24
        with patch.dict(sys.modules, {'sample_decay_dose.isotopes': fake}):
            rho_hdpe = utils.get_rho_from_atom_density(sd.ADENS_HDPE_COLD)
            expected_hdpe = (12.0 * 3.992647e-02 + 13.003355 * 4.318339e-04 + 1.007825 * 8.071660e-02) * 1.660539
            self.assertAlmostEqual(rho_hdpe, expected_hdpe, places=6)
            self.assertAlmostEqual(rho_hdpe, 0.94, delta=0.02)  # HDPE is about 0.94 g/cm3
            rho = utils.get_rho_from_atom_density({'sm': 1e-3, 'am-242m': 1e-4})
            self.assertAlmostEqual(rho, (150.36 * 1e-3 + 242.0595 * 1e-4) * 1.660539, places=6)

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_atom_dens_minimum_is_read_at_call_time(self, mock_run):
        from sample_decay_dose import constants
        mock_run.return_value.returncode = 0
        mock_run.return_value.stderr = b""
        mock_run.return_value.stdout = b"case,1\nU235,1.0e-3\nPu239,1.0e-8\n"
        with patch.object(constants, 'ATOM_DENS_MINIMUM', 1e-6):
            self.assertEqual(list(utils.get_burned_nuclide_atom_dens('x.f71', 1)), ['u-235'])
        self.assertEqual(list(utils.get_burned_nuclide_atom_dens('x.f71', 1)), ['u-235', 'pu-239'])

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_get_F33_num_sets_matches_basename_exactly(self, mock_run):
        mock_run.return_value.returncode = 0
        mock_run.return_value.stderr = b""
        mock_run.return_value.stdout = (
            b"               file        dataType numSets  fileFormat version           directory\n"
            b" EIRENE.mix0002.f33 Origen::Library       2 bof(binary)   6.3.1 examples/irradiator\n")
        self.assertEqual(utils.get_F33_num_sets('examples/irradiator/EIRENE.mix0002.f33'), 2)
        with self.assertRaises(RuntimeError):  # 'EIRENE.mix0002.f33' is a substring of this path
            utils.get_F33_num_sets('examples/irradiator/xEIRENE.mix0002.f33')

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_get_f71_volume(self, mock_run):
        mock_run.return_value.returncode = 0
        mock_run.return_value.stderr = b""
        mock_run.return_value.stdout = (b"          case,             1,             2\n"
                                        b"        volume, 1.1293750E+05, 2.0000000E+00\n"
                                        b"      Co060, 5.0000E-04, 5.0000E-04\n")
        self.assertAlmostEqual(utils.get_f71_volume('x.f71', 1), 1.129375e5)
        self.assertAlmostEqual(utils.get_f71_volume('x.f71', 2), 2.0)
        with self.assertRaises(RuntimeError):
            utils.get_f71_volume('x.f71', 5)

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_get_f71_positions_index(self, mock_run):
        # Tokens: id time power flux fluence energy initialhm libpos case step DCGNAB
        mock_run.return_value.returncode = 0
        mock_run.return_value.stderr = b""
        mock_run.return_value.stdout = (
            b"1 0.0 0.0 0.0 0.0 0.0 0.0 1 1 1 1\n"
            b"2 1.0 1.0 1.0 1.0 1.0 1.0 2 1 2 1\n"
            b"state definition present\n"
        )
        idx = utils.get_f71_positions_index('x.f71')
        self.assertIn(1, idx)
        self.assertIn(2, idx)
        # case is token 8 (zero-based) -> '1' for the second line
        self.assertEqual(idx[2]['case'], '1')
        # optional: also check libpos for clarity
        self.assertEqual(idx[2]['libpos'], '2')

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_get_f71_positions_index_raises_on_failed_subprocess(self, mock_run):
        proc = MagicMock()
        proc.returncode = 2
        proc.stdout = b""
        proc.stderr = b"obiwan failed"
        mock_run.return_value = proc
        with self.assertRaises(RuntimeError):
            utils.get_f71_positions_index('x.f71')

    @patch('sample_decay_dose.utils.subprocess.run')
    def test_get_burned_nuclide_atom_dens(self, mock_run):
        mock_run.return_value.returncode = 0
        mock_run.return_value.stderr = b""
        mock_run.return_value.stdout = (
            b"case,1,2,3,4,5\n"
            b"U235,0,0,0,0,1.20e-3\n"
            b"Pu239,0,0,0,0,5.00e-5\n"
        )
        dens = utils.get_burned_nuclide_atom_dens('x.f71', 5)
        self.assertAlmostEqual(dens['u-235'], 1.20e-3)
        self.assertAlmostEqual(dens['pu-239'], 5.00e-5)

    def test_as_nuclide_name_rejects_embedded_garbage(self):
        self.assertIsNone(utils._as_nuclide_name("prefix u235 suffix"))
        self.assertEqual(utils._as_nuclide_name("U-235"), "u-235")
        self.assertEqual(utils._as_nuclide_name("Xe135m"), "xe-135m")


def expand_origen_times(time_line: str) -> list[float]:
    """ Expands an ORIGEN 't=[nL a b]' or 't=[a ...]' list into times [days] """
    body = time_line.strip()[len('t=['):-1].split()
    if len(body) == 3 and body[0].endswith('L'):
        n = int(body[0][:-1])
        a, b = float(body[1]), float(body[2])
        return [a] + list(np.geomspace(a, b, n + 2)[1:-1]) + [b]
    return [float(x) for x in body]


class TestDecayTimeGrid(unittest.TestCase):

    def test_grids_are_strictly_increasing(self):
        for days in (30.0, 1.0, 0.0005, 5e-5, 1e-7):
            for t_first in (1e-4, 1e-3):
                times = [0.0] + expand_origen_times(sd.decay_time_grid(days, 12, t_first))
                self.assertEqual(len(times), 12)  # start=0 position plus 11 time points
                self.assertTrue(all(b > a for a, b in zip(times, times[1:])), (days, t_first, times))
                self.assertAlmostEqual(times[-1], days)
                self.assertEqual(sd.decayed_f71_position(days, 12), 12)

    def test_zero_day_reads_start_position(self):
        line = sd.decay_time_grid(0.0, 12)
        times = expand_origen_times(line)
        self.assertEqual(times, [sd.ZERO_DAY_STEP_DAYS])
        self.assertEqual(sd.decayed_f71_position(0.0, 12), 1)

    def test_invalid_grids_raise(self):
        with self.assertRaises(ValueError):
            sd.decay_time_grid(-1.0, 12)
        with self.assertRaises(ValueError):
            sd.decay_time_grid(30.0, 3)


class TestOrigenDecks(unittest.TestCase):

    def _time_line(self, deck: str, occurrence: int = -1) -> str:
        lines = [ln.strip() for ln in deck.splitlines() if ln.strip().startswith('t=[')]
        return lines[occurrence]

    def _npos(self, deck: str) -> list[int]:
        return [int(x) for x in re.findall(r'npos=(\d+) end', deck)]

    @patch('sample_decay_dose.SampleDose.get_f71_positions_index', return_value={})
    def test_origen_from_triton_deck(self, _mock_idx):
        o = sd.OrigenFromTriton(_f71='core.f71', _mass=1.0)
        o.debug = 0
        o.sample_volume = 0.5
        for days, npos in ((30.0, 12), (5e-5, 12), (0.0, 1)):
            o.set_decay_days(days)
            deck = o.origen_deck()
            times = expand_origen_times(self._time_line(deck))
            self.assertTrue(all(b > a for a, b in zip([0.0] + times, times)))
            self.assertEqual(self._npos(deck), [npos] * 3)
            self.assertEqual(o.decayed_F71_position, npos)
            self.assertIn('start=0', deck)
            self.assertIn('volume=0.5', deck)
            self.assertEqual(deck.count('=opus'), 3)

    @patch('sample_decay_dose.SampleDose.get_F33_num_sets', return_value=2)
    def test_origen_irradiation_deck(self, _mock_sets):
        o = sd.OrigenIrradiation(_f33='libs/EIRENE.mix0002.f33', _mass=1.0)
        o.sample_density = 8.0
        for days, npos in ((1.0 / 24.0, 12), (5e-5, 12), (0.0, 1)):
            o.set_decay_days(days)
            deck = o.origen_deck()
            self.assertIn('file="EIRENE.mix0002.f33"', deck)
            self.assertIn('pos=2', deck)
            self.assertIn(f'cp -r {o.cwd}/libs/EIRENE.mix0002.f33 .', deck)
            irr_times = expand_origen_times(self._time_line(deck, 0))
            self.assertEqual(len(irr_times), o.irradiate_steps)
            self.assertIn(f'flux=[{o.irradiate_steps}R {o.irradiate_flux}]', deck)
            decay_times = expand_origen_times(self._time_line(deck, 1))
            self.assertTrue(all(b > a for a, b in zip([0.0] + decay_times, decay_times)))
            self.assertEqual(self._npos(deck), [npos] * 3)
        o.irradiate_days = 0.0
        with self.assertRaises(ValueError):
            o.origen_deck()

    def test_origen_decay_box_deck_and_case_dir(self):
        o = sd.OrigenDecayBox({'co-60': 1e-3}, 2.0)
        for days, npos in ((30.0, 12), (0.0005, 12), (0.0, 1)):
            o.set_decay_days(days)
            deck = o.origen_deck()
            times = expand_origen_times(self._time_line(deck))
            self.assertTrue(all(b > a for a, b in zip([0.0] + times, times)), times)
            self.assertEqual(self._npos(deck), [npos] * 3)
            self.assertIn('volume=2.0', deck)
        # The case directory tracks volume, decay time, and composition
        o.set_decay_days(30.0)
        dir_30 = o.case_dir
        o.set_decay_days(60.0)
        self.assertNotEqual(o.case_dir, dir_30)
        self.assertIn('_cm3-', o.case_dir)
        other = sd.OrigenDecayBox({'cs-137': 1e-3}, 2.0)
        other.set_decay_days(60.0)
        self.assertNotEqual(other.case_dir, o.case_dir)

    @patch('os.path.isfile', return_value=True)
    @patch('sample_decay_dose.SampleDose.get_burned_material_total_mass_dens', return_value=7.5)
    @patch('sample_decay_dose.SampleDose.get_f71_positions_index', return_value={1: {}, 2: {}, 31: {}})
    @patch('sample_decay_dose.SampleDose.get_F33_num_sets', return_value=2)
    def test_read_irradiated_material_density(self, _sets, _idx, mock_rho, _isfile):
        o = sd.OrigenIrradiation(_f33='x.f33', _mass=15.0)
        o.read_irradiated_material_density()
        mock_rho.assert_called_with(os.path.join(o.cwd, o.case_dir, 'irradiate.f71'), 31)  # end of irradiation
        self.assertAlmostEqual(o.sample_volume, 2.0)


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

    def test_zero_day_deck_keeps_later_decks_valid(self):
        self.o.set_decay_days(0.0)
        deck0 = self.o.origen_deck()
        self.assertIn(f't=[{sd.ZERO_DAY_STEP_DAYS}]', deck0)
        self.assertEqual(re.findall(r'npos=(\d+) end', deck0), ['1'] * 3)
        self.assertEqual(self.o.SAMPLE_F71_position, 12)  # configuration is unchanged
        self.o.set_decay_days(30.0)
        deck30 = self.o.origen_deck()  # raised "Too few time steps" before
        self.assertIn('t=[9L 0.0001 30.0]', deck30)
        self.assertEqual(self.o.decayed_F71_position, 12)

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
        mock_run.assert_called_with(self.o.ORIGEN_input_file_name, 1)
        self.assertIn('u-235', self.o.decayed_atom_dens)

    @patch('os.path.exists', return_value=False)
    @patch('os.mkdir')
    @patch('os.chdir')
    @patch('builtins.open', new_callable=mock_open)
    @patch('sample_decay_dose.SampleDose.get_burned_nuclide_atom_dens', return_value={'u-235': 1e-3})
    def test_run_decay_sample_raises_on_scale_failure(self, mock_get, mock_file, mock_chdir,
                                                      mock_mkdir, mock_exists):
        self.o.burned_atom_dens = {'u-238': 0.1}
        with patch('sample_decay_dose.SampleDose.run_scale', return_value=False):
            with self.assertRaises(RuntimeError):
                self.o.run_decay_sample()


class TestF71PositionSelection(unittest.TestCase):

    @patch('sample_decay_dose.SampleDose.get_f71_positions_index', return_value={
        10: {'case': '1', 'time': '100.0'},
        20: {'case': '1', 'time': '200.0'},
        30: {'case': '1', 'time': '300.0'},
    })
    def test_set_f71_pos_uses_f71_position_key(self, _mock_idx):
        o = sd.OrigenFromTriton(_f71='core.f71', _mass=1.0)
        o.set_f71_pos(260.0, case='1')
        self.assertEqual(o.BURNED_MATERIAL_F71_position, 30)

    @patch('sample_decay_dose.SampleDose.get_f71_positions_index', return_value={
        10: {'case': '1', 'time': '100.0'},
        20: {'case': '1', 'time': '200.0'},
        30: {'case': '1', 'time': '300.0'},
        31: {'case': '2', 'time': '300.0'},
        32: {'case': '2', 'time': '300.0'},
        33: {'case': '2', 'time': '400.0'},
    })
    def test_set_f71_pos_picks_closest_time(self, _mock_idx):
        o = sd.OrigenFromTriton(_f71='core.f71', _mass=1.0)
        for t, expected in ((201.0, 20), (250.0, 20), (251.0, 30), (100.0, 10), (50.0, 10), (900.0, 30)):
            o.set_f71_pos(t, case='1')
            self.assertEqual(o.BURNED_MATERIAL_F71_position, expected, t)
        o.set_f71_pos(300.0, case='2')  # equal times: the lower position
        self.assertEqual(o.BURNED_MATERIAL_F71_position, 31)


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
        mock_run.assert_called_with(self.d.MAVRIC_input_file_name, 1)

    @patch('os.path.isfile', return_value=True)
    @patch('os.chdir')
    def test_get_responses(self, mock_chdir, mock_isfile):
        mock_text = mavric_tally_summary({'1': ('Neutron', 'neutron', '1', 1.0e-3, 1.0e-4),
                                          '2': ('Photon', 'photon', '2', 2.0e-3, 2.0e-4)})
        m = mock_open(read_data=mock_text)
        with patch('builtins.open', m):
            self.d.get_responses()
        self.assertAlmostEqual(self.d.responses['1']['value'], 1.0e-3)
        self.assertAlmostEqual(self.d.responses['1']['stdev'], 1.0e-4)
        self.assertAlmostEqual(self.d.responses['2']['value'], 2.0e-3)
        self.assertAlmostEqual(self.d.responses['3']['value'], 0.5 * 2.0e-3)
        self.assertAlmostEqual(self.d.responses['3']['stdev'], 0.5 * 2.0e-4)

    @patch('os.path.isfile', return_value=True)
    @patch('os.chdir')
    def test_get_responses_raises_without_neutron_tally(self, mock_chdir, mock_isfile):
        mock_text = mavric_tally_summary({'2': ('Photon', 'photon', '2', 2.0e-3, 2.0e-4)})
        with patch('builtins.open', mock_open(read_data=mock_text)):
            with self.assertRaises(RuntimeError):
                self.d.get_responses()

    def test_total_dose(self):
        self.d.responses = {
            '1': {'value': 1.0, 'stdev': 0.1},
            '2': {'value': 2.0, 'stdev': 0.2},
            '3': {'value': 3.0, 'stdev': 0.3},
        }
        total = self.d.total_dose
        self.assertAlmostEqual(total['value'], 6.0)
        # Beta '3' is k * gamma, so sigma_2 and sigma_3 add linearly: sqrt(0.1^2 + (0.2 + 0.3)^2) = sqrt(0.26)
        self.assertAlmostEqual(total['stdev'], 0.5099019513592785)

    def test_total_dose_review_example(self):
        # n = 1 +/- 0.1 and gamma = 100 +/- 1 without beta: sigma = sqrt(0.01 + 1) = 1.005
        self.d.responses = {'1': {'value': 1.0, 'stdev': 0.1}, '2': {'value': 100.0, 'stdev': 1.0},
                            '3': {'value': 0.0, 'stdev': 0.0}}
        self.assertAlmostEqual(self.d.total_dose['value'], 101.0)
        self.assertAlmostEqual(self.d.total_dose['stdev'], 1.0049875621120890)

    def test_total_dose_requires_responses(self):
        self.d.responses = {}
        with self.assertRaises(RuntimeError):
            _ = self.d.total_dose


class TestDoseEstimatorWithoutOrigenRun(unittest.TestCase):
    """ ORIGEN case not run in this process: decayed_atom_dens is empty """

    def _mock_origen(self, beta=None, neutron=None):
        mock_o = MagicMock(spec=sd.Origen)
        mock_o.sample_weight = 10.0
        mock_o.sample_density = 7.8
        mock_o.sample_volume = 10.0 / 7.8
        mock_o.SAMPLE_F71_file_name = 'decay.f71'
        mock_o.SAMPLE_F71_position = 12
        mock_o.SAMPLE_DECAY_days = 30.0
        mock_o.decayed_atom_dens = {}
        mock_o.get_beta_to_gamma.side_effect = FileNotFoundError('no plt') if beta is None else None
        mock_o.get_beta_to_gamma.return_value = beta
        mock_o.get_neutron_integral.side_effect = FileNotFoundError('no plt') if neutron is None else None
        mock_o.get_neutron_integral.return_value = neutron
        mock_o.case_dir = 'run_10.0_g-30_days'
        mock_o.cwd = '/tmp'
        return mock_o

    def test_reads_earlier_opus_spectra(self):
        d = sd.DoseEstimator(self._mock_origen(beta=0.25, neutron=0.0))
        self.assertEqual(d.beta_over_gamma, 0.25)
        self.assertEqual(d.neutron_intensity, 0.0)
        d.debug = 0
        deck = d.mavric_deck()
        self.assertNotIn('src 1', deck)  # a known zero neutron intensity omits the neutron source

    def test_unknown_spectra(self):
        d = sd.DoseEstimator(self._mock_origen())
        d.debug = 0
        self.assertIsNone(d.beta_over_gamma)
        self.assertIsNone(d.neutron_intensity)
        deck = d.mavric_deck()
        self.assertIn('src 1', deck)  # unknown neutron intensity includes the neutron source
        self.assertIn('distribution 1', deck)
        with patch('os.path.isfile', return_value=True), patch('os.chdir'), \
                patch('builtins.open', mock_open(read_data=mavric_tally_summary(
                    {'1': ('Neutron', 'neutron', '1', 1e-3, 1e-4), '2': ('Photon', 'photon', '2', 2e-3, 2e-4)}))):
            with self.assertRaises(RuntimeError):  # bare sample: beta applies but the ratio is unknown
                d.get_responses()
            d.beta_over_gamma = 0.5
            d.get_responses()
        self.assertAlmostEqual(d.responses['3']['value'], 1e-3)

    def test_shielded_estimator_reports_zero_beta(self):
        d = sd.DoseEstimatorGenericTank(self._mock_origen())
        d.debug = 0
        self.assertFalse(d.beta_applies)
        with patch('os.path.isfile', return_value=True), patch('os.chdir'), \
                patch('builtins.open', mock_open(read_data=mavric_tally_summary(
                    {'1': ('Neutron', 'neutron', '1', 1e-3, 1e-4), '2': ('Photon', 'photon', '2', 2e-3, 2e-4)}))):
            d.get_responses()
        self.assertEqual(d.responses['3'], {'value': 0.0, 'stdev': 0.0})
        d.layers_mats, d.layers_thicknesses, d.layers_temperature_K = [], [], []
        self.assertTrue(d.beta_applies)


class TestHandlingContactResponses(unittest.TestCase):

    def setUp(self):
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
        self.o = mock_o
        self.text = mavric_tally_summary({
            '1': ('Neutron', 'neutron', '1', 1.0, 0.1),
            '2': ('Photon', 'photon', '2', 100.0, 1.0),
            '5': ('Neutron', 'neutron', '1', 0.1, 0.01),
            '6': ('Photon', 'photon', '2', 10.0, 0.2),
        })

    def _read(self, est, text=None):
        est.debug = 0
        with patch('os.path.isfile', return_value=True), \
                patch('builtins.open', mock_open(read_data=self.text if text is None else text)):
            est.get_responses()

    def test_contact_and_handling_doses(self):
        for cls in (sd.HandlingContactDoseEstimatorGenericTank, sd.MHATank):
            est = cls(self.o)
            self._read(est)
            self.assertEqual(est.responses['5']['particle'], 'neutron')
            self.assertEqual(est.responses['5']['pid'], '1')
            self.assertEqual(est.responses['3']['value'], 0.0)  # shielded by the default layers
            self.assertAlmostEqual(est.contact_dose['value'], 101.0)
            self.assertAlmostEqual(est.contact_dose['stdev'], np.sqrt(0.1 ** 2 + 1.0 ** 2))
            self.assertAlmostEqual(est.handling_dose['value'], 10.1)
            self.assertAlmostEqual(est.handling_dose['stdev'], np.sqrt(0.01 ** 2 + 0.2 ** 2))
            self.assertEqual(est.total_dose, est.contact_dose)

    def test_bare_sample_beta(self):
        est = sd.HandlingContactDoseEstimatorGenericTank(self.o)
        est.layers_mats, est.layers_thicknesses, est.layers_temperature_K = [], [], []
        self._read(est)
        self.assertAlmostEqual(est.responses['3']['value'], 50.0)
        self.assertAlmostEqual(est.responses['7']['value'], 5.0)
        self.assertAlmostEqual(est.contact_dose['value'], 151.0)
        self.assertAlmostEqual(est.contact_dose['stdev'], np.sqrt(0.1 ** 2 + (1.0 + 0.5) ** 2))
        self.assertAlmostEqual(est.handling_dose['stdev'], np.sqrt(0.01 ** 2 + (0.2 + 0.1) ** 2))

    def test_missing_detector_raises(self):
        est = sd.HandlingContactDoseEstimatorGenericTank(self.o)
        text = mavric_tally_summary({'1': ('Neutron', 'neutron', '1', 1.0, 0.1),
                                     '2': ('Photon', 'photon', '2', 100.0, 1.0)})
        with self.assertRaises(RuntimeError):
            self._read(est, text)
        with self.assertRaises(RuntimeError):
            self._read(est, 'no tallies here\n')

    def test_hot_cell_inherits_parser(self):
        import sample_decay_dose.HotCell as hc
        cell = hc.HotCellDoses(self.o)
        self._read(cell)
        self.assertAlmostEqual(cell.handling_dose['value'], 10.1)


class TestReadOpus(unittest.TestCase):

    def _write_plt(self, tmpdir, lines):
        fname = os.path.join(tmpdir, 'spectrum.plt')
        with open(fname, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return fname

    def test_integrate_opus(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            # OPUS format: 6 header lines, then (energy, intensity) pairs;
            # a bin is terminated when the next line repeats the previous intensity.
            lines = [
                'OPUS plot', 'title', '1', '2', '3', '4',
                '    1.0 5.0',   # bin 1 starts
                '    2.0 5.0',   # same y -> integral += 5 * (2 - 1)
                '    2.0 3.0',   # y changes -> bin 2 starts
                '    4.0 3.0',   # same y -> integral += 3 * (4 - 2)
                'tail line without two numbers stops parsing',
            ]
            fname = self._write_plt(tmp, lines)
            self.assertAlmostEqual(integrate_opus(fname), 5.0 * 1.0 + 3.0 * 2.0)

    def test_integrate_opus_unterminated_point_is_zero(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            for point in ('   1.000e+00   7.000e+00', '1.0 7.0'):  # parsed with and without leading spaces
                fname = self._write_plt(tmp, ['h'] * 6 + [point])
                self.assertEqual(integrate_opus(fname), 0.0)
                fname = self._write_plt(tmp, ['h'] * 6 + [point, '   3.000e+00   7.000e+00'])
                self.assertAlmostEqual(integrate_opus(fname), 14.0)

    def test_integrate_opus_adjacent_equal_bins(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            # bins of 5, 5, and 1 over unit widths; the two 5 bins are separate bins
            lines = ['h'] * 6 + ['   0.0 5.0', '   1.0 5.0', '   1.0 5.0', '   2.0 5.0', '   2.0 1.0', '   3.0 1.0',
                                 '0 years']
            fname = self._write_plt(tmp, lines)
            self.assertAlmostEqual(integrate_opus(fname), 11.0)

    def test_integrate_opus_stops_at_next_spectrum(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            lines = ['h'] * 6 + ['   1.0 2.0', '   3.0 2.0', '2.738e-03 years', '   1.0 9.0', '   3.0 9.0']
            fname = self._write_plt(tmp, lines)
            self.assertAlmostEqual(integrate_opus(fname), 4.0)

    def test_integrate_opus_rejects_mismatched_bin(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            fname = self._write_plt(tmp, ['h'] * 6 + ['   1.0 2.0', '   3.0 4.0'])
            with self.assertRaises(ValueError):
                integrate_opus(fname)


if __name__ == '__main__':
    ensure_dummy_read_opus()
    unittest.main()
