import os
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from leaky_box_origen.LeakyBox import (
    LEAKED_ELEMENTS,
    LeakyBox,
    _sorted_steps_by_time,
    _average_rates,
    _time_average_rates,
    _add_analytic_columns,
    _inventory_to_release_rate,
    plot_results,
)


def _make_box_a(index: dict | None = None):
    from leaky_box_origen.LeakyBox import DecayBoxA
    if index is None:
        index = {16: {'case': '20', 'time': '1.0e8'}}
    with patch('leaky_box_origen.LeakyBox.get_f71_positions_index', return_value=index):
        return DecayBoxA('dummy.f71', 500e3)


class TestLeakyBoxHelpers(unittest.TestCase):

    def test_sorted_steps_by_time(self):
        steps = {
            3: {"time": "30"},
            1: {"time": "10"},
            2: {"time": "20"},
        }
        sorted_steps = _sorted_steps_by_time(steps)
        times = [float(v["time"]) for _, v in sorted_steps]
        self.assertEqual(times, [10.0, 20.0, 30.0])

    def test_average_rates_union_and_volume(self):
        prev = {"a": 1.0, "b": 2.0}
        cur = {"b": 4.0, "c": 6.0}
        avg = _average_rates(prev, cur, volume_ratio=2.0)
        self.assertAlmostEqual(avg["a"], 1.0)  # (1 + 0)/2 * 2
        self.assertAlmostEqual(avg["b"], 6.0)  # (2 + 4)/2 * 2
        self.assertAlmostEqual(avg["c"], 6.0)  # (0 + 6)/2 * 2

    def test_time_average_rates(self):
        leak_rates = {
            1: {"time": "0", "rate": {"x": 1.0}},
            2: {"time": "10", "rate": {"x": 3.0}},
            3: {"time": "20", "rate": {"x": 5.0}},
        }
        avg = _time_average_rates(leak_rates, volume_ratio=1.0)
        # Trapezoidal average across two equal intervals -> 3.0
        self.assertAlmostEqual(avg["x"], 3.0)

    def test_add_analytic_columns(self):
        isotope = "xe-136"
        times = np.array([1.0, 2.0])
        eps_a = 0.1
        eps_b = 0.2
        n0_density = 2.0
        n0_total = 10.0

        pd_A = pd.DataFrame({"time [s]": times, isotope: np.nan})
        pd_B = pd.DataFrame({"time [s]": times, "total": np.nan})
        pd_C = pd.DataFrame({"time [s]": times, "total": np.nan})

        # Populate ORIGEN columns with analytic values so diff ~ 0.
        pd_A[isotope] = n0_density * np.exp(-eps_a * times)
        pd_B["total"] = (
            n0_total * eps_a * (np.exp(-eps_a * times) - np.exp(-eps_b * times)) / (eps_b - eps_a)
        )
        pd_C["total"] = (
            n0_total * (-eps_b * np.exp(-eps_a * times) + eps_a * (np.exp(-eps_b * times) - 1.0) + eps_b)
            / (eps_b - eps_a)
        )

        _add_analytic_columns(pd_A, pd_B, pd_C, isotope, n0_density, n0_total, eps_a, eps_b, None)

        self.assertIn("analytic", pd_A.columns)
        self.assertIn("analytic", pd_B.columns)
        self.assertIn("analytic", pd_C.columns)
        self.assertTrue(np.allclose(pd_A["diff [%]"].values, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(pd_B["diff [%]"].values, 0.0, atol=1e-12))
        self.assertTrue(np.allclose(pd_C["diff [%]"].values, 0.0, atol=1e-12))

    def test_add_analytic_columns_equal_rates_is_stable(self):
        isotope = "xe-136"
        times = np.array([0.0, 1.0, 2.0])
        eps = 0.1
        n0_density = 1.0
        n0_total = 1.0

        pd_A = pd.DataFrame({"time [s]": times, isotope: np.ones_like(times)})
        pd_B = pd.DataFrame({"time [s]": times, "total": np.ones_like(times)})
        pd_C = pd.DataFrame({"time [s]": times, "total": np.ones_like(times)})

        _add_analytic_columns(pd_A, pd_B, pd_C, isotope, n0_density, n0_total, eps, eps, None)
        self.assertTrue(np.isfinite(pd_B["analytic"]).all())
        self.assertTrue(np.isfinite(pd_C["analytic"]).all())

    def test_inventory_to_release_rate(self):
        df = pd.DataFrame({
            "time [s]": [0.0, 1.0],
            "time [d]": [0.0, 1.0 / 86400.0],
            "xe-135": [10.0, 20.0],
            "i-131": [4.0, 8.0],
            "cs-137": [100.0, 100.0],  # retained daughter, never removed by ORIGEN
            "rb-88": [7.0, 7.0],  # retained daughter of kr-88
        })
        out = _inventory_to_release_rate(df, removal_rate_s=0.25)
        self.assertAlmostEqual(out.loc[0, "xe-135"], 2.5)
        self.assertAlmostEqual(out.loc[1, "xe-135"], 5.0)
        self.assertAlmostEqual(out.loc[1, "i-131"], 2.0)
        self.assertEqual(out.loc[0, "cs-137"], 0.0)
        self.assertEqual(out.loc[1, "rb-88"], 0.0)
        self.assertAlmostEqual(out.loc[1, "time [s]"], 1.0)
        # Explicit element list overrides the default LEAKED_ELEMENTS.
        out_cs = _inventory_to_release_rate(df, removal_rate_s=0.25, released_elements=['cs'])
        self.assertAlmostEqual(out_cs.loc[0, "cs-137"], 25.0)
        self.assertEqual(out_cs.loc[0, "xe-135"], 0.0)
        with self.assertRaises(ValueError):
            _inventory_to_release_rate(df, removal_rate_s=0.0)

    def test_leaked_elements_match_origen_removal(self):
        box = LeakyBox()
        base = _make_box_a()
        box.setup_cases(base)
        self.assertEqual(tuple(box.nuclide_removal_rates.keys()), LEAKED_ELEMENTS)
        deck = base.origen_deck()
        for ele in LEAKED_ELEMENTS:
            self.assertIn(f'ele=[{ele}]', deck)
        self.assertNotIn('ele=[cs]', deck)

    def test_setup_cases_appends_leak_once(self):
        base = _make_box_a()
        LeakyBox().setup_cases(base)
        LeakyBox().setup_cases(base)
        self.assertEqual(base.case_dir, '_box_A_leak')

    def test_plot_results_creates_file(self):
        isotope = "xe-136"
        times = np.array([1.0, 2.0])
        pd_A = pd.DataFrame({"time [d]": times, isotope: [1.0, 0.5], "analytic": [1.0, 0.5]})
        pd_B = pd.DataFrame({"time [d]": times, "total": [1.0, 0.5], "analytic": [1.0, 0.5]})
        pd_C = pd.DataFrame({"time [d]": times, "total": [0.2, 0.1], "analytic": [0.2, 0.1]})

        with tempfile.TemporaryDirectory() as tmpdir:
            cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                out = plot_results(pd_A, pd_B, pd_C, isotope, out_prefix="testplot", logy=False)
                self.assertTrue(os.path.isfile(out))
            finally:
                os.chdir(cwd)


class TestRunDecaySampleCwdRestore(unittest.TestCase):
    """ A failed SCALE run inside run_decay_sample must not leak the case directory as cwd """

    def _make_box(self):
        from leaky_box_origen.LeakyBox import DecayBoxA
        with patch('leaky_box_origen.LeakyBox.get_f71_positions_index',
                   return_value={16: {'case': '20', 'time': '1.0e8'}}):
            return DecayBoxA('dummy.f71', 500e3)

    def test_cwd_restored_on_scale_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            start_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                box = self._make_box()  # constructed inside tmpdir so box.cwd == tmpdir
                box.atom_dens = {'xe-136': 1.0}
                box.debug = 0
                with patch('leaky_box_origen.LeakyBox.run_scale', return_value=False):
                    with self.assertRaises(RuntimeError):
                        box.run_decay_sample()
                self.assertEqual(os.getcwd(), tmpdir)
            finally:
                os.chdir(start_cwd)

    def test_cwd_restored_on_success(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            start_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                box = self._make_box()
                box.atom_dens = {'xe-136': 1.0}
                box.debug = 0
                with patch('leaky_box_origen.LeakyBox.run_scale', return_value=True), \
                     patch('leaky_box_origen.LeakyBox.get_burned_nuclide_atom_dens',
                           return_value={'xe-136': 0.5}):
                    box.run_decay_sample()
                self.assertEqual(os.getcwd(), tmpdir)
                self.assertEqual(box.final_atom_dens, {'xe-136': 0.5})
            finally:
                os.chdir(start_cwd)

    def test_restore_goes_to_callers_cwd(self):
        # The box is built in one directory and run from another; cwd must return to the caller's.
        with tempfile.TemporaryDirectory() as build_dir, tempfile.TemporaryDirectory() as run_dir:
            start_cwd = os.getcwd()
            try:
                os.chdir(build_dir)
                box = self._make_box()
                box.atom_dens = {'xe-136': 1.0}
                box.debug = 0
                os.chdir(run_dir)
                with patch('leaky_box_origen.LeakyBox.run_scale', return_value=True), \
                     patch('leaky_box_origen.LeakyBox.get_burned_nuclide_atom_dens',
                           return_value={'xe-136': 0.5}):
                    box.run_decay_sample()
                self.assertEqual(os.getcwd(), os.path.realpath(run_dir))
                self.assertTrue(os.path.isdir(os.path.join(run_dir, box.case_dir)))
                self.assertEqual(box.case_path, os.path.join(os.path.realpath(run_dir), box.case_dir))
            finally:
                os.chdir(start_cwd)


class TestDecayBoxBCwdRestore(unittest.TestCase):
    """ Same cwd guarantees for the box B/C decay class """

    @staticmethod
    def _make_box():
        from leaky_box_origen.LeakyBox import DecayBoxB
        box = DecayBoxB({'xe-136': 1.0}, 1.0)
        box.debug = 0
        box.DECAY_steps = 8
        return box

    def test_cwd_restored_on_scale_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            start_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                box = self._make_box()
                with patch('leaky_box_origen.LeakyBox.run_scale', return_value=False):
                    with self.assertRaises(RuntimeError):
                        box.run_decay_sample()
                self.assertEqual(os.getcwd(), os.path.realpath(tmpdir))
            finally:
                os.chdir(start_cwd)

    def test_cwd_restored_on_success_to_callers_cwd(self):
        with tempfile.TemporaryDirectory() as build_dir, tempfile.TemporaryDirectory() as run_dir:
            start_cwd = os.getcwd()
            try:
                os.chdir(build_dir)
                box = self._make_box()
                os.chdir(run_dir)
                with patch('leaky_box_origen.LeakyBox.run_scale', return_value=True), \
                     patch('leaky_box_origen.LeakyBox.get_burned_nuclide_atom_dens',
                           return_value={'xe-136': 0.25}):
                    box.run_decay_sample()
                self.assertEqual(os.getcwd(), os.path.realpath(run_dir))
                self.assertEqual(box.final_atom_dens, {'xe-136': 0.25})
                self.assertTrue(os.path.isfile(os.path.join(run_dir, '_box_B', 'origen.inp')))
            finally:
                os.chdir(start_cwd)

    def test_cwd_restored_when_skip_file_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            start_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                box = self._make_box()
                box.skip_calculation = True
                with self.assertRaises(ValueError):
                    box.run_decay_sample()
                self.assertEqual(os.getcwd(), os.path.realpath(tmpdir))
            finally:
                os.chdir(start_cwd)


def _sequential_parallel(n_jobs=None):
    # Stand-in for joblib.Parallel that runs delayed calls in-process (patches stay active).
    def run(gen):
        return [f(*a, **kw) for f, a, kw in gen]
    return run


class TestLeakRateAbsolutePath(unittest.TestCase):
    """ Workers must get an absolute F71 path, independent of their cwd """

    def test_leak_rate_uses_absolute_f71_path(self):
        index = {1: {'case': '1', 'time': '0.0'}, 2: {'case': '1', 'time': '100.0'}}
        seen_paths = []

        def fake_adens(path, pos):
            seen_paths.append(path)
            return {'xe-135': 2.0 * pos, 'cs-135': 1.0}

        with tempfile.TemporaryDirectory() as run_dir, tempfile.TemporaryDirectory() as other_dir:
            start_cwd = os.getcwd()
            try:
                os.chdir(run_dir)
                base = _make_box_a()
                base.debug = 0
                box = LeakyBox()
                box.removal_rate = 0.5
                box.setup_cases(base)
                with patch('leaky_box_origen.LeakyBox.run_scale', return_value=True), \
                     patch('leaky_box_origen.LeakyBox.get_burned_nuclide_atom_dens', return_value={}):
                    box.run_case()
                os.chdir(other_dir)  # e.g. a reused worker that still sits in another run directory
                with patch('leaky_box_origen.LeakyBox.get_f71_positions_index', return_value=index) as idx, \
                     patch('leaky_box_origen.LeakyBox.get_burned_nuclide_atom_dens', side_effect=fake_adens), \
                     patch('joblib.Parallel', _sequential_parallel):
                    box.get_leak_rate_parallel()
                expected = os.path.join(os.path.realpath(run_dir), '_box_A_leak', base.F71_file_name)
                self.assertEqual(idx.call_args[0][0], expected)
                self.assertEqual(seen_paths, [expected, expected])
                self.assertAlmostEqual(box.leak_rates[2]['rate']['xe-135'], 2.0)  # 4.0 * 0.5
                self.assertNotIn('cs-135', box.leak_rates[2]['rate'])
            finally:
                os.chdir(start_cwd)


class TestSetF71Pos(unittest.TestCase):

    @staticmethod
    def _box():
        index = {
            1: {'case': '20', 'time': '100.0'},
            2: {'case': '20', 'time': '200.0'},
            3: {'case': '20', 'time': '300.0'},
            4: {'case': '21', 'time': '201.0'},
        }
        return _make_box_a(index)

    def test_closest_slot(self):
        box = self._box()
        for t, pos in [(201.0, 2), (250.0, 2), (260.0, 3), (200.0, 2), (100.0, 1), (149.0, 1), (151.0, 2)]:
            box.set_f71_pos(t, case='20')
            self.assertEqual(box.BURNED_MATERIAL_F71_position, pos, msg=f't={t}')

    def test_out_of_range_clamps(self):
        box = self._box()
        box.set_f71_pos(1.0, case='20')
        self.assertEqual(box.BURNED_MATERIAL_F71_position, 1)
        box.set_f71_pos(1e9, case='20')
        self.assertEqual(box.BURNED_MATERIAL_F71_position, 3)

    def test_int_case_accepted(self):
        box = self._box()
        box.set_f71_pos(201.0, case=20)
        self.assertEqual(box.BURNED_MATERIAL_F71_position, 2)
        box.set_f71_pos(201.0, case=21)
        self.assertEqual(box.BURNED_MATERIAL_F71_position, 4)
        with self.assertRaises(ValueError):
            box.set_f71_pos(201.0, case=99)


class TestChiQ(unittest.TestCase):

    def test_scalar_types(self):
        from leaky_box_origen.LeakyBox import _chi_q_at_time
        for val in (2, 2.0, np.int64(2), np.int32(2), np.float32(2.0), np.float64(2.0)):
            self.assertEqual(_chi_q_at_time(10.0, val), 2.0, msg=repr(type(val)))

    def test_schedule_lookup(self):
        from leaky_box_origen.LeakyBox import _chi_q_at_time, chi_q_schedule_rg145_400m
        sched = [(10.0, 3.0), (20.0, 2.0), (30.0, 1.0)]
        self.assertEqual(_chi_q_at_time(0.0, sched), 3.0)
        self.assertEqual(_chi_q_at_time(10.0, sched), 3.0)  # t_end is inclusive
        self.assertEqual(_chi_q_at_time(10.5, sched), 2.0)
        self.assertEqual(_chi_q_at_time(np.float32(25.0), sched), 1.0)
        self.assertEqual(_chi_q_at_time(1e9, sched), 1.0)  # beyond the last entry keeps the last value
        rg = chi_q_schedule_rg145_400m()
        self.assertAlmostEqual(_chi_q_at_time(3600.0, rg), 3.35e-3)
        self.assertAlmostEqual(_chi_q_at_time(4.0 * 3600.0, rg), (1.91e-3 * 8.0 - 3.35e-3 * 2.0) / 6.0)
        self.assertAlmostEqual(_chi_q_at_time(400.0 * 86400.0, rg), 1.60e-4)


class TestComputeDoseTimeseries(unittest.TestCase):

    def test_units_and_integration(self):
        from leaky_box_origen.LeakyBox import compute_dose_timeseries
        df = pd.DataFrame({
            'time [s]': [0.0, 10.0, 30.0],
            'i-131': [1.0e3, 1.0e3, 3.0e3],  # Bq/s
            'xe-133': [1.0e6, 1.0e6, 1.0e6],  # Bq/s
        })
        chi_q = 1.0e-3  # s/m^3
        br = 3.3e-4  # m^3/s
        out = compute_dose_timeseries(df, {'i-131': 7.4e-9}, chi_q, br, {'xe-133': 1.2e-10})
        # Sv/s = s/m^3 * m^3/s * Bq/s * Sv/Bq
        inh0 = chi_q * br * 1.0e3 * 7.4e-9
        imm0 = chi_q * 1.0e6 * 1.2e-10 / 86400.0
        self.assertAlmostEqual(out.loc[0, 'dose_rate_inhalation [Sv/s]'], inh0, delta=1e-12 * inh0)
        self.assertAlmostEqual(out.loc[0, 'dose_rate_immersion [Sv/s]'], imm0, delta=1e-12 * imm0)
        # Trapezoid: inhalation rate r, r, 3r over dt 10, 20 -> 10 r + 20 * 2 r = 50 r
        self.assertAlmostEqual(out['dose_inhalation [Sv]'].iloc[-1] / inh0, 50.0)
        self.assertAlmostEqual(out['dose_immersion [Sv]'].iloc[-1] / imm0, 30.0)
        self.assertAlmostEqual(out['dose [rem]'].iloc[-1], 100.0 * out['dose [Sv]'].iloc[-1])
        self.assertTrue((out['activity_fraction_without_dcf [-]'] == 0.0).all())

    def test_non_range_index_and_unsorted(self):
        from leaky_box_origen.LeakyBox import compute_dose_timeseries
        df = pd.DataFrame({'time [s]': [20.0, 0.0, 10.0], 'i-131': [1.0, 1.0, 1.0]},
                          index=[7, 3, 5])
        out = compute_dose_timeseries(df, {'i-131': 1.0}, 1.0, 1.0)
        self.assertEqual(out['time [s]'].tolist(), [0.0, 10.0, 20.0])
        self.assertEqual(out['dose [Sv]'].tolist(), [0.0, 10.0, 20.0])

    def test_missing_dcf_fraction_and_warning(self):
        from leaky_box_origen.LeakyBox import compute_dose_timeseries
        df = pd.DataFrame({
            'time [s]': [0.0, 10.0],
            'i-131': [3.0, 1.0],
            'kr-85': [0.0, 1.0],  # immersion only, counts as covered
            'h-3': [1.0, 2.0],  # no DCF in either map
            'cs-137': [0.0, 0.0],  # no DCF but no release
        })
        with self.assertWarns(UserWarning) as cm:
            out = compute_dose_timeseries(df, {'i-131': 1e-9}, 1e-3, 3.3e-4, {'kr-85': 2.2e-11},
                                          max_missing_dcf_fraction=None)
        self.assertIn('h-3', str(cm.warning))
        self.assertNotIn('cs-137', str(cm.warning))
        self.assertEqual(out['activity_fraction_without_dcf [-]'].tolist(), [0.25, 0.5])
        with self.assertRaises(ValueError):
            compute_dose_timeseries(df, {'i-131': 1e-9}, 1e-3, 3.3e-4, {'kr-85': 2.2e-11},
                                    max_missing_dcf_fraction=0.4)

    def test_default_missing_fraction_is_fail_closed(self):
        from leaky_box_origen.LeakyBox import (compute_dose_timeseries, DEFAULT_MAX_MISSING_DCF_FRACTION,
                                               DEFAULT_ICRP72_INHALATION_CSV, LEAKY_BOX_DATA_DIR)
        import os
        self.assertEqual(DEFAULT_MAX_MISSING_DCF_FRACTION, 0.05)
        # 25% of the release without a DCF exceeds the default: raises without opt-out.
        df = pd.DataFrame({'time [s]': [0.0, 10.0], 'i-131': [3.0, 1.0], 'h-3': [1.0, 2.0]})
        with self.assertRaises(ValueError):
            compute_dose_timeseries(df, {'i-131': 1e-9}, 1e-3, 3.3e-4)
        # Full coverage passes under the default.
        ok = pd.DataFrame({'time [s]': [0.0, 10.0], 'i-131': [1.0, 1.0]})
        out = compute_dose_timeseries(ok, {'i-131': 1e-9}, 1e-3, 3.3e-4)
        self.assertTrue((out['activity_fraction_without_dcf [-]'] == 0.0).all())
        # Paper default is the conservative max-across-types table, which ships with the package.
        self.assertEqual(DEFAULT_ICRP72_INHALATION_CSV, 'dcf_icrp72_inhalation_adult_max.csv')
        self.assertTrue(os.path.isfile(os.path.join(str(LEAKY_BOX_DATA_DIR),
                                                    DEFAULT_ICRP72_INHALATION_CSV)))


class TestDcfLoaders(unittest.TestCase):

    def _write(self, tmpdir, name, text):
        path = os.path.join(tmpdir, name)
        with open(path, 'w') as f:
            f.write(text)
        return path

    def test_inhalation_loader_caps_and_bounds(self):
        from leaky_box_origen.LeakyBox import _load_dcf_csv
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, 'inh.csv',
                               'nuclide,dcf_sv_bq\n'
                               'Cs-137 ,4.6e-09\n'
                               'pu-238,1.06e-4\n'
                               'zn-69,1.06e-97\n'  # garbled exponent, below lower bound
                               'ar-39,0.013\n'  # above the upper cap
                               'x-1,-1.0\n'
                               'y-1,abc\n')
            dcf = _load_dcf_csv(path)
        self.assertEqual(dcf, {'cs-137': 4.6e-09, 'pu-238': 1.06e-4})

    def test_inhalation_loader_rem_scale(self):
        from leaky_box_origen.LeakyBox import _load_dcf_csv
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, 'inh_rem.csv', 'nuclide,dcf_rem_bq\ncs-137,4.6e-07\nzz-1,1e-19\n')
            dcf = _load_dcf_csv(path)
        self.assertAlmostEqual(dcf['cs-137'], 4.6e-09)
        self.assertNotIn('zz-1', dcf)  # 1e-21 Sv/Bq after scaling

    def test_immersion_loader_caps_and_bounds(self):
        from leaky_box_origen.LeakyBox import _load_dcf_immersion_csv
        with tempfile.TemporaryDirectory() as tmpdir:
            path = self._write(tmpdir, 'imm.csv',
                               'nuclide,dcf_sv_per_bq_m3_day\n'
                               'kr-85,2.2e-11\n'
                               'ar-39,0.013296\n'
                               'kr-83m,1e-30\n')
            dcf = _load_dcf_immersion_csv(path)
        self.assertEqual(dcf, {'kr-85': 2.2e-11})

    def test_missing_columns_raise(self):
        from leaky_box_origen.LeakyBox import _load_dcf_csv, _load_dcf_immersion_csv
        with tempfile.TemporaryDirectory() as tmpdir:
            no_nuc = self._write(tmpdir, 'a.csv', 'name,dcf_sv_bq\ncs-137,1e-9\n')
            no_val = self._write(tmpdir, 'b.csv', 'nuclide,value\ncs-137,1e-9\n')
            with self.assertRaises(ValueError):
                _load_dcf_csv(no_nuc)
            with self.assertRaises(ValueError):
                _load_dcf_csv(no_val)
            with self.assertRaises(ValueError):
                _load_dcf_immersion_csv(no_val)

    def test_missing_fgr_file_hint(self):
        from leaky_box_origen.LeakyBox import _load_dcf_csv
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(FileNotFoundError) as cm:
                _load_dcf_csv(os.path.join(tmpdir, 'dcf_fgr11_missing.csv'))
        msg = str(cm.exception)
        self.assertIn('--pdf', msg)
        self.assertIn('not distributed', msg)
        self.assertNotIn('PDF/', msg)


class TestCaseDirNaming(unittest.TestCase):

    def test_case_dir_names(self):
        from leaky_box_origen.LeakyBox import _case_dir_name, _run_prefix_from_box_json, _step_case_dirs
        self.assertEqual(_case_dir_name('box_A', None), '_box_A')
        self.assertEqual(_case_dir_name('box_A', 'xe135'), '_box_A_xe135')
        self.assertEqual(_case_dir_name('box_B', 'xe135', 3), '_box_B_xe135_0003')
        self.assertEqual(_case_dir_name('box_C', None, 12), '_box_C_0012')
        self.assertEqual(_run_prefix_from_box_json('/x/boxB_xe135.json5'), 'xe135')
        self.assertIsNone(_run_prefix_from_box_json('/x/boxB.json5'))
        self.assertEqual(_step_case_dirs('box_B', 2, 'xe136'),
                         ['_box_B_xe136_0002', '_box_B_xe136_0002_leak', '_box_B_0002', '_box_B_0002_leak'])

    def test_activity_reader_prefers_prefixed_dirs(self):
        import json5
        from leaky_box_origen.LeakyBox import _activity_timeseries_per_nuclide_from_box_json
        with tempfile.TemporaryDirectory() as tmpdir:
            for d in ('_box_B_xe135_0002_leak', '_box_B_xe136_0002_leak', '_box_B_0003_leak'):
                os.makedirs(os.path.join(tmpdir, d))
                open(os.path.join(tmpdir, d, 'origen.f71'), 'w').close()
            box_json = os.path.join(tmpdir, 'boxB_xe135.json5')
            with open(box_json, 'w') as f:
                json5.dump({'2': {'time': 10.0, 'adens': {}}, '3': {'time': 20.0, 'adens': {}}}, f)

            def fake_data(path, pos, f71units='becq'):
                return {'xe-135': float(len(os.path.basename(os.path.dirname(path))))}

            with patch('leaky_box_origen.LeakyBox.get_burned_nuclide_data', side_effect=fake_data) as data:
                df = _activity_timeseries_per_nuclide_from_box_json(box_json, 'box_B')
            used = [os.path.basename(os.path.dirname(c[0][0])) for c in data.call_args_list]
        self.assertEqual(used, ['_box_B_xe135_0002_leak', '_box_B_0003_leak'])  # legacy fallback for step 3
        self.assertEqual(df['time [s]'].tolist(), [10.0, 20.0])


class TestComputeDoseFromBox(unittest.TestCase):

    def test_activity_scale_and_released_elements(self):
        import json5
        from leaky_box_origen.LeakyBox import compute_dose_from_box
        with tempfile.TemporaryDirectory() as tmpdir:
            for k in (1, 2):
                d = os.path.join(tmpdir, f'_box_B_xe135_{k:04d}_leak')
                os.makedirs(d)
                open(os.path.join(d, 'origen.f71'), 'w').close()
            box_json = os.path.join(tmpdir, 'boxB_xe135.json5')
            with open(box_json, 'w') as f:
                json5.dump({'1': {'time': 0.0, 'adens': {}}, '2': {'time': 100.0, 'adens': {}}}, f)
            dcf_csv = os.path.join(tmpdir, 'inh.csv')
            with open(dcf_csv, 'w') as f:
                f.write('nuclide,dcf_sv_bq\ni-131,1e-8\ncs-137,1e-8\n')

            def fake_data(path, pos, f71units='becq'):
                return {'i-131': 1.0e6, 'cs-137': 1.0e9}  # Bq in box B (per cm^3 of box A)

            kwargs = dict(chi_q_schedule=1e-3, breathing_rate_m3_s=3.3e-4, out_prefix='t',
                          activity_representation='inventory_bq', removal_rate_s=1e-5)
            with patch('leaky_box_origen.LeakyBox.get_burned_nuclide_data', side_effect=fake_data), \
                 patch('leaky_box_origen.LeakyBox.plot_dose_timeseries', return_value='x.png'):
                per_cm3 = compute_dose_from_box(box_json, 'box_B', dcf_csv, **kwargs)
                scaled = compute_dose_from_box(box_json, 'box_B', dcf_csv, activity_scale=5.0e5, **kwargs)
            self.assertTrue(os.path.isfile(os.path.join(tmpdir, 'leaky_boxes_dose_t.csv')))
        # Only iodine is released; retained cs-137 contributes nothing.
        expected_rate = 1e-3 * 3.3e-4 * (1.0e6 * 1e-5) * 1e-8
        self.assertAlmostEqual(per_cm3.loc[0, 'dose_rate [Sv/s]'] / expected_rate, 1.0)
        self.assertEqual(per_cm3['activity_scale [-]'].tolist(), [1.0, 1.0])
        self.assertAlmostEqual(scaled.loc[1, 'dose [Sv]'] / per_cm3.loc[1, 'dose [Sv]'], 5.0e5)
        self.assertEqual(scaled['activity_scale [-]'].tolist(), [5.0e5, 5.0e5])


class TestRunSimulationCaseDirs(unittest.TestCase):
    """ Runs with an output prefix must not share ORIGEN case directories """

    def _run(self, out_prefix):
        from leaky_box_origen.LeakyBox import DecayBoxB, _run_simulation
        case_dirs = []

        def fake_run_case(self):
            case_dirs.append(self.decay_leaks.case_dir)
            self.decay_leaks.final_atom_dens = {'xe-135': 0.5}

        def fake_leak(self):
            for i, t in ((1, 0.0), (2, 100.0)):
                self.leak_rates[i] = {'time': t, 'rate': {'xe-135': 1e-6}}
                self.adens[i] = {'time': t, 'adens': {'xe-135': 1.0}}

        def fake_decay_c(self):
            case_dirs.append(self.case_dir)
            self.final_atom_dens = {'xe-135': 0.1}

        with tempfile.TemporaryDirectory() as tmpdir:
            start_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                box_a = _make_box_a()
                with patch.object(LeakyBox, 'run_case', fake_run_case), \
                     patch.object(LeakyBox, 'get_leak_rate_parallel', fake_leak), \
                     patch.object(DecayBoxB, 'run_decay_sample', fake_decay_c):
                    _run_simulation(box_a, False, False, out_prefix)
                written = sorted(os.listdir(tmpdir))
            finally:
                os.chdir(start_cwd)
        return case_dirs, written

    def test_prefixed_case_dirs(self):
        case_dirs, written = self._run('xe135')
        self.assertEqual(case_dirs, ['_box_A_xe135_leak', '_box_B_xe135_0002_leak', '_box_C_xe135_0002'])
        self.assertIn('boxB_xe135.json5', written)

    def test_legacy_case_dirs_without_prefix(self):
        case_dirs, written = self._run(None)
        self.assertEqual(case_dirs, ['_box_A_leak', '_box_B_0002_leak', '_box_C_0002'])
        self.assertIn('boxB.json5', written)
