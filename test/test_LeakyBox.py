import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from leaky_box_origen.LeakyBox import (
    _sorted_steps_by_time,
    _average_rates,
    _time_average_rates,
    _add_analytic_columns,
    _inventory_to_release_rate,
    plot_results,
)


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
        })
        out = _inventory_to_release_rate(df, removal_rate_s=0.25)
        self.assertAlmostEqual(out.loc[0, "xe-135"], 2.5)
        self.assertAlmostEqual(out.loc[1, "xe-135"], 5.0)
        self.assertAlmostEqual(out.loc[1, "time [s]"], 1.0)
        with self.assertRaises(ValueError):
            _inventory_to_release_rate(df, removal_rate_s=0.0)

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
