import json5
import numpy as np
import pytest

from leaky_box_origen import AnalysisRerun as ar
from leaky_box_origen.LeakyBox import PCTperDAY


def _analytic(times, eps_a, eps_b, lam):
    # Independent restatement of the single-isotope A -> B -> C solution (n0 = 1, volume = 1).
    a = np.exp(-(eps_a + lam) * times)
    b = eps_a * (np.exp(-(eps_a + lam) * times) - np.exp(-(eps_b + lam) * times)) / (eps_b - eps_a)
    c = np.exp(-lam * times) * (-eps_b * np.exp(-eps_a * times) + eps_a * (np.exp(-eps_b * times) - 1.0)
                                + eps_b) / (eps_b - eps_a)
    return a, b, c


def _write_box_jsons(root, prefix, isotope, lam):
    times = np.array([43200.0 * k for k in range(1, 6)])
    a, b, c = _analytic(times, PCTperDAY, 0.1 * PCTperDAY, lam or 0.0)
    for box, vals in (("A", a), ("B", b), ("C", c)):
        data = {str(k + 1): {"time": float(t), "adens": {isotope: float(v)}}
                for k, (t, v) in enumerate(zip(times, vals))}
        name = f"box{box}_{prefix}.json5" if prefix else f"box{box}.json5"
        with open(root / name, "w") as f:
            json5.dump(data, f)


@pytest.mark.parametrize("prefix,isotope", [("xe135", "xe-135"), ("xe136", "xe-136")])
def test_rerun_matches_exact_analytic_inputs(tmp_path, prefix, isotope):
    _write_box_jsons(tmp_path, prefix, isotope, ar.TEST_CASES[prefix][1])
    pd_A, pd_B, pd_C = ar.main(run_dir=str(tmp_path), prefix=prefix, write_excel=False)
    for df in (pd_A, pd_B, pd_C):
        assert np.nanmax(np.abs(df["diff [%]"].to_numpy())) < 1e-9


def test_rerun_falls_back_to_unprefixed_names_and_writes_excel(tmp_path):
    _write_box_jsons(tmp_path, None, "xe-135", ar.TEST_CASES["xe135"][1])
    pd_A, _, _ = ar.main(run_dir=str(tmp_path), prefix="xe135")
    assert np.nanmax(np.abs(pd_A["diff [%]"].to_numpy())) < 1e-9
    assert (tmp_path / "leaky_boxes_rerun_xe135.xlsx").is_file()


def test_rerun_missing_inputs_and_unknown_prefix(tmp_path):
    with pytest.raises(FileNotFoundError):
        ar.main(run_dir=str(tmp_path), prefix="xe135", write_excel=False)
    with pytest.raises(ValueError):
        ar.main(run_dir=str(tmp_path), prefix="kr85", write_excel=False)


def test_resolve_run_dir_skips_runs_without_prefix(tmp_path, monkeypatch):
    from pathlib import Path
    base = Path(ar.__file__).resolve().parent
    run_old = base / "run_2000-01-01-residual-test"
    run_new = base / "run_2000-01-02-residual-test"
    try:
        run_old.mkdir(exist_ok=True)
        run_new.mkdir(exist_ok=True)
        _write_box_jsons(run_old, "xe135", "xe-135", ar.TEST_CASES["xe135"][1])
        (run_new / "boxA_xe136.json5").write_text("{}")
        (run_new / "boxB_xe136.json5").write_text("{}")
        (run_new / "boxC_xe136.json5").write_text("{}")
        monkeypatch.chdir(tmp_path)  # cwd holds no box JSONs, so the search falls to run_* dirs
        resolved = ar._resolve_run_dir(None, "xe135")
        assert resolved == run_old
        with pytest.raises(FileNotFoundError):
            ar._resolve_run_dir(None, "kr85")
    finally:
        import shutil
        shutil.rmtree(run_old, ignore_errors=True)
        shutil.rmtree(run_new, ignore_errors=True)
