import argparse
import subprocess
from pathlib import Path

import pytest

from sample_decay_dose.gw_exe import scalerte_ssh_agent as agent
from sample_decay_dose.gw_exe.scalerte_ssh_client import ScalerteSshClient


def _write_meta(job_dir: Path, worker: str, state: str = "RUNNING") -> None:
    (job_dir / "meta.json").write_text(
        (
            "{\n"
            f'  "job_id": "{job_dir.name}",\n'
            f'  "worker": "{worker}",\n'
            f'  "state": "{state}"\n'
            "}\n"
        ),
        encoding="utf-8",
    )


def test_build_submit_command_from_cmd() -> None:
    args = argparse.Namespace(cmd="scalerte run.inp", command=[])
    assert agent._build_submit_command(args) == "scalerte run.inp"


def test_build_submit_command_from_tokens() -> None:
    args = argparse.Namespace(cmd=None, command=["--", "scalerte", "run file.inp"])
    assert agent._build_submit_command(args) == "scalerte 'run file.inp'"


def test_compute_state_completed(tmp_path: Path) -> None:
    job_dir = tmp_path / "job1"
    job_dir.mkdir()
    _write_meta(job_dir, "c0801")
    (job_dir / "exit_code").write_text("0\n", encoding="utf-8")
    state, rc = agent._compute_state(job_dir)
    assert state == "COMPLETED"
    assert rc == 0


def test_pick_worker_least_busy(tmp_path: Path) -> None:
    jobs_root = tmp_path / "jobs"
    jobs_root.mkdir()

    job_a = jobs_root / "job-a"
    job_a.mkdir()
    _write_meta(job_a, "c0801")

    job_b = jobs_root / "job-b"
    job_b.mkdir()
    _write_meta(job_b, "c0802")
    (job_b / "exit_code").write_text("0\n", encoding="utf-8")

    picked = agent._pick_worker(jobs_root, ["c0801", "c0802", "c0803"])
    assert picked == "c0802"


def test_client_submit_builds_expected_ssh_command(monkeypatch: pytest.MonkeyPatch) -> None:
    recorded: list[list[str]] = []

    def fake_run(cmd, capture_output, text, check):  # noqa: ANN001
        recorded.append(cmd)
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout='{"ok": true, "job_id": "j123", "worker": "c0801"}',
            stderr="",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)

    client = ScalerteSshClient(host="cl", ssh_options=("-o", "BatchMode=yes"))
    result = client.submit(["scalerte", "input.inp"], workdir="/home/u/case")

    assert result["job_id"] == "j123"
    assert recorded, "subprocess.run was not called"
    cmd = recorded[0]
    assert cmd[:5] == ["ssh", "-o", "BatchMode=yes", "cl", "scalerte-ssh-agent"]
    assert "submit" in cmd
    assert "--workdir" in cmd
    assert "--" in cmd
    assert "scalerte" in cmd


def test_client_wait_until_complete(monkeypatch: pytest.MonkeyPatch) -> None:
    client = ScalerteSshClient()
    states = iter(
        [
            {"state": "RUNNING", "job_id": "j1"},
            {"state": "COMPLETED", "job_id": "j1", "exit_code": 0},
        ]
    )

    monkeypatch.setattr(client, "status", lambda job_id: next(states))
    monkeypatch.setattr("time.sleep", lambda _seconds: None)

    final = client.wait("j1", poll_seconds=0.01, timeout_seconds=1.0)
    assert final["state"] == "COMPLETED"
