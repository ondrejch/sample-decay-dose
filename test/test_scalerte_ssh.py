import argparse
import json
import shlex
import subprocess
from pathlib import Path

import pytest

from sample_decay_dose.gw_exe import scalerte_ssh_agent as agent
from sample_decay_dose.gw_exe import scalerte_ssh_client as client_mod
from sample_decay_dose.gw_exe.scalerte_ssh_client import ScalerteSshClient

JOB_ID = "20260927T120000Z-1a2b3c4d"


def _remote_argv(cmd: list[str]) -> list[str]:
    """Simulate OpenSSH: join the words after the host with spaces, then split them as the remote shell does."""
    host_index = cmd.index("--") + 1
    return shlex.split(" ".join(cmd[host_index + 1:]))


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
    assert cmd[:5] == ["ssh", "-o", "BatchMode=yes", "--", "cl"]
    assert len(cmd) == 6, "the agent call must be a single remote command argument"
    argv = _remote_argv(cmd)
    assert argv[:2] == ["scalerte-ssh-agent", "submit"]
    args = agent._build_parser().parse_args(argv[1:])
    assert args.workdir == "/home/u/case"
    assert agent._build_submit_command(args) == "scalerte input.inp"


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


def _recording_client(monkeypatch: pytest.MonkeyPatch, **kwargs) -> tuple[ScalerteSshClient, list[list[str]]]:
    recorded: list[list[str]] = []

    def fake_run(cmd, **_kwargs):  # noqa: ANN001
        recorded.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout='{"ok": true, "job_id": "j1", "job": {}}', stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    return ScalerteSshClient(**kwargs), recorded


@pytest.mark.parametrize(
    "command",
    ["scalerte input.inp", "cd sub && scalerte 'my input.inp'", "echo $HOME; scalerte x.inp | tee log"],
)
def test_cmd_string_survives_openssh_join(monkeypatch: pytest.MonkeyPatch, command: str) -> None:
    client, recorded = _recording_client(monkeypatch, jobs_root="~/jobs dir")
    client.submit(command, workdir="/home/u/case A", name="my job")
    argv = _remote_argv(recorded[0])
    assert argv[0] == "scalerte-ssh-agent"
    args = agent._build_parser().parse_args(argv[1:])
    assert args.cmd == command
    assert args.workdir == "/home/u/case A"
    assert args.jobs_root == "~/jobs dir"
    assert args.name == "my job"


def test_command_tokens_survive_openssh_join(monkeypatch: pytest.MonkeyPatch) -> None:
    client, recorded = _recording_client(monkeypatch)
    client.submit(["scalerte", "run file.inp", "a&&b"], workdir="/w")
    args = agent._build_parser().parse_args(_remote_argv(recorded[0])[1:])
    assert agent._build_submit_command(args) == "scalerte 'run file.inp' 'a&&b'"


def test_multiword_agent_command(monkeypatch: pytest.MonkeyPatch) -> None:
    client, recorded = _recording_client(
        monkeypatch, agent_command="python3 -m sample_decay_dose.gw_exe.scalerte_ssh_agent"
    )
    client.status(JOB_ID)
    argv = _remote_argv(recorded[0])
    assert argv[:3] == ["python3", "-m", "sample_decay_dose.gw_exe.scalerte_ssh_agent"]
    args = agent._build_parser().parse_args(argv[3:])
    assert (args.subcommand, args.job_id) == ("status", JOB_ID)


def test_readme_client_example(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture) -> None:
    recorded: list[list[str]] = []

    def fake_run(cmd, **_kwargs):  # noqa: ANN001
        recorded.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout='{"ok": true, "job_id": "j1"}', stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    rc = client_mod.main(
        ["--host", "cl", "submit", "--workdir", "/home/you/caseA", "--cmd", "scalerte input.inp"]
    )
    assert rc == 0
    assert recorded[0][:3] == ["ssh", "--", "cl"]
    args = agent._build_parser().parse_args(_remote_argv(recorded[0])[1:])
    assert args.cmd == "scalerte input.inp"
    assert json.loads(capsys.readouterr().out)["job_id"] == "j1"


def test_launch_uses_option_terminator_and_null_stdin(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    recorded: list[tuple[list[str], dict]] = []

    def fake_run(cmd, **kwargs):  # noqa: ANN001
        recorded.append((cmd, kwargs))
        return subprocess.CompletedProcess(cmd, 0, stdout="4242\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    pid = agent._launch_remote_job(
        ssh_bin="ssh", ssh_options=["-o", "BatchMode=yes"], worker="c0801", run_script_path=tmp_path / "run.sh"
    )
    assert pid == "4242"
    cmd, kwargs = recorded[0]
    assert cmd[:5] == ["ssh", "-o", "BatchMode=yes", "--", "c0801"]
    assert "</dev/null" in cmd[5]
    assert kwargs["stdin"] is subprocess.DEVNULL


@pytest.mark.parametrize("workers", ["-oProxyCommand=evil", "c0801,-x", "c08 01"])
def test_invalid_worker_names_rejected(workers: str) -> None:
    with pytest.raises(ValueError, match="Invalid worker"):
        agent._parse_workers(workers)


@pytest.mark.parametrize("job_id", ["../etc", "job1", "", JOB_ID + "/.."])
def test_invalid_job_id_rejected(tmp_path: Path, capsys: pytest.CaptureFixture, job_id: str) -> None:
    rc = agent.main(["--jobs-root", str(tmp_path), "status", job_id])
    assert rc == 2
    assert "Invalid job_id" in json.loads(capsys.readouterr().err)["error"]


def _submit_job(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> tuple[Path, str]:
    jobs_root = tmp_path / "jobs"

    def fake_run(cmd, **_kwargs):  # noqa: ANN001
        return subprocess.CompletedProcess(cmd, 0, stdout="4242\n", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.chdir(tmp_path)
    rc = agent.main(
        ["--jobs-root", str(jobs_root), "submit", "--workers", "c0801", "--ssh-option=-oBatchMode=yes",
         "--cmd", "scalerte input.inp"]
    )
    assert rc == 0
    (job_dir,) = list(jobs_root.iterdir())
    return jobs_root, job_dir.name


def test_submit_records_absolute_workdir_and_ssh_options(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    jobs_root, job_id = _submit_job(monkeypatch, tmp_path)
    assert agent.JOB_ID_RE.fullmatch(job_id)
    meta = json.loads((jobs_root / job_id / "meta.json").read_text(encoding="utf-8"))
    assert meta["workdir"] == str(tmp_path)
    assert meta["ssh_options"] == ["-oBatchMode=yes"]
    assert meta["state"] == "RUNNING" and meta["launcher_pid"] == "4242"
    assert not [p for p in (jobs_root / job_id).iterdir() if p.name.endswith(".tmp")]


def _fake_probe(monkeypatch: pytest.MonkeyPatch, returncode: int, stdout: str) -> list[list[str]]:
    recorded: list[list[str]] = []

    def fake_run(cmd, **_kwargs):  # noqa: ANN001
        recorded.append(cmd)
        return subprocess.CompletedProcess(cmd, returncode, stdout=stdout, stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    return recorded


def test_dead_worker_job_becomes_lost_after_grace(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    jobs_root, job_id = _submit_job(monkeypatch, tmp_path)
    job_dir = jobs_root / job_id
    probes = _fake_probe(monkeypatch, 1, "")  # ps finds no such PID
    agent._refresh_liveness(job_dir, grace_seconds=60.0, now=1000.0)
    assert agent._compute_state(job_dir)[0] == "RUNNING"
    assert probes[0][:4] == ["ssh", "-oBatchMode=yes", "--", "c0801"]
    assert probes[0][4] == "ps -ww -o args= -p 4242"
    agent._refresh_liveness(job_dir, grace_seconds=60.0, now=1030.0)
    assert agent._compute_state(job_dir)[0] == "RUNNING"
    agent._refresh_liveness(job_dir, grace_seconds=60.0, now=1061.0)
    assert agent._compute_state(job_dir)[0] == "LOST"
    assert "LOST" in agent.TERMINAL_STATES
    # A LOST job no longer counts against its worker.
    assert agent._active_jobs_by_worker(jobs_root) == {}
    # An exit code that shows up late still wins.
    (job_dir / "exit_code").write_text("0\n", encoding="utf-8")
    assert agent._compute_state(job_dir)[0] == "COMPLETED"


def test_live_or_unknown_worker_job_stays_running(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    jobs_root, job_id = _submit_job(monkeypatch, tmp_path)
    job_dir = jobs_root / job_id
    _fake_probe(monkeypatch, 1, "")
    agent._refresh_liveness(job_dir, grace_seconds=60.0, now=1000.0)
    assert (job_dir / "missing_since").exists()
    _fake_probe(monkeypatch, 0, f"bash {job_dir}/run.sh\n")
    agent._refresh_liveness(job_dir, grace_seconds=0.0)
    assert not (job_dir / "missing_since").exists()
    _fake_probe(monkeypatch, 255, "")  # ssh connection failure: state unknown
    agent._refresh_liveness(job_dir, grace_seconds=0.0)
    assert agent._compute_state(job_dir)[0] == "RUNNING"


def test_reused_pid_counts_as_dead(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    jobs_root, job_id = _submit_job(monkeypatch, tmp_path)
    _fake_probe(monkeypatch, 0, "/usr/sbin/sshd -D\n")
    agent._refresh_liveness(jobs_root / job_id, grace_seconds=0.0)
    assert agent._compute_state(jobs_root / job_id)[0] == "LOST"


def test_status_command_reports_lost(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture
) -> None:
    jobs_root, job_id = _submit_job(monkeypatch, tmp_path)
    capsys.readouterr()
    _fake_probe(monkeypatch, 1, "")
    monkeypatch.setattr(agent, "LOST_GRACE_SECONDS", 0.0)
    assert agent.main(["--jobs-root", str(jobs_root), "status", job_id]) == 0
    job = json.loads(capsys.readouterr().out)["job"]
    assert job["state"] == "LOST" and job["lost_at"]


def test_client_wait_returns_on_lost_and_times_out(monkeypatch: pytest.MonkeyPatch) -> None:
    client = ScalerteSshClient()
    monkeypatch.setattr("time.sleep", lambda _seconds: None)
    monkeypatch.setattr(client, "status", lambda job_id: {"state": "LOST", "job_id": job_id})
    assert client.wait("j1")["state"] == "LOST"
    monkeypatch.setattr(client, "status", lambda job_id: {"state": "RUNNING", "job_id": job_id})
    with pytest.raises(TimeoutError):
        client.wait("j1", poll_seconds=0.0, timeout_seconds=0.0)


def test_meta_write_is_atomic(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    path = tmp_path / "meta.json"
    agent._write_json(path, {"state": "SUBMITTED"})

    def failing_replace(src, dst):  # noqa: ANN001
        raise OSError("disk full")

    monkeypatch.setattr(agent.os, "replace", failing_replace)
    with pytest.raises(OSError):
        agent._write_json(path, {"state": "RUNNING"})
    assert json.loads(path.read_text(encoding="utf-8")) == {"state": "SUBMITTED"}
    assert [p.name for p in tmp_path.iterdir()] == ["meta.json"]
