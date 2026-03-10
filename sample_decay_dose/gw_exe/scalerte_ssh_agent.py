"""SSH-native job agent for dispatching scalerte commands to worker nodes.

This module is intended to run on the login node (for example ``cl``). It
stores job metadata in a shared filesystem directory and executes workload
commands on worker nodes via SSH.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import time
import uuid
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


DEFAULT_WORKERS = ("c0801", "c0802", "c0803", "c0804")
TERMINAL_STATES = {"COMPLETED", "FAILED", "CANCELED", "SUBMIT_FAILED"}


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _emit_json(payload: dict[str, Any], *, stream: Any = sys.stdout) -> None:
    json.dump(payload, stream, indent=2, sort_keys=True)
    stream.write("\n")
    stream.flush()


def _job_dir(jobs_root: Path, job_id: str) -> Path:
    return jobs_root / job_id


def _optional_text(path: Path) -> str | None:
    if not path.exists():
        return None
    return path.read_text(encoding="utf-8", errors="replace").strip()


def _optional_int(path: Path) -> int | None:
    raw = _optional_text(path)
    if raw is None or raw == "":
        return None
    try:
        return int(raw)
    except ValueError:
        return None


def _compute_state(job_directory: Path) -> tuple[str, int | None]:
    if (job_directory / "canceled_at").exists():
        return "CANCELED", None

    exit_code = _optional_int(job_directory / "exit_code")
    if exit_code is not None:
        return ("COMPLETED" if exit_code == 0 else "FAILED"), exit_code

    meta_path = job_directory / "meta.json"
    if meta_path.exists():
        meta = _read_json(meta_path)
        if meta.get("state") == "SUBMIT_FAILED":
            return "SUBMIT_FAILED", None

    return "RUNNING", None


def _tail_lines(path: Path, line_count: int) -> list[str]:
    if line_count <= 0 or not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return list(deque(handle, maxlen=line_count))


def _active_jobs_by_worker(jobs_root: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    if not jobs_root.exists():
        return counts

    for child in jobs_root.iterdir():
        if not child.is_dir():
            continue
        meta_path = child / "meta.json"
        if not meta_path.exists():
            continue
        try:
            meta = _read_json(meta_path)
        except json.JSONDecodeError:
            continue
        worker = meta.get("worker")
        if not worker:
            continue
        state, _ = _compute_state(child)
        if state not in TERMINAL_STATES:
            counts[worker] = counts.get(worker, 0) + 1
    return counts


def _parse_workers(worker_csv: str | None) -> list[str]:
    if not worker_csv:
        return list(DEFAULT_WORKERS)
    workers = [item.strip() for item in worker_csv.split(",") if item.strip()]
    if not workers:
        raise ValueError("No workers configured.")
    return workers


def _pick_worker(jobs_root: Path, workers: list[str]) -> str:
    counts = _active_jobs_by_worker(jobs_root)
    return min(workers, key=lambda worker: (counts.get(worker, 0), worker))


def _build_submit_command(args: argparse.Namespace) -> str:
    if args.cmd and args.command:
        raise ValueError("Use either --cmd or positional command tokens, not both.")

    if args.cmd:
        command = args.cmd.strip()
        if not command:
            raise ValueError("--cmd cannot be empty.")
        return command

    tokens = args.command
    if tokens and tokens[0] == "--":
        tokens = tokens[1:]
    if not tokens:
        raise ValueError("No command provided. Use --cmd or append command tokens after '--'.")
    return " ".join(shlex.quote(token) for token in tokens)


def _build_run_script(job_directory: Path, workdir: Path, command: str) -> str:
    return f"""#!/usr/bin/env bash
set -u

JOB_DIR={shlex.quote(str(job_directory))}
WORKDIR={shlex.quote(str(workdir))}
STDOUT_FILE={shlex.quote(str(job_directory / "stdout.log"))}
STDERR_FILE={shlex.quote(str(job_directory / "stderr.log"))}
EXIT_FILE={shlex.quote(str(job_directory / "exit_code"))}
START_FILE={shlex.quote(str(job_directory / "started_at"))}
END_FILE={shlex.quote(str(job_directory / "finished_at"))}
CMD={shlex.quote(command)}

umask 077
date -Is > "$START_FILE"
echo "${{HOSTNAME:-unknown}}" > "$JOB_DIR/hostname"

cd "$WORKDIR" || {{
  echo "Failed to cd to $WORKDIR" > "$STDERR_FILE"
  echo 200 > "$EXIT_FILE"
  date -Is > "$END_FILE"
  exit 200
}}

bash -lc "$CMD" > "$STDOUT_FILE" 2> "$STDERR_FILE"
rc=$?
echo "$rc" > "$EXIT_FILE"
date -Is > "$END_FILE"
exit "$rc"
"""


def _launch_remote_job(
    *,
    ssh_bin: str,
    ssh_options: list[str],
    worker: str,
    run_script_path: Path,
) -> str:
    remote_command = f"nohup bash {shlex.quote(str(run_script_path))} >/dev/null 2>&1 & echo $!"
    cmd = [ssh_bin, *ssh_options, worker, remote_command]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"SSH launch failed on {worker}: rc={proc.returncode}; "
            f"stdout={proc.stdout.strip()!r}; stderr={proc.stderr.strip()!r}"
        )
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError(f"SSH launch on {worker} returned no PID.")
    return lines[-1]


def _load_job_view(jobs_root: Path, job_id: str) -> dict[str, Any]:
    job_directory = _job_dir(jobs_root, job_id)
    meta_path = job_directory / "meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(f"Unknown job_id: {job_id}")
    meta = _read_json(meta_path)

    state, exit_code = _compute_state(job_directory)
    stdout_path = job_directory / "stdout.log"
    stderr_path = job_directory / "stderr.log"

    view: dict[str, Any] = dict(meta)
    view["state"] = state
    view["exit_code"] = exit_code
    view["started_at"] = _optional_text(job_directory / "started_at")
    view["finished_at"] = _optional_text(job_directory / "finished_at")
    view["hostname"] = _optional_text(job_directory / "hostname")
    view["stdout_path"] = str(stdout_path)
    view["stderr_path"] = str(stderr_path)
    view["stdout_bytes"] = stdout_path.stat().st_size if stdout_path.exists() else 0
    view["stderr_bytes"] = stderr_path.stat().st_size if stderr_path.exists() else 0
    return view


def _stream_log_until_complete(job_directory: Path, log_path: Path, poll_seconds: float) -> None:
    position = 0
    while True:
        if log_path.exists():
            with log_path.open("r", encoding="utf-8", errors="replace") as handle:
                handle.seek(position)
                chunk = handle.read()
                if chunk:
                    sys.stdout.write(chunk)
                    sys.stdout.flush()
                position = handle.tell()
        state, _ = _compute_state(job_directory)
        if state in TERMINAL_STATES:
            if log_path.exists():
                with log_path.open("r", encoding="utf-8", errors="replace") as handle:
                    handle.seek(position)
                    chunk = handle.read()
                    if chunk:
                        sys.stdout.write(chunk)
                        sys.stdout.flush()
            break
        time.sleep(max(poll_seconds, 0.05))


def _command_submit(args: argparse.Namespace) -> int:
    jobs_root = Path(args.jobs_root).expanduser().resolve()
    jobs_root.mkdir(parents=True, exist_ok=True)

    workers = _parse_workers(args.workers)
    command = _build_submit_command(args)
    chosen_worker = args.worker or _pick_worker(jobs_root, workers)
    if chosen_worker not in workers:
        raise ValueError(f"Worker {chosen_worker!r} not in allowed workers: {workers}")

    job_id = f"{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}-{uuid.uuid4().hex[:8]}"
    job_directory = _job_dir(jobs_root, job_id)
    job_directory.mkdir(parents=False, exist_ok=False)
    run_script_path = job_directory / "run.sh"
    run_script_path.write_text(
        _build_run_script(job_directory, Path(args.workdir).expanduser(), command),
        encoding="utf-8",
    )
    os.chmod(run_script_path, 0o700)

    meta = {
        "job_id": job_id,
        "name": args.name,
        "worker": chosen_worker,
        "workers": workers,
        "command": command,
        "workdir": str(Path(args.workdir).expanduser()),
        "jobs_root": str(jobs_root),
        "submitted_at": _utc_now_iso(),
        "state": "SUBMITTED",
    }
    _write_json(job_directory / "meta.json", meta)
    (job_directory / "command.sh").write_text(command + "\n", encoding="utf-8")

    try:
        launcher_pid = _launch_remote_job(
            ssh_bin=args.ssh_bin,
            ssh_options=list(args.ssh_option or []),
            worker=chosen_worker,
            run_script_path=run_script_path,
        )
        meta["launcher_pid"] = launcher_pid
        meta["launched_at"] = _utc_now_iso()
        meta["state"] = "RUNNING"
        _write_json(job_directory / "meta.json", meta)
    except Exception as exc:  # noqa: BLE001
        meta["state"] = "SUBMIT_FAILED"
        meta["submit_error"] = str(exc)
        _write_json(job_directory / "meta.json", meta)
        _emit_json({"ok": False, "job_id": job_id, "error": str(exc)}, stream=sys.stderr)
        return 2

    _emit_json(
        {
            "ok": True,
            "job_id": job_id,
            "worker": chosen_worker,
            "launcher_pid": launcher_pid,
            "job_dir": str(job_directory),
        }
    )
    return 0


def _command_status(args: argparse.Namespace) -> int:
    jobs_root = Path(args.jobs_root).expanduser().resolve()
    try:
        payload = _load_job_view(jobs_root, args.job_id)
    except FileNotFoundError as exc:
        _emit_json({"ok": False, "error": str(exc)}, stream=sys.stderr)
        return 1
    _emit_json({"ok": True, "job": payload})
    return 0


def _command_result(args: argparse.Namespace) -> int:
    jobs_root = Path(args.jobs_root).expanduser().resolve()
    try:
        payload = _load_job_view(jobs_root, args.job_id)
    except FileNotFoundError as exc:
        _emit_json({"ok": False, "error": str(exc)}, stream=sys.stderr)
        return 1

    job_directory = _job_dir(jobs_root, args.job_id)
    payload["stdout_tail"] = _tail_lines(job_directory / "stdout.log", args.tail_lines)
    payload["stderr_tail"] = _tail_lines(job_directory / "stderr.log", args.tail_lines)
    _emit_json({"ok": True, "job": payload})
    return 0


def _command_list(args: argparse.Namespace) -> int:
    jobs_root = Path(args.jobs_root).expanduser().resolve()
    jobs: list[dict[str, Any]] = []
    if jobs_root.exists():
        for child in sorted(jobs_root.iterdir(), key=lambda path: path.name, reverse=True):
            if not child.is_dir():
                continue
            meta_path = child / "meta.json"
            if not meta_path.exists():
                continue
            try:
                jobs.append(_load_job_view(jobs_root, child.name))
            except Exception:  # noqa: BLE001
                continue
            if len(jobs) >= args.limit:
                break
    _emit_json({"ok": True, "jobs": jobs})
    return 0


def _command_logs(args: argparse.Namespace) -> int:
    jobs_root = Path(args.jobs_root).expanduser().resolve()
    job_directory = _job_dir(jobs_root, args.job_id)
    if not job_directory.exists():
        print(f"Unknown job_id: {args.job_id}", file=sys.stderr)
        return 1

    log_name = "stdout.log" if args.stream == "stdout" else "stderr.log"
    log_path = job_directory / log_name

    if not args.follow:
        for line in _tail_lines(log_path, args.lines):
            sys.stdout.write(line)
        sys.stdout.flush()
        return 0

    for line in _tail_lines(log_path, args.lines):
        sys.stdout.write(line)
    sys.stdout.flush()
    _stream_log_until_complete(job_directory, log_path, poll_seconds=args.poll_seconds)
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scalerte-ssh-agent",
        description=(
            "Run scalerte-related commands on worker nodes via SSH and track job state "
            "in a shared filesystem."
        ),
    )
    parser.add_argument(
        "--jobs-root",
        default="~/.scalerte_ssh/jobs",
        help="Directory for job metadata and logs (default: ~/.scalerte_ssh/jobs).",
    )

    sub = parser.add_subparsers(dest="subcommand", required=True)

    submit = sub.add_parser("submit", help="Submit a command to a worker node.")
    submit.add_argument("--workdir", default=".", help="Working directory for the command.")
    submit.add_argument("--workers", default=",".join(DEFAULT_WORKERS), help="CSV list of worker hostnames.")
    submit.add_argument("--worker", default=None, help="Force a specific worker.")
    submit.add_argument("--name", default="", help="Optional human-readable label.")
    submit.add_argument("--cmd", default=None, help="Command string run with bash -lc.")
    submit.add_argument("--ssh-bin", default="ssh", help="SSH executable (default: ssh).")
    submit.add_argument(
        "--ssh-option",
        action="append",
        default=[],
        help="Additional option passed to SSH. Repeat for multiple values.",
    )
    submit.add_argument("command", nargs=argparse.REMAINDER, help="Command tokens, optionally after '--'.")
    submit.set_defaults(handler=_command_submit)

    status = sub.add_parser("status", help="Get job status.")
    status.add_argument("job_id", help="Job ID from submit output.")
    status.set_defaults(handler=_command_status)

    result = sub.add_parser("result", help="Get status and output tail.")
    result.add_argument("job_id", help="Job ID from submit output.")
    result.add_argument("--tail-lines", type=int, default=40, help="Number of tail lines for stdout/stderr.")
    result.set_defaults(handler=_command_result)

    list_cmd = sub.add_parser("list", help="List recent jobs.")
    list_cmd.add_argument("--limit", type=int, default=20, help="Max jobs to return.")
    list_cmd.set_defaults(handler=_command_list)

    logs = sub.add_parser("logs", help="Print job logs.")
    logs.add_argument("job_id", help="Job ID from submit output.")
    logs.add_argument("--stream", choices=("stdout", "stderr"), default="stdout", help="Which log to print.")
    logs.add_argument("--lines", type=int, default=100, help="Initial tail lines.")
    logs.add_argument("--follow", action="store_true", help="Follow log stream until job completes.")
    logs.add_argument(
        "--poll-seconds",
        type=float,
        default=0.5,
        help="Polling interval while following logs.",
    )
    logs.set_defaults(handler=_command_logs)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except Exception as exc:  # noqa: BLE001
        _emit_json({"ok": False, "error": str(exc)}, stream=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
