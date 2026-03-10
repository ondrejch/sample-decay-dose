"""Client helpers for invoking ``scalerte-ssh-agent`` over SSH."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Any, Iterable


TERMINAL_STATES = {"COMPLETED", "FAILED", "CANCELED", "SUBMIT_FAILED"}


class ScalerteSshError(RuntimeError):
    """Raised when communication with the remote agent fails."""


@dataclass
class ScalerteSshClient:
    """Thin Python API around ``ssh <host> scalerte-ssh-agent ...``."""

    host: str = "cl"
    agent_command: str = "scalerte-ssh-agent"
    ssh_bin: str = "ssh"
    ssh_options: tuple[str, ...] = ()
    jobs_root: str | None = None

    def _ssh_base_command(self) -> list[str]:
        return [self.ssh_bin, *self.ssh_options, self.host, self.agent_command]

    def _run_json(self, args: list[str]) -> dict[str, Any]:
        cmd = self._ssh_base_command() + args
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise ScalerteSshError(
                "Remote agent call failed: "
                f"rc={proc.returncode}, cmd={cmd!r}, stdout={proc.stdout.strip()!r}, stderr={proc.stderr.strip()!r}"
            )
        try:
            payload = json.loads(proc.stdout)
        except json.JSONDecodeError as exc:
            raise ScalerteSshError(
                f"Remote agent returned non-JSON output for cmd={cmd!r}: {proc.stdout!r}"
            ) from exc
        if not payload.get("ok", False):
            raise ScalerteSshError(f"Remote agent reported failure: {payload!r}")
        return payload

    def _run_text(self, args: list[str], *, capture_output: bool = True) -> str | int:
        cmd = self._ssh_base_command() + args
        if capture_output:
            proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if proc.returncode != 0:
                raise ScalerteSshError(
                    "Remote logs call failed: "
                    f"rc={proc.returncode}, cmd={cmd!r}, stdout={proc.stdout.strip()!r}, stderr={proc.stderr.strip()!r}"
                )
            return proc.stdout
        proc = subprocess.run(cmd, check=False)
        return proc.returncode

    def _agent_args(self, *items: str) -> list[str]:
        args: list[str] = []
        if self.jobs_root:
            args.extend(["--jobs-root", self.jobs_root])
        args.extend(items)
        return args

    def submit(
        self,
        command: str | Iterable[str],
        *,
        workdir: str = ".",
        workers: Iterable[str] | None = None,
        worker: str | None = None,
        name: str | None = None,
    ) -> dict[str, Any]:
        args = self._agent_args("submit", "--workdir", workdir)
        if workers:
            args.extend(["--workers", ",".join(workers)])
        if worker:
            args.extend(["--worker", worker])
        if name:
            args.extend(["--name", name])

        if isinstance(command, str):
            args.extend(["--cmd", command])
        else:
            tokens = list(command)
            if not tokens:
                raise ValueError("command cannot be empty")
            args.append("--")
            args.extend(tokens)
        payload = self._run_json(args)
        return payload["job"] if "job" in payload else payload

    def status(self, job_id: str) -> dict[str, Any]:
        payload = self._run_json(self._agent_args("status", job_id))
        return payload["job"]

    def result(self, job_id: str, *, tail_lines: int = 40) -> dict[str, Any]:
        payload = self._run_json(self._agent_args("result", job_id, "--tail-lines", str(tail_lines)))
        return payload["job"]

    def list(self, *, limit: int = 20) -> list[dict[str, Any]]:
        payload = self._run_json(self._agent_args("list", "--limit", str(limit)))
        return payload["jobs"]

    def logs(
        self,
        job_id: str,
        *,
        stream: str = "stdout",
        lines: int = 100,
        follow: bool = False,
        poll_seconds: float = 0.5,
    ) -> str | int:
        args = self._agent_args(
            "logs",
            job_id,
            "--stream",
            stream,
            "--lines",
            str(lines),
            "--poll-seconds",
            str(poll_seconds),
        )
        if follow:
            args.append("--follow")
            return self._run_text(args, capture_output=False)
        return str(self._run_text(args, capture_output=True))

    def wait(self, job_id: str, *, poll_seconds: float = 5.0, timeout_seconds: float | None = None) -> dict[str, Any]:
        start = time.monotonic()
        while True:
            job = self.status(job_id)
            if job.get("state") in TERMINAL_STATES:
                return job
            if timeout_seconds is not None and (time.monotonic() - start) > timeout_seconds:
                raise TimeoutError(f"Timed out waiting for job {job_id}")
            time.sleep(max(poll_seconds, 0.05))


def _print_json(payload: dict[str, Any]) -> None:
    json.dump(payload, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="scalerte-ssh",
        description="Local SSH client for scalerte-ssh-agent running on a login node.",
    )
    parser.add_argument("--host", default="cl", help="SSH host alias where scalerte-ssh-agent runs.")
    parser.add_argument("--agent-command", default="scalerte-ssh-agent", help="Remote agent command.")
    parser.add_argument("--ssh-bin", default="ssh", help="SSH executable.")
    parser.add_argument(
        "--ssh-option",
        action="append",
        default=[],
        help="Additional option passed to SSH. Repeat for multiple values.",
    )
    parser.add_argument("--jobs-root", default=None, help="Override remote jobs root.")

    sub = parser.add_subparsers(dest="subcommand", required=True)

    submit = sub.add_parser("submit", help="Submit command to remote worker pool.")
    submit.add_argument("--workdir", default=".", help="Working directory on remote side.")
    submit.add_argument("--workers", default=None, help="CSV list of workers.")
    submit.add_argument("--worker", default=None, help="Force specific worker.")
    submit.add_argument("--name", default=None, help="Optional job label.")
    submit.add_argument("--cmd", default=None, help="Command string to run remotely.")
    submit.add_argument("command", nargs=argparse.REMAINDER, help="Command tokens, optionally after '--'.")

    status = sub.add_parser("status", help="Get job status.")
    status.add_argument("job_id")

    result = sub.add_parser("result", help="Get job result and tail.")
    result.add_argument("job_id")
    result.add_argument("--tail-lines", type=int, default=40)

    list_cmd = sub.add_parser("list", help="List recent jobs.")
    list_cmd.add_argument("--limit", type=int, default=20)

    logs = sub.add_parser("logs", help="Print logs for a job.")
    logs.add_argument("job_id")
    logs.add_argument("--stream", choices=("stdout", "stderr"), default="stdout")
    logs.add_argument("--lines", type=int, default=100)
    logs.add_argument("--follow", action="store_true")
    logs.add_argument("--poll-seconds", type=float, default=0.5)

    wait = sub.add_parser("wait", help="Wait until job finishes.")
    wait.add_argument("job_id")
    wait.add_argument("--poll-seconds", type=float, default=5.0)
    wait.add_argument("--timeout-seconds", type=float, default=None)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    client = ScalerteSshClient(
        host=args.host,
        agent_command=args.agent_command,
        ssh_bin=args.ssh_bin,
        ssh_options=tuple(args.ssh_option or []),
        jobs_root=args.jobs_root,
    )

    try:
        if args.subcommand == "submit":
            workers = [w.strip() for w in args.workers.split(",") if w.strip()] if args.workers else None
            if args.cmd and args.command:
                raise ValueError("Use either --cmd or positional command tokens, not both.")
            if args.cmd:
                job = client.submit(
                    args.cmd,
                    workdir=args.workdir,
                    workers=workers,
                    worker=args.worker,
                    name=args.name,
                )
            else:
                tokens = args.command
                if tokens and tokens[0] == "--":
                    tokens = tokens[1:]
                if not tokens:
                    raise ValueError("No command provided.")
                job = client.submit(
                    tokens,
                    workdir=args.workdir,
                    workers=workers,
                    worker=args.worker,
                    name=args.name,
                )
            _print_json(job)
            return 0

        if args.subcommand == "status":
            _print_json(client.status(args.job_id))
            return 0

        if args.subcommand == "result":
            _print_json(client.result(args.job_id, tail_lines=args.tail_lines))
            return 0

        if args.subcommand == "list":
            _print_json({"jobs": client.list(limit=args.limit)})
            return 0

        if args.subcommand == "wait":
            _print_json(
                client.wait(
                    args.job_id,
                    poll_seconds=args.poll_seconds,
                    timeout_seconds=args.timeout_seconds,
                )
            )
            return 0

        if args.subcommand == "logs":
            output = client.logs(
                args.job_id,
                stream=args.stream,
                lines=args.lines,
                follow=args.follow,
                poll_seconds=args.poll_seconds,
            )
            if args.follow:
                return int(output)
            sys.stdout.write(str(output))
            return 0

    except Exception as exc:  # noqa: BLE001
        print(str(exc), file=sys.stderr)
        return 2

    parser.print_help(sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
