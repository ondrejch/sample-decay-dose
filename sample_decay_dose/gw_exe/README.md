# scalerte SSH gateway executables

This folder contains the SSH-native gateway tools for dispatching scalerte
jobs from `cl` to worker nodes (`c0801,c0802,c0803,c0804`) without opening
network ports.

## Files

- `scalerte_ssh_agent.py`: Run on `cl`, schedules and tracks jobs.
- `scalerte_ssh_client.py`: Run locally, calls the remote agent over SSH.

## Installation

The console scripts are defined in the repository's `pyproject.toml` under
`[project.scripts]`. `pip install ./` installs two commands:

- `scalerte-ssh-agent` runs `sample_decay_dose.gw_exe.scalerte_ssh_agent:main`.
- `scalerte-ssh` runs `sample_decay_dose.gw_exe.scalerte_ssh_client:main`.

The client calls `scalerte-ssh-agent` on the login node by default, so install
the package there as well. Without that install, point the client at the module
with `--agent-command "python3 -m sample_decay_dose.gw_exe.scalerte_ssh_agent"`.

## Typical usage

From local machine:

```bash
scalerte-ssh --host cl submit --workdir /home/you/caseA --cmd "scalerte input.inp"
scalerte-ssh --host cl wait <job_id> --timeout-seconds 86400
```

Global options such as `--host` go before the subcommand. The module form
`python -m sample_decay_dose.gw_exe.scalerte_ssh_client` takes the same
arguments.

The client quotes the agent command and all its arguments into one remote
command string. A `--cmd` value therefore reaches the worker unchanged,
including spaces, quotes, `&&` and pipes. Command tokens after `--` are quoted
one by one.

Give `--workdir` as an absolute path. A relative path is resolved on the login
node against the agent's working directory, which is `$HOME` for an SSH call.

SSH options whose value starts with `-` need the `=` form, for example
`--ssh-option=-oBatchMode=yes`.

From `cl` (agent side):

```bash
scalerte-ssh-agent list
```

Python API import path:

```python
from sample_decay_dose.gw_exe.scalerte_ssh_client import ScalerteSshClient
```

## Job states

A job is `RUNNING` until its run script writes `exit_code`. It then becomes
`COMPLETED` (exit code 0) or `FAILED`. A launch error gives `SUBMIT_FAILED`.

The agent records the worker and the PID of the job's launcher process.
`status`, `result` and `logs --follow` check over SSH that this process still
runs `run.sh` for the job. A job whose launcher has been gone for 60 s
(`LOST_GRACE_SECONDS`) without an exit code becomes `LOST`. This covers node
reboots, the OOM killer and `kill -9`. The grace period allows for NFS caching
of a freshly written `exit_code`. `COMPLETED`, `FAILED`, `SUBMIT_FAILED` and
`LOST` are terminal, so `wait` returns for all of them. `list` and worker
selection read the recorded state without probing.

Job metadata is written to a temporary file and renamed into place, so readers
on the shared filesystem never see a partial `meta.json`.
