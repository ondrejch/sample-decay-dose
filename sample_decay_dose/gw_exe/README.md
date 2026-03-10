# scalerte SSH gateway executables

This folder contains the SSH-native gateway tools for dispatching scalerte
jobs from `cl` to worker nodes (`c0801,c0802,c0803,c0804`) without opening
network ports.

## Files

- `scalerte_ssh_agent.py`: Run on `cl`, schedules and tracks jobs.
- `scalerte_ssh_client.py`: Run locally, calls the remote agent over SSH.
- `setup.cfg`: Console-script entry point definitions for these executables.

## Typical usage

From local machine:

```bash
python -m sample_decay_dose.gw_exe.scalerte_ssh_client submit --host cl --workdir /home/you/caseA --cmd "scalerte input.inp"
```

From `cl` (agent side):

```bash
python -m sample_decay_dose.gw_exe.scalerte_ssh_agent list
```

Python API import path:

```python
from sample_decay_dose.gw_exe.scalerte_ssh_client import ScalerteSshClient
```
