# Tests

Unit tests for the core library and leaky-box utilities.

**Run** (from the repository root)
- `pytest -q`

Nine of the modules hold pytest-style test functions, which `python -m unittest discover` does not collect,
so the unittest runner covers only the `TestCase` modules.
