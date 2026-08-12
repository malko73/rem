"""Pytest root conftest: make src/ and analysis/ importable.

The test suite imports `rem4_numerical` / `quotient_geometry` /
`quotient_optimizer` (src/) and `gamma_F` / `d1_hamiltonians` (analysis/).
Without this file, plain `pytest tests/` fails at collection with
ModuleNotFoundError both locally and in CI (the test job does not set
PYTHONPATH). Root-cause fix: register both directories on sys.path here.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
for _p in (ROOT / "src", ROOT / "analysis"):
    _s = str(_p)
    if _s not in sys.path:
        sys.path.insert(0, _s)
