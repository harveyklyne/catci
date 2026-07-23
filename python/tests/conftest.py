"""Shared fixtures: load the frozen R oracle from ``<repo>/fixtures``.

See ``fixtures/README.md`` for the JSON conventions (row-major matrices, 1-based
indices, all arrays are JSON arrays). Reconstruction helpers live here so every
test parses the oracle the same way.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

# python/tests/conftest.py -> repo root is parents[2]
FIXTURE_PATH = Path(__file__).resolve().parents[2] / "fixtures" / "catci_fixtures.json"


def mat(m: dict) -> np.ndarray:
    """Reconstruct a row-major ``{nrow, ncol, data}`` matrix."""
    return np.array(m["data"], dtype=float).reshape(m["nrow"], m["ncol"])


def vec(v) -> np.ndarray:
    return np.asarray(v, dtype=float)


@pytest.fixture(scope="session")
def oracle() -> dict:
    if not FIXTURE_PATH.exists():
        pytest.skip(f"fixtures not found at {FIXTURE_PATH}")
    with open(FIXTURE_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="session")
def shared_TS(oracle) -> tuple[np.ndarray, np.ndarray]:
    """The (T, Sigma) shared by the search / scalar / update fixtures."""
    sp = oracle["search_paths"]
    return vec(sp["T_vector"]), mat(sp["Sigma"])
