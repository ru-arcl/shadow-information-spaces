"""Shared pytest fixtures: seeded random instances (generators live in :mod:`shadowinfo.testing`)."""

from __future__ import annotations

import json
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from shadowinfo.events import ShadowSequence  # noqa: E402
from shadowinfo.testing import random_instance  # noqa: E402

N_RANDOM = 320


@pytest.fixture(scope="session")
def instances():
    """``N_RANDOM`` consistent random instances (FOV events, ranges, inf bounds, totals)."""
    return [random_instance(random.Random(1000 + i)) for i in range(N_RANDOM)]


@pytest.fixture(scope="session")
def fig11():
    d = json.loads((ROOT / "tests" / "fixtures" / "tor_fig11.json").read_text())
    d["seq"] = ShadowSequence.from_dict(d["sequence"])
    return d
