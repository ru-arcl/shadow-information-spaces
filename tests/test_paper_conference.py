"""Conference-paper fixtures (ICRA 2008 Fig. 4 and Sec. VI tasks, ICRA 2010 Fig. 5 structure).

Expected values in the fixtures were computed with an independent integer
program (see each fixture's ``notes``); ``hi = null`` means +inf.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from shadowinfo.bipartite import BipartiteIState, vertex_name
from shadowinfo.events import ShadowSequence
from shadowinfo.lp import lp_all_bounds
from shadowinfo.nondeterministic import CombinatorialFilter, evader_status

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def _load(name):
    return json.loads((FIXTURES / f"{name}.json").read_text())


def _b(x):
    return (x[0], math.inf if x[1] is None else x[1])


def _bmap(d):
    return {int(k): _b(v) for k, v in d.items()}


def _edges(seq):
    st = BipartiteIState.from_sequence(seq, fov="naive")
    return sorted((int(vertex_name(u)), int(vertex_name(r[1]))) for u, r in st.edges())


FIG4 = _load("icra08_fig4")
TASKS = _load("icra08_tasks")["cases"]
FIG5 = _load("icra10_fig5_structure")["configs"]


@pytest.mark.parametrize("fov", ["batch", "naive"])
def test_icra08_fig4_bounds(fov):
    seq = ShadowSequence.from_dict(FIG4["sequence"])
    f = CombinatorialFilter.from_sequence(seq, fov=fov)
    want = _bmap(FIG4["expected_bounds"])
    assert f.all_bounds() == want == lp_all_bounds(seq)
    sets = FIG4["expected_bounds_sets"]
    assert f.bounds([13, 15]) == _b(sets["13+15"])
    assert f.bounds(f.alive()) == _b(sets["all_final"])


def test_icra08_fig4_bipartite_edges():
    seq = ShadowSequence.from_dict(FIG4["sequence"])
    assert _edges(seq) == sorted(map(tuple, FIG4["bipartite"]["edges"]))


@pytest.mark.parametrize("case", TASKS, ids=[c["name"] for c in TASKS])
@pytest.mark.parametrize("fov", ["batch", "naive"])
def test_icra08_tasks(case, fov):
    seq = ShadowSequence.from_dict(case["sequence"])
    total = None if case["total"] is None else _b(case["total"])
    f = CombinatorialFilter.from_sequence(seq, total=total, fov=fov)
    final = _bmap(case["expected_final_bounds"])
    assert f.all_bounds() == final == lp_all_bounds(seq, total)
    want = {**_bmap(case["expected_refined_initial"]), **_bmap(case["expected_refined_appear"])}
    assert f.refine_initial_bounds() == want
    assert f.initial_total_bounds() == _b(case["expected_initial_total"])
    if case["name"].startswith("pursuit_evasion"):
        status = {s: "clear" if hi == 0 else "evader" if lo >= 1 else "contaminated"
                  for s, (lo, hi) in final.items()}
        assert evader_status(f) == status


@pytest.mark.parametrize("cfg", FIG5, ids=[c["name"] for c in FIG5])
def test_icra10_fig5_structure(cfg):
    seq = ShadowSequence.from_dict(cfg["sequence"])
    f = CombinatorialFilter.from_sequence(seq)
    want = _bmap(cfg["expected_bounds_without_fov"])
    assert f.all_bounds() == want == lp_all_bounds(seq)
    assert _edges(seq) == sorted(map(tuple, cfg["bipartite"]["edges"]))
