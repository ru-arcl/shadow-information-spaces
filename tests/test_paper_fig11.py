"""Fig. 11 max-flow example (T-RO 2012, Sec. V-C): s19 has between 10 and 24 targets."""

from __future__ import annotations

from collections import Counter

from shadowinfo.lp import lp_all_bounds, lp_bounds
from shadowinfo.nondeterministic import CombinatorialFilter
from shadowinfo.testing import replay_witness


def test_bounds_match_paper_and_lp(fig11):
    seq = fig11["seq"]
    q = fig11["query"]
    assert lp_bounds(seq, q) == (fig11["lp_lower"], fig11["lp_upper"])
    for fov in ("batch", "naive"):
        f = CombinatorialFilter.from_sequence(seq, fov=fov)
        assert f.bounds(q) == (fig11["paper_lower"], fig11["paper_upper"]) == (10, 24)
        expect = {int(s): tuple(b) for s, b in fig11["lp_all_bounds"].items()}
        assert f.all_bounds() == expect == lp_all_bounds(seq)


def test_literal_paper_recipe_reproduces_fig11(fig11):
    f = CombinatorialFilter.from_sequence(fig11["seq"])
    assert f.bounds_paper(fig11["query"]) == (fig11["paper_lower"], fig11["paper_upper"])
    assert f.bounds_paper(fig11["query"], corrected=True) == (10, 24)


def test_witnesses_replay(fig11):
    seq = fig11["seq"]
    f = CombinatorialFilter.from_sequence(seq, fov="naive")
    for sense, value in (("max", 24), ("min", 10)):
        w = f.witness([19], sense)
        assert w.value == value
        assert replay_witness(seq, w.supply, w.flow)[19] == value


def test_fixture_semantic_witnesses_are_feasible(fig11):
    seq = fig11["seq"]
    gone = {9: (2, 3), 14: (2, 4)}
    for key, value in (("semantic_witness_upper", 24), ("semantic_witness_lower", 10)):
        w = fig11[key]
        supply = {int(k): v for k, v in w["supply"].items()}
        flow = {(u, ("gone", v) if v in gone else ("alive", v)): k for u, v, k in w["flow"]}
        assert replay_witness(seq, supply, flow)[19] == value
        into = Counter()
        for u, v, k in w["flow"]:
            into[v] += k
        assert all(lo <= into[v] <= hi for v, (lo, hi) in gone.items())


def test_refined_initial_bounds_are_already_tight(fig11):
    f = CombinatorialFilter.from_sequence(fig11["seq"])
    refined = f.refine_initial_bounds()
    initial = dict(fig11["seq"].initial)
    assert {s: refined[s] for s in initial} == initial
