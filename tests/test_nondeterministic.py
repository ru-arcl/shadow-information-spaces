"""Combinatorial filter: exactness vs. the LP, ground truth, tightness, FOV batching, extensions."""

from __future__ import annotations

import itertools
import math
import random
import time

import pytest

from shadowinfo.events import Appear, Disappear, Enter, Exit, Merge, ShadowSequence, Split
from shadowinfo.lp import ShadowLP, lp_bounds, lp_feasible
from shadowinfo.maxflow import InfeasibleError
from shadowinfo.nondeterministic import (CombinatorialFilter, MultiTeamFilter, RedOrBlueFilter,
                                         counting_filter, evader_status, pursuit_evasion_filter)
from shadowinfo.testing import growing_instance, perturb, random_instance, replay_witness

INF = math.inf


def _subsets(rng, alive, n=3):
    out = []
    for _ in range(n):
        if len(alive) >= 2:
            out.append(rng.sample(alive, rng.randint(2, len(alive))))
    return out


# -- exactness and soundness on random instances ------------------------------------------
def test_coverage_of_random_instances(instances):
    kinds = {type(e).__name__ for inst in instances for e in inst.seq.events}
    assert kinds == {"Appear", "Disappear", "Split", "Merge", "Enter", "Exit"}
    assert sum(inst.total is not None for inst in instances) > 50
    assert sum(any(hi == INF for _, hi in inst.seq.initial.values()) for inst in instances) > 50
    assert sum(any(isinstance(e, Disappear) and e.lo < e.hi for e in inst.seq.events) for inst in instances) > 50
    assert sum(any(isinstance(e, Appear) and e.lo < e.hi for e in inst.seq.events) for inst in instances) > 50


def test_bounds_equal_lp_and_contain_truth(instances):
    for i, inst in enumerate(instances):
        rng = random.Random(i)
        seq, total = inst.seq, inst.total
        f = CombinatorialFilter.from_sequence(seq, total)
        lp = ShadowLP(seq, total)
        ab = f.all_bounds()
        assert set(ab) == set(seq.alive_at_end())
        for s, b in ab.items():
            assert b == lp.extremes([lp.last[s]]), (i, s)
            assert b[0] <= inst.truth_final[s] <= b[1]
        for q in _subsets(rng, list(ab)):
            b = f.bounds(q)
            assert b == lp.extremes([lp.last[s] for s in q]), (i, q)
            assert b[0] <= sum(inst.truth_final[s] for s in q) <= b[1]
        for v, b in f.refine_initial_bounds().items():
            assert b == lp.extremes([lp.first[v]]), (i, v)
            assert b[0] <= inst.truth_left[v] <= b[1]
        tb = f.initial_total_bounds()
        assert tb == lp.extremes([lp.first[s] for s in seq.initial])
        assert tb[0] <= sum(inst.truth_left[s] for s in seq.initial) <= tb[1]
        if len(seq.initial) >= 2:
            q = rng.sample(list(seq.initial), 2)
            assert f.left_bounds(q) == lp.extremes([lp.first[s] for s in q])


def test_batched_fov_equals_naive(instances):
    for inst in instances:
        a = CombinatorialFilter.from_sequence(inst.seq, inst.total, fov="batch")
        b = CombinatorialFilter.from_sequence(inst.seq, inst.total, fov="naive")
        assert a.all_bounds() == b.all_bounds()
        assert a.refine_initial_bounds() == b.refine_initial_bounds()
        assert a.initial_total_bounds() == b.initial_total_bounds()
        assert len(a.istate.left) <= len(b.istate.left)


def test_bounds_are_tight_witnesses_replay(instances):
    checked = 0
    for i, inst in enumerate(instances[:200]):
        rng = random.Random(i)
        f = CombinatorialFilter.from_sequence(inst.seq, inst.total, fov="naive")
        queries = [[s] for s in f.alive()] + _subsets(rng, f.alive(), 1)
        for q in queries:
            lo, hi = f.bounds(q)
            for sense, value in (("min", lo), ("max", hi)):
                if value == INF:
                    with pytest.raises(ValueError):
                        f.witness(q, sense)
                    continue
                w = f.witness(q, sense)
                assert w.value == value
                final = replay_witness(inst.seq, w.supply, w.flow, inst.total)
                assert sum(final[s] for s in q) == value, (i, q, sense)
                checked += 1
    assert checked > 500


def test_incremental_queries_match_lp_on_prefixes(instances):
    for inst in instances[:40]:
        f = CombinatorialFilter(inst.seq.initial, inst.total)
        for k, e in enumerate(inst.seq.events):
            f.apply(e)
            prefix = ShadowSequence(inst.seq.initial, inst.seq.events[:k + 1])
            lp = ShadowLP(prefix, inst.total)
            assert f.all_bounds() == {s: lp.extremes([lp.last[s]]) for s in prefix.alive_at_end()}


def test_paper_recipe_on_random_instances():
    n = 0
    for i in range(300):
        inst = random_instance(random.Random(5000 + i), p_total=0)
        f = CombinatorialFilter.from_sequence(inst.seq)
        for s, (lo, hi) in f.all_bounds().items():
            plo, phi = f.bounds_paper([s])
            assert phi == hi
            assert plo >= lo
            assert f.bounds_paper([s], corrected=True) == (lo, hi)
            n += 1
    assert n > 500


# -- infeasibility -----------------------------------------------------------------------------
@pytest.mark.parametrize("initial,events,total", [
    ({1: (2, 2)}, [Disappear(1, 0, 1)], None),
    ({1: (0, 2)}, [Exit(1, 3)], None),
    ({1: (0, 1), 2: (0, 1)}, [Merge(1, 2, 3), Disappear(3, 3, 3)], None),
    ({1: (0, 2), 2: (0, 2)}, [], (5, 6)),
    ({1: (1, 1)}, [Split(1, 2, 3), Disappear(2, 0, 0), Disappear(3, 0, 0)], None),
    ({1: (0, 5)}, [Exit(1, 2), Enter(1, 1), Disappear(1, 5, 5)], None),
])
def test_inconsistent_observations_raise(initial, events, total):
    seq = ShadowSequence(initial, events)
    assert not lp_feasible(seq, total)
    for fov in ("batch", "naive"):
        f = CombinatorialFilter.from_sequence(seq, total, fov)
        assert not f.feasible()
        with pytest.raises(InfeasibleError):
            f.all_bounds()
        with pytest.raises(InfeasibleError):
            f.refine_initial_bounds()
        if f.alive():
            with pytest.raises(InfeasibleError):
                f.bounds(f.alive())


def test_random_perturbations_feasibility_matches_lp():
    n_infeasible = 0
    for i in range(300):
        rng = random.Random(9000 + i)
        inst = random_instance(rng)
        seq = perturb(rng, inst.seq)
        ok = lp_feasible(seq, inst.total)
        for fov in ("batch", "naive"):
            f = CombinatorialFilter.from_sequence(seq, inst.total, fov)
            assert f.feasible() == ok, i
            if not ok:
                with pytest.raises(InfeasibleError):
                    f.all_bounds()
        if ok:
            lp = ShadowLP(seq, inst.total)
            assert f.all_bounds() == {s: lp.extremes([lp.last[s]]) for s in seq.alive_at_end()}
        n_infeasible += not ok
    assert n_infeasible > 30


def test_query_errors():
    f = CombinatorialFilter({1: (0, 1)})
    f.apply(Split(1, 2, 3))
    with pytest.raises(ValueError):
        f.bounds([1])
    with pytest.raises(ValueError):
        f.bounds([])
    f.apply(Appear(4, 1, 1))
    with pytest.raises(ValueError):
        f.left_bounds([1, 4])
    for bad in ([], [999], [4, 999], [999, 998], [2]):  # 2 is alive, not a left vertex
        with pytest.raises(ValueError):
            f.left_bounds(bad)
    with pytest.raises(ValueError):
        f.refine_initial_bounds([1, 999])
    assert f.refine_initial_bounds([]) == {}
    assert f.left_bounds([4]) == (1, 1)


def test_group_queries_warm_started_from_cached_circulation():
    rng = random.Random(5)
    for _ in range(30):
        inst = random_instance(rng)
        f = CombinatorialFilter.from_sequence(inst.seq, inst.total)
        if not f.feasible() or len(f.alive()) < 2:
            continue
        fresh = CombinatorialFilter.from_sequence(inst.seq, inst.total)
        q = rng.sample(f.alive(), 2)
        n = f._feasible_net(right_group=q)
        assert n.net.is_circulation()
        assert _ints_of(n.net.extremes(n.group)) == f.bounds(q)
        ref = fresh._build(right_group=q)
        ref.net.find_feasible()
        assert f.bounds(q) == _ints_of(ref.net.extremes(ref.group))
        left = list(inst.seq.initial)[:2]
        if len(left) == 2:
            n = f._feasible_net(left_group=left)
            assert n.net.is_circulation()


def _ints_of(b):
    return int(b[0]), (INF if b[1] == INF else int(b[1]))


# -- Sec. V-E helpers -------------------------------------------------------------------------
def test_counting():
    f = counting_filter([1, 2])
    assert f.initial_total_bounds() == (0, INF)
    f.apply(Disappear(1, 3, 3))
    assert f.initial_total_bounds() == (3, INF)
    f.apply(Split(2, 3, 4))
    f.apply(Disappear(3, 1, 1))
    f.apply(Disappear(4, 0, 0))
    assert f.initial_total_bounds() == (4, 4)
    assert f.refine_initial_bounds() == {1: (3, 3), 2: (1, 1)}


def test_pursuit_evasion():
    f = pursuit_evasion_filter([1, 2, 3])
    assert evader_status(f) == {1: "contaminated", 2: "contaminated", 3: "contaminated"}
    f.apply(Disappear(1, 0, 0))
    f.apply(Merge(2, 3, 4))
    f.apply(Split(4, 5, 6))
    f.apply(Appear(7))
    assert evader_status(f) == {5: "contaminated", 6: "contaminated", 7: "clear"}
    f.apply(Disappear(5, 0, 0))
    assert evader_status(f) == {6: "evader", 7: "clear"}
    assert f.bounds([6, 7]) == (1, 1)


def test_evader_status_by_thresholds():
    # An Enter with no earlier Exit breaks the single-evader model but stays feasible.
    f = pursuit_evasion_filter([1, 2])
    f.apply(Enter(2, 1))
    assert f.feasible() and f.all_bounds() == {1: (0, 1), 2: (1, 2)}
    assert evader_status(f) == {1: "contaminated", 2: "evader"}


# -- Sec. V-F distinguishability ------------------------------------------------------------
def test_multi_team_independent():
    init = {"red": {1: (1, 2), 2: (0, 1)}, "blue": {1: (0, 3), 3: (2, 2)}}
    m = MultiTeamFilter(init)
    m.apply(Merge(1, 2, 4))
    m.apply(Disappear(3, 0, 0), {"blue": (2, 2)})
    m.apply(Exit(4, 1), {"red": 1})
    m.apply(Appear(5), None)
    m.apply(Enter(5, 2), {"blue": 2})
    assert m.bounds([4], ["red"]) == (0, 2)
    assert m.bounds([4], ["blue"]) == (0, 3)
    assert m.bounds([4]) == (0, 5)
    assert m.bounds([5]) == (2, 2)
    red = CombinatorialFilter({1: (1, 2), 2: (0, 1), 3: (0, 0)})
    red.extend([Merge(1, 2, 4), Disappear(3, 0, 0), Exit(4, 1), Appear(5)])
    assert m.all_bounds()["red"] == red.all_bounds()
    with pytest.raises(ValueError):
        m.apply(Enter(4))
    with pytest.raises(ValueError):
        m.apply(Disappear(4, 1, 1))


def _red_or_blue_exact_all(pure, mixed, events, queries):
    """Brute force over every red/blue split of every mixed shadow (tiny instances).

    ``queries`` is a list of ``(q, teams)``; returns the exact ``(lo, hi)`` of each.
    """
    ranges = []
    labels = sorted(mixed)
    for s in labels:
        lo, hi = mixed[s]
        ranges.append([(t, r) for t in range(lo, hi + 1) for r in range(t + 1)])
    best = [[INF, -INF] for _ in queries]
    for choice in itertools.product(*ranges):
        init = {t: dict(pure.get(t, {})) for t in ("red", "blue")}
        for s, (t, r) in zip(labels, choice):
            for team, k in (("red", r), ("blue", t - r)):
                plo, phi = init[team].get(s, (0, 0))
                init[team][s] = (plo + k, phi + k)
        m = MultiTeamFilter(init)
        try:
            for e, counts in events:
                m.apply(e, counts)
            for f in m.filters.values():  # every team must be feasible
                f.all_bounds()
        except InfeasibleError:
            continue
        for b, (q, teams) in zip(best, queries):
            lo, hi = m.bounds(q, teams)
            b[0], b[1] = min(b[0], lo), max(b[1], hi)
    return [tuple(b) for b in best]


def _red_or_blue_exact(pure, mixed, events, q, teams=None):
    return _red_or_blue_exact_all(pure, mixed, events, [(q, teams)])[0]


def _red_or_blue(pure, mixed, events):
    f = RedOrBlueFilter(("red", "blue"), pure, mixed)
    for e, counts in events:
        f.apply(e, counts)
    return f


def test_red_or_blue_two_pass():
    pure = {"red": {1: (0, 1)}, "blue": {2: (1, 1)}}
    mixed = {1: (1, 2), 2: (0, 1), 3: (2, 2)}
    events = [(Split(1, 4, 5), None), (Merge(5, 2, 6), None), (Disappear(4, 0, 0), {"red": (1, 1)}),
              (Merge(6, 3, 7), None)]
    f = _red_or_blue(pure, mixed, events)
    r = f.two_pass([7])
    exact = _red_or_blue_exact(pure, mixed, events, [7])
    assert f.bounds([7]) == exact
    assert r.passes["red"] is not None and r.passes["blue"] is not None
    assert exact[0] <= r.lower and r.upper <= exact[1]
    assert r.first_team_lower == (r.passes["blue"]["red"][0], r.passes["red"]["red"][0])
    for team in ("red", "blue"):
        assert f.bounds([7], [team]) == _red_or_blue_exact(pure, mixed, events, [7], [team])


def test_red_or_blue_mixed_group_held_both_colours():
    # Shadow 1 holds 1 red + 1 blue, known only as "2, red or blue"; both go to 2, which
    # disappears showing 1 blue (red not observed).  The two-pass procedure is wrong here.
    mixed = {1: (2, 2)}
    events = [(Split(1, 2, 3), None), (Disappear(2, 0, 0), {"red": (0, INF), "blue": (1, 1)})]
    f = _red_or_blue({}, mixed, events)
    assert f.bounds([3]) == _red_or_blue_exact({}, mixed, events, [3]) == (0, 1)
    for team in ("red", "blue"):
        assert f.bounds([3], [team]) == _red_or_blue_exact({}, mixed, events, [3], [team]) == (0, 1)
    r = f.two_pass([3])
    assert r.passes["red"] is None and (r.lower, r.upper) == (1, 1)  # excludes the truth (0)

    # Both passes contradict the observations, yet the input is consistent.
    pure = {"red": {2: (1, 1)}}
    events = [(Split(1, 3, 4), None), (Disappear(3, 0, 0), {"red": (1, 1), "blue": (1, 1)}),
              (Merge(2, 4, 5), None)]
    f = _red_or_blue(pure, mixed, events)
    assert f.feasible()
    assert f.bounds([5]) == _red_or_blue_exact(pure, mixed, events, [5]) == (1, 1)
    r = f.two_pass([5])
    assert r.passes == {"red": None, "blue": None} and r.lower is None and r.upper is None

    f.apply(Exit(5), {"blue": 1})  # shadow 5 holds no blue target
    assert not f.feasible()
    with pytest.raises(InfeasibleError):
        f.bounds([5])
    with pytest.raises(InfeasibleError):
        f.two_pass([5])
    with pytest.raises(ValueError):
        _red_or_blue({}, mixed, []).bounds([1], ["green"])
    with pytest.raises(ValueError):
        _red_or_blue({}, mixed, []).bounds([])


def _random_red_or_blue(rng):
    """Consistent random red/blue instance with ground truth (component and FOV events)."""
    labels = [1, 2, 3][:rng.randint(2, 3)]
    cnt = {s: [rng.randint(0, 2), rng.randint(0, 2)] for s in labels}
    pure, mixed = {"red": {}, "blue": {}}, {}
    for s in labels:
        c = sum(cnt[s])
        if rng.random() < 0.6:
            mixed[s] = (max(0, c - rng.randint(0, 1)), c + rng.randint(0, 1))
        else:
            pure["red"][s], pure["blue"][s] = (cnt[s][0],) * 2, (cnt[s][1],) * 2
    nxt, events = 10, []
    for _ in range(rng.randint(2, 6)):
        alive = list(cnt)
        k = rng.choice(["split", "merge", "dis", "appear", "fov"] if len(alive) > 1 else ["split", "fov"])
        if k == "split":
            s = rng.choice(alive)
            r, b = cnt.pop(s)
            ra, ba = rng.randint(0, r), rng.randint(0, b)
            cnt[nxt], cnt[nxt + 1] = [ra, ba], [r - ra, b - ba]
            events.append((Split(s, nxt, nxt + 1), None))
            nxt += 2
        elif k == "merge":
            a, b = rng.sample(alive, 2)
            (ra, ba), (rb, bb) = cnt.pop(a), cnt.pop(b)
            cnt[nxt] = [ra + rb, ba + bb]
            events.append((Merge(a, b, nxt), None))
            nxt += 1
        elif k == "appear":
            c = [rng.randint(0, 1), rng.randint(0, 1)]
            cnt[nxt] = c
            events.append((Appear(nxt), {"red": (c[0], c[0]), "blue": (c[1], c[1])}))
            nxt += 1
        elif k == "dis":
            s = rng.choice(alive)
            r, b = cnt.pop(s)
            events.append((Disappear(s), {"red": (r, r) if rng.random() < .7 else (0, INF), "blue": (b, b)}))
        else:
            s, i = rng.choice(alive), rng.randint(0, 1)
            if cnt[s][i] and rng.random() < 0.5:
                cnt[s][i] -= 1
                events.append((Exit(s), {("red", "blue")[i]: 1}))
            else:
                cnt[s][i] += 1
                events.append((Enter(s), {("red", "blue")[i]: 1}))
    return pure, mixed, events, cnt


def test_red_or_blue_exact_and_sound_on_random_instances():
    rng = random.Random(11)
    two_pass_wrong = 0
    for _ in range(300):
        pure, mixed, events, cnt = _random_red_or_blue(rng)
        if not mixed:
            continue
        f = _red_or_blue(pure, mixed, events)
        assert f.feasible()
        alive = f.alive()
        assert sorted(alive) == sorted(cnt)
        qs = [[s] for s in alive] + ([rng.sample(alive, 2)] if len(alive) > 1 else [])
        queries = [(q, teams) for q in qs for teams in (None, ["red"], ["blue"])]
        for (q, teams), exact in zip(queries, _red_or_blue_exact_all(pure, mixed, events, queries)):
            b = f.bounds(q, teams)
            assert b == exact
            truth = sum(c for s in q for t, c in zip(("red", "blue"), cnt[s]) if teams is None or t in teams)
            assert b[0] <= truth <= b[1]
        for q in qs:
            r = f.two_pass(q)
            two_pass_wrong += r.lower is None or (r.lower, r.upper) != f.bounds(q)
    assert two_pass_wrong > 0  # the paper's procedure is not exact on these instances


def test_red_or_blue_without_team_observations_equals_union():
    rng = random.Random(3)
    vacuous = {"red": (0, INF), "blue": (0, INF)}
    for i in range(40):
        inst = random_instance(rng, p_fov=0, p_total=0, p_inf=0)
        events = []
        for e in inst.seq.events:
            if isinstance(e, Appear):
                e = Appear(e.s)
            elif isinstance(e, Disappear):
                e = Disappear(e.s, 0, INF)
            events.append(e)
        seq = ShadowSequence(inst.seq.initial, events)
        union = CombinatorialFilter.from_sequence(seq)
        f = RedOrBlueFilter(("red", "blue"), {}, seq.initial)
        for e in events:
            f.apply(e, vacuous if isinstance(e, Disappear) else None)
        for s in seq.alive_at_end():
            assert f.bounds([s]) == union.bounds([s])
            r = f.two_pass([s])
            assert (r.lower, r.upper) == union.bounds([s])
            assert r.first_team_lower == (0, r.lower)
            assert r.first_team_upper == (0, r.upper)


# -- performance ---------------------------------------------------------------------------------
def test_all_bounds_performance():
    inst = growing_instance(random.Random(7), n_component=500, target_alive=100)
    n_comp = sum(not isinstance(e, (Enter, Exit)) for e in inst.seq.events)
    assert n_comp == 500
    f = CombinatorialFilter.from_sequence(inst.seq)
    t = time.perf_counter()
    ab = f.all_bounds()
    elapsed = time.perf_counter() - t
    assert len(ab) >= 80
    assert elapsed < 2.0, elapsed
    for s, v in inst.truth_final.items():
        assert ab[s][0] <= v <= ab[s][1]
    sample = sorted(ab)[:5]
    assert {s: ab[s] for s in sample} == {s: lp_bounds(inst.seq, [s]) for s in sample}
