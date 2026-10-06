"""Probabilistic filter (T-RO 2012 Sec. VI): paper regressions, brute force,
Monte Carlo agreement and truncation heuristics."""

from __future__ import annotations

import itertools
import json
import math
import random
from fractions import Fraction
from pathlib import Path

import pytest

from shadowinfo.events import Appear, Disappear, Enter, Exit, Merge, ShadowSequence, Split
from shadowinfo.probabilistic import (
    BinomialSplit, ExactFilter, FovObservation, InconsistentObservationError, ProbAppear,
    ProbDisappear, ProbSequence, ProbSplit, TruncatedFilter, TruncationFailure, monte_carlo,
    observation_from_dict, observation_to_dict, split_rule_from_dict, symmetric_obs_model,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def load(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text())


def table(joint) -> dict:
    return {tuple(k): p for k, p in joint}


def assert_table(f: ExactFilter, labels, expected: dict, tol: float = 1e-9) -> None:
    assert f.labels == tuple(labels)
    assert set(f.table) == set(expected)
    for k, p in expected.items():
        assert f.table[k] == pytest.approx(p, abs=tol), k


def assert_normalized(f) -> None:
    assert math.isclose(f.total(), 1.0, abs_tol=1e-9)
    assert all(p > 0 for p in f.table.values())


@pytest.fixture(scope="module")
def table3():
    fx = load("tor_fig12_table3.json")
    seq = ProbSequence.from_dict(fx["prob_sequence"])
    kw = dict(split_rule=split_rule_from_dict(fx["split_rule"]), obs_model=fx["obs_model"])
    return fx, seq, kw


# -- Table III / Fig. 12 (Sec. VI-D) ---------------------------------------------

def test_table3_every_step(table3):
    fx, seq, kw = table3
    f = ExactFilter.from_sequence(seq, **kw)
    steps = fx["steps"]
    assert_table(f, steps[0]["labels"], table([e[:2] for e in steps[0]["joint_paper_order"]]))
    for t, o in enumerate(seq.observations):
        f.apply(o)
        step = steps[t + 1]
        assert step["after"] == t
        assert_table(f, step["labels"], table([e[:2] for e in step["joint_paper_order"]]))
        for k, p, frac in step["joint_paper_order"]:
            assert f.table[tuple(k)] == pytest.approx(float(Fraction(frac)), abs=1e-12)
        assert_normalized(f)
        if t + 1 < len(seq.observations) and isinstance(seq.observations[t + 1], Merge):
            assert f.entries == fx["entries_before_merge"] == 10


def test_table3_final_matches_paper(table3):
    fx, seq, kw = table3
    f = ExactFilter.from_sequence(seq, **kw).run(seq.observations)
    m = f.marginal(4)
    for x, p in fx["expected_final"]["paper_rounded"].items():
        assert m[int(x)] == pytest.approx(p, abs=1e-4)
    wrong = fx["expected_final_if_global_renormalize"]["marginal_s4"]
    assert abs(m[2] - wrong["2"]) > 0.01
    assert f.expected_counts() == pytest.approx({4: 1 * 2 / 13 + 2 * 10 / 13})


def test_table3_exact_rationals():
    half = Fraction(1, 2)
    p = Fraction(9, 10)
    model = symmetric_obs_model(p)
    obs = [FovObservation(1, "exit"), Split(2, 3, 4), FovObservation(3, "enter"), Merge(1, 3, 5),
           ProbDisappear.of(5, {1: half, 2: half})]
    f = ExactFilter.from_counts({1: 2, 2: 2}, split_rule=BinomialSplit(half), obs_model=model)
    f.run(obs)
    assert f.marginal(4) == {0: Fraction(1, 13), 1: Fraction(2, 13), 2: Fraction(10, 13)}


def test_icra10_fixture_agrees():
    fx = load("icra10_fig4_table3.json")
    seq = ProbSequence.from_dict(fx["prob_sequence"])
    f = ExactFilter.from_sequence(seq, split_rule=split_rule_from_dict(fx["split_rule"]),
                                  obs_model=fx["obs_model"])
    for t, o in enumerate(seq.observations):
        f.apply(o)
        step = next((s for s in fx["steps"] if s["after_event"] == t), None)
        if step is not None:
            assert_table(f, step["labels"], table(step["joint"]))
    for x, p in fx["expected_final_marginal"]["p"].items():
        assert f.marginal(4)[int(x)] == pytest.approx(p, abs=1e-9)


def test_table1_merge():
    fx = load("icra10_table1_merge.json")
    seq = ProbSequence.from_dict(fx["prob_sequence"])
    f = ExactFilter.from_sequence(seq).run(seq.observations)
    assert_table(f, fx["expected_final"]["labels"], table(fx["expected_final"]["joint"]))


def test_table3_truncating_merge_output_loses_s4_zero(table3):
    """Sec. VI-E: dropping P(s4=0, s5=2) = 0.0225 loses P(s4=0) = 0.0769."""
    fx, seq, kw = table3
    f = ExactFilter.from_sequence(seq, **kw)
    for o in seq.observations[:3]:
        f.apply(o)
    g = TruncatedFilter(f.labels, f.table, max_entries=6, mode="TR", **kw)
    g.run(seq.observations[3:])
    assert g.truncations == 1
    expected = fx["truncation_example"]["expected_final_marginal_s4"]
    m = g.marginal(4)
    for x, p in expected.items():
        assert m.get(int(x), 0.0) == pytest.approx(p, abs=1e-12)


def test_table3_monte_carlo(table3):
    fx, seq, kw = table3
    exact = ExactFilter.from_sequence(seq, **kw).run(seq.observations).marginal(4)
    n = 10000
    mc = monte_carlo(seq, trials=n, seed=7, **kw)
    assert mc.successes == n and mc.attempts > n
    assert mc.labels == (4,)
    assert_normalized(mc)
    for x, p in exact.items():
        assert abs(mc.marginal(4).get(x, 0) - p) < 4 * math.sqrt(p * (1 - p) / n) + 1e-3
    wrong = fx["expected_final_if_global_renormalize"]["marginal_s4"]
    assert abs(mc.marginal(4)[2] - exact[2]) < abs(mc.marginal(4)[2] - wrong["2"])
    paper = fx["monte_carlo_paper"]
    small = monte_carlo(seq, trials=paper["trials_successful"], seed=1, **kw).marginal(4)
    for x, p in paper["marginal_s4"].items():
        assert abs(small.get(int(x), 0) - p) < 0.04


def test_monte_carlo_deterministic(table3):
    _, seq, kw = table3
    a = monte_carlo(seq, trials=300, seed=11, **kw)
    b = monte_carlo(seq, trials=300, seed=11, **kw)
    assert a.table == b.table and a.attempts == b.attempts


# -- the longer Fig. 12 sequence (Sec. VI-E, Fig. 14) -------------------------------

@pytest.mark.parametrize("variant", ["derived_figure_order", "derived_s4_obs_before_merge"])
def test_fig12_complex_sequence(variant):
    fx = load("tor_fig12_complex.json")
    seq = ProbSequence.from_dict(fx["prob_sequence"])
    v = fx[variant]
    order = v.get("order", range(len(seq.observations)))
    obs = [seq.observations[i] for i in order]
    f = ExactFilter.from_sequence(seq, split_rule=split_rule_from_dict(fx["split_rule"]),
                                  obs_model=fx["obs_model"])
    for t, o in enumerate(obs):
        if isinstance(o, Merge):
            bm = v["before_merge"]
            assert f.labels == tuple(bm["labels"])
            assert f.entries == bm["num_entries"]
            for s, e in bm["expected"].items():
                assert f.expected_counts()[int(s)] == pytest.approx(e, abs=1e-9)
        f.apply(o)
        assert_normalized(f)
    end = v["end"]
    assert f.labels == tuple(end["labels"]) and f.entries == end["num_entries"]
    for s, e in end["expected_exact"].items():
        assert f.expected_counts()[int(s)] == pytest.approx(float(Fraction(e)), abs=1e-9)
    for x, p in end["marginal_s4"].items():
        assert f.marginal(4)[int(x)] == pytest.approx(p, abs=1e-12)
    for x, p in end["marginal_s5"].items():
        assert f.marginal(5)[int(x)] == pytest.approx(p, abs=1e-12)


def test_fig12_complex_paper_entry_count():
    """The paper's 135 entries before the merge (s4 observations first)."""
    fx = load("tor_fig12_complex.json")
    derived = fx["derived_s4_obs_before_merge"]["before_merge"]["num_entries"]
    assert fx["paper"]["entries_before_merge"] == derived == 135


# -- Fig. 16 (Sec. VII-B), regression without the unpublished FOV events -----------

@pytest.fixture(scope="module")
def fig16():
    fx = load("tor_fig16_sequence.json")["table_iv"]
    return fx, ProbSequence.from_dict(fx["prob_sequence"])


def test_fig16_exact_regression(fig16):
    fx, seq = fig16
    ShadowSequence({s: (n, n) for s, n in zip(seq.labels, next(iter(seq.joint)))}, seq.observations)
    f = ExactFilter.from_sequence(seq).run(seq.observations)
    d = fx["derived_no_fov_exact"]
    assert f.labels == (17, 18, 19)
    assert f.peak_entries == d["peak_entries"] and f.entries == d["final_entries"]
    assert f.expected_counts() == pytest.approx({int(s): e for s, e in d["expected"].items()}, abs=1e-9)
    assert sum(f.expected_counts().values()) == pytest.approx(10 + 7 + 8 + 9 - 6 - 9 - 4)


@pytest.mark.parametrize("mode", ["TR", "RT", "RT-LA"])
def test_truncation_large_budget_is_exact(fig16, mode):
    _, seq = fig16
    exact = ExactFilter.from_sequence(seq).run(seq.observations)
    g = TruncatedFilter(seq.labels, seq.joint, exact.peak_entries, mode=mode, seed=5).run(seq.observations)
    assert g.truncations == 0 and g.table == exact.table


@pytest.mark.parametrize("mode", ["TR", "RT"])
def test_truncation_small_budget(fig16, mode):
    _, seq = fig16
    exact = ExactFilter.from_sequence(seq).run(seq.observations).expected_counts()
    with pytest.raises(TruncationFailure):
        TruncatedFilter(seq.labels, seq.joint, 2000, mode=mode, seed=5).run(seq.observations)
    g = TruncatedFilter(seq.labels, seq.joint, 20000, mode=mode, seed=5)
    for t, o in enumerate(seq.observations):
        g.apply(o, seq.observations[t + 1:])
        assert_normalized(g)
        assert g.entries <= 20000
    assert g.truncations >= 1
    for s, e in exact.items():
        assert g.expected_counts()[s] == pytest.approx(e, abs=0.2)


def _noisy_sequence():
    F = FovObservation
    return [Split(1, 3, 4), F(3, "enter"), F(4, "exit"), Split(2, 5, 6), F(5, "null"),
            F(6, "enter"), F(3, "exit"), F(4, "enter"), F(6, "exit"), F(5, "enter"), F(3, "null"),
            Merge(3, 5, 7), F(7, "exit"), F(6, "null"), F(4, "exit"), F(7, "enter"),
            F(6, "enter"), Disappear(4, 2, 3)]


@pytest.mark.parametrize("mode", ["RT", "RT-LA"])
def test_randomized_truncation_is_seeded(mode):
    kw = dict(obs_model=symmetric_obs_model(0.8, null=True))
    obs = _noisy_sequence()
    runs = [TruncatedFilter((1, 2), {(6, 6): 1.0}, 100, mode=mode, seed=s, **kw).run(obs)
            for s in (21, 21, 22)]
    assert runs[0].table == runs[1].table
    assert runs[0].truncations > 0
    assert runs[0].table != runs[2].table
    exact = ExactFilter.from_counts({1: 6, 2: 6}, **kw).run(obs)
    assert exact.peak_entries > 100
    for s, e in exact.expected_counts().items():
        assert runs[0].expected_counts()[s] == pytest.approx(e, abs=0.5)


def test_rt_la_truncates_less_than_rt():
    kw = dict(obs_model=symmetric_obs_model(0.8, null=True))
    obs = _noisy_sequence()
    rt = TruncatedFilter((1, 2), {(6, 6): 1.0}, 100, mode="RT", seed=1, **kw).run(obs)
    la = TruncatedFilter((1, 2), {(6, 6): 1.0}, 100, mode="RT-LA", seed=1, **kw).run(obs)
    assert 0 < la.truncations < rt.truncations


# -- RT-LA lookahead ------------------------------------------------------------

def _wide(n: int = 3):
    return ExactFilter.from_independent({1: {x: 1 for x in range(n)}, 2: {x: 1 for x in range(n)}})


@pytest.mark.parametrize("gap,kind,counts_fov,truncated", [
    (3, "disappear", True, False), (4, "disappear", True, True),
    (1, "merge", True, False), (2, "merge", True, True),
    (4, "disappear", False, False), (2, "merge", False, False),
])
def test_rt_la_lookahead(gap, kind, counts_fov, truncated):
    base = _wide()
    f = TruncatedFilter(base.labels, base.table, 4, mode="RT-LA", seed=3,
                        lookahead_counts_fov=counts_fov)
    future = [Enter(1)] * gap + [Disappear(1, 0, 10) if kind == "disappear" else Merge(1, 2, 3)]
    f.apply(Enter(2), future)
    assert (f.truncations == 1) is truncated
    assert f.entries == (4 if truncated else 9)


def test_rt_la_without_future_is_rt():
    base = _wide()
    a = TruncatedFilter(base.labels, base.table, 4, mode="RT-LA", seed=9)
    b = TruncatedFilter(base.labels, base.table, 4, mode="RT", seed=9)
    a.apply(Enter(1))
    b.apply(Enter(1))
    assert a.table == b.table and a.entries == 4


# -- individual updates -------------------------------------------------------------

def test_fov_zero_exit_renormalizes_within_entry():
    f = ExactFilter.from_independent({1: {0: 0.5, 2: 0.5}}, obs_model=symmetric_obs_model(0.9))
    f.apply(FovObservation(1, "exit"))
    assert f.as_dict() == pytest.approx({(1,): 0.45 + 0.5, (3,): 0.05})


def test_fov_no_feasible_branch_drops_entry():
    f = ExactFilter.from_independent({1: {0: 0.5, 2: 0.5}})
    f.apply(FovObservation(1, "exit"))
    assert f.as_dict() == {(1,): 1.0}
    g = ExactFilter.from_counts({1: 0})
    with pytest.raises(InconsistentObservationError):
        g.apply(FovObservation(1, "exit"))


def test_null_observation_with_null_events():
    model = symmetric_obs_model(0.8, null=True)
    f = ExactFilter.from_counts({1: 1}, obs_model=model)
    f.apply(FovObservation(1, "null"))
    assert f.as_dict() == pytest.approx({(0,): 0.1, (1,): 0.8, (2,): 0.1})


def test_noise_free_enter_exit():
    f = ExactFilter.from_independent({1: {0: 0.2, 1: 0.3, 3: 0.5}, 2: {4: 1}})
    f.apply(Exit(1, 1))
    assert f.as_dict() == pytest.approx({(0, 4): 0.375, (2, 4): 0.625})
    f.apply(Enter(2, 2))
    assert f.marginal(2) == {6: 1.0}
    with pytest.raises(InconsistentObservationError):
        f.apply(Exit(1, 3))


def test_appear_disappear():
    f = ExactFilter.from_counts({1: 2})
    f.apply(ProbAppear.of(2, {0: 0.25, 3: 0.75}))
    f.apply(Appear(3, 1, 2))
    assert f.marginal(3) == pytest.approx({1: 0.5, 2: 0.5})
    assert f.probability({2: 3, 3: 1}) == pytest.approx(0.375)
    f.apply(Merge(2, 3, 4))
    f.apply(ProbDisappear.of(4, {1: 0.5, 4: 0.5}))
    assert f.labels == (1,) and f.as_dict() == {(2,): 1.0}
    f.apply(Split(1, 5, 6))
    with pytest.raises(InconsistentObservationError):
        f.apply(Disappear(5, 3, 3))


def test_per_event_split_probability():
    f = ExactFilter.from_counts({1: 3})
    f.apply(ProbSplit(1, 2, 3, 0.25))
    assert f.marginal(2) == pytest.approx({a: math.comb(3, a) * 0.25 ** a * 0.75 ** (3 - a)
                                           for a in range(4)})
    custom = ExactFilter.from_counts({1: 2}, split_rule=lambda n: {(n, 0): 1})
    custom.apply(Split(1, 2, 3))
    assert custom.as_dict() == {(2, 0): 1}


def test_serialization_roundtrip():
    obs = [FovObservation(1, "null"), ProbSplit(1, 2, 3, 0.3), Split(2, 4, 5), Merge(3, 4, 6),
           ProbAppear.of(7, {1: 0.5, 2: 0.5}), Appear(8, 0, 2), Enter(6, 2), Exit(7),
           ProbDisappear.of(5, {0: 1.0}), Disappear(8, 1, 1)]
    seq = ProbSequence((1,), {(4,): 1.0}, obs)
    back = ProbSequence.from_dict(json.loads(json.dumps(seq.to_dict())))
    assert back == seq
    assert observation_from_dict({"type": "fov", "s": 1, "y": "exit", "paper_font": "bold"}) == \
        FovObservation(1, "exit")
    d = observation_to_dict(ProbAppear.of(3, {2: 1.0}))
    assert d == {"type": "appear", "s": 3, "lo": 2, "hi": 2, "dist": {"2": 1.0}}
    ShadowSequence({1: (4, 4)}, [o for o in obs if not isinstance(o, FovObservation)])


def test_serialization_fractions_and_duplicate_keys():
    half = Fraction(1, 2)
    seq = ProbSequence((1,), {(4,): Fraction(1, 3), (5,): Fraction(2, 3)},
                       [ProbSplit(1, 2, 3, half), ProbDisappear.of(2, {1: half, 2: half})])
    d = json.loads(json.dumps(seq.to_dict()))
    assert d["joint"] == [[[4], 1 / 3], [[5], 2 / 3]]
    assert d["events"][0]["p"] == 0.5 and d["events"][1]["dist"] == {"1": 0.5, "2": 0.5}
    dup = ProbSequence.from_dict({"labels": [1], "joint": [[[1], 0.5], [[1], 0.25], [[2], 0.25]],
                                  "events": []})
    assert dup.joint == {(1,): 0.75, (2,): 0.25}


# -- malformed models and labels -------------------------------------------------

def test_fov_impossible_observation_raises():
    model = {"enter": {"enter": 0.9, "exit": 0.1}, "exit": {"enter": 0.1, "exit": 0.9}}
    seq = ProbSequence.from_counts({1: 2, 2: 1}, [FovObservation(1, "null"), Merge(1, 2, 3)])
    with pytest.raises(InconsistentObservationError):
        ExactFilter.from_sequence(seq, obs_model=model).run(seq.observations)
    with pytest.raises(InconsistentObservationError):
        monte_carlo(seq, trials=10, obs_model=model, max_attempts=100)
    for bad in ({"exit": {"Exit": 1.0}}, {"Exit": {"exit": 1.0}}, {"exit": {"exit": -0.5, "enter": 1.5}}):
        with pytest.raises(ValueError):
            ExactFilter.from_counts({1: 2}, obs_model=bad)


def test_unnormalised_obs_row_and_split_rule_are_normalised():
    prior = {1: {0: 0.5, 1: 0.5}}
    raw = ExactFilter.from_independent(prior, obs_model={"exit": {"exit": 1.8, "enter": 0.2}})
    raw.apply(FovObservation(1, "exit"))
    ref = ExactFilter.from_independent(prior, obs_model={"exit": {"exit": 0.9, "enter": 0.1}})
    ref.apply(FovObservation(1, "exit"))
    assert raw.as_dict() == pytest.approx({(0,): 0.45, (1,): 0.5, (2,): 0.05}) == ref.as_dict()
    uniform = ExactFilter.from_independent({1: {0: 0.5, 2: 0.5}},
                                           split_rule=lambda n: {(a, n - a): 1 for a in range(n + 1)})
    uniform.run([Split(1, 2, 3), Merge(2, 3, 4)])
    assert uniform.as_dict() == pytest.approx({(0,): 0.5, (2,): 0.5})
    assert uniform.total() == pytest.approx(1)
    with pytest.raises(ValueError):
        ExactFilter.from_counts({1: 2}, split_rule=lambda n: {(n, 1): 1}).apply(Split(1, 2, 3))


def test_binomial_split_of_many_targets():
    n = 1100
    exact = BinomialSplit(Fraction(3, 10))(n)
    for p in (0.5, 0.3):
        dist = BinomialSplit(p)(n)
        assert math.fsum(dist.values()) == pytest.approx(1, abs=1e-12)
        assert all(v > 0 for v in dist.values())
    dist = BinomialSplit(0.3)(n)
    for k in [(330, 770), (300, 800), (400, 700)]:
        assert dist[k] == pytest.approx(float(exact[k]), rel=1e-9)
    assert BinomialSplit(1.0)(n) == {(n, 0): 1.0} and BinomialSplit(0.0)(n) == {(0, n): 1.0}
    f = ExactFilter.from_counts({1: n})
    f.apply(Split(1, 2, 3))
    assert f.expected_counts() == pytest.approx({2: n / 2, 3: n / 2})


def test_split_merge_appear_into_alive_label_raise():
    for o in (Split(1, 2, 3), Split(1, 4, 4), Merge(1, 3, 2)):
        with pytest.raises(KeyError):
            ExactFilter.from_counts({1: 2, 2: 1, 3: 0}).apply(o)
    f = ExactFilter.from_counts({1: 2, 2: 1})
    f.apply(Merge(1, 2, 1))  # the merged label may reuse a consumed one
    assert f.labels == (1,) and f.as_dict() == {(3,): 1}
    for obs in ([Split(1, 2, 3)], [Merge(1, 3, 2)], [Appear(1, 0, 1)]):
        with pytest.raises(KeyError):
            monte_carlo(ProbSequence.from_counts({1: 2, 2: 1, 3: 0}, obs), trials=5)


def test_underflowed_mass_is_not_stored():
    tiny = 1e-200
    model = {"enter": {"enter": 1 - tiny, "exit": tiny}, "exit": {"exit": 1 - tiny, "enter": tiny},
             "null": {"null": 1}}
    f = ExactFilter.from_counts({1: 3}, obs_model=model)
    for _ in range(3):
        f.apply(FovObservation(1, "enter"))
        assert_normalized(f)
    assert f.as_dict() == pytest.approx({(4,): tiny, (6,): 1.0}, rel=1e-9)
    g = ExactFilter.from_independent({1: {0: tiny, 3: 1.0}, 2: {0: tiny, 1: 1.0}})
    g.apply(Merge(1, 2, 3))
    g.apply(ProbAppear.of(4, {0: tiny, 1: 1.0}))
    assert_normalized(g)


def test_truncation_failure_on_inconsistent_input():
    """TruncationFailure only says the mass ran out after a truncation; the
    exact filter tells whether the input was inconsistent anyway."""
    obs = [Split(1, 3, 4), Disappear(2, 5, 5)]
    with pytest.raises(InconsistentObservationError) as exact:
        ExactFilter((1, 2), {(3, 1): 1.0}).run(obs)
    assert not isinstance(exact.value, TruncationFailure)
    with pytest.raises(TruncationFailure, match=r"no probability mass left after 1 truncation"):
        TruncatedFilter((1, 2), {(3, 1): 1.0}, 2, mode="TR").run(obs)


def test_monte_carlo_labels_reordered_to_exact():
    seq = ProbSequence((2, 1), {(1, 1): 1.0}, [Split(1, 5, 3)])
    exact = ExactFilter.from_sequence(seq).run(seq.observations)
    mc = monte_carlo(seq, trials=2000, seed=3)
    assert exact.labels == (2, 5, 3) and mc.labels == (2, 3, 5)
    got = mc.reordered(exact.labels)
    assert got.labels == exact.labels and set(got.table) == set(exact.table)
    for k, p in exact.table.items():
        assert got.table[k] == pytest.approx(p, abs=0.05)
    with pytest.raises(KeyError):
        mc.reordered((2, 3))


# -- brute force: per-target enumeration --------------------------------------------

def brute_force(labels, joint, observations, split_p, model):
    """Enumerate every history of distinguishable targets and normalise at the end."""
    worlds = []
    for key, p in joint.items():
        targets = tuple(s for s, x in zip(labels, key) for _ in range(x))
        worlds.append((targets, frozenset(labels), p))
    for o in observations:
        out = []
        for targets, alive, w in worlds:
            if isinstance(o, FovObservation):
                inside = [i for i, s in enumerate(targets) if s == o.s]
                weights = {e: model[o.y][e] for e in ("enter", "exit", "null")}
                if not inside:
                    weights["exit"] = 0
                z = sum(weights.values())
                if not z:
                    continue
                if weights["enter"]:
                    out.append((targets + (o.s,), alive, w * weights["enter"] / z))
                if weights["null"]:
                    out.append((targets, alive, w * weights["null"] / z))
                if weights["exit"]:
                    for i in inside:
                        out.append((targets[:i] + targets[i + 1:], alive,
                                    w * weights["exit"] / z / len(inside)))
            elif isinstance(o, Enter):
                out.append((targets + (o.s,) * o.k, alive, w))
            elif isinstance(o, Exit):
                inside = [i for i, s in enumerate(targets) if s == o.s]
                if len(inside) < o.k:
                    continue
                for gone in itertools.combinations(inside, o.k):
                    t = tuple(s for i, s in enumerate(targets) if i not in gone)
                    out.append((t, alive, w / math.comb(len(inside), o.k)))
            elif isinstance(o, Split):
                p = o.p if isinstance(o, ProbSplit) else split_p
                inside = [i for i, s in enumerate(targets) if s == o.s]
                for route in itertools.product((o.a, o.b), repeat=len(inside)):
                    t = list(targets)
                    q = w
                    for i, dest in zip(inside, route):
                        t[i] = dest
                        q *= p if dest == o.a else 1 - p
                    out.append((tuple(t), (alive - {o.s}) | {o.a, o.b}, q))
            elif isinstance(o, Merge):
                t = tuple(o.s if s in (o.a, o.b) else s for s in targets)
                out.append((t, (alive - {o.a, o.b}) | {o.s}, w))
            elif isinstance(o, ProbAppear):
                for n, q in o.dist:
                    out.append((targets + (o.s,) * n, alive | {o.s}, w * q))
            elif isinstance(o, ProbDisappear):
                n = sum(s == o.s for s in targets)
                q = dict(o.dist).get(n, 0)
                if q:
                    out.append((tuple(s for s in targets if s != o.s), alive - {o.s}, w * q))
        worlds = out
    if not worlds:
        return None, None
    final = {}
    for targets, alive, w in worlds:
        final_labels = tuple(sorted(alive))
        key = tuple(sum(s == a for s in targets) for a in final_labels)
        final[key] = final.get(key, 0) + w
    z = sum(final.values())
    return final_labels, {k: v / z for k, v in final.items() if v}


def random_instance(rng: random.Random):
    nshadows = rng.randint(1, 2)
    labels = tuple(range(1, nshadows + 1))
    joint = {}
    for _ in range(rng.randint(1, 3)):
        key = tuple(rng.randint(0, 1) for _ in labels)
        if sum(key) <= 3:
            joint[key] = joint.get(key, 0) + Fraction(rng.randint(1, 4))
    if not joint:
        joint = {tuple(0 for _ in labels): Fraction(1)}
    alive = list(labels)
    nxt = nshadows + 1
    obs = []
    for _ in range(rng.randint(1, 6)):
        kind = rng.choice(["fov", "fov", "split", "merge", "appear", "disappear", "enter", "exit"])
        s = rng.choice(alive) if alive else None
        if kind == "appear" or not alive:
            dist = {n: Fraction(rng.randint(1, 3), 6) for n in rng.sample(range(3), 2)}
            obs.append(ProbAppear.of(nxt, dist))
            alive.append(nxt)
            nxt += 1
        elif kind == "fov":
            obs.append(FovObservation(s, rng.choice(["enter", "exit", "null"])))
        elif kind == "split":
            p = Fraction(rng.randint(1, 4), 5)
            obs.append(ProbSplit(s, nxt, nxt + 1, p) if rng.random() < 0.5 else Split(s, nxt, nxt + 1))
            alive.remove(s)
            alive += [nxt, nxt + 1]
            nxt += 2
        elif kind == "merge" and len(alive) >= 2:
            a, b = rng.sample(alive, 2)
            obs.append(Merge(a, b, nxt))
            alive = [x for x in alive if x not in (a, b)] + [nxt]
            nxt += 1
        elif kind == "disappear":
            obs.append(ProbDisappear.of(s, {n: Fraction(1, 2) for n in rng.sample(range(4), 2)}))
            alive.remove(s)
        elif kind == "enter":
            obs.append(Enter(s))
        elif kind == "exit":
            obs.append(Exit(s))
    return labels, joint, obs


def random_model(rng: random.Random):
    model = {}
    for y in ("enter", "exit", "null"):
        w = [Fraction(rng.randint(0, 3)) for _ in range(3)]
        if not any(w):
            w[rng.randrange(3)] = Fraction(1)
        model[y] = {e: x / sum(w) for e, x in zip(("enter", "exit", "null"), w)}
    return model


@pytest.mark.parametrize("seed", range(150))
def test_exact_equals_brute_force(seed):
    rng = random.Random(seed)
    labels, joint, obs = random_instance(rng)
    model = random_model(rng)
    split_p = Fraction(rng.randint(1, 4), 5)
    final_labels, expected = brute_force(labels, joint, obs, split_p, model)
    f = ExactFilter(labels, joint, split_rule=BinomialSplit(split_p), obs_model=model)
    if expected is None:
        with pytest.raises(InconsistentObservationError):
            f.run(obs)
        return
    f.run(obs)
    order = [f.labels.index(s) for s in final_labels]
    got = {tuple(k[i] for i in order): p for k, p in f.table.items()}
    assert got == expected
    assert f.total() == 1
    t = TruncatedFilter(labels, joint, max(f.peak_entries, 1), mode="RT-LA", seed=seed,
                        split_rule=BinomialSplit(split_p), obs_model=model).run(obs)
    assert t.table == f.table and t.truncations == 0


@pytest.mark.parametrize("seed", [3, 8, 13])
def test_monte_carlo_matches_exact_on_random_instances(seed):
    rng = random.Random(seed)
    while True:
        labels, joint, obs = random_instance(rng)
        model = random_model(rng)
        try:
            exact = ExactFilter(labels, joint, obs_model=model).run(obs)
        except InconsistentObservationError:
            continue
        if exact.labels and len(obs) >= 3:
            break
    n = 4000
    floats = {y: {e: float(p) for e, p in row.items()} for y, row in model.items()}
    mc = monte_carlo(ProbSequence(labels, joint, obs), trials=n, seed=seed, obs_model=floats)
    order = [mc.labels.index(s) for s in exact.labels]
    got = {tuple(k[i] for i in order): p for k, p in mc.table.items()}
    for k in set(got) | set(exact.table):
        p = float(exact.table.get(k, 0))
        assert abs(got.get(k, 0) - p) < 4.5 * math.sqrt(p * (1 - p) / n) + 2e-3, k


@pytest.mark.parametrize("mode", ["TR", "RT", "RT-LA"])
def test_mass_sums_to_one_under_truncation(mode):
    rng = random.Random(99)
    checked = 0
    for _ in range(60):
        labels, joint, obs = random_instance(rng)
        model = {y: {e: float(p) for e, p in r.items()} for y, r in random_model(rng).items()}
        f = TruncatedFilter(labels, {k: float(p) for k, p in joint.items()}, 2, mode=mode, seed=4,
                            obs_model=model)
        try:
            for t, o in enumerate(obs):
                f.apply(o, obs[t + 1:])
                assert_normalized(f)
                checked += 1
        except InconsistentObservationError:
            pass
    assert checked > 50
