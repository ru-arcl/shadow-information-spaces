"""Oracle, equations, bipartite I-state and bounds (shadowinfo.polygon.oracle / stagent / __main__) vs. the Java.

Golden data: tests/fixtures/java (tools/java_reference; docs/notes/original_java.md §5.2 items 7, 8):

7. ``oracle_runs`` -- the original ``SingleTypeAgentOracle`` (with its merge bug B7) and the harness's
   ``MergeFixedOracle``, seeded through ``Math.random()``: the event list, ``toString()``, the final hidden
   counts, the equations and the bipartite graph of ``deriveShadowInfoState``;
8. ``java_bounds_outcomes`` -- every distinct ``deriveShadowBounds`` output over 24 identity-hash orders.

``JavaRandom`` + ``SingleTypeAgentOracle(compat="java")`` must reproduce item 7 exactly, and so must
``compat="fixed"`` on the same (Java-shaped) histories.  ``java_bounds(order="sorted")`` must equal the Java
output where the Java output is order-independent; the exact bounds must contain the true counts of every
consistent oracle and detect the infeasible observations of the merge bug (all recorded T-RO Fig. 15(b) runs).
The command line must print what the ``ProjectPanel*`` constructors print.
"""

from __future__ import annotations

import io
import json
import math
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path as FsPath

import pytest

from shadowinfo.events import Appear, Disappear, Merge, ShadowSequence, Split
from shadowinfo.maxflow import InfeasibleError
from shadowinfo.nondeterministic import CombinatorialFilter
from shadowinfo.polygon import (Equation, EventType, GapTrackingError, JavaRandom, OracleError, SingleTypeAgentEvent,
                                SingleTypeAgentOracle, bounds_lines, derive_gap_evolving_equations,
                                derive_shadow_info_state, exact_bounds, gap_history_to_events, get_gaps,
                                java_bounds, load_path, path_labels, path_record, project_panel_lines, sample_path,
                                transition_events)
from shadowinfo.polygon import __main__ as cli
from shadowinfo.polygon.gaps import Gap, GapHistory
from shadowinfo.polygon.oracle import strict_log

ROOT = FsPath(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "java"
sys.path.insert(0, str(ROOT / "tools" / "java_reference"))

import convert  # noqa: E402

BIG = {5, 13, 14}


def _load(name):
    return json.loads((FIX / f"{name}.json").read_text())


RUNS = sorted(p.stem for p in FIX.glob("*_*.json"))
FIXTURES = {name: _load(name) for name in RUNS}
OK_RUNS = [r for r in RUNS if FIXTURES[r]["get_gaps"]["ok"]]
ORACLE_CASES = [(r, k) for r in OK_RUNS for k in range(len(FIXTURES[r]["oracle_runs"]))]
# all 5 seeds of the T-RO Fig. 15(b) panel run and the digitised Fig. 15(b) path (README of the harness)
INFEASIBLE = {("14_PP5", "original", s) for s in range(1, 6)} | {("14_fig_TRO-Fig15b", "original", 1)}


def _slow(name):
    return int(name.split("_")[0]) in BIG


def _param(names):
    return [pytest.param(n, marks=pytest.mark.slow) if _slow(n) else n for n in names]


def _param_cases(cases):
    return [pytest.param(*c, marks=pytest.mark.slow, id=f"{c[0]}-{c[1]}") if _slow(c[0])
            else pytest.param(*c, id=f"{c[0]}-{c[1]}") for c in cases]


_SAMPLES: dict = {}
_HIST: dict = {}


def _samples(label):
    """Sampled critical points of a stored path (mode independent: no exact 45-degree segment)."""
    label = cli.path_label(label)  # one cache entry per path, whatever name is used
    if label not in _SAMPLES:
        poly, path = load_path(label)
        s = sample_path(poly, path, "java")
        assert s.mode_independent
        _SAMPLES[label] = s
    return _SAMPLES[label]


def history(label, compat="java"):
    label = cli.path_label(label)
    if (label, compat) not in _HIST:
        s = _samples(label)
        _HIST[label, compat] = get_gaps(s.polygon, s.path, compat, samples=s)
    return _HIST[label, compat]


def fixture_history(name):
    poly_n, label = name.split("_", 1)
    return history(f"{label}@{poly_n}")


def java_oracle(name, k, compat="java"):
    orun = FIXTURES[name]["oracle_runs"][k]
    o = SingleTypeAgentOracle(orun["n_targets"], orun["seed"], merge_bug=orun["variant"] == "original",
                              compat=compat)
    return o.initialize(fixture_history(name)), orun


def outcomes(orun):
    return [{int(s): tuple(b) for s, b in oc["bounds"].items()} for oc in orun["java_bounds_outcomes"]]


# --------------------------------------------------------------------------- java.util.Random

# OpenJDK 17.0.20.1: new Random(seed); 5 x nextDouble, 5 x nextInt(), nextInt(b) for the bounds below,
# 3 x nextLong, 3 x nextBoolean, 3 x nextFloat, 6 x nextGaussian (generated with the JDK of the harness).
BOUNDS = [1, 2, 3, 7, 10, 16, 100, 1000000, 1 << 30, 2147483647, 1073741825]
JDK17 = {
    42: "0.7275636800328681 0.6832234717598454 0.30871945533265976 0.27707849007413665 0.6655489517945736 "
        "-415012931 1938135004 1583910553 1639144584 1184328952 0 0 0 0 9 14 63 357226 472231456 1610411243 "
        "830146030 -3720955133854962410 2777110785106071796 6240679879620383709 false true false 0.1722179 "
        "0.15865451 0.5874274 -1.3341828825971889 0.8485163405804642 1.0157756767477628 -0.6240830275297751 "
        "0.05931978782536116 -0.09229029794439661",
    -7: "0.2691218093260761 0.6731811902798038 0.5214896372918163 0.8518782305056931 0.08574465212415472 "
        "1976714354 -2105189454 -430535443 -100754905 740122818 0 1 0 5 7 5 17 717690 935790419 701801676 "
        "16831854 4909513586090054363 693266197946214462 3675442807801244372 false true false 0.74494654 "
        "0.66148126 0.1327824 -0.5637057298271061 0.6205894041770273 -0.3156321051884695 -1.5102040510739227 "
        "0.7213789004340653 -0.46379413490177135",
    0x5DEECE66D: "1.4683476656784933E-11 0.0416310008860028 0.36460223930076674 0.09229765068318274 "
                 "0.5267502733235571 1001493746 -724596254 -293211053 -1855169837 -1906560325 0 1 2 2 3 14 74 "
                 "245324 553439984 1740383999 404629406 -2097128362344201561 1416240747525818885 "
                 "-278711534114550735 true true false 0.25331062 0.019841611 0.37837744 0.7999714132183818 "
                 "1.1301488277057372 1.5051474995431822 1.4221884786544783 -1.3653364099486198 "
                 "0.041380253987791736",
}


@pytest.mark.parametrize("seed", sorted(JDK17))
def test_java_random_equals_jdk17(seed):
    import numpy as np

    want = JDK17[seed].split()
    r = JavaRandom(seed)
    got = [r.next_double() for _ in range(5)] + [r.next_int() for _ in range(5)] + [r.next_int(b) for b in BOUNDS]
    got += [r.next_long() for _ in range(3)] + [r.next_boolean() for _ in range(3)]
    floats = [r.next_float() for _ in range(3)]
    gauss = [r.next_gaussian() for _ in range(6)]
    assert [float(x) for x in want[:5]] == got[:5]
    assert [int(x) for x in want[5:24]] == got[5:24]
    assert [x == "true" for x in want[24:27]] == got[24:27]
    assert [float(np.float32(x)) for x in want[27:30]] == floats  # Float.toString is shortest-float32
    assert [float(x) for x in want[30:]] == gauss  # bit for bit thanks to strict_log (fdlibm)


def test_java_random_index_check_and_reference_clone():
    chk = json.loads((FIX / "index.json").read_text())["java_random_check"]
    r = JavaRandom(chk["seed"])
    assert [r.next_double() for _ in chk["math_random_after_seeding"]] == chk["math_random_after_seeding"]
    r = JavaRandom(42)
    assert [int(r.next_double() * 100) for _ in range(6)] == chk["new_Random_42_int_nextDouble_x100"]
    r = JavaRandom(12345)
    assert [r.random() for _ in range(50)] == convert.java_random_doubles(12345, 50)
    r.set_seed(12345)
    assert r.next_double() == convert.java_random_doubles(12345, 1)[0]
    with pytest.raises(ValueError):
        r.next_int(0)
    assert JavaRandom().next_double() < 1.0  # unseeded, like Math.random()


def test_strict_log():
    xs = [1.0, 2.0, 0.5, math.e, 1e-300, 5e-324, 1.0 + 2 ** -30, 1.0 - 2 ** -40, 0.999, 123456.789, 1e300]
    for x in xs:
        assert abs(strict_log(x) - math.log(x)) <= 2 * math.ulp(math.log(x)) + 5e-324, x
    assert strict_log(1.0) == 0.0 and strict_log(0.0) == -math.inf and math.isnan(strict_log(-1.0))
    assert strict_log(math.inf) == math.inf


# --------------------------------------------------------------------------- events and equations


def test_event_and_equation_strings():
    ev = [SingleTypeAgentEvent(EventType.INIT, 0.0, None, [1, 2, 3, 4], 36, 0),
          SingleTypeAgentEvent(EventType.APPEAR, 0.1, None, [5], 2, 34),
          SingleTypeAgentEvent(EventType.DISAPPEAR, 0.2, [14], None, 37, 36),
          SingleTypeAgentEvent(EventType.SPLIT, 0.3, [1], [6, 7]),
          SingleTypeAgentEvent(EventType.MERGE, 0.4, [7, 2], [15]),
          SingleTypeAgentEvent(EventType.AGENT_APPEAR, 0.5, [3], None, 0, 2),
          SingleTypeAgentEvent(EventType.AGENT_DISAPPEAR, 0.6, None, [4], 0, 1)]
    assert [str(e) for e in ev] == ["INIT EVENT 1 2 3 4", "APPEAR EVENT 5 <- 34", "DISAPPEAR EVENT 14 -> 36",
                                    "SPLIT EVENT 1 -> 6 7", "MERGE EVENT 7 2 -> 15",
                                    "AGENT_APPEAR: 2 appeared from gap 3", "AGENT_DISAPPEAR: 1 disappeared from gap 4"]
    assert ev[1].shadow_event() == Appear(5, 34, 34) and ev[2].shadow_event() == Disappear(14, 36, 36)
    assert ev[3].shadow_event() == Split(1, 6, 7) and ev[4].shadow_event() == Merge(7, 2, 15)
    assert ev[0].shadow_event() is None
    with pytest.raises(OracleError) as ei:
        str(SingleTypeAgentEvent())  # null type: NPE in the switch
    assert ei.value.java_exception == "java.lang.NullPointerException"
    # Equation.main prints "5 x3 + 8 x4 = 120"; negative coefficients are joined by a space only
    e = Equation()
    e.add_term(3, 5)
    e.add_term(4, 8)
    e.constant = 120
    assert str(e) == "5 x3 + 8 x4 = 120"
    e = Equation({11: 1, 10: -1, 4: -1})
    assert str(e) == "-1 x4 -1 x10 + 1 x11 = 0" and e.evaluate({4: 1, 10: 2, 11: 3}) == 0
    e.add_term(11, 7)  # overwrites
    assert str(e) == "-1 x4 -1 x10 + 7 x11 = 0" and str(Equation()) == " = 0"


def test_distribute_one_batch():
    o = SingleTypeAgentOracle(10, 3, compat="java")
    for t, g in [(0, 3), (5, 1), (10, 2), (7, 7), (1000000, 5)]:
        b = o.distribute_one_batch(t, g)
        assert len(b) == g and sum(b) == t and min(b) >= 0
        if t and g > 1:
            assert b[-1] >= 1  # B8: the last part is never empty
    assert SingleTypeAgentOracle(10, 1, compat="java").distribute_one_batch(1, 2) == [0, 1]  # B8
    with pytest.raises(OracleError, match="never terminates"):
        o.distribute_one_batch(2, 5)  # Java loops forever
    f = SingleTypeAgentOracle(10, 3)
    b = f.distribute_one_batch(2, 5)
    assert len(b) == 5 and sum(b) == 2


# --------------------------------------------------------------------------- item 7: the oracle runs


@pytest.mark.parametrize("name,k", _param_cases(ORACLE_CASES))
def test_oracle_runs_equal_java(name, k):
    h = fixture_history(name)
    for compat in ("java", "fixed"):
        o, orun = java_oracle(name, k, compat)
        assert o.merge_bug == (orun["variant"] == "original")
        assert [e.to_row() for e in o.events] == [r[:6] for r in orun["events"]], compat
        assert o.number_of_events == orun["n_events"]
        if len(orun["events"][0]) > 6:
            assert o.event_strings() == [r[6] for r in orun["events"]]
        assert o.truth() == orun["truth_final"] == o.truth(h.final_ids)
        eqs = derive_gap_evolving_equations(h, o)
        assert [str(e) for e in eqs] == orun["equations"], compat
        bg = derive_shadow_info_state(h, o)
        d = bg.to_dict()
        assert d == {key: orun["bg"][key] for key in d}, compat
        if "bg_text_example" in orun:  # Java prints the adjacency in identity-hash order
            def norm(text):
                return [re.sub(r"goes to: (.*)", lambda m: "goes to: " + " ".join(sorted(m.group(1).split(), key=int)),
                               ln) for ln in text.split("\n")]

            assert norm(str(bg) + "\n") == norm(orun["bg_text_example"])
    # the reference replica of the harness agrees as well
    assert [r[:6] for r in orun["events"]] == convert.replay_oracle(
        {"get_gaps": FIXTURES[name]["get_gaps"]}, orun["seed"], orun["n_targets"], orun["variant"] == "original")


def test_oracle_accessors():
    """``Oracle.getTotalAgentNumber``, ``getNumberOfEvents``, ``getEvent(i)`` and ``getEventByTime(t)`` (the lookup
    ``deriveGapEvolvingEquations`` uses: ``gapss[k][0].relativeTime``) on the ProjectPanel run, seed 1."""
    h = fixture_history("12_PP")
    o, orun = java_oracle("12_PP", 0)
    assert o.get_total_agent_number() == orun["n_targets"] == 100
    assert o.number_of_events == len(h) == 15
    assert [o.get_event(i).to_row() for i in range(o.number_of_events)] == [r[:6] for r in orun["events"]]
    assert [o.event_by_time(gs[0].relative_time) for gs in h] == o.events
    assert o.event_by_time(0.5) is None


@pytest.mark.parametrize("name,k", _param_cases(ORACLE_CASES))
def test_ground_truth_and_equations_hold(name, k):
    """Every equation holds for the oracle's own counts while targets are conserved (merge-fixed oracle);
    the merge bug breaks exactly the merge equations."""
    o, orun = java_oracle(name, k)
    h = fixture_history(name)
    ids = {g.id for gs in h for g in gs}
    counts = o.label_counts
    assert set(counts) == ids and all(v is not None and v >= 0 for v in counts.values())
    visible = [e.visible for e in o.events]
    for kk, c in enumerate(o.counts):
        total = sum(c.values()) + visible[kk]
        if orun["variant"] == "merge_fixed":
            assert total == orun["n_targets"]
    eqs = derive_gap_evolving_equations(h, o)
    assert len(eqs) == len(o.events)  # Java-shaped: one event per set, distinct times
    broken = [e.type for e, eq in zip(o.events, eqs) if eq.evaluate(counts) != 0]
    if orun["variant"] == "merge_fixed":
        assert not broken
    else:
        assert set(broken) <= {EventType.MERGE}  # B7 breaks exactly the merge equations


# --------------------------------------------------------------------------- item 8: bounds


@pytest.mark.parametrize("name,k", _param_cases(ORACLE_CASES))
def test_java_bounds_and_exact_bounds(name, k):
    o, orun = java_oracle(name, k)
    bg = derive_shadow_info_state(fixture_history(name), o)
    outs = outcomes(orun)
    jb = java_bounds(bg)
    if not orun["java_bounds_order_dependent"]:
        assert dict(jb) == outs[0]  # Java's answer does not depend on the hash order: we must print it
    assert dict(jb) in outs  # observed on every golden run, though Java does not guarantee it
    assert jb.lines == [f"{w} flow in shadow {s} is {jb[s][i]}" for s in jb for w, i in (("Max", 1), ("Min", 0))]
    assert bg.to_dict() == derive_shadow_info_state(fixture_history(name), o).to_dict()  # bg not mutated
    seq = o.sequence()
    key = (name, orun["variant"], orun["seed"])
    if key in INFEASIBLE:
        with pytest.raises(InfeasibleError):
            exact_bounds(o)
        assert not CombinatorialFilter.from_sequence(seq, total=(o.hidden_total,) * 2).feasible()
        return
    ex = exact_bounds(o)
    assert ex == convert.exact_bounds(seq, o.hidden_total)
    assert set(ex) == set(fixture_history(name).final_ids)
    assert set(jb) == {s for s in ex if s not in {e.s for e in seq.events if isinstance(e, Appear)}}  # B11
    truth = o.final_counts
    if orun["variant"] == "merge_fixed":
        assert all(ex[s][0] <= truth[s] <= ex[s][1] for s in ex)
    if any(dict(jb) == oc for oc in outs) and all(oc == ex for oc in outs):
        assert dict(jb) == {s: b for s, b in ex.items() if s in jb}


def test_infeasible_runs_are_exactly_the_merge_bug_fig15b_runs():
    seen = set()
    for name in OK_RUNS:
        for orun in FIXTURES[name]["oracle_runs"]:
            seq, hidden, _t, _p = convert.oracle_sequence(orun)
            if convert.exact_bounds(seq, hidden) is None:
                seen.add((name, orun["variant"], orun["seed"]))
    assert seen == INFEASIBLE


def test_random_hash_orders_reproduce_every_java_outcome_on_p12():
    for k, orun in enumerate(FIXTURES["12_PP"]["oracle_runs"]):
        o, _ = java_oracle("12_PP", k)
        bg = derive_shadow_info_state(None, o)
        ours = {json.dumps(sorted(java_bounds(bg, seed).items())) for seed in range(200)}
        theirs = {json.dumps(sorted(oc.items())) for oc in outcomes(orun)}
        assert theirs <= ours, orun["seed"]


def test_b10_counterexample():
    """docs/notes/original_java.md B10: over 40 identity-hash orders Java printed 5:[0,0] 8:[1,1] (31x),
    5:[0,1] 8:[0,1] (7x) and the wrong 5:[1,1] 8:[0,0] (2x); the truth is s5 = 0, s8 = 1."""
    inf = math.inf
    seq = ShadowSequence({1: (0, inf), 2: (0, inf)}, [
        Appear(3, 1, 1), Split(3, 4, 5), Merge(1, 4, 6), Appear(7, 1, 1), Split(7, 8, 9), Disappear(6, 1, 1),
        Disappear(2, 1, 1), Disappear(9, 0, 0)])
    o = SingleTypeAgentOracle.from_sequence(seq, n_targets=3, visible0=2)
    assert o.hidden_total == 1 and o.sequence() == seq
    assert [str(e) for e in derive_gap_evolving_equations(None, o)][:2] == ["1 x1 + 1 x2 = 1", "1 x3 = 1"]
    bg = derive_shadow_info_state(None, o)
    assert bg.edges() == [(0, 2), (0, 6), (3, 5), (3, 6), (7, 8), (7, 9)]
    assert exact_bounds(o) == {5: (0, 0), 8: (1, 1)}
    seen = Counter(tuple(sorted(java_bounds(bg, seed).items())) for seed in range(200))
    assert set(seen) == {((5, (0, 0)), (8, (1, 1))), ((5, (0, 1)), (8, (0, 1))), ((5, (1, 1)), (8, (0, 0)))}
    assert seen.most_common(1)[0][0] == ((5, (0, 0)), (8, (1, 1)))
    assert dict(java_bounds(bg)) == {5: (0, 0), 8: (1, 1)}
    import random

    assert dict(java_bounds(bg, random.Random(3))) == dict(java_bounds(bg, 3))
    with pytest.raises(ValueError):
        java_bounds(bg, "unsorted")


# --------------------------------------------------------------------------- compat="fixed" end to end


FIXED_PATHS = [pytest.param(lab, marks=pytest.mark.slow) if path_record(lab)["polygon"] in BIG else lab
               for lab in path_labels()]


@pytest.mark.parametrize("label", FIXED_PATHS)
def test_fixed_pipeline_on_every_path(label):
    """Default mode on every stored path, including the four on which Java's getGaps throws: the oracle
    reports exactly the component events of the history, conserves targets, all its equations hold and
    the exact bounds contain the true counts."""
    h = history(label, "fixed")
    n = path_record(label)["harness"]["n_targets"]
    ge = gap_history_to_events(h)
    for seed in (1, 2):
        o = SingleTypeAgentOracle(n, seed).initialize(h)
        assert len(o.events) == 1 + len(ge.events)
        for oe, e in zip(o.events[1:], ge.events):
            se = oe.shadow_event()
            if isinstance(e, (Appear, Disappear)):
                assert type(se) is type(e) and se.s == e.s
            else:
                assert se == e
        assert [e.time for e in o.events[1:]] == ge.times
        assert len(o.counts) == len(h) and all(set(c) == {g.id for g in gs} for c, gs in zip(o.counts, h))
        assert sum(o.final_counts.values()) + o.events[-1].visible == n
        eqs = derive_gap_evolving_equations(h, o)
        assert len(eqs) == len(o.events) and all(eq.evaluate(o.label_counts) == 0 for eq in eqs)
        assert set(o.label_counts) >= {g.id for gs in h for g in gs}
        ex = exact_bounds(o)
        truth = o.final_counts
        assert set(ex) == set(truth) and all(ex[s][0] <= truth[s] <= ex[s][1] for s in ex)
    # replaying the original merge bug on this history makes the T-RO Fig. 15(b) observations infeasible
    if label in ("P14b", "fig_TRO-Fig15b"):
        with pytest.raises(InfeasibleError):
            exact_bounds(SingleTypeAgentOracle(n, 1, merge_bug=True).initialize(h))


def test_fixed_oracle_handles_multi_event_transitions():
    """P13+ leaves polygon 13: its fixed history has transitions with several events and >2-way components,
    which the Java oracle cannot represent (B9)."""
    h = history("P13+", "fixed")
    assert max(len(x) for x in transition_events(h)) > 1
    o = SingleTypeAgentOracle(100, 7).initialize(h)
    assert o.sequence().events == [e.shadow_event() for e in o.events[1:]]
    assert sum(o.final_counts.values()) + o.events[-1].visible == 100
    with pytest.raises(GapTrackingError):  # Java's getGaps already throws on this path
        history("P13+", "java")


def _java_like_history(sets, times):
    """A hand-made gap history: ``sets`` = [[(id, [links])...]...]."""
    h = GapHistory()
    for gs, t in zip(sets, times):
        row = []
        for gid, links in gs:
            g = Gap(None, gid, 0)
            for x in links:
                g.add_gap(x)
            row.append(g)
        row[0].relative_time = t
        h.append(row)
    h.times = list(times)
    return h


def test_java_quirks_untyped_events_and_equal_times():
    # B9: a transition without an event leaves the type null; equations and graph throw NPE in Java
    h = _java_like_history([[(1, []), (2, [])], [(1, []), (2, [])]], [0.0, 0.5])
    o = SingleTypeAgentOracle(10, 1, compat="java").initialize(h)
    assert o.events[1].type is None
    with pytest.raises(OracleError) as ei:
        derive_gap_evolving_equations(h, o)
    assert ei.value.java_exception == "java.lang.NullPointerException"
    with pytest.raises(OracleError):
        derive_shadow_info_state(h, o)
    assert len(SingleTypeAgentOracle(10, 1).initialize(h).events) == 1  # fixed: no event, nothing invented
    # TreeMap<Double, Event>: two sets with the same time keep only the later event
    h = _java_like_history([[(1, [2, 3])], [(2, []), (3, [])], [(2, []), (3, []), (4, [])]], [0.0, 0.5, 0.5])
    o = SingleTypeAgentOracle(10, 1, compat="java").initialize(h)
    assert [e.type for e in o.events] == [EventType.INIT, EventType.APPEAR]
    assert [str(e) for e in derive_gap_evolving_equations(h, o)][1:] == ["1 x4 = %d" % o.events[1].moved] * 2
    f = SingleTypeAgentOracle(10, 1).initialize(h)
    assert [e.type for e in f.events] == [EventType.INIT, EventType.SPLIT, EventType.APPEAR]
    # a merge parent without a state in Java (gap missing from the next set): NullPointerException
    h = _java_like_history([[(1, [3]), (2, [9])], [(3, [])]], [0.0, 0.5])
    with pytest.raises(OracleError):
        SingleTypeAgentOracle(10, 1, compat="java").initialize(h)
    # merge and disappear in one transition: the disappearance between the two merging parents replaces
    # fromGap with a 1-element array, and the second parent's e.getFromGap()[1] = id throws (live Java: AIOOBE)
    h = _java_like_history([[(1, [4]), (2, []), (3, [4])], [(4, [])]], [0.0, 0.5])
    with pytest.raises(OracleError) as ei:
        SingleTypeAgentOracle(100, 1, compat="java").initialize(h)
    assert ei.value.java_exception == "java.lang.ArrayIndexOutOfBoundsException"
    assert len(SingleTypeAgentOracle(100, 1).initialize(h).counts) == 2  # fixed: no crash


def test_merge_bug_in_java_matches_the_original_line():
    """B7 on a minimal history: the second parent's count is doubled and the first one's is lost."""
    h = _java_like_history([[(1, [3]), (2, [3])], [(3, [])]], [0.0, 0.5])
    for bug in (True, False):
        o = SingleTypeAgentOracle(100, 5, merge_bug=bug, compat="java").initialize(h)
        c1, c2 = o.counts[0][1], o.counts[0][2]
        assert o.counts[1][3] == (2 * c2 if bug else c1 + c2)
        f = SingleTypeAgentOracle(100, 5, merge_bug=bug).initialize(h)
        assert [e.to_row() for e in f.events] == [e.to_row() for e in o.events] and f.counts == o.counts


# --------------------------------------------------------------------------- the command line


def _cli(*argv):
    buf = io.StringIO()
    rc = cli.main(list(argv), out=buf)
    return rc, buf.getvalue()


def _panel_expected(name, k, bounds):
    run = FIXTURES[name]
    orun = run["oracle_runs"][k]
    text = "".join(s["line"] + "\n" for s in run["get_gaps"]["sets"]) + "\n\n"
    text += orun["bg_text_example"]
    text += "".join(f"Max flow in shadow {s} is {bounds[s][1]}\nMin flow in shadow {s} is {bounds[s][0]}\n"
                    for s in sorted(bounds))
    text += "\n\nComputation time: <ms>\n" + "".join(e + "\n" for e in orun["equations"]) + "\n"
    text += "".join(r[6] + "\n" for r in orun["events"])
    return text


def _normalize(text):
    text = re.sub(r"Computation time: \d+", "Computation time: <ms>", text)
    return re.sub(r"goes to: ([^\n]*)", lambda m: "goes to: " + " ".join(sorted(m.group(1).split(), key=int)) + " ",
                  text)


def _printed_bounds(text):
    b: dict = {}
    for w, s, v in re.findall(r"(Max|Min) flow in shadow (-?\d+) is (-?\d+)", text):
        b.setdefault(int(s), [None, None])[0 if w == "Min" else 1] = int(v)
    return {s: tuple(v) for s, v in b.items()}


def test_cli_project_panel_java_equals_golden_output():
    """``run ProjectPanel --compat java --seed 1`` prints the ProjectPanel stdout of the seeded Java run."""
    rc, out = _cli("run", "ProjectPanel", "--compat", "java", "--seed", "1")
    assert rc == 0
    orun = FIXTURES["12_PP"]["oracle_runs"][0]
    assert orun["seed"] == 1 and orun["variant"] == "original"
    printed = _printed_bounds(out)
    assert printed in outcomes(orun)  # order-dependent in Java: one of the outcomes it can print
    assert _normalize(out) == _normalize(_panel_expected("12_PP", 0, printed))
    assert out.startswith("0.0 [1, 0] [2, 11] [3, 16] [4, 23] \n0.056234042553191504 [1, 0 | 6 7 ] ")


@pytest.mark.slow
def test_cli_project_panel5_java_equals_golden_output():
    """T-RO Fig. 15(b) panel, seed 1: Java's bounds do not depend on the hash order there, so every line
    except the computation time equals the Java output (adjacency lists compared as sets)."""
    rc, out = _cli("run", "ProjectPanel5", "--compat", "java", "--seed", "1")
    assert rc == 0
    orun = FIXTURES["14_PP5"]["oracle_runs"][0]
    assert not orun["java_bounds_order_dependent"]
    assert _normalize(out) == _normalize(_panel_expected("14_PP5", 0, outcomes(orun)[0]))


def test_cli_nse_panel_and_drawing_panel():
    rc, out = _cli("run", "ProjectPanel2", "--compat", "java")
    assert rc == 0
    run = FIXTURES["1_PP2"]
    want = []
    for s, st in zip(run["get_gaps"]["sets"], run["nse"]["states"]):
        want.append("".join(f"[{gid}, {ie}| {c}" + (" | " + "".join(f"{x} " for x in to) if to else "") + "] "
                            for (gid, ie, to), c in zip(s["gaps"], st)))
    assert out == "".join(w + "\n" for w in want) + "\n\n"
    assert _cli("run", "ProjectPanel2")[1] == out  # fixed mode: same on a Java-shaped history
    assert _cli("run", "ProjectPanel4", "--compat", "java") == (0, "")  # drawing only


def test_cli_fixed_mode_and_errors(capsys):
    rc, out = _cli("run", "ProjectPanel", "--seed", "1")
    assert rc == 0
    h = history("P12", "fixed")
    o = SingleTypeAgentOracle(100, 1).initialize(h)
    assert _printed_bounds(out) == exact_bounds(o) and 18 in _printed_bounds(out)  # B11 fixed
    assert "\n".join(bounds_lines(exact_bounds(o))) in out
    lines = list(project_panel_lines(h, SingleTypeAgentOracle(100, 1)))
    assert _normalize("".join(x + "\n" for x in lines)) == _normalize(out)
    # same by map number and by path label
    assert _normalize(_cli("run", "12", "--seed", "1")[1]) == _normalize(out)
    assert _normalize(_cli("run", "P12", "--seed", "1", "--targets", "100")[1]) == _normalize(out)
    # Java's getGaps throws on P13+: nothing printed, the exception on stderr, exit status 1
    rc, out = _cli("run", "P13+", "--compat", "java", "--seed", "1")
    assert rc == 1 and out == ""
    assert "ArrayIndexOutOfBoundsException" in capsys.readouterr().err
    rc, out = _cli("run", "13", "--path", "P13+", "--seed", "1")
    assert rc == 0 and "Computation time" in out
    assert "warning: path P13+:" in capsys.readouterr().err  # fixed mode tracks it but says it leaves the map
    _cli("run", "P13", "--seed", "1")
    assert capsys.readouterr().err == ""
    with pytest.raises(SystemExit):
        _cli("run", "12", "--path", "P13")
    with pytest.raises(SystemExit):
        _cli("run", "nonsense")


def test_cli_cuts_vis_list():
    poly = FIXTURES["12_PP"]["polygon"]
    counts = json.loads((FIX / f"poly{poly}.json").read_text())["counts"]
    rc, out = _cli("cuts", "12", "--summary")
    assert rc == 0
    assert (f"single tangents {counts['single_tangent']}, non-general inflections {counts['inflection_nongeneral']}, "
            f"general inflections {counts['general_inflection_cut']}, bitangent rays {counts['bitangent_cut']}") in out
    rc, out = _cli("cuts", "12")
    total = sum(counts[k] for k in ("single_tangent", "inflection_nongeneral", "general_inflection_cut",
                                    "bitangent_cut"))
    body = out.splitlines()[2:]
    assert len(body) == total and body[0].startswith("0 SINGLETANGENT ")
    assert any("fromLine=" in x for x in body) and any("this=" in x and "opp=" in x for x in body)
    vis = json.loads((FIX / "poly12.json").read_text())["visibility"][0]
    rc, out = _cli("vis", "12", str(vis["q"][0]), str(vis["q"][1]), "--compat", "java")
    assert rc == 0 and "inside" in out
    sp = vis["physical_gaps"][0]
    assert f"[{sp[0]}, {sp[1]}, ({sp[2][0]}, {sp[2][1]}), null]" in out
    assert out.endswith("".join(f"  {x} {y}\n" for x, y in vis["visibility_polygon"]))
    rc, out = _cli("list")
    assert rc == 0 and all(lab in out for lab in path_labels())
    assert _cli("cuts", "15") == (2, "")


def test_cli_entry_point():
    r = subprocess.run([sys.executable, "-m", "shadowinfo.polygon", "run", "ProjectPanel", "--compat", "java",
                        "--seed", "1"], capture_output=True, text=True, cwd=ROOT, timeout=120)
    assert r.returncode == 0
    assert _normalize(r.stdout) == _normalize(_cli("run", "ProjectPanel", "--compat", "java", "--seed", "1")[1])
