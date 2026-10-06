"""Gap tracking (shadowinfo.polygon.gaps) and gap-history events (shadowinfo.polygon.events) vs. the Java.

Golden data: tests/fixtures/java (tools/java_reference; docs/notes/original_java.md §5.2 items 5, 6, 9):

5. ``getGaps`` -- gap sets with IDs, iEdges, links (Java HashSet order), relative times and the printed
   ``ProjectPanel`` lines on every run where Java succeeds; the exception, critical point and source
   line on the four runs where it throws;
6. never-see-evader states (``processGapHistoryInformation(NEVER_SEE_EVADER)``);
9. the paper runs: polygon 12 = ICRA'08 Fig. 4 / T-RO Fig. 11, polygon 13 = T-RO Fig. 15(a) (85 events,
   18 bounded final shadows), polygon 14 = T-RO Fig. 15(b) (385 events, 491 IDs, 12 bounded final shadows).

``compat="java"`` must reproduce all of it.  ``compat="fixed"`` must give valid histories everywhere: by
default (``java_matching=True``) identical to Java wherever Java succeeds; with the opt-in
``java_matching=False`` identical except at the steps where Java's ``samePhysicalGap`` matching rotates the
labels (pinned below); and on the crashing runs a history whose shadows agree with the physical gaps at
every sample.
"""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from pathlib import Path as FsPath

import pytest

from shadowinfo.bipartite import BipartiteIState
from shadowinfo.events import Appear, Disappear, Merge, ShadowSequence, Split, event_labels
from shadowinfo.polygon import (CutType, Gap, GapTrackingError, NSEState, PhysicalGap, chain_graph, format_gap_set,
                                gap_history_to_events, get_gaps, java_hashset_order, load_path, load_polygon,
                                never_see_evader, never_see_evader_via_filter, sample_path, transition_events)
from shadowinfo.polygon import gaps as G
from shadowinfo.polygon.events import nse_string

ROOT = FsPath(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "java"
sys.path.insert(0, str(ROOT / "tools" / "java_reference"))

import convert  # noqa: E402

RUNS = sorted(p.stem for p in FIX.glob("*_*.json"))
BIG = {5, 13, 14}


def _load(name):
    return json.loads((FIX / f"{name}.json").read_text())


FIXTURES = {name: _load(name) for name in RUNS}
OK_RUNS = [r for r in RUNS if FIXTURES[r]["get_gaps"]["ok"]]
CRASH_RUNS = [r for r in RUNS if not FIXTURES[r]["get_gaps"]["ok"]]

# Steps where Java's samePhysicalGap match of pgs[k][0] is ambiguous (two adjacent small gaps within one
# edge), the last candidate wins and every label moves one gap along the boundary (quirk B3).  The opt-in
# get_gaps(compat="fixed", java_matching=False) (the simulator's tracker) re-derives exactly these steps from
# the hidden chains; the default keeps Java's labels.  Not rare: on random valid paths in the 14 maps about
# 18% of the runs have such a step, and about 12% end with different final shadow IDs.
ROTATED = {
    "14_PP5": [124, 125, 126, 127, 128, 129, 130, 219, 220, 221, 222, 223],
    "14_fig_TRO-Fig15b": [124, 125, 126, 127, 128, 129, 221, 222, 223, 224, 225],
    "5_PP3": [123],
    "5_c_stub": [30],
}


def _param(names):
    return [pytest.param(n, marks=pytest.mark.slow) if int(n.split("_")[0]) in BIG else n for n in names]


_SAMPLES: dict = {}


def samples(name):
    """Critical points and physical gaps of a fixture run (mode independent: no 45-degree segments)."""
    if name not in _SAMPLES:
        poly_n, label = name.split("_", 1)
        poly, path = load_path(f"{label}@{poly_n}")
        s = sample_path(poly, path, "java")
        assert s.mode_independent and not s.errors
        _SAMPLES[name] = s
    return _SAMPLES[name]


def history(name, compat="java", **kw):
    s = samples(name)
    return get_gaps(s.polygon, s.path, compat, samples=s, **kw)


def as_lists(h):
    return [[g.to_list() for g in gs] for gs in h]


def java_pairs(h, i):
    """Java's step i-1 -> i as (index in sample i-1, index in sample i) pairs: continuations and links."""
    links = {g.id: set(g.to_gaps) for gs in h for g in gs}
    a, b = h.sample_ids[i - 1], h.sample_ids[i]
    return {(j, k) for j, x in enumerate(a) for k, y in enumerate(b) if x == y or y in links.get(x, ())}


def check_samples_agree(h):
    """Geometric sanity: replaying the events, the shadows alive at every sample are exactly the IDs
    given to the physical gaps there (so their number is the number of physical gaps)."""
    ge = gap_history_to_events(h)
    ge.sequence()  # validates
    per = transition_events(h)
    assert sum(map(len, per)) == len(ge.events)
    alive = set(ge.initial)
    k = 0
    for i, ids in enumerate(h.sample_ids):
        while k < len(h.sample_indices) and h.sample_indices[k] < i:
            for e in per[k]:
                gone, new = event_labels(e)
                alive = (alive - set(gone)) | set(new)
            k += 1
        assert len(ids) == len(set(ids)) == len(h.samples.physical[i]), i
        assert set(ids) == alive, (i, sorted(ids), sorted(alive))
    for k, gs in enumerate(h):  # each set is the last sample before an event
        assert [g.id for g in gs] == h.sample_ids[h.sample_indices[k]]
        assert [g.physical for g in gs] == h.samples.physical[h.sample_indices[k]]


# --------------------------------------------------------------------------- Java helpers


# java.util.HashSet<Integer> iteration orders printed by JDK 17 (scratch check: 405 random insertion
# sequences, all equal to JavaIntHashSet)
JAVA_HASHSETS = [
    ([15, 16], [16, 15]),
    ([31, 17, 1], [17, 1, 31]),
    ([65539, 3], [65539, 3]),
    ([-1, 5, -17], [-1, -17, 5]),
    # one bucket: the 9th entry resizes the table (treeifyBin below capacity 64), then size > 0.75 cap
    ([20, 4, 36, 52, 68, 84, 100, 116, 132, 148, 164, 180, 196],
     [4, 36, 68, 100, 132, 164, 196, 20, 52, 84, 116, 148, 180]),
]


@pytest.mark.parametrize("ids, order", JAVA_HASHSETS)
def test_java_hashset_order(ids, order):
    assert java_hashset_order(ids) == order


def test_java_hashset_basics():
    assert java_hashset_order([6, 7]) == [6, 7]
    assert java_hashset_order(range(12, 0, -1)) == list(range(1, 13))
    assert java_hashset_order(range(40, 0, -1)) == list(range(1, 41))  # 64 buckets
    assert java_hashset_order([5, 5, 21]) == [5, 21]
    s = G.JavaIntHashSet([21, 5])
    assert list(s) == [21, 5] and 5 in s and len(s) == 2 and not s.add(21)


def test_link_order_in_fixtures_is_java_hashset_order():
    """Every toGapSet Java printed equals the HashSet order of its IDs inserted in ascending order."""
    n = 0
    for name in OK_RUNS:
        for st in FIXTURES[name]["get_gaps"]["sets"]:
            for _gid, _ie, to in st["gaps"]:
                if to:
                    assert to == java_hashset_order(sorted(to)), (name, to)
                    n += 1
    assert n > 300


def test_gap_helpers():
    n = 10
    g = PhysicalGap(3, 6, (0.0, 0.0), None)  # from a point on edge 3 to vertex 7
    assert (g.full_start_edge(n), g.full_end_edge(n)) == (4, 6)
    assert [G.vertex_in_gap(v, g, n) for v in range(n)] == [v in (4, 5, 6) for v in range(n)]
    assert [G.line_in_gap(e, g, n) for e in range(n)] == [e in (4, 5, 6) for e in range(n)]
    w = PhysicalGap(8, 1, None, (0.0, 0.0))  # wraps: from vertex 8 to a point on edge 1
    assert (w.full_start_edge(n), w.full_end_edge(n)) == (8, 0)
    assert [G.line_in_gap(e, w, n) for e in range(n)] == [e in (8, 9, 0) for e in range(n)]
    assert [G.vertex_in_gap(v, w, n) for v in range(n)] == [v in (9, 0, 1) for v in range(n)]
    # B4: Java's wrap test accepts any start/end edges summing to n-1
    a, b = PhysicalGap(2, 4, None, (0, 0)), PhysicalGap(7, 5, None, (0, 0))
    assert G.same_physical_gap(a, b, n, "java") and not G.same_physical_gap(a, b, n, "fixed")
    c, d = PhysicalGap(0, 2, None, (0, 0)), PhysicalGap(9, 3, None, (0, 0))
    assert G.same_physical_gap(c, d, n, "java") and G.same_physical_gap(c, d, n, "fixed")


def test_gap_to_string():
    g = Gap(PhysicalGap(0, 1, None, (0, 0)), 19, 0)
    assert str(g) == "[19, 0]"
    g.add_gap(16)
    g.add_gap(15)
    assert str(g) == "[19, 0 | 16 15 ]" and g.to_gaps == (16, 15) and g.to_list() == [19, 0, [16, 15]]
    g.state = NSEState.CONTAMINATED
    assert str(g) == "[19, 0| 1 | 16 15 ]"
    h = Gap(None, 3, 7)
    h.state = NSEState.CLEAR
    assert format_gap_set([h], with_time=False) == "[3, 7| 0] "
    assert format_gap_set([h], time=1e-4) == "1.0E-4 [3, 7| 0] "


def test_collapse_physical_gaps_crashes_like_java():
    n = 20

    def gs(*edges):
        return [Gap(PhysicalGap(s, e, None, (0.0, 0.0))) for s, e in edges]

    a = gs((1, 3), (6, 8), (12, 14))
    b = gs((6, 8), (12, 14), (1, 3))
    gen = G.IDGenerator()
    out = G.collapse_physical_gaps([a, b], gen, n)
    assert [g.id for g in a] == [1, 2, 3] and [g.id for g in out] == [2, 3, 1] and out is b
    assert [g.iedge for g in out] == [6, 12, 1]  # (fullStart + fullEnd) // 2, the end edge is partial
    with pytest.raises(G._JavaThrow) as ex:  # fewer gaps in the next sample, no match before the end
        G.collapse_physical_gaps([a, gs((16, 17), (9, 10))], gen, n)
    assert (ex.value.exc, ex.value.line) == ("java.lang.ArrayIndexOutOfBoundsException", 380)
    assert ex.value.message == "Index 2 out of bounds for length 2"
    with pytest.raises(G._JavaThrow) as ex:  # fewer gaps, a match, then the shift runs off the end
        G.collapse_physical_gaps([a, gs((1, 3), (6, 8))], gen, n)
    assert (ex.value.exc, ex.value.line) == ("java.lang.ArrayIndexOutOfBoundsException", 390)
    with pytest.raises(G._JavaThrow) as ex:  # one gap more: it gets no ID
        G.collapse_physical_gaps([a, gs((1, 3), (6, 8), (12, 14), (16, 17))], gen, n)
    assert (ex.value.exc, ex.value.line) == ("java.lang.NullPointerException", 460)
    with pytest.raises(G._JavaThrow) as ex:  # nothing matches
        G.collapse_physical_gaps([a, gs((4, 5), (9, 10), (16, 17))], gen, n)
    assert (ex.value.exc, ex.value.line) == ("java.lang.NullPointerException", 411)
    with pytest.raises(G._JavaThrow) as ex:  # no IDs and no generator
        G.collapse_physical_gaps([gs((1, 3))], None, n)
    assert ex.value.line == 371
    assert G.collapse_physical_gaps([[], b], gen, n) == []


# --------------------------------------------------------------------------- compat="java" vs. the Java


@pytest.mark.parametrize("name", _param(OK_RUNS))
def test_get_gaps_java_matches_golden(name):
    gg = FIXTURES[name]["get_gaps"]
    h = history(name)
    assert as_lists(h) == [st["gaps"] for st in gg["sets"]]
    assert h.times == [st["t"] for st in gg["sets"]]
    assert [gs[0].relative_time for gs in h] == h.times  # only element 0 carries the time
    assert all(g.relative_time == 0.0 for gs in h for g in gs[1:])
    assert h.sample_indices == [st["sample_index"] for st in gg["sets"]]
    assert h.lines() == [st["line"] for st in gg["sets"]]  # ProjectPanel / ProjectPanel5 stdout
    assert (h.final_ids, h.max_id, len(h)) == (gg["final_ids"], gg["max_id"], gg["n_sets"])
    assert len(h) == 1 + FIXTURES[name]["n_gi_bt_crossings"]
    check_samples_agree(h)


@pytest.mark.parametrize("name", _param(CRASH_RUNS))
def test_get_gaps_java_throws_where_java_throws(name):
    gg = FIXTURES[name]["get_gaps"]
    with pytest.raises(GapTrackingError) as ex:
        history(name)
    e = ex.value
    fl = gg["failure"]
    assert (e.critical_index, e.critical_type.name, e.getgaps_line) == \
        (fl["critical_index"], fl["critical_type"], fl["getGaps_line"])
    assert e.java_exception == gg["exception"]["class"]
    assert f"collapsePhysicalGaps:{e.java_line}" == fl["innermost_frame"].split(".")[-1]
    if "Index" in gg["exception"]["message"]:
        assert gg["exception"]["message"] in str(e)
    assert e.phase.startswith(fl["phase"].split(" of ")[0])


def test_java_mode_reports_b1_and_empty_sets():
    poly, _ = load_path("P12")
    path45 = G.Path([(100.0, 460.0), (150.0, 510.0)])  # slope exactly 1: purturbPointAlongSeg NPE (B1)
    with pytest.raises(GapTrackingError) as ex:
        get_gaps(poly, path45, "java")
    assert ex.value.critical_index == 0 and ex.value.java_exception.endswith("NullPointerException")
    h = get_gaps(poly, path45, "fixed")
    check_samples_agree(h)
    square = G.Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])  # no gaps: previousGaps[0] throws in Java
    p = G.Path([(2.0, 2.0), (7.0, 3.0)])
    with pytest.raises(GapTrackingError) as ex:
        get_gaps(square, p, "java")
    assert ex.value.getgaps_line == 307 and "out of bounds for length 0" in str(ex.value)
    h = get_gaps(square, p)
    assert list(h) == [[]] and h.times == [0.0] and gap_history_to_events(h).events == []
    with pytest.raises(ValueError):
        get_gaps(poly, path45, "java", samples=sample_path(poly, path45, "fixed"))


# --------------------------------------------------------------------------- compat="fixed"


@pytest.mark.parametrize("name", _param(OK_RUNS))
def test_fixed_is_identical_to_java(name):
    """Settled requirement: the default compat="fixed" gives Java's gap history wherever Java succeeds."""
    hj = history(name)
    for h in (history(name, "fixed"), history(name, "fixed", java_matching=True)):
        assert h.inferred == []
        assert as_lists(h) == as_lists(hj) and h.times == hj.times and h.lines() == hj.lines()
        assert h.sample_indices == hj.sample_indices and h.sample_ids == hj.sample_ids
        assert gap_history_to_events(h) == gap_history_to_events(hj)


# Live Java (tools/java_reference GoldenDump) on a random valid path in map 7, where the chain matching
# rotates Java's labels at critical point 1 and changes the final shadow IDs ({4, 7, 8} vs {5, 7, 8}).
RANDOM_7 = [(368.24264464516, 517.7757563437503), (43.90993696296563, 77.31204259753078),
            (111.56240747077494, 40.59598986746982)]
RANDOM_7_JAVA = [
    "0.0 [1, 0] [2, 2] [3, 7] [4, 10] [5, 14] [6, 17] ",
    "0.02721444906993155 [5, 0] [1, 2] [2, 7] [3, 10] [4, 14] ",
    "0.06603092785001376 [5, 0 | 7 ] [1, 2 | 7 ] [3, 11] [4, 14] ",
    "0.28960638696324437 [3, 11] [4, 14] [7, 1] ",
    "0.34924228623071846 [4, 14] [7, 1] ",
    "0.39865361225062945 [8, 6] [4, 14] [7, 1] ",
]


def test_fixed_default_matches_live_java_on_a_random_path():
    poly = load_polygon(7)
    path = G.Path(RANDOM_7)
    s = sample_path(poly, path)
    h = get_gaps(poly, path, samples=s)
    assert h.lines() == RANDOM_7_JAVA == get_gaps(poly, path, "java", samples=s).lines()
    assert h.inferred == [] and sorted(h.final_ids) == [4, 7, 8]
    assert gap_history_to_events(h).events == [Disappear(6, 0, 0), Disappear(2, 0, 0), Merge(5, 1, 7),
                                               Disappear(3, 0, 0), Appear(8, 0, 0)]
    hc = get_gaps(poly, path, samples=s, java_matching=False)  # opt-in: labels stay on their shadows
    assert [i for i, _ in hc.inferred] == [1] and "disagrees with the hidden chains" in hc.inferred[0][1]
    assert sorted(hc.final_ids) == [5, 7, 8]
    check_samples_agree(hc)


# Synthetic polygon 124 of the review differential (live Java, GoldenDump): the path retraces part of its
# previous segment, so Java's ptDistFromStart dates the crossings on the third segment by the second (B17).
POLY_124 = [
    (968.4732354806313, 517.809972188041), (808.7094571156526, 694.5376721993525),
    (793.176421133058, 859.295800489348), (430.4646928512673, 657.351907972535),
    (453.13880512756026, 602.0426731376547), (451.28623229840406, 591.6752637032869),
    (450.46560801086775, 589.8516836593719), (248.4914457531979, 717.0399325008848),
    (315.1932933263803, 640.0907545535463), (170.44775683157388, 671.7758558421967),
    (319.8107192125494, 540.6955271330625), (288.9002935827283, 544.2826586116165),
    (417.80035676885814, 509.68180908837047), (347.00238069117546, 499.86795430741194),
    (197.6014950722273, 493.5869673835119), (179.6751215217945, 402.67776597000574),
    (200.07711113266345, 392.5733219678857), (212.08520298102536, 323.6547128577861),
    (346.0672680925676, 360.6740498137657), (246.2863710388393, 224.28281991089005),
    (458.00776438910606, 358.652218080578), (460.05442462397457, 52.10586651154779),
    (515.1738085451282, 227.93700891540118), (546.7281811252379, 275.0496084997412),
    (612.2140573701993, 98.47935538948519), (628.1639233441944, 56.10566613636098),
    (634.5856068435317, 42.823520363035016), (762.169222568371, 457.2461411524895),
    (759.2143359443986, 462.90490428322664),
]
RETRACE_JAVA_TIMES = [0.0, 0.0026595690367123346, 0.06513413591174319, 0.1588247752603057, 0.17017564467364893,
                      0.3654847611458886, 0.9605408788210644, 0.7985971433774186, 0.7508973810878307,
                      0.708409049602836, 0.708409049602836, 0.7508973810878307, 0.7985971433774186,
                      0.9605408788210644]


def test_fixed_equals_java_on_a_retracing_path_except_the_b17_times():
    poly = G.Polygon(POLY_124)
    path = G.Path([(888.0, 521.0), (853.0, 591.0), (463.0, 607.0), (658.0, 599.0)])
    hj, h = get_gaps(poly, path, "java"), get_gaps(poly, path)
    assert hj.times == RETRACE_JAVA_TIMES  # repeated and decreasing: Java's TreeMap-keyed oracle drops events (B9)
    assert h.times == sorted(h.times) and len(set(h.times)) == len(h.times) and h.times[:6] == hj.times[:6]
    assert as_lists(h) == as_lists(hj) and h.sample_ids == hj.sample_ids and h.inferred == []
    assert gap_history_to_events(h).events == gap_history_to_events(hj).events


@pytest.mark.parametrize("name", _param(OK_RUNS))
def test_chain_matching_equals_java_except_where_java_rotates_labels(name):
    hj, h = history(name), history(name, "fixed", java_matching=False)
    rotated = ROTATED.get(name, [])
    assert [i for i, _ in h.inferred] == rotated
    assert all("disagrees with the hidden chains" in why for _, why in h.inferred)
    # the hidden-chain graph is exactly Java's step everywhere else
    for i in range(1, len(hj.sample_ids)):
        same = java_pairs(hj, i) == set(chain_graph(hj.samples.polygon, hj.samples.physical[i - 1],
                                                    hj.samples.physical[i]))
        assert same == (i not in rotated), i
    if not rotated:
        assert as_lists(h) == as_lists(hj) and h.lines() == hj.lines()
    else:
        first = rotated[0]
        before = [k for k, i in enumerate(hj.sample_indices) if i < first - 1]
        assert [as_lists(h)[k] for k in before] == [as_lists(hj)[k] for k in before]
        assert h.sample_ids[:first] == hj.sample_ids[:first] and h.sample_ids[first] != hj.sample_ids[first]
        for i in rotated:
            assert _ambiguous_java_match(h, i), i
    # same structure: number of sets, set times and sizes, event types, IDs used, final shadows
    assert h.times == hj.times and h.sample_indices == hj.sample_indices
    assert [len(gs) for gs in h] == [len(gs) for gs in hj]
    ej, ef = gap_history_to_events(hj), gap_history_to_events(h)
    assert [type(e) for e in ef.events] == [type(e) for e in ej.events] and ef.times == ej.times
    assert (h.max_id, sorted(h.final_ids)) == (hj.max_id, sorted(hj.final_ids))
    check_samples_agree(h)


def _ambiguous_java_match(h, i):
    """Redo Java's step i-1 -> i on the fixed history and report whether a ``samePhysicalGap`` match of
    ``pgs[k][0]`` had several candidates (Java keeps the last one)."""
    seen = []
    orig = G.collapse_physical_gaps

    def spy(pgs, id_gen, n):
        for a, b in zip(pgs, pgs[1:]):
            if a:
                seen.append(sum(G.same_physical_gap(a[0], g, n, "java") for g in b[:len(pgs[0])]))
        return orig(pgs, id_gen, n)

    s = h.samples
    prev = [Gap(p, gid) for p, gid in zip(s.physical[i - 1], h.sample_ids[i - 1])]
    G.collapse_physical_gaps = spy
    try:
        G._java_proposal(s.critical_points[i], i, s, prev, s.polygon.n, [s.polygon.edge(k) for k in range(s.polygon.n)],
                         10 ** 6)
    finally:
        G.collapse_physical_gaps = orig
    return any(c >= 2 for c in seen)


@pytest.mark.parametrize("name", _param(CRASH_RUNS))
@pytest.mark.parametrize("java_matching", [False, True])
def test_fixed_repairs_the_crashing_runs(name, java_matching):
    crash = FIXTURES[name]["get_gaps"]["failure"]["critical_index"]
    h = history(name, "fixed", java_matching=java_matching)
    steps = [i for i, _ in h.inferred]
    assert crash in steps
    assert any(i == crash and why.startswith("Java throws") for i, why in h.inferred)
    check_samples_agree(h)
    ge = gap_history_to_events(h)
    seq = ge.sequence()
    assert seq.alive_at_end() == sorted(h.final_ids)
    assert len(ge.times) == len(ge.events) and ge.times == sorted(ge.times)
    nse = never_see_evader(h)
    assert [nse_string(x) for x in nse] == [nse_string(x) for x in never_see_evader_via_filter(h)]


def test_fixed_repairs_known_crashes_physically():
    """At 14_PP4's critical point 654 the waypoint (633.93, 714.29) lies on a bitangent: in one step two
    shadows merge and one appears (ProjectPanel5 moved the waypoint by 2 units).  In the digitised
    ICRA'08 Fig. 8(b) path the 0.005 perturbation at critical point 596 overshoots the bitangent
    crossing 0.0013 further on, so the sample shows the merge and the appear together."""
    h = history("14_PP4", "fixed", java_matching=True)
    assert h.inferred == [(654, "Java throws; from the hidden chains: 476 -> [484], 470 -> [484], appear 483")]
    k = h.sample_indices.index(653)
    assert transition_events(h)[k] == [Merge(476, 470, 484), Appear(483, 0, 0)]
    h = history("14_fig_ICRA08-Fig8b", "fixed", java_matching=True)
    assert [i for i, _ in h.inferred] == [596, 597]
    k = h.sample_indices.index(595)
    assert transition_events(h)[k] == [Merge(448, 446, 451), Appear(452, 0, 0)]
    # the bitangent at 597 was already crossed by the sample of 596: no event, no extra gap set
    assert h.sample_indices[k + 1] == 597 and 596 not in h.sample_indices


def test_fixed_relative_time_ignores_far_segments():
    """B17: ``ptDistFromStart`` keeps the last segment within eps, which on a self-crossing path dates a
    crossing on the first leg by the last leg; ``compat="fixed"`` keeps the segment it was found on."""
    path = G.Path([(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (5.0, 10.0), (5.0, -5.0)])
    p = (5.0, 0.0)
    cp = G.CriticalPoint(5.0, p, CutType.SINGLETANGENT, 0, 0, path.segments[0])
    length = path.length()
    assert G._relative_time(path, cp, length, "java", 5e-5) == 35.0 / length
    assert G._relative_time(path, cp, length, "fixed", 5e-5) == 5.0 / length
    end = G.CriticalPoint(10.0, (10.0, 0.0), CutType.NONE, -1, 0, path.segments[0])  # a waypoint
    assert G._relative_time(path, end, length, "java", 5e-5) == G._relative_time(path, end, length, "fixed", 5e-5)


# --------------------------------------------------------------------------- events and NSE


@pytest.mark.parametrize("name", _param(OK_RUNS))
def test_events_and_nse_match_java(name):
    run = FIXTURES[name]
    hj = history(name)
    ge = gap_history_to_events(hj)
    init, per, anomalies = convert.gap_history_events(run)
    assert anomalies == []
    assert ge.initial == init and ge.events == [e for evs in per for e in evs]
    assert transition_events(hj) == per and all(len(evs) == 1 for evs in per)
    assert ge.times == hj.times[1:]
    seq = ge.sequence()
    assert seq.alive_at_end() == sorted(run["get_gaps"]["final_ids"])
    states = never_see_evader(hj, "java")
    assert [nse_string(s) for s in states] == run["nse"]["states"]
    assert [nse_string(g.state for g in gs) for gs in hj] == run["nse"]["states"]
    assert [nse_string(s) for s in never_see_evader(hj)] == run["nse"]["states"]
    assert [nse_string(s) for s in never_see_evader_via_filter(hj)] == run["nse"]["states"]
    # ProjectPanel2/3 print the sets after NSE: [id, iEdge| s | links ]
    want = ["".join(f"[{i}, {ie}| {c}" + (" | " + "".join(f"{x} " for x in to) if to else "") + "] "
                    for (i, ie, to), c in zip(st["gaps"], cs))
            for st, cs in zip(run["get_gaps"]["sets"], run["nse"]["states"])]
    for gs in hj:  # Java's NSE does not reset states
        for g in gs:
            g.state = None
    never_see_evader(hj, "java")
    assert hj.lines(with_time=False) == want
    hf = history(name, "fixed")
    assert [nse_string(s) for s in never_see_evader(hf)] == [nse_string(s) for s in never_see_evader_via_filter(hf)]


def _hist(sets):
    """A gap history from ``[[(id, [links]), ...], ...]``."""
    out = []
    for st in sets:
        gs = []
        for gid, to in st:
            g = Gap(None, gid, 0)
            for x in to:
                g.add_gap(x)
            gs.append(g)
        out.append(gs)
    return out


def test_transition_events_general_components():
    h = _hist([[(1, [4, 5, 6]), (2, []), (3, [])],
               [(4, [9]), (5, [9]), (6, [9]), (2, [7, 8]), (3, [10, 11]), (12, [10, 11])],
               [(9, []), (7, []), (8, []), (10, []), (11, []), (13, [])]])
    per = transition_events(h)
    assert per[0] == [Split(1, 4, 14), Split(14, 5, 6), Appear(12, 0, 0)]
    assert per[1] == [Merge(4, 5, 15), Merge(15, 6, 9), Split(2, 7, 8), Merge(3, 12, 16), Split(16, 10, 11),
                      Appear(13, 0, 0)]
    ge = gap_history_to_events(h)
    seq = ge.sequence()
    assert seq.alive_at_end() == [7, 8, 9, 10, 11, 13]
    # gap 5 joins two link groups that looked separate: one component, merged then split
    h = _hist([[(1, [3]), (2, [4]), (5, [3, 4])], [(3, []), (4, [])]])
    assert transition_events(h) == [[Merge(1, 2, 6), Merge(6, 5, 7), Split(7, 3, 4)]]
    with pytest.raises(ValueError, match="alone"):
        transition_events(_hist([[(1, [2])], [(2, [])]]))
    with pytest.raises(ValueError, match="survives"):
        transition_events(_hist([[(1, [2, 3])], [(1, []), (2, []), (3, [])]]))
    with pytest.raises(ValueError, match="not in set"):
        transition_events(_hist([[(1, [2, 3])], [(2, [])]]))


def test_never_see_evader_b6_and_states():
    # 1 splits into three, then 4 and 5 merge.  Java ignores |toGapSet| > 2 (B6), so 2, 3, 4 become clear.
    h = _hist([[(1, [2, 3, 4]), (5, [])], [(2, []), (3, []), (4, [6]), (5, [6])], [(2, []), (3, []), (6, [])]])
    assert [nse_string(s) for s in never_see_evader(h, "java")] == ["11", "0001", "001"]
    assert format_gap_set(h[1], with_time=False) == "[2, 0| 0] [3, 0| 0] [4, 0| 0 | 6 ] [5, 0| 1 | 6 ] "
    fixed = [nse_string(s) for s in never_see_evader(h)]
    assert fixed == ["11", "1111", "111"]
    assert fixed == [nse_string(s) for s in never_see_evader_via_filter(h)]
    # a clear part does not clear a merge with a contaminated part; a new gap is clear
    h = _hist([[(1, [])], [(1, [3]), (2, [3])], [(3, [])]])
    assert [nse_string(s) for s in never_see_evader(h, "java")] == ["1", "10", "1"]
    assert [nse_string(s) for s in never_see_evader(h)] == ["1", "10", "1"]
    assert str(NSEState.CLEAR) == "| 0" and str(NSEState.CONTAMINATED) == "| 1"
    assert never_see_evader([]) == [] and nse_string([None, NSEState.CLEAR]) == "-0"


def test_never_see_evader_with_an_empty_start_set():
    """The start point of an L room sees no gap (Java: getGaps throws at line 307), the fixed tracker gives
    an empty set 0 and then gap 1 appears.  Nothing can be contaminated: the fixed labels are clear, as the
    counting filter says; Java's pursuerNeverSeeEvader returns at once and leaves the states unset."""
    room = G.Polygon([(0, 0), (10, 0), (10, 4), (4, 4), (4, 10), (0, 10)])
    path = G.Path([(2.0, 2.0), (9.0, 2.0)])
    with pytest.raises(GapTrackingError) as ex:
        get_gaps(room, path, "java")
    assert ex.value.getgaps_line == 307
    h = get_gaps(room, path)
    assert h.lines() == ["0.0 ", "0.2857142857142857 [1, 3] "]
    for g in h[1]:
        g.state = NSEState.CONTAMINATED  # a stale state from an earlier call must be reset
    assert [nse_string(s) for s in never_see_evader(h)] == ["", "0"]
    assert [nse_string(s) for s in never_see_evader_via_filter(h)] == ["", "0"]
    assert format_gap_set(h[1], with_time=False) == "[1, 3| 0] "
    for g in h[1]:
        g.state = None
    assert [nse_string(s) for s in never_see_evader(h, "java")] == ["", "-"]


# --------------------------------------------------------------------------- the paper runs (item 9)


def _bounded_finals(h):
    """Final shadows the Java prints bounds for: all but those that appeared and never split or merged (B11)."""
    appeared = {e.s for e in gap_history_to_events(h).events if isinstance(e, Appear)}
    return [s for s in h.final_ids if s not in appeared]


@pytest.mark.parametrize("compat", ["java", "fixed"])
def test_paper_icra08_fig4_and_tro_fig11(compat):
    """P12 is ICRA'08 Fig. 2/4: labels 1-19.  T-RO Fig. 11 swaps 6 and 7 and draws s5 as initial."""
    poly, path = load_path("P12")
    h = get_gaps(poly, path, compat)
    ge = gap_history_to_events(h)
    assert ge.initial == [1, 2, 3, 4] and ge.events[0] == Appear(5, 0, 0)
    swap = {6: 7, 7: 6}
    r = lambda x: swap.get(x, x)  # noqa: E731
    ev = []
    for e in ge.events[1:]:
        if isinstance(e, Split):
            ev.append(Split(r(e.s), r(e.a), r(e.b)))
        elif isinstance(e, Merge):
            ev.append(Merge(r(e.a), r(e.b), r(e.s)))
        else:
            ev.append(type(e)(r(e.s), 0, 0))
    fig = json.loads((ROOT / "tests" / "fixtures" / "tor_fig11.json").read_text())
    paper = ShadowSequence.from_dict(fig["sequence"])

    def shape(e):
        if isinstance(e, Split):
            return ("split", e.s, frozenset((e.a, e.b)))
        if isinstance(e, Merge):
            return ("merge", frozenset((e.a, e.b)), e.s)
        return (type(e).__name__, e.s)

    assert Counter(map(shape, ev)) == Counter(map(shape, paper.events))
    seq = ShadowSequence({s: (0, math.inf) for s in [1, 2, 3, 4, 5]}, ev)
    edges = sorted([u, v] for u, (_t, v) in BipartiteIState.from_sequence(seq).edges())
    assert edges == sorted(fig["bipartite"]["edges"])
    assert sorted(h.final_ids) == [13, 15, 18, 19] and _bounded_finals(h) == [19, 15, 13]


@pytest.mark.slow
@pytest.mark.parametrize("compat", ["java", "fixed"])
def test_paper_tro_fig15(compat):
    """T-RO Fig. 15(a): 85 component events, 18 final shadows with bounds.  Fig. 15(b): 385 events, 491
    shadows in total, 12 final shadows with bounds (ProjectPanel5 path, both modes)."""
    a = history("13_c_first", compat)  # the stored path P13 (commented '// polygon 13' block, 9 points)
    assert a.samples.path.points == load_path("P13")[1].points
    ea = gap_history_to_events(a)
    assert len(ea.events) == len(a) - 1 == 85
    assert len(a.final_ids) == 20 and len(_bounded_finals(a)) == 18
    b = history("14_PP5", compat)
    assert b.samples.path.points == load_path("P14b")[1].points
    eb = gap_history_to_events(b)
    assert len(eb.events) == len(b) - 1 == 385 and b.max_id == 491
    assert len(b.final_ids) == 14 and len(_bounded_finals(b)) == 12
    for h, e in ((a, ea), (b, eb)):  # every ID is an initial shadow or created by exactly one event
        assert h.max_id == len(e.initial) + sum(len(event_labels(x)[1]) for x in e.events)
        assert h.times[1:] == e.times  # one event per set transition, at the set's time
        opened = [h.samples.critical_points[i + 1].cut_type for i in h.sample_indices[:-1]]
        assert set(opened) == {CutType.GENERAL_INFLECTION, CutType.BITANGENT}
