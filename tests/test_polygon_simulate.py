"""Polygon-world simulator (shadowinfo.polygon.simulate): robot on a path, random-walk targets, pockets, events.

Checked here:

* the JS-mirrorable predicates (scalar == vectorised; closed visibility on hand-made cases);
* shadow pockets: inside the polygon, disjoint, and with the visibility polygon they tile it (area);
  a point is hidden iff it lies in exactly one pocket;
* the default demo path of every map (``data/paths.json`` section ``demo``): valid, and it produces every
  kind of component event;
* the component events and labels of a run are exactly those of
  ``get_gaps(compat="fixed", java_matching=False)`` on the same path;
* tick semantics (docs/DESIGN.md §5): replaying a run's events on the true counts reproduces the ground
  truth after every tick, on all 14 maps and many seeds; the event stream is a valid
  :class:`~shadowinfo.events.ShadowSequence` and the combinatorial filter's bounds contain the truth;
* loops, back-and-forth runs and click-to-move (``go_to`` / ``follow``).
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pytest

from shadowinfo.events import Appear, Disappear, Enter, Exit, Merge, ShadowSequence, Split, event_from_dict
from shadowinfo.nondeterministic import CombinatorialFilter
from shadowinfo.polygon import (EPSILON, PolygonSimulator, Shadow, demo_paths, gap_history_to_events, get_gaps,
                                load_demo_path, load_path, load_paths, load_polygon, path_record, physical_gaps,
                                sample_path, seg_in_polygon, shadow_pockets, simulate, validate_path,
                                visibility_polygon)
from shadowinfo.polygon.geometry import Path, Polygon, pt_seg_dist
from shadowinfo.polygon.simulate import (_on_boundary, _track, point_in_polygon, points_in_polygon, polygon_area,
                                         segment_in_polygon, segment_touches_boundary, segments_in_polygon)
from shadowinfo.rng import Rng

MAPS = list(range(1, 15))
BIG = {5, 13, 14}
LONG = {2, 4, 9, 10, 11}  # demo tours of 4900-9000 units
COMPONENT = ("appear", "disappear", "split", "merge")


def _maps(slow_maps=BIG):
    return [pytest.param(n, marks=pytest.mark.slow) if n in slow_maps else n for n in MAPS]


def _random_points(P, k, rng):
    x0, x1, y0, y1 = P.vx.min(), P.vx.max(), P.vy.min(), P.vy.max()
    X, Y = [], []
    while len(X) < k:
        x, y = x0 + rng.random() * (x1 - x0), y0 + rng.random() * (y1 - y0)
        if point_in_polygon(P, (x, y)):
            X.append(x)
            Y.append(y)
    return np.array(X), np.array(Y)


# --------------------------------------------------------------------------- predicates


def test_predicates_scalar_equals_vectorised():
    rng = Rng(3)
    for n in (2, 5, 12, 14):
        P = load_polygon(n)
        X = np.array([rng.random() * 1000 for _ in range(400)])
        Y = np.array([rng.random() * 1000 for _ in range(400)])
        inside = points_in_polygon(P, X, Y)
        assert inside.tolist() == [point_in_polygon(P, (x, y)) for x, y in zip(X, Y)]
        assert inside.any() and not inside.all()
        Xi, Yi = X[inside], Y[inside]
        for q in [(Xi[0], Yi[0]), P.vertices[3], ((P.vertices[0][0] + P.vertices[1][0]) / 2, P.vertices[0][1])]:
            assert segments_in_polygon(P, q, Xi, Yi).tolist() == [segment_in_polygon(P, q, (x, y))
                                                                  for x, y in zip(Xi, Yi)]


def test_closed_visibility_cases():
    # an L-shaped room (CCW): reflex vertex at (10, 10)
    L = Polygon([(0, 0), (20, 0), (20, 10), (10, 10), (10, 20), (0, 20)], name="L")
    assert L.signed_area() > 0
    assert segment_in_polygon(L, (5, 5), (15, 5))
    assert segment_in_polygon(L, (15, 5), (5, 15))  # grazes the reflex vertex (10, 10): visible (closed)
    assert not segment_in_polygon(L, (15, 6), (6, 15))  # cuts the corner outside
    assert segment_in_polygon(L, (12, 10), (18, 10))  # along an edge
    assert segment_in_polygon(L, (0, 5), (0, 15))  # along the outer wall
    # through a convex vertex the segment leaves the polygon
    T = Polygon([(0, 0), (10, 0), (5, 5)], name="T")
    assert not segment_in_polygon(T, (2, 1), (8, 7))
    # segment_in_polygon agrees with Java's segInPolygon away from grazing contacts
    rng = Rng(11)
    for n in (1, 12, 13):
        P = load_polygon(n)
        X, Y = _random_points(P, 120, rng)
        for k in range(0, 120, 2):
            a, b = (X[k], Y[k]), (X[k + 1], Y[k + 1])
            assert segment_in_polygon(P, a, b) == seg_in_polygon(P, a + b)
    assert segment_touches_boundary(L, (5, 5), (10, 15)) and not segment_touches_boundary(L, (5, 5), (8, 8))


def test_validate_path():
    P13 = load_polygon(13)
    assert validate_path(P13, path_record("P13")["waypoints"])[0] == (195.0, 620.0)
    with pytest.raises(ValueError, match="waypoint 9"):  # (69, 418) is outside polygon 13
        validate_path(P13, path_record("P13+")["waypoints"])
    with pytest.raises(ValueError, match="at least 2"):
        validate_path(P13, [(195.0, 620.0)])
    with pytest.raises(ValueError, match="zero length"):
        validate_path(P13, [(195.0, 620.0), (195.0, 620.0)])
    P12 = load_polygon(12)
    with pytest.raises(ValueError, match="leaves polygon"):
        validate_path(P12, [(90.0, 450.0), (580.0, 840.0)])
    with pytest.raises(ValueError, match="strictly inside"):
        validate_path(P12, [P12.vertices[0], (90.0, 450.0)])
    for label, rec in load_paths().items():  # every stored path that the Java runs inside the polygon
        if label not in ("P13+", "P13+tail"):
            validate_path(load_polygon(rec["polygon"]), rec["waypoints"])
    # a segment along an edge is refused (its inflection cuts are collinear with it: no ordered crossings)
    with pytest.raises(ValueError, match="runs along edge 4"):
        validate_path(COMB, [(5, 25), (5, 10), (45, 10), (45, 25)])
    validate_path(COMB, [(5, 25), (5, 9.99), (45, 9.99), (45, 25)])
    validate_path(Polygon(L_ROOM), [(15, 5), (5, 15)])  # through a reflex vertex: allowed


COMB = Polygon([(0, 0), (50, 0), (50, 30), (40, 30), (40, 10), (30, 10), (30, 30), (20, 30), (20, 10), (10, 10),
                (10, 30), (0, 30)], name="comb")
L_ROOM = [(0, 0), (20, 0), (20, 10), (10, 10), (10, 20), (0, 20)]


def _truth_mismatches(sim, ticks):
    """Target-ticks whose label disagrees with an independent visibility test from the robot."""
    bad = 0
    for _ in range(ticks):
        sim.step()
        bad += sum(segment_in_polygon(sim.poly, sim.robot, p) != (lab == 0)
                   for p, lab in zip(sim.targets, sim.target_labels))
    return bad


def test_ground_truth_near_walls_and_on_vertex_lines():
    """Regressions: (a) a path 0.0069 from the comb's wall (y = 9.99) has the right ground truth; (b) on a path
    through the reflex vertex of an L room, along the line through two vertices, the pockets are not phantoms
    (Java's scan reports the visible part as a shadow there); (c) a go_to that ends within PERTURB * sqrt(2)
    of a wall: Java's end sample is outside the polygon, where its scan sees no gaps, so every target would be
    labelled visible -- the fixed sample is pulled back inside."""
    sim = PolygonSimulator(COMB, [(5, 25), (5, 9.99), (45, 9.99), (45, 25)], 40, 5, speed=0.5, target_step=0.0)
    assert _truth_mismatches(sim, sim.ticks_per_lap()) == 0
    assert sim.stats["membership_fallbacks"] == 0
    L = Polygon(L_ROOM)
    h = get_gaps(L, Path([(15, 5), (5, 15)]))
    ge = gap_history_to_events(h)
    assert ge.initial == [1] and ge.events == [Appear(2, 0, 0), Disappear(1, 0, 0)]
    sim = PolygonSimulator(L, [(15, 5), (5, 15)], 40, 5, speed=0.5, target_step=0.0)
    assert _truth_mismatches(sim, sim.ticks_per_lap()) == 0
    assert sim.stats["membership_fallbacks"] == 0 and sim.stats["sense_fallbacks"] == 0
    sim = PolygonSimulator(12, [(300.0, 100.0), (400.0, 100.0)], 40, 3, speed=5.0)
    for _ in range(25):
        sim.step()
    sim.go_to((428.567, 100.0))
    for _ in range(15):
        sim.step()
    assert sim.poly.contains(sim.sense_point) and sim.stats["membership_fallbacks"] == 0
    hidden = sum(not segment_in_polygon(sim.poly, sim.robot, p) for p in sim.targets)
    assert hidden == sum(1 for x in sim.target_labels if x) > 0


def test_samples_past_a_waypoint_near_a_wall_stay_inside():
    """get_gaps on a path whose end lies 0.002 from a wall: Java samples 0.005 past the end, outside the
    polygon, and every shadow 'disappears' there; the fixed sample is pulled back inside."""
    L = Polygon([(0, 0), (20, 0), (20, 10), (10, 10), (10, 20), (0, 20)])
    for x in (19.998, 19.9):
        h = get_gaps(L, Path([(5, 15), (5, 2), (x, 2)]))
        assert h.final_ids == [2] and gap_history_to_events(h).events == [Disappear(1, 0, 0), Appear(2, 0, 0)]
        assert all(L.contains(p) for p in h.samples.points)
        assert (4 in h.samples.moved) == (x == 19.998) and h.samples.mode_independent == (x == 19.9)
    sj = sample_path(L, Path([(5, 15), (5, 2), (19.998, 2)]), "java")
    assert sj.points[-1] == (20.003, 2.0) and sj.physical[-1] == []  # Java: outside, sees nothing


# --------------------------------------------------------------------------- demo paths


def test_demo_paths_records():
    recs = demo_paths()
    assert sorted(recs, key=int) == [str(n) for n in MAPS]
    stored = load_paths()
    for k, r in recs.items():
        assert r["polygon"] == int(k) and r["provenance"] == "demo" and r["note"]
        assert r["closed"] == (r["waypoints"][0] == r["waypoints"][-1])
        if r["source"] is None:  # a patrol tour designed for the demo
            assert r["closed"] and r["at_end"] == "loop" and r["design"]["vias"]
        else:  # a path of the Java code, run back and forth
            assert stored[r["source"]]["polygon"] == int(k) and stored[r["source"]]["kind"] == "code"
            assert stored[r["source"]]["waypoints"] == r["waypoints"] and r["at_end"] == "reverse"
        path, at_end = load_demo_path(int(k))
        assert at_end == r["at_end"] and list(path.points) == [tuple(p) for p in r["waypoints"]]
        validate_path(load_polygon(int(k)), r["waypoints"])
    assert {recs[str(n)]["source"] for n in (1, 5, 12, 13, 14)} == {"P1", "P5", "P12", "P13", "P14b"}
    with pytest.raises(KeyError):
        load_demo_path(15)


@pytest.mark.parametrize("n", _maps())
def test_demo_path_event_mix(n):
    P = load_polygon(n)
    path, _ = load_demo_path(n)
    ev = gap_history_to_events(_track(P, path).history)  # get_gaps(P, path, java_matching=False), cached below
    c = Counter(type(e).__name__ for e in ev.events)
    assert all(c[t] >= 2 for t in ("Appear", "Disappear", "Split", "Merge")), c
    assert len(ev.events) >= 14


# --------------------------------------------------------------------------- pockets


def _boundary_dist(P, p):
    return min(pt_seg_dist(e, p) for e in P.edges)


@pytest.mark.parametrize("n", MAPS)
def test_pockets_tile_the_polygon(n):
    P = load_polygon(n)
    path, _ = load_demo_path(n)
    rng = Rng(100 + n)
    area = polygon_area(P.vertices)
    # points along the path and the waypoints themselves: some waypoints lie exactly on a line through two
    # vertices (e.g. P14b's (633.93, 712.29) on a general inflection), where Java's scan is degenerate; the
    # fixed physical gaps (shadow_pockets' default) are not
    pts = [(a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]))
           for a, b in zip(path.points[:10], path.points[1:11]) for f in (0.0, 0.31, 0.67)]
    X, Y = _random_points(P, 250, rng)
    for q in pts:
        sh = shadow_pockets(P, q)
        assert [s.label for s in sh] == list(range(1, len(sh) + 1))
        assert len(sh) == len(physical_gaps(P, q))
        V = visibility_polygon(P, q)
        assert abs(polygon_area(V) + sum(s.area for s in sh) - area) <= 1e-9 * area
        for s in sh:
            # counter-clockwise like the polygon; zero area only when q lies on the extension of an edge (P14b's
            # waypoint (616.07, 196.43) is below the vertical edge 63-64: Java's scan reports a zero-width gap)
            assert isinstance(s, Shadow) and s.area >= 0
            # chain vertices + free end on the boundary (the free end within Java's EPSILON: getIntersect is
            # slope-intercept, so near-vertical edges lose a few digits)
            assert all(_boundary_dist(P, p) < EPSILON for p in s.pocket)
            anchor, free = s.window
            assert anchor in P.vertices and free in s.pocket
            inner = (anchor[0] + (1 - 1e-9) * (free[0] - anchor[0]), anchor[1] + (1 - 1e-9) * (free[1] - anchor[1]))
            assert segment_in_polygon(P, anchor, inner)  # the window lies in the polygon
            assert all(point_in_polygon(P, ((s.pocket[i][0] + s.pocket[i + 1][0]) / 2,
                                            (s.pocket[i][1] + s.pocket[i + 1][1]) / 2)) or
                       _boundary_dist(P, ((s.pocket[i][0] + s.pocket[i + 1][0]) / 2,
                                          (s.pocket[i][1] + s.pocket[i + 1][1]) / 2)) < EPSILON
                       for i in range(len(s.pocket) - 1))
        vis = segments_in_polygon(P, q, X, Y)
        hits = sum(points_in_polygon(s.pocket, X, Y).astype(int) for s in sh) if sh else np.zeros(len(X), int)
        assert (hits <= 1).all()  # disjoint
        assert ((hits == 1) == ~vis).all()  # hidden iff in a pocket


def test_shadow_to_dict():
    P = load_polygon(12)
    sh = shadow_pockets(P, (90.0, 450.0), labels=[7, 8, 9, 10][:len(physical_gaps(P, (90.0, 450.0)))])
    d = sh[0].to_dict()
    assert d["label"] == 7 and len(d["window"]) == 2 and d["pocket"][0] == list(sh[0].pocket[0])


# --------------------------------------------------------------------------- the run equals get_gaps


def _component(frames):
    out = []
    for f in frames[1:]:
        for d in f["events"]:
            if d["type"] in COMPONENT:
                e = event_from_dict(d)
                if isinstance(e, Appear):
                    e = Appear(e.s, 0, 0)
                elif isinstance(e, Disappear):
                    e = Disappear(e.s, 0, 0)
                out.append(e)
    return out


@pytest.mark.parametrize("n", _maps())
def test_run_equals_get_gaps(n):
    P = load_polygon(n)
    path, _ = load_demo_path(n)
    h = get_gaps(P, path, "fixed", java_matching=False)  # the hidden-chain labels the simulator uses
    ge = gap_history_to_events(h)
    assert _track(P, path).history.lines() == h.lines()  # what the simulator consumes
    sim = PolygonSimulator(P, path, 6, n, speed=9.0, at_end="stop")
    assert sim.initial_labels == ge.initial == h.sample_ids[0]
    frames = [sim.frame()]
    while sim.leg.s < sim.leg.length or sim.leg.next_cp < len(sim.leg.track.distances):
        frames.append(sim.step())
        leg = sim.leg
        if sim.sense_point == sim.robot:  # labels are those get_gaps gives at the last critical point crossed
            assert sim.labels == h.sample_ids[leg.next_cp - 1]
        else:  # numerical fallback: the sensor is at that critical point's sample
            assert sim.sense_point == h.samples.points[leg.next_cp - 1]
    assert _component(frames) == ge.events
    assert sim.labels == h.final_ids
    assert sim.stats["crossings"] == len(h.samples.critical_points) - 1
    assert sim.stats["sense_fallbacks"] <= 3 and sim.stats["membership_fallbacks"] == 0
    # stopped at the end: the robot stays put and only the targets move
    f = sim.step()
    assert f["robot"] == list(path.points[-1]) and not [d for d in f["events"] if d["type"] in COMPONENT]


# --------------------------------------------------------------------------- tick semantics and ground truth


def _descendants(label, events):
    out = {label}
    for d in events:
        if d["type"] == "split" and d["s"] in out:
            out |= {d["a"], d["b"]}
        elif d["type"] == "merge" and (d["a"] in out or d["b"] in out):
            out.add(d["s"])
    return out


def _replay(counts, frame):
    """Exact replay of one tick's events on true counts (as for the grid, tests/test_grid.py)."""
    created = dict(map(tuple, frame["created"]))
    for d in frame["events"]:
        typ, s = d["type"], d["s"]
        if typ == "exit":
            counts[s] -= d["k"]
            assert counts[s] >= 0
        elif typ == "enter":
            counts[s] += d["k"]
        elif typ == "appear":
            assert d["lo"] == d["hi"] and s not in counts
            counts[s] = d["lo"]
        elif typ == "disappear":
            assert d["lo"] == d["hi"] == counts.pop(s)
        elif typ == "split":
            n = counts.pop(s)
            counts[d["a"]], counts[d["b"]] = created[d["a"]], created[d["b"]]
            assert counts[d["a"]] + counts[d["b"]] == n
        elif typ == "merge":
            counts[s] = counts.pop(d["a"]) + counts.pop(d["b"])
            assert counts[s] == created[s]
    return counts


def _check_run(sim, ticks):
    counts = dict(sim.initial_counts)
    assert sorted(counts) == sorted(sim.labels)
    covered = revealed = 0
    for _ in range(ticks):
        before = list(sim.target_labels)
        f = sim.step()
        counts = _replay(counts, f)
        assert sorted(counts.items()) == [tuple(p) for p in f["counts"]]
        assert sorted(counts) == sorted(f["labels"])
        assert f["target_labels"] == sim.target_labels
        for a, b in zip(before, sim.target_labels):  # a hidden target changing shadow is explained by events
            if a and b and b not in _descendants(a, f["events"]):
                types = {(d["type"], d["s"]) for d in f["events"]}
                assert ("enter", b) in types or ("appear", b) in types
        covered += sum(d["lo"] for d in f["events"] if d["type"] == "appear")
        revealed += sum(d["lo"] for d in f["events"] if d["type"] == "disappear")
    for mode in ("exact", "unknown", "evader"):
        seq = sim.sequence(mode)  # validates labels (alive, never reused)
        assert seq.alive_at_end() == sorted(sim.labels)
    return covered, revealed


@pytest.mark.parametrize("n", _maps())
def test_ground_truth_replay_many_seeds(n):
    seeds = (1, 2) if n in BIG | LONG else (1, 2, 3, 4, 5)
    tot = Counter()
    for seed in seeds:
        sim = PolygonSimulator(n, None, 12, seed, speed=12.0, target_step=15.0)
        cov, rev = _check_run(sim, sim.ticks_per_lap() + 3)
        tot["covered"] += cov
        tot["revealed"] += rev
        tot.update(d.__class__.__name__ for d in sim.history)
        assert sim.stats["membership_fallbacks"] == 0
    # shadows appear and disappear with zero area (the sample is PERTURB after the crossing), so Appear /
    # Disappear counts are almost always 0 here; test_sensor_step_attribution covers nonzero ones
    assert tot["Enter"] and tot["Exit"] and tot["Split"] and tot["Merge"], tot


@pytest.mark.parametrize("n", [3, 12])
def test_filter_bounds_contain_truth(n):
    sim = PolygonSimulator(n, None, 15, 4, speed=15.0, target_step=20.0)
    for _ in range(sim.ticks_per_lap()):
        sim.step()
    truth = sim.counts()
    for mode in ("exact", "unknown"):
        filt = CombinatorialFilter(sim.initial_condition(mode))
        filt.extend(sim.history)
        bounds = filt.all_bounds()
        assert set(bounds) == set(truth)
        for s, k in truth.items():
            assert bounds[s][0] <= k <= bounds[s][1]


def test_sensor_step_attribution():
    """One sensor step: per-target labels before / after plus the step's component events -> events."""
    sim = PolygonSimulator(12, None, 0, 1)
    evs = [Split(1, 10, 11), Merge(2, 3, 12), Appear(13, 0, 0), Disappear(4, 0, 0)]
    old = [0, 0, 4, 4, 1, 1, 1, 2, 3, 5, 5, 4]
    new = [13, 0, 0, 5, 10, 11, 0, 12, 12, 0, 12, 13]
    out, created = sim._sensor_events(old, new, evs)
    assert out == [Exit(1, 1), Exit(5, 2),  # 1 -> 0; 5 -> 0 and 5 -> 12 (12 is no descendant of 5)
                   Split(1, 10, 11), Merge(2, 3, 12), Appear(13, 2, 2), Disappear(4, 3, 3),
                   Enter(5, 1), Enter(12, 1)]  # 4 -> 5 (4 disappeared: counted by Disappear, then Enter)
    assert created == {10: 1, 11: 1, 12: 2}
    # chains through intermediate labels (transition_events for 3-way splits / merges)
    evs = [Merge(1, 2, 20), Split(20, 3, 21), Split(21, 4, 5)]
    out, created = sim._sensor_events([1, 1, 2, 2, 0], [3, 4, 5, 0, 4], evs)
    assert out == [Exit(2, 1)] + evs + [Enter(4, 1)]
    assert created == {20: 3, 3: 1, 21: 2, 4: 1, 5: 1}


def test_frames_and_determinism():
    a = PolygonSimulator(12, None, 5, 9, speed=6.0)
    b = PolygonSimulator(12, None, 5, 9, speed=6.0)
    fa, fb = a.run(60), b.run(60)
    assert fa == fb
    events = [event_from_dict(d) for f in fa[1:] for d in f["events"]]
    assert events == a.history
    assert fa[0]["t"] == 0 and fa[-1]["t"] == 60
    c = PolygonSimulator(12, None, 5, 10, speed=6.0)
    assert c.run(60) != fa
    g = a.frame(geometry=True)
    assert [s["label"] for s in g["shadows"]] == g["labels"] and len(g["visibility"]) >= 3
    # the targets: uniform in the polygon from shadowinfo.rng.Rng (two draws per placement attempt)
    rng = Rng(9)
    P = load_polygon(12)
    x0, x1, y0, y1 = (float(P.vx.min()), float(P.vx.max()), float(P.vy.min()), float(P.vy.max()))
    first = None
    while first is None:
        p = (x0 + rng.random() * (x1 - x0), y0 + rng.random() * (y1 - y0))
        if point_in_polygon(P, p) and not _on_boundary(P, p):
            first = p
    assert PolygonSimulator(12, None, 5, 9).targets[0] == first


def test_targets_stay_inside():
    sim = PolygonSimulator(7, None, 30, 3, speed=20.0, target_step=40.0)
    moved = 0
    for _ in range(80):
        before = list(sim.targets)
        sim.step()
        moved += sum(p != q for p, q in zip(before, sim.targets))
        assert all(point_in_polygon(sim.poly, p) for p in sim.targets)
        for p, q in zip(before, sim.targets):
            if p != q:
                assert segment_in_polygon(sim.poly, p, q) and not segment_touches_boundary(sim.poly, p, q)
    assert moved > 30 * 80 // 2


# --------------------------------------------------------------------------- laps, back and forth, click-to-move


@pytest.mark.parametrize("n, at_end", [(3, "loop"), (12, "reverse"), pytest.param(13, "reverse", marks=pytest.mark.slow)])
def test_several_laps(n, at_end):
    sim = PolygonSimulator(n, None, 10, 5, speed=11.0)
    assert sim.at_end == at_end
    lap = sim.ticks_per_lap()
    _check_run(sim, int(2.6 * lap))
    assert sim.stats["legs"] == 3 and len(sim.legs) == 3
    if at_end == "loop":  # same path again: the same events, relabelled by a bijection
        ev = gap_history_to_events(sim.legs[0].track.history).events
        assert len(sim.legs[1].track.events_at) == len(sim.legs[0].track.events_at)
        assert sim.legs[1].track is sim.legs[0].track and len(ev) > 0
    else:
        assert sim.legs[1].track.path.points == tuple(reversed(sim.legs[0].track.path.points))


def test_go_to_and_follow():
    sim = PolygonSimulator(12, None, 10, 2, speed=8.0)
    for _ in range(30):
        sim.step()
    with pytest.raises(ValueError, match="polygon 12"):
        sim.go_to((900.0, 100.0))
    target = (330.0, 470.0)
    validate_path(sim.poly, [sim.robot, target])
    sim.go_to(target)
    assert sim.at_end == "stop" and sim.path.points == (sim.robot, target)
    counts = dict(sim.counts())
    for _ in range(sim.ticks_per_lap() + 2):
        counts = _replay(counts, sim.step())
    assert sim.robot == target and sorted(counts.items()) == sorted(sim.counts().items())
    sim.follow([(560.0, 462.0), (575.0, 800.0)])
    for _ in range(sim.ticks_per_lap() + 2):
        counts = _replay(counts, sim.step())
    assert sim.robot == (575.0, 800.0)
    assert sorted(counts.items()) == sorted(sim.counts().items())
    assert ShadowSequence(sim.initial_condition(), sim.history).alive_at_end() == sorted(sim.labels)
    assert sim.stats["legs"] == 3


def test_constructor_options():
    P1 = load_polygon(1)
    sim = PolygonSimulator(P1, "P1s", 3, 1)
    assert sim.path.points[0] == (886.0, 281.0) and sim.at_end == "stop"
    with pytest.raises(ValueError, match="belongs to polygon"):
        PolygonSimulator(12, "P1s")
    with pytest.raises(ValueError, match="closed path"):
        PolygonSimulator(12, None, at_end="loop")
    with pytest.raises(ValueError, match="at_end"):
        PolygonSimulator(12, None, at_end="bounce")
    with pytest.raises(ValueError, match="speed"):
        PolygonSimulator(12, None, speed=0)
    with pytest.raises(ValueError):
        PolygonSimulator(load_path("P13+")[0], load_path("P13+")[1])  # leaves the polygon
    s2 = PolygonSimulator(13, load_path("P13+")[1], 3, 1, validate=False)  # allowed if asked: get_gaps copes
    assert s2.labels
    sim, frames = simulate(12, None, 4, 3, ticks=5)
    assert len(frames) == 6 and isinstance(sim, PolygonSimulator)
    sim, frames = simulate(12, Path([(90.0, 450.0), (580.0, 450.0)]), 4, 3, speed=50.0)
    assert len(frames) == 1 + math.ceil(490.0 / 50.0) and sim.robot == (580.0, 450.0)


def test_split_merge_counts_are_true_counts():
    """``created`` holds the number of targets each new shadow received (some nonzero over a run)."""
    sim = PolygonSimulator(14, None, 40, 1, speed=12.0, target_step=12.0)
    nonzero = 0
    for _ in range(120):
        f = sim.step()
        made = {d["a"] for d in f["events"] if d["type"] == "split"} | {d["b"] for d in f["events"] if d["type"] == "split"}
        made |= {d["s"] for d in f["events"] if d["type"] == "merge"}
        assert set(dict(map(tuple, f["created"]))) == made
        nonzero += sum(1 for _, k in f["created"] if k)
    assert nonzero > 0
