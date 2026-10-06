"""Grid world: RNG, visibility, shadow tracking, simulator self-checks, JS parity."""

from __future__ import annotations

import importlib.util
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import pytest

from shadowinfo.events import ShadowSequence, event_from_dict
from shadowinfo.grid import (
    GridMap, GridSimulator, ShadowTracker, expand_path, initial_bounds, label_components, load_office,
)
from shadowinfo.rng import Rng

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "parity_office.json"


@pytest.fixture(scope="module")
def office():
    return load_office()


# -- RNG ---------------------------------------------------------------------

def test_rng_reference_values():
    r = Rng(1)
    assert [r.next_u32() for _ in range(3)] == [2693262067, 11749833, 2265367787]
    r = Rng(42)
    xs = [r.random() for _ in range(1000)]
    assert all(0 <= x < 1 for x in xs)
    assert abs(sum(xs) / len(xs) - 0.5) < 0.05
    r = Rng(3)
    assert Counter(r.randint(5) for _ in range(5000)).keys() == set(range(5))
    with pytest.raises(ValueError):
        r.randint(0)


# -- visibility --------------------------------------------------------------

def _exact_los(m: GridMap, a, b) -> bool:
    """Reference: segment between centres meets no open obstacle square, and
    does not slip through a corner flanked by two obstacles."""
    (r0, c0), (r1, c1) = a, b
    y0, x0 = Fraction(2 * r0 + 1, 2), Fraction(2 * c0 + 1, 2)
    dy, dx = r1 - r0, c1 - c0

    def interval(p0, d, lo):
        if d == 0:
            return (Fraction(-1), Fraction(2)) if lo < p0 < lo + 1 else None
        t1, t2 = (lo - p0) / d, (lo + 1 - p0) / d
        return (min(t1, t2), max(t1, t2))

    for r in range(min(r0, r1), max(r0, r1) + 1):
        for c in range(min(c0, c1), max(c0, c1) + 1):
            if not m.blocked[m.index(r, c)]:
                continue
            iy, ix = interval(y0, dy, r), interval(x0, dx, c)
            if iy is None or ix is None:
                continue
            lo, hi = max(iy[0], ix[0], Fraction(0)), min(iy[1], ix[1], Fraction(1))
            if lo < hi:
                return False
    if dx and dy:
        sy, sx = (1 if dy > 0 else -1), (1 if dx > 0 else -1)
        for xc in range(min(c0, c1) + 1, max(c0, c1) + 1):
            t = (xc - x0) / dx
            yc = y0 + t * dy
            if 0 < t < 1 and yc.denominator == 1:
                yc = int(yc)
                rb, cb = (yc - 1 if sy > 0 else yc), (xc - 1 if sx > 0 else xc)
                if m.blocked[m.index(rb, cb + sx)] and m.blocked[m.index(rb + sy, cb)]:
                    return False
    return True


def test_corner_rule():
    m = GridMap([".#", "#."])
    assert not m.line_of_sight(0, 0, 1, 1)
    assert not m.line_of_sight(1, 1, 0, 0)
    m = GridMap([".#", ".."])
    assert m.line_of_sight(0, 0, 1, 1)
    m = GridMap(["...", ".#.", "..."])
    assert not m.line_of_sight(0, 0, 2, 2)
    assert m.line_of_sight(0, 0, 0, 2) and m.line_of_sight(0, 0, 2, 0)


def test_line_of_sight_matches_exact_geometry(office):
    m, _ = office
    rng = Rng(7)
    free = m.free_cells
    for _ in range(3000):
        a, b = m.cell(free[rng.randint(len(free))]), m.cell(free[rng.randint(len(free))])
        los = m.line_of_sight(*a, *b)
        assert los == m.line_of_sight(*b, *a)
        assert los == _exact_los(m, a, b), (a, b)


def test_radius_and_cache(office):
    m, path = office
    v = m.visibility(path[0], 5)
    r0, c0 = path[0]
    assert all((r - r0) ** 2 + (c - c0) ** 2 <= 30 for r, c in map(m.cell, (i for i, x in enumerate(v) if x)))
    assert m.visibility(path[0], 5) is v
    full = m.visibility(path[0])
    assert all(full[i] >= v[i] for i in range(len(v)))
    assert v[m.index(*path[0])] == 1


def test_radius_is_round_disc():
    # Range is "centre within R + 1/2": dr^2 + dc^2 <= R^2 + R, with no 1-cell
    # nubs at the axis extremes (the row at offset R is at least 3 cells wide).
    n = 35
    m = GridMap(["#" * n] + ["#" + "." * (n - 2) + "#"] * (n - 2) + ["#" * n])
    r0 = c0 = n // 2
    for radius in range(1, 16):
        v = m.visibility((r0, c0), radius)
        got = {m.cell(i) for i, x in enumerate(v) if x}
        want = {(r, c) for r in range(1, n - 1) for c in range(1, n - 1)
                if (r - r0) ** 2 + (c - c0) ** 2 <= radius * radius + radius}
        assert got == want, radius
        assert {(r0 + radius, c0 - 1), (r0 + radius, c0 + 1)} <= got, radius


def test_radius_no_wall_nub_shadow(office):
    # Regression: with the old <= R^2 disc, robot (13, 2) at R = 12 saw the nub
    # (25, 2) but not the wall-side cell (25, 1), which split off as a spurious
    # 1-cell shadow.
    m, path = office
    assert path[13] == (13, 2)
    v = m.visibility((13, 2), 12)
    assert v[m.index(25, 1)] == 1 and v[m.index(25, 2)] == 1


# -- map and path ------------------------------------------------------------

def test_office_map_and_path(office):
    m, path = office
    assert 48 <= m.cols <= 64 and 32 <= m.rows <= 40
    assert all(m.is_free(r, c) for r, c in path)
    for (r0, c0), (r1, c1) in zip(path, path[1:] + path[:1]):
        assert abs(r0 - r1) + abs(c0 - c1) == 1
    _, cells = label_components([1 - b for b in m.blocked], m.rows, m.cols)
    assert len(cells) == 1, "free space must be connected"


def test_maps_ship_inside_package():
    # load_office() must work from a non-editable install: the maps live in
    # the package and are declared as package data.
    import shadowinfo.grid as grid_pkg
    assert grid_pkg.MAPS_DIR.parent == Path(grid_pkg.__file__).resolve().parent
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert '"shadowinfo.grid" = ["maps/*.txt", "maps/*.json"]' in pyproject
    assert {p.suffix for p in grid_pkg.MAPS_DIR.iterdir()} <= {".txt", ".json"}


def test_map_trailing_whitespace():
    # Trailing spaces, tabs and CR are ignored (as in docs/js/grid.js); short
    # rows are padded with obstacles.
    for m in (GridMap.from_text("####\n#..#  \n####\t\n"), GridMap(["####\r\n", "#..#\r\n", "####\r\n"])):
        assert (m.rows, m.cols, m.free_cells) == (3, 4, [5, 6])
    assert GridMap.from_text("####\n#.\n####\n").free_cells == [5]
    with pytest.raises(ValueError):
        GridMap.from_text("####\n#. #\n####\n")


def test_expand_path():
    p = expand_path([[0, 0], [0, 2], [2, 2], [2, 0]])
    assert p == [(0, 0), (0, 1), (0, 2), (1, 2), (2, 2), (2, 1), (2, 0), (1, 0)]
    assert expand_path([[0, 0], [0, 2]], loop=False) == [(0, 0), (0, 1), (0, 2)]
    with pytest.raises(ValueError):
        expand_path([[0, 0], [1, 1]])


# -- shadow tracking ---------------------------------------------------------

def _mask(rows):
    return [1 if ch == "1" else 0 for row in rows for ch in row]


def _types(events):
    return [type(e).__name__ for e in events]


def test_label_components_scan_order():
    comp, cells = label_components(_mask(["1.1", "1.1", "..1"]), 3, 3)
    assert cells == [[0, 3], [2, 5, 8]]
    assert comp[8] == 1 and comp[1] == -1


def test_tracker_basic_events():
    t = ShadowTracker(1, 7)
    assert t.reset(_mask(["1111111"])) == [1]
    tr = t.update(_mask(["111.111"]))
    assert [repr(e) for e in tr.events] == ["Split(s=1, a=2, b=3)"]
    assert t.labels == [2, 2, 2, 0, 3, 3, 3]
    tr = t.update(_mask(["11..111"]))
    assert tr.events == [] and t.alive() == [2, 3]
    tr = t.update(_mask(["1111111"]))
    assert [repr(e) for e in tr.events] == ["Merge(a=2, b=3, s=4)"]
    tr = t.update(_mask(["......."]))
    assert _types(tr.events) == ["Disappear"] and tr.disappeared == [4]
    tr = t.update(_mask(["1.....1"]))
    assert [(type(e).__name__, e.s) for e in tr.events] == [("Appear", 5), ("Appear", 6)]


def test_tracker_general_group():
    t = ShadowTracker(5, 5)
    t.reset(_mask(["1...1"] * 5))
    tr = t.update(_mask(["11111", ".....", ".....", ".....", "11111"]))
    assert [repr(e) for e in tr.events] == ["Merge(a=1, b=2, s=3)", "Split(s=3, a=4, b=5)"]
    assert tr.merged_from == {3: [1, 2]} and tr.covers == {4: [0], 5: [1]}
    t = ShadowTracker(1, 9)
    t.reset(_mask(["111111111"]))
    tr = t.update(_mask(["1.1.1.1.1"]))
    assert [repr(e) for e in tr.events] == [
        "Split(s=1, a=2, b=3)", "Split(s=3, a=4, b=5)", "Split(s=5, a=6, b=7)", "Split(s=7, a=8, b=9)"]
    assert t.alive() == [2, 4, 6, 8, 9]


# -- simulator self-checks ---------------------------------------------------

def _descendants(label, events):
    out = {label}
    for d in events:
        if d["type"] == "split" and d["s"] in out:
            out |= {d["a"], d["b"]}
        elif d["type"] == "merge" and (d["a"] in out or d["b"] in out):
            out.add(d["s"])
    return out


def _replay(counts, frame):
    """Exact replay of one tick's events on true counts."""
    created = dict(map(tuple, frame["created"]))
    for d in frame["events"]:
        typ, s = d["type"], d["s"]
        if typ == "exit":
            counts[s] -= d["k"]
            assert counts[s] >= 0
        elif typ == "enter":
            counts[s] += d["k"]
        elif typ == "appear":
            assert d["lo"] == d["hi"]
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


@pytest.mark.parametrize("seed,radius", [(1, None), (2, None), (3, None), (4, 12), (5, 6)])
def test_simulator_self_consistency(office, seed, radius):
    m, path = office
    sim = GridSimulator(m, path, 30, seed, radius)
    counts = dict(sim.initial_counts)
    for _ in range(len(path)):
        before_labels, before_targets = list(sim.labels), list(sim.targets)
        f = sim.step()
        counts = _replay(counts, f)
        assert sorted(counts.items()) == [tuple(p) for p in f["counts"]]
        assert sorted(counts) == sim.tracker.alive()
        for i0, i1 in zip(before_targets, sim.targets):
            a, b = before_labels[i0], sim.labels[i1]
            if a and b and a != b and b not in _descendants(a, f["events"]):
                types = {(d["type"], d["s"]) for d in f["events"]}
                assert ("enter", b) in types and any(
                    ("exit", s) in types for s in _descendants(a, f["events"])), (f["t"], a, b)
    for mode in ("exact", "unknown", "evader"):
        seq = sim.sequence(mode)
        assert seq.alive_at_end() == sim.tracker.alive()
        assert ShadowSequence.from_json(seq.to_json()).events == seq.events
    assert sim.robot == path[0]


def test_frames_events_roundtrip(office):
    m, path = office
    sim = GridSimulator(m, path, 5, 9)
    frames = sim.run(40)
    events = [event_from_dict(d) for f in frames[1:] for d in f["events"]]
    assert events == sim.history
    assert frames[0]["t"] == 0 and frames[-1]["t"] == 40


def test_initial_bounds():
    assert initial_bounds({1: 2, 4: 0}) == {1: (2, 2), 4: (0, 0)}
    assert initial_bounds({1: 2}, "unknown") == {1: (0, math.inf)}
    assert initial_bounds({1: 2}, "evader") == {1: (0, 1)}
    with pytest.raises(ValueError):
        initial_bounds({1: 2}, "bogus")


def test_office_lap_event_mix(office):
    m, path = office
    sim = GridSimulator(m, path, 30, 1)
    frames = sim.run(len(path))
    c = Counter(d["type"] for f in frames[1:] for d in f["events"])
    for typ in ("appear", "disappear", "split", "merge", "enter", "exit"):
        assert c[typ] >= 3, c
    comp = c["appear"] + c["disappear"] + c["split"] + c["merge"]
    assert 50 <= comp <= 300, c
    assert any(f["created"] for f in frames[1:])


# -- Python/JS parity fixture --------------------------------------------------

def _load_tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_parity_fixture_up_to_date():
    fixture = json.loads(FIXTURE.read_text())
    fresh = json.loads(json.dumps(_load_tool("gen_parity_fixture").build_fixture()))
    assert fresh.keys() == fixture.keys()
    for key in fresh:
        if key != "runs":
            assert fresh[key] == fixture[key], key
    assert len(fresh["runs"]) == len(fixture["runs"])

    def header(run):
        return {k: v for k, v in run.items() if k != "frames"}

    for run_new, run_old in zip(fresh["runs"], fixture["runs"]):
        assert header(run_new) == header(run_old)
        for f_new, f_old in zip(run_new["frames"], run_old["frames"]):
            assert f_new == f_old, (run_new["seed"], f_new["t"])
        assert len(run_new["frames"]) == len(run_old["frames"])


def test_js_map_up_to_date():
    expected = _load_tool("gen_js_map").render("office")
    assert (ROOT / "docs" / "js" / "map_office.js").read_text() == expected
