"""Parity of shadowinfo.polygon (io, geometry, cuts, visibility) with the original Java implementation.

Golden data: tests/fixtures/java (tools/java_reference, docs/notes/original_java.md §5.2 items 1-4):

1. cuts of all 14 polygons (single tangents, inflections, general-inflection cuts, bitangent cuts,
   getCuts layout) -- compared exactly, double for double;
2. visibility: physical gaps and visibility polygons at lattice points, waypoints and the perturbed
   samples of every run -- exact;
3. pointInPolygon (float32 Path2D, non-zero rule) and segInPolygon(q -> v_j) -- exact;
4. critical points of every path (distance, point, cut type, cut index, path segment) and the
   perturbed sample points -- exact.

The scalar reference implementations are also checked against the vectorised ones on random input.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from pathlib import Path as FsPath

import numpy as np
import pytest

from shadowinfo.polygon import (Bitangent, CutType, GeneralInflection, GeometryError, InflectionType, Path,
                                PhysicalGap, Polygon, all_critical_points, all_cuts, bitangent_cuts,
                                bitangent_lines, general_inflection_cuts, inflections, java_double_str, load_path,
                                load_paths, load_polygon, path_labels, path_record, physical_gaps,
                                purturb_point_along_seg, seg_in_polygon, single_tangent_cuts, vertex_visibility,
                                visibility_polygon)
from shadowinfo.polygon import geometry as G
from shadowinfo.polygon import io as pio
from shadowinfo.polygon.cuts import Cut
from shadowinfo.polygon.simulate import pocket_polygon
from shadowinfo.polygon.visibility import _scan, _scan_reference

ROOT = FsPath(__file__).resolve().parents[1]
FIX = ROOT / "tests" / "fixtures" / "java"
POLYS = range(1, 15)
RUNS = sorted(p.stem for p in FIX.glob("*_*.json"))
BIG = {5, 13, 14}  # polygons with n >= 140: per-sample visibility takes ~0.5-2 s per run


def load(name):
    return json.loads((FIX / f"{name}.json").read_text())


def bits(xs) -> str:
    return "".join("1" if x else "0" for x in xs)


def gaps_list(gs):
    return [g.to_list() for g in gs]


def pts(ps):
    return [list(p) for p in ps]


_RUN_CACHE: dict = {}


def run_fixture(name):
    if name not in _RUN_CACHE:
        _RUN_CACHE[name] = load(name)
    return _RUN_CACHE[name]


# --------------------------------------------------------------------------- io and data


def test_dat_files_and_provenance():
    meta = json.loads((pio.DATA_DIR / "paths.json").read_text())
    assert meta["provenance"]["repository"] == "https://github.com/arc-l/shadow-information-space"
    assert "BSD-3-Clause" in meta["provenance"]["license"]
    for n in POLYS:
        raw = (pio.DATA_DIR / f"{n}.dat").read_bytes()
        rec = meta["polygons"][str(n)]
        assert hashlib.sha256(raw).hexdigest() == rec["sha256"]
        assert len(raw.decode().splitlines()) == rec["n_vertices"]
    assert {k for k, v in meta["polygons"].items() if "figure" in v} == {"12", "13", "14"}


@pytest.mark.parametrize("n", POLYS)
def test_load_polygon_matches_java(n):
    d = load(f"poly{n}")
    P = load_polygon(n)
    assert P is load_polygon(n)  # cached
    assert P.n == d["n"] == len(P.edges)
    assert pts(P.vertices) == d["vertices"]  # no transform: y-up .dat frame, as loaded by Java
    assert P.signed_area() > 0  # counter-clockwise
    assert list(P.turn) == d["turn"]
    assert P.reflex == d["reflex"]
    assert all(P.is_reflex(i) == (i in d["reflex"]) for i in range(P.n))
    for i, e in enumerate(P.edges):
        assert e == P.vertices[i] + P.vertices[(i + 1) % P.n]


def test_parse_dat_fixes_blank_lines_and_spacing(tmp_path):
    assert pio.parse_dat("1 2\n\n3  4\r\n 5 6 \n   \n7\t8\n") == [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0), (7.0, 8.0)]
    with pytest.raises(ValueError):
        pio.parse_dat("1 2\n3\n")
    f = tmp_path / "sq.dat"
    f.write_text("0 0\n10 0\n10 10\n0 10\n\n")
    P = load_polygon(f)
    assert P.n == 4 and P.name == "sq"
    with pytest.raises(ValueError):
        load_polygon(15)


def test_paths_json_has_every_path():
    paths = load_paths()
    code = {"P12", "P1", "P5", "P14a", "P14b", "P1s", "P13", "P13+", "P13+tail", "P5s", "P1L"}
    figs = {"fig_ICRA08-Fig2a", "fig_ICRA08-Fig3a", "fig_TRO-Fig11a", "fig_ICRA08-Fig8a", "fig_TRO-Fig15a",
            "fig_ICRA08-Fig8b", "fig_TRO-Fig15b"}
    assert set(paths) == code | figs and path_labels() == list(paths)
    # every row of the harness path table, with the same waypoints
    rows = [ln.split("\t") for ln in (ROOT / "tools/java_reference/paths.tsv").read_text().splitlines()
            if ln.strip() and not ln.startswith("#")]
    assert len(rows) == len(paths)
    for name, poly, ntg, _, wps in rows:
        rec = path_record(f"{name}@{poly}")
        assert rec["polygon"] == int(poly) and rec["harness"]["n_targets"] == int(ntg)
        assert rec["waypoints"] == [[float(c) for c in xy.split(",")] for xy in wps.split()]
        assert (ROOT / rec["harness"]["fixture"]).exists()
        fx = run_fixture(FsPath(rec["harness"]["fixture"]).stem)
        assert fx["waypoints"] == rec["waypoints"] and rec["harness"]["java_get_gaps_ok"] == fx["get_gaps"]["ok"]
    for k in code:
        assert paths[k]["kind"] == "code" and paths[k]["panels"] and "provenance" in paths[k]
    for k in figs:
        r = paths[k]
        assert r["kind"] == "digitised" and r["figure"] and r["paper_reported"]
        assert r["provenance"]["closest_code_path"] in code
        assert r["provenance"]["inside_polygon_check"]["ok"]
    assert paths["P13+"]["waypoints"][:9] == paths["P13"]["waypoints"]
    assert paths["P14a"]["waypoints"][:19] == paths["P14b"]["waypoints"][:19]
    assert path_record("ProjectPanel5") == paths["P14b"] and pio.PANELS["ProjectPanel"] == "P12"
    P, path = load_path("P12")
    assert P is load_polygon(12) and path.points[0] == (90.0, 450.0)
    with pytest.raises(KeyError):
        path_record("nope")


@pytest.mark.parametrize("x, s", [
    (0.0, "0.0"), (-0.0, "-0.0"), (1.0, "1.0"), (100.0, "100.0"), (-2.5, "-2.5"), (0.001, "0.001"),
    (0.0009, "9.0E-4"), (1e-4, "1.0E-4"), (1e7, "1.0E7"), (9999999.0, "9999999.0"), (1234567.125, "1234567.125"),
    (12345678.9, "1.23456789E7"), (0.05623, "0.05623"), (1.5e-10, "1.5E-10"), (6.02e23, "6.02E23"),
    (math.nan, "NaN"), (math.inf, "Infinity"), (-math.inf, "-Infinity")
])
def test_java_double_str(x, s):
    assert java_double_str(x) == s


def test_java_double_str_matches_gap_set_lines():
    n = 0
    for name in RUNS:
        gg = run_fixture(name)["get_gaps"]
        for s in gg.get("sets", []):
            assert s["line"].startswith(java_double_str(s["t"]) + " ")
            n += 1
    assert n > 1000


# --------------------------------------------------------------------------- JDK primitives


def test_relative_ccw_semantics():
    ab = (0.0, 0.0, 10.0, 0.0)
    assert G.relative_ccw(ab, (5.0, -1.0)) == 1  # right of A->B in y-up
    assert G.relative_ccw(ab, (5.0, 1.0)) == -1
    assert G.relative_ccw(ab, (-1.0, 0.0)) == -1  # collinear, before A
    assert G.relative_ccw(ab, (11.0, 0.0)) == 1  # collinear, beyond B
    assert G.relative_ccw(ab, (5.0, 0.0)) == 0 and G.relative_ccw(ab, (10.0, 0.0)) == 0
    assert G.lines_intersect(ab, (10.0, 0.0, 10.0, 5.0))  # touching counts
    assert not G.lines_intersect(ab, (11.0, 0.0, 12.0, 5.0))


def test_numpy_primitives_equal_scalar():
    rng = random.Random(3)
    vals = [0.0, 1.0, 100.0, 285.71, 1000.0, 1e-5, 3e-5]
    for _ in range(3000):
        c = [rng.choice(vals) if rng.random() < 0.3 else rng.uniform(-50, 1050) for _ in range(8)]
        l1, l2 = tuple(c[:4]), tuple(c[4:])
        p = (c[4], c[5])
        assert int(G._np_rccw(*c[:6])) == G.relative_ccw_xy(*c[:6])
        assert float(G._np_pt_seg_dist(*l1, *p)) == G.pt_seg_dist(l1, p)
        x, y, ok = G._np_get_intersect(*l1, *l2)
        ref = G.get_intersect(l1, l2)
        assert bool(ok) == (ref is not None)
        if ref is not None:
            assert (float(x), float(y)) == ref
        assert bool(G._np_on_extension(*l1, *p)) == G.on_extension(l1, p)
        assert bool(G._np_on_reverse_extension(*l1, *p)) == G.on_reverse_extension(l1, p)


def test_get_intersect_quirks():
    assert G.get_intersect((0, 0, 1, 1), (0, 1, 1, 2)) is None  # parallel
    assert G.get_intersect((0, 0, 0, 1), (5, 0, 5, 1)) is None  # both vertical
    assert G.get_intersect((0, 0, 0, 1), (-1, 3, 1, 3)) == (0, 3.0)
    assert G.get_point_on_line_from_x((2, 0, 2, 5), 3.0) == (3.0, 0.0)  # meaningless, as in Java
    assert G.get_point_on_line_from_x((2, 0, 2, 5), 2.0) is None
    assert G.get_point_on_line_from_y((0, 2, 5, 2), 2.0) is None
    assert G.on_extension((0, 0, 10, 10), (11, 0))  # per-axis OR, not a projection
    assert G.get_segment_intersect((0, 0, 10, 0), (5, -1, 5, 1)) == (5, 0.0)
    assert G.get_segment_intersect((0, 0, 10, 0), (15, -1, 15, 1)) is None
    assert G.get_x_intercept((0, 1, 2, 5)) == G.get_y_intercept((0, 1, 2, 5)) == 1.0


def test_small_primitives_match_java():
    """``epsilonEqual``, ``getSlope``, ``pointOnSegment`` and the buggy, unused ``segmentsOnSameLine``; the expected
    values were printed by the original ``Algorithm`` class (JDK 17), not derived from the Python."""
    assert [G.epsilon_equal(*x) for x in [(1, 1.00004), (1, 1.00005), (1, 1.00006), (-3, -3.00005000001)]] == \
        [True, False, False, False]  # 1.00005 - 1 > 5e-5 in doubles
    assert [G.get_slope(ln) for ln in [(0, 0, 0.00004, 10), (0, 0, 0.00006, 10), (1, 2, 3, 7), (3, 7, 1, 2),
                                       (0, 5, 10, 5)]] == [math.inf, 166666.66666666666, 2.5, 2.5, 0.0]
    assert [G.point_on_segment(q, (0, 5, 10, 5)) for q in [(5, 5.00004), (5, 5.0001), (10.00004, 5), (10.0001, 5),
                                                          (-0.00003, 5)]] == [True, False, True, False, True]
    # Java's segmentsOnSameLine compares the second intercept at line1.P2: parallel lines count as "same line"
    assert all(G.segments_on_same_line(a, b) for a, b in [((0, 0, 1, 1), (2, 2, 3, 3)), ((0, 0, 1, 1), (0, 1, 1, 2)),
                                                          ((0, 0, 1, 2), (5, 10, 6, 12)), ((0, 0, 0, 1), (0, 2, 0, 3)),
                                                          ((1, 0, 2, 1), (0, 0, 3, 3))])
    sq = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    assert sq.edge(0) == (0, 0, 10, 0) and sq.edge(-1) == sq.edge(3) == (0, 10, 0, 0)  # getLineSegmentArray
    assert isinstance(sq.general_path, G.GeneralPath) and sq.general_path.contains(5, 5)  # getGeneralPath


def test_path_get_intersect_point_matches_java():
    """``Path.getIntersectPoint(Line2D)``: the first path segment the line crosses (values from the JDK 17 run)."""
    p = Path([(90, 450), (580, 450), (580, 840), (520, 840)])
    got = [p.intersect_point(ln) for ln in [(500, 0, 500, 1000), (0, 600, 1000, 600), (550, 800, 600, 900),
                                             (0, 0, 10, 10), (600, 400, 560, 860)]]
    assert got == [(500.0, 450.0), (580.0, 600.0), (570.0, 840.0), None, (580.0, 630.0)]


def test_dat_writer_and_random_polygons_match_java(tmp_path):
    """``Polygon.write`` (byte-identical to the Java output for map 12, trailing blank line included),
    ``randomPolygonFromData`` and ``randomPolygon`` after seeding ``Math.random()`` (values from the JDK 17 run)."""
    txt = pio.format_dat(load_polygon(12))
    assert hashlib.sha256(txt.encode()).hexdigest()[:16] == "f1f1703a2a4b499a"  # sha256 of the Java output
    assert txt.startswith("0.0 0.0\n428.57 0.0\n428.57 250.0\n571.43 250.0\n") and txt.endswith("\n\n")
    for n in POLYS:  # parse_dat skips the blank line Java's own read() throws on (B16): a round trip
        P = load_polygon(n)
        pio.write_dat(P, tmp_path / f"{n}.dat")
        assert pio.read_dat(tmp_path / f"{n}.dat") == list(P.vertices)
    assert (pio.NUMBER_OF_POLYGON_DATA_FILES, pio.DATA_FILE_POSTFIX) == (7, ".dat")
    # Math.random() seeded with s -> polygon file (int) (r * 7 + 1); small seeds all draw r ~ 0.73 (map 6)
    seeds = {1: "6", 12: "6", 1000: "5", 2024: "5", 31337: "1", 123456789: "5", 987654321: "3", -5: "2"}
    assert {s: pio.random_polygon_from_data(s).name for s in seeds} == seeds
    r = pio.random_polygon(5, 1000, 600, rng=7)  # Java: randomPolygon(5, new Dimension(1000, 600), dc)
    assert r.vertices == ((730.6990420600421, 449.50176188017986), (348.30970303125696, 538.3662856452628),
                          (708.1771577767972, 211.1488665877841), (120.73605139050846, 509.94786050836353),
                          (83.21971724462152, 557.2488436935114))


def test_purturb_point_along_seg():
    p = (5.0, 5.0)
    assert G.purturb_point_along_seg(p, (0, 0, 10, 0)) == (5.005, 5.0)
    assert G.purturb_point_along_seg(p, (0, 0, 10, 0), far=False) == (4.995, 5.0)
    assert G.purturb_point_along_seg(p, (5, 10, 5, 0)) == (5.0, 4.995)
    s = G.purturb_point_along_seg(p, (0, 0, 10, 4))  # |slope| < 1: step 0.005 in x
    assert s[0] == 5.005
    s = G.purturb_point_along_seg(p, (0, 0, 2, 10))  # |slope| > 1: step 0.005 in y
    assert s[1] == 5.005
    with pytest.raises(GeometryError):  # B1: exact 45 degrees throws in Java
        G.purturb_point_along_seg(p, (0, 0, 10, 10), compat="java")
    assert G.purturb_point_along_seg(p, (0, 0, 10, 10)) == (5.005, 5.005)
    with pytest.raises(ValueError):
        G.purturb_point_along_seg(p, (0, 0, 10, 0), compat="nope")


def test_general_path_contains_float32_rule():
    sq = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    assert sq.contains((5, 5)) and sq.contains((0, 5)) and not sq.contains((10, 5))  # half-open
    assert sq.contains((5, 0)) and not sq.contains((5, 10))
    assert not sq.contains((math.nan, 1)) and not sq.contains((math.inf, 1))
    tiny = Polygon([(0, 0), (1, 0), (1, 0.1 + 1e-9), (0, 0.1 + 1e-9)])  # 0.1+1e-9 rounds to float32 0.1
    y = float(np.float32(0.1 + 1e-9))
    assert not tiny.contains((0.5, y)) and tiny.contains((0.5, np.nextafter(y, 0)))
    rng = random.Random(5)
    for n in (1, 12, 14):
        P = load_polygon(n)
        xs = [rng.uniform(-10, 1010) for _ in range(300)] + [v[0] for v in P.vertices]
        ys = [rng.uniform(-10, 1010) for _ in range(300)] + [v[1] for v in P.vertices]
        assert P.contains_many(xs, ys).tolist() == [P.general_path.contains_scalar(x, y) for x, y in zip(xs, ys)]


@pytest.mark.parametrize("n", [1, 3, 12, 13])
def test_seg_in_polygon_vectorised_equals_scalar(n):
    P = load_polygon(n)
    rng = random.Random(n)
    V = P.vertices
    segs = []
    for _ in range(150):
        a = V[rng.randrange(P.n)] if rng.random() < 0.5 else (rng.uniform(0, 1000), rng.uniform(0, 1000))
        b = V[rng.randrange(P.n)] if rng.random() < 0.5 else (rng.uniform(0, 1000), rng.uniform(0, 1000))
        segs.append(a + b)
    a = np.asarray(segs)
    many = G.seg_in_polygon_many(P, a[:, 0], a[:, 1], a[:, 2], a[:, 3])
    assert many.tolist() == [G.seg_in_polygon_scalar(P, s) for s in segs] == [seg_in_polygon(P, s) for s in segs]
    assert any(many) and not all(many)


# --------------------------------------------------------------------------- cuts (fixture item 1)


@pytest.mark.parametrize("n", POLYS)
def test_cuts_match_java(n):
    d = load(f"poly{n}")
    P = load_polygon(n)
    c = d["cuts"]
    st = single_tangent_cuts(P)
    assert all(x.type is CutType.SINGLETANGENT for x in st)
    assert [list(x.line) for x in st] == c["single_tangent"]
    assert pts(inflections(P, InflectionType.NONGENERAL)) == c["inflection_nongeneral"]
    assert pts(inflections(P, "ALL")) == c["inflection_all"]
    gi = general_inflection_cuts(P)
    assert [[list(g.line), g.from_line, g.counter_clockwise] for g in gi] == c["general_inflection_cut"]
    assert c["inflection_general_equals_general_inflection_cut_lines"]
    assert inflections(P, InflectionType.GENERAL) == [g.line for g in gi]
    bt = bitangent_cuts(P)
    assert [[list(b.line), b.this_point, b.opposite_point, list(b.opposite_segment), list(b.curve_to_point)]
            for b in bt] == c["bitangent_cut"]
    assert c["bitangent_lines_equal_bitangent_cut_lines"] and bitangent_lines(P) == [b.line for b in bt]
    cuts = all_cuts(P)
    layout, cur = [], None
    for x in cuts:
        if x.type.value != cur:
            cur = x.type.value
            layout.append([cur, 0])
        layout[-1][1] += 1
    assert layout == c["get_cuts_layout"] and len(cuts) == c["get_cuts_total"]
    assert cuts == st + [Cut(ln, CutType.NONGENERAL_INFLECTION) for ln in inflections(P, "NONGENERAL")] + gi + bt
    assert all(isinstance(x, GeneralInflection) for x in gi) and all(isinstance(x, Bitangent) for x in bt)
    counts = d["counts"]
    assert (len(st), len(gi), len(bt)) == (counts["single_tangent"], counts["general_inflection_cut"],
                                          counts["bitangent_cut"])
    assert all(P.is_reflex(g.from_line + 1 if g.counter_clockwise else g.from_line) for g in gi)


def test_convex_polygon_has_no_cuts_or_gaps():
    sq = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
    assert all_cuts(sq) == [] and inflections(sq) == [] and sq.reflex == []
    assert physical_gaps(sq, (3.0, 4.0)) == []
    assert vertex_visibility(sq, (3.0, 4.0)).all()
    assert visibility_polygon(sq, (3.0, 4.0)) == list(sq.vertices)
    assert len(visibility_polygon(sq, (3.0, 4.0), compat="java")) == 8  # each vertex twice (B15)


L_ROOM = [(0, 0), (20, 0), (20, 10), (10, 10), (10, 20), (0, 20)]


def test_fixed_gaps_on_a_line_through_two_vertices():
    """(12, 8) lies on the line through the reflex vertex (10, 10) and the vertex (0, 20): Java's scan reports
    the one shadow (area 50) three times, once as its complement (the visible part, area 250)."""
    L = Polygon(L_ROOM)
    java = physical_gaps(L, (12, 8), compat="java")
    assert len(java) == 3 and sorted(round(G_area(L, g), 6) for g in java) == [50.0, 50.0, 250.0]
    assert physical_gaps(L, (12, 8)) == [PhysicalGap(3, 5, None, (0.0, 20.0))]
    assert physical_gaps(L, (12, 8.0001)) == physical_gaps(L, (12, 8.0001), compat="java")  # off the line
    for q in [(12, 8), (15, 5), (11, 9), (12, 8.0001), (12, 7.9999)]:
        assert_tiles(L, q)
    # map 14 at the P14a waypoint (633.93, 714.29): 19 Java pockets covering 1.99 times the polygon
    P = load_polygon(14)
    q = (633.93, 714.29)
    assert len(physical_gaps(P, q, compat="java")) == 19 and len(physical_gaps(P, q)) == 18
    assert_tiles(P, q)


def G_area(P, g):
    from shadowinfo.polygon.simulate import polygon_area

    return polygon_area(pocket_polygon(P, g))


def test_fixed_gaps_tile_at_waypoints_and_on_vertex_lines():
    """Every stored waypoint inside its polygon, and points on lines through a reflex vertex and another vertex
    of every map: the fixed pockets tile the polygon (Java's do not at some of them)."""
    rng = random.Random(5)
    for label in load_paths():
        P, path = load_path(label)
        for w in path.points:
            if pio_status(P, w) == "inside":
                assert_tiles(P, w, 800)
    for n in range(1, 15):
        P = load_polygon(n)
        V = P.vertices
        done = 0
        while done < 6:
            r, v = rng.choice(P.reflex), rng.randrange(P.n)
            t = rng.choice([0.25, 0.5, 0.75, 1.5, 2.0, -0.5])
            q = (V[r][0] + t * (V[r][0] - V[v][0]), V[r][1] + t * (V[r][1] - V[v][1]))
            if v == r or pio_status(P, q) != "inside" or min(G.pt_seg_dist(e, q) for e in P.edges) < 1e-3:
                continue
            done += 1
            assert_tiles(P, q, 800)


def pio_status(P, q):
    from shadowinfo.polygon.visibility import point_status

    return point_status(P, q)


def test_fixed_visibility_refuses_points_outside():
    L = Polygon(L_ROOM)
    for q in [(15.0, 15.0), (25.0, 5.0)]:
        assert physical_gaps(L, q, compat="java") == []
        with pytest.raises(GeometryError):
            physical_gaps(L, q)
        with pytest.raises(GeometryError):
            visibility_polygon(L, q)
    assert physical_gaps(L, (10.0, 15.0)) == physical_gaps(L, (10.0, 15.0), compat="java")  # on an edge: allowed


def test_polygon_must_be_simple_and_counter_clockwise():
    with pytest.raises(ValueError, match="counter-clockwise"):
        Polygon(L_ROOM[::-1])
    with pytest.raises(ValueError, match="not simple"):
        Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (5, -5)])  # edges 0 and 3 cross
    with pytest.raises(ValueError, match="not simple"):
        Polygon([(0, 0), (10, 0), (5, 0), (5, 5)])  # edge 1 folds back onto edge 0
    with pytest.raises(ValueError, match="not simple"):
        Polygon([(0, 0), (4, 0), (4, 4), (2, 0), (0, 4)])  # vertex 3 touches edge 0
    with pytest.raises(ValueError, match="not simple"):
        Polygon([(0, 0), (1, 0), (1, 0), (0, 1)])  # zero-length edge
    with pytest.raises(ValueError, match="finite"):
        Polygon([(0, 0), (1, 0), (float("nan"), 1)])
    Polygon([(0, 0), (5, 0), (10, 0), (10, 10), (0, 10)])  # a straight vertex is fine
    assert Polygon(L_ROOM[::-1], validate=False).signed_area() < 0  # Java's behaviour, on request
    for n in range(1, 15):
        assert load_polygon(n).signed_area() > 0


# --------------------------------------------------------------------------- pointInPolygon / segInPolygon / visibility


@pytest.mark.parametrize("n", POLYS)
def test_point_in_polygon_matches_java(n):
    pip = load(f"poly{n}")["point_in_polygon"]
    P = load_polygon(n)
    xs = [50.0 * xi for xi in range(1, 20)]
    rows = [bits(P.contains_many(xs, [50.0 * yi] * 19)) for yi in range(1, 20)]
    assert rows == pip["rows"]
    assert bits(P.contains(v) for v in P.vertices) == pip["at_vertices"]
    assert bits(P.contains(((e[0] + e[2]) / 2, (e[1] + e[3]) / 2)) for e in P.edges) == pip["at_edge_midpoints"]


# Lattice points of the golden visibility queries that lie on a line through a reflex vertex and another
# vertex (or along an edge): Java's scan misses, duplicates or inverts a gap there, the fixed scan does not.
DEGENERATE_LATTICE = {
    11: [(100.0, 500.0), (400.0, 500.0), (500.0, 500.0), (600.0, 500.0)],
    12: [(500.0, 300.0), (500.0, 400.0), (500.0, 500.0)],
    13: [(400.0, 400.0), (500.0, 500.0), (900.0, 900.0)],
    14: [(100.0, 600.0), (900.0, 600.0), (800.0, 700.0)],
}


def assert_tiles(P, q, n_samples=3000, seed=0):
    """The fixed shadow pockets of ``q`` and its visibility polygon tile ``P``: the areas add up, and every
    sample point of ``P`` is visible (closed visibility, exact predicates) or in exactly one pocket."""
    from shadowinfo.polygon.simulate import points_in_polygon, polygon_area, segments_in_polygon

    gaps = physical_gaps(P, q)
    pockets = [pocket_polygon(P, g) for g in gaps]
    area = sum(polygon_area(pk) for pk in pockets) + polygon_area(visibility_polygon(P, q))
    assert area == pytest.approx(P.signed_area(), rel=1e-9), q
    rng = np.random.default_rng(seed)
    px = rng.uniform(P.vx.min(), P.vx.max(), n_samples)
    py = rng.uniform(P.vy.min(), P.vy.max(), n_samples)
    inside = points_in_polygon(P, px, py)
    px, py = px[inside], py[inside]
    vis = segments_in_polygon(P, q, px, py)
    cnt = sum((points_in_polygon(pk, px, py).astype(int) for pk in pockets), np.zeros(len(px), dtype=int))
    assert np.all(cnt == (~vis).astype(int)), (q, int(np.sum(cnt != (~vis).astype(int))))


@pytest.mark.parametrize("n", POLYS)
def test_visibility_matches_java(n):
    d = load(f"poly{n}")
    P = load_polygon(n)
    n_vp = 0
    degenerate = []
    for v in d["visibility"]:
        q = tuple(v["q"])
        assert bits(vertex_visibility(P, q)) == v["seg_in_polygon_to_vertices"], q
        gaps = physical_gaps(P, q, compat="java")
        assert gaps_list(gaps) == v["physical_gaps"], q
        for g in gaps:
            assert (g.start_point is None) != (g.end_point is None)
            fs, fe = g.full_start_edge(P.n), g.full_end_edge(P.n)
            assert 0 <= fs < P.n and 0 <= fe < P.n
        fixed_gaps = physical_gaps(P, q)
        if fixed_gaps != gaps:
            degenerate.append(q)
            assert_tiles(P, q)
        if "visibility_polygon" in v:
            n_vp += 1
            java = visibility_polygon(P, q, compat="java")
            assert pts(java) == v["visibility_polygon"], q
            fixed = visibility_polygon(P, q)
            assert all(a != b for a, b in zip(fixed, fixed[1:] + fixed[:1]))
            if fixed_gaps == gaps:
                assert [p for k, p in enumerate(java) if k == 0 or p != java[k - 1]][:len(fixed)] == fixed
    assert n_vp > 0
    # lattice points on a line through a reflex vertex and another vertex, where Java's scan is degenerate
    assert degenerate == DEGENERATE_LATTICE.get(n, [])


@pytest.mark.parametrize("n", POLYS)
def test_seg_in_polygon_to_vertices_scalar(n):
    """segInPolygon(q -> v_j) through the generic and the scalar code paths (first lattice points)."""
    d = load(f"poly{n}")
    P = load_polygon(n)
    for v in d["visibility"][:3 if n in BIG else 8]:
        q = tuple(v["q"])
        segs = [q + p for p in P.vertices]
        assert bits(seg_in_polygon(P, s) for s in segs) == v["seg_in_polygon_to_vertices"]
        if P.n <= 70:
            assert bits(G.seg_in_polygon_scalar(P, s) for s in segs) == v["seg_in_polygon_to_vertices"]


@pytest.mark.parametrize("n", [1, 3, 6, 7, 8, 12, 13])
def test_scan_equals_literal_java_loop(n):
    """The vectorised scan equals a literal scalar transcription of getPhysicalGaps' loop."""
    P = load_polygon(n)
    rng = random.Random(100 + n)
    qs = []
    while len(qs) < (4 if n == 13 else 25):
        q = (rng.uniform(0, 1000), rng.uniform(0, 1000))
        if P.contains(q) or rng.random() < 0.1:
            qs.append(q)
    qs += [P.vertices[0], ((P.vertices[0][0] + P.vertices[1][0]) / 2, (P.vertices[0][1] + P.vertices[1][1]) / 2)]
    for q in qs:
        assert _scan(P, q) == _scan_reference(P, q), q


# --------------------------------------------------------------------------- paths and critical points (item 4)


@pytest.mark.parametrize("name", RUNS)
def test_critical_points_match_java(name):
    r = run_fixture(name)
    P = load_polygon(r["polygon"])
    path = Path(r["waypoints"])
    assert path.length() == r["path_length"]
    assert bits(P.contains(w) for w in path.points) == r["waypoints_inside"]
    cuts = all_cuts(P)
    cps = all_critical_points(path, cuts)
    assert cps == path.all_critical_points(cuts)
    assert len(cps) == r["n_critical_points"]
    assert sum(cp.cut_type in (CutType.GENERAL_INFLECTION, CutType.BITANGENT) for cp in cps) == r["n_gi_bt_crossings"]
    for cp, row in zip(cps, r["critical_points"]):
        assert [cp.distance, list(cp.point), cp.cut_type.value, cp.cut_index, cp.seg_index] == row[:5]
        assert cp.seg == path.segments[cp.seg_index]
        assert (cp.cut is None) == (cp.cut_index == -1) and (cp.cut is None or cp.cut is cuts[cp.cut_index])
        assert list(purturb_point_along_seg(cp.point, cp.seg, compat="java")) == row[5]
    assert cps[0].distance == 0.0 and cps[0].cut_type is CutType.NONE
    assert cps[-1].cut_type is CutType.NONE and cps[-1].point == path.points[-1]
    pi = r["panel_intersect_points"]
    assert pts(path.intersect_points(bitangent_lines(P))) == pi["bitangent_lines"]
    assert pts(path.intersect_points(inflections(P, InflectionType.GENERAL))) == pi["general_inflection_lines"]


@pytest.mark.parametrize("name", RUNS)
def test_run_visibility_matches_java(name):
    r = run_fixture(name)
    P = load_polygon(r["polygon"])
    for v in r["visibility"]:
        q = tuple(v["q"])
        assert gaps_list(physical_gaps(P, q, compat="java")) == v["physical_gaps"], v["where"]
        assert pts(visibility_polygon(P, q, compat="java")) == v["visibility_polygon"], v["where"]


def _sample_gaps(name):
    r = run_fixture(name)
    P = load_polygon(r["polygon"])
    rows = [row for row in r["critical_points"] if len(row) > 6]
    assert rows, name
    for k, row in enumerate(rows):
        assert gaps_list(physical_gaps(P, tuple(row[5]), compat="java")) == row[6], (name, k)


SAMPLE_RUNS = [n for n in RUNS if not n.split("_", 1)[1].startswith("fig_")]


@pytest.mark.parametrize("name", [pytest.param(n, marks=pytest.mark.slow) if int(n.split("_")[0]) in BIG else n
                                  for n in SAMPLE_RUNS])
def test_sample_physical_gaps_match_java(name):
    """getPhysicalGaps at every perturbed critical point (the samples getGaps uses)."""
    _sample_gaps(name)


def test_path_distance_quirks():
    path = Path([(0, 0), (10, 0), (10, 10)])
    assert path.length() == 20.0
    assert path.pt_dist_from_start((5, 0)) == 5.0
    assert path.pt_dist_from_start((10, 0)) == 10.0  # on both segments: the last one wins (B17)
    assert path.pt_dist_from_start((10.0, 1e-5)) == 10.00001
    assert math.isnan(path.pt_dist_from_start((5, 5)))
    assert path.pt_dist_on_segment((10.0, 1e-5), 0) != path.pt_dist_from_start((10.0, 1e-5))
    assert path.point_on_path((3, 0)) and not path.point_on_path((3, 1))
    with pytest.raises(ValueError):
        Path([(0, 0)])


@pytest.mark.parametrize("label, n_java, n_all, digest", [("P12", 28, 30, "e0010a1b519ab938"),
                                                          ("P13", 190, 192, "e103b17e872a741e")])
def test_cut_intersect_points_match_java(label, n_java, n_all, digest):
    """``Path.getIntersectPoints(getCuts(poly))``: the Java printed the count, ``getAllCriticalPoints`` length and
    one ``distance cutType`` line per crossing (JDK 17); the digest is the sha256 of that text."""
    from shadowinfo.polygon import cut_intersect_points

    poly, path = load_path(label)
    ip = cut_intersect_points(path, all_cuts(poly))
    txt = f"{len(ip)} {len(all_critical_points(path, all_cuts(poly)))}\n" + "".join(
        f"{java_double_str(c.distance)} {c.cut_type.name}\n" for c in ip)
    assert (len(ip), hashlib.sha256(txt.encode()).hexdigest()[:16]) == (n_java, digest)
    allc = all_critical_points(path, all_cuts(poly))
    assert len(allc) == n_all and [c for c in allc if c.cut_index >= 0] == ip  # start/end added, nothing lost


def test_critical_points_last_writer_wins():
    path = Path([(0, 0), (10, 0), (10, 10)])
    a = Cut((5, -1, 5, 1), CutType.SINGLETANGENT)
    b = Cut((5, 1, 5, -1), CutType.BITANGENT)  # same crossing distance, later in the list
    start = Cut((0, -1, 0, 1), CutType.GENERAL_INFLECTION)  # crosses at distance 0: overwritten by start
    corner = Cut((9, -1, 11, 1), CutType.NONGENERAL_INFLECTION)  # crosses at the waypoint (10, 0)
    cps = all_critical_points(path, [a, b, start, corner])
    assert [(cp.distance, cp.cut_type, cp.cut_index) for cp in cps] == [
        (0.0, CutType.NONE, -1), (5.0, CutType.BITANGENT, 1), (10.0, CutType.NONGENERAL_INFLECTION, 3),
        (20.0, CutType.NONE, -1)]
    assert cps[2].seg_index == 1  # found on both segments with key 10.0: the second put wins
    assert path.intersect_points([a.line, corner.line]) == [(5, 0.0), (10, 0.0)]


def test_physical_gap_full_edges():
    g = PhysicalGap(5, 9, (1.0, 2.0), None)
    assert (g.full_start_edge(10), g.full_end_edge(10)) == (6, 9)
    g = PhysicalGap(9, 0, (1.0, 2.0), None)
    assert g.full_start_edge(10) == 0
    g = PhysicalGap(3, 0, None, (1.0, 2.0))
    assert (g.full_start_edge(10), g.full_end_edge(10)) == (3, 9)
    assert g.to_list() == [3, 0, None, [1.0, 2.0]]
