"""Critical lines ("cuts") of a polygon and their crossings by a robot path.

Port of ``Algorithm.getInflection``, ``getGeneralInflectionCut``, ``getBitangentCut``, ``getBitangent``,
``getSingletangentCut``, ``getCuts`` and ``Path.getAllCriticalPoints`` (``docs/notes/original_java.md``
§2.5-§2.6).  Paper: T-RO 2012 Sec. II-A (a shadow component changes only when the sensor crosses an
inflection or a bitangent ray) and Sec. V-C (the visibility cell decomposition behind the sampling).

* **Inflection rays** start at a reflex vertex and extend a boundary edge to the nearest boundary
  hit.  *General* ones (the neighbouring vertex is convex) make a gap appear or disappear;
  *non-general* ones (neighbour reflex) make it jump to another vertex.
* **Bitangent rays** continue the segment between two mutually visible reflex vertices that are
  tangent on the same side; crossing one splits or merges gaps.
* **Single-tangent rays** (tangent at one vertex, crossing at the other) change only the vertex a gap
  hangs on.  Crossing them triggers a resample in ``getGaps`` but no component event.

``all_cuts`` concatenates single tangents, non-general inflections, general inflections and
bitangents in this order; the order matters because equal crossing distances keep the *last* cut
(quirk B2).  All results are bit-identical to Java (golden fixtures ``tests/fixtures/java/poly*.json``).
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import (EPSILON, GeometryError, Line, Path, Point, Polygon, _np_get_intersect, _np_on_extension,
                       _np_on_reverse_extension, _np_pt_seg_dist, _np_rccw, distance_sq, relative_ccw,
                       seg_in_polygon_many)

__all__ = [
    "CutType", "InflectionType", "Cut", "GeneralInflection", "Bitangent", "CriticalPoint", "inflections",
    "general_inflection_cuts", "bitangent_cuts", "bitangent_lines", "single_tangent_cuts", "all_cuts",
    "all_critical_points", "cut_intersect_points",
]


class CutType(enum.Enum):
    """``Algorithm.CUT_TYPE``."""

    GENERAL_INFLECTION = "GENERAL_INFLECTION"
    NONGENERAL_INFLECTION = "NONGENERAL_INFLECTION"
    SINGLETANGENT = "SINGLETANGENT"
    BITANGENT = "BITANGENT"
    NONE = "NONE"


class InflectionType(enum.Enum):
    """``Algorithm.INFLECTION_TYPE``."""

    GENERAL = "GENERAL"
    NONGENERAL = "NONGENERAL"
    ALL = "ALL"


@dataclass(frozen=True)
class Cut:
    """``cut.Cut``: a critical ray ``line`` (``P1`` = the tangent/reflex vertex, ``P2`` = boundary hit)."""

    line: Line
    type: CutType


@dataclass(frozen=True)
class GeneralInflection(Cut):
    """``cut.GeneralInflection``: extension of edge ``from_line`` past a reflex vertex.

    ``counter_clockwise`` is True for the extension of the incoming edge ``i-1`` beyond ``v_i`` and
    False for the backward extension of the outgoing edge ``i``.
    """

    type: CutType = field(default=CutType.GENERAL_INFLECTION, init=False)
    from_line: int = -1
    counter_clockwise: bool = False


@dataclass(frozen=True)
class Bitangent(Cut):
    """``cut.Bitangent``: the ray leaving the bitangent segment beyond vertex ``this_point``.

    ``opposite_point`` is the other tangent vertex, ``opposite_segment`` the companion ray beyond it,
    and ``curve_to_point`` is ``v_{this_point + 1}`` (decides split vs. merge in ``getGaps``).
    """

    type: CutType = field(default=CutType.BITANGENT, init=False)
    this_point: int = -1
    opposite_point: int = -1
    opposite_segment: Optional[Line] = None
    curve_to_point: Optional[Point] = None


@dataclass(frozen=True)
class CriticalPoint:
    """``cut.PathCutIntersectPoint``: where the path crosses a cut (or its start/end, type ``NONE``).

    ``cut_index`` indexes the cut list given to :func:`all_critical_points` (-1 for start/end) and
    ``seg_index`` the path segment the crossing was found on.
    """

    distance: float
    point: Point
    cut_type: CutType
    cut_index: int
    seg_index: int
    seg: Line
    cut: Optional[Cut] = None


# --------------------------------------------------------------------------- shared vectorised scans


def _nearest_hits(poly: Polygon, rays: np.ndarray, exclude: np.ndarray, mode: str, eps: float):
    """Nearest boundary hits of a batch of lines, exactly as the Java ``for (j ...)`` loops pick them.

    ``rays``: ``(B, 4)`` lines; ``exclude``: ``(B, n)`` bool, edges skipped by identity tests.

    * ``mode="forward"`` (inflection from ``prevLine``): ``p = getIntersect(ray, e_j)``, keep if
      ``e_j.ptSegDist(p) < eps`` and ``onExtension(ray, p)``; distance ``ray.ptSegDist(p)``.
    * ``mode="backward"`` (inflection from ``nextLine``): same with ``onReverseExtension``.
    * ``mode="both"`` (tangent lines): ``p = getIntersect(e_j, ray)``, keep if ``e_j.ptSegDist(p) < eps``,
      split by ``onExtension(ray, p)`` into a ccw and a cw candidate set.

    Each set keeps the first strict minimum in edge order (``if (dist < best)``).  Returns, per set, a
    list of ``(x, y)`` or ``None``.
    """
    ex1, ey1, ex2, ey2 = poly.ex1, poly.ey1, poly.ex2, poly.ey2
    rx1, ry1, rx2, ry2 = (rays[:, k][:, None] for k in range(4))
    if mode == "both":
        x, y, valid = _np_get_intersect(ex1, ey1, ex2, ey2, rx1, ry1, rx2, ry2, eps)
    else:
        x, y, valid = _np_get_intersect(rx1, ry1, rx2, ry2, ex1, ey1, ex2, ey2, eps)
    with np.errstate(invalid="ignore", over="ignore"):
        ok = valid & ~exclude & (_np_pt_seg_dist(ex1, ey1, ex2, ey2, x, y) < eps)
        d = _np_pt_seg_dist(rx1, ry1, rx2, ry2, x, y)
        if mode == "forward":
            sets = [ok & _np_on_extension(rx1, ry1, rx2, ry2, x, y)]
        elif mode == "backward":
            sets = [ok & _np_on_reverse_extension(rx1, ry1, rx2, ry2, x, y)]
        else:
            ext = _np_on_extension(rx1, ry1, rx2, ry2, x, y)
            sets = [ok & ext, ok & ~ext]
    out = []
    for mask in sets:
        dd = np.where(mask, d, np.inf)
        best = np.argmin(dd, axis=1)
        rows = np.arange(len(best))
        has = mask[rows, best]
        out.append([(float(x[r, b]), float(y[r, b])) if h else None for r, b, h in zip(rows, best, has)])
    return out


def _exclude(n: int, idx: Sequence[Sequence[int]]) -> np.ndarray:
    ex = np.zeros((len(idx), n), dtype=bool)
    for r, ks in enumerate(idx):
        ex[r, [k % n for k in ks]] = True
    return ex


# --------------------------------------------------------------------------- inflections


def _inflection_rays(poly: Polygon, eps: float):
    """Every inflection ray of every reflex vertex, in Java output order (cached).

    Records ``(i, from_prev, kind, hit)``: ``from_prev`` is True for case (a) (``prevLine`` = edge
    ``i-1`` extended beyond ``v_i``), False for case (b) (edge ``i`` extended backwards); ``kind`` is
    the turn of the neighbouring vertex (``-1`` general, ``1`` non-general, ``0`` collinear).
    """
    key = ("inflection_rays", eps)
    if key in poly._cache:
        return poly._cache[key]
    n, E, V = poly.n, poly.edges, poly.vertices
    recs = []
    for i in range(n):
        if not poly.is_reflex(i):
            continue
        poi = V[i]
        pp_line, prev_line, next_line, nn_line = E[(i - 2) % n], E[(i - 1) % n], E[i], E[(i + 1) % n]
        recs.append([i, True, relative_ccw(pp_line, poi), prev_line, ((i - 2), (i - 1), i)])
        recs.append([i, False, relative_ccw(nn_line, poi), next_line, ((i + 1), (i - 1), i)])
    out = []
    for from_prev, mode in ((True, "forward"), (False, "backward")):
        sub = [r for r in recs if r[1] == from_prev]
        if not sub:
            continue
        hits = _nearest_hits(poly, np.asarray([r[3] for r in sub], dtype=np.float64),
                             _exclude(n, [r[4] for r in sub]), mode, eps)[0]
        for r, h in zip(sub, hits):
            r.append(h)
    for r in recs:
        out.append((r[0], r[1], r[2], r[5]))
    poly._cache[key] = out
    return out


def _ray_line(poly: Polygon, i: int, hit: Optional[Point]) -> Line:
    if hit is None:
        raise GeometryError(f"inflection ray at vertex {i} of {poly!r} hits no edge (Java NullPointerException)")
    v = poly.vertices[i]
    return (v[0], v[1], hit[0], hit[1])


def inflections(poly: Polygon, kind: InflectionType = InflectionType.ALL, eps: float = EPSILON) -> List[Line]:
    """``Algorithm.getInflection(poly, type)``: the inflection rays as lines ``v_i -> hit``.

    For each reflex ``v_i``, (a) the extension of edge ``i-1`` beyond ``v_i`` -- GENERAL when ``v_{i-1}``
    is convex (``relativeCCW(e_{i-2}, v_i) == -1``), NONGENERAL when it is reflex -- then (b) the
    backward extension of edge ``i``, typed by ``relativeCCW(e_{i+1}, v_i)``.  ``ALL`` keeps every ray,
    collinear neighbours included.  The hit is the nearest edge point (``ptSegDist`` to the extended
    edge) over all edges except the two/three adjacent ones.
    """
    kind = InflectionType(kind)
    want = {InflectionType.GENERAL: (-1,), InflectionType.NONGENERAL: (1,), InflectionType.ALL: (-1, 0, 1)}[kind]
    return [_ray_line(poly, i, hit) for i, _, k, hit in _inflection_rays(poly, eps) if k in want]


def general_inflection_cuts(poly: Polygon, eps: float = EPSILON) -> List[GeneralInflection]:
    """``Algorithm.getGeneralInflectionCut``: the GENERAL rays with ``from_line`` and ``counter_clockwise``.

    ``from_line`` is edge ``i-1`` for case (a) (``counter_clockwise=True``) and edge ``i`` for (b).
    """
    n = poly.n
    return [GeneralInflection(_ray_line(poly, i, hit), from_line=(i - 1) % n if prev else i, counter_clockwise=prev)
            for i, prev, k, hit in _inflection_rays(poly, eps) if k == -1]


# --------------------------------------------------------------------------- tangents


def _pairs(poly: Polygon, i: int) -> np.ndarray:
    """``for (int j = i + 2; j < lines.length && j != i; j++)``."""
    return np.arange(i + 2, poly.n)


def _side_products(poly: Polygon, i: int, js: np.ndarray):
    """``rccw(seg12, v_{i-1}) * rccw(seg12, v_{i+1})`` and the same at ``v_j`` for ``seg12 = v_i -> v_j``."""
    n, vx, vy = poly.n, poly.vx, poly.vy
    x1, y1 = vx[i], vy[i]
    x2, y2 = vx[js], vy[js]
    from_p1 = (_np_rccw(x1, y1, x2, y2, vx[(i - 1) % n], vy[(i - 1) % n])
               * _np_rccw(x1, y1, x2, y2, vx[(i + 1) % n], vy[(i + 1) % n]))
    jm, jp = (js - 1) % n, (js + 1) % n
    from_p2 = _np_rccw(x1, y1, x2, y2, vx[jm], vy[jm]) * _np_rccw(x1, y1, x2, y2, vx[jp], vy[jp])
    return from_p1.astype(np.int64), from_p2.astype(np.int64)


def _tangent_hits(poly: Polygon, pairs: List[Tuple[int, int]], eps: float):
    """``ccwPoint`` / ``cwPoint`` of the segments ``v_i -> v_j`` (edges ``i-1, i, j-1, j`` excluded)."""
    if not pairs:
        return [], []
    V = poly.vertices
    rays = np.asarray([(V[i][0], V[i][1], V[j][0], V[j][1]) for i, j in pairs], dtype=np.float64)
    ex = _exclude(poly.n, [(i - 1, i, j - 1, j) for i, j in pairs])
    ccw, cw = _nearest_hits(poly, rays, ex, "both", eps)
    return ccw, cw


def _bitangent_data(poly: Polygon, eps: float):
    key = ("bitangent", eps)
    if key in poly._cache:
        return poly._cache[key]
    n, V = poly.n, poly.vertices
    ex1, ey1, ex2, ey2 = poly.ex1, poly.ey1, poly.ex2, poly.ey2
    turn = np.asarray(poly.turn)
    pairs: List[Tuple[int, int]] = []
    for i in range(n):
        if turn[i] != 1:
            continue
        js = _pairs(poly, i)
        js = js[turn[js] == 1]
        if js.size == 0:
            continue
        fp1, fp2 = _side_products(poly, i, js)
        js = js[(fp1 == 1) & (fp2 == 1)]
        if js.size == 0:
            continue
        x1, y1 = poly.vx[i], poly.vy[i]
        x2, y2 = poly.vx[js][:, None], poly.vy[js][:, None]
        # seg.intersectsLine(seg12) for every edge other than l11, l12, l21, l22
        r_seg = _np_rccw(x1, y1, x2, y2, ex1, ey1)
        hit = ((r_seg * np.roll(r_seg, -1, axis=1) <= 0)
               & (_np_rccw(ex1, ey1, ex2, ey2, x1, y1) * _np_rccw(ex1, ey1, ex2, ey2, x2, y2) <= 0))
        hit &= ~_exclude(n, [(i - 1, i, j - 1, j) for j in js])
        pairs.extend((i, int(j)) for j in js[~hit.any(axis=1)])
    ccw, cw = _tangent_hits(poly, pairs, eps)
    out: List[Bitangent] = []
    for (i, j), c, w in zip(pairs, ccw, cw):
        if c is None or w is None:
            raise GeometryError(f"bitangent {i}-{j} of {poly!r} hits no edge (Java NullPointerException)")
        p1, p2 = V[i], V[j]
        to_i = V[(i + 1) % n]
        to_j = V[(j + 1) % n]
        if distance_sq(p1, c) < distance_sq(p2, c):
            a, b = (p1[0], p1[1], c[0], c[1]), (p2[0], p2[1], w[0], w[1])
        else:
            a, b = (p1[0], p1[1], w[0], w[1]), (p2[0], p2[1], c[0], c[1])
        out.append(Bitangent(a, this_point=i, opposite_point=j, opposite_segment=b, curve_to_point=to_i))
        out.append(Bitangent(b, this_point=j, opposite_point=i, opposite_segment=a, curve_to_point=to_j))
    out_t = tuple(out)
    poly._cache[key] = out_t
    return out_t


def bitangent_cuts(poly: Polygon, eps: float = EPSILON) -> List[Bitangent]:
    """``Algorithm.getBitangentCut``: two :class:`Bitangent` rays per bitangent pair.

    Pairs ``i < j`` with ``j >= i+2`` (no cyclic wrap) of reflex vertices whose neighbours lie strictly
    on the same side of ``v_i v_j`` at both ends, and whose segment crosses no edge other than the
    four incident ones.  For the line through them the nearest boundary hit past ``v_j``
    (``onExtension``, "ccw point") and the nearest other hit ("cw point") are found.  Emitted, in this
    order: the ray from ``v_i`` (``this=i, opp=j``) and the ray from ``v_j`` (``this=j, opp=i``).
    """
    return list(_bitangent_data(poly, eps))


def bitangent_lines(poly: Polygon, eps: float = EPSILON) -> List[Line]:
    """``Algorithm.getBitangent``: the same rays as :func:`bitangent_cuts` without metadata (panels draw these)."""
    return [b.line for b in _bitangent_data(poly, eps)]


def single_tangent_cuts(poly: Polygon, eps: float = EPSILON) -> List[Cut]:
    """``Algorithm.getSingletangentCut``.

    Pairs ``i < j``, ``j >= i+2``, with at least one reflex end, where the line ``v_i v_j`` is tangent
    at exactly one end and crosses at the other (``fromP1 * fromP2 == -1``) and ``segInPolygon(v_i v_j)``.
    The cut is the ray from the tangent vertex away from the other one, to the nearest boundary hit.
    """
    key = ("single_tangent", eps)
    if key in poly._cache:
        return list(poly._cache[key])
    n, V = poly.n, poly.vertices
    turn = np.asarray(poly.turn)
    cand: List[Tuple[int, int, int]] = []
    for i in range(n):
        js = _pairs(poly, i)
        if turn[i] != 1:
            js = js[turn[js] == 1]
        if js.size == 0:
            continue
        fp1, fp2 = _side_products(poly, i, js)
        keep = fp1 * fp2 == -1
        js, fp1 = js[keep], fp1[keep]
        if js.size == 0:
            continue
        inside = seg_in_polygon_many(poly, np.full(js.size, poly.vx[i]), np.full(js.size, poly.vy[i]),
                                     poly.vx[js], poly.vy[js], eps)
        cand.extend((i, int(j), int(f)) for j, f in zip(js[inside], fp1[inside]))
    ccw, cw = _tangent_hits(poly, [(i, j) for i, j, _ in cand], eps)
    out: List[Cut] = []
    for (i, j, f1), c, w in zip(cand, ccw, cw):
        p1, p2 = V[i], V[j]
        a, b = (p1, p2) if f1 == 1 else (p2, p1)
        end = c if (c is not None and distance_sq(a, c) < distance_sq(b, c)) else w
        if end is None:
            raise GeometryError(f"single tangent {i}-{j} of {poly!r} hits no edge (Java NullPointerException)")
        out.append(Cut((a[0], a[1], end[0], end[1]), CutType.SINGLETANGENT))
    poly._cache[key] = tuple(out)
    return out


def all_cuts(poly: Polygon, eps: float = EPSILON) -> List[Cut]:
    """``Algorithm.getCuts``: single tangents + non-general inflections + general inflections + bitangents."""
    key = ("all_cuts", eps)
    if key not in poly._cache:
        cuts: List[Cut] = list(single_tangent_cuts(poly, eps))
        cuts += [Cut(line, CutType.NONGENERAL_INFLECTION) for line in inflections(poly, InflectionType.NONGENERAL, eps)]
        cuts += general_inflection_cuts(poly, eps)
        cuts += bitangent_cuts(poly, eps)
        poly._cache[key] = tuple(cuts)
    return list(poly._cache[key])  # type: ignore[arg-type]


# --------------------------------------------------------------------------- path crossings


def _cut_crossings(path: Path, cuts: Sequence[Cut], eps: float) -> Dict[float, CriticalPoint]:
    m: Dict[float, CriticalPoint] = {}
    segs = path.segments
    for i, j, p in path._hits([c.line for c in cuts], eps):
        d = path.pt_dist_on_segment(p, i)
        m[d] = CriticalPoint(d, p, cuts[j].type, j, i, segs[i], cuts[j])
    return m


def cut_intersect_points(path: Path, cuts: Sequence[Cut], eps: float = EPSILON) -> List[CriticalPoint]:
    """``Path.getIntersectPoints(Cut[])``: the crossings of the path with the cuts, by arc length, without the
    start and end points of :func:`all_critical_points` (same ``TreeMap`` keys, so B2 applies)."""
    m = _cut_crossings(path, cuts, eps)
    return [m[k] for k in sorted(m)]


def all_critical_points(path: Path, cuts: Sequence[Cut], eps: float = EPSILON) -> List[CriticalPoint]:
    """``Path.getAllCriticalPoints(cuts)``: crossings of the path with the cuts, plus start and end.

    For each path segment ``i`` (outer loop) and cut ``j`` (inner loop), ``getSegmentIntersect`` gives
    the crossing, keyed by ``ptDistFromStart(p, lines, i)`` in a ``TreeMap<Double, ...>``.  Quirk B2
    kept: equal keys keep the **last** put, and ``put(0.0, start)`` / ``put(ptDistFromStart(end), end)``
    (type ``NONE``) overwrite a cut crossing exactly there.  A crossing at a waypoint can be found from
    both adjacent segments with keys differing in the last bits; both survive.
    """
    m = _cut_crossings(path, cuts, eps)
    segs = path.segments
    m[0.0] = CriticalPoint(0.0, path.points[0], CutType.NONE, -1, 0, segs[0])
    end = path.points[-1]
    d = path.pt_dist_from_start(end, eps)
    m[d] = CriticalPoint(d, end, CutType.NONE, -1, len(segs) - 1, segs[-1])
    return [m[k] for k in sorted(m)]
