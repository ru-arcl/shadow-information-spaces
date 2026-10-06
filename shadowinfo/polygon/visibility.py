"""Visibility from a point inside a polygon: physical gaps and the visibility polygon.

Port of ``Algorithm.getPhysicalGaps`` and ``Algorithm.getVisibilityPolygon``
(``docs/notes/original_java.md`` §2.4).  A *gap* (T-RO 2012 Sec. II-A) is a discontinuity of the
visibility polygon; the boundary chain behind it bounds one shadow component.  A
:class:`PhysicalGap` records that hidden chain as a cyclic interval of boundary edges plus the point
where the gap edge meets the boundary.

Both functions scan the edges ``i = 0, 1, ...``.  If ``v_i`` or ``v_{i+1}`` is not visible
(``segInPolygon(q -> v)``), each other vertex ``v_j`` (``j != i, i+1``) is a candidate blocker: the
ray ``q -> v_j`` continued past ``v_j`` must hit edge ``i`` *on* the edge, with ``q -> v_j`` and
``v_j -> hit`` both inside.  The hit nearest to the hidden ``v_i`` is ``ip1`` and the hit nearest to
the hidden ``v_{i+1}`` is ``ip2``; the scan stops as soon as it has what it needs and then jumps past
the edges hidden behind the ``ip2`` blocker.  Every rule (including the early break, the identity
test ``ip1 != ip2`` and the jump) is kept, so the output equals Java's in content and order.

Speed: the vertex visibilities are computed once per query point (Java recomputes them for every
candidate) and the cheap per-candidate tests run as one numpy matrix; the scan itself, with its
order-dependent accept/reject, stays scalar.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Optional, Tuple

import numpy as np

from .geometry import (EPSILON, COMPAT_MODES, GeometryError, Point, Polygon, _check_compat, _np_get_intersect,
                       _np_on_extension, _np_on_reverse_extension, _np_pt_seg_dist, _np_rccw, _unblocked, distance,
                       get_intersect, on_extension, on_reverse_extension, seg_in_polygon_many, seg_in_polygon_scalar)

__all__ = ["PhysicalGap", "vertex_visibility", "physical_gaps", "repair_gaps", "point_status", "visibility_polygon",
           "COMPAT_MODES"]


@dataclass(frozen=True)
class PhysicalGap:
    """``gap.PhysicalGap``: a hidden boundary chain seen from one point.

    Exactly one of ``start_point`` / ``end_point`` is set:

    * ``PhysicalGap(j, i, None, ip)``: hidden chain from vertex ``v_j`` counter-clockwise to the point
      ``ip`` on edge ``i`` (the gap edge leaves ``q`` through ``v_j``);
    * ``PhysicalGap(i, j-1, ip, None)``: hidden chain from ``ip`` on edge ``i`` to vertex ``v_j``.

    (Only the geometry; the ID/link bookkeeping of ``gap.Gap`` lives in the gap tracker.)
    """

    start_edge: int
    end_edge: int
    start_point: Optional[Point] = None
    end_point: Optional[Point] = None

    def full_start_edge(self, n: int) -> int:
        """``getFullStartEdge``: first edge that is *entirely* hidden."""
        if self.start_point is None:
            return self.start_edge
        return 0 if self.start_edge + 1 == n else self.start_edge + 1

    def full_end_edge(self, n: int) -> int:
        """``getFullEndEdge``: last edge that is entirely hidden."""
        if self.end_point is None:
            return self.end_edge
        return n - 1 if self.end_edge == 0 else self.end_edge - 1

    def to_list(self) -> list:
        """``[startEdge, endEdge, startPoint|None, endPoint|None]`` (golden-fixture format)."""
        sp = None if self.start_point is None else list(self.start_point)
        ep = None if self.end_point is None else list(self.end_point)
        return [self.start_edge, self.end_edge, sp, ep]


def vertex_visibility(poly: Polygon, q: Point, eps: float = EPSILON) -> np.ndarray:
    """``vis[j] = segInPolygon(poly, q -> v_j)`` for every vertex, as a bool array.

    Equal to :func:`~shadowinfo.polygon.geometry.seg_in_polygon` on each segment; the parts that
    depend only on the polygon (vertex-on-boundary, ``contains(v_j)``, ``relativeCCW(e_k, v_j)``) are
    cached on the polygon.
    """
    n = poly.n
    qx, qy = float(q[0]), float(q[1])
    vx, vy = poly.vx, poly.vy
    out = np.zeros(n, dtype=bool)
    q_on = bool(np.any(_np_pt_seg_dist(poly.ex1, poly.ey1, poly.ex2, poly.ey2, qx, qy) < eps))
    if not (q_on or poly.contains((qx, qy))):
        return out
    v_on, v_in, rccw_ev = poly._vertex_tables(eps)
    mid_in = poly.contains_many((qx + vx) / 2, (qy + vy) / 2)
    rows = np.flatnonzero(mid_in & (v_on | v_in))
    if rows.size == 0:
        return out
    a1 = np.full(rows.size, qx)
    b1 = np.full(rows.size, qy)
    a2, b2 = vx[rows], vy[rows]
    rccw_seg = _np_rccw(qx, qy, a2[:, None], b2[:, None], vx[None, :], vy[None, :])
    r_q = _np_rccw(poly.ex1, poly.ey1, poly.ex2, poly.ey2, qx, qy)
    rccw_edge = r_q[None, :] * rccw_ev[:, rows].T
    out[rows] = _unblocked(poly, a1, b1, a2, b2, eps, rccw_seg, rccw_edge)
    return out


# One visited edge of the scan: (i, v_i visible, v_{i+1} visible, ip1, ip1n, ip2, ip2n).
_Record = Tuple[int, bool, bool, Optional[Point], int, Optional[Point], int]


def _scan(poly: Polygon, q: Point, eps: float = EPSILON) -> List[_Record]:
    """The edge scan shared by ``getPhysicalGaps`` and ``getVisibilityPolygon``."""
    n = poly.n
    qx, qy = float(q[0]), float(q[1])
    vis = vertex_visibility(poly, q, eps)
    vis_next = np.roll(vis, -1)
    rows = np.flatnonzero(~(vis & vis_next))  # edges not entirely visible
    cols = np.flatnonzero(vis)  # only visible v_j can pass segInPolygon(q -> v_j)
    # Accepted blockers per edge, in j order: (j, ip) with segInPolygon(v_j -> ip) true.  Candidates
    # that fail a test change no state, and the break test after them is false, so they can be dropped.
    accepted: Dict[int, List[Tuple[int, Point]]] = {}
    if rows.size and cols.size:
        cx, cy = poly.vx[cols][None, :], poly.vy[cols][None, :]
        sx1, sy1 = poly.ex1[rows][:, None], poly.ey1[rows][:, None]
        sx2, sy2 = poly.ex2[rows][:, None], poly.ey2[rows][:, None]
        # ip = getIntersect(pp = q -> v_j, seg = e_i); keep if onExtension(pp, ip) and ip lies on e_i
        X, Y, valid = _np_get_intersect(qx, qy, cx, cy, sx1, sy1, sx2, sy2, eps)
        with np.errstate(invalid="ignore", over="ignore"):
            cand = (valid & _np_on_extension(qx, qy, cx, cy, X, Y)
                    & ~_np_on_extension(sx1, sy1, sx2, sy2, X, Y) & ~_np_on_reverse_extension(sx1, sy1, sx2, sy2, X, Y))
        rr, cc = np.nonzero(cand)  # row-major: edges ascending, then j ascending
        ii, js = rows[rr], cols[cc]
        keep = (js != ii) & (js != (ii + 1) % n)
        rr, cc, ii, js = rr[keep], cc[keep], ii[keep], js[keep]
        if js.size:
            ext_ok = seg_in_polygon_many(poly, poly.vx[js], poly.vy[js], X[rr, cc], Y[rr, cc], eps)
            for i_, j_, x_, y_ in zip(ii[ext_ok].tolist(), js[ext_ok].tolist(), X[rr, cc][ext_ok].tolist(),
                                      Y[rr, cc][ext_ok].tolist()):
                accepted.setdefault(i_, []).append((j_, (x_, y_)))
    records: List[_Record] = []
    vis_l = vis.tolist()
    vis_next_l = vis_next.tolist()
    big = float(np.finfo(np.float64).max)
    i = 0
    while i < n:
        pp1_in = vis_l[i]
        pp2_in = vis_next_l[i]
        if pp1_in and pp2_in:
            records.append((i, True, True, None, -1, None, -1))
            i += 1
            continue
        p1 = poly.vertices[i]
        p2 = poly.vertices[(i + 1) % n]
        ip1: Optional[Point] = None
        ip2: Optional[Point] = None
        ip1n = ip2n = -1
        dist1 = dist2 = big
        ni = i
        last_i = i
        for j, ip in accepted.get(i, ()):
            d = distance(ip, p1)
            if d < dist1 and not pp1_in:
                dist1, ip1, ip1n = d, ip, j
            d = distance(ip, p2)
            if d < dist2 and not pp2_in:
                dist2, ip2, ip2n, ni = d, ip, j, j
            # ``ip1 != ip2`` is an object identity test: the same object iff set by the same j
            if (pp1_in and ip2 is not None) or (pp2_in and ip1 is not None) or (
                    ip1 is not None and ip2 is not None and ip1n != ip2n):
                if ni > i:
                    i = ni - 1
                elif ip2 is not None:
                    i = n - 1
                break
        records.append((last_i, pp1_in, pp2_in, ip1, ip1n, ip2, ip2n))
        i += 1
    return records


def _scan_reference(poly: Polygon, q: Point, eps: float = EPSILON) -> List[_Record]:
    """Literal scalar transcription of the Java loop (slow; used by the tests to check :func:`_scan`)."""
    n = poly.n
    lines = poly.edges
    records: List[_Record] = []
    nj = 0
    i = 0
    while i < n:
        seg = lines[i]
        p1, p2 = (seg[0], seg[1]), (seg[2], seg[3])
        pp1_in = seg_in_polygon_scalar(poly, (q[0], q[1], p1[0], p1[1]), eps)
        pp2_in = seg_in_polygon_scalar(poly, (q[0], q[1], p2[0], p2[1]), eps)
        ip1 = ip2 = None
        ip1n = ip2n = -1
        if not pp1_in or not pp2_in:
            dist1 = dist2 = np.finfo(np.float64).max
            ni = i
            nj = i
            last_i = i
            for j in range(n):
                if j == i or (i < n - 1 and j == i + 1) or (i == n - 1 and j == 0) or (nj < j < i):
                    continue
                p = poly.vertices[j]
                pp = (q[0], q[1], p[0], p[1])
                ip = get_intersect(pp, seg, eps)
                if ip is not None and on_extension(pp, ip):
                    if on_extension(seg, ip) or on_reverse_extension(seg, ip):
                        continue
                    if (seg_in_polygon_scalar(poly, pp, eps)
                            and seg_in_polygon_scalar(poly, (p[0], p[1], ip[0], ip[1]), eps)):
                        d = distance(ip, p1)
                        if d < dist1 and not pp1_in:
                            dist1, ip1, ip1n = d, ip, j
                        d = distance(ip, p2)
                        if d < dist2 and not pp2_in:
                            dist2, ip2, ip2n, ni = d, ip, j, j
                    if (pp1_in and ip2 is not None) or (pp2_in and ip1 is not None) or (
                            ip1 is not None and ip2 is not None and ip1n != ip2n):
                        if ni > i:
                            i = ni - 1
                        elif ip2 is not None:
                            i = n - 1
                        break
            records.append((last_i, pp1_in, pp2_in, ip1, ip1n, ip2, ip2n))
        else:
            records.append((i, True, True, None, -1, None, -1))
        i += 1
    return records


def _gaps_from_records(n: int, records: List[_Record]) -> List[PhysicalGap]:
    out: List[PhysicalGap] = []
    for i, _, _, ip1, ip1n, ip2, ip2n in records:
        if ip1 is not None:
            out.append(PhysicalGap(ip1n, i, None, ip1))
        if ip2 is not None:
            out.append(PhysicalGap(i, ip2n - 1 if ip2n != 0 else n - 1, ip2, None))
    return out


def physical_gaps(poly: Polygon, q: Point, eps: float = EPSILON, compat: str = "fixed") -> List[PhysicalGap]:
    """``Algorithm.getPhysicalGaps(poly, q)``: the gaps seen from ``q``, in Java output order.

    Order is counter-clockwise by the edge on which the gap's free end lies, starting at edge 0.
    For every gap at least one edge is fully hidden (the scan skips ``j == i`` and ``j == i+1``).

    ``compat="java"``: Java's scan as is.  It is degenerate where the ray from ``q`` through a reflex
    vertex hits another vertex (exactly or within rounding: every point on the line through two vertices,
    common on the integer maps): the hit is the end of one edge and the start of the next, so the same
    hidden chain can be reported twice, the scan can go the wrong way around the blocker (a phantom gap
    whose "hidden" chain is the visible part of the boundary), or a rounded hit can fall off both edges
    (a gap is missed).  It returns ``[]`` for a point outside the polygon.

    ``compat="fixed"`` (default) raises :class:`~shadowinfo.polygon.geometry.GeometryError` for a point
    strictly outside the polygon, repairs the scan (:func:`_clean_gaps`: a phantom is replaced by its
    complement, duplicates are dropped) and checks the result against the gap structure computed from exact
    orientation signs (:func:`_exact_gaps`; which vertex anchors which gap, on which side, ending on which
    edge).  Where they differ the exact gaps are returned.  At points off those lines the structures agree
    and the result is Java's, bit for bit.
    """
    _check_compat(compat)
    if compat == "java":
        return _gaps_from_records(poly.n, _scan(poly, q, eps))
    on_boundary = _check_inside(poly, q, eps)
    return repair_gaps(poly, q, _gaps_from_records(poly.n, _scan(poly, q, eps)), on_boundary)


def repair_gaps(poly: Polygon, q: Point, raw: List[PhysicalGap], on_boundary: bool = False) -> List[PhysicalGap]:
    """``physical_gaps(poly, q, compat="fixed")`` from ``raw = physical_gaps(poly, q, compat="java")`` (``q`` in
    the closed polygon; ``on_boundary``: within eps of an edge, where only the scan repair applies)."""
    gaps = _clean_gaps(poly, q, raw)
    if not on_boundary and (gaps != raw or _near_vertex_line(poly, q)):
        exact = _exact_gaps(poly, q)
        if exact is not None and _structure(exact) != _structure(gaps):
            return exact
    return gaps


def point_status(poly: Polygon, q: Point, eps: float = EPSILON) -> str:
    """``"inside"``, ``"boundary"`` (within eps of an edge) or ``"outside"`` (where the Java scan returns ``[]``)."""
    try:
        return "boundary" if _check_inside(poly, q, eps) else "inside"
    except GeometryError:
        return "outside"


def _check_inside(poly: Polygon, q: Point, eps: float) -> bool:
    """Raise for a point strictly outside the polygon; return whether it is within eps of the boundary."""
    qx, qy = float(q[0]), float(q[1])
    if not (math.isfinite(qx) and math.isfinite(qy)):
        raise GeometryError(f"point {q} is not finite")
    on = bool(np.any(_np_pt_seg_dist(poly.ex1, poly.ey1, poly.ex2, poly.ey2, qx, qy) < eps))
    if on or poly.contains((qx, qy)):
        return on
    raise GeometryError(f"point {q} is outside polygon {poly.name}")


def _core_point(poly: Polygon, g: PhysicalGap) -> Point:
    """A point just inside the polygon at the middle of the edge at the gap's anchor vertex, which the gap
    claims to hide entirely (the midpoint itself, rounded, may lie just outside the edge)."""
    n = poly.n
    e = g.start_edge if g.end_point is not None else g.end_edge
    x1, y1, x2, y2 = poly.edges[e % n]
    d = 1e-6  # inward normal of a counter-clockwise edge: (-dy, dx), scaled by d * |edge|
    return (0.5 * (x1 + x2) - d * (y2 - y1), 0.5 * (y1 + y2) + d * (x2 - x1))


def _full_edges(g: PhysicalGap, n: int) -> int:
    return (g.full_end_edge(n) - g.full_start_edge(n)) % n + 1


def _clean_gaps(poly: Polygon, q: Point, gaps: List[PhysicalGap]) -> List[PhysicalGap]:
    """The ``compat="fixed"`` repair of a degenerate scan (see :func:`physical_gaps`).

    A genuine gap hides its whole boundary chain, in particular the edge at its anchor vertex, so a gap
    whose anchor edge is visible from ``q`` (closed visibility, exact orientation predicates) is a phantom;
    if its complement (same anchor and free end, the rest of the boundary) hides its anchor edge, that is
    the real gap, otherwise the phantom is dropped.
    Two gaps with the same hidden chain are one gap described twice (the hit is the vertex shared by two
    edges); the description that counts the shared vertex's edge as fully hidden is kept.
    """
    if not gaps:
        return gaps
    from .gaps import hidden_chain  # (gaps imports this module)
    from .simulate import segments_in_polygon

    n = poly.n

    def hidden(gs: List[PhysicalGap]) -> List[bool]:
        cores = [_core_point(poly, g) for g in gs]
        return (~segments_in_polygon(poly, q, [c[0] for c in cores], [c[1] for c in cores])).tolist()

    ok = hidden(gaps)
    if all(ok) and len(gaps) < 2:
        return gaps
    # a phantom describes the complement of a real hidden chain (the scan went the wrong way around the
    # blocker): try the complementary description, with the same anchor and the same free end
    flips = {k: _complement(g, n) for k, g in enumerate(gaps) if not ok[k]}
    flips = {k: g for k, g in flips.items() if g is not None}
    flip_ok = dict(zip(flips, hidden(list(flips.values())))) if flips else {}
    keep: List[PhysicalGap] = []
    for k, g in enumerate(gaps):
        if ok[k]:
            keep.append(g)
        elif flip_ok.get(k):
            keep.append(flips[k])  # type: ignore[arg-type]
    if len(keep) < 2:
        return keep
    per = sum(math.hypot(e[2] - e[0], e[3] - e[1]) for e in poly.edges)
    tol = 1e-9 * max(per, 1.0)
    chains = [hidden_chain(poly, g) for g in keep]
    out: List[PhysicalGap] = []
    used = [False] * len(keep)
    for k, g in enumerate(keep):
        if used[k]:
            continue
        best = k
        for m in range(k + 1, len(keep)):
            if not used[m] and _same_chain(chains[k], chains[m], per, tol):
                used[m] = True
                if _full_edges(keep[m], n) > _full_edges(keep[best], n):
                    best = m
        out.append(keep[best])
    return out


def _complement(g: PhysicalGap, n: int) -> Optional[PhysicalGap]:
    """The gap with the same anchor and free end that hides the rest of the boundary, if it is well formed."""
    if g.end_point is not None:  # chain v_j -> ip on edge i; complement: ip on edge i -> v_j
        i, j = g.end_edge, g.start_edge
        c = PhysicalGap(i, (j - 1) % n, g.end_point, None)
    else:  # chain ip on edge i -> v_{e+1}; complement: v_{e+1} -> ip on edge i
        i, j = g.start_edge, (g.end_edge + 1) % n
        c = PhysicalGap(j, i, None, g.start_point)
    if j == i or j == (i + 1) % n:  # Java never reports these (no fully hidden edge)
        return None
    return c


def _same_chain(a: Tuple[float, float], b: Tuple[float, float], per: float, tol: float) -> bool:
    d = abs(a[0] - b[0]) % per
    return min(d, per - d) <= tol and abs(a[1] - b[1]) <= tol


def visibility_polygon(poly: Polygon, q: Point, compat: str = "fixed", eps: float = EPSILON) -> List[Point]:
    """``Algorithm.getVisibilityPolygon(poly, q, dc)``: vertex list of the visibility polygon of ``q``.

    Per visited edge: ``v_i`` if visible, then ``ip1``, then ``ip2``, then ``v_{i+1}`` if visible.
    Quirk B15: Java emits consecutive duplicates (``v_{i+1}`` of edge ``i`` is ``v_i`` of edge
    ``i+1``).  ``compat="java"`` keeps them; ``compat="fixed"`` (default) drops consecutive repeats,
    including a last vertex equal to the first, and raises :class:`~shadowinfo.polygon.geometry.GeometryError`
    for a point strictly outside the polygon (Java returns an empty or meaningless list there).
    """
    _check_compat(compat)
    if compat == "fixed":
        _check_inside(poly, q, eps)
        return _visibility_from_gaps(poly, physical_gaps(poly, q, eps, "fixed"))
    n = poly.n
    pts: List[Point] = []
    for i, in1, in2, ip1, _, ip2, _ in _scan(poly, q, eps):
        if in1:
            pts.append(poly.vertices[i])
        if ip1 is not None:
            pts.append(ip1)
        if ip2 is not None:
            pts.append(ip2)
        if in2:
            pts.append(poly.vertices[(i + 1) % n])
    if compat == "java":
        return pts
    out: List[Point] = []
    for p in pts:
        if not out or out[-1] != p:
            out.append(p)
    while len(out) > 1 and out[-1] == out[0]:
        out.pop()
    return out


def _visibility_from_gaps(poly: Polygon, gaps: List[PhysicalGap]) -> List[Point]:
    """The visibility polygon as the boundary with every hidden chain replaced by its window.

    Vertices not strictly inside a hidden chain, plus both ends of every chain (the anchor vertex and the
    point ``ip``), in boundary order from ``v_0`` (edge index, then distance from the edge's start, which is
    Java's order: ``v_i``, ``ip1``, ``ip2``, ``v_{i+1}``), consecutive repeats dropped.  With the repaired
    gaps of :func:`physical_gaps` this is also right where Java's scan is degenerate.
    """
    from .gaps import _perimeter, hidden_chain

    n = poly.n
    V = poly.vertices
    cum = _perimeter(poly)
    per = cum[-1]
    chains = [hidden_chain(poly, g) for g in gaps]
    keyed: List[Tuple[int, float, Point]] = []
    for k in range(n):
        if not any(0.0 < (cum[k] - a) % per < ln for a, ln in chains):
            keyed.append((k, 0.0, V[k]))
    for g in gaps:
        if g.end_point is not None:  # chain v_start -> ip on edge end_edge
            keyed.append((g.start_edge, 0.0, V[g.start_edge]))
            ep = tuple(g.end_point)  # type: ignore[arg-type]
            keyed.append((g.end_edge, distance(V[g.end_edge], ep), ep))  # type: ignore[arg-type]
        else:  # chain ip on edge start_edge -> v_{end_edge + 1}
            sp = tuple(g.start_point)  # type: ignore[arg-type]
            keyed.append((g.start_edge, distance(V[g.start_edge], sp), sp))  # type: ignore[arg-type]
            e = (g.end_edge + 1) % n
            keyed.append((e, 0.0, V[e]))
    keyed.sort(key=lambda t: (t[0], t[1]))
    out: List[Point] = []
    for _, _, p in keyed:
        if not out or out[-1] != p:
            out.append(p)
    while len(out) > 1 and out[-1] == out[0]:
        out.pop()
    return out


# --------------------------------------------------------------------------- compat="fixed": exact gap structure


_NEAR = 1e-3  # distance from a line through a reflex vertex and another vertex that triggers the exact check


def _near_vertex_line(poly: Polygon, q: Point) -> bool:
    """``q`` lies within :data:`_NEAR` of a line through a reflex vertex and another vertex, where Java's scan
    can go wrong (elsewhere it never does: the scan only fails when a line of sight through a blocker meets
    another vertex)."""
    key = ("vertex_lines",)
    t = poly._cache.get(key)
    if t is None:
        refl = np.array([k for k, s in enumerate(_exact_turns(poly)) if s < 0], dtype=np.int64)
        rx, ry = poly.vx[refl][:, None], poly.vy[refl][:, None]
        ux, uy = poly.vx[None, :] - rx, poly.vy[None, :] - ry
        ln = np.hypot(ux, uy)
        keep = ln > 0
        t = (np.broadcast_to(rx, ux.shape)[keep], np.broadcast_to(ry, ux.shape)[keep], ux[keep] / ln[keep],
             uy[keep] / ln[keep])
        poly._cache[key] = t
    rx, ry, ux, uy = t  # type: ignore[misc]
    if rx.size == 0:
        return False
    return bool(np.any(np.abs(ux * (q[1] - ry) - uy * (q[0] - rx)) < _NEAR))


def _cross_sign(ax: float, ay: float, bx: float, by: float, cx: float, cy: float, dx: float, dy: float) -> int:
    """``sign((b - a) x (d - c))``, exact on the float inputs (float filter, rational fallback)."""
    v = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx)
    if abs(v) > 1e-9 * (abs(bx - ax) + abs(by - ay)) * (abs(dx - cx) + abs(dy - cy)):
        return 1 if v > 0 else -1
    f = Fraction
    w = (f(bx) - f(ax)) * (f(dy) - f(cy)) - (f(by) - f(ay)) * (f(dx) - f(cx))
    return (w > 0) - (w < 0)


def _orient_signs(ax: float, ay: float, bx: float, by: float, xs: np.ndarray, ys: np.ndarray,
                  skip: int = -1) -> np.ndarray:
    """``sign((b - a) x (p - a))`` for every point ``p``, exact (see :func:`_cross_sign`); ``p = xs[skip]`` is
    ``b`` itself (sign 0)."""
    o = (bx - ax) * (ys - ay) - (by - ay) * (xs - ax)
    bound = 1e-9 * (abs(bx - ax) + abs(by - ay)) * (np.abs(xs - ax) + np.abs(ys - ay))
    sg = np.sign(o).astype(np.int64)
    for k in np.flatnonzero(np.abs(o) <= bound).tolist():
        sg[k] = 0 if k == skip else _cross_sign(ax, ay, bx, by, ax, ay, float(xs[k]), float(ys[k]))
    return sg


def _exact_turns(poly: Polygon) -> List[int]:
    """``sign((v_k - v_{k-1}) x (v_{k+1} - v_k))``: +1 convex, -1 reflex, 0 straight (exact)."""
    key = ("exact_turns",)
    t = poly._cache.get(key)
    if t is None:
        V, n = poly.vertices, poly.n
        t = [_cross_sign(*V[k - 1], *V[k], *V[k], *V[(k + 1) % n]) for k in range(n)]
        poly._cache[key] = t
    return t  # type: ignore[return-value]


def _cone(poly: Polygon, k: int, fx: float, fy: float, tx: float, ty: float) -> bool:
    """The direction ``f -> t`` points into the closed interior angle at ``v_k`` (exact)."""
    n = poly.n
    u, v, w = poly.vertices[k - 1], poly.vertices[k], poly.vertices[(k + 1) % n]
    c1 = _cross_sign(*u, *v, fx, fy, tx, ty)
    c2 = _cross_sign(*v, *w, fx, fy, tx, ty)
    t = _exact_turns(poly)[k]
    if t > 0:
        return c1 >= 0 and c2 >= 0
    if t < 0:
        return c1 >= 0 or c2 >= 0
    return c2 >= 0


def _exact_gaps(poly: Polygon, q: Point) -> Optional[List[PhysicalGap]]:
    """The gaps seen from a point strictly inside the polygon, from exact orientation signs (closed visibility).

    A vertex ``r`` visible from ``q`` (no edge crosses ``q r`` properly, and at every vertex on it both
    directions point into the closed interior angle) anchors a gap iff the line of sight continues past ``r``
    into the polygon with the boundary at ``r`` on one side ``s`` of it (both neighbours on side ``s``, or one
    on side ``s`` and the other on the line behind ``r``).  The hidden chain leaves ``r`` along the neighbour
    closer in angle to the line of sight, counter-clockwise (``PhysicalGap(r, i, None, ip)``) or clockwise
    (``PhysicalGap(i, r-1, ip, None)``), and ends at the first point past ``r`` where the boundary meets the
    line of sight from side ``s``: a proper crossing of edge ``i`` (``ip = getIntersect(q -> r, e_i)``, Java's
    arithmetic) or a vertex ``v_k`` touched from side ``s`` (``ip = v_k``, on the edge that keeps the hidden
    chain's full edges complete).  Returns None if the structure is not consistent (the caller then keeps
    the repaired Java scan).
    """
    n = poly.n
    V = poly.vertices
    xs, ys = poly.vx, poly.vy
    qx, qy = float(q[0]), float(q[1])
    out: List[Tuple[int, float, PhysicalGap]] = []
    # cheap prefilter: an anchor is not convex, and its neighbours are not clearly on opposite sides
    dxs, dys = xs - qx, ys - qy
    xa, ya, xb, yb = np.roll(xs, 1), np.roll(ys, 1), np.roll(xs, -1), np.roll(ys, -1)
    oa = dxs * (ya - qy) - dys * (xa - qx)
    ob = dxs * (yb - qy) - dys * (xb - qx)
    sz = np.abs(dxs) + np.abs(dys)
    clear = (np.abs(oa) > 1e-9 * sz * (np.abs(xa - qx) + np.abs(ya - qy))) & (
        np.abs(ob) > 1e-9 * sz * (np.abs(xb - qx) + np.abs(yb - qy)))
    turns = np.array(_exact_turns(poly))
    cands = np.flatnonzero((turns <= 0) & ~(clear & (oa * ob < 0))).tolist()
    for r in cands:
        rx, ry = V[r]
        dx, dy = rx - qx, ry - qy
        if dx == 0.0 and dy == 0.0:
            continue
        sg = _orient_signs(qx, qy, rx, ry, xs, ys, r)
        a, b = (r - 1) % n, (r + 1) % n
        ahead = (xs - rx) * dx + (ys - ry) * dy  # > 0: past r along the line of sight
        sides = []
        for u in (a, b):
            if sg[u] != 0:
                sides.append(int(sg[u]))
            elif ahead[u] > 0:
                sides = []  # an edge runs on along the line of sight: the gap (if any) is at its far end
                break
        if not sides or any(x != sides[0] for x in sides) or not _cone(poly, r, qx, qy, rx, ry):
            continue
        if not _sees(poly, q, r, sg):
            continue
        side = sides[0]
        # the hidden neighbour: the one on side ``side`` closer in angle to the forward direction
        cand = [u for u in (a, b) if sg[u] == side]
        if len(cand) == 2:
            ua, ub = cand
            da = (xs[ua] - rx) * dx + (ys[ua] - ry) * dy
            db = (xs[ub] - rx) * dx + (ys[ub] - ry) * dy
            ca = abs(dx * (ys[ua] - ry) - dy * (xs[ua] - rx))
            cb = abs(dx * (ys[ub] - ry) - dy * (xs[ub] - rx))
            hid = ua if da * cb > db * ca else ub
        else:
            hid = cand[0]
        ccw = hid == b
        # first contact past r from side ``side``: a vertex on the line touched from that side, or a proper
        # crossing of an edge (no edge crosses between q and r, which is visible)
        best_t = math.inf
        best: Optional[Tuple[str, int]] = None
        dd = dx * dx + dy * dy
        nx = _next_index(n)
        sg_next = sg[nx]
        sg_prev = np.empty_like(sg)
        sg_prev[nx] = sg
        touch = (sg == 0) & (ahead > 0) & ((sg_prev == side) | (sg_next == side))
        touch[r] = False
        if touch.any():
            ks = np.flatnonzero(touch)
            tv = ((xs[ks] - qx) * dx + (ys[ks] - qy) * dy) / dd
            m = int(np.argmin(tv))
            best_t, best = float(tv[m]), ("v", int(ks[m]))
        cross = sg * sg_next < 0
        cross[r] = False
        cross[(r - 1) % n] = False
        if cross.any():
            ks = np.flatnonzero(cross)
            ex, ey = poly.ex2[ks] - xs[ks], poly.ey2[ks] - ys[ks]
            den = dx * ey - dy * ex
            with np.errstate(divide="ignore", invalid="ignore"):
                te = ((xs[ks] - qx) * ey - (ys[ks] - qy) * ex) / den
            for m in np.argsort(te).tolist():
                t = float(te[m])
                if not (t > 0.0) or t >= best_t:
                    continue
                k = int(ks[m])
                k1 = (k + 1) % n
                if _cross_sign(*V[k], *V[k1], *V[k], qx, qy) * _cross_sign(*V[k], *V[k1], *V[k], rx, ry) < 0:
                    continue  # between q and r (cannot happen for a visible r, up to rounding of t)
                best_t, best = t, ("e", k)
                break
        if best is None:
            return None
        kind, k = best
        if kind == "e":
            ip = get_intersect((qx, qy, rx, ry), poly.edges[k])
            x1, y1, x2, y2 = poly.edges[k]
            if ip is None or not (min(x1, x2) <= ip[0] <= max(x1, x2) and min(y1, y2) <= ip[1] <= max(y1, y2)):
                ip = (qx + best_t * dx, qy + best_t * dy)
            e_end, e_start = k, k
        else:
            ip = V[k]
            e_end, e_start = k, (k - 1) % n
        if ccw:
            if e_end == r or e_end == (r - 1) % n:
                return None
            g = PhysicalGap(r, e_end, None, ip)
            out.append((e_end, distance(V[e_end], ip), g))
        else:
            if (r - 1) % n == e_start or r == e_start:
                return None
            g = PhysicalGap(e_start, (r - 1) % n, ip, None)
            out.append((e_start, distance(V[e_start], ip), g))
    out.sort(key=lambda t: (t[0], t[1]))
    return [g for _, _, g in out]


def _signs_vs_edges(poly: Polygon, ks: np.ndarray, px: float, py: float) -> np.ndarray:
    """``sign((v_{k+1} - v_k) x (p - v_k))`` for the edges ``ks``, exact (see :func:`_cross_sign`)."""
    x1, y1, x2, y2 = poly.ex1[ks], poly.ey1[ks], poly.ex2[ks], poly.ey2[ks]
    o = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)
    bound = 1e-9 * (np.abs(x2 - x1) + np.abs(y2 - y1)) * (np.abs(px - x1) + np.abs(py - y1))
    sg = np.sign(o).astype(np.int64)
    for m in np.flatnonzero(np.abs(o) <= bound).tolist():
        k = int(ks[m])
        sg[m] = _cross_sign(*poly.edges[k][:2], *poly.edges[k][2:], *poly.edges[k][:2], px, py)
    return sg


def _sees(poly: Polygon, q: Point, r: int, sg: np.ndarray) -> bool:
    """Closed visibility of ``v_r`` from ``q`` (exact): ``sg`` are the orientation signs of the vertices
    with respect to the line ``q -> v_r``."""
    n = poly.n
    qx, qy = float(q[0]), float(q[1])
    rx, ry = poly.vertices[r]
    dx, dy = rx - qx, ry - qy
    ks = np.flatnonzero(sg * sg[_next_index(n)] < 0)  # edges straddling the line
    ks = ks[(ks != r) & (ks != (r - 1) % n)]
    if ks.size:
        if np.any(_signs_vs_edges(poly, ks, qx, qy) * _signs_vs_edges(poly, ks, rx, ry) < 0):
            return False
    along = (poly.vx - qx) * dx + (poly.vy - qy) * dy
    for k in np.flatnonzero((sg == 0) & (along > 0) & (along < dx * dx + dy * dy)).tolist():  # on the open segment
        if k != r and not (_cone(poly, k, qx, qy, rx, ry) and _cone(poly, k, rx, ry, qx, qy)):
            return False
    return True


def _next_index(n: int) -> np.ndarray:
    a = _NEXT.get(n)
    if a is None:
        a = _NEXT[n] = (np.arange(n) + 1) % n
    return a


_NEXT: Dict[int, np.ndarray] = {}


def _structure(gaps: List[PhysicalGap]) -> List[Tuple[int, int, bool]]:
    return sorted((g.start_edge, g.end_edge, g.end_point is not None) for g in gaps)
