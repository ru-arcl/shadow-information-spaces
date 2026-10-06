"""Polygon world: a robot on a path, random-walk targets, shadow pockets, critical events and ground truth.

New functionality built on the ported engine (the Java code has no geometric targets and no FOV
events, ``docs/notes/original_java.md`` §1 and B13).  It is the polygon counterpart of
:mod:`shadowinfo.grid.simulate` and uses the same tick semantics (``docs/DESIGN.md`` §5) and the same
event types (:mod:`shadowinfo.events`), so its output feeds every filter of the package unchanged.

**Robot and shadows.**  The robot moves along a polyline path by a fixed arc length per tick
(``speed``).  The shadows are the gaps of :func:`~shadowinfo.polygon.gaps.get_gaps` with
``compat="fixed", java_matching=False`` on that path (every step checked against the hidden boundary
chains, so a label never moves to another physical shadow; the ground truth below relies on that): the
critical points crossed during a tick are processed in order, each with the component events and the
labels ``get_gaps`` gives at that critical point (T-RO 2012 Sec. II-A), so the component events and
labels of a run are exactly those of ``get_gaps(..., java_matching=False)`` on the same path.  At the
robot position itself the gaps are recomputed (:func:`~shadowinfo.polygon.visibility.physical_gaps`)
and take their labels from the last critical point crossed.  When the path ends the robot stops, runs
the path again (closed paths) or runs it backwards; :meth:`PolygonSimulator.go_to` / :meth:`~PolygonSimulator.follow` start a new path from
the current position (click-to-move).  Each new path is tracked by ``get_gaps`` on its own and joined
to the current shadows through the hidden boundary chains (:func:`~shadowinfo.polygon.gaps.chain_graph`).

**Shadow pockets.**  The pocket of a gap is the part of the polygon cut off by the gap's window
(the segment from the gap's anchor vertex to the point where the gap edge meets the boundary), on
the side away from the robot: the hidden boundary chain closed by the window
(:func:`pocket_polygon`).  The pockets are disjoint, lie in the polygon and, with the visibility
polygon, tile it.

**Targets.**  ``n_targets`` points placed uniformly in the polygon and doing a random walk: each
tick, each target draws ``dx, dy`` uniform in ``[-target_step, target_step)`` from
:class:`shadowinfo.rng.Rng` (two draws per target and tick, accepted or not) and moves unless the
step leaves the polygon or touches its boundary.  A target is *visible* iff the segment from the
sensing point to it lies in the closed polygon (:func:`segment_in_polygon`); a hidden target belongs
to the pocket that contains it.  These predicates use only ``+ - * /`` and comparisons on doubles, so
a JS port gives the same bits.

**One tick** (``docs/DESIGN.md`` §5):

1. *sensor substep* -- the robot advances.  For every critical point crossed, then for the final
   position: targets that were in a shadow and are now visible (or now in a shadow that is not a
   descendant) emit ``Exit(old)`` before the component events, unless their shadow disappears, in
   which case they are counted by ``Disappear(s, k, k)``; then the component events; then targets
   that were visible (or in an unrelated shadow) and are now in a shadow emit ``Enter(new)``, unless
   the shadow just appeared, in which case they are counted by ``Appear(s, k, k)``.
2. *target substep* -- the targets move; shadow -> visible emits ``Exit``, visible -> shadow
   ``Enter`` (and shadow -> another shadow both).

FOV events are aggregated per shadow (``Exit(s, k)``): within each step all exits sorted by label,
then the component events, then all enters sorted by label.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ..events import Appear, Bound, Disappear, Enter, Event, Exit, Merge, ShadowSequence, Split, event_to_dict
from ..grid.simulate import initial_bounds
from ..rng import Rng
from .events import transition_events
from .gaps import Gap, GapHistory, _components, _from_graph, chain_graph, get_gaps
from .geometry import Path, Point, Polygon
from .io import load_demo_path
from .visibility import PhysicalGap, physical_gaps, visibility_polygon

__all__ = [
    "AT_END", "Shadow", "PolygonSimulator", "simulate", "orient", "point_in_polygon", "points_in_polygon",
    "segment_in_polygon", "segments_in_polygon", "segment_touches_boundary", "segment_along_edge", "validate_path",
    "pocket_polygon", "gap_window", "shadow_pockets", "polygon_area",
]

AT_END = ("stop", "loop", "reverse")
"""What the robot does at the end of its path: stay there, run the (closed) path again, or run it backwards."""


# --------------------------------------------------------------------------- predicates (JS-mirrorable)


def orient(ax: float, ay: float, bx: float, by: float, cx: float, cy: float) -> float:
    """``(b - a) x (c - a)``: positive when ``c`` is to the left of ``a -> b`` (y-up frame)."""
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def _coords(poly_or_pts) -> Tuple[np.ndarray, np.ndarray]:
    if isinstance(poly_or_pts, Polygon):
        return poly_or_pts.vx, poly_or_pts.vy
    a = np.asarray(poly_or_pts, dtype=np.float64).reshape(-1, 2)
    return a[:, 0].copy(), a[:, 1].copy()


def point_in_polygon(poly_or_pts, p: Point) -> bool:
    """Even-odd crossing test (the classic ``(yi > py) != (yj > py)`` rule); ``poly_or_pts`` is a
    :class:`~shadowinfo.polygon.geometry.Polygon` or a vertex list.  Points on the boundary may go
    either way; the simulator never relies on them."""
    xs, ys = _coords(poly_or_pts)
    px, py = float(p[0]), float(p[1])
    inside = False
    n = len(xs)
    j = n - 1
    for i in range(n):
        xi, yi, xj, yj = float(xs[i]), float(ys[i]), float(xs[j]), float(ys[j])
        if (yi > py) != (yj > py) and px < (xj - xi) * (py - yi) / (yj - yi) + xi:
            inside = not inside
        j = i
    return inside


def points_in_polygon(poly_or_pts, px, py) -> np.ndarray:
    """Vectorised :func:`point_in_polygon` (same operations elementwise, so the same answers)."""
    xs, ys = _coords(poly_or_pts)
    return _pip_arrays(xs, ys, np.asarray(px, dtype=np.float64), np.asarray(py, dtype=np.float64))


def _pip_arrays(xs: np.ndarray, ys: np.ndarray, px: np.ndarray, py: np.ndarray) -> np.ndarray:
    px = px.reshape(-1, 1)
    py = py.reshape(-1, 1)
    xj, yj = np.roll(xs, 1), np.roll(ys, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        hit = ((ys > py) != (yj > py)) & (px < (xj - xs) * (py - ys) / (yj - ys) + xs)
    return hit.sum(axis=1) % 2 == 1


def _turns(poly: Polygon) -> np.ndarray:
    """``orient(v_{k-1}, v_k, v_{k+1})`` sign per vertex: ``+1`` convex, ``-1`` reflex, ``0`` straight."""
    key = ("sim_turns",)
    t = poly._cache.get(key)
    if t is None:
        n = poly.n
        V = poly.vertices
        t = np.array([np.sign(orient(V[k - 1][0], V[k - 1][1], V[k][0], V[k][1], V[(k + 1) % n][0],
                                     V[(k + 1) % n][1])) for k in range(n)], dtype=np.int8)
        poly._cache[key] = t
    return t  # type: ignore[return-value]


def _in_cone(poly: Polygon, k: int, dx: float, dy: float) -> bool:
    """Direction ``(dx, dy)`` points into the closed interior angle of the (CCW) polygon at ``v_k``."""
    n = poly.n
    ux, uy = poly.vertices[k - 1]
    vx, vy = poly.vertices[k]
    wx, wy = poly.vertices[(k + 1) % n]
    c1 = (vx - ux) * dy - (vy - uy) * dx
    c2 = (wx - vx) * dy - (wy - vy) * dx
    t = _turns(poly)[k]
    if t > 0:
        return c1 >= 0 and c2 >= 0
    if t < 0:
        return c1 >= 0 or c2 >= 0
    return c2 >= 0


def _through_vertex_ok(poly: Polygon, ax: float, ay: float, bx: float, by: float, k: int) -> bool:
    """For ``v_k`` on the line ``a b``: False iff ``v_k`` is strictly inside the segment and the segment
    leaves the polygon there (one of its two directions points outside at ``v_k``)."""
    vx, vy = poly.vertices[k]
    dx, dy = bx - ax, by - ay
    if not ((vx - ax) * dx + (vy - ay) * dy > 0 and (vx - bx) * dx + (vy - by) * dy < 0):
        return True
    return _in_cone(poly, k, dx, dy) and _in_cone(poly, k, -dx, -dy)


def segment_in_polygon(poly: Polygon, a: Point, b: Point) -> bool:
    """The segment ``a b`` lies in the *closed* polygon (closed visibility), for ``a`` and ``b`` in it.

    No edge crosses the segment properly (all four orientations nonzero with opposite signs), and at
    every vertex strictly inside the segment both directions of the segment point into the polygon's
    closed interior angle.  Grazing a reflex vertex or running along an edge is allowed.
    """
    ax, ay, bx, by = float(a[0]), float(a[1]), float(b[0]), float(b[1])
    n = poly.n
    V = poly.vertices
    o = [orient(ax, ay, bx, by, V[k][0], V[k][1]) for k in range(n)]
    for k in range(n):
        o1, o2 = o[k], o[(k + 1) % n]
        if (o1 < 0 < o2) or (o2 < 0 < o1):
            px, py = V[k]
            rx, ry = V[(k + 1) % n]
            o3 = orient(px, py, rx, ry, ax, ay)
            o4 = orient(px, py, rx, ry, bx, by)
            if (o3 < 0 < o4) or (o4 < 0 < o3):
                return False
    for k in range(n):
        if o[k] == 0 and not _through_vertex_ok(poly, ax, ay, bx, by, k):
            return False
    return True


def segments_in_polygon(poly: Polygon, a: Point, bx, by) -> np.ndarray:
    """Vectorised :func:`segment_in_polygon` from one point ``a`` to many points ``(bx[j], by[j])``."""
    ax, ay = float(a[0]), float(a[1])
    bx = np.asarray(bx, dtype=np.float64).reshape(-1, 1)
    by = np.asarray(by, dtype=np.float64).reshape(-1, 1)
    vx, vy, ex2, ey2 = poly.vx, poly.vy, poly.ex2, poly.ey2
    o = (bx - ax) * (vy - ay) - (by - ay) * (vx - ax)  # orient(a, b, v_k)
    o2 = np.roll(o, -1, axis=1)
    o3 = (ex2 - vx) * (ay - vy) - (ey2 - vy) * (ax - vx)  # orient(v_k, v_k+1, a)
    o4 = (ex2 - vx) * (by - vy) - (ey2 - vy) * (bx - vx)  # orient(v_k, v_k+1, b)
    proper = (((o < 0) & (o2 > 0)) | ((o2 < 0) & (o > 0))) & (((o3 < 0) & (o4 > 0)) | ((o4 < 0) & (o3 > 0)))
    ok = ~proper.any(axis=1)
    rows, ks = np.nonzero((o == 0) & ok[:, None])
    for r, k in zip(rows.tolist(), ks.tolist()):
        if ok[r] and not _through_vertex_ok(poly, ax, ay, float(bx[r, 0]), float(by[r, 0]), k):
            ok[r] = False
    return ok


def segment_touches_boundary(poly: Polygon, a: Point, b: Point) -> bool:
    """Closed segment-edge intersection test against every edge (touching counts; collinear edges
    count even when disjoint, which only makes the caller more cautious)."""
    return bool(_touch_many(poly, np.array([a[0]]), np.array([a[1]]), np.array([b[0]]), np.array([b[1]]))[0])


def _touch_many(poly: Polygon, ax, ay, bx, by) -> np.ndarray:
    ax, ay, bx, by = (np.asarray(v, dtype=np.float64).reshape(-1, 1) for v in (ax, ay, bx, by))
    vx, vy, ex2, ey2 = poly.vx, poly.vy, poly.ex2, poly.ey2
    o1 = (bx - ax) * (vy - ay) - (by - ay) * (vx - ax)
    o2 = (bx - ax) * (ey2 - ay) - (by - ay) * (ex2 - ax)
    o3 = (ex2 - vx) * (ay - vy) - (ey2 - vy) * (ax - vx)
    o4 = (ex2 - vx) * (by - vy) - (ey2 - vy) * (bx - vx)
    s12 = ((o1 <= 0) & (o2 >= 0)) | ((o1 >= 0) & (o2 <= 0))
    s34 = ((o3 <= 0) & (o4 >= 0)) | ((o3 >= 0) & (o4 <= 0))
    return (s12 & s34).any(axis=1)


def polygon_area(pts: Sequence[Point]) -> float:
    """Signed shoelace area (positive counter-clockwise)."""
    s = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return 0.5 * s


def validate_path(poly: Polygon, waypoints: Sequence[Sequence[float]]) -> List[Point]:
    """Check a robot path: at least 2 waypoints, each strictly inside the polygon (even-odd test and not
    on an edge), each segment inside the closed polygon (:func:`segment_in_polygon`) and not running along
    an edge (:func:`segment_along_edge`).  Returns the waypoints as tuples; raises ``ValueError`` naming the
    first violation.

    A segment may still pass through a reflex vertex (closed visibility).  A segment along an edge is
    refused because every critical line of that edge (its inflection cuts) is collinear with it, so the
    crossings cannot be ordered along the segment and every point of it is on the boundary, where the
    physical gaps are degenerate: the tracker would compress the events into wrong instants and the
    simulator's ground truth would be wrong (move the segment off the wall by any amount instead)."""
    pts = [(float(p[0]), float(p[1])) for p in waypoints]
    if len(pts) < 2:
        raise ValueError("a path needs at least 2 waypoints")
    for i, p in enumerate(pts):
        if not all(math.isfinite(c) for c in p):
            raise ValueError(f"waypoint {i} {p} is not finite")
        if not point_in_polygon(poly, p) or _on_boundary(poly, p):
            raise ValueError(f"waypoint {i} {p} is not strictly inside polygon {poly.name}")
    for i in range(len(pts) - 1):
        if pts[i] == pts[i + 1]:
            raise ValueError(f"segment {i} has zero length")
        if not segment_in_polygon(poly, pts[i], pts[i + 1]):
            raise ValueError(f"segment {i} {pts[i]} -> {pts[i + 1]} leaves polygon {poly.name}")
        k = segment_along_edge(poly, pts[i], pts[i + 1])
        if k is not None:
            raise ValueError(f"segment {i} {pts[i]} -> {pts[i + 1]} runs along edge {k} of polygon {poly.name}")
    return pts


def segment_along_edge(poly: Polygon, a: Point, b: Point) -> Optional[int]:
    """The first edge that the segment ``a b`` overlaps in more than a point (collinear), or None."""
    ax, ay, bx, by = float(a[0]), float(a[1]), float(b[0]), float(b[1])
    dx, dy = bx - ax, by - ay
    dd = dx * dx + dy * dy
    V = poly.vertices
    n = poly.n
    for k in range(n):
        (px, py), (rx, ry) = V[k], V[(k + 1) % n]
        if orient(ax, ay, bx, by, px, py) != 0 or orient(ax, ay, bx, by, rx, ry) != 0:
            continue
        t1 = (px - ax) * dx + (py - ay) * dy
        t2 = (rx - ax) * dx + (ry - ay) * dy
        if min(max(t1, t2), dd) > max(min(t1, t2), 0.0):
            return k
    return None


def _on_boundary(poly: Polygon, p: Point) -> bool:
    px, py = p
    vx, vy, ex2, ey2 = poly.vx, poly.vy, poly.ex2, poly.ey2
    o = (ex2 - vx) * (py - vy) - (ey2 - vy) * (px - vx)
    within = ((px - vx) * (ex2 - vx) + (py - vy) * (ey2 - vy) >= 0) & ((px - ex2) * (vx - ex2) + (py - ey2) * (vy - ey2) >= 0)
    return bool(np.any((o == 0) & within))


# --------------------------------------------------------------------------- pockets


def _cyclic(a: int, b: int, n: int) -> List[int]:
    out = [a]
    while out[-1] != b:
        out.append((out[-1] + 1) % n)
    return out


def pocket_polygon(poly: Polygon, g: PhysicalGap) -> List[Point]:
    """The shadow pocket of a gap: its hidden boundary chain closed by the window, counter-clockwise.

    ``PhysicalGap(j, i, None, ip)``: ``v_j, v_{j+1}, ..., v_i, ip`` (``ip`` on edge ``i``, anchor ``v_j``);
    ``PhysicalGap(i, e, ip, None)``: ``ip, v_{i+1}, ..., v_{e+1}`` (``ip`` on edge ``i``, anchor ``v_{e+1}``).
    """
    n = poly.n
    V = poly.vertices
    if g.end_point is not None:
        return [V[k] for k in _cyclic(g.start_edge, g.end_edge, n)] + [tuple(g.end_point)]  # type: ignore[list-item]
    return [tuple(g.start_point)] + [V[k] for k in _cyclic((g.start_edge + 1) % n, (g.end_edge + 1) % n, n)]  # type: ignore[arg-type,list-item]


def gap_window(poly: Polygon, g: PhysicalGap) -> Tuple[Point, Point]:
    """``(anchor, free end)`` of a gap: the window segment from the occluding reflex vertex to the point
    where the gap edge meets the boundary."""
    n = poly.n
    if g.end_point is not None:
        return poly.vertices[g.start_edge], tuple(g.end_point)  # type: ignore[return-value]
    return poly.vertices[(g.end_edge + 1) % n], tuple(g.start_point)  # type: ignore[return-value]


@dataclass(frozen=True)
class Shadow:
    """One shadow component seen from a point: its label, gap, pocket polygon and window."""

    label: int
    gap: PhysicalGap
    pocket: Tuple[Point, ...]
    window: Tuple[Point, Point]

    @property
    def area(self) -> float:
        return polygon_area(self.pocket)

    def to_dict(self) -> dict:
        return {"label": self.label, "window": [list(p) for p in self.window],
                "pocket": [list(p) for p in self.pocket]}


def shadow_pockets(poly: Polygon, q: Point, gaps: Optional[Sequence[PhysicalGap]] = None,
                   labels: Optional[Sequence[int]] = None) -> List[Shadow]:
    """The shadows seen from ``q`` (``physical_gaps`` order), labelled ``labels`` or ``1, 2, ...``."""
    gaps = physical_gaps(poly, q) if gaps is None else list(gaps)
    labels = list(range(1, len(gaps) + 1)) if labels is None else list(labels)
    return [Shadow(s, g, tuple(pocket_polygon(poly, g)), gap_window(poly, g)) for s, g in zip(labels, gaps)]


# --------------------------------------------------------------------------- legs and configurations


@dataclass
class _Track:
    """``get_gaps(compat="fixed", java_matching=False)`` of one path, as the simulator consumes it (cached per
    polygon+path).  The hidden-chain matching keeps every label on its physical shadow, which the target ground
    truth needs; Java's own label rotations (quirk B3) would move targets between labels without an event."""

    path: Path
    history: GapHistory
    distances: List[float]
    events_at: Dict[int, List[Event]]  # critical index -> component events (local labels)
    m: int  # number of shadows at the start (local labels 1..m)
    max_local: int


def _track(poly: Polygon, path: Path) -> _Track:
    key = ("sim_track", path.points)
    tr = poly._cache.get(key)
    if tr is None:
        h = get_gaps(poly, path, "fixed", java_matching=False)
        per = transition_events(h)
        events_at: Dict[int, List[Event]] = {}
        max_local = h.max_id
        for k, evs in enumerate(per):
            events_at[h.sample_indices[k] + 1] = list(evs)
            for e in evs:
                max_local = max(max_local, *[v for kk, v in e.__dict__.items() if kk in ("s", "a", "b")])
        m = len(h.sample_ids[0])
        assert h.sample_ids[0] == list(range(1, m + 1))
        tr = _Track(path, h, [cp.distance for cp in h.samples.critical_points], events_at, m, max_local)  # type: ignore[union-attr]
        poly._cache[key] = tr
    return tr  # type: ignore[return-value]


@dataclass
class _Leg:
    track: _Track
    init_map: Dict[int, int]
    base: int
    next_cp: int = 1
    s: float = 0.0

    def g(self, local: int) -> int:
        return self.init_map[local] if local <= self.track.m else local - self.track.m - 1 + self.base

    def events(self, i: int) -> List[Event]:
        return [_relabel(e, self.g) for e in self.track.events_at.get(i, ())]

    def labels(self, i: int) -> List[int]:
        return [self.g(x) for x in self.track.history.sample_ids[i]]

    @property
    def length(self) -> float:
        return self.track.path.length()

    def point_at(self, s: float) -> Point:
        path = self.track.path
        pre = path._prefix
        if s >= pre[-1]:
            return path.points[-1]
        i = 0
        while i + 1 < len(pre) - 1 and pre[i + 1] <= s:
            i += 1
        x1, y1, x2, y2 = path.segments[i]
        f = (s - pre[i]) / (pre[i + 1] - pre[i])
        return (x1 + f * (x2 - x1), y1 + f * (y2 - y1))


def _relabel(e: Event, g) -> Event:
    if isinstance(e, Appear):
        return Appear(g(e.s), e.lo, e.hi)
    if isinstance(e, Disappear):
        return Disappear(g(e.s), e.lo, e.hi)
    if isinstance(e, Split):
        return Split(g(e.s), g(e.a), g(e.b))
    if isinstance(e, Merge):
        return Merge(g(e.a), g(e.b), g(e.s))
    raise TypeError(e)


@dataclass
class _Config:
    """Where the sensor is and what it sees: the gaps with their (global) labels."""

    point: Point
    gaps: List[PhysicalGap]
    labels: List[int]
    critical_index: int = -1  # the critical point whose sample this is (-1: the robot position)
    _shadows: Optional[List[Shadow]] = field(default=None, repr=False)
    _arrays: Optional[list] = field(default=None, repr=False)

    def shadows(self, poly: Polygon) -> List[Shadow]:
        if self._shadows is None:
            self._shadows = shadow_pockets(poly, self.point, self.gaps, self.labels)
        return self._shadows

    def pocket_arrays(self, poly: Polygon):
        if self._arrays is None:
            arr = []
            for sh in self.shadows(poly):
                a = np.asarray(sh.pocket, dtype=np.float64)
                arr.append((sh.label, a[:, 0].copy(), a[:, 1].copy(), a[:, 0].min(), a[:, 0].max(),
                            a[:, 1].min(), a[:, 1].max()))
            self._arrays = arr
        return self._arrays


def _structure(gaps: Sequence[PhysicalGap]) -> List[Tuple[int, int, bool]]:
    return [(g.start_edge, g.end_edge, g.start_point is None) for g in gaps]


def _bijection(poly: Polygon, a: Sequence[PhysicalGap], b: Sequence[PhysicalGap]) -> Optional[List[int]]:
    """``perm[k]`` = index in ``a`` of the gap ``b[k]`` continues, if the hidden chains pair them 1-1."""
    if len(a) != len(b):
        return None
    comps = _components(len(a), len(b), chain_graph(poly, a, b))
    perm = [-1] * len(b)
    for ja, kb in comps:
        if len(ja) != 1 or len(kb) != 1:
            return None
        perm[kb[0]] = ja[0]
    return perm


def _dist_to_polyline(xs: np.ndarray, ys: np.ndarray, px: float, py: float) -> float:
    x1, y1 = xs, ys
    x2, y2 = np.roll(xs, -1), np.roll(ys, -1)
    dx, dy = x2 - x1, y2 - y1
    L = dx * dx + dy * dy
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.clip(np.where(L > 0, ((px - x1) * dx + (py - y1) * dy) / L, 0.0), 0.0, 1.0)
    cx, cy = x1 + t * dx, y1 + t * dy
    return float(np.min((cx - px) ** 2 + (cy - py) ** 2))


def _fov_events(exits: Dict[int, int], enters: Dict[int, int]) -> List[Event]:
    return [Exit(s, k) for s, k in sorted(exits.items())] + [Enter(s, k) for s, k in sorted(enters.items())]


def _inc(d: Dict[int, int], k: int) -> None:
    d[k] = d.get(k, 0) + 1


# --------------------------------------------------------------------------- the simulator


class PolygonSimulator:
    """Deterministic polygon-world simulator (robot on a path, random-walk targets).

    ``poly``: a :class:`~shadowinfo.polygon.geometry.Polygon` or a map number 1..14.  ``path``: a
    :class:`~shadowinfo.polygon.geometry.Path`, a waypoint list, a stored path label
    (:func:`~shadowinfo.polygon.io.path_record`) or ``None`` for the map's default demo path
    (:func:`~shadowinfo.polygon.io.load_demo_path`).  Waypoint paths are validated
    (:func:`validate_path`) unless ``validate=False``.  ``speed`` is the robot's arc length per tick,
    ``target_step`` the half-width of a target step.  ``at_end`` in :data:`AT_END`, or ``"auto"``:
    ``"loop"`` for a closed path, the demo record's choice for a demo path, ``"stop"`` otherwise.

    ``frame()`` describes the current state as a JSON-ready dict; ``step()`` advances one tick and
    returns the frame with that tick's ``events`` and ``created`` (``[[label, n], ...]``: the true count
    of every shadow created by a split or merge this tick, at the moment it was created).

    ``stats`` counts the rare numerical fallbacks: ``"sense_fallbacks"`` (the gaps at the robot
    position do not continue those of the last critical point crossed -- the robot sits on a critical
    line, or ``get_gaps``' sample overshot the next crossing -- so the sensor keeps the last sample
    point for that tick, reported as ``frame()["sense"]``) and ``"membership_fallbacks"`` (a hidden
    target lies in no pocket by a rounding error and is assigned to the nearest one).
    """

    def __init__(self, poly: Union[Polygon, int], path: Union[Path, Sequence[Sequence[float]], str, None] = None,
                 n_targets: int = 10, seed: int = 1, speed: float = 4.0, *, target_step: float = 8.0,
                 at_end: str = "auto", validate: bool = True) -> None:
        from .io import load_path, load_polygon

        if not isinstance(poly, Polygon):
            poly = load_polygon(int(poly))
        self.poly: Polygon = poly
        demo_at_end = None
        if path is None:
            if poly.name is None or not poly.name.isdigit():
                raise ValueError("path=None needs a numbered map (1..14)")
            path, demo_at_end = load_demo_path(int(poly.name))
        elif isinstance(path, str):
            label = path
            p2, path = load_path(label)
            if p2 is not poly and p2.vertices != poly.vertices:
                raise ValueError(f"path {label!r} belongs to polygon {p2.name}, not {poly.name}")
        if not isinstance(path, Path):
            path = Path(path)
        if validate:
            validate_path(poly, path.points)
        if speed <= 0 or not math.isfinite(speed):
            raise ValueError("speed must be a positive number")
        if at_end == "auto":
            at_end = demo_at_end or ("loop" if path.points[0] == path.points[-1] else "stop")
        self._check_at_end(at_end, path)
        self.path: Path = path
        self.at_end = at_end
        self.speed = float(speed)
        self.target_step = float(target_step)
        self.rng = Rng(seed)
        self.t = 0
        self.stats: Dict[str, int] = {"sense_fallbacks": 0, "membership_fallbacks": 0, "crossings": 0, "legs": 1}
        self.history: List[Event] = []
        self._pending: List[Tuple[List[Event], _Config]] = []  # steps queued by follow() for the next tick
        self._next_label = 1
        self.legs: List[_Leg] = []
        # first leg: labels are get_gaps' own
        tr = _track(poly, path)
        leg = _Leg(tr, {k: k for k in range(1, tr.m + 1)}, tr.m + 1)
        self._next_label = tr.max_local + 1
        self.legs.append(leg)
        first = _Config(tr.history.samples.points[0], list(tr.history.samples.physical[0]), leg.labels(0), 0)  # type: ignore[union-attr,arg-type]
        self.robot: Point = path.points[0]
        self._cfg = self._config_at(self.robot, first)
        self.targets: List[Point] = self._place_targets(int(n_targets))
        self._states = self._target_states(self._cfg, self.targets)
        self.initial_counts = self.counts()
        self.initial_labels = list(self._cfg.labels)

    # ------------------------------------------------------------------ public state

    @property
    def leg(self) -> _Leg:
        return self.legs[-1]

    @property
    def labels(self) -> List[int]:
        """Labels of the shadows alive now, in boundary order (``physical_gaps`` order)."""
        return list(self._cfg.labels)

    @property
    def sense_point(self) -> Point:
        """Where the sensor is: the robot position, except for a ``sense_fallbacks`` tick."""
        return self._cfg.point

    @property
    def target_labels(self) -> List[int]:
        """Shadow of each target (0 = visible)."""
        return list(self._states)

    def shadows(self) -> List[Shadow]:
        """The shadows now, with their pockets."""
        return list(self._cfg.shadows(self.poly))

    def counts(self) -> Dict[int, int]:
        """True number of targets in each alive shadow."""
        c = {s: 0 for s in self._cfg.labels}
        for s in self._states:
            if s:
                c[s] += 1
        return c

    def initial_condition(self, mode: str = "exact") -> Dict[int, Bound]:
        return initial_bounds(self.initial_counts, mode)

    def sequence(self, mode: str = "exact") -> ShadowSequence:
        """All events so far, with the initial condition for ``mode`` (``exact``/``unknown``/``evader``)."""
        return ShadowSequence(self.initial_condition(mode), list(self.history))

    def visibility_polygon(self) -> List[Point]:
        """The visibility polygon of the sensing point (``visibility_polygon(compat="fixed")``)."""
        return visibility_polygon(self.poly, self.sense_point)

    def frame(self, geometry: bool = False) -> dict:
        """JSON-ready state.  ``geometry=True`` adds the pockets and the visibility polygon."""
        f = {
            "t": self.t,
            "robot": list(self.robot),
            "sense": list(self._cfg.point),
            "leg": len(self.legs) - 1,
            "s": self.leg.s,
            "targets": [list(p) for p in self.targets],
            "target_labels": list(self._states),
            "counts": [[s, n] for s, n in sorted(self.counts().items())],
            "labels": list(self._cfg.labels),
        }
        if geometry:
            f["shadows"] = [sh.to_dict() for sh in self.shadows()]
            f["visibility"] = [list(p) for p in self.visibility_polygon()]
        return f

    # ------------------------------------------------------------------ routing

    @staticmethod
    def _check_at_end(at_end: str, path: Path) -> None:
        if at_end not in AT_END:
            raise ValueError(f"at_end must be one of {AT_END} or 'auto', got {at_end!r}")
        if at_end == "loop" and path.points[0] != path.points[-1]:
            raise ValueError("at_end='loop' needs a closed path (first waypoint == last waypoint)")

    def follow(self, waypoints: Sequence[Sequence[float]], at_end: str = "stop", validate: bool = True) -> None:
        """Leave the current path: from the robot position, follow ``waypoints`` (validated)."""
        pts = [self.robot] + [(float(p[0]), float(p[1])) for p in waypoints]
        pts = [p for i, p in enumerate(pts) if i == 0 or p != pts[i - 1]]
        if len(pts) < 2:
            return
        if validate:
            validate_path(self.poly, pts)
        path = Path(pts)
        self._check_at_end(at_end, path)
        self.path = path
        self.at_end = at_end
        self._start_leg(path, self._pending)

    def go_to(self, point: Sequence[float]) -> None:
        """Click-to-move: drive straight to ``point`` (the segment must lie in the polygon) and stop."""
        self.follow([point], "stop")

    def _start_leg(self, path: Path, steps: List[Tuple[List[Event], _Config]]) -> None:
        """Begin tracking ``path`` (starting at the robot position) and join its shadows to the current ones."""
        tr = _track(self.poly, path)
        cur = steps[-1][1] if steps else self._cfg
        new0 = list(tr.history.samples.physical[0])  # type: ignore[union-attr,arg-type]
        point0 = tr.history.samples.points[0]  # type: ignore[union-attr]
        if _structure(cur.gaps) == _structure(new0):
            ids, links, evs = list(cur.labels), {}, []
        else:
            perm = _bijection(self.poly, cur.gaps, new0)
            if perm is not None:
                ids, links, evs = [cur.labels[j] for j in perm], {}, []
            else:
                prev = [Gap(g, lab) for g, lab in zip(cur.gaps, cur.labels)]
                ids, links, last = _from_graph(prev, len(new0), chain_graph(self.poly, cur.gaps, new0),
                                               self._next_label - 1)
                self._next_label = last + 1
                evs = self._link_events(cur.labels, ids, links)
        base = self._next_label
        leg = _Leg(tr, {k + 1: ids[k] for k in range(tr.m)}, base)
        self._next_label = base + tr.max_local - tr.m
        self.legs.append(leg)
        self.stats["legs"] += 1
        steps.append((evs, _Config(point0, new0, ids, 0)))  # type: ignore[arg-type]

    def _link_events(self, prev_labels: List[int], ids: List[int], links: Dict[int, List[int]]) -> List[Event]:
        prev = [Gap(None, lab) for lab in prev_labels]
        for j, cs in links.items():
            for c in cs:
                prev[j].add_gap(c)
        nxt = [Gap(None, x) for x in ids]
        known = set(prev_labels) | set(ids)
        fresh: Dict[int, int] = {}

        def g(x: int) -> int:
            if x in known:
                return x
            if x not in fresh:
                fresh[x] = self._next_label
                self._next_label += 1
            return fresh[x]

        return [_relabel(e, g) for e in transition_events([prev, nxt])[0]]

    def _next_path(self) -> Optional[Path]:
        if self.at_end == "stop":
            return None
        if self.at_end == "loop":
            return self.path
        return Path(list(reversed(self.leg.track.path.points)))

    # ------------------------------------------------------------------ sensing

    def _config_at(self, q: Point, last: _Config) -> _Config:
        """The gaps at ``q`` labelled as continuations of ``last`` (no critical point in between)."""
        if q == last.point:
            return last
        gaps = physical_gaps(self.poly, q)
        if _structure(gaps) == _structure(last.gaps):
            return _Config(q, gaps, list(last.labels))
        perm = _bijection(self.poly, last.gaps, gaps)
        if perm is not None:
            return _Config(q, gaps, [last.labels[j] for j in perm])
        self.stats["sense_fallbacks"] += 1
        return last

    def _target_states(self, cfg: _Config, pts: Sequence[Point], idx: Optional[Sequence[int]] = None) -> List[int]:
        """Shadow label of each point (0 = visible from ``cfg.point``)."""
        if not pts:
            return []
        X = np.array([p[0] for p in pts], dtype=np.float64)
        Y = np.array([p[1] for p in pts], dtype=np.float64)
        vis = segments_in_polygon(self.poly, cfg.point, X, Y)
        out = [0] * len(pts)
        hidden = np.flatnonzero(~vis)
        if hidden.size == 0:
            return out
        hx, hy = X[hidden], Y[hidden]
        found = np.zeros(hidden.size, dtype=np.int64)
        nfound = np.zeros(hidden.size, dtype=np.int64)
        for label, xs, ys, x0, x1, y0, y1 in cfg.pocket_arrays(self.poly):
            box = (hx >= x0) & (hx <= x1) & (hy >= y0) & (hy <= y1)
            if not box.any():
                continue
            sel = np.flatnonzero(box)
            inside = _pip_arrays(xs, ys, hx[sel], hy[sel])
            hit = sel[inside]
            nfound[hit] += 1
            found[hit] = np.where(found[hit] == 0, label, found[hit])
        for r in range(hidden.size):
            if nfound[r] != 1:
                self.stats["membership_fallbacks"] += 1
                best, lab = math.inf, 0
                for label, xs, ys, *_ in cfg.pocket_arrays(self.poly):
                    d = _dist_to_polyline(xs, ys, float(hx[r]), float(hy[r]))
                    if d < best:
                        best, lab = d, label
                found[r] = lab
            out[int(hidden[r])] = int(found[r])
        return out

    # ------------------------------------------------------------------ targets

    def _place_targets(self, n: int) -> List[Point]:
        x0, x1 = float(self.poly.vx.min()), float(self.poly.vx.max())
        y0, y1 = float(self.poly.vy.min()), float(self.poly.vy.max())
        out: List[Point] = []
        while len(out) < n:
            p = (x0 + self.rng.random() * (x1 - x0), y0 + self.rng.random() * (y1 - y0))
            if point_in_polygon(self.poly, p) and not _on_boundary(self.poly, p):
                out.append(p)
        return out

    def _move_targets(self) -> List[int]:
        """One random-walk step for every target; returns the indices of the targets that moved."""
        r = self.target_step
        cand = []
        for (x, y) in self.targets:
            dx = (2.0 * self.rng.random() - 1.0) * r
            dy = (2.0 * self.rng.random() - 1.0) * r
            cand.append((x + dx, y + dy))
        if not cand:
            return []
        cx = np.array([c[0] for c in cand])
        cy = np.array([c[1] for c in cand])
        ox = np.array([p[0] for p in self.targets])
        oy = np.array([p[1] for p in self.targets])
        ok = _pip_arrays(self.poly.vx, self.poly.vy, cx, cy) & ~_touch_many(self.poly, ox, oy, cx, cy)
        moved = np.flatnonzero(ok).tolist()
        for j in moved:
            self.targets[j] = cand[j]
        return moved

    # ------------------------------------------------------------------ one tick

    def _advance(self, steps: List[Tuple[List[Event], _Config]]) -> None:
        remaining = self.speed
        while remaining > 0:
            leg = self.leg
            L = leg.length
            if leg.s >= L:
                nxt = self._next_path()
                if nxt is None:
                    break
                self._start_leg(nxt, steps)
                continue
            s_new = leg.s + remaining
            if s_new >= L:
                remaining = s_new - L
                s_new = L
            else:
                remaining = 0.0
            dist = leg.track.distances
            samples = leg.track.history.samples
            while leg.next_cp < len(dist) and (dist[leg.next_cp] <= s_new or s_new >= L):
                i = leg.next_cp
                steps.append((leg.events(i), _Config(samples.points[i], list(samples.physical[i]),  # type: ignore[union-attr,arg-type]
                                                     leg.labels(i), i)))
                self.stats["crossings"] += 1
                leg.next_cp += 1
            leg.s = s_new
        self.robot = self.leg.point_at(self.leg.s)

    def _sensor_events(self, old: List[int], new: List[int], evs: List[Event]) -> Tuple[List[Event], Dict[int, int]]:
        """FOV + component events of one step from per-target labels before/after (see module docstring)."""
        appeared = {e.s for e in evs if isinstance(e, Appear)}
        disappeared = {e.s for e in evs if isinstance(e, Disappear)}
        children: Dict[int, List[int]] = {}
        for e in evs:
            if isinstance(e, Split):
                children.setdefault(e.s, []).extend([e.a, e.b])
            elif isinstance(e, Merge):
                children.setdefault(e.a, []).append(e.s)
                children.setdefault(e.b, []).append(e.s)
        desc_cache: Dict[int, set] = {}

        def desc(a: int) -> set:
            if a not in desc_cache:
                out, todo = set(), [a]
                while todo:
                    for c in children.get(todo.pop(), ()):
                        if c not in out:
                            out.add(c)
                            todo.append(c)
                desc_cache[a] = out
            return desc_cache[a]

        exits: Dict[int, int] = {}
        enters: Dict[int, int] = {}
        revealed: Dict[int, int] = {}
        covered: Dict[int, int] = {}
        stay_from: Dict[int, int] = {}
        stay_to: Dict[int, int] = {}
        for a, b in zip(old, new):
            if a and b and (a == b or b in desc(a)):
                if a != b:
                    _inc(stay_from, a)
                    _inc(stay_to, b)
                continue
            if a:
                _inc(revealed if a in disappeared else exits, a)
            if b:
                _inc(covered if b in appeared else enters, b)
        made = set()
        for e in evs:
            if isinstance(e, Split):
                made.update((e.a, e.b))
            elif isinstance(e, Merge):
                made.add(e.s)
        val = {x: stay_to.get(x, 0) for x in made}
        for e in reversed(evs):  # split chains: an intermediate label holds what its children receive
            if isinstance(e, Split) and e.s in made:
                val[e.s] = val[e.a] + val[e.b]
        cnt = dict(stay_from)
        created: Dict[int, int] = {}
        comp: List[Event] = []
        for e in evs:
            if isinstance(e, Appear):
                k = covered.get(e.s, 0)
                e = Appear(e.s, k, k)
            elif isinstance(e, Disappear):
                k = revealed.get(e.s, 0)
                e = Disappear(e.s, k, k)
            elif isinstance(e, Merge):
                cnt[e.s] = created[e.s] = cnt.pop(e.a, 0) + cnt.pop(e.b, 0)
            elif isinstance(e, Split):
                n = cnt.pop(e.s, 0)
                cnt[e.a] = created[e.a] = val[e.a]
                cnt[e.b] = created[e.b] = val[e.b]
                if n != val[e.a] + val[e.b]:  # pragma: no cover - would be a bookkeeping bug
                    raise AssertionError(f"split {e}: {n} targets but children receive {val[e.a]} + {val[e.b]}")
            comp.append(e)
        return _fov_events(exits, {}) + comp + _fov_events({}, enters), created

    def step(self) -> dict:
        """Advance one tick (sensor substep, then target substep) and return the frame."""
        self.t += 1
        steps: List[Tuple[List[Event], _Config]] = list(self._pending)
        self._pending.clear()
        self._advance(steps)
        last = steps[-1][1] if steps else self._cfg
        final = self._config_at(self.robot, last)
        if final is not last or not steps:
            steps.append(([], final))
        events: List[Event] = []
        created: Dict[int, int] = {}
        states = self._states
        for evs, cfg in steps:
            if cfg is self._cfg and not evs:
                continue
            new_states = self._target_states(cfg, self.targets)
            ev, cr = self._sensor_events(states, new_states, evs)
            events += ev
            created.update(cr)
            states = new_states
        self._cfg = final
        # target substep
        moved = self._move_targets()
        new_states = list(states)
        if moved:
            ms = self._target_states(final, [self.targets[j] for j in moved])
            for j, s in zip(moved, ms):
                new_states[j] = s
        exits: Dict[int, int] = {}
        enters: Dict[int, int] = {}
        for a, b in zip(states, new_states):
            if a != b:
                if a:
                    _inc(exits, a)
                if b:
                    _inc(enters, b)
        events += _fov_events(exits, enters)
        self._states = new_states
        self.history.extend(events)
        f = self.frame()
        f["events"] = [event_to_dict(e) for e in events]
        f["created"] = [[s, n] for s, n in sorted(created.items())]
        return f

    def run(self, ticks: int) -> List[dict]:
        """Frame 0 followed by ``ticks`` stepped frames."""
        frames = [self.frame()]
        for _ in range(ticks):
            frames.append(self.step())
        return frames

    def ticks_per_lap(self) -> int:
        """Ticks needed to run the current path once (``ceil(length / speed)``)."""
        return int(math.ceil(self.path.length() / self.speed))


def simulate(poly: Union[Polygon, int], path=None, n_targets: int = 10, seed: int = 1, ticks: Optional[int] = None,
             **kw) -> Tuple[PolygonSimulator, List[dict]]:
    """Build a :class:`PolygonSimulator` and run it ``ticks`` ticks (default: one pass of the path)."""
    sim = PolygonSimulator(poly, path, n_targets, seed, **kw)
    return sim, sim.run(sim.ticks_per_lap() if ticks is None else ticks)
