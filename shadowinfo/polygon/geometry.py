"""Geometry primitives of the original Java implementation, reproduced bit for bit.

Port of ``geometry.Algorithm`` (primitives), ``drawable.Polygon`` and ``drawable.Path`` together with
the ``java.awt.geom`` library semantics they depend on (``docs/notes/original_java.md`` §2.1-§2.3 and
§2.6).  Every formula copies the Java operation order, so that cuts, critical points and visibility
come out as the *same doubles* as in Java (quirk B14):

* ``Line2D.relativeCCW``, ``linesIntersect``, ``ptSegDist(Sq)`` and ``Point2D.distance`` are the JDK
  formulas (``sqrt(dx*dx + dy*dy)``, never ``hypot``);
* ``Polygon.pointInPolygon`` is ``Path2D.Float.contains``: vertices rounded to float32, non-zero
  winding, ``Curve.pointCrossingsForLine`` (half-open in y);
* the ``Algorithm`` helpers (``getIntersect`` in slope-intercept form, ``onExtension`` as a per-axis
  OR, absolute ``EPSILON``) keep their quirks.

Coordinates are the y-up ``.dat`` frame.  A point is a tuple ``(x, y)`` and a line/segment a tuple
``(x1, y1, x2, y2)`` (``Line2D.Double``).  Hot loops are vectorised with numpy.  IEEE-754 ``+ - * /``
and ``sqrt`` are correctly rounded both in Java and in numpy, so vectorising over independent
elements does not change any bit; only integer sums (crossing counts) and boolean ``any`` are
reduced.  The scalar functions are the reference; ``tests/test_polygon_geometry.py`` checks that the
vectorised paths agree with them.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .io import EPSILON, PERTURB

if TYPE_CHECKING:  # pragma: no cover
    from .cuts import CriticalPoint, Cut

Point = Tuple[float, float]
Line = Tuple[float, float, float, float]

COMPAT_MODES = ("java", "fixed")

__all__ = [
    "Point", "Line", "COMPAT_MODES", "GeometryError", "epsilon_equal", "relative_ccw", "relative_ccw_xy",
    "lines_intersect", "pt_seg_dist_sq", "pt_seg_dist", "distance", "distance_sq", "get_slope",
    "get_y_intercept", "get_x_intercept", "get_point_on_line_from_x", "get_point_on_line_from_y",
    "get_intersect", "get_segment_intersect", "segments_on_same_line", "on_extension", "on_reverse_extension",
    "point_on_segment", "purturb_point_along_seg", "point_crossings_for_line", "GeneralPath", "Polygon",
    "seg_in_polygon", "seg_in_polygon_scalar", "seg_in_polygon_many", "Path",
]


class GeometryError(ValueError):
    """A geometric computation that throws in Java (``NullPointerException``) has no answer.

    Raised e.g. by :func:`purturb_point_along_seg` at slope exactly +-1 in ``compat="java"`` (B1), or by
    the cut functions when a ray never hits the boundary (cannot happen in a closed simple polygon).
    """


def _check_compat(compat: str) -> None:
    if compat not in COMPAT_MODES:
        raise ValueError(f"compat must be one of {COMPAT_MODES}, got {compat!r}")


# --------------------------------------------------------------------------- scalar JDK primitives


def epsilon_equal(n1: float, n2: float, eps: float = EPSILON) -> bool:
    """``Algorithm.epsilonEqual``: ``|n1 - n2| <= eps``."""
    return abs(n1 - n2) <= eps


def relative_ccw_xy(x1: float, y1: float, x2: float, y2: float, px: float, py: float) -> int:
    """``Line2D.relativeCCW(x1, y1, x2, y2, px, py)`` (JDK formula, operation for operation).

    In the y-up frame it returns ``+1`` when ``P`` is to the *right* of ``A->B`` and ``-1`` when it is
    to the left.  For collinear points: ``-1`` before ``A``, ``+1`` beyond ``B``, ``0`` on the segment.
    """
    x2 -= x1
    y2 -= y1
    px -= x1
    py -= y1
    ccw = px * y2 - py * x2
    if ccw == 0.0:
        ccw = px * x2 + py * y2
        if ccw > 0.0:
            px -= x2
            py -= y2
            ccw = px * x2 + py * y2
            if ccw < 0.0:
                ccw = 0.0
    return -1 if ccw < 0.0 else (1 if ccw > 0.0 else 0)


def relative_ccw(line: Line, p: Point) -> int:
    """``line.relativeCCW(p)``.  Reflex test of the Java code: ``relative_ccw(e_{i-1}, v_{i+1}) == 1``."""
    return relative_ccw_xy(line[0], line[1], line[2], line[3], p[0], p[1])


def lines_intersect(l1: Line, l2: Line) -> bool:
    """``Line2D.intersectsLine`` / ``linesIntersect``: closed segment test (touching counts)."""
    x1, y1, x2, y2 = l1
    x3, y3, x4, y4 = l2
    return (relative_ccw_xy(x1, y1, x2, y2, x3, y3) * relative_ccw_xy(x1, y1, x2, y2, x4, y4) <= 0
            and relative_ccw_xy(x3, y3, x4, y4, x1, y1) * relative_ccw_xy(x3, y3, x4, y4, x2, y2) <= 0)


def pt_seg_dist_sq(line: Line, p: Point) -> float:
    """``Line2D.ptSegDistSq`` (JDK projection formula)."""
    x1, y1, x2, y2 = line
    px, py = p
    x2 -= x1
    y2 -= y1
    px -= x1
    py -= y1
    dotprod = px * x2 + py * y2
    if dotprod <= 0.0:
        projlen_sq = 0.0
    else:
        px = x2 - px
        py = y2 - py
        dotprod = px * x2 + py * y2
        if dotprod <= 0.0:
            projlen_sq = 0.0
        else:
            projlen_sq = dotprod * dotprod / (x2 * x2 + y2 * y2)
    len_sq = px * px + py * py - projlen_sq
    if len_sq < 0:
        len_sq = 0.0
    return len_sq


def pt_seg_dist(line: Line, p: Point) -> float:
    """``Line2D.ptSegDist`` = ``sqrt(ptSegDistSq)``."""
    return math.sqrt(pt_seg_dist_sq(line, p))


def distance(p: Point, q: Point) -> float:
    """``p.distance(q)`` = ``Math.sqrt(dx*dx + dy*dy)`` with ``dx = q.x - p.x`` (not ``hypot``)."""
    dx = q[0] - p[0]
    dy = q[1] - p[1]
    return math.sqrt(dx * dx + dy * dy)


def distance_sq(p: Point, q: Point) -> float:
    """``p.distanceSq(q)``."""
    dx = q[0] - p[0]
    dy = q[1] - p[1]
    return dx * dx + dy * dy


def get_slope(line: Line, eps: float = EPSILON) -> float:
    """``Algorithm.getSlope``: ``+inf`` when ``|x1 - x2| <= eps``, else ``dy/dx``."""
    if epsilon_equal(line[0], line[2], eps):
        return math.inf
    return (line[3] - line[1]) / (line[2] - line[0])


def get_y_intercept(line_or_k, p: Optional[Point] = None, eps: float = EPSILON) -> float:
    """``Algorithm.getYIntercept(line)`` or, with ``p``, ``getYIntercept(k, p)``.

    With an infinite slope the "intercept" is the slope itself (``+inf``).
    """
    if p is None:
        line = line_or_k
        k = get_slope(line, eps)
        return k if math.isinf(k) else line[1] - line[0] * k
    k = line_or_k
    return k if math.isinf(k) else p[1] - p[0] * k


def get_x_intercept(line: Line, eps: float = EPSILON) -> float:
    """``Algorithm.getXIntercept``.  Bug kept: it is identical to :func:`get_y_intercept` (unused in Java)."""
    return get_y_intercept(line, eps=eps)


def get_point_on_line_from_x(line: Line, x: float, eps: float = EPSILON) -> Optional[Point]:
    """``Algorithm.getPointOnLineFromX``: ``(x, k*x + b)``.

    Quirk kept: for an eps-vertical line it returns the meaningless ``(x, 0)`` when ``x != x1`` and
    ``None`` when ``x == x1``.
    """
    if epsilon_equal(line[0], line[2], eps):
        return (x, 0.0) if x != line[0] else None
    k = get_slope(line, eps)
    intercept = get_y_intercept(line, eps=eps)
    return (x, k * x + intercept)


def get_point_on_line_from_y(line: Line, y: float, eps: float = EPSILON) -> Optional[Point]:
    """``Algorithm.getPointOnLineFromY``: ``((y - b)/k, y)``; ``(0, y)`` or ``None`` for eps-horizontal lines."""
    if epsilon_equal(line[1], line[3], eps):
        return (0.0, y) if y != line[1] else None
    k = get_slope(line, eps)
    intercept = get_y_intercept(line, eps=eps)
    return ((y - intercept) / k, y)


def get_intersect(line1: Line, line2: Line, eps: float = EPSILON) -> Optional[Point]:
    """``Algorithm.getIntersect``: intersection of the two *infinite* lines, slope-intercept form.

    If exactly one line is eps-vertical the other is evaluated at its x; both vertical, or slopes
    eps-equal, give ``None``.  Quirk kept: the second intercept uses ``line2``'s *second* point.
    """
    x11, y11, x12, y12 = line1
    x21, y21, x22, y22 = line2
    v1 = epsilon_equal(x12, x11, eps)
    v2 = epsilon_equal(x22, x21, eps)
    if v1 or v2:
        if not v1:
            return get_point_on_line_from_x(line1, x21, eps)
        if not v2:
            return get_point_on_line_from_x(line2, x11, eps)
        return None
    k1 = get_slope(line1, eps)
    k2 = get_slope(line2, eps)
    i1 = get_y_intercept(k1, (x11, y11))
    i2 = get_y_intercept(k2, (x22, y22))
    if epsilon_equal(k1, k2, eps):
        return None
    x = (i2 - i1) / (k1 - k2)
    y = x * k1 + i1
    return (x, y)


def get_segment_intersect(line1: Line, line2: Line, eps: float = EPSILON) -> Optional[Point]:
    """``Algorithm.getSegmentIntersect``: :func:`get_intersect` within ``eps`` (``ptSegDist``) of both segments."""
    p = get_intersect(line1, line2, eps)
    if p is not None and pt_seg_dist(line1, p) < eps and pt_seg_dist(line2, p) < eps:
        return p
    return None


def segments_on_same_line(line1: Line, line2: Line, eps: float = EPSILON) -> bool:
    """``Algorithm.segmentsOnSameLine``.  Bug kept: the second intercept uses ``line1.P2`` with ``k2``
    (unused in Java)."""
    if (epsilon_equal(line1[2], line1[0], eps) and epsilon_equal(line2[2], line2[0], eps)
            and epsilon_equal(line1[2], line2[0], eps)):
        return True
    k1 = get_slope(line1, eps)
    k2 = get_slope(line2, eps)
    i1 = get_y_intercept(k1, (line1[0], line1[1]))
    i2 = get_y_intercept(k2, (line1[2], line1[3]))
    return epsilon_equal(k1, k2, eps) and epsilon_equal(i1, i2, eps)


def on_extension(line: Line, p: Point) -> bool:
    """``Algorithm.onExtension``: ``p`` lies past ``P2`` -- a per-axis OR, not a projection."""
    x1, y1, x2, y2 = line
    return (p[0] - x2) * (x2 - x1) > 0 or (p[1] - y2) * (y2 - y1) > 0


def on_reverse_extension(line: Line, p: Point) -> bool:
    """``Algorithm.onReverseExtension``: :func:`on_extension` of the reversed line (``p`` before ``P1``)."""
    x1, y1, x2, y2 = line
    return (p[0] - x1) * (x1 - x2) > 0 or (p[1] - y1) * (y1 - y2) > 0


def point_on_segment(p: Point, line: Line, eps: float = EPSILON) -> bool:
    """``Algorithm.pointOnSegment``: ``ptSegDist < eps``."""
    return pt_seg_dist(line, p) < eps


def purturb_point_along_seg(p: Point, seg: Line, amount: float = PERTURB, far: bool = True,
                            compat: str = "fixed", eps: float = EPSILON) -> Point:
    """``Algorithm.purturbPointAlongSeg``: move ``p`` along ``seg`` by ``amount`` in x or in y.

    The step is ``amount`` in x when ``|slope| < 1`` (or horizontal) and in y when ``|slope| > 1`` (or
    vertical), so its length is between ``amount`` and ``amount*sqrt(2)``.  Of the two candidates it
    returns the one farther from (``far``) or nearer to ``seg``'s first point.

    Quirk B1: at ``|slope| == 1`` exactly Java dereferences ``null`` and throws; ``compat="java"``
    raises :class:`GeometryError` there, ``compat="fixed"`` (default) uses the x step.
    """
    _check_compat(compat)
    x, y = p
    s = get_slope(seg, eps)
    if math.isinf(s):
        p1, p2 = (x, y + amount), (x, y - amount)
    elif s == 0:
        p1, p2 = (x + amount, y), (x - amount, y)
    elif abs(s) > 1:
        p1 = get_point_on_line_from_y(seg, y + amount, eps)
        p2 = get_point_on_line_from_y(seg, y - amount, eps)
    elif abs(s) < 1 or compat == "fixed":
        p1 = get_point_on_line_from_x(seg, x + amount, eps)
        p2 = get_point_on_line_from_x(seg, x - amount, eps)
    else:
        raise GeometryError(f"purturbPointAlongSeg: slope {s!r} of {seg} is exactly +-1 (Java NullPointerException)")
    if p1 is None or p2 is None:
        raise GeometryError(f"purturbPointAlongSeg: no candidate point on {seg} (Java NullPointerException)")
    a = (seg[0], seg[1])
    dist1 = distance(p1, a)
    dist2 = distance(p2, a)
    if (dist1 > dist2 and far) or (dist1 < dist2 and not far):
        return p1
    return p2


def point_crossings_for_line(px: float, py: float, x0: float, y0: float, x1: float, y1: float) -> int:
    """``sun.awt.geom.Curve.pointCrossingsForLine``: crossings of the ray to the right of ``(px, py)``."""
    if py < y0 and py < y1:
        return 0
    if py >= y0 and py >= y1:
        return 0
    if px >= x0 and px >= x1:
        return 0
    if px < x0 and px < x1:
        return 1 if y0 < y1 else -1
    xintercept = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
    if px >= xintercept:
        return 0
    return 1 if y0 < y1 else -1


# --------------------------------------------------------------------------- numpy versions
# Each mirrors the scalar function above operation for operation, elementwise with broadcasting.


def _np_rccw(x1, y1, x2, y2, px, py):
    """Vectorised :func:`relative_ccw_xy` (int8 result)."""
    x2 = x2 - x1
    y2 = y2 - y1
    px = px - x1
    py = py - y1
    ccw = px * y2 - py * x2
    c2 = px * x2 + py * y2
    c3 = (px - x2) * x2 + (py - y2) * y2
    c3 = np.where(c3 < 0.0, 0.0, c3)
    ccw = np.where(ccw == 0.0, np.where(c2 > 0.0, c3, c2), ccw)
    return np.where(ccw < 0.0, -1, np.where(ccw > 0.0, 1, 0)).astype(np.int8)


def _np_pt_seg_dist(x1, y1, x2, y2, px, py):
    """Vectorised :func:`pt_seg_dist`."""
    x1, y1, x2, y2, px, py = (np.asarray(a, dtype=np.float64) for a in (x1, y1, x2, y2, px, py))
    x2 = x2 - x1
    y2 = y2 - y1
    px = px - x1
    py = py - y1
    dot1 = px * x2 + py * y2
    bx = x2 - px
    by = y2 - py
    dot2 = bx * x2 + by * y2
    first = dot1 <= 0.0
    mid = ~first & ~(dot2 <= 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        proj = np.where(mid, dot2 * dot2 / (x2 * x2 + y2 * y2), 0.0)
    fx = np.where(first, px, bx)
    fy = np.where(first, py, by)
    len_sq = fx * fx + fy * fy - proj
    len_sq = np.where(len_sq < 0, 0.0, len_sq)
    return np.sqrt(len_sq)


def _np_distance(ax, ay, bx, by):
    dx = bx - ax
    dy = by - ay
    return np.sqrt(dx * dx + dy * dy)


def _np_get_intersect(ax1, ay1, ax2, ay2, bx1, by1, bx2, by2, eps: float = EPSILON):
    """Vectorised :func:`get_intersect` of line ``a`` (line1) and line ``b`` (line2): ``(x, y, valid)``."""
    ax1, ay1, ax2, ay2, bx1, by1, bx2, by2 = (np.asarray(a, dtype=np.float64)
                                              for a in (ax1, ay1, ax2, ay2, bx1, by1, bx2, by2))
    with np.errstate(all="ignore"):
        v1 = np.abs(ax2 - ax1) <= eps
        v2 = np.abs(bx2 - bx1) <= eps
        k1 = (ay2 - ay1) / (ax2 - ax1)
        k2 = (by2 - by1) / (bx2 - bx1)
        i1 = ay1 - ax1 * k1  # getYIntercept(k1, line1.P1) == getYIntercept(line1)
        i2b = by1 - bx1 * k2  # getYIntercept(line2) (uses P1)
        y_a = k1 * bx1 + i1  # getPointOnLineFromX(line1, line2.x1)
        y_b = k2 * ax1 + i2b  # getPointOnLineFromX(line2, line1.x1)
        i2 = by2 - bx2 * k2  # getYIntercept(k2, line2.P2)
        xg = (i2 - i1) / (k1 - k2)
        yg = xg * k1 + i1
        anyv = v1 | v2
        x = np.where(anyv, np.where(v1, ax1, bx1), xg)
        y = np.where(anyv, np.where(v1, y_b, y_a), yg)
        valid = np.where(anyv, ~(v1 & v2), ~(np.abs(k1 - k2) <= eps))
    return x, y, valid


def _np_on_extension(x1, y1, x2, y2, px, py):
    return ((px - x2) * (x2 - x1) > 0) | ((py - y2) * (y2 - y1) > 0)


def _np_on_reverse_extension(x1, y1, x2, y2, px, py):
    return ((px - x1) * (x1 - x2) > 0) | ((py - y1) * (y1 - y2) > 0)


# --------------------------------------------------------------------------- Path2D.Float


class GeneralPath:
    """``java.awt.geom.GeneralPath`` as built by ``Polygon.done()``: a float32 ``Path2D.Float``.

    ``done()`` issues ``moveTo(v0)``, ``lineTo(v0)``, ``lineTo(v1)``, ..., ``lineTo(v_{n-1})`` with
    every coordinate cast to ``float``; ``contains`` closes the path implicitly and applies the
    ``WIND_NON_ZERO`` rule.  The zero-length ``v0 -> v0`` and a horizontal closing edge never cross, so
    the crossing count is the sum over the ``n`` closed-polygon edges.
    """

    def __init__(self, vertices: Sequence[Point]):
        f = np.asarray(vertices, dtype=np.float64).astype(np.float32).astype(np.float64)
        self.fx0 = f[:, 0].copy()
        self.fy0 = f[:, 1].copy()
        self.fx1 = np.roll(self.fx0, -1)
        self.fy1 = np.roll(self.fy0, -1)
        self._dir = np.where(self.fy0 < self.fy1, 1, -1).astype(np.int32)
        self._coords = [(float(x), float(y)) for x, y in f]

    def point_crossings(self, px: float, py: float) -> int:
        """``Path2D.Float.pointCrossings`` (scalar reference)."""
        c = self._coords
        movx, movy = c[0]
        curx, cury = movx, movy
        crossings = 0
        for endx, endy in [c[0]] + c[1:]:  # lineTo(v0), lineTo(v1), ...
            crossings += point_crossings_for_line(px, py, curx, cury, endx, endy)
            curx, cury = endx, endy
        if cury != movy:
            crossings += point_crossings_for_line(px, py, curx, cury, movx, movy)
        return crossings

    def contains_scalar(self, x: float, y: float) -> bool:
        """``Path2D.contains(x, y)`` (scalar reference)."""
        if not (x * 0.0 + y * 0.0 == 0.0):
            return False
        return self.point_crossings(x, y) != 0

    def contains_many(self, px, py) -> np.ndarray:
        """Vectorised ``contains`` for arrays of points (any shape); returns a bool array of that shape."""
        px = np.asarray(px, dtype=np.float64)
        py = np.asarray(py, dtype=np.float64)
        if px.shape != py.shape:
            px, py = np.broadcast_arrays(px, py)
        shape = px.shape
        px = px.reshape(-1, 1)
        py = py.reshape(-1, 1)
        x0, y0, x1, y1 = self.fx0, self.fy0, self.fx1, self.fy1
        with np.errstate(all="ignore"):
            finite = (px[:, 0] * 0.0 + py[:, 0] * 0.0) == 0.0
            live = ~((py < y0) & (py < y1)) & ~((py >= y0) & (py >= y1)) & ~((px >= x0) & (px >= x1))
            left = (px < x0) & (px < x1)
            xint = x0 + (py - y0) * (x1 - x0) / (y1 - y0)
            cross = live & (left | ~(px >= xint))
        crossings = (cross * self._dir).sum(axis=1)
        return (finite & (crossings != 0)).reshape(shape)

    def contains(self, x: float, y: float) -> bool:
        """``Path2D.contains(x, y)``."""
        return bool(self.contains_many(np.float64(x), np.float64(y)))


# --------------------------------------------------------------------------- Polygon


class Polygon:
    """A simple polygon as loaded by ``drawable.Polygon`` (vertices in the y-up ``.dat`` frame, CCW).

    ``edges[i]`` is ``v_i -> v_{i+1 mod n}`` (``getLineSegmentArray``).  Instances are treated as
    immutable; derived data (turns, cuts, visibility tables) is cached on the instance.

    The vertices must be finite, counter-clockwise (positive signed area) and form a simple polygon
    (no zero-length edge, no two edges touching except neighbours at their shared vertex); otherwise
    ``ValueError`` is raised at once.  ``validate=False`` skips these checks (Java has none).
    """

    def __init__(self, vertices: Iterable[Sequence[float]], name: Optional[str] = None, *, validate: bool = True):
        vs = tuple((float(v[0]), float(v[1])) for v in vertices)
        if len(vs) < 3:
            raise ValueError("a polygon needs at least 3 vertices")
        self.vertices: Tuple[Point, ...] = vs
        self.n = len(vs)
        self.name = name
        n = self.n
        self.edges: Tuple[Line, ...] = tuple((vs[i][0], vs[i][1], vs[(i + 1) % n][0], vs[(i + 1) % n][1])
                                             for i in range(n))
        a = np.asarray(vs, dtype=np.float64)
        self.vx = a[:, 0].copy()
        self.vy = a[:, 1].copy()
        self.ex1, self.ey1 = self.vx, self.vy
        self.ex2 = np.roll(self.vx, -1)
        self.ey2 = np.roll(self.vy, -1)
        self.general_path = GeneralPath(vs)
        self._cache: Dict[object, object] = {}
        if validate:
            self._validate()

    def _validate(self) -> None:
        """Every algorithm here (the reflex test ``relativeCCW``, the cut directions, the visibility scan) assumes
        a simple polygon listed counter-clockwise in the y-up frame, as the ``.dat`` files are; Java checks
        neither and fails late (an inflection ray that hits no edge, a NullPointerException)."""
        if not all(math.isfinite(c) for v in self.vertices for c in v):
            raise ValueError("polygon vertices must be finite")
        area = self.signed_area()
        if not area > 0:
            raise ValueError(f"vertices must be listed counter-clockwise in the y-up frame (got signed area "
                             f"{area!r}); reverse the vertex list")
        bad = _self_intersection(self)
        if bad is not None:
            raise ValueError(f"polygon is not simple: edges {bad[0]} and {bad[1]} intersect")

    def __repr__(self) -> str:
        return f"Polygon(name={self.name!r}, n={self.n})"

    def __len__(self) -> int:
        return self.n

    def edge(self, i: int) -> Line:
        """``getLineSegmentArray()[i mod n]``."""
        return self.edges[i % self.n]

    @property
    def turn(self) -> Tuple[int, ...]:
        """``turn[i] = relativeCCW(e_{i-1}, v_{i+1})``: ``+1`` reflex, ``-1`` convex, ``0`` collinear."""
        t = self._cache.get("turn")
        if t is None:
            n = self.n
            t = tuple(relative_ccw(self.edges[i - 1], self.edges[i][2:]) for i in range(n))
            self._cache["turn"] = t
        return t  # type: ignore[return-value]

    @property
    def reflex(self) -> List[int]:
        """Indices of the reflex vertices (``turn == 1``)."""
        return [i for i, t in enumerate(self.turn) if t == 1]

    def is_reflex(self, i: int) -> bool:
        """The reflex test used by every cut function: ``e_{i-1}.relativeCCW(v_{i+1}) == 1``."""
        return self.turn[i % self.n] == 1

    def signed_area(self) -> float:
        """Shoelace area (positive for the counter-clockwise, y-up ``.dat`` files)."""
        return 0.5 * float(np.sum(self.ex1 * self.ey2 - self.ex2 * self.ey1))

    def contains(self, p: Point) -> bool:
        """``Polygon.pointInPolygon`` = ``GeneralPath.contains`` (float32 vertices, non-zero rule)."""
        return self.general_path.contains(p[0], p[1])

    point_in_polygon = contains

    def contains_many(self, px, py) -> np.ndarray:
        """Vectorised :meth:`contains`."""
        return self.general_path.contains_many(px, py)

    def seg_in_polygon(self, seg: Line) -> bool:
        """``Algorithm.segInPolygon(this, seg)`` (see :func:`seg_in_polygon`)."""
        return seg_in_polygon(self, seg)

    def _vertex_tables(self, eps: float = EPSILON):
        """Per-polygon tables for segments ending at a vertex (cached).

        ``on[j]``: ``v_j`` is within eps of some edge; ``inside[j]``: ``contains(v_j)``;
        ``rccw[k, j] = relativeCCW(e_k, v_j)``.
        """
        t = self._cache.get(("vertex_tables", eps))
        if t is None:
            vx, vy = self.vx, self.vy
            d = _np_pt_seg_dist(self.ex1[:, None], self.ey1[:, None], self.ex2[:, None], self.ey2[:, None],
                                vx[None, :], vy[None, :])
            on = np.any(d < eps, axis=0)
            inside = self.contains_many(vx, vy)
            rccw = _np_rccw(self.ex1[:, None], self.ey1[:, None], self.ex2[:, None], self.ey2[:, None],
                            vx[None, :], vy[None, :])
            t = (on, inside, rccw)
            self._cache[("vertex_tables", eps)] = t
        return t


def _self_intersection(poly: Polygon) -> Optional[Tuple[int, int]]:
    """First pair of edges ``(i, j)`` that intersect although they should not, or None for a simple polygon.

    Non-adjacent edges must be disjoint (touching counts); adjacent edges may only share their common
    vertex (they must not fold back onto each other); no edge may have zero length.
    """
    n = poly.n
    x1, y1, x2, y2 = poly.ex1, poly.ey1, poly.ex2, poly.ey2
    zero = np.flatnonzero((x1 == x2) & (y1 == y2))
    if zero.size:
        k = int(zero[0])
        return k, (k + 1) % n
    X1, Y1, X2, Y2 = x1[:, None], y1[:, None], x2[:, None], y2[:, None]
    o1 = (X2 - X1) * (y1[None, :] - Y1) - (Y2 - Y1) * (x1[None, :] - X1)  # orient(e_i, start of e_j)
    o2 = (X2 - X1) * (y2[None, :] - Y1) - (Y2 - Y1) * (x2[None, :] - X1)  # orient(e_i, end of e_j)
    s1, s2 = np.sign(o1), np.sign(o2)
    straddle = s1 * s2 <= 0  # e_j touches or crosses the line of e_i
    hit = straddle & straddle.T
    # all four orientations zero: collinear, intersect iff the bounding boxes overlap
    col = (s1 == 0) & (s2 == 0)
    col = col & col.T
    box = ((np.maximum(X1, X2) >= np.minimum(x1, x2)[None, :]) & (np.maximum(x1, x2)[None, :] >= np.minimum(X1, X2))
           & (np.maximum(Y1, Y2) >= np.minimum(y1, y2)[None, :]) & (np.maximum(y1, y2)[None, :] >= np.minimum(Y1, Y2)))
    hit = np.where(col, box, hit)
    idx = np.arange(n)
    nxt = (idx + 1) % n
    hit[idx, idx] = False
    # adjacent edges: allowed to meet at their shared vertex only; a fold-back is collinear and overlapping
    # with the next edge pointing back along this one
    dot = (x2 - x1) * (x2[nxt] - x1[nxt]) + (y2 - y1) * (y2[nxt] - y1[nxt])
    fold = col[idx, nxt] & (dot < 0)
    hit[idx, nxt] = fold
    hit[nxt, idx] = fold
    if n == 3:
        return None if not hit.any() else tuple(int(v) for v in np.argwhere(hit)[0])  # type: ignore[return-value]
    ii, jj = np.nonzero(np.triu(hit))
    if ii.size == 0:
        return None
    return int(ii[0]), int(jj[0])


# --------------------------------------------------------------------------- segInPolygon


def seg_in_polygon_scalar(poly: Polygon, seg: Line, eps: float = EPSILON) -> bool:
    """``Algorithm.segInPolygon``, scalar reference (see :func:`seg_in_polygon`)."""
    gp = poly.general_path
    x1, y1, x2, y2 = seg
    p1, p2 = (x1, y1), (x2, y2)
    if not gp.contains_scalar((x1 + x2) / 2, (y1 + y2) / 2):
        return False
    p1_on = p2_on = False
    for e in poly.edges:
        if pt_seg_dist(e, p1) < eps:
            p1_on = True
        if pt_seg_dist(e, p2) < eps:
            p2_on = True
    if (p1_on or gp.contains_scalar(x1, y1)) and (p2_on or gp.contains_scalar(x2, y2)):
        for e in poly.edges:
            if lines_intersect(seg, e):
                p = get_intersect(seg, e, eps)
                if p is not None and distance(p, p1) > eps and distance(p, p2) > eps:
                    return False
        return True
    return False


def seg_in_polygon_many(poly: Polygon, x1, y1, x2, y2, eps: float = EPSILON) -> np.ndarray:
    """Vectorised :func:`seg_in_polygon` for ``m`` segments given as four 1-D arrays."""
    x1 = np.asarray(x1, dtype=np.float64).reshape(-1)
    y1 = np.asarray(y1, dtype=np.float64).reshape(-1)
    x2 = np.asarray(x2, dtype=np.float64).reshape(-1)
    y2 = np.asarray(y2, dtype=np.float64).reshape(-1)
    out = np.zeros(x1.shape[0], dtype=bool)
    idx = np.flatnonzero(poly.contains_many((x1 + x2) / 2, (y1 + y2) / 2))
    if idx.size == 0:
        return out
    ex1, ey1, ex2, ey2 = poly.ex1, poly.ey1, poly.ex2, poly.ey2
    a1, b1, a2, b2 = x1[idx], y1[idx], x2[idx], y2[idx]
    on1 = np.any(_np_pt_seg_dist(ex1, ey1, ex2, ey2, a1[:, None], b1[:, None]) < eps, axis=1)
    on2 = np.any(_np_pt_seg_dist(ex1, ey1, ex2, ey2, a2[:, None], b2[:, None]) < eps, axis=1)
    ok = (on1 | poly.contains_many(a1, b1)) & (on2 | poly.contains_many(a2, b2))
    idx = idx[ok]
    if idx.size == 0:
        return out
    a1, b1, a2, b2 = a1[ok], b1[ok], a2[ok], b2[ok]
    out[idx] = _unblocked(poly, a1, b1, a2, b2, eps)
    return out


def _unblocked(poly: Polygon, a1, b1, a2, b2, eps: float, rccw_seg=None, rccw_edge=None) -> np.ndarray:
    """Step (3) of ``segInPolygon`` for segments whose end/mid tests passed: no edge crosses the segment
    at a point farther than eps from both ends.  ``rccw_seg[r, k] = relativeCCW(seg_r, v_k)`` and
    ``rccw_edge[r, k] = relativeCCW(e_k, P1_r) * relativeCCW(e_k, P2_r)`` may be precomputed."""
    ex1, ey1, ex2, ey2 = poly.ex1, poly.ey1, poly.ex2, poly.ey2
    if rccw_seg is None:
        rccw_seg = _np_rccw(a1[:, None], b1[:, None], a2[:, None], b2[:, None], ex1, ey1)
    s_prod = rccw_seg * np.roll(rccw_seg, -1, axis=1)  # relativeCCW(seg, e.P1) * relativeCCW(seg, e.P2)
    if rccw_edge is None:
        rccw_edge = (_np_rccw(ex1, ey1, ex2, ey2, a1[:, None], b1[:, None])
                     * _np_rccw(ex1, ey1, ex2, ey2, a2[:, None], b2[:, None]))
    r, k = np.nonzero((s_prod <= 0) & (rccw_edge <= 0))
    res = np.ones(a1.shape[0], dtype=bool)
    if r.size:
        px, py, valid = _np_get_intersect(a1[r], b1[r], a2[r], b2[r], ex1[k], ey1[k], ex2[k], ey2[k], eps)
        with np.errstate(invalid="ignore"):
            blk = valid & (_np_distance(px, py, a1[r], b1[r]) > eps) & (_np_distance(px, py, a2[r], b2[r]) > eps)
        res[r[blk]] = False
    return res


def seg_in_polygon(poly: Polygon, seg: Line, eps: float = EPSILON) -> bool:
    """``Algorithm.segInPolygon(polygon, seg)``.

    (1) the midpoint must satisfy ``GeneralPath.contains``; (2) each endpoint must be within eps of an
    edge or be contained; (3) no edge may ``intersectsLine`` the segment at a point (``getIntersect``)
    farther than eps from both endpoints.  A segment that grazes a vertex in its interior therefore
    counts as blocked; parallel edges are ignored because ``getIntersect`` returns ``None``.
    """
    return bool(seg_in_polygon_many(poly, [seg[0]], [seg[1]], [seg[2]], [seg[3]], eps)[0])


# --------------------------------------------------------------------------- Path


class Path:
    """``drawable.Path``: a polyline robot path (``addPoint`` order).

    ``segments[i]`` is ``points[i] -> points[i+1]``.  Arc lengths accumulate segment lengths
    left to right exactly as ``Path.getLength`` / ``ptDistFromStart`` do.
    """

    def __init__(self, points: Iterable[Sequence[float]]):
        pts = tuple((float(p[0]), float(p[1])) for p in points)
        if len(pts) < 2:
            raise ValueError("a path needs at least 2 points")
        self.points: Tuple[Point, ...] = pts
        self.segments: Tuple[Line, ...] = tuple((pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1])
                                                for i in range(len(pts) - 1))
        prefix = [0.0]
        d = 0.0
        for s in self.segments:
            d += distance((s[0], s[1]), (s[2], s[3]))
            prefix.append(d)
        self._prefix = tuple(prefix)

    def __repr__(self) -> str:
        return f"Path({list(self.points)!r})"

    def length(self) -> float:
        """``Path.getLength``: sum of the segment lengths, left to right."""
        return self._prefix[-1]

    def pt_dist_on_segment(self, p: Point, i: int) -> float:
        """Private ``ptDistFromStart(p, lines, i)``: length of segments ``< i`` plus ``|p - P1_i|``."""
        s = self.segments[i]
        return self._prefix[i] + distance((s[0], s[1]), p)

    def pt_dist_from_start(self, p: Point, eps: float = EPSILON) -> float:
        """``Path.ptDistFromStart(p)``: arc length to ``p``; ``nan`` if ``p`` is on no segment.

        Quirk B17 kept: all segments are scanned and the **last** one within eps wins (wrong for
        self-crossing paths; near a waypoint the next segment wins).  Use :meth:`pt_dist_on_segment`
        when the segment is known.
        """
        d = math.nan
        for i, s in enumerate(self.segments):
            if pt_seg_dist(s, p) < eps:
                d = self.pt_dist_on_segment(p, i)
        return d

    def point_on_path(self, p: Point, eps: float = EPSILON) -> bool:
        """``Path.pointOnPath``."""
        return any(pt_seg_dist(s, p) < eps for s in self.segments)

    def intersect_point(self, line: Line, eps: float = EPSILON) -> Optional[Point]:
        """``Path.getIntersectPoint(Line2D)``: the crossing of ``line`` with the **first** path segment it
        crosses (segment order, not arc-length order among several crossings of one segment), or None."""
        for s in self.segments:
            p = get_segment_intersect(s, line, eps)
            if p is not None:
                return p
        return None

    def _hits(self, lines: Sequence[Line], eps: float):
        """``(seg_index, line_index, point)`` of every ``getSegmentIntersect(segment, line)``, in Java loop order."""
        if not lines:
            return
        L = np.asarray(lines, dtype=np.float64).reshape(-1, 4)
        cx1, cy1, cx2, cy2 = L[:, 0], L[:, 1], L[:, 2], L[:, 3]
        for i, (sx1, sy1, sx2, sy2) in enumerate(self.segments):
            x, y, valid = _np_get_intersect(sx1, sy1, sx2, sy2, cx1, cy1, cx2, cy2, eps)
            with np.errstate(invalid="ignore"):
                ok = (valid & (_np_pt_seg_dist(sx1, sy1, sx2, sy2, x, y) < eps)
                      & (_np_pt_seg_dist(cx1, cy1, cx2, cy2, x, y) < eps))
            for j in np.flatnonzero(ok):
                yield i, int(j), (float(x[j]), float(y[j]))

    def intersect_points(self, lines: Sequence[Line], eps: float = EPSILON) -> List[Point]:
        """``Path.getIntersectPoints(Line2D[])``: crossings ordered by arc length.

        Keys are exact doubles in a ``TreeMap``: equal distances keep only the **last** crossing (B2).
        """
        m: Dict[float, Point] = {}
        for i, _, p in self._hits(lines, eps):
            m[self.pt_dist_on_segment(p, i)] = p
        return [m[k] for k in sorted(m)]

    def all_critical_points(self, cuts: Sequence["Cut"], eps: float = EPSILON) -> List["CriticalPoint"]:
        """``Path.getAllCriticalPoints(cuts)``, see :func:`shadowinfo.polygon.cuts.all_critical_points`."""
        from .cuts import all_critical_points

        return all_critical_points(self, cuts, eps)
