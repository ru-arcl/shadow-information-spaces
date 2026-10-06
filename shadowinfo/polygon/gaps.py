"""Gap tracking along a robot path: shadow components, their IDs and their links.

Port of ``Algorithm.getGaps`` and its helpers (``collapsePhysicalGaps``, ``samePhysicalGap``,
``vertexInGap``, ``lineInGap``) together with ``gap.Gap`` and ``IDGenerator``
(``docs/notes/original_java.md`` §2.7).  The robot walks the path; at every critical point (a
crossing of a cut, plus the start and the end, :func:`~shadowinfo.polygon.cuts.all_critical_points`)
the physical gaps are sampled just *after* the crossing.  Crossings of general inflections and of
bitangents are the component events of T-RO 2012 Sec. II-A (appear / disappear, split / merge); the
other crossings only change how a gap is described.  The result is ``Gap[][]``: one gap set per
interval between component events, each gap with an integer ID (never reused), the IDs it turns
into at the next event (``toGapSet``) and, on element 0 only, the relative time of the set.

Two modes (``compat``):

* ``"java"`` reproduces the Java bookkeeping literally, including the quirks B3-B5 and the crashes:
  :class:`GapTrackingError` is raised at the same critical point, with the same Java exception class
  and source line, as the original (``tests/fixtures/java``: ``14_PP4``, ``13_c_second``,
  ``13_c_tail``, ``14_fig_ICRA08-Fig8b``).
* ``"fixed"`` (default) does every step Java's way (same rules for appear / disappear / split /
  merge, same ID order, same times and iEdges) and keeps Java's step whenever it is consistent: every
  gap gets one ID, the new IDs are exactly the ones handed out, no linked gap survives and every link
  group is a split (1 -> 2) or a merge (2 -> 1).  Otherwise -- Java throws, or its step is inconsistent
  -- the step is read off the boundary chains hidden behind the gaps (:func:`chain_graph`), which also
  finds events the critical points missed (a crossing dropped by the last-writer-wins rule B2, a
  perturbation that overshoots the next crossing).  It never throws, and at every sample the shadows
  alive are exactly the physical gaps.  On every run where Java does not throw the history is Java's,
  set for set and link for link: the golden runs, and a differential test against the live Java on
  about 1260 runs where it succeeds (random valid paths in the 14 maps, perturbed paper paths, 30
  synthetic polygons).  The exceptions are where Java's own answer is wrong: the set times on a path
  that crosses or retraces itself (B17, below), a sample that Java's perturbation puts outside the
  polygon (:func:`sample_path`), and a sample on a line through two vertices where Java's visibility
  scan is degenerate (:func:`~shadowinfo.polygon.visibility.physical_gaps`).  Of these, only B17 changed
  a history in those runs (one path that retraces a segment).

  ``java_matching=False`` (opt-in) also re-derives from the hidden chains every step where Java's
  consistent step *disagrees* with them.  That happens when Java's ``samePhysicalGap`` match is
  ambiguous: two adjacent small gaps both lie within one edge of ``pgs[k][0]``, the last one wins and
  every label moves one gap along the boundary (B3).  The labels then stay on their physical shadows
  (the simulator, :mod:`~shadowinfo.polygon.simulate`, needs that for its ground truth), but the
  history is no longer Java's: on the golden runs ``P5`` / ``P5s`` (one step) and the T-RO Fig. 15(b)
  runs ``P14b`` / ``fig_TRO-Fig15b`` (critical points 124-130 and 219-225) 5, 5 and 29 of the events
  name different shadows, and on random valid paths in the 14 maps about 18% of the runs differ, about
  12% with different final shadow IDs.  Event types, set times and the number of IDs are unchanged.

  B17: ``compat="fixed"`` dates each set by the arc length of its critical point on the segment it
  was found on, where Java's ``ptDistFromStart`` takes the *last* segment within eps.  On a path that
  crosses or retraces itself Java's times can repeat or go backwards (and its ``TreeMap``-keyed oracle
  then drops events, B9); the fixed times are increasing, so on such paths the times (only) differ.

Java ``HashSet<Integer>`` iteration order is reproduced (:func:`java_hashset_order`): it decides
which split child is printed first and which one the oracle fills first.
"""

from __future__ import annotations

import math
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .cuts import Bitangent, CriticalPoint, Cut, CutType, GeneralInflection, all_critical_points, all_cuts
from .geometry import (COMPAT_MODES, EPSILON, PERTURB, GeometryError, Path, Point, Polygon, _check_compat,
                       distance, get_slope, pt_seg_dist, purturb_point_along_seg, relative_ccw)
from .io import java_double_str
from .visibility import PhysicalGap, physical_gaps, point_status, repair_gaps

__all__ = [
    "GapTrackingError", "JavaIntHashSet", "java_hashset_order", "IDGenerator", "Gap", "vertex_in_gap",
    "line_in_gap", "same_physical_gap", "collapse_physical_gaps", "PathSamples", "sample_path", "GapHistory",
    "get_gaps", "format_gap", "format_gap_set", "format_gap_sets", "hidden_chain", "chain_graph", "COMPAT_MODES",
]


# --------------------------------------------------------------------------- errors


class GapTrackingError(RuntimeError):
    """``getGaps`` throws: the Java gap bookkeeping is inconsistent at a critical point (quirk B3).

    Attributes:

    * ``critical_index`` -- index into the critical points (``getAllCriticalPoints``) of the step that
      throws, i.e. the loop variable ``i`` of ``getGaps``;
    * ``critical_type`` -- the :class:`~shadowinfo.polygon.cuts.CutType` of that critical point;
    * ``java_exception`` -- ``"java.lang.NullPointerException"`` or
      ``"java.lang.ArrayIndexOutOfBoundsException"``;
    * ``java_line`` -- the line of ``Algorithm.java`` that throws (``collapsePhysicalGaps`` lines
      371-460, or ``getGaps`` itself);
    * ``getgaps_line`` -- the line of ``getGaps`` on the stack (the ``collapsePhysicalGaps`` call, or
      the throwing statement itself);
    * ``phase`` -- which branch of ``getGaps`` was running.
    """

    def __init__(self, message: str, *, critical_index: int = -1, critical_type: Optional[CutType] = None,
                 java_exception: str = "", java_line: int = -1, getgaps_line: int = -1, phase: str = ""):
        super().__init__(message)
        self.critical_index = critical_index
        self.critical_type = critical_type
        self.java_exception = java_exception
        self.java_line = java_line
        self.getgaps_line = getgaps_line
        self.phase = phase


class _JavaThrow(Exception):
    """Internal: a Java runtime exception at ``line`` of ``Algorithm.java``."""

    def __init__(self, exc: str, line: int, message: str = ""):
        super().__init__(f"{exc} at Algorithm.java:{line}: {message}")
        self.exc = exc
        self.line = line
        self.message = message


_NPE = "java.lang.NullPointerException"
_AIOOBE = "java.lang.ArrayIndexOutOfBoundsException"


def _aioobe(line: int, index: int, length: int) -> _JavaThrow:
    return _JavaThrow(_AIOOBE, line, f"Index {index} out of bounds for length {length}")


# --------------------------------------------------------------------------- java.util.HashSet<Integer>


def _spread(h: int) -> int:
    """``HashMap.hash(Integer)``: ``h ^ (h >>> 16)`` on the 32-bit ``Integer.hashCode``."""
    h &= 0xFFFFFFFF
    return h ^ (h >> 16)


class JavaIntHashSet:
    """``java.util.HashSet<Integer>`` with Java's iteration order (OpenJDK 8+ ``HashMap``).

    Iteration is by bucket ``hash & (capacity - 1)`` (capacity 16, doubled when the size exceeds
    ``0.75 * capacity`` or when a bucket reaches 8 entries while the capacity is below 64), and by
    insertion order inside a bucket (resizing splits buckets without reordering them).  Elements are
    never removed (``Gap.toGapSet`` only grows).  Tree bins (8+ colliding entries at capacity >= 64)
    are not modelled and raise ``NotImplementedError``; gap IDs never get near that.
    """

    __slots__ = ("_items", "_cap")

    def __init__(self, items: Iterable[int] = ()):
        self._items: List[int] = []
        self._cap = 16
        for x in items:
            self.add(x)

    def add(self, x: int) -> bool:
        x = int(x)
        if x in self._items:
            return False
        mask = self._cap - 1
        b = _spread(x) & mask
        chain = sum(1 for y in self._items if _spread(y) & mask == b)
        self._items.append(x)
        if chain >= 8:  # treeifyBin: resize instead while the table is small
            if self._cap < 64:
                self._cap *= 2
            else:  # pragma: no cover - needs 9 IDs colliding in one bucket of a 64-table
                raise NotImplementedError("HashSet tree bins are not modelled")
        if len(self._items) > 0.75 * self._cap:
            self._cap *= 2
        return True

    def __contains__(self, x: object) -> bool:
        return x in self._items

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        mask = self._cap - 1
        return iter(sorted(self._items, key=lambda y: _spread(y) & mask))  # stable: insertion order in a bucket

    def to_list(self) -> List[int]:
        """``toArray()``: the elements in iteration order."""
        return list(self)

    def __repr__(self) -> str:
        return f"JavaIntHashSet({self.to_list()!r})"


def java_hashset_order(ids: Iterable[int]) -> List[int]:
    """Iteration order of a ``new HashSet<Integer>()`` after adding ``ids`` in this order.

    For fewer than 13 non-negative IDs below 65 536 this is "ascending by ``id mod 16``, ties in
    insertion order"; e.g. ``[16, 15]`` for the insertions ``15, 16``.
    """
    return JavaIntHashSet(ids).to_list()


# --------------------------------------------------------------------------- Gap, IDGenerator


class IDGenerator:
    """``geometry.IDGenerator``: ``next_id()`` returns 1, 2, 3, ... (one generator per ``getGaps``)."""

    def __init__(self) -> None:
        self.last = 0

    def next_id(self) -> int:
        self.last += 1
        return self.last


class Gap:
    """``gap.Gap`` / ``gap.PhysicalGap``: one shadow component in a gap set.

    * ``id`` -- the shadow label (``-1`` = not assigned yet);
    * ``iedge`` -- Java's "invariant edge", a cosmetic boundary-edge label printed by ``toString``
      (quirk B5: ad-hoc clamping, ``0`` for gaps whose hidden chain wraps past vertex 0);
    * ``relative_time`` -- arc-length fraction of the event that opened the set; Java sets it on
      element 0 of each set only (the others keep ``0.0``);
    * ``state`` -- scenario state, e.g. :class:`~shadowinfo.polygon.events.NSEState` (``None`` = unset);
    * ``to_gaps`` -- the IDs this gap turns into at the next event (split: 2 children, merge: the
      merged gap), in Java ``HashSet`` iteration order;
    * ``physical`` -- the sampled :class:`~shadowinfo.polygon.visibility.PhysicalGap` geometry.
    """

    __slots__ = ("id", "iedge", "relative_time", "state", "_to", "physical")

    def __init__(self, physical: Optional[PhysicalGap] = None, id: int = -1, iedge: int = -1):
        self.id = id
        self.iedge = iedge
        self.relative_time = 0.0
        self.state: object = None
        self._to = JavaIntHashSet()
        self.physical = physical

    # PhysicalGap accessors used by the tracker
    @property
    def start_edge(self) -> int:
        return self.physical.start_edge  # type: ignore[union-attr]

    @property
    def end_edge(self) -> int:
        return self.physical.end_edge  # type: ignore[union-attr]

    def full_start_edge(self, n: int) -> int:
        return self.physical.full_start_edge(n)  # type: ignore[union-attr]

    def full_end_edge(self, n: int) -> int:
        return self.physical.full_end_edge(n)  # type: ignore[union-attr]

    def add_gap(self, gid: int) -> None:
        """``Gap.addGap``: link this gap to ``gid`` in the next gap set."""
        self._to.add(gid)

    @property
    def to_gaps(self) -> Tuple[int, ...]:
        """``toGapSet`` in Java ``HashSet`` iteration order."""
        return tuple(self._to)

    def to_list(self) -> list:
        """``[id, iEdge, toGapSet]`` (golden-fixture format)."""
        return [self.id, self.iedge, list(self.to_gaps)]

    def __str__(self) -> str:
        return format_gap(self)

    def __repr__(self) -> str:
        return f"Gap(id={self.id}, iedge={self.iedge}, to={list(self.to_gaps)}, physical={self.physical})"


def format_gap(g: Gap) -> str:
    """``Gap.toString()``: ``[id, iEdge]``, ``[id, iEdge | a b ]``; a state is printed after iEdge
    (``[19, 0| 1]`` for a contaminated gap in the never-see-evader scenario)."""
    s = f"[{g.id}, {g.iedge}"
    if g.state is not None:
        s += str(g.state)
    to = g.to_gaps
    if to:
        s += " | " + "".join(f"{x} " for x in to)
    return s + "]"


def format_gap_set(gaps: Sequence[Gap], with_time: bool = True, time: Optional[float] = None) -> str:
    """One printed gap set.

    ``with_time=True``: the ``ProjectPanel`` / ``ProjectPanel5`` line, ``Double.toString(relativeTime)``
    of element 0 (or ``time``) followed by ``gap.toString() + " "`` for each gap.  ``with_time=False``:
    the ``ProjectPanel2`` / ``ProjectPanel3`` line (gaps only, normally after never-see-evader states
    were set).
    """
    body = "".join(f"{format_gap(g)} " for g in gaps)
    if not with_time:
        return body
    if time is None:
        time = gaps[0].relative_time if gaps else 0.0
    return f"{java_double_str(time)} {body}"


def format_gap_sets(history: Sequence[Sequence[Gap]], with_time: bool = True) -> List[str]:
    """:func:`format_gap_set` for every set of a gap history (times from ``history.times`` if present)."""
    times = getattr(history, "times", None)
    return [format_gap_set(gs, with_time, None if times is None else times[k]) for k, gs in enumerate(history)]


# --------------------------------------------------------------------------- Java helpers


def vertex_in_gap(iv: int, gap, n: int) -> bool:
    """``Algorithm.vertexInGap``: ``start < iv <= end`` cyclically, on the *raw* start/end edges."""
    s, e = gap.start_edge, gap.end_edge
    return (s < iv <= e) or (e < s and ((s < iv < n) or (0 <= iv <= e)))


def line_in_gap(e: int, gap, n: int) -> bool:
    """``Algorithm.lineInGap``: ``fullStart <= e <= fullEnd`` cyclically (edge ``e`` entirely hidden)."""
    fs, fe = gap.full_start_edge(n), gap.full_end_edge(n)
    return (fs <= e <= fe) or (fe < fs and ((fs <= e < n) or (0 <= e <= fe)))


def same_physical_gap(a, b, n: int, compat: str = "fixed") -> bool:
    """``Algorithm.samePhysicalGap``: start edges and end edges each differ by at most one.

    Quirk B4: Java's wrap-around test ``s0 + s1 == n - 1`` accepts *any* pair of edges whose indices
    sum to ``n - 1``.  ``compat="java"`` keeps it; ``"fixed"`` tests cyclic distance ``<= 1``.
    """
    s0, s1, e0, e1 = a.start_edge, b.start_edge, a.end_edge, b.end_edge
    if compat == "java":
        return (abs(s0 - s1) <= 1 or s0 + s1 == n - 1) and (abs(e0 - e1) <= 1 or e0 + e1 == n - 1)
    return (s0 - s1) % n in (0, 1, n - 1) and (e0 - e1) % n in (0, 1, n - 1)


def _iedges(pgs: Sequence[Sequence[Gap]], n: int) -> Dict[int, int]:
    """The "invariant edge" per ID (``collapsePhysicalGaps`` lines 396-454, quirk B5).

    Raises :class:`_JavaThrow` where Java dereferences a missing map entry.
    """
    num = len(pgs[0])
    cse: Dict[int, int] = {}
    cee: Dict[int, int] = {}
    for g in pgs[0]:
        cse[g.id] = g.full_start_edge(n)
        cee[g.id] = g.full_end_edge(n)
    for i in range(1, len(pgs)):
        for j in range(num):
            if j >= len(pgs[i]):
                raise _aioobe(411, j, len(pgs[i]))
            g = pgs[i][j]
            if g.id not in cse:
                raise _JavaThrow(_NPE, 411, f"no gap {g.id} in the first sample")
            c = cse[g.id]
            se = g.full_start_edge(n)
            if se != c:
                if se > c:
                    if se == c + 1:
                        c = se
                elif se + c == n:
                    c = se
            cse[g.id] = c
            c = cee[g.id]
            ee = g.full_end_edge(n)
            if ee != c:
                if ee > c:
                    if ee + c == n:
                        c = ee
                elif ee + 1 == c:
                    c = ee
            cee[g.id] = c
    return {gid: ((cse[gid] + cee[gid]) // 2 if cee[gid] >= cse[gid] else 0) for gid in cse}


def collapse_physical_gaps(pgs: Sequence[List[Gap]], id_gen: Optional[IDGenerator], n: int) -> List[Gap]:
    """``Algorithm.collapsePhysicalGaps`` (literal): propagate IDs through consecutive samples.

    Sample ``k + 1`` takes its IDs from sample ``k`` by a cyclic shift ``j``: the last ``j`` with
    ``samePhysicalGap(pgs[k][0], pgs[k+1][j])`` (Java wrap test, B4).  If ``pgs[0][0]`` has no ID,
    all of ``pgs[0]`` get fresh IDs from ``id_gen``.  Then every gap of the last sample gets its
    ``iEdge`` and the last sample is returned (an empty first sample returns a new empty list).

    Quirk B3: every sample must have ``len(pgs[0])`` gaps in the same cyclic order.  Otherwise Java
    throws ``ArrayIndexOutOfBoundsException`` (lines 380, 390) or ``NullPointerException`` (lines
    371, 411, 460); this function raises the internal ``_JavaThrow`` with that line.
    """
    num = len(pgs[0])
    if num == 0:
        return []
    if pgs[0][0].id == -1:
        if id_gen is None:
            raise _JavaThrow(_NPE, 371, "idGen is null")
        for g in pgs[0]:
            g.id = id_gen.next_id()
    for i in range(len(pgs) - 1):
        pg0 = pgs[i][0]
        cur, nxt = pgs[i], pgs[i + 1]
        for j in range(num):
            if j >= len(nxt):
                raise _aioobe(380, j, len(nxt))
            if same_physical_gap(pg0, nxt[j], n, "java"):
                for k in range(num):
                    t = k + j
                    if t >= num:
                        t -= num
                    if t >= len(nxt):
                        raise _aioobe(390, t, len(nxt))
                    nxt[t].id = cur[k].id
    ce = _iedges(pgs, n)
    last = pgs[-1]
    for g in last:
        if g.id not in ce:
            raise _JavaThrow(_NPE, 460, f"no invariant edge for gap {g.id}")
        g.iedge = ce[g.id]
    return last


# --------------------------------------------------------------------------- sampling the path


@dataclass
class PathSamples:
    """Everything ``getGaps`` computes before its bookkeeping, one entry per critical point.

    ``points[i]`` is ``purturbPointAlongSeg(critical_points[i].point, seg, PERTURB, far=true)`` and
    ``physical[i]`` the physical gaps there (:func:`~shadowinfo.polygon.visibility.physical_gaps` of the
    mode); if the perturbation has no answer (B1 in ``compat="java"``), both are ``None`` and ``errors[i]``
    holds the :class:`GeometryError`.  ``compat="fixed"`` moves a sample that Java's perturbation puts
    outside the polygon (past a waypoint or the path end that lies within ``PERTURB * sqrt(2)`` of a wall,
    where Java's scan silently returns no gaps) back along the same direction, halving the step until it
    is inside, or onto the critical point itself.
    """

    polygon: Polygon
    path: Path
    cuts: List[Cut]
    critical_points: List[CriticalPoint]
    points: List[Optional[Point]]
    physical: List[Optional[List[PhysicalGap]]]
    errors: Dict[int, GeometryError] = field(default_factory=dict)
    compat: str = "fixed"
    mode_independent: bool = True
    """True when the samples are the same in both modes and may be reused for either: no critical point
    lies on a segment of slope exactly +-1 (B1), no sample is outside the polygon and the fixed physical
    gaps equal Java's at every sample (they differ only at points on a line through two vertices)."""
    moved: Dict[int, Point] = field(default_factory=dict)
    """``compat="fixed"``: Java's sample point of every critical point whose sample was moved inside."""


def _inside_sample(poly: Polygon, cp: CriticalPoint, p: Point, eps: float) -> Point:
    """A sample on the ray from the critical point through ``p`` (Java's perturbed point, outside the polygon)
    that is in the closed polygon: ``p`` pulled back by halving, or the critical point itself."""
    x0, y0 = cp.point
    dx, dy = p[0] - x0, p[1] - y0
    for k in range(1, 40):
        f = 0.5 ** k
        c = (x0 + f * dx, y0 + f * dy)
        if point_status(poly, c, eps) == "inside":
            return c
    return (x0, y0)


_SAMPLE_CACHE_SIZE = 32


def sample_path(poly: Polygon, path: Path, compat: str = "fixed", eps: float = EPSILON,
                purturb: float = PERTURB, cuts: Optional[Sequence[Cut]] = None) -> PathSamples:
    """Critical points of ``path`` and the physical gaps just after each (``getGaps`` lines 76-94).

    This is the expensive part of :func:`get_gaps` (one visibility scan per critical point: about 2 s for the
    671 critical points of ``P14b``).  With the polygon's own cuts (``cuts=None``) the result is cached on the
    polygon (the last :data:`_SAMPLE_CACHE_SIZE` paths per mode), so treat it as read-only.
    """
    _check_compat(compat)
    if cuts is None:
        key = (path.points, compat, eps, purturb)
        lru = poly._cache.setdefault("sample_path", OrderedDict())
        hit = lru.get(key)  # type: ignore[attr-defined]
        if hit is None:
            hit = _sample_path(poly, path, compat, eps, purturb, list(all_cuts(poly, eps)))
            lru[key] = hit  # type: ignore[index]
            while len(lru) > _SAMPLE_CACHE_SIZE:  # type: ignore[arg-type]
                lru.popitem(last=False)  # type: ignore[attr-defined]
        else:
            lru.move_to_end(key)  # type: ignore[attr-defined]
        return hit
    return _sample_path(poly, path, compat, eps, purturb, list(cuts))


def _sample_path(poly: Polygon, path: Path, compat: str, eps: float, purturb: float,
                 cuts: List[Cut]) -> PathSamples:
    cps = all_critical_points(path, cuts, eps)
    pts: List[Optional[Point]] = []
    phys: List[Optional[List[PhysicalGap]]] = []
    errors: Dict[int, GeometryError] = {}
    moved: Dict[int, Point] = {}
    independent = all(abs(get_slope(cp.seg, eps)) != 1.0 for cp in cps)
    for i, cp in enumerate(cps):
        try:
            p = purturb_point_along_seg(cp.point, cp.seg, purturb, True, compat=compat, eps=eps)
        except GeometryError as ex:
            errors[i] = ex
            pts.append(None)
            phys.append(None)
            continue
        raw = physical_gaps(poly, p, eps, "java")
        status = point_status(poly, p, eps)
        if status == "outside":
            independent = False
            if compat == "fixed":
                moved[i] = p
                p = _inside_sample(poly, cp, p, eps)
                status = point_status(poly, p, eps)
                raw = physical_gaps(poly, p, eps, "java")
        fixed = repair_gaps(poly, p, raw, status == "boundary") if status != "outside" else raw
        if fixed != raw:
            independent = False
        pts.append(p)
        phys.append(fixed if compat == "fixed" else raw)
    return PathSamples(poly, path, cuts, cps, pts, phys, errors, compat, independent, moved)


# --------------------------------------------------------------------------- the gap history


class GapHistory(list):
    """``Gap[][]`` returned by :func:`get_gaps`: a list of gap sets, plus what produced them.

    * ``times[k]`` -- relative time of set ``k`` (Java stores it on ``set[0]`` only, so an empty set
      has none there);
    * ``sample_indices[k]`` -- the critical point whose sample set ``k`` is (the last sample before
      the next event);
    * ``sample_ids[i]`` -- the IDs of the physical gaps at every critical point ``i`` (in
      ``getPhysicalGaps`` order), after the bookkeeping;
    * ``samples`` -- the :class:`PathSamples`; ``compat`` -- the mode;
    * ``inferred`` -- ``compat="fixed"`` only: ``(critical_index, description)`` of every step that
      was read off the hidden chains instead of done Java's way (Java throws there, its step is
      inconsistent, e.g. an event the critical points missed, or -- with ``java_matching=False`` only --
      it disagrees with the chains).
    """

    def __init__(self, sets: Iterable[List[Gap]] = ()):
        super().__init__(sets)
        self.times: List[float] = []
        self.sample_indices: List[int] = []
        self.sample_ids: List[List[int]] = []
        self.samples: Optional[PathSamples] = None
        self.compat = "fixed"
        self.inferred: List[Tuple[int, str]] = []

    @property
    def max_id(self) -> int:
        """Largest ID handed out (= number of shadow labels, IDs start at 1)."""
        return max((g.id for gs in self for g in gs), default=0)

    @property
    def final_ids(self) -> List[int]:
        """IDs of the last gap set (the shadows alive at the end of the path), in boundary order."""
        return [g.id for g in self[-1]] if self else []

    def lines(self, with_time: bool = True) -> List[str]:
        """Printed gap sets (:func:`format_gap_sets`)."""
        return format_gap_sets(self, with_time)


def _relative_time(path: Path, cp: CriticalPoint, length: float, compat: str, eps: float) -> float:
    """``path.ptDistFromStart(point) / pathLength`` (getGaps line 295).

    ``ptDistFromStart`` keeps the *last* segment within eps of the point (B17).  ``compat="fixed"``
    only considers the segments whose arc length at the point agrees with the critical point's own
    distance (within ``1e-6`` of the path length), which is the same answer unless the path crosses
    itself.
    """
    p = cp.point
    if compat == "java":
        return path.pt_dist_from_start(p, eps) / length
    d = math.nan
    tol = 1e-6 * max(length, 1.0)
    for i, s in enumerate(path.segments):
        if pt_seg_dist(s, p) < eps:
            di = path.pt_dist_on_segment(p, i)
            if abs(di - cp.distance) <= tol:
                d = di
    if math.isnan(d):
        d = cp.distance
    return d / length


def get_gaps(poly: Polygon, path: Path, compat: str = "fixed", *, java_matching: bool = True,
             eps: float = EPSILON, purturb: float = PERTURB, samples: Optional[PathSamples] = None) -> GapHistory:
    """``Algorithm.getGaps(polygon, path)``: the gap history of a path, as gap sets with IDs and links.

    1. ``cuts = getCuts(poly)``; critical points sorted by arc length (start and end included).
    2. At each critical point ``i`` the gaps are sampled at the perturbed point just after it.
    3. A general-inflection crossing is an *appear* if the previous critical point lies on the
       interior side of the inflection's edge (``relativeCCW == -1``; the gap containing that edge as
       a full edge gets a fresh ID), otherwise a *disappear* (the previous gap containing the edge is
       dropped).  A bitangent crossing is a *split* if the previous sample and ``v_{this+1}`` lie on
       the same side of the cut (the previous gap containing ``v_opp`` is the parent; every new gap
       containing one of the edges ``this-1, this, opp-1, opp`` gets a fresh ID and a link from the
       parent), otherwise a *merge* (every previous gap containing one of those edges links to the
       new ID of the gap containing ``v_opp``).  The other gaps keep their IDs.
    4. Between events, IDs are propagated sample to sample; the last sample before each event (and
       the end) becomes a gap set, its element 0 carrying the relative time of the event (or ``0``).

    ``compat="java"`` raises :class:`GapTrackingError` exactly where the Java code throws (see the
    module docstring).  ``compat="fixed"`` (default) never throws, also on paths that leave the polygon.
    ``samples`` may pass a precomputed :func:`sample_path` (of either mode when ``mode_independent``).

    ``java_matching`` (``compat="fixed"`` only, ignored otherwise): True (default) keeps Java's own
    step whenever it is consistent, which reproduces the Java history on every run where Java does not
    throw (up to the B17 set times on a self-crossing path).  False also re-derives from the hidden
    chains the steps where Java's consistent step disagrees with them -- its label rotations: when two
    adjacent small gaps are both within one edge of ``pgs[k][0]``, ``samePhysicalGap`` matches both and
    the last match wins, which shifts every ID one gap along the boundary (quirk B3; on the T-RO Fig.
    15(b) run ``P14b`` at critical points 124-130 and 219-223, on ``P5`` / ``P5s``, and on about 18% of
    random valid paths).  IDs then stay with their physical shadows (what the simulator uses), at the
    price of a history that differs from Java's.
    """
    _check_compat(compat)
    if samples is None:
        lru = poly._cache.get("sample_path")
        other = lru.get((path.points, "java" if compat == "fixed" else "fixed", eps, purturb)) if lru else None
        samples = other if other is not None and other.mode_independent else sample_path(poly, path, compat, eps,
                                                                                            purturb)
    elif samples.compat != compat and not samples.mode_independent:
        raise ValueError(f"samples were computed with compat={samples.compat!r} and differ in compat={compat!r} (B1)")
    if compat == "java":
        return _get_gaps_java(poly, path, samples, eps)
    return _get_gaps_fixed(poly, path, samples, eps, java_matching)


def _event_step(cp: CriticalPoint, i: int, samples: PathSamples, previous: Optional[List[Gap]], pg: List[Gap],
                n: int, edges: Sequence[Tuple[float, float, float, float]], id_gen: IDGenerator) -> None:
    """The event branches of ``getGaps`` (lines 102-290) at critical point ``i``, literally.

    Assigns fresh IDs to the gaps of ``pg`` that an appear / split / merge creates, links the
    split or merging gaps of ``previous`` to them, and collapses the remaining gaps of both to carry
    their IDs over.  Raises :class:`_JavaThrow` (with ``gline``, the ``getGaps`` line on the stack)
    where Java throws.
    """
    cps = samples.critical_points
    gline = 107
    try:
        if cp.cut_type is CutType.GENERAL_INFLECTION:
            cut = cp.cut
            assert isinstance(cut, GeneralInflection)
            if i == 0:
                raise _aioobe(107, -1, len(cps))
            last = cps[i - 1].point
            if previous is None:  # pragma: no cover - a set is always closed before an event
                raise _JavaThrow(_NPE, 147, "previousGaps is null")
            if relative_ccw(edges[cut.from_line], last) == -1:  # a gap appears
                gline = 138
                rest = []
                for g in pg:
                    if line_in_gap(cut.from_line, g, n):
                        g.id = id_gen.next_id()
                    else:
                        rest.append(g)
                collapse_physical_gaps([previous, rest], None, n)
            else:  # a gap disappears
                gline = 162
                rest = [g for g in previous if not line_in_gap(cut.from_line, g, n)]
                collapse_physical_gaps([rest, pg], None, n)
        elif cp.cut_type is CutType.BITANGENT:
            cut = cp.cut
            assert isinstance(cut, Bitangent)
            gline = 174
            if i == 0:
                raise _aioobe(174, -1, len(cps))
            if (i - 1) in samples.errors:
                raise _JavaThrow(_NPE, 174, str(samples.errors[i - 1]))
            last = samples.points[i - 1]
            if previous is None:  # pragma: no cover
                raise _JavaThrow(_NPE, 198, "previousGaps is null")
            prod = relative_ccw(cut.line, cut.curve_to_point) * relative_ccw(cut.line, last)
            itp, iop = cut.this_point, cut.opposite_point
            four = (n - 1 if itp == 0 else itp - 1, itp, n - 1 if iop == 0 else iop - 1, iop)
            if prod == 1:  # same side as v_{this+1}: a gap splits
                gline = 231
                parent = None
                old = []
                for g in previous:
                    if vertex_in_gap(iop, g, n):
                        parent = g
                    else:
                        old.append(g)
                new = []
                for g in pg:
                    if any(line_in_gap(e, g, n) for e in four):
                        g.id = id_gen.next_id()
                        if parent is None:
                            gline = 220
                            raise _JavaThrow(_NPE, 220, "no previous gap contains the opposite vertex")
                        parent.add_gap(g.id)
                    else:
                        new.append(g)
                collapse_physical_gaps([old, new], None, n)
            else:  # gaps merge
                gline = 286
                parents, old = [], []
                for g in previous:
                    (parents if any(line_in_gap(e, g, n) for e in four) else old).append(g)
                new = []
                for g in pg:
                    if vertex_in_gap(iop, g, n):
                        g.id = id_gen.next_id()
                        for p in parents:
                            p.add_gap(g.id)
                    else:
                        new.append(g)
                collapse_physical_gaps([old, new], None, n)
    except _JavaThrow as ex:
        ex.gline = gline  # type: ignore[attr-defined]
        raise


_EVENT_TYPES = (CutType.GENERAL_INFLECTION, CutType.BITANGENT)
_PHASES = {CutType.GENERAL_INFLECTION: "general-inflection branch", CutType.BITANGENT: "bitangent branch"}


def _get_gaps_java(poly: Polygon, path: Path, samples: PathSamples, eps: float) -> GapHistory:
    """The literal port of ``getGaps`` (``compat="java"``)."""
    n = poly.n
    edges = [poly.edge(i) for i in range(n)]
    id_gen = IDGenerator()
    length = path.length()
    rel_time = 0.0
    cps = samples.critical_points
    out = GapHistory()
    out.samples = samples
    out.compat = "java"
    sample_gaps: List[List[Gap]] = []
    pg_vec: List[List[Gap]] = []
    previous: Optional[List[Gap]] = None
    for i, cp in enumerate(cps):
        phase, gline = "sampling", 92
        try:
            if i in samples.errors:
                raise _JavaThrow(_NPE, 92, str(samples.errors[i]))
            pg = [Gap(p) for p in samples.physical[i]]  # type: ignore[union-attr]
            sample_gaps.append(pg)
            pg_vec.append(pg)
            if cp.cut_type in _EVENT_TYPES:
                phase = _PHASES[cp.cut_type]
                try:
                    _event_step(cp, i, samples, previous, pg, n, edges, id_gen)
                except _JavaThrow as ex:
                    gline = ex.gline  # type: ignore[attr-defined]
                    raise
            if i == 0 or cp.cut_type in _EVENT_TYPES:
                rel_time = path.pt_dist_from_start(cp.point, eps) / length
            if i == len(cps) - 1 or cps[i + 1].cut_type in _EVENT_TYPES:
                phase, gline = "closing the gap set", 305
                previous = collapse_physical_gaps(pg_vec, id_gen, n)
                if not previous:
                    gline = 307
                    raise _aioobe(307, 0, 0)
                previous[0].relative_time = rel_time
                out.append(previous)
                out.times.append(rel_time)
                out.sample_indices.append(i)
                pg_vec = []
        except _JavaThrow as ex:
            raise GapTrackingError(
                f"getGaps throws {ex.exc} at critical point {i} ({cp.cut_type.name}, {phase}; "
                f"Algorithm.java:{ex.line}): {ex.message}",
                critical_index=i, critical_type=cp.cut_type, java_exception=ex.exc, java_line=ex.line,
                getgaps_line=gline, phase=phase) from None
    out.sample_ids = [[g.id for g in gs] for gs in sample_gaps]
    return out


# --------------------------------------------------------------------------- compat="fixed"


def _perimeter(poly: Polygon) -> Tuple[float, ...]:
    """Cumulative boundary arc length at each vertex (``cum[n]`` = perimeter)."""
    key = ("perimeter",)
    if key not in poly._cache:
        cum = [0.0]
        for i in range(poly.n):
            s = poly.edge(i)
            cum.append(cum[-1] + distance((s[0], s[1]), (s[2], s[3])))
        poly._cache[key] = tuple(cum)
    return poly._cache[key]  # type: ignore[return-value]


def hidden_chain(poly: Polygon, g: PhysicalGap) -> Tuple[float, float]:
    """The boundary chain hidden behind a gap, as ``(start, length)`` in boundary arc length.

    The chain runs counter-clockwise from ``v_start`` to ``end_point`` (gap edge through
    ``v_start``) or from ``start_point`` to ``v_{end+1}`` (gap edge through ``v_{end+1}``); the vertex
    the gap edge touches is the gap's *anchor*.
    """
    cum = _perimeter(poly)
    per = cum[-1]
    v = poly.vertices
    if g.end_point is not None:
        a = cum[g.start_edge]
        b = cum[g.end_edge] + distance(v[g.end_edge], g.end_point)
    else:
        a = cum[g.start_edge] + distance(v[g.start_edge], g.start_point)  # type: ignore[arg-type]
        b = cum[g.end_edge + 1]
    return a, (b - a) % per


def _core(poly: Polygon, g: PhysicalGap) -> float:
    """Midpoint (boundary arc length) of the fully hidden edge at the gap's anchor vertex."""
    cum = _perimeter(poly)
    e = g.start_edge if g.end_point is not None else g.end_edge
    return 0.5 * (cum[e] + cum[e + 1])


def _inside(x: float, chain: Tuple[float, float], per: float) -> bool:
    d = (x - chain[0]) % per
    return 0.0 < d < chain[1]


def chain_graph(poly: Polygon, a: Sequence[PhysicalGap], b: Sequence[PhysicalGap]) -> List[Tuple[int, int]]:
    """Which gaps of two consecutive samples are the same shadow or its parent / child.

    ``(j, k)`` is an edge when the *core* of one gap (the middle of the fully hidden edge at its
    anchor vertex) lies inside the hidden chain of the other.  Between samples with no component
    event the anchor of a gap stays put or moves to the next reflex vertex along its own chain
    (non-general inflection), while the free end may jump (single tangent), so continuing gaps are
    joined and neighbouring gaps are not.  A splitting gap is joined to both children (each child's
    anchor was hidden by the parent), merging gaps to the merged gap, and an appearing gap (hiding
    boundary that was visible) or a disappearing one to nothing.
    """
    per = _perimeter(poly)[-1]
    ca = [hidden_chain(poly, g) for g in a]
    cb = [hidden_chain(poly, g) for g in b]
    ka = [_core(poly, g) for g in a]
    kb = [_core(poly, g) for g in b]
    return [(j, k) for j in range(len(a)) for k in range(len(b))
            if _inside(kb[k], ca[j], per) or _inside(ka[j], cb[k], per)]


def _java_proposal(cp: CriticalPoint, i: int, samples: PathSamples, prev: List[Gap], n: int,
                   edges: Sequence[Tuple[float, float, float, float]], last_id: int):
    """What Java's bookkeeping does for the step ``i-1 -> i``, on copies.

    Returns ``(ids, links, last_id)`` (IDs of the gaps of sample ``i``, ``{index in prev: [new IDs in
    insertion order]}``, last ID handed out) or ``None`` where Java would throw.
    """
    a = [Gap(g.physical, g.id) for g in prev]
    b = [Gap(p) for p in samples.physical[i]]  # type: ignore[union-attr]
    gen = IDGenerator()
    gen.last = last_id
    try:
        if cp.cut_type in _EVENT_TYPES:
            _event_step(cp, i, samples, a, b, n, edges, gen)
        else:
            collapse_physical_gaps([a, b], gen, n)
    except _JavaThrow:
        return None
    links = {j: list(g._to._items) for j, g in enumerate(a) if len(g._to)}
    return [g.id for g in b], links, gen.last


def _components(na: int, nb: int, pairs: Iterable[Tuple[int, int]]) -> List[Tuple[List[int], List[int]]]:
    """Connected components of a bipartite graph, as ``(sorted a-nodes, sorted b-nodes)``,
    ordered by their first node in boundary order (``a`` nodes first, then unmatched ``b`` nodes)."""
    parent = list(range(na + nb))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for j, k in pairs:
        ra, rb = find(j), find(na + k)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    comps: Dict[int, Tuple[List[int], List[int]]] = {}
    for x in range(na + nb):
        c = comps.setdefault(find(x), ([], []))
        (c[0] if x < na else c[1]).append(x if x < na else x - na)
    return [comps[r] for r in sorted(comps)]


def _proposal_problem(prev: List[Gap], ids: List[int], links: Dict[int, List[int]], new_ids: Sequence[int],
                      graph: Optional[set]) -> Optional[str]:
    """Why Java's step is not kept, or None.  It must give every gap a unique ID, hand out exactly the
    new IDs, not let a linked gap survive, make every link group a split (1 -> 2) or a merge (2 -> 1)
    and, if ``graph`` is given, have exactly the edges of :func:`chain_graph` as its continuations and
    links."""
    if -1 in ids:
        return "Java leaves a gap without ID"
    if len(set(ids)) != len(ids):
        return "Java gives two gaps the same ID"
    prev_ids = [g.id for g in prev]
    if set(ids) - set(prev_ids) != set(new_ids):
        return "Java's new IDs are inconsistent"
    pos = {gid: k for k, gid in enumerate(ids)}
    pairs = set()
    for j, gid in enumerate(prev_ids):
        if gid in pos:
            if j in links:
                return "a linked gap survives"
            pairs.add((j, pos[gid]))
        for c in links.get(j, ()):
            if c not in pos:  # pragma: no cover - links only point to new IDs of this sample
                return "a link points nowhere"
            pairs.add((j, pos[c]))
    for ja, kb in _components(len(prev), len(ids), [p for p in pairs if prev_ids[p[0]] != ids[p[1]]]):
        if ja and kb and (len(ja), len(kb)) not in ((1, 2), (2, 1)):
            return f"Java links {len(ja)} gap(s) to {len(kb)}"
    if graph is not None and pairs != graph:
        return "Java's step disagrees with the hidden chains"
    return None


def _get_gaps_fixed(poly: Polygon, path: Path, samples: PathSamples, eps: float, java_matching: bool) -> GapHistory:
    """``compat="fixed"``: Java's rules where they are consistent, hidden-chain matching elsewhere.

    Every step ``i-1 -> i`` between consecutive samples is first done Java's way (the event branch at
    a general-inflection / bitangent crossing, the ``samePhysicalGap`` shift otherwise) on copies.
    The result is kept if it is consistent (:func:`_proposal_problem`; with ``java_matching=False`` it
    must also agree with :func:`chain_graph`).  Otherwise the step is read off the chain graph: each
    connected component with ``p`` previous and ``q`` new gaps is a continuation (1-1, the ID carries
    over), a disappearance (1-0), an appearance (0-1), a split (1-q), a merge (p-1) or, for ``p, q >=
    2``, a merge followed by a split (every previous gap linked to every new one).  New IDs are handed
    out in boundary order.  A gap set is closed before every step that has a component event.
    """
    if samples.errors:
        raise ValueError("compat='fixed' needs samples computed with compat='fixed'")
    n = poly.n
    edges = [poly.edge(i) for i in range(n)]
    cps = samples.critical_points
    length = path.length()
    out = GapHistory()
    out.samples = samples
    out.compat = "fixed"
    last_id = 0
    first = [Gap(p) for p in samples.physical[0]]  # type: ignore[union-attr]
    for g in first:
        last_id += 1
        g.id = last_id
    group: List[List[Gap]] = [first]
    group_start = 0
    sample_gaps: List[List[Gap]] = [first]

    def close(end: int) -> None:
        last = group[-1]
        ce = _iedges(group, n)
        for g in last:
            g.iedge = ce[g.id]
        t = _relative_time(path, cps[group_start], length, "fixed", eps)
        if last:
            last[0].relative_time = t
        out.append(last)
        out.times.append(t)
        out.sample_indices.append(end)

    for i in range(1, len(cps)):
        cp = cps[i]
        prev = group[-1]
        phys = samples.physical[i]
        assert phys is not None
        graph = chain_graph(poly, [g.physical for g in prev], phys)  # type: ignore[misc]
        prop = _java_proposal(cp, i, samples, prev, n, edges, last_id)
        if prop is None:
            problem: Optional[str] = "Java throws"
        else:
            ids, links, new_last = prop
            problem = _proposal_problem(prev, ids, links, range(last_id + 1, new_last + 1),
                                        None if java_matching else set(graph))
        if problem is not None:
            ids, links, new_last = _from_graph(prev, len(phys), graph, last_id)
            out.inferred.append((i, f"{problem}; {_describe(prev, ids, links, last_id)}"))
        cur = [Gap(p, gid) for p, gid in zip(phys, ids)]
        event = new_last > last_id or bool(set(g.id for g in prev) - set(ids))
        if event:
            for j, children in links.items():
                for c in children:
                    prev[j].add_gap(c)
            close(i - 1)
            group = [cur]
            group_start = i
        else:
            group.append(cur)
        last_id = new_last
        sample_gaps.append(cur)
    close(len(cps) - 1)
    out.sample_ids = [[g.id for g in gs] for gs in sample_gaps]
    return out


def _from_graph(prev: List[Gap], nb: int, graph: Sequence[Tuple[int, int]], last_id: int):
    """IDs and links of a step read off the chain graph (see :func:`_get_gaps_fixed`)."""
    ids = [-1] * nb
    comps = _components(len(prev), nb, graph)
    for ja, kb in comps:
        if len(ja) == 1 and len(kb) == 1:  # continuation
            ids[kb[0]] = prev[ja[0]].id
    for k in range(nb):  # new IDs in boundary order
        if ids[k] == -1:
            last_id += 1
            ids[k] = last_id
    links: Dict[int, List[int]] = {}
    for ja, kb in comps:
        if ja and kb and not (len(ja) == 1 and len(kb) == 1):
            for j in ja:
                links[j] = [ids[k] for k in kb]
    return ids, links, last_id


def _describe(prev: List[Gap], ids: List[int], links: Dict[int, List[int]], last_id: int) -> str:
    prev_ids = [g.id for g in prev]
    linked = {c for cs in links.values() for c in cs}
    parts = [f"{prev_ids[j]} -> {cs}" for j, cs in links.items()]
    parts += [f"appear {g}" for g in ids if g > last_id and g not in linked]
    parts += [f"disappear {g}" for j, g in enumerate(prev_ids) if g not in ids and j not in links]
    return "from the hidden chains: " + (", ".join(parts) if parts else "no event, IDs carried over")
