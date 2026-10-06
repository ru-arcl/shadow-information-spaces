"""The single-type-agent filter of the original Java: equations, bipartite I-state and bounds.

Port of ``pe.is.singleTypeAgent.SingleTypeAgentAlgorithm``, ``pe.is.Equation``,
``graph.BipartiteGraph`` / ``Vertex`` / ``Edge`` and ``graph.algorithm.maxflow.MaxFlow``
(``docs/notes/original_java.md`` §2.10).  Given a gap history and the observations of a
:class:`~shadowinfo.polygon.oracle.SingleTypeAgentOracle`:

* :func:`derive_gap_evolving_equations` -- the linear equations of T-RO 2012 Sec. IV-B (Eq. 6 style),
  printed by :class:`Equation` exactly as Java does (``-1 x1 + 1 x6 + 1 x7 = 0``);
* :func:`derive_shadow_info_state` -- the bipartite I-state of Sec. V-B (Fig. 9) in Java's shape: all
  initial shadows share ONE left vertex ``0`` holding the known hidden total (B12), an appeared shadow
  is a single vertex that is its own left and right vertex (a self-loop) until it splits or merges;
* :func:`java_bounds` -- ``deriveShadowBounds``: the max-flow recipe of Sec. V-C as coded, on a literal
  port of the Java Edmonds-Karp (``MaxFlow``).  Its result depends on ``HashMap<Vertex, ...>``
  iteration order, i.e. on JVM identity hashes, and is not a correct bound computation (B10, B11):
  kept for comparison only, with a documented deterministic order;
* :func:`exact_bounds` -- the correct answer: :class:`~shadowinfo.nondeterministic.CombinatorialFilter`
  with the oracle's hidden total, for every shadow alive at the end.  It raises
  :class:`~shadowinfo.maxflow.InfeasibleError` when the observations are inconsistent (as the merge
  bug B7 of the original oracle makes them, e.g. on all five recorded T-RO Fig. 15(b) runs).

The graph built here equals :class:`shadowinfo.bipartite.BipartiteIState` once Java's pooling of the
initial shadows and its self-loop vertices are allowed for (``tools/java_reference/convert.py``,
``compare_bipartite``); the exact filter uses the latter.
"""

from __future__ import annotations

import math
import random as _random
import time as _time
from collections import deque
from typing import Callable, Dict, Iterator, List, Optional, Sequence, Tuple, Union

from ..events import ShadowSequence
from ..maxflow import InfeasibleError
from ..nondeterministic import CombinatorialFilter
from .gaps import Gap, format_gap_sets
from .oracle import EventType, OracleError, SingleTypeAgentEvent, SingleTypeAgentOracle, java_int

__all__ = [
    "Equation", "derive_gap_evolving_equations", "BGVertex", "BGEdge", "JavaBipartiteGraph",
    "derive_shadow_info_state", "JavaBounds", "java_bounds", "exact_bounds", "bounds_lines", "project_panel_lines",
]

_NPE = "java.lang.NullPointerException"


# --------------------------------------------------------------------------- equations


class Equation:
    """``pe.is.Equation``: ``sum(coeff * x_id) = constant`` over shadow IDs.

    Terms are kept sorted by ID (``TreeMap``); :meth:`add_term` overwrites an existing coefficient.
    ``str()`` is Java's ``toString``: ``"c xN"`` joined by ``" + "`` when the next coefficient is
    ``>= 0`` and by ``" "`` otherwise, then ``" = constant"`` (``-1 x4 -1 x10 + 1 x11 = 0``).
    """

    __slots__ = ("coeffs", "constant")

    def __init__(self, coeffs: Optional[Dict[int, int]] = None, constant: int = 0):
        self.coeffs: Dict[int, int] = dict(coeffs or {})
        self.constant = constant

    def add_term(self, unknown: int, coeff: int) -> None:
        """``addTerm``: set the coefficient of ``x_unknown`` (overwrites)."""
        self.coeffs[unknown] = coeff

    def terms(self) -> List[Tuple[int, int]]:
        """``[(id, coeff)]`` sorted by ID."""
        return sorted(self.coeffs.items())

    def __str__(self) -> str:
        ts = self.terms()
        buf = ""
        for i, (k, c) in enumerate(ts):
            buf += f"{c} x{k}"
            if i != len(ts) - 1:
                buf += " + " if ts[i + 1][1] >= 0 else " "
        return buf + f" = {self.constant}"

    def __repr__(self) -> str:
        return f"Equation({self.terms()!r}, {self.constant})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Equation) and self.coeffs == other.coeffs and self.constant == other.constant

    def evaluate(self, values: Dict[int, int]) -> int:
        """Left-hand side minus constant for the given unknowns (0 iff the equation holds)."""
        return sum(c * values[k] for k, c in self.coeffs.items()) - self.constant


def _event_equation(e: SingleTypeAgentEvent) -> Equation:
    eq = Equation()
    t = e.type
    try:
        if t is EventType.APPEAR:
            eq.add_term(e.to_gaps[0], 1)  # type: ignore[index]
            eq.constant = e.moved
        elif t is EventType.DISAPPEAR:
            eq.add_term(e.from_gaps[0], 1)  # type: ignore[index]
            eq.constant = e.moved
        elif t is EventType.SPLIT:
            eq.add_term(e.to_gaps[0], 1)  # type: ignore[index]
            eq.add_term(e.to_gaps[1], 1)  # type: ignore[index]
            eq.add_term(e.from_gaps[0], -1)  # type: ignore[index]
        elif t is EventType.MERGE:
            eq.add_term(e.to_gaps[0], 1)  # type: ignore[index]
            eq.add_term(e.from_gaps[0], -1)  # type: ignore[index]
            eq.add_term(e.from_gaps[1], -1)  # type: ignore[index]
        elif t is None:
            raise OracleError(f"deriveGapEvolvingEquations: event at t={e.time} has no type (B9)",
                              java_exception=_NPE)
    except TypeError as ex:
        raise OracleError(f"deriveGapEvolvingEquations: {ex}", java_exception=_NPE) from None
    return eq


def derive_gap_evolving_equations(history: Optional[Sequence[Sequence[Gap]]],
                                  oracle: SingleTypeAgentOracle) -> List[Equation]:
    """``SingleTypeAgentAlgorithm.deriveGapEvolvingEquations(gapss, oracle)`` (T-RO Sec. IV-B).

    First ``sum_{g initial} x_g = N - visible0``; then one equation per event: APPEAR / DISAPPEAR
    ``x_s = moved``, SPLIT ``x_a + x_b - x_s = 0``, MERGE ``x_s - x_a - x_b = 0``.

    ``oracle.compat == "java"``: literal -- one equation per gap set, with the event looked up by the
    set's relative time (``getEventByTime(gapss[i][0].getRelativeTime())``); a set without an event
    type raises :class:`~shadowinfo.polygon.oracle.OracleError` (Java NPE, B9), INIT and FOV events give
    the empty equation ``" = 0"``.  ``"fixed"``: one equation per oracle event (``history`` unused), which
    is the same list on Java-shaped histories.
    """
    n = oracle.get_total_agent_number()
    if oracle.compat == "fixed":
        if not oracle.events:
            return []
        init = oracle.events[0]
        out = [Equation({g: 1 for g in init.to_gaps or []}, n - init.visible)]
        return out + [_event_equation(e) for e in oracle.events[1:]]
    if history is None:
        raise ValueError("compat='java' needs the gap history (one equation per gap set)")
    out = []
    for i, gaps in enumerate(history):
        if not gaps:
            raise OracleError(f"deriveGapEvolvingEquations: gap set {i} is empty",
                              java_exception="java.lang.ArrayIndexOutOfBoundsException")
        e = oracle.event_by_time(gaps[0].relative_time)
        if e is None:
            raise OracleError(f"deriveGapEvolvingEquations: no event at t={gaps[0].relative_time}",
                              java_exception=_NPE)
        if i == 0:
            eq = Equation()
            for g in gaps:
                eq.add_term(g.id, 1)
            eq.constant = java_int(n - e.visible)
        else:
            eq = _event_equation(e)
        out.append(eq)
    return out


# --------------------------------------------------------------------------- bipartite graph


class BGVertex:
    """``graph.Vertex``: ``id``, ``min_weight`` / ``max_weight`` (``-1`` = unknown, an alive shadow) and the
    out-edges ``ve`` (``veMap``: target vertex -> :class:`BGEdge`)."""

    __slots__ = ("id", "min_weight", "max_weight", "ve")

    def __init__(self, id: int):
        self.id = id
        self.min_weight = -1
        self.max_weight = -1
        self.ve: Dict["BGVertex", "BGEdge"] = {}

    def __repr__(self) -> str:
        return f"BGVertex({self.id}, W={self.max_weight})"


class BGEdge:
    """``graph.Edge``: ``from_v -> to_v`` with ``weight`` (flow) and ``capacity``."""

    __slots__ = ("from_v", "to_v", "weight", "capacity")

    def __init__(self, from_v: BGVertex, to_v: BGVertex, weight: int = 0, capacity: int = 0):
        self.from_v = from_v
        self.to_v = to_v
        self.weight = weight
        self.capacity = capacity


class JavaBipartiteGraph:
    """``graph.BipartiteGraph`` as built by ``deriveShadowInfoState``.

    * ``left`` / ``right`` -- ``leftVertexMap`` / ``rightVertexMap`` (``TreeMap<Integer, Vertex>``; iterate
      with :meth:`left_vertices` / :meth:`right_vertices`, sorted by ID);
    * ``right_edges`` -- ``rightVertexEdgeMap``: the in-edges of every right vertex.

    ``str()`` is ``BipartiteGraph.toString()`` except that each ``goes to:`` list is sorted by ID (Java
    prints it in ``HashMap`` identity-hash order).  :meth:`to_dict` is the golden fixtures' ``bg`` record.
    """

    def __init__(self) -> None:
        self.left: Dict[int, BGVertex] = {}
        self.right: Dict[int, BGVertex] = {}
        self.right_edges: Dict[BGVertex, List[BGEdge]] = {}

    def left_vertices(self) -> List[BGVertex]:
        return [self.left[k] for k in sorted(self.left)]

    def right_vertices(self) -> List[BGVertex]:
        return [self.right[k] for k in sorted(self.right)]

    def edges(self) -> List[Tuple[int, int]]:
        """``(left id, right id)`` of every edge, sorted."""
        return sorted((v.id, w.id) for v in self.left.values() for w in v.ve)

    def __str__(self) -> str:
        buf = ["Left vertices\n"]
        for v in self.left_vertices():
            buf.append(f"Vertex: {v.id} W: {v.max_weight}, goes to: ")
            buf.extend(f"{i} " for i in sorted(w.id for w in v.ve))
            buf.append("\n")
        buf.append("Right vertices\n")
        buf.extend(f"Vertex: {v.id} W: {v.max_weight}\n" for v in self.right_vertices())
        return "".join(buf)

    def to_dict(self) -> dict:
        """``{"n_left", "n_right", "n_edges", "left": [[id, min, max, sorted targets]], "right": [[id, min, max]]}``."""
        left = [[v.id, v.min_weight, v.max_weight, sorted(w.id for w in v.ve)] for v in self.left_vertices()]
        return {"n_left": len(self.left), "n_right": len(self.right), "n_edges": sum(len(x[3]) for x in left),
                "left": left, "right": [[v.id, v.min_weight, v.max_weight] for v in self.right_vertices()]}


def _need(v: Optional[BGVertex], what: str, e: SingleTypeAgentEvent) -> BGVertex:
    if v is None:
        raise OracleError(f"deriveShadowInfoState: {what} at t={e.time} is not a right vertex", java_exception=_NPE)
    return v


def derive_shadow_info_state(history: Optional[Sequence[Sequence[Gap]]],
                             oracle: SingleTypeAgentOracle) -> JavaBipartiteGraph:
    """``SingleTypeAgentAlgorithm.deriveShadowInfoState(gapss, oracle)``: the bipartite I-state (Fig. 9).

    Walks the oracle's events (``getEvent(i)``):

    * event 0 (if the first gap set is not empty): ONE left vertex ``0`` with weight ``N - visible0`` and an
      edge to a right vertex per initial shadow (B12; the per-shadow variant is commented out in Java);
    * APPEAR ``s <- k``: one vertex ``s`` of weight ``k`` that is both left and right, with a self-loop;
    * DISAPPEAR ``s -> k``: the right vertex ``s`` gets weight ``k`` (it stays, as a known count);
    * SPLIT ``s -> a b``: the right vertex ``s`` is replaced by ``a`` and ``b``; every in-edge is duplicated;
    * MERGE ``a b -> s``: the right vertices ``a``, ``b`` become ``s``; in-edges are deduplicated by source.

    ``history`` gives the initial gaps (``gapss[0]``); if ``None`` the INIT event's shadows are used.
    Raises :class:`~shadowinfo.polygon.oracle.OracleError` where Java throws (an untyped event, B9; a
    merge with one parent; an unknown shadow).  Java prints the graph; here ``str(graph)`` does.
    """
    bg = JavaBipartiteGraph()
    gaps0 = ([g.id for g in history[0]] if history is not None and len(history)
             else (list(oracle.events[0].to_gaps or []) if oracle.events else []))
    for i, e in enumerate(oracle.events):
        if i == 0 and gaps0:
            vl = BGVertex(0)
            vl.min_weight = vl.max_weight = java_int(oracle.get_total_agent_number() - e.visible)
            bg.left[vl.id] = vl
            for gid in gaps0:
                v = BGVertex(gid)
                ed = BGEdge(vl, v)
                vl.ve[v] = ed
                bg.right[v.id] = v
                bg.right_edges[v] = [ed]
            continue
        t = e.type
        if t is None:
            raise OracleError(f"deriveShadowInfoState: event at t={e.time} has no type (B9)", java_exception=_NPE)
        if t is EventType.APPEAR:
            v = BGVertex(e.to_gaps[0])  # type: ignore[index]
            v.min_weight = v.max_weight = e.moved
            ed = BGEdge(v, v)
            v.ve[v] = ed
            bg.right_edges[v] = [ed]
            bg.left[v.id] = v
            bg.right[v.id] = v
        elif t is EventType.DISAPPEAR:
            v = _need(bg.right.get(e.from_gaps[0]), f"disappearing shadow {e.from_gaps[0]}", e)  # type: ignore[index]
            v.max_weight = v.min_weight = e.moved
        elif t is EventType.SPLIT:
            frm = _need(bg.right.get(e.from_gaps[0]), f"splitting shadow {e.from_gaps[0]}", e)  # type: ignore[index]
            edges = bg.right_edges.get(frm)
            if edges is None:
                raise OracleError(f"deriveShadowInfoState: no edges into {frm.id}", java_exception=_NPE)
            v1, v2 = BGVertex(e.to_gaps[0]), BGVertex(e.to_gaps[1])  # type: ignore[index]
            del bg.right[frm.id]
            bg.right[v1.id] = v1
            bg.right[v2.id] = v2
            del bg.right_edges[frm]
            es1: List[BGEdge] = []
            es2: List[BGEdge] = []
            bg.right_edges[v1] = es1
            bg.right_edges[v2] = es2
            for ce in edges:
                ne1, ne2 = BGEdge(ce.from_v, v1), BGEdge(ce.from_v, v2)
                ce.from_v.ve.pop(frm, None)
                ce.from_v.ve[v1] = ne1
                ce.from_v.ve[v2] = ne2
                es1.append(ne1)
                es2.append(ne2)
        elif t is EventType.MERGE:
            a, b = e.from_gaps[0], e.from_gaps[1]  # type: ignore[index]
            v1 = _need(bg.right.get(a), f"merging shadow {a}", e)
            v2 = _need(bg.right.get(b), f"merging shadow {b}", e)
            del bg.right[v1.id]
            bg.right.pop(v2.id, None)
            v = BGVertex(e.to_gaps[0])  # type: ignore[index]
            bg.right[v.id] = v
            es1, es2 = bg.right_edges.pop(v1, None), bg.right_edges.pop(v2, None)
            if es1 is None or es2 is None:
                raise OracleError(f"deriveShadowInfoState: no edges into {a} or {b}", java_exception=_NPE)
            sources: Dict[BGVertex, None] = {}
            for ed in list(dict.fromkeys(es1 + es2)):
                sources[ed.from_v] = None
                ed.from_v.ve.pop(ed.to_v, None)
            es = []
            for s in sources:
                ed = BGEdge(s, v)
                s.ve[v] = ed
                es.append(ed)
            bg.right_edges[v] = es
        # INIT (with no initial gaps), AGENT_APPEAR, AGENT_DISAPPEAR: nothing (stubbed out in Java, B13)
    return bg


# --------------------------------------------------------------------------- Java max-flow


class _FVertex:
    """``graph.Vertex`` of the flow network: ``visited`` persists across calls, as in Java."""

    __slots__ = ("id", "rank", "visited", "max_weight", "ve", "vre")

    def __init__(self, id: int, rank: Tuple, order: "_Order"):
        self.id = id
        self.rank = rank
        self.visited = False
        self.max_weight = -1
        self.ve = _JMap(order)
        self.vre = _JMap(order)


class _FEdge:
    __slots__ = ("from_v", "to_v", "weight", "capacity")

    def __init__(self, from_v: _FVertex, to_v: _FVertex, weight: int = 0, capacity: int = 0):
        self.from_v = from_v
        self.to_v = to_v
        self.weight = weight
        self.capacity = capacity


class _Order:
    """Iteration order of ``HashMap<Vertex, ...>``.

    ``"sorted"`` / ``"reverse"``: by vertex ID.  ``"hash"`` (a seeded random identity hash per vertex):
    Java's ``HashMap`` order -- by bucket ``(h ^ h >>> 16) & (capacity - 1)``, then insertion order in the
    bucket, capacity 16 doubled when the size exceeds 0.75 of it (never shrinks; tree bins not modelled).
    """

    def __init__(self, order: Union[str, int, _random.Random]):
        if isinstance(order, str):
            if order not in ("sorted", "reverse"):
                raise ValueError(f"order must be 'sorted', 'reverse', an int seed or a random.Random, got {order!r}")
            self.kind = order
        elif isinstance(order, _random.Random):
            self.kind, self.rng = "hash", order
        elif isinstance(order, int) and not isinstance(order, bool):
            self.kind, self.rng = "hash", _random.Random(order)
        else:
            raise ValueError(f"order must be 'sorted', 'reverse', an int seed or a random.Random, got {order!r}")

    def rank(self, vid: int) -> Tuple:
        if self.kind == "sorted":
            return (vid,)
        if self.kind == "reverse":
            return (-vid,)
        h = self.rng.getrandbits(31)  # HotSpot identity hashes are 31-bit
        return (h ^ (h >> 16),)


class _JMap:
    """``HashMap<Vertex, Edge>`` with a modelled iteration order (see :class:`_Order`)."""

    __slots__ = ("d", "cap", "seq", "order")

    def __init__(self, order: _Order):
        self.d: Dict[_FVertex, Tuple[_FEdge, int]] = {}
        self.cap = 16
        self.seq = 0
        self.order = order

    def put(self, k: _FVertex, v: _FEdge) -> None:
        if k in self.d:
            self.d[k] = (v, self.d[k][1])
            return
        self.seq += 1
        self.d[k] = (v, self.seq)
        if len(self.d) > 0.75 * self.cap:
            self.cap *= 2

    def get(self, k: _FVertex) -> Optional[_FEdge]:
        x = self.d.get(k)
        return None if x is None else x[0]

    def remove(self, k: _FVertex) -> None:
        self.d.pop(k, None)

    def clear(self) -> None:
        self.d.clear()  # HashMap.clear keeps the table size

    def values(self) -> List[_FEdge]:
        if self.order.kind == "hash":
            m = self.cap - 1
            items = sorted(self.d.items(), key=lambda kv: (kv[0].rank[0] & m, kv[1][1]))
        else:
            items = sorted(self.d.items(), key=lambda kv: kv[0].rank)
        return [x[1][0] for x in items]


def _reset_visited(source: _FVertex) -> None:
    """``MaxFlow.resetGraphVisitedState``: BFS over ``veMap`` through *visited* vertices only."""
    source.visited = False
    q = deque([source])
    while q:
        fv = q.popleft()
        for e in fv.ve.values():
            if e.to_v.visited:
                e.to_v.visited = False
                q.append(e.to_v)


def _build_residue_graph(source: _FVertex) -> None:
    source.visited = True
    source.vre.clear()
    q = deque([source])
    while q:
        fv = q.popleft()
        fv.vre.clear()
        for e in fv.ve.values():
            to = e.to_v
            if not to.visited:
                to.visited = True
                q.append(to)
            fv.vre.put(to, _FEdge(e.from_v, e.to_v, e.capacity, 0))  # createResidueEdge
    _reset_visited(source)


def _search_for_path(source: _FVertex, sink: _FVertex) -> Optional[List[_FEdge]]:
    q: deque = deque([(source, None, None)])
    source.visited = True
    path = None
    while q:
        node = q.popleft()
        for e in node[0].vre.values():
            to = e.to_v
            if not to.visited:
                to.visited = True
                temp = (to, e, node)
                if to is sink:
                    path = temp
                    break
                q.append(temp)
        if path is not None:
            break
    _reset_visited(source)
    if path is None:
        return None
    out: List[_FEdge] = []
    while path is not None and path[2] is not None:
        out.insert(0, path[1])
        path = path[2]
    return out


def _search_all_augment_paths(source: _FVertex, sink: _FVertex) -> None:
    while True:
        p = _search_for_path(source, sink)
        if p is None:
            return
        m = 2147483647
        for e in p:
            if m > e.weight:
                m = e.weight
        for e in p:
            if m == e.weight:
                e.from_v.vre.remove(e.to_v)
            else:
                e.weight = java_int(e.weight - m)
            re = e.to_v.vre.get(e.from_v)
            if re is None:
                re = _FEdge(e.to_v, e.from_v, 0, 0)
                e.to_v.vre.put(e.from_v, re)
            re.weight = java_int(re.weight + m)


def _construct_flow_graph(source: _FVertex) -> None:
    q = deque([source])
    source.visited = True
    while q:
        fv = q.popleft()
        for e in fv.ve.values():
            to = e.to_v
            if not to.visited:
                to.visited = True
                q.append(to)
            re = to.vre.get(fv)
            e.weight = re.weight if re is not None else 0
    _reset_visited(source)


def _get_max_flow(source: _FVertex, sink: _FVertex) -> None:
    """``MaxFlow.getMaxFlow``: Edmonds-Karp on the residual maps; flows are written to ``edge.weight``.

    Literal quirks: zero-capacity residual edges are traversed and "augmented" by 0, a residual edge
    whose weight equals the bottleneck is removed, and ``visited`` flags are reset only along ``veMap``
    paths through visited vertices, so a vertex reached through a reverse residual edge can stay
    visited and be skipped by later searches.
    """
    _build_residue_graph(source)
    _search_all_augment_paths(source, sink)
    _construct_flow_graph(source)


class JavaBounds(Dict[int, Tuple[int, int]]):
    """``{shadow: (min, max)}`` printed by ``deriveShadowBounds``; ``lines`` are the printed lines."""

    lines: List[str]


def java_bounds(bg: JavaBipartiteGraph, order: Union[str, int, _random.Random] = "sorted") -> JavaBounds:
    """``SingleTypeAgentAlgorithm.deriveShadowBounds(bg)``: Java's per-shadow "bounds" (comparison only).

    1. Left vertices with a self-loop (appeared, never split or merged) are dropped (B11: no bound for
       them).  ``source -> L`` and every ``L -> R`` get capacity ``w_L``; ``max = sum w_L``.
    2. A disappeared right vertex gets ``R -> sink2`` of capacity ``w_R`` (``sure = sum w_R``), an alive one
       ``R -> sink1``; ``sink1 -> sink2`` has capacity ``max - sure``.
    3. For each alive shadow ``j`` in ID order: with only ``j -> sink1`` open (capacity ``max - sure``),
       ``Max flow in shadow j is f(j -> sink1)``; with ``j`` closed and the others open,
       ``Min flow in shadow j is max - (sure + sum_i f(i -> sink1))``.

    The flows come from a literal port of ``MaxFlow`` (see :func:`_get_max_flow`), whose result depends on
    ``HashMap<Vertex, ...>`` iteration order, i.e. on JVM identity hashes (B10): ``order="sorted"`` (by
    vertex ID; source ``-1``, sinks ``-2``, ``-3``), ``"reverse"``, or an int seed / ``random.Random`` for
    random identity hashes in Java's ``HashMap`` bucket order.  It equals the Java output only where the
    Java output does not depend on the order.  Lower bounds of disappeared shadows are not enforced and
    no feasibility check is made, so the numbers can be wrong; use :func:`exact_bounds`.

    ``bg`` is not modified (Java mutates it).  Returns a :class:`JavaBounds` dict with the printed lines in
    ``.lines`` (without the final ``println("\\n")``).
    """
    od = _Order(order)
    left = bg.left_vertices()
    right = bg.right_vertices()
    # identity hashes in order of creation: graph vertices (left then right, by ID), then source/sinks
    fv: Dict[BGVertex, _FVertex] = {}
    for v in left + [r for r in right if r not in left]:
        if v not in fv:
            fv[v] = _FVertex(v.id, od.rank(v.id), od)
            fv[v].max_weight = v.max_weight
    for v in left:  # veMap of the graph's left vertices (insertion by target ID)
        for w in sorted(v.ve, key=lambda x: x.id):
            if w not in fv:
                fv[w] = _FVertex(w.id, od.rank(w.id), od)
                fv[w].max_weight = w.max_weight
            fv[v].ve.put(fv[w], _FEdge(fv[v], fv[w]))
    source = _FVertex(-1, od.rank(-1), od)
    sink1 = _FVertex(-2, od.rank(-2), od)
    sink2 = _FVertex(-3, od.rank(-3), od)
    total = 0
    dropped = set()
    for v in left:
        f = fv[v]
        if v in v.ve:  # self-loop: an appeared shadow that never split or merged
            dropped.add(v)
            continue
        e = _FEdge(source, f, 0, v.max_weight)
        source.ve.put(f, e)
        total = java_int(total + e.capacity)
        for ed in f.ve.values():
            ed.capacity = e.capacity
    sure = 0
    current: List[_FEdge] = []
    for v in right:
        if v in dropped:
            continue
        f = fv[v]
        if v.max_weight == -1:
            e = _FEdge(f, sink1, 0, 0)
            f.ve.put(sink1, e)
            current.append(e)
        else:
            e = _FEdge(f, sink2, 0, 0)
            f.ve.put(sink2, e)
            e.capacity = v.max_weight
            sure = java_int(sure + v.max_weight)
    sink1.ve.put(sink2, _FEdge(sink1, sink2, 0, java_int(total - sure)))
    free = java_int(total - sure)
    out = JavaBounds()
    lines: List[str] = []
    for j, ej in enumerate(current):
        for i, ei in enumerate(current):
            ei.capacity = free if i == j else 0
        _get_max_flow(source, sink2)
        hi = ej.weight
        lines.append(f"Max flow in shadow {ej.from_v.id} is {hi}")
        for i, ei in enumerate(current):
            ei.capacity = 0 if i == j else free
        _get_max_flow(source, sink2)
        subtotal = sure
        for ei in current:
            subtotal = java_int(subtotal + ei.weight)
        lo = java_int(total - subtotal)
        lines.append(f"Min flow in shadow {ej.from_v.id} is {lo}")
        out[ej.from_v.id] = (lo, hi)
    out.lines = lines
    return out


# --------------------------------------------------------------------------- exact bounds


def exact_bounds(source: Union[SingleTypeAgentOracle, ShadowSequence], hidden_total: Optional[int] = None,
                 ) -> Dict[int, Tuple[int, float]]:
    """Exact bounds on every shadow alive at the end: ``{s: (lo, hi)}`` (T-RO Sec. V-C, done right).

    ``source`` is an initialised oracle (its :meth:`~shadowinfo.polygon.oracle.SingleTypeAgentOracle.sequence`
    and ``hidden_total`` are used) or a :class:`~shadowinfo.events.ShadowSequence` together with
    ``hidden_total``.  Computed by :class:`~shadowinfo.nondeterministic.CombinatorialFilter` with
    ``total=(hidden_total, hidden_total)``: every target distribution consistent with all observations is
    a feasible circulation, so the bounds are tight.  Unlike :func:`java_bounds` it also bounds appeared
    shadows (B11 fixed) and raises :class:`~shadowinfo.maxflow.InfeasibleError` when the observations are
    inconsistent.
    """
    if isinstance(source, SingleTypeAgentOracle):
        seq = source.sequence()
        h = source.hidden_total if hidden_total is None else hidden_total
    else:
        if hidden_total is None:
            raise ValueError("hidden_total is required with a ShadowSequence")
        seq, h = source, hidden_total
    f = CombinatorialFilter.from_sequence(seq, total=(h, h))
    return {s: (int(lo), hi if hi == math.inf else int(hi)) for s, (lo, hi) in sorted(f.all_bounds().items())}


def bounds_lines(bounds: Dict[int, Tuple[int, float]]) -> List[str]:
    """Bounds in the ``deriveShadowBounds`` print format, by shadow ID: ``Max flow in shadow s is hi`` then
    ``Min flow in shadow s is lo`` (``inf`` prints as ``Infinity``)."""
    out: List[str] = []
    for s in sorted(bounds):
        lo, hi = bounds[s]
        out.append(f"Max flow in shadow {s} is {'Infinity' if hi == math.inf else int(hi)}")
        out.append(f"Min flow in shadow {s} is {int(lo)}")
    return out


# --------------------------------------------------------------------------- the panel output


def project_panel_lines(history: Sequence[Sequence[Gap]], oracle: SingleTypeAgentOracle,
                        order: Union[str, int, _random.Random] = "sorted",
                        clock: Callable[[], float] = _time.perf_counter) -> Iterator[str]:
    """What the ``ProjectPanel`` / ``ProjectPanel5`` constructors print to stdout (§2.11), line by line.

    1. every gap set: ``relativeTime`` and ``Gap.toString() + " "`` per gap; two empty lines;
    2. ``oracle.initialize(history)``; equations, bipartite graph and bounds are computed;
    3. the bipartite graph (``str(graph)``; Java's ``println`` adds an empty line);
    4. the bounds: :func:`java_bounds` (``order``) when ``oracle.compat == "java"``, else the exact bounds of
       every alive shadow (:func:`exact_bounds`, same line format; a line ``Infeasible observations: ...``
       if there are none); then two empty lines (``println("\\n")``);
    5. ``Computation time: <ms>`` (steps 2-4 without printing, measured with ``clock``);
    6. every equation; an empty line; every event's ``toString()``.

    A generator: when a step raises (``compat="java"`` reproduces Java's exceptions), the lines printed
    before it have been yielded, as on the Java console.
    """
    yield from format_gap_sets(history, True)
    yield ""
    yield ""
    oracle.initialize(history)
    t0 = clock()
    eqs = derive_gap_evolving_equations(history, oracle)
    bg = derive_shadow_info_state(history, oracle)
    elapsed = clock() - t0
    yield from str(bg).split("\n")  # toString ends with a newline, println adds the empty line
    t0 = clock()
    if oracle.compat == "java":
        blines = java_bounds(bg, order).lines
    else:
        try:
            blines = bounds_lines(exact_bounds(oracle))
        except InfeasibleError as ex:
            blines = [f"Infeasible observations: {ex}"]
    elapsed += clock() - t0
    yield from blines
    yield ""
    yield ""
    yield f"Computation time: {int(elapsed * 1000)}"
    for eq in eqs:
        yield str(eq)
    yield ""
    for e in oracle.events:
        es = str(e)
        if es:
            yield es
