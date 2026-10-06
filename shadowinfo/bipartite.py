"""Incremental bipartite I-state (T-RO 2012, Sec. V-B, Fig. 9).

Left vertices are the shadows at ``t0`` and the appeared shadows; right
vertices are the disappeared shadows and the shadows alive now.  An edge
``L_i -> R_j`` means a target hidden in ``s_i`` may end up in ``s_j``.  Each
alive shadow keeps its *reach set* (the left vertices it is connected to):

* appear ``s``      -> ``reach[s] = {L_s}``                         (Fig. 9(a))
* split ``s -> a, b`` -> both children copy ``reach[s]``             (Fig. 9(b))
* merge ``a, b -> s`` -> ``reach[s] = reach[a] | reach[b]``          (Fig. 9(c))
* disappear ``s``   -> right vertex ``R_s`` with edges from ``reach[s]`` (Fig. 9(d))

FOV events (Sec. V-D) add pseudo-shadows: an enter is an appear of ``k`` targets
merged into ``s`` (a left :class:`Pseudo` vertex), an exit is a split of ``s``
whose ``k``-target part disappears (a right :class:`Pseudo` vertex).  With
``fov="batch"`` the events are accumulated per shadow by a
:class:`~shadowinfo.fov.FovBatch` and flushed before the shadow takes part in a
component event or a query.

Left and right vertices live in separate namespaces: a shadow that appears and
is still alive (``s18`` in Fig. 11) is both ``left[18]`` and an alive right
vertex ``18``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Tuple, Union

from .events import (Appear, Bound, Disappear, Enter, Event, Exit, InvalidSequenceError, Merge,
                     ShadowSequence, Split)
from .fov import FovBatch


@dataclass(frozen=True)
class Pseudo:
    """Synthetic shadow created by an FOV event; never equal to an integer label.

    ``kind`` is ``"enter"`` (left vertex) or ``"exit"`` (right vertex), ``s`` the
    shadow the event happened to and ``n`` a running counter of pseudo-shadows.
    """

    kind: str
    s: int
    n: int

    def __str__(self) -> str:
        return f"{'+' if self.kind == 'enter' else '-'}{self.s}#{self.n}"


Vertex = Union[int, Pseudo]


def vertex_name(v: Vertex) -> str:
    """JSON-friendly name: ``"7"`` for shadow 7, ``"+7#3"``/``"-7#3"`` for pseudo-shadows."""
    return str(v)


class BipartiteIState:
    """Bipartite I-state maintained one event at a time (Fig. 9).

    ``left`` and ``disappeared`` map vertices to ``(lo, hi)`` bounds in creation
    order; ``reach`` maps each alive shadow to its frozenset of left vertices and
    ``disappeared_reach`` does the same for disappeared right vertices.
    """

    def __init__(self, initial: Dict[int, Bound], fov: str = "batch") -> None:
        if fov not in ("batch", "naive"):
            raise ValueError("fov must be 'batch' or 'naive'")
        self.fov = fov
        self.left: Dict[Vertex, Bound] = {}
        self.left_order: Dict[Vertex, int] = {}
        self.initial: List[int] = []
        self.disappeared: Dict[Vertex, Bound] = {}
        self.disappeared_reach: Dict[Vertex, FrozenSet[Vertex]] = {}
        self.reach: Dict[int, FrozenSet[Vertex]] = {}
        self.pending: Dict[int, FovBatch] = {}
        self.seen: set = set()
        self.n_pseudo = 0
        for s, (lo, hi) in initial.items():
            self._appear(int(s), lo, hi)
            self.initial.append(int(s))

    @classmethod
    def from_sequence(cls, seq: ShadowSequence, fov: str = "batch") -> "BipartiteIState":
        st = cls(seq.initial, fov)
        for e in seq.events:
            st.apply(e)
        return st

    # -- events ---------------------------------------------------------------
    def apply(self, e: Event) -> None:
        if isinstance(e, (Enter, Exit)):
            self._require_alive(e.s, e)
            if e.k < 1:
                raise InvalidSequenceError(f"{e}: k must be >= 1")
            if self.fov == "naive":
                (self._enter if isinstance(e, Enter) else self._exit)(e.s, e.k)
            else:
                self.pending.setdefault(e.s, FovBatch()).add(e)
        elif isinstance(e, Appear):
            self._appear(e.s, e.lo, e.hi)
        elif isinstance(e, Disappear):
            self._require_alive(e.s, e)
            _check_bound(e.lo, e.hi, e)
            self.flush(e.s)
            self.disappeared[e.s] = (e.lo, e.hi)
            self.disappeared_reach[e.s] = self.reach.pop(e.s)
        elif isinstance(e, Split):
            if e.a == e.b:
                raise InvalidSequenceError(f"{e}: repeated label")
            self._require_alive(e.s, e)
            self._require_fresh((e.a, e.b), e)
            self.flush(e.s)
            r = self.reach.pop(e.s)
            self.reach[e.a] = self.reach[e.b] = r
            self.seen.update((e.a, e.b))
        elif isinstance(e, Merge):
            if e.a == e.b:
                raise InvalidSequenceError(f"{e}: repeated label")
            self._require_alive(e.a, e)
            self._require_alive(e.b, e)
            self._require_fresh((e.s,), e)
            self.flush(e.a)
            self.flush(e.b)
            self.reach[e.s] = self.reach.pop(e.a) | self.reach.pop(e.b)
            self.seen.add(e.s)
        else:
            raise TypeError(f"unknown event {e!r}")

    def extend(self, events: Iterable[Event]) -> None:
        for e in events:
            self.apply(e)

    def flush(self, s: Union[int, None] = None) -> None:
        """Turn pending FOV batches into pseudo-shadows (for ``s``, or all shadows)."""
        for t in ([s] if s is not None else list(self.pending)):
            b = self.pending.pop(t, None)
            if b is None:
                continue
            for e in b.batch_events(t):
                (self._enter if isinstance(e, Enter) else self._exit)(t, e.k)

    def _appear(self, s: int, lo: int, hi: float) -> None:
        _check_bound(lo, hi, s)
        self._require_fresh((s,), s)
        self.seen.add(s)
        self._add_left(s, (lo, hi))
        self.reach[s] = frozenset((s,))

    def _enter(self, s: int, k: int) -> None:
        p = self._pseudo("enter", s)
        self._add_left(p, (k, k))
        self.reach[s] = self.reach[s] | {p}

    def _exit(self, s: int, k: int) -> None:
        p = self._pseudo("exit", s)
        self.disappeared[p] = (k, k)
        self.disappeared_reach[p] = self.reach[s]

    def _pseudo(self, kind: str, s: int) -> Pseudo:
        self.n_pseudo += 1
        return Pseudo(kind, s, self.n_pseudo)

    def _add_left(self, v: Vertex, b: Bound) -> None:
        self.left_order[v] = len(self.left_order)
        self.left[v] = b

    def _require_alive(self, s: int, e: object) -> None:
        if s not in self.reach:
            raise InvalidSequenceError(f"{e}: shadow {s} not alive")

    def _require_fresh(self, labels: Tuple[int, ...], e: object) -> None:
        for s in labels:
            if s in self.seen:
                raise InvalidSequenceError(f"{e}: label {s} reused")

    # -- queries ----------------------------------------------------------------
    def alive(self) -> List[int]:
        return list(self.reach)

    def right(self) -> Dict[Tuple[str, Vertex], Bound]:
        """Right vertices as ``{("gone", v) | ("alive", s): (lo, hi)}`` (alive ones are ``[0, inf)``).

        Pending FOV batches are not reflected; call :meth:`flush` first.
        """
        out: Dict[Tuple[str, Vertex], Bound] = {("gone", v): b for v, b in self.disappeared.items()}
        out.update({("alive", s): (0, math.inf) for s in self.reach})
        return out

    def edges(self) -> List[Tuple[Vertex, Tuple[str, Vertex]]]:
        """Edges ``(left, right)`` with right vertices tagged as in :meth:`right`."""
        out = []
        for tag, src in (("gone", self.disappeared_reach), ("alive", self.reach)):
            for v, r in src.items():
                out += [(u, (tag, v)) for u in self.sorted_left(r)]
        return out

    def sorted_left(self, vertices: Iterable[Vertex]) -> List[Vertex]:
        return sorted(vertices, key=self.left_order.__getitem__)

    def real_left(self) -> List[int]:
        return [v for v in self.left if not isinstance(v, Pseudo)]

    def to_dict(self) -> dict:
        """JSON-friendly snapshot (Fig. 11(c) style); flushes pending FOV batches."""
        self.flush()
        return {
            "left": [[vertex_name(v), lo, None if hi == math.inf else hi] for v, (lo, hi) in self.left.items()],
            "right": [[vertex_name(v), lo, None if hi == math.inf else hi, tag]
                      for (tag, v), (lo, hi) in self.right().items()],
            "edges": [[vertex_name(u), vertex_name(v), tag] for u, (tag, v) in self.edges()],
        }


def _check_bound(lo: int, hi: float, where: object) -> None:
    if lo < 0 or hi < lo:
        raise InvalidSequenceError(f"{where}: bad bound {(lo, hi)}")


__all__ = ["BipartiteIState", "Pseudo", "Vertex", "vertex_name"]
