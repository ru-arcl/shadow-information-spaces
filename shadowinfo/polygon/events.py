"""Gap histories as shadow component events, and the never-see-evader scenario.

:func:`gap_history_to_events` reads the component events (T-RO 2012 Sec. II-A: appear, disappear,
split, merge) off a gap history from :func:`~shadowinfo.polygon.gaps.get_gaps`, as the existing
:mod:`shadowinfo.events` types, so the polygon front end feeds every filter of the package unchanged
(combinatorial, LP, probabilistic).  This is how the Java oracle and filter read ``Gap[][]``
(``docs/notes/original_java.md`` §2.7, last paragraph; ``SingleTypeAgentOracle.distributeAgents``).

:func:`never_see_evader` ports ``pe.Algorithm.processGapHistoryInformation(gapss, NEVER_SEE_EVADER)``
(§2.8; pursuit-evasion as a special case, T-RO Sec. V-E): every shadow at the start is contaminated,
a new shadow is clear, continuing and split shadows keep their label and a merged shadow is
contaminated if any part was.  :func:`never_see_evader_via_filter` derives the same labels from the
counting filter (:func:`shadowinfo.nondeterministic.counting_filter`): a shadow is clear iff its
upper bound is 0.
"""

from __future__ import annotations

import math
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from ..events import Appear, Bound, Disappear, Event, Merge, ShadowSequence, Split
from ..nondeterministic import counting_filter
from .gaps import Gap, GapTrackingError, _check_compat

__all__ = [
    "NSEState", "GapEvents", "transition_events", "gap_history_to_events", "never_see_evader",
    "never_see_evader_via_filter", "nse_string",
]


class NSEState:
    """``pe.is.basic.NSEState``: a gap is clear (no evader can be in it) or contaminated.

    Two singletons, :attr:`CLEAR` and :attr:`CONTAMINATED`; ``str()`` is Java's ``"| 0"`` / ``"| 1"``,
    which ``Gap.toString`` prints after the iEdge.
    """

    CLEAR: "NSEState"
    CONTAMINATED: "NSEState"

    __slots__ = ("clear",)

    def __init__(self, clear: bool):
        self.clear = clear

    def is_clear(self) -> bool:
        return self.clear

    def __str__(self) -> str:
        return "| 0" if self.clear else "| 1"

    def __repr__(self) -> str:
        return "NSEState.CLEAR" if self.clear else "NSEState.CONTAMINATED"


NSEState.CLEAR = NSEState(True)
NSEState.CONTAMINATED = NSEState(False)


class GapEvents(NamedTuple):
    """Component events of a gap history: ``(initial, events, times)``.

    ``initial`` are the IDs of the first gap set (the shadows at ``t0``), ``events`` the
    :class:`~shadowinfo.events.Appear` / ``Disappear`` / ``Split`` / ``Merge`` events in order and
    ``times[i]`` the relative time (arc-length fraction) of ``events[i]``.  Appear and disappear carry
    the counts ``(0, 0)``: the gap history alone says nothing about targets.
    """

    initial: List[int]
    events: List[Event]
    times: List[float]

    def sequence(self, initial_bound: Bound = (0, math.inf)) -> ShadowSequence:
        """A :class:`~shadowinfo.events.ShadowSequence` with ``initial_bound`` for every initial shadow."""
        return ShadowSequence({s: initial_bound for s in self.initial}, list(self.events))


def _fresh_labels(history: Sequence[Sequence[Gap]]):
    nxt = max((g.id for gs in history for g in gs), default=0)
    while True:
        nxt += 1
        yield nxt


def transition_events(history: Sequence[Sequence[Gap]]) -> List[List[Event]]:
    """The events between gap sets ``k`` and ``k + 1``, for every ``k``.

    * A gap of set ``k`` linked to 2 IDs (``toGapSet``) splits: ``Split(s, a, b)`` with ``a, b`` in
      Java ``HashSet`` order (the order the Java oracle fills them).
    * Gaps of set ``k`` linked to the same single ID merge: ``Merge(a, b, s)``, parents in set order.
    * An ID of set ``k + 1`` that is neither in set ``k`` nor linked appears; an ID of set ``k`` that
      is not in set ``k + 1`` and not linked disappears.

    Java histories have exactly one event per transition.  ``compat="fixed"`` histories may have
    several (a missed crossing) and, after a repair, a component with ``p`` parents and ``q`` children
    beyond 1-2 / 2-1; it becomes a chain of binary merges then binary splits through fresh labels
    (above every ID of the history), which only these events use.  Events of one transition are
    independent; they are listed component by component in boundary order of set ``k``, then appears,
    then disappears.  Raises ``ValueError`` if the links are not a history (a linked gap survives, a
    link points to an ID missing from set ``k + 1``, or a gap is linked to a single ID alone).
    """
    fresh = _fresh_labels(history)
    out: List[List[Event]] = []
    for k in range(len(history) - 1):
        cur, nxt = history[k], history[k + 1]
        cur_ids = [g.id for g in cur]
        nxt_ids = [g.id for g in nxt]
        nxt_set = set(nxt_ids)
        # link components: union parents sharing a child
        comps: List[Tuple[List[int], List[int]]] = []
        child_comp: Dict[int, int] = {}
        for g in cur:
            to = list(g.to_gaps)
            if not to:
                continue
            if g.id in nxt_set:
                raise ValueError(f"set {k}: linked gap {g.id} survives into set {k + 1}")
            for c in to:
                if c not in nxt_set:
                    raise ValueError(f"set {k}: gap {g.id} links to {c}, which is not in set {k + 1}")
            hit = sorted({child_comp[c] for c in to if c in child_comp})
            if hit:
                ci = hit[0]
                for other in hit[1:]:  # join components
                    ps, cs = comps[other]
                    comps[ci][0].extend(ps)
                    comps[ci][1].extend(c for c in cs if c not in comps[ci][1])
                    for c in cs:
                        child_comp[c] = ci
                    comps[other] = ([], [])
            else:
                ci = len(comps)
                comps.append(([], []))
            comps[ci][0].append(g.id)
            for c in to:
                if c not in comps[ci][1]:
                    comps[ci][1].append(c)
                child_comp[c] = ci
        evs: List[Event] = []
        linked = set(child_comp)
        pos = {gid: j for j, gid in enumerate(cur_ids)}
        for parents, children in sorted((c for c in comps if c[0]), key=lambda c: min(pos[p] for p in c[0])):
            parents = sorted(parents, key=pos.__getitem__)
            if len(parents) == 1 and len(children) == 1:
                raise ValueError(f"set {k}: gap {parents[0]} is linked to {children[0]} alone (not an event)")
            s = parents[0]
            if len(children) == 1:  # merge (a chain of binary merges for more than 2 parents)
                for p in parents[1:-1]:
                    t = next(fresh)
                    evs.append(Merge(s, p, t))
                    s = t
                evs.append(Merge(s, parents[-1], children[0]))
                continue
            for p in parents[1:]:  # p > 1 parents, q > 1 children: merge first, then split
                t = next(fresh)
                evs.append(Merge(s, p, t))
                s = t
            for c in children[:-2]:
                t = next(fresh)
                evs.append(Split(s, c, t))
                s = t
            evs.append(Split(s, children[-2], children[-1]))
        cur_set = set(cur_ids)
        evs += [Appear(g, 0, 0) for g in nxt_ids if g not in cur_set and g not in linked]
        evs += [Disappear(g.id, 0, 0) for g in cur if not g.to_gaps and g.id not in nxt_set]
        out.append(evs)
    return out


def gap_history_to_events(history: Sequence[Sequence[Gap]]) -> GapEvents:
    """Component events of a gap history (:func:`transition_events`, flattened) with their times.

    The time of an event is the relative time of the set it leads to (``history.times`` when the
    history has them, else ``set[0].relative_time`` as Java stores it).
    """
    if not history:
        return GapEvents([], [], [])
    times = getattr(history, "times", None) or [gs[0].relative_time if gs else 0.0 for gs in history]
    events: List[Event] = []
    ev_times: List[float] = []
    for k, evs in enumerate(transition_events(history)):
        events += evs
        ev_times += [times[k + 1]] * len(evs)
    return GapEvents([g.id for g in history[0]], events, ev_times)


def never_see_evader(history: Sequence[Sequence[Gap]], compat: str = "fixed") -> List[List[NSEState]]:
    """``pe.Algorithm.pursuerNeverSeeEvader``: clear/contaminated label of every gap of every set.

    Sets ``gap.state`` on the history (as Java does, so that :func:`~shadowinfo.polygon.gaps.format_gap_set`
    prints ``[id, iEdge| 1]``, the ``ProjectPanel2/3`` output) and returns the states per set.

    ``compat="java"``: literal.  All gaps of set 0 are contaminated; a gap still without a state is
    set clear (and its links are not followed); otherwise its state goes to the same ID in the next set
    (no link), the merged gap becomes contaminated if this part is (one link), or both children copy
    it (two links).  Quirk B6: a gap linked to more than 2 IDs is ignored, so its children end up
    clear.  Java returns at once when set 0 is empty, leaving every state untouched (``None``).
    ``compat="fixed"`` (default) first clears every state, then labels a gap of set ``k + 1``
    contaminated iff it continues a contaminated gap or any gap linked to it is contaminated; this is
    the same on every Java-shaped history, and an empty set 0 (a start point that sees no gap, which
    ``get_gaps(compat="fixed")`` allows where Java throws) makes every later gap clear.
    """
    _check_compat(compat)
    if not history:
        return []
    if compat == "fixed":
        for gs in history:
            for g in gs:
                g.state = None
        for g in history[0]:
            g.state = NSEState.CONTAMINATED
        for k in range(len(history) - 1):
            cur, nxt = history[k], history[k + 1]
            dirty = {g.id for g in cur if not g.state.is_clear()}  # type: ignore[union-attr]
            for g in cur:
                if g.id in dirty:
                    dirty.update(g.to_gaps)
            for g in nxt:
                g.state = NSEState.CONTAMINATED if g.id in dirty else NSEState.CLEAR
        return [[g.state for g in gs] for gs in history]  # type: ignore[misc]
    if not history[0]:  # Java returns at once (gapss[0].length == 0) and leaves every state as it was
        return [[g.state for g in gs] for gs in history]  # type: ignore[misc]
    for g in history[0]:
        g.state = NSEState.CONTAMINATED
    for k, gaps in enumerate(history):
        nxt_map: Dict[int, Gap] = {g.id: g for g in history[k + 1]} if k + 1 < len(history) else {}
        for g in gaps:
            to = g.to_gaps
            if g.state is None:
                g.state = NSEState.CLEAR
            elif len(to) == 0:
                ng = nxt_map.get(g.id)
                if ng is not None:
                    ng.state = g.state
            elif len(to) == 1:
                if not g.state.is_clear():  # type: ignore[attr-defined]
                    _java_get(nxt_map, to[0], k).state = NSEState.CONTAMINATED
            elif len(to) == 2:
                _java_get(nxt_map, to[0], k).state = g.state
                _java_get(nxt_map, to[1], k).state = g.state
    return [[g.state for g in gs] for gs in history]  # type: ignore[misc]


def _java_get(m: Dict[int, Gap], gid: int, k: int) -> Gap:
    if gid not in m:
        raise GapTrackingError(f"pursuerNeverSeeEvader: gap {gid} linked from set {k} is not in set {k + 1} "
                               "(Java NullPointerException)", java_exception="java.lang.NullPointerException")
    return m[gid]


def never_see_evader_via_filter(history: Sequence[Sequence[Gap]]) -> List[List[NSEState]]:
    """The never-see-evader labels from the counting filter, as a cross-check of :func:`never_see_evader`.

    Every shadow at ``t0`` starts with ``[0, inf)`` targets (:func:`~shadowinfo.nondeterministic.counting_filter`),
    appearing and disappearing shadows hold 0 (``Appear(s, 0, 0)``, ``Disappear(s, 0, 0)``).  A gap is
    clear iff the upper bound of its shadow is 0.  (``pursuit_evasion_filter`` adds ``total = (1, 1)``,
    a different model.)
    """
    if not history:
        return []
    per = transition_events(history)
    filt = counting_filter([g.id for g in history[0]])
    out: List[List[NSEState]] = []
    for k, gs in enumerate(history):
        if k:
            filt.extend(per[k - 1])
        b = filt.all_bounds()
        out.append([NSEState.CLEAR if b[g.id][1] == 0 else NSEState.CONTAMINATED for g in gs])
    return out


def nse_string(states: Sequence[Optional[NSEState]]) -> str:
    """One character per gap: ``0`` clear, ``1`` contaminated, ``-`` no state (golden-fixture format)."""
    return "".join("-" if s is None else "0" if s.is_clear() else "1" for s in states)
