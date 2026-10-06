"""Combinatorial filter for nondeterministically moving targets (T-RO 2012, Sec. V).

The filter keeps the bipartite I-state (:mod:`shadowinfo.bipartite`) and answers
queries with flows on the augmented graph of Fig. 11(d)::

    T -> S          [0, inf)                    closes the circulation
    S -> I          [N_lo, N_hi] or [0, inf)    total targets in the shadows at t0
    I -> L_i        [l_i, u_i]                  shadow i at t0
    S -> L_i        [l_i, u_i]                  appeared shadow / FOV enter pseudo-shadow
    L_i -> R_j      [0, inf)                    bipartite edges
    R_j -> T        [l_j, u_j]                  disappeared shadow / FOV exit pseudo-shadow
    R_j -> T        [0, inf)                    shadow alive now

Every feasible circulation is a target distribution consistent with all
observations and vice versa (Prop. 5), so the bounds on a set ``Q`` of alive
shadows are the min/max of the flow into ``Q`` over feasible circulations.
These are computed exactly with :class:`~shadowinfo.maxflow.FlowNetwork`
(flow with lower bounds, then augmentation around the query edge).  The
paper's own recipe, Eqs. (7)/(8), is available as :meth:`CombinatorialFilter.bounds_paper`.

Inconsistent observations make the circulation infeasible; queries then raise
:class:`~shadowinfo.maxflow.InfeasibleError`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Hashable, Iterable, List, Mapping, Optional, Sequence, Tuple

from .bipartite import BipartiteIState, Vertex
from .events import Appear, Bound, Disappear, Event, Merge, ShadowSequence, Split
from .maxflow import INF, FlowNetwork, InfeasibleError

RightKey = Tuple[str, Vertex]  # ("gone", v) or ("alive", s)

_S, _T, _I, _G = ("S",), ("T",), ("I",), ("G",)


@dataclass
class _Net:
    net: FlowNetwork
    supply: Dict[Vertex, int]
    demand: Dict[RightKey, int]
    inner: Dict[Tuple[Vertex, RightKey], int]
    total: int
    group: int = -1


@dataclass
class Witness:
    """A feasible target distribution attaining a bound.

    ``supply[v]`` targets start in left vertex ``v`` and ``flow[(v, r)]`` of them
    end in right vertex ``r`` (``("gone", v)`` or ``("alive", s)``).
    """

    value: int
    supply: Dict[Vertex, int]
    flow: Dict[Tuple[Vertex, RightKey], int]


class CombinatorialFilter:
    """Incremental combinatorial filter over the bipartite I-state (Sec. V-B..E).

    ``initial`` gives ``(lo, hi)`` targets for each shadow at ``t0`` (``hi`` may
    be ``math.inf``); ``total=(lo, hi)`` optionally bounds the total number of
    targets in those shadows.  ``fov`` selects the naive or batched (Sec. V-D)
    treatment of FOV events; both give identical bounds.
    """

    def __init__(self, initial: Mapping[int, Bound], total: Optional[Bound] = None,
                 fov: str = "batch") -> None:
        self.istate = BipartiteIState(dict(initial), fov)
        self.total = total
        self._cache: Optional[_Net] = None

    @classmethod
    def from_sequence(cls, seq: ShadowSequence, total: Optional[Bound] = None,
                      fov: str = "batch") -> "CombinatorialFilter":
        f = cls(seq.initial, total, fov)
        f.extend(seq.events)
        return f

    def apply(self, e: Event) -> None:
        self.istate.apply(e)
        self._cache = None

    def extend(self, events: Iterable[Event]) -> None:
        for e in events:
            self.apply(e)

    def alive(self) -> List[int]:
        return self.istate.alive()

    # -- network construction ---------------------------------------------------
    def _build(self, right_group: Sequence[int] = (), left_group: Sequence[int] = ()) -> _Net:
        """Augmented graph of Fig. 11(d); a query group is routed through node ``G``."""
        st = self.istate
        st.flush()
        net = FlowNetwork()
        for k in (_S, _T, _I):
            net.add_node(k)
        net.add_edge(_T, _S)
        total = net.add_edge(_S, _I, *(self.total or (0, INF)))
        initial = set(st.initial)
        group = -1
        lg = set(left_group)
        if lg:
            kinds = {v in initial for v in lg}
            if len(kinds) != 1:
                raise ValueError("a left query set must contain only initial or only appeared shadows")
            group = net.add_edge(_I if kinds.pop() else _S, _G)
        supply = _add_left(net, st, lambda v: _G if v in lg else (_I if v in initial else _S))
        rg = set(right_group)
        if rg:
            group = net.add_edge(_G, _T)
        demand, inner = _add_right(net, st, lambda tag, v: _G if tag == "alive" and v in rg else _T)
        return _Net(net, supply, demand, inner, total, group)

    def _feasible_net(self, right_group: Sequence[int] = (), left_group: Sequence[int] = ()) -> _Net:
        """The cached feasible network, or a group network warm-started from its circulation."""
        if self._cache is None:
            n = self._build()
            n.net.find_feasible()
            self._cache = n
        if not right_group and not left_group:
            return self._cache
        n = self._build(right_group, left_group)
        assert n.group >= 0
        # Same feasible set as the cached network: reuse its circulation edge by edge.
        c, f = self._cache, n.net.flow
        cf = c.net.flow
        f[0] = cf[0]  # T -> S
        f[n.total] = cf[c.total]
        for key in ("supply", "demand", "inner"):
            new, old = getattr(n, key), getattr(c, key)
            for k, e in new.items():
                f[e] = cf[old[k]]
        f[n.group] = (sum(cf[c.demand[("alive", s)]] for s in right_group) if right_group
                      else sum(cf[c.supply[v]] for v in left_group))
        return n

    def feasible(self) -> bool:
        """``True`` iff the observations so far admit some target distribution."""
        try:
            self._feasible_net()
        except InfeasibleError:
            return False
        return True

    def _check_alive(self, shadows: Iterable[int]) -> List[int]:
        q = list(dict.fromkeys(shadows))
        for s in q:
            if s not in self.istate.reach:
                raise ValueError(f"shadow {s} is not alive")
        if not q:
            raise ValueError("empty query")
        return q

    # -- bounds -----------------------------------------------------------------------
    def bounds(self, shadows: Iterable[int]) -> Tuple[int, float]:
        """Tight ``(lo, hi)`` on the total number of targets in a set of alive shadows."""
        q = self._check_alive(shadows)
        if len(q) == 1:
            n = self._feasible_net()
            return _ints(n.net.extremes(n.demand[("alive", q[0])]))
        n = self._feasible_net(right_group=q)
        return _ints(n.net.extremes(n.group))

    def all_bounds(self) -> Dict[int, Tuple[int, float]]:
        """``{s: (lo, hi)}`` for every alive shadow (one feasible flow, two augmentations each)."""
        n = self._feasible_net()
        return {s: _ints(n.net.extremes(n.demand[("alive", s)])) for s in self.istate.reach}

    def refine_initial_bounds(self, labels: Optional[Iterable[int]] = None) -> Dict[int, Tuple[int, float]]:
        """Tightened ``(lo, hi)`` for shadows at ``t0`` and appeared shadows (Sec. V-E, Eqs. (9), (10)).

        Computed exactly as the min/max flow on the supply edge ``S -> L_i``.
        """
        vs = self.istate.real_left() if labels is None else self._check_left(labels, allow_empty=True)
        n = self._feasible_net()
        return {v: _ints(n.net.extremes(n.supply[v])) for v in vs}

    def _check_left(self, labels: Iterable[int], allow_empty: bool = False) -> List[int]:
        lg = list(dict.fromkeys(labels))
        for v in lg:
            if v not in self.istate.left:
                raise ValueError(f"shadow {v} is not a left vertex (initial or appeared)")
        if not lg and not allow_empty:
            raise ValueError("empty query")
        return lg

    def left_bounds(self, labels: Iterable[int]) -> Tuple[int, float]:
        """Tight ``(lo, hi)`` on the total number of targets that were in a set of left shadows."""
        lg = self._check_left(labels)
        if len(lg) == 1:
            return self.refine_initial_bounds(lg)[lg[0]]
        n = self._feasible_net(left_group=lg)
        return _ints(n.net.extremes(n.group))

    def initial_total_bounds(self) -> Tuple[int, float]:
        """Tight ``(lo, hi)`` on the number of targets in the shadows at ``t0`` (counting, Sec. V-E)."""
        n = self._feasible_net()
        return _ints(n.net.extremes(n.total))

    def witness(self, shadows: Iterable[int], sense: str = "max") -> Witness:
        """A feasible distribution attaining the ``"max"`` or ``"min"`` bound of ``shadows``."""
        q = self._check_alive(shadows)
        n = self._build(right_group=q)
        n.net.find_feasible()
        value = n.net.maximize(n.group) if sense == "max" else n.net.minimize(n.group)
        if value == INF:
            raise ValueError("unbounded: no finite witness")
        f = n.net.flow
        return Witness(int(value), {v: int(f[e]) for v, e in n.supply.items()},
                       {k: int(f[e]) for k, e in n.inner.items() if f[e]})

    # -- the paper's literal recipe ---------------------------------------------------
    def bounds_paper(self, shadows: Iterable[int], corrected: bool = False) -> Tuple[int, float]:
        """Bounds via the max-flow recipe of Sec. V-C, Eqs. (7) and (8), taken literally.

        Upper bound (Eq. 7): ``c(S,i) = u_i``, ``c(i,T) = l_i`` for disappeared
        shadows, ``c(i,T) = 0`` for other alive shadows, ``inf`` for the query;
        result ``f(Q,T) + sum_i f(i,T) - sum_i c(i,T)`` over disappeared ``i``.
        Lower bound (Eq. 8): ``c(S,i) = l_i``, ``c(i,T) = l_i`` for disappeared
        shadows, ``inf`` for other alive shadows, ``0`` for the query; result
        ``sum_i c(S,i) - sum_j f(j,T)``.

        The recipe does not check feasibility and ignores ``total``.  Eq. (7) is
        exact on feasible inputs.  The literal Eq. (8) is not when disappeared
        shadows have ``l_i < u_i``: capping them at ``l_i`` makes it *overstate*
        the lower bound.  ``corrected=True`` uses ``c(i,T) = u_i`` for them in the
        lower-bound run, which agrees with :meth:`bounds` on all tested inputs.
        """
        q = set(self._check_alive(shadows))
        st = self.istate
        st.flush()

        def run(left_cap, gone_cap, alive_cap):
            net = FlowNetwork()
            net.add_node(_S)
            net.add_node(_T)
            out = {}
            for v, b in st.left.items():
                net.add_edge(_S, ("L", v), 0, left_cap(b))
            for tag, src in (("gone", st.disappeared_reach), ("alive", st.reach)):
                for v, reach in src.items():
                    r = (tag, v)
                    for u in st.sorted_left(reach):
                        net.add_edge(("L", u), ("R", r))
                    cap = gone_cap(st.disappeared[v]) if tag == "gone" else alive_cap(v)
                    out[r] = net.add_edge(("R", r), _T, 0, cap)
            value = net.max_flow(_S, _T)
            return value, {r: net.flow[e] for r, e in out.items()}

        value, f = run(lambda b: b[1], lambda b: b[0], lambda s: INF if s in q else 0)
        lsum = sum(lo for lo, _ in st.disappeared.values())
        if value == INF:
            upper = INF
        else:
            upper = sum(f[("alive", s)] for s in q) + sum(f[("gone", v)] for v in st.disappeared) - lsum
        value, _ = run(lambda b: b[0], (lambda b: b[1]) if corrected else (lambda b: b[0]),
                       lambda s: 0 if s in q else INF)
        lower = sum(lo for lo, _ in st.left.values()) - value
        return _ints((lower, upper))


def _add_left(net: FlowNetwork, st: BipartiteIState, parent, tag: tuple = ()) -> Dict[Vertex, int]:
    """Supply edges ``parent(v) -> L_v`` with the left bounds; returns ``{v: edge}``."""
    return {v: net.add_edge(parent(v), ("L", v) + tag, lo, hi) for v, (lo, hi) in st.left.items()}


def _add_right(net: FlowNetwork, st: BipartiteIState, sink, tag: tuple = ()):
    """Bipartite edges ``L_u -> R_r`` and demand edges ``R_r -> sink(kind, v)``; returns ``(demand, inner)``."""
    demand, inner = {}, {}
    for kind, src in (("gone", st.disappeared_reach), ("alive", st.reach)):
        for v, reach in src.items():
            r = (kind, v)
            for u in st.sorted_left(reach):
                inner[(u, r)] = net.add_edge(("L", u) + tag, ("R", r) + tag)
            lo, hi = st.disappeared[v] if kind == "gone" else (0, INF)
            demand[r] = net.add_edge(("R", r) + tag, sink(kind, v), lo, hi)
    return demand, inner


def _ints(b: Tuple[float, float]) -> Tuple[int, float]:
    return int(b[0]), (INF if b[1] == INF else int(b[1]))


# -- thin helpers for Sec. V-E tasks ---------------------------------------------------
def counting_filter(shadows: Iterable[int], fov: str = "batch") -> CombinatorialFilter:
    """Counting: the number of targets is unknown, every shadow at ``t0`` starts at ``[0, inf)``.

    ``initial_total_bounds()`` then bounds the count; it is determined once ``lo == hi``.
    """
    return CombinatorialFilter({s: (0, INF) for s in shadows}, fov=fov)


def pursuit_evasion_filter(shadows: Iterable[int], fov: str = "batch") -> CombinatorialFilter:
    """Passive pursuit-evasion: one evader hidden in one of the shadows at ``t0`` (``[0, 1]`` each)."""
    return CombinatorialFilter({s: (0, 1) for s in shadows}, total=(1, 1), fov=fov)


def evader_status(filt: CombinatorialFilter) -> Dict[int, str]:
    """``"clear"`` (``[0,0]``), ``"evader"`` (``[1,1]``) or ``"contaminated"`` (``[0,1]``) per alive shadow.

    Classified by thresholds (``hi == 0`` clear, ``lo >= 1`` evader, otherwise
    contaminated), as in the JS port.  ``total=(1, 1)`` only counts targets in
    the shadows at ``t0``: an ``Enter`` not matched by an earlier ``Exit``
    contradicts the single-evader model but stays feasible, giving bounds such
    as ``(1, 2)``, which are reported as ``"evader"``.
    """
    return {s: "clear" if hi == 0 else "evader" if lo >= 1 else "contaminated"
            for s, (lo, hi) in filt.all_bounds().items()}


# -- distinguishability (Sec. V-F) ----------------------------------------------------
class MultiTeamFilter:
    """One :class:`CombinatorialFilter` per team when attributes are not intertwined (Sec. V-F).

    ``initial[team][s] = (lo, hi)``; shadows missing for a team start at
    ``(0, 0)``.  Teams are independent, so bounds for a union of teams are the
    sums of per-team bounds.  Split and merge go to every team; count-carrying
    events need per-team counts (see :meth:`apply`).
    """

    def __init__(self, initial: Mapping[Hashable, Mapping[int, Bound]],
                 total: Optional[Mapping[Hashable, Bound]] = None, fov: str = "batch") -> None:
        labels = sorted({s for b in initial.values() for s in b})
        self.teams = list(initial)
        self.filters = {t: CombinatorialFilter({s: initial[t].get(s, (0, 0)) for s in labels},
                                               (total or {}).get(t), fov) for t in self.teams}

    def apply(self, e: Event, counts: Optional[Mapping[Hashable, object]] = None) -> None:
        """Apply ``e`` to every team.

        ``Appear``/``Disappear``: ``counts[team] = (lo, hi)``; teams not listed get
        ``(0, 0)`` (``counts`` may be omitted only when ``e.lo == e.hi == 0``).
        ``Enter``/``Exit``: ``counts[team] = k``; the event goes to listed teams
        with ``k > 0``.  ``Split``/``Merge``: ``counts`` is ignored.
        """
        if isinstance(e, (Split, Merge)):
            for f in self.filters.values():
                f.apply(e)
            return
        if counts is None:
            if isinstance(e, (Appear, Disappear)) and e.lo == e.hi == 0:
                counts = {}
            else:
                raise ValueError(f"{e}: per-team counts required")
        unknown = set(counts) - set(self.filters)
        if unknown:
            raise ValueError(f"unknown teams {unknown}")
        for t, f in self.filters.items():
            if isinstance(e, (Appear, Disappear)):
                lo, hi = counts.get(t, (0, 0))
                f.apply(type(e)(e.s, lo, hi))
            elif counts.get(t, 0) > 0:
                f.apply(type(e)(e.s, counts[t]))

    def bounds(self, shadows: Iterable[int], teams: Optional[Iterable[Hashable]] = None) -> Tuple[int, float]:
        """Bounds on the number of targets of the given teams (default: all) in ``shadows``."""
        q = list(shadows)
        bs = [self.filters[t].bounds(q) for t in (self.teams if teams is None else teams)]
        return sum(b[0] for b in bs), sum(b[1] for b in bs)

    def all_bounds(self) -> Dict[Hashable, Dict[int, Tuple[int, float]]]:
        return {t: f.all_bounds() for t, f in self.filters.items()}


@dataclass
class TwoPassBounds:
    """Result of the paper's literal two-pass procedure of Sec. V-F (:meth:`RedOrBlueFilter.two_pass`).

    ``passes[t]`` holds per-team ``(lo, hi)`` when every mixed target is assigned
    to team ``t`` (``None`` if that assignment contradicts the observations).
    In the paper's terms (``a`` = red): ``lower = l_r1 + l_b1`` targets of either
    team, with between ``first_team_lower[0] = l_r1`` (mixed targets blue) and
    ``first_team_lower[1] = l_r2`` (mixed targets red) of them red; the upper
    bound is reported the same way.  When the two passes disagree on the total,
    ``lower`` is the smaller and ``upper`` the larger of the two; both are
    ``None`` when neither pass is feasible.

    Not sound in general: see :meth:`RedOrBlueFilter.two_pass`.
    """

    passes: Dict[Hashable, Optional[Dict[Hashable, Tuple[int, float]]]]
    lower: Optional[int]
    upper: Optional[float]
    first_team_lower: Tuple[Optional[int], Optional[int]]
    first_team_upper: Tuple[Optional[float], Optional[float]]


class RedOrBlueFilter:
    """Two teams with some targets known only as "``a`` or ``b``" (Sec. V-F, second case).

    ``pure[team][s]`` are single-team initial bounds and ``mixed[s]`` bounds on
    targets in ``s`` that may belong to either team.  Events go to both teams as
    in :meth:`MultiTeamFilter.apply`.

    Teams mix only in the initial condition, so the problem is still a single
    flow: :meth:`bounds` builds one network holding one copy of the augmented
    graph per team (shared ``S`` and ``T``).  A mixed shadow ``s`` gets a node
    ``M_s`` fed by ``S -> M_s [lo, hi]`` that sends into each team's ``L_s``
    (the pure part stays ``S -> L_s [pl, ph]``).  Feasible integral circulations
    are exactly the team-labelled distributions consistent with the
    observations, so the bounds are exact over all red/blue assignments of the
    mixed targets.

    :meth:`two_pass` keeps the paper's literal procedure for comparison.
    """

    def __init__(self, teams: Tuple[Hashable, Hashable], pure: Mapping[Hashable, Mapping[int, Bound]],
                 mixed: Mapping[int, Bound], fov: str = "batch") -> None:
        self.teams = tuple(teams)
        if len(self.teams) != 2 or self.teams[0] == self.teams[1]:
            raise ValueError("RedOrBlueFilter needs two distinct teams")
        unknown = set(pure) - set(self.teams)
        if unknown:
            raise ValueError(f"unknown teams {unknown}")
        self.mixed = {int(s): b for s, b in mixed.items()}
        base = {t: dict(pure.get(t, {})) for t in self.teams}
        for t in self.teams:
            for s in self.mixed:
                base[t].setdefault(s, (0, 0))
        self.joint = MultiTeamFilter(base, fov=fov)
        self.passes = {}
        for owner in self.teams:
            init = {t: dict(base[t]) for t in self.teams}
            for s, (lo, hi) in self.mixed.items():
                plo, phi = init[owner][s]
                init[owner][s] = (plo + lo, phi + hi)
            self.passes[owner] = MultiTeamFilter(init, fov=fov)

    def apply(self, e: Event, counts: Optional[Mapping[Hashable, object]] = None) -> None:
        self.joint.apply(e, counts)
        for p in self.passes.values():
            p.apply(e, counts)

    def alive(self) -> List[int]:
        return self.joint.filters[self.teams[0]].alive()

    def _build(self, shadows: Sequence[int], teams: Sequence[Hashable]) -> Tuple[FlowNetwork, int]:
        net = FlowNetwork()
        net.add_edge(_T, _S)
        for s, (lo, hi) in self.mixed.items():
            net.add_edge(_S, ("M", s), lo, hi)
        q, qt = set(shadows), set(teams)
        for t in self.teams:
            st = self.joint.filters[t].istate
            st.flush()
            initial = set(st.initial)
            _add_left(net, st, lambda v: _S, (t,))
            for s in self.mixed:
                if s in initial:
                    net.add_edge(("M", s), ("L", s, t))
            into = t in qt
            _add_right(net, st, lambda kind, v: _G if into and kind == "alive" and v in q else _T, (t,))
        return net, net.add_edge(_G, _T)

    def feasible(self) -> bool:
        """``True`` iff some red/blue assignment of the mixed targets fits the observations."""
        net, _ = self._build((), ())
        try:
            net.find_feasible()
        except InfeasibleError:
            return False
        return True

    def bounds(self, shadows: Iterable[int], teams: Optional[Iterable[Hashable]] = None) -> Tuple[int, float]:
        """Exact ``(lo, hi)`` on the number of targets of ``teams`` (default: both) in ``shadows``.

        The min/max is over all assignments of the mixed targets to teams;
        ``bounds(q, [a])`` is the paper's "between ``l_r1`` and ``l_r2`` red".
        Raises :class:`InfeasibleError` if no assignment fits the observations.
        """
        q = self.joint.filters[self.teams[0]]._check_alive(shadows)
        ts = list(self.teams if teams is None else dict.fromkeys(teams))
        unknown = set(ts) - set(self.teams)
        if unknown:
            raise ValueError(f"unknown teams {unknown}")
        net, g = self._build(q, ts)
        net.find_feasible()
        return _ints(net.extremes(g))

    def two_pass(self, shadows: Iterable[int]) -> TwoPassBounds:
        """The paper's literal two-pass procedure (Sec. V-F), kept for comparison with :meth:`bounds`.

        Each pass gives every mixed target to one team and runs one filter per
        team.  This is **not sound** in general: the paper assumes both passes
        agree (``l_r1 + l_b1 = l_r2 + l_b2``), but when per-team observations
        show that a mixed group held both colours, a pass (or both) contradicts
        them and the surviving pass can report a range that excludes the true
        count.  A failed pass is reported as ``None``; this raises
        :class:`InfeasibleError` only when the observations themselves are
        inconsistent (no assignment fits, as in :meth:`bounds`).
        """
        q = self.joint.filters[self.teams[0]]._check_alive(shadows)
        if not self.feasible():
            raise InfeasibleError("no assignment of the mixed targets fits the observations")
        a, b = self.teams
        res: Dict[Hashable, Optional[Dict[Hashable, Tuple[int, float]]]] = {}
        for owner, p in self.passes.items():
            try:
                res[owner] = {t: p.filters[t].bounds(q) for t in self.teams}
            except InfeasibleError:
                res[owner] = None
        ok = [r for r in res.values() if r is not None]
        pick = [res[b], res[a]]
        return TwoPassBounds(
            res,
            min((r[a][0] + r[b][0] for r in ok), default=None),
            max((r[a][1] + r[b][1] for r in ok), default=None),
            tuple(None if r is None else r[a][0] for r in pick),
            tuple(None if r is None else r[a][1] for r in pick),
        )


__all__ = [
    "CombinatorialFilter", "Witness", "InfeasibleError", "counting_filter", "pursuit_evasion_filter",
    "evader_status", "MultiTeamFilter", "RedOrBlueFilter", "TwoPassBounds",
]
