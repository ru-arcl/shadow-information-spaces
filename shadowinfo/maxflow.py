"""Edmonds–Karp max-flow and flows with lower bounds (T-RO 2012, Sec. V-C).

The paper reduces target counting over the bipartite I-state to max-flow
(Fig. 11(d)) and solves it with Edmonds–Karp.  :class:`FlowNetwork` provides

* :meth:`FlowNetwork.max_flow` -- classic Edmonds–Karp from the current flow
  (used for the paper's literal recipes, Eqs. (7)–(10));
* :meth:`FlowNetwork.find_feasible` -- a feasible circulation respecting lower
  bounds ``lo <= f(e) <= hi`` via the super source/sink reduction;
* :meth:`FlowNetwork.maximize` / :meth:`FlowNetwork.minimize` -- the extreme
  flow on a designated edge over *all* feasible circulations, obtained from a
  feasible circulation by augmenting around cycles through that edge.

Capacities are integers or ``math.inf``.  Everything is pure Python and
deterministic (arcs are scanned in insertion order).
"""

from __future__ import annotations

import math
from collections import deque
from typing import Dict, Hashable, Iterable, List, Optional, Sequence, Tuple

INF = math.inf


class InfeasibleError(ValueError):
    """Raised when no flow satisfies all lower/upper bounds (inconsistent observations)."""


class FlowNetwork:
    """Directed network with integer lower/upper bounds on edges.

    Edge ``e`` owns arcs ``2e`` (forward, residual ``hi - f``) and ``2e + 1``
    (backward, residual ``f - lo``).  Flows start at ``0`` for :meth:`max_flow`
    and at a feasible circulation after :meth:`find_feasible`.
    """

    def __init__(self) -> None:
        self.keys: List[Hashable] = []
        self.index: Dict[Hashable, int] = {}
        self.adj: List[List[int]] = []
        self.head: List[int] = []
        self.lo: List[float] = []
        self.hi: List[float] = []
        self.flow: List[float] = []

    # -- construction ----------------------------------------------------
    def add_node(self, key: Hashable) -> int:
        if key in self.index:
            return self.index[key]
        self.index[key] = len(self.keys)
        self.keys.append(key)
        self.adj.append([])
        return self.index[key]

    def add_edge(self, u: Hashable, v: Hashable, lo: float = 0, hi: float = INF) -> int:
        """Add edge ``u -> v`` with ``lo <= f <= hi``; returns the edge id."""
        if lo < 0 or lo == INF:
            raise ValueError(f"bad lower bound {lo}")
        iu, iv = self.add_node(u), self.add_node(v)
        e = len(self.lo)
        self.lo.append(lo)
        self.hi.append(hi)
        self.flow.append(0)
        self.head += [iv, iu]
        self.adj[iu].append(2 * e)
        self.adj[iv].append(2 * e + 1)
        return e

    @property
    def num_edges(self) -> int:
        return len(self.lo)

    def tail_of(self, e: int) -> int:
        return self.head[2 * e + 1]

    def head_of(self, e: int) -> int:
        return self.head[2 * e]

    def snapshot(self) -> List[float]:
        return list(self.flow)

    def restore(self, flows: Sequence[float]) -> None:
        self.flow[:] = flows

    # -- residual graph ----------------------------------------------------
    def _residual(self, a: int) -> float:
        e = a >> 1
        return self.hi[e] - self.flow[e] if not a & 1 else self.flow[e] - self.lo[e]

    def _augmenting_path(self, s: int, t: int, skip: int) -> Optional[Tuple[List[int], float]]:
        """Shortest residual path ``s -> t`` (BFS) avoiding edge ``skip``."""
        parent = [-1] * len(self.keys)
        parent[s] = -2
        queue = deque([s])
        head, flow, lo, hi = self.head, self.flow, self.lo, self.hi
        while queue:
            u = queue.popleft()
            for a in self.adj[u]:
                e = a >> 1
                if e == skip:
                    continue
                v = head[a]
                if parent[v] != -1:
                    continue
                if (hi[e] - flow[e] if not a & 1 else flow[e] - lo[e]) <= 0:
                    continue
                parent[v] = a
                if v == t:
                    path = []
                    bottleneck = INF
                    while v != s:
                        a = parent[v]
                        path.append(a)
                        bottleneck = min(bottleneck, self._residual(a))
                        v = head[a ^ 1]
                    return path, bottleneck
                queue.append(v)
        return None

    def _push(self, path: Iterable[int], amount: float) -> None:
        for a in path:
            if a & 1:
                self.flow[a >> 1] -= amount
            else:
                self.flow[a >> 1] += amount

    def _augment(self, s: int, t: int, limit: float = INF, skip: int = -1) -> float:
        """Edmonds–Karp from the current flow; returns the amount pushed (``inf`` if unbounded).

        When an augmenting path of infinite capacity exists and ``limit`` is
        infinite, ``inf`` is returned and the flow is left unchanged (finite
        pushes made before that path was found are undone).
        """
        total = 0
        pushed: List[Tuple[List[int], float]] = []
        while total < limit:
            found = self._augmenting_path(s, t, skip)
            if found is None:
                break
            path, bottleneck = found
            amount = min(bottleneck, limit - total)
            if amount == INF:
                for p, a in reversed(pushed):
                    self._push(p, -a)
                return INF
            self._push(path, amount)
            pushed.append((path, amount))
            total += amount
        return total

    # -- algorithms ----------------------------------------------------------
    def max_flow(self, s: Hashable, t: Hashable) -> float:
        """Edmonds–Karp max-flow value from ``s`` to ``t`` (lower bounds must be 0)."""
        return self._augment(self.index[s], self.index[t])

    def find_feasible(self) -> None:
        """Find a circulation with ``lo <= f <= hi`` on every edge.

        Super source/sink reduction: start from ``f = lo``, route each node's
        excess from a super source to the deficient nodes via a super sink.
        Raises :class:`InfeasibleError` if no such circulation exists.
        """
        for e in range(self.num_edges):
            if self.lo[e] > self.hi[e]:
                raise InfeasibleError(f"edge {self.keys[self.tail_of(e)]}->{self.keys[self.head_of(e)]}: "
                                      f"lo {self.lo[e]} > hi {self.hi[e]}")
        self.flow[:] = self.lo
        excess = [0] * len(self.keys)
        for e in range(self.num_edges):
            excess[self.head_of(e)] += self.lo[e]
            excess[self.tail_of(e)] -= self.lo[e]
        n_nodes, n_edges = len(self.keys), self.num_edges
        ss, tt = self.add_node(object()), self.add_node(object())
        need = 0
        for v in range(n_nodes):
            if excess[v] > 0:
                self.add_edge(self.keys[ss], self.keys[v], 0, excess[v])
                need += excess[v]
            elif excess[v] < 0:
                self.add_edge(self.keys[v], self.keys[tt], 0, -excess[v])
        got = self._augment(ss, tt)
        self._truncate(n_nodes, n_edges)
        if got < need:
            raise InfeasibleError("no flow satisfies the bounds")

    def _truncate(self, n_nodes: int, n_edges: int) -> None:
        for k in self.keys[n_nodes:]:
            del self.index[k]
        del self.keys[n_nodes:], self.adj[n_nodes:]
        del self.lo[n_edges:], self.hi[n_edges:], self.flow[n_edges:], self.head[2 * n_edges:]
        for lst in self.adj:
            while lst and lst[-1] >> 1 >= n_edges:
                lst.pop()

    def maximize(self, e: int) -> float:
        """Raise the flow on edge ``e`` as far as feasibility allows; returns ``f(e)``.

        Requires a feasible circulation (:meth:`find_feasible`).  Any other
        feasible circulation differs from the current one by a circulation in
        the residual graph, so the extra flow on ``e = (x, y)`` is the max-flow
        ``y -> x`` in the residual graph without ``e``, capped by ``hi(e) - f(e)``.
        Returns ``inf`` (flow unchanged) if unbounded.
        """
        x, y = self.tail_of(e), self.head_of(e)
        extra = self._augment(y, x, self.hi[e] - self.flow[e], skip=e)
        if extra == INF:
            return INF
        self.flow[e] += extra
        return self.flow[e]

    def minimize(self, e: int) -> float:
        """Lower the flow on edge ``e`` as far as feasibility allows; returns ``f(e)``."""
        x, y = self.tail_of(e), self.head_of(e)
        self.flow[e] -= self._augment(x, y, self.flow[e] - self.lo[e], skip=e)
        return self.flow[e]

    def extremes(self, e: int) -> Tuple[float, float]:
        """``(min f(e), max f(e))`` over feasible circulations; the flow is restored."""
        saved = self.snapshot()
        lo = self.minimize(e)
        self.restore(saved)
        hi = self.maximize(e)
        self.restore(saved)
        return lo, hi

    def is_circulation(self) -> bool:
        """Check bounds and conservation of the current flow (for tests)."""
        bal = [0] * len(self.keys)
        for e in range(self.num_edges):
            if not self.lo[e] <= self.flow[e] <= self.hi[e]:
                return False
            bal[self.head_of(e)] += self.flow[e]
            bal[self.tail_of(e)] -= self.flow[e]
        return all(b == 0 for b in bal)
