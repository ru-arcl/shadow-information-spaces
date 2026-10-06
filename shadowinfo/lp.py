"""Linear-programming baseline on the shadow sequence (T-RO 2012, Sec. IV).

One unknown ``x_i`` per shadow; constraints (2)–(4):

* initial shadow, appear, disappear: ``l_i <= x_i <= u_i``;
* split ``s -> a, b``: ``x_s = x_a + x_b``; merge ``a, b -> s``: ``x_a + x_b = x_s``;
* ``x_i >= 0``.

FOV events become component events with their own unknowns (Sec. V-D):
``Enter(s, k)`` is an appear ``x_t = k`` plus a merge ``x_s + x_t = x_s'`` and
``Exit(s, k)`` a split ``x_s = x_s' + x_t`` plus a disappear ``x_t = k``.  An
optional ``total=(lo, hi)`` bounds the sum over the shadows at ``t0``.

The constraint matrix is totally unimodular (Lemma 3), so the LP optimum is
integral (Prop. 4).  This module does not use the bipartite I-state or any
flow code and serves as an independent oracle in the tests.
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from scipy.optimize import linprog

from .events import Appear, Bound, Disappear, Enter, Exit, Merge, ShadowSequence, Split
from .maxflow import InfeasibleError


class ShadowLP:
    """Constraints (2)–(4) for a shadow sequence, ready for :func:`scipy.optimize.linprog`."""

    def __init__(self, seq: ShadowSequence, total: Optional[Bound] = None) -> None:
        self.lo: List[float] = []
        self.hi: List[float] = []
        self.eq: List[Dict[int, int]] = []
        self.first: Dict[int, int] = {}
        self.last: Dict[int, int] = {}
        for s, (lo, hi) in seq.initial.items():
            self._var(s, lo, hi)
        self.initial = [self.first[s] for s in seq.initial]
        for e in seq.events:
            if isinstance(e, Appear):
                self._var(e.s, e.lo, e.hi)
            elif isinstance(e, Disappear):
                v = self.last[e.s]
                self.lo[v] = max(self.lo[v], e.lo)
                self.hi[v] = min(self.hi[v], e.hi)
            elif isinstance(e, Split):
                s = self.last[e.s]
                self.eq.append({s: 1, self._var(e.a): -1, self._var(e.b): -1})
            elif isinstance(e, Merge):
                a, b = self.last[e.a], self.last[e.b]
                self.eq.append({a: 1, b: 1, self._var(e.s): -1})
            elif isinstance(e, Enter):
                old, t = self.last[e.s], self._var(None, e.k, e.k)
                self.eq.append({old: 1, t: 1, self._var(e.s): -1})
            elif isinstance(e, Exit):
                old, t = self.last[e.s], self._var(None, e.k, e.k)
                self.eq.append({old: 1, self._var(e.s): -1, t: -1})
        self.total = total

    def _var(self, s: Optional[int], lo: float = 0, hi: float = math.inf) -> int:
        v = len(self.lo)
        self.lo.append(lo)
        self.hi.append(hi)
        if s is not None:
            self.first.setdefault(s, v)
            self.last[s] = v
        return v

    def _solve(self, c: np.ndarray):
        n = len(self.lo)
        a_eq = np.zeros((len(self.eq), n))
        for r, row in enumerate(self.eq):
            for v, coef in row.items():
                a_eq[r, v] = coef
        a_ub = b_ub = None
        if self.total is not None:
            row = np.zeros(n)
            row[self.initial] = 1
            lo, hi = self.total
            a_ub, b_ub = np.array([-row]), np.array([-lo])
            if hi != math.inf:
                a_ub, b_ub = np.vstack([a_ub, row]), np.append(b_ub, hi)
        bounds = [(lo, None if hi == math.inf else hi) for lo, hi in zip(self.lo, self.hi)]
        return linprog(c, A_ub=a_ub, b_ub=b_ub, A_eq=a_eq if self.eq else None,
                       b_eq=np.zeros(len(self.eq)) if self.eq else None, bounds=bounds, method="highs")

    def feasible(self) -> bool:
        if any(lo > hi for lo, hi in zip(self.lo, self.hi)):
            return False
        return self._solve(np.zeros(len(self.lo))).status == 0

    def extremes(self, variables: Iterable[int]) -> Tuple[int, float]:
        """``(min, max)`` of the sum of the given unknowns; raises :class:`InfeasibleError`."""
        if not self.feasible():
            raise InfeasibleError("LP infeasible")
        c = np.zeros(len(self.lo))
        c[list(variables)] = 1
        lo = self._solve(c)
        hi = self._solve(-c)
        if lo.status != 0:
            raise RuntimeError(f"LP min failed: {lo.message}")
        upper = math.inf if hi.status == 3 else _integral(-hi.fun, hi)
        return _integral(lo.fun, lo), upper


def _integral(x: float, res) -> int:
    if res.status != 0 or abs(x - round(x)) > 1e-6:
        raise RuntimeError(f"LP result not integral/optimal: {x} ({res.message})")
    return int(round(x))


def lp_bounds(seq: ShadowSequence, shadows: Iterable[int], total: Optional[Bound] = None,
              at: str = "end") -> Tuple[int, float]:
    """LP ``(lo, hi)`` on the number of targets in ``shadows``.

    ``at="end"`` refers to each shadow's last unknown (its count now / when it
    disappeared); ``at="start"`` to its first one (at ``t0`` or when it appeared).
    """
    lp = ShadowLP(seq, total)
    pick = lp.last if at == "end" else lp.first
    return lp.extremes(pick[s] for s in dict.fromkeys(shadows))


def lp_all_bounds(seq: ShadowSequence, total: Optional[Bound] = None) -> Dict[int, Tuple[int, float]]:
    """LP bounds for every shadow alive at the end of ``seq``."""
    lp = ShadowLP(seq, total)
    return {s: lp.extremes([lp.last[s]]) for s in seq.alive_at_end()}


def lp_feasible(seq: ShadowSequence, total: Optional[Bound] = None) -> bool:
    return ShadowLP(seq, total).feasible()


__all__ = ["ShadowLP", "lp_bounds", "lp_all_bounds", "lp_feasible"]
