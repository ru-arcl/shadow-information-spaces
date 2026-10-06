"""Edmonds–Karp, flows with lower bounds, and extreme flows on an edge (vs. scipy LP)."""

from __future__ import annotations

import math
import random

import numpy as np
import pytest
from scipy.optimize import linprog

from shadowinfo.maxflow import FlowNetwork, InfeasibleError

INF = math.inf


def test_max_flow_textbook():
    net = FlowNetwork()
    for u, v, c in [("s", "a", 16), ("s", "b", 13), ("a", "b", 10), ("b", "a", 4), ("a", "c", 12),
                    ("c", "b", 9), ("b", "d", 14), ("d", "c", 7), ("c", "t", 20), ("d", "t", 4)]:
        net.add_edge(u, v, 0, c)
    assert net.max_flow("s", "t") == 23


def test_max_flow_infinite_path_is_unbounded():
    net = FlowNetwork()
    net.add_edge("s", "a")
    net.add_edge("a", "t")
    net.add_edge("s", "t", 0, 3)
    assert net.max_flow("s", "t") == INF


def test_max_flow_infinite_edges_finite_cut():
    net = FlowNetwork()
    net.add_edge("s", "a", 0, 5)
    net.add_edge("a", "b")
    net.add_edge("b", "t")
    net.add_edge("s", "b", 0, 2)
    assert net.max_flow("s", "t") == 7


def test_feasible_with_lower_bounds_and_extremes():
    net = FlowNetwork()
    net.add_edge("t", "s")
    e1 = net.add_edge("s", "a", 2, 4)
    e2 = net.add_edge("a", "x", 0)
    e3 = net.add_edge("a", "y", 0)
    net.add_edge("x", "t", 1, 1)
    ey = net.add_edge("y", "t", 0)
    net.find_feasible()
    assert net.is_circulation()
    assert net.extremes(ey) == (1, 3)
    assert net.extremes(e1) == (2, 4)
    assert net.extremes(e2) == (1, 1)
    assert net.is_circulation()
    assert net.extremes(e3) == (1, 3)


def test_extremes_unbounded():
    net = FlowNetwork()
    net.add_edge("t", "s")
    net.add_edge("s", "a", 1, INF)
    e = net.add_edge("a", "t", 0, INF)
    net.find_feasible()
    assert net.extremes(e) == (1, INF)


def test_unbounded_maximize_leaves_flow_unchanged():
    # BFS finds the short finite path y->a->x before the infinite y->b->c->x.
    net = FlowNetwork()
    e = net.add_edge("x", "y", 0, INF)
    net.add_edge("y", "a", 0, 1)
    net.add_edge("a", "x", 0, 1)
    net.add_edge("y", "b")
    net.add_edge("b", "c")
    net.add_edge("c", "x")
    net.find_feasible()
    before = net.snapshot()
    assert net.maximize(e) == INF
    assert net.flow == before and net.is_circulation()


def test_infeasible_lower_bounds():
    net = FlowNetwork()
    net.add_edge("t", "s")
    net.add_edge("s", "a", 3, 3)
    net.add_edge("a", "t", 0, 2)
    with pytest.raises(InfeasibleError):
        net.find_feasible()
    bad = FlowNetwork()
    bad.add_edge("a", "b", 3, 2)
    with pytest.raises(InfeasibleError):
        bad.find_feasible()


def test_find_feasible_leaves_network_unchanged():
    net = FlowNetwork()
    net.add_edge("t", "s")
    net.add_edge("s", "a", 1, 2)
    net.add_edge("a", "t", 1, 5)
    n_nodes, n_edges, adj = len(net.keys), net.num_edges, [list(a) for a in net.adj]
    net.find_feasible()
    assert (len(net.keys), net.num_edges, net.adj) == (n_nodes, n_edges, adj)
    assert len(net.index) == n_nodes


def _lp_extremes(nodes, edges, e):
    m = len(edges)
    a_eq = np.zeros((len(nodes), m))
    for k, (u, v, _, _) in enumerate(edges):
        a_eq[v, k] += 1
        a_eq[u, k] -= 1
    bounds = [(lo, None if hi == INF else hi) for _, _, lo, hi in edges]
    c = np.zeros(m)
    c[e] = 1
    r1 = linprog(c, A_eq=a_eq, b_eq=np.zeros(len(nodes)), bounds=bounds, method="highs")
    if r1.status == 2:
        return None
    r2 = linprog(-c, A_eq=a_eq, b_eq=np.zeros(len(nodes)), bounds=bounds, method="highs")
    return round(r1.fun), INF if r2.status == 3 else round(-r2.fun)


def test_random_circulations_match_lp():
    n_feasible = 0
    for seed in range(250):
        rng = random.Random(seed)
        n = rng.randint(3, 7)
        edges = []
        for _ in range(rng.randint(n, 3 * n)):
            u, v = rng.sample(range(n), 2)
            lo = rng.choice([0, 0, 0, 1, 2])
            hi = INF if rng.random() < 0.2 else lo + rng.randint(0, 4)
            edges.append((u, v, lo, hi))
        net = FlowNetwork()
        for i in range(n):
            net.add_node(i)
        for u, v, lo, hi in edges:
            net.add_edge(u, v, lo, hi)
        expect = [_lp_extremes(range(n), edges, e) for e in range(len(edges))]
        if expect[0] is None:
            with pytest.raises(InfeasibleError):
                net.find_feasible()
            continue
        n_feasible += 1
        net.find_feasible()
        assert net.is_circulation()
        for e in range(len(edges)):
            assert net.extremes(e) == expect[e], (seed, e)
            assert net.is_circulation()
    assert n_feasible > 50
