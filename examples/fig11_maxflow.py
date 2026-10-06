"""T-RO 2012 Fig. 11: bounds on the targets in shadow s19 by max-flow (Sec. V-C).

The shadow sequence of Fig. 11(b) (``tests/fixtures/tor_fig11.json``) is fed
into the incremental bipartite I-state (Fig. 11(c)).  The bounds on ``s19`` are
then computed three independent ways:

* exact feasible flow with lower bounds on the augmented graph of Fig. 11(d)
  (:meth:`CombinatorialFilter.bounds`);
* the paper's literal recipe, Eqs. (7) and (8) (:meth:`CombinatorialFilter.bounds_paper`);
* the LP relaxation of Sec. IV, Eqs. (2)-(4), solved with HiGHS (:func:`lp_bounds`);
  its optimum is integral because the constraint matrix is totally unimodular.

The paper reports ``10 <= s19 <= 24``.

The second part runs the geometry that produced the figure: map 12 of the original Java code and the
path of ``ProjectPanel`` (``P12``), through :mod:`shadowinfo.polygon`.  The gap tracker yields the
shadow labels 1-19 of ICRA'08 Fig. 4; T-RO Fig. 11(b) swaps labels 6 and 7 and draws s5 (which
appears at t = 0.056) as an initial shadow (``docs/notes/original_java.md`` §1).  After that
relabelling the events and the bipartite I-state are those of Fig. 11(b)-(c), and with the paper's
counts attached the filter bounds s19 to [10, 24] again.

Usage::

    python examples/fig11_maxflow.py
"""

from __future__ import annotations

import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo import (Appear, BipartiteIState, CombinatorialFilter, Disappear, Merge,  # noqa: E402
                        ShadowSequence, Split, lp_all_bounds, lp_bounds, vertex_name)
from shadowinfo.polygon import all_critical_points, all_cuts, gap_history_to_events, get_gaps, load_path  # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "tor_fig11.json"


def fmt(b) -> str:
    lo, hi = b
    return f"[{lo}, {'inf' if hi == math.inf else hi}]"


TRO_LABELS = {6: 7, 7: 6}
"""ICRA'08 Fig. 4 label -> T-RO Fig. 11(b) label (the two children of the first split are swapped)."""


def relabel(e):
    r = lambda x: TRO_LABELS.get(x, x)  # noqa: E731
    if isinstance(e, Split):
        return Split(r(e.s), r(e.a), r(e.b))
    if isinstance(e, Merge):
        return Merge(r(e.a), r(e.b), r(e.s))
    return type(e)(r(e.s), e.lo, e.hi)


def shape(e):
    """An event up to the order of the two children of a split / parents of a merge."""
    if isinstance(e, Split):
        return ("split", e.s, frozenset((e.a, e.b)))
    if isinstance(e, Merge):
        return ("merge", frozenset((e.a, e.b)), e.s)
    return (type(e).__name__.lower(), e.s)


def geometric_run(fx: dict) -> bool:
    """Fig. 11(a): map 12 + the ``ProjectPanel`` path through the polygon front end."""
    poly, path = load_path("P12")
    history = get_gaps(poly, path)  # identical to compat="java" on this path
    ge = gap_history_to_events(history)
    n_cp = len(all_critical_points(path, all_cuts(poly)))
    print(f"\nFig. 11(a) geometry: map {poly.name} ({poly.n} vertices), ProjectPanel path P12 "
          f"{' '.join(f'({x:g},{y:g})' for x, y in path.points)}")
    print(f"  {n_cp} critical points, {len(ge.events)} component events, labels 1-{history.max_id}; "
          f"initial shadows {ge.initial}")
    paper = ShadowSequence.from_dict(fx["sequence"])
    first, rest = ge.events[0], ge.events[1:]
    print(f"  t = {ge.times[0]:.3f}  {first}   (T-RO Fig. 11(b) draws s5 as initial)")
    events = [relabel(e) for e in rest]
    for t, e, raw in zip(ge.times[1:], events, rest):
        note = "   (Java labels: " + str(raw) + ")" if e != raw else ""
        print(f"  t = {t:.3f}  {e}{note}")
    same_events = Counter(map(shape, events)) == Counter(map(shape, paper.events))
    order = [shape(e) for e in events] == [shape(e) for e in paper.events]
    print(f"  the same {len(events)} events as Fig. 11(b) after swapping labels 6 and 7: {'yes' if same_events else 'NO'}")
    if same_events and not order:  # list the pairs of events the figure draws in the other order
        pos = {shape(e): i for i, e in enumerate(paper.events)}
        swapped = [f"{events[i]!s} before {events[j]!s}" for i in range(len(events)) for j in range(i + 1, len(events))
                   if pos[shape(events[i])] > pos[shape(events[j])]]
        print("  in a different order from the fixture (events on different lineages; the I-state is the same):")
        for x in swapped:
            print(f"    {x}")

    # attach the paper's counts (initial bounds, appear / disappear bounds) to the geometric events
    counts = {e.s: (e.lo, e.hi) for e in paper.events if isinstance(e, (Appear, Disappear))}
    events = [type(e)(e.s, *counts[e.s]) if isinstance(e, (Appear, Disappear)) else e for e in events]
    seq = ShadowSequence(dict(paper.initial), events)
    edges = sorted([u, v] for u, (_t, v) in BipartiteIState.from_sequence(seq).edges())
    same_bg = edges == sorted(fx["bipartite"]["edges"])
    print(f"  bipartite I-state = Fig. 11(c) ({len(edges)} edges): {'yes' if same_bg else 'NO'}")
    q = fx["query"]
    b = CombinatorialFilter.from_sequence(seq).bounds(q)
    print(f"  with the paper's counts attached, bounds on s{q[0]}: {fmt(b)}; final shadows "
          f"{sorted(history.final_ids)}")
    return same_events and same_bg and b == (fx["paper_lower"], fx["paper_upper"])


def main() -> None:
    fx = json.loads(FIXTURE.read_text())
    seq = ShadowSequence.from_dict(fx["sequence"])
    filt = CombinatorialFilter.from_sequence(seq)
    st = filt.istate

    print("Fig. 11(b) shadow sequence:", len(seq.events), "component events")
    print("  initial bounds:", ", ".join(f"s{s}{fmt(b)}" for s, b in seq.initial.items()))
    print("Fig. 11(c) bipartite I-state")
    print("  left :", ", ".join(f"{vertex_name(v)}{fmt(b)}" for v, b in st.left.items()))
    print("  right:", ", ".join(f"{vertex_name(v)}{fmt(b)} disappeared" if tag == "gone" else f"{vertex_name(v)} alive"
                               for (tag, v), b in st.right().items()))
    print("  edges:", " ".join(f"{vertex_name(u)}-{vertex_name(v)}" for u, (_, v) in st.edges()))

    q = fx["query"]
    flow = filt.bounds(q)
    paper = filt.bounds_paper(q)
    lp = lp_bounds(seq, q)
    print(f"\nBounds on s{q[0]} (paper: [{fx['paper_lower']}, {fx['paper_upper']}])")
    print(f"  feasible flow with lower bounds : {fmt(flow)}")
    print(f"  paper's recipe, Eqs. (7)/(8)    : {fmt(paper)}")
    print(f"  LP relaxation, Eqs. (2)-(4)     : {fmt(lp)}")

    w_hi, w_lo = filt.witness(q, "max"), filt.witness(q, "min")
    for name, w in (("max", w_hi), ("min", w_lo)):
        start = ", ".join(f"{vertex_name(v)}:{n}" for v, n in w.supply.items())
        print(f"  witness ({name}, s{q[0]}={w.value}): start counts {start}")

    flows, lps = filt.all_bounds(), lp_all_bounds(seq)
    print("\nAll final shadows   flow        LP")
    for s in sorted(flows):
        print(f"  s{s:<4}            {fmt(flows[s]):<10}  {fmt(lps[s])}")
    print("  refined initial bounds (Eqs. 9, 10):",
          ", ".join(f"s{s}{fmt(b)}" for s, b in filt.refine_initial_bounds().items()))

    ok = flow == paper == lp == (fx["paper_lower"], fx["paper_upper"]) and flows == lps
    print("\nall methods agree with the paper:", "yes" if ok else "NO")
    geo = geometric_run(fx)
    print("\nthe geometric run reproduces Fig. 11(b)-(c) and the bounds:", "yes" if geo else "NO")
    if not (ok and geo):
        sys.exit(1)


if __name__ == "__main__":
    main()
