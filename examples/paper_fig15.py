"""T-RO 2012 Fig. 15 (= ICRA'08 Fig. 8): the two large runs of the original Java, through the Python port.

* Fig. 15(a): map 13 and the commented-out ``// polygon 13`` path of ``ProjectPanel`` (first 9 points,
  ``P13``), 100 targets;
* Fig. 15(b): map 14 and the ``ProjectPanel5`` path (``P14b``), 10^6 targets.

For each run the script computes the critical lines and critical points (``getCuts``,
``getAllCriticalPoints``), tracks the gaps (``getGaps``), distributes the targets with the
single-type-agent oracle with the merge bug fixed (``docs/notes/original_java.md`` B7), builds the
bipartite I-state (Sec. V-B) and computes exact bounds on every final shadow (Sec. V-C) with
:class:`~shadowinfo.nondeterministic.CombinatorialFilter`.  It does this twice:

* ``compat="java"`` -- the Java gap history, bit for bit (the paper's numbers come from it);
* ``compat="fixed"`` (the default of the public API) -- gap tracking that never crashes; wherever Java
  does not crash (as on both runs here) it gives Java's history, so the two columns agree.

It also reports the opt-in ``get_gaps(..., java_matching=False)``, which keeps each ID on its physical
shadow where Java's ``samePhysicalGap`` matching shifts IDs along the boundary (B3); the simulator uses it.

It then prints the counts next to the numbers the paper reports and explains every difference.  The
paper counts are collected in ``docs/notes/original_maps.md`` ("Figure -> map table").

Usage::

    python examples/paper_fig15.py [--seed 1] [--paths P13,P14b]
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo import BipartiteIState, InfeasibleError  # noqa: E402
from shadowinfo.polygon import (SingleTypeAgentOracle, all_critical_points, all_cuts,  # noqa: E402
                                derive_shadow_info_state, exact_bounds, gap_history_to_events, get_gaps,
                                java_bounds, load_path, sample_path)

PAPER = {
    "P13": dict(figure="T-RO Fig. 15(a) = ICRA'08 Fig. 8(a)", n_targets=100, events=85, labels=None, bounded=18,
                bg_vertices=41, bg_edges=60, time="0.1 s for one team"),
    "P14b": dict(figure="T-RO Fig. 15(b) = ICRA'08 Fig. 8(b)", n_targets=1_000_000, events=385, labels=491,
                 bounded=12, bg_vertices=124, bg_edges=339, time="under 1 s (T-RO); 2.5 s (ICRA'08)"),
}
"""What the papers report for the two runs (``docs/notes/original_maps.md``)."""

KINDS = ("Split", "Merge", "Appear", "Disappear")


def analyse(history, n_targets: int, seed: int) -> dict:
    """Events, bipartite graphs, oracle observations and bounds of one gap history."""
    r: dict = {"history": history}
    ge = gap_history_to_events(history)
    r["events"] = len(ge.events)
    r["kinds"] = Counter(type(e).__name__ for e in ge.events)
    r["labels"] = history.max_id
    r["alive"] = len(history.final_ids)
    r["initial"] = len(ge.initial)
    st = BipartiteIState.from_sequence(ge.sequence())
    r["si_bg"] = (len(st.left), len(st.right()), len(st.edges()))

    t0 = time.perf_counter()
    oracle = SingleTypeAgentOracle(n_targets, rng=seed, merge_bug=False, compat=history.compat).initialize(history)
    t1 = time.perf_counter()
    bg = derive_shadow_info_state(history, oracle)
    t2 = time.perf_counter()
    r["bounds"] = exact_bounds(oracle)
    t3 = time.perf_counter()
    r["t_oracle"], r["t_bg"], r["t_bounds"] = t1 - t0, t2 - t1, t3 - t2
    d = bg.to_dict()
    r["java_bg"] = (d["n_left"], d["n_right"], d["n_edges"])
    r["oracle"] = oracle
    r["truth"] = oracle.final_counts
    jb = java_bounds(bg, "sorted")
    r["java_bounded"] = len(jb)
    r["java_agree"] = sum(jb[s] == r["bounds"][s] for s in jb)
    buggy = SingleTypeAgentOracle(n_targets, rng=seed, merge_bug=True, compat=history.compat).initialize(history)
    try:
        exact_bounds(buggy)
        r["original_oracle"] = "consistent"
    except InfeasibleError:
        r["original_oracle"] = "infeasible (B7)"
    return r


def fmt_b(b) -> str:
    lo, hi = b
    return f"[{lo}, {'inf' if hi == math.inf else hi}]"


def report(label: str, seed: int) -> dict:
    paper = PAPER[label]
    poly, path = load_path(label)
    t0 = time.perf_counter()
    cuts = all_cuts(poly)
    t1 = time.perf_counter()
    samples = sample_path(poly, path, "fixed")
    t2 = time.perf_counter()
    fixed = get_gaps(poly, path, "fixed", samples=samples)
    t3 = time.perf_counter()
    java = get_gaps(poly, path, "java", samples=samples if samples.mode_independent else None)
    n_cp = len(all_critical_points(path, cuts))
    rj = analyse(java, paper["n_targets"], seed)
    rf = analyse(fixed, paper["n_targets"], seed)

    print(f"{paper['figure']}: map {poly.name} ({poly.n} vertices, {len(cuts)} critical lines), path {label} "
          f"({len(path.points)} waypoints, {n_cp} critical points), N = {paper['n_targets']:,} targets, seed {seed}")
    print(f"  {'':34s} {'paper':>9}  {'compat=java':>13}  {'default':>13}")

    def row(name, p, j, f, note=""):
        p = "-" if p is None else p
        print(f"  {name:34s} {p!s:>9}  {j!s:>13}  {f!s:>13}  {note}")

    def kinds(r):
        return ", ".join(f"{k.lower()} {r['kinds'][k]}" for k in KINDS)

    row("component events", paper["events"], rj["events"], rf["events"], kinds(rf))
    row("shadow labels (IDs handed out)", paper["labels"], rj["labels"], rf["labels"])
    row("shadows alive at the end", None, rj["alive"], rf["alive"])
    row("  of which Java prints bounds for", paper["bounded"], rj["java_bounded"], rf["java_bounded"],
        "appeared and never split or merged: no bound (B11)")
    jv = lambda r: f"{r['java_bg'][0]}+{r['java_bg'][1]}={r['java_bg'][0] + r['java_bg'][1]}"  # noqa: E731
    sv = lambda r: f"{r['si_bg'][0]}+{r['si_bg'][1]}={r['si_bg'][0] + r['si_bg'][1]}"  # noqa: E731
    row("BG vertices, Java (pooled t0)", paper["bg_vertices"], jv(rj), jv(rf), "left + right")
    row("BG edges, Java (pooled t0)", paper["bg_edges"], rj["java_bg"][2], rf["java_bg"][2])
    row("BG vertices, one per t0 shadow", paper["bg_vertices"], sv(rj), sv(rf), "BipartiteIState")
    row("BG edges, one per t0 shadow", paper["bg_edges"], rj["si_bg"][2], rf["si_bg"][2])

    hidden = rf["oracle"].hidden_total
    print(f"  exact bounds, default mode (merge-fixed oracle: {hidden:,} hidden at t0; true count in parentheses):")
    items = [f"s{s} {fmt_b(b)} ({rf['truth'].get(s)})" for s, b in rf["bounds"].items()]
    for i in range(0, len(items), 4):
        print("    " + "   ".join(items[i:i + 4]))
    for name, r in (("compat=java", rj), ("default", rf)):
        inside = all(b[0] <= r["truth"][s] <= b[1] for s, b in r["bounds"].items())
        print(f"  {name:12s} truth inside the exact bounds: {'yes' if inside else 'NO'}; Java's max-flow bounds "
              f"(sorted order) equal them on {r['java_agree']}/{r['java_bounded']} shadows; "
              f"original oracle (merge bug): {r['original_oracle']}")
    print(f"  time: cuts {t1 - t0:.2f} s, visibility at the critical points {t2 - t1:.2f} s, gap tracking "
          f"{t3 - t2:.3f} s; oracle {rf['t_oracle']:.3f} s, bipartite graph {rf['t_bg']:.3f} s, exact bounds "
          f"{rf['t_bounds']:.3f} s (paper: {paper['time']})")
    if fixed.inferred:
        print(f"  default mode: {len(fixed.inferred)} critical points re-derived from the hidden chains, where Java's "
              f"step throws or is inconsistent (critical points {fixed.inferred[0][0]}..{fixed.inferred[-1][0]})")
    own = get_gaps(poly, path, "fixed", java_matching=False, samples=samples)
    ev_f, ev_o = gap_history_to_events(fixed).events, gap_history_to_events(own).events
    n_diff = sum(x != y for x, y in zip(ev_f, ev_o)) + abs(len(ev_f) - len(ev_o))
    if n_diff:
        e_own = len(BipartiteIState.from_sequence(gap_history_to_events(own).sequence()).edges())
        print(f"  get_gaps(..., java_matching=False) (opt-in, IDs stay on their physical shadows): {n_diff} of "
              f"{len(ev_f)} events name different shadows (critical points {own.inferred[0][0]}.."
              f"{own.inferred[-1][0]}); BG edges, one per t0 shadow: {e_own}")
    else:
        print("  get_gaps(..., java_matching=False) (opt-in): the same history")
    print()
    return {"java": rj, "fixed": rf}


NOTES = """\
How the paper's numbers relate to the port (docs/notes/original_java.md §1, original_maps.md):
  * events, labels: identical; 385 events / 491 labels on Fig. 15(b), 85 events on Fig. 15(a).
  * "final shadows" (18, 12) are the shadows Java prints bounds for: deriveShadowBounds drops a shadow
    that appeared and never split or merged (B11), although its count is known; 2 such shadows are alive
    at the end of each run.  The exact bounds above include them ([k, k]).
  * bipartite graph: Java pools all shadows at t0 into ONE left vertex (B12).  Fig. 15(a)'s 41 vertices
    are the pooled count (12 + 29); Fig. 15(b)'s 124 vertices are the per-shadow count (61 + 63) while its
    339 edges are the pooled count.  No convention found gives Fig. 15(a)'s 60 edges (pooled: 80).
  * the default mode gives Java's gap history wherever Java does not crash, so its columns equal
    compat=java.  The opt-in java_matching=False keeps IDs on their shadows where Java shifts them (B3):
    on Fig. 15(b) some events then connect different shadows; the counts of events, labels and final
    shadows do not change.
  * the original oracle double-counts merges (B7), which makes the observations inconsistent (all five
    recorded Fig. 15(b) runs are infeasible); the oracle here conserves targets, so the exact bounds always contain the truth.
  * "10^6 targets in 5 teams" (Fig. 15(b)): the Java code has a single team; the run uses one team of 10^6.
  * times: the Java panels print "Computation time" for the equations, bipartite graph and bounds only,
    not the geometry; the paper's times presumably are that line, on the hardware of 2008-2012.  Here
    the same steps take milliseconds and the geometry (pure Python) about two seconds on Fig. 15(b).
Reproduce the Java console output: python -m shadowinfo.polygon run ProjectPanel5 --compat java --seed 1"""


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=1, help="seed of the oracle's java.util.Random")
    ap.add_argument("--paths", default="P13,P14b", help="comma-separated: P13 (Fig. 15a), P14b (Fig. 15b)")
    a = ap.parse_args(argv)
    for label in a.paths.split(","):
        if label not in PAPER:
            ap.error(f"unknown path {label!r}; choose from {', '.join(PAPER)}")
        report(label, a.seed)
    print(NOTES)


if __name__ == "__main__":
    main()
