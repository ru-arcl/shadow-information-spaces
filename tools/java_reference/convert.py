"""Convert the Java golden dumps (``tests/fixtures/java/*.json``) to ``shadowinfo`` objects and cross-check them.

The dumps are written by ``tools/java_reference/src/GoldenDump.java``, which runs
the original Java implementation (see ``README.md`` next to this file). This
module:

* reads component events off the Java gap history, ``Gap[][]`` from
  ``Algorithm.getGaps`` (:func:`gap_history_events`);
* converts the Java oracle's event list into a :class:`shadowinfo.events.ShadowSequence`
  (:func:`oracle_sequence`) and checks that it validates;
* checks that the Java bipartite graph (``deriveShadowInfoState``) equals
  :class:`shadowinfo.bipartite.BipartiteIState` built from that sequence, with
  the Java "pooled" initial vertex 0 standing for all initial shadows
  (:func:`compare_bipartite`);
* compares the Java max-flow bounds (``deriveShadowBounds``, every outcome the
  harness saw over several identity-hash orders) with the exact bounds of
  :class:`shadowinfo.nondeterministic.CombinatorialFilter`
  (:func:`compare_bounds`);
* checks the Java never-see-evader labels against the reach sets of the
  bipartite I-state (:func:`compare_nse`);
* counts what the paper reports for Fig. 15(a)/(b) (:func:`paper_counts`).

Usage::

    python tools/java_reference/convert.py               # report on every dump
    python tools/java_reference/convert.py --json out.json
    python tools/java_reference/convert.py 12_PP 14_PP5  # selected dumps
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from shadowinfo.bipartite import BipartiteIState  # noqa: E402
from shadowinfo.events import Appear, Disappear, Merge, ShadowSequence, Split  # noqa: E402
from shadowinfo.maxflow import InfeasibleError  # noqa: E402
from shadowinfo.nondeterministic import CombinatorialFilter  # noqa: E402

FIXTURES = ROOT / "tests" / "fixtures" / "java"
INF = math.inf

# Paper numbers (T-RO 2012 Sec. VI / Fig. 15; ICRA 2008 Fig. 8), see docs/notes/original_maps.md.
PAPER = {
    "fig15a": {"events": 85, "targets": 100, "bg_vertices": 41, "bg_edges": 60, "final_shadows": 18},
    "fig15b": {"events": 385, "shadows": 491, "bg_vertices": 124, "bg_edges": 339, "final_shadows": 12,
               "targets": 1000000},
}


def load(name: str) -> dict:
    p = Path(name)
    if not p.suffix:
        p = FIXTURES / f"{name}.json"
    return json.loads(p.read_text())


def run_files() -> List[Path]:
    return sorted(p for p in FIXTURES.glob("*_*.json"))


# ----------------------------------------------------------------------------- java.util.Random
def java_random_doubles(seed: int, n: int) -> List[float]:
    """``new java.util.Random(seed).nextDouble()`` x n, i.e. ``Math.random()`` after the harness seeds it.

    Reference for the Python port's ``JavaRandom`` clone; checked against ``index.json``.
    """
    mask = (1 << 48) - 1
    s = (seed ^ 0x5DEECE66D) & mask

    def nxt(bits: int) -> int:
        nonlocal s
        s = (s * 0x5DEECE66D + 0xB) & mask
        return s >> (48 - bits)

    return [((nxt(26) << 27) + nxt(27)) * (1.0 / (1 << 53)) for _ in range(n)]


def replay_oracle(run: dict, seed: int, n_targets: int, merge_bug: bool = True) -> List[list]:
    """Re-run ``SingleTypeAgentOracle.distributeAgents`` on the dumped gap history (reference replica).

    Mirrors the Java statement by statement, including the merge bug (``merge_bug=True``:
    the second merging parent sets the target to ``count2 + count2``). Returns rows like the
    dump's ``events``: ``[t, type, fromGap, toGap, movedAgents, visibleAgents]``. The test suite
    checks that this reproduces every Java oracle run, so the dumps plus ``java_random_doubles``
    are enough to port the oracle exactly.
    """
    rand = _JavaRandomStream(seed).next_double

    def dob(t: int, g: int) -> List[int]:
        if g == 1:
            return [t]
        if t == 0:
            return [0] * g
        sp = set()
        while len(sp) < g - 1:
            sp.add(int(rand() * t))
        sp.add(t)
        r = sorted(sp)
        return [r[0]] + [r[i] - r[i - 1] for i in range(1, g)]

    sets = run["get_gaps"]["sets"]
    counts: List[Dict[int, Optional[int]]] = []
    rows = []
    prev_vis = 0
    for i, st in enumerate(sets):
        gaps = st["gaps"]
        cnt: Dict[int, Optional[int]] = {g[0]: None for g in gaps}
        typ, frm, to, moved, vis = None, None, None, 0, 0
        if i == 0:
            prev_vis = int(rand() * 0.5 * n_targets)
            vis, typ = prev_vis, "INIT"
            b = dob(n_targets - prev_vis, len(gaps))
            for j, g in enumerate(gaps):
                cnt[g[0]] = b[j]
            to = [g[0] for g in gaps]
        else:
            pgaps, pcnt = sets[i - 1]["gaps"], counts[i - 1]
            if len(gaps) > len(pgaps):
                gm = set(cnt)
                for pid, _ie, pto in pgaps:
                    if len(pto) == 2:
                        typ = "SPLIT"
                        b = dob(pcnt[pid], 2)
                        cnt[pto[0]], cnt[pto[1]] = b[0], b[1]
                        frm, to = [pid], [pto[0], pto[1]]
                    else:
                        if pid not in cnt:
                            raise KeyError(f"set {i}: NullPointerException (gap {pid} missing)")
                        cnt[pid] = pcnt[pid]
                    gm.discard(pid)
                if len(gm) == 1:
                    typ = "APPEAR"
                    gn = next(iter(gm))
                    b = dob(prev_vis, 2)
                    cnt[gn] = b[0]
                    prev_vis = b[1]
                    vis, to, moved = prev_vis, [gn], b[0]
                else:
                    vis = prev_vis
            else:
                for pid, _ie, pto in pgaps:
                    if len(pto) == 1:
                        t = pcnt[pid]
                        tgt = pto[0]
                        if cnt[tgt] is None:
                            typ, vis = "MERGE", prev_vis
                            cnt[tgt] = t
                            to, frm = [tgt], [pid, 0]
                        else:
                            cnt[tgt] = (t + t) if merge_bug else (cnt[tgt] + t)
                            frm[1] = pid
                    elif pid not in cnt:
                        typ = "DISAPPEAR"
                        ng = pcnt[pid]
                        prev_vis += ng
                        vis, frm, moved = prev_vis, [pid], ng
                    else:
                        cnt[pid] = pcnt[pid]
        counts.append(cnt)
        rows.append([st["t"], typ, frm, to, moved, vis])
    # TreeMap<Double, event>: sorted by time, equal times overwrite
    by_t: Dict[float, list] = {}
    for r in rows:
        by_t[r[0]] = r
    return [by_t[t] for t in sorted(by_t)]


def format_event(row: list) -> str:
    """``SingleTypeAgentEvent.toString()`` for a dumped event row."""
    _t, typ, frm, to, moved, _vis = row[:6]
    if typ == "INIT":
        return "INIT EVENT" + "".join(f" {g}" for g in to)
    if typ == "APPEAR":
        return f"APPEAR EVENT {to[0]} <- {moved}"
    if typ == "DISAPPEAR":
        return f"DISAPPEAR EVENT {frm[0]} -> {moved}"
    if typ == "SPLIT":
        return f"SPLIT EVENT {frm[0]} -> {to[0]} {to[1]}"
    if typ == "MERGE":
        return f"MERGE EVENT {frm[0]} {frm[1]} -> {to[0]}"
    raise ValueError(f"event type {typ}")  # Java: NullPointerException in the switch


def format_equations(orun: dict) -> List[str]:
    """``deriveGapEvolvingEquations`` + ``Equation.toString()`` from a dumped oracle run."""
    out = []
    for i, (_t, typ, frm, to, moved, vis) in enumerate(e[:6] for e in orun["events"]):
        if i == 0:
            terms, const = {g: 1 for g in to}, orun["n_targets"] - vis
        elif typ in ("APPEAR", "DISAPPEAR"):
            terms, const = {(to if typ == "APPEAR" else frm)[0]: 1}, moved
        elif typ == "SPLIT":
            terms, const = {to[0]: 1, to[1]: 1}, 0
            terms[frm[0]] = -1
        elif typ == "MERGE":
            terms, const = {to[0]: 1}, 0
            terms[frm[0]] = -1
            terms[frm[1]] = -1
        else:
            raise ValueError(f"event type {typ}")
        ks = sorted(terms)
        buf = ""
        for j, k in enumerate(ks):
            buf += f"{terms[k]} x{k}"
            if j != len(ks) - 1:
                buf += " + " if terms[ks[j + 1]] >= 0 else " "
        out.append(buf + f" = {const}")
    return out


class _JavaRandomStream:
    """Incremental ``java.util.Random`` (48-bit LCG) producing ``nextDouble()``."""

    def __init__(self, seed: int) -> None:
        self.s = (seed ^ 0x5DEECE66D) & ((1 << 48) - 1)

    def _next(self, bits: int) -> int:
        self.s = (self.s * 0x5DEECE66D + 0xB) & ((1 << 48) - 1)
        return self.s >> (48 - bits)

    def next_double(self) -> float:
        return ((self._next(26) << 27) + self._next(27)) * (1.0 / (1 << 53))


# ----------------------------------------------------------------------------- events
def gap_history_events(run: dict) -> Tuple[List[int], List[List[object]], List[str]]:
    """Component events read off the Java gap history, as described in docs/notes/original_java.md §2.7.

    Returns ``(initial_ids, events_per_transition, anomalies)``. ``events_per_transition[k]``
    lists the events between gap set ``k`` and ``k + 1``: Split / Merge, plus Appear /
    Disappear with zero counts, since the gap history carries no counts. The Java code
    assumes exactly one event per transition; anything else is reported as an anomaly.
    """
    sets = run["get_gaps"]["sets"]
    ids = [[g[0] for g in s["gaps"]] for s in sets]
    out: List[List[object]] = []
    anomalies: List[str] = []
    for k in range(len(sets) - 1):
        cur, nxt = sets[k]["gaps"], set(ids[k + 1])
        evs: List[object] = []
        targets = set()
        merges: Dict[int, List[int]] = defaultdict(list)
        for gid, _iedge, to in cur:
            if len(to) == 2:
                evs.append(Split(gid, to[0], to[1]))
                targets.update(to)
            elif len(to) == 1:
                merges[to[0]].append(gid)
                targets.add(to[0])
            elif len(to) > 2:
                anomalies.append(f"set {k}: gap {gid} links to {to} (>2 children)")
                targets.update(to)
        for s, parents in merges.items():
            if len(parents) == 2:
                evs.append(Merge(parents[0], parents[1], s))
            else:
                anomalies.append(f"set {k}: {parents} -> {s} ({len(parents)} merging parents)")
        cur_ids = set(ids[k])
        for gid in ids[k + 1]:
            if gid not in cur_ids and gid not in targets:
                evs.append(Appear(gid, 0, 0))
        for gid, _iedge, to in cur:
            if not to and gid not in nxt:
                evs.append(Disappear(gid, 0, 0))
            if to and gid in nxt:
                anomalies.append(f"set {k}: linked gap {gid} survives")
        if len(evs) != 1:
            anomalies.append(f"set {k}->{k + 1}: {len(evs)} events {evs}")
        out.append(evs)
    return ids[0], out, anomalies


def _shape(e) -> tuple:
    """Event without counts, with split children / merge parents unordered."""
    if isinstance(e, Split):
        return ("split", e.s, frozenset((e.a, e.b)))
    if isinstance(e, Merge):
        return ("merge", frozenset((e.a, e.b)), e.s)
    if isinstance(e, Appear):
        return ("appear", e.s)
    if isinstance(e, Disappear):
        return ("disappear", e.s)
    raise TypeError(e)


def oracle_sequence(orun: dict) -> Tuple[ShadowSequence, int, List[float], List[str]]:
    """``ShadowSequence`` for one Java oracle run, plus the known hidden total ``N - visible0``.

    Initial shadows get ``[0, inf)``; the oracle's observations are exact:
    ``Appear(s, moved, moved)``, ``Disappear(s, moved, moved)``. The total is passed to
    :class:`CombinatorialFilter` separately (``total=(h, h)``), which is how the Java
    code's pooled left vertex 0 encodes it.
    """
    evs = orun["events"]
    t0, typ, _frm, to, _moved, vis = evs[0][:6]
    assert typ == "INIT", evs[0]
    hidden = orun["n_targets"] - vis
    events, times, problems = [], [], []
    for t, typ, frm, to, moved, _vis in (e[:6] for e in evs[1:]):
        times.append(t)
        if typ == "APPEAR":
            events.append(Appear(to[0], moved, moved))
        elif typ == "DISAPPEAR":
            events.append(Disappear(frm[0], moved, moved))
        elif typ == "SPLIT":
            events.append(Split(frm[0], to[0], to[1]))
        elif typ == "MERGE":
            events.append(Merge(frm[0], frm[1], to[0]))
        else:
            problems.append(f"event at t={t}: type {typ}")
    seq = ShadowSequence({g: (0, INF) for g in to_initial(evs[0])}, events)
    return seq, hidden, times, problems


def to_initial(init_event: list) -> List[int]:
    return list(init_event[3])


# ----------------------------------------------------------------------------- bipartite graph
def java_bg_edges(bg: dict) -> set:
    return {(lid, rid) for lid, _lo, _hi, tos in bg["left"] for rid in tos}


def shadowinfo_bg(seq: ShadowSequence) -> BipartiteIState:
    return BipartiteIState.from_sequence(seq)


def compare_bipartite(seq: ShadowSequence, hidden: int, bg: dict) -> dict:
    """Java BG (pooled initial vertex 0) vs ``BipartiteIState`` (one left vertex per initial shadow).

    Java quirk modelled here: an appeared shadow is a single ``Vertex`` object that is both
    its own left and right vertex (self-loop). While it never splits or merges, the right
    vertex therefore shows the appear count instead of -1 (alive/unknown), and a disappear
    overwrites the shared weight with the disappear count.
    """
    st = shadowinfo_bg(seq)
    init = set(st.initial)
    pooled = lambda u: 0 if u in init else u  # noqa: E731
    si_edges = {(pooled(u), v) for u, (_tag, v) in st.edges()}
    jv_edges = java_bg_edges(bg)
    left_w = {lid: hi for lid, _lo, hi, _ in bg["left"]}
    right_w = {rid: hi for rid, _lo, hi in bg["right"]}
    si_left_w = {0: hidden}
    si_left_w.update({v: (st.disappeared[v][1] if v in st.disappeared else b[1])
                      for v, b in st.left.items() if v not in init})
    si_right_w = {v: b[1] for v, b in st.disappeared.items()}
    si_right_w.update({s: (st.left[s][1] if s in st.left and s not in init else -1) for s in st.reach})
    res = {
        "edges_equal": si_edges == jv_edges,
        "left_weights_equal": left_w == si_left_w,
        "right_weights_equal": right_w == si_right_w,
        "java": {"n_left": bg["n_left"], "n_right": bg["n_right"], "n_edges": bg["n_edges"]},
        "shadowinfo": {"n_left": len(st.left), "n_right": len(st.disappeared) + len(st.reach),
                       "n_edges": len(st.edges()), "n_edges_pooled": len(si_edges)},
    }
    if not res["edges_equal"]:
        res["only_java"] = sorted(jv_edges - si_edges)
        res["only_shadowinfo"] = sorted(si_edges - jv_edges)
    if not res["left_weights_equal"]:
        res["left_weight_diff"] = {k: (left_w.get(k), si_left_w.get(k)) for k in set(left_w) | set(si_left_w)
                                   if left_w.get(k) != si_left_w.get(k)}
    if not res["right_weights_equal"]:
        res["right_weight_diff"] = {k: (right_w.get(k), si_right_w.get(k)) for k in set(right_w) | set(si_right_w)
                                    if right_w.get(k) != si_right_w.get(k)}
    return res


# ----------------------------------------------------------------------------- bounds
def exact_bounds(seq: ShadowSequence, hidden: int):
    """Exact bounds for every final shadow, or ``None`` if the observations are infeasible."""
    f = CombinatorialFilter.from_sequence(seq, total=(hidden, hidden))
    try:
        return {s: (int(lo), hi if hi == INF else int(hi)) for s, (lo, hi) in f.all_bounds().items()}
    except InfeasibleError:
        return None


def compare_bounds(seq: ShadowSequence, hidden: int, orun: dict, final_ids: List[int]) -> dict:
    exact = exact_bounds(seq, hidden)
    truth = dict(zip(final_ids, orun["truth_final"]))
    res: dict = {"feasible": exact is not None}
    outs = orun.get("java_bounds_outcomes", [])
    res["n_java_outcomes"] = len(outs)
    if exact is None:
        res["java_outcomes"] = [o["bounds"] for o in outs]
        return res
    res["exact"] = {str(s): list(b) for s, b in sorted(exact.items())}
    res["truth_in_exact"] = all(exact[s][0] <= truth[s] <= exact[s][1] for s in exact)
    res["truth_not_in_exact"] = {str(s): (truth[s], exact[s]) for s in exact if not exact[s][0] <= truth[s] <= exact[s][1]}
    bounded = set()
    per_outcome = []
    for o in outs:
        jb = {int(k): tuple(v) for k, v in o["bounds"].items()}
        bounded |= set(jb)
        diff = {str(s): {"java": list(jb[s]), "exact": list(exact[s])} for s in sorted(jb) if jb[s] != exact.get(s)}
        kinds = Counter(_classify(jb[s], exact[s]) for s in jb)
        per_outcome.append({"burnins": o["burnins"], "equal_exact": not diff, "diff": diff, "kinds": dict(kinds),
                            "unsound": sorted(s for s in jb if not (jb[s][0] <= truth[s] <= jb[s][1]))})
    res["outcomes"] = per_outcome
    res["any_outcome_equals_exact"] = any(p["equal_exact"] for p in per_outcome)
    res["all_outcomes_equal_exact"] = bool(per_outcome) and all(p["equal_exact"] for p in per_outcome)
    res["shadows_without_java_bound"] = sorted(set(exact) - bounded)
    return res


def _classify(j: Tuple[int, int], e: Tuple[int, float]) -> str:
    """How a Java [min, max] relates to the exact [lo, hi]."""
    if tuple(j) == tuple(e):
        return "equal"
    if j[0] > j[1]:
        return "min>max"
    if j[0] < 0:
        return "negative_min"
    if j[0] <= e[0] and j[1] >= e[1]:
        return "looser"
    return "excludes_feasible_values"


# ----------------------------------------------------------------------------- NSE
def compare_nse(run: dict) -> dict:
    """Java NSE labels vs 'contaminated iff the reach set contains an initial shadow'.

    The latter is the counting filter: initial shadows ``[0, inf)``, appeared and
    disappeared shadows hold 0, so the upper bound of a shadow is ``inf`` iff it can
    contain a target from ``t0``.
    """
    init, per, anomalies = gap_history_events(run)
    st = BipartiteIState({g: (0, INF) for g in init})
    states = run["nse"]["states"]
    sets = run["get_gaps"]["sets"]
    mism = []
    for k, s in enumerate(sets):
        ids = [g[0] for g in s["gaps"]]
        want = "".join("1" if (st.reach[g] & set(st.initial)) else "0" for g in ids)
        if want != states[k]:
            mism.append({"set": k, "java": states[k], "reach_rule": want, "ids": ids})
        if k < len(per):
            for e in per[k]:
                st.apply(e)
    return {"equal": not mism, "mismatches": mism[:5], "n_mismatched_sets": len(mism)}


# ----------------------------------------------------------------------------- paper
def paper_counts(run: dict, orun: dict) -> dict:
    bg = orun["bg"]
    seq, hidden, _t, _p = oracle_sequence(orun)
    st = shadowinfo_bg(seq)
    self_loop = [lid for lid, _lo, _hi, tos in bg["left"] if lid in tos]
    final = run["get_gaps"]["final_ids"]
    return {
        "component_events": len(run["get_gaps"]["sets"]) - 1,
        "shadow_ids": run["get_gaps"]["max_id"],
        "final_shadows": len(final),
        "final_shadows_bounded_by_java": len([g for g in final if g not in self_loop]),
        "java_bg_left": bg["n_left"], "java_bg_right": bg["n_right"],
        "java_bg_vertices": bg["n_left"] + bg["n_right"], "java_bg_edges": bg["n_edges"],
        "java_bg_vertices_distinct_ids": len({v[0] for v in bg["left"]} | {v[0] for v in bg["right"]}),
        "shadowinfo_bg_vertices": len(st.left) + len(st.disappeared) + len(st.reach),
        "shadowinfo_bg_edges": len(st.edges()),
        "initial_shadows": len(st.initial),
        "appeared": sum(isinstance(e, Appear) for e in seq.events),
        "disappeared": sum(isinstance(e, Disappear) for e in seq.events),
        "splits": sum(isinstance(e, Split) for e in seq.events),
        "merges": sum(isinstance(e, Merge) for e in seq.events),
    }


# ----------------------------------------------------------------------------- per run
def check_run(run: dict) -> dict:
    rep: dict = {"name": run["name"], "polygon": run["polygon"], "get_gaps_ok": run["get_gaps"]["ok"],
                 "n_critical_points": run["n_critical_points"], "n_gi_bt": run["n_gi_bt_crossings"]}
    if not run["get_gaps"]["ok"]:
        rep["failure"] = {"exception": run["get_gaps"]["exception"]["class"], **run["get_gaps"]["failure"]}
        return rep
    gg = run["get_gaps"]
    rep["n_sets"] = gg["n_sets"]
    rep["sets_equal_samples"] = gg["sets_equal_samples"]
    rep["getGaps_repeatable"] = gg["repeat_calls_identical"]
    init, per, anomalies = gap_history_events(run)
    rep["gap_history_anomalies"] = anomalies
    rep["nse"] = compare_nse(run)
    flat = [e for evs in per for e in evs]
    try:
        ShadowSequence({g: (0, INF) for g in init}, flat)
        rep["gap_history_sequence_valid"] = True
    except Exception as ex:  # noqa: BLE001
        rep["gap_history_sequence_valid"] = repr(ex)
    oruns = []
    bg_shapes = set()
    for orun in run["oracle_runs"]:
        o: dict = {"variant": orun["variant"], "seed": orun["seed"]}
        if "initialize_exception" in orun:
            o["initialize_exception"] = orun["initialize_exception"]["class"]
            oruns.append(o)
            continue
        seq, hidden, times, problems = oracle_sequence(orun)
        o["hidden_total"] = hidden
        o["problems"] = problems
        o["n_events"] = orun["n_events"]
        o["events_equal_gap_history"] = ([_shape(e) for e in seq.events] == [_shape(e) for e in flat]
                                         and to_initial(orun["events"][0]) == init)
        o["times_equal_set_times"] = times == [s["t"] for s in gg["sets"][1:]]
        o["bg"] = compare_bipartite(seq, hidden, orun["bg"])
        bg_shapes.add(json.dumps(sorted(java_bg_edges(orun["bg"]))))
        o["bounds"] = compare_bounds(seq, hidden, orun, gg["final_ids"])
        oruns.append(o)
    rep["oracle_runs"] = oruns
    rep["bg_structure_seed_independent"] = len(bg_shapes) <= 1
    if run["oracle_runs"] and "bg" in run["oracle_runs"][0]:
        rep["counts"] = paper_counts(run, run["oracle_runs"][0])
    return rep


def icra08_fig4_check(run: dict) -> dict:
    """The polygon-12 ProjectPanel run against tests/fixtures/icra08_fig4.json (labels only).

    The paper treats s5 as initial; in Java s5 appears at t = 0.056. The paper also
    orders two pairs of independent events differently, so the comparison is on the
    final bipartite graph rather than on the event order.
    """
    fig = json.loads((ROOT / "tests" / "fixtures" / "icra08_fig4.json").read_text())
    seq, _h, _t, _p = oracle_sequence(run["oracle_runs"][0])
    ev = [e for e in seq.events if not (isinstance(e, Appear) and e.s == 5)]
    s2 = ShadowSequence({g: (0, INF) for g in list(seq.initial) + [5]}, ev)
    st = BipartiteIState.from_sequence(s2)
    edges = sorted([u, v] for u, (_t2, v) in st.edges())
    paper = fig["sequence"]["events"]
    ours = [_shape(e) for e in ev]
    theirs = [_shape(e) for e in ShadowSequence.from_dict(fig["sequence"]).events]
    return {"bipartite_edges_equal": edges == sorted(fig["bipartite"]["edges"]),
            "same_event_multiset": Counter(ours) == Counter(theirs),
            "same_event_order": ours == theirs, "n_paper_events": len(paper)}


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="*", help="fixture names (default: all runs)")
    ap.add_argument("--json", help="write the full report here")
    a = ap.parse_args(argv)
    files = [FIXTURES / f"{r}.json" for r in a.runs] if a.runs else run_files()
    report = {}
    for f in files:
        run = json.loads(f.read_text())
        rep = check_run(run)
        if run["name"] == "PP" and run["polygon"] == 12 and run["get_gaps"]["ok"]:
            rep["icra08_fig4"] = icra08_fig4_check(run)
        report[f.stem] = rep
        _print(f.stem, rep)
    if a.json:
        Path(a.json).write_text(json.dumps(report, indent=1, default=list))
    return 0


def _print(name: str, rep: dict) -> None:
    if not rep["get_gaps_ok"]:
        fl = rep["failure"]
        print(f"{name}: getGaps FAILS ({fl['exception']}) at critical point {fl.get('critical_index')} "
              f"({fl.get('critical_type')}; {fl.get('phase')})")
        return
    print(f"{name}: {rep['n_sets'] - 1} component events, sets==samples {rep['sets_equal_samples']}, "
          f"history anomalies {len(rep['gap_history_anomalies'])}, sequence valid {rep['gap_history_sequence_valid']}, "
          f"NSE==reach-rule {rep['nse']['equal']}")
    if "counts" in rep:
        print(f"   counts {rep['counts']}")
    for o in rep["oracle_runs"]:
        if "initialize_exception" in o:
            print(f"   {o['variant']:11s} seed {o['seed']}: oracle.initialize throws {o['initialize_exception']}")
            continue
        b = o["bounds"]
        bgok = o["bg"]["edges_equal"] and o["bg"]["left_weights_equal"] and o["bg"]["right_weights_equal"]
        line = (f"   {o['variant']:11s} seed {o['seed']}: events==history {o['events_equal_gap_history']}, "
                f"BG==BipartiteIState {bgok}, feasible {b['feasible']}, java outcomes {b['n_java_outcomes']}")
        if b["feasible"]:
            line += (f", any==exact {b['any_outcome_equals_exact']}, all==exact {b['all_outcomes_equal_exact']}, "
                     f"unsound outcomes {sum(1 for p in b['outcomes'] if p['unsound'])}, "
                     f"no java bound for {b['shadows_without_java_bound']}")
        print(line)
    if "icra08_fig4" in rep:
        print(f"   ICRA'08 Fig. 4: {rep['icra08_fig4']}")


if __name__ == "__main__":
    raise SystemExit(main())
