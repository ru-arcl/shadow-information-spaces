"""Command line for the polygon front end: what the original Java panels print, without the applets.

Usage::

    python -m shadowinfo.polygon list                       # the 14 maps and every stored robot path
    python -m shadowinfo.polygon cuts N [--summary]         # critical lines of map N (getCuts order)
    python -m shadowinfo.polygon vis N X Y [--compat ...]   # physical gaps and visibility polygon at (X, Y)
    python -m shadowinfo.polygon show [N] [--at X Y] [--save FILE]   # live viewer: the robot follows the mouse
    python -m shadowinfo.polygon run TARGET [--path NAME] [--seed S] [--targets N] [--compat java|fixed]
                                            [--merge-bug yes|no] [--bounds-order sorted|reverse|INT]
                                            [--report panel|nse]

``run`` reproduces the stdout of the ``ProjectPanel*`` constructors (``docs/notes/original_java.md``
§2.11).  ``TARGET`` is a panel class (``ProjectPanel`` = path ``P12`` with 100 targets, ``ProjectPanel2``
/ ``ProjectPanel3`` = never-see-evader on ``P1`` / ``P5``, ``ProjectPanel4`` prints nothing,
``ProjectPanel5`` = ``P14b`` with 10^6 targets), a path label (``P13``, ``fig_TRO-Fig15b``, ... see
``list``) or a map number (with ``--path``; default: the panel path on that map).

* ``--compat java`` reproduces the Java output: gap tracking with its crashes, the oracle with the merge
  bug (B7) and Java's max-flow "bounds" (``--bounds-order``: ``HashMap`` order, B10).  With ``--seed S`` the
  oracle draws from ``new java.util.Random(S)``, which is what the Java prints after ``Math.random()`` is
  seeded with ``S`` (the golden fixtures' ``oracle_runs``).
* ``--compat fixed`` (default) tracks gaps robustly, conserves targets and prints the exact bounds of every
  shadow alive at the end in the same line format.  A path that leaves the polygon (``P13+``) is still tracked,
  with a warning on stderr (:func:`path_warning`).
* Without ``--seed`` the run is not reproducible, like ``Math.random()``.

A Java exception (``--compat java``) ends the output where the Java console output ends; the error goes
to stderr and the exit status is 1.
"""

from __future__ import annotations

import argparse
import sys
from typing import Iterable, Optional, Sequence, TextIO, Tuple

from .cuts import (Bitangent, GeneralInflection, InflectionType, all_cuts, bitangent_cuts, general_inflection_cuts,
                   inflections, single_tangent_cuts)
from .events import never_see_evader
from .gaps import GapTrackingError, format_gap_sets, get_gaps
from .geometry import COMPAT_MODES, GeometryError
from .io import N_POLYGONS, PANELS, java_double_str, load_path, load_paths, load_polygon, path_record
from .oracle import OracleError, SingleTypeAgentOracle
from .stagent import project_panel_lines
from .visibility import physical_gaps, visibility_polygon

PANEL_TARGETS = {"ProjectPanel": 100, "ProjectPanel2": 100, "ProjectPanel3": 100, "ProjectPanel4": 1000000,
                 "ProjectPanel5": 1000000}
"""``SingleTypeAgentOracle(N)`` of each panel (only ``ProjectPanel`` and ``ProjectPanel5`` run an oracle)."""

PANEL_REPORTS = {"ProjectPanel": "panel", "ProjectPanel2": "nse", "ProjectPanel3": "nse", "ProjectPanel4": "none",
                 "ProjectPanel5": "panel"}
"""What each panel prints: the oracle/filter report, never-see-evader gap sets, or nothing (drawing only)."""


def _fmt_point(p) -> str:
    return "null" if p is None else f"({java_double_str(p[0])}, {java_double_str(p[1])})"


def resolve_target(target: str, path: Optional[str] = None) -> Tuple[str, str, int, str]:
    """``(panel or "", path label, default number of targets, report)`` for a ``run`` target."""
    labels = load_paths()
    if target in PANELS:
        label = path or PANELS[target]
        return target, path_label(label, labels), PANEL_TARGETS[target], PANEL_REPORTS[target]
    if target.isdigit():
        n = int(target)
        if not 1 <= n <= N_POLYGONS:
            raise SystemExit(f"map number must be 1..{N_POLYGONS}")
        if path is None:
            for panel, lab in PANELS.items():
                if labels[lab]["polygon"] == n and PANEL_REPORTS[panel] != "none":
                    path = lab
                    break
            else:
                cands = [k for k, r in labels.items() if r["polygon"] == n]
                if not cands:
                    raise SystemExit(f"no stored path on map {n}; give --path (see 'list')")
                path = cands[0]
        label = path_label(path, labels)
        if labels[label]["polygon"] != n:
            raise SystemExit(f"path {label} is on map {labels[label]['polygon']}, not {n}")
    else:
        label = path_label(target if path is None else path, labels)
    rec = labels[label]
    return "", label, rec["harness"]["n_targets"], "panel"


def path_label(name: str, labels: Optional[dict] = None) -> str:
    labels = labels if labels is not None else load_paths()
    if name in labels:
        return name
    try:
        rec = path_record(name)
    except KeyError as ex:
        raise SystemExit(str(ex.args[0])) from None
    for k, r in labels.items():
        if r == rec:
            return k
    raise SystemExit(f"unknown path {name!r}")  # pragma: no cover


def _parse_order(s: str):
    if s in ("sorted", "reverse"):
        return s
    try:
        return int(s)
    except ValueError:
        raise argparse.ArgumentTypeError("--bounds-order must be 'sorted', 'reverse' or an integer seed") from None


def run_lines(target: str, path: Optional[str] = None, seed: Optional[int] = None, targets: Optional[int] = None,
              compat: str = "fixed", merge_bug: Optional[bool] = None, bounds_order="sorted",
              report: Optional[str] = None) -> Iterable[str]:
    """The lines ``run`` prints (a generator; exceptions propagate after the lines printed before them)."""
    _panel, label, n_default, rep = resolve_target(target, path)
    rep = report or rep
    if rep == "none":
        return
    poly, pth = load_path(label)
    history = get_gaps(poly, pth, compat)
    if rep == "nse":
        never_see_evader(history, compat)
        yield from format_gap_sets(history, with_time=False)
        yield ""
        yield ""
        return
    oracle = SingleTypeAgentOracle(n_default if targets is None else targets, seed, merge_bug=merge_bug,
                                   compat=compat)
    yield from project_panel_lines(history, oracle, bounds_order)


def path_warning(target: str, path: Optional[str] = None) -> Optional[str]:
    """Why the ``run`` path is not a valid robot path (:func:`~shadowinfo.polygon.simulate.validate_path`),
    or None.  ``compat="fixed"`` still tracks such a path (samples outside are pulled back inside), so
    ``run`` prints this to stderr as a warning; Java's ``getGaps`` usually throws on it (``P13+``)."""
    from .simulate import validate_path

    _panel, label, _n, rep = resolve_target(target, path)
    if rep == "none":
        return None
    poly, pth = load_path(label)
    try:
        validate_path(poly, pth.points)
    except ValueError as ex:
        return f"path {label}: {ex}"
    return None


def cmd_run(a: argparse.Namespace, out: TextIO) -> int:
    mb = None if a.merge_bug is None else a.merge_bug == "yes"
    if a.compat == "fixed":
        warning = path_warning(a.target, a.path)
        if warning:
            print(f"warning: {warning}; the output is not a valid run", file=sys.stderr)
    try:
        for line in run_lines(a.target, a.path, a.seed, a.targets, a.compat, mb, a.bounds_order, a.report):
            out.write(line + "\n")
    except (GapTrackingError, OracleError) as ex:
        out.flush()
        exc = getattr(ex, "java_exception", "") or type(ex).__name__
        print(f"{exc}: {ex}", file=sys.stderr)
        return 1
    return 0


def cmd_list(a: argparse.Namespace, out: TextIO) -> int:
    out.write("maps (shadowinfo/polygon/data/N.dat; y-up, counter-clockwise):\n")
    for n in range(1, N_POLYGONS + 1):
        p = load_polygon(n)
        uses = sorted(k for k, r in load_paths().items() if r["polygon"] == n)
        out.write(f"  {n:2d}  vertices {p.n:3d}  reflex {len(p.reflex):3d}  paths: {', '.join(uses) or '-'}\n")
    out.write("paths (label  map  kind  waypoints  panels  figure):\n")
    for k, r in load_paths().items():
        panels = ",".join(r.get("panels") or []) or "-"
        fig = r.get("figure") or "-"
        if len(fig) > 60:
            fig = fig[:57] + "..."
        ok = "" if r["harness"].get("java_get_gaps_ok", True) else "  [Java getGaps throws]"
        out.write(f"  {k:18s} {r['polygon']:2d}  {r['kind']:9s} {len(r['waypoints']):2d}  {panels:22s} {fig}{ok}\n")
    out.write("panels: " + ", ".join(f"{p} = {lab}" for p, lab in PANELS.items()) + "\n")
    return 0


def _cut_extra(c) -> str:
    if isinstance(c, GeneralInflection):
        return f" fromLine={c.from_line} ccw={str(c.counter_clockwise).lower()}"
    if isinstance(c, Bitangent):
        return f" this={c.this_point} opp={c.opposite_point} curveTo={_fmt_point(c.curve_to_point)}"
    return ""


def cmd_cuts(a: argparse.Namespace, out: TextIO) -> int:
    poly = load_polygon(a.map)
    st = single_tangent_cuts(poly)
    ng = inflections(poly, InflectionType.NONGENERAL)
    gi = general_inflection_cuts(poly)
    bt = bitangent_cuts(poly)
    out.write(f"map {a.map}: {poly.n} vertices, {len(poly.reflex)} reflex\n")
    out.write(f"single tangents {len(st)}, non-general inflections {len(ng)}, general inflections {len(gi)}, "
              f"bitangent rays {len(bt)}; getCuts {len(st) + len(ng) + len(gi) + len(bt)}\n")
    if a.summary:
        return 0
    for i, c in enumerate(all_cuts(poly)):
        x1, y1, x2, y2 = c.line
        out.write(f"{i} {c.type.name} {java_double_str(x1)} {java_double_str(y1)} {java_double_str(x2)} "
                  f"{java_double_str(y2)}{_cut_extra(c)}\n")
    return 0


def cmd_vis(a: argparse.Namespace, out: TextIO) -> int:
    poly = load_polygon(a.map)
    q = (a.x, a.y)
    out.write(f"point {_fmt_point(q)} in map {a.map}: {'inside' if poly.contains(q) else 'OUTSIDE'}\n")
    try:
        gaps = physical_gaps(poly, q, compat=a.compat)
    except GeometryError as ex:  # compat="fixed" refuses a point outside the polygon
        print(f"{type(ex).__name__}: {ex}", file=sys.stderr)
        return 1
    out.write(f"physical gaps {len(gaps)} (startEdge, endEdge, startPoint, endPoint):\n")
    for g in gaps:
        out.write(f"  [{g.start_edge}, {g.end_edge}, {_fmt_point(g.start_point)}, {_fmt_point(g.end_point)}]\n")
    vp = visibility_polygon(poly, q, a.compat)
    out.write(f"visibility polygon {len(vp)} vertices:\n")
    for p in vp:
        out.write(f"  {java_double_str(p[0])} {java_double_str(p[1])}\n")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="python -m shadowinfo.polygon",
                                 description="Polygon front end: the original Java implementation without applets.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="print what a ProjectPanel prints (gap sets, bipartite graph, bounds, ...)")
    r.add_argument("target", help="panel class (ProjectPanel..ProjectPanel5), path label or map number")
    r.add_argument("--path", help="path label (see 'list'); replaces the panel's path")
    r.add_argument("--seed", type=int, help="seed of java.util.Random (default: unseeded, not reproducible)")
    r.add_argument("--targets", type=int, help="number of targets N (default: the panel's)")
    r.add_argument("--compat", choices=COMPAT_MODES, default="fixed")
    r.add_argument("--merge-bug", choices=("yes", "no"), help="oracle merge bug B7 (default: yes for java)")
    r.add_argument("--bounds-order", type=_parse_order, default="sorted",
                   help="HashMap order of Java's max-flow (compat java): sorted, reverse or an int hash seed")
    r.add_argument("--report", choices=("panel", "nse", "none"),
                   help="oracle report (ProjectPanel/5), never-see-evader sets (ProjectPanel2/3) or nothing")
    sub.add_parser("list", help="the maps and the stored robot paths")
    c = sub.add_parser("cuts", help="critical lines of a map")
    c.add_argument("map", type=int)
    c.add_argument("--summary", action="store_true", help="counts only")
    v = sub.add_parser("vis", help="physical gaps and visibility polygon at a point")
    v.add_argument("map", type=int)
    v.add_argument("x", type=float)
    v.add_argument("y", type=float)
    v.add_argument("--compat", choices=COMPAT_MODES, default="fixed",
                   help="java keeps the duplicate vertices of getVisibilityPolygon (B15)")
    s = sub.add_parser("show", help="live viewer: visibility region follows the mouse (needs matplotlib)")
    s.add_argument("map", type=int, nargs="?", default=12)
    s.add_argument("--at", type=float, nargs=2, metavar=("X", "Y"), help="initial robot position")
    s.add_argument("--save", metavar="FILE", help="write one frame to FILE instead of opening a window")
    return ap


def cmd_show(a: argparse.Namespace, out: TextIO) -> int:
    try:
        from .viewer import show
    except ImportError as e:  # matplotlib missing
        print(f"show needs matplotlib ({e}); pip install matplotlib", file=sys.stderr)
        return 2
    try:
        show(a.map, tuple(a.at) if a.at else None, a.save)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 2
    return 0


def main(argv: Optional[Sequence[str]] = None, out: Optional[TextIO] = None) -> int:
    a = build_parser().parse_args(argv)
    out = out if out is not None else sys.stdout
    cmds = {"run": cmd_run, "list": cmd_list, "cuts": cmd_cuts, "vis": cmd_vis, "show": cmd_show}
    if a.cmd in ("cuts", "vis", "show") and not 1 <= a.map <= N_POLYGONS:
        print(f"map number must be 1..{N_POLYGONS}", file=sys.stderr)
        return 2
    return cmds[a.cmd](a, out)


if __name__ == "__main__":
    sys.exit(main())
