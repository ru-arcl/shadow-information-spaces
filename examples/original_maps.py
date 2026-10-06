"""The 14 maps of the original Java implementation with every stored robot path and its provenance.

The maps are ``shadowinfo/polygon/data/{1..14}.dat``, byte-identical copies of ``source/polygons/N.dat``
of https://github.com/arc-l/shadow-information-space (``docs/notes/original_java.md`` §2.2).  The paths
come from ``data/paths.json`` (``docs/notes/original_maps.md``):

* **code** -- waypoints copied from the ``ProjectPanel*`` classes (active or commented out);
* **digitised** -- the path drawn in a paper figure, digitised in the ``.dat`` frame;
* **demo** -- a closed patrol tour designed for the simulator and the browser demo, on the maps the
  Java code never drove a robot through.

Everything is drawn in the ``.dat`` frame, y-up, exactly as the Java window shows it; nothing is
flipped.  The paper draws maps 12 and 14 mirrored (y-down), map 13 y-up.

The script prints one table of the maps (vertices, reflex vertices, critical lines from ``getCuts``,
paths, paper figure) and one of the paths, and writes ``original_maps.png``.

Usage::

    python examples/original_maps.py [--out docs/assets] [--no-plot]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo.polygon import all_cuts, demo_paths, load_paths, load_polygon  # noqa: E402
from shadowinfo.polygon.io import N_POLYGONS  # noqa: E402

FIGURES = {12: "T-RO Fig. 11 · ICRA'08 Figs. 2–4", 13: "T-RO Fig. 15(a) · ICRA'08 Fig. 8(a)",
           14: "T-RO Fig. 15(b) · ICRA'08 Fig. 8(b)"}
"""Paper figures drawn from a map (``original_maps.md``, "Figure → map table")."""
SHORT_FIGURES = {12: "T-RO Fig. 11", 13: "T-RO Fig. 15(a)", 14: "T-RO Fig. 15(b)"}

DRAWN_DIGITISED = ("fig_TRO-Fig11a", "fig_TRO-Fig15a", "fig_TRO-Fig15b")
"""Digitised paths drawn in the figure (the ICRA'08 ones coincide with them within a line width)."""

LEAVES_MAP = ("P13+", "P13+tail")
"""Commented-out code paths that leave their map (``(69, 418)`` is outside 13.dat); not drawn."""

# reference palette (dataviz skill), first three categorical slots: they validate all-pairs
KIND_STYLE = {
    "code": dict(color="#2a78d6", ls="-", label="Java code path"),
    "digitised": dict(color="#eb6834", ls=(0, (4, 2)), label="digitised from the paper"),
    "demo": dict(color="#1baf7a", ls="-", label="demo tour"),
}
SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
FLOOR, WALL = "#f0efec", "#46453f"


def map_rows():
    """One record per map: vertices, reflex vertices, cut counts, stored paths, demo path, figure."""
    paths = load_paths()
    demo = demo_paths()
    rows = []
    for n in range(1, N_POLYGONS + 1):
        poly = load_polygon(n)
        cuts = all_cuts(poly)
        kinds = {}
        for c in cuts:
            kinds[c.type.name] = kinds.get(c.type.name, 0) + 1
        rows.append(dict(n=n, poly=poly, cuts=len(cuts), kinds=kinds,
                         paths=[k for k, r in paths.items() if r["polygon"] == n],
                         demo=demo[str(n)], figure=FIGURES.get(n, "")))
    return rows


def drawn_paths(row, paths):
    """``(label, kind, waypoints)`` drawn on one map: every code path that stays inside the map, the
    T-RO digitised paths and, where the Java code has no path, the demo tour."""
    out = []
    for k in row["paths"]:
        r = paths[k]
        if r["kind"] == "code" and k not in LEAVES_MAP:
            out.append((k, "code", r["waypoints"]))
        elif k in DRAWN_DIGITISED:
            out.append((k, "digitised", r["waypoints"]))
    if row["demo"]["source"] is None:
        out.append(("demo", "demo", row["demo"]["waypoints"]))
    return out


def _path_figure(r: dict) -> str:
    if r["kind"] == "digitised":
        return r["figure"]
    fig = r.get("figure") or ""
    for key, short in (("Fig. 15(b)", "T-RO Fig. 15(b)"), ("Fig. 15(a)", "T-RO Fig. 15(a)"),
                       ("Fig. 4", "ICRA'08 Fig. 4, T-RO Fig. 11(b)")):
        if key in fig:
            return short + (" (v0-v19 drawn)" if fig.startswith("Vertices") else "")
    return "-"


def print_tables(rows) -> None:
    paths = load_paths()
    print("Maps (shadowinfo/polygon/data/N.dat; .dat frame, y-up, counter-clockwise)")
    print(f"{'map':>3}  {'vert':>4}  {'reflex':>6}  {'cuts':>5}  {'ST':>4} {'NGI':>4} {'GI':>4} {'BT':>4}  "
          f"{'code paths':<20}  {'digitised':>9}  figure")
    for r in rows:
        k = r["kinds"]
        code = ", ".join(x for x in r["paths"] if paths[x]["kind"] == "code") or "- (demo tour)"
        dig = sum(paths[x]["kind"] == "digitised" for x in r["paths"])
        print(f"{r['n']:>3}  {r['poly'].n:>4}  {len(r['poly'].reflex):>6}  {r['cuts']:>5}  "
              f"{k.get('SINGLETANGENT', 0):>4} {k.get('NONGENERAL_INFLECTION', 0):>4} "
              f"{k.get('GENERAL_INFLECTION', 0):>4} {k.get('BITANGENT', 0):>4}  {code:<20}  {dig:>9}  "
              f"{r['figure'] or '-'}")
    print("  cuts = getCuts: ST single tangents, NGI non-general inflections, GI general inflections, "
          "BT bitangent rays")
    print()
    print("Paths (data/paths.json; Java getGaps = outcome of the original code on the path)")
    print(f"{'label':<17} {'map':>3}  {'kind':<9} {'pts':>3}  {'source':<53} {'getGaps':<7}  figure")
    for lab, r in paths.items():
        if r["kind"] == "code":
            src = f"{r['status']}: {', '.join(r['panels'])}"
        else:
            pr = r["provenance"]
            src = f"{pr['source_image'].split('.pdf')[0]}.pdf, drawn {pr['orientation_in_figure']}"
        ok = "ok" if r["harness"].get("java_get_gaps_ok", True) else "throws"
        print(f"{lab:<17} {r['polygon']:>3}  {r['kind']:<9} {len(r['waypoints']):>3}  {src:<53} {ok:<7}  "
              f"{_path_figure(r)}")
    tours = [n for n, r in sorted(demo_paths().items(), key=lambda kv: int(kv[0])) if r["source"] is None]
    print(f"demo: closed tours on maps {', '.join(tours)} (provenance 'demo'); the other maps use their code path")


def plot(rows, out: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Polygon as MplPolygon

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 8,
                         "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE, "text.color": INK})
    paths = load_paths()
    ncol, nrow = 5, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(12.5, 8.6))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.01, wspace=0.06, hspace=0.32)
    for ax in axes.flat:
        ax.set_axis_off()
        ax.set_aspect("equal")
        ax.set_xlim(-15, 1015)
        ax.set_ylim(-15, 1015)
    for r in rows:
        ax = axes.flat[r["n"] - 1]
        poly = r["poly"]
        ax.add_patch(MplPolygon(poly.vertices, closed=True, fc=FLOOR, ec=WALL, lw=0.9, joinstyle="round"))
        for lab, kind, wps in drawn_paths(r, paths):
            st = KIND_STYLE[kind]
            xs, ys = [p[0] for p in wps], [p[1] for p in wps]
            ax.plot(xs, ys, color=st["color"], ls=st["ls"], lw=1.6, solid_capstyle="round", zorder=3)
            if kind != "demo":  # start dot and an arrowhead at the end
                ax.plot(xs[0], ys[0], "o", ms=4, color=st["color"], mec=SURFACE, mew=0.8, zorder=4)
                ax.annotate("", xy=(xs[-1], ys[-1]), xytext=(xs[-2], ys[-2]), zorder=4,
                            arrowprops=dict(arrowstyle="-|>", color=st["color"], lw=0, mutation_scale=9,
                                            shrinkA=0, shrinkB=0))
        drawn = [lab for lab, kind, _ in drawn_paths(r, paths) if kind == "code"]
        sub = f"{poly.n} vertices · " + (", ".join(drawn) if drawn else "demo tour")
        ax.set_title(f"Map {r['n']}", loc="left", fontsize=10, fontweight="bold", color=INK, pad=14)
        ax.text(0.0, 1.0, sub, transform=ax.transAxes, fontsize=7.2, color=INK2, va="bottom", ha="left")
        if r["figure"]:
            ax.set_title(SHORT_FIGURES[r["n"]], loc="right", fontsize=8, fontweight="bold", color=INK2, pad=14)
    legend_ax = axes.flat[N_POLYGONS]
    handles = [Line2D([], [], color=s["color"], ls=s["ls"], lw=1.8, label=s["label"]) for s in KIND_STYLE.values()]
    legend_ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=8, handlelength=2.6,
                     borderaxespad=0.2)
    note = ("Dot = start, arrow = end. Demo tours\n"
            "are closed loops, drawn on the maps\n"
            "where the Java code has no path.\n\n"
            "Drawn in the .dat frame (y-up), as in\n"
            "the Java window; the papers draw\n"
            "maps 12 and 14 mirrored. P13+ leaves\n"
            "map 13 and is not drawn.")
    legend_ax.text(0.03, 0.66, note, transform=legend_ax.transAxes, fontsize=7.6, color=INK2, va="top", ha="left",
                   linespacing=1.45)
    fig.suptitle("The 14 maps of the original Java implementation and their robot paths", x=0.01, ha="left",
                 fontsize=12, fontweight="bold", color=INK, y=0.985)
    fig.text(0.01, 0.945, "arc-l/shadow-information-space, source/polygons/{1..14}.dat (BSD-3-Clause); "
                          "paths from shadowinfo/polygon/data/paths.json", fontsize=8, color=MUTED, ha="left")
    out.mkdir(parents=True, exist_ok=True)
    dest = out / "original_maps.png"
    fig.savefig(dest, dpi=130)
    plt.close(fig)
    return dest


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "docs" / "assets", help="folder for original_maps.png")
    ap.add_argument("--no-plot", action="store_true", help="print the tables only")
    a = ap.parse_args(argv)
    rows = map_rows()
    print_tables(rows)
    if not a.no_plot:
        print(f"\nwrote {plot(rows, a.out)}")


if __name__ == "__main__":
    main()
