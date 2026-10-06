"""README figure: the four shadow component events -- appear, disappear, split and merge -- on real geometry.

The robot follows the ``ProjectPanel`` path ``P12`` on map 12 of the original Java code (the running
example of the README and T-RO 2012 Fig. 11).  The gap tracker (:func:`shadowinfo.polygon.get_gaps`)
labels the shadows and turns the changes of the gap set into component events
(:func:`shadowinfo.polygon.transition_events`).  For one event of each type the figure shows the robot
just *before* and just *after* the critical line it crosses (dashed: the inflection ray of a reflex
vertex for appear / disappear, the bitangent for split / merge): free space, the visible region, the
shadow pockets (:func:`shadowinfo.polygon.pocket_polygon`) with their labels, the robot and a piece of
its path.  The shadows taking part in the event are coloured; the others are grey.

Usage::

    python examples/event_illustrations.py [--out docs/assets/events.png] [--path P12]
        [--events appear=2 split=4 disappear=17 merge=21]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402
from matplotlib.patches import Polygon as MplPolygon  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo import Appear, Disappear, Merge, Split  # noqa: E402
from shadowinfo.polygon import (Bitangent, get_gaps, load_path, physical_gaps, pocket_polygon,  # noqa: E402
                                transition_events, visibility_polygon)
from shadowinfo.polygon.viewer import gap_anchor  # noqa: E402

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
WALL, FREE, VISIBLE, OTHER = "#8f8d85", "#ffffff", "#fff1b8", "#dcdbd5"
ROBOT, CUT = "#1a1a1a", "#c0392b"
ROLE = {"old": "#2a78d6", "new": "#eb6834", "new2": "#1baf7a"}
DPI = 150
ASPECT = 1.0
KINDS = ("appear", "disappear", "split", "merge")
DEFAULT_EVENTS = {"appear": 2, "disappear": 15, "split": 4, "merge": 21}
"""Critical point (index into the path's critical points) of the event shown for each type on ``P12``."""

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
                     "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE})


class Run:
    """The gap history of a path, with the arc position of every critical point and labels anywhere."""

    def __init__(self, name: str):
        self.path_name = name
        self.poly, self.path = load_path(name)
        self.history = get_gaps(self.poly, self.path)
        self.samples = self.history.samples
        pre = self.path._prefix
        self.arc = [pre[cp.seg_index] + math.dist(cp.seg[:2], cp.point) for cp in self.samples.critical_points]
        # events at the critical point that starts gap set k + 1
        self.events = {}
        for k, evs in enumerate(transition_events(self.history)):
            i = self.history.sample_indices[k] + 1
            assert self.history.sample_ids[i] == [g.id for g in self.history[k + 1]]
            self.events[i] = evs

    def point_at(self, s: float):
        pts, pre = self.path.points, self.path._prefix
        s = min(max(s, 0.0), pre[-1])
        k = max(j for j in range(len(pts) - 1) if pre[j] <= s)
        (x0, y0), (x1, y1) = pts[k], pts[k + 1]
        t = (s - pre[k]) / max(pre[k + 1] - pre[k], 1e-12)
        return (x0 + t * (x1 - x0), y0 + t * (y1 - y0))

    def shadows(self, s: float):
        """``[(label, pocket)]`` seen from the path point at arc length ``s`` (between critical points)."""
        j = max(i for i, a in enumerate(self.arc) if a <= s)
        by_anchor = {gap_anchor(g): gid for g, gid in zip(self.samples.physical[j], self.history.sample_ids[j])}
        q = self.point_at(s)
        return [(by_anchor[gap_anchor(g)], pocket_polygon(self.poly, g)) for g in physical_gaps(self.poly, q)]

    def kind(self, i: int) -> str:
        (e,) = self.events[i]
        return {Appear: "appear", Disappear: "disappear", Split: "split", Merge: "merge"}[type(e)]

    def window(self, i: int, frac: float = 0.6):
        """Arc positions just before and after event ``i``, at most ``frac`` of the way to the neighbouring events."""
        ev = sorted(self.events)
        k = ev.index(i)
        lo = self.arc[ev[k - 1]] if k else 0.0
        hi = self.arc[ev[k + 1]] if k + 1 < len(ev) else self.path.length()
        a = self.arc[i]
        d = min(60.0, frac * (a - lo), frac * (hi - a))
        return a - d, a + d


def describe(e) -> tuple:
    """Caption and ``{label: role}`` of an event."""
    if isinstance(e, Appear):
        return f"appear: s{e.s} appears", {e.s: "new"}
    if isinstance(e, Disappear):
        return f"disappear: s{e.s} disappears", {e.s: "old"}
    if isinstance(e, Split):
        return f"split: s{e.s} → s{e.a}, s{e.b}", {e.s: "old", e.a: "new", e.b: "new2"}
    return f"merge: s{e.a}, s{e.b} → s{e.s}", {e.a: "old", e.b: "new2", e.s: "new"}


def cut_segments(run: Run, i: int):
    cut = run.samples.critical_points[i].cut
    (x1, y1, x2, y2) = cut.line
    segs = [((x1, y1), (x2, y2))]
    if isinstance(cut, Bitangent):
        segs.append((run.poly.vertices[cut.opposite_point], (x1, y1)))
    return segs


def inner_point(pocket, view):
    """Point of ``pocket`` ∩ view farthest from their boundary, on a coarse grid (no extra dependencies)."""
    from matplotlib.path import Path as MPath
    x0, x1, y0, y1 = view
    mp = MPath(pocket)
    best, bd = None, -1.0
    n = 40
    edges = list(zip(pocket, pocket[1:] + pocket[:1]))
    for a in range(1, n):
        for b in range(1, n):
            p = (x0 + (x1 - x0) * a / n, y0 + (y1 - y0) * b / n)
            if not mp.contains_point(p):
                continue
            d = min(p[0] - x0, x1 - p[0], (p[1] - y0) * 1.0, y1 - p[1])
            for (ax_, ay), (bx, by) in edges:
                vx, vy = bx - ax_, by - ay
                L = vx * vx + vy * vy or 1e-12
                t = max(0.0, min(1.0, ((p[0] - ax_) * vx + (p[1] - ay) * vy) / L))
                d = min(d, math.hypot(p[0] - ax_ - t * vx, p[1] - ay - t * vy))
            if d > bd:
                best, bd = p, d
    return best, bd


def view_box(run: Run, i: int, sa: float, sb: float, roles, aspect: float):
    pts = [run.point_at(sa), run.point_at(sb), run.samples.critical_points[i].point]
    for s in (sa, sb):
        for lab, pk in run.shadows(s):
            if lab in roles:
                pts += pk
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    w, h = x1 - x0, y1 - y0
    pad = 0.08 * max(w, h) + 25
    x0, x1, y0, y1 = x0 - pad, x1 + pad, y0 - pad, y1 + pad
    w, h = x1 - x0, y1 - y0
    if w / h < aspect:
        cx, w = (x0 + x1) / 2, h * aspect
        x0, x1 = cx - w / 2, cx + w / 2
    else:
        cy, h = (y0 + y1) / 2, w / aspect
        y0, y1 = cy - h / 2, cy + h / 2
    return x0, x1, y0, y1


def draw_panel(ax, run: Run, i: int, s: float, sa: float, sb: float, roles, view, title: str):
    poly = run.poly
    x0, x1, y0, y1 = view
    ax.set_facecolor(WALL)
    ax.add_patch(MplPolygon(poly.vertices, closed=True, fc=FREE, ec="none", zorder=0))
    q = run.point_at(s)
    ax.add_patch(MplPolygon(visibility_polygon(poly, q), closed=True, fc=VISIBLE, ec="none", zorder=1))
    span = x1 - x0
    for lab, pk in run.shadows(s):
        role = roles.get(lab)
        color = ROLE[role] if role else OTHER
        ax.add_patch(MplPolygon(pk, closed=True, fc=color, alpha=0.55 if role else 0.9, ec=color if role else "none",
                                lw=1.0, zorder=2))
        p, d = inner_point(pk, view)
        if p is not None and d > 0.012 * span:
            ax.text(*p, f"s{lab}", ha="center", va="center", zorder=7, fontsize=10 if role else 8.5,
                    fontweight="bold" if role else "normal", color=INK if role else MUTED,
                    bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec=color if role else "none", lw=1.2,
                              alpha=0.92) if role else None)
    ax.add_patch(MplPolygon(poly.vertices, closed=True, fc="none", ec=INK2, lw=1.4, zorder=3))
    for a, b in cut_segments(run, i):
        ax.plot([a[0], b[0]], [a[1], b[1]], ls=(0, (4, 3)), color=CUT, lw=1.1, zorder=4)
    # the stretch of path around the event, the heading as an arrow and (after) where the robot was before
    step = 0.09 * span
    ps = [run.point_at(sa - 0.45 * span + k * (sb - sa + 0.65 * span) / 200) for k in range(201)]
    ax.plot([p[0] for p in ps], [p[1] for p in ps], color=INK2, lw=1.2, zorder=5, solid_capstyle="round",
            solid_joinstyle="round")
    if s == sb:
        ax.scatter(*run.point_at(sa), s=70, facecolors="none", edgecolors=INK2, linewidths=1.2, ls=(0, (2, 1.5)),
                   zorder=6)
    a, b = run.point_at(s + 0.6 * step), run.point_at(s + 1.6 * step)
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14, color=INK, lw=1.6, zorder=8,
                                 shrinkA=0, shrinkB=0))
    ax.scatter(*q, s=110, color=ROBOT, edgecolors=SURFACE, linewidths=1.8, zorder=8)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(MUTED)
        sp.set_linewidth(0.8)
    ax.set_title(title, fontsize=9, color=INK2, loc="left", pad=3)


def plot(run: Run, chosen: dict, out: Path) -> Path:
    W, side, gap, mid = 10.67, 0.12, 0.08, 0.32      # inches: width, margins, gap in a pair, gap between pairs
    pw = (W - 2 * side - 2 * gap - mid) / 4          # panel width
    ph = pw / ASPECT
    top, head, sub, rowgap, foot = 0.62, 0.42, 0.2, 0.22, 0.1
    H = top + 2 * (head + sub + ph) + rowgap + foot
    fig = plt.figure(figsize=(W, H))
    for n, kind in enumerate(KINDS):
        i = chosen[kind]
        assert run.kind(i) == kind, (kind, i, run.events[i])
        (e,) = run.events[i]
        caption, roles = describe(e)
        sa, sb = run.window(i)
        view = view_box(run, i, sa, sb, roles, ASPECT)
        r, c = divmod(n, 2)
        x = side + c * (2 * pw + gap + mid)
        y = H - top - r * (head + sub + ph + rowgap) - head - sub - ph
        axes = [fig.add_axes([(x + k * (pw + gap)) / W, y / H, pw / W, ph / H]) for k in range(2)]
        draw_panel(axes[0], run, i, sa, sa, sb, roles, view, "before")
        draw_panel(axes[1], run, i, sb, sa, sb, roles, view, "after")
        ty = (y + ph + sub + 0.06) / H
        fig.text(x / W, ty, caption, fontsize=12.5, fontweight="bold", color=INK, va="bottom")
        line = "a bitangent" if isinstance(run.samples.critical_points[i].cut, Bitangent) else "an inflection ray"
        fig.text((x + 2 * pw + gap) / W, ty, f"robot crosses {line}", fontsize=8.5, color=MUTED, va="bottom",
                 ha="right")
    fig.text(side / W, 1 - 0.14 / H, f"Shadow component events along path {run.path_name} on map {run.poly.name}",
             ha="left", va="top", fontsize=13.5, color=INK)
    fig.text(side / W, 1 - 0.42 / H, "yellow: visible region · coloured: shadows in the event · light grey: other "
             "shadows · dark grey: obstacles · dashed red: the critical line crossed · ring: robot's position before",
             ha="left", va="top", fontsize=8.5, color=INK2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI)
    plt.close(fig)
    return out


def save(out: Path, run=None, chosen=None):
    """Write the figure to ``out`` for ``run`` (default: ``P12``); return ``(out, {kind: critical point})``.

    ``chosen`` overrides the event shown for some types; a type left out gets :data:`DEFAULT_EVENTS` on ``P12``
    and the first event of that type on any other path.
    """
    run = Run("P12") if run is None else run
    chosen = {**(DEFAULT_EVENTS if run.path_name == "P12" else {}), **(chosen or {})}
    for kind in KINDS:
        if kind not in chosen:
            chosen[kind] = next(i for i in sorted(run.events) if len(run.events[i]) == 1 and run.kind(i) == kind)
    return plot(run, chosen, out), chosen


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(ROOT / "docs" / "assets" / "events.png"))
    ap.add_argument("--path", default="P12")
    ap.add_argument("--events", nargs="*", default=[], metavar="KIND=CP",
                    help="critical point index of the event shown for a type (default: P12's choice)")
    ap.add_argument("--list", action="store_true", help="list the path's component events and exit")
    args = ap.parse_args()
    run = Run(args.path)
    if args.list:
        for i, evs in sorted(run.events.items()):
            print(f"cp {i:3d}  arc {run.arc[i]:8.1f}  {run.samples.critical_points[i].cut_type.name:18s} "
                  f"{'; '.join(map(str, evs))}")
        return
    chosen = dict(kv.split("=") for kv in args.events)
    out, chosen = save(Path(args.out), run, {k: int(v) for k, v in chosen.items()})
    print(f"wrote {out}: " + ", ".join(f"{k} at critical point {chosen[k]}" for k in KINDS))


if __name__ == "__main__":
    main()
