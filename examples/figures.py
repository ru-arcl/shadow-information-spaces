"""README figures for the office map: snapshot, shadow sequence and bipartite I-state (plus ``events.png``).

All three describe one window of ticks ``(t0, t1]`` of an office-map run:

* ``office_snapshot.png`` -- the map at ``t1``: visible region, shadows
  coloured and labelled with the combinatorial bounds, the true count and the
  expected count, the robot, its path and the targets (ground truth);
* ``shadow_sequence.png`` -- the shadow sequence of the window in the style
  of T-RO Figs. 5 and 7 (appear, split, merge, disappear and FOV arrows);
* ``bipartite.png`` -- the bipartite I-state of the window (Fig. 11(c)).

It also writes ``events.png``, the appear / disappear / split / merge illustrations on map 12
(:mod:`event_illustrations`).

The snapshot shows the bounds of the filter over the whole run.  The
bipartite I-state is that of a window filter which starts from the
full-history bounds at ``t0`` plus the bound on their sum; its bounds are
sound and, for the default window, equal to the full filter's.

Usage::

    python examples/figures.py [--out docs/assets] [--targets 12] [--seed 1] [--t0 148] [--t1 165]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]

import event_illustrations  # noqa: E402
from office_grid import expected_counts, run  # noqa: E402
from shadowinfo import (Appear, CombinatorialFilter, Disappear, Enter, Exit, Merge, Pseudo,  # noqa: E402
                        Split, event_labels)

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
WALL = "#46453f"
DPI = 150

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
})


def fmt(b) -> str:
    return f"[{b[0]}, {'∞' if b[1] == math.inf else b[1]}]"


def mix(c, w: float):
    return tuple(w * x + (1 - w) * y for x, y in zip(to_rgb(c), to_rgb(SURFACE)))


def window_filter(r, t0: int, t1: int):
    """Full-run filter after ``t1`` ticks and a filter over ticks ``(t0, t1]`` seeded with its bounds at ``t0``."""
    full = CombinatorialFilter(r["sim"].initial_condition(r["mode"]))
    for tick in r["events"][:t0]:
        full.extend(tick)
    b0 = full.all_bounds()
    win = CombinatorialFilter(b0, total=full.bounds(list(b0)))
    for tick in r["events"][t0:t1]:
        full.extend(tick)
        win.extend(tick)
    return full, win


def shadow_colors(r, t1: int):
    alive = sorted(set(r["grids"][t1]) - {0})
    return {s: SERIES[i % len(SERIES)] for i, s in enumerate(alive)}


# -- map snapshot -----------------------------------------------------------------------------
def plot_snapshot(r, t1: int, filt, colors, path: Path) -> None:
    gm, grid, frame = r["map"], r["grids"][t1], r["frames"][t1]
    rows, cols = gm.rows, gm.cols
    img = [[to_rgb(WALL) if gm.blocked[i * cols + j] else
            (mix(colors[grid[i * cols + j]], 0.42) if grid[i * cols + j] else to_rgb(SURFACE))
            for j in range(cols)] for i in range(rows)]
    fig, ax = plt.subplots(figsize=(9.2, 6.6))
    ax.imshow(img, interpolation="nearest", extent=(-0.5, cols - 0.5, rows - 0.5, -0.5))
    loop = r["path"] + r["path"][:1]
    ax.plot([c for _, c in loop], [rr for rr, _ in loop], color=MUTED, lw=1, alpha=0.8, zorder=2)
    counts, bounds = dict(frame["counts"]), filt.all_bounds()
    exp = expected_counts(r, t1)[0] if r["mode"] == "exact" else {}
    for tr, tc in frame["targets"]:
        hidden = grid[tr * cols + tc] != 0
        ax.scatter(tc, tr, s=30 if hidden else 48, zorder=4, linewidths=1.4 if hidden else 0.8,
                   edgecolors=INK if hidden else SURFACE,
                   facecolors=SURFACE if hidden else INK)
    ax.scatter(frame["robot"][1], frame["robot"][0], s=150, marker="o", color=INK, edgecolors=SURFACE,
               linewidths=2, zorder=5)
    for s, color in colors.items():
        cells = [divmod(i, cols) for i, v in enumerate(grid) if v == s]
        cr, cc = sum(c[0] for c in cells) / len(cells), sum(c[1] for c in cells) / len(cells)
        pr, pc = min(cells, key=lambda c: (c[0] - cr) ** 2 + (c[1] - cc) ** 2)
        text = f"s{s}  {fmt(bounds[s])}\ntrue {counts[s]}" + (f" · E {exp[s]:.1f}" if s in exp else "")
        ax.annotate(text, (pc, pr), ha="center", va="center", fontsize=8.5, color=INK, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.35", fc=SURFACE, ec=color, lw=1.6, alpha=0.95))
    ax.set_axis_off()
    handles = [Line2D([], [], marker="o", ls="", ms=10, mfc=INK, mec=SURFACE, mew=2, label="robot"),
               Line2D([], [], color=MUTED, lw=1, label="patrol path"),
               Line2D([], [], marker="s", ls="", ms=9, mfc=SURFACE, mec=AXIS, label="visible region"),
               Line2D([], [], marker="o", ls="", ms=6, mfc=INK, mec=SURFACE, label="visible target"),
               Line2D([], [], marker="o", ls="", ms=5, mfc=SURFACE, mec=INK, mew=1.4,
                      label="hidden target (ground truth)"),
               Line2D([], [], marker="s", ls="", ms=9, mfc=mix(SERIES[0], 0.42), mec="none",
                      label="shadow s#: [lower, upper], true count, E = expected")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False,
              fontsize=8.5, labelcolor=INK2, handletextpad=0.4, columnspacing=1.4)
    ax.set_title(f"Office map at tick {t1}: shadows with the combinatorial filter's bounds", loc="left",
                 fontsize=11, color=INK, pad=8)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


# -- shadow sequence (Figs. 5, 7) -------------------------------------------------------------
def layout_sequence(r, t0: int, t1: int):
    """x of every event (one unit each, a gap between ticks with events) and the lane of every label."""
    alive0 = sorted(set(r["grids"][t0]) - {0})
    xs, ticks, x = [], [], 0.0
    for t in range(t0 + 1, t1 + 1):
        tick = r["events"][t - 1]
        if not tick:
            continue
        ticks.append((t, x + 0.5, x + len(tick) + 0.5))
        for e in tick:
            x += 1
            xs.append((x, t, e))
        x += 0.7
    end = x + 0.5
    birth = {s: 0.0 for s in alive0}
    death = {}
    for x, _, e in xs:
        consumed, created = event_labels(e)
        for s in consumed:
            death[s] = x
        for s in created:
            birth[s] = x
    lane = {}

    def free(y, s):
        a, b = birth[s], death.get(s, end)
        return all(not (lane[o] == y and birth[o] < b and a < death.get(o, end)) for o in lane)

    def nearest(y0, s):
        for d in range(0, 64):
            for y in (y0 + d, y0 - d):
                if free(y, s):
                    return y
        raise RuntimeError("no free lane")

    for i, s in enumerate(alive0):
        lane[s] = 2 * i
    for x, _, e in xs:
        if isinstance(e, Split):
            lane[e.a] = nearest(lane[e.s], e.a)
            lane[e.b] = nearest(lane[e.s] + 1, e.b)
        elif isinstance(e, Merge):
            lane[e.s] = nearest(min(lane[e.a], lane[e.b]), e.s)
        elif isinstance(e, Appear):
            lane[e.s] = nearest(max(lane.values()) + 1, e.s)
    rank = {y: k for k, y in enumerate(sorted(set(lane.values())))}
    return xs, ticks, end, birth, death, {s: rank[y] for s, y in lane.items()}


def plot_sequence(r, t0: int, t1: int, colors, path: Path) -> None:
    xs, ticks, end, birth, death, lane = layout_sequence(r, t0, t1)
    n = max(lane.values()) + 1
    fig, ax = plt.subplots(figsize=(min(13, 2.2 + 0.36 * end), 1.5 + 0.42 * n))
    for k, (t, a, b) in enumerate(ticks):
        if k % 2:
            ax.axvspan(a, b, color=GRID, alpha=0.45, lw=0, zorder=0)
    for s, y in lane.items():
        x0, x1 = birth[s], death.get(s, end)
        final = s not in death
        ax.plot([x0, x1], [y, y], color=colors[s] if final else INK2, lw=2.4 if final else 1.6,
                solid_capstyle="butt", zorder=1)
    for x, _, e in xs:
        if isinstance(e, (Split, Merge)):
            ys = [lane[v] for v in (e.s, e.a, e.b)]
            ax.plot([x, x], [min(ys), max(ys)], color=INK2, lw=1.6, zorder=1)
        elif isinstance(e, Disappear):
            ax.plot([x, x], [lane[e.s] - 0.28, lane[e.s] + 0.28], color=INK, lw=2.4, zorder=2)
            ax.text(x + 0.15, lane[e.s] - 0.3, f"{e.lo}", fontsize=7.5, color=INK2, va="bottom")
        elif isinstance(e, (Enter, Exit)):
            y, up = lane[e.s], isinstance(e, Exit)
            ax.annotate("", (x, y - (0.62 if up else 0.12)), (x, y - (0.12 if up else 0.62)),
                        arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.1, mutation_scale=8))
            tag = ("$e_x$" if up else "$e_e$") + (f"×{e.k}" if e.k > 1 else "")
            ax.text(x + 0.1, y - 0.5, tag, fontsize=8, color=INK2, va="center")
    for s, y in lane.items():
        final = s not in death
        spots = []
        if birth[s] == 0:
            spots.append((0.0, SURFACE))
        elif not final or end - birth[s] > 2:
            spots.append((birth[s] + (0 if s in _appeared(xs) else 0.5), SURFACE))
        if final:
            spots.append((end, colors[s]))
        for px, fc in spots:
            ax.scatter(px, y, s=330, color=fc, edgecolors=INK2 if fc == SURFACE else SURFACE, lw=1.2, zorder=3)
            ax.text(px, y, str(s), ha="center", va="center", fontsize=6.8, zorder=4,
                    color=SURFACE if fc != SURFACE else INK)
    ax.set_xticks([(a + b) / 2 for _, a, b in ticks], [str(t) for t, _, _ in ticks])
    ax.tick_params(axis="x", length=0, labelsize=8)
    ax.set_xlabel("tick (ticks without events omitted; events of one tick in order)", fontsize=8.5)
    ax.set_yticks([])
    ax.set_ylim(n - 0.4, -0.95)
    ax.set_xlim(-0.7, end + 0.7)
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.set_title(f"Shadow sequence, ticks {t0 + 1}–{t1}", loc="left", fontsize=11, pad=24)
    ax.text(0, 1.015, "vertical joins: split / merge · bar: disappear (revealed count) · arrows: FOV enter $e_e$ / "
            "exit $e_x$ · filled circles: shadows alive at the end", transform=ax.transAxes, fontsize=8.5,
            color=INK2, va="bottom")
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def _appeared(xs):
    return {e.s for _, _, e in xs if isinstance(e, Appear)}


# -- bipartite I-state (Fig. 11(c)) -----------------------------------------------------------
def vname(v) -> str:
    if isinstance(v, Pseudo):
        return f"{'+' if v.kind == 'enter' else '−'} s{v.s}"
    return f"s{v}"


def plot_bipartite(r, t0: int, t1: int, filt, colors, path: Path) -> None:
    st = filt.istate
    st.flush()
    left, right, edges = list(st.left), list(st.right()), st.edges()
    li = {v: k for k, v in enumerate(left)}
    nbr = {rv: [li[u] for u, w in edges if w == rv] for rv in right}
    right.sort(key=lambda rv: (rv[0] == "alive", sum(nbr[rv]) / max(1, len(nbr[rv]))))
    ri = {rv: k for k, rv in enumerate(right)}
    bounds, counts = filt.all_bounds(), dict(r["frames"][t1]["counts"])
    h = max(len(left), len(right))
    ly = lambda k: -k - (h - len(left)) / 2  # noqa: E731
    ry = lambda k: -k - (h - len(right)) / 2  # noqa: E731
    fig, ax = plt.subplots(figsize=(8.6, 1.2 + 0.4 * h))
    for u, rv in edges:
        alive = rv[0] == "alive"
        ax.plot([0, 1], [ly(li[u]), ry(ri[rv])], color=colors[rv[1]] if alive else AXIS, lw=1.1,
                alpha=0.8 if alive else 1, zorder=1)
    for v in left:
        y = ly(li[v])
        kind = "FOV enter" if isinstance(v, Pseudo) else ("at t0" if v in st.initial else "appeared")
        ax.scatter(0, y, s=60, color=INK2 if kind == "at t0" else SURFACE, edgecolors=INK2, lw=1.4,
                   marker="D" if kind == "FOV enter" else "o", zorder=2)
        ax.text(-0.05, y, f"{vname(v)}  {fmt(st.left[v])}", ha="right", va="center", fontsize=8.5)
        ax.text(-0.75, y, kind, ha="left", va="center", fontsize=7.5, color=MUTED)
    for rv in right:
        (tag, v), y = rv, ry(ri[rv])
        alive, pseudo = tag == "alive", isinstance(v, Pseudo)
        ax.scatter(1, y, s=90 if alive else 60, color=colors[v] if alive else SURFACE,
                   edgecolors=SURFACE if alive else INK2, lw=1.4, marker="D" if pseudo else "o", zorder=2)
        n = st.disappeared[v][0] if not alive else None
        text = (f"{vname(v)}  {fmt(bounds[v])}   true {counts[v]}" if alive
                else f"{vname(v)}  {n} {'exited' if pseudo else 'revealed'}")
        ax.text(1.05, y, text, ha="left", va="center", fontsize=8.5)
        ax.text(1.95, y, "alive at end" if alive else ("FOV exit" if pseudo else "disappeared"),
                ha="right", va="center", fontsize=7.5, color=MUTED)
    top = 1.1
    ax.text(-0.05, top, "sources (left vertices)", ha="right", fontsize=8.5, color=INK2)
    ax.text(1.05, top, "sinks (right vertices)", ha="left", fontsize=8.5, color=INK2)
    ax.set_xlim(-0.78, 1.97)
    ax.set_ylim(-h + 0.4, top + 0.5)
    ax.set_axis_off()
    ax.set_title(f"Bipartite I-state, ticks {t0 + 1}–{t1}: a target may move along any edge; bounds by max-flow",
                 loc="left", fontsize=10.5, pad=6)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def save_all(r, out: Path, tick: int, window: int = 17, t0=None):
    """Write the three figures for the window ending at ``tick``; return their paths."""
    out.mkdir(parents=True, exist_ok=True)
    t0 = max(0, tick - window) if t0 is None else t0
    full, win = window_filter(r, t0, tick)
    colors = shadow_colors(r, tick)
    paths = [out / "office_snapshot.png", out / "shadow_sequence.png", out / "bipartite.png"]
    plot_snapshot(r, tick, full, colors, paths[0])
    plot_sequence(r, t0, tick, colors, paths[1])
    plot_bipartite(r, t0, tick, win, colors, paths[2])
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(HERE.parent / "docs" / "assets"))
    ap.add_argument("--targets", type=int, default=12)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--t0", type=int, default=148)
    ap.add_argument("--t1", type=int, default=165)
    args = ap.parse_args()
    r = run(args.targets, args.seed, max(196, args.t1))
    paths = save_all(r, Path(args.out), args.t1, t0=args.t0)
    paths.append(event_illustrations.save(Path(args.out) / "events.png")[0])
    print("saved:", *paths, sep="\n  ")


if __name__ == "__main__":
    main()
