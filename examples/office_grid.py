"""One patrol lap on the fixed office map: events, true counts, bounds and expected counts.

The grid simulator (DESIGN §5) moves the robot once around its looped path
while ``--targets`` targets random-walk; every shadow and FOV event is
recorded (Secs. II-III).  The event stream is fed to

* the combinatorial filter (Sec. V), whose bounds are checked against the
  true counts after *every* tick; and
* the probabilistic filter (Sec. VI) with area-proportional binomial splits
  and noise-free FOV events, truncated by RT-LA at ``--max-entries``
  (exact whenever the table stays smaller).

``--mode`` sets the initial knowledge: ``exact`` counts, ``unknown`` counts
``[0, inf)`` or ``evader`` (one target, ``[0, 1]`` per shadow, hidden in
total iff it is not visible at ``t0``); expected counts need ``exact``.
``--plot`` saves the snapshot, shadow-sequence and bipartite figures of the
final ticks (see ``examples/figures.py``) into ``--out``.

Usage::

    python examples/office_grid.py [--targets 10] [--seed 1] [--mode exact] [--plot --out DIR]
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo import CombinatorialFilter, Split  # noqa: E402
from shadowinfo.grid import GridSimulator, load_office  # noqa: E402
from shadowinfo.probabilistic import ProbSplit, TruncatedFilter  # noqa: E402


def area_splits(events, labels):
    """Replace each ``Split`` of one tick by a ``ProbSplit`` with ``p = |a| / (|a| + |b|)``."""
    area = Counter(v for v in labels if v)
    children = {e.s: (e.a, e.b) for e in events if isinstance(e, Split)}

    def size(s: int) -> int:
        if s in area or s not in children:
            return area.get(s, 0)
        return sum(size(c) for c in children[s])

    out = []
    for e in events:
        if isinstance(e, Split):
            a, b = size(e.a), size(e.b)
            e = ProbSplit(e.s, e.a, e.b, a / (a + b) if a + b else 0.5)
        out.append(e)
    return out


def run(targets: int = 10, seed: int = 1, ticks: int = 0, mode: str = "exact", radius=None):
    """Simulate; return the simulator, its frames, per-tick event lists and probabilistic queues."""
    gridmap, path = load_office()
    sim = GridSimulator(gridmap, path, targets, seed, radius)
    frames, events, queue, grids = [sim.frame()], [], [], [list(sim.labels)]
    for _ in range(ticks or len(path)):
        n = len(sim.history)
        frames.append(sim.step())
        tick = sim.history[n:]
        events.append(tick)
        queue.append(area_splits(tick, sim.labels))
        grids.append(list(sim.labels))
    return dict(sim=sim, map=gridmap, path=path, frames=frames, events=events, queue=queue, grids=grids, mode=mode)


def expected_counts(r, upto: int, max_entries: int = 20000):
    """Expected counts after ``upto`` ticks (RT-LA truncation, exact if never triggered)."""
    init = r["sim"].initial_counts
    labels = tuple(sorted(init))
    q = [o for tick in r["queue"][:upto] for o in tick]
    f = TruncatedFilter(labels, {tuple(init[s] for s in labels): 1}, max_entries, mode="RT-LA").run(q)
    return f.expected_counts(), f


def fmt(b) -> str:
    return f"[{b[0]}, {'inf' if b[1] == math.inf else b[1]}]"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--targets", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ticks", type=int, default=0, help="ticks to simulate (default: one lap)")
    ap.add_argument("--mode", choices=("exact", "unknown", "evader"), default="exact")
    ap.add_argument("--radius", type=int, default=None, help="sensing radius in cells")
    ap.add_argument("--max-entries", type=int, default=20000)
    ap.add_argument("--plot", action="store_true", help="save figures of the final ticks")
    ap.add_argument("--out", default="figures", help="output directory for --plot")
    args = ap.parse_args()

    if args.mode == "evader" and args.targets != 1:
        ap.error("--mode evader needs --targets 1")
    r = run(args.targets, args.seed, args.ticks, args.mode, args.radius)
    sim = r["sim"]
    hidden = sum(sim.initial_counts.values())
    total = (hidden, hidden) if args.mode == "evader" else None
    filt = CombinatorialFilter(sim.initial_condition(args.mode), total=total)
    violations = 0
    for tick, frame in zip(r["events"], r["frames"][1:]):
        filt.extend(tick)
        bounds = filt.all_bounds()
        violations += sum(not bounds[s][0] <= n <= bounds[s][1] for s, n in frame["counts"])
    kinds = Counter(type(e).__name__.lower() for tick in r["events"] for e in tick)
    print(f"office map {r['map'].rows}x{r['map'].cols}, {len(r['events'])} ticks, "
          f"{args.targets} targets, seed {args.seed}, initial knowledge '{args.mode}'")
    print("events:", ", ".join(f"{k} {kinds[k]}" for k in ("appear", "disappear", "split", "merge", "enter", "exit")),
          f"(total {sum(kinds.values())})")
    print("truth outside the bounds at any tick:", violations)

    counts = dict(r["frames"][-1]["counts"])
    area = Counter(v for v in sim.labels if v)
    exp = expected_counts(r, len(r["events"]), args.max_entries)[0] if args.mode == "exact" else {}
    print(f"\n{'shadow':>7} {'cells':>6} {'true':>5} {'bounds':>10} {'E[count]':>9}")
    for s, b in filt.all_bounds().items():
        e = f"{exp[s]:.2f}" if s in exp else "-"
        print(f"{'s%d' % s:>7} {area[s]:>6} {counts[s]:>5} {fmt(b):>10} {e:>9}")
    if args.mode == "unknown":
        print("bounds on the initial number of hidden targets:", fmt(filt.initial_total_bounds()))
    if args.plot:
        import figures
        print("saved:", *figures.save_all(r, Path(args.out), tick=len(r["events"])), sep="\n  ")
    if violations:
        sys.exit(1)


if __name__ == "__main__":
    main()
