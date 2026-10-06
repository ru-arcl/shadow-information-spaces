"""Counting and passive pursuit-evasion on the office map (T-RO 2012, Sec. V-E).

*Counting.*  The number of targets is unknown: every shadow at ``t0`` starts
at ``[0, inf)``.  As the robot patrols, the filter refines the bound on the
number of targets that were hidden at ``t0``; adding the targets seen at
``t0`` bounds the total ``n``.  The lower bound verifies a known ``n`` once it
reaches it; the upper bound stays infinite unless every initial shadow is
eventually cleared ("the free space is not completely explored").

*Pursuit-evasion.*  One evader, ``[0, 1]`` per shadow at ``t0`` and a total of
exactly one.  Each shadow is *clear* ``[0, 0]``, *evader* ``[1, 1]`` or
*contaminated* ``[0, 1]``.  Two runs are shown: a worst-case evader that is
never seen (only component events, all counts zero), reported as the
contaminated fraction of the shadow area; and a random-walk evader, which is
caught the first time it leaves a shadow.

Usage::

    python examples/counting_and_pursuit.py [--targets 10] [--seed 1] [--laps 2] [--every 28]
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from office_grid import fmt, run  # noqa: E402
from shadowinfo import (Appear, Disappear, Enter, Exit, Merge, Split, counting_filter,  # noqa: E402
                        evader_status, pursuit_evasion_filter)


def unseen(tick):
    """The component events of a tick as seen by an evader that never shows itself."""
    out = []
    for e in tick:
        if isinstance(e, (Appear, Disappear)):
            out.append(type(e)(e.s, 0, 0))
        elif isinstance(e, (Split, Merge)):
            out.append(e)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--targets", type=int, default=10, help="targets for the counting task")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--evader-seed", type=int, default=3)
    ap.add_argument("--laps", type=int, default=2)
    ap.add_argument("--every", type=int, default=28, help="print a row every this many ticks")
    args = ap.parse_args()

    ticks = 196 * args.laps
    cnt = run(args.targets, args.seed, ticks)
    ev = run(1, args.evader_seed, ticks)
    seen0 = args.targets - sum(cnt["sim"].initial_counts.values())
    fc = counting_filter(cnt["sim"].initial_counts)
    worst = pursuit_evasion_filter(ev["sim"].initial_counts)
    rand = pursuit_evasion_filter(ev["sim"].initial_counts)
    caught = verified = None

    print(f"counting: {args.targets} targets, {seen0} visible at t0 | pursuit-evasion: one evader")
    print(f"\n{'tick':>5} {'n (counting)':>14}   {'unseen evader: contaminated':>28}   {'random evader'}")
    for t in range(1, ticks + 1):
        fc.extend(cnt["events"][t - 1])
        worst.extend(unseen(ev["events"][t - 1]))
        if caught is None:
            tick = ev["events"][t - 1]
            if any(isinstance(e, Exit) for e in tick):
                caught = t
            else:
                rand.extend(e for e in tick if not isinstance(e, Enter))
        lo, hi = fc.initial_total_bounds()
        if verified is None and seen0 + lo == args.targets:
            verified = t
        if t % args.every and t != ticks:
            continue
        area = Counter(v for v in ev["grids"][t] if v)
        st = evader_status(worst)
        bad = [s for s, v in st.items() if v != "clear"]
        frac = sum(area[s] for s in bad) / max(1, sum(area.values()))
        if caught is None:
            r = Counter(evader_status(rand).values())
            rs = f"{r['contaminated']} contaminated, {r['evader']} certain"
        else:
            rs = f"seen at tick {caught}"
        n = fmt((seen0 + lo, hi + seen0 if hi != math.inf else hi))
        print(f"{t:>5} {n:>14}   {len(bad):>3}/{len(st)} shadows, {frac:>4.0%} of area   {rs}")

    print(f"\ncounting: true n = {args.targets};",
          f"lower bound reached it at tick {verified}" if verified else "lower bound never reached it",
          "| upper bound:", "finite" if fc.initial_total_bounds()[1] != math.inf else "inf (not all explored)")


if __name__ == "__main__":
    main()
