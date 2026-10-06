"""T-RO 2012 Fig. 12 / Table III: the exact probabilistic filter, step by step (Sec. VI-D).

Two targets start in each of ``s1`` and ``s2``.  The observation queue is
``y_x`` on s1, split s2 -> s3, s4, ``y_e`` on s3, merge s1, s3 -> s5, and
disappear s5 revealing 1 or 2 targets with probability 1/2 each.  Splits are
binomial(1/2) and FOV observations have true-positive rate 0.9 (no null
events).  The joint pmf is propagated with exact rational arithmetic
(Algorithms 1 and 2), giving ``P(s4 = 0, 1, 2) = 1/13, 2/13, 10/13``; the paper
prints 0.0769 / 0.1538 / 0.7692.

The rejection-sampling Monte Carlo baseline of Sec. VI-E is then run for
several seeds (the paper reports 0.079 / 0.154 / 0.767 with 1000 trials), and
the truncation remark of Sec. VI-E is reproduced: dropping the two smallest
entries after the merge loses ``P(s4 = 0)`` altogether.

Usage::

    python examples/table3_probabilistic.py [--trials 1000] [--seeds 5]
"""

from __future__ import annotations

import argparse
import sys
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo import Merge, Split  # noqa: E402
from shadowinfo.probabilistic import (BinomialSplit, ExactFilter, FovObservation, ProbDisappear,  # noqa: E402
                                      ProbSequence, monte_carlo, symmetric_obs_model)

PAPER_MC = {0: 0.079, 1: 0.154, 2: 0.767}
SHORT = {"enter": "e", "exit": "x", "null": "n"}


def describe(o) -> str:
    if isinstance(o, FovObservation):
        return f"y_{SHORT[o.y]}, s{o.s}"
    if isinstance(o, Split):
        return f"split s{o.s} -> s{o.a}, s{o.b}"
    if isinstance(o, Merge):
        return f"merge s{o.a}, s{o.b} -> s{o.s}"
    return f"disappear s{o.s}"


def show(filt: ExactFilter, title: str) -> None:
    print(f"{title}  ({filt.entries} entries)")
    for key, p in sorted(filt.table.items()):
        cells = ", ".join(f"s{s}={x}" for s, x in zip(filt.labels, key))
        print(f"    P({cells}) = {str(p):>9} = {float(p):.4f}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--trials", type=int, default=1000, help="successful Monte Carlo trials per seed")
    ap.add_argument("--seeds", type=int, default=5, help="number of Monte Carlo seeds")
    args = ap.parse_args()

    queue = [FovObservation(1, "exit"), Split(2, 3, 4), FovObservation(3, "enter"), Merge(1, 3, 5),
             ProbDisappear.of(5, {1: F(1, 2), 2: F(1, 2)})]
    model = symmetric_obs_model(F(9, 10))
    kw = dict(split_rule=BinomialSplit(F(1, 2)), obs_model=model)
    filt = ExactFilter.from_counts({1: 2, 2: 2}, **kw)

    print("Table III: exact joint pmf after each observation\n")
    show(filt, "initial")
    snapshots = []
    for t, o in enumerate(queue):
        snapshots.append(filt.copy())
        filt.apply(o, queue[t + 1:])
        show(filt, describe(o))
    exact = filt.marginal(4)

    print("\nMonte Carlo (Sec. VI-E), P(s4 = 0, 1, 2):")
    seq = ProbSequence.from_counts({1: 2, 2: 2}, queue)
    print(f"    {'exact':<16}" + "".join(f"{float(exact.get(n, 0)):>8.4f}" for n in range(3)))
    print(f"    {'paper MC (1000)':<16}" + "".join(f"{PAPER_MC[n]:>8.3f}" for n in range(3)))
    for seed in range(1, args.seeds + 1):
        mc = monte_carlo(seq, trials=args.trials, seed=seed, split_rule=BinomialSplit(0.5),
                         obs_model=symmetric_obs_model(0.9))
        m = mc.marginal(4)
        print(f"    {f'seed {seed}':<16}" + "".join(f"{m.get(n, 0):>8.3f}" for n in range(3))
              + f"   ({mc.successes}/{mc.attempts} trials kept)")

    trunc = snapshots[-1]
    trunc.table = dict(sorted(trunc.table.items(), key=lambda kp: -kp[1])[:6])
    trunc.apply(queue[-1])
    m = trunc.marginal(4)
    print("\nTruncation remark (Sec. VI-E): keep the 6 largest entries after the merge")
    print("    P(s4 = 0, 1, 2) =", ", ".join(str(m.get(n, 0)) for n in range(3)))

    ok = exact == {0: F(1, 13), 1: F(2, 13), 2: F(10, 13)}
    print("\nexact result equals 1/13, 2/13, 10/13:", "yes" if ok else "NO")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
