"""Truncation heuristics vs exact filtering vs Monte Carlo (T-RO 2012, Sec. VII-B, Table IV).

The instance is "Fig. 16-like": the 14 component events of Fig. 16 with the
Table IV counts (10, 7, 8 targets in s1, s2, s3; 9 appear in s11), plus
``--fov`` noisy FOV observations scattered along the sequence.  The paper does
not publish its 32 FOV observations, so they are generated here from a
simulated ground truth (binomial(1/2) splits, true events enter w.p. 0.6,
sensor true-positive rate 0.9); the counts revealed by the three disappear
events are the simulated ones.  The numbers therefore differ from Table IV,
but the comparison is the same: exact (Algorithms 1, 2) against TR-X, RT-X and
RT-LA-X (mean (std) over ``--runs`` seeds) and 1000-trial rejection sampling.

The default of 24 observations peaks at about 3.8M exact entries (~1 GB,
~10 s); ``--fov 32`` matches the paper's count but needs ~2 GB for the exact
filter (``--no-exact`` skips it).

Usage::

    python examples/truncation_heuristics.py [--fov 24] [--max-entries 500 2000 10000 50000] [--runs 3]
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo import Appear, Disappear, Merge, Split  # noqa: E402
from shadowinfo.probabilistic import (ExactFilter, FovObservation, ProbDisappear, ProbSequence,  # noqa: E402
                                      TruncatedFilter, TruncationFailure, monte_carlo, symmetric_obs_model)
from shadowinfo.rng import Rng  # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "tor_fig16_sequence.json"
MODEL = symmetric_obs_model(0.9)


def fig16_like(n_fov: int, seed: int, p_enter: float = 0.6, p_true: float = 0.9):
    """Fig. 16 component events plus ``n_fov`` FOV observations of a simulated run."""
    base = ProbSequence.from_dict(json.loads(FIXTURE.read_text())["table_iv"]["prob_sequence"])
    rng = Rng(seed)
    x = dict(zip(base.labels, next(iter(base.joint))))
    slots = sorted(rng.randint(len(base.observations) + 1) for _ in range(n_fov))
    queue = []
    for t, o in enumerate(base.observations + [None]):
        while slots and slots[0] == t:
            slots.pop(0)
            s = sorted(x)[rng.randint(len(x))]
            e = "enter" if x[s] == 0 or rng.random() < p_enter else "exit"
            x[s] += 1 if e == "enter" else -1
            flip = {"enter": "exit", "exit": "enter"}
            queue.append(FovObservation(s, e if rng.random() < p_true else flip[e]))
        if isinstance(o, Split):
            n = x.pop(o.s)
            x[o.a] = sum(rng.random() < 0.5 for _ in range(n))
            x[o.b] = n - x[o.a]
        elif isinstance(o, Merge):
            x[o.s] = x.pop(o.a) + x.pop(o.b)
        elif isinstance(o, Appear):
            x[o.s] = o.lo
        elif isinstance(o, Disappear):
            o = ProbDisappear.of(o.s, {x.pop(o.s): 1})
        if o is not None:
            queue.append(o)
    return ProbSequence(base.labels, base.joint, queue), x


def row(name: str, means, secs: float, extra: str = "") -> None:
    print(f"  {name:<14}" + "".join(f"{m:>15}" for m in means) + f"{secs:>9.2f}  {extra}".rstrip())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fov", type=int, default=24, help="number of FOV observations")
    ap.add_argument("--seed", type=int, default=1, help="seed of the simulated ground truth")
    ap.add_argument("--max-entries", type=int, nargs="+", default=[500, 2000, 10000, 50000])
    ap.add_argument("--runs", type=int, default=3, help="seeds per randomized heuristic")
    ap.add_argument("--trials", type=int, default=1000, help="successful Monte Carlo trials")
    ap.add_argument("--no-exact", action="store_true", help="skip the exact filter (memory heavy)")
    args = ap.parse_args()

    seq, truth = fig16_like(args.fov, args.seed)
    final = sorted(truth)
    q = seq.observations
    print(f"Fig. 16-like instance: {len(q) - args.fov} component events + {args.fov} FOV observations")
    print("  revealed counts:", ", ".join(f"s{o.s}={o.lo}" for o in q if isinstance(o, Disappear)),
          "| true final counts:", ", ".join(f"s{s}={truth[s]}" for s in final))
    print(f"\n  {'heuristic':<14}" + "".join(f"{'E[s%d]' % s:>15}" for s in final) + f"{'t (s)':>9}")

    if not args.no_exact:
        t = time.perf_counter()
        f = ExactFilter(seq.labels, seq.joint, obs_model=MODEL).run(q)
        e = f.expected_counts()
        row("exact", [f"{e[s]:.2f}" for s in final], time.perf_counter() - t, f"peak {f.peak_entries} entries")

    for mode in ("TR", "RT", "RT-LA"):
        for x in args.max_entries:
            seeds = [1] if mode == "TR" else range(1, args.runs + 1)
            results, fails, t = [], 0, time.perf_counter()
            for seed in seeds:
                try:
                    g = TruncatedFilter(seq.labels, seq.joint, x, mode=mode, seed=seed, obs_model=MODEL).run(q)
                    results.append(g.expected_counts())
                except TruncationFailure:
                    fails += 1
            secs = (time.perf_counter() - t) / len(seeds)
            name = f"{mode}-{x}"
            if 3 * fails > len(seeds):
                print(f"  {name:<14}{'failure' if not results else 'frequent failure':>15}")
                continue
            cells = []
            for s in final:
                v = [r[s] for r in results]
                sd = f" ({statistics.pstdev(v):.2f})" if len(seeds) > 1 else ""
                cells.append(f"{statistics.fmean(v):.2f}{sd}")
            row(name, cells, secs, f"{fails}/{len(seeds)} failed" if fails else "")

    t = time.perf_counter()
    mc = monte_carlo(seq, trials=args.trials, seed=1, obs_model=MODEL, max_attempts=10 ** 7)
    e = mc.expected_counts()
    row("Monte Carlo", [f"{e[s]:.2f}" for s in final], time.perf_counter() - t,
        f"{mc.successes}/{mc.attempts} trials kept")


if __name__ == "__main__":
    main()
