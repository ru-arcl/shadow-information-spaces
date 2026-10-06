"""Generate ``tests/fixtures/js_filter_cases.json`` for the Python/JS filter parity tests.

The JS ports in ``docs/js/{maxflow,bipartite,filter,probabilistic}.js`` must
reproduce :mod:`shadowinfo` exactly: integer bounds bit for bit, probabilities
to ``1e-9``.  This script records the Python answers on

* random consistent shadow sequences (:func:`shadowinfo.testing.random_instance`)
  with FOV events, ``inf`` bounds and total constraints: ``all_bounds()``, set
  bounds, refined initial bounds, the initial total, the literal and corrected
  paper recipe (Eqs. 7, 8), witnesses, the bipartite I-state and bounds after
  every prefix of the events;
* perturbed (often infeasible) sequences (:func:`shadowinfo.testing.perturb`);
* counting and pursuit-evasion initial conditions (Sec. V-E);
* a large instance (:func:`shadowinfo.testing.growing_instance`) for timing;
* random probabilistic instances (Sec. VI): exact joint pmf, marginals and
  expected counts, TR / RT / RT-LA truncation and Monte Carlo runs.

Run from the repository root: ``python3 tools/gen_js_filter_fixtures.py``.
"""

from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from shadowinfo.bipartite import vertex_name  # noqa: E402
from shadowinfo.events import Appear, Disappear, Enter, Exit, Merge, ShadowSequence, Split  # noqa: E402
from shadowinfo.maxflow import InfeasibleError  # noqa: E402
from shadowinfo.nondeterministic import CombinatorialFilter  # noqa: E402
from shadowinfo.probabilistic import (  # noqa: E402
    BinomialSplit, ExactFilter, FovObservation, InconsistentObservationError, ProbAppear,
    ProbDisappear, ProbSequence, ProbSplit, TruncatedFilter, TruncationFailure, monte_carlo,
)
from shadowinfo.testing import growing_instance, perturb, random_instance  # noqa: E402

OUT = ROOT / "tests" / "fixtures" / "js_filter_cases.json"


def hi_json(hi: float) -> Optional[int]:
    return None if hi == math.inf else int(hi)


def bound_json(b) -> list:
    return [int(b[0]), hi_json(b[1])]


# -- nondeterministic -----------------------------------------------------------

def record_combinatorial(name: str, seq: ShadowSequence, total, rng: random.Random,
                         prefixes: bool = False, extras: bool = True) -> dict:
    case = {"name": name, "sequence": seq.to_dict(), "total": None if total is None else bound_json(total)}
    f = CombinatorialFilter.from_sequence(seq, total)
    case["feasible"] = f.feasible()
    if case["feasible"]:
        case["all_bounds"] = {str(s): bound_json(b) for s, b in f.all_bounds().items()}
        if extras:
            alive = f.alive()
            sets = []
            for _ in range(3 if len(alive) >= 2 else 0):
                q = rng.sample(alive, rng.randint(2, len(alive)))
                sets.append({"shadows": q, "bounds": bound_json(f.bounds(q))})
            case["set_bounds"] = sets
            case["refine_initial"] = {str(s): bound_json(b) for s, b in f.refine_initial_bounds().items()}
            case["initial_total"] = bound_json(f.initial_total_bounds())
            if alive:
                q = [rng.choice(alive)]
                case["paper"] = {"shadows": q, "literal": bound_json(f.bounds_paper(q)),
                                 "corrected": bound_json(f.bounds_paper(q, corrected=True))}
                ws = []
                for sense in ("max", "min"):
                    try:
                        w = f.witness(q, sense)
                    except ValueError:
                        continue
                    ws.append({"shadows": q, "sense": sense, "value": w.value,
                               "supply": [[vertex_name(v), n] for v, n in w.supply.items()],
                               "flow": [[vertex_name(u), r[0], vertex_name(r[1]), n]
                                        for (u, r), n in w.flow.items()]})
                case["witness"] = ws
        case["bipartite"] = f.istate.to_dict()
    if prefixes:
        g = CombinatorialFilter(seq.initial, total)
        out: List[Optional[dict]] = []
        for e in seq.events:
            g.apply(e)
            try:
                out.append({str(s): bound_json(b) for s, b in g.all_bounds().items()})
            except InfeasibleError:
                out.append(None)
        case["prefix_bounds"] = out
    return case


def combinatorial_cases() -> List[dict]:
    cases = []
    rng = random.Random(20120424)
    for i in range(240):
        kw = {}
        if i % 4 == 1:
            kw = dict(n_events=(20, 60), max_alive=12, p_fov=0.5)
        elif i % 4 == 2:
            kw = dict(p_inf=0.4, p_total=0.6)
        inst = random_instance(rng, **kw)
        cases.append(record_combinatorial(f"random-{i}", inst.seq, inst.total, rng, prefixes=i % 8 == 0))
    for i in range(80):
        inst = random_instance(rng, n_events=(3, 25))
        seq = perturb(rng, inst.seq)
        cases.append(record_combinatorial(f"perturbed-{i}", seq, inst.total, rng))
    for i in range(20):
        inst = random_instance(rng, p_fov=0.5)
        initial = {s: (0, math.inf) for s in inst.seq.initial}
        seq = ShadowSequence(initial, inst.seq.events)
        cases.append(record_combinatorial(f"counting-{i}", seq, None, rng))
    for i in range(20):
        inst = random_instance(rng, max_count=1, n_initial=(2, 5), p_fov=0.3)
        initial = {s: (0, 1) for s in inst.seq.initial}
        seq = ShadowSequence(initial, inst.seq.events)
        cases.append(record_combinatorial(f"pursuit-{i}", seq, (1, 1), rng))
    return cases


def performance_case() -> dict:
    inst = growing_instance(random.Random(7), n_component=500, target_alive=100)
    f = CombinatorialFilter.from_sequence(inst.seq)
    n_comp = sum(isinstance(e, (Appear, Disappear, Split, Merge)) for e in inst.seq.events)
    return {"name": "growing-500", "sequence": inst.seq.to_dict(), "total": None,
            "component_events": n_comp, "alive": len(f.alive()),
            "all_bounds": {str(s): bound_json(b) for s, b in f.all_bounds().items()}}


# -- probabilistic ------------------------------------------------------------------

def random_prob_instance(rng: random.Random, big: bool):
    n_shadows = rng.randint(1, 3)
    labels = list(range(1, n_shadows + 1))
    joint: Dict[tuple, float] = {}
    for _ in range(rng.randint(1, 4)):
        key = tuple(rng.randint(0, 4 if big else 2) for _ in labels)
        joint[key] = joint.get(key, 0.0) + rng.randint(1, 5) / 7
    alive = list(labels)
    nxt = n_shadows + 1
    obs = []
    for _ in range(rng.randint(8, 22) if big else rng.randint(2, 9)):
        kind = rng.choice(["fov", "fov", "fov", "split", "split", "merge", "appear", "disappear",
                           "enter", "exit"])
        s = rng.choice(alive) if alive else None
        if kind == "appear" or not alive:
            if rng.random() < 0.6:
                dist = {n: rng.randint(1, 4) / 5 for n in rng.sample(range(4), rng.randint(1, 3))}
                obs.append(ProbAppear.of(nxt, dist))
            else:
                lo = rng.randint(0, 2)
                obs.append(Appear(nxt, lo, lo + rng.randint(0, 2)))
            alive.append(nxt)
            nxt += 1
        elif kind == "fov":
            obs.append(FovObservation(s, rng.choice(["enter", "exit", "null"])))
        elif kind == "split" and len(alive) < (6 if big else 4):
            if rng.random() < 0.5:
                obs.append(ProbSplit(s, nxt, nxt + 1, rng.choice([0.2, 0.3, 0.5, 0.75, 0.9])))
            else:
                obs.append(Split(s, nxt, nxt + 1))
            alive.remove(s)
            alive += [nxt, nxt + 1]
            nxt += 2
        elif kind == "merge" and len(alive) >= 2:
            a, b = rng.sample(alive, 2)
            obs.append(Merge(a, b, nxt))
            alive = [x for x in alive if x not in (a, b)] + [nxt]
            nxt += 1
        elif kind == "disappear" and len(alive) >= 2:
            if rng.random() < 0.6:
                support = rng.sample(range(5 if big else 4), rng.randint(1, 3))
                obs.append(ProbDisappear.of(s, {n: rng.randint(1, 3) / 4 for n in support}))
            else:
                lo = rng.randint(0, 2)
                obs.append(Disappear(s, lo, lo + rng.randint(0, 3)))
            alive.remove(s)
        elif kind == "enter":
            obs.append(Enter(s, rng.randint(1, 2)))
        elif kind == "exit":
            obs.append(Exit(s, 1))
    return labels, joint, obs


def random_obs_model(rng: random.Random) -> dict:
    model = {}
    for y in ("enter", "exit", "null"):
        w = [rng.randint(0, 4) for _ in range(3)]
        if not any(w):
            w[rng.randrange(3)] = 1
        w[("enter", "exit", "null").index(y)] += 3
        model[y] = {e: x / sum(w) for e, x in zip(("enter", "exit", "null"), w)}
    return model


MAX_STORED_ENTRIES = 300


def pmf_json(f) -> dict:
    """Labels, entry count and marginals; the full joint only when it is small."""
    out = {"labels": list(f.labels), "entries": f.entries,
           "marginals": {str(s): {str(n): float(p) for n, p in f.marginal(s).items()} for s in f.labels}}
    if f.entries <= MAX_STORED_ENTRIES:
        out["joint"] = [[list(k), float(p)] for k, p in sorted(f.table.items())]
    return out


def run_exact(labels, joint, obs, split_p, model) -> dict:
    f = ExactFilter(labels, joint, split_rule=BinomialSplit(split_p), obs_model=model)
    try:
        f.run(obs)
    except InconsistentObservationError as err:
        return {"error": type(err).__name__}
    return {**pmf_json(f), "peak_entries": f.peak_entries,
            "expected": {str(s): float(v) for s, v in f.expected_counts().items()}}


def run_truncated(labels, joint, obs, split_p, model, mode, max_entries, seed, counts_fov=True) -> dict:
    f = TruncatedFilter(labels, joint, max_entries, mode=mode, seed=seed, split_rule=BinomialSplit(split_p),
                        obs_model=model, lookahead_counts_fov=counts_fov)
    spec = {"mode": mode, "max_entries": max_entries, "seed": seed, "lookahead_counts_fov": counts_fov}
    try:
        f.run(obs)
    except InconsistentObservationError as err:
        return {**spec, "error": type(err).__name__, "truncations": f.truncations}
    return {**spec, **pmf_json(f), "truncations": f.truncations, "truncated_mass": float(f.truncated_mass),
            "peak_entries": f.peak_entries,
            "expected": {str(s): float(v) for s, v in f.expected_counts().items()}}


def run_mc(labels, joint, obs, split_p, model, trials, seed) -> dict:
    spec = {"trials": trials, "seed": seed}
    try:
        r = monte_carlo(ProbSequence(tuple(labels), dict(joint), list(obs)), trials=trials, seed=seed,
                        split_rule=BinomialSplit(split_p), obs_model=model, max_attempts=50 * trials)
    except InconsistentObservationError as err:
        return {**spec, "error": type(err).__name__}
    return {**spec, "labels": list(r.labels), "joint": [[list(k), p] for k, p in sorted(r.table.items())],
            "successes": r.successes, "attempts": r.attempts}


def prob_case(name, labels, joint, obs, split_p, model, rng, truncate=True, mc=True) -> dict:
    seq = ProbSequence(tuple(labels), dict(joint), list(obs)).to_dict()
    case = {"name": name, "sequence": seq, "split_rule": {"type": "binomial", "p": split_p},
            "obs_model": model, "exact": run_exact(labels, joint, obs, split_p, model)}
    if truncate:
        peak = case["exact"].get("peak_entries", 50)
        budgets = sorted({max(1, peak // 4), max(1, peak // 2), max(2, peak - 1)})
        runs = []
        for mode in ("TR", "RT", "RT-LA"):
            for m in budgets:
                runs.append(run_truncated(labels, joint, obs, split_p, model, mode, m, rng.randint(1, 10 ** 6)))
        runs.append(run_truncated(labels, joint, obs, split_p, model, "RT-LA", budgets[0], 5, counts_fov=False))
        case["truncated"] = runs
    if mc:
        case["monte_carlo"] = run_mc(labels, joint, obs, split_p, model, 200, rng.randint(1, 10 ** 6))
    return case


def probabilistic_cases() -> List[dict]:
    rng = random.Random(2010)
    cases = []
    for i in range(40):
        big = i % 2 == 1
        labels, joint, obs = random_prob_instance(rng, big)
        model = random_obs_model(rng)
        split_p = rng.choice([0.5, 0.4, 0.7])
        cases.append(prob_case(f"prob-random-{i}", labels, joint, obs, split_p, model, rng))
    fx = json.loads((ROOT / "tests" / "fixtures" / "tor_fig12_table3.json").read_text())
    seq = ProbSequence.from_dict(fx["prob_sequence"])
    cases.append(prob_case("tor-table3", list(seq.labels), seq.joint, seq.observations, 0.5,
                           fx["obs_model"], rng))
    fx = json.loads((ROOT / "tests" / "fixtures" / "tor_fig16_sequence.json").read_text())
    seq = ProbSequence.from_dict(fx["table_iv"]["prob_sequence"])
    case = prob_case("tor-fig16-table4", list(seq.labels), seq.joint, seq.observations, 0.5, None, rng,
                     truncate=False, mc=False)
    case["truncated"] = [
        run_truncated(list(seq.labels), seq.joint, seq.observations, 0.5, None, mode, m, 11)
        for mode, m in (("TR", 2000), ("RT", 2000), ("RT-LA", 2000), ("TR", 20000), ("RT", 20000),
                        ("RT-LA", 20000))
    ]
    cases.append(case)
    return cases


def main() -> None:
    data = {
        "generator": "tools/gen_js_filter_fixtures.py",
        "combinatorial": combinatorial_cases(),
        "performance": performance_case(),
        "probabilistic": probabilistic_cases(),
    }
    OUT.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    comb = data["combinatorial"]
    prob = data["probabilistic"]
    print(f"wrote {OUT.relative_to(ROOT)}: {len(comb)} combinatorial "
          f"({sum(not c['feasible'] for c in comb)} infeasible), {len(prob)} probabilistic "
          f"({sum('error' in c['exact'] for c in prob)} inconsistent), "
          f"performance {data['performance']['component_events']} component events / "
          f"{data['performance']['alive']} alive")


if __name__ == "__main__":
    main()
