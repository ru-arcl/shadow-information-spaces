"""Random shadow sequences with ground truth, and witness replay (for tests and fixtures).

:func:`random_instance` simulates targets first and then derives observations
that contain the truth, so every generated instance is consistent:

* each shadow's true count is tracked through split/merge (random split of the
  targets), enter/exit (``k`` targets cross the FOV boundary);
* initial / appear / disappear observations are ranges around the true count
  (sometimes exact, sometimes with ``hi = inf``);
* an optional ``total`` range contains the true number of targets at ``t0``.

:func:`replay_witness` pushes a bipartite witness (from a filter with
``fov="naive"``) through the original events target by target and checks every
observation, which certifies that a reported bound is attained.
"""

from __future__ import annotations

import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Tuple

from .bipartite import Pseudo
from .events import (Appear, Bound, Disappear, Enter, Event, Exit, Merge, ShadowSequence, Split)


@dataclass
class Instance:
    """A consistent shadow sequence plus the target counts that generated it."""

    seq: ShadowSequence
    total: Optional[Bound]
    truth_left: Dict[int, int]
    truth_final: Dict[int, int]
    truth_gone: Dict[int, int] = field(default_factory=dict)


def _noisy(rng: random.Random, c: int, p_exact: float, p_inf: float, spread: int = 2) -> Bound:
    if rng.random() < p_exact:
        return c, c
    lo = max(0, c - rng.randint(0, spread))
    hi = math.inf if rng.random() < p_inf else c + rng.randint(0, spread)
    return lo, hi


def random_instance(rng: random.Random, n_initial: Tuple[int, int] = (1, 5),
                    n_events: Tuple[int, int] = (0, 25), p_fov: float = 0.3, p_inf: float = 0.1,
                    p_total: float = 0.3, max_alive: int = 8, max_count: int = 4,
                    p_exact: float = 0.4) -> Instance:
    """Random consistent instance (see module docstring)."""
    label = rng.randint(1, 3)
    counts: Dict[int, int] = {}
    initial: Dict[int, Bound] = {}
    truth_left: Dict[int, int] = {}
    for _ in range(rng.randint(*n_initial)):
        c = rng.randint(0, max_count)
        counts[label] = truth_left[label] = c
        initial[label] = _noisy(rng, c, p_exact, p_inf)
        label += rng.randint(1, 2)

    def fresh() -> int:
        nonlocal label
        label += rng.randint(1, 2)
        return label - 1

    events: List[Event] = []
    gone: Dict[int, int] = {}
    for _ in range(rng.randint(*n_events)):
        alive = list(counts)
        options = [("appear", 1.0)]
        if alive:
            options += [("disappear", 1.0 if len(alive) > 1 else 0.3)]
            if len(alive) < max_alive:
                options += [("split", 1.5)]
            options += [("enter", 2 * p_fov)]
            if any(counts.values()):
                options += [("exit", 2 * p_fov)]
        if len(alive) >= 2:
            options += [("merge", 1.5)]
        kind = rng.choices([o for o, _ in options], [w for _, w in options])[0]
        if kind == "appear":
            s, c = fresh(), rng.randint(0, max_count - 1)
            counts[s] = truth_left[s] = c
            lo, hi = _noisy(rng, c, p_exact, p_inf)
            events.append(Appear(s, lo, hi))
        elif kind == "disappear":
            s = rng.choice(alive)
            c = gone[s] = counts.pop(s)
            lo, hi = _noisy(rng, c, 0.5, 0.0)
            events.append(Disappear(s, lo, hi))
        elif kind == "split":
            s = rng.choice(alive)
            c = counts.pop(s)
            a, b = fresh(), fresh()
            counts[a] = rng.randint(0, c)
            counts[b] = c - counts[a]
            events.append(Split(s, a, b))
        elif kind == "merge":
            a, b = rng.sample(alive, 2)
            s = fresh()
            counts[s] = counts.pop(a) + counts.pop(b)
            events.append(Merge(a, b, s))
        elif kind == "enter":
            s, k = rng.choice(alive), rng.randint(1, 2)
            counts[s] += k
            events.append(Enter(s, k))
        else:
            s = rng.choice([t for t in alive if counts[t] > 0])
            k = rng.randint(1, counts[s])
            counts[s] -= k
            events.append(Exit(s, k))
    total = None
    if rng.random() < p_total:
        n = sum(truth_left[s] for s in initial)
        total = (max(0, n - rng.randint(0, 2)), math.inf if rng.random() < p_inf else n + rng.randint(0, 2))
    return Instance(ShadowSequence(initial, events), total, truth_left, dict(counts), gone)


def growing_instance(rng: random.Random, n_component: int = 500, target_alive: int = 100,
                     n_initial: int = 20, p_fov: float = 0.5) -> Instance:
    """Large instance: about ``n_component`` component events with ~``target_alive`` alive shadows."""
    label = 1
    counts: Dict[int, int] = {}
    initial: Dict[int, Bound] = {}
    truth_left: Dict[int, int] = {}
    for _ in range(n_initial):
        c = rng.randint(0, 5)
        counts[label] = truth_left[label] = c
        initial[label] = _noisy(rng, c, 0.3, 0.05, 3)
        label += 1
    events: List[Event] = []
    n = 0
    while n < n_component:
        alive = list(counts)
        grow = len(alive) < target_alive
        r = rng.random()
        if rng.random() < p_fov and alive:
            s = rng.choice(alive)
            if counts[s] > 0 and rng.random() < 0.5:
                k = rng.randint(1, counts[s])
                counts[s] -= k
                events.append(Exit(s, k))
            else:
                counts[s] += 1
                events.append(Enter(s, 1))
            continue
        n += 1
        if r < (0.35 if grow else 0.15) or len(alive) < 2:
            c = rng.randint(0, 3)
            counts[label] = truth_left[label] = c
            events.append(Appear(label, *_noisy(rng, c, 0.5, 0.05)))
            label += 1
        elif r < (0.7 if grow else 0.45):
            s = rng.choice(alive)
            c = counts.pop(s)
            counts[label] = rng.randint(0, c)
            counts[label + 1] = c - counts[label]
            events.append(Split(s, label, label + 1))
            label += 2
        elif r < 0.85:
            a, b = rng.sample(alive, 2)
            counts[label] = counts.pop(a) + counts.pop(b)
            events.append(Merge(a, b, label))
            label += 1
        else:
            s = rng.choice(alive)
            c = counts.pop(s)
            events.append(Disappear(s, *_noisy(rng, c, 0.6, 0.0)))
    return Instance(ShadowSequence(initial, events), None, truth_left, dict(counts))


def perturb(rng: random.Random, seq: ShadowSequence) -> ShadowSequence:
    """Shift one observation (an initial range, appear/disappear range or exit size)."""
    d = seq.to_dict()
    sites = [("initial", k) for k in d["initial"]]
    sites += [("event", i) for i, e in enumerate(d["events"]) if e["type"] in ("appear", "disappear", "exit")]
    kind, key = rng.choice(sites)
    shift = rng.choice([-3, -2, -1, 1, 2, 3])
    if kind == "initial":
        lo, hi = d["initial"][key]
        lo = max(0, lo + shift)
        d["initial"][key] = [lo, None if hi is None else max(lo, hi + shift)]
    else:
        e = d["events"][key]
        if e["type"] == "exit":
            e["k"] = max(1, e["k"] + abs(shift))
        else:
            e["lo"] = max(0, e["lo"] + shift)
            e["hi"] = None if e["hi"] is None else max(e["lo"], e["hi"] + shift)
    return ShadowSequence.from_dict(d)


def _descendants(seq: ShadowSequence) -> Dict[int, set]:
    """Right vertices (as tagged keys of a ``fov="naive"`` filter) reachable from each label."""
    desc: Dict[int, set] = defaultdict(set)
    for s in seq.alive_at_end():
        desc[s].add(("alive", s))
    n = sum(isinstance(e, (Enter, Exit)) for e in seq.events) + 1
    for e in reversed(seq.events):
        if isinstance(e, (Enter, Exit)):
            n -= 1
            if isinstance(e, Exit):
                desc[e.s].add(("gone", Pseudo("exit", e.s, n)))
        elif isinstance(e, Disappear):
            desc[e.s].add(("gone", e.s))
        elif isinstance(e, Split):
            desc[e.s] |= desc[e.a] | desc[e.b]
        elif isinstance(e, Merge):
            desc[e.a] |= desc[e.s]
            desc[e.b] |= desc[e.s]
    return desc


class ReplayError(AssertionError):
    """A witness does not describe a valid evolution of targets."""


def replay_witness(seq: ShadowSequence, supply: Mapping, flow: Mapping,
                   total: Optional[Bound] = None) -> Dict[int, int]:
    """Replay a ``fov="naive"`` witness through the events; return final counts per alive shadow.

    Each unit of ``flow[(left, right)]`` is a target that starts in ``left`` and
    must end in ``right``; at a split it follows the child from which ``right``
    is still reachable.  Every observation is checked along the way.
    """
    out_of: Dict[object, Counter] = defaultdict(Counter)
    for (u, r), k in flow.items():
        if k < 0:
            raise ReplayError(f"negative flow on {(u, r)}")
        out_of[u][r] += k

    def start(v: object, lo: float, hi: float) -> Counter:
        c = Counter(out_of.pop(v, Counter()))
        n = sum(c.values())
        if n != supply.get(v, 0) or not lo <= n <= hi:
            raise ReplayError(f"left vertex {v}: {n} targets, supply {supply.get(v)}, bounds {(lo, hi)}")
        return c

    desc = _descendants(seq)
    content: Dict[int, Counter] = {s: start(s, lo, hi) for s, (lo, hi) in seq.initial.items()}
    if total is not None and not total[0] <= sum(sum(c.values()) for c in content.values()) <= total[1]:
        raise ReplayError("total constraint violated")
    n = 0
    for e in seq.events:
        if isinstance(e, Appear):
            content[e.s] = start(e.s, e.lo, e.hi)
        elif isinstance(e, Disappear):
            c = content.pop(e.s)
            if set(c) - {("gone", e.s)} or not e.lo <= c[("gone", e.s)] <= e.hi:
                raise ReplayError(f"{e}: content {dict(c)}")
        elif isinstance(e, Split):
            c = content.pop(e.s)
            ca, cb = Counter(), Counter()
            for r, k in c.items():
                if r in desc[e.a]:
                    ca[r] += k
                elif r in desc[e.b]:
                    cb[r] += k
                else:
                    raise ReplayError(f"{e}: target bound for {r} has no route")
            content[e.a], content[e.b] = ca, cb
        elif isinstance(e, Merge):
            content[e.s] = content.pop(e.a) + content.pop(e.b)
        elif isinstance(e, Enter):
            n += 1
            content[e.s] += start(Pseudo("enter", e.s, n), e.k, e.k)
        elif isinstance(e, Exit):
            n += 1
            key = ("gone", Pseudo("exit", e.s, n))
            if content[e.s][key] != e.k:
                raise ReplayError(f"{e}: {content[e.s][key]} targets leave")
            del content[e.s][key]
    if out_of and any(sum(c.values()) for c in out_of.values()):
        raise ReplayError(f"flow from unknown left vertices {list(out_of)}")
    final = {}
    for s, c in content.items():
        if set(c) - {("alive", s)}:
            raise ReplayError(f"shadow {s} ends with misrouted targets {dict(c)}")
        final[s] = c[("alive", s)]
    return final


__all__ = ["Instance", "random_instance", "growing_instance", "perturb", "replay_witness", "ReplayError"]
