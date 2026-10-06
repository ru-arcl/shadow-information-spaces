"""Robot + random-walk targets in the grid world (DESIGN §5).

Produces the critical-event stream of T-RO 2012, Secs. II-III, together with
the ground truth needed to check the filters: per-shadow target counts after
every tick and the true counts of every shadow created by a split or merge.

One tick:

1. *sensor substep* -- the robot advances one cell along its looped path.
   Targets that were in a shadow and are now visible emit ``Exit(old)``
   before the component events (unless their shadow disappears, in which
   case the count is revealed by ``Disappear``); then the component events;
   then targets that were visible and are now in a shadow emit
   ``Enter(new)`` (unless the shadow just appeared, in which case the count
   is carried by ``Appear``).
2. *target substep* -- each target draws ``randint(5)``: stay, up, down,
   left or right (blocked moves stay).  shadow -> visible emits ``Exit``,
   visible -> shadow emits ``Enter``; a 4-neighbour move between two shadow
   cells never changes the shadow.

FOV events of one substep are aggregated per shadow (``Exit(s, k)``), all
exits sorted by label first, then all enters sorted by label.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

from ..events import Appear, Bound, Disappear, Enter, Event, Exit, ShadowSequence, event_to_dict
from ..rng import Rng
from .gridmap import Cell, GridMap
from .shadows import ShadowTracker

MOVES = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))
MODES = ("exact", "unknown", "evader")


def label_hash(labels: Sequence[int]) -> int:
    """FNV-1a style 32-bit hash of a label grid (one step per cell)."""
    h = 2166136261
    for v in labels:
        h = ((h ^ v) * 16777619) & 0xFFFFFFFF
    return h


def initial_bounds(counts: Dict[int, int], mode: str = "exact") -> Dict[int, Bound]:
    """Initial condition for a :class:`ShadowSequence` from true counts.

    ``"exact"``: ``lo = hi =`` true count; ``"unknown"``: ``[0, inf)``
    (counting, Sec. V-F); ``"evader"``: ``[0, 1]`` per shadow (pursuit-evasion).
    """
    if mode == "exact":
        return {s: (n, n) for s, n in counts.items()}
    if mode == "unknown":
        return {s: (0, math.inf) for s in counts}
    if mode == "evader":
        return {s: (0, 1) for s in counts}
    raise ValueError(f"unknown mode {mode!r}; expected one of {MODES}")


def _tally(labels: Sequence[int], targets: Sequence[int]) -> Dict[int, int]:
    counts = {s: 0 for s in set(labels) - {0}}
    for i in targets:
        if labels[i]:
            counts[labels[i]] += 1
    return counts


def _fov_events(exits: Dict[int, int], enters: Dict[int, int]) -> List[Event]:
    return [Exit(s, k) for s, k in sorted(exits.items())] + [Enter(s, k) for s, k in sorted(enters.items())]


class GridSimulator:
    """Deterministic simulator; identical results to ``SI.GridSimulator`` in JS.

    ``frame()`` describes the current state; ``step()`` advances one tick and
    returns the frame, which additionally holds that tick's ``events`` and
    ``created``: ``[[label, n], ...]`` for every shadow created by a split or
    merge this tick, ``n`` being its true count at the moment the component
    events fire (after the exits, before the enters of the sensor substep).
    Frames are plain JSON-ready dicts; Python and JS produce identical ones.
    """

    def __init__(self, gridmap: GridMap, path: Sequence[Cell], n_targets: int = 10,
                 seed: int = 1, radius: Optional[int] = None) -> None:
        if not path:
            raise ValueError("empty path")
        self.map = gridmap
        self.path = [tuple(p) for p in path]
        self.radius = radius
        self.rng = Rng(seed)
        self.t = 0
        self.robot = self.path[0]
        self.tracker = ShadowTracker(gridmap.rows, gridmap.cols)
        self.tracker.reset(gridmap.shadow_mask(self.robot, radius))
        free = gridmap.free_cells
        self.targets: List[int] = [free[self.rng.randint(len(free))] for _ in range(n_targets)]
        self.initial_counts = _tally(self.tracker.labels, self.targets)
        self.history: List[Event] = []

    @property
    def labels(self) -> List[int]:
        return self.tracker.labels

    def counts(self) -> Dict[int, int]:
        return _tally(self.tracker.labels, self.targets)

    def initial_condition(self, mode: str = "exact") -> Dict[int, Bound]:
        return initial_bounds(self.initial_counts, mode)

    def sequence(self, mode: str = "exact") -> ShadowSequence:
        """All events so far, with the initial condition for ``mode``."""
        return ShadowSequence(self.initial_condition(mode), list(self.history))

    def frame(self) -> dict:
        cols = self.map.cols
        return {
            "t": self.t,
            "robot": list(self.robot),
            "targets": [list(divmod(i, cols)) for i in self.targets],
            "counts": [[s, n] for s, n in sorted(self.counts().items())],
            "labels_hash": label_hash(self.tracker.labels),
        }

    def step(self) -> dict:
        self.t += 1
        self.robot = self.path[self.t % len(self.path)]
        old = self.tracker.labels
        tr = self.tracker.update(self.map.shadow_mask(self.robot, self.radius))
        new = tr.labels
        appeared, disappeared = set(tr.appeared), set(tr.disappeared)

        exits: Dict[int, int] = {}
        enters: Dict[int, int] = {}
        revealed: Dict[int, int] = {}
        hidden: Dict[int, int] = {}
        for i in self.targets:
            a, b = old[i], new[i]
            if a and not b:
                if a in disappeared:
                    revealed[a] = revealed.get(a, 0) + 1
                else:
                    exits[a] = exits.get(a, 0) + 1
            elif b and not a:
                if b in appeared:
                    hidden[b] = hidden.get(b, 0) + 1
                else:
                    enters[b] = enters.get(b, 0) + 1

        comp_count = [0] * len(tr.cells)
        comp_of = {s: k for k, s in enumerate(tr.comp_labels)}
        stay: Dict[int, int] = {}
        for i in self.targets:
            if old[i] and new[i]:
                comp_count[comp_of[new[i]]] += 1
                stay[old[i]] = stay.get(old[i], 0) + 1
        created = {s: sum(comp_count[k] for k in ks) for s, ks in tr.covers.items()}
        for s, olds in tr.merged_from.items():
            created[s] = sum(stay.get(o, 0) for o in olds)

        comp_events: List[Event] = []
        for e in tr.events:
            if isinstance(e, Appear):
                n = hidden.get(e.s, 0)
                e = Appear(e.s, n, n)
            elif isinstance(e, Disappear):
                n = revealed.get(e.s, 0)
                e = Disappear(e.s, n, n)
            comp_events.append(e)
        events = _fov_events(exits, {}) + comp_events + _fov_events({}, enters)

        exits, enters = {}, {}
        rows, cols, blocked = self.map.rows, self.map.cols, self.map.blocked
        for j, i in enumerate(self.targets):
            dr, dc = MOVES[self.rng.randint(5)]
            r, c = divmod(i, cols)
            r, c = r + dr, c + dc
            if not (0 <= r < rows and 0 <= c < cols) or blocked[r * cols + c]:
                continue
            k = r * cols + c
            a, b = new[i], new[k]
            if a and not b:
                exits[a] = exits.get(a, 0) + 1
            elif b and not a:
                enters[b] = enters.get(b, 0) + 1
            elif a != b:
                raise AssertionError(f"target moved between shadows {a} and {b} without an event")
            self.targets[j] = k
        events += _fov_events(exits, enters)
        self.history.extend(events)

        f = self.frame()
        f["events"] = [event_to_dict(e) for e in events]
        f["created"] = [[s, n] for s, n in sorted(created.items())]
        return f

    def run(self, ticks: int) -> List[dict]:
        """Frame 0 followed by ``ticks`` stepped frames."""
        frames = [self.frame()]
        for _ in range(ticks):
            frames.append(self.step())
        return frames


def simulate(gridmap: GridMap, path: Sequence[Cell], n_targets: int, seed: int, ticks: int,
             radius: Optional[int] = None) -> Tuple[GridSimulator, List[dict]]:
    """Convenience wrapper: build a simulator and run it for ``ticks`` ticks."""
    sim = GridSimulator(gridmap, path, n_targets, seed, radius)
    return sim, sim.run(ticks)


__all__ = ["GridSimulator", "simulate", "initial_bounds", "label_hash", "MOVES", "MODES"]
