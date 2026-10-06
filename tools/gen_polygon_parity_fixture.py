"""Write ``tests/fixtures/polygon_parity.json`` for the Python/JS parity tests of the polygon world.

For each of the 14 maps and its default demo path (``data/paths.json`` section ``demo``) the fixture holds:

* ``histories`` -- ``get_gaps(compat="fixed", java_matching=False)`` of the path (the hidden-chain labels the
  simulator tracks with, :func:`~shadowinfo.polygon.simulate._track`): every gap set (``[id, iEdge, toGapSet]``),
  the set times and sample indices, the IDs at every critical point, the repaired steps and the
  component events with their times (:func:`~shadowinfo.polygon.events.gap_history_to_events`);
* ``runs`` -- :class:`~shadowinfo.polygon.simulate.PolygonSimulator` runs with two seeds: seed 1 runs one
  full pass of the path in about 160 ticks and a few ticks into the next leg (loop or reverse), seed 2
  runs 120 ticks at the default speed with more targets, and on four maps a third run drives the robot
  off its path (``go_to`` / ``follow``).  Every frame is stored (robot, sensing point, labels, per-target
  ground truth, events, created counts), compacted by :func:`compact_frame` (target positions only every
  ``TARGET_EVERY`` ticks and on the last frame -- any divergence of the random walk shows up in the
  per-tick ground truth anyway --, the sensing point only where it is not the robot, no ``counts``,
  which follow from ``labels`` and ``target_labels``, and no empty ``events`` / ``created``);
* ``geometry`` -- the shadow pockets, windows and visibility polygon at a few ticks of the seed-1 runs.

``js-tests/polygon_parity.test.js`` recomputes all of it in Node and deep-compares (exact doubles).
``inputs_sha256`` (:func:`input_fingerprint`) hashes everything the numbers are computed from (the
``shadowinfo`` sources, the polygon data and this script); ``tests/test_polygon_js.py`` regenerates the
fixture and compares it whenever the fingerprint differs, so the slow regeneration runs only after a change.

Usage: ``python tools/gen_polygon_parity_fixture.py``.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo.events import event_to_dict  # noqa: E402
from shadowinfo.polygon import (PolygonSimulator, gap_history_to_events, load_demo_path, load_polygon,  # noqa: E402
                                validate_path)
from shadowinfo.polygon.simulate import _track  # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "polygon_parity.json"
MAPS = tuple(range(1, 15))
SEEDS = (1, 2)
FAST_TICKS = 160  # seed 1: speed chosen so that one pass takes about this many ticks
EXTRA_TICKS = 12  # ... plus this many into the next leg
SLOW_TICKS = 120  # seed 2: default speed
N_TARGETS = {1: 10, 2: 15}
TARGET_EVERY = 25
GEOMETRY_TICKS = (0, 70, 140)
SCRIPTED = (3, 8, 12, 14)  # maps with a go_to / follow run


def input_files() -> list:
    """What the fixture is computed from: ``shadowinfo/**/*.py``, ``shadowinfo/polygon/data/*`` and this script."""
    data = ROOT / "shadowinfo" / "polygon" / "data"
    return (sorted((ROOT / "shadowinfo").rglob("*.py")) + sorted(p for p in data.iterdir() if p.is_file())
            + [Path(__file__).resolve()])


def input_fingerprint() -> str:
    """sha256 over the repository paths and bytes of :func:`input_files`."""
    h = hashlib.sha256()
    for f in input_files():
        h.update(f.relative_to(ROOT).as_posix().encode() + b"\0" + f.read_bytes() + b"\0")
    return h.hexdigest()


def _bounds_json(b: dict) -> dict:
    return {str(s): [lo, None if hi == math.inf else hi] for s, (lo, hi) in b.items()}


def compact_frame(f: dict, last: int) -> dict:
    """The stored form of a frame (``js-tests/polygon_parity.test.js`` applies the same rules)."""
    g = dict(f)
    if f["t"] % TARGET_EVERY and f["t"] != last:
        del g["targets"]
    if g["sense"] == g["robot"]:
        del g["sense"]
    del g["counts"]
    for k in ("events", "created"):
        if k in g and not g[k]:
            del g[k]
    return g


def _history(n: int) -> dict:
    poly = load_polygon(n)
    path, at_end = load_demo_path(n)
    h = _track(poly, path).history  # = get_gaps(poly, path, java_matching=False), cached like the simulator's
    ge = gap_history_to_events(h)
    return {
        "map": n,
        "waypoints": [list(p) for p in path.points],
        "at_end": at_end,
        "n_critical_points": len(h.sample_ids),
        "critical_distances": [cp.distance for cp in h.samples.critical_points],
        "sets": [[g.to_list() for g in gs] for gs in h],
        "times": h.times,
        "sample_indices": h.sample_indices,
        "sample_ids": h.sample_ids,
        "inferred": [list(x) for x in h.inferred],
        "initial": ge.initial,
        "events": [event_to_dict(e) for e in ge.events],
        "event_times": ge.times,
    }


def _first_valid(poly, robot, candidates):
    for c in candidates:
        c = (float(c[0]), float(c[1]))
        if c == tuple(robot):
            continue
        try:
            validate_path(poly, [robot, c])
        except ValueError:
            continue
        return c
    raise RuntimeError("no valid go_to target")


def _run(n: int, seed: int, scripted: bool = False) -> dict:
    poly = load_polygon(n)
    path, _ = load_demo_path(n)
    if seed == 1 and not scripted:
        speed = math.ceil(path.length() / FAST_TICKS * 10) / 10
        ticks = math.ceil(path.length() / speed) + EXTRA_TICKS
    else:
        speed, ticks = 4.0, SLOW_TICKS
    nt = N_TARGETS[seed]
    sim = PolygonSimulator(poly, None, nt, seed, speed)
    run = {"map": n, "seed": seed, "n_targets": nt, "speed": speed, "ticks": ticks, "script": [],
           "initial": {m: _bounds_json(sim.initial_condition(m)) for m in ("exact", "unknown", "evader")},
           "initial_labels": sim.initial_labels}
    frames = [sim.frame()]
    geometry = {}
    if seed == 1 and not scripted:
        geometry["0"] = {k: v for k, v in sim.frame(geometry=True).items() if k in ("shadows", "visibility")}
    wps = list(path.points)
    for t in range(1, ticks + 1):
        if scripted and t == 40:  # click-to-move to the first reachable waypoint (from the end of the path)
            c = _first_valid(poly, sim.robot, wps[::-1])
            sim.go_to(c)
            run["script"].append([t, "go_to", list(c)])
        if scripted and t == 80:  # then follow two waypoints
            a = _first_valid(poly, sim.robot, wps)
            b = _first_valid(poly, a, wps[::-1])
            sim.follow([a, b], "reverse")
            run["script"].append([t, "follow", [list(a), list(b)], "reverse"])
        frames.append(sim.step())
        if t in GEOMETRY_TICKS and seed == 1 and not scripted:
            geometry[str(t)] = {k: v for k, v in sim.frame(geometry=True).items() if k in ("shadows", "visibility")}
    run["frames"] = [compact_frame(f, ticks) for f in frames]
    run["stats"] = dict(sim.stats)
    run["n_history"] = len(sim.history)
    if geometry:
        run["geometry"] = geometry
    return run


def build_fixture() -> dict:
    runs = []
    for n in MAPS:
        for seed in SEEDS:
            runs.append(_run(n, seed))
        if n in SCRIPTED:
            runs.append(_run(n, 2, scripted=True))
    return {
        "_doc": "Python/JS parity of the polygon world (tools/gen_polygon_parity_fixture.py); "
                "frames are compacted by compact_frame: 'targets' only every %d ticks and on the last frame, "
                "'sense' only where it differs from 'robot', no 'counts', no empty 'events'/'created'." % TARGET_EVERY,
        "inputs_sha256": input_fingerprint(),
        "histories": [_history(n) for n in MAPS],
        "runs": runs,
    }


def main() -> None:
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text(json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False) + "\n")
    print(f"wrote {FIXTURE.relative_to(ROOT)} ({FIXTURE.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
