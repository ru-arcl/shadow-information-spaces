"""Write ``tests/fixtures/parity_office.json`` for the Python/JS parity tests.

The fixture holds, for the fixed office map and robot path, the RNG streams,
a few full visibility masks and label grids, and every simulator frame of one
full lap for several seeds (``js-tests/parity.test.js`` recomputes all of it
in Node and deep-compares; ``tests/test_grid.py`` recomputes it in Python).

Usage: ``python tools/gen_parity_fixture.py``.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from shadowinfo.grid import GridSimulator, ShadowTracker, load_office  # noqa: E402
from shadowinfo.rng import Rng  # noqa: E402

FIXTURE = ROOT / "tests" / "fixtures" / "parity_office.json"
SEEDS = (1, 2, 3)
N_TARGETS = 30
RADIUS_RUNS = ((1, 12), (7, 6))  # (seed, radius)
SNAPSHOT_TICKS = (0, 47, 94, 141)


def _bounds_json(b: dict) -> dict:
    return {str(s): [lo, None if hi == math.inf else int(hi)] for s, (lo, hi) in sorted(b.items())}


def _run(gmap, path, seed: int, radius, n_targets: int = N_TARGETS) -> dict:
    sim = GridSimulator(gmap, path, n_targets, seed, radius)
    initial = {m: _bounds_json(sim.initial_condition(m)) for m in ("exact", "unknown", "evader")}
    frames = sim.run(len(path))
    return {"seed": seed, "radius": radius, "n_targets": n_targets, "initial": initial, "frames": frames}


def build_fixture() -> dict:
    gmap, path = load_office()
    rng_streams = []
    for seed in (0, 1, 2, 3, 12345, 0xFFFFFFFF):
        r = Rng(seed)
        rng_streams.append({
            "seed": seed,
            "u32": [r.next_u32() for _ in range(8)],
            "randint7": [r.randint(7) for _ in range(8)],
            "randint1000": [r.randint(1000) for _ in range(8)],
        })
    tracker = ShadowTracker(gmap.rows, gmap.cols)
    tracker.reset(gmap.shadow_mask(path[0]))
    snapshots = []
    for t in range(max(SNAPSHOT_TICKS) + 1):
        if t:
            tracker.update(gmap.shadow_mask(path[t % len(path)]))
        if t in SNAPSHOT_TICKS:
            robot = path[t % len(path)]
            snapshots.append({
                "t": t,
                "robot": list(robot),
                "visible": "".join(map(str, gmap.visibility(robot))),
                "visible_r12": "".join(map(str, gmap.visibility(robot, 12))),
                "labels": list(tracker.labels),
            })
    runs = [_run(gmap, path, s, None) for s in SEEDS]
    runs += [_run(gmap, path, s, r) for s, r in RADIUS_RUNS]
    runs.append(_run(gmap, path, 1, None, 0))  # no targets: component events only
    return {
        "map": "office",
        "rows": gmap.rows,
        "cols": gmap.cols,
        "path": [list(p) for p in path],
        "rng": rng_streams,
        "snapshots": snapshots,
        "runs": runs,
    }


def main() -> None:
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    FIXTURE.write_text(json.dumps(build_fixture(), separators=(",", ":")) + "\n")
    print(f"wrote {FIXTURE.relative_to(ROOT)} ({FIXTURE.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
