/* Python/JS parity: recompute tests/fixtures/parity_office.json in Node
 * (written by tools/gen_parity_fixture.py) and deep-compare. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

require("../docs/js/rng.js");
require("../docs/js/grid.js");
const SI = require("../docs/js/map_office.js");

const fixture = JSON.parse(
  fs.readFileSync(path.join(__dirname, "..", "tests", "fixtures", "parity_office.json"), "utf8"));
const spec = SI.maps.office;
const gmap = SI.GridMap.fromText(spec.text);
const lap = SI.expandPath(spec.waypoints, spec.loop);
const roundtrip = (x) => JSON.parse(JSON.stringify(x));

test("embedded map and path match shadowinfo/grid/maps/ (via fixture)", () => {
  assert.equal(gmap.rows, fixture.rows);
  assert.equal(gmap.cols, fixture.cols);
  assert.deepEqual(lap, fixture.path);
});

test("RNG streams", () => {
  for (const s of fixture.rng) {
    const r = new SI.Rng(s.seed);
    const u32 = Array.from({ length: 8 }, () => r.nextU32());
    const r7 = Array.from({ length: 8 }, () => r.randint(7));
    const r1000 = Array.from({ length: 8 }, () => r.randint(1000));
    assert.deepEqual({ seed: s.seed, u32, randint7: r7, randint1000: r1000 }, s);
  }
});

test("visibility masks and label grids", () => {
  const tracker = new SI.ShadowTracker(gmap.rows, gmap.cols);
  tracker.reset(gmap.shadowMask(lap[0]));
  const ticks = fixture.snapshots.map((s) => s.t);
  const last = Math.max(...ticks);
  let k = 0;
  for (let t = 0; t <= last; t++) {
    if (t) tracker.update(gmap.shadowMask(lap[t % lap.length]));
    if (t !== ticks[k]) continue;
    const snap = fixture.snapshots[k++];
    const robot = lap[t % lap.length];
    assert.deepEqual(robot, snap.robot);
    assert.equal(gmap.visibility(robot).join(""), snap.visible, `visibility at t=${t}`);
    assert.equal(gmap.visibility(robot, 12).join(""), snap.visible_r12, `visibility r=12 at t=${t}`);
    assert.deepEqual(tracker.labels, snap.labels, `labels at t=${t}`);
  }
});

for (const run of fixture.runs) {
  test(`simulator frames: seed=${run.seed} radius=${run.radius} n=${run.n_targets}`, () => {
    const sim = new SI.GridSimulator(gmap, lap, run.n_targets, run.seed, run.radius);
    const initial = {};
    for (const m of ["exact", "unknown", "evader"]) initial[m] = roundtrip(sim.initialCondition(m));
    assert.deepEqual(initial, run.initial);
    const frames = roundtrip(sim.run(lap.length));
    assert.equal(frames.length, run.frames.length);
    for (let t = 0; t < frames.length; t++) assert.deepEqual(frames[t], run.frames[t], `frame ${t}`);
  });
}
