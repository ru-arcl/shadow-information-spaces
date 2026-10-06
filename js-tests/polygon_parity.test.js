/* Python/JS parity of the polygon world: recompute tests/fixtures/polygon_parity.json (written by
 * tools/gen_polygon_parity_fixture.py) in Node and deep-compare, exact doubles throughout.
 * - histories: get_gaps(compat="fixed", java_matching=False) (what the simulator tracks with) of the
 *   default path of every map (gap sets, times, IDs at every critical point, repaired steps,
 *   component events);
 * - runs: PolygonSimulator frames tick by tick (robot, sensing point, labels, ground truth, events,
 *   created counts, target positions), two seeds per map plus click-to-move runs;
 * - geometry: shadow pockets, windows and visibility polygons at a few ticks. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

require("../docs/js/rng.js");
require("../docs/js/maps_polygon.js");
require("../docs/js/polygon.js");
const SI = require("../docs/js/polygon_sim.js");
const P = SI.polygon;

const fixture = JSON.parse(
  fs.readFileSync(path.join(__dirname, "..", "tests", "fixtures", "polygon_parity.json"), "utf8"));
const roundtrip = (x) => JSON.parse(JSON.stringify(x));
const TARGET_EVERY = Number(/every (\d+) ticks/.exec(fixture._doc)[1]);

/* gen_polygon_parity_fixture.compact_frame */
function compactFrame(f, last) {
  const g = Object.assign({}, f);
  if (f.t % TARGET_EVERY && f.t !== last) delete g.targets;
  if (g.sense[0] === g.robot[0] && g.sense[1] === g.robot[1]) delete g.sense;
  delete g.counts;
  for (const k of ["events", "created"]) if (k in g && !g[k].length) delete g[k];
  return g;
}

for (const hx of fixture.histories) {
  test(`map ${hx.map}: gap history of the default path`, () => {
    const rec = SI.polygonMaps.demo[String(hx.map)];
    assert.deepEqual(rec.waypoints, hx.waypoints);
    assert.equal(rec.at_end, hx.at_end);
    const poly = P.loadPolygon(hx.map);
    const h = P.getGaps(poly, new P.Path(rec.waypoints), { javaMatching: false });
    assert.deepEqual(h.samples.criticalPoints.map((cp) => cp.distance), hx.critical_distances);
    assert.deepEqual(h.sets.map((gs) => gs.map((g) => g.toList())), hx.sets);
    assert.deepEqual(h.times, hx.times);
    assert.deepEqual(h.sampleIndices, hx.sample_indices);
    assert.deepEqual(h.sampleIds, hx.sample_ids);
    assert.deepEqual(h.inferred, hx.inferred);
    const ge = P.gapHistoryToEvents(h);
    assert.deepEqual(ge.initial, hx.initial);
    assert.deepEqual(roundtrip(ge.events), hx.events);
    assert.deepEqual(ge.times, hx.event_times);
  });
}

for (const run of fixture.runs) {
  const what = run.script.length ? "click-to-move" : `speed ${run.speed}`;
  test(`map ${run.map} seed ${run.seed} (${run.n_targets} targets, ${what}): simulator frames`, () => {
    const sim = new SI.PolygonSimulator(run.map, null, run.n_targets, run.seed, run.speed);
    const initial = {};
    for (const m of ["exact", "unknown", "evader"]) initial[m] = roundtrip(sim.initialCondition(m));
    assert.deepEqual(initial, run.initial);
    assert.deepEqual(sim.initialLabels, run.initial_labels);
    const check = (f) => {
      const counts = new Map(f.labels.map((s) => [s, 0]));
      f.target_labels.forEach((s) => { if (s) counts.set(s, counts.get(s) + 1); });
      assert.deepEqual(f.counts, [...counts.entries()].sort((a, b) => a[0] - b[0]), `counts at t=${f.t}`);
      assert.deepEqual(roundtrip(compactFrame(f, run.ticks)), run.frames[f.t], `frame ${f.t}`);
      if (run.geometry && run.geometry[String(f.t)]) {
        const g = sim.frame(true);
        assert.deepEqual(roundtrip({ shadows: g.shadows, visibility: g.visibility }), run.geometry[String(f.t)], `geometry ${f.t}`);
      }
    };
    check(sim.frame());
    for (let t = 1; t <= run.ticks; t++) {
      for (const s of run.script) {
        if (s[0] !== t) continue;
        if (s[1] === "go_to") sim.goTo(s[2]);
        else sim.follow(s[2], s[3]);
      }
      check(sim.step());
    }
    assert.deepEqual(sim.stats, run.stats);
    assert.equal(sim.history.length, run.n_history);
  });
}

test("simulator API: path validation, stored paths, click-to-move", () => {
  const S = SI.polygonSim;
  const poly = P.loadPolygon(12);
  assert.throws(() => S.validatePath(poly, [[90, 450]]), /at least 2 waypoints/);
  assert.throws(() => S.validatePath(poly, [[90, 450], [90, 450]]), /zero length/);
  assert.throws(() => S.validatePath(poly, [[90, 450], [-5, 450]]), /not strictly inside/);
  assert.throws(() => S.validatePath(poly, [[90, 450], [700, 100]]), /leaves polygon/);
  assert.deepEqual(S.validatePath(poly, [[90, 450], [580, 450]]), [[90, 450], [580, 450]]);
  const sim = new SI.PolygonSimulator(12, "P12", 5, 3);
  assert.equal(sim.atEnd, "stop");
  assert.throws(() => new SI.PolygonSimulator(13, "P12"), /belongs to polygon 12/);
  assert.throws(() => sim.goTo([700, 100]), /leaves polygon/);
  const frames = sim.run(sim.ticksPerLap() + 3);
  const last = frames[frames.length - 1];
  assert.deepEqual(last.robot, [520, 840]); /* stopped at the end of the path */
  sim.goTo([580, 840]);
  sim.step();
  assert.equal(sim.legs.length, 2);
  assert.deepEqual(sim.sequence("unknown").initial, { 1: [0, Infinity], 2: [0, Infinity], 3: [0, Infinity], 4: [0, Infinity] });
});

test("pockets and the visibility polygon tile every map", () => {
  const S = SI.polygonSim;
  for (let n = 1; n <= 14; n++) {
    const sim = new SI.PolygonSimulator(n, null, 0, 1, 37);
    for (let t = 0; t < 3; t++) {
      const area = S.polygonArea(P.loadPolygon(n).vertices);
      const parts = sim.shadows().reduce((acc, sh) => acc + S.polygonArea(sh.pocket), S.polygonArea(sim.visibilityPolygon()));
      assert.ok(Math.abs(parts - area) <= 1e-9 * area, `map ${n} t=${t}: ${parts} vs ${area}`);
      sim.shadows().forEach((sh) => assert.ok(S.polygonArea(sh.pocket) >= 0));
      sim.step();
    }
  }
});
