/* Smoke test for the demo controller (docs/js/controller.js): the same
 * simulation + filter loop the page runs, checked for soundness, on the grid
 * office map and on the 14 polygon maps of the original Java code. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const SI = require("../docs/js/controller.js");

const LAP = SI.expandPath(SI.maps.office.waypoints, SI.maps.office.loop).length;
const MODES = ["exact", "unknown", "evader"];
const SEEDS = [1, 2, 3];

/* One tick's soundness checks: bounds contain the truth, E[n] inside the bounds and summing to the
 * hidden total, the I-state compacted, and the bipartite view consistent with the shadows. */
function checkTick(c, tag) {
  const mode = c.opts.mode;
  for (const sh of c.shadows) {
    assert.ok(sh.lo <= sh.truth && sh.truth <= sh.hi, `${tag} t=${c.t} s${sh.label}: ${sh.truth} not in [${sh.lo}, ${sh.hi}]`);
    if (sh.expected !== null) {
      assert.ok(sh.expected >= sh.lo - 1e-9 && sh.expected <= sh.hi + 1e-9, `${tag} t=${c.t} s${sh.label}: E=${sh.expected} outside bounds`);
    }
  }
  const h = c.hidden;
  assert.ok(h.lo <= h.truth && h.truth <= h.hi, `${tag} t=${c.t}: hidden ${h.truth} not in [${h.lo}, ${h.hi}]`);
  if (mode !== "unknown") assert.equal(h.lo, h.hi, `${tag}: hidden total is known exactly in this mode`);
  if (c.prob) {
    const sum = c.shadows.reduce((acc, sh) => acc + sh.expected, 0);
    assert.ok(Math.abs(sum - h.truth) < 1e-6, `${tag} t=${c.t}: sum of expectations ${sum} != ${h.truth}`);
  }
  assert.ok(c.edges <= c.opts.hardEdges + 400, `${tag}: I-state not compacted: ${c.edges} edges`);
  const view = c.bipartiteView();
  assert.ok(view.right.filter((r) => r.kind === "alive").length === c.shadows.length);
  for (const [li, ri] of view.edges) assert.ok(view.left[li] && view.right[ri]);
}

for (const mode of MODES) {
  test(`3 laps x 3 seeds, mode ${mode}: bounds contain the truth`, () => {
    for (const seed of SEEDS) {
      const c = new SI.Controller({ nTargets: 10, seed, mode });
      for (let t = 0; t < 3 * LAP; t++) {
        c.step();
        for (const sh of c.shadows) {
          assert.ok(sh.lo <= sh.truth && sh.truth <= sh.hi,
                    `seed ${seed} t=${c.t} s${sh.label}: ${sh.truth} not in [${sh.lo}, ${sh.hi}]`);
          if (sh.expected !== null) {
            assert.ok(sh.expected >= sh.lo - 1e-9 && sh.expected <= sh.hi + 1e-9,
                      `seed ${seed} t=${c.t} s${sh.label}: E=${sh.expected} outside bounds`);
          }
        }
        const h = c.hidden;
        assert.ok(h.lo <= h.truth && h.truth <= h.hi, `seed ${seed} t=${c.t}: hidden ${h.truth} not in [${h.lo}, ${h.hi}]`);
        if (mode !== "unknown") assert.equal(h.lo, h.hi, "hidden total is known exactly in this mode");
        if (c.prob) {
          const sum = c.shadows.reduce((acc, sh) => acc + sh.expected, 0);
          assert.ok(Math.abs(sum - h.truth) < 1e-6, `t=${c.t}: sum of expectations ${sum} != ${h.truth}`);
        }
        assert.ok(c.edges <= c.opts.hardEdges + 400, `I-state not compacted: ${c.edges} edges`);
        const view = c.bipartiteView();
        assert.ok(view.right.filter((r) => r.kind === "alive").length === c.shadows.length);
        for (const [li, ri] of view.edges) assert.ok(view.left[li] && view.right[ri]);
      }
      assert.deepEqual(c.violations, []);
      assert.deepEqual(c.errors, []);
      assert.equal(c.t, 3 * LAP);
      if (mode === "unknown") assert.equal(c.prob, null);
      if (mode === "evader") assert.equal(c.sim.targets.length, 1);
    }
  });
}

test("exact re-rooting keeps the bounds of the uncompacted filter", () => {
  for (const mode of MODES) {
    const c = new SI.Controller({ nTargets: 12, seed: 4, mode, compactEdges: 300, checkEvery: 3 });
    const sim = SI.grid.simulator("office", c.sim.targets.length, 4, null);
    const ref = new SI.CombinatorialFilter(sim.initialCondition(mode), c.total);
    for (let t = 0; t < 2 * LAP; t++) {
      c.step();
      ref.extend(sim.step().events);
      if (c.compactions.relaxed) continue;
      const rb = ref.allBounds();
      for (const sh of c.shadows) assert.deepEqual([sh.lo, sh.hi], rb[sh.label], `${mode} t=${c.t} s${sh.label}`);
    }
    assert.ok(c.compactions.exact > 0, `${mode}: no re-root happened`);
  }
});

test("forced relaxed re-root stays sound", () => {
  const c = new SI.Controller({ nTargets: 15, seed: 5, mode: "unknown", compactEdges: 50, hardEdges: 60, maxCheckShadows: 0 });
  for (let t = 0; t < LAP; t++) {
    c.step();
    for (const sh of c.shadows) assert.ok(sh.lo <= sh.truth && sh.truth <= sh.hi);
  }
  assert.ok(c.compactions.relaxed > 0);
  assert.deepEqual(c.violations, []);
});

/* Polygon maps: the default path of each map in every mode.  Seed 1 drives at the simulator's
 * default 4 units per tick (up to 300 ticks); seed 2 at a stride that makes one pass about 200
 * ticks, run 10 ticks past the end of the pass (the turnaround of an open path, or the next lap). */
for (let map = 1; map <= 14; map++) {
  test(`map ${map}, default path, 3 modes x 2 seeds: bounds contain the truth`, () => {
    for (const mode of MODES) {
      for (const seed of [1, 2]) {
        const length = new SI.polygon.Path(SI.polygonMaps.demo[String(map)].waypoints).length();
        const stepLength = seed === 1 ? 4 : Math.max(4, Math.ceil(length / 200));
        const c = seed === 1 ? new SI.Controller({ map, mode, seed }) : new SI.Controller({ map, mode, seed, stepLength });
        assert.equal(c.world, "polygon");
        assert.equal(c.opts.radius, null);
        const ticks = seed === 1 ? Math.min(c.lapTicks + 10, 300) : c.lapTicks + 10;
        const tag = `map ${map} ${mode} seed ${seed}`;
        let events = 0;
        checkTick(c, tag);
        for (let t = 0; t < ticks; t++) {
          c.step();
          events += c.lastEvents.length;
          checkTick(c, tag);
          const truth = c.sim.counts();
          assert.deepEqual(c.shadows.map((sh) => sh.label), Array.from(truth.keys()).sort((a, b) => a - b), `${tag} t=${c.t}: labels`);
        }
        assert.deepEqual(c.violations, [], tag);
        assert.deepEqual(c.errors, [], tag);
        assert.ok(events > 0, `${tag}: no events`);
        if (mode === "unknown") assert.equal(c.prob, null);
        if (mode === "evader") assert.equal(c.sim.targets.length, 1);
      }
    }
  });
}

test("map list: the office grid plus the 14 original maps, paper figures noted", () => {
  const maps = SI.demo.mapList();
  assert.deepEqual(maps.map((m) => m.id), ["office"].concat(Array.from({ length: 14 }, (_, i) => String(i + 1))));
  assert.equal(maps[0].name, "Office (grid)");
  assert.equal(maps[12].name, "Map 12 \u2014 T-RO Fig. 11(a) / ICRA'08 Fig. 4");
  assert.equal(maps[13].name, "Map 13 \u2014 T-RO Fig. 15(a)");
  assert.equal(maps[14].name, "Map 14 \u2014 T-RO Fig. 15(b)");
  assert.equal(maps[3].name, "Map 3");
  for (const m of maps.slice(1)) assert.equal(m.world, "polygon");
});

test("path options: the demo path first, then every stored path of the map with its provenance", () => {
  const stored = SI.polygonMaps.paths;
  let total = 0;
  for (let map = 1; map <= 14; map++) {
    const opts = SI.demo.pathOptions(map);
    const demo = SI.polygonMaps.demo[String(map)];
    assert.equal(opts[0].id, "default");
    assert.deepEqual(opts[0].waypoints, demo.waypoints);
    const labels = Object.keys(stored).filter((k) => stored[k].polygon === map && k !== demo.source);
    assert.deepEqual(opts.slice(1).map((o) => o.id), labels);
    total += opts.length - 1 + (demo.source ? 1 : 0);
    for (const o of opts) {
      assert.ok(o.info.length > 20, `map ${map} ${o.id}: no provenance`);
      if (!o.valid) {
        assert.ok(["P13+", "P13+tail"].includes(o.id), `map ${map} ${o.id} should be valid: ${o.error}`);
        continue;
      }
      const c = new SI.Controller({ map, path: o.id, seed: 3 });
      assert.deepEqual(c.path, o.waypoints);
      for (let t = 0; t < 40; t++) c.step();
      assert.deepEqual(c.violations, [], `map ${map} ${o.id}`);
    }
  }
  assert.equal(total, Object.keys(stored).length);
  assert.deepEqual(SI.demo.pathOptions("office"), []);
  assert.match(SI.demo.pathOptions(14).find((o) => o.id === "fig_TRO-Fig15b").info, /digitised from T-RO Fig\. 15\(b\)/);
  assert.match(SI.demo.pathOptions(14)[0].info, /ProjectPanel5/);
  assert.throws(() => new SI.Controller({ map: 13, path: "P12" }), /not on map 13/);
  assert.throws(() => new SI.Controller({ map: 13, path: "P13+" }), /leaves polygon 13|not strictly inside/);
});

test("drawn paths: waypoint checks, back and forth on open paths, loops on closed ones", () => {
  const poly = SI.polygon.loadPolygon(12);
  const check = SI.demo.checkWaypoint;
  assert.equal(check(poly, null, [100, 450]), null);
  assert.equal(check(poly, null, [300, 800]), "outside"); /* inside the wall block of 12.dat */
  assert.equal(check(poly, null, [-5, 450]), "outside");
  assert.equal(check(poly, null, [0, 450]), "outside"); /* on the boundary */
  assert.equal(check(poly, [100, 450], [100, 450]), "repeat");
  assert.equal(check(poly, [560, 450], [150, 950]), "leaves");
  assert.equal(check(poly, [100, 450], [560, 450]), null);
  const comb = new SI.polygon.Polygon([[0, 0], [50, 0], [50, 30], [40, 30], [40, 10], [30, 10], [30, 30], [20, 30],
    [20, 10], [10, 10], [10, 30], [0, 30]], "comb");
  assert.equal(check(comb, [5, 10], [45, 10]), "along"); /* runs along the walls at y = 10 */
  assert.equal(check(comb, [5, 9.9], [45, 9.9]), null);

  const open = [[100, 450], [560, 450], [560, 800]];
  const c = new SI.Controller({ map: 12, path: open, seed: 2, mode: "unknown" });
  assert.equal(c.sim.atEnd, "reverse");
  assert.deepEqual(c.path, open);
  for (let t = 0; t < 2 * c.lapTicks + 5; t++) {
    c.step();
    checkTick(c, "drawn open");
  }
  assert.deepEqual(c.violations, []);

  const closed = [[100, 450], [560, 450], [560, 800], [560, 450], [100, 450]];
  const d = new SI.Controller({ map: 12, path: closed, seed: 2 });
  assert.equal(d.sim.atEnd, "loop");
  for (let t = 0; t < d.lapTicks + 5; t++) {
    d.step();
    checkTick(d, "drawn closed");
  }
  assert.throws(() => new SI.Controller({ map: 12, path: [[100, 450], [300, 800]] }), /not strictly inside/);
});

test("prefetch computes the reversed track once, only for back-and-forth paths", () => {
  const c = new SI.Controller({ map: 1, seed: 5 });
  assert.equal(c.prefetch(), true);
  assert.equal(c.prefetch(), false);
  c.reset();
  assert.equal(c.prefetch(), true);
  assert.equal(new SI.Controller({ map: 2 }).prefetch(), false); /* closed tour: loops */
  assert.equal(new SI.Controller({}).prefetch(), false); /* grid */
});

test("exact re-rooting keeps the bounds on a polygon map", () => {
  for (const mode of MODES) {
    const c = new SI.Controller({ map: 3, nTargets: 12, seed: 4, mode, compactEdges: 150, checkEvery: 3 });
    const sim = new SI.PolygonSimulator(3, null, c.sim.targets.length, 4, 4, { targetStep: 8 });
    const ref = new SI.CombinatorialFilter(sim.initialCondition(mode), c.total);
    for (let t = 0; t < c.lapTicks; t++) {
      c.step();
      ref.extend(sim.step().events);
      if (c.compactions.relaxed) continue;
      const rb = ref.allBounds();
      for (const sh of c.shadows) assert.deepEqual([sh.lo, sh.hi], rb[sh.label], `${mode} t=${c.t} s${sh.label}`);
    }
    assert.ok(c.compactions.exact > 0, `${mode}: no re-root happened`);
  }
});

test("the exactness check is spread over ticks and re-roots at its snapshot", () => {
  /* checkBudgetMs: 0 compares one subset per tick: the re-root lands ticks after its snapshot, and the
   * replayed filter keeps the bounds of the uncompacted one */
  for (const mode of MODES) {
    const c = new SI.Controller({ map: 3, nTargets: 12, seed: 4, mode, compactEdges: 150, checkEvery: 3, checkBudgetMs: 0 });
    const sim = new SI.PolygonSimulator(3, null, c.sim.targets.length, 4, 4, { targetStep: 8 });
    const ref = new SI.CombinatorialFilter(sim.initialCondition(mode), c.total);
    let spread = 0;
    for (let t = 0; t < c.lapTicks; t++) {
      const before = c.compactions.exact;
      c.step();
      ref.extend(sim.step().events);
      if (c.compactions.exact > before && c.compactions.lastT < c.t) spread += 1;
      if (c.compactions.relaxed) continue;
      const rb = ref.allBounds();
      for (const sh of c.shadows) assert.deepEqual([sh.lo, sh.hi], rb[sh.label], `${mode} t=${c.t} s${sh.label}`);
      assert.ok(c.edges <= c.opts.hardEdges);
    }
    assert.ok(c.compactions.exact > 0, `${mode}: no re-root happened`);
    /* (with one evader the re-roots are decided at once: a single hidden shadow, or all bounds fixed) */
    if (mode !== "evader") assert.ok(spread > 0, `${mode}: no re-root after a multi-tick check`);
    assert.deepEqual(c.violations, []);
  }
});

test("split probabilities follow shadow areas", () => {
  const labels = [0, 7, 7, 7, 8, 9, 9, 0];
  const p = SI.demo.splitProbabilities([{ type: "split", s: 1, a: 7, b: 6 }, { type: "split", s: 6, a: 8, b: 9 }], labels);
  assert.equal(p.get(1), 1 / 3);
  assert.equal(p.get(0), 3 / 6);
  const q = SI.demo.splitProbabilities([{ type: "split", s: 1, a: 7, b: 6 }, { type: "split", s: 6, a: 8, b: 9 }],
    new Map([[7, 3], [8, 1], [9, 2]]));
  assert.deepEqual([q.get(0), q.get(1)], [3 / 6, 1 / 3]);
  const c = new SI.Controller({ map: 12 });
  const areas = c.areas(), poly = c.sim.poly;
  let hidden = 0;
  areas.forEach((a) => { hidden += a; });
  const total = Math.abs(SI.polygonSim.polygonArea(poly.vertices));
  const visible = Math.abs(SI.polygonSim.polygonArea(c.sim.visibilityPolygon()));
  assert.ok(Math.abs(hidden + visible - total) < 1e-6 * total, "pockets and the visible region tile the map");
});

test("box-with-total bounds", () => {
  const b = { 1: [0, 2], 2: [1, Infinity], 3: [0, 1] };
  assert.deepEqual(SI.demo.boxTotalBounds([1, 3], [2], b, [2, 4]), [0, 3]);
  assert.deepEqual(SI.demo.boxTotalBounds([1, 2], [3], b, [0, Infinity]), [1, Infinity]);
  assert.deepEqual(SI.demo.boxTotalBounds([1, 3], [2], b, [5, 5]), [0, 3]);
});

test("demo scripts load as classic scripts in index.html order", () => {
  const html = fs.readFileSync(path.join(__dirname, "..", "docs", "index.html"), "utf8");
  const srcs = [...html.matchAll(/<script src="([^"]+)"><\/script>/g)].map((m) => m[1]);
  assert.ok(srcs.includes("js/controller.js") && srcs.includes("js/app.js"));
  assert.ok(!/type="module"/.test(html));
  const ctx = vm.createContext({ console });
  ctx.globalThis = ctx;
  for (const src of srcs.filter((s) => s !== "js/app.js")) {
    vm.runInContext(fs.readFileSync(path.join(__dirname, "..", "docs", src), "utf8"), ctx, { filename: src });
  }
  const c = new ctx.SI.Controller({ seed: 2 });
  for (let t = 0; t < 20; t++) c.step();
  assert.equal(c.violations.length, 0);
  const p = new ctx.SI.Controller({ map: "12", seed: 2 });
  for (let t = 0; t < 20; t++) p.step();
  assert.equal(p.world, "polygon");
  assert.equal(p.violations.length, 0);
});

test("bipartite view counts the disappeared vertices it leaves out and reports the re-root time", () => {
  const c = new SI.Controller({ nTargets: 20, seed: 4 });
  let hidden = 0, rooted = false;
  for (let t = 0; t < 2 * LAP; t++) {
    c.step();
    let gone = 0;
    c.filter.istate.disappeared.forEach((bd, v) => { if (!SI.bipartite.isPseudo(v)) gone += 1; });
    const view = c.bipartiteView(14, 5);
    const shown = view.right.filter((r) => r.kind === "gone").length;
    assert.ok(shown <= 5);
    assert.equal(shown + view.hiddenGone, gone, `t=${c.t}`);
    if (view.hiddenGone) hidden += 1;
    assert.equal(view.rootT, c.compactions.lastT);
    if (view.rootT !== null) rooted = true;
  }
  assert.ok(hidden > 0, "no tick hid a disappeared vertex");
  assert.ok(rooted, "no re-root happened");
});
