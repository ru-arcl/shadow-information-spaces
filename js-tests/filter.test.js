/* Parity tests: docs/js/{events,maxflow,bipartite,filter}.js vs shadowinfo (Python).
 * Fixture: tests/fixtures/js_filter_cases.json from tools/gen_js_filter_fixtures.py. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

require("../docs/js/events.js");
require("../docs/js/maxflow.js");
require("../docs/js/bipartite.js");
const SI = require("../docs/js/filter.js");

const FIXTURES = path.join(__dirname, "..", "tests", "fixtures");
const DATA = JSON.parse(fs.readFileSync(path.join(FIXTURES, "js_filter_cases.json"), "utf8"));

const hi = (x) => (x === null ? Infinity : x);
const bound = (b) => [b[0], hi(b[1])];
const boundsMap = (obj) => Object.fromEntries(Object.entries(obj).map(([k, b]) => [k, bound(b)]));

function build(c, fov) {
  const seq = SI.ShadowSequence.fromDict(c.sequence);
  return SI.CombinatorialFilter.fromSequence(seq, c.total, fov);
}

test("fixture size", () => {
  assert.ok(DATA.combinatorial.length >= 200);
  assert.ok(DATA.combinatorial.some((c) => !c.feasible));
  assert.ok(DATA.combinatorial.some((c) => c.total !== null));
  assert.ok(DATA.combinatorial.some((c) => c.sequence.events.some((e) => e.type === "enter" || e.type === "exit")));
  assert.ok(DATA.combinatorial.some((c) => c.all_bounds && Object.values(c.all_bounds).some((b) => b[1] === null)));
});

test("ShadowSequence JSON round trip and validation", () => {
  for (const c of DATA.combinatorial.slice(0, 50)) {
    const seq = SI.ShadowSequence.fromDict(c.sequence);
    assert.deepEqual(JSON.parse(JSON.stringify(seq)), c.sequence);
  }
  const ev = SI.ev;
  assert.throws(() => new SI.ShadowSequence({ 1: [0, 1] }, [ev.split(2, 3, 4)]), SI.InvalidSequenceError);
  assert.throws(() => new SI.ShadowSequence({ 1: [0, 1] }, [ev.split(1, 2, 3), ev.appear(2, 0, 0)]),
                SI.InvalidSequenceError);
  assert.throws(() => new SI.ShadowSequence({ 1: [0, 1] }, [ev.exit(1, 0)]), SI.InvalidSequenceError);
  const s = new SI.ShadowSequence({ 1: [0, null], 2: [1, 1] }, [ev.merge(1, 2, 3), ev.split(3, 4, 5)]);
  assert.deepEqual(s.aliveAtEnd(), [4, 5]);
  assert.deepEqual(s.allLabels(), [1, 2, 3, 4, 5]);
  assert.equal(s.initial[1][1], Infinity);
});

test("bounds, refinement, paper recipe and witnesses match Python", () => {
  let checked = 0;
  for (const c of DATA.combinatorial) {
    const f = build(c);
    assert.equal(f.feasible(), c.feasible, c.name);
    if (!c.feasible) {
      assert.throws(() => f.allBounds(), SI.InfeasibleError, c.name);
      continue;
    }
    assert.deepEqual(f.allBounds(), boundsMap(c.all_bounds), c.name);
    for (const q of c.set_bounds) assert.deepEqual(f.bounds(q.shadows), bound(q.bounds), c.name);
    assert.deepEqual(f.refineInitialBounds(), boundsMap(c.refine_initial), c.name);
    assert.deepEqual(f.initialTotalBounds(), bound(c.initial_total), c.name);
    if (c.paper) {
      assert.deepEqual(f.boundsPaper(c.paper.shadows), bound(c.paper.literal), c.name);
      assert.deepEqual(f.boundsPaper(c.paper.shadows, true), bound(c.paper.corrected), c.name);
    }
    for (const w of c.witness) {
      const got = f.witness(w.shadows, w.sense);
      assert.equal(got.value, w.value, c.name);
      const supply = Array.from(got.supply, ([v, n]) => [String(v), n]);
      assert.deepEqual(supply, w.supply, c.name);
      assert.deepEqual(got.flow.map(([u, tag, r, n]) => [String(u), tag, String(r), n]), w.flow, c.name);
    }
    assert.deepEqual(f.istate.toDict(), c.bipartite, c.name);
    checked++;
  }
  assert.ok(checked >= 200);
});

test("naive FOV conversion gives the same bounds as batching", () => {
  for (const c of DATA.combinatorial) {
    if (!c.feasible) continue;
    assert.deepEqual(build(c, "naive").allBounds(), boundsMap(c.all_bounds), c.name);
  }
});

test("incremental apply: bounds after every prefix", () => {
  let n = 0;
  for (const c of DATA.combinatorial) {
    if (!c.prefix_bounds) continue;
    const seq = SI.ShadowSequence.fromDict(c.sequence);
    const f = new SI.CombinatorialFilter(seq.initial, c.total);
    seq.events.forEach((e, i) => {
      f.apply(e);
      const want = c.prefix_bounds[i];
      if (want === null) assert.throws(() => f.allBounds(), SI.InfeasibleError);
      else assert.deepEqual(f.allBounds(), boundsMap(want), `${c.name} prefix ${i}`);
      n++;
    });
  }
  assert.ok(n > 100);
});

test("batchFov reproduces the four cases of Sec. V-D", () => {
  const ev = SI.ev;
  const out = SI.batchFov([ev.exit(1, 2), ev.enter(1, 3), ev.exit(1, 1), ev.split(1, 2, 3), ev.enter(2, 1)]);
  assert.deepEqual(out, [ev.exit(1, 2), ev.enter(1, 2), ev.split(1, 2, 3), ev.enter(2, 1)]);
  assert.deepEqual(SI.batchFov([ev.enter(1, 1), ev.exit(1, 1)]), []);
});

test("paper Fig. 11: s19 in [10, 24]", () => {
  const fx = JSON.parse(fs.readFileSync(path.join(FIXTURES, "tor_fig11.json"), "utf8"));
  const f = SI.CombinatorialFilter.fromSequence(SI.ShadowSequence.fromDict(fx.sequence));
  assert.deepEqual(f.bounds(fx.query), [fx.lp_lower, fx.lp_upper]);
  assert.deepEqual(f.boundsPaper(fx.query), [fx.paper_lower, fx.paper_upper]);
  assert.deepEqual(f.allBounds(), boundsMap(fx.lp_all_bounds));
  const edges = f.istate.edges().map(([u, , r]) => [u, r]);
  assert.deepEqual(edges.sort((a, b) => a[0] - b[0] || a[1] - b[1]), fx.bipartite.edges);
});

test("counting and pursuit-evasion helpers", () => {
  const ev = SI.ev;
  const pe = SI.pursuitEvasionFilter([1, 2]);
  pe.apply(ev.split(1, 3, 4));
  pe.apply(ev.disappear(3, 0, 0));
  pe.apply(ev.disappear(2, 0, 0));
  assert.deepEqual(SI.evaderStatus(pe), { 4: "evader" });
  const cf = SI.countingFilter([1]);
  cf.apply(ev.exit(1, 2));
  cf.apply(ev.disappear(1, 1, 1));
  assert.deepEqual(cf.initialTotalBounds(), [3, 3]);
  const pe2 = SI.pursuitEvasionFilter([1, 2]);
  pe2.apply(ev.enter(2, 1));
  assert.deepEqual(pe2.allBounds(), { 1: [0, 1], 2: [1, 2] });
  assert.deepEqual(SI.evaderStatus(pe2), { 1: "contaminated", 2: "evader" });
});

test("leftBounds / refineInitialBounds reject empty and unknown labels", () => {
  const ev = SI.ev;
  const f = new SI.CombinatorialFilter({ 1: [2, 2], 2: [0, 5] });
  f.extend([ev.appear(3, 1, 4), ev.disappear(2, 3, 3)]);
  for (const bad of [[], [999], [3, 999], [999, 998]]) assert.throws(() => f.leftBounds(bad), RangeError);
  assert.throws(() => f.refineInitialBounds([1, 999]), RangeError);
  assert.deepEqual(f.refineInitialBounds([]), {});
  assert.deepEqual(f.leftBounds([1, 2]), [5, 5]);
  // group networks are warm-started from the cached circulation
  assert.ok(f._feasibleNet(null, [1, 2]).net.isCirculation());
  assert.ok(f._feasibleNet([1, 3]).net.isCirculation());
  assert.deepEqual(f.bounds([1, 3]), [3, 6]);
});

test("maxflow: an unbounded maximize leaves the flow unchanged", () => {
  const net = new SI.FlowNetwork();
  const e = net.addEdge("x", "y", 0, Infinity);
  net.addEdge("y", "a", 0, 1);
  net.addEdge("a", "x", 0, 1);
  net.addEdge("y", "b");
  net.addEdge("b", "c");
  net.addEdge("c", "x");
  net.findFeasible();
  const before = net.flow.slice();
  assert.equal(net.maximize(e), Infinity);
  assert.deepEqual(net.flow, before);
  assert.ok(net.isCirculation());
});

test("maxflow: lower bounds, infeasibility and unbounded edges", () => {
  const net = new SI.FlowNetwork();
  net.addEdge("t", "s");
  const a = net.addEdge("s", "a", 2, 5);
  net.addEdge("a", "t", 0, 3);
  net.findFeasible();
  assert.ok(net.isCirculation());
  assert.deepEqual(net.extremes(a), [2, 3]);
  const bad = new SI.FlowNetwork();
  bad.addEdge("t", "s", 0, 1);
  bad.addEdge("s", "t", 2, 4);
  assert.throws(() => bad.findFeasible(), SI.InfeasibleError);
  const inf = new SI.FlowNetwork();
  inf.addEdge("t", "s");
  const e = inf.addEdge("s", "t", 1, Infinity);
  inf.findFeasible();
  assert.deepEqual(inf.extremes(e), [1, Infinity]);
  assert.equal(inf.numEdges, 2);
});

test("performance: 500 component events / ~100 alive shadows, allBounds < 200 ms", () => {
  const c = DATA.performance;
  const seq = SI.ShadowSequence.fromDict(c.sequence);
  build(c).allBounds();
  let best = Infinity;
  for (let rep = 0; rep < 3; rep++) {
    const t0 = performance.now();
    const f = SI.CombinatorialFilter.fromSequence(seq, c.total);
    const b = f.allBounds();
    best = Math.min(best, performance.now() - t0);
    assert.deepEqual(b, boundsMap(c.all_bounds));
  }
  assert.ok(c.component_events >= 500 && c.alive >= 100);
  assert.ok(best < 200, `allBounds took ${best.toFixed(1)} ms`);
});
