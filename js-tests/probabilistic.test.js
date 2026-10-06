/* Parity tests: docs/js/probabilistic.js vs shadowinfo.probabilistic (Python).
 * Fixture: tests/fixtures/js_filter_cases.json from tools/gen_js_filter_fixtures.py;
 * paper examples: tests/fixtures/tor_fig12_table3.json. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

require("../docs/js/rng.js");
require("../docs/js/events.js");
const SI = require("../docs/js/probabilistic.js");

const FIXTURES = path.join(__dirname, "..", "tests", "fixtures");
const load = (name) => JSON.parse(fs.readFileSync(path.join(FIXTURES, name), "utf8"));
const DATA = load("js_filter_cases.json");
const TOL = 1e-9;

function close(got, want, what) {
  assert.ok(Math.abs(got - want) <= TOL, `${what}: ${got} vs ${want}`);
}

function closeMap(got, want, what) {
  assert.deepEqual(Object.keys(got).sort(), Object.keys(want).sort(), what);
  for (const k of Object.keys(want)) close(got[k], want[k], `${what}[${k}]`);
}

function closeJoint(got, want, what) {
  assert.equal(got.length, want.length, `${what}: entry count`);
  got.forEach(([k, p], i) => {
    assert.deepEqual(k, want[i][0], `${what}: key ${i}`);
    close(p, want[i][1], `${what}: p${JSON.stringify(k)}`);
  });
}

function checkPmf(f, want, what) {
  assert.deepEqual(f.labels, want.labels, `${what}: labels`);
  assert.equal(f.entries, want.entries, `${what}: entries`);
  for (const s of f.labels) closeMap(f.marginal(s), want.marginals[s], `${what}: marginal ${s}`);
  if (want.joint) closeJoint(f.asList(), want.joint, what);
  close(f.total(), 1, `${what}: mass`);
}

function setup(c) {
  const seq = SI.prob.probSequenceFromDict(c.sequence);
  const opts = { splitRule: SI.prob.splitRuleFromDict(c.split_rule) };
  if (c.obs_model) opts.obsModel = c.obs_model;
  return { seq, opts };
}

test("fixture size", () => {
  const cases = DATA.probabilistic;
  assert.ok(cases.length >= 20);
  assert.ok(cases.some((c) => c.exact.error));
  assert.ok(cases.some((c) => (c.truncated || []).some((t) => t.error === "TruncationFailure")));
  const types = new Set(cases.flatMap((c) => c.sequence.events.map((e) => e.type + (e.dist ? "*" : "") + (e.p ? "p" : ""))));
  for (const t of ["fov", "split", "splitp", "merge", "appear", "appear*", "disappear", "disappear*", "enter", "exit"])
    assert.ok(types.has(t), t);
});

test("exact filter matches Python ExactFilter", () => {
  let ok = 0;
  for (const c of DATA.probabilistic) {
    const { seq, opts } = setup(c);
    const f = new SI.ExactFilter(seq.labels, seq.joint, opts);
    if (c.exact.error) {
      assert.throws(() => f.run(seq.observations), SI.InconsistentObservationError, c.name);
      continue;
    }
    f.run(seq.observations);
    checkPmf(f, c.exact, c.name);
    closeMap(f.expectedCounts(), c.exact.expected, `${c.name}: expected`);
    assert.equal(f.peakEntries, c.exact.peak_entries, c.name);
    ok++;
  }
  assert.ok(ok >= 20);
});

test("incremental apply equals run when no truncation", () => {
  for (const c of DATA.probabilistic.slice(0, 20)) {
    if (c.exact.error) continue;
    const { seq, opts } = setup(c);
    const f = new SI.ExactFilter(seq.labels, seq.joint, opts);
    for (const o of seq.observations) f.apply(o);
    checkPmf(f, c.exact, c.name);
  }
});

test("TR / RT / RT-LA truncation matches Python TruncatedFilter", () => {
  let n = 0;
  for (const c of DATA.probabilistic) {
    const { seq, opts } = setup(c);
    for (const t of c.truncated || []) {
      const what = `${c.name} ${t.mode}-${t.max_entries} seed ${t.seed}`;
      const f = new SI.TruncatedFilter(seq.labels, seq.joint, Object.assign({}, opts, {
        maxEntries: t.max_entries, mode: t.mode, seed: t.seed, lookaheadCountsFov: t.lookahead_counts_fov,
      }));
      if (t.error) {
        assert.throws(() => f.run(seq.observations), (err) => err.name === t.error, what);
        assert.equal(f.truncations, t.truncations, what);
        continue;
      }
      f.run(seq.observations);
      checkPmf(f, t, what);
      assert.equal(f.truncations, t.truncations, what);
      assert.equal(f.peakEntries, t.peak_entries, what);
      close(f.truncatedMass, t.truncated_mass, `${what}: truncated mass`);
      closeMap(f.expectedCounts(), t.expected, `${what}: expected`);
      n++;
    }
  }
  assert.ok(n >= 100);
});

test("Monte Carlo matches Python sample for sample", () => {
  for (const c of DATA.probabilistic) {
    if (!c.monte_carlo) continue;
    const { seq, opts } = setup(c);
    const m = c.monte_carlo;
    const run = () => SI.monteCarlo(seq, Object.assign({}, opts, { trials: m.trials, seed: m.seed, maxAttempts: 50 * m.trials }));
    if (m.error) {
      assert.throws(run, SI.InconsistentObservationError, c.name);
      continue;
    }
    const r = run();
    assert.deepEqual(r.labels, m.labels, c.name);
    assert.equal(r.successes, m.successes, c.name);
    assert.equal(r.attempts, m.attempts, c.name);
    closeJoint(r.asList(), m.joint, c.name);
  }
});

test("paper Table III: every step and P(s4) = 1/13, 2/13, 10/13", () => {
  const fx = load("tor_fig12_table3.json");
  const seq = SI.prob.probSequenceFromDict(fx.prob_sequence);
  const f = new SI.ExactFilter(seq.labels, seq.joint, {
    splitRule: SI.prob.splitRuleFromDict(fx.split_rule), obsModel: fx.obs_model,
  });
  const steps = fx.steps;
  seq.observations.forEach((o, t) => {
    f.apply(o);
    const want = steps[t + 1];
    assert.deepEqual(f.labels, want.labels, want.title);
    const got = new Map(f.asList().map(([k, p]) => [k.join(","), p]));
    assert.equal(got.size, want.joint_paper_order.length, want.title);
    for (const [k, p] of want.joint_paper_order) close(got.get(k.join(",")), p, `${want.title} ${k}`);
  });
  closeMap(f.marginal(4), fx.expected_final.marginal_s4, "final s4");
  close(f.expectedCounts()[4], 0 / 13 + 2 / 13 + 20 / 13, "E[s4]");
});

test("paper truncation example: keeping the 6 largest merge-output entries loses P(s4=0)", () => {
  const fx = load("tor_fig12_table3.json");
  const ex = fx.truncation_example;
  const seq = SI.prob.probSequenceFromDict(fx.prob_sequence);
  const f = new SI.ExactFilter(seq.labels, seq.joint, { obsModel: fx.obs_model });
  seq.observations.slice(0, 4).forEach((o) => f.apply(o));
  const ranked = f.asList().sort((x, y) => y[1] - x[1]);
  const dropped = ranked.slice(6).sort((x, y) => x[0][0] - y[0][0] || x[0][1] - y[0][1]);
  assert.deepEqual(dropped.map((kp) => kp[0]), ex.dropped.map((kp) => kp[0]));
  const g = new SI.ExactFilter(f.labels, ranked.slice(0, 6), { obsModel: fx.obs_model });
  g.apply(seq.observations[4]);
  const want = Object.fromEntries(Object.entries(ex.expected_final_marginal_s4).filter(([, p]) => p));
  closeMap(g.marginal(4), want, "s4 after truncation");
});

test("observation JSON round trip and models", () => {
  for (const c of DATA.probabilistic.slice(0, 10)) {
    const seq = SI.prob.probSequenceFromDict(c.sequence);
    assert.deepEqual(seq.observations.map(SI.prob.observationToDict), c.sequence.events, c.name);
  }
  const m = SI.prob.symmetricObsModel(0.9);
  assert.equal(m.enter.exit, 1 - 0.9);
  assert.deepEqual(SI.prob.perfectObsModel().null, { enter: 0, exit: 0, null: 1 });
  const rule = SI.prob.binomialSplit(0.5);
  assert.deepEqual(rule(2), [[0, 2, 0.25], [1, 1, 0.5], [2, 0, 0.25]]);
});

test("noise-free enter/exit and inconsistent disappear", () => {
  const f = SI.ExactFilter.fromCounts({ 1: 2 });
  f.apply({ type: "enter", s: 1, k: 3 });
  f.apply({ type: "exit", s: 1, k: 1 });
  assert.deepEqual(f.asList(), [[[4], 1]]);
  f.apply({ type: "split", s: 1, a: 2, b: 3, p: 1 });
  assert.deepEqual(f.asList(), [[[4, 0], 1]]);
  assert.throws(() => f.apply({ type: "disappear", s: 3, lo: 1, hi: 2 }), SI.InconsistentObservationError);
  const g = SI.ExactFilter.fromIndependent({ 1: { 0: 0.5, 1: 0.5 } }, { obsModel: SI.prob.symmetricObsModel(0.9) });
  g.apply({ type: "fov", s: 1, y: "exit" });
  closeMap(g.marginal(1), { 0: 0.45, 1: 0.5, 2: 0.05 }, "zero-exit rescaled within entry");
  const h = SI.ExactFilter.fromIndependent({ 1: { 0: 0.5, 1: 0.5 } });
  h.apply({ type: "fov", s: 1, y: "exit" });
  closeMap(h.marginal(1), { 0: 1 }, "no feasible branch: entry dropped");
});

test("fig. 16 (Table IV setting) exact run is fast enough for the demo", () => {
  const c = DATA.probabilistic.find((x) => x.name === "tor-fig16-table4");
  const { seq, opts } = setup(c);
  const t0 = performance.now();
  const f = new SI.ExactFilter(seq.labels, seq.joint, opts).run(seq.observations);
  const ms = performance.now() - t0;
  assert.equal(f.peakEntries, 28600);
  assert.ok(ms < 5000, `${ms.toFixed(0)} ms`);
});

/* ---- malformed models and labels (mirrors tests/test_probabilistic.py) ---------------- */

test("impossible FOV observation throws; unknown model keys and negative weights rejected", () => {
  const model = { enter: { enter: 0.9, exit: 0.1 }, exit: { enter: 0.1, exit: 0.9 } };
  const obs = [{ type: "fov", s: 1, y: "null" }, { type: "merge", a: 1, b: 2, s: 3 }];
  assert.throws(() => SI.ExactFilter.fromCounts({ 1: 2, 2: 1 }, { obsModel: model }).run(obs),
                SI.InconsistentObservationError);
  assert.throws(() => SI.monteCarlo({ labels: [1, 2], joint: [[[2, 1], 1]], observations: obs },
                                    { trials: 10, obsModel: model, maxAttempts: 100 }),
                SI.InconsistentObservationError);
  for (const bad of [{ exit: { Exit: 1 } }, { Exit: { exit: 1 } }, { exit: { exit: -0.5, enter: 1.5 } }]) {
    assert.throws(() => SI.ExactFilter.fromCounts({ 1: 2 }, { obsModel: bad }), RangeError, JSON.stringify(bad));
  }
});

test("unnormalised observation rows and custom split rules are normalised", () => {
  const prior = { 1: { 0: 0.5, 1: 0.5 } };
  const raw = SI.ExactFilter.fromIndependent(prior, { obsModel: { exit: { exit: 1.8, enter: 0.2 } } });
  raw.apply({ type: "fov", s: 1, y: "exit" });
  closeMap(raw.marginal(1), { 0: 0.45, 1: 0.5, 2: 0.05 }, "rescaled row");
  const uniform = (n) => Array.from({ length: n + 1 }, (_, a) => [a, n - a, 1]);
  const u = SI.ExactFilter.fromIndependent({ 1: { 0: 0.5, 2: 0.5 } }, { splitRule: uniform });
  u.run([{ type: "split", s: 1, a: 2, b: 3 }, { type: "merge", a: 2, b: 3, s: 4 }]);
  closeMap(u.marginal(4), { 0: 0.5, 2: 0.5 }, "uniform split rule");
  close(u.total(), 1, "total");
  assert.throws(() => SI.ExactFilter.fromCounts({ 1: 2 }, { splitRule: (n) => [[n, 1, 1]] })
    .apply({ type: "split", s: 1, a: 2, b: 3 }), RangeError);
});

test("binomial split of more than 1030 targets stays finite and normalised", () => {
  const n = 1100;
  for (const p of [0.5, 0.3]) {
    const d = SI.prob.binomialSplit(p)(n);
    const z = d.reduce((acc, r) => acc + r[2], 0);
    assert.ok(Math.abs(z - 1) < 1e-12, `p=${p}: sum ${z}`);
    assert.ok(d.every((r) => r[2] > 0 && Number.isFinite(r[2])));
  }
  /* log P(330 | 1100, 0.3) from the exact Fraction value (Python). */
  const d = new Map(SI.prob.binomialSplit(0.3)(n).map((r) => [r[0], r[2]]));
  const lc = (k) => { let s = 0; for (let j = 1; j <= k; j++) s += Math.log(j); return s; };
  const ref = Math.exp(lc(n) - lc(330) - lc(770) + 330 * Math.log(0.3) + 770 * Math.log(0.7));
  assert.ok(Math.abs(d.get(330) / ref - 1) < 1e-9);
  assert.deepEqual(SI.prob.binomialSplit(1)(n), [[n, 0, 1]]);
  assert.deepEqual(SI.prob.binomialSplit(0)(n), [[0, n, 1]]);
  const f = SI.ExactFilter.fromCounts({ 1: n });
  f.apply({ type: "split", s: 1, a: 2, b: 3 });
  const e = f.expectedCounts();
  assert.ok(Math.abs(e[2] - n / 2) < 1e-6 && Math.abs(e[3] - n / 2) < 1e-6);
});

test("split, merge and appear into an alive label throw", () => {
  for (const o of [{ type: "split", s: 1, a: 2, b: 3 }, { type: "split", s: 1, a: 4, b: 4 }, { type: "merge", a: 1, b: 3, s: 2 }]) {
    assert.throws(() => SI.ExactFilter.fromCounts({ 1: 2, 2: 1, 3: 0 }).apply(o), RangeError, JSON.stringify(o));
  }
  const f = SI.ExactFilter.fromCounts({ 1: 2, 2: 1 });
  f.apply({ type: "merge", a: 1, b: 2, s: 1 });
  assert.deepEqual([f.labels, f.asList()], [[1], [[[3], 1]]]);
  for (const o of [{ type: "split", s: 1, a: 2, b: 3 }, { type: "merge", a: 1, b: 3, s: 2 }, { type: "appear", s: 1, lo: 0, hi: 1 }]) {
    assert.throws(() => SI.monteCarlo({ labels: [1, 2, 3], joint: [[[2, 1, 0], 1]], observations: [o] }, { trials: 5 }),
                  RangeError, JSON.stringify(o));
  }
});

test("underflowed mass is not stored", () => {
  const tiny = 1e-200;
  const model = { enter: { enter: 1 - tiny, exit: tiny }, exit: { exit: 1 - tiny, enter: tiny }, null: { null: 1 } };
  const f = SI.ExactFilter.fromCounts({ 1: 3 }, { obsModel: model });
  for (let j = 0; j < 3; j++) {
    f.apply({ type: "fov", s: 1, y: "enter" });
    assert.ok(f.asList().every((kp) => kp[1] > 0));
  }
  assert.deepEqual(f.asList().map((kp) => kp[0]), [[4], [6]]);
  const g = SI.ExactFilter.fromIndependent({ 1: { 0: tiny, 3: 1 }, 2: { 0: tiny, 1: 1 } });
  g.apply({ type: "merge", a: 1, b: 2, s: 3 });
  g.apply({ type: "appear", s: 4, dist: { 0: tiny, 1: 1 } });
  assert.ok(g.asList().every((kp) => kp[1] > 0));
  close(g.total(), 1, "total");
});

test("duplicate joint keys are summed; truncation failure message; reordered()", () => {
  const seq = SI.prob.probSequenceFromDict({ labels: [1], joint: [[[1], 0.5], [[1], 0.25], [[2], 0.25]], events: [] });
  assert.deepEqual(seq.joint, [[[1], 0.75], [[2], 0.25]]);
  const obs = [{ type: "split", s: 1, a: 3, b: 4 }, { type: "disappear", s: 2, lo: 5, hi: 5 }];
  assert.throws(() => new SI.ExactFilter([1, 2], [[[3, 1], 1]]).run(obs),
                (err) => err instanceof SI.InconsistentObservationError && !(err instanceof SI.TruncationFailure));
  assert.throws(() => new SI.TruncatedFilter([1, 2], [[[3, 1], 1]], { maxEntries: 2, mode: "TR" }).run(obs),
                (err) => err instanceof SI.TruncationFailure && /no probability mass left after 1 truncation/.test(err.message));
  const s2 = { labels: [2, 1], joint: [[[1, 1], 1]], observations: [{ type: "split", s: 1, a: 5, b: 3 }] };
  const exact = new SI.ExactFilter(s2.labels, s2.joint).run(s2.observations);
  const mc = SI.monteCarlo(s2, { trials: 2000, seed: 3 });
  assert.deepEqual([exact.labels, mc.labels], [[2, 5, 3], [2, 3, 5]]);
  const got = mc.reordered(exact.labels);
  assert.deepEqual(got.labels, exact.labels);
  const want = new Map(exact.asList().map(([k, p]) => [k.join(), p]));
  assert.deepEqual(got.asList().map(([k]) => k.join()).sort(), [...want.keys()].sort());
  for (const [k, p] of got.asList()) assert.ok(Math.abs(p - want.get(k.join())) < 0.05);
  assert.throws(() => mc.reordered([2, 3]), RangeError);
});
