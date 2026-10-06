/* Unit tests for docs/js/{rng,grid,map_office}.js (mirrors tests/test_grid.py). */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");

require("../docs/js/rng.js");
require("../docs/js/grid.js");
const SI = require("../docs/js/map_office.js");

const ROOT = path.join(__dirname, "..");
const MAPS = path.join(ROOT, "shadowinfo", "grid", "maps");
const office = () => SI.GridMap.fromText(SI.maps.office.text);
const lap = () => SI.expandPath(SI.maps.office.waypoints, SI.maps.office.loop);
const mask = (rows) => rows.join("").split("").map((ch) => (ch === "1" ? 1 : 0));

test("rng reference values", () => {
  const r = new SI.Rng(1);
  assert.deepEqual([r.nextU32(), r.next_u32(), r.nextU32()], [2693262067, 11749833, 2265367787]);
  const q = new SI.Rng(42);
  for (let i = 0; i < 1000; i++) {
    const x = q.random();
    assert.ok(x >= 0 && x < 1);
  }
  assert.throws(() => q.randint(0), RangeError);
});

test("map_office.js matches shadowinfo/grid/maps/", () => {
  assert.equal(SI.maps.office.text, fs.readFileSync(path.join(MAPS, "office.txt"), "utf8"));
  const spec = JSON.parse(fs.readFileSync(path.join(MAPS, "office_path.json"), "utf8"));
  assert.deepEqual(SI.maps.office.waypoints, spec.waypoints);
  const m = office();
  assert.equal(m.toText(), SI.maps.office.text);
  const p = lap();
  for (let k = 0; k < p.length; k++) {
    const a = p[k], b = p[(k + 1) % p.length];
    assert.ok(m.isFree(a[0], a[1]));
    assert.equal(Math.abs(a[0] - b[0]) + Math.abs(a[1] - b[1]), 1);
  }
});

test("map parser ignores trailing whitespace like Python", () => {
  for (const m of [SI.GridMap.fromText("####\n#..#  \n####\t\n"), new SI.GridMap(["####\r\n", "#..#\r\n", "####\r\n"])]) {
    assert.deepEqual([m.rows, m.cols, m.freeCells], [3, 4, [5, 6]]);
  }
  assert.deepEqual(SI.GridMap.fromText("####\n#.\n####\n").freeCells, [5]);
  assert.throws(() => SI.GridMap.fromText("####\n#. #\n####\n"));
});

test("sensing radius is a round disc (dr^2 + dc^2 <= R^2 + R)", () => {
  const n = 35, rows = ["#".repeat(n)];
  for (let k = 0; k < n - 2; k++) rows.push("#" + ".".repeat(n - 2) + "#");
  rows.push("#".repeat(n));
  const m = new SI.GridMap(rows), r0 = n >> 1, c0 = n >> 1;
  for (let R = 1; R <= 15; R++) {
    const v = m.visibility([r0, c0], R);
    for (let r = 1; r < n - 1; r++) {
      for (let c = 1; c < n - 1; c++) {
        const inside = (r - r0) * (r - r0) + (c - c0) * (c - c0) <= R * R + R;
        assert.equal(v[m.index(r, c)], inside ? 1 : 0, `R=${R} (${r},${c})`);
      }
    }
    assert.equal(v[m.index(r0 + R, c0 + 1)], 1);
  }
  const o = office();
  const w = o.visibility([13, 2], 12);
  assert.equal(w[o.index(25, 1)], 1, "no 1-cell wall nub shadow at (25, 1)");
});

test("corner rule and symmetry", () => {
  assert.equal(new SI.GridMap([".#", "#."]).lineOfSight(0, 0, 1, 1), false);
  assert.equal(new SI.GridMap([".#", ".."]).lineOfSight(0, 0, 1, 1), true);
  const m = office(), r = new SI.Rng(7), free = m.freeCells;
  for (let k = 0; k < 2000; k++) {
    const a = m.cell(free[r.randint(free.length)]), b = m.cell(free[r.randint(free.length)]);
    assert.equal(m.lineOfSight(a[0], a[1], b[0], b[1]), m.lineOfSight(b[0], b[1], a[0], a[1]));
  }
});

test("expandPath", () => {
  assert.deepEqual(SI.expandPath([[0, 0], [0, 2], [2, 2], [2, 0]]),
    [[0, 0], [0, 1], [0, 2], [1, 2], [2, 2], [2, 1], [2, 0], [1, 0]]);
  assert.deepEqual(SI.expandPath([[0, 0], [0, 2]], false), [[0, 0], [0, 1], [0, 2]]);
  assert.throws(() => SI.expandPath([[0, 0], [1, 1]]));
});

test("tracker: split, merge, disappear, appear, general group", () => {
  const t = new SI.ShadowTracker(1, 7);
  assert.deepEqual(t.reset(mask(["1111111"])), [1]);
  assert.deepEqual(t.update(mask(["111.111"])).events, [{ type: "split", s: 1, a: 2, b: 3 }]);
  assert.deepEqual(t.labels, [2, 2, 2, 0, 3, 3, 3]);
  assert.deepEqual(t.update(mask(["1111111"])).events, [{ type: "merge", a: 2, b: 3, s: 4 }]);
  assert.deepEqual(t.update(mask(["......."])).events, [{ type: "disappear", s: 4, lo: 0, hi: 0 }]);
  assert.deepEqual(t.update(mask(["1.....1"])).events.map((e) => [e.type, e.s]), [["appear", 5], ["appear", 6]]);
  const g = new SI.ShadowTracker(5, 5);
  g.reset(mask(Array(5).fill("1...1")));
  const tr = g.update(mask(["11111", ".....", ".....", ".....", "11111"]));
  assert.deepEqual(tr.events, [{ type: "merge", a: 1, b: 2, s: 3 }, { type: "split", s: 3, a: 4, b: 5 }]);
});

function replay(counts, frame) {
  const created = new Map(frame.created);
  for (const d of frame.events) {
    const s = d.s;
    if (d.type === "exit") {
      counts.set(s, counts.get(s) - d.k);
      assert.ok(counts.get(s) >= 0);
    } else if (d.type === "enter") counts.set(s, counts.get(s) + d.k);
    else if (d.type === "appear") {
      assert.equal(d.lo, d.hi);
      counts.set(s, d.lo);
    } else if (d.type === "disappear") {
      assert.equal(d.lo, counts.get(s));
      assert.equal(d.hi, counts.get(s));
      counts.delete(s);
    } else if (d.type === "split") {
      const n = counts.get(s);
      counts.delete(s);
      counts.set(d.a, created.get(d.a));
      counts.set(d.b, created.get(d.b));
      assert.equal(created.get(d.a) + created.get(d.b), n);
    } else if (d.type === "merge") {
      const n = counts.get(d.a) + counts.get(d.b);
      counts.delete(d.a);
      counts.delete(d.b);
      counts.set(s, n);
      assert.equal(created.get(s), n);
    }
  }
}

for (const [seed, radius] of [[1, null], [5, 6]]) {
  test(`simulator replay self-check seed=${seed} radius=${radius}`, () => {
    const p = lap();
    const sim = new SI.GridSimulator(office(), p, 30, seed, radius);
    const counts = new Map(sim.initialCounts);
    const kinds = {};
    for (let t = 0; t < p.length; t++) {
      const f = sim.step();
      f.events.forEach((e) => { kinds[e.type] = (kinds[e.type] || 0) + 1; });
      replay(counts, f);
      assert.deepEqual(Array.from(counts.entries()).sort((x, y) => x[0] - y[0]), f.counts);
    }
    assert.deepEqual(sim.robot, p[0]);
    if (radius === null) for (const k of ["appear", "disappear", "split", "merge", "enter", "exit"]) assert.ok(kinds[k] > 0, k);
    const seq = JSON.parse(JSON.stringify(sim.sequence("unknown")));
    assert.ok(Object.values(seq.initial).every((b) => b[0] === 0 && b[1] === null));
    assert.equal(seq.events.length, sim.history.length);
  });
}

test("SI.grid.simulator helper", () => {
  const sim = SI.grid.simulator("office", 4, 2);
  assert.equal(sim.targets.length, 4);
  assert.equal(sim.path.length, lap().length);
  assert.deepEqual(sim.initialCondition("evader")[1], [0, 1]);
});
