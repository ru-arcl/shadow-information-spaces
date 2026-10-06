/* docs/js/polygon.js against the original Java (golden fixtures tests/fixtures/java, written by
 * tools/java_reference): cuts, pointInPolygon, segInPolygon, physical gaps and visibility polygons,
 * critical points and their samples, and gap histories.  Exact equality of doubles throughout.
 *
 * Physical gaps / visibility polygons: {java: true} is Java's scan; the default (compat="fixed") equals it
 * except at the lattice points on a line through a reflex vertex and another vertex (DEGENERATE_LATTICE),
 * where the repaired pockets tile the polygon.
 *
 * Gap histories: by default the JS tracker must reproduce the Java history (IDs, iEdges, links in
 * HashSet order, times, printed lines) on every run where Java does not throw; with {javaMatching: false}
 * it re-derives exactly the steps where Java's samePhysicalGap match rotates the labels (quirk B3).  On
 * the runs where Java throws it repairs the crash, with the shadows alive at every sample equal to the
 * physical gaps there. */
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

const FIX = path.join(__dirname, "..", "tests", "fixtures", "java");
const load = (name) => JSON.parse(fs.readFileSync(path.join(FIX, `${name}.json`), "utf8"));
const RUNS = fs.readdirSync(FIX).filter((f) => /^\d+_.*\.json$/.test(f)).map((f) => f.slice(0, -5)).sort();
const bits = (xs) => Array.from(xs, (b) => (b ? "1" : "0")).join("");
const gapsList = (gs) => gs.map((g) => g.toList());

/* Lattice points of the golden visibility queries where Java's scan is degenerate
 * (tests/test_polygon_geometry.py DEGENERATE_LATTICE). */
const DEGENERATE_LATTICE = {
  11: [[100, 500], [400, 500], [500, 500], [600, 500]],
  12: [[500, 300], [500, 400], [500, 500]],
  13: [[400, 400], [500, 500], [900, 900]],
  14: [[100, 600], [900, 600], [800, 700]],
};

/* The fixed shadow pockets of q and its visibility polygon tile the polygon (areas add up), and a grid
 * of sample points is visible (closed visibility) or in exactly one pocket. */
function assertTiles(poly, q) {
  const S = SI.polygonSim;
  const gaps = P.physicalGaps(poly, q);
  const pockets = gaps.map((g) => S.pocketPolygon(poly, g));
  const area = pockets.reduce((a, pk) => a + S.polygonArea(pk), S.polygonArea(P.visibilityPolygon(poly, q)));
  assert.ok(Math.abs(area - poly.signedArea()) <= 1e-9 * poly.signedArea(), `${q}: ${area} vs ${poly.signedArea()}`);
  const xs = Array.from(poly.vx), ys = Array.from(poly.vy);
  const [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
  for (let i = 1; i < 40; i++) {
    for (let j = 1; j < 40; j++) {
      const p = [x0 + ((x1 - x0) * (i + 0.37)) / 40, y0 + ((y1 - y0) * (j + 0.61)) / 40];
      if (!S.pointInPolygon(poly, p)) continue;
      const cnt = pockets.filter((pk) => S.pointInPolygon(pk, p)).length;
      assert.equal(cnt, S.segmentInPolygon(poly, q, p) ? 0 : 1, `${q} sample ${p}`);
    }
  }
}

/* Steps where Java's label propagation rotates the IDs (tests/test_polygon_gaps.py ROTATED). */
const ROTATED = {
  "14_PP5": [124, 125, 126, 127, 128, 129, 130, 219, 220, 221, 222, 223],
  "14_fig_TRO-Fig15b": [124, 125, 126, 127, 128, 129, 221, 222, 223, 224, 225],
  "5_PP3": [123],
  "5_c_stub": [30],
};

test("primitives: Python float modulo, relativeCCW, Double.toString layout", () => {
  assert.equal(P.pyMod(-1, 10), 9);
  assert.equal(P.pyMod(7.5, 2), 1.5);
  assert.ok(Object.is(P.pyMod(-4, 2), 0));
  assert.equal(P.pyMod(-1e-20, 10), 10);
  assert.equal(P.relativeCCW(0, 0, 1, 0, 0.5, -1), 1); /* right of A->B in the y-up frame */
  assert.equal(P.relativeCCW(0, 0, 1, 0, 0.5, 1), -1);
  assert.equal(P.relativeCCW(0, 0, 1, 0, -1, 0), -1);
  assert.equal(P.relativeCCW(0, 0, 1, 0, 2, 0), 1);
  assert.equal(P.relativeCCW(0, 0, 1, 0, 0.5, 0), 0);
  const cases = [[0, "0.0"], [-0, "-0.0"], [100, "100.0"], [0.05623, "0.05623"], [1e-3, "0.001"], [1e-4, "1.0E-4"],
    [12345678, "1.2345678E7"], [9999999.5, "9999999.5"], [1e21, "1.0E21"], [1.5e-7, "1.5E-7"],
    [0.056234042553191504, "0.056234042553191504"], [-2.5, "-2.5"], [123.0, "123.0"]];
  for (const [x, s] of cases) assert.equal(P.javaDoubleStr(x), s, String(x));
});

test("java.util.HashSet<Integer> iteration order", () => {
  const cases = [
    [[15, 16], [16, 15]],
    [[31, 17, 1], [17, 1, 31]],
    [[65539, 3], [65539, 3]],
    [[-1, 5, -17], [-1, -17, 5]],
    [[20, 4, 36, 52, 68, 84, 100, 116, 132, 148, 164, 180, 196], [4, 36, 68, 100, 132, 164, 196, 20, 52, 84, 116, 148, 180]],
  ];
  for (const [ids, order] of cases) assert.deepEqual(P.javaHashsetOrder(ids), order);
  assert.deepEqual(P.javaHashsetOrder(Array.from({ length: 40 }, (_, i) => 40 - i)), Array.from({ length: 40 }, (_, i) => i + 1));
});

for (let n = 1; n <= 14; n++) {
  test(`polygon ${n}: vertices, cuts, pointInPolygon, visibility vs. Java`, () => {
    const d = load(`poly${n}`);
    const poly = P.loadPolygon(n);
    assert.deepEqual(poly.vertices, d.vertices);
    assert.deepEqual(Array.from(poly.turn), d.turn);
    assert.deepEqual(poly.reflex(), d.reflex);
    const c = d.cuts;
    const st = P.singleTangentCuts(poly);
    assert.deepEqual(st.map((x) => x.line), c.single_tangent);
    assert.deepEqual(P.inflections(poly, "NONGENERAL"), c.inflection_nongeneral);
    assert.deepEqual(P.inflections(poly, "ALL"), c.inflection_all);
    const gi = P.generalInflectionCuts(poly);
    assert.deepEqual(gi.map((g) => [g.line, g.fromLine, g.counterClockwise]), c.general_inflection_cut);
    assert.deepEqual(P.inflections(poly, "GENERAL"), gi.map((g) => g.line));
    const bt = P.bitangentCuts(poly);
    assert.deepEqual(bt.map((b) => [b.line, b.thisPoint, b.oppositePoint, b.oppositeSegment, b.curveToPoint]), c.bitangent_cut);
    const layout = [];
    for (const x of P.allCuts(poly)) {
      if (!layout.length || layout[layout.length - 1][0] !== x.type) layout.push([x.type, 0]);
      layout[layout.length - 1][1] += 1;
    }
    assert.deepEqual(layout, c.get_cuts_layout);

    const pip = d.point_in_polygon;
    const rows = [];
    for (let yi = 1; yi < 20; yi++) {
      const r = [];
      for (let xi = 1; xi < 20; xi++) r.push(poly.contains(50 * xi, 50 * yi));
      rows.push(bits(r));
    }
    assert.deepEqual(rows, pip.rows);
    assert.equal(bits(poly.vertices.map((v) => poly.contains(v[0], v[1]))), pip.at_vertices);
    assert.equal(bits(poly.edges.map((e) => poly.contains((e[0] + e[2]) / 2, (e[1] + e[3]) / 2))), pip.at_edge_midpoints);

    let nvp = 0;
    const degenerate = [];
    for (const v of d.visibility) {
      assert.equal(bits(P.vertexVisibility(poly, v.q)), v.seg_in_polygon_to_vertices, String(v.q));
      assert.equal(bits(poly.vertices.map((p) => P.segInPolygon(poly, v.q[0], v.q[1], p[0], p[1]))),
        v.seg_in_polygon_to_vertices);
      const java = P.physicalGaps(poly, v.q, { java: true });
      assert.deepEqual(gapsList(java), v.physical_gaps, String(v.q));
      const fixedGaps = P.physicalGaps(poly, v.q);
      const same = JSON.stringify(gapsList(fixedGaps)) === JSON.stringify(gapsList(java));
      if (!same) {
        degenerate.push(v.q);
        assertTiles(poly, v.q);
      }
      if (v.visibility_polygon) {
        nvp++;
        assert.deepEqual(P.visibilityPolygon(poly, v.q, { java: true }), v.visibility_polygon);
        const fixed = P.visibilityPolygon(poly, v.q);
        fixed.forEach((p, k) => assert.notDeepEqual(p, fixed[(k + 1) % fixed.length]));
        if (same) {
          const dedup = v.visibility_polygon.filter((p, k, a) => k === 0 || p[0] !== a[k - 1][0] || p[1] !== a[k - 1][1]);
          assert.deepEqual(dedup.slice(0, fixed.length), fixed, String(v.q));
        }
      }
    }
    assert.ok(nvp > 0);
    assert.deepEqual(degenerate, DEGENERATE_LATTICE[n] || []);
  });
}

/* Replay the events: the shadows alive at every sample are exactly the IDs of its physical gaps. */
function checkSamplesAgree(h) {
  const per = P.transitionEvents(h.sets);
  let alive = new Set(h.sets[0].map((g) => g.id));
  let k = 0;
  h.sampleIds.forEach((ids, i) => {
    while (k < h.sampleIndices.length && h.sampleIndices[k] < i) {
      for (const e of per[k]) {
        if (e.type === "appear") alive.add(e.s);
        else if (e.type === "disappear") alive.delete(e.s);
        else if (e.type === "split") {
          alive.delete(e.s);
          alive.add(e.a);
          alive.add(e.b);
        } else {
          alive.delete(e.a);
          alive.delete(e.b);
          alive.add(e.s);
        }
      }
      k++;
    }
    assert.equal(new Set(ids).size, ids.length, `sample ${i}`);
    assert.equal(ids.length, h.samples.physical[i].length);
    assert.deepEqual([...ids].sort((a, b) => a - b), [...alive].sort((a, b) => a - b), `sample ${i}`);
  });
  h.sets.forEach((gs, k2) => assert.deepEqual(gs.map((g) => g.id), h.sampleIds[h.sampleIndices[k2]]));
}

for (const name of RUNS) {
  test(`run ${name}: critical points, samples, visibility and gap history vs. Java`, () => {
    const r = load(name);
    const poly = P.loadPolygon(r.polygon);
    const pth = new P.Path(r.waypoints);
    assert.equal(pth.length(), r.path_length);
    const samples = P.samplePath(poly, pth);
    const cps = samples.criticalPoints;
    assert.equal(cps.length, r.n_critical_points);
    cps.forEach((cp, i) => {
      const row = r.critical_points[i];
      const javaPoint = i in samples.moved ? samples.moved[i] : samples.points[i];
      assert.deepEqual([cp.distance, cp.point, cp.cutType, cp.cutIndex, cp.segIndex, javaPoint], row.slice(0, 6), `cp ${i}`);
      if (row.length > 6) assert.deepEqual(gapsList(P.physicalGaps(poly, row[5], { java: true })), row[6], `gaps at cp ${i}`);
    });
    for (const v of r.visibility) {
      assert.deepEqual(gapsList(P.physicalGaps(poly, v.q, { java: true })), v.physical_gaps, v.where);
      assert.deepEqual(P.visibilityPolygon(poly, v.q, { java: true }), v.visibility_polygon, v.where);
    }

    const gg = r.get_gaps;
    if (gg.ok) {
      const h = P.getGaps(poly, pth, { samples });
      assert.deepEqual(h.inferred, []);
      assert.deepEqual(h.sets.map((gs) => gs.map((g) => g.toList())), gg.sets.map((s) => s.gaps));
      assert.deepEqual(h.times, gg.sets.map((s) => s.t));
      assert.deepEqual(h.sampleIndices, gg.sets.map((s) => s.sample_index));
      assert.deepEqual(P.historyLines(h), gg.sets.map((s) => s.line));
      assert.deepEqual(h.sets[h.sets.length - 1].map((g) => g.id), gg.final_ids);
      assert.equal(P.historyMaxId(h.sets), gg.max_id);
      const states = P.neverSeeEvader(h).map((ss) => ss.map((s) => (s === "clear" ? "0" : "1")).join(""));
      assert.deepEqual(states, r.nse.states);
      checkSamplesAgree(h);
    }
    const hf = P.getGaps(poly, pth, { samples, javaMatching: false });
    if (gg.ok) {
      assert.deepEqual(hf.inferred.map((x) => x[0]), ROTATED[name] || []);
      for (const [, why] of hf.inferred) assert.match(why, /^Java's step disagrees with the hidden chains; /);
      assert.equal(hf.sets.length, gg.n_sets);
      assert.equal(P.historyMaxId(hf.sets), gg.max_id);
      assert.deepEqual(hf.times, gg.sets.map((s) => s.t));
    } else {
      const crash = gg.failure.critical_index;
      assert.ok(hf.inferred.some(([i, why]) => i === crash && why.startsWith("Java throws")), JSON.stringify(hf.inferred));
    }
    checkSamplesAgree(hf);
    const ge = P.gapHistoryToEvents(hf);
    assert.equal(ge.events.length, ge.times.length);
    assert.deepEqual(ge.initial, hf.sets[0].map((g) => g.id));
  });
}

test("gap events of the paper runs (T-RO Fig. 15): P13 85 events, P14b 385 events and 491 IDs", () => {
  for (const [label, nev, nid] of [["P13", 85, 123], ["P14b", 385, 491]]) {
    const [poly, pth] = P.loadPath(label);
    const h = P.getGaps(poly, pth);
    const ge = P.gapHistoryToEvents(h);
    assert.equal(ge.events.length, nev, label);
    assert.equal(P.historyMaxId(h.sets), nid, label);
  }
});

test("fixed tracker repairs 14_PP4 physically: a merge and an appear in one step", () => {
  const r = load("14_PP4");
  const h = P.getGaps(P.loadPolygon(14), new P.Path(r.waypoints));
  assert.deepEqual(h.inferred, [[654, "Java throws; from the hidden chains: 476 -> [484], 470 -> [484], appear 483"]]);
  const k = h.sampleIndices.indexOf(653);
  assert.deepEqual(P.transitionEvents(h.sets)[k], [{ type: "merge", a: 476, b: 470, s: 484 }, { type: "appear", s: 483, lo: 0, hi: 0 }]);
});
