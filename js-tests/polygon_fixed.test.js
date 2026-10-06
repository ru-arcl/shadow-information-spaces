/* compat="fixed" behaviour of docs/js/polygon.js and polygon_sim.js that the Java golden fixtures cannot
 * show (ports of the regressions in tests/test_polygon_geometry.py, test_polygon_gaps.py and
 * test_polygon_simulate.py): polygon validation, the repaired physical gaps on lines through two
 * vertices, points outside the polygon, samples pulled back inside, never-see-evader with an empty
 * first set, and paths along an edge. */
"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");

require("../docs/js/rng.js");
require("../docs/js/maps_polygon.js");
require("../docs/js/polygon.js");
const SI = require("../docs/js/polygon_sim.js");
const P = SI.polygon;
const S = SI.polygonSim;

const L_ROOM = [[0, 0], [20, 0], [20, 10], [10, 10], [10, 20], [0, 20]];
const COMB = new P.Polygon([[0, 0], [50, 0], [50, 30], [40, 30], [40, 10], [30, 10], [30, 30], [20, 30], [20, 10], [10, 10],
  [10, 30], [0, 30]], "comb");
const gapsList = (gs) => gs.map((g) => g.toList());

/* Pocket areas plus the visibility polygon's area equal the polygon's, and every grid point is visible
 * (closed visibility) or in exactly one pocket. */
function assertTiles(poly, q) {
  const pockets = P.physicalGaps(poly, q).map((g) => S.pocketPolygon(poly, g));
  const area = pockets.reduce((a, pk) => a + S.polygonArea(pk), S.polygonArea(P.visibilityPolygon(poly, q)));
  assert.ok(Math.abs(area - poly.signedArea()) <= 1e-9 * poly.signedArea(), `${q}: ${area} vs ${poly.signedArea()}`);
  const xs = Array.from(poly.vx), ys = Array.from(poly.vy);
  const [x0, x1, y0, y1] = [Math.min(...xs), Math.max(...xs), Math.min(...ys), Math.max(...ys)];
  for (let i = 0; i < 50; i++) {
    for (let j = 0; j < 50; j++) {
      const p = [x0 + ((x1 - x0) * (i + 0.37)) / 50, y0 + ((y1 - y0) * (j + 0.61)) / 50];
      if (!S.pointInPolygon(poly, p)) continue;
      const cnt = pockets.filter((pk) => S.pointInPolygon(pk, p)).length;
      assert.equal(cnt, S.segmentInPolygon(poly, q, p) ? 0 : 1, `${q} sample ${p}`);
    }
  }
}

test("Polygon refuses clockwise, non-simple and non-finite vertices unless {validate: false}", () => {
  assert.throws(() => new P.Polygon(L_ROOM.slice().reverse()), /counter-clockwise/);
  assert.throws(() => new P.Polygon([[0, 0], [10, 0], [10, 10], [0, 10], [5, -5]]), /not simple: edges 0 and 3/);
  assert.throws(() => new P.Polygon([[0, 0], [10, 0], [5, 0], [5, 5]]), /not simple/); /* fold-back */
  assert.throws(() => new P.Polygon([[0, 0], [4, 0], [4, 4], [2, 0], [0, 4]]), /not simple/); /* vertex on an edge */
  assert.throws(() => new P.Polygon([[0, 0], [1, 0], [1, 0], [0, 1]]), /not simple/); /* zero-length edge */
  assert.throws(() => new P.Polygon([[0, 0], [1, 0], [NaN, 1]]), /finite/);
  new P.Polygon([[0, 0], [5, 0], [10, 0], [10, 10], [0, 10]]); /* a straight vertex is fine */
  assert.ok(new P.Polygon(L_ROOM.slice().reverse(), null, { validate: false }).signedArea() < 0);
  for (let n = 1; n <= 14; n++) assert.ok(P.loadPolygon(n).signedArea() > 0);
});

test("fixed physical gaps on a line through two vertices; points outside are refused", () => {
  const L = new P.Polygon(L_ROOM);
  /* (12, 8) is on the line through the reflex vertex (10, 10) and (0, 20): Java reports the one shadow
   * three times, once as its complement (the visible part) */
  assert.equal(P.physicalGaps(L, [12, 8], { java: true }).length, 3);
  assert.deepEqual(gapsList(P.physicalGaps(L, [12, 8])), [[3, 5, null, [0, 20]]]);
  assert.deepEqual(gapsList(P.physicalGaps(L, [12, 8.0001])), gapsList(P.physicalGaps(L, [12, 8.0001], { java: true })));
  for (const q of [[12, 8], [15, 5], [11, 9], [12, 8.0001], [12, 7.9999]]) assertTiles(L, q);
  /* map 14 at the P14a waypoint (633.93, 714.29): 19 Java pockets, 18 real ones */
  const P14 = P.loadPolygon(14);
  assert.equal(P.physicalGaps(P14, [633.93, 714.29], { java: true }).length, 19);
  assert.equal(P.physicalGaps(P14, [633.93, 714.29]).length, 18);
  assertTiles(P14, [633.93, 714.29]);
  for (const q of [[15, 15], [25, 5]]) {
    assert.deepEqual(P.physicalGaps(L, q, { java: true }), []);
    assert.throws(() => P.physicalGaps(L, q), P.GeometryError);
    assert.throws(() => P.visibilityPolygon(L, q), P.GeometryError);
    assert.equal(P.pointStatus(L, q), "outside");
  }
  assert.deepEqual(P.physicalGaps(L, [10, 15]), P.physicalGaps(L, [10, 15], { java: true })); /* on an edge */
  assert.equal(P.pointStatus(L, [10, 15]), "boundary");
});

test("exact orientation signs: crossSign equals Python's Fraction fallback (visibility._cross_sign)", () => {
  const cases = [ /* [a, b, c, d, sign((b - a) x (d - c))] as computed by Python */
    [[0, 0, 1, 1, 0, 0, 3, 3], 0],
    [[0, 0, 0.1, 0.2, 0, 0, 0.30000000000000004, 0.6], -1],
    [[0, 0, 0.1, 0.2, 0, 0, 0.3, 0.6000000000000001], 1],
    [[1e300, 0, 1e300 + 1e284, 1e-300, 1e300, 0, 1e300 + 2e284, 2e-300], 1],
    [[5e-324, 0, 1e-323, 5e-324, 0, 0, 1, 1], 0],
    [[0, 0, 0.1, 0.2, 0, 0, 0.3, 0.6], 0],
    [[0.1, 0.7, 0.3, 1.1, 0.2, 0.9, 0.5, 1.5], -1],
  ];
  for (const [a, s] of cases) assert.equal(P.crossSign(...a), s, String(a));
});

test("fixed tiling at every stored waypoint and on vertex lines of every map", () => {
  for (const label of Object.keys(SI.polygonMaps.paths)) {
    const rec = SI.polygonMaps.paths[label], poly = P.loadPolygon(rec.polygon);
    for (const w of rec.waypoints) if (P.pointStatus(poly, w) === "inside") assertTiles(poly, w);
  }
  let seed = 5;
  const rand = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
  for (let n = 1; n <= 14; n++) {
    const poly = P.loadPolygon(n), V = poly.vertices, refl = poly.reflex();
    let done = 0;
    while (done < 4) {
      const r = refl[Math.floor(rand() * refl.length)], v = Math.floor(rand() * poly.n);
      const t = [0.25, 0.5, 0.75, 1.5, 2.0, -0.5][Math.floor(rand() * 6)];
      const q = [V[r][0] + t * (V[r][0] - V[v][0]), V[r][1] + t * (V[r][1] - V[v][1])];
      if (v === r || P.pointStatus(poly, q) !== "inside" || poly.edges.some((e) => P.ptSegDist(e[0], e[1], e[2], e[3], q[0], q[1]) < 1e-3)) continue;
      done++;
      assertTiles(poly, q);
    }
  }
});

test("a sample past a path end near a wall is pulled back inside", () => {
  const L = new P.Polygon(L_ROOM);
  for (const x of [19.998, 19.9]) {
    const h = P.getGaps(L, new P.Path([[5, 15], [5, 2], [x, 2]]));
    assert.deepEqual(h.sets[h.sets.length - 1].map((g) => g.id), [2]);
    assert.deepEqual(P.gapHistoryToEvents(h).events, [{ type: "disappear", s: 1, lo: 0, hi: 0 }, { type: "appear", s: 2, lo: 0, hi: 0 }]);
    assert.ok(h.samples.points.every((p) => L.contains(p[0], p[1])));
    assert.equal(4 in h.samples.moved, x === 19.998);
  }
  const h = P.getGaps(L, new P.Path([[5, 15], [5, 2], [19.998, 2]]));
  assert.deepEqual(h.samples.moved[4], [20.003, 2]); /* Java's sample, outside the polygon */
});

test("never-see-evader: an empty first gap set makes every later gap clear", () => {
  const room = new P.Polygon([[0, 0], [10, 0], [10, 4], [4, 4], [4, 10], [0, 10]], "L");
  const h = P.getGaps(room, new P.Path([[2, 2], [9, 2]]));
  assert.deepEqual(P.historyLines(h, true), ["0.0 ", "0.2857142857142857 [1, 3] "]);
  h.sets[1].forEach((g) => { g.state = "contaminated"; }); /* a stale state is reset */
  assert.deepEqual(P.neverSeeEvader(h), [[], ["clear"]]);
  assert.equal(P.formatGapSet(h.sets[1], false), "[1, 3| 0] ");
  assert.deepEqual(P.neverSeeEvader([]), []);
});

test("validatePath refuses a segment along an edge; the simulator's ground truth near walls", () => {
  assert.throws(() => S.validatePath(COMB, [[5, 25], [5, 10], [45, 10], [45, 25]]), /runs along edge 4/);
  S.validatePath(COMB, [[5, 25], [5, 9.99], [45, 9.99], [45, 25]]);
  S.validatePath(new P.Polygon(L_ROOM), [[15, 5], [5, 15]]); /* through a reflex vertex: allowed */
  assert.equal(S.segmentAlongEdge(COMB, [5, 10], [45, 10]), 4);
  assert.equal(S.segmentAlongEdge(COMB, [5, 9.99], [45, 9.99]), null);
  const mismatches = (sim, ticks) => {
    let bad = 0;
    for (let t = 0; t < ticks; t++) {
      sim.step();
      sim.targets.forEach((p, k) => { if (S.segmentInPolygon(sim.poly, sim.robot, p) !== (sim.targetLabels[k] === 0)) bad++; });
    }
    return bad;
  };
  let sim = new SI.PolygonSimulator(COMB, [[5, 25], [5, 9.99], [45, 9.99], [45, 25]], 40, 5, 0.5, { targetStep: 0 });
  assert.equal(mismatches(sim, sim.ticksPerLap()), 0);
  assert.equal(sim.stats.membership_fallbacks, 0);
  const L = new P.Polygon(L_ROOM);
  const ge = P.gapHistoryToEvents(P.getGaps(L, new P.Path([[15, 5], [5, 15]])));
  assert.deepEqual(ge.initial, [1]);
  assert.deepEqual(ge.events, [{ type: "appear", s: 2, lo: 0, hi: 0 }, { type: "disappear", s: 1, lo: 0, hi: 0 }]);
  sim = new SI.PolygonSimulator(L, [[15, 5], [5, 15]], 40, 5, 0.5, { targetStep: 0 });
  assert.equal(mismatches(sim, sim.ticksPerLap()), 0);
  assert.equal(sim.stats.membership_fallbacks + sim.stats.sense_fallbacks, 0);
  sim = new SI.PolygonSimulator(12, [[300, 100], [400, 100]], 40, 3, 5);
  for (let t = 0; t < 25; t++) sim.step();
  sim.goTo([428.567, 100]);
  for (let t = 0; t < 15; t++) sim.step();
  assert.ok(sim.poly.contains(sim.sensePoint[0], sim.sensePoint[1]));
  assert.equal(sim.stats.membership_fallbacks, 0);
  const hidden = sim.targets.filter((p) => !S.segmentInPolygon(sim.poly, sim.robot, p)).length;
  assert.ok(hidden > 0);
  assert.equal(sim.targetLabels.filter((x) => x).length, hidden);
});
