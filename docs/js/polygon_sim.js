/* Polygon world: a robot on a path, random-walk targets, shadow pockets, critical events and ground truth.
 *
 * Bit-for-bit port of shadowinfo/polygon/simulate.py (PolygonSimulator), on top of docs/js/polygon.js
 * (gap tracking = Python get_gaps(compat="fixed")) and SI.Rng (= shadowinfo.rng.Rng).  Same tick
 * semantics as the grid world (DESIGN §5): per sensor step, Exit events, then the component events
 * (Appear/Disappear carrying the covered/revealed counts), then Enter events; then the target substep.
 * Frames are the JSON-ready objects of the Python frame()/step() (events in the form of DESIGN §2).
 * The predicates use only + - * / and comparisons on doubles, in the Python operation order.
 * Browser load order: rng.js, maps_polygon.js, polygon.js, polygon_sim.js (Node requires them itself). */
(function (SI) {
  "use strict";

  if (typeof require === "function") {
    if (!SI.Rng) require("./rng.js");
    if (!SI.polygonMaps) require("./maps_polygon.js");
    if (!SI.polygon) require("./polygon.js");
  }
  var P = SI.polygon;
  var AT_END = ["stop", "loop", "reverse"];

  /* ---- predicates ---------------------------------------------------------- */

  /* (b - a) x (c - a): positive when c is left of a -> b (y-up). */
  var orient = P.orient;

  /* Even-odd crossing test on vertex arrays xs, ys. */
  function pipXY(xs, ys, px, py) {
    var inside = false, n = xs.length, j = n - 1;
    for (var i = 0; i < n; i++) {
      var xi = xs[i], yi = ys[i], xj = xs[j], yj = ys[j];
      if (yi > py !== yj > py && px < ((xj - xi) * (py - yi)) / (yj - yi) + xi) inside = !inside;
      j = i;
    }
    return inside;
  }

  /* point_in_polygon(poly or [[x, y], ...], p) */
  function pointInPolygon(polyOrPts, p) {
    if (polyOrPts instanceof P.Polygon) return pipXY(polyOrPts.vx, polyOrPts.vy, p[0], p[1]);
    return pipXY(polyOrPts.map(function (q) { return q[0]; }), polyOrPts.map(function (q) { return q[1]; }), p[0], p[1]);
  }

  /* The segment a b lies in the closed polygon (closed visibility): no proper edge crossing, and at
   * every vertex strictly inside the segment both directions point into the interior angle
   * (simulate.segment_in_polygon; shared with the gap repair of polygon.js). */
  var segmentInPolygon = P.closedSegmentInPolygon;

  /* Closed segment-edge intersection test against every edge (touching counts). */
  function segmentTouchesBoundary(poly, a, b) {
    var ax = a[0], ay = a[1], bx = b[0], by = b[1], n = poly.n, vx = poly.vx, vy = poly.vy;
    for (var k = 0; k < n; k++) {
      var k1 = k + 1 === n ? 0 : k + 1, x0 = vx[k], y0 = vy[k], x1 = vx[k1], y1 = vy[k1];
      var o1 = (bx - ax) * (y0 - ay) - (by - ay) * (x0 - ax);
      var o2 = (bx - ax) * (y1 - ay) - (by - ay) * (x1 - ax);
      var o3 = (x1 - x0) * (ay - y0) - (y1 - y0) * (ax - x0);
      var o4 = (x1 - x0) * (by - y0) - (y1 - y0) * (bx - x0);
      var s12 = (o1 <= 0 && o2 >= 0) || (o1 >= 0 && o2 <= 0);
      var s34 = (o3 <= 0 && o4 >= 0) || (o3 >= 0 && o4 <= 0);
      if (s12 && s34) return true;
    }
    return false;
  }

  function onBoundary(poly, p) {
    var px = p[0], py = p[1], n = poly.n, vx = poly.vx, vy = poly.vy;
    for (var k = 0; k < n; k++) {
      var k1 = k + 1 === n ? 0 : k + 1, x0 = vx[k], y0 = vy[k], x1 = vx[k1], y1 = vy[k1];
      var o = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0);
      if (o === 0 && (px - x0) * (x1 - x0) + (py - y0) * (y1 - y0) >= 0 && (px - x1) * (x0 - x1) + (py - y1) * (y0 - y1) >= 0) {
        return true;
      }
    }
    return false;
  }

  /* Signed shoelace area (positive counter-clockwise). */
  function polygonArea(pts) {
    var s = 0, n = pts.length;
    for (var i = 0; i < n; i++) {
      var a = pts[i], b = pts[(i + 1) % n];
      s += a[0] * b[1] - b[0] * a[1];
    }
    return 0.5 * s;
  }

  function fmtPt(p) {
    return "(" + p[0] + ", " + p[1] + ")";
  }

  /* At least 2 waypoints, each strictly inside, each segment in the closed polygon and not running along
   * an edge (whose critical lines are all collinear with it, so the events on it cannot be ordered and
   * every point of it is on the boundary); throws an Error naming the first violation otherwise. */
  function validatePath(poly, waypoints) {
    var pts = waypoints.map(function (p) { return [+p[0], +p[1]]; }), i;
    if (pts.length < 2) throw new Error("a path needs at least 2 waypoints");
    for (i = 0; i < pts.length; i++) {
      if (!isFinite(pts[i][0]) || !isFinite(pts[i][1])) throw new Error("waypoint " + i + " " + fmtPt(pts[i]) + " is not finite");
      if (!pointInPolygon(poly, pts[i]) || onBoundary(poly, pts[i])) {
        throw new Error("waypoint " + i + " " + fmtPt(pts[i]) + " is not strictly inside polygon " + poly.name);
      }
    }
    for (i = 0; i + 1 < pts.length; i++) {
      if (pts[i][0] === pts[i + 1][0] && pts[i][1] === pts[i + 1][1]) throw new Error("segment " + i + " has zero length");
      if (!segmentInPolygon(poly, pts[i], pts[i + 1])) {
        throw new Error("segment " + i + " " + fmtPt(pts[i]) + " -> " + fmtPt(pts[i + 1]) + " leaves polygon " + poly.name);
      }
      var k = segmentAlongEdge(poly, pts[i], pts[i + 1]);
      if (k !== null) {
        throw new Error("segment " + i + " " + fmtPt(pts[i]) + " -> " + fmtPt(pts[i + 1]) + " runs along edge " + k + " of polygon " + poly.name);
      }
    }
    return pts;
  }

  /* The first edge that the segment a b overlaps in more than a point (collinear), or null. */
  function segmentAlongEdge(poly, a, b) {
    var ax = a[0], ay = a[1], bx = b[0], by = b[1], dx = bx - ax, dy = by - ay, dd = dx * dx + dy * dy, n = poly.n;
    for (var k = 0; k < n; k++) {
      var p = poly.vertices[k], r = poly.vertices[(k + 1) % n];
      if (orient(ax, ay, bx, by, p[0], p[1]) !== 0 || orient(ax, ay, bx, by, r[0], r[1]) !== 0) continue;
      var t1 = (p[0] - ax) * dx + (p[1] - ay) * dy, t2 = (r[0] - ax) * dx + (r[1] - ay) * dy;
      if (Math.min(Math.max(t1, t2), dd) > Math.max(Math.min(t1, t2), 0.0)) return k;
    }
    return null;
  }

  /* ---- pockets -------------------------------------------------------------- */

  function cyclic(a, b, n) {
    var out = [a];
    while (out[out.length - 1] !== b) out.push((out[out.length - 1] + 1) % n);
    return out;
  }

  /* The hidden boundary chain of a gap closed by its window, counter-clockwise. */
  function pocketPolygon(poly, g) {
    var n = poly.n, V = poly.vertices;
    if (g.endPoint !== null) {
      return cyclic(g.startEdge, g.endEdge, n).map(function (k) { return V[k].slice(); }).concat([g.endPoint.slice()]);
    }
    return [g.startPoint.slice()].concat(cyclic((g.startEdge + 1) % n, (g.endEdge + 1) % n, n).map(function (k) { return V[k].slice(); }));
  }

  /* [anchor, free end]: the window from the occluding reflex vertex to where the gap edge meets the boundary. */
  function gapWindow(poly, g) {
    if (g.endPoint !== null) return [poly.vertices[g.startEdge].slice(), g.endPoint.slice()];
    return [poly.vertices[(g.endEdge + 1) % poly.n].slice(), g.startPoint.slice()];
  }

  /* The shadows seen from q: [{label, gap, pocket, window}] in physicalGaps order. */
  function shadowPockets(poly, q, gaps, labels) {
    gaps = gaps || P.physicalGaps(poly, q);
    return gaps.map(function (g, k) {
      return { label: labels ? labels[k] : k + 1, gap: g, pocket: pocketPolygon(poly, g), window: gapWindow(poly, g) };
    });
  }

  function shadowToDict(sh) {
    return { label: sh.label, window: sh.window.map(function (p) { return p.slice(); }),
      pocket: sh.pocket.map(function (p) { return p.slice(); }) };
  }

  /* ---- tracks, legs, configurations ------------------------------------------ */

  function eventLabels(e) {
    return e.type === "merge" ? [e.a, e.b, e.s] : e.type === "split" ? [e.s, e.a, e.b] : [e.s];
  }

  /* getGaps(compat="fixed", javaMatching: false) of one path as the simulator consumes it (cached per
   * polygon and path): the hidden-chain matching keeps every label on its physical shadow, which the
   * target ground truth needs (Java's label rotations, quirk B3, would move targets between labels). */
  function track(poly, path) {
    var key = "sim_track" + JSON.stringify(path.points), tr = poly._cache.get(key);
    if (tr) return tr;
    var h = P.getGaps(poly, path, { javaMatching: false }), per = P.transitionEvents(h.sets), eventsAt = new Map();
    var maxLocal = P.historyMaxId(h.sets);
    per.forEach(function (evs, k) {
      eventsAt.set(h.sampleIndices[k] + 1, evs);
      evs.forEach(function (e) { eventLabels(e).forEach(function (x) { if (x > maxLocal) maxLocal = x; }); });
    });
    var m = h.sampleIds[0].length;
    h.sampleIds[0].forEach(function (x, k) { if (x !== k + 1) throw new Error("unexpected initial labels"); });
    tr = { path: path, history: h, distances: h.samples.criticalPoints.map(function (cp) { return cp.distance; }),
      eventsAt: eventsAt, m: m, maxLocal: maxLocal };
    poly._cache.set(key, tr);
    return tr;
  }

  function relabel(e, g) {
    var s, a, b;
    if (e.type === "appear" || e.type === "disappear") return { type: e.type, s: g(e.s), lo: e.lo, hi: e.hi };
    if (e.type === "split") {
      s = g(e.s);
      a = g(e.a);
      b = g(e.b);
      return { type: "split", s: s, a: a, b: b };
    }
    if (e.type === "merge") {
      a = g(e.a);
      b = g(e.b);
      s = g(e.s);
      return { type: "merge", a: a, b: b, s: s };
    }
    throw new Error("not a component event: " + JSON.stringify(e));
  }

  function Leg(tr, initMap, base) {
    this.track = tr;
    this.initMap = initMap;
    this.base = base;
    this.nextCp = 1;
    this.s = 0;
  }

  Leg.prototype.g = function (local) {
    return local <= this.track.m ? this.initMap.get(local) : local - this.track.m - 1 + this.base;
  };

  Leg.prototype.events = function (i) {
    var self = this;
    return (this.track.eventsAt.get(i) || []).map(function (e) { return relabel(e, function (x) { return self.g(x); }); });
  };

  Leg.prototype.labels = function (i) {
    var self = this;
    return this.track.history.sampleIds[i].map(function (x) { return self.g(x); });
  };

  Object.defineProperty(Leg.prototype, "length", { get: function () { return this.track.path.length(); } });

  Leg.prototype.pointAt = function (s) {
    var path = this.track.path, pre = path.prefix;
    if (s >= pre[pre.length - 1]) return path.points[path.points.length - 1].slice();
    var i = 0;
    while (i + 1 < pre.length - 1 && pre[i + 1] <= s) i++;
    var seg = path.segments[i], f = (s - pre[i]) / (pre[i + 1] - pre[i]);
    return [seg[0] + f * (seg[2] - seg[0]), seg[1] + f * (seg[3] - seg[1])];
  };

  /* Where the sensor is and what it sees: the gaps with their (global) labels. */
  function Config(point, gaps, labels, criticalIndex) {
    this.point = point;
    this.gaps = gaps;
    this.labels = labels;
    this.criticalIndex = criticalIndex === undefined ? -1 : criticalIndex;
    this._shadows = null;
    this._arrays = null;
  }

  Config.prototype.shadows = function (poly) {
    if (!this._shadows) this._shadows = shadowPockets(poly, this.point, this.gaps, this.labels);
    return this._shadows;
  };

  Config.prototype.pocketArrays = function (poly) {
    if (!this._arrays) {
      this._arrays = this.shadows(poly).map(function (sh) {
        var xs = sh.pocket.map(function (p) { return p[0]; }), ys = sh.pocket.map(function (p) { return p[1]; });
        return { label: sh.label, xs: xs, ys: ys, x0: Math.min.apply(null, xs), x1: Math.max.apply(null, xs),
          y0: Math.min.apply(null, ys), y1: Math.max.apply(null, ys) };
      });
    }
    return this._arrays;
  };

  function sameStructure(a, b) {
    if (a.length !== b.length) return false;
    for (var k = 0; k < a.length; k++) {
      if (a[k].startEdge !== b[k].startEdge || a[k].endEdge !== b[k].endEdge || (a[k].startPoint === null) !== (b[k].startPoint === null)) {
        return false;
      }
    }
    return true;
  }

  /* perm[k] = index in a of the gap b[k] continues, if the hidden chains pair them 1-1. */
  function bijection(poly, a, b) {
    if (a.length !== b.length) return null;
    var comps = P.components(a.length, b.length, P.chainGraph(poly, a, b)), perm = b.map(function () { return -1; });
    for (var c = 0; c < comps.length; c++) {
      if (comps[c][0].length !== 1 || comps[c][1].length !== 1) return null;
      perm[comps[c][1][0]] = comps[c][0][0];
    }
    return perm;
  }

  function distToPolyline(xs, ys, px, py) {
    var n = xs.length, best = Infinity;
    for (var i = 0; i < n; i++) {
      var x1 = xs[i], y1 = ys[i], x2 = xs[(i + 1) % n], y2 = ys[(i + 1) % n];
      var dx = x2 - x1, dy = y2 - y1, L = dx * dx + dy * dy;
      var t = L > 0 ? ((px - x1) * dx + (py - y1) * dy) / L : 0;
      t = t < 0 ? 0 : t > 1 ? 1 : t;
      var cx = x1 + t * dx, cy = y1 + t * dy, ex = cx - px, ey = cy - py, d = ex * ex + ey * ey;
      if (d < best) best = d;
    }
    return best;
  }

  function inc(m, k) {
    m.set(k, (m.get(k) || 0) + 1);
  }

  function sortedPairs(m) {
    return Array.from(m.entries()).sort(function (p, q) { return p[0] - q[0]; });
  }

  function fovEvents(exits, enters) {
    var out = [];
    sortedPairs(exits).forEach(function (p) { out.push({ type: "exit", s: p[0], k: p[1] }); });
    sortedPairs(enters).forEach(function (p) { out.push({ type: "enter", s: p[0], k: p[1] }); });
    return out;
  }

  function initialBounds(counts, mode) {
    mode = mode || "exact";
    var out = {};
    counts.forEach(function (n, s) {
      if (mode === "exact") out[s] = [n, n];
      else if (mode === "unknown") out[s] = [0, Infinity];
      else if (mode === "evader") out[s] = [0, 1];
      else throw new Error("unknown mode " + JSON.stringify(mode) + "; expected one of exact,unknown,evader");
    });
    return out;
  }

  /* ---- the simulator --------------------------------------------------------- */

  /* new PolygonSimulator(poly, path, nTargets, seed, speed, opts)
   * poly: Polygon or map number 1..14; path: Path, waypoints, a stored path label (SI.polygonMaps.paths)
   * or null for the map's demo path; opts: {targetStep (8), atEnd ("auto"), validate (true)}. */
  function PolygonSimulator(poly, path, nTargets, seed, speed, opts) {
    opts = opts || {};
    if (!(poly instanceof P.Polygon)) poly = P.loadPolygon(+poly);
    this.poly = poly;
    var demoAtEnd = null;
    if (path === null || path === undefined) {
      if (poly.name === null || !/^\d+$/.test(poly.name)) throw new Error("path=null needs a numbered map (1..14)");
      var rec = SI.polygonMaps.demo[poly.name];
      if (!rec) throw new Error("no demo path for polygon " + poly.name);
      path = new P.Path(rec.waypoints);
      demoAtEnd = rec.at_end;
    } else if (typeof path === "string") {
      var label = path, lp = P.loadPath(label), p2 = lp[0];
      if (p2 !== poly && JSON.stringify(p2.vertices) !== JSON.stringify(poly.vertices)) {
        throw new Error("path " + JSON.stringify(label) + " belongs to polygon " + p2.name + ", not " + poly.name);
      }
      path = lp[1];
    }
    if (!(path instanceof P.Path)) path = new P.Path(path);
    if (opts.validate !== false) validatePath(poly, path.points);
    speed = speed === undefined ? 4 : +speed;
    if (!(speed > 0) || !isFinite(speed)) throw new Error("speed must be a positive number");
    var atEnd = opts.atEnd === undefined ? "auto" : opts.atEnd;
    if (atEnd === "auto") atEnd = demoAtEnd || (samePoint(path.points[0], path.points[path.points.length - 1]) ? "loop" : "stop");
    checkAtEnd(atEnd, path);
    this.path = path;
    this.atEnd = atEnd;
    this.speed = speed;
    this.targetStep = opts.targetStep === undefined ? 8 : +opts.targetStep;
    this.rng = new SI.Rng(seed === undefined ? 1 : seed);
    this.t = 0;
    this.stats = { sense_fallbacks: 0, membership_fallbacks: 0, crossings: 0, legs: 1 };
    this.history = [];
    this._pending = [];
    this._nextLabel = 1;
    this.legs = [];
    var tr = track(poly, path), init = new Map();
    for (var k = 1; k <= tr.m; k++) init.set(k, k);
    var leg = new Leg(tr, init, tr.m + 1);
    this._nextLabel = tr.maxLocal + 1;
    this.legs.push(leg);
    var smp = tr.history.samples;
    var first = new Config(smp.points[0], smp.physical[0].slice(), leg.labels(0), 0);
    this.robot = path.points[0].slice();
    this._cfg = this._configAt(this.robot, first);
    this.targets = this._placeTargets(nTargets === undefined ? 10 : nTargets | 0);
    this._states = this._targetStates(this._cfg, this.targets);
    this.initialCounts = this.counts();
    this.initialLabels = this._cfg.labels.slice();
  }

  function samePoint(a, b) {
    return a[0] === b[0] && a[1] === b[1];
  }

  function checkAtEnd(atEnd, path) {
    if (AT_END.indexOf(atEnd) < 0) throw new Error("atEnd must be one of stop,loop,reverse or 'auto', got " + JSON.stringify(atEnd));
    if (atEnd === "loop" && !samePoint(path.points[0], path.points[path.points.length - 1])) {
      throw new Error("atEnd='loop' needs a closed path (first waypoint == last waypoint)");
    }
  }

  Object.defineProperty(PolygonSimulator.prototype, "leg", { get: function () { return this.legs[this.legs.length - 1]; } });
  /* Labels of the shadows alive now, in boundary order. */
  Object.defineProperty(PolygonSimulator.prototype, "labels", { get: function () { return this._cfg.labels.slice(); } });
  /* Where the sensor is: the robot position, except for a sense-fallback tick. */
  Object.defineProperty(PolygonSimulator.prototype, "sensePoint", { get: function () { return this._cfg.point; } });
  /* Shadow of each target (0 = visible). */
  Object.defineProperty(PolygonSimulator.prototype, "targetLabels", { get: function () { return this._states.slice(); } });

  PolygonSimulator.prototype.shadows = function () {
    return this._cfg.shadows(this.poly).slice();
  };

  /* True number of targets in each alive shadow (Map label -> count, in label order of the boundary). */
  PolygonSimulator.prototype.counts = function () {
    var c = new Map();
    this._cfg.labels.forEach(function (s) { c.set(s, 0); });
    this._states.forEach(function (s) { if (s) c.set(s, c.get(s) + 1); });
    return c;
  };

  PolygonSimulator.prototype.initialCondition = function (mode) {
    return initialBounds(this.initialCounts, mode);
  };

  /* {initial, events}: all events so far with the initial condition for mode (exact/unknown/evader). */
  PolygonSimulator.prototype.sequence = function (mode) {
    return { initial: this.initialCondition(mode), events: this.history.slice() };
  };

  PolygonSimulator.prototype.visibilityPolygon = function () {
    return P.visibilityPolygon(this.poly, this.sensePoint);
  };

  /* JSON-ready state; geometry adds the shadow pockets and the visibility polygon. */
  PolygonSimulator.prototype.frame = function (geometry) {
    var f = {
      t: this.t,
      robot: this.robot.slice(),
      sense: this._cfg.point.slice(),
      leg: this.legs.length - 1,
      s: this.leg.s,
      targets: this.targets.map(function (p) { return p.slice(); }),
      target_labels: this._states.slice(),
      counts: sortedPairs(this.counts()),
      labels: this._cfg.labels.slice(),
    };
    if (geometry) {
      f.shadows = this.shadows().map(shadowToDict);
      f.visibility = this.visibilityPolygon();
    }
    return f;
  };

  /* Leave the current path: from the robot position, follow waypoints (validated). */
  PolygonSimulator.prototype.follow = function (waypoints, atEnd, validate) {
    var pts = [this.robot].concat(waypoints.map(function (p) { return [+p[0], +p[1]]; }));
    pts = pts.filter(function (p, i) { return i === 0 || !samePoint(p, pts[i - 1]); });
    if (pts.length < 2) return;
    if (validate !== false) validatePath(this.poly, pts);
    var path = new P.Path(pts);
    atEnd = atEnd || "stop";
    checkAtEnd(atEnd, path);
    this.path = path;
    this.atEnd = atEnd;
    this._startLeg(path, this._pending);
  };

  /* Click-to-move: drive straight to point (the segment must lie in the polygon) and stop. */
  PolygonSimulator.prototype.goTo = function (point) {
    this.follow([point], "stop");
  };

  PolygonSimulator.prototype._startLeg = function (path, steps) {
    var tr = track(this.poly, path), cur = steps.length ? steps[steps.length - 1][1] : this._cfg;
    var new0 = tr.history.samples.physical[0].slice(), point0 = tr.history.samples.points[0];
    var ids, evs = [];
    if (sameStructure(cur.gaps, new0)) {
      ids = cur.labels.slice();
    } else {
      var perm = bijection(this.poly, cur.gaps, new0);
      if (perm !== null) {
        ids = perm.map(function (j) { return cur.labels[j]; });
      } else {
        var prev = cur.gaps.map(function (g, k) { return new P.Gap(g, cur.labels[k]); });
        var fg = P.fromGraph(prev, new0.length, P.chainGraph(this.poly, cur.gaps, new0), this._nextLabel - 1);
        ids = fg.ids;
        this._nextLabel = fg.last + 1;
        evs = this._linkEvents(cur.labels, ids, fg.links);
      }
    }
    var base = this._nextLabel, init = new Map();
    for (var k = 0; k < tr.m; k++) init.set(k + 1, ids[k]);
    this._nextLabel = base + tr.maxLocal - tr.m;
    this.legs.push(new Leg(tr, init, base));
    this.stats.legs += 1;
    steps.push([evs, new Config(point0, new0, ids, 0)]);
  };

  PolygonSimulator.prototype._linkEvents = function (prevLabels, ids, links) {
    var prev = prevLabels.map(function (lab) { return new P.Gap(null, lab); });
    links.forEach(function (l) { l[1].forEach(function (c) { prev[l[0]].addGap(c); }); });
    var nxt = ids.map(function (x) { return new P.Gap(null, x); });
    var known = new Set(prevLabels.concat(ids)), fresh = new Map(), self = this;
    function g(x) {
      if (known.has(x)) return x;
      if (!fresh.has(x)) {
        fresh.set(x, self._nextLabel);
        self._nextLabel += 1;
      }
      return fresh.get(x);
    }
    return P.transitionEvents([prev, nxt])[0].map(function (e) { return relabel(e, g); });
  };

  PolygonSimulator.prototype._nextPath = function () {
    if (this.atEnd === "stop") return null;
    if (this.atEnd === "loop") return this.path;
    return new P.Path(this.leg.track.path.points.slice().reverse());
  };

  /* The gaps at q labelled as continuations of last (no critical point in between). */
  PolygonSimulator.prototype._configAt = function (q, last) {
    if (samePoint(q, last.point)) return last;
    var gaps = P.physicalGaps(this.poly, q);
    if (sameStructure(gaps, last.gaps)) return new Config(q, gaps, last.labels.slice());
    var perm = bijection(this.poly, last.gaps, gaps);
    if (perm !== null) return new Config(q, gaps, perm.map(function (j) { return last.labels[j]; }));
    this.stats.sense_fallbacks += 1;
    return last;
  };

  /* Shadow label of each point (0 = visible from cfg.point). */
  PolygonSimulator.prototype._targetStates = function (cfg, pts) {
    var poly = this.poly, out = [], arrays = null;
    for (var r = 0; r < pts.length; r++) {
      var p = pts[r];
      if (segmentInPolygon(poly, cfg.point, p)) {
        out.push(0);
        continue;
      }
      arrays = arrays || cfg.pocketArrays(poly);
      var found = 0, nfound = 0;
      for (var k = 0; k < arrays.length; k++) {
        var a = arrays[k];
        if (!(p[0] >= a.x0 && p[0] <= a.x1 && p[1] >= a.y0 && p[1] <= a.y1)) continue;
        if (pipXY(a.xs, a.ys, p[0], p[1])) {
          nfound++;
          if (found === 0) found = a.label;
        }
      }
      if (nfound !== 1) {
        this.stats.membership_fallbacks += 1;
        var best = Infinity, lab = 0;
        for (k = 0; k < arrays.length; k++) {
          var d = distToPolyline(arrays[k].xs, arrays[k].ys, p[0], p[1]);
          if (d < best) {
            best = d;
            lab = arrays[k].label;
          }
        }
        found = lab;
      }
      out.push(found);
    }
    return out;
  };

  PolygonSimulator.prototype._placeTargets = function (n) {
    var poly = this.poly;
    var x0 = Math.min.apply(null, poly.vx), x1 = Math.max.apply(null, poly.vx);
    var y0 = Math.min.apply(null, poly.vy), y1 = Math.max.apply(null, poly.vy);
    var out = [];
    while (out.length < n) {
      var px = x0 + this.rng.random() * (x1 - x0);
      var p = [px, y0 + this.rng.random() * (y1 - y0)];
      if (pointInPolygon(poly, p) && !onBoundary(poly, p)) out.push(p);
    }
    return out;
  };

  /* One random-walk step per target (two draws each, accepted or not); returns the moved indices. */
  PolygonSimulator.prototype._moveTargets = function () {
    var r = this.targetStep, cand = [], j;
    for (j = 0; j < this.targets.length; j++) {
      var dx = (2.0 * this.rng.random() - 1.0) * r;
      var dy = (2.0 * this.rng.random() - 1.0) * r;
      cand.push([this.targets[j][0] + dx, this.targets[j][1] + dy]);
    }
    var moved = [];
    for (j = 0; j < cand.length; j++) {
      if (pipXY(this.poly.vx, this.poly.vy, cand[j][0], cand[j][1]) && !segmentTouchesBoundary(this.poly, this.targets[j], cand[j])) {
        moved.push(j);
      }
    }
    for (j = 0; j < moved.length; j++) this.targets[moved[j]] = cand[moved[j]];
    return moved;
  };

  PolygonSimulator.prototype._advance = function (steps) {
    var remaining = this.speed;
    while (remaining > 0) {
      var leg = this.leg, L = leg.length;
      if (leg.s >= L) {
        var nxt = this._nextPath();
        if (nxt === null) break;
        this._startLeg(nxt, steps);
        continue;
      }
      var sNew = leg.s + remaining;
      if (sNew >= L) {
        remaining = sNew - L;
        sNew = L;
      } else {
        remaining = 0;
      }
      var dist = leg.track.distances, smp = leg.track.history.samples;
      while (leg.nextCp < dist.length && (dist[leg.nextCp] <= sNew || sNew >= L)) {
        var i = leg.nextCp;
        steps.push([leg.events(i), new Config(smp.points[i], smp.physical[i].slice(), leg.labels(i), i)]);
        this.stats.crossings += 1;
        leg.nextCp += 1;
      }
      leg.s = sNew;
    }
    this.robot = this.leg.pointAt(this.leg.s);
  };

  /* FOV + component events of one sensor step from the per-target labels before/after. */
  PolygonSimulator.prototype._sensorEvents = function (old, nw, evs) {
    var appeared = new Set(), disappeared = new Set(), children = new Map(), made = new Set();
    evs.forEach(function (e) {
      if (e.type === "appear") appeared.add(e.s);
      else if (e.type === "disappear") disappeared.add(e.s);
      else if (e.type === "split") {
        if (!children.has(e.s)) children.set(e.s, []);
        children.get(e.s).push(e.a, e.b);
        made.add(e.a);
        made.add(e.b);
      } else if (e.type === "merge") {
        if (!children.has(e.a)) children.set(e.a, []);
        children.get(e.a).push(e.s);
        if (!children.has(e.b)) children.set(e.b, []);
        children.get(e.b).push(e.s);
        made.add(e.s);
      }
    });
    var descCache = new Map();
    function desc(a) {
      if (!descCache.has(a)) {
        var out = new Set(), todo = [a];
        while (todo.length) {
          (children.get(todo.pop()) || []).forEach(function (c) {
            if (!out.has(c)) {
              out.add(c);
              todo.push(c);
            }
          });
        }
        descCache.set(a, out);
      }
      return descCache.get(a);
    }
    var exits = new Map(), enters = new Map(), revealed = new Map(), covered = new Map();
    var stayFrom = new Map(), stayTo = new Map();
    for (var j = 0; j < old.length; j++) {
      var a = old[j], b = nw[j];
      if (a && b && (a === b || desc(a).has(b))) {
        if (a !== b) {
          inc(stayFrom, a);
          inc(stayTo, b);
        }
        continue;
      }
      if (a) inc(disappeared.has(a) ? revealed : exits, a);
      if (b) inc(appeared.has(b) ? covered : enters, b);
    }
    var val = new Map();
    made.forEach(function (x) { val.set(x, stayTo.get(x) || 0); });
    for (var k = evs.length - 1; k >= 0; k--) { /* split chains: an intermediate label holds what its children receive */
      var e = evs[k];
      if (e.type === "split" && made.has(e.s)) val.set(e.s, val.get(e.a) + val.get(e.b));
    }
    var cnt = new Map(stayFrom), created = new Map(), comp = [];
    function pop(x) {
      var v = cnt.has(x) ? cnt.get(x) : 0;
      cnt.delete(x);
      return v;
    }
    evs.forEach(function (ev) {
      var x;
      if (ev.type === "appear") {
        x = covered.get(ev.s) || 0;
        ev = { type: "appear", s: ev.s, lo: x, hi: x };
      } else if (ev.type === "disappear") {
        x = revealed.get(ev.s) || 0;
        ev = { type: "disappear", s: ev.s, lo: x, hi: x };
      } else if (ev.type === "merge") {
        x = pop(ev.a);
        x += pop(ev.b);
        cnt.set(ev.s, x);
        created.set(ev.s, x);
      } else if (ev.type === "split") {
        var total = pop(ev.s);
        cnt.set(ev.a, val.get(ev.a));
        created.set(ev.a, val.get(ev.a));
        cnt.set(ev.b, val.get(ev.b));
        created.set(ev.b, val.get(ev.b));
        if (total !== val.get(ev.a) + val.get(ev.b)) throw new Error("split bookkeeping: " + JSON.stringify(ev));
      }
      comp.push(ev);
    });
    return [fovEvents(exits, new Map()).concat(comp, fovEvents(new Map(), enters)), created];
  };

  /* Advance one tick (sensor substep, then target substep) and return the frame. */
  PolygonSimulator.prototype.step = function () {
    this.t += 1;
    var steps = this._pending.slice(), self = this;
    this._pending = [];
    this._advance(steps);
    var last = steps.length ? steps[steps.length - 1][1] : this._cfg;
    var fin = this._configAt(this.robot, last);
    if (fin !== last || !steps.length) steps.push([[], fin]);
    var events = [], created = new Map(), states = this._states;
    steps.forEach(function (st) {
      var evs = st[0], cfg = st[1];
      if (cfg === self._cfg && !evs.length) return;
      var newStates = self._targetStates(cfg, self.targets);
      var res = self._sensorEvents(states, newStates, evs);
      Array.prototype.push.apply(events, res[0]);
      res[1].forEach(function (v, k) { created.set(k, v); });
      states = newStates;
    });
    this._cfg = fin;
    var moved = this._moveTargets(), newStates = states.slice();
    if (moved.length) {
      var ms = this._targetStates(fin, moved.map(function (j) { return self.targets[j]; }));
      moved.forEach(function (j, k) { newStates[j] = ms[k]; });
    }
    var exits = new Map(), enters = new Map();
    for (var j = 0; j < states.length; j++) {
      var a = states[j], b = newStates[j];
      if (a !== b) {
        if (a) inc(exits, a);
        if (b) inc(enters, b);
      }
    }
    Array.prototype.push.apply(events, fovEvents(exits, enters));
    this._states = newStates;
    Array.prototype.push.apply(this.history, events);
    var f = this.frame();
    f.events = events;
    f.created = sortedPairs(created);
    return f;
  };

  /* Frame 0 followed by ticks stepped frames. */
  PolygonSimulator.prototype.run = function (ticks) {
    var frames = [this.frame()];
    for (var t = 0; t < ticks; t++) frames.push(this.step());
    return frames;
  };

  /* Ticks needed to run the current path once. */
  PolygonSimulator.prototype.ticksPerLap = function () {
    return Math.ceil(this.path.length() / this.speed);
  };

  /* [sim, frames]: a simulator run for ticks ticks (default: one pass of the path). */
  function simulate(poly, path, nTargets, seed, ticks, speed, opts) {
    var sim = new PolygonSimulator(poly, path, nTargets, seed, speed, opts);
    return [sim, sim.run(ticks === undefined || ticks === null ? sim.ticksPerLap() : ticks)];
  }

  SI.PolygonSimulator = PolygonSimulator;
  SI.polygonSim = {
    AT_END: AT_END,
    orient: orient,
    pointInPolygon: pointInPolygon,
    segmentInPolygon: segmentInPolygon,
    segmentTouchesBoundary: segmentTouchesBoundary,
    onBoundary: onBoundary,
    polygonArea: polygonArea,
    validatePath: validatePath,
    segmentAlongEdge: segmentAlongEdge,
    pocketPolygon: pocketPolygon,
    gapWindow: gapWindow,
    shadowPockets: shadowPockets,
    shadowToDict: shadowToDict,
    initialBounds: initialBounds,
    track: track,
    simulate: simulate,
    PolygonSimulator: PolygonSimulator,
  };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
