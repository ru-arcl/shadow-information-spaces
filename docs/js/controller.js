/* DOM-free controller for the demo (DESIGN §6): a world simulation + combinatorial
 * filter (T-RO 2012, Sec. V) + truncated probabilistic filter (Sec. VI).
 *
 * Two worlds feed the same filters (events of DESIGN §2): the grid office map
 * (docs/js/grid.js, map "office") and the 14 polygons of the original Java
 * implementation (docs/js/polygon_sim.js, maps "1".."14"; a path is the map's demo
 * path, a stored path label of SI.polygonMaps.paths or a list of waypoints).
 *
 * The bipartite I-state grows with every event (each FOV batch adds a
 * pseudo-shadow, Sec. V-D), so a long run is kept small by re-rooting: the
 * history is replaced by a fresh filter whose shadows at t0 are the alive
 * shadows with their current bounds plus the bound on their sum.  The
 * feasible count vectors of the alive shadows are terminal demand vectors of a
 * bounded network, an integral g-polymatroid (Hoffman; DESIGN §6), so the
 * re-rooted filter is exact iff its
 * tight bounds agree with the old ones on every subset of alive shadows;
 * this is checked before an "exact" re-root.  If the graph keeps growing
 * without an exact re-root, a "relaxed" one is forced; it is an outer
 * approximation, so the bounds stay sound (they still contain the truth). */
(function (SI) {
  "use strict";

  if (typeof require === "function") {
    if (!SI.Rng) require("./rng.js");
    if (!SI.events) require("./events.js");
    if (!SI.CombinatorialFilter) require("./filter.js");
    if (!SI.TruncatedFilter) require("./probabilistic.js");
    if (!SI.GridSimulator) require("./grid.js");
    if (!SI.maps || !SI.maps.office) require("./map_office.js");
    if (!SI.PolygonSimulator) require("./polygon_sim.js");
  }

  var INF = Infinity;
  var DEFAULTS = {
    map: "office",
    path: null,
    stepLength: 4,
    targetStep: 8,
    nTargets: 10,
    seed: 1,
    mode: "exact",
    radius: null,
    prob: true,
    maxEntries: 1500,
    truncation: "RT-LA",
    probBudgetMs: 400,
    compactEdges: 1200,
    hardEdges: 5000,
    checkEvery: 8,
    maxCheckShadows: 9,
    checkBudgetMs: 6,
    forcedCheckMs: 40,
  };

  function byNum(a, b) {
    return a - b;
  }

  function addHi(a, b) {
    return a === INF || b === INF ? INF : a + b;
  }

  /* Number of grid cells per shadow label. */
  function shadowAreas(labels) {
    var area = new Map();
    for (var i = 0; i < labels.length; i++) if (labels[i]) area.set(labels[i], (area.get(labels[i]) || 0) + 1);
    return area;
  }

  /* Area of every pocket of a polygon-world simulator (Map label -> square units). */
  function pocketAreas(sim) {
    var area = new Map();
    sim.shadows().forEach(function (sh) { area.set(sh.label, Math.abs(SI.polygonSim.polygonArea(sh.pocket))); });
    return area;
  }

  /* Area-proportional split parameter p = |a| / (|a| + |b|) for every split of one tick,
   * from the final shadow areas (a Map label -> area, or a grid label array) by a backward
   * pass over the tick's events. */
  function splitProbabilities(events, areas) {
    var area = areas instanceof Map ? new Map(areas) : shadowAreas(areas), p = new Map();
    for (var i = events.length - 1; i >= 0; i--) {
      var e = events[i];
      if (e.type === "split") {
        var a = area.get(e.a) || 0, b = area.get(e.b) || 0;
        p.set(i, a + b > 0 ? a / (a + b) : 0.5);
        area.set(e.s, a + b);
      }
    }
    return p;
  }

  /* Edges of the bipartite I-state (Fig. 9) without building them. */
  function edgeCount(istate) {
    var n = 0;
    istate.disappearedReach.forEach(function (r) { n += r.size; });
    istate.reach.forEach(function (r) { n += r.size; });
    return n;
  }

  /* Tight bounds of the box {lo_s <= x_s <= hi_s} intersected with {tlo <= sum x <= thi} on subset A. */
  function boxTotalBounds(A, rest, b, tot) {
    var lo = 0, hi = 0, loRest = 0, hiRest = 0;
    A.forEach(function (s) {
      lo += b[s][0];
      hi = addHi(hi, b[s][1]);
    });
    rest.forEach(function (s) {
      loRest += b[s][0];
      hiRest = addHi(hiRest, b[s][1]);
    });
    return [Math.max(lo, hiRest === INF ? -INF : tot[0] - hiRest), Math.min(hi, tot[1] === INF ? INF : tot[1] - loRest)];
  }

  /* ---- maps and paths ------------------------------------------------------ */

  /* Paper figures drawn from the original maps (docs/notes/original_maps.md). */
  var FIGURES = { 12: "T-RO Fig. 11(a) / ICRA'08 Fig. 4", 13: "T-RO Fig. 15(a)", 14: "T-RO Fig. 15(b)" };

  function isPolygonMap(map) {
    return /^([1-9]|1[0-4])$/.test(String(map));
  }

  /* [{id, name, world, figure, title}]: the grid office map, then the 14 original maps. */
  function mapList() {
    var out = [{ id: "office", name: "Office (grid)", world: "grid", figure: null, title: "56 x 36 grid office floor plan" }];
    var polys = SI.polygonMaps.polygons;
    for (var k = 1; k <= 14; k++) {
      var rec = polys[String(k)];
      out.push({ id: String(k), name: "Map " + k + (FIGURES[k] ? " \u2014 " + FIGURES[k] : ""), world: "polygon",
                 figure: FIGURES[k] || null, title: rec.file + ", " + rec.n_vertices + " vertices" +
                   (rec.figure ? "; " + rec.figure : "") });
    }
    return out;
  }

  function sentence(s) {
    s = String(s).trim();
    return /[.!?)]$/.test(s) ? s : s + ".";
  }

  /* One-paragraph provenance of a stored path record or the demo record of a map. */
  function pathInfo(rec, label) {
    var parts = [];
    if (rec.provenance === "demo") {
      parts.push(rec.source ? "Default path: " + rec.note : "Default path: a closed patrol tour designed for this demo " +
        "(hand-picked via points joined by shortest routes keeping clear of the walls).");
    } else if (rec.kind === "code") {
      parts.push(label + ": waypoints copied from the Java code (" + rec.source.join(", ") + ", " + rec.status + ").");
      if (rec.figure) parts.push(sentence(rec.figure));
      if (rec.notes) parts.push(sentence(rec.notes));
    } else {
      parts.push(label + ": digitised from " + rec.figure + " (" + (rec.provenance.source_image || "paper figure") + ").");
      if (rec.paper_reported) parts.push("Paper: " + sentence(rec.paper_reported));
    }
    return parts.join(" ");
  }

  /* Path choices of a polygon map: [{id, name, kind, valid, error, info, waypoints}].
   * "default" is the map's demo path; the others are every stored path on the map (code
   * paths of the Java ProjectPanel classes and paths digitised from the paper figures). */
  function pathOptions(map) {
    if (!isPolygonMap(map)) return [];
    var maps = SI.polygonMaps, demo = maps.demo[String(map)], poly = SI.polygon.loadPolygon(+map);
    var out = [{ id: "default", name: (demo.source ? demo.source + " (default)" : "Patrol tour (default)"),
                 kind: "demo", valid: true, error: null, info: pathInfo(demo), waypoints: demo.waypoints }];
    Object.keys(maps.paths).forEach(function (label) {
      var rec = maps.paths[label];
      if (rec.polygon !== +map || label === demo.source) return;
      var error = null;
      try {
        SI.polygonSim.validatePath(poly, rec.waypoints);
      } catch (err) {
        error = err.message;
      }
      var kind = rec.kind === "code" ? "Java code" : "digitised, " + rec.figure;
      out.push({ id: label, name: label + " (" + kind + ")", kind: rec.kind, valid: error === null, error: error,
                 info: pathInfo(rec, label), waypoints: rec.waypoints });
    });
    return out;
  }

  /* Why a drawn waypoint b cannot follow a (null when it can): b strictly inside the
   * polygon, the segment a b in the closed polygon and not along an edge
   * (SI.polygonSim.validatePath). */
  function checkWaypoint(poly, a, b) {
    var S = SI.polygonSim;
    if (!S.pointInPolygon(poly, b) || S.onBoundary(poly, b)) return "outside";
    if (a && a[0] === b[0] && a[1] === b[1]) return "repeat";
    if (a && !S.segmentInPolygon(poly, a, b)) return "leaves";
    if (a && S.segmentAlongEdge(poly, a, b) !== null) return "along";
    return null;
  }

  /* ---- the controller ------------------------------------------------------- */

  function Controller(opts) {
    this.opts = Object.assign({}, DEFAULTS, opts || {});
    if (this.opts.mode === "evader") this.opts.nTargets = 1;
    if (SI.grid.MODES.indexOf(this.opts.mode) < 0) throw new RangeError("unknown mode " + this.opts.mode);
    this.opts.map = String(this.opts.map);
    this.world = isPolygonMap(this.opts.map) ? "polygon" : "grid";
    if (this.world === "grid" && !SI.maps[this.opts.map]) throw new RangeError("unknown map " + this.opts.map);
    if (this.world === "polygon") this.opts.radius = null;
    this.reset();
  }

  Controller.DEFAULTS = DEFAULTS;

  /* The polygon simulator for opts.path: null/"default" (demo path), a stored label or waypoints. */
  Controller.prototype._polygonSim = function () {
    var o = this.opts, poly = SI.polygon.loadPolygon(+o.map), path = o.path, atEnd;
    if (path === null || path === undefined || path === "default") {
      path = null;
    } else if (typeof path === "string") {
      var rec = SI.polygonMaps.paths[path];
      if (!rec || rec.polygon !== +o.map) throw new RangeError("path " + path + " is not on map " + o.map);
      path = rec.waypoints;
    }
    if (path !== null) {
      var first = path[0], last = path[path.length - 1];
      atEnd = first[0] === last[0] && first[1] === last[1] ? "loop" : "reverse";
    }
    return new SI.PolygonSimulator(poly, path, o.nTargets, o.seed, o.stepLength,
                                   { targetStep: o.targetStep, atEnd: atEnd });
  };

  Controller.prototype.reset = function () {
    var o = this.opts;
    if (this.world === "polygon") {
      this.sim = this._polygonSim();
      this.map = this.sim.poly;
      this.path = this.sim.path.points;
      this.lapTicks = this.sim.ticksPerLap();
      this._prefetched = false;
    } else {
      this.sim = SI.grid.simulator(o.map, o.nTargets, o.seed, o.radius);
      this.map = this.sim.map;
      this.path = this.sim.path;
      this.lapTicks = this.path.length;
    }
    this.t = 0;
    this.lastEvents = [];
    this.violations = [];
    this.errors = [];
    this.compactions = { exact: 0, relaxed: 0, lastT: null, lastKind: null };
    this._lastCheck = -INF;
    this._check = null;
    var counts = this.sim.initialCounts;
    var hidden0 = 0;
    counts.forEach(function (n) { hidden0 += n; });
    this.total = o.mode === "evader" ? [hidden0, hidden0] : null;
    this._setRoot(this.sim.initialCondition(o.mode), this.total, []);
    this.prob = null;
    this.probNote = null;
    this.probMs = 0;
    if (!o.prob) {
      this.probNote = "probabilistic filter off";
    } else if (o.mode === "unknown") {
      this.probNote = "needs a prior; not defined for unknown counts";
    } else {
      this._initProb(counts, hidden0);
    }
    this._update();
  };

  Controller.prototype._initProb = function (counts, hidden0) {
    var o = this.opts;
    var labels = Array.from(counts.keys()).sort(byNum);
    var joint;
    if (o.mode === "exact" || hidden0 === 0) {
      joint = [[labels.map(function (s) { return counts.get(s); }), 1]];
    } else {
      var area = this.areas(), z = 0;
      labels.forEach(function (s) { z += area.get(s) || 0; });
      joint = labels.map(function (s, i) {
        return [labels.map(function (_, j) { return j === i ? 1 : 0; }), z > 0 ? (area.get(s) || 0) / z : 1 / labels.length];
      });
    }
    var popts = { maxEntries: o.maxEntries, mode: o.truncation, seed: o.seed };
    this.prob = labels.length ? new SI.TruncatedFilter(labels, joint, popts) : new SI.TruncatedFilter([], [[[], 1]], popts);
  };

  /* Area of every alive shadow: grid cells, or square units of the pocket polygons. */
  Controller.prototype.areas = function () {
    return this.world === "polygon" ? pocketAreas(this.sim) : shadowAreas(this.sim.labels);
  };

  /* Compute the gap track of the next leg (the reversed path) ahead of time, so the robot does
   * not stall when it turns around; returns true when it did any work. */
  Controller.prototype.prefetch = function () {
    if (this.world !== "polygon" || this.sim.atEnd !== "reverse" || this._prefetched) return false;
    this._prefetched = true;
    var pts = this.sim.path.points.slice().reverse();
    SI.polygonSim.track(this.sim.poly, new SI.polygon.Path(pts));
    return true;
  };

  /* Advance one tick; returns the frame of the simulator's step() (events included). */
  Controller.prototype.step = function () {
    var frame = this.sim.step();
    this.t = frame.t;
    this.lastEvents = frame.events;
    this.lastFrame = frame;
    try {
      this.filter.extend(frame.events);
      for (var i = 0; i < frame.events.length; i++) this._log.push(frame.events[i]);
    } catch (err) {
      this.errors.push("t=" + this.t + ": " + err.message);
      throw err;
    }
    if (this.prob) this._stepProb(frame.events);
    this._update();
    this._maybeCompact();
    return frame;
  };

  Controller.prototype._stepProb = function (events) {
    var ps = splitProbabilities(events, this.areas());
    var t0 = Date.now();
    try {
      for (var i = 0; i < events.length; i++) {
        var e = events[i];
        if (e.type === "split") e = Object.assign({}, e, { p: ps.get(i) });
        this.prob.apply(e, events.slice(i + 1));
      }
    } catch (err) {
      if (!(err instanceof SI.InconsistentObservationError)) throw err;
      this.prob = null;
      this.probNote = (err instanceof SI.TruncationFailure ? "lost after truncation" : "inconsistent") +
        " at t=" + this.t + "; reset to restart";
      return;
    }
    this.probMs = Date.now() - t0;
    if (this.probMs > this.opts.probBudgetMs) {
      this.prob = null;
      this.probNote = "too expensive (" + this.probMs + " ms/tick); disabled";
    }
  };

  /* Recompute bounds, truth and the containment check. */
  Controller.prototype._update = function () {
    var filt = this.filter;
    var alive = filt.alive().sort(byNum);
    var bounds = filt.allBounds();
    var counts = this.sim.counts();
    var expected = this.prob ? this.prob.expectedCounts() : null;
    var area = this.areas();
    var self = this, hiddenTrue = 0;
    this.shadows = alive.map(function (s) {
      var truth = counts.get(s) || 0;
      hiddenTrue += truth;
      var b = bounds[s];
      var ok = b[0] <= truth && truth <= b[1];
      if (!ok) self.violations.push({ t: self.t, label: s, truth: truth, lo: b[0], hi: b[1] });
      return { label: s, truth: truth, lo: b[0], hi: b[1], ok: ok, area: area.get(s) || 0,
               expected: expected && s in expected ? expected[s] : null };
    });
    var hb = alive.length ? filt.bounds(alive) : [0, 0];
    var hok = hb[0] <= hiddenTrue && hiddenTrue <= hb[1];
    if (!hok) this.violations.push({ t: this.t, label: "all", truth: hiddenTrue, lo: hb[0], hi: hb[1] });
    this.hidden = { lo: hb[0], hi: hb[1], truth: hiddenTrue, ok: hok, total: this.sim.targets.length };
    this.edges = edgeCount(filt.istate);
  };

  /* A new filter rooted at init (bounds per label) and total, extended by events (the log of
   * everything applied since the root, which a snapshot for the exactness check replays). */
  Controller.prototype._setRoot = function (init, total, events) {
    this._rootInit = init;
    this._rootTotal = total;
    this._log = events.slice();
    this.filter = new SI.CombinatorialFilter(init, total);
    this.filter.extend(events);
  };

  /* Re-rooting is checked without blocking the page: once the I-state exceeds compactEdges (at most
   * every checkEvery ticks), _startCheck snapshots the filter and the alive shadows' bounds, and each
   * tick _advanceCheck compares at most checkBudgetMs worth of subset bounds on the snapshot.  When
   * every subset agrees, the filter re-rooted at the snapshot is extended by the events since then
   * (the same feasible set, so the same bounds as without re-rooting).  Above hardEdges a re-root is
   * forced at once (relaxed, unless a check within forcedCheckMs proves it exact). */
  Controller.prototype._maybeCompact = function () {
    var o = this.opts;
    if (this.edges > o.hardEdges) {
      this._check = null;
      this._lastCheck = this.t;
      var exact = this.reroot(true, Date.now() + o.forcedCheckMs);
      if (exact !== null) this._update();
      return;
    }
    if (!this._check) {
      if (this.edges <= o.compactEdges || this.t - this._lastCheck < o.checkEvery) return;
      this._lastCheck = this.t;
      if (!this._startCheck()) return;
    }
    var done = this._advanceCheck(Date.now() + o.checkBudgetMs);
    if (done === null) return;
    var c = this._check;
    this._check = null;
    this._lastCheck = this.t;
    if (!done) return;
    var init = {};
    c.alive.forEach(function (s) { init[s] = c.b[s].slice(); });
    this._setRoot(init, c.alive.length ? c.tot : null, this._log.slice(c.logIndex));
    this._rooted(true, c.t);
    this._update();
  };

  Controller.prototype._rooted = function (exact, t) {
    this.compactions[exact ? "exact" : "relaxed"] += 1;
    this.compactions.lastT = t;
    this.compactions.lastKind = exact ? "exact" : "relaxed";
  };

  /* Snapshot for an exactness check; false when the answer is known at once (and acted on). */
  Controller.prototype._startCheck = function () {
    var alive = this.filter.alive().sort(byNum), b = {};
    this.shadows.forEach(function (sh) { b[sh.label] = [sh.lo, sh.hi]; });
    var tot = [this.hidden.lo, this.hidden.hi], k = alive.length;
    if (k <= 1 || alive.every(function (s) { return b[s][0] === b[s][1]; })) {
      this.reroot(false);
      this._update();
      return false;
    }
    if (k > this.opts.maxCheckShadows) return false;
    var snap = new SI.CombinatorialFilter(this._rootInit, this._rootTotal);
    snap.extend(this._log);
    this._check = { t: this.t, alive: alive, b: b, tot: tot, filter: snap, logIndex: this._log.length,
                    masks: subsetMasks(k), next: 0 };
    return true;
  };

  /* Compare subset bounds on the snapshot until the deadline (at least one subset): true (exact),
   * false (a subset differs) or null (not finished). */
  Controller.prototype._advanceCheck = function (deadline) {
    var c = this._check;
    while (c.next < c.masks.length) {
      if (!subsetAgrees(c.filter, c.alive, c.masks[c.next], c.b, c.tot)) return false;
      c.next += 1;
      if (Date.now() >= deadline && c.next < c.masks.length) return null;
    }
    return true;
  };

  /* Replace the I-state by one rooted at the alive shadows now.  Returns true (exact), false
   * (relaxed, only when force) or null (not exact and not forced).  deadline (ms since the epoch,
   * optional) bounds the exactness check; running out counts as not exact. */
  Controller.prototype.reroot = function (force, deadline) {
    var filt = this.filter;
    var alive = filt.alive().sort(byNum);
    var b = filt.allBounds();
    var tot = alive.length ? filt.bounds(alive) : [0, 0];
    var exact = this._rerootIsExact(alive, b, tot, deadline);
    if (!exact && !force) return null;
    var init = {};
    alive.forEach(function (s) { init[s] = b[s].slice(); });
    this._check = null;
    this._setRoot(init, alive.length ? tot : null, []);
    this._rooted(exact, this.t);
    return exact;
  };

  /* Masks of the subsets 2 <= |A| < k of k shadows, by size then value. */
  function subsetMasks(k) {
    var masks = [];
    var pop = function (x) {
      var c = 0;
      for (; x; x &= x - 1) c++;
      return c;
    };
    for (var m = 1; m < (1 << k) - 1; m++) if (pop(m) >= 2) masks.push(m);
    masks.sort(function (x, y) { return pop(x) - pop(y) || x - y; });
    return masks;
  }

  function subsetAgrees(filter, alive, mask, b, tot) {
    var A = [], rest = [];
    for (var j = 0; j < alive.length; j++) (mask & (1 << j) ? A : rest).push(alive[j]);
    var want = boxTotalBounds(A, rest, b, tot), got = filter.bounds(A);
    return got[0] === want[0] && got[1] === want[1];
  }

  /* True iff the alive shadows' feasible set equals box(bounds) ∩ {sum in tot}, i.e. the
   * tight bounds agree on every subset (both sides are integral g-polymatroids).  False also when
   * the deadline (optional) passes first. */
  Controller.prototype._rerootIsExact = function (alive, b, tot, deadline) {
    var k = alive.length;
    if (k <= 1) return true;
    if (alive.every(function (s) { return b[s][0] === b[s][1]; })) return true;
    if (k > this.opts.maxCheckShadows) return false;
    var masks = subsetMasks(k);
    for (var i = 0; i < masks.length; i++) {
      if (!subsetAgrees(this.filter, alive, masks[i], b, tot)) return false;
      if (deadline !== undefined && Date.now() >= deadline && i + 1 < masks.length) return false;
    }
    return true;
  };

  /* Compact view of the bipartite I-state for drawing (Fig. 11(c) style):
   * {left: [{id, lo, hi, kind}], right: [{id, lo, hi, kind}], edges: [[li, ri]], hiddenLeft, hiddenGone, rootT}.
   * FOV pseudo-shadows are aggregated into one "enter" (left) and one "exit" (right) node.
   * Only the latest maxGone disappeared vertices are kept (hiddenGone counts the rest, whose
   * edges and left-only vertices are omitted too).  "initial" vertices are the shadows at t0,
   * or at the last re-root when rootT (its time) is not null. */
  Controller.prototype.bipartiteView = function (maxLeft, maxGone) {
    maxLeft = maxLeft || 14;
    maxGone = maxGone || 5;
    var st = this.filter.istate;
    st.flush();
    var isPseudo = SI.bipartite.isPseudo;
    var initial = new Set(st.initial);
    var right = [], rightReach = [];
    var gone = [];
    st.disappeared.forEach(function (bd, v) { if (!isPseudo(v)) gone.push(v); });
    var hiddenGone = Math.max(0, gone.length - maxGone);
    gone = gone.slice(hiddenGone);
    gone.forEach(function (v) {
      var bd = st.disappeared.get(v);
      right.push({ id: v, lo: bd[0], hi: bd[1], kind: "gone" });
      rightReach.push(st.disappearedReach.get(v));
    });
    var bounds = {};
    this.shadows.forEach(function (sh) { bounds[sh.label] = sh; });
    Array.from(st.reach.keys()).sort(byNum).forEach(function (s) {
      var sh = bounds[s];
      right.push({ id: s, lo: sh ? sh.lo : 0, hi: sh ? sh.hi : INF, kind: "alive" });
      rightReach.push(st.reach.get(s));
    });
    var exitK = 0, exitReach = new Set();
    st.disappeared.forEach(function (bd, v) {
      if (!isPseudo(v)) return;
      exitK += bd[0];
      st.disappearedReach.get(v).forEach(function (u) { exitReach.add(u); });
    });
    var used = new Set();
    rightReach.forEach(function (r) { r.forEach(function (u) { used.add(u); }); });
    exitReach.forEach(function (u) { used.add(u); });
    var realLeft = st.sortedLeft(Array.from(used).filter(function (u) { return !isPseudo(u); }));
    var hiddenLeft = Math.max(0, realLeft.length - maxLeft);
    realLeft = realLeft.slice(hiddenLeft);
    var left = realLeft.map(function (u) {
      var bd = st.left.get(u);
      return { id: u, lo: bd[0], hi: bd[1], kind: initial.has(u) ? "initial" : "appeared" };
    });
    var enterK = 0, hasEnter = false;
    st.left.forEach(function (bd, u) {
      if (isPseudo(u) && u[0] === "+") {
        enterK += bd[0];
        hasEnter = true;
      }
    });
    if (hasEnter) left.push({ id: "enter", lo: enterK, hi: enterK, kind: "fov" });
    if (exitK) {
      right.push({ id: "exit", lo: exitK, hi: exitK, kind: "fov" });
      rightReach.push(exitReach);
    }
    var edges = [];
    right.forEach(function (rv, ri) {
      var r = rightReach[ri], fovSeen = false;
      left.forEach(function (lv, li) {
        if (lv.kind !== "fov" && r.has(lv.id)) edges.push([li, ri]);
      });
      if (hasEnter) {
        r.forEach(function (u) { if (isPseudo(u) && u[0] === "+") fovSeen = true; });
        if (fovSeen) edges.push([left.length - 1, ri]);
      }
    });
    return { left: left, right: right, edges: edges, hiddenLeft: hiddenLeft, hiddenGone: hiddenGone,
             rootT: this.compactions.lastT };
  };

  /* Human-readable one-line description of an event. */
  function describeEvent(e) {
    switch (e.type) {
      case "appear":
        return "s" + e.s + " appears" + (e.hi ? " hiding " + e.lo : "");
      case "disappear":
        return "s" + e.s + " disappears, revealing " + e.lo;
      case "split":
        return "s" + e.s + " splits into s" + e.a + ", s" + e.b;
      case "merge":
        return "s" + e.a + " + s" + e.b + " merge into s" + e.s;
      case "enter":
        return (e.k > 1 ? e.k + " targets enter" : "a target enters") + " s" + e.s;
      case "exit":
        return (e.k > 1 ? e.k + " targets exit" : "a target exits") + " s" + e.s;
      default:
        return JSON.stringify(e);
    }
  }

  SI.Controller = Controller;
  SI.demo = {
    FIGURES: FIGURES,
    isPolygonMap: isPolygonMap,
    mapList: mapList,
    pathOptions: pathOptions,
    pathInfo: pathInfo,
    checkWaypoint: checkWaypoint,
    describeEvent: describeEvent,
    splitProbabilities: splitProbabilities,
    shadowAreas: shadowAreas,
    pocketAreas: pocketAreas,
    edgeCount: edgeCount,
    boxTotalBounds: boxTotalBounds,
  };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
