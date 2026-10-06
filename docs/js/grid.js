/* Grid world: map, visibility, shadow tracking and simulation.
 * Bit-for-bit port of shadowinfo/grid/{gridmap,shadows,simulate}.py (DESIGN §5).
 * Events are plain objects in the JSON form of DESIGN §2, e.g.
 * {type: "split", s: 1, a: 2, b: 3}; an unbounded hi is Infinity. */
(function (SI) {
  "use strict";

  if (!SI.Rng && typeof require === "function") require("./rng.js");

  /* ---- gridmap ---------------------------------------------------------- */

  function GridMap(lines) {
    var rows = [];
    for (var i = 0; i < lines.length; i++) {
      var line = lines[i].replace(/\s+$/, "");
      if (line.trim().length) rows.push(line);
    }
    if (!rows.length) throw new Error("empty map");
    this.rows = rows.length;
    this.cols = 0;
    for (i = 0; i < rows.length; i++) this.cols = Math.max(this.cols, rows[i].length);
    this.blocked = [];
    for (i = 0; i < rows.length; i++) {
      for (var c = 0; c < this.cols; c++) {
        var ch = c < rows[i].length ? rows[i][c] : "#";
        if (ch !== "#" && ch !== ".") throw new Error("bad map character " + JSON.stringify(ch));
        this.blocked.push(ch === "#" ? 1 : 0);
      }
    }
    this.freeCells = [];
    for (i = 0; i < this.blocked.length; i++) if (!this.blocked[i]) this.freeCells.push(i);
    this._visCache = new Map();
  }

  GridMap.fromText = function (text) {
    return new GridMap(text.split(/\r?\n/));
  };

  GridMap.prototype.toText = function () {
    var out = [];
    for (var r = 0; r < this.rows; r++) {
      var s = "";
      for (var c = 0; c < this.cols; c++) s += this.blocked[r * this.cols + c] ? "#" : ".";
      out.push(s);
    }
    return out.join("\n") + "\n";
  };

  GridMap.prototype.index = function (r, c) {
    return r * this.cols + c;
  };

  GridMap.prototype.cell = function (i) {
    return [Math.floor(i / this.cols), i % this.cols];
  };

  GridMap.prototype.isFree = function (r, c) {
    return r >= 0 && r < this.rows && c >= 0 && c < this.cols && !this.blocked[r * this.cols + c];
  };

  /* Integer supercover walk between cell centres; a corner pass is blocked
   * only when both side cells are obstacles. */
  GridMap.prototype.lineOfSight = function (r0, c0, r1, c1) {
    var nx = Math.abs(c1 - c0), ny = Math.abs(r1 - r0);
    var sx = c1 > c0 ? 1 : -1, sy = r1 > r0 ? 1 : -1;
    var cols = this.cols, blocked = this.blocked;
    var r = r0, c = c0, ix = 0, iy = 0;
    while (ix < nx || iy < ny) {
      var decision = (1 + 2 * ix) * ny - (1 + 2 * iy) * nx;
      if (decision === 0) {
        if (blocked[r * cols + c + sx] && blocked[(r + sy) * cols + c]) return false;
        r += sy;
        c += sx;
        ix++;
        iy++;
      } else if (decision < 0) {
        c += sx;
        ix++;
      } else {
        r += sy;
        iy++;
      }
      if (blocked[r * cols + c]) return false;
    }
    return true;
  };

  /* Flat 0/1 mask of free cells visible from robot [r, c] (cached); a radius
   * keeps cells with dr^2 + dc^2 <= radius^2 + radius (centre within R + 1/2). */
  GridMap.prototype.visibility = function (robot, radius) {
    var hasR = radius !== undefined && radius !== null;
    var key = this.index(robot[0], robot[1]) + "|" + (hasR ? radius : "");
    var mask = this._visCache.get(key);
    if (mask) return mask;
    var r0 = robot[0], c0 = robot[1];
    if (!this.isFree(r0, c0)) throw new Error("robot cell " + robot + " is not free");
    var r2 = hasR ? radius * radius + radius : 0;
    mask = new Array(this.rows * this.cols).fill(0);
    for (var j = 0; j < this.freeCells.length; j++) {
      var i = this.freeCells[j];
      var r = Math.floor(i / this.cols), c = i % this.cols;
      if (hasR && (r - r0) * (r - r0) + (c - c0) * (c - c0) > r2) continue;
      if (this.lineOfSight(r0, c0, r, c)) mask[i] = 1;
    }
    this._visCache.set(key, mask);
    return mask;
  };

  GridMap.prototype.shadowMask = function (robot, radius) {
    var vis = this.visibility(robot, radius);
    var out = new Array(vis.length);
    for (var i = 0; i < vis.length; i++) out[i] = this.blocked[i] || vis[i] ? 0 : 1;
    return out;
  };

  function expandPath(waypoints, loop) {
    if (loop === undefined) loop = true;
    var pts = waypoints.map(function (p) { return [p[0] | 0, p[1] | 0]; });
    var same = function (p, q) { return p[0] === q[0] && p[1] === q[1]; };
    if (loop && pts.length > 1 && !same(pts[pts.length - 1], pts[0])) pts.push(pts[0]);
    var path = [pts[0].slice()];
    for (var k = 0; k + 1 < pts.length; k++) {
      var r = pts[k][0], c = pts[k][1], r1 = pts[k + 1][0], c1 = pts[k + 1][1];
      if (r !== r1 && c !== c1) throw new Error("segment " + pts[k] + "->" + pts[k + 1] + " is not axis aligned");
      var dr = Math.sign(r1 - r), dc = Math.sign(c1 - c);
      while (r !== r1 || c !== c1) {
        r += dr;
        c += dc;
        path.push([r, c]);
      }
    }
    if (loop && path.length > 1 && same(path[path.length - 1], path[0])) path.pop();
    return path;
  }

  /* ---- shadows ---------------------------------------------------------- */

  function labelComponents(mask, rows, cols) {
    var n = rows * cols;
    var comp = new Array(n).fill(-1);
    var cells = [];
    for (var start = 0; start < n; start++) {
      if (!mask[start] || comp[start] >= 0) continue;
      var k = cells.length;
      comp[start] = k;
      var members = [start], stack = [start];
      while (stack.length) {
        var i = stack.pop();
        var r = Math.floor(i / cols), c = i % cols;
        var nb = [[i - cols, r > 0], [i + cols, r < rows - 1], [i - 1, c > 0], [i + 1, c < cols - 1]];
        for (var q = 0; q < 4; q++) {
          var j = nb[q][0];
          if (nb[q][1] && mask[j] && comp[j] < 0) {
            comp[j] = k;
            members.push(j);
            stack.push(j);
          }
        }
      }
      members.sort(function (a, b) { return a - b; });
      cells.push(members);
    }
    return { comp: comp, cells: cells };
  }

  function ShadowTracker(rows, cols) {
    this.rows = rows;
    this.cols = cols;
    this.labels = new Array(rows * cols).fill(0);
    this.nextLabel = 1;
  }

  ShadowTracker.prototype._newLabel = function () {
    return this.nextLabel++;
  };

  ShadowTracker.prototype.reset = function (mask) {
    var lc = labelComponents(mask, this.rows, this.cols);
    this.nextLabel = 1;
    var compLabels = [];
    for (var k = 0; k < lc.cells.length; k++) compLabels.push(this._newLabel());
    this.labels = lc.comp.map(function (k) { return k >= 0 ? compLabels[k] : 0; });
    return compLabels;
  };

  ShadowTracker.prototype.alive = function () {
    return uniqueLabels(this.labels);
  };

  function uniqueLabels(labels) {
    var seen = new Set();
    for (var i = 0; i < labels.length; i++) if (labels[i]) seen.add(labels[i]);
    return Array.from(seen).sort(function (a, b) { return a - b; });
  }

  /* Returns {events, labels, cells, compLabels, appeared, disappeared, covers,
   * mergedFrom} (see Transition in shadows.py);
   * appear/disappear events carry lo = hi = 0 (the simulator fills counts). */
  ShadowTracker.prototype.update = function (mask) {
    var old = this.labels;
    var lc = labelComponents(mask, this.rows, this.cols);
    var comp = lc.comp, cells = lc.cells, nNew = cells.length;
    var oldLabels = uniqueLabels(old);
    var nodeOf = new Map();
    oldLabels.forEach(function (lab, j) { nodeOf.set(lab, nNew + j); });
    var parent = [];
    for (var x = 0; x < nNew + oldLabels.length; x++) parent.push(x);
    function find(x) {
      while (parent[x] !== x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
      }
      return x;
    }
    for (var i = 0; i < old.length; i++) {
      if (old[i] && comp[i] >= 0) {
        var rx = find(comp[i]), ry = find(nodeOf.get(old[i]));
        if (rx !== ry) parent[Math.max(rx, ry)] = Math.min(rx, ry);
      }
    }
    var newsOf = new Map(), oldsOf = new Map();
    for (var k = 0; k < nNew; k++) {
      var root = find(k);
      if (!newsOf.has(root)) newsOf.set(root, []);
      newsOf.get(root).push(k);
    }
    oldLabels.forEach(function (lab) {
      var root = find(nodeOf.get(lab));
      if (!oldsOf.has(root)) oldsOf.set(root, []);
      oldsOf.get(root).push(lab);
    });
    var byNum = function (a, b) { return a - b; };
    var events = [], compLabels = new Array(nNew).fill(0);
    var tr = { events: events, labels: null, cells: cells, compLabels: compLabels,
               appeared: [], disappeared: [], covers: new Map(), mergedFrom: new Map() };
    Array.from(oldsOf.keys()).sort(byNum).forEach(function (root) {
      if (!newsOf.has(root)) {
        var lab = oldsOf.get(root)[0];
        events.push({ type: "disappear", s: lab, lo: 0, hi: 0 });
        tr.disappeared.push(lab);
      }
    });
    var self = this;
    Array.from(newsOf.keys()).sort(byNum).forEach(function (root) {
      var news = newsOf.get(root), olds = oldsOf.get(root) || [];
      if (!olds.length) {
        var s = self._newLabel();
        events.push({ type: "appear", s: s, lo: 0, hi: 0 });
        tr.appeared.push(s);
        compLabels[news[0]] = s;
        return;
      }
      var cur = olds[0];
      for (var q = 1; q < olds.length; q++) {
        var merged = self._newLabel();
        events.push({ type: "merge", a: cur, b: olds[q], s: merged });
        tr.mergedFrom.set(merged, olds.slice(0, q + 1));
        cur = merged;
      }
      for (var j = 0; j + 1 < news.length; j++) {
        var a = self._newLabel(), b = self._newLabel();
        events.push({ type: "split", s: cur, a: a, b: b });
        compLabels[news[j]] = a;
        tr.covers.set(a, [news[j]]);
        tr.covers.set(b, news.slice(j + 1));
        cur = b;
      }
      compLabels[news[news.length - 1]] = cur;
    });
    this.labels = comp.map(function (k) { return k >= 0 ? compLabels[k] : 0; });
    tr.labels = this.labels;
    return tr;
  };

  /* ---- simulate --------------------------------------------------------- */

  var MOVES = [[0, 0], [-1, 0], [1, 0], [0, -1], [0, 1]];
  var MODES = ["exact", "unknown", "evader"];

  function labelHash(labels) {
    var h = 2166136261;
    for (var i = 0; i < labels.length; i++) h = Math.imul(h ^ labels[i], 16777619) >>> 0;
    return h;
  }

  /* counts: Map or [[label, n], ...]; returns {label: [lo, hi]}. */
  function initialBounds(counts, mode) {
    mode = mode || "exact";
    var entries = counts instanceof Map ? Array.from(counts.entries()) : counts;
    var out = {};
    entries.forEach(function (p) {
      if (mode === "exact") out[p[0]] = [p[1], p[1]];
      else if (mode === "unknown") out[p[0]] = [0, Infinity];
      else if (mode === "evader") out[p[0]] = [0, 1];
      else throw new Error("unknown mode " + JSON.stringify(mode) + "; expected one of " + MODES);
    });
    return out;
  }

  function tally(labels, targets) {
    var counts = new Map();
    uniqueLabels(labels).forEach(function (s) { counts.set(s, 0); });
    targets.forEach(function (i) { if (labels[i]) counts.set(labels[i], counts.get(labels[i]) + 1); });
    return counts;
  }

  function bump(m, key) {
    m.set(key, (m.get(key) || 0) + 1);
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

  function GridSimulator(map, path, nTargets, seed, radius) {
    if (!path.length) throw new Error("empty path");
    this.map = map;
    this.path = path.map(function (p) { return [p[0], p[1]]; });
    this.radius = radius === undefined ? null : radius;
    this.rng = new SI.Rng(seed === undefined ? 1 : seed);
    this.t = 0;
    this.robot = this.path[0];
    this.tracker = new ShadowTracker(map.rows, map.cols);
    this.tracker.reset(map.shadowMask(this.robot, this.radius));
    var free = map.freeCells;
    this.targets = [];
    for (var j = 0; j < (nTargets === undefined ? 10 : nTargets); j++) {
      this.targets.push(free[this.rng.randint(free.length)]);
    }
    this.initialCounts = tally(this.tracker.labels, this.targets);
    this.history = [];
  }

  Object.defineProperty(GridSimulator.prototype, "labels", {
    get: function () { return this.tracker.labels; },
  });

  GridSimulator.prototype.counts = function () {
    return tally(this.tracker.labels, this.targets);
  };

  GridSimulator.prototype.initialCondition = function (mode) {
    return initialBounds(this.initialCounts, mode);
  };

  /* {initial, events} with all events so far (DESIGN §2 JSON shape). */
  GridSimulator.prototype.sequence = function (mode) {
    return { initial: this.initialCondition(mode), events: this.history.slice() };
  };

  GridSimulator.prototype.frame = function () {
    var cols = this.map.cols;
    return {
      t: this.t,
      robot: [this.robot[0], this.robot[1]],
      targets: this.targets.map(function (i) { return [Math.floor(i / cols), i % cols]; }),
      counts: sortedPairs(this.counts()),
      labels_hash: labelHash(this.tracker.labels),
    };
  };

  GridSimulator.prototype.step = function () {
    this.t += 1;
    this.robot = this.path[this.t % this.path.length];
    var old = this.tracker.labels;
    var tr = this.tracker.update(this.map.shadowMask(this.robot, this.radius));
    var nw = tr.labels;
    var appeared = new Set(tr.appeared), disappeared = new Set(tr.disappeared);
    var exits = new Map(), enters = new Map(), revealed = new Map(), hidden = new Map();
    var targets = this.targets, j, i, a, b;
    for (j = 0; j < targets.length; j++) {
      i = targets[j];
      a = old[i];
      b = nw[i];
      if (a && !b) bump(disappeared.has(a) ? revealed : exits, a);
      else if (b && !a) bump(appeared.has(b) ? hidden : enters, b);
    }
    var compCount = new Array(tr.cells.length).fill(0), compOf = new Map();
    tr.compLabels.forEach(function (s, k) { compOf.set(s, k); });
    var stay = new Map();
    for (j = 0; j < targets.length; j++) {
      i = targets[j];
      if (old[i] && nw[i]) {
        compCount[compOf.get(nw[i])] += 1;
        bump(stay, old[i]);
      }
    }
    var created = new Map();
    tr.covers.forEach(function (ks, s) {
      created.set(s, ks.reduce(function (acc, k) { return acc + compCount[k]; }, 0));
    });
    tr.mergedFrom.forEach(function (olds, s) {
      created.set(s, olds.reduce(function (acc, o) { return acc + (stay.get(o) || 0); }, 0));
    });
    var compEvents = tr.events.map(function (e) {
      if (e.type === "appear") {
        var n = hidden.get(e.s) || 0;
        return { type: "appear", s: e.s, lo: n, hi: n };
      }
      if (e.type === "disappear") {
        var m = revealed.get(e.s) || 0;
        return { type: "disappear", s: e.s, lo: m, hi: m };
      }
      return e;
    });
    var events = fovEvents(exits, new Map()).concat(compEvents, fovEvents(new Map(), enters));

    exits = new Map();
    enters = new Map();
    var rows = this.map.rows, cols = this.map.cols, blocked = this.map.blocked;
    for (j = 0; j < targets.length; j++) {
      var mv = MOVES[this.rng.randint(5)];
      i = targets[j];
      var r = Math.floor(i / cols) + mv[0], c = (i % cols) + mv[1];
      if (r < 0 || r >= rows || c < 0 || c >= cols || blocked[r * cols + c]) continue;
      var k = r * cols + c;
      a = nw[i];
      b = nw[k];
      if (a && !b) bump(exits, a);
      else if (b && !a) bump(enters, b);
      else if (a !== b) throw new Error("target moved between shadows " + a + " and " + b + " without an event");
      targets[j] = k;
    }
    events = events.concat(fovEvents(exits, enters));
    Array.prototype.push.apply(this.history, events);

    var f = this.frame();
    f.events = events;
    f.created = sortedPairs(created);
    return f;
  };

  GridSimulator.prototype.run = function (ticks) {
    var frames = [this.frame()];
    for (var t = 0; t < ticks; t++) frames.push(this.step());
    return frames;
  };

  SI.GridMap = GridMap;
  SI.expandPath = expandPath;
  SI.labelComponents = labelComponents;
  SI.ShadowTracker = ShadowTracker;
  SI.GridSimulator = GridSimulator;
  SI.grid = {
    MOVES: MOVES,
    MODES: MODES,
    labelHash: labelHash,
    initialBounds: initialBounds,
    /* Build a simulator on a map from SI.maps (e.g. "office"). */
    simulator: function (name, nTargets, seed, radius) {
      var spec = SI.maps[name];
      var map = GridMap.fromText(spec.text);
      return new GridSimulator(map, expandPath(spec.waypoints, spec.loop), nTargets, seed, radius);
    },
  };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
