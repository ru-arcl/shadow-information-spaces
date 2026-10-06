/* Probabilistic shadow information spaces (T-RO 2012, Sec. VI; ICRA 2010).
 * Port of shadowinfo/probabilistic.py.  The state is a joint pmf over the
 * number of targets in each alive shadow: an ordered label array plus a table
 * Map("x1,x2,..." -> {k: [x1, x2, ...], p}) (Sec. VI-A).  ExactFilter runs
 * Algorithm 1 (PROCESSPROBABILITYMASS) and Algorithm 2 (PROCESSFOVEVENT);
 * TruncatedFilter adds the TR / RT / RT-LA heuristics (Sec. VI-E, VII-B);
 * monteCarlo is the rejection-sampling baseline (Sec. VI-E).
 *
 * Observations are plain objects: the component/FOV events of events.js,
 * {type: "fov", s, y} with y in "enter" | "exit" | "null", "dist" on
 * appear/disappear ({count: p} or [[count, p], ...]) and "p" on split
 * (per-event binomial parameter for child a). */
(function (SI) {
  "use strict";

  if (typeof require === "function") {
    if (!SI.events) require("./events.js");
    if (!SI.Rng) {
      try {
        require("./rng.js");
      } catch (err) {
        /* fall back to the local generator below */
      }
    }
  }

  var Rng = SI.Rng || (function () {
    function Mulberry32(seed) {
      this.state = (seed === undefined ? 1 : seed) >>> 0;
    }
    Mulberry32.prototype.nextU32 = function () {
      this.state = (this.state + 0x6d2b79f5) >>> 0;
      var t = this.state;
      t = Math.imul(t ^ (t >>> 15), t | 1) >>> 0;
      t = (t ^ ((t + (Math.imul(t ^ (t >>> 7), t | 61) >>> 0)) >>> 0)) >>> 0;
      return (t ^ (t >>> 14)) >>> 0;
    };
    Mulberry32.prototype.random = function () {
      return this.nextU32() / 4294967296;
    };
    Mulberry32.prototype.randint = function (n) {
      return Math.floor((this.nextU32() * n) / 4294967296);
    };
    return Mulberry32;
  })();

  var FOV_EVENTS = ["enter", "exit", "null"];
  var DELTA = { enter: 1, null: 0, exit: -1 };
  var TRUNCATION_MODES = ["TR", "RT", "RT-LA"];

  class InconsistentObservationError extends Error {}
  InconsistentObservationError.prototype.name = "InconsistentObservationError";

  /* A truncated filter was left with no mass after it had truncated (Tables IV, V "failure").
   * This does not prove truncation was the cause: observations that are inconsistent even for
   * ExactFilter raise it too once a truncation has happened; run ExactFilter to tell them apart. */
  class TruncationFailure extends InconsistentObservationError {}
  TruncationFailure.prototype.name = "TruncationFailure";

  function keyStr(k) {
    return k.join(",");
  }

  function cmpKeys(a, b) {
    for (var i = 0; i < a.length && i < b.length; i++) if (a[i] !== b[i]) return a[i] - b[i];
    return a.length - b.length;
  }

  function addTo(table, k, p) {
    var ks = keyStr(k), e = table.get(ks);
    if (e === undefined) table.set(ks, { k: k, p: 0 + p });
    else e.p += p;
  }

  /* Count distribution as sorted [[n, p], ...] with zero weights removed. */
  function distItems(dist) {
    var items = Array.isArray(dist) ? dist.map(function (x) { return [+x[0], x[1]]; })
      : Object.keys(dist).map(function (n) { return [+n, dist[n]]; });
    items = items.filter(function (x) { return x[1]; }).sort(function (a, b) { return a[0] - b[0] || a[1] - b[1]; });
    if (!items.length || items.some(function (x) { return x[0] < 0 || x[1] < 0; }))
      throw new RangeError("bad count distribution " + JSON.stringify(dist));
    return items;
  }

  var combCache = [];

  function comb(n, a) {
    var row = combCache[n];
    if (row === undefined) {
      row = combCache[n] = [];
      var c = 1n;
      for (var j = 0; j <= n; j++) {
        row.push(Number(c));
        c = (c * BigInt(n - j)) / BigInt(j + 1);
      }
    }
    return row[a];
  }

  var logFactCache = [0];

  function logFact(n) {
    for (var j = logFactCache.length; j <= n; j++) logFactCache.push(logFactCache[j - 1] + Math.log(j));
    return logFactCache[n];
  }

  /* Binomial pmf in log space, renormalised; used once comb(n, a) overflows (n >= 1030). */
  function binomialLogPmf(n, p) {
    if (p === 0 || p === 1) return [p ? [n, 0, 1] : [0, n, 1]];
    var lp = Math.log(p), lq = Math.log1p(-p), c = logFact(n), out = [], z = 0, a;
    for (a = 0; a <= n; a++) {
      var v = Math.exp(c - logFact(a) - logFact(n - a) + a * lp + (n - a) * lq);
      z += v;
      if (v) out.push([a, n - a, v]);
    }
    for (a = 0; a < out.length; a++) out[a][2] /= z;
    return out;
  }

  /* Split rule: each target independently enters the first child with probability p.
   * rule(n) -> [[nA, nB, prob], ...]. */
  function binomialSplit(p) {
    if (p === undefined) p = 0.5;
    if (!(p >= 0 && p <= 1)) throw new RangeError("p must lie in [0, 1]");
    var cache = new Map();
    var rule = function (n) {
      var out = cache.get(n);
      if (out === undefined) {
        var q = 1 - p;
        out = [];
        for (var a = 0; a <= n; a++) {
          var c = comb(n, a);
          if (c === Infinity) {
            out = binomialLogPmf(n, p);
            break;
          }
          var v = c * Math.pow(p, a) * Math.pow(q, n - a);
          if (v) out.push([a, n - a, v]);
        }
        cache.set(n, out);
      }
      return out;
    };
    rule.p = p;
    rule.binomial = true;
    return rule;
  }

  /* P(e | y) = 1 iff e == y: noise-free FOV sensing. */
  function perfectObsModel() {
    var m = {};
    FOV_EVENTS.forEach(function (y) {
      m[y] = {};
      FOV_EVENTS.forEach(function (e) { m[y][e] = e === y ? 1 : 0; });
    });
    return m;
  }

  /* True-positive rate p (Sec. VI-D uses 0.9); with nullEvents the error mass is shared. */
  function symmetricObsModel(p, nullEvents) {
    if (p === undefined) p = 0.9;
    if (!nullEvents) {
      return {
        enter: { enter: p, exit: 1 - p, null: 0 },
        exit: { enter: 1 - p, exit: p, null: 0 },
        null: { enter: 0, exit: 0, null: 1 },
      };
    }
    var q = (1 - p) / 2, m = {};
    FOV_EVENTS.forEach(function (y) {
      m[y] = {};
      FOV_EVENTS.forEach(function (e) { m[y][e] = e === y ? p : q; });
    });
    return m;
  }

  var SUM_TOL = 1e-9;

  /* Weights {key: w} rescaled to sum to one.  Negative weights throw; all-zero weights and
   * weights that already sum to one up to rounding (SUM_TOL) are returned unchanged. */
  function normalized(weights, what) {
    var keys = Object.keys(weights), z = 0;
    keys.forEach(function (k) {
      if (weights[k] < 0) throw new RangeError(what + ": negative weight in " + JSON.stringify(weights));
      z += weights[k];
    });
    if (!z || Math.abs(z - 1) <= SUM_TOL) return weights;
    var out = {};
    keys.forEach(function (k) { out[k] = weights[k] / z; });
    return out;
  }

  /* obsModel[y][e] = P(e | y), zeros filled in.  Unknown y/e keys and negative weights throw;
   * a row that does not sum to one is normalised; an omitted or all-zero row makes y impossible. */
  function obsModelFromDict(d) {
    var unknown = [];
    Object.keys(d).forEach(function (y) {
      if (FOV_EVENTS.indexOf(y) < 0) {
        unknown.push(JSON.stringify(y));
        return;
      }
      Object.keys(d[y] || {}).forEach(function (e) {
        if (FOV_EVENTS.indexOf(e) < 0) unknown.push(y + "/" + JSON.stringify(e));
      });
    });
    if (unknown.length)
      throw new RangeError("observation model: unknown FOV event(s) " + unknown.join(", ") + "; expected " + FOV_EVENTS);
    var m = {};
    FOV_EVENTS.forEach(function (y) {
      var row = {};
      FOV_EVENTS.forEach(function (e) { row[e] = (d[y] && d[y][e]) || 0; });
      m[y] = normalized(row, "observation model row " + JSON.stringify(y));
    });
    return m;
  }

  /* rule(n); a custom (non-binomial) rule's output is checked and normalised. */
  function splitDist(rule, n) {
    var dist = rule(n);
    if (rule.binomial) return dist;
    var w = {}, any = false;
    dist.forEach(function (r) {
      if (r[0] < 0 || r[1] < 0 || r[0] + r[1] !== n)
        throw new RangeError("split rule: " + JSON.stringify(dist) + " is not a split of " + n + " targets");
      w[r[0]] = (w[r[0]] || 0) + r[2];
      if (r[2]) any = true;
    });
    if (!any) throw new RangeError("split rule: no probability mass for " + n + " targets");
    w = normalized(w, "split rule output for " + n + " targets");
    return Object.keys(w).map(Number).sort(function (a, b) { return a - b; })
      .map(function (a) { return [a, n - a, w[a]]; });
  }

  function splitRuleFromDict(d) {
    if ((d.type || "binomial") !== "binomial") throw new RangeError("unknown split rule " + JSON.stringify(d));
    return binomialSplit(d.p === undefined ? 0.5 : d.p);
  }

  /* One queue item from JSON (fixture format of shadowinfo.probabilistic.ProbSequence). */
  function observationFromDict(d) {
    if (d.type === "fov") {
      if (FOV_EVENTS.indexOf(d.y) < 0) throw new RangeError("bad FOV observation " + JSON.stringify(d.y));
      return { type: "fov", s: d.s, y: d.y };
    }
    if ((d.type === "appear" || d.type === "disappear") && d.dist) {
      var items = distItems(d.dist);
      return { type: d.type, s: d.s, lo: items[0][0], hi: items[items.length - 1][0], dist: items };
    }
    var e = SI.events.eventFromDict(d);
    if (d.type === "split" && d.p !== undefined) e.p = d.p;
    return e;
  }

  function observationToDict(o) {
    if (o.type === "fov") return { type: "fov", s: o.s, y: o.y };
    var d = SI.events.eventToDict(o);
    if (o.dist) {
      d.dist = {};
      distItems(o.dist).forEach(function (x) { d.dist[x[0]] = x[1]; });
    }
    if (o.type === "split" && o.p !== undefined) d.p = o.p;
    return d;
  }

  /* {labels, joint: [[key, p], ...], observations} from the JSON format; repeated joint keys are summed. */
  function probSequenceFromDict(d) {
    var joint = new Map();
    d.joint.forEach(function (kp) {
      var ks = keyStr(kp[0]), e = joint.get(ks);
      if (e === undefined) joint.set(ks, [kp[0].slice(), kp[1]]);
      else e[1] += kp[1];
    });
    return {
      labels: d.labels.slice(),
      joint: Array.from(joint.values()),
      observations: d.events.map(observationFromDict),
    };
  }

  /* ---- joint pmf ------------------------------------------------------------ */

  function JointPmf(labels, joint) {
    this.labels = labels.slice();
    this.table = new Map();
    var self = this;
    var items = joint instanceof Map ? Array.from(joint.values()).map(function (e) { return [e.k, e.p]; }) : joint;
    items.forEach(function (kp) {
      if (!kp[1]) return;
      if (kp[0].length !== self.labels.length) throw new RangeError("joint table keys must match the labels");
      if (kp[1] < 0) throw new RangeError("negative probability mass");
      self.table.set(keyStr(kp[0]), { k: kp[0].slice(), p: kp[1] });
    });
  }

  Object.defineProperty(JointPmf.prototype, "entries", {
    get: function () { return this.table.size; },
  });

  JointPmf.prototype.total = function () {
    var z = 0;
    this.table.forEach(function (e) { z += e.p; });
    return z;
  };

  JointPmf.prototype.index = function (s) {
    var i = this.labels.indexOf(s);
    if (i < 0) throw new RangeError("shadow " + s + " is not alive");
    return i;
  };

  /* {count: P(s = count)} for counts with positive mass. */
  JointPmf.prototype.marginal = function (s) {
    var i = this.index(s), m = new Map();
    this.table.forEach(function (e) { m.set(e.k[i], (m.get(e.k[i]) || 0) + e.p); });
    var out = {};
    Array.from(m.keys()).sort(function (a, b) { return a - b; }).forEach(function (x) { out[x] = m.get(x); });
    return out;
  };

  /* {label: E[s]} for every alive shadow. */
  JointPmf.prototype.expectedCounts = function () {
    var labels = this.labels, sums = labels.map(function () { return 0; });
    this.table.forEach(function (e) {
      for (var i = 0; i < labels.length; i++) sums[i] += e.k[i] * e.p;
    });
    var out = {};
    labels.forEach(function (s, i) { out[s] = sums[i]; });
    return out;
  };

  /* Probability that every shadow in assignment {label: count} holds that count. */
  JointPmf.prototype.probability = function (assignment) {
    var self = this;
    var idx = Object.keys(assignment).map(function (s) { return [self.index(+s), assignment[s]]; });
    var z = 0;
    this.table.forEach(function (e) {
      if (idx.every(function (ix) { return e.k[ix[0]] === ix[1]; })) z += e.p;
    });
    return z;
  };

  /* Entries sorted by key: [[key, p], ...]. */
  JointPmf.prototype.asList = function () {
    return Array.from(this.table.values()).sort(function (a, b) { return cmpKeys(a.k, b.k); })
      .map(function (e) { return [e.k.slice(), e.p]; });
  };

  /* The same pmf with its key columns in the order of labels, a permutation of this.labels
   * (e.g. to compare a monteCarlo result, whose labels are sorted, with an ExactFilter). */
  JointPmf.prototype.reordered = function (labels) {
    var a = labels.slice().sort(function (x, y) { return x - y; });
    var b = this.labels.slice().sort(function (x, y) { return x - y; });
    if (a.length !== b.length || a.some(function (x, j) { return x !== b[j]; }))
      throw new RangeError("labels " + JSON.stringify(labels) + " are not a permutation of " + JSON.stringify(this.labels));
    var idx = labels.map(function (s) { return this.index(s); }, this);
    return new JointPmf(labels, Array.from(this.table.values()).map(function (e) {
      return [idx.map(function (i) { return e.k[i]; }), e.p];
    }));
  };

  JointPmf.prototype.toJSON = function () {
    return { labels: this.labels.slice(), joint: this.asList() };
  };

  JointPmf.prototype._normalize = function () {
    var z = this.total();
    if (!z) throw new InconsistentObservationError("no probability mass left");
    if (z !== 1) {
      var table = this.table;
      table.forEach(function (e, ks) {
        e.p = e.p / z;
        if (!e.p) table.delete(ks);
      });
    }
  };

  /* ---- exact filter (Algorithms 1 and 2) -------------------------------------- */

  /* opts: {splitRule: rule(n) -> [[nA, nB, p], ...] (default binomialSplit(0.5)),
   *        obsModel: {y: {e: P(e | y)}} (default noise-free)}. */
  function ExactFilter(labels, joint, opts) {
    opts = opts || {};
    JointPmf.call(this, labels, joint);
    this._normalize();
    this.splitRule = opts.splitRule || binomialSplit(0.5);
    this.obsModel = opts.obsModel ? obsModelFromDict(opts.obsModel) : perfectObsModel();
    this.peakEntries = this.entries;
    this._splitRules = new Map();
  }

  ExactFilter.prototype = Object.create(JointPmf.prototype);
  ExactFilter.prototype.constructor = ExactFilter;

  /* Known initial counts {label: n}: P = 1 on that single entry. */
  ExactFilter.fromCounts = function (counts, opts) {
    var labels = Object.keys(counts).map(Number);
    return new ExactFilter(labels, [[labels.map(function (s) { return counts[s]; }), 1]], opts);
  };

  /* Independent initial marginals {label: {count: p}}. */
  ExactFilter.fromIndependent = function (dists, opts) {
    var labels = Object.keys(dists).map(Number);
    var rows = [[[], 1]];
    labels.forEach(function (s) {
      var next = [];
      rows.forEach(function (r) {
        Object.keys(dists[s]).forEach(function (n) {
          next.push([r[0].concat([+n]), r[1] * dists[s][n]]);
        });
      });
      rows = next;
    });
    return new ExactFilter(labels, rows.filter(function (r) { return r[1]; }), opts);
  };

  ExactFilter.prototype.copy = function () {
    var other = Object.create(Object.getPrototypeOf(this));
    Object.assign(other, this);
    other.labels = this.labels.slice();
    other.table = new Map();
    this.table.forEach(function (e, ks) { other.table.set(ks, { k: e.k, p: e.p }); });
    if (this.rng) {
      other.rng = new Rng();
      other.rng.state = this.rng.state;
    }
    return other;
  };

  /* Process the whole queue (Algorithm 1) and return this. */
  ExactFilter.prototype.run = function (observations) {
    for (var t = 0; t < observations.length; t++) this._apply(observations[t], observations, t + 1);
    return this;
  };

  /* Process one queue item; upcoming (the rest of the queue) is used only by RT-LA. */
  ExactFilter.prototype.apply = function (o, upcoming) {
    this._apply(o, upcoming || null, 0);
  };

  ExactFilter.prototype._apply = function (o, queue, start) {
    switch (o.type) {
      case "fov":
        this._fov(o.s, o.y);
        break;
      case "split":
        this._split(o);
        break;
      case "merge":
        this._merge(o.a, o.b, o.s);
        break;
      case "appear":
        this._appear(o);
        break;
      case "disappear":
        this._disappear(o);
        break;
      case "enter":
        this._shift(o.s, o.k === undefined ? 1 : o.k);
        break;
      case "exit":
        this._shift(o.s, -(o.k === undefined ? 1 : o.k));
        break;
      default:
        throw new TypeError("unknown observation " + JSON.stringify(o));
    }
    this.peakEntries = Math.max(this.peakEntries, this.entries);
    this._afterUpdate(queue, start);
  };

  ExactFilter.prototype._afterUpdate = function () {};

  ExactFilter.prototype._appear = function (e) {
    if (this.labels.indexOf(e.s) >= 0) throw new RangeError("shadow " + e.s + " already alive");
    var dist;
    if (e.dist) {
      dist = distItems(e.dist);
    } else {
      var hi = SI.events.hiFromJson(e.hi === undefined ? 0 : e.hi), lo = e.lo || 0;
      if (hi === Infinity) throw new RangeError("an appear event needs a finite upper bound or a dist");
      var w = hi > lo ? 1 / (hi - lo + 1) : 1;
      dist = [];
      for (var n = lo; n <= hi; n++) dist.push([n, w]);
    }
    var out = new Map();
    this.table.forEach(function (en) {
      for (var j = 0; j < dist.length; j++) {
        var v = en.p * dist[j][1];
        if (!v) continue;
        var nk = en.k.concat([dist[j][0]]);
        out.set(keyStr(nk), { k: nk, p: v });
      }
    });
    this.table = out;
    this.labels.push(e.s);
    this._normalize();
  };

  ExactFilter.prototype._disappear = function (e) {
    var i = this.index(e.s), like;
    if (e.dist) {
      var dm = new Map(distItems(e.dist));
      like = function (x) { return dm.get(x); };
    } else {
      var lo = e.lo || 0, hi = SI.events.hiFromJson(e.hi === undefined ? 0 : e.hi);
      like = function (x) { return lo <= x && x <= hi ? 1 : 0; };
    }
    var out = new Map();
    this.table.forEach(function (en) {
      var w = like(en.k[i]), v = w ? en.p * w : 0;
      if (v) addTo(out, en.k.slice(0, i).concat(en.k.slice(i + 1)), v);
    });
    this.table = out;
    this.labels.splice(i, 1);
    this._normalizeOrFail("disappear of shadow " + e.s);
  };

  ExactFilter.prototype._split = function (e) {
    var i = this.index(e.s), rule;
    var rest = this.labels.slice(0, i).concat(this.labels.slice(i + 1));
    if (e.a === e.b || rest.indexOf(e.a) >= 0 || rest.indexOf(e.b) >= 0)
      throw new RangeError("split of shadow " + e.s + ": children " + e.a + ", " + e.b + " must be distinct and not alive");
    if (e.p !== undefined) {
      rule = this._splitRules.get(e.p);
      if (rule === undefined) {
        rule = binomialSplit(e.p);
        this._splitRules.set(e.p, rule);
      }
    } else {
      rule = this.splitRule;
    }
    var out = new Map(), dists = new Map();
    this.table.forEach(function (en) {
      var restK = en.k.slice(0, i).concat(en.k.slice(i + 1));
      var parts = dists.get(en.k[i]);
      if (parts === undefined) dists.set(en.k[i], (parts = splitDist(rule, en.k[i])));
      for (var j = 0; j < parts.length; j++) {
        var v = en.p * parts[j][2];
        if (v) addTo(out, restK.concat([parts[j][0], parts[j][1]]), v);
      }
    });
    this.table = out;
    this.labels = rest.concat([e.a, e.b]);
  };

  ExactFilter.prototype._merge = function (a, b, s) {
    var i = this.index(a), j = this.index(b);
    if (i === j) throw new RangeError("cannot merge a shadow with itself");
    if (this.labels.indexOf(s) >= 0 && s !== a && s !== b) throw new RangeError("merge into shadow " + s + ": already alive");
    var keep = [];
    for (var q = 0; q < this.labels.length; q++) if (q !== i && q !== j) keep.push(q);
    var out = new Map();
    this.table.forEach(function (en) {
      var nk = keep.map(function (q) { return en.k[q]; });
      nk.push(en.k[i] + en.k[j]);
      addTo(out, nk, en.p);
    });
    this.table = out;
    var labels = this.labels;
    this.labels = keep.map(function (q) { return labels[q]; });
    this.labels.push(s);
  };

  /* Algorithm 2.  When x_s = 0 the exit branch is impossible and the enter/null branches of
   * that entry are rescaled to its mass; an entry with no feasible branch is dropped.  An
   * observation whose model row is all zero is impossible. */
  ExactFilter.prototype._fov = function (s, y) {
    var i = this.index(s), model = this.obsModel[y];
    if (!model) throw new RangeError("bad FOV observation " + JSON.stringify(y));
    var branches = [];
    FOV_EVENTS.forEach(function (e) { if (model[e]) branches.push([DELTA[e], model[e]]); });
    var z0 = 0;
    branches.forEach(function (b) { if (b[0] >= 0) z0 += b[1]; });
    var dropped = false, out = new Map();
    this.table.forEach(function (en) {
      var x = en.k[i];
      if (x === 0 && !z0) {
        dropped = true;
        return;
      }
      for (var j = 0; j < branches.length; j++) {
        var d = branches[j][0], w = branches[j][1];
        if (x === 0) {
          if (d < 0) continue;
          w = w / z0;
        }
        var v = en.p * w;
        if (!v) continue;
        var nk = en.k.slice();
        nk[i] = x + d;
        addTo(out, nk, v);
      }
    });
    this.table = out;
    if (dropped || !out.size) this._normalizeOrFail("FOV observation " + JSON.stringify(y) + " on shadow " + s);
  };

  /* Noise-free FOV events: Enter adds k targets to every entry; Exit removes k, ruling out x_s < k. */
  ExactFilter.prototype._shift = function (s, k) {
    var i = this.index(s), out = new Map();
    this.table.forEach(function (en) {
      var x = en.k[i] + k;
      if (x >= 0) {
        var nk = en.k.slice();
        nk[i] = x;
        out.set(keyStr(nk), { k: nk, p: en.p });
      }
    });
    var dropped = out.size < this.table.size;
    this.table = out;
    if (dropped) this._normalizeOrFail("exit of " + -k + " target(s) from shadow " + s);
  };

  ExactFilter.prototype._normalizeOrFail = function (what) {
    if (!this.table.size) throw new InconsistentObservationError(what + ": no probability mass left");
    this._normalize();
  };

  /* ---- truncation heuristics (Sec. VI-E, VII-B) ------------------------------- */

  /* opts: ExactFilter opts plus {maxEntries, mode: "TR" | "RT" | "RT-LA", seed (default 1),
   * lookahead: [4, 2], lookaheadCountsFov: true}.  RT-LA needs the future queue: use run(),
   * or pass upcoming to apply(); without it RT-LA behaves as RT. */
  function TruncatedFilter(labels, joint, opts) {
    opts = opts || {};
    var mode = opts.mode || "TR";
    if (TRUNCATION_MODES.indexOf(mode) < 0) throw new RangeError("mode must be one of " + TRUNCATION_MODES);
    if (!(opts.maxEntries >= 1)) throw new RangeError("maxEntries must be positive");
    ExactFilter.call(this, labels, joint, opts);
    this.maxEntries = opts.maxEntries;
    this.mode = mode;
    this.rng = new Rng(opts.seed === undefined ? 1 : opts.seed);
    this.lookahead = opts.lookahead || [4, 2];
    this.lookaheadCountsFov = opts.lookaheadCountsFov === undefined ? true : opts.lookaheadCountsFov;
    this.truncations = 0;
    this.truncatedMass = 0;
  }

  TruncatedFilter.prototype = Object.create(ExactFilter.prototype);
  TruncatedFilter.prototype.constructor = TruncatedFilter;

  TruncatedFilter.prototype._skipForLookahead = function (queue, start) {
    if (this.mode !== "RT-LA" || !queue || start >= queue.length) return false;
    var nDis = this.lookahead[0], nMer = this.lookahead[1], window = Math.max(nDis, nMer), pos = 0;
    for (var t = start; t < queue.length && pos < window; t++) {
      var o = queue[t];
      if (!this.lookaheadCountsFov && !SI.events.isComponent(o)) continue;
      pos += 1;
      if ((o.type === "disappear" && pos <= nDis) || (o.type === "merge" && pos <= nMer)) return true;
    }
    return false;
  };

  TruncatedFilter.prototype._afterUpdate = function (queue, start) {
    if (this.entries <= this.maxEntries || this._skipForLookahead(queue, start)) return;
    var items = Array.from(this.table.values()).sort(function (a, b) { return cmpKeys(a.k, b.k); });
    var scores = this.mode === "TR" ? items.map(function (e) { return e.p; })
      : items.map(function (e) { return e.p * this.rng.random(); }, this);
    var order = items.map(function (_, j) { return j; });
    order.sort(function (a, b) { return scores[b] - scores[a] || a - b; });
    var before = this.total();
    var kept = new Map();
    for (var j = 0; j < this.maxEntries; j++) {
      var e = items[order[j]];
      kept.set(keyStr(e.k), e);
    }
    this.table = kept;
    this.truncatedMass += before - this.total();
    this.truncations += 1;
    this._normalize();
  };

  TruncatedFilter.prototype._normalizeOrFail = function (what) {
    if (!this.table.size && this.truncations)
      throw new TruncationFailure(what + ": no probability mass left after " + this.truncations + " truncation(s)");
    ExactFilter.prototype._normalizeOrFail.call(this, what);
  };

  /* ---- Monte Carlo baseline (Sec. VI-E) ------------------------------------------ */

  function sample(rng, items) {
    var total = 0, j;
    for (j = 0; j < items.length; j++) total += items[j][1];
    var u = rng.random() * total, acc = 0;
    for (j = 0; j < items.length; j++) {
      acc += items[j][1];
      if (u < acc) return items[j][0];
    }
    return items[items.length - 1][0];
  }

  /* Rejection sampling until `trials` consistent runs; seq = {labels, joint, observations}.
   * The result's labels are sorted (ExactFilter keeps event order): use reordered() to compare.
   * opts: {trials (1000), seed (1), splitRule, obsModel, maxAttempts (1000 * trials)}.
   * Returns a JointPmf with successes and attempts. */
  function monteCarlo(seq, opts) {
    opts = opts || {};
    var trials = opts.trials === undefined ? 1000 : opts.trials;
    var splitRule = opts.splitRule || binomialSplit(0.5);
    var model = opts.obsModel ? obsModelFromDict(opts.obsModel) : perfectObsModel();
    var maxAttempts = opts.maxAttempts === undefined ? 1000 * trials : opts.maxAttempts;
    var rng = new Rng(opts.seed === undefined ? 1 : opts.seed);
    var init = seq.joint.filter(function (kp) { return kp[1]; })
      .sort(function (a, b) { return cmpKeys(a[0], b[0]) || a[1] - b[1]; });
    var fovItems = {};
    FOV_EVENTS.forEach(function (y) {
      fovItems[y] = FOV_EVENTS.filter(function (e) { return model[y][e]; }).map(function (e) { return [e, model[y][e]]; });
    });
    var obs = seq.observations.map(function (o) {
      return o.dist ? Object.assign({}, o, { dist: distItems(o.dist) }) : o;
    });
    var finalLabels = null, counts = new Map(), successes = 0, attempts = 0;
    while (successes < trials) {
      if (attempts >= maxAttempts)
        throw new InconsistentObservationError("only " + successes + " of " + attempts +
                                               " Monte Carlo trials were consistent");
      attempts += 1;
      var x = new Map(), start = sample(rng, init), ok = true, n, i;
      seq.labels.forEach(function (s, j) { x.set(s, start[j]); });
      for (i = 0; i < obs.length; i++) {
        var o = obs[i];
        if (o.type === "fov") {
          n = x.get(o.s);
          var items = fovItems[o.y].filter(function (ew) { return n > 0 || ew[0] !== "exit"; });
          if (!items.length) {
            ok = false;
            break;
          }
          x.set(o.s, n + DELTA[sample(rng, items)]);
        } else if (o.type === "split") {
          n = x.get(o.s);
          x.delete(o.s);
          if (o.a === o.b || x.has(o.a) || x.has(o.b))
            throw new RangeError("split of shadow " + o.s + ": children " + o.a + ", " + o.b + " must be distinct and not alive");
          var rule = o.p !== undefined ? binomialSplit(o.p) : splitRule, na = 0;
          if (rule.binomial) {
            for (var t = 0; t < n; t++) if (rng.random() < rule.p) na += 1;
          } else {
            na = sample(rng, splitDist(rule, n).map(function (r) { return [r[0], r[2]]; }));
          }
          x.set(o.a, na);
          x.set(o.b, n - na);
        } else if (o.type === "merge") {
          if (o.a === o.b) throw new RangeError("cannot merge a shadow with itself");
          var m = x.get(o.a) + x.get(o.b);
          x.delete(o.a);
          x.delete(o.b);
          if (x.has(o.s)) throw new RangeError("merge into shadow " + o.s + ": already alive");
          x.set(o.s, m);
        } else if (o.type === "appear") {
          if (x.has(o.s)) throw new RangeError("shadow " + o.s + " already alive");
          if (o.dist) {
            x.set(o.s, sample(rng, o.dist));
          } else {
            var hi = SI.events.hiFromJson(o.hi === undefined ? 0 : o.hi), lo = o.lo || 0;
            if (hi === Infinity) throw new RangeError("an appear event needs a finite upper bound or a dist");
            x.set(o.s, lo + rng.randint(hi - lo + 1));
          }
        } else if (o.type === "disappear") {
          n = x.get(o.s);
          x.delete(o.s);
          if (o.dist) ok = sample(rng, o.dist) === n;
          else ok = (o.lo || 0) <= n && n <= SI.events.hiFromJson(o.hi === undefined ? 0 : o.hi);
          if (!ok) break;
        } else if (o.type === "enter" || o.type === "exit") {
          var k = o.k === undefined ? 1 : o.k;
          x.set(o.s, x.get(o.s) + (o.type === "enter" ? k : -k));
          if (x.get(o.s) < 0) {
            ok = false;
            break;
          }
        } else {
          throw new TypeError("unknown observation " + JSON.stringify(o));
        }
      }
      if (!ok) continue;
      if (finalLabels === null) finalLabels = Array.from(x.keys()).sort(function (a, b) { return a - b; });
      var key = finalLabels.map(function (s) { return x.get(s); }), ks = keyStr(key);
      var c = counts.get(ks);
      if (c === undefined) counts.set(ks, { k: key, c: 1 });
      else c.c += 1;
      successes += 1;
    }
    var table = Array.from(counts.values()).map(function (e) { return [e.k, e.c / successes]; });
    var res = new JointPmf(finalLabels || [], table);
    res.successes = successes;
    res.attempts = attempts;
    return res;
  }

  SI.InconsistentObservationError = InconsistentObservationError;
  SI.TruncationFailure = TruncationFailure;
  SI.JointPmf = JointPmf;
  SI.ExactFilter = ExactFilter;
  SI.TruncatedFilter = TruncatedFilter;
  SI.monteCarlo = monteCarlo;
  SI.prob = {
    FOV_EVENTS: FOV_EVENTS,
    TRUNCATION_MODES: TRUNCATION_MODES,
    binomialSplit: binomialSplit,
    perfectObsModel: perfectObsModel,
    symmetricObsModel: symmetricObsModel,
    obsModelFromDict: obsModelFromDict,
    splitRuleFromDict: splitRuleFromDict,
    observationFromDict: observationFromDict,
    observationToDict: observationToDict,
    probSequenceFromDict: probSequenceFromDict,
    distItems: distItems,
    Rng: Rng,
  };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
