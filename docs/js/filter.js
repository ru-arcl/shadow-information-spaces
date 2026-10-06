/* Combinatorial filter for nondeterministically moving targets (T-RO 2012, Sec. V).
 * Port of shadowinfo/nondeterministic.py.  Queries are flows on the augmented
 * graph of Fig. 11(d):
 *
 *   T -> S      [0, inf)                 closes the circulation
 *   S -> I      [N_lo, N_hi] or [0, inf) total targets in the shadows at t0
 *   I -> L_i    [l_i, u_i]               shadow i at t0
 *   S -> L_i    [l_i, u_i]               appeared shadow / FOV enter pseudo-shadow
 *   L_i -> R_j  [0, inf)                 bipartite edges
 *   R_j -> T    [l_j, u_j]               disappeared shadow / FOV exit pseudo-shadow
 *   R_j -> T    [0, inf)                 shadow alive now
 *
 * Feasible circulations are exactly the target distributions consistent with
 * the observations (Prop. 5), so bounds are min/max flows into the query set.
 * Inconsistent observations make queries throw SI.InfeasibleError.  Bounds are
 * [lo, hi] arrays; hi may be Infinity. */
(function (SI) {
  "use strict";

  if (!SI.events && typeof require === "function") require("./events.js");
  if (!SI.FlowNetwork && typeof require === "function") require("./maxflow.js");
  if (!SI.BipartiteIState && typeof require === "function") require("./bipartite.js");

  var INF = Infinity;
  var FlowNetwork = SI.FlowNetwork, InfeasibleError = SI.InfeasibleError;
  var S = "S", T = "T", I = "I", G = "G";

  function L(v) {
    return "L:" + v;
  }

  function R(tag, v) {
    return "R:" + tag + ":" + v;
  }

  function rightKey(tag, v) {
    return tag + ":" + v;
  }

  function uniq(xs) {
    return Array.from(new Set(xs));
  }

  /* initial: {label: [lo, hi]} (hi may be Infinity or null); total: optional [lo, hi] on the
   * number of targets in the shadows at t0; fov: "batch" (default) or "naive". */
  function CombinatorialFilter(initial, total, fov) {
    this.istate = new SI.BipartiteIState(initial, fov);
    this.total = total ? [total[0], SI.events.hiFromJson(total[1])] : null;
    this._cache = null;
  }

  CombinatorialFilter.fromSequence = function (seq, total, fov) {
    var f = new CombinatorialFilter(seq.initial, total, fov);
    f.extend(seq.events);
    return f;
  };

  CombinatorialFilter.prototype.apply = function (e) {
    this.istate.apply(e);
    this._cache = null;
  };

  CombinatorialFilter.prototype.extend = function (events) {
    for (var i = 0; i < events.length; i++) this.apply(events[i]);
  };

  CombinatorialFilter.prototype.alive = function () {
    return this.istate.alive();
  };

  /* Augmented graph of Fig. 11(d); a query group is routed through node G. */
  CombinatorialFilter.prototype._build = function (rightGroup, leftGroup) {
    var st = this.istate;
    st.flush();
    var net = new FlowNetwork();
    net.addNode(S);
    net.addNode(T);
    net.addNode(I);
    net.addEdge(T, S);
    var total = this.total ? net.addEdge(S, I, this.total[0], this.total[1]) : net.addEdge(S, I);
    var initial = new Set(st.initial);
    var group = -1;
    var lg = new Set(leftGroup || []);
    if (lg.size) {
      var kinds = new Set(Array.from(lg).map(function (v) { return initial.has(v); }));
      if (kinds.size !== 1)
        throw new RangeError("a left query set must contain only initial or only appeared shadows");
      group = net.addEdge(kinds.has(true) ? I : S, G);
    }
    var supply = new Map();
    st.left.forEach(function (b, v) {
      var parent = lg.has(v) ? G : initial.has(v) ? I : S;
      supply.set(v, net.addEdge(parent, L(v), b[0], b[1]));
    });
    var rg = new Set(rightGroup || []);
    if (rg.size) group = net.addEdge(G, T);
    var demand = new Map(), inner = [];
    [["gone", st.disappearedReach], ["alive", st.reach]].forEach(function (ts) {
      var tag = ts[0];
      ts[1].forEach(function (reach, v) {
        var rn = R(tag, v);
        st.sortedLeft(reach).forEach(function (u) {
          inner.push([u, tag, v, net.addEdge(L(u), rn)]);
        });
        var b = tag === "gone" ? st.disappeared.get(v) : [0, INF];
        demand.set(rightKey(tag, v), net.addEdge(rn, tag === "alive" && rg.has(v) ? G : T, b[0], b[1]));
      });
    });
    return { net: net, supply: supply, demand: demand, inner: inner, total: total, group: group };
  };

  /* The cached feasible network, or a group network warm-started from its circulation. */
  CombinatorialFilter.prototype._feasibleNet = function (rightGroup, leftGroup) {
    var rg = rightGroup && rightGroup.length ? rightGroup : null;
    var lg = leftGroup && leftGroup.length ? leftGroup : null;
    if (this._cache === null) {
      var base = this._build();
      base.net.findFeasible();
      this._cache = base;
    }
    if (!rg && !lg) return this._cache;
    var n = this._build(rg, lg);
    if (n.group < 0) throw new Error("group query without a group edge");
    // Same feasible set as the cached network: reuse its circulation edge by edge.
    var c = this._cache, f = n.net.flow, cf = c.net.flow, g = 0;
    f[0] = cf[0]; // T -> S
    f[n.total] = cf[c.total];
    n.supply.forEach(function (e, v) { f[e] = cf[c.supply.get(v)]; });
    n.demand.forEach(function (e, r) { f[e] = cf[c.demand.get(r)]; });
    n.inner.forEach(function (x, i) { f[x[3]] = cf[c.inner[i][3]]; });
    if (rg) rg.forEach(function (s) { g += cf[c.demand.get(rightKey("alive", s))]; });
    else lg.forEach(function (v) { g += cf[c.supply.get(v)]; });
    f[n.group] = g;
    return n;
  };

  /* true iff the observations so far admit some target distribution. */
  CombinatorialFilter.prototype.feasible = function () {
    try {
      this._feasibleNet();
    } catch (err) {
      if (err instanceof InfeasibleError) return false;
      throw err;
    }
    return true;
  };

  CombinatorialFilter.prototype._checkAlive = function (shadows) {
    var q = uniq(shadows);
    for (var i = 0; i < q.length; i++)
      if (!this.istate.reach.has(q[i])) throw new RangeError("shadow " + q[i] + " is not alive");
    if (!q.length) throw new RangeError("empty query");
    return q;
  };

  /* Tight [lo, hi] on the total number of targets in a set of alive shadows. */
  CombinatorialFilter.prototype.bounds = function (shadows) {
    var q = this._checkAlive(typeof shadows === "number" ? [shadows] : shadows);
    var n;
    if (q.length === 1) {
      n = this._feasibleNet();
      return n.net.extremes(n.demand.get(rightKey("alive", q[0])));
    }
    n = this._feasibleNet(q);
    return n.net.extremes(n.group);
  };

  /* {label: [lo, hi]} for every alive shadow (one feasible flow, two augmentations each). */
  CombinatorialFilter.prototype.allBounds = function () {
    var n = this._feasibleNet(), out = {};
    this.istate.reach.forEach(function (_, s) {
      out[s] = n.net.extremes(n.demand.get(rightKey("alive", s)));
    });
    return out;
  };

  /* Tightened [lo, hi] for shadows at t0 and appeared shadows (Sec. V-E, Eqs. (9), (10)),
   * computed exactly as min/max flow on the supply edge. */
  CombinatorialFilter.prototype.refineInitialBounds = function (labels) {
    var vs = labels === undefined || labels === null ? this.istate.realLeft() : this._checkLeft(labels, true);
    var n = this._feasibleNet(), out = {};
    vs.forEach(function (v) { out[v] = n.net.extremes(n.supply.get(v)); });
    return out;
  };

  CombinatorialFilter.prototype._checkLeft = function (labels, allowEmpty) {
    var lg = uniq(labels);
    for (var i = 0; i < lg.length; i++)
      if (!this.istate.left.has(lg[i]))
        throw new RangeError("shadow " + lg[i] + " is not a left vertex (initial or appeared)");
    if (!lg.length && !allowEmpty) throw new RangeError("empty query");
    return lg;
  };

  /* Tight [lo, hi] on the total number of targets that were in a set of left shadows. */
  CombinatorialFilter.prototype.leftBounds = function (labels) {
    var lg = this._checkLeft(labels);
    if (lg.length === 1) return this.refineInitialBounds(lg)[lg[0]];
    var n = this._feasibleNet(null, lg);
    return n.net.extremes(n.group);
  };

  /* Tight [lo, hi] on the number of targets in the shadows at t0 (counting, Sec. V-E). */
  CombinatorialFilter.prototype.initialTotalBounds = function () {
    var n = this._feasibleNet();
    return n.net.extremes(n.total);
  };

  /* A feasible distribution attaining the "max" or "min" bound of shadows:
   * {value, supply: Map(left -> n), flow: [[left, tag, right, n], ...] (nonzero only)}. */
  CombinatorialFilter.prototype.witness = function (shadows, sense) {
    var q = this._checkAlive(shadows);
    var n = this._build(q);
    n.net.findFeasible();
    var value = sense === "min" ? n.net.minimize(n.group) : n.net.maximize(n.group);
    if (value === INF) throw new RangeError("unbounded: no finite witness");
    var f = n.net.flow, supply = new Map();
    n.supply.forEach(function (e, v) { supply.set(v, f[e]); });
    var flow = n.inner.filter(function (x) { return f[x[3]]; }).map(function (x) {
      return [x[0], x[1], x[2], f[x[3]]];
    });
    return { value: value, supply: supply, flow: flow };
  };

  /* Bounds via the literal max-flow recipe of Sec. V-C, Eqs. (7) and (8) (see
   * shadowinfo CombinatorialFilter.bounds_paper; corrected=true fixes Eq. (8)). */
  CombinatorialFilter.prototype.boundsPaper = function (shadows, corrected) {
    var q = new Set(this._checkAlive(shadows));
    var st = this.istate;
    st.flush();

    function run(leftCap, goneCap, aliveCap) {
      var net = new FlowNetwork(), out = new Map();
      net.addNode(S);
      net.addNode(T);
      st.left.forEach(function (b, v) { net.addEdge(S, L(v), 0, leftCap(b)); });
      [["gone", st.disappearedReach], ["alive", st.reach]].forEach(function (ts) {
        var tag = ts[0];
        ts[1].forEach(function (reach, v) {
          var rn = R(tag, v);
          st.sortedLeft(reach).forEach(function (u) { net.addEdge(L(u), rn); });
          var cap = tag === "gone" ? goneCap(st.disappeared.get(v)) : aliveCap(v);
          out.set(rightKey(tag, v), net.addEdge(rn, T, 0, cap));
        });
      });
      var value = net.maxFlow(S, T), f = new Map();
      out.forEach(function (e, r) { f.set(r, net.flow[e]); });
      return [value, f];
    }

    var lo = function (b) { return b[0]; }, hi = function (b) { return b[1]; };
    var res = run(hi, lo, function (s) { return q.has(s) ? INF : 0; });
    var lsum = 0;
    st.disappeared.forEach(function (b) { lsum += b[0]; });
    var upper;
    if (res[0] === INF) {
      upper = INF;
    } else {
      upper = -lsum;
      q.forEach(function (s) { upper += res[1].get(rightKey("alive", s)); });
      st.disappeared.forEach(function (_, v) { upper += res[1].get(rightKey("gone", v)); });
    }
    res = run(lo, corrected ? hi : lo, function (s) { return q.has(s) ? 0 : INF; });
    var lower = -res[0];
    st.left.forEach(function (b) { lower += b[0]; });
    return [lower, upper];
  };

  /* Counting: unknown number of targets, every shadow at t0 starts at [0, inf). */
  function countingFilter(shadows, fov) {
    var init = {};
    shadows.forEach(function (s) { init[s] = [0, INF]; });
    return new CombinatorialFilter(init, null, fov);
  }

  /* Passive pursuit-evasion: one evader hidden in one of the shadows at t0. */
  function pursuitEvasionFilter(shadows, fov) {
    var init = {};
    shadows.forEach(function (s) { init[s] = [0, 1]; });
    return new CombinatorialFilter(init, [1, 1], fov);
  }

  /* "clear" ([0,0]), "evader" ([1,1]) or "contaminated" ([0,1]) per alive shadow. */
  function evaderStatus(filt) {
    var b = filt.allBounds(), out = {};
    for (var s in b) out[s] = b[s][1] === 0 ? "clear" : b[s][0] >= 1 ? "evader" : "contaminated";
    return out;
  }

  SI.CombinatorialFilter = CombinatorialFilter;
  SI.countingFilter = countingFilter;
  SI.pursuitEvasionFilter = pursuitEvasionFilter;
  SI.evaderStatus = evaderStatus;
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
