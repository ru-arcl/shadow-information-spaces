/* Edmonds–Karp max-flow and flows with lower bounds (T-RO 2012, Sec. V-C).
 * Port of shadowinfo/maxflow.py.  Edge e owns arcs 2e (forward, residual
 * hi - f) and 2e + 1 (backward, residual f - lo).  Capacities are integers or
 * Infinity; arcs are scanned in insertion order, so results are deterministic
 * and identical to the Python implementation. */
(function (SI) {
  "use strict";

  var INF = Infinity;

  class InfeasibleError extends Error {}
  InfeasibleError.prototype.name = "InfeasibleError";

  function FlowNetwork() {
    this.keys = [];
    this.index = new Map();
    this.adj = [];
    this.head = [];
    this.lo = [];
    this.hi = [];
    this.flow = [];
    this._uid = 0;
  }

  FlowNetwork.prototype.addNode = function (key) {
    var i = this.index.get(key);
    if (i !== undefined) return i;
    i = this.keys.length;
    this.index.set(key, i);
    this.keys.push(key);
    this.adj.push([]);
    return i;
  };

  /* Edge u -> v with lo <= f <= hi; returns the edge id. */
  FlowNetwork.prototype.addEdge = function (u, v, lo, hi) {
    if (lo === undefined) lo = 0;
    if (hi === undefined || hi === null) hi = INF;
    if (lo < 0 || lo === INF) throw new RangeError("bad lower bound " + lo);
    var iu = this.addNode(u), iv = this.addNode(v);
    var e = this.lo.length;
    this.lo.push(lo);
    this.hi.push(hi);
    this.flow.push(0);
    this.head.push(iv, iu);
    this.adj[iu].push(2 * e);
    this.adj[iv].push(2 * e + 1);
    return e;
  };

  Object.defineProperty(FlowNetwork.prototype, "numEdges", {
    get: function () { return this.lo.length; },
  });

  FlowNetwork.prototype.tailOf = function (e) {
    return this.head[2 * e + 1];
  };

  FlowNetwork.prototype.headOf = function (e) {
    return this.head[2 * e];
  };

  FlowNetwork.prototype.snapshot = function () {
    return this.flow.slice();
  };

  FlowNetwork.prototype.restore = function (flows) {
    for (var i = 0; i < flows.length; i++) this.flow[i] = flows[i];
  };

  FlowNetwork.prototype._residual = function (a) {
    var e = a >> 1;
    return a & 1 ? this.flow[e] - this.lo[e] : this.hi[e] - this.flow[e];
  };

  /* Shortest residual path s -> t (BFS) avoiding edge skip: [arcs t..s, bottleneck] or null. */
  FlowNetwork.prototype._augmentingPath = function (s, t, skip) {
    var n = this.keys.length;
    var parent = new Int32Array(n).fill(-1);
    parent[s] = -2;
    var queue = new Int32Array(n), qh = 0, qt = 0;
    queue[qt++] = s;
    var head = this.head, flow = this.flow, lo = this.lo, hi = this.hi, adj = this.adj;
    while (qh < qt) {
      var u = queue[qh++], arcs = adj[u];
      for (var j = 0; j < arcs.length; j++) {
        var a = arcs[j], e = a >> 1;
        if (e === skip) continue;
        var v = head[a];
        if (parent[v] !== -1) continue;
        if ((a & 1 ? flow[e] - lo[e] : hi[e] - flow[e]) <= 0) continue;
        parent[v] = a;
        if (v === t) {
          var path = [], bottleneck = INF;
          while (v !== s) {
            a = parent[v];
            path.push(a);
            bottleneck = Math.min(bottleneck, this._residual(a));
            v = head[a ^ 1];
          }
          return [path, bottleneck];
        }
        queue[qt++] = v;
      }
    }
    return null;
  };

  FlowNetwork.prototype._push = function (path, amount) {
    for (var i = 0; i < path.length; i++) {
      var a = path[i];
      if (a & 1) this.flow[a >> 1] -= amount;
      else this.flow[a >> 1] += amount;
    }
  };

  /* Edmonds–Karp from the current flow; returns the amount pushed (Infinity if unbounded,
   * in which case the flow is left unchanged when limit is infinite: finite pushes made
   * before the infinite path was found are undone). */
  FlowNetwork.prototype._augment = function (s, t, limit, skip) {
    if (limit === undefined) limit = INF;
    if (skip === undefined) skip = -1;
    var total = 0, pushed = [];
    while (total < limit) {
      var found = this._augmentingPath(s, t, skip);
      if (found === null) break;
      var amount = Math.min(found[1], limit - total);
      if (amount === INF) {
        for (var i = pushed.length - 1; i >= 0; i--) this._push(pushed[i][0], -pushed[i][1]);
        return INF;
      }
      this._push(found[0], amount);
      pushed.push([found[0], amount]);
      total += amount;
    }
    return total;
  };

  /* Edmonds–Karp max-flow value from s to t (lower bounds must be 0). */
  FlowNetwork.prototype.maxFlow = function (s, t) {
    return this._augment(this.index.get(s), this.index.get(t));
  };

  /* Circulation with lo <= f <= hi on every edge via the super source/sink reduction;
   * throws InfeasibleError if none exists. */
  FlowNetwork.prototype.findFeasible = function () {
    var m = this.numEdges, e;
    for (e = 0; e < m; e++) {
      if (this.lo[e] > this.hi[e])
        throw new InfeasibleError("edge " + this.keys[this.tailOf(e)] + "->" + this.keys[this.headOf(e)] +
                                  ": lo " + this.lo[e] + " > hi " + this.hi[e]);
    }
    for (e = 0; e < m; e++) this.flow[e] = this.lo[e];
    var nNodes = this.keys.length;
    var excess = new Array(nNodes).fill(0);
    for (e = 0; e < m; e++) {
      excess[this.headOf(e)] += this.lo[e];
      excess[this.tailOf(e)] -= this.lo[e];
    }
    var tag = "\u0000super" + this._uid++;
    var ss = this.addNode(tag + "s"), tt = this.addNode(tag + "t");
    var need = 0;
    for (var v = 0; v < nNodes; v++) {
      if (excess[v] > 0) {
        this.addEdge(this.keys[ss], this.keys[v], 0, excess[v]);
        need += excess[v];
      } else if (excess[v] < 0) {
        this.addEdge(this.keys[v], this.keys[tt], 0, -excess[v]);
      }
    }
    var got = this._augment(ss, tt);
    this._truncate(nNodes, m);
    if (got < need) throw new InfeasibleError("no flow satisfies the bounds");
  };

  FlowNetwork.prototype._truncate = function (nNodes, nEdges) {
    for (var i = nNodes; i < this.keys.length; i++) this.index.delete(this.keys[i]);
    this.keys.length = nNodes;
    this.adj.length = nNodes;
    this.lo.length = this.hi.length = this.flow.length = nEdges;
    this.head.length = 2 * nEdges;
    for (var u = 0; u < nNodes; u++) {
      var lst = this.adj[u];
      while (lst.length && lst[lst.length - 1] >> 1 >= nEdges) lst.pop();
    }
  };

  /* Raise the flow on edge e as far as feasibility allows; returns f(e) (Infinity if
   * unbounded, flow unchanged).  Requires a feasible circulation. */
  FlowNetwork.prototype.maximize = function (e) {
    var x = this.tailOf(e), y = this.headOf(e);
    var extra = this._augment(y, x, this.hi[e] - this.flow[e], e);
    if (extra === INF) return INF;
    this.flow[e] += extra;
    return this.flow[e];
  };

  /* Lower the flow on edge e as far as feasibility allows; returns f(e). */
  FlowNetwork.prototype.minimize = function (e) {
    var x = this.tailOf(e), y = this.headOf(e);
    this.flow[e] -= this._augment(x, y, this.flow[e] - this.lo[e], e);
    return this.flow[e];
  };

  /* [min f(e), max f(e)] over feasible circulations; the flow is restored. */
  FlowNetwork.prototype.extremes = function (e) {
    var saved = this.snapshot();
    var lo = this.minimize(e);
    this.restore(saved);
    var hi = this.maximize(e);
    this.restore(saved);
    return [lo, hi];
  };

  FlowNetwork.prototype.isCirculation = function () {
    var bal = new Array(this.keys.length).fill(0);
    for (var e = 0; e < this.numEdges; e++) {
      if (!(this.lo[e] <= this.flow[e] && this.flow[e] <= this.hi[e])) return false;
      bal[this.headOf(e)] += this.flow[e];
      bal[this.tailOf(e)] -= this.flow[e];
    }
    return bal.every(function (b) { return b === 0; });
  };

  SI.INF = INF;
  SI.InfeasibleError = InfeasibleError;
  SI.FlowNetwork = FlowNetwork;
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
