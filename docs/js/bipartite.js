/* Incremental bipartite I-state (T-RO 2012, Sec. V-B, Fig. 9) and FOV batching (Sec. V-D).
 * Port of shadowinfo/bipartite.py and shadowinfo/fov.py.
 *
 * Left vertices: shadows at t0, appeared shadows and FOV-enter pseudo-shadows;
 * right vertices: disappeared shadows, FOV-exit pseudo-shadows and alive shadows.
 * A real shadow is its integer label; a pseudo-shadow is the string "+s#n"
 * (enter) or "-s#n" (exit), so it can never equal a label.  Reach sets are
 * treated as immutable Sets and shared between split children. */
(function (SI) {
  "use strict";

  if (!SI.events && typeof require === "function") require("./events.js");

  var InvalidSequenceError = SI.InvalidSequenceError;
  var describe = SI.events.describe;

  /* Running (dTot, dMin) counter for the FOV events of one shadow (Sec. V-D). */
  function FovBatch() {
    this.dTot = 0;
    this.dMin = 0;
  }

  FovBatch.prototype.add = function (e) {
    if (e.type === "enter") {
      this.dTot += e.k;
    } else if (e.type === "exit") {
      this.dTot -= e.k;
      this.dMin = Math.min(this.dMin, this.dTot);
    } else {
      throw new TypeError("not a FOV event: " + describe(e));
    }
  };

  /* At most one batch exit then one batch enter, equivalent to the run (four cases). */
  FovBatch.prototype.batchEvents = function (s) {
    var out = [];
    if (this.dMin < 0) out.push({ type: "exit", s: s, k: -this.dMin });
    if (this.dTot - this.dMin > 0) out.push({ type: "enter", s: s, k: this.dTot - this.dMin });
    return out;
  };

  /* Replace every run of FOV events per shadow by its batch form. */
  function batchFov(events) {
    var pending = new Map(), out = [];
    function flush(s) {
      var b = pending.get(s);
      if (b === undefined) return;
      pending.delete(s);
      Array.prototype.push.apply(out, b.batchEvents(s));
    }
    events.forEach(function (e) {
      if (e.type === "enter" || e.type === "exit") {
        if (!pending.has(e.s)) pending.set(e.s, new FovBatch());
        pending.get(e.s).add(e);
        return;
      }
      if (e.type === "split" || e.type === "disappear") flush(e.s);
      else if (e.type === "merge") {
        flush(e.a);
        flush(e.b);
      }
      out.push(e);
    });
    Array.from(pending.keys()).forEach(flush);
    return out;
  }

  function isPseudo(v) {
    return typeof v === "string";
  }

  function vertexName(v) {
    return String(v);
  }

  function checkBound(lo, hi, where) {
    if (lo < 0 || hi < lo) throw new InvalidSequenceError(where + ": bad bound [" + lo + ", " + hi + "]");
  }

  /* Bipartite I-state maintained one event at a time (Fig. 9).
   * initial: {label: [lo, hi]}; fov: "batch" (default) or "naive". */
  function BipartiteIState(initial, fov) {
    fov = fov || "batch";
    if (fov !== "batch" && fov !== "naive") throw new RangeError("fov must be 'batch' or 'naive'");
    this.fov = fov;
    this.left = new Map();
    this.leftOrder = new Map();
    this.initial = [];
    this.disappeared = new Map();
    this.disappearedReach = new Map();
    this.reach = new Map();
    this.pending = new Map();
    this.seen = new Set();
    this.nPseudo = 0;
    var labels = initial instanceof Map ? Array.from(initial.keys()) : Object.keys(initial).map(Number);
    for (var i = 0; i < labels.length; i++) {
      var s = labels[i], b = initial instanceof Map ? initial.get(s) : initial[s];
      this._appear(s, b[0], SI.events.hiFromJson(b[1]));
      this.initial.push(s);
    }
  }

  BipartiteIState.fromSequence = function (seq, fov) {
    var st = new BipartiteIState(seq.initial, fov);
    st.extend(seq.events);
    return st;
  };

  BipartiteIState.prototype.apply = function (e) {
    var where = describe(e), r;
    switch (e.type) {
      case "enter":
      case "exit":
        this._requireAlive(e.s, where);
        if (!(e.k >= 1)) throw new InvalidSequenceError(where + ": k must be >= 1");
        if (this.fov === "naive") {
          if (e.type === "enter") this._enter(e.s, e.k);
          else this._exit(e.s, e.k);
        } else {
          if (!this.pending.has(e.s)) this.pending.set(e.s, new FovBatch());
          this.pending.get(e.s).add(e);
        }
        break;
      case "appear":
        this._appear(e.s, e.lo, e.hi);
        break;
      case "disappear":
        this._requireAlive(e.s, where);
        checkBound(e.lo, e.hi, where);
        this.flush(e.s);
        this.disappeared.set(e.s, [e.lo, e.hi]);
        this.disappearedReach.set(e.s, this.reach.get(e.s));
        this.reach.delete(e.s);
        break;
      case "split":
        this._requireAlive(e.s, where);
        this._requireFresh([e.a, e.b], where);
        if (e.a === e.b) throw new InvalidSequenceError(where + ": repeated label");
        this.flush(e.s);
        r = this.reach.get(e.s);
        this.reach.delete(e.s);
        this.reach.set(e.a, r);
        this.reach.set(e.b, r);
        this.seen.add(e.a);
        this.seen.add(e.b);
        break;
      case "merge":
        if (e.a === e.b) throw new InvalidSequenceError(where + ": repeated label");
        this._requireAlive(e.a, where);
        this._requireAlive(e.b, where);
        this._requireFresh([e.s], where);
        this.flush(e.a);
        this.flush(e.b);
        r = new Set(this.reach.get(e.a));
        this.reach.get(e.b).forEach(function (v) { r.add(v); });
        this.reach.delete(e.a);
        this.reach.delete(e.b);
        this.reach.set(e.s, r);
        this.seen.add(e.s);
        break;
      default:
        throw new TypeError("unknown event " + where);
    }
  };

  BipartiteIState.prototype.extend = function (events) {
    for (var i = 0; i < events.length; i++) this.apply(events[i]);
  };

  /* Turn pending FOV batches into pseudo-shadows (for s, or for all shadows). */
  BipartiteIState.prototype.flush = function (s) {
    var targets = s === undefined || s === null ? Array.from(this.pending.keys()) : [s];
    for (var i = 0; i < targets.length; i++) {
      var t = targets[i], b = this.pending.get(t);
      if (b === undefined) continue;
      this.pending.delete(t);
      var batch = b.batchEvents(t);
      for (var j = 0; j < batch.length; j++) {
        if (batch[j].type === "enter") this._enter(t, batch[j].k);
        else this._exit(t, batch[j].k);
      }
    }
  };

  BipartiteIState.prototype._appear = function (s, lo, hi) {
    checkBound(lo, hi, "shadow " + s);
    this._requireFresh([s], "shadow " + s);
    this.seen.add(s);
    this._addLeft(s, [lo, hi]);
    this.reach.set(s, new Set([s]));
  };

  BipartiteIState.prototype._enter = function (s, k) {
    var p = this._pseudo("+", s);
    this._addLeft(p, [k, k]);
    var r = new Set(this.reach.get(s));
    r.add(p);
    this.reach.set(s, r);
  };

  BipartiteIState.prototype._exit = function (s, k) {
    var p = this._pseudo("-", s);
    this.disappeared.set(p, [k, k]);
    this.disappearedReach.set(p, this.reach.get(s));
  };

  BipartiteIState.prototype._pseudo = function (sign, s) {
    this.nPseudo += 1;
    return sign + s + "#" + this.nPseudo;
  };

  BipartiteIState.prototype._addLeft = function (v, b) {
    this.leftOrder.set(v, this.leftOrder.size);
    this.left.set(v, b);
  };

  BipartiteIState.prototype._requireAlive = function (s, where) {
    if (!this.reach.has(s)) throw new InvalidSequenceError(where + ": shadow " + s + " not alive");
  };

  BipartiteIState.prototype._requireFresh = function (labels, where) {
    for (var i = 0; i < labels.length; i++)
      if (this.seen.has(labels[i])) throw new InvalidSequenceError(where + ": label " + labels[i] + " reused");
  };

  BipartiteIState.prototype.alive = function () {
    return Array.from(this.reach.keys());
  };

  BipartiteIState.prototype.sortedLeft = function (vertices) {
    var order = this.leftOrder;
    return Array.from(vertices).sort(function (a, b) { return order.get(a) - order.get(b); });
  };

  BipartiteIState.prototype.realLeft = function () {
    return Array.from(this.left.keys()).filter(function (v) { return !isPseudo(v); });
  };

  /* Right vertices as [tag, vertex, lo, hi] with tag "gone" or "alive" (pending batches not flushed). */
  BipartiteIState.prototype.right = function () {
    var out = [];
    this.disappeared.forEach(function (b, v) { out.push(["gone", v, b[0], b[1]]); });
    this.reach.forEach(function (_, s) { out.push(["alive", s, 0, Infinity]); });
    return out;
  };

  /* Edges [left, tag, right] in the order of shadowinfo BipartiteIState.edges(). */
  BipartiteIState.prototype.edges = function () {
    var out = [], self = this;
    [["gone", this.disappearedReach], ["alive", this.reach]].forEach(function (ts) {
      ts[1].forEach(function (r, v) {
        self.sortedLeft(r).forEach(function (u) { out.push([u, ts[0], v]); });
      });
    });
    return out;
  };

  /* JSON snapshot (Fig. 11(c) style), same format as Python to_dict(); flushes pending batches. */
  BipartiteIState.prototype.toDict = function () {
    this.flush();
    var hj = SI.events.hiToJson;
    var left = [];
    this.left.forEach(function (b, v) { left.push([vertexName(v), b[0], hj(b[1])]); });
    return {
      left: left,
      right: this.right().map(function (r) { return [vertexName(r[1]), r[2], hj(r[3]), r[0]]; }),
      edges: this.edges().map(function (e) { return [vertexName(e[0]), vertexName(e[2]), e[1]]; }),
    };
  };

  SI.FovBatch = FovBatch;
  SI.batchFov = batchFov;
  SI.BipartiteIState = BipartiteIState;
  SI.bipartite = { isPseudo: isPseudo, vertexName: vertexName };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
