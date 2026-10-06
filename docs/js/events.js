/* Critical events over shadows (T-RO 2012, Secs. II-III).
 * Port of shadowinfo/events.py.  Events are plain objects in the JSON form of
 * DESIGN §2, e.g. {type: "split", s: 1, a: 2, b: 3}; an unbounded hi is
 * Infinity in memory and null in JSON.  Shadow labels are integers, never reused. */
(function (SI) {
  "use strict";

  var FIELDS = {
    appear: ["s", "lo", "hi"],
    disappear: ["s", "lo", "hi"],
    split: ["s", "a", "b"],
    merge: ["a", "b", "s"],
    enter: ["s", "k"],
    exit: ["s", "k"],
  };
  var COMPONENT = { appear: 1, disappear: 1, split: 1, merge: 1 };
  var FOV = { enter: 1, exit: 1 };

  class InvalidSequenceError extends Error {}
  InvalidSequenceError.prototype.name = "InvalidSequenceError";

  var ev = {
    appear: function (s, lo, hi) {
      return { type: "appear", s: s, lo: lo || 0, hi: hi === undefined ? lo || 0 : hi };
    },
    disappear: function (s, lo, hi) {
      return { type: "disappear", s: s, lo: lo || 0, hi: hi === undefined ? lo || 0 : hi };
    },
    split: function (s, a, b) {
      return { type: "split", s: s, a: a, b: b };
    },
    merge: function (a, b, s) {
      return { type: "merge", a: a, b: b, s: s };
    },
    enter: function (s, k) {
      return { type: "enter", s: s, k: k === undefined ? 1 : k };
    },
    exit: function (s, k) {
      return { type: "exit", s: s, k: k === undefined ? 1 : k };
    },
  };

  function hiToJson(hi) {
    return hi === Infinity ? null : hi;
  }

  function hiFromJson(hi) {
    return hi === null || hi === undefined || hi === Infinity ? Infinity : hi;
  }

  function boundFromJson(b) {
    return [b[0], hiFromJson(b[1])];
  }

  /* Event in JSON form -> event object (null hi -> Infinity, defaults filled, extra keys dropped). */
  function eventFromDict(d) {
    var f = FIELDS[d.type];
    if (!f) throw new InvalidSequenceError("unknown event type " + JSON.stringify(d.type));
    var e = { type: d.type };
    for (var i = 0; i < f.length; i++) e[f[i]] = d[f[i]];
    if (d.type === "appear" || d.type === "disappear") {
      if (e.lo === undefined) e.lo = 0;
      e.hi = "hi" in d ? hiFromJson(d.hi) : 0;
    } else if (FOV[d.type] && e.k === undefined) {
      e.k = 1;
    }
    return e;
  }

  function eventToDict(e) {
    var f = FIELDS[e.type];
    var d = { type: e.type };
    for (var i = 0; i < f.length; i++) d[f[i]] = f[i] === "hi" ? hiToJson(e.hi) : e[f[i]];
    return d;
  }

  /* [consumed, created] labels; FOV events neither consume nor create labels. */
  function eventLabels(e) {
    switch (e.type) {
      case "appear":
        return [[], [e.s]];
      case "disappear":
        return [[e.s], []];
      case "split":
        return [[e.s], [e.a, e.b]];
      case "merge":
        return [[e.a, e.b], [e.s]];
      default:
        return [[], []];
    }
  }

  function isComponent(e) {
    return COMPONENT[e.type] === 1;
  }

  function isFov(e) {
    return FOV[e.type] === 1;
  }

  function describe(e) {
    return JSON.stringify(eventToDict(e));
  }

  /* Initial condition {label: [lo, hi]} plus an ordered, validated event list. */
  function ShadowSequence(initial, events) {
    this.initial = {};
    for (var k in initial) this.initial[+k] = boundFromJson(initial[k]);
    this.events = (events || []).slice();
    this.validate();
  }

  ShadowSequence.prototype.validate = function () {
    var alive = new Set(), seen = new Set(), k;
    for (k in this.initial) {
      var b = this.initial[k];
      if (b[0] < 0 || b[1] < b[0]) throw new InvalidSequenceError("bad initial bound " + JSON.stringify(b));
      alive.add(+k);
      seen.add(+k);
    }
    for (var i = 0; i < this.events.length; i++) {
      var e = this.events[i], where = "event " + i + " " + describe(e);
      if (!FIELDS[e.type]) throw new InvalidSequenceError(where + ": unknown type");
      if (isFov(e)) {
        if (!alive.has(e.s)) throw new InvalidSequenceError(where + ": shadow " + e.s + " not alive");
        if (!(e.k >= 1)) throw new InvalidSequenceError(where + ": k must be >= 1");
        continue;
      }
      if ((e.type === "appear" || e.type === "disappear") && (e.lo < 0 || e.hi < e.lo))
        throw new InvalidSequenceError(where + ": bad bound");
      var cc = eventLabels(e), consumed = cc[0], created = cc[1];
      if (new Set(consumed).size !== consumed.length || new Set(created).size !== created.length)
        throw new InvalidSequenceError(where + ": repeated label");
      consumed.forEach(function (s) {
        if (!alive.has(s)) throw new InvalidSequenceError(where + ": shadow " + s + " not alive");
      });
      created.forEach(function (s) {
        if (seen.has(s)) throw new InvalidSequenceError(where + ": label " + s + " reused");
      });
      consumed.forEach(function (s) { alive.delete(s); });
      created.forEach(function (s) { alive.add(s); seen.add(s); });
    }
  };

  ShadowSequence.prototype.aliveAtEnd = function () {
    var alive = Object.keys(this.initial).map(Number);
    this.events.forEach(function (e) {
      var cc = eventLabels(e);
      alive = alive.filter(function (s) { return cc[0].indexOf(s) < 0; }).concat(cc[1]);
    });
    return alive.sort(function (a, b) { return a - b; });
  };

  ShadowSequence.prototype.allLabels = function () {
    var labels = new Set(Object.keys(this.initial).map(Number));
    this.events.forEach(function (e) {
      eventLabels(e)[1].forEach(function (s) { labels.add(s); });
    });
    return Array.from(labels).sort(function (a, b) { return a - b; });
  };

  ShadowSequence.prototype.append = function (e) {
    this.events.push(e);
    try {
      this.validate();
    } catch (err) {
      this.events.pop();
      throw err;
    }
  };

  ShadowSequence.prototype.toDict = function () {
    var initial = {};
    for (var k in this.initial) initial[k] = [this.initial[k][0], hiToJson(this.initial[k][1])];
    return { initial: initial, events: this.events.map(eventToDict) };
  };

  ShadowSequence.prototype.toJSON = ShadowSequence.prototype.toDict;

  ShadowSequence.fromDict = function (d) {
    return new ShadowSequence(d.initial, d.events.map(eventFromDict));
  };

  ShadowSequence.fromJSON = function (text) {
    return ShadowSequence.fromDict(JSON.parse(text));
  };

  SI.ev = ev;
  SI.InvalidSequenceError = InvalidSequenceError;
  SI.ShadowSequence = ShadowSequence;
  SI.events = {
    FIELDS: FIELDS,
    hiToJson: hiToJson,
    hiFromJson: hiFromJson,
    boundFromJson: boundFromJson,
    eventFromDict: eventFromDict,
    eventToDict: eventToDict,
    eventLabels: eventLabels,
    isComponent: isComponent,
    isFov: isFov,
    describe: describe,
  };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
