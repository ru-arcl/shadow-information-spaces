/* Polygon engine: geometry, critical lines (cuts), visibility, gap tracking and component events.
 *
 * Bit-for-bit port of shadowinfo/polygon/{geometry,cuts,visibility,gaps,events}.py, which port the
 * original Java implementation (docs/notes/original_java.md §2.3-§2.7; T-RO 2012 Secs. II-A, V-C).
 * Gap tracking is the Python compat="fixed" mode: Java's rules where they are consistent, the hidden
 * boundary chains elsewhere, which reproduces the Java gap history on every run where Java does not
 * throw; {javaMatching: false} is the opt-in chain relabelling the simulator uses.  physicalGaps and
 * visibilityPolygon are the fixed ones too (they refuse a point outside the polygon and repair Java's
 * degenerate scan on lines through two vertices); {java: true} gives Java's raw scan.  Polygon refuses
 * clockwise, non-simple or non-finite vertices unless {validate: false}.
 *
 * Every formula keeps the Python (= Java) operation order on doubles, so cuts, critical points,
 * samples, physical gaps and gap histories are the same numbers as in Python; only + - * / sqrt and
 * comparisons are used (no trig).  Polygons are y-up and counter-clockwise (the .dat frame), never
 * flipped.  A point is [x, y], a line or segment [x1, y1, x2, y2]; a physical gap is
 * {startEdge, endEdge, startPoint, endPoint} with exactly one of the two points set; component
 * events are plain objects in the JSON form of DESIGN §2. */
(function (SI) {
  "use strict";

  var EPSILON = 5e-5; /* Algorithm.epsilon */
  var PERTURB = 0.005; /* Algorithm.purturb */
  var CUT = {
    GENERAL_INFLECTION: "GENERAL_INFLECTION",
    NONGENERAL_INFLECTION: "NONGENERAL_INFLECTION",
    SINGLETANGENT: "SINGLETANGENT",
    BITANGENT: "BITANGENT",
    NONE: "NONE",
  };

  function imod(a, n) {
    var r = a % n;
    return r < 0 ? r + n : r;
  }

  /* Python float ``x % y`` (CPython float_rem: the result has the sign of y). */
  function pyMod(x, y) {
    var m = x % y;
    if (m !== 0) {
      if (y < 0 !== m < 0) m += y;
    } else {
      m = y < 0 ? -0 : 0;
    }
    return m;
  }

  /* ---- JDK / Algorithm primitives (geometry.py) --------------------------- */

  /* Line2D.relativeCCW: +1 when P is right of A->B (y-up), -1 left; collinear: -1 before A, +1 past B, 0 on. */
  function relativeCCW(x1, y1, x2, y2, px, py) {
    x2 -= x1;
    y2 -= y1;
    px -= x1;
    py -= y1;
    var ccw = px * y2 - py * x2;
    if (ccw === 0) {
      ccw = px * x2 + py * y2;
      if (ccw > 0) {
        px -= x2;
        py -= y2;
        ccw = px * x2 + py * y2;
        if (ccw < 0) ccw = 0;
      }
    }
    return ccw < 0 ? -1 : ccw > 0 ? 1 : 0;
  }

  function rccwLine(l, p) {
    return relativeCCW(l[0], l[1], l[2], l[3], p[0], p[1]);
  }

  /* Line2D.ptSegDistSq (JDK projection formula). */
  function ptSegDistSq(x1, y1, x2, y2, px, py) {
    x2 -= x1;
    y2 -= y1;
    px -= x1;
    py -= y1;
    var dot = px * x2 + py * y2, proj;
    if (dot <= 0) {
      proj = 0;
    } else {
      px = x2 - px;
      py = y2 - py;
      dot = px * x2 + py * y2;
      proj = dot <= 0 ? 0 : (dot * dot) / (x2 * x2 + y2 * y2);
    }
    var len = px * px + py * py - proj;
    return len < 0 ? 0 : len;
  }

  function ptSegDist(x1, y1, x2, y2, px, py) {
    return Math.sqrt(ptSegDistSq(x1, y1, x2, y2, px, py));
  }

  /* Point2D.distance: sqrt(dx*dx + dy*dy), never hypot. */
  function distance(ax, ay, bx, by) {
    var dx = bx - ax, dy = by - ay;
    return Math.sqrt(dx * dx + dy * dy);
  }

  function distanceSq(ax, ay, bx, by) {
    var dx = bx - ax, dy = by - ay;
    return dx * dx + dy * dy;
  }

  /* Algorithm.getPointOnLineFromX: (x, k*x + b); for an eps-vertical line (x, 0) or null. */
  function pointOnLineFromX(x1, y1, x2, y2, x, eps) {
    if (Math.abs(x1 - x2) <= eps) return x !== x1 ? [x, 0] : null;
    var k = (y2 - y1) / (x2 - x1);
    var b = y1 - x1 * k;
    return [x, k * x + b];
  }

  /* Algorithm.getPointOnLineFromY: ((y - b)/k, y); for an eps-horizontal line (0, y) or null. */
  function pointOnLineFromY(x1, y1, x2, y2, y, eps) {
    if (Math.abs(y1 - y2) <= eps) return y !== y1 ? [0, y] : null;
    var k = Math.abs(x1 - x2) <= eps ? Infinity : (y2 - y1) / (x2 - x1);
    var b = k === Infinity ? k : y1 - x1 * k;
    return [(y - b) / k, y];
  }

  /* Algorithm.getIntersect: the infinite lines, slope-intercept form (line 2's intercept from its P2). */
  function getIntersect(x11, y11, x12, y12, x21, y21, x22, y22, eps) {
    var v1 = Math.abs(x12 - x11) <= eps, v2 = Math.abs(x22 - x21) <= eps;
    if (v1 || v2) {
      if (!v1) return pointOnLineFromX(x11, y11, x12, y12, x21, eps);
      if (!v2) return pointOnLineFromX(x21, y21, x22, y22, x11, eps);
      return null;
    }
    var k1 = (y12 - y11) / (x12 - x11);
    var k2 = (y22 - y21) / (x22 - x21);
    var i1 = y11 - x11 * k1;
    var i2 = y22 - x22 * k2;
    if (Math.abs(k1 - k2) <= eps) return null;
    var x = (i2 - i1) / (k1 - k2);
    return [x, x * k1 + i1];
  }

  function onExtension(x1, y1, x2, y2, px, py) {
    return (px - x2) * (x2 - x1) > 0 || (py - y2) * (y2 - y1) > 0;
  }

  function onReverseExtension(x1, y1, x2, y2, px, py) {
    return (px - x1) * (x1 - x2) > 0 || (py - y1) * (y1 - y2) > 0;
  }

  /* Algorithm.getSlope */
  function getSlope(l, eps) {
    if (Math.abs(l[0] - l[2]) <= eps) return Infinity;
    return (l[3] - l[1]) / (l[2] - l[0]);
  }

  /* Algorithm.purturbPointAlongSeg (compat="fixed": the x step at slope exactly +-1, quirk B1). */
  function purturbPointAlongSeg(p, seg, amount, far, eps) {
    amount = amount === undefined ? PERTURB : amount;
    far = far === undefined ? true : far;
    eps = eps === undefined ? EPSILON : eps;
    var x = p[0], y = p[1], s = getSlope(seg, eps), p1, p2;
    if (s === Infinity) {
      p1 = [x, y + amount];
      p2 = [x, y - amount];
    } else if (s === 0) {
      p1 = [x + amount, y];
      p2 = [x - amount, y];
    } else if (Math.abs(s) > 1) {
      p1 = pointOnLineFromY(seg[0], seg[1], seg[2], seg[3], y + amount, eps);
      p2 = pointOnLineFromY(seg[0], seg[1], seg[2], seg[3], y - amount, eps);
    } else {
      p1 = pointOnLineFromX(seg[0], seg[1], seg[2], seg[3], x + amount, eps);
      p2 = pointOnLineFromX(seg[0], seg[1], seg[2], seg[3], x - amount, eps);
    }
    if (p1 === null || p2 === null) throw new Error("purturbPointAlongSeg: no candidate point on " + seg);
    var d1 = distance(p1[0], p1[1], seg[0], seg[1]);
    var d2 = distance(p2[0], p2[1], seg[0], seg[1]);
    if ((d1 > d2 && far) || (d1 < d2 && !far)) return p1;
    return p2;
  }

  /* ---- Polygon ------------------------------------------------------------ */

  /* A simple polygon (drawable.Polygon): vertices y-up CCW; edges[i] = v_i -> v_{i+1 mod n};
   * contains() is Path2D.Float.contains (float32 vertices, non-zero rule, half-open in y). */
  function Polygon(vertices, name, opts) {
    var vs = vertices.map(function (v) { return [+v[0], +v[1]]; });
    if (vs.length < 3) throw new Error("a polygon needs at least 3 vertices");
    var n = vs.length, i;
    this.vertices = vs;
    this.n = n;
    this.name = name === undefined ? null : name;
    this.edges = [];
    this.vx = new Float64Array(n);
    this.vy = new Float64Array(n);
    this.fx = new Float64Array(n);
    this.fy = new Float64Array(n);
    for (i = 0; i < n; i++) {
      var a = vs[i], b = vs[(i + 1) % n];
      this.edges.push([a[0], a[1], b[0], b[1]]);
      this.vx[i] = a[0];
      this.vy[i] = a[1];
      this.fx[i] = Math.fround(a[0]);
      this.fy[i] = Math.fround(a[1]);
    }
    this.turn = new Int8Array(n);
    for (i = 0; i < n; i++) {
      var e = this.edges[imod(i - 1, n)], v = this.edges[i];
      this.turn[i] = relativeCCW(e[0], e[1], e[2], e[3], v[2], v[3]);
    }
    this._cache = new Map();
    if (!(opts && opts.validate === false)) validatePolygon(this);
  }

  /* geometry.Polygon._validate: every algorithm assumes a simple polygon listed counter-clockwise in the
   * y-up frame (Java checks neither and fails late); throws Error otherwise. */
  function validatePolygon(poly) {
    if (!poly.vertices.every(function (v) { return isFinite(v[0]) && isFinite(v[1]); })) throw new Error("polygon vertices must be finite");
    var area = poly.signedArea();
    if (!(area > 0)) {
      throw new Error("vertices must be listed counter-clockwise in the y-up frame (got signed area " + area + "); reverse the vertex list");
    }
    var bad = selfIntersection(poly);
    if (bad !== null) throw new Error("polygon is not simple: edges " + bad[0] + " and " + bad[1] + " intersect");
  }

  /* geometry._self_intersection: the first pair of edges [i, j] (i < j) that intersect although they should
   * not, or null.  Non-adjacent edges must be disjoint (touching counts); adjacent edges may only share
   * their common vertex (no fold-back); no edge may have zero length. */
  function selfIntersection(poly) {
    var n = poly.n, E = poly.edges, i, j;
    for (i = 0; i < n; i++) if (E[i][0] === E[i][2] && E[i][1] === E[i][3]) return [i, (i + 1) % n];
    var sgn = function (x) { return x > 0 ? 1 : x < 0 ? -1 : 0; };
    var s1 = new Int8Array(n * n), s2 = new Int8Array(n * n);
    for (i = 0; i < n; i++) {
      var a = E[i];
      for (j = 0; j < n; j++) {
        var b = E[j];
        s1[i * n + j] = sgn((a[2] - a[0]) * (b[1] - a[1]) - (a[3] - a[1]) * (b[0] - a[0]));
        s2[i * n + j] = sgn((a[2] - a[0]) * (b[3] - a[1]) - (a[3] - a[1]) * (b[2] - a[0]));
      }
    }
    var straddle = function (p, q) { return s1[p * n + q] * s2[p * n + q] <= 0; };
    var col1 = function (p, q) { return s1[p * n + q] === 0 && s2[p * n + q] === 0; };
    var col = function (p, q) { return col1(p, q) && col1(q, p); };
    var fold = function (p) {
      var q = (p + 1) % n, a = E[p], b = E[q];
      return col(p, q) && (a[2] - a[0]) * (b[2] - b[0]) + (a[3] - a[1]) * (b[3] - b[1]) < 0;
    };
    for (i = 0; i < n; i++) {
      for (j = i + 1; j < n; j++) {
        var hit;
        if (j === (i + 1) % n || i === (j + 1) % n) {
          hit = (j === (i + 1) % n && fold(i)) || (i === (j + 1) % n && fold(j));
        } else if (col(i, j)) {
          var a2 = E[i], b2 = E[j];
          hit = Math.max(a2[0], a2[2]) >= Math.min(b2[0], b2[2]) && Math.max(b2[0], b2[2]) >= Math.min(a2[0], a2[2]) &&
            Math.max(a2[1], a2[3]) >= Math.min(b2[1], b2[3]) && Math.max(b2[1], b2[3]) >= Math.min(a2[1], a2[3]);
        } else {
          hit = straddle(i, j) && straddle(j, i);
        }
        if (hit) return [i, j];
      }
    }
    return null;
  }

  /* Shoelace area (positive for the counter-clockwise, y-up .dat files). */
  Polygon.prototype.signedArea = function () {
    var s = 0;
    for (var i = 0; i < this.n; i++) {
      var e = this.edges[i];
      s += e[0] * e[3] - e[2] * e[1];
    }
    return 0.5 * s;
  };

  Polygon.prototype.edge = function (i) {
    return this.edges[imod(i, this.n)];
  };

  Polygon.prototype.isReflex = function (i) {
    return this.turn[imod(i, this.n)] === 1;
  };

  Polygon.prototype.reflex = function () {
    var out = [];
    for (var i = 0; i < this.n; i++) if (this.turn[i] === 1) out.push(i);
    return out;
  };

  /* Polygon.pointInPolygon = GeneralPath.contains (Curve.pointCrossingsForLine per edge). */
  Polygon.prototype.contains = function (px, py) {
    if (!(px * 0 + py * 0 === 0)) return false;
    var fx = this.fx, fy = this.fy, n = this.n, crossings = 0;
    for (var k = 0; k < n; k++) {
      var x0 = fx[k], y0 = fy[k], k1 = k + 1 === n ? 0 : k + 1, x1 = fx[k1], y1 = fy[k1];
      if (py < y0 && py < y1) continue;
      if (py >= y0 && py >= y1) continue;
      if (px >= x0 && px >= x1) continue;
      if (!(px < x0 && px < x1)) {
        var xint = x0 + ((py - y0) * (x1 - x0)) / (y1 - y0);
        if (px >= xint) continue;
      }
      crossings += y0 < y1 ? 1 : -1;
    }
    return crossings !== 0;
  };

  /* Some edge within eps (ptSegDist) of the point. */
  Polygon.prototype.nearBoundary = function (px, py, eps) {
    var E = this.edges;
    for (var k = 0; k < this.n; k++) {
      var e = E[k];
      if (ptSegDist(e[0], e[1], e[2], e[3], px, py) < eps) return true;
    }
    return false;
  };

  /* Per-polygon tables for segments ending at a vertex: on[j] (v_j within eps of an edge),
   * inside[j] (contains(v_j)), rccw[k*n + j] = relativeCCW(e_k, v_j). */
  Polygon.prototype._vertexTables = function (eps) {
    var key = "vt" + eps, t = this._cache.get(key);
    if (t) return t;
    var n = this.n, on = new Uint8Array(n), inside = new Uint8Array(n), rccw = new Int8Array(n * n);
    for (var j = 0; j < n; j++) {
      on[j] = this.nearBoundary(this.vx[j], this.vy[j], eps) ? 1 : 0;
      inside[j] = this.contains(this.vx[j], this.vy[j]) ? 1 : 0;
      for (var k = 0; k < n; k++) {
        var e = this.edges[k];
        rccw[k * n + j] = relativeCCW(e[0], e[1], e[2], e[3], this.vx[j], this.vy[j]);
      }
    }
    t = { on: on, inside: inside, rccw: rccw };
    this._cache.set(key, t);
    return t;
  };

  /* Step (3) of segInPolygon: no edge crosses the segment at a point farther than eps from both ends. */
  function unblocked(poly, x1, y1, x2, y2, eps, rq1, rq2) {
    var n = poly.n, E = poly.edges, vx = poly.vx, vy = poly.vy;
    var r0 = relativeCCW(x1, y1, x2, y2, vx[0], vy[0]), rk = r0;
    for (var k = 0; k < n; k++) {
      var rn = k + 1 === n ? r0 : relativeCCW(x1, y1, x2, y2, vx[k + 1], vy[k + 1]);
      if (rk * rn <= 0) {
        var e = E[k];
        var re = rq1 ? rq1[k] * rq2(k) : relativeCCW(e[0], e[1], e[2], e[3], x1, y1) * relativeCCW(e[0], e[1], e[2], e[3], x2, y2);
        if (re <= 0) {
          var p = getIntersect(x1, y1, x2, y2, e[0], e[1], e[2], e[3], eps);
          if (p !== null && distance(p[0], p[1], x1, y1) > eps && distance(p[0], p[1], x2, y2) > eps) return false;
        }
      }
      rk = rn;
    }
    return true;
  }

  /* Algorithm.segInPolygon: midpoint contained, each end within eps of an edge or contained, and no
   * edge intersecting the segment (closed test) at a point farther than eps from both ends. */
  function segInPolygon(poly, x1, y1, x2, y2, eps) {
    eps = eps === undefined ? EPSILON : eps;
    if (!poly.contains((x1 + x2) / 2, (y1 + y2) / 2)) return false;
    if (!(poly.nearBoundary(x1, y1, eps) || poly.contains(x1, y1))) return false;
    if (!(poly.nearBoundary(x2, y2, eps) || poly.contains(x2, y2))) return false;
    return unblocked(poly, x1, y1, x2, y2, eps);
  }

  /* vis[j] = segInPolygon(q -> v_j) for every vertex (Uint8Array). */
  function vertexVisibility(poly, q, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var n = poly.n, qx = +q[0], qy = +q[1], vx = poly.vx, vy = poly.vy, out = new Uint8Array(n);
    if (!(poly.nearBoundary(qx, qy, eps) || poly.contains(qx, qy))) return out;
    var t = poly._vertexTables(eps), E = poly.edges, rq = new Int8Array(n), k, j;
    for (k = 0; k < n; k++) rq[k] = relativeCCW(E[k][0], E[k][1], E[k][2], E[k][3], qx, qy);
    var rccw = t.rccw;
    for (j = 0; j < n; j++) {
      if (!poly.contains((qx + vx[j]) / 2, (qy + vy[j]) / 2)) continue;
      if (!(t.on[j] || t.inside[j])) continue;
      var jj = j;
      out[j] = unblocked(poly, qx, qy, vx[j], vy[j], eps, rq, function (kk) { return rccw[kk * n + jj]; }) ? 1 : 0;
    }
    return out;
  }

  /* ---- visibility (visibility.py) ------------------------------------------ */

  /* The edge scan of getPhysicalGaps / getVisibilityPolygon: one record per visited edge,
   * [i, v_i visible, v_{i+1} visible, ip1, ip1n, ip2, ip2n]. */
  function scan(poly, q, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var n = poly.n, qx = +q[0], qy = +q[1], vx = poly.vx, vy = poly.vy, V = poly.vertices;
    var vis = vertexVisibility(poly, q, eps);
    var records = [];
    var i = 0;
    while (i < n) {
      var i1 = i + 1 === n ? 0 : i + 1;
      var in1 = !!vis[i], in2 = !!vis[i1];
      if (in1 && in2) {
        records.push([i, true, true, null, -1, null, -1]);
        i++;
        continue;
      }
      var p1 = V[i], p2 = V[i1], e = poly.edges[i];
      var ip1 = null, ip2 = null, ip1n = -1, ip2n = -1, dist1 = Number.MAX_VALUE, dist2 = Number.MAX_VALUE;
      var ni = i, lastI = i;
      for (var j = 0; j < n; j++) {
        if (!vis[j] || j === i || j === i1) continue;
        var ip = getIntersect(qx, qy, vx[j], vy[j], e[0], e[1], e[2], e[3], eps);
        if (ip === null || !onExtension(qx, qy, vx[j], vy[j], ip[0], ip[1])) continue;
        if (onExtension(e[0], e[1], e[2], e[3], ip[0], ip[1]) || onReverseExtension(e[0], e[1], e[2], e[3], ip[0], ip[1])) continue;
        if (!segInPolygon(poly, vx[j], vy[j], ip[0], ip[1], eps)) continue;
        var d = distance(ip[0], ip[1], p1[0], p1[1]);
        if (d < dist1 && !in1) {
          dist1 = d;
          ip1 = ip;
          ip1n = j;
        }
        d = distance(ip[0], ip[1], p2[0], p2[1]);
        if (d < dist2 && !in2) {
          dist2 = d;
          ip2 = ip;
          ip2n = j;
          ni = j;
        }
        /* Java's ip1 != ip2 is an identity test: the same object iff set by the same j */
        if ((in1 && ip2 !== null) || (in2 && ip1 !== null) || (ip1 !== null && ip2 !== null && ip1n !== ip2n)) {
          if (ni > i) i = ni - 1;
          else if (ip2 !== null) i = n - 1;
          break;
        }
      }
      records.push([lastI, in1, in2, ip1, ip1n, ip2, ip2n]);
      i++;
    }
    return records;
  }

  function PhysicalGap(startEdge, endEdge, startPoint, endPoint) {
    this.startEdge = startEdge;
    this.endEdge = endEdge;
    this.startPoint = startPoint || null;
    this.endPoint = endPoint || null;
  }

  /* getFullStartEdge: first edge that is entirely hidden. */
  PhysicalGap.prototype.fullStartEdge = function (n) {
    if (this.startPoint === null) return this.startEdge;
    return this.startEdge + 1 === n ? 0 : this.startEdge + 1;
  };

  /* getFullEndEdge: last edge that is entirely hidden. */
  PhysicalGap.prototype.fullEndEdge = function (n) {
    if (this.endPoint === null) return this.endEdge;
    return this.endEdge === 0 ? n - 1 : this.endEdge - 1;
  };

  /* [startEdge, endEdge, startPoint|null, endPoint|null] (golden-fixture format). */
  PhysicalGap.prototype.toList = function () {
    return [this.startEdge, this.endEdge, this.startPoint && this.startPoint.slice(), this.endPoint && this.endPoint.slice()];
  };

  /* Algorithm.getPhysicalGaps(poly, q).  {java: true}: Java's scan as is (Java output order; [] for a
   * point outside).  Default (compat="fixed"): throws GeometryError for a point strictly outside, and
   * repairs the scan where it is degenerate (on a line through a reflex vertex and another vertex):
   * repairGaps.  Elsewhere the result is Java's, bit for bit (visibility.physical_gaps). */
  function physicalGaps(poly, q, opts) {
    if (typeof opts === "number") opts = { eps: opts };
    var eps = (opts && opts.eps) || EPSILON, raw = scanGaps(poly, q, eps);
    if (opts && opts.java) return raw;
    var onBoundary = checkInside(poly, q, eps);
    return repairGaps(poly, q, raw, onBoundary);
  }

  function scanGaps(poly, q, eps) {
    var out = [], n = poly.n;
    scan(poly, q, eps).forEach(function (r) {
      if (r[3] !== null) out.push(new PhysicalGap(r[4], r[0], null, r[3]));
      if (r[5] !== null) out.push(new PhysicalGap(r[0], r[6] !== 0 ? r[6] - 1 : n - 1, r[5], null));
    });
    return out;
  }

  /* Algorithm.getVisibilityPolygon.  {java: true}: Java's vertex list (consecutive duplicates, B15).
   * Default (compat="fixed"): throws GeometryError for a point strictly outside; the boundary with every
   * hidden chain of the repaired gaps replaced by its window, consecutive repeats dropped. */
  function visibilityPolygon(poly, q, opts) {
    var eps = (opts && opts.eps) || EPSILON, n = poly.n, pts = [];
    if (!(opts && opts.java)) {
      checkInside(poly, q, eps);
      return visibilityFromGaps(poly, physicalGaps(poly, q, { eps: eps }));
    }
    scan(poly, q, eps).forEach(function (r) {
      if (r[1]) pts.push(poly.vertices[r[0]]);
      if (r[3] !== null) pts.push(r[3]);
      if (r[5] !== null) pts.push(r[5]);
      if (r[2]) pts.push(poly.vertices[(r[0] + 1) % n]);
    });
    return pts.map(function (p) { return p.slice(); });
  }

  /* Vertices not strictly inside a hidden chain plus both ends of every chain, in boundary order from v_0
   * (edge, then distance from its start), consecutive repeats dropped (visibility._visibility_from_gaps). */
  function visibilityFromGaps(poly, gaps) {
    var n = poly.n, V = poly.vertices, cum = perimeter(poly), per = cum[n], keyed = [], k;
    var chains = gaps.map(function (g) { return hiddenChain(poly, g); });
    for (k = 0; k < n; k++) {
      var ck = cum[k];
      if (!chains.some(function (c) { var d = pyMod(ck - c[0], per); return 0 < d && d < c[1]; })) keyed.push([k, 0, V[k]]);
    }
    gaps.forEach(function (g) {
      if (g.endPoint !== null) {
        keyed.push([g.startEdge, 0, V[g.startEdge]]);
        keyed.push([g.endEdge, distance(V[g.endEdge][0], V[g.endEdge][1], g.endPoint[0], g.endPoint[1]), g.endPoint]);
      } else {
        keyed.push([g.startEdge, distance(V[g.startEdge][0], V[g.startEdge][1], g.startPoint[0], g.startPoint[1]), g.startPoint]);
        var e = (g.endEdge + 1) % n;
        keyed.push([e, 0, V[e]]);
      }
    });
    keyed.sort(function (a, b) { return a[0] - b[0] || a[1] - b[1]; });
    var out = [];
    keyed.forEach(function (t) {
      var p = t[2], l = out[out.length - 1];
      if (!l || l[0] !== p[0] || l[1] !== p[1]) out.push(p);
    });
    while (out.length > 1 && out[out.length - 1][0] === out[0][0] && out[out.length - 1][1] === out[0][1]) out.pop();
    return out.map(function (p) { return p.slice(); });
  }

  /* ---- compat="fixed": point status and the repair of degenerate scans (visibility.py) ---- */

  /* geometry.GeometryError: a query the fixed mode refuses (a point outside the polygon). */
  class GeometryError extends Error {
    constructor(message) {
      super(message);
      this.name = "GeometryError";
    }
  }

  /* Throws GeometryError for a point strictly outside; returns whether it is within eps of an edge. */
  function checkInside(poly, q, eps) {
    var qx = +q[0], qy = +q[1];
    if (!(isFinite(qx) && isFinite(qy))) throw new GeometryError("point (" + q[0] + ", " + q[1] + ") is not finite");
    var on = poly.nearBoundary(qx, qy, eps);
    if (on || poly.contains(qx, qy)) return on;
    throw new GeometryError("point (" + qx + ", " + qy + ") is outside polygon " + poly.name);
  }

  /* "inside", "boundary" (within eps of an edge) or "outside" (where Java's scan returns []). */
  function pointStatus(poly, q, eps) {
    try {
      return checkInside(poly, q, eps === undefined ? EPSILON : eps) ? "boundary" : "inside";
    } catch (ex) {
      if (ex instanceof GeometryError) return "outside";
      throw ex;
    }
  }

  /* The fixed-mode repair of Java's scan raw at q (visibility.repair_gaps): cleanGaps, then, at a point
   * off the boundary where the scan changed or that lies near a vertex line, the exact gap structure
   * wins where it differs. */
  function repairGaps(poly, q, raw, onBoundary) {
    var gaps = cleanGaps(poly, q, raw);
    if (!onBoundary && (!sameGaps(gaps, raw) || nearVertexLine(poly, q))) {
      var exact = exactGaps(poly, q);
      if (exact !== null && structureKey(exact) !== structureKey(gaps)) return exact;
    }
    return gaps;
  }

  function samePoint(a, b) {
    return a === b || (a !== null && b !== null && a[0] === b[0] && a[1] === b[1]);
  }

  function sameGaps(a, b) {
    if (a.length !== b.length) return false;
    for (var k = 0; k < a.length; k++) {
      var g = a[k], h = b[k];
      if (g.startEdge !== h.startEdge || g.endEdge !== h.endEdge || !samePoint(g.startPoint, h.startPoint) ||
        !samePoint(g.endPoint, h.endPoint)) return false;
    }
    return true;
  }

  /* visibility._structure: sorted (startEdge, endEdge, endPoint set), as a string. */
  function structureKey(gaps) {
    return JSON.stringify(gaps.map(function (g) { return [g.startEdge, g.endEdge, g.endPoint !== null ? 1 : 0]; })
      .sort(function (a, b) { return a[0] - b[0] || a[1] - b[1] || a[2] - b[2]; }));
  }

  /* Closed visibility (simulate.segment_in_polygon), the predicate of the simulator's ground truth. */
  function orient(ax, ay, bx, by, cx, cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
  }

  /* orient(v_{k-1}, v_k, v_{k+1}) sign per vertex: +1 convex, -1 reflex, 0 straight. */
  function simTurns(poly) {
    var t = poly._cache.get("sim_turns");
    if (!t) {
      var n = poly.n, V = poly.vertices;
      t = new Int8Array(n);
      for (var k = 0; k < n; k++) {
        var a = V[imod(k - 1, n)], b = V[k], c = V[(k + 1) % n];
        var o = orient(a[0], a[1], b[0], b[1], c[0], c[1]);
        t[k] = o > 0 ? 1 : o < 0 ? -1 : 0;
      }
      poly._cache.set("sim_turns", t);
    }
    return t;
  }

  /* Direction (dx, dy) points into the closed interior angle of the (CCW) polygon at v_k. */
  function inCone(poly, k, dx, dy) {
    var n = poly.n, u = poly.vertices[imod(k - 1, n)], v = poly.vertices[k], w = poly.vertices[(k + 1) % n];
    var c1 = (v[0] - u[0]) * dy - (v[1] - u[1]) * dx;
    var c2 = (w[0] - v[0]) * dy - (w[1] - v[1]) * dx;
    var t = simTurns(poly)[k];
    if (t > 0) return c1 >= 0 && c2 >= 0;
    if (t < 0) return c1 >= 0 || c2 >= 0;
    return c2 >= 0;
  }

  function throughVertexOk(poly, ax, ay, bx, by, k) {
    var v = poly.vertices[k], dx = bx - ax, dy = by - ay;
    if (!((v[0] - ax) * dx + (v[1] - ay) * dy > 0 && (v[0] - bx) * dx + (v[1] - by) * dy < 0)) return true;
    return inCone(poly, k, dx, dy) && inCone(poly, k, -dx, -dy);
  }

  /* The segment a b lies in the closed polygon: no proper edge crossing, and at every vertex strictly
   * inside the segment both directions point into the interior angle. */
  function closedSegmentInPolygon(poly, a, b) {
    var ax = a[0], ay = a[1], bx = b[0], by = b[1], n = poly.n, vx = poly.vx, vy = poly.vy, k;
    var o = new Float64Array(n);
    for (k = 0; k < n; k++) o[k] = orient(ax, ay, bx, by, vx[k], vy[k]);
    for (k = 0; k < n; k++) {
      var k1 = k + 1 === n ? 0 : k + 1, o1 = o[k], o2 = o[k1];
      if ((o1 < 0 && 0 < o2) || (o2 < 0 && 0 < o1)) {
        var o3 = orient(vx[k], vy[k], vx[k1], vy[k1], ax, ay);
        var o4 = orient(vx[k], vy[k], vx[k1], vy[k1], bx, by);
        if ((o3 < 0 && 0 < o4) || (o4 < 0 && 0 < o3)) return false;
      }
    }
    for (k = 0; k < n; k++) if (o[k] === 0 && !throughVertexOk(poly, ax, ay, bx, by, k)) return false;
    return true;
  }

  /* A point just inside the polygon at the middle of the edge at the gap's anchor vertex. */
  function corePoint(poly, g) {
    var e = poly.edges[imod(g.endPoint !== null ? g.startEdge : g.endEdge, poly.n)];
    var x1 = e[0], y1 = e[1], x2 = e[2], y2 = e[3], d = 1e-6;
    return [0.5 * (x1 + x2) - d * (y2 - y1), 0.5 * (y1 + y2) + d * (x2 - x1)];
  }

  function anchorHidden(poly, q, g) {
    return !closedSegmentInPolygon(poly, [+q[0], +q[1]], corePoint(poly, g));
  }

  function fullEdges(g, n) {
    return imod(g.fullEndEdge(n) - g.fullStartEdge(n), n) + 1;
  }

  /* The gap with the same anchor and free end that hides the rest of the boundary, or null. */
  function complementGap(g, n) {
    var i, j, c;
    if (g.endPoint !== null) {
      i = g.endEdge;
      j = g.startEdge;
      c = new PhysicalGap(i, imod(j - 1, n), g.endPoint, null);
    } else {
      i = g.startEdge;
      j = (g.endEdge + 1) % n;
      c = new PhysicalGap(j, i, null, g.startPoint);
    }
    if (j === i || j === (i + 1) % n) return null;
    return c;
  }

  /* visibility._clean_gaps: a gap whose anchor edge is visible is a phantom, replaced by its complement
   * if that hides its anchor edge (else dropped); two gaps with the same hidden chain are one, the
   * description with more fully hidden edges is kept. */
  function cleanGaps(poly, q, gaps) {
    if (!gaps.length) return gaps;
    var n = poly.n, ok = gaps.map(function (g) { return anchorHidden(poly, q, g); });
    if (ok.every(Boolean) && gaps.length < 2) return gaps;
    var keep = [], k, m;
    gaps.forEach(function (g, i) {
      if (ok[i]) {
        keep.push(g);
      } else {
        var c = complementGap(g, n);
        if (c !== null && anchorHidden(poly, q, c)) keep.push(c);
      }
    });
    if (keep.length < 2) return keep;
    var per = 0;
    poly.edges.forEach(function (e) { per += distance(e[0], e[1], e[2], e[3]); });
    var tol = 1e-9 * Math.max(per, 1.0);
    var chains = keep.map(function (g) { return hiddenChain(poly, g); }), out = [], used = keep.map(function () { return false; });
    for (k = 0; k < keep.length; k++) {
      if (used[k]) continue;
      var best = k;
      for (m = k + 1; m < keep.length; m++) {
        if (!used[m] && sameChain(chains[k], chains[m], per, tol)) {
          used[m] = true;
          if (fullEdges(keep[m], n) > fullEdges(keep[best], n)) best = m;
        }
      }
      out.push(keep[best]);
    }
    return out;
  }

  function sameChain(a, b, per, tol) {
    var d = Math.abs(a[0] - b[0]) % per;
    return Math.min(d, per - d) <= tol && Math.abs(a[1] - b[1]) <= tol;
  }

  var NEAR = 1e-3; /* distance from a line through a reflex vertex and another vertex that triggers the exact check */

  /* q within NEAR of a line through a reflex vertex and another vertex (where Java's scan can go wrong). */
  function nearVertexLine(poly, q) {
    var t = poly._cache.get("vertex_lines"), k;
    if (!t) {
      var turns = exactTurns(poly), n = poly.n, rx = [], ry = [], ux = [], uy = [];
      for (var r = 0; r < n; r++) {
        if (turns[r] >= 0) continue;
        for (var v = 0; v < n; v++) {
          var dx = poly.vx[v] - poly.vx[r], dy = poly.vy[v] - poly.vy[r], ln = Math.sqrt(dx * dx + dy * dy);
          if (!(ln > 0)) continue;
          rx.push(poly.vx[r]);
          ry.push(poly.vy[r]);
          ux.push(dx / ln);
          uy.push(dy / ln);
        }
      }
      t = { rx: rx, ry: ry, ux: ux, uy: uy };
      poly._cache.set("vertex_lines", t);
    }
    var qx = +q[0], qy = +q[1];
    for (k = 0; k < t.rx.length; k++) {
      if (Math.abs(t.ux[k] * (qy - t.ry[k]) - t.uy[k] * (qx - t.rx[k])) < NEAR) return true;
    }
    return false;
  }

  /* Exact arithmetic on doubles: x = m * 2^e with an integer (BigInt) m. */
  var F64VIEW = new DataView(new ArrayBuffer(8));
  function decompose(x) {
    if (x === 0) return [BigInt(0), 0];
    F64VIEW.setFloat64(0, x);
    var hi = F64VIEW.getUint32(0), lo = F64VIEW.getUint32(4), e = (hi >>> 20) & 0x7ff;
    var m = BigInt(hi & 0xfffff) * BigInt(4294967296) + BigInt(lo);
    if (e === 0) e = 1;
    else m += BigInt(1) << BigInt(52);
    return [hi >>> 31 ? -m : m, e - 1075];
  }

  /* sign((b - a) x (d - c)), exact on the double inputs (float filter, exact fallback). */
  function crossSign(ax, ay, bx, by, cx, cy, dx, dy) {
    var v = (bx - ax) * (dy - cy) - (by - ay) * (dx - cx);
    if (Math.abs(v) > 1e-9 * (Math.abs(bx - ax) + Math.abs(by - ay)) * (Math.abs(dx - cx) + Math.abs(dy - cy))) return v > 0 ? 1 : -1;
    var ds = [ax, ay, bx, by, cx, cy, dx, dy].map(decompose), emin = Infinity;
    ds.forEach(function (d) { if (d[0] !== BigInt(0) && d[1] < emin) emin = d[1]; });
    if (emin === Infinity) return 0;
    var b = ds.map(function (d) { return d[0] << BigInt(d[1] - emin); });
    var w = (b[2] - b[0]) * (b[7] - b[5]) - (b[3] - b[1]) * (b[6] - b[4]);
    return w > BigInt(0) ? 1 : w < BigInt(0) ? -1 : 0;
  }

  /* sign((b - a) x (v_k - a)) for every vertex, exact; v_skip is b itself (sign 0). */
  function orientSigns(poly, ax, ay, bx, by, skip) {
    var n = poly.n, xs = poly.vx, ys = poly.vy, sg = new Int8Array(n), A = Math.abs(bx - ax) + Math.abs(by - ay);
    for (var k = 0; k < n; k++) {
      var o = (bx - ax) * (ys[k] - ay) - (by - ay) * (xs[k] - ax);
      var bound = 1e-9 * A * (Math.abs(xs[k] - ax) + Math.abs(ys[k] - ay));
      if (Math.abs(o) <= bound) sg[k] = k === skip ? 0 : crossSign(ax, ay, bx, by, ax, ay, xs[k], ys[k]);
      else sg[k] = o > 0 ? 1 : -1;
    }
    return sg;
  }

  /* sign((v_k - v_{k-1}) x (v_{k+1} - v_k)): +1 convex, -1 reflex, 0 straight (exact). */
  function exactTurns(poly) {
    var t = poly._cache.get("exact_turns");
    if (!t) {
      var V = poly.vertices, n = poly.n;
      t = new Int8Array(n);
      for (var k = 0; k < n; k++) {
        var u = V[imod(k - 1, n)], v = V[k], w = V[(k + 1) % n];
        t[k] = crossSign(u[0], u[1], v[0], v[1], v[0], v[1], w[0], w[1]);
      }
      poly._cache.set("exact_turns", t);
    }
    return t;
  }

  /* The direction f -> t points into the closed interior angle at v_k (exact). */
  function cone(poly, k, fx, fy, tx, ty) {
    var n = poly.n, u = poly.vertices[imod(k - 1, n)], v = poly.vertices[k], w = poly.vertices[(k + 1) % n];
    var c1 = crossSign(u[0], u[1], v[0], v[1], fx, fy, tx, ty);
    var c2 = crossSign(v[0], v[1], w[0], w[1], fx, fy, tx, ty);
    var t = exactTurns(poly)[k];
    if (t > 0) return c1 >= 0 && c2 >= 0;
    if (t < 0) return c1 >= 0 || c2 >= 0;
    return c2 >= 0;
  }

  /* sign((v_{k+1} - v_k) x (p - v_k)), exact. */
  function edgeSign(poly, k, px, py) {
    var e = poly.edges[k], x1 = e[0], y1 = e[1], x2 = e[2], y2 = e[3];
    var o = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1);
    var bound = 1e-9 * (Math.abs(x2 - x1) + Math.abs(y2 - y1)) * (Math.abs(px - x1) + Math.abs(py - y1));
    if (Math.abs(o) <= bound) return crossSign(x1, y1, x2, y2, x1, y1, px, py);
    return o > 0 ? 1 : -1;
  }

  /* Closed visibility of v_r from q (exact); sg: orientation signs of the vertices w.r.t. q -> v_r. */
  function sees(poly, qx, qy, r, sg) {
    var n = poly.n, rx = poly.vx[r], ry = poly.vy[r], dx = rx - qx, dy = ry - qy, rm = imod(r - 1, n), k;
    for (k = 0; k < n; k++) {
      if (k === r || k === rm || !(sg[k] * sg[(k + 1) % n] < 0)) continue;
      if (edgeSign(poly, k, qx, qy) * edgeSign(poly, k, rx, ry) < 0) return false;
    }
    var dd = dx * dx + dy * dy;
    for (k = 0; k < n; k++) {
      var along = (poly.vx[k] - qx) * dx + (poly.vy[k] - qy) * dy;
      if (sg[k] === 0 && along > 0 && along < dd && k !== r &&
        !(cone(poly, k, qx, qy, rx, ry) && cone(poly, k, rx, ry, qx, qy))) return false;
    }
    return true;
  }

  /* visibility._exact_gaps: the gaps seen from a point strictly inside, from exact orientation signs
   * (closed visibility): a visible reflex (or straight) vertex r anchors a gap iff the line of sight
   * continues past r into the polygon with the boundary at r on one side; the hidden chain leaves r along
   * the neighbour closer in angle to the line of sight and ends at the first proper crossing of an edge
   * (ip = getIntersect) or vertex touched from that side past r.  null if inconsistent. */
  function exactGaps(poly, q) {
    var n = poly.n, V = poly.vertices, xs = poly.vx, ys = poly.vy, qx = +q[0], qy = +q[1];
    var turns = exactTurns(poly), out = [], k, r, ahead = new Float64Array(n);
    for (r = 0; r < n; r++) {
      if (turns[r] > 0) continue;
      var a = imod(r - 1, n), b = (r + 1) % n;
      var dxr = xs[r] - qx, dyr = ys[r] - qy, sz = Math.abs(dxr) + Math.abs(dyr);
      var oa = dxr * (ys[a] - qy) - dyr * (xs[a] - qx), ob = dxr * (ys[b] - qy) - dyr * (xs[b] - qx);
      var clear = Math.abs(oa) > 1e-9 * sz * (Math.abs(xs[a] - qx) + Math.abs(ys[a] - qy)) &&
        Math.abs(ob) > 1e-9 * sz * (Math.abs(xs[b] - qx) + Math.abs(ys[b] - qy));
      if (clear && oa * ob < 0) continue;
      var rx = V[r][0], ry = V[r][1], dx = rx - qx, dy = ry - qy;
      if (dx === 0 && dy === 0) continue;
      var sg = orientSigns(poly, qx, qy, rx, ry, r);
      for (k = 0; k < n; k++) ahead[k] = (xs[k] - rx) * dx + (ys[k] - ry) * dy;
      var sides = [];
      for (var u of [a, b]) {
        if (sg[u] !== 0) {
          sides.push(sg[u]);
        } else if (ahead[u] > 0) {
          sides = [];
          break;
        }
      }
      if (!sides.length || sides.some(function (x) { return x !== sides[0]; }) || !cone(poly, r, qx, qy, rx, ry)) continue;
      if (!sees(poly, qx, qy, r, sg)) continue;
      var side = sides[0], hid;
      var cand = [a, b].filter(function (x) { return sg[x] === side; });
      if (cand.length === 2) {
        var ua = cand[0], ub = cand[1];
        var da = (xs[ua] - rx) * dx + (ys[ua] - ry) * dy, db = (xs[ub] - rx) * dx + (ys[ub] - ry) * dy;
        var ca = Math.abs(dx * (ys[ua] - ry) - dy * (xs[ua] - rx)), cb = Math.abs(dx * (ys[ub] - ry) - dy * (xs[ub] - rx));
        hid = da * cb > db * ca ? ua : ub;
      } else {
        hid = cand[0];
      }
      var ccw = hid === b, bestT = Infinity, best = null, dd = dx * dx + dy * dy;
      /* first contact past r from side `side`: a vertex on the line touched from that side ... */
      for (k = 0; k < n; k++) {
        if (k === r || sg[k] !== 0 || !(ahead[k] > 0) || !(sg[imod(k - 1, n)] === side || sg[(k + 1) % n] === side)) continue;
        var tv = ((xs[k] - qx) * dx + (ys[k] - qy) * dy) / dd;
        if (best === null || tv < bestT) {
          bestT = tv;
          best = ["v", k];
        }
      }
      /* ... or a proper crossing of an edge */
      var crossings = [];
      for (k = 0; k < n; k++) {
        if (k === r || k === a || !(sg[k] * sg[(k + 1) % n] < 0)) continue;
        var ex = poly.edges[k][2] - xs[k], ey = poly.edges[k][3] - ys[k], den = dx * ey - dy * ex;
        crossings.push([((xs[k] - qx) * ey - (ys[k] - qy) * ex) / den, k]);
      }
      crossings.sort(function (p1, p2) {
        var n1 = p1[0] !== p1[0], n2 = p2[0] !== p2[0];
        if (n1 || n2) return n1 === n2 ? p1[1] - p2[1] : n1 ? 1 : -1;
        return p1[0] < p2[0] ? -1 : p1[0] > p2[0] ? 1 : p1[1] - p2[1];
      });
      for (var c = 0; c < crossings.length; c++) {
        var t = crossings[c][0];
        if (!(t > 0) || t >= bestT) continue;
        k = crossings[c][1];
        var k1 = (k + 1) % n;
        if (crossSign(V[k][0], V[k][1], V[k1][0], V[k1][1], V[k][0], V[k][1], qx, qy) *
          crossSign(V[k][0], V[k][1], V[k1][0], V[k1][1], V[k][0], V[k][1], rx, ry) < 0) continue;
        bestT = t;
        best = ["e", k];
        break;
      }
      if (best === null) return null;
      var ip, eEnd, eStart;
      k = best[1];
      if (best[0] === "e") {
        var e = poly.edges[k];
        ip = getIntersect(qx, qy, rx, ry, e[0], e[1], e[2], e[3], EPSILON);
        if (ip === null || !(Math.min(e[0], e[2]) <= ip[0] && ip[0] <= Math.max(e[0], e[2]) &&
          Math.min(e[1], e[3]) <= ip[1] && ip[1] <= Math.max(e[1], e[3]))) ip = [qx + bestT * dx, qy + bestT * dy];
        eEnd = eStart = k;
      } else {
        ip = V[k].slice();
        eEnd = k;
        eStart = imod(k - 1, n);
      }
      if (ccw) {
        if (eEnd === r || eEnd === a) return null;
        out.push([eEnd, distance(V[eEnd][0], V[eEnd][1], ip[0], ip[1]), new PhysicalGap(r, eEnd, null, ip)]);
      } else {
        if (a === eStart || r === eStart) return null;
        out.push([eStart, distance(V[eStart][0], V[eStart][1], ip[0], ip[1]), new PhysicalGap(eStart, a, ip, null)]);
      }
    }
    out.sort(function (x, y) { return x[0] - y[0] || x[1] - y[1]; });
    return out.map(function (x) { return x[2]; });
  }

  /* ---- cuts (cuts.py) ------------------------------------------------------ */

  /* Nearest boundary hits of a line as the Java loops pick them (first strict minimum in edge order).
   * mode "forward"/"backward": getIntersect(ray, e_j) past P2 (onExtension) / before P1; returns [hit].
   * mode "both": getIntersect(e_j, ray) split by onExtension(ray): returns [ccwPoint, cwPoint]. */
  function nearestHits(poly, ray, excl, mode, eps) {
    var n = poly.n, E = poly.edges, skip = new Uint8Array(n);
    excl.forEach(function (k) { skip[imod(k, n)] = 1; });
    var rx1 = ray[0], ry1 = ray[1], rx2 = ray[2], ry2 = ray[3];
    var best1 = Infinity, best2 = Infinity, h1 = null, h2 = null;
    for (var j = 0; j < n; j++) {
      if (skip[j]) continue;
      var e = E[j], p;
      if (mode === "both") p = getIntersect(e[0], e[1], e[2], e[3], rx1, ry1, rx2, ry2, eps);
      else p = getIntersect(rx1, ry1, rx2, ry2, e[0], e[1], e[2], e[3], eps);
      if (p === null || !(ptSegDist(e[0], e[1], e[2], e[3], p[0], p[1]) < eps)) continue;
      var d = ptSegDist(rx1, ry1, rx2, ry2, p[0], p[1]);
      if (mode === "forward" || mode === "backward") {
        var ok = mode === "forward" ? onExtension(rx1, ry1, rx2, ry2, p[0], p[1]) : onReverseExtension(rx1, ry1, rx2, ry2, p[0], p[1]);
        if (ok && d < best1) {
          best1 = d;
          h1 = p;
        }
      } else if (onExtension(rx1, ry1, rx2, ry2, p[0], p[1])) {
        if (d < best1) {
          best1 = d;
          h1 = p;
        }
      } else if (d < best2) {
        best2 = d;
        h2 = p;
      }
    }
    return [h1, h2];
  }

  /* Every inflection ray, in Java output order: {i, fromPrev, kind (-1 general, 1 non-general,
   * 0 collinear), hit}. */
  function inflectionRays(poly, eps) {
    var key = "infl" + eps, out = poly._cache.get(key);
    if (out) return out;
    var n = poly.n, E = poly.edges, V = poly.vertices;
    out = [];
    for (var i = 0; i < n; i++) {
      if (!poly.isReflex(i)) continue;
      var poi = V[i];
      var pp = E[imod(i - 2, n)], prev = E[imod(i - 1, n)], next = E[i], nn = E[imod(i + 1, n)];
      out.push({ i: i, fromPrev: true, kind: rccwLine(pp, poi),
        hit: nearestHits(poly, prev, [i - 2, i - 1, i], "forward", eps)[0] });
      out.push({ i: i, fromPrev: false, kind: rccwLine(nn, poi),
        hit: nearestHits(poly, next, [i + 1, i - 1, i], "backward", eps)[0] });
    }
    poly._cache.set(key, out);
    return out;
  }

  function rayLine(poly, r) {
    if (r.hit === null) throw new Error("inflection ray at vertex " + r.i + " hits no edge (Java NullPointerException)");
    var v = poly.vertices[r.i];
    return [v[0], v[1], r.hit[0], r.hit[1]];
  }

  /* Algorithm.getInflection(poly, kind): kind "GENERAL", "NONGENERAL" or "ALL". */
  function inflections(poly, kind, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var want = { GENERAL: [-1], NONGENERAL: [1], ALL: [-1, 0, 1] }[kind || "ALL"];
    return inflectionRays(poly, eps).filter(function (r) { return want.indexOf(r.kind) >= 0; })
      .map(function (r) { return rayLine(poly, r); });
  }

  /* Algorithm.getGeneralInflectionCut. */
  function generalInflectionCuts(poly, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var n = poly.n;
    return inflectionRays(poly, eps).filter(function (r) { return r.kind === -1; }).map(function (r) {
      return { line: rayLine(poly, r), type: CUT.GENERAL_INFLECTION, fromLine: r.fromPrev ? imod(r.i - 1, n) : r.i,
        counterClockwise: r.fromPrev };
    });
  }

  function sideProducts(poly, i, j) {
    var n = poly.n, vx = poly.vx, vy = poly.vy, x1 = vx[i], y1 = vy[i], x2 = vx[j], y2 = vy[j];
    var im = imod(i - 1, n), ip = (i + 1) % n, jm = imod(j - 1, n), jp = (j + 1) % n;
    return [relativeCCW(x1, y1, x2, y2, vx[im], vy[im]) * relativeCCW(x1, y1, x2, y2, vx[ip], vy[ip]),
      relativeCCW(x1, y1, x2, y2, vx[jm], vy[jm]) * relativeCCW(x1, y1, x2, y2, vx[jp], vy[jp])];
  }

  function tangentHits(poly, i, j, eps) {
    var V = poly.vertices;
    return nearestHits(poly, [V[i][0], V[i][1], V[j][0], V[j][1]], [i - 1, i, j - 1, j], "both", eps);
  }

  /* Algorithm.getBitangentCut: two rays per bitangent pair (this=i, opp=j), then (this=j, opp=i). */
  function bitangentCuts(poly, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var key = "bt" + eps, out = poly._cache.get(key);
    if (out) return out.slice();
    var n = poly.n, V = poly.vertices, vx = poly.vx, vy = poly.vy, E = poly.edges, turn = poly.turn;
    out = [];
    for (var i = 0; i < n; i++) {
      if (turn[i] !== 1) continue;
      for (var j = i + 2; j < n; j++) {
        if (turn[j] !== 1) continue;
        var sp = sideProducts(poly, i, j);
        if (sp[0] !== 1 || sp[1] !== 1) continue;
        var x1 = vx[i], y1 = vy[i], x2 = vx[j], y2 = vy[j];
        var skip = [imod(i - 1, n), i, imod(j - 1, n), j], hit = false;
        for (var k = 0; k < n && !hit; k++) {
          if (skip.indexOf(k) >= 0) continue;
          var e = E[k];
          hit = relativeCCW(x1, y1, x2, y2, e[0], e[1]) * relativeCCW(x1, y1, x2, y2, e[2], e[3]) <= 0 &&
            relativeCCW(e[0], e[1], e[2], e[3], x1, y1) * relativeCCW(e[0], e[1], e[2], e[3], x2, y2) <= 0;
        }
        if (hit) continue;
        var h = tangentHits(poly, i, j, eps), c = h[0], w = h[1];
        if (c === null || w === null) throw new Error("bitangent " + i + "-" + j + " hits no edge (Java NullPointerException)");
        var p1 = V[i], p2 = V[j], a, b;
        if (distanceSq(p1[0], p1[1], c[0], c[1]) < distanceSq(p2[0], p2[1], c[0], c[1])) {
          a = [p1[0], p1[1], c[0], c[1]];
          b = [p2[0], p2[1], w[0], w[1]];
        } else {
          a = [p1[0], p1[1], w[0], w[1]];
          b = [p2[0], p2[1], c[0], c[1]];
        }
        out.push({ line: a, type: CUT.BITANGENT, thisPoint: i, oppositePoint: j, oppositeSegment: b,
          curveToPoint: V[(i + 1) % n] });
        out.push({ line: b, type: CUT.BITANGENT, thisPoint: j, oppositePoint: i, oppositeSegment: a,
          curveToPoint: V[(j + 1) % n] });
      }
    }
    poly._cache.set(key, out);
    return out.slice();
  }

  /* Algorithm.getSingletangentCut: tangent at one end, crossing at the other, segInPolygon(v_i v_j). */
  function singleTangentCuts(poly, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var key = "st" + eps, out = poly._cache.get(key);
    if (out) return out.slice();
    var n = poly.n, V = poly.vertices, vx = poly.vx, vy = poly.vy, turn = poly.turn;
    out = [];
    for (var i = 0; i < n; i++) {
      for (var j = i + 2; j < n; j++) {
        if (turn[i] !== 1 && turn[j] !== 1) continue;
        var sp = sideProducts(poly, i, j);
        if (sp[0] * sp[1] !== -1) continue;
        if (!segInPolygon(poly, vx[i], vy[i], vx[j], vy[j], eps)) continue;
        var h = tangentHits(poly, i, j, eps), c = h[0], w = h[1];
        var a = sp[0] === 1 ? V[i] : V[j], b = sp[0] === 1 ? V[j] : V[i];
        var end = c !== null && distanceSq(a[0], a[1], c[0], c[1]) < distanceSq(b[0], b[1], c[0], c[1]) ? c : w;
        if (end === null) throw new Error("single tangent " + i + "-" + j + " hits no edge (Java NullPointerException)");
        out.push({ line: [a[0], a[1], end[0], end[1]], type: CUT.SINGLETANGENT });
      }
    }
    poly._cache.set(key, out);
    return out.slice();
  }

  /* Algorithm.getCuts: single tangents + non-general inflections + general inflections + bitangents. */
  function allCuts(poly, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var key = "cuts" + eps, out = poly._cache.get(key);
    if (!out) {
      out = singleTangentCuts(poly, eps);
      inflections(poly, "NONGENERAL", eps).forEach(function (l) { out.push({ line: l, type: CUT.NONGENERAL_INFLECTION }); });
      out = out.concat(generalInflectionCuts(poly, eps), bitangentCuts(poly, eps));
      poly._cache.set(key, out);
    }
    return out.slice();
  }

  /* ---- Path and critical points ------------------------------------------- */

  /* drawable.Path: a polyline; segments[i] = points[i] -> points[i+1]; prefix arc lengths. */
  function Path(points) {
    var pts = points.map(function (p) { return [+p[0], +p[1]]; });
    if (pts.length < 2) throw new Error("a path needs at least 2 points");
    this.points = pts;
    this.segments = [];
    this.prefix = [0];
    var d = 0;
    for (var i = 0; i + 1 < pts.length; i++) {
      var s = [pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1]];
      this.segments.push(s);
      d += distance(s[0], s[1], s[2], s[3]);
      this.prefix.push(d);
    }
  }

  Path.prototype.length = function () {
    return this.prefix[this.prefix.length - 1];
  };

  /* ptDistFromStart(p, lines, i): length of segments < i plus |p - P1_i|. */
  Path.prototype.ptDistOnSegment = function (p, i) {
    var s = this.segments[i];
    return this.prefix[i] + distance(s[0], s[1], p[0], p[1]);
  };

  /* Path.ptDistFromStart(p): the last segment within eps wins (quirk B17); NaN if on none. */
  Path.prototype.ptDistFromStart = function (p, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var d = NaN;
    for (var i = 0; i < this.segments.length; i++) {
      var s = this.segments[i];
      if (ptSegDist(s[0], s[1], s[2], s[3], p[0], p[1]) < eps) d = this.ptDistOnSegment(p, i);
    }
    return d;
  };

  /* Path.getAllCriticalPoints(cuts): crossings sorted by arc length (equal distances keep the last
   * one, B2), plus start and end (type NONE). */
  function allCriticalPoints(path, cuts, eps) {
    eps = eps === undefined ? EPSILON : eps;
    var m = new Map(), segs = path.segments;
    for (var i = 0; i < segs.length; i++) {
      var s = segs[i];
      for (var j = 0; j < cuts.length; j++) {
        var c = cuts[j].line;
        var p = getIntersect(s[0], s[1], s[2], s[3], c[0], c[1], c[2], c[3], eps);
        if (p === null || !(ptSegDist(s[0], s[1], s[2], s[3], p[0], p[1]) < eps) ||
            !(ptSegDist(c[0], c[1], c[2], c[3], p[0], p[1]) < eps)) continue;
        var d = path.ptDistOnSegment(p, i);
        m.set(d, { distance: d, point: p, cutType: cuts[j].type, cutIndex: j, segIndex: i, seg: s, cut: cuts[j] });
      }
    }
    m.set(0, { distance: 0, point: path.points[0], cutType: CUT.NONE, cutIndex: -1, segIndex: 0, seg: segs[0], cut: null });
    var end = path.points[path.points.length - 1], de = path.ptDistFromStart(end, eps);
    m.set(de, { distance: de, point: end, cutType: CUT.NONE, cutIndex: -1, segIndex: segs.length - 1,
      seg: segs[segs.length - 1], cut: null });
    var keys = Array.from(m.keys()).sort(function (a, b) { return a - b; });
    return keys.map(function (k) { return m.get(k); });
  }

  /* ---- java.util.HashSet<Integer> order, Gap, helpers (gaps.py) ---------- */

  function spread(h) {
    h = h >>> 0;
    return (h ^ (h >>> 16)) >>> 0;
  }

  /* HashSet<Integer> with Java's iteration order (bucket, then insertion order inside a bucket). */
  function JavaIntHashSet(items) {
    this.items = [];
    this.cap = 16;
    var self = this;
    (items || []).forEach(function (x) { self.add(x); });
  }

  JavaIntHashSet.prototype.add = function (x) {
    if (this.items.indexOf(x) >= 0) return false;
    var mask = this.cap - 1, b = spread(x) & mask, chain = 0;
    for (var k = 0; k < this.items.length; k++) if ((spread(this.items[k]) & mask) === b) chain++;
    this.items.push(x);
    if (chain >= 8) {
      if (this.cap < 64) this.cap *= 2;
      else throw new Error("HashSet tree bins are not modelled");
    }
    if (this.items.length > 0.75 * this.cap) this.cap *= 2;
    return true;
  };

  JavaIntHashSet.prototype.toList = function () {
    var mask = this.cap - 1;
    return this.items.slice().sort(function (a, b) { return (spread(a) & mask) - (spread(b) & mask); });
  };

  JavaIntHashSet.prototype.has = function (x) {
    return this.items.indexOf(x) >= 0;
  };

  Object.defineProperty(JavaIntHashSet.prototype, "size", { get: function () { return this.items.length; } });

  function javaHashsetOrder(ids) {
    return new JavaIntHashSet(ids).toList();
  }

  /* gap.Gap: a shadow component in a gap set (id, iEdge, relativeTime, links to the next set). */
  function Gap(physical, id, iedge) {
    this.id = id === undefined ? -1 : id;
    this.iedge = iedge === undefined ? -1 : iedge;
    this.relativeTime = 0;
    this.state = null;
    this.to = new JavaIntHashSet();
    this.physical = physical || null;
  }

  Object.defineProperty(Gap.prototype, "startEdge", { get: function () { return this.physical.startEdge; } });
  Object.defineProperty(Gap.prototype, "endEdge", { get: function () { return this.physical.endEdge; } });
  Gap.prototype.fullStartEdge = function (n) { return this.physical.fullStartEdge(n); };
  Gap.prototype.fullEndEdge = function (n) { return this.physical.fullEndEdge(n); };
  Gap.prototype.addGap = function (gid) { this.to.add(gid); };
  /* toGapSet in Java HashSet iteration order. */
  Gap.prototype.toGaps = function () { return this.to.toList(); };
  /* [id, iEdge, toGapSet] (golden-fixture format). */
  Gap.prototype.toList = function () { return [this.id, this.iedge, this.toGaps()]; };
  Gap.prototype.toString = function () { return formatGap(this); };

  /* Java Double.toString layout of the shortest round-trip digits (what JDK >= 19 prints). */
  function javaDoubleStr(x) {
    if (x !== x) return "NaN";
    if (x === Infinity) return "Infinity";
    if (x === -Infinity) return "-Infinity";
    var sign = x < 0 || Object.is(x, -0) ? "-" : "";
    var ax = Math.abs(x);
    if (ax === 0) return sign + "0.0";
    var s = String(ax), exp = 0, k = s.indexOf("e");
    if (k >= 0) {
      exp = parseInt(s.slice(k + 1), 10);
      s = s.slice(0, k);
    }
    var dot = s.indexOf("."), ip = dot >= 0 ? s.slice(0, dot) : s, fp = dot >= 0 ? s.slice(dot + 1) : "";
    var digits = ip + fp, e10 = ip.length - 1 + exp;
    var lead = digits.search(/[1-9]/);
    e10 -= lead;
    digits = digits.slice(lead).replace(/0+$/, "") || "0";
    if (ax >= 1e-3 && ax < 1e7) {
      if (e10 >= 0) {
        var a = digits.slice(0, e10 + 1);
        while (a.length < e10 + 1) a += "0";
        return sign + a + "." + (digits.slice(e10 + 1) || "0");
      }
      return sign + "0." + "0".repeat(-e10 - 1) + digits;
    }
    return sign + digits[0] + "." + (digits.slice(1) || "0") + "E" + e10;
  }

  /* Gap.toString: [id, iEdge], [id, iEdge | a b ]; a never-see-evader state is printed after iEdge. */
  function formatGap(g) {
    var s = "[" + g.id + ", " + g.iedge;
    if (g.state !== null && g.state !== undefined) s += g.state === "clear" ? "| 0" : "| 1";
    var to = g.toGaps();
    if (to.length) s += " | " + to.map(function (x) { return x + " "; }).join("");
    return s + "]";
  }

  /* One printed gap set (ProjectPanel line: Double.toString(time) + each Gap.toString() + " "). */
  function formatGapSet(gaps, withTime, time) {
    var body = gaps.map(function (g) { return formatGap(g) + " "; }).join("");
    if (withTime === false) return body;
    if (time === undefined) time = gaps.length ? gaps[0].relativeTime : 0;
    return javaDoubleStr(time) + " " + body;
  }

  function vertexInGap(iv, g, n) {
    var s = g.startEdge, e = g.endEdge;
    return (s < iv && iv <= e) || (e < s && ((s < iv && iv < n) || (0 <= iv && iv <= e)));
  }

  function lineInGap(e, g, n) {
    var fs = g.fullStartEdge(n), fe = g.fullEndEdge(n);
    return (fs <= e && e <= fe) || (fe < fs && ((fs <= e && e < n) || (0 <= e && e <= fe)));
  }

  /* Algorithm.samePhysicalGap with Java's wrap test (B4), as collapsePhysicalGaps uses it. */
  function samePhysicalGapJava(a, b, n) {
    var s0 = a.startEdge, s1 = b.startEdge, e0 = a.endEdge, e1 = b.endEdge;
    return (Math.abs(s0 - s1) <= 1 || s0 + s1 === n - 1) && (Math.abs(e0 - e1) <= 1 || e0 + e1 === n - 1);
  }

  /* A Java runtime exception of the literal bookkeeping (where getGaps would throw). */
  function JavaThrow(exc, line, message) {
    this.exc = exc;
    this.line = line;
    this.message = message || "";
  }

  function IDGenerator(last) {
    this.last = last || 0;
  }

  IDGenerator.prototype.nextId = function () {
    this.last += 1;
    return this.last;
  };

  /* The "invariant edge" per ID (collapsePhysicalGaps lines 396-454, quirk B5). */
  function iedges(pgs, n) {
    var num = pgs[0].length, cse = new Map(), cee = new Map(), i, j, g, c;
    pgs[0].forEach(function (g0) {
      cse.set(g0.id, g0.fullStartEdge(n));
      cee.set(g0.id, g0.fullEndEdge(n));
    });
    for (i = 1; i < pgs.length; i++) {
      for (j = 0; j < num; j++) {
        if (j >= pgs[i].length) throw new JavaThrow("AIOOBE", 411);
        g = pgs[i][j];
        if (!cse.has(g.id)) throw new JavaThrow("NPE", 411);
        c = cse.get(g.id);
        var se = g.fullStartEdge(n);
        if (se !== c) {
          if (se > c) {
            if (se === c + 1) c = se;
          } else if (se + c === n) c = se;
        }
        cse.set(g.id, c);
        c = cee.get(g.id);
        var ee = g.fullEndEdge(n);
        if (ee !== c) {
          if (ee > c) {
            if (ee + c === n) c = ee;
          } else if (ee + 1 === c) c = ee;
        }
        cee.set(g.id, c);
      }
    }
    var out = new Map();
    cse.forEach(function (s, gid) {
      var e = cee.get(gid);
      out.set(gid, e >= s ? Math.floor((s + e) / 2) : 0);
    });
    return out;
  }

  /* Algorithm.collapsePhysicalGaps (literal): IDs propagate by the cyclic shift of the last
   * samePhysicalGap match of pgs[k][0]; then iEdges; returns the last sample. */
  function collapsePhysicalGaps(pgs, idGen, n) {
    var num = pgs[0].length, i, j, k;
    if (num === 0) return [];
    if (pgs[0][0].id === -1) {
      if (!idGen) throw new JavaThrow("NPE", 371);
      pgs[0].forEach(function (g) { g.id = idGen.nextId(); });
    }
    for (i = 0; i + 1 < pgs.length; i++) {
      var pg0 = pgs[i][0], cur = pgs[i], nxt = pgs[i + 1];
      for (j = 0; j < num; j++) {
        if (j >= nxt.length) throw new JavaThrow("AIOOBE", 380);
        if (samePhysicalGapJava(pg0, nxt[j], n)) {
          for (k = 0; k < num; k++) {
            var t = k + j;
            if (t >= num) t -= num;
            if (t >= nxt.length) throw new JavaThrow("AIOOBE", 390);
            nxt[t].id = cur[k].id;
          }
        }
      }
    }
    var ce = iedges(pgs, n), last = pgs[pgs.length - 1];
    last.forEach(function (g) {
      if (!ce.has(g.id)) throw new JavaThrow("NPE", 460);
      g.iedge = ce.get(g.id);
    });
    return last;
  }

  /* ---- sampling and gap tracking ------------------------------------------ */

  /* Critical points of the path and the physical gaps just after each (getGaps lines 76-94), with the
   * compat="fixed" repairs (gaps.sample_path): a sample that Java's perturbation puts outside the
   * polygon (past a waypoint or path end within PERTURB * sqrt(2) of a wall, where Java's scan sees no
   * gaps) is pulled back along the same direction (moved[i] keeps Java's point), and the gaps are
   * repaired (repairGaps).  {samples}.physical[i] are the fixed gaps at points[i]. */
  function samplePath(poly, path, opts) {
    var eps = (opts && opts.eps) || EPSILON, purturb = (opts && opts.purturb) || PERTURB;
    var cuts = (opts && opts.cuts) || allCuts(poly, eps);
    var cps = allCriticalPoints(path, cuts, eps), points = [], physical = [], moved = {};
    cps.forEach(function (cp, i) {
      var p = purturbPointAlongSeg(cp.point, cp.seg, purturb, true, eps);
      var raw = physicalGaps(poly, p, { eps: eps, java: true }), status = pointStatus(poly, p, eps);
      if (status === "outside") {
        moved[i] = p;
        p = insideSample(poly, cp, p, eps);
        status = pointStatus(poly, p, eps);
        raw = physicalGaps(poly, p, { eps: eps, java: true });
      }
      points.push(p);
      physical.push(status !== "outside" ? repairGaps(poly, p, raw, status === "boundary") : raw);
    });
    return { polygon: poly, path: path, cuts: cuts, criticalPoints: cps, points: points, physical: physical, moved: moved };
  }

  /* A sample on the ray from the critical point through p (outside the polygon) that is inside it: p
   * pulled back by halving, or the critical point itself (gaps._inside_sample). */
  function insideSample(poly, cp, p, eps) {
    var x0 = cp.point[0], y0 = cp.point[1], dx = p[0] - x0, dy = p[1] - y0;
    for (var k = 1; k < 40; k++) {
      var f = Math.pow(0.5, k), c = [x0 + f * dx, y0 + f * dy];
      if (pointStatus(poly, c, eps) === "inside") return c;
    }
    return [x0, y0];
  }

  /* The event branches of getGaps (lines 102-290) at critical point i, literally. */
  function eventStep(cp, i, samples, previous, pg, n, edges, idGen) {
    var cps = samples.criticalPoints, cut = cp.cut, rest, g, k;
    if (cp.cutType === CUT.GENERAL_INFLECTION) {
      if (i === 0) throw new JavaThrow("AIOOBE", 107);
      var last = cps[i - 1].point;
      if (rccwLine(edges[cut.fromLine], last) === -1) { /* a gap appears */
        rest = [];
        pg.forEach(function (x) {
          if (lineInGap(cut.fromLine, x, n)) x.id = idGen.nextId();
          else rest.push(x);
        });
        collapsePhysicalGaps([previous, rest], null, n);
      } else { /* a gap disappears */
        rest = previous.filter(function (x) { return !lineInGap(cut.fromLine, x, n); });
        collapsePhysicalGaps([rest, pg], null, n);
      }
    } else if (cp.cutType === CUT.BITANGENT) {
      if (i === 0) throw new JavaThrow("AIOOBE", 174);
      var lastS = samples.points[i - 1];
      var prod = rccwLine(cut.line, cut.curveToPoint) * rccwLine(cut.line, lastS);
      var itp = cut.thisPoint, iop = cut.oppositePoint;
      var four = [itp === 0 ? n - 1 : itp - 1, itp, iop === 0 ? n - 1 : iop - 1, iop];
      var inFour = function (x) { return four.some(function (e) { return lineInGap(e, x, n); }); };
      var old = [], nw = [];
      if (prod === 1) { /* same side as v_{this+1}: a gap splits */
        var parent = null;
        previous.forEach(function (x) {
          if (vertexInGap(iop, x, n)) parent = x;
          else old.push(x);
        });
        for (k = 0; k < pg.length; k++) {
          g = pg[k];
          if (inFour(g)) {
            g.id = idGen.nextId();
            if (parent === null) throw new JavaThrow("NPE", 220);
            parent.addGap(g.id);
          } else {
            nw.push(g);
          }
        }
        collapsePhysicalGaps([old, nw], null, n);
      } else { /* gaps merge */
        var parents = [];
        previous.forEach(function (x) { (inFour(x) ? parents : old).push(x); });
        pg.forEach(function (x) {
          if (vertexInGap(iop, x, n)) {
            x.id = idGen.nextId();
            parents.forEach(function (p) { p.addGap(x.id); });
          } else {
            nw.push(x);
          }
        });
        collapsePhysicalGaps([old, nw], null, n);
      }
    }
  }

  function isEventType(t) {
    return t === CUT.GENERAL_INFLECTION || t === CUT.BITANGENT;
  }

  /* Java's step i-1 -> i on copies: {ids, links: [[index in prev, [new IDs in insertion order]], ...],
   * last} or null where Java throws. */
  function javaProposal(cp, i, samples, prev, n, edges, lastId) {
    var a = prev.map(function (g) { return new Gap(g.physical, g.id); });
    var b = samples.physical[i].map(function (p) { return new Gap(p); });
    var gen = new IDGenerator(lastId);
    try {
      if (isEventType(cp.cutType)) eventStep(cp, i, samples, a, b, n, edges, gen);
      else collapsePhysicalGaps([a, b], gen, n);
    } catch (ex) {
      if (ex instanceof JavaThrow) return null;
      throw ex;
    }
    var links = [];
    a.forEach(function (g, j) { if (g.to.size) links.push([j, g.to.items.slice()]); });
    return { ids: b.map(function (g) { return g.id; }), links: links, last: gen.last };
  }

  function perimeter(poly) {
    var cum = poly._cache.get("perimeter");
    if (!cum) {
      cum = [0];
      for (var i = 0; i < poly.n; i++) {
        var s = poly.edges[i];
        cum.push(cum[i] + distance(s[0], s[1], s[2], s[3]));
      }
      poly._cache.set("perimeter", cum);
    }
    return cum;
  }

  /* The boundary chain hidden behind a gap, as [start, length] in boundary arc length. */
  function hiddenChain(poly, g) {
    var cum = perimeter(poly), per = cum[poly.n], v = poly.vertices, a, b;
    if (g.endPoint !== null) {
      a = cum[g.startEdge];
      b = cum[g.endEdge] + distance(v[g.endEdge][0], v[g.endEdge][1], g.endPoint[0], g.endPoint[1]);
    } else {
      a = cum[g.startEdge] + distance(v[g.startEdge][0], v[g.startEdge][1], g.startPoint[0], g.startPoint[1]);
      b = cum[g.endEdge + 1];
    }
    return [a, pyMod(b - a, per)];
  }

  /* Midpoint (boundary arc length) of the fully hidden edge at the gap's anchor vertex. */
  function core(poly, g) {
    var cum = perimeter(poly), e = g.endPoint !== null ? g.startEdge : g.endEdge;
    return 0.5 * (cum[e] + cum[e + 1]);
  }

  function insideChain(x, chain, per) {
    var d = pyMod(x - chain[0], per);
    return 0 < d && d < chain[1];
  }

  /* Pairs [j, k] (gap j of a, gap k of b) whose core lies in the other's hidden chain: continuations,
   * parent/child of a split or merge (gaps.chain_graph). */
  function chainGraph(poly, a, b) {
    var per = perimeter(poly)[poly.n];
    var ca = a.map(function (g) { return hiddenChain(poly, g); }), cb = b.map(function (g) { return hiddenChain(poly, g); });
    var ka = a.map(function (g) { return core(poly, g); }), kb = b.map(function (g) { return core(poly, g); });
    var out = [];
    for (var j = 0; j < a.length; j++) {
      for (var k = 0; k < b.length; k++) {
        if (insideChain(kb[k], ca[j], per) || insideChain(ka[j], cb[k], per)) out.push([j, k]);
      }
    }
    return out;
  }

  /* Connected components of a bipartite graph: [[sorted a-nodes], [sorted b-nodes]] ordered by their
   * smallest node (a nodes first, then b nodes). */
  function components(na, nb, pairs) {
    var parent = [], x;
    for (x = 0; x < na + nb; x++) parent.push(x);
    function find(y) {
      while (parent[y] !== y) {
        parent[y] = parent[parent[y]];
        y = parent[y];
      }
      return y;
    }
    pairs.forEach(function (p) {
      var ra = find(p[0]), rb = find(na + p[1]);
      if (ra !== rb) parent[Math.max(ra, rb)] = Math.min(ra, rb);
    });
    var comps = new Map();
    for (x = 0; x < na + nb; x++) {
      var r = find(x);
      if (!comps.has(r)) comps.set(r, [[], []]);
      if (x < na) comps.get(r)[0].push(x);
      else comps.get(r)[1].push(x - na);
    }
    return Array.from(comps.keys()).sort(function (p, q) { return p - q; }).map(function (r) { return comps.get(r); });
  }

  function pairKey(j, k) {
    return j + "," + k;
  }

  /* Why Java's step is not kept, or null (gaps._proposal_problem). */
  function proposalProblem(prev, ids, links, lastId, newLast, graph) {
    if (ids.indexOf(-1) >= 0) return "Java leaves a gap without ID";
    var idSet = new Set(ids);
    if (idSet.size !== ids.length) return "Java gives two gaps the same ID";
    var prevIds = prev.map(function (g) { return g.id; }), prevSet = new Set(prevIds), x;
    var fresh = [];
    idSet.forEach(function (v) { if (!prevSet.has(v)) fresh.push(v); });
    if (fresh.length !== newLast - lastId || fresh.some(function (v) { return !(v > lastId && v <= newLast); })) {
      return "Java's new IDs are inconsistent";
    }
    var pos = new Map(), linkOf = new Map(), pairs = new Map();
    ids.forEach(function (gid, k) { pos.set(gid, k); });
    links.forEach(function (l) { linkOf.set(l[0], l[1]); });
    for (var j = 0; j < prevIds.length; j++) {
      var gid = prevIds[j];
      if (pos.has(gid)) {
        if (linkOf.has(j)) return "a linked gap survives";
        pairs.set(pairKey(j, pos.get(gid)), [j, pos.get(gid)]);
      }
      var cs = linkOf.get(j) || [];
      for (x = 0; x < cs.length; x++) {
        if (!pos.has(cs[x])) return "a link points nowhere";
        pairs.set(pairKey(j, pos.get(cs[x])), [j, pos.get(cs[x])]);
      }
    }
    var changed = Array.from(pairs.values()).filter(function (p) { return prevIds[p[0]] !== ids[p[1]]; });
    var comps = components(prev.length, ids.length, changed);
    for (x = 0; x < comps.length; x++) {
      var la = comps[x][0].length, lb = comps[x][1].length;
      if (la && lb && !((la === 1 && lb === 2) || (la === 2 && lb === 1))) {
        return "Java links " + la + " gap(s) to " + lb;
      }
    }
    if (graph !== null) {
      if (graph.length !== pairs.size || graph.some(function (p) { return !pairs.has(pairKey(p[0], p[1])); })) {
        return "Java's step disagrees with the hidden chains";
      }
    }
    return null;
  }

  /* IDs and links of a step read off the chain graph (gaps._from_graph). */
  function fromGraph(prev, nb, graph, lastId) {
    var ids = [], k, comps = components(prev.length, nb, graph);
    for (k = 0; k < nb; k++) ids.push(-1);
    comps.forEach(function (c) {
      if (c[0].length === 1 && c[1].length === 1) ids[c[1][0]] = prev[c[0][0]].id;
    });
    for (k = 0; k < nb; k++) {
      if (ids[k] === -1) {
        lastId += 1;
        ids[k] = lastId;
      }
    }
    var links = [];
    comps.forEach(function (c) {
      if (c[0].length && c[1].length && !(c[0].length === 1 && c[1].length === 1)) {
        var ch = c[1].map(function (kk) { return ids[kk]; });
        c[0].forEach(function (j) { links.push([j, ch]); });
      }
    });
    return { ids: ids, links: links, last: lastId };
  }

  function describe(prev, ids, links, lastId) {
    var prevIds = prev.map(function (g) { return g.id; }), linked = new Set(), linkedJ = new Set();
    links.forEach(function (l) {
      linkedJ.add(l[0]);
      l[1].forEach(function (c) { linked.add(c); });
    });
    var parts = links.map(function (l) { return prevIds[l[0]] + " -> [" + l[1].join(", ") + "]"; });
    ids.forEach(function (g) { if (g > lastId && !linked.has(g)) parts.push("appear " + g); });
    prevIds.forEach(function (g, j) { if (ids.indexOf(g) < 0 && !linkedJ.has(j)) parts.push("disappear " + g); });
    return "from the hidden chains: " + (parts.length ? parts.join(", ") : "no event, IDs carried over");
  }

  /* path.ptDistFromStart(point) / length, keeping only the segments whose arc length agrees with the
   * critical point's own distance (B17 fixed). */
  function relativeTime(path, cp, length, eps) {
    var p = cp.point, d = NaN, tol = 1e-6 * Math.max(length, 1.0);
    for (var i = 0; i < path.segments.length; i++) {
      var s = path.segments[i];
      if (ptSegDist(s[0], s[1], s[2], s[3], p[0], p[1]) < eps) {
        var di = path.ptDistOnSegment(p, i);
        if (Math.abs(di - cp.distance) <= tol) d = di;
      }
    }
    if (d !== d) d = cp.distance;
    return d / length;
  }

  /* Algorithm.getGaps(poly, path) with the Python compat="fixed" bookkeeping: every step is first done
   * Java's way and kept when it is consistent; otherwise (Java throws, or its step is inconsistent) it is
   * read off the chain graph.  This reproduces the Java history on every run where Java does not throw
   * (up to the B17 set times on a self-crossing or retracing path).  {javaMatching: false} also
   * re-derives the consistent steps that disagree with the hidden chains (Java's label rotations, quirk
   * B3), so that every ID stays on its physical shadow (what the simulator uses).  Returns {sets, times,
   * sampleIndices, sampleIds, samples, inferred: [[criticalIndex, why], ...]}; sets[k] is an array of Gap. */
  function getGaps(poly, path, opts) {
    opts = opts || {};
    var eps = opts.eps || EPSILON, javaMatching = opts.javaMatching === undefined ? true : !!opts.javaMatching;
    var samples = opts.samples || samplePath(poly, path, opts);
    var n = poly.n, edges = poly.edges, cps = samples.criticalPoints, length = path.length();
    var out = { sets: [], times: [], sampleIndices: [], sampleIds: [], samples: samples, inferred: [], compat: "fixed" };
    var lastId = 0;
    var first = samples.physical[0].map(function (p) { return new Gap(p, ++lastId); });
    var group = [first], groupStart = 0, sampleGaps = [first];

    function close(end) {
      var last = group[group.length - 1], ce = iedges(group, n);
      last.forEach(function (g) { g.iedge = ce.get(g.id); });
      var t = relativeTime(path, cps[groupStart], length, eps);
      if (last.length) last[0].relativeTime = t;
      out.sets.push(last);
      out.times.push(t);
      out.sampleIndices.push(end);
    }

    for (var i = 1; i < cps.length; i++) {
      var cp = cps[i], prev = group[group.length - 1], phys = samples.physical[i];
      var graph = chainGraph(poly, prev.map(function (g) { return g.physical; }), phys);
      var prop = javaProposal(cp, i, samples, prev, n, edges, lastId), problem, ids, links, newLast;
      if (prop === null) {
        problem = "Java throws";
      } else {
        ids = prop.ids;
        links = prop.links;
        newLast = prop.last;
        problem = proposalProblem(prev, ids, links, lastId, newLast, javaMatching ? null : graph);
      }
      if (problem !== null) {
        var fg = fromGraph(prev, phys.length, graph, lastId);
        ids = fg.ids;
        links = fg.links;
        newLast = fg.last;
        out.inferred.push([i, problem + "; " + describe(prev, ids, links, lastId)]);
      }
      var cur = phys.map(function (p, k) { return new Gap(p, ids[k]); });
      var idSet = new Set(ids);
      var event = newLast > lastId || prev.some(function (g) { return !idSet.has(g.id); });
      if (event) {
        links.forEach(function (l) { l[1].forEach(function (c) { prev[l[0]].addGap(c); }); });
        close(i - 1);
        group = [cur];
        groupStart = i;
      } else {
        group.push(cur);
      }
      lastId = newLast;
      sampleGaps.push(cur);
    }
    close(cps.length - 1);
    out.sampleIds = sampleGaps.map(function (gs) { return gs.map(function (g) { return g.id; }); });
    return out;
  }

  function historyMaxId(sets) {
    var m = 0;
    sets.forEach(function (gs) { gs.forEach(function (g) { if (g.id > m) m = g.id; }); });
    return m;
  }

  /* Printed gap sets of a history (ProjectPanel lines). */
  function historyLines(h, withTime) {
    return h.sets.map(function (gs, k) { return formatGapSet(gs, withTime, h.times[k]); });
  }

  /* ---- events (events.py) -------------------------------------------------- */

  /* The component events between gap sets k and k+1, for every k (events.transition_events): splits
   * (children in Java HashSet order), merges, appears, disappears; components beyond 1-2 / 2-1 become
   * chains of binary events through fresh labels above every ID of the history. */
  function transitionEvents(sets) {
    var fresh = historyMaxId(sets), out = [];
    function nextFresh() {
      fresh += 1;
      return fresh;
    }
    for (var k = 0; k + 1 < sets.length; k++) {
      var cur = sets[k], nxt = sets[k + 1];
      var curIds = cur.map(function (g) { return g.id; }), nxtIds = nxt.map(function (g) { return g.id; });
      var nxtSet = new Set(nxtIds), comps = [], childComp = new Map();
      cur.forEach(function (g) {
        var to = g.toGaps();
        if (!to.length) return;
        if (nxtSet.has(g.id)) throw new Error("set " + k + ": linked gap " + g.id + " survives into set " + (k + 1));
        to.forEach(function (c) {
          if (!nxtSet.has(c)) throw new Error("set " + k + ": gap " + g.id + " links to " + c + ", which is not in set " + (k + 1));
        });
        var hit = Array.from(new Set(to.filter(function (c) { return childComp.has(c); })
          .map(function (c) { return childComp.get(c); }))).sort(function (p, q) { return p - q; });
        var ci;
        if (hit.length) {
          ci = hit[0];
          hit.slice(1).forEach(function (other) {
            var ps = comps[other][0], cs = comps[other][1];
            Array.prototype.push.apply(comps[ci][0], ps);
            cs.forEach(function (c) { if (comps[ci][1].indexOf(c) < 0) comps[ci][1].push(c); });
            cs.forEach(function (c) { childComp.set(c, ci); });
            comps[other] = [[], []];
          });
        } else {
          ci = comps.length;
          comps.push([[], []]);
        }
        comps[ci][0].push(g.id);
        to.forEach(function (c) {
          if (comps[ci][1].indexOf(c) < 0) comps[ci][1].push(c);
          childComp.set(c, ci);
        });
      });
      var evs = [], pos = new Map();
      curIds.forEach(function (gid, j) { pos.set(gid, j); });
      var minPos = function (c) { return Math.min.apply(null, c[0].map(function (p) { return pos.get(p); })); };
      comps.filter(function (c) { return c[0].length; })
        .map(function (c, idx) { return [minPos(c), idx, c]; })
        .sort(function (p, q) { return p[0] - q[0] || p[1] - q[1]; })
        .forEach(function (item) {
          var parents = item[2][0].slice().sort(function (p, q) { return pos.get(p) - pos.get(q); });
          var children = item[2][1];
          if (parents.length === 1 && children.length === 1) {
            throw new Error("set " + k + ": gap " + parents[0] + " is linked to " + children[0] + " alone (not an event)");
          }
          var s = parents[0], t, x;
          if (children.length === 1) {
            for (x = 1; x < parents.length - 1; x++) {
              t = nextFresh();
              evs.push({ type: "merge", a: s, b: parents[x], s: t });
              s = t;
            }
            evs.push({ type: "merge", a: s, b: parents[parents.length - 1], s: children[0] });
            return;
          }
          for (x = 1; x < parents.length; x++) {
            t = nextFresh();
            evs.push({ type: "merge", a: s, b: parents[x], s: t });
            s = t;
          }
          for (x = 0; x < children.length - 2; x++) {
            t = nextFresh();
            evs.push({ type: "split", s: s, a: children[x], b: t });
            s = t;
          }
          evs.push({ type: "split", s: s, a: children[children.length - 2], b: children[children.length - 1] });
        });
      var curSet = new Set(curIds);
      nxtIds.forEach(function (g) {
        if (!curSet.has(g) && !childComp.has(g)) evs.push({ type: "appear", s: g, lo: 0, hi: 0 });
      });
      cur.forEach(function (g) {
        if (!g.to.size && !nxtSet.has(g.id)) evs.push({ type: "disappear", s: g.id, lo: 0, hi: 0 });
      });
      out.push(evs);
    }
    return out;
  }

  /* {initial, events, times}: the component events of a gap history with their relative times. */
  function gapHistoryToEvents(h) {
    var sets = h.sets || h;
    if (!sets.length) return { initial: [], events: [], times: [] };
    var times = h.times || sets.map(function (gs) { return gs.length ? gs[0].relativeTime : 0; });
    var events = [], evTimes = [];
    transitionEvents(sets).forEach(function (evs, k) {
      evs.forEach(function (e) {
        events.push(e);
        evTimes.push(times[k + 1]);
      });
    });
    return { initial: sets[0].map(function (g) { return g.id; }), events: events, times: evTimes };
  }

  /* Never-see-evader labels ("clear" / "contaminated") of every gap of every set (compat="fixed"):
   * all shadows at t0 are contaminated; a gap of set k+1 is contaminated iff it continues a
   * contaminated gap or a contaminated gap is linked to it.  An empty first set (a start that sees no
   * gap) makes every later gap clear (Java returns at once there).  Sets gap.state, returns the states. */
  function neverSeeEvader(h) {
    var sets = h.sets || h;
    if (!sets.length) return [];
    sets.forEach(function (gs) { gs.forEach(function (g) { g.state = null; }); });
    sets[0].forEach(function (g) { g.state = "contaminated"; });
    for (var k = 0; k + 1 < sets.length; k++) {
      var dirty = new Set();
      sets[k].forEach(function (g) { if (g.state !== "clear") dirty.add(g.id); });
      sets[k].forEach(function (g) { if (dirty.has(g.id)) g.toGaps().forEach(function (c) { dirty.add(c); }); });
      sets[k + 1].forEach(function (g) { g.state = dirty.has(g.id) ? "contaminated" : "clear"; });
    }
    return sets.map(function (gs) { return gs.map(function (g) { return g.state; }); });
  }

  /* ---- maps (docs/js/maps_polygon.js) -------------------------------------- */

  var polyCache = new Map();

  /* The original map n (1..14) from SI.polygonMaps, cached. */
  function loadPolygon(n) {
    if (polyCache.has(n)) return polyCache.get(n);
    var maps = SI.polygonMaps;
    if (!maps || !maps.polygons[String(n)]) throw new Error("unknown polygon " + n + " (load docs/js/maps_polygon.js)");
    var p = new Polygon(maps.polygons[String(n)].vertices, String(n));
    polyCache.set(n, p);
    return p;
  }

  /* A stored path record (paths.json "paths") by label, e.g. "P14b" or "fig_TRO-Fig15b". */
  function pathRecord(label) {
    var rec = SI.polygonMaps && SI.polygonMaps.paths[label];
    if (!rec) throw new Error("unknown path " + JSON.stringify(label));
    return rec;
  }

  /* [polygon, Path] of a stored path. */
  function loadPath(label) {
    var rec = pathRecord(label);
    return [loadPolygon(rec.polygon), new Path(rec.waypoints)];
  }

  SI.polygon = {
    EPSILON: EPSILON,
    PERTURB: PERTURB,
    CUT: CUT,
    pyMod: pyMod,
    relativeCCW: relativeCCW,
    ptSegDist: ptSegDist,
    ptSegDistSq: ptSegDistSq,
    distance: distance,
    getIntersect: getIntersect,
    onExtension: onExtension,
    onReverseExtension: onReverseExtension,
    purturbPointAlongSeg: purturbPointAlongSeg,
    Polygon: Polygon,
    Path: Path,
    PhysicalGap: PhysicalGap,
    segInPolygon: segInPolygon,
    vertexVisibility: vertexVisibility,
    physicalGaps: physicalGaps,
    visibilityPolygon: visibilityPolygon,
    GeometryError: GeometryError,
    pointStatus: pointStatus,
    repairGaps: repairGaps,
    exactGaps: exactGaps,
    nearVertexLine: nearVertexLine,
    crossSign: crossSign,
    orient: orient,
    simTurns: simTurns,
    inCone: inCone,
    closedSegmentInPolygon: closedSegmentInPolygon,
    inflections: inflections,
    generalInflectionCuts: generalInflectionCuts,
    bitangentCuts: bitangentCuts,
    singleTangentCuts: singleTangentCuts,
    allCuts: allCuts,
    allCriticalPoints: allCriticalPoints,
    JavaIntHashSet: JavaIntHashSet,
    javaHashsetOrder: javaHashsetOrder,
    javaDoubleStr: javaDoubleStr,
    Gap: Gap,
    formatGap: formatGap,
    formatGapSet: formatGapSet,
    historyLines: historyLines,
    historyMaxId: historyMaxId,
    vertexInGap: vertexInGap,
    lineInGap: lineInGap,
    samplePath: samplePath,
    hiddenChain: hiddenChain,
    chainGraph: chainGraph,
    components: components,
    fromGraph: fromGraph,
    getGaps: getGaps,
    transitionEvents: transitionEvents,
    gapHistoryToEvents: gapHistoryToEvents,
    neverSeeEvader: neverSeeEvader,
    loadPolygon: loadPolygon,
    pathRecord: pathRecord,
    loadPath: loadPath,
  };
})((globalThis.SI = globalThis.SI || {}));

if (typeof module !== "undefined" && module.exports) module.exports = globalThis.SI;
