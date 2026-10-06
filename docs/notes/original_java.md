# The original Java implementation (2008–2012) and what a port must reproduce

Source: the original Java code by J. Yu. It was cloned outside this repository (`orig/source/`: package
`geometry`, 14 polygons in `polygons/{1..14}.dat`, config in `config/*.properties`).
It compiles with JDK 17 plus `log4j-1.2.12.jar`. Every non-UI file was read in full. The
`pe/ui` panels were read only to see which functions they call, with which inputs.
A headless harness was run on the real classes. It calls the same functions the panels call,
with the panels' polygons and paths, and without AWT windows. Every run result quoted below
comes from that harness, not from reasoning about the code.

The owner's requirement is that "the Python implementation must realize the same function as
the original Java implementation (we don't need the applets)". §1 says what that function is.
§2 is the per-function inventory, with its status in `shadowinfo`. §3 lists the robot paths.
§4 lists the quirks and bugs, and §5 is the porting plan with its parity tests.

---

## 1. What the program can compute (applet/GUI excluded)

The program is a geometric front end for the combinatorial filter. It drives a single
omnidirectional, infinite-range robot along a polyline path inside a **simple polygon**.
It then:

1. **Loads a polygon** from `polygons/N.dat`. The format is one vertex `x y` per line.
   Coordinates lie in `[0,1000]²`, are **y-up**, and run **counter-clockwise** (see §2.1).
2. **Computes the critical lines ("cuts")** of the polygon:
   - **inflection rays** at reflex vertices, which are extensions of a boundary edge. They come in two
     kinds: *general*, which make a gap appear or disappear, and *non-general*, which make a gap jump
     to another vertex;
   - **bitangent rays** through pairs of mutually visible reflex vertices, which make gaps split or merge;
   - **single-tangent rays** through a reflex vertex and a vertex where the line crosses the
     boundary, which change gap bookkeeping only.
3. **Computes visibility** from a point: the visibility polygon, and the list of
   *physical gaps*. Each gap is a cyclic interval of hidden boundary edges, plus the point where
   the gap edge meets the boundary.
4. **Tracks gaps (shadow components) along the path.** It finds where the path crosses the cuts
   (sorted by arc length) and samples the gaps just after each crossing. It then propagates
   integer gap IDs and turns crossings of general inflections and bitangents into
   **appear / disappear / split / merge** component events. The output is
   `Gap[][]`: one gap set per inter-event interval, with IDs, a "relative time"
   (arc-length fraction) and successor links (`toGapSet`).
5. **Pursuit–evasion where the evader is never seen** (`NEVER_SEE_EVADER`). It propagates a
   clear/contaminated label for each gap through the gap history.
6. **Simulates targets with an oracle.** `SingleTypeAgentOracle(N)` randomly distributes `N`
   indistinguishable targets as the gaps evolve. This yields the observations: the initial
   hidden total, the counts entering newly appeared gaps, and the counts revealed by
   disappearing gaps.
7. **Runs the combinatorial filter.** It prints the event list and the linear equations
   (Sec. IV-B, Eq. 6 style). It builds the **bipartite I-state** (Sec. V-B, Fig. 9) and
   computes **per-shadow lower and upper bounds** for every final shadow with an
   Edmonds–Karp **max-flow** (Sec. V-C). It prints these as text.

The code does **not** contain multiple teams or attributes (Sec. V-F), FOV events (Sec. V-D; the
hooks exist but are commented out), anything probabilistic (Sec. VI), polygons with holes, the
LP baseline, or initial bounds other than "the hidden total is known exactly".

**Which paper figures this code produced** (checked against the harness output):

| run | reproduces |
|---|---|
| polygon 12 + `ProjectPanel` path, oracle(100) | **ICRA'08 Fig. 2/4** exactly: labels 1–19, split 1→6,7; 6→8,9; 4+10→11; 3→13,14; 7+2→15; 11+5→16; 8+16→17; 17+12→19; appear 5,10,12,18; disappear 9,14; final 13,15,18,19. Also **T-RO Fig. 11(b)**, where T-RO swaps labels 6↔7 and draws s5 as initial. In Java, s5 appears at t = 0.056. The figure renders the map y-down. |
| polygon 13 + commented "polygon 13" path (first 9 points), oracle(100) | **T-RO Fig. 15(a) / ICRA'08 Fig. 8(a)**: 85 component events, 18 final shadows get bounds. The printed BG has 12 left + 29 right map entries = 41 (paper: "41 vertices"; appeared shadows with a self-loop are counted in both maps), but 80 edges (paper: 60). |
| polygon 14 + `ProjectPanel5` path, oracle(1 000 000) | **T-RO Fig. 15(b) / ICRA'08 Fig. 8(b)**: 385 component events, 491 shadow IDs, 12 final shadows get bounds, 339 BG edges (paper: 339). The printed BG has 51 left + 63 right map entries = 114, against 124 in the paper. 124 = 114 − 1 + 11 is what you get if the 11 initial gaps each had their own left vertex (the commented-out "general case" code): `shadowinfo.BipartiteIState`, which has one left vertex per initial shadow, has exactly 61 + 63 = 124 vertices on this run (and 396 edges), so the paper mixes the two counts (`examples/paper_fig15.py`). The figure renders the map y-down. "Five teams" is not in this code. |

---

## 2. Function-by-function inventory

Notation: `n` = number of polygon vertices. Edge `i` is `v_i → v_{i+1 mod n}`, so
`lines[i].P1 = v_i`. `ε = Algorithm.epsilon = 5e-5` (absolute, in map units).
`δ = Algorithm.purturb = 0.005`.
"Status" names the Python counterpart in `shadowinfo.polygon` (the port is complete; the modules are listed in §5.1, and §6 summarises the result). Before the port, all polygon geometry was missing; `shadowinfo` only had the grid-world analogue (`shadowinfo.grid`).

### 2.1 Conventions and Java library semantics the port must copy

| item | semantics |
|---|---|
| Coordinates | Math coordinates, **y-up**. `Polygon.draw` and the mouse handler flip with `Y_MAX − y`, so what the Java GUI shows is the y-up picture. The published figures are not consistent: ICRA'08 Fig. 2a and T-RO Figs. 11(a) and 15(b) are y-down renderings, while T-RO Fig. 15(a) is y-up. **Do not flip the data.** |
| Orientation | All 14 files have positive shoelace area in y-up (CCW). The reflex test only works for this orientation. |
| `Line2D.relativeCCW(A→B, P)` | `ccw = (P−A)×'(B−A)` computed as `px*y2 − py*x2`. In y-up terms it returns **+1 when P is to the right** of A→B and −1 when P is to the left. If `ccw == 0.0` exactly, it returns −1 when P projects before A, +1 when P projects beyond B, and 0 on the segment. **Reflex test:** `relativeCCW(e_{i-1}, v_{i+1}) == 1`. **Convex:** `== −1`. |
| `Line2D.intersectsLine` / `linesIntersect` | Closed test: touching counts. `rccw(a,b,c)*rccw(a,b,d) <= 0 && rccw(c,d,a)*rccw(c,d,b) <= 0`. |
| `Line2D.ptSegDist` | `sqrt(ptSegDistSq)` using the JDK's projection formula. Copy it literally for bit parity. |
| `Point2D.distance` | `Math.sqrt(dx*dx+dy*dy)`, **not** `hypot`. |
| `GeneralPath.contains` | `Polygon.done()` builds a `Path2D.Float`: coordinates are **rounded to float32**, the path is implicitly closed, the rule is `WIND_NON_ZERO`, and the crossing rule is `Curve.pointCrossingsForLine`. That rule is half-open in y (`py < y0 && py < y1` → 0, `py >= y0 && py >= y1` → 0) and needs `px < xintercept`. |
| `epsilonEqual(a,b)` | `|a−b| <= ε`. |
| `getSlope` | `+∞` when `|x1−x2| <= ε`, else `dy/dx`. |
| `getIntersect(l1,l2)` | Intersection of the **infinite lines** in slope–intercept form. If exactly one line is ε-vertical, it evaluates the other line at that x. If both are vertical, or the slopes are ε-equal, it returns null. Note that `i2` uses `line2.P2`, and that near-vertical lines lose precision. |
| `getSegmentIntersect` | `getIntersect` plus `ptSegDist < ε` to **both** segments. |
| `onExtension(line, p)` | `(p.x−x2)(x2−x1) > 0 || (p.y−y2)(y2−y1) > 0`, meaning p is past P2. This is an **OR** of per-axis tests, not a projection. `onReverseExtension` applies the same test to the reversed line. |
| `TreeMap<Double,…>` | Sorted by exact `double`. **On equal keys the last `put` wins.** |
| `HashSet<Integer>` iteration (`toGapSet`) | Bucket order is `(h ^ h>>>16) & (cap−1)` with `cap = 16` (≤ 12 elements). Inside a bucket, insertion order applies. For IDs below 65 536 the order is "ascending by `id mod 16`, ties by insertion". That is why `[6, 10 | 16 15 ]` prints 16 before 15. This decides which split child is `toGap[0]`. |
| `HashMap<Vertex,…>` iteration | Keyed by **identity hash**. HotSpot derives identity hashes from a per-thread xorshift, so the order is repeatable for one JVM and allocation sequence but otherwise arbitrary. It affects BG printing and max-flow augmenting order (see §4, B10). |

### 2.2 Loading and configuration

| Java | semantics | paper | shadowinfo |
|---|---|---|---|
| `Polygon.getPolygon(k, dc)` / `read(BufferedReader)` | Reads `polygons/k.dat` relative to `FileHelper.baseUrl`. Each line is `line.trim().split(" ")`, giving x and y as doubles. It stops at EOF. The `line.equals("\n")` guard can never fire, so a **blank line throws** `NumberFormatException`. Two spaces between the numbers also break parsing. There is no validation. It calls `done()` to build the float `GeneralPath`. | – | `io.load_polygon` / `io.parse_dat` (skips blank lines, B16) |
| `Polygon.write` | Writes `"x y\n"` per vertex (Java `Double.toString`) and then a **blank line**. `read` cannot read that file back. | – | `io.format_dat` / `io.write_dat` (byte-identical; `parse_dat` reads it back) |
| `Polygon.randomPolygonFromData` | Picks `1 + floor(rand·number)`. `polygons.properties: number = 7`, so only files 1–7 are ever picked; 8–14 are reachable only by explicit number. `postfix = .dat`; `prefix` is commented out. | – | `io.random_polygon_from_data(rng)`, `io.NUMBER_OF_POLYGON_DATA_FILES`, `io.DATA_FILE_POSTFIX` |
| `Polygon.randomPolygon` | Random points. These are generally **not** a simple polygon. Unused. | – | `io.random_polygon(k, w, h, rng)` (`validate=False`) |
| `Polygon.getLineSegmentArray` | Edges `i: v_i→v_{i+1}`, cached. | – | `geometry.Polygon.edge(i)` |
| `Polygon.pointInPolygon` | `GeneralPath.contains`. See §2.1. | – | `geometry.Polygon.contains`, `geometry.GeneralPath` |
| `config/geometry.properties` | `scaling-factor 1000` (world size), `canvas-width/height 550`, `canvas-margin 8`, `auto-resize 1`, `minimum-width/height 20`. Drawing only. | – | n/a |
| `config/projects.properties` | `demo-project = Empty Demo Project`. Unused. | – | n/a |
| `FileHelper.*` | Base-URL resource loading (applet code base or jar folder), properties loading, and file write/append/copy helpers. Only `loadUrlProfile` and `getStreamFromUrl` are used by the algorithms. | – | n/a (use package data) |
| `Geometry` (static) | Reads `geometry.properties` and sets `DEFAULT_DC = DrawingContext(1000,1000,0,0)`. | – | n/a |
| `IDGenerator` | Counter. `getNextId()` returns 1, 2, 3, … One instance per `getGaps` call. | – | trivial |
| `Triple` | Generic 3-tuple (BFS back-pointers). | – | n/a |

### 2.3 Geometry primitives (`geometry.Algorithm`)

| Java | semantics / quirks | shadowinfo |
|---|---|---|
| `epsilonEqual`, `getSlope`, `getYIntercept(line)`, `getYIntercept(k,p)` | As in §2.1. With infinite slope, the intercept returns `k` (+∞). | `geometry.epsilon_equal`, `get_slope`, `get_y_intercept` |
| `getXIntercept(line)` | **Bug:** identical to `getYIntercept`. Unused. | – |
| `getPointOnLineFromX(line,x)` | If the line is ε-vertical and `x != x1` it returns `(x, 0)`, a meaningless point (doc says "return null"). If `x == x1` it returns null. Otherwise `(x, kx+b)`. | `geometry.get_point_on_line_from_x` |
| `getPointOnLineFromY(line,y)` | Mirror image of the above. If ε-horizontal it returns `(0, y)` or null. Otherwise `((y−b)/k, y)`. | `geometry.get_point_on_line_from_y` |
| `getIntersect`, `getSegmentIntersect`, `onExtension`, `onReverseExtension` | §2.1. | `geometry.get_intersect`, `get_segment_intersect`, `on_extension`, `on_reverse_extension` |
| `pointOnSegment(p,l)` | `ptSegDist < ε`. | `geometry.point_on_segment` |
| `segmentsOnSameLine` | **Bug:** the second intercept uses `line1.P2` with `k2`. Unused. | – |
| `purturbPointAlongSeg(p, seg, δ, far)` | Moves `p` along `seg` by **δ in x** (`|slope| < 1` or horizontal) or **δ in y** (`|slope| > 1` or vertical), so the step length is between δ and δ√2. It returns whichever of the two candidates is farther from (`far`) or nearer to `seg.P1`. **Bug:** if `|slope| == 1` exactly, both candidates are null and it throws an NPE. None of the panel paths has an exact 45° segment. | `geometry.purturb_point_along_seg` (NPE → `GeometryError` in java mode, x branch in fixed mode) |
| `segInPolygon(poly, seg)` | (1) The midpoint must satisfy `GeneralPath.contains`. (2) Each endpoint must be within ε of some edge, or be contained. (3) No edge may `intersectsLine` the segment at a point farther than ε from both endpoints. A segment that grazes a vertex in its interior therefore counts as **blocked**. Parallel or overlapping edges are ignored, because `getIntersect` returns null for them. | `geometry.seg_in_polygon` (scalar and vectorised) |

### 2.4 Visibility

| Java | semantics | paper | shadowinfo |
|---|---|---|---|
| `getPhysicalGaps(poly, q)` → `PhysicalGap[]` | For each edge `i`: `vis(v) = segInPolygon(q→v)`. If both endpoints are visible, there is no gap on this edge. Otherwise, scan `j` (skipping `j == i` and `j == i+1`). For blocker vertex `v_j = lines[j].P1`, intersect the line `q→v_j` with the line of edge `i`. Keep the hit `ip` only if it is past `v_j`, lies **on** edge `i` (neither `onExtension` nor `onReverseExtension`), and both `q→v_j` and `v_j→ip` pass `segInPolygon`. When `v_i` is hidden, keep the hit nearest to `v_i` as `ip1` (blocker `ip1n`). When `v_{i+1}` is hidden, keep the hit nearest to `v_{i+1}` as `ip2` (blocker `ip2n`). Break as soon as the needed hit(s) are found; `ip1 != ip2` is an *object* identity test. Then jump `i ← ip2n − 1` (skip edges hidden behind the blocker), or end the loop if `ip2n < i`. Emit `PhysicalGap(ip1n, i, null, ip1)`: hidden chain from `v_{ip1n}` ccw to `ip1` on edge `i`. Emit `PhysicalGap(i, ip2n−1 (or n−1 if ip2n == 0), ip2, null)`: hidden chain from `ip2` on edge `i` to `v_{ip2n}`. Output order is ccw by the edge on which the gap's free end lies, starting from edge 0. (`nj` is assigned `i` before it is used, so the `j < i && j > nj` skip never fires.) | Sec. II-A (gaps = shadow-region boundaries) | `visibility.physical_gaps` |
| `PhysicalGap.getFullStartEdge/EndEdge` | First/last edge that is **entirely** hidden: `startEdge+1` if `startPoint != null`, and `endEdge−1` if `endPoint != null`, cyclically. Because `j == i` and `j == i+1` are skipped, every gap has at least one full edge. | – | `visibility.PhysicalGap.full_start_edge/full_end_edge` |
| `getVisibilityPolygon(poly, q, dc)` → `Polygon` | Same scan, emitting vertices per edge: `v_i` if visible, then `ip1`, then `ip2`, then `v_{i+1}` if visible. **Consecutive duplicates** occur, because `v_{i+1}` of edge `i` is re-added as `v_i` of edge `i+1`. Only the GUI calls it, on mouse move. | – | `visibility.visibility_polygon` (duplicates kept with `compat="java"`, B15) |

### 2.5 Cuts (critical lines)

`getCuts(poly)` returns the following lists **concatenated in this order**:
single tangents, non-general inflections (as `Cut(line, NONGENERAL_INFLECTION)`), general
inflections (`GeneralInflection`), and bitangents (`Bitangent`).
The order matters because ties in the critical-point `TreeMap` are last-writer-wins.

| Java | semantics | paper | shadowinfo |
|---|---|---|---|
| `getInflection(poly, type)` → `Line2D[]` | For each **reflex** `v_i` (`rccw(e_{i-1}, v_{i+1}) == 1`): **(a)** "from prevLine" (edge `i−1`), where `v_{i−1}`'s turn decides the type: `rccw(e_{i−2}, v_i) == −1` (convex) → GENERAL, `== 1` (reflex) → NONGENERAL, 0 (collinear) → dropped. The ray goes from `v_i` along the extension of edge `i−1` beyond `v_i` to the nearest boundary hit. "Nearest" means `ptSegDist` to edge `i−1`, over edges other than `i−2`, `i−1`, `i`, with the hit `ptSegDist < ε` on that edge and `onExtension`. **(b)** "from nextLine" (edge `i`), typed by `v_{i+1}`'s turn (`rccw(e_{i+1}, v_i)`); the ray goes from `v_i` along the reverse extension of edge `i`. Excluded edges are `i−1`, `i`, `i+1`. `ALL` keeps both kinds and also collinear ones. If nothing is hit, the result is `Line2D(v_i, null)` and throws an NPE. That never happens in a closed polygon. Output order: by `i`, (a) before (b). | Sec. II-A, inflection crossings ([33] aspect graphs; ICRA'08 refs [7], [13]) | `cuts.inflections` |
| `getGeneralInflectionCut(poly)` → `GeneralInflection[]` | Same as GENERAL, carrying `fromLine` (`i−1` for (a), `i` for (b)) and `counterClockWise` (true for (a), false for (b)). Used by `getGaps` for appear and disappear. | same | `cuts.general_inflection_cuts` |
| `getBitangentCut(poly)` → `Bitangent[]` | For pairs `i < j` with `j >= i+2` (no cyclic wrap, so the pair `(0, n−1)` is examined but always rejected) where **both** `v_i` and `v_j` are reflex. Both pairs of neighbours must lie strictly on the **same side** of the line `v_i v_j` (`rccw` product == 1 at each end). The segment `v_i v_j` must not `intersectsLine` any edge other than the four incident ones. For the line through them, find the nearest boundary hit past `v_j` (`ccwPoint`, `onExtension`) and the nearest on the other side (`cwPoint`, which in practice is before `v_i`). It emits **two** cuts: `Bitangent(v_i→cwPoint, this=i, opp=j, oppSeg=v_j→ccwPoint, curveTo=v_{i+1})` and `Bitangent(v_j→ccwPoint, this=j, opp=i, oppSeg=v_i→cwPoint, curveTo=v_{j+1})`. These are the two rays leaving the bitangent segment beyond each tangent vertex. The `if (p1.distSq(ccw) < p2.distSq(ccw))` branch is the swapped variant and is effectively dead code. | Sec. II-A, bitangent crossings | `cuts.bitangent_cuts` |
| `getBitangent(poly)` → `Line2D[]` | The same rays without metadata. Drawn by the panels. | same | `cuts.bitangent_lines` |
| `getSingletangentCut(poly)` → `Cut[]` | Pairs `i < j`, `j >= i+2`, where **at least one** of `v_i`, `v_j` is reflex. The line must be tangent at exactly one end and cross at the other (`fromP1 * fromP2 == −1`, where `from = rccw(seg, prev) * rccw(seg, next)`). It requires `segInPolygon(v_i v_j)`. The cut is the ray from the tangent vertex away from the other vertex, out to the nearest boundary hit. These are the critical lines where the vertex a gap emanates from changes without a component event. | implicit (visibility cell decomposition, Sec. V-C "[6]") | `cuts.single_tangent_cuts` |
| Cut classes | `Cut(line, CUT_TYPE)`. `GeneralInflection(+fromLine, +counterClockWise)`. `Bitangent(+thisPoint, oppositePoint, oppositeSegment, curveToPoint)`. `PathCutIntersectPoint(point, distance, seg, cut \| cutType)`. `CUT_TYPE ∈ {GENERAL_INFLECTION, NONGENERAL_INFLECTION, SINGLETANGENT, BITANGENT, NONE}`. | – | `cuts.Cut`, `GeneralInflection`, `Bitangent`, `CriticalPoint`, `CutType`; `all_cuts` = `getCuts` |

Cut counts from the harness (golden values):

| poly | n | reflex | single-tangent | non-general inf. | general inf. | bitangent rays |
|---|---|---|---|---|---|---|
| 1 | 31 | 13 | 41 | 14 | 12 | 10 |
| 2 | 68 | 32 | 104 | 34 | 30 | 80 |
| 3 | 18 | 9 | 17 | 10 | 8 | 20 |
| 4 | 124 | 58 | 181 | 70 | 46 | 116 |
| 5 | 253 | 127 | 403 | 124 | 130 | 394 |
| 6 | 18 | 10 | 16 | 12 | 8 | 12 |
| 7 | 22 | 8 | 31 | 6 | 10 | 18 |
| 8 | 30 | 12 | 35 | 12 | 12 | 24 |
| 9 | 33 | 13 | 44 | 12 | 14 | 20 |
| 10 | 54 | 25 | 87 | 30 | 20 | 56 |
| 11 | 45 | 21 | 58 | 24 | 18 | 56 |
| 12 | 34 | 15 | 47 | 18 | 12 | 40 |
| 13 | 140 | 68 | 407 | 74 | 62 | 364 |
| 14 | 267 | 134 | 727 | 78 | 190 | 882 |

### 2.6 Paths and path–cut intersections (`drawable.Path`)

| Java | semantics | shadowinfo |
|---|---|---|
| `addPoint`, `getLineArray`, `getLength` | Polyline. The length is the sum of segment lengths (`Point2D.distance`). | `geometry.Path`, `Path.length` |
| `ptDistFromStart(p)` | Arc length to `p`. It scans **all** segments and keeps the **last** one with `ptSegDist < ε`, so the answer is wrong for self-intersecting paths. None of the panel paths self-intersects. It returns NaN if `p` is on no segment. | `geometry.Path.pt_dist_from_start` (B17 fixed in `get_gaps(compat="fixed")`) |
| `getIntersectPoint(Line2D)`, `pointOnPath(p)` | The crossing of one line with the **first** path segment it crosses (`getSegmentIntersect`), or null; whether `p` is within ε of a segment. Neither is called by the algorithms. | `geometry.Path.intersect_point`, `Path.point_on_path` |
| `getIntersectPoints(Line2D[])`, `getIntersectPoints(Cut[])` | For each path segment `i` and each segment/cut `j`, `getSegmentIntersect`, keyed by `ptDistFromStart(p, lines, i)` = (length of segments before `i`) + `|p − seg_i.P1|`. Collected in a `TreeMap`, so **duplicate distances are dropped (last put wins)**. A crossing at a waypoint can be found from both adjacent segments with keys that differ in the last bits; both then survive. | `geometry.Path.intersect_points` (`Line2D[]`), `cuts.cut_intersect_points` (`Cut[]`) |
| `getAllCriticalPoints(cuts)` | The same, then `put(0.0, start, NONE, seg 0)` and `put(len, end, NONE, last seg)`. These overwrite any cut that crosses exactly there. | `cuts.all_critical_points` (`Path.all_critical_points`) |

### 2.7 Gap tracking along a path: `Algorithm.getGaps(poly, path)` → `Gap[][]`

Algorithm (Sec. II-A component events; Sec. V-C "visibility cell decomposition"):

1. `cuts = getCuts(poly)`. `pcips = path.getAllCriticalPoints(cuts)`, sorted by arc length.
2. For every critical point `i`: `pp = purturbPointAlongSeg(point_i, seg_i, δ, far=true)`.
   This is a point just **after** the crossing (or just past the end point).
   Then `pg = getPhysicalGaps(poly, pp)`, appended to the current group `pgVec`.
3. **Event handling.** Only GENERAL_INFLECTION and BITANGENT crossings are component events.
   NONGENERAL and SINGLETANGENT crossings are resampling points only.
   - *General inflection* (cut `c`, `fromLine = f`): `last = point_{i−1}` (**not** perturbed).
     If `rccw(edge f, last) == −1` (the robot was on the interior side of edge `f`), it is
     an **appear**. The gap in `pg` whose full-edge interval contains `f` (`lineInGap`) gets a
     fresh ID. The other gaps take their IDs from `previousGaps` via
     `collapsePhysicalGaps([prev, rest], null)`. Otherwise it is a **disappear**. The gap of
     `previousGaps` containing `f` is dropped; no link is recorded. The remaining previous gaps
     pass their IDs to `pg`.
   - *Bitangent* (`this = a`, `opp = b`, `curveTo = v_{a+1}`):
     `last = perturb(point_{i−1}, seg_{i−1}, far)`.
     It is a **split** when `rccw(cut, curveTo) * rccw(cut, last) == 1` (same side), and a
     **merge** otherwise (including 0).
     - Split: the previous gap with `vertexInGap(b)` is the parent. Every gap in `pg` that
       contains, as a full edge, one of edges `a−1, a, b−1, b` gets a fresh ID and is added to
       `parent.toGapSet`.
     - Merge: every previous gap containing one of those four edges is a merging parent. The
       gap in `pg` with `vertexInGap(b)` gets a fresh ID, which is added to each parent's
       `toGapSet`.
     - In both cases the remaining gaps pass their IDs by collapsing.
4. `relativeTime` = `ptDistFromStart(point_i)/length`. It is set at `i = 0` and at every event.
5. When the next critical point is an event, or `i` is the last point, the group is closed:
   `previousGaps = collapsePhysicalGaps(pgVec, idGen)`.
   - The returned set is the **last** sample of the group. Its IDs are propagated from the
     group's first sample. On the very first group, the first sample gets fresh IDs `1, 2, …`
     in boundary order.
   - `previousGaps[0].relativeTime` = the event time. **Only element 0 carries the time.**
   - The set is appended to the output, and `pgVec` is cleared.
6. Output: `Gap[][]`, one set per interval (1 + number of GI/BT crossings). Gaps are in ccw
   boundary order. A gap continues under the same ID. `toGapSet` links a set `k` to set `k+1`.

Helpers:

| Java | semantics |
|---|---|
| `collapsePhysicalGaps(pgs[][], idGen, n)` | Assumes **every** sample in `pgs` has `numGaps = |pgs[0]|` gaps in the same cyclic order. For consecutive samples it takes `pg0 = pgs[k][0]` and, for each `j`, if `samePhysicalGap(pg0, pgs[k+1][j])`, copies all IDs with a cyclic shift `j`. Several matches mean **the last one wins**. If there is no match, IDs stay −1 (or stale). Then it computes an "invariant edge" `iEdge` per ID by clamping each ID's full start/end edges across samples: the start may only move +1 or wrap, the end only −1 or wrap. `iEdge = (cse+cee)/2` if `cee >= cse`, else **0**. `iEdge` is written onto the last sample. If the sample counts differ, the result is `ArrayIndexOutOfBounds` or an NPE (see §4, B3). |
| `samePhysicalGap(a,b,n)` | `(|Δstart| <= 1 || sa+sb == n−1) && (|Δend| <= 1 || ea+eb == n−1)`. The wrap test is wrong (§4, B4). |
| `vertexInGap(v, g, n)` | `start < v <= end`, cyclic, using the **raw** start/end edges. |
| `lineInGap(e, g, n)` | `fullStart <= e <= fullEnd`, cyclic. |
| `Gap` | `id` (−1 = unset), `iEdge` (a cosmetic label, only printed), `relativeTime`, `state` (`GapState`), `toGapSet: HashSet<Integer>`. `toString()` = `"[" + id + ", " + iEdge + (state?) + (" | " + each id + " ")? + "]"`. Examples: `[1, 0 | 6 7 ]`, `[19, 0| 1]`. `getGapMap` builds an id→gap map. |
| `PhysicalGap(se, ee, sp, ep)` | Exactly one of `sp`/`ep` is non-null. |

**How the component events are read off `Gap[][]`** (this is what the oracle does, §2.9):
for consecutive sets `k`, `k+1`:
- **split**: a gap in `k` whose `|toGapSet| == 2`;
- **merge**: gaps in `k` whose `toGapSet` is the same singleton;
- **appear**: a set with more gaps and no split, so the leftover new ID appeared;
- **disappear**: a set with fewer or equal gaps, no link, and an ID of `k` missing from `k+1`.

Each set transition is assumed to contain **exactly one** event.

shadowinfo: `gaps.get_gaps` → `GapHistory`, `gaps.collapse_physical_gaps`, `same_physical_gap`, `vertex_in_gap`, `line_in_gap`, `format_gap_sets`; `events.gap_history_to_events` turns the history into `shadowinfo.events`. The grid analogue is `grid.shadows.ShadowTracker.update` → `Transition`,
which derives events from component overlap rather than from cut crossings.
The event types themselves exist in `events.py`.

### 2.8 Gap-history processing: `pe.Algorithm.processGapHistoryInformation(gapss, SCENARIO, data)`

There is only one scenario, `NEVER_SEE_EVADER` → `pursuerNeverSeeEvader(gapss)`
(Sec. V-E, pursuit–evasion as a special case):

- If `gapss` is empty or `gapss[0]` is empty, return. Otherwise every gap of `gapss[0]` becomes `NSEState.CONTAMINATED`.
- For each set `k`, in order:
  - a gap with `state == null` (a new gap) becomes **CLEAR**. Note that in this branch its own links are not followed.
  - `|toGapSet| == 0`: copy the state to the same ID in `k+1`, if that ID is present (a disappearance just drops it).
  - `|toGapSet| == 1` (merge): if this gap is not clear, the target becomes CONTAMINATED. Otherwise the target is untouched and becomes CLEAR when visited if still null.
  - `|toGapSet| == 2` (split): both children copy the state.
  - Larger sets are silently ignored, so the children become CLEAR (§4, B6).
- `NSEState.toString()`: `"| 0"` is clear and `"| 1"` is contaminated.
- `additionalData` is unused.

shadowinfo: `events.never_see_evader` (Java semantics with `compat="java"`; `compat="fixed"` propagates contamination along every link, so a split into more than 2 gaps is handled, B6) and `events.never_see_evader_via_filter`, which uses `counting_filter`, with initial
`[0, ∞)` for each initial gap, `Appear(s,0,0)` and `Disappear(s,0,0)`. A gap is CLEAR ⇔ its upper bound is 0.
(`pursuit_evasion_filter` / `evader_status` add `total=(1,1)`, which is a *different*, stronger
model. Do not use it for parity.)

### 2.9 Oracle: `SingleTypeAgentOracle(N)` (`pe.is.Oracle`, `Event`)

`Oracle` is abstract: `initialize(gapss)`, `getTotalAgentNumber()`, `getNumberOfEvents()`,
`getEvent(i)`, `getEventByTime(t)`. `EVENT_TYPE = {INIT, SPLIT, MERGE, APPEAR, DISAPPEAR,
AGENT_APPEAR, AGENT_DISAPPEAR}`. `Event` only has `getEventRelativeTime()`.

`initialize` → `distributeAgents(gapss)`. One event per gap set, keyed by
`gapss[k][0].relativeTime` in a `TreeMap<Double,…>`, so equal times overwrite. All randomness is
`Math.random()`, which is **unseeded and cannot be reproduced**.

- `k = 0`: `visible = (int)(rand·0.5·N)`. INIT: `distributeOneBatch(N − visible, |gaps|)`
  gives each initial gap a `SingleTypeAgentState(count)`. `toGap` = the initial IDs.
- `|gaps_k| > |gaps_{k−1}|`:
  - each parent with `|toGapSet| == 2` (iteration order) is a **SPLIT**:
    `distributeOneBatch(count, 2)` goes to `toGaps[0]` and `toGaps[1]`, in HashSet order;
  - the other parents copy their state forward;
  - if exactly one ID is left over, it is an **APPEAR**:
    `b = distributeOneBatch(visible, 2)`, the new gap gets `b[0]` (= `movedAgents`), and
    `visible = b[1]`.
- Otherwise:
  - each parent with `|toGapSet| == 1` is a **MERGE**. The first parent seen sets
    `target = count_1`; the second sets `target = count_2 + count_2` (§4, B7 doubles it);
    `fromGap = {first, second}`, `toGap = {target}`;
  - a parent whose ID is missing from `k` is a **DISAPPEAR**: `moved = count` and
    `visible += count`;
  - other parents copy their state forward.
- `distributeOneBatch(t, g)`: if `g == 1`, return `{t}`; if `t == 0`, return zeros. Otherwise draw
  `g−1` **distinct** integers in `[0, t)` until there are enough, add `t`, sort, and take the
  differences. The last part is therefore always ≥ 1. If `t < g−1` the loop **never ends**.
- `addRandomVisibilityEvents` (AGENT_APPEAR/AGENT_DISAPPEAR, the FOV events) exists but is
  **never called**, because the call is commented out. It also has a bug: AGENT_DISAPPEAR
  *sets* the gap count to `deltaA` instead of adding.
- `SingleTypeAgentEvent.toString()`:
  - `INIT EVENT 1 2 3 4`
  - `APPEAR EVENT 5 <- 2` (new gap ← moved)
  - `DISAPPEAR EVENT 14 -> 36`
  - `SPLIT EVENT 1 -> 6 7`
  - `MERGE EVENT 7 2 -> 15`
  - `AGENT_APPEAR: k appeared from gap g` and `AGENT_DISAPPEAR: k disappeared from gap g`.
  An event with a null type throws an NPE in `switch`.

shadowinfo: `oracle.SingleTypeAgentOracle` and `oracle.JavaRandom` (bit-exact `java.util.Random`, incl. `nextGaussian`). Before the port: `grid.simulate.GridSimulator` is a different model: targets do a random
walk. `rng.Rng` (Mulberry32) exists, but Java's `Math.random` is `java.util.Random`'s 48-bit LCG,
and `rng.Rng` does not reproduce it.

### 2.10 Filter: `SingleTypeAgentAlgorithm`, `Equation`, `BipartiteGraph`, `MaxFlow`

| Java | semantics | paper | shadowinfo |
|---|---|---|---|
| `Equation` | A `TreeMap<unknown, coeff>` (sorted by ID) plus an int constant. `addTerm` **overwrites** an existing coefficient. `toString`: `"c xN"` terms joined by `" + "` if the next coefficient is ≥ 0, else `" "`; then `" = const"`. Example: `-1 x1 + 1 x6 + 1 x7 = 0`. | Sec. IV-B (Eq. 6) | `stagent.Equation` (Java printer); also `lp.ShadowLP` |
| `deriveGapEvolvingEquations(gapss, oracle)` | One equation per set, using `oracle.getEventByTime(gapss[k][0].relativeTime)`. k = 0: `Σ x_g = N − visible`. APPEAR: `x_new = moved`. DISAPPEAR: `x_old = moved`. SPLIT: `x_a + x_b − x_s = 0`. MERGE: `x_s − x_a − x_b = 0`. | Sec. IV-B | `stagent.derive_gap_evolving_equations` |
| `deriveShadowInfoState(gapss, oracle)` → `BipartiteGraph` | **All initial gaps share ONE left vertex `0`** with weight `N − visible` and edges to each initial gap. This is the "known hidden total" initial condition; per-gap left vertices are commented out as "for more general case". APPEAR: one `Vertex(id)` is **both left and right**, with a self-loop edge and weight = moved. DISAPPEAR: the right vertex's weight = moved (it stays in the right map). SPLIT: the right vertex is replaced by two, and every in-edge is duplicated. MERGE: two right vertices become one, with in-edges deduplicated by source. A self-loop source becomes an ordinary left vertex once its gap splits or merges. AGENT_* events are no-ops. It prints `bg.toString()`: `Left vertices` / `Vertex: id W: w, goes to: ids…` / `Right vertices` / `Vertex: id W: w` (−1 = unknown, i.e. alive at the end). | Sec. V-B, Fig. 9 | present: `bipartite.BipartiteIState` (incremental reach sets, one left vertex per initial shadow). **Equivalent** with `initial = {g: (0, ∞)}` and `total = (N−v, N−v)`. |
| `deriveShadowBounds(bg)` | See the next paragraph. It returns `{source, sink2}` and **mutates `bg`**. | Sec. V-C, Eqs. 7, 8 | **differs**: `CombinatorialFilter.all_bounds()` is exact (lower-bounded flow) and `bounds_paper()` is the literal Eqs. 7/8. The Java recipe is `stagent.java_bounds(bg, order)` (comparison only); `stagent.exact_bounds` is the default. |
| `MaxFlow.getMaxFlow(source, sink)` | **Edmonds–Karp.** The residual graph is built by BFS from the source; residual capacity = `capacity`. Each round BFS-searches the residual maps (HashMap order) for the sink and pushes the bottleneck. A residual edge whose weight equals the bottleneck is **removed**. Reverse edges are created lazily. Zero-capacity edges *are* in the residual graph: they are "augmented" by 0 and then removed. It ends when no path is found. The flow on an edge is the weight of its reverse residual edge, written into `edge.weight`. Two-way edges are not supported. `getIncrementalMaxFlow` (no rebuild) is unused. `main` builds a 6-vertex test graph and prints nothing. | Sec. V-C, [10] | present: `maxflow.FlowNetwork.max_flow` (EK, deterministic insertion order, plus lower bounds) |
| `BipartiteGraph` | `leftVertexMap`, `rightVertexMap` (TreeMaps by ID), `rightVertexEdgeMap: Vertex → Set<Edge>`. `Vertex{id, visited, minWeight, maxWeight (−1 = unknown), veMap, vreMap}`. `Edge{from, to, weight (flow or residual), capacity}`. `createResidueEdge` = `Edge(from, to, capacity, 0)`. | Fig. 9 | `BipartiteIState.to_dict()` |
| `NSEState`, `SingleTypeAgentState`, `GapState` (two identical interfaces, in `gap` and `gap.state`) | Per-gap state holders. `SingleTypeAgentState` has no `toString`, so printing it gives `Object@hash`. | – | n/a |

**`deriveShadowBounds` in detail.**

1. Every left vertex with a self-loop is **deleted** from both maps. This is a shadow that
   appeared and never split or merged, whether it is still alive or already disappeared.
   Consequence: a final shadow that appeared with a known count gets **no bound printed**.
   For example s18 in Fig. 11 is dropped, although its count, 22, is known exactly.
2. `S → L` has capacity `w_L`. Every `L → R` edge also has capacity `w_L`. Let
   `max = Σ w_L`.
3. A right vertex with a known weight (disappeared) gets an edge `R → sink2` with capacity
   `w_R`. Let `sure = Σ w_R`. A right vertex with an unknown weight (alive) gets `R → sink1`.
   The edge `sink1 → sink2` has capacity `free = max − sure`.
4. For each alive shadow `j`, in ID order:
   - **Upper bound:** set `cap(j → sink1) = free` and all other alive edges to 0. Run
     max-flow and print `Max flow in shadow j is f(j→sink1)`.
   - **Lower bound:** set `cap(j) = 0` and the others to `free`. Run max-flow and print
     `Min flow in shadow j is max − (sure + Σ_{i≠j} f(i→sink1))`.
5. Finally it prints a blank line.

Differences from the paper and from `shadowinfo`:
- **Upper bound vs Eq. 7.** Java reports `f(j,T)` without the term
  `+ Σ_i f(i,T) − Σ_i c(i,T)` over disappeared shadows. The result is too high whenever the max
  flow it happens to find leaves a disappeared edge unsaturated.
- **Lower bound vs Eq. 8.** Java subtracts the disappeared **capacities** (`sure`) rather than the
  actual flows. The result is too low or too high depending on which max flow is found.
- **Lower bounds are never enforced** (each disappeared shadow had exactly `w_R` targets), so the
  printed numbers depend on which maximum flow EK returns. That depends on `HashMap` identity
  order (see §4, B10).
- No feasibility check is made.
- Bounds are only for **single** shadows (no group queries).
- All weights are exact (`l = u`), because the oracle knows the counts.

`shadowinfo.CombinatorialFilter` with `total=(N−v, N−v)` computes the exact bounds.
On the ICRA'08 Fig. 4 run it gives exactly Java's numbers for 13, 15 and 19:
13: [0, 38], 15: [0, 38], 19: [24, 62]. It also gives 18: [22, 22], which Java does not print.

### 2.11 Output formats (what the panels print to stdout)

Each `ProjectPanel*` constructor prints the following. These lines are the program's "output":

- `ProjectPanel` / `ProjectPanel5`: for each set, `relativeTime + " "` followed by each
  `gap.toString() + " "`. The time uses Java `Double.toString`, e.g. `0.0`, `0.05623…`; it
  differs from Python `repr` for exponents (`1.0E-4`).
  Then two blank lines. Then the BG (from `deriveShadowInfoState`). Then the Max/Min lines.
  Then `Computation time: <ms>`. Then each equation. Then a blank line and each event's `toString()`.
- `ProjectPanel2`: the gap sets **after** NSE processing (so each gap shows `| 0`/`| 1`), with
  no times. Then two blank lines.
- `ProjectPanel3`: the same as `ProjectPanel2` but for polygon 5.
- `ProjectPanel4`: prints nothing. It draws only the general inflections and bitangents.
- All panels: on mouse move, `getPhysicalGaps` and `getVisibilityPolygon` at the cursor
  (drawn, not printed).

---

## 3. Robot paths in `ProjectPanel*.java`

Coordinates are y-up, exactly as written in the source. "Active" means the path is used by that
panel; "commented" means it is commented out. The harness labels are used in §5.

| label | polygon | panel(s), state | waypoints |
|---|---|---|---|
| **P12** | 12 | `ProjectPanel` (active; `getGaps`, oracle `N=100`, equations, BG, bounds) | (90,450) (580,450) (580,840) (520,840) |
| **P1** | 1 | `ProjectPanel2` (active; `getGaps` + NSE); `ProjectPanel3` (commented) | (930,886) (913,700) (746,583) (75,565) (136,215) (375,69) (788,27) (905,95) |
| **P5** | 5 | `ProjectPanel3` (active; `getGaps` + NSE). Comment `// The one` on (245,286). | (946,501) (723,591) (670,511) (596,643) (606,930) (580,919) (588,672) (250,611) (266,535) (115,433) (93,358) (245,286) (98,288) (93,123) |
| **P14a** | 14 | `ProjectPanel4` (active, drawing only: GI + bitangent lines, no `getGaps`) | (89.29,625.00) (178.57,517.86) (330.36,446.43) (375.00,357.14) (178.57,321.43) (107.14,241.07) (133.93,142.86) (196.43,285.71) (392.86,321.43) (455.36,151.79) (598.21,133.93) (660.71,187.50) (616.07,196.43) (526.79,169.64) (455.36,383.93) (580.36,366.07) (544.64,437.50) (383.93,437.50) (437.50,625.00) **(633.93,714.29) (678.57,687.50)** |
| **P14b** | 14 | `ProjectPanel5` (active; `getGaps`, oracle `N=1 000 000`, equations, BG, bounds) | as P14a but with the last two waypoints changed to **(633.93,712.29) (688.57,657.50)** |
| **P1s** | 1 | commented, labelled `// polygon 1`, in `ProjectPanel`, `ProjectPanel2` and `ProjectPanel3` | (886,281) (886,48) (556,36) |
| **P13** | 13 | commented, labelled `// polygon 13`, in `ProjectPanel` | (195,620) (195,510) (620,510) (620,660) (760,660) (760,470) (590,470) (590,440) (460,440) |
| **P13+** | 13 | same comment block, as `// Second piece of the path`; it continues P13 | P13 + (69,418) (60,575) (356,620). A doubly-commented tail follows: (490,858) (675,876). Note that (69,418) is **outside** polygon 13. |
| **P1L** | 1 | commented in `ProjectPanel` (inside the "polygon 13" block, but it fits polygon 1), `ProjectPanel2` and `ProjectPanel3` | (930,886) (913,700) (746,583) (471,580) (469,780) (665,860) (650,965) (415,786) (295,568) (75,565) (136,215) (375,69) (788,27) (905,95) |
| **P5s** | 5 | commented, labelled `// polygon 5`, in `ProjectPanel` and `ProjectPanel2` (it is a sub-path of P5) | (580,919) (588,672) (250,611) |

The 14-polygon paths have decimals like `89.29 = 50·1000/560`. They were apparently digitised on a
560-px canvas. No path has an exact 45° segment, and no path crosses itself.

Harness results, used as golden data:

| run | critical pts | GI+BT events | gap sets | final gaps | Java outcome |
|---|---|---|---|---|---|
| P12 on 12 | 30 | 14 | 15 | 4 (IDs 13,15,18,19) | ok; = ICRA'08 Fig. 4 |
| P1 on 1 | 49 | 14 | 15 | 1 | ok |
| P5 on 5 | 306 | 174 | 175 | 7 | ok |
| P14a on 14 | 664 | 382 | – | – | **`getGaps` throws an NPE** in `collapsePhysicalGaps` (line 460, reached from the bitangent-merge branch). The panel never calls `getGaps`, so the bug was never seen. |
| P14b on 14 | 671 | 385 | 386 | 14 (12 bounded) | ok; = T-RO Fig. 15(b) |
| P1s on 1 | 9 | 2 | 3 | 3 | ok |
| P13 on 13 | 192 | 85 | 86 | 20 (18 bounded) | ok; = T-RO Fig. 15(a) |
| P13+ on 13 | 300 | 131 | – | – | **`ArrayIndexOutOfBounds`** in `collapsePhysicalGaps` (GI-appear branch); the path leaves the polygon |
| P5s on 5 | 84 | 47 | 48 | 6 | ok |
| P1L on 1 | 73 | 20 | 21 | 1 | ok |

---

## 4. Bugs and quirks (decide: match or fix, and document)

Each item gives its effect, followed by a recommendation for the port: **[match]** reproduce
the Java behaviour exactly in a `java_compat` mode, or **[fix]** use the correct behaviour in
the default mode.

- **B1. `purturbPointAlongSeg` NPE at slope exactly ±1.** Separately, the step is δ in x or in y, not δ of arc length. [match the step rule; fix the NPE by treating |s|==1 as the x branch]
- **B2. Critical points use last-writer-wins on exact distance.** Two cuts crossing the path at the same arc length keep only the later one, in `getCuts` order. The start and end points overwrite cuts at 0 and at the full length. A crossing at a waypoint may appear twice with almost equal keys. [match]
- **B3. `collapsePhysicalGaps` assumes every sample has the same number of gaps, in the same cyclic order, and matches only on `pgs[k][0]`.** When the count changes inside a group, the result is `ArrayIndexOutOfBounds` or, if no match is found, `id = −1` followed by an NPE. Examples: P14a, P13+, or a numerical miss of a crossing. [match in compat mode by raising a typed `GapTrackingError`; the default mode could match gaps by blocker vertex instead]
- **B4. `samePhysicalGap` wrap test `s0 + s1 == n − 1`** accepts any pair of edges whose indices sum to `n−1`, not just `{0, n−1}`. Spurious matches can relabel gaps, because the last match wins. [match in compat; fix as `(s0−s1) mod n ∈ {−1,0,1}`]
- **B5. `iEdge`** comes from ad hoc clamping, the wrap tests `se+cse == n` and `ee+cee == n` are wrong, and it is 0 for wrapping gaps. It is only used in `toString`. [match, for printing only]
- **B6. NSE ignores a split into more than 2 gaps** (the children become CLEAR, which is unsound). New gaps are set CLEAR without following their own links. That shortcut is harmless, because the children of a clear gap default to CLEAR anyway. [fix by treating any `|toGapSet| ≥ 2` as a split; matching is the same on all panel runs]
- **B7. The oracle double-counts merges.** The second merging parent sets `target = count₂ + count₂`, and `count₁` is lost. Target conservation breaks, so later DISAPPEAR counts of the merged lineage (and `visible`) are inconsistent with the initial total. The filter can then receive **infeasible** observations, and Java does not detect it. [fix; keep a `java_merge_bug=True` switch only to replay old runs]
- **B8. `distributeOneBatch(t, g)` loops forever if `t < g−1`.** It is also biased: the last part is ≥ 1 whenever t > 0, so an APPEAR never takes *all* visible targets; with `visible == 1` it always moves 0. [fix the loop; match the distribution in compat mode]
- **B9. Events are keyed by relative time in a `TreeMap`.** Equal times overwrite. `deriveGapEvolvingEquations` looks events up by `gapss[k][0].relativeTime`. A set with no event (equal gap counts and no change) leaves the event type null, and the equations step throws an NPE. [match with an error]
- **B10. `deriveShadowBounds` is not a correct bound computation** (§2.10). Crafted counterexample, run on the real Java classes with an injected oracle. Events: INIT gaps {1,2} with hidden total 1; APPEAR 3 (1 target); SPLIT 3→4,5; MERGE 1+4→6; APPEAR 7 (1); SPLIT 7→8,9; DISAPPEAR 6 (1), 2 (1), 9 (0). The truth (= `shadowinfo`) is s5 ∈ [0,0] and s8 ∈ [1,1]. Over 40 different identity-hash orders Java printed `5:[0,0] 8:[1,1]` 31 times, `5:[0,1] 8:[0,1]` 7 times, and `5:[1,1] 8:[0,0]` twice. The last answer is **wrong**: neither shadow is inside its printed bounds. The output depends on `HashMap<Vertex>` identity order. Within one JVM it is repeatable, but not across allocation sequences. [fix: the default is `CombinatorialFilter` (exact). Provide `java_bounds()` with a fixed, documented order for comparisons only.]
- **B11. Shadows that appeared and never split or merged are deleted before the bounds step.** No bound is printed for a final shadow that appeared, even though its count is known. [fix: report `[k,k]`]
- **B12. Initial condition.** All initial shadows are pooled into one left vertex with an exact total. There are no per-shadow initial bounds (`minWeight` is never used). [shadowinfo supports both; use `total`]
- **B13. FOV events (AGENT_APPEAR/DISAPPEAR) are stubbed out** in both the oracle and the filter, and `addRandomVisibilityEvents` also has a set-instead-of-add bug. [not needed for parity; shadowinfo has `Enter`/`Exit`]
- **B14. Floating-point dependence.** `GeneralPath` uses float32 vertices, while all other tests use double. Absolute `ε = 5e-5` is compared against slope–intercept intersections, which lose precision near vertical lines. `segInPolygon` treats grazing a vertex as blocked. To reproduce the cut and gap lists bit for bit, the port must copy these formulas literally, in the same operation order, with `math.sqrt(dx*dx+dy*dy)` and no `hypot`. [match in the geometry layer]
- **B15. `getVisibilityPolygon` emits duplicate consecutive vertices.** [match in compat; dedupe in default]
- **B16. Loader issues.** A blank line inside a `.dat` file throws, and `Polygon.write` produces files that `read` cannot load. `randomPolygonFromData` only draws from files 1–7. [fix the loader to skip blank lines]
- **B17. `ptDistFromStart` takes the last matching segment.** This is wrong for self-crossing paths. None are used. [fix by using the segment index already known from the intersection]
- **B18. Nondeterminism.** Everything after `getGaps` uses unseeded `Math.random`. [port the oracle with a seedable `java.util.Random` clone, see §5]

---

## 5. Porting plan

### 5.1 New modules (`shadowinfo/polygon/`)

**Status:** done. The table below lists the modules as they were built, with their actual public
names and signatures. The module docstrings are the reference; they cite the Java method for each
function and the paper section it implements.

All modules are pure Python plus numpy, with no new dependencies. Each function name follows
the Java name it mirrors. Every function whose behaviour differs takes `compat=`:
`compat="fixed"` (the default of the public API) applies the [fix] items of §4, and
`compat="java"` (used by the parity tests) keeps the quirks marked [match].

| module | contents | Java source |
|---|---|---|
| `polygon/data/{1..14}.dat` | Copy of the 14 maps as package data (byte-identical), plus `paths.json` holding the §3 table, the digitised figure paths and the provenance of each. | `polygons/` |
| `polygon/io.py` | `load_polygon(n_or_path)` (skips blank lines), `load_path(label)`, `load_paths()`, `path_labels()`, `path_record(label)`, `load_demo_path(n)`, `parse_dat`/`read_dat`, `format_dat`/`write_dat`, `random_polygon_from_data`, `random_polygon`, `PANELS`, `EPSILON = 5e-5`, `PERTURB = 0.005`, `java_double_str(x)` (Java `Double.toString` formatting for golden-text comparison). | `Polygon.read/getPolygon`, properties |
| `polygon/geometry.py` | `relative_ccw`, `lines_intersect`, `pt_seg_dist(_sq)`, `distance`, `epsilon_equal`, `get_slope`, `get_y_intercept`, `get_point_on_line_from_x/y`, `get_intersect`, `get_segment_intersect`, `on_extension`, `on_reverse_extension`, `point_on_segment`, `purturb_point_along_seg(p, seg, ..., compat=)`, `GeneralPath` (float32 non-zero crossing rule), `seg_in_polygon` (plus the vectorised `seg_in_polygon_many`), `Polygon(vertices, name=None, *, validate=True)` (`edge`, `contains`, `is_reflex(i)`, `signed_area`; refuses clockwise, non-simple or non-finite input unless `validate=False`), `Path` (`length`, `pt_dist_from_start`, `point_on_path`, `intersect_point`, `intersect_points`, `all_critical_points`), `GeometryError`. Scalar code that copies the Java operation order; numpy only where bit-exactness is not affected. | `Algorithm` (primitives), `drawable.Polygon/Path` |
| `polygon/visibility.py` | `PhysicalGap` (with `full_start_edge/full_end_edge`), `physical_gaps(poly, q, eps, compat="fixed")`, `repair_gaps`, `point_status`, `visibility_polygon(poly, q, compat="fixed")`. | `getPhysicalGaps`, `getVisibilityPolygon` |
| `polygon/cuts.py` | `CutType`, `InflectionType`, `Cut`, `GeneralInflection`, `Bitangent`, `CriticalPoint`, `inflections(poly, kind)`, `general_inflection_cuts`, `bitangent_cuts`, `bitangent_lines`, `single_tangent_cuts`, `all_cuts` (Java order), `cut_intersect_points(path, cuts)`, `all_critical_points(path, cuts)`. | `getInflection`, `getGeneralInflectionCut`, `getBitangentCut`, `getBitangent`, `getSingletangentCut`, `getCuts`, `getAllCriticalPoints` |
| `polygon/gaps.py` | `Gap` (id, iedge, relative_time, state, to_gaps as an **ordered tuple in Java HashSet order**), `java_hashset_order(ids)`, `vertex_in_gap`, `line_in_gap`, `same_physical_gap(..., compat=)`, `collapse_physical_gaps`, `sample_path(poly, path, compat=)` → `PathSamples`, `get_gaps(poly, path, compat="fixed", *, java_matching=True)` → `GapHistory`, `format_gap_sets(history, with_time=True)`, `GapTrackingError`. | `getGaps` and helpers, `Gap.toString` |
| `polygon/events.py` | `gap_history_to_events(history)` → `GapEvents` (initial labels, `Appear/Disappear/Split/Merge` steps, times). This uses the existing `shadowinfo.events` types, so the polygon world feeds every existing filter (combinatorial, LP, probabilistic). It also has `never_see_evader(history, compat=)` (Java NSE, which sets `NSEState`) and `never_see_evader_via_filter(history)` (a cross-check through `counting_filter`). | `pe.Algorithm` |
| `polygon/oracle.py` | `JavaRandom` (an exact clone of `java.util.Random`: 48-bit LCG and `nextDouble`; `Math.random()` is one of these), `SingleTypeAgentOracle(n_targets, rng=None, merge_bug=None, compat="fixed")` (`merge_bug` defaults to `compat == "java"`) with `distribute_one_batch` and events matching `SingleTypeAgentEvent` (`__str__` in Java format), `OracleError`. `oracle.sequence()` converts to a `ShadowSequence` of `Appear/Disappear/Split/Merge` with the hidden total as `total`. | `SingleTypeAgentOracle`, `SingleTypeAgentEvent`, `Oracle`, `Event` |
| `polygon/stagent.py` | `Equation` (Java `toString`), `derive_gap_evolving_equations(history, oracle)`, `derive_shadow_info_state(history, oracle)` → `JavaBipartiteGraph` (Java-shaped BG with pooled vertex 0 and self-loops; `__str__` with sorted adjacency), `java_bounds(bg, order="sorted")` (the §2.10 recipe on `maxflow.FlowNetwork`; prints `Max/Min flow in shadow …`), and `exact_bounds(source, hidden_total=None)` (`source` is an oracle or a `ShadowSequence`) → `CombinatorialFilter(..., total=...)` bounds. | `SingleTypeAgentAlgorithm`, `Equation`, `BipartiteGraph`, `MaxFlow` |
| `polygon/simulate.py` | New, with no Java counterpart: `PolygonSimulator`, `simulate`, `validate_path`, `shadow_pockets` (§6). | – |
| `polygon/__main__.py` | `python -m shadowinfo.polygon run P12` reproduces a panel's stdout (§2.11) with a seed. It also has `list`, `cuts N` and `vis N x y`. | `ProjectPanel*` constructors |

Reuse from existing code:
- `events.Appear/Disappear/Split/Merge` and `ShadowSequence`.
- `bipartite.BipartiteIState`, which is semantically the same as `deriveShadowInfoState`, without the pooling quirk.
- `nondeterministic.CombinatorialFilter(total=...)` for exact bounds; `bounds_paper` for Eqs. 7/8.
- `counting_filter` for NSE cross-checks.
- `maxflow.FlowNetwork` (Edmonds–Karp) as the engine for `java_bounds`.
- `lp.ShadowLP` / `lp_bounds` as a third opinion.
- `probabilistic.ExactFilter`, which can now run on polygon-generated sequences, something the Java never did.

`grid/` stays as it is. The polygon pipeline is a second front end that produces the same event
JSON (§2 of `DESIGN.md`).

Performance: the Java `getPhysicalGaps` is O(n²)–O(n³) per sample. Java takes about 3.5 s
for all panel runs, mostly P14b: n = 267 and 671 samples. The port precomputes edge arrays,
caches `seg_in_polygon(q→v_j)` per sample (Java recomputes it for every `i`), vectorises the
exclusion tests with numpy and keeps the scalar code for the final accept/reject, so results stay
bit-identical. `sample_path` results are cached on the polygon (the last 32 paths per mode). The
full P14b pipeline (sampling, gap tracking, oracle, exact bounds) takes about 2.6 s. The slowest
tests carry the `slow` marker and still run by default.

### 5.2 Golden outputs to generate from the Java (parity fixtures)

**Status:** this is done. The harness is `tools/java_reference/`, with `GoldenDump.java`,
`run.sh`, `convert.py` and a README. The fixtures are in `tests/fixtures/java/` and are
checked by `tests/test_java_reference.py`. The README records the results. The most
important ones:

- Java's `deriveShadowBounds` matches the exact bounds in only about 74% of identity-hash
  orders on the real runs, not just on crafted cases.
- All 5 seeds of the original oracle on the Fig. 15(b) path give infeasible observations,
  because of B7.
- `getGaps` also crashes on the digitised ICRA'08 Fig. 8(b) path, which lies entirely inside
  the polygon.

The plan as originally written follows.

Generate these with a headless harness (it landed as `tools/java_reference/src/GoldenDump.java`,
which also does the identity-hash burn-in probes of item 8). The harness must call `FileHelper.setBaseUrl("file:<orig>/source/")`
and run with `-Djava.awt.headless=true`.
Write the outputs to `tests/fixtures/java/*.json`, using `Double.toString` text for every double
(it round-trips exactly in Python via `float()`).

1. **Cuts, per polygon 1–14:** the full lists from `getSingletangentCut`,
   `getInflection(NONGENERAL)`, `getInflection(GENERAL)`, `getGeneralInflectionCut` (with
   `fromLine` and `ccw`), `getBitangentCut` (with `this`, `opp`, `curveTo` and both segments),
   `getBitangent`, and `getCuts` (with type and order). The counts are in §2.5. Compare exactly,
   or within 1e-9 if a deliberate non-bit-exact path is chosen.
2. **Visibility, per polygon:** `getPhysicalGaps(p, q)` and `getVisibilityPolygon(p, q)` for
   (a) every waypoint of every §3 path that lies inside the polygon, (b) every perturbed
   critical point that `getGaps` samples on the runs in item 4, and (c) a 25×25 lattice of
   points with `pointInPolygon`. Each case records `(startEdge, endEdge, startPoint, endPoint)`
   and the vertex list.
3. **`pointInPolygon` and `segInPolygon`** on the same lattice, plus segments from lattice
   points to every vertex (this checks the float32 path rule and the grazing rule).
4. **Critical points** `path.getAllCriticalPoints(getCuts(p))` as `(distance, type, point)`
   for P12, P1, P5, P14a, P14b, P1s, P13, P13+, P5s and P1L. The counts are in §3.
5. **`getGaps`**, as the printed set lines (`relativeTime` + `Gap.toString`) for P12, P1, P5,
   P14b, P1s, P13, P5s and P1L. Also the **expected exceptions**: P14a gives an NPE in
   `collapsePhysicalGaps`, and P13+ gives `ArrayIndexOutOfBounds`. In compat mode Python must
   raise `GapTrackingError` at the same critical-point index.
6. **NSE**, from `processGapHistoryInformation(gapss, NEVER_SEE_EVADER)`: the printed sets for
   all runs in item 5. P1 and P5 are what `ProjectPanel2/3` print.
7. **Oracle pipeline.** `Math.random` is unseeded, so the harness seeds it: it gets
   `Math.random`'s static `java.util.Random` by reflection, calls `setSeed(seed)` on it and runs
   the original, unmodified `SingleTypeAgentOracle` (see `tools/java_reference/README.md`; the plan first proposed a `SeededSingleTypeAgentOracle` copy). Record,
   for P12 (`N=100`), P13 (`N=100`), P14b (`N=1 000 000`), P1 and P5 with seeds 1–5: the event
   strings, the equations, and the BG as sets (adjacency order is identity-hash dependent, so
   compare sorted). Python's `JavaRandom` + oracle must reproduce the events and equations
   exactly, including B7.
8. **Bounds.** Record the Java `Max/Min flow` lines for the runs in item 7, together with
   `shadowinfo` exact bounds for the same events. The test asserts:
   - `exact_bounds` contains the true counts the oracle used, whenever the run has no merge (B7);
   - `java_bounds(order="sorted")` equals the Java lines **only where the Java result does not
     depend on order** (verified by re-running the Java with several identity-hash burn-ins);
   - the B10 counterexample is checked against the three Java outcomes and the exact `[0,0]` / `[1,1]`.
9. **Paper cross-checks**:
   - P12 = ICRA'08 Fig. 4. Its event list maps onto `tests/fixtures/tor_fig11.json` by
     swapping 6↔7 and treating 5 as initial.
   - P13 gives 85 events and 18 bounded final shadows (T-RO Fig. 15a).
   - P14b gives 385 events, 491 IDs and 12 bounded final shadows (T-RO Fig. 15b).

### 5.3 Order of work

**Status:** steps 1–7 are done (§6).

1. `io` + `geometry`, with fixture tests 1 and 3.
2. `cuts`, with tests 1 and 4.
3. `visibility`, with test 2.
4. `gaps` + `events`, with tests 5, 6 and 9.
5. `oracle` + `stagent` (with `JavaRandom`), with tests 7 and 8.
6. `__main__` and a README section: "Polygon front end (port of the original Java)".
7. Optional: a JS port for the demo, reusing the same fixtures, as for `grid`.

---

## 6. Port status

Every function of §1 is ported to `shadowinfo.polygon` (module table in §5.1; the "shadowinfo"
column of the §2 tables names each counterpart), with the CLI `python -m shadowinfo.polygon`
standing in for the panels' stdout (§2.11). Parity with the Java is tested against the §5.2
fixtures: `tests/test_java_reference.py`, `tests/test_polygon_geometry.py`, `test_polygon_gaps.py`
and `test_polygon_stagent.py`.

- `compat="java"` reproduces cuts, critical points, physical gaps, visibility polygons, gap
  histories, NSE, oracle runs (seeded `java.util.Random`, both the original and the merge-fixed
  oracle), equations and the bipartite graph exactly, and raises `GapTrackingError` at the same
  critical point where Java throws (P14a, P13+, P13+tail, `fig_ICRA08-Fig8b`).
  `java_bounds(order="sorted")` gives Java's printed bounds wherever they do not depend on the
  identity-hash order, and one of Java's outcomes on every golden run.
- `compat="fixed"` (default) applies the [fix] items of §4 and never throws. Gap tracking
  (`get_gaps`, default `java_matching=True`) keeps Java's step whenever it is consistent and reads
  the step off the hidden boundary chains only where Java throws or its step is inconsistent. On every
  run where Java does not throw, the history is Java's, set for set and link for link: all golden runs,
  and a differential against the live Java on about 1260 such runs (random valid paths in the 14
  maps, perturbed paper paths, synthetic polygons). The one difference is B17: on a path that crosses
  or retraces itself the fixed set times increase, where Java's repeat or go backwards (one retracing
  run, times only).
- `get_gaps(..., java_matching=False)` (opt-in) also re-derives from the hidden chains every step
  where Java's consistent step disagrees with them: an ambiguous `samePhysicalGap` match (B3) that
  shifts every ID one gap along the boundary. IDs then stay on their physical shadows, which the
  simulator needs for its ground truth, but the history is no longer Java's: 5 of 174 events differ
  on P5, 5 of 47 on P5s and 29 of 385 on P14b (critical points 124–130 and 219–223); on random valid
  paths about 18% of runs differ, about 12% with different final shadow IDs. Event types, set times
  and the number of IDs are unchanged.
- Fixed-mode geometry repairs what Java gets wrong at degenerate points. `physical_gaps` and
  `visibility_polygon` replace a phantom gap by its complement, drop duplicate gaps and, within 1e-3
  of a line through two vertices, check the scan against exact orientation signs; elsewhere they are
  Java's bit for bit. Both raise `GeometryError` for a point outside the polygon. `sample_path` pulls
  a sample that Java's perturbation puts outside the polygon back inside along the same direction
  (recorded in `PathSamples.moved`). `python -m shadowinfo.polygon run` in fixed mode warns on stderr when the path fails
  `simulate.validate_path` (`P13+` and `P13+tail` leave map 13). `Polygon` refuses clockwise, non-simple or non-finite vertices
  (`validate=False` skips the check). NSE with an empty first gap set labels later gaps CLEAR (Java
  returns early). The oracle raises `OracleError` where Java throws.
- New, with no Java counterpart: `simulate.PolygonSimulator` (random-walk targets, shadow pockets,
  FOV events, ground truth; it tracks with `java_matching=False`), `validate_path` (refuses paths
  that leave the polygon or run along an edge) and the browser demo on all 14 maps (`docs/js/polygon*.js`).
- Paper runs (`examples/`): `fig11_maxflow.py` re-derives T-RO Fig. 11(b)–(c) and the bounds
  [10, 24] from map 12 and the `ProjectPanel` path; `paper_fig15.py` compares the Fig. 15 runs with
  every number the paper reports. Fig. 15(b)'s 124 bipartite vertices are the per-shadow count
  (61 + 63, `BipartiteIState`) while its 339 edges are Java's pooled graph (§1);
  `original_maps.py` draws the 14 maps with all paths (`docs/assets/original_maps.png`).

### 6.1 Every non-UI public Java function, its Python counterpart and its test

This table covers every `public` method of every non-UI class in `source/geometry` (listed from the source,
not from §2), plus the private helpers that carry the algorithm. Each "test" compares with output of the
real Java classes (golden fixtures of `tools/java_reference/`, or values printed by the original classes
on JDK 17 and quoted in the test), unless it says otherwise. Test files: **G** `tests/test_polygon_geometry.py`, **Gp**
`tests/test_polygon_gaps.py`, **S** `tests/test_polygon_stagent.py`, **J** `tests/test_java_reference.py`.
Python names are relative to `shadowinfo.polygon`.

| Java (class.method) | Python counterpart | test |
|---|---|---|
| **`geometry.Algorithm`** | | |
| `setEpsilon`, `setPurturb` (static fields `epsilon`, `purturb`) | `io.EPSILON`, `io.PERTURB`; `eps=` / `purturb=` keyword of every function that uses them | defaults: every parity test; other values: not tested |
| `INFLECTION_TYPE`, `CUT_TYPE` | `cuts.InflectionType`, `cuts.CutType` | G `test_cuts_match_java` |
| `epsilonEqual`, `getSlope` | `geometry.epsilon_equal`, `get_slope` | G `test_small_primitives_match_java` |
| `getYIntercept(line)`, `getYIntercept(k, p)`, `getXIntercept` (bug: = y-intercept) | `geometry.get_y_intercept`, `get_x_intercept` | G `test_get_intersect_quirks` |
| `getPointOnLineFromX`, `getPointOnLineFromY` | `geometry.get_point_on_line_from_x`, `_from_y` | G `test_get_intersect_quirks` |
| `getIntersect`, `getSegmentIntersect` | `geometry.get_intersect`, `get_segment_intersect` | G `test_get_intersect_quirks`, `test_numpy_primitives_equal_scalar`; through every cut test |
| `onExtension`, `onReverseExtension`, `pointOnSegment` | `geometry.on_extension`, `on_reverse_extension`, `point_on_segment` | G `test_get_intersect_quirks`, `test_numpy_primitives_equal_scalar`, `test_small_primitives_match_java` |
| `segmentsOnSameLine` (bug; unused) | `geometry.segments_on_same_line` (bug kept) | G `test_small_primitives_match_java` |
| `purturbPointAlongSeg` | `geometry.purturb_point_along_seg(compat=)` | G `test_purturb_point_along_seg`, `test_critical_points_match_java` |
| `segInPolygon` | `geometry.seg_in_polygon` (scalar `seg_in_polygon_scalar`, vectorised `seg_in_polygon_many`) | G `test_visibility_matches_java`, `test_seg_in_polygon_to_vertices_scalar`, `test_seg_in_polygon_vectorised_equals_scalar` |
| `getPhysicalGaps` | `visibility.physical_gaps(compat=)` | G `test_visibility_matches_java`, `test_run_visibility_matches_java`, `test_sample_physical_gaps_match_java`, `test_scan_equals_literal_java_loop` |
| `getVisibilityPolygon` | `visibility.visibility_polygon(compat=)` | G `test_visibility_matches_java`, `test_run_visibility_matches_java` |
| `getCuts` | `cuts.all_cuts` | G `test_cuts_match_java`, J `test_polygon_dump` |
| `getSingletangentCut` | `cuts.single_tangent_cuts` | G `test_cuts_match_java` |
| `getInflection(poly, type)` | `cuts.inflections(poly, kind)` | G `test_cuts_match_java` |
| `getGeneralInflectionCut` | `cuts.general_inflection_cuts` | G `test_cuts_match_java` |
| `getBitangentCut`, `getBitangent` | `cuts.bitangent_cuts`, `bitangent_lines` | G `test_cuts_match_java` |
| `getGaps` | `gaps.get_gaps(compat=)` → `GapHistory`; sampling in `gaps.sample_path` | Gp `test_get_gaps_java_matches_golden`, `test_get_gaps_java_throws_where_java_throws`, `test_fixed_is_identical_to_java`, `test_fixed_default_matches_live_java_on_a_random_path`; J `test_get_gaps_failures` |
| private `collapsePhysicalGaps`, `samePhysicalGap`, `vertexInGap`, `lineInGap` | `gaps.collapse_physical_gaps`, `same_physical_gap(compat=)`, `vertex_in_gap`, `line_in_gap` | Gp `test_collapse_physical_gaps_crashes_like_java`, `test_gap_helpers` |
| **`geometry.cut.*`** | | |
| `Cut(line, type)`, `getCut`, `getCutType` | `cuts.Cut` (`line`, `type`) | G `test_cuts_match_java` |
| `GeneralInflection`: `getFromLine`, `isCounterClockWise` | `cuts.GeneralInflection` (`from_line`, `counter_clockwise`) | G `test_cuts_match_java` |
| `Bitangent`: `getThisPoint`, `getOppositePoint`, `getOppositeSegment`, `getCurveToPoint` | `cuts.Bitangent` (`this_point`, `opposite_point`, `opposite_segment`, `curve_to_point`) | G `test_cuts_match_java` |
| `PathCutIntersectPoint`: `getPoint`, `getDistance`, `getSeg`, `getCut`, `getCutType` | `cuts.CriticalPoint` (`point`, `distance`, `seg`, `cut`, `cut_type`, plus `cut_index`, `seg_index`) | G `test_critical_points_match_java` |
| **`geometry.drawable.Path`** (non-drawing part) | | |
| `Path(list)`, `addPoint`, `getPointArray`, `getLineArray` | `geometry.Path(points)` (`points`, `segments`) | G `test_critical_points_match_java` |
| `getLength`, `ptDistFromStart`, `pointOnPath` | `Path.length`, `pt_dist_from_start` (B17 kept; `pt_dist_on_segment` is the private overload), `point_on_path` | G `test_path_distance_quirks` |
| `getIntersectPoint(Line2D)` (unused) | `Path.intersect_point` | G `test_path_get_intersect_point_matches_java` |
| `getIntersectPoints(Line2D[])` | `Path.intersect_points` | G `test_critical_points_match_java` (the panels' `panel_intersect_points`), `test_critical_points_last_writer_wins` |
| `getIntersectPoints(Cut[])` | `cuts.cut_intersect_points` | G `test_cut_intersect_points_match_java` |
| `getAllCriticalPoints` | `cuts.all_critical_points` (= `Path.all_critical_points`) | G `test_critical_points_match_java`, `test_critical_points_last_writer_wins` |
| **`geometry.drawable.Polygon`** (non-drawing part) | | |
| `Polygon(dc)`, `addVertex`, `done`, `getNumberOfVertices`, `getPointArray`, `getGeneralPath` | `geometry.Polygon(vertices)` (`n`, `vertices`, `general_path`) | G `test_load_polygon_matches_java`, `test_small_primitives_match_java` |
| `getLineSegmentArray` | `Polygon.edges`, `Polygon.edge(i)` | G `test_small_primitives_match_java`, every cut test |
| `pointInPolygon` (float32 `GeneralPath.contains`) | `Polygon.contains`, `geometry.GeneralPath` | G `test_point_in_polygon_matches_java`, `test_general_path_contains_float32_rule` |
| `read`, `getPolygon(n, dc)` | `io.parse_dat` / `read_dat` (B16 fixed), `io.load_polygon(n)` | G `test_load_polygon_matches_java`, `test_parse_dat_fixes_blank_lines_and_spacing`, `test_dat_files_and_provenance` |
| `write` | `io.format_dat`, `io.write_dat` | G `test_dat_writer_and_random_polygons_match_java` |
| `randomPolygonFromData`, `getNumberOfPolygonDataFiles`, `getDataFilePostfix` | `io.random_polygon_from_data(rng)`, `io.NUMBER_OF_POLYGON_DATA_FILES`, `io.DATA_FILE_POSTFIX` | G `test_dat_writer_and_random_polygons_match_java` |
| `randomPolygon` (unused) | `io.random_polygon` | G `test_dat_writer_and_random_polygons_match_java` |
| **`geometry.gap.*`, `IDGenerator`** | | |
| `Gap(id, ie)`, `get/setId`, `get/setIEdge`, `get/setRelativeTime`, `get/setState`, `getToGapSet`, `addGap` | `gaps.Gap` (`id`, `iedge`, `relative_time`, `state`, `to_gaps` in Java `HashSet` order, `add_gap`) | Gp `test_get_gaps_java_matches_golden`, `test_link_order_in_fixtures_is_java_hashset_order`, `test_java_hashset_order` |
| `Gap.toString` | `gaps.format_gap`, `format_gap_set(s)`, `GapHistory.lines` | Gp `test_gap_to_string`; G `test_java_double_str_matches_gap_set_lines`; J `test_oracle_replica_and_text_formats` |
| `Gap.getGapMap` | `{g.id: g for g in gaps}` inside the oracle (last wins, as in Java) | through S `test_oracle_runs_equal_java` |
| `PhysicalGap(se, ee, sp, ep)`, getters, `getFullStartEdge`, `getFullEndEdge` | `visibility.PhysicalGap` (`start_edge`, `end_edge`, `start_point`, `end_point`, `full_start_edge`, `full_end_edge`) | G `test_physical_gap_full_edges`, `test_visibility_matches_java` |
| `GapState` (both interfaces) | `Gap.state` (any object) | – |
| `IDGenerator.getNextId` | `gaps.IDGenerator.next_id` | Gp `test_collapse_physical_gaps_crashes_like_java`; through `get_gaps` |
| **`geometry.pe`** | | |
| `pe.Algorithm.processGapHistoryInformation(gapss, NEVER_SEE_EVADER, data)`, `pursuerNeverSeeEvader`, `SCENARIO` | `events.never_see_evader(history, compat=)` (the only scenario); cross-check `never_see_evader_via_filter` | Gp `test_events_and_nse_match_java`, `test_never_see_evader_b6_and_states`; J `test_gap_history_events_bg_nse` |
| `NSEState.CLEAR`, `CONTAMINATED`, `isClear`, `toString` | `events.NSEState` | Gp `test_gap_to_string`, `test_never_see_evader_b6_and_states` |
| `Equation`: `addTerm`, `get/setCoeffMap`, `get/setConstant`, `toString`, `main` | `stagent.Equation` (`add_term`, `coeffs`, `constant`, `str`) | S `test_event_and_equation_strings` (includes the `main` output) |
| `Event.getEventRelativeTime`; `SingleTypeAgentEvent` getters/setters, `toString` | `oracle.SingleTypeAgentEvent` (`type`, `time`, `from_gaps`, `to_gaps`, `moved`, `visible`, `str`) | S `test_event_and_equation_strings`, `test_oracle_runs_equal_java` |
| `Oracle.EVENT_TYPE` | `oracle.EventType` | S `test_event_and_equation_strings` |
| `SingleTypeAgentOracle(N)`, `initialize` (private `distributeAgents`) | `oracle.SingleTypeAgentOracle(n, rng, merge_bug=, compat=)`, `.initialize(history)` | S `test_oracle_runs_equal_java` (all 76 golden oracle runs), `test_merge_bug_in_java_matches_the_original_line`; J `test_oracle_replica_and_text_formats` |
| `getTotalAgentNumber`, `getNumberOfEvents`, `getEvent`, `getEventByTime` | `get_total_agent_number`, `number_of_events`, `get_event`, `event_by_time` | S `test_oracle_accessors` |
| private `distributeOneBatch` | `SingleTypeAgentOracle.distribute_one_batch` | S `test_distribute_one_batch` |
| private `addRandomVisibilityEvents` (never called, B13) | not ported; FOV events come from `simulate.PolygonSimulator` (`Enter` / `Exit`) | – |
| `SingleTypeAgentState(n)`, `get/setNumberOfAgents` | `SingleTypeAgentOracle.counts` / `final_counts` / `truth()` | S `test_oracle_runs_equal_java` (`truth_final`), `test_ground_truth_and_equations_hold` |
| `Math.random()` (`java.util.Random`) | `oracle.JavaRandom` | S `test_java_random_equals_jdk17`, `test_java_random_index_check_and_reference_clone`; J `test_index_and_random_clone` |
| `SingleTypeAgentAlgorithm.deriveGapEvolvingEquations` | `stagent.derive_gap_evolving_equations` | S `test_oracle_runs_equal_java` |
| `SingleTypeAgentAlgorithm.deriveShadowInfoState` | `stagent.derive_shadow_info_state` → `JavaBipartiteGraph` | S `test_oracle_runs_equal_java` (structure and printed text); J `test_gap_history_events_bg_nse` (= `BipartiteIState`) |
| `SingleTypeAgentAlgorithm.deriveShadowBounds` | `stagent.java_bounds(bg, order)` (comparison only); default `stagent.exact_bounds` | S `test_java_bounds_and_exact_bounds`, `test_random_hash_orders_reproduce_every_java_outcome_on_p12`, `test_b10_counterexample`, `test_infeasible_runs_are_exactly_the_merge_bug_fig15b_runs`; J `test_java_bounds_are_order_dependent_and_b7_breaks_consistency` |
| **`geometry.graph.*`** | | |
| `BipartiteGraph` (`leftVertexMap`, `rightVertexMap`, `rightVertexEdgeMap`, `toString`) | `stagent.JavaBipartiteGraph` (`left_vertices`, `right_vertices`, `edges`, `to_dict`, `str`) | S `test_oracle_runs_equal_java` |
| `Vertex(id)`, `Edge(from, to[, w, c])`, `createResidueEdge` | `stagent.BGVertex`, `BGEdge`; flow copies `_FVertex`, `_FEdge` (private) | through the `java_bounds` tests |
| `MaxFlow.getMaxFlow` | `stagent._get_max_flow` (literal Edmonds–Karp port with `HashMap` order as a parameter, the engine of `java_bounds`) | through the `java_bounds` tests |
| `MaxFlow.getIncrementalMaxFlow` (unused), `MaxFlow.main` (prints nothing) | not exposed: they are `_search_all_augment_paths` + `_construct_flow_graph` without the rebuild | – |
| **panel output** | | |
| `ProjectPanel`…`ProjectPanel5` constructors' stdout (§2.11) | `python -m shadowinfo.polygon run …` (`stagent.project_panel_lines`) | S `test_cli_project_panel_java_equals_golden_output`, `test_cli_project_panel5_java_equals_golden_output`, `test_cli_nse_panel_and_drawing_panel`, `test_cli_entry_point` |

UI and infrastructure, deliberately not ported: `AbstractDrawable`, `Drawable`, `DrawingContext`, `LineSegment`, `Point`,
`PointAgent`, every `draw*`/`paint`/colour accessor, the `pe/ui` applets and panels, `FileHelper` (applet/base-URL
resource loading, replaced by package data in `io`), `Geometry` (reads `geometry.properties` for `DEFAULT_DC`) and
`Triple` (a 3-tuple used by the BFS).

How this was verified (2026-10): the full pytest and `npm test` suites; the golden fixtures regenerated
from the original classes (`tools/java_reference/run.sh`) are byte-identical to the committed ones; the
functions the panels never call (`write`, `randomPolygon*`, `getIntersectPoint(s)`, `segmentsOnSameLine`,
the oracle accessors) were checked against values printed by the original classes on JDK 17.
