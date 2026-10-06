# Design: `shadowinfo` (Python) and the JS demo

Reference: J. Yu and S. M. LaValle, *Shadow Information Spaces: Combinatorial
Filters for Tracking Targets*, IEEE T-RO 28(2), 2012 (section numbers below
refer to it), plus ICRA 2008 / ICRA 2010 conference versions.

## 1. Repository layout

```
shadowinfo/                Python package (pure Python + numpy/scipy)
  __init__.py              re-exports the public API below
  events.py                event types, ShadowSequence (validated event list)
  bipartite.py             incremental bipartite I-state (Fig. 9)
  maxflow.py               Edmonds–Karp + flow with lower bounds (no deps)
  nondeterministic.py      combinatorial filter: bounds, refine, counting, PE, teams
  fov.py                   FOV-event batching (Sec. V-D) and FOV→component conversion
  lp.py                    LP baseline (Sec. IV) via scipy.optimize.linprog
                           (integral by total unimodularity; asserted)
  probabilistic.py         exact joint-distribution filter (Sec. VI-B..D),
                           truncation heuristics TR / RT / RT-LA (VI-E), Monte Carlo
  rng.py                   Mulberry32 PRNG shared bit-for-bit with JS
  testing.py               random consistent instances + witness replay (tests, fixtures)
  grid/
    __init__.py
    gridmap.py             ASCII map loading, cell visibility (shared algorithm with JS)
    shadows.py             shadow components + frame-to-frame component events
    simulate.py            robot path + random-walk targets → events + ground truth
    maps/office.txt        the fixed map (ASCII) shared by Python and JS (package data)
    maps/office_path.json  robot waypoints for the fixed map
  polygon/                 port of the original Java (no applets) + polygon simulator (§8)
    io.py                  maps, paths, Java constants and Double.toString
    geometry.py            Algorithm primitives, java.awt.geom semantics, Polygon, Path
    cuts.py                inflection / bitangent / single-tangent cuts, critical points
    visibility.py          physical gaps, visibility polygon
    gaps.py                getGaps: gap tracking along a path (compat java | fixed)
    events.py              gap history -> shadowinfo.events, never-see-evader
    oracle.py              java.util.Random clone, SingleTypeAgentOracle
    stagent.py             equations, Java bipartite graph, java_bounds, exact_bounds
    simulate.py            PolygonSimulator: random-walk targets, pockets, FOV events
    __main__.py            CLI: python -m shadowinfo.polygon list | run | cuts | vis
    data/                  {1..14}.dat (byte-identical originals), paths.json
examples/                  runnable scripts reproducing paper examples + figures
tests/                     pytest; tests/fixtures/ holds paper and parity fixtures,
                           tests/fixtures/java/ the golden outputs of the original Java
docs/index.html            JS demo (static, no build), served by GitHub Pages
docs/js/*.js               classic scripts on a global `SI` namespace (see §6); polygon.js,
                           polygon_sim.js and maps_polygon.js port shadowinfo/polygon
                           (gap tracking, simulator, the 14 maps and all paths)
docs/css/style.css         demo styles (light/dark)
docs/assets/               README figures (examples/figures.py) + demo screenshots
docs/notes/                paper extraction notes, conference-vs-journal deltas
js-tests/                  node --test unit tests + Python/JS parity tests
tools/                     generators for JS map data and parity fixtures
tools/java_reference/      headless harness that runs the original Java classes and writes
                           tests/fixtures/java/ (GoldenDump.java, run.sh, convert.py, README)
package.json               `npm test` = `node --test js-tests/*.test.js`
```

## 2. Event model (`events.py`)

Shadows are identified by integer labels, never reused. Events (dataclasses,
frozen):

| event | fields | meaning |
|---|---|---|
| `Appear(s, lo=0, hi=0)` | new shadow `s` holding `lo..hi` targets | component appear (Eq. 3, `d_i`) |
| `Disappear(s, lo, hi)` | shadow `s` revealed holding `lo..hi` (normally `lo==hi`) | component disappear |
| `Split(s, a, b)` | `s` → `a`, `b` | evolve/split |
| `Merge(a, b, s)` | `a`, `b` → `s` | evolve/merge |
| `Enter(s, k=1)` | `k` targets enter shadow `s` from FOV | FOV enter |
| `Exit(s, k=1)` | `k` targets exit shadow `s` into FOV | FOV exit |

`ShadowSequence(initial: dict[label, (lo, hi)], events: list)` validates
(labels alive / fresh, no reuse) and exposes `alive_at_end()`, `all_labels()`,
`to_json()/from_json()` (JSON format shared with JS:
`{"initial": {"1": [0, 3], ...}, "events": [{"type": "split", "s": 1, "a": 2, "b": 3}, ...]}`;
`hi` may be `null` meaning +∞).

Invalid sequences raise `InvalidSequenceError` (e.g. a dead or reused label, or
`Split(s, a, a)`).

## 3. Nondeterministic filter

### 3.1 Bipartite I-state (`bipartite.py`, Sec. V-B, Fig. 9)
Left vertices = initial shadows + appeared shadows (+ FOV-enter pseudo-shadows),
each with `[lo, hi]`. Right vertices = disappeared shadows (with `[lo, hi]`) +
shadows alive at the end. Maintain `reach[s] = set(left vertices)` for each
alive shadow: appear → `{L_s}`; split → both children copy; merge → union;
disappear → right vertex with edges from `reach[s]`. Incremental: process events
one at a time (`BipartiteIState.apply(event)`).

FOV events (Sec. V-D): `Enter(s,k)` ≡ `Appear(t, k, k)` + `Merge(s, t, s')`;
`Exit(s,k)` ≡ `Split(s, s', t)` + `Disappear(t, k, k)`. Implement both the naive
per-event conversion and the batched form (`d_min`, `d_tot` counter per shadow,
at most one batch exit then one batch enter, flushed whenever the shadow takes
part in a component event). Both must give identical bounds.

### 3.2 Bounds (`nondeterministic.py`)
Primary algorithm: **feasible flow with lower bounds** on the network
`S → L_i [lo_i, hi_i]`, `L_i → R_j [0, ∞)` for edges of the bipartite graph,
`R_j → T [lo_j, hi_j]` for disappeared shadows and `[0, ∞)` for final shadows,
optional total-count constraint `[N_lo, N_hi]` on the number of targets hidden
in the **initial** shadows (an edge from a super source feeding the initial
shadows' source edges). It deliberately does not count FOV enters: a target that
exits and re-enters would otherwise be counted twice, making e.g. a seen-once
single evader infeasible. Python and JS use the same semantics. Upper bound of a query set
`Q` of final shadows = max total flow into `Q` subject to feasibility; lower
bound = min. Implemented via the standard reduction to max-flow (super
source/sink) on our own Edmonds–Karp (`maxflow.py`), ∞ represented by `math.inf`.
Queries raise `InfeasibleError` if the observations are inconsistent, and
`ValueError` for an empty query or a label that is not alive (for
`left_bounds`: not an initial or appeared shadow). Group queries warm-start
from the cached feasible circulation.

Also implement the paper's literal recipe (Eqs. 7, 8) as `bounds_paper(...)`.
Finding: Eq. (7) (upper bound) matched the exact optimum on every random case,
but the literal Eq. (8) (lower bound, with `c(i,T) = l_i` for disappeared
shadows) can overstate the lower bound when a disappeared shadow has `l_i < u_i`;
using `c(i,T) = u_i` in the lower-bound run (`bounds_paper(..., corrected=True)`)
matches the LP. Both give [10, 24] on Fig. 11. Implement
the refinement of initial bounds (Eqs. 9, 10) as `refine_initial_bounds(...)`.
Every result must agree with the LP baseline in `lp.py` (TU ⇒ LP optimum is
integral) — tests check this on random sequences.

Public API:
```python
filt = CombinatorialFilter(initial={1: (2, 4), 2: (0, 3)}, total=None)   # total=(lo,hi) optional
filt.apply(event)                    # incremental
filt.bounds(shadows=[19])            # -> (lo, hi) for the set (sum), hi may be math.inf
filt.all_bounds()                    # -> {label: (lo, hi)} for each alive shadow
filt.refine_initial_bounds()         # -> {label: (lo, hi)} tightened initial bounds
filt.left_bounds([1, 2])             # -> (lo, hi) for a set of initial/appeared shadows
filt.initial_total_bounds()          # -> (lo, hi) on the number of targets at t0
```
Teams/attributes (Sec. V-F): `MultiTeamFilter` runs one `CombinatorialFilter`
per attribute. For "red-or-blue" initial conditions, `RedOrBlueFilter.bounds(q,
teams=None)` solves one joint network (one copy of the augmented graph per
team, shared `S`/`T`, a node `M_s` per mixed shadow feeding both teams' `L_s`)
and is exact over all red/blue assignments. The paper's two-pass procedure is
kept as `two_pass()` (returns `TwoPassBounds`) for comparison only: it is not
sound when per-team observations show that a mixed group held both colours. Counting
(`counting_filter`, initial `[0, ∞)`) and pursuit–evasion
(`pursuit_evasion_filter`, `[0,1]` per shadow with `total=(1,1)`;
`evader_status` classifies `hi == 0` clear, `lo >= 1` evader, else
contaminated) are thin helpers.

## 4. Probabilistic filter (`probabilistic.py`, Sec. VI)

State: ordered label tuple + `dict[tuple[int,...], float]` (joint pmf).
- appear: product with independent `P(s_k)`;
- disappear: Bayes update with observed `P(s_k)` then marginalise `s_k`, renormalise;
- merge: sum `x_i + x_j`;
- split: split rule `rule(n) -> {(n_a, n_b): p}`; default binomial with `p_a`
  (paper example: 0.5; demo: area-proportional);
- FOV observation `y` on shadow `s`: branch into enter/exit/null weighted by
  `P(e|y)`. When `x_s = 0` the exit branch is impossible and, following
  Algorithm 2, the remaining branches are rescaled **within that entry** (the
  entry keeps its mass) — global renormalisation gives 1/12, 1/6, 3/4 instead of
  Table III's 1/13, 2/13, 10/13. Monte Carlo correspondingly samples each FOV
  event among the events feasible at the current count and rejects trials only
  at disappear events.
- Observations: `FovObservation(s, y)` with `y ∈ {"enter","exit","null"}` and
  model `obs_model[y][e] = P(e | y)`; component events with a pmf are
  `ProbAppear`/`ProbDisappear` (`dist`) and `ProbSplit` (own binomial `p`).
  A plain `Appear(s, lo, hi)` is uniform over `lo..hi`; a plain
  `Disappear(s, lo, hi)` is the indicator likelihood `1[lo <= x <= hi]`;
  noise-free `Enter`/`Exit` add/remove `k` (ruling out `x_s < k`).
- Observation-model rows and custom split rules are rescaled to sum to one
  (negative weights and unknown keys raise `ValueError`). An observation that
  leaves no mass raises `InconsistentObservationError`. Split, merge or appear
  into an alive label raises `KeyError`.
- JSON extensions for probabilistic sequences: `{"type": "fov", "s", "y"}`,
  `"dist": {count: prob}` on appear/disappear, `"p"` on split
  (`ProbSequence.to_dict/from_dict`; masses written as floats, repeated joint
  keys summed).
- `ExactFilter`, and `TruncatedFilter` with `max_entries` and
  `mode in {"TR","RT","RT-LA"}` (RT: keep by `p * U(0,1)`; RT-LA: no truncation
  if a disappear is within the next 4 events or a merge within the next 2;
  whether FOV observations count toward the window is a parameter,
  `lookahead_counts_fov`, default true). If no mass is left after a
  truncation it raises `TruncationFailure` (the input may also be inconsistent;
  the exact filter tells the two apart).
- `monte_carlo(sequence, trials=1000, seed=...)` rejection-sampling baseline;
  its labels are sorted (`JointPmf.reordered(labels)` aligns it with
  `ExactFilter`, which keeps event order).
- `expected_counts()`, `marginal(label)`.

Regression targets from the paper: Table III / Fig. 12 example
(final `P(s4=0) ≈ 0.0769`, Monte Carlo ≈ 0.079 / 0.154 / 0.767), Fig. 11 max-flow
example (upper bound of `s19` = 24, lower bound = 10).

## 5. Grid world (`grid/`, mirrored in `docs/js/grid.js`)

**Bit-for-bit parity with JS is a requirement** — keep the algorithms simple
and integer-based where possible.

- Map: ASCII, `#` obstacle, `.` free. Cell `(r, c)`, row-major.
- Robot: sits at a free cell center. Visibility: cell `(r,c)` is visible iff the
  segment between the robot cell center and `(r,c)` center does not pass
  through the interior of any obstacle cell — an integer supercover traversal
  between cell centres; when the segment passes exactly through a lattice
  corner it is blocked only if *both* side cells are obstacles (the side cells
  are otherwise not tested) and the traversal steps diagonally. The predicate
  is symmetric. Sampling cell centres means that a visible wedge narrower than
  one cell (an oblique view through a doorway) yields isolated visible cells;
  they are treated as ordinary visible cells.
  Optional sensing radius `R` (Euclidean, in cells): a cell is in range iff its
  centre is within `R + 1/2` of the robot's, i.e. the integer compare
  `dr² + dc² <= R² + R` (a round disc; `<= R²` leaves 1-cell nubs on the axes
  that create spurious 1-cell shadows along walls).
- Shadow region = free ∧ ¬visible. Components: 4-connectivity, labelled by
  scanning row-major.
- Tracking: overlap bipartite graph between consecutive frames' components;
  each connected group with `m` old and `n` new components gives: (0,1) appear,
  (1,0) disappear, (1,1) live on (label kept), (1,2) split, (2,1) merge, general
  (m,n) → chain of merges then chain of splits. New labels are allocated by a
  global counter in a deterministic order (sorted by component scan order).
- Targets: `N` targets on free cells, random 4-neighbour walk (or stay) using
  `rng.py` (Mulberry32 — same integers in Python and JS).
- One tick = (1) robot advances one cell along its path; targets that were in
  shadow and are now visible emit `Exit(old_label)` *before* the component
  events; component events; targets that were visible and are now in shadow
  emit `Enter(new_label)` after them; (2) each target moves one step; shadow →
  visible emits `Exit`, visible → shadow emits `Enter`.
- A target revealed by a disappearing shadow is reported only through
  `Disappear(s, k, k)` (no `Exit`); a target covered by a newly appearing shadow
  only through `Appear(s, k, k)` (no `Enter`). FOV events are aggregated per
  shadow per substep (`k >= 1`).
- Output: event list (same JSON as §2), per-tick ground truth counts per
  shadow, frames for rendering.

## 6. JS demo (`docs/`)

Vanilla JS, no build step, no external deps. Files are **classic scripts** (not
ES modules, which Chrome blocks on `file://`) that attach to a global `SI`
namespace and also `module.exports = SI` under Node, so the page works when
opened by double-click, on GitHub Pages, and in `node --test`. Canvas rendering of the
fixed map: obstacles, visible region, shadows coloured by label with label ids,
robot + path, targets (ground truth toggle). Side panel: per-shadow true count
vs. filter bounds `[lo, hi]` (and expected count from the probabilistic filter
with truncation), event log, live bipartite-graph view (Fig. 11c style).
Controls: play/pause/step/reset, speed, #targets, seed, initial-knowledge mode
(exact counts / unknown count `[0,∞)` / single evader), sensing range (grid only).
Respects light/dark mode; usable at phone width.

**Map selector.** "Office (grid)" plus the 14 original maps ("Map N", with the
paper figure where there is one: 12 = T-RO Fig. 11(a) / ICRA'08 Fig. 4,
13 = T-RO Fig. 15(a), 14 = T-RO Fig. 15(b)). The page opens on map 12 with its
default path. A polygon map is drawn in the `.dat` frame, y up like the Java
window, never flipped. The page draws the polygon, the visibility polygon from
the sensor point, each shadow pocket (the hidden boundary chain closed by its
window, the dashed edge) tinted by label with a label pill (a split child inherits its parent's
colour; on every tick with events the 12 largest shadows are given distinct colours, smaller
ones keep theirs), the path, the robot
and the targets. The side panel is the same as for the grid: bounds from the
`CombinatorialFilter`, `E[n]` from the truncated probabilistic filter (split
parameter = pocket-area ratio), the event log and the bipartite view. The
polygon simulator (`SI.PolygonSimulator`) emits the same events as the grid,
so the filters need no changes.

**Paths.** Each map offers its demo path (`paths.json` "demo": a Java code path
driven back and forth, or a closed patrol tour) and every stored path on the map
(Java code paths, paths digitised from the paper figures) with its provenance.
Stored paths that leave the polygon (P13+, P13+tail) are listed but disabled.
**Draw your own path**: click waypoints. A waypoint must lie strictly inside the
polygon and each leg inside the closed polygon without running along a wall (the
`validatePath` predicates); a rejected click flashes the offending leg and explains
why. Double-click, Enter or Start runs the path from t = 0, back and forth. With 3 or
more waypoints, a click within 9 px of the start closes the path and starts it as a
loop.
URL parameters: `map`, `path` (a stored label, `default`, or `x1,y1;x2,y2;...`),
`n`, `seed`, `mode`, `radius`, `gt`, `speed`, `ticks`, `autoplay`, `theme`.

The gap track of a path is computed once per polygon and path (about 0.5-1 s for
P14b in Node 22). The reversed path's track is prefetched while the page is idle, so the
turnaround of an open path does not stall the animation.

The DOM-free `docs/js/controller.js` drives sim + filters (tested in Node, on
all 14 maps × 3 modes × 2 seeds in `js-tests/demo_smoke.test.js`). It also lists
the maps (`SI.demo.mapList`) and their paths (`SI.demo.pathOptions`), and checks
drawn waypoints (`SI.demo.checkWaypoint`).

The bipartite I-state grows with every FOV batch, so the controller
**re-roots** long runs: it replaces the history by a filter whose initial
shadows are the alive shadows with their current bounds plus the bound on
their sum. Each alive count `x_s` is the demand of terminal `R_s` (its net
inflow; every `R_s` drains into `T`), and the terminal demand vectors of a
network with lower/upper bounds form an integral base polyhedron (Hoffman);
dropping `T`'s coordinate gives an integral g-polymatroid, so the feasible
set is determined by its subset min/max. The re-rooted filter's set (box ∩
total) is one too, and singletons and the full set match by construction, so
a re-root is declared exact only after checking that bounds agree on every
subset `2 ≤ |A| < k` of the `k` alive shadows (`k ≤ maxCheckShadows`). The check
runs on a snapshot of the filter, spread over several ticks (at most `checkBudgetMs`
= 6 ms of subset checks per tick, so the animation does not stall); once it proves
the re-root exact, the filter re-rooted at the snapshot replays the events applied
since then, so the bounds are the same as without re-rooting. Otherwise a *relaxed*
re-root (an outer approximation — still sound) is forced only above a hard edge
limit (`hardEdges`), after a last exactness check capped at `forcedCheckMs` = 40 ms.
How many subsets a tick checks depends on the machine, so the tick at which a
re-root lands can differ between machines; the bounds cannot. In evader mode the total is `[h, h]` with `h` the number of
targets hidden at t0 (the evader may start in view). After a re-root the
bipartite view tags the root vertices `root` (header `root (t = …)`), and it
shows only the latest disappeared vertices, with a "+ N earlier disappeared not
shown" note.

## 7. Testing

- unit tests per module; paper regression examples;
- randomized: random valid shadow sequences + ground-truth target trajectories
  → true counts always within bounds; bounds equal LP optimum; bounds are tight
  (witness flows realised);
- batched FOV == naive FOV;
- probabilistic exact vs Monte Carlo (statistical tolerance) and vs brute-force
  per-target enumeration on tiny instances;
- parity: Python and Node produce identical visibility masks, labels and
  event JSON for the fixed map, path and seed.

## 8. Polygon front end (`polygon/`)

A port of the original Java implementation
([arc-l/shadow-information-space](https://github.com/arc-l/shadow-information-space), 2008-2012)
without its applets: what the `ProjectPanel*` classes compute, function by function
(`docs/notes/original_java.md` has the inventory, the quirks B1-B18 and the porting plan;
`docs/notes/original_maps.md` matches maps and paths to the paper figures). It is a second front end
next to `grid/`: both produce the §2 events, so every filter runs on either.

- **Data.** `data/{1..14}.dat` are byte-identical copies of the original maps: `x y` per line,
  `[0, 1000]^2`, **y-up, counter-clockwise**. Nothing is ever flipped (the papers draw maps 12 and 14
  mirrored). `data/paths.json` holds every path of the Java panels (`P12`, `P13`, `P14b`, ...), the
  paths digitised from the paper figures (`fig_*`) and one demo path per map, each with its
  provenance. `examples/original_maps.py` draws them (`docs/assets/original_maps.png`).
- **Pipeline** (Sec. II-A, IV-B, V-B, V-C): `all_cuts(poly)` (critical lines, Java order) ->
  `all_critical_points(path, cuts)` -> `get_gaps(poly, path)` (physical gaps sampled just after each
  crossing; general-inflection crossings give appear/disappear, bitangent crossings split/merge;
  a `GapHistory` of gap sets with IDs and links) -> `gap_history_to_events` (§2 events) ->
  `SingleTypeAgentOracle(N, rng).initialize(history)` (simulated counts) ->
  `derive_gap_evolving_equations`, `derive_shadow_info_state` (Java's bipartite graph) and
  `exact_bounds` (`CombinatorialFilter` with `total=(hidden, hidden)`).
- **Two modes.** `compat="java"` copies the Java operation by operation (float32 `GeneralPath`,
  slope-intercept intersections, last-writer-wins `TreeMap`s, `HashSet` iteration order,
  `java.util.Random`, the merge bug) and raises `GapTrackingError` / `OracleError` where Java throws,
  at the same critical point. Golden fixtures of the real Java classes (`tools/java_reference/`,
  `tests/fixtures/java/`) check it. `compat="fixed"` (default) fixes the bugs and never throws. Gap
  tracking keeps Java's step whenever it is consistent and reads the step off the hidden boundary
  chains only where Java throws or is inconsistent, so on every run where Java succeeds the history is
  Java's (all golden runs and about 1260 live-Java differential runs; only the set times of a
  self-crossing or retracing path differ, B17). `get_gaps(..., java_matching=False)` (opt-in, used by
  the simulator) also overrides Java's ambiguous `samePhysicalGap` matches (B3), so IDs stay on their
  physical shadows; that history differs from Java's on about 18% of random runs (P14b: 29 of 385
  events). Fixed-mode `physical_gaps` / `visibility_polygon` repair Java's scan at degenerate points
  (phantom and duplicate gaps, points on a line through two vertices) and raise `GeometryError` for a
  point outside the polygon; a path sample that Java's perturbation puts outside is pulled back in;
  `Polygon` refuses clockwise or non-simple input. The oracle conserves targets (B7) and cannot loop
  forever (B8); bounds are exact and cover appeared shadows (B10, B11). Java's max-flow
  recipe is kept as `java_bounds(bg, order)` for comparison only; its result depends on `HashMap`
  identity order, which `order` simulates.
- **Speed.** Java's `getPhysicalGaps` is O(n^2)-O(n^3) per sample. The port caches the vertex
  visibility per sample and vectorises the exclusion tests with numpy, keeping Java's scalar
  accept/reject for bit-exactness, and caches `sample_path` results on the polygon (treat them as
  read-only): the full P14b pipeline (n = 267, 671 critical points) takes about 2.6 s, of which the
  oracle, graph and bounds take milliseconds. Tests on maps 5, 13 and 14 carry the `slow` marker
  (registered in `pyproject.toml`, run by default; the whole suite takes about 2 min).
- **Simulator** (`simulate.py`, no Java counterpart). `PolygonSimulator` moves the robot a fixed arc
  length per tick along a path (stop, loop or reverse at the end, `go_to` / `follow` for new paths);
  its component events and labels are exactly those of `get_gaps(..., java_matching=False)` on the
  same path. `validate_path` refuses waypoints outside the polygon and legs that leave it or run
  along an edge. Shadow pockets
  are the hidden boundary chains closed by their windows; targets do a random walk with `rng.Rng`,
  and the tick semantics (Exit before, Enter after the component events) are those of §5. The
  predicates use only `+ - * /`, so `docs/js/polygon_sim.js` reproduces the frames bit for bit.
- **CLI.** `python -m shadowinfo.polygon run ProjectPanel5 --compat java --seed 1` prints what the
  panel prints (gap sets, bipartite graph, bounds, equations, events); `list`, `cuts N`, `vis N x y`.
- **Paper runs.** `P12` gives ICRA'08 Fig. 4 / T-RO Fig. 11 (`examples/fig11_maxflow.py`); `P13`
  and `P14b` give T-RO Fig. 15(a)/(b): 85 and 385 component events, 491 labels, 18 and 12 bounded
  final shadows (`examples/paper_fig15.py` compares every reported number).

