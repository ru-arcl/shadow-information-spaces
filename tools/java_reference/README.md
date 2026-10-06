# Java reference harness (golden outputs of the original implementation)

This folder runs our **original Java implementation** of *Shadow Information Spaces*
(T-RO 2012, ICRA 2008) without its applet UI. It writes golden JSON fixtures to
`tests/fixtures/java/`, so the Python port can be tested against the original. The original
sources are **not modified**. The harness only calls their public methods, plus two small,
documented additions (see "What the harness adds").

| file | purpose |
|---|---|
| `src/GoldenDump.java` | The harness. It loads polygons and paths, calls the original classes, and writes JSON. |
| `src/Json.java` | Minimal JSON writer, so no extra dependency is needed. |
| `src/MergeFixedOracle.java` | Copy of `SingleTypeAgentOracle` with the merge double-count bug fixed. Used only for the `merge_fixed` variant. |
| `paths.tsv` | Every robot path: the active and commented paths from `ProjectPanel*.java`, and the paths digitised from the paper figures (`docs/notes/original_maps.md`). |
| `run.sh` | Builds the original code and the harness, then runs it. |
| `convert.py` | Converts the dumps to `shadowinfo` objects and cross-checks them. It also contains reference replicas of the oracle and the Java string formats. |

The background for all of this is in `docs/notes/original_java.md`, which inventories the
Java code and lists its quirks B1–B18, and in `docs/notes/original_maps.md`, which matches
the maps and paths to the paper figures.

## Regenerating the fixtures

1. Get the original sources from <https://github.com/arc-l/shadow-information-space> at
   commit `1a52234199b03749ca9eae6630f7e2f02494630d`. The committed fixtures, `index.json` and
   the byte-identical `.dat` copies in `shadowinfo/polygon/data/` were all made from that commit
   (it is also recorded in the `provenance` of `paths.json`). The harness needs the `source/`
   folder, which contains `geometry/`, `polygons/` and `config/`.

   ```sh
   git clone https://github.com/arc-l/shadow-information-space
   git -C shadow-information-space checkout 1a52234199b03749ca9eae6630f7e2f02494630d
   ```

2. Get `log4j-1.2.12.jar`, for example from Maven Central at
   `https://repo1.maven.org/maven2/log4j/log4j/1.2.12/log4j-1.2.12.jar`.
3. Use a JDK 11 or later. The fixtures in the repo were made with Temurin 17.0.20.1.

```sh
ORIG_SRC=/path/to/shadow-information-space/source \
LOG4J_JAR=/path/to/log4j-1.2.12.jar \
JAVA_HOME=/path/to/jdk-17 \
tools/java_reference/run.sh                      # all 14 polygons + all paths, about 35–45 s
python tools/java_reference/convert.py           # cross-check report (add --json report.json)
```

`run.sh` runs these commands:

```sh
javac -nowarn -encoding ISO-8859-1 -cp log4j-1.2.12.jar -d build/orig $(find $ORIG_SRC/geometry -name '*.java')
javac -encoding UTF-8 -cp build/orig:log4j-1.2.12.jar -d build/harness tools/java_reference/src/*.java
java -Djava.awt.headless=true -XX:-OmitStackTraceInFastThrow \
     --add-opens java.base/java.lang=ALL-UNNAMED -Xss16m \
     -cp build/harness:build/orig:log4j-1.2.12.jar \
     GoldenDump $ORIG_SRC tools/java_reference/paths.tsv tests/fixtures/java
```

- **`-Djava.awt.headless=true`** keeps AWT from opening a window. Only `java.awt.geom` is
  used.
- **`--add-opens`** is needed to seed `Math.random()` (see below).
- **`-XX:-OmitStackTraceInFastThrow`** keeps stack traces on exceptions that are thrown
  repeatedly, so that crashes can be located.

The log4j "no appenders" warnings on stderr are harmless.

**`GoldenDump` options:**
- `--polys 1,5`
- `--runs PP@12,PP5@14`
- `--seeds 1,2,3,4,5`
- `--fixed-seeds 1,2,3`
- `--burnins 0,1,...`
- `--no-polys` and `--no-runs`

Set the environment variable `OUT=` to write somewhere else. `index.json` is written only by
a full run.

Two full runs on the same JVM produce **byte-identical** files (checked). See
"Nondeterminism" for the one part that can differ on another JVM, or when only some runs are
selected.

## What the harness adds (everything else is the original code)

1. **Seeded `Math.random()`.** The original oracle uses unseeded `Math.random()`.
   `Math.random()` delegates to a single static `java.util.Random`, the field
   `java.lang.Math$RandomNumberGeneratorHolder.randomNumberGenerator`. The harness gets that
   object by reflection and calls `setSeed(seed)`. The **unmodified**
   `SingleTypeAgentOracle` then draws exactly `new java.util.Random(seed).nextDouble()`,
   …. `index.json → java_random_check` asserts this. `convert.java_random_doubles` is a
   48-bit-LCG clone of it in Python.
2. **Per-sample physical gaps.** `Algorithm.getGaps` keeps the gaps it computes at each
   sample to itself. The harness repeats the sampling loop, which uses only public methods:
   `path.getAllCriticalPoints(getCuts(poly))`, then
   `purturbPointAlongSeg(point, seg, purturb, true)`, then `getPhysicalGaps(poly, sample)`.
   For the ID propagation it calls `getGaps` itself. The JSON also checks that gap set *k*
   of `getGaps` is exactly the physical-gap array of sample `sample_index`
   (`sets_equal_samples`, true in every run).
3. **`RecordingPath`**, a subclass of `drawable.Path`, logs the points that `getGaps` passes
   to `ptDistFromStart`. `getGaps` does this at `i = 0` and at every GI/BT crossing.
   Together with the stack trace's `getGaps` line number, the log locates the critical point
   at which a crashing `getGaps` fails (`get_gaps.failure`).
4. **Identity-hash "burn-ins" around `deriveShadowBounds`.** Its result depends on
   `HashMap<Vertex,…>` iteration order. `Vertex` has no `hashCode`, so that order follows
   the JVM identity hashes. The harness rebuilds the bipartite graph 24 times, each after
   `k` extra `System.identityHashCode(new Object())` calls, and records every distinct
   bounds outcome.
5. **`MergeFixedOracle`**, the `merge_fixed` variant. It is the original oracle with one line
   changed: the second merging parent adds to the first parent's count instead of doubling
   its own (bug B7). This gives consistent observations, so the bounds comparison is
   meaningful. The `original` runs always use the unmodified class.

## Fixture layout (`tests/fixtures/java/`)

**Format conventions:**
- Every double is written with Java `Double.toString`. Python's `float()` parses it back to
  the identical double. The writer checks the round trip, and `non_round_trip_doubles` is 0.
- Points are `[x, y]` in the `.dat` frame, which is y-up and counter-clockwise. **Do not flip
  the data.**
- Lines are `[x1, y1, x2, y2]`.
- A physical gap is `[startEdge, endEdge, startPoint|null, endPoint|null]`.

**`index.json`** holds the JVM, ε, δ, the seeds, the burn-ins, the `java.util.Random` check,
and the list of runs.

**`poly<N>.json`, one file per polygon 1–14:**
- `vertices` as loaded. No transform is applied: the Java scales only for drawing.
- `turn[i]`: `relativeCCW(e_{i-1}, v_{i+1})`. A reflex vertex is one with +1.
- `reflex`: the list of reflex vertex indices.
- `counts`.
- `cuts`:
  - `single_tangent`
  - `inflection_nongeneral`
  - `inflection_all`
  - `general_inflection_cut`, with `fromLine` and `ccw`
  - `bitangent_cut`, with `this`, `opp`, `oppSeg` and `curveTo`
  - `get_cuts_layout`: `getCuts` is the concatenation ST + NGI + GI + BT, in this order.
  - Two flags confirming that `getInflection(GENERAL)` and `getBitangent()` repeat the
    lines of the corresponding `*Cut` functions.
- `point_in_polygon`, which tests the `Path2D.Float` non-zero rule:
  - a 19×19 lattice (step 50);
  - every vertex;
  - every edge midpoint.
- `visibility`, for every interior point of a 9×9 lattice (step 100):
  - `segInPolygon(q→v_j)` for all `j`;
  - `getPhysicalGaps(q)`;
  - for about 12 of these points, `getVisibilityPolygon(q)` with Java's duplicate
    vertices.

**`<polygon>_<path>.json`, one file per path in `paths.tsv`:**
- `critical_points`: `[distance, point, cutType, index into getCuts, path segment, perturbed sample, physical gaps at the sample]`.
  The last element is omitted for the `fig_*` runs to save space.
- `panel_intersect_points`: what the panels draw.
- `visibility`: physical gaps and visibility polygon at the waypoints and at five samples.
- `get_gaps`, if it succeeded:
  - one entry per gap set: `t` (the relative time, carried by `gaps[0]` only),
    `sample_index`, and `gaps = [id, iEdge, toGapSet in Java HashSet order]`;
  - `line`: the exact `ProjectPanel` stdout line;
  - `final_ids` and `max_id`.
- `get_gaps`, if it failed: `ok: false`, the exception, and
  `failure.critical_index` / `phase`.
- `nse`: `processGapHistoryInformation(NEVER_SEE_EVADER)` states per set. `0` = clear,
  `1` = contaminated.
- `oracle_runs`, one per seed and variant. Panel and code paths get original seeds 1–5 and
  merge-fixed seeds 1–3. `fig_*` paths get one seed of each. Each run has:
  - `events`: `[t, type, from, to, moved, visible(, toString)]`;
  - `truth_final`: the oracle's hidden counts at the end;
  - `equations` (`Equation.toString`);
  - `bg`: the structure after `deriveShadowInfoState` (left `[id, min, max, sorted targets]`,
    right `[id, min, max]`), plus `bg_text_example`;
  - `java_bounds_outcomes`: every distinct `{shadow: [min, max]}` printed by
    `deriveShadowBounds`, with the burn-ins that produced it.

Total size is about 5.2 MB. To stay near 5 MB, a few things are left out of the dumps.
Each can be recomputed from what is kept:

- the per-gap edge and point data in `get_gaps`, which equal the physical gaps of the
  sample given by `sample_index`;
- the NSE text lines, which are `Gap.toString` with the state;
- the oracle's per-set hidden counts. Only the final counts are kept; `convert.replay_oracle`
  recomputes the others.
- `toString()` of events, except in the first oracle run of each file;
- the per-sample physical gaps for the `fig_*` runs.

## Cross-checks (`convert.py`, `tests/test_java_reference.py`)

The checks run on every dump. All of them pass unless the result says otherwise.

- **Gap history:**
  - The gap history read as component events always gives exactly one event per set
    transition.
  - It validates as a `ShadowSequence`.
  - It equals the oracle's event list.
- **`replay_oracle`.** This Python replica of `SingleTypeAgentOracle.distributeAgents`,
  driven by the `java.util.Random` clone, reproduces **all 76 oracle runs** exactly. That
  covers both variants, and it covers 1 000 000 targets.
- **Text formats.** `format_event` and `format_equations` reproduce Java's strings exactly.
  The gap-set `line` is `repr(t)` followed by the gap strings, for all 1336 sets.
- **Bipartite graph.** The Java graph equals `BipartiteIState` in every run once two things
  are allowed for:
  - Java pools all initial shadows into a single left vertex `0`.
  - An appeared shadow is a single `Vertex` that is its own left and right vertex. An alive
    one therefore shows its appear count instead of −1.
- **NSE.** Every run equals "contaminated ⇔ the reach set contains an initial shadow", which
  is the counting filter.
- **ICRA'08 Fig. 4.** The polygon-12 `ProjectPanel` run, with s5 treated as initial, gives
  exactly the bipartite edges of `tests/fixtures/icra08_fig4.json` and the same event
  multiset. Two independent pairs of events are in a different order.

The bounds comparison is not a pass: see "Results the port must know about" below.

## Results the port must know about

**Java `deriveShadowBounds` is order-dependent and often wrong.**
- In every feasible run, at least one hash order printed exactly the
  `CombinatorialFilter(total=(N−v₀, N−v₀))` bounds.
- Over all runs × 24 hash orders, Java equals the exact bounds in 1242 of 1680 cases
  (about 74%).
- In the other cases the Java bound has one of these faults:
  - `min > max`;
  - a negative min;
  - looser than the exact bound;
  - it excludes feasible values, so it is wrong.
- In the `merge_fixed` runs the true count is always inside the exact bounds. It is outside
  the printed Java bound in many hash orders.
- Shadows that appeared and never split or merged get no Java bound (B11).

**The original oracle's merge bug (B7) makes observations inconsistent.**
- In 30 of 40 feasible original runs, the oracle's own final counts lie outside the exact
  bounds.
- Six original runs are infeasible: all 5 seeds of `14_PP5`, and `14_fig_TRO-Fig15b`.
  `shadowinfo` raises `InfeasibleError` on them; Java prints numbers anyway.

**`getGaps` crashes:**

| run | exception | `getGaps` line | `collapsePhysicalGaps` line | critical point |
|---|---|---|---|---|
| `14_PP4` | `NullPointerException` | 286 (bitangent merge) | 460 | 654 (BT) |
| `13_c_second`, `13_c_tail` | `ArrayIndexOutOfBounds` | 138 (GI appear) | 380 | 239 (GI) |
| `14_fig_ICRA08-Fig8b` | `ArrayIndexOutOfBounds` | 138 | 390 | 596 (GI) |

The `14_fig_ICRA08-Fig8b` path lies entirely inside the polygon. Valid paths can make the
Java gap tracking crash.

**Nondeterminism:**
- `getGaps`, NSE, the oracle (once seeded) and the equations are deterministic.
- Gap `toGapSet` order is `HashSet<Integer>` order. This is deterministic, but it is not
  sorted.
- The only JVM-dependent output is two things that follow identity hashes:
  - `java_bounds_outcomes`, meaning which outcomes were sampled and their burn-in lists;
  - the adjacency order inside `bg_text_example`.

  Compare `bg` structurally. Treat the Java bounds as a set of possible outcomes.

## Paper cross-check (T-RO Fig. 15, ICRA'08 Fig. 8)

| quantity | paper 15(a) | Java `13_c_first` (code path) | paper 15(b) | Java `14_PP5` |
|---|---|---|---|---|
| component events | 85 | **85** | 385 | **385** |
| shadow IDs | – | 123 | 491 | **491** |
| final shadows (with a bound) | 18 | 20 (**18**) | 12 | 14 (**12**) |
| BG vertices | 41 | **41** (12 L + 29 R, pooled) | 124 | 114 (51 L + 63 R, pooled); **124** with one left vertex per initial shadow (`BipartiteIState`) |
| BG edges | 60 | 80 pooled / 156 per-initial | 339 | **339** pooled / 396 per-initial |

The paths digitised from the figures give 83 / 22 (`13_fig_TRO-Fig15a`) and 381 / 18
(`14_fig_TRO-Fig15b`), as events / final shadows. The code paths are the ones that reproduce
the paper's numbers. The `14_fig_ICRA08-Fig8b` digitisation crashes `getGaps`.
