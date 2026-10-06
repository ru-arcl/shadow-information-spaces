# Original Java maps and robot paths vs. the paper figures

This note matches the 14 polygons in the original Java code
(`source/polygons/{1..14}.dat`) and the robot paths hard-coded in
`source/geometry/pe/ui/ProjectPanel*.java` to the figures in the three papers.
The applets are out of scope.

Scratch artefacts, which are not part of the repo:

- `digitised_paths.json` has the digitised paths, the fit parameters and every code path with its validity check.
- `maps_vs_figures.png` is a contact sheet showing each paper crop, the matched `.dat` with the code path and the digitised path, and an overlay.
- `orig_maps_yup.png` and `orig_maps_ydown.png` show all 14 maps with their code paths.

In the repo, `examples/original_maps.py` draws all 14 maps with every stored path (code, digitised,
demo) in the `.dat` frame and prints the map and path tables (`docs/assets/original_maps.png`). The
paths and their provenance are in `shadowinfo/polygon/data/paths.json`.

## Coordinate frame

- A `.dat` file holds one vertex `x y` per line in `[0,1000]^2`. `Polygon.read` does no scaling.
- `Polygon.draw` and `Path.draw` paint `screen_y = Y_MAX - y`, with `Y_MAX = scaling-factor = 1000` from `config/geometry.properties`. So the `.dat` frame is **y-up**, and the Java window shows the map y-up. `orig/shadow.png` confirms this: it is a Java screenshot of 1.dat with the ProjectPanel2 path.
- Every robot path in the code uses the same `.dat` frame. All digitised paths below are given in that frame.
- "Orientation" in the table describes how the paper figure shows the `.dat`:
  - `yup` means the figure looks like the Java window.
  - `ydown` means the figure is a vertical mirror of the Java window, with the raw `.dat` y drawn downwards.
  - Mirroring or transposing the map does not change the shadow topology, but it reverses the left/right labels. Use the `.dat` frame as the source of truth.

## Robot paths in the Java code

| class | polygon | path | status | inside the polygon? |
|---|---|---|---|---|
| `ProjectPanel` | 12 | `(90,450) (580,450) (580,840) (520,840)` | active | yes, clearance ≥ 20 |
| `ProjectPanel` | 13 | 9 points `(195,620) … (460,440)` | commented | yes, clearance ≥ 34 |
| `ProjectPanel` | 13 | "Second piece" `(69,418) (60,575) (356,620)` plus a doubly commented tail | commented | **no**: `(69,418)` is outside 13.dat. The piece fits inside 1.dat, 7.dat and 9.dat, so it is a leftover from another map |
| `ProjectPanel*` | 1 | `(886,281) (886,48) (556,36)` | commented stub | yes |
| `ProjectPanel2` | 1 | 8 points `(930,886) … (905,95)` | active | yes |
| `ProjectPanel2` | 1 | 14 points, a longer variant | commented | yes |
| `ProjectPanel3` | 5 | 14 points `(946,501) … (93,123)` | active | yes, clearance ≥ 8.8 |
| `ProjectPanel4` | 14 | 21 points `(89.29,625) … (678.57,687.5)` | active | yes, clearance ≥ 20 |
| `ProjectPanel5` | 14 | same as PP4, except v19 = `(633.93,712.29)` and v20 = `(688.57,657.5)` | active | yes, clearance ≥ 20 |

`ProjectPanel5` builds `SingleTypeAgentOracle(1000000)`. That is the "1 million targets" of the Fig. 15(b) run.

Polygons 2, 3, 4, 6, 7, 8, 9, 10 and 11 have no path in the code and appear in no paper figure.

## Figure → map table

Method:

1. Render each map in the 8 dihedral orientations.
2. Fit it to the black-pixel distance map of the paper crop with a per-axis scale and offset (Nelder–Mead).
3. Brute-force all 14 maps × 8 orientations for each figure.

The matching map/orientation was unique every time. The runner-up always had a forward RMS error ≥ 7 units and a reverse 95th-percentile error ≥ 31 units.

The "fit RMS" column is the RMS distance from the `.dat` boundary to the drawn walls, in `.dat` units. Vector figures from ICRA08 were rendered at 600 dpi. T-RO figures are 150 ppi JPEGs and were read at native size.

| figure | map | orientation | fit RMS | path source | paper-reported numbers for this run |
|---|---|---|---|---|---|
| T-RO Fig. 11(a) = ICRA08 Fig. 2(a) | 12.dat | ydown (T-RO squeezes x to 0.58) | 3.2 (T-RO), 6.3 (ICRA08)\* | Same route as the active `ProjectPanel` polygon-12 path. **Digitised**, because the drawn coordinates differ (see below) | 5 shadows at `t0`. Appear 10, 12, 18; split 1→6,7, 7→9,8, 3→13,14; merge 2+6→15, 4+10→11, 11+5→16, 8+16→17, 17+12→19; disappear 9, 14; final 13, 15, 18, 19. Bipartite: 8 left (1,2,3,4,5,10,12,18) + 6 right (9,14,13,15,18,19), 11 edges. Bounds on s19: [10, 24] (`fig11_extraction.md`) |
| ICRA08 Fig. 3(a)–(d) | 12.dat | ydown | 6.3\* | Same as Fig. 2(a). Digitised from 3(a); the start x is ≈115 instead of ≈80 | Illustrates split, appear, merge and disappear only. No counts |
| T-RO Fig. 15(a) = ICRA08 Fig. 8(a) | 13.dat | **yup** | 0.0 (ICRA08), 0.9 (T-RO) | Commented `ProjectPanel` polygon-13 path, **first 9 points only**. **Digitised**: vertices 3–4 and 7–8 are 18–24 units off | 85 component events; 100 targets placed by the oracle; bipartite graph 41 vertices / 60 edges; 18 final shadows; 0.1 s for one team |
| T-RO Fig. 15(b) = ICRA08 Fig. 8(b) | 14.dat | ydown | 0.2 (ICRA08), 0.9 (T-RO) | **Code**: `ProjectPanel4`/`5` polygon-14 path, vertices 0–19 (median deviation 0.7 units on T-RO, 2.8 on ICRA08; the latter is a systematic offset under 1 line width). The final vertex v20 is not drawn | 385 component events; 491 total shadows; bipartite graph 124 vertices / 339 edges; 10⁶ targets in 5 teams; 12 final shadows; under 1 s (T-RO) / 2.5 s (ICRA08) |
| T-RO Fig. 1, 3, 4; ICRA08 Fig. 1; ICRA10 Fig. 1 | none | n/a | n/a | n/a | Schematics: a rectangle with a slot, spotlight discs, rectangle obstacles, photos |
| T-RO Fig. 6 = ICRA10 Fig. 2 | none | n/a | n/a | n/a | FOV-event schematic |
| T-RO Fig. 12 = ICRA10 Fig. 4 | none | n/a | n/a | n/a | Hand-drawn polygon **with a hole**. No `.dat` has a hole. The event sequence is given in the figure itself, on the left |

\* The 12.dat residual comes from a single wall. The paper drawing puts vertices 21–22 (`(285.71,571.43)–(285.71,1000)`) at **x ≈ 250**: measured 249.8 on ICRA08 Fig. 2(a), with the same offset in Fig. 3 and T-RO Fig. 11(a). Every other wall lies on the 1000/7 and 1000/28 grid within about 1 unit. The figure was redrawn by hand rather than taken from the Java window. The orientation also shows this: the figure is ydown, while the Java window draws the map y-up.

## Digitised paths (`.dat` frame, start → arrow tip)

All digitised polylines pass a scratch check: every vertex is inside the polygon (ray casting), no segment crosses a polygon edge, and 41 points sampled per segment are all inside.

| figure | digitised polyline | min clearance |
|---|---|---|
| ICRA08 Fig. 2(a) | `(79.8,467.2) (611.8,466.9) (612.8,926.8)` | 31 |
| T-RO Fig. 11(a) | `(86.3,464.5) (612.1,464.5) (607.9,933.0)` | 31 |
| ICRA08 Fig. 3(a) | `(114.8,467.5) (612.0,467.1) (612.6,926.2)` | 31 |
| ICRA08 Fig. 8(a) | `(196.9,614.2) (196.6,508.4) (615.4,508.5) (614.8,642.0) (757.5,642.3) (757.5,473.1) (578.5,473.3) (579.9,419.1) (448.7,420.0)` | 37 |
| T-RO Fig. 15(a) | `(196.7,610.8) (196.0,507.7) (617.6,507.9) (616.9,642.2) (757.4,642.8) (758.2,473.5) (580.1,473.4) (579.8,418.7) (457.2,420.5)` | 37 |
| T-RO Fig. 15(b) / ICRA08 8(b) | 20 vertices, equal to PP4 v0–v19 within line width (see JSON) | 19 |

Residuals of the green pixels to the fitted polyline: 95th percentile ≤ 0.8 px on the T-RO JPEGs and ≤ 6.6 px at 600 dpi. That is about 8 `.dat` units or less, roughly the drawn line width.

Comparison with the code paths:

- **12.dat.** The code goes `(90,450)→(580,450)→(580,840)→(520,840)`. The figure goes `(≈80,≈466)→(≈612,≈466)→(≈612,≈927)` with no final hook. The horizontal leg is 15–17 units lower, the vertical leg 32 units further right, and it ends 87 units further along. The route is the same corridor sequence.
- **13.dat.** The code y-values 660 and 440 are drawn at ≈642 and ≈420. Everything else agrees within about 6–12 units.
- **14.dat.** The figure is the code path. The arrow tip `(629.6,714.8)` is PP4 v19 `(633.93,714.29)`.

## Implications for the Python port

- To reproduce the Fig. 15 runs, load 13.dat and 14.dat unchanged in the `.dat` frame and use the code paths: 13 with the first 9 points (`P13`); 14 with the `ProjectPanel5` path (`P14b`). `examples/paper_fig15.py` does this and compares every count with the numbers above: events, labels and bounded final shadows match exactly; the bipartite-graph sizes follow two different conventions (see `original_java.md` §1).
- The paper's 13.dat figure is y-up, while the 12.dat and 14.dat figures are mirrored. Do not flip the polygons or paths to match a figure.
- For Fig. 11, the code path on 12.dat is the reference. The paper's polygon has one wall shifted (x = 250 instead of 285.71) and its path is drawn by hand. `examples/fig11_maxflow.py` runs it and recovers the Fig. 11(b) events, the Fig. 11(c) bipartite I-state and the bounds [10, 24] on s19.
