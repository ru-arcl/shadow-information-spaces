# Conference papers vs. the journal: what is extra or different

Sources:

- **ICRA08**: J. Yu, S. M. LaValle, *Tracking Hidden Agents Through Shadow
  Information Spaces*, ICRA 2008, pp. 2331–2338. Read in full.
- **ICRA10**: J. Yu, S. M. LaValle, *Probabilistic Shadow Information Spaces*,
  ICRA 2010, pp. 3543–3549. Read in full.
- **T-RO**: J. Yu, S. M. LaValle, *Shadow Information Spaces: Combinatorial
  Filters for Tracking Targets*, IEEE T-RO 28(2), 2012, pp. 440–456.

T-RO is the reference for this repository (`docs/DESIGN.md`). This note lists
what the two conference papers contain that T-RO Secs. II–VII leave out or state
differently. "Implication" says what that means for the code.

New fixtures from this pass are in `tests/fixtures/`, listed in §6. Expected
values were computed with an independent integer program (scipy `milp`) or an
independent exact filter. These scripts live outside the repo, so the numbers
are not produced by `shadowinfo`.

---

## 1. Plain-language summary

A robot moving through an environment with sensors sees only part of it at any
moment. The part it cannot see, the *shadow region*, splits into connected
pieces (*shadows*). Shadows appear, disappear, split and merge as the sensors
move. We observed something about targets that move unpredictably
and are seen only when they cross into view: to track them, you only need the
order of these shadow events and the times targets enter or leave a shadow
across the edge of the field of view. Nothing else in the sensor history
matters. This compressed record is called a *shadow information space*, and
filters that work on it are *combinatorial filters*.

**Nondeterministic case.** Here any target motion consistent with the
observations is possible. The possible numbers of hidden targets per shadow
form an integer program whose constraint matrix is totally unimodular. So exact
lower and upper bounds can be computed in polynomial time. In practice they
come from a maximum-flow computation on a small bipartite graph that links
where targets could have come from to where they could be now.

**Other tasks.** The same machinery handles:

- counting targets whose total number is unknown;
- keeping the status of a passive pursuit–evasion search;
- tracking teams of partially distinguishable targets.

**Probabilistic case.** When target motion and field-of-view observations are
probabilistic, the filter carries the exact joint distribution of target counts
over the current shadows through each event. That distribution can grow
exponentially, so truncation heuristics similar to particle filtering are used
for large instances.

**The three papers.** ICRA 2008 introduced the combinatorial filter. ICRA 2010
added the probabilistic version. The 2012 IEEE Transactions on Robotics article
unifies and extends both.

---

## 2. BibTeX

Each DOI was checked against Crossref. Pages and venue were checked against the
PDF headers and footers.

```bibtex
@inproceedings{YuLav08ICRA,
  author    = {Yu, Jingjin and LaValle, Steven M.},
  title     = {Tracking Hidden Agents Through Shadow Information Spaces},
  booktitle = {Proceedings of the IEEE International Conference on Robotics and Automation (ICRA)},
  address   = {Pasadena, CA, USA},
  month     = may,
  year      = {2008},
  pages     = {2331--2338},
  doi       = {10.1109/ROBOT.2008.4543562}
}

@inproceedings{YuLav10ICRA,
  author    = {Yu, Jingjin and LaValle, Steven M.},
  title     = {Probabilistic Shadow Information Spaces},
  booktitle = {Proceedings of the IEEE International Conference on Robotics and Automation (ICRA)},
  address   = {Anchorage, AK, USA},
  month     = may,
  year      = {2010},
  pages     = {3543--3549},
  doi       = {10.1109/ROBOT.2010.5509588}
}

@article{YuLav12TRO,
  author    = {Yu, Jingjin and LaValle, Steven M.},
  title     = {Shadow Information Spaces: Combinatorial Filters for Tracking Targets},
  journal   = {IEEE Transactions on Robotics},
  volume    = {28},
  number    = {2},
  pages     = {440--456},
  month     = apr,
  year      = {2012},
  doi       = {10.1109/TRO.2011.2174494}
}
```

Dates: ICRA08 ran May 19–23, 2008. ICRA10 ran May 3–8, 2010. The T-RO article
appeared online on Dec. 23, 2011, in the April 2012 issue. ICRA10 footnote 1
points to a longer version at `http://msl.cs.uiuc.edu/~jyu18/icra10/shadow.pdf`,
which is probably no longer online.

---

## 3. ICRA 2008: deltas

### 3.1 Formulation

| # | ICRA08 | T-RO | Implication |
|---|---|---|---|
| 1 | Bodies are called *agents*. Agent `a` belongs to team `l(a) ∈ T`, `|T| = m`. Agents can be told apart only if they are in different teams. Example: 30 agents, 10 red, 15 green, 5 blue. | Bodies are called *targets*. *Location* (binary detector vs. counter) and *identity* (fully, partially or not distinguishable) are treated as separate attributes. | Teams are the only identity model we need. `MultiTeamFilter` = one filter per team. |
| 2 | The sensor gives `y ∈ ℕ^{m+1}`. Entries 1..m count visible agents per team. Entry m+1 is the label of the shadow entered or exited at a FOV event, and 0 otherwise. | No explicit observation vector; FOV events are attached to shadows. | Gives a concrete per-tick format for the grid simulator: per-team visible counts plus `(shadow, enter/exit)` pairs. |
| 3 | Initial condition per shadow: `{(T_1,l_1,u_1),…,(T_k,l_k,u_k)}`. Each `T_i` is a set of teams and the `T_i` are disjoint. Examples: "at least 2 and at most 5 red" → `({red},2,5)`; "exactly 7, teams unknown" → `(T,7,7)`; "completely unknown" → `(T,0,∞)`. | The same in Eq. (1), with attribute sets `a_i`. Example `{(green,6,9),(blue or red,5,5)}`. | None (these are extra examples). |
| 4 | General-position footnote: "three or more components involved in a change" is excluded. | "Four or more" components in an evolve event are excluded. | ICRA08's wording is a slip, since a split already involves 3 components. Follow T-RO. |
| 5 | Appearing shadows must have `v_i = (0,0)`: "a new shadow component must be clear of agents". | Appear and disappear carry an arbitrary `d_i` (Eq. 3). ICRA10 justifies this: targets may slip in during the appear, and a sub-team may sweep the shadow. | `Appear(s, lo, hi)` with `lo, hi > 0` is the general case. For geometric simulators, `(0,0)` is the natural default. |

### 3.2 Extra tasks (ICRA08 Sec. II-D, intro, Sec. VI)

ICRA08 lists these queries. T-RO Sec. V-E keeps only the bound, counting and
pursuit–evasion tasks.

- Minimum and maximum number of red (or any colour) agents in a final shadow.
  T-RO has this too.
- Total number of agents in `F` (counting). T-RO has this too.
- Whether a final shadow contains any agent at all, i.e. whether its lower
  bound is > 0 or its upper bound is 0.
- **Whether a final shadow holds more blue agents than red ones.** With per-team
  bounds, "certainly more blue" is `lo_blue > hi_red`. "Possibly more" needs a
  joint check, because the two teams' flows are independent only if no
  constraint couples them.
- **Whether blue and red agents are separated into different shadows.**
  Certainly separated if no shadow has `hi_blue > 0` and `hi_red > 0` at the
  same time.
- **Monitoring team movement** "to ensure that sufficient numbers from each team
  are present in critical parts of an environment" (intro, item 3). This is the
  only concrete reading of the "herding" task that the T-RO abstract mentions.
  The word "herding" appears nowhere else in the three papers.
- **Recognising task completion.** Example: in a wildlife preserve with known
  `n`, the census is done once the lower bound on the total reaches `n`. T-RO
  folds this into "counting".

Implication: these are thin helpers on top of per-team `all_bounds()`. They
could go into `nondeterministic.py` (`certainly_more`, `certainly_separated`) if
the coordinator wants them. They are optional; T-RO does not require them.

### 3.3 A FOV shortcut that T-RO omits (ICRA08 Sec. IV)

ICRA08 handles two kinds of FOV events without extra component events:

1. If `k` agents enter a shadow *right after it appeared*, set its appear bound
   to `v_i = (k,k)`.
2. If `k` agents leave a shadow and *the next event is its disappearance*, set
   its disappear bound to `v_i = (k,k)`.

All other FOV events are converted to component events (ICRA08 Sec. V-C). An
agent appearing in view becomes split followed by disappear. An agent vanishing
from view becomes appear followed by merge. T-RO Sec. V-D uses the same
conversion. T-RO then adds per-shadow batching (`d_min`, `d_tot`), which
ICRA08 does not have.

Implication: the shortcut fits the grid simulator well. A sensor sweep that
reveals a shadow typically produces a run of `Exit` events followed by
`Disappear`. Folding the run into `Disappear(s, k, k)` keeps the bipartite
graph small. It must give the same bounds as the batched form, which makes it a
good extra equivalence test.

### 3.4 Max-flow network as ICRA08 states it (Sec. V-B). Superseded.

- **Vertex names.** Incoming = initial or appeared (left side). Inactive =
  disappeared (right side). Active = alive at the end (right side).
- **Two sinks.** Source `S0` feeds the incoming vertices. Inactive vertices go to
  sink `S1`. Active vertices and `S1` go to a second sink `S2`.
- **Exact weights first.** The derivation assumes `l_i = u_i = w_i`.
  - *Minimum of `s_m`:* set capacity 0 on `m`'s edge to the sink, run max-flow,
    and read off the deficit, Eq. (1):
    `min f(m) = Σ_i f(S0,i) − Σ_j f(j,S1) − f(S1,S2)`.
  - *Maximum of `s_m`:* give `m` infinite capacity, give the other active
    vertices 0, and read `f(m,·)`.
- **Intervals.** ICRA08 then says: "replace `v_i` by `l_i` for minima and `u_i`
  for maxima".

Problems:

- Text and figure disagree. The text sends inactive vertices to `S1` and active
  ones to `S2`. Fig. 7(b) draws active → `S1` → `S2` and inactive → `S2`.
- Fig. 7(b) merges initial shadows 1–4 into one source vertex "1-4", "because we
  can equivalently assume that they are split from same `v_i`". That only holds
  if the initial knowledge is a single bound on `x1+x2+x3+x4`. With separate
  per-shadow bounds it loses information. T-RO Fig. 11(d) keeps them separate.
  **Do not implement the merging.**
- The "`l` for minima, `u` for maxima" rule breaks when the interval bounds are
  not tight. Example: `s1 ∈ [2,10]` splits into `a, b`, and `b` disappears
  revealing exactly 3. Using `l` gives source 2 and sink 3, so the network
  cannot saturate and Eq. (1) gives `2 − 3 = −1`. The true minimum of `a` is 0.
  T-RO Sec. IV-A works around this by assuming tight initial conditions
  temporarily. The sibling note `fig11_extraction.md` shows that the literal
  T-RO lower-bound recipe has a similar problem. The flow-with-lower-bounds
  formulation in DESIGN §3.2 is the right reference.

### 3.5 Refining initial bounds (ICRA08 Sec. VI, Eqs. 2–3, with a worked example)

**Worked example (text only, no figure numbers).** `s1` starts with a lower
bound of 4. Later, 6 agents are revealed when `s9` disappears, and `s9`
descends only from `s1`. So `s1` must have started with at least 6. This is now
the `refine_lower_bound` case in `tests/fixtures/icra08_tasks.json`, built on
the Fig. 4 sequence (§6). Expected refined `s1 = [6, 8]`.

**ICRA08 recipe for the lower bound:**

- `c(S0,m) = l_m`;
- every other `c(S0,i)` = that vertex's weight;
- inactive edges weighted as before, finals infinite;
- `c(S1,S2) = Σ c(S0,i) − Σ c(j,S1)`;
- result: `l'_m = l_m + Σ_j c(j,S1) − Σ_j f(j,S1)`, with `j` over inactive
  vertices.

**ICRA08 recipe for the upper bound:**

- `c(S0,m) = u_m` and `c(S0,i) = l_i` for `i ≠ m`;
- result: `u'_m = Σ_j f(m,j)`, the total flow *out of* `m`.

**T-RO (Eqs. 9–10).** For the lower bound, the other sources and the
disappearing shadows are set to `u_i` and the rest of the sinks to 0. For the
upper bound, the result is `u'_1 = f(S,1)`. Neither paper defines "weight" for
an interval-bounded vertex.

Implication: implement `refine_initial_bounds` as min/max of `x_i` under the
feasible-flow model, as DESIGN says. The fixture gives the expected values.

### 3.6 Counting and pursuit–evasion (ICRA08 Sec. VI)

**Counting.** Set every initial shadow to `(0, ∞)`. If refinement makes
`l_i = u_i` for all of them, `n` is known; otherwise Eqs. (2)/(3) bound `n`.
The upper bound stays ∞ while any part of `F` is unexplored. T-RO says the
same. Fixture case: `counting`.

**Pursuit–evasion.** Each initial shadow is `(0,1)`. Each final shadow ends in
one of three states: `(0,0)` cleared, `(1,1)` evader certainly there, `(0,1)`
unknown. Neither paper states the single-evader total constraint
`Σ = 1` explicitly; DESIGN adds it as `total=(1,1)`. Fixture cases:
`pursuit_evasion_not_found` and `pursuit_evasion_found`. In the second, the
evader is revealed when `s14` disappears. Every final shadow is then cleared and
the evader's initial shadow is pinned to `s3`.

### 3.7 Teams with cross-team constraints (ICRA08 Sec. VI, last paragraph)

ICRA08 gives a two-run recipe. For a constraint like "3–7 agents that are red
or blue": "set red capacity from `S0` … to zero, and flow the reds … to obtain a
minimum of the red. We then do the same for blue. Adding those two numbers up
gives the answer."

The recipe contradicts itself: it zeroes the red capacity and then flows reds.
T-RO Sec. V-F replaces it with four max-flow runs:

1. red = 0, blue = `l_i`;
2. then the other way round;
3. check `lr1 + lb1 = lr2 + lb2`;
4. report the red count as lying between `lr1` and `lr2`.

Implication: ICRA08 adds nothing to the T-RO version. The T-RO procedure is
available as `RedOrBlueFilter.two_pass()`, but it is not sound: when per-team
observations show that a mixed group held both colours, a pass contradicts
them and the surviving pass can report a range that misses the true count.
`RedOrBlueFilter.bounds()` instead solves one joint flow network over all
red/blue assignments and is exact (DESIGN §3.2).

### 3.8 Implementation numbers (ICRA08 Sec. VII; T-RO Sec. VII-A)

| | ICRA08 | T-RO |
|---|---|---|
| Platform | Java 1.5, Intel U2500 1.2 GHz, 1 GB | Java 1.6, Core 2 Quad 3.0 GHz, 1.5 GB JVM |
| Fig. 8(a) = T-RO Fig. 15(a) | 85 component events; 100 agents; bipartite 41 vertices / 60 edges; 18 final shadows, single team: 0.1 s | same, 0.1 s |
| Fig. 8(b) = T-RO Fig. 15(b) | 385 component events; 491 shadows; bipartite 124 vertices / 339 edges; 10⁶ agents in 5 teams; 12 final shadows: **2.5 s** | same instance: **under 1 s** |
| Max-flow | Edmonds–Karp, `O(VE²)`. `V` = number of initial, appearing, disappearing and final shadows. | Same, plus a complexity analysis (`O(n⁵)` EK, `O(n³)` FIFO push-relabel, `O(n³ m)` for `m` teams). |

ICRA08 is clearer about one point. The simulator did *not* simulate agent
motion. An "oracle randomly distributed 100 agents in the free space as the
component events occur". Our grid simulator, which uses random walks, is a
different and more concrete ground-truth model. That difference is fine.

### 3.9 Open problems ICRA08 lists that T-RO dropped

ICRA08 names weaker sensor models as open questions:

1. a binary sensor that reports only whether at least one agent is visible;
2. a counter that reports only the number of visible agents, without the
   shadow involved;
3. uncertainty about *which* shadow an agent entered or exited.

T-RO Sec. III-A mentions binary vs. counting location sensors but does not
solve these cases. The probabilistic FOV model in T-RO Sec. VI covers only
confusion between enter, exit and null on a *known* shadow. Possible "future
work" items for the README.

---

## 4. ICRA 2010: deltas

### 4.1 Formulation

- **Component events are exact.** ICRA10 assumes the components and component
  events "are either provided or can be efficiently obtained, with 100%
  accuracy". Imprecise component events are left to future work. T-RO VI-F
  sketches the fix: keep a distribution over consistent shadow sequences.
- **Timing of null events.** "It may appear that time information associated
  with null observations is lost; this can be fixed by recording the time
  duration of a null event and introducing time parameter in the sensors'
  statistical model." *Not in T-RO.* This is a possible extension: a
  duration-dependent `P(e | y_n, Δt)`.
- **Three regimes** (Sec. III-D, *not in T-RO*):
  1. Few agents and few events: run the exact filter.
  2. Few agents and many events: the exact filter works only with very reliable
     split and sensor statistics. Otherwise the nondeterministic filter "can be
     a better alternative".
  3. Many agents: use heuristics.

  Scale quoted for the non-optimised Java code: the exact filter handles "tens
  of agents and events combined" before running out of memory (10+ min). The
  heuristics handle "up to a thousand" agents and events. This framing is good
  for the README and the demo: the demo should show the nondeterministic bounds
  first and the probabilistic expectation as an extra.
- **Exact vs. distributional counts.** ICRA10 Table II/III carry *exact*
  counts:
  - an appear sets `o.s_e = n_e`;
  - a disappear removes every entry with `o.s_v ≠ o.n_v`.

  The prose (Sec. III-B), the Fig. 4 example (`a5 ∈ {1,2}`), and T-RO Table II
  use *distributions* `P(s_e = n_e)` and `P(s_v = n_v)`. Implement the
  distributional form; the exact form is the special case.

### 4.2 FOV update: per-entry renormalisation (important)

ICRA10 Table IV, `PROCESS_FOV_EVENT`, matches T-RO Algorithm 2. For each joint
entry `p_j` with `x_i = 0` (exit impossible), the `else` branch reads:
"normalize `P'`, `P''` such that `P' + P'' = p_j`". So the entry keeps its full
mass, split between the enter and null outcomes. No mass is dropped globally.

The prose in both papers ("the affected probability mass needs to be removed
and the remaining renormalized") reads more like global dropping. That is not
what produces the published numbers:

| `P(s4 = 0, 1, 2)` after the Fig. 4 / Table III sequence | value |
|---|---|
| T-RO Table III (exact) | 0.0769 / 0.1538 / 0.7692 |
| per-entry renormalisation (Alg. 2) | **0.0769 / 0.1538 / 0.7692**; every intermediate row of Table III matches |
| drop impossible mass, renormalise globally | 0.0833 / 0.1667 / 0.7500 |
| paper Monte Carlo, 1000 trials | 0.079 / 0.154 / 0.767 |
| our rejection MC (discard trial if an exit is sampled on an empty shadow), 400k trials | 0.0838 / 0.1671 / 0.7491 |
| our MC that resamples enter/null instead, 400k trials | 0.0774 / 0.1543 / 0.7684 |

DESIGN §4 adopts per-entry renormalisation, and `monte_carlo` samples each FOV
event among the events feasible at the current count, so exact and Monte Carlo
agree on instances with FOV exits.

Typos in ICRA10 Table IV that T-RO fixes:

- `P'''` uses `P(e = e_e | y)`; it should be `e_x`.
- Both papers' prose say "If an *enter* event is not possible … the two entries
  left need to be renormalized"; it should say *exit*.

### 4.3 Worked numbers that also appear in T-RO

These are now fixtures (§6):

- **Table I merge.** Five equally likely `(s1,s2,s3)` entries; `s2, s3` merge
  into `s4`, giving `P(s1=1, s4=5) = 0.6` and `P(s1=2, s4=4) = 0.4`. T-RO
  Table I is identical.
- **Fig. 4 example.** 2 agents in each of `s1, s2`; binomial(½) split; no null
  observations; true-positive rate 0.9; `a5 ∈ {1,2}` with probability ½ each.
  - The sequence (bold observations only) is: `y_x(s1)`, split `s2→s3,s4`,
    `y_e(s3)`, merge `s1+s3→s5`, disappear `s5`.
  - The "original sequence has only 10 entries" before the merge. Reproduced:
    10.
- **"Slightly more complicated" variant.** 5 agents each in `s1, s2`, plus six
  lightened observations, reported as "135 joint probability table entries
  before the merge step". T-RO Fig. 14 gives the axes as `s1∈[1,9]`,
  `s3∈[0,7]`, `s4∈[0,5]`. **Reproduced** once both `y_e(s4)` observations are
  processed before the merge; see `probabilistic_extraction.md` §6 and
  `tor_fig12_complex.json` (tested in `tests/test_probabilistic.py`).
  - Reading the six extra observations from the figure in vertical order:
    - before the split: `y_e(s1)`, `y_x(s2)`;
    - after `y_e(s3)`: `y_x(s1)`, `y_e(s1)`, `y_e(s4)`;
    - after the merge: `y_e(s4)`.
  - Under the no-null 0.9 model, in strict figure order, this gives **110**
    entries with ranges `[1,9] × [0,7] × [0,7]`. Moving the last `y_e(s4)`
    before the merge (it commutes with the merge) gives exactly **135**.
  - The caption's `s4∈[0,5]` is still unexplained: `s3 ≤ 7` needs `s2 = 6`,
    which allows `s4` up to 6. Do not test the `s4` range.

### 4.4 Simulation (ICRA10 Sec. V-C) vs. T-RO Sec. VII-B

- **Sequence.** ICRA10 Fig. 5 = T-RO Fig. 16, decoded below. It has 14
  component events and 20 shadows. Both papers add 32 FOV observations "not
  marked in the figure". **The FOV events and the split rule used are given in
  neither paper**, so the expected counts cannot be reproduced.
- **Counts.** ICRA10 uses the large configuration only: shadows
  1, 2, 3, 11, 13, 16, 20 hold 25, 22, 23, 9, 8, 15, 9. T-RO adds a small
  configuration (10, 7, 8, 9, 6, 9, 4) for its Table IV and reuses the large
  one for Table V.
- **Results.** ICRA10 Table V:
  - exact: out of memory after 10 min, with more than 2×10⁷ entries at peak;
  - Monte Carlo: 18.26 / 22.25 / 11.85 in 43.2 s;
  - RT-100000: 17.78 / 21.83 / 12.34 in 66.3 s.

  T-RO Table V repeats these numbers, adds standard deviations (0.03 / 0.03 /
  0.02), and adds RT-LA-50000 (18.24 / 21.99 / 12.57, 40.6 s).
- **Naming clash.** In ICRA10, "RT-X" means *keep the X largest entries*
  (deterministic). That is what T-RO calls **TR** (basic truncation). In T-RO,
  RT means *random* truncation: keep by `p · U(0,1)`. Yet T-RO labels ICRA10's
  RT-100000 row as RT-100000 and gives it a standard deviation. Either ICRA10's
  run was in fact randomised, or T-RO relabelled it. For our code, use T-RO's
  definitions (TR, RT, RT-LA as in DESIGN §4). Do not treat the ICRA10 row as a
  TR regression target.

Decoded Fig. 5 / Fig. 16 (event order follows vertical position):

| kind | shadows |
|---|---|
| initial | 1, 2, 3 |
| splits | 2 → 4, 5; 1 → 6, 7; 7 → 8, 9; 8 → 13, 14; 10 → 15, 16; 12 → 19, 20 |
| merges | 5 + 3 → 10; 9 + 4 → 12; 11 + 6 → 17; 14 + 15 → 18 (this line hops over the 12 → 19/20 line) |
| appear | 11 |
| disappear | 13, 16, 20 |
| final | 17, 18, 19 |

Bipartite I-state edges: 1–13, 1–17, 1–18, 1–19, 1–20, 2–16, 2–18, 2–19, 2–20,
3–16, 3–18, 11–17.

Nondeterministic bounds without the unknown FOV events:

- small configuration: `s17 ∈ [9,13]`, `s18 ∈ [0,6]`, `s19 ∈ [0,6]`, total 15;
- large configuration: `[9,26]`, `[8,38]`, `[0,30]`, total 47.

Fixture: `icra10_fig5_structure.json`. It is useful as a structural test and as
a benchmark input for the TR / RT / RT-LA heuristics. Supply our own split rule
and FOV events.

---

## 5. Figures worth redrawing for the README

All three papers are © IEEE, so redraw these rather than copy them. Most can be
generated from our own code on the fixed grid map.

| Figure | Content | Why it is useful |
|---|---|---|
| ICRA08 Fig. 2 | Pipeline: environment + robot path → shadow-sequence graph → bipartite graph → max-flow graph | The best one-picture summary. Can be regenerated from the grid simulator (map, event DAG, bipartite, flow network). |
| ICRA08 Fig. 3 | Four panels, one per component event (split, appear, merge, disappear), each caused by the robot crossing a bitangent or inflection ray in a polygon | Explains the event types. The grid demo can capture one frame per event type. |
| ICRA08 Fig. 5 | Three panels: after a merge, a hidden agent can move into the former neighbouring shadow | Shows why merges lose information. |
| ICRA08 Fig. 1 / T-RO Fig. 4 | Application scenarios: multi-robot corridors, sensor-network holes, urban helicopter view, satellite + clouds | Motivation. Text only, or a small collage of our own drawings. |
| ICRA10 Fig. 1 / T-RO Fig. 3 | "Watch lights": disc-shaped FOVs whose motion creates appear/disappear and merge/split | A clear abstract illustration. Could also be a second demo mode with disc sensors and no obstacles. |
| ICRA10 Fig. 3 / T-RO Fig. 7 | Shadow sequence with FOV arrows (`y_e`, `y_x`) and a disappear box | Shows what the filter consumes. |
| ICRA10 Fig. 4 / T-RO Figs. 12–13 | Small polygon with a hole + 5-event sequence; joint pmf drawn as balls | Probabilistic filter example. Our fixture reproduces it exactly. |

---

## 6. Fixtures written in this pass

All live in `tests/fixtures/`. Nondeterministic sequences use the DESIGN §2
JSON (`ShadowSequence.to_dict()`, `null` = +∞). All of them load with
`ShadowSequence.from_dict`.

| File | Content | Key expected values |
|---|---|---|
| `icra08_fig4.json` | ICRA08 Fig. 4 sequence in ICRA08's labelling (children of `s1` named opposite to T-RO Fig. 11) with T-RO's numeric bounds; expected bipartite graph (= ICRA08 Fig. 7(a)); the ILP rows quoted in the text | `s13 [1,3]`, `s15 [0,5]`, `s18 [5,7]`, `s19 [10,24]`; `s13+s15 [1,8]`; all finals `[16,37]` |
| `icra08_tasks.json` | Four Sec. VI cases on the same sequence: `refine_lower_bound`, `counting`, `pursuit_evasion_not_found`, `pursuit_evasion_found` (with `total`) | refined `s1 = [6,8]`; counting `s1 ≥ 3`, `s3 ≥ 1`, `s19 ∈ [7,∞)`, initial total `[4,∞)`; pursuit–evasion statuses |
| `icra10_fig5_structure.json` | ICRA10 Fig. 5 / T-RO Fig. 16 sequence in both count configurations; bipartite graph; bounds without FOV; paper expectations kept for reference only | see §4.4 |
| `icra10_table1_merge.json` | Probabilistic merge example (Table I) | `{(1,5): 0.6, (2,4): 0.4}` |
| `icra10_fig4_table3.json` | Probabilistic Fig. 4 example: every intermediate joint table (matches T-RO Table III row by row), final marginal, paper's MC, the global-renormalisation alternative, and a note on Monte Carlo semantics | `P(s4) = 0.0769 / 0.1538 / 0.7692` |

**Probabilistic fixture format.** Adopted as DESIGN §4's JSON extension and
loaded by `ProbSequence.from_dict`, `split_rule_from_dict` and
`obs_model_from_dict`:

```json
{"prob_sequence": {"labels": [1, 2], "joint": [[[2, 2], 1.0]],
  "events": [{"type": "fov", "s": 1, "y": "exit"},
             {"type": "split", "s": 2, "a": 3, "b": 4},
             {"type": "merge", "a": 1, "b": 3, "s": 5},
             {"type": "disappear", "s": 5, "dist": {"1": 0.5, "2": 0.5}}]},
 "split_rule": {"type": "binomial", "p": 0.5},
 "obs_model": {"enter": {"enter": 0.9, "exit": 0.1, "null": 0.0}, "...": "..."}}
```

- `obs_model[y][e] = P(e | y)`.
- Joint entries are `[values in label order, probability]`.
- An appear with a distribution is `{"type": "appear", "s": k, "dist": {...}}`.

`tests/test_paper_conference.py` checks the three nondeterministic fixtures;
`tests/test_probabilistic.py` checks the probabilistic ones.

---

## 7. Resolved decisions

1. The FOV update uses per-entry renormalisation, and Monte Carlo matches it
   (§4.2, DESIGN §4).
2. ICRA08's interval-bound flow recipes and its merging of shadows 1–4 into one
   vertex are not exact (§3.4). They are not ported; the flow-with-lower-bounds
   reference stands.
3. The "135 entries" claim is reproduced when both `y_e(s4)` observations come
   before the merge (§4.3); the Fig. 14 `s4` range is not.
4. The TR/RT naming differs between ICRA10 and T-RO (§4.4).
5. Not implemented (optional): the ICRA08 helpers "certainly more blue than
   red", "certainly separated", "census complete" (§3.2), and the FOV shortcut
   that folds exits right before a disappear into the disappear bound (§3.3).
