# Probabilistic filter: extraction from the papers

Sources: T-RO 2012 Sec. VI (pp. 450–455) and Sec. VII-B; ICRA 2010
("Probabilistic Shadow Information Spaces"), the conference version of
Sec. VI. ICRA 2008 has no probabilistic material.

Fixtures:

| file | content |
|---|---|
| `tests/fixtures/tor_fig12_table3.json` | Fig. 12 simple sequence, every Table III step, final pmf, Monte Carlo numbers |
| `tests/fixtures/tor_fig12_complex.json` | Fig. 12 longer sequence (5 + 5 targets), the paper's 135-entry count, our derived values |
| `tests/fixtures/tor_fig16_sequence.json` | Fig. 16 component events, both initial-count settings, all of Tables IV and V, our no-FOV exact values |

All Table III numbers were recomputed exactly with rational arithmetic. They
match the paper to every printed digit once the FOV rule of Sec. 4 below
(renormalise within each entry) is used. The verification scripts were
throwaway code and are not in the repository; the fixtures carry the results,
and `tests/test_probabilistic.py` re-derives them with `shadowinfo`.

---

## 1. Problem formulation (T-RO VI-A)

Notation. `s_i` is both shadow *i* and the random variable "number of targets in
shadow *i*". The state is a joint pmf `P(s_1, …, s_n)`, stored as a table of
entries `P(s_1 = x_1, …, s_n = x_n) = p_j`. The set of variables changes with
component events.

Assumptions:
1. Component events are observed without error.
2. Targets are indistinguishable. The initial condition is a joint pmf
   `P(s_1, …, s_n)` at `t0`.
3. A **split rule** says how the targets of a splitting shadow are distributed
   between the two children.
4. FOV observations follow a known `P(e = e | y = y)`, with `e ∈ E_FOV = {e_e, e_x, e_n}`
   (enter, exit, null) and `y ∈ Y_FOV = {y_e, y_x, y_n}`. The sensor may report an
   enter event as enter, exit or null, and likewise for the other two events.

Output: `P(s_1, …, s_m)` over the shadows alive at `tf`. Expected counts are
read off its marginals.

## 2. Component events (T-RO VI-B; Algorithm 1)

- **Appear** `s_k` with distribution `P(s_k)` (independent of the rest):
  `P'(s_1=x_1,…,s_n=x_n, s_k=x_k) = P(s_1=x_1,…,s_n=x_n) · P(s_k=x_k)`.
  The ICRA 2010 version only allows a fixed count `n_e`, i.e. `P(s_k = n_e) = 1`.
- **Disappear** `s_k`, revealing a count with distribution `P(s_k)`:
  `P'(…, s_{k-1}=x_{k-1}, s_{k+1}=x_{k+1}, …) ∝ Σ_{x_k} P(…, s_k = x_k, …) · P(s_k = x_k)`,
  followed by normalisation. Entries with weight 0 ("stale entries") are removed.
  In ICRA 2010 the count is deterministic: remove every entry with `s_v ≠ n_v`,
  then renormalise. If the total mass is 0, the observations are inconsistent
  with the retained entries. A truncation heuristic counts this as a
  **failure** (see Sec. 7).
- **Split** `s_s → s_s1, s_s2`: replace each entry by the entries the split rule
  produces: `P'(…, s_s1 = a, s_s2 = n − a) = p_j · rule(n)(a, n − a)`, where
  `n = x_s`. In the paper example each target independently enters each child
  with probability 0.5, so `rule(n)(a, n−a) = C(n, a) 0.5^n`. The text also
  mentions an "area-proportional" rule. For one split of a shadow holding `n`
  targets, the number of entries can grow by a factor of up to `n + 1`.
- **Merge** `s_i, s_j → s_k`:
  `P'(…, s_k = x_k) = Σ_{x_i + x_j = x_k} P(…, s_i = x_i, …, s_j = x_j, …)`.

### Table I (merge example), verbatim

Before the merge the shadows are `s1, s2, s3`; then `s2, s3` merge into `s4`.

| | entry |
|---|---|
| before merge | `P(s1=1, s2=1, s3=4) = 0.2` |
| | `P(s1=1, s2=2, s3=3) = 0.2` |
| | `P(s1=1, s2=3, s3=2) = 0.2` |
| | `P(s1=2, s2=1, s3=3) = 0.2` |
| | `P(s1=2, s2=2, s3=2) = 0.2` |
| after merge | `P(s1=1, s4=5) = 0.2 + 0.2 + 0.2 = 0.6` |
| | `P(s1=2, s4=4) = 0.2 + 0.2 = 0.4` |

(The fixture is `tests/fixtures/icra10_table1_merge.json`, written by the
ICRA-2010 extraction.)

## 3. Table II: the `observation` data structure used by Algorithm 1

| field | meaning |
|---|---|
| `event` | event type: one of the component events *appear, disappear, split, merge* or the FOV events *enter, exit, null* |
| `s_s` | the originating shadow in a split event |
| `s_s1` | the first new shadow after a split event |
| `s_s2` | the second new shadow in a split event |
| `s_m1` | the first shadow in a merge event |
| `s_m2` | the second shadow in a merge event |
| `s_m` | the newly merged shadow |
| `s_e` | the newly appeared shadow from an appear event |
| `P(s_e = n_e)` | probability that `s_e` contains `n_e` targets |
| `s_v` | the disappearing shadow in a disappear event |
| `P(s_v = n_v)` | probability that `s_v` contains `n_v` targets |

ICRA 2010 Table II calls the type field `t` and has `n_e` and `n_v` ("the number
of agents in an appear event" and "the number revealed in a disappear event")
where T-RO has distributions.

Mapping to our types. `Split(s, a, b)` ↔ `(s_s, s_s1, s_s2)`.
`Merge(a, b, s)` ↔ `(s_m1, s_m2, s_m)`. `Appear(s, lo, hi)` and
`Disappear(s, lo, hi)` cover `s_e` and `s_v`. For the probabilistic filter they
also need a pmf over the count. The fixtures write it as
`"dist": {"<count>": prob}`. The FOV observation is
`{"type": "fov", "s": <label>, "y": "enter" | "exit" | "null"}`, i.e.
`FovObservation(s, y)` in `probabilistic.py` (DESIGN §4).
`events.event_from_dict` handles neither `"dist"` nor `"fov"`;
`ProbSequence.from_dict` parses both.

## 4. FOV events and observations (T-RO VI-C; Algorithm 2)

The three FOV events act on a single entry with `x_s = x`:
- enter `e_e`: the entry moves to `x + 1`;
- exit `e_x`: the entry moves to `x − 1`. This is **impossible if `x = 0`**;
- null `e_n`: nothing changes.

An observation `y` on shadow `s` splits every entry into up to three pieces,
weighted by `P(e_e|y)`, `P(e_n|y)` and `P(e_x|y)`.

### Algorithm 1: PROCESSPROBABILITYMASS (verbatim, T-RO p. 451)

```
Input:  P(s1,...,sn), the initial target distribution
        Q, the queue of observation sequences
        a split rule
        P(e | y), the sensor statistics
Output: the target distribution after all observations

foreach event observation o in Q
  switch(o.event)
  case appear:
    update all P(s1 = x1,...,sn = xn) = pj entries to
      P(s1 = x1,...,sn = xn, o.se = ne) = pj * P(o.se = ne)
  case disappear:
    set P(s1 = x1,...,sn = xn) to
      Σ P(s1 = x1,...,o.sv = nv,...,sn = xn) * P(o.sv = nv)
    remove stale entries and renormalize the probability masses
  case split:
    add two new shadows o.ss1, o.ss2
    split prob. mass in o.ss into o.ss1, o.ss2 by split rule
  case merge:
    add a new shadow o.sm and set P(...,o.sm = n) to
      Σ_{n1+n2=n} P(...,o.sm1 = n1,...,o.sm2 = n2,...)
  case enter, exit, null:
    call PROCESSFOVEVENT
return the updated target distribution
```

### Algorithm 2: PROCESSFOVEVENT (verbatim, T-RO p. 451)

```
Input:  P(s1,...,sn), the target distribution
        P(e | y), the sensor statistics
        y ∈ {ye, yx, yn}, the FOV observation
        si, the affected shadow
Output: the target distribution after the observation

foreach P(..., si = xi, ...) = pj entry in the distribution
  let P'  (..., si = xi + 1, ...) = pj * P(e = ee | y = y)
  let P'' (..., si = xi, ...)     = pj * P(e = en | y = y)
  if xi > 0
    let P'''(..., si = xi − 1, ...) = pj * P(e = ex | y = y)
  else
    normalize P', P'' such that P' + P'' = pj
  remove P(..., si = j, ...) = pj entry
  store entries P', P'' and P''' if applicable
return the updated target distribution
```

ICRA 2010 Table IV is the same algorithm with a typo. Its `P'''` line uses
`P(e = e_e | y = y)` where it should use `e_x`. The T-RO prose has a related
slip: "If an **enter** event is not possible for the observation, the two
remaining entries are renormalized" should say **exit**, as the `x_i > 0` test
in Algorithm 2 shows.

### The renormalisation rule, verified

There are two ways to read "the affected probability mass needs to be removed
and the remaining values renormalized":

| rule | meaning | Table III final `P(s4 = 0,1,2)` |
|---|---|---|
| **(A) per entry**, Algorithm 2 (`zero_exit_rule = "renormalize_within_entry"`) | at `x = 0`, rescale `P', P''` so they sum to `p_j`. Each FOV observation keeps the total mass at 1 and carries no evidence about `x`. | **1/13, 2/13, 10/13 = 0.0769, 0.1538, 0.7692** ✔ |
| (B) global | drop `p_j · P(e_x|y)` at `x = 0`, then renormalise the whole table | 1/12, 1/6, 3/4 = 0.0833, 0.1667, 0.75 ✘ |

Rule **(A)** is correct. It reproduces every intermediate entry in Table III.
For example, after `split` the entry `P(s1=1, s3=0, s4=2) = 0.225` becomes
`P(s1=1, s3=1, s4=2) = 0.225`, not `0.225 · 0.9`. The text says
"Renormalization is performed in the third step for the first and sixth
entries", which are exactly the two `s3 = 0` entries.

DESIGN.md §4 adopts the **per-entry** renormalisation (A); a global
renormalisation (B) gives different numbers. In probability terms, (A) samples the event from `P(e | y, x) ∝ P(e | y) · 1[e feasible at x]`.

Degenerate case (not covered by the paper). If `x = 0` and
`P(e_e|y) = P(e_n|y) = 0`, the entry has no feasible branch. Drop it and
renormalise globally. If no mass is left, the observations are inconsistent.

The entry count can at most triple per FOV observation (enter, exit and null
pieces). Identical keys produced by different entries must be **summed**: the
paper's counts (10 entries, 135 entries) assume this.

## 5. Fig. 12 / Table III worked example (T-RO VI-D; ICRA 2010 Fig. 4)

Assumptions, verbatim in substance:
1. initially there are two targets each in `s1` and `s2`, so `P(s1=2, s2=2) = 1`;
2. the split rule is that each target goes into each of the two split shadows
   with probability 0.5 (binomial, `p = 0.5`);
3. there is no null event or observation, and the true-positive rate of every
   observation is `p = 0.9`:
   `P(e_e|y_e) = 0.9, P(e_x|y_e) = 0.1, P(e_x|y_x) = 0.9, P(e_e|y_x) = 0.1`, and every `P(e_n|·) = 0`;
4. `a5 = 1` with probability 0.5 and `a5 = 2` with probability 0.5.

### Fig. 12 event sequence

Time runs downwards. In the figure, **bold** observations belong to the simple
sequence. *Light* ones are the six extra observations of the "slightly more
complicated" sequence. The order comes from the vertical position in the
figure (500-dpi render). The figure is a raster image, but the order is
unambiguous.

| # | event | in simple sequence | in complex sequence |
|---|---|---|---|
| 1 | `y_x` on s1 (exit observation) | **bold** | yes |
| 2 | `y_e` on s1 | – | light |
| 3 | `y_x` on s2 | – | light |
| 4 | split s2 → s3, s4 | yes | yes |
| 5 | `y_e` on s3 | **bold** | yes |
| 6 | `y_x` on s1 | – | light |
| 7 | `y_e` on s1 | – | light |
| 8 | `y_e` on s4 | – | light |
| 9 | merge s1, s3 → s5 | yes | yes |
| 10 | `y_e` on s4 | – | light |
| 11 | disappear s5, revealing `a5` | yes | yes (no `a5` pmf given) |

At `tf` only shadow `s4` is alive (green circle). The paper notes that making
s5 disappear needs extra resources, e.g. a sub-search team.

### Table III, verbatim, in the paper's entry order

| observation | probability masses |
|---|---|
| *initial* | `P(s1=2, s2=2) = 1` |
| `y_x`, s1 | `P(s1=1, s2=2) = 0.9` |
| | `P(s1=3, s2=2) = 0.1` |
| *split*, s2 → s3, s4 | `P(s1=1, s3=0, s4=2) = 0.9*0.25 = 0.225` |
| | `P(s1=1, s3=1, s4=1) = 0.9*0.5 = 0.45` |
| | `P(s1=1, s3=2, s4=0) = 0.9*0.25 = 0.225` |
| | `P(s1=3, s3=0, s4=2) = 0.1*0.25 = 0.025` |
| | `P(s1=3, s3=1, s4=1) = 0.1*0.5 = 0.05` |
| | `P(s1=3, s3=2, s4=0) = 0.1*0.25 = 0.025` |
| `y_e`, s3 | `P(s1=1, s3=1, s4=2) = 0.225` |
| | `P(s1=1, s3=0, s4=1) = 0.45*0.1 = 0.045` |
| | `P(s1=1, s3=2, s4=1) = 0.45*0.9 = 0.405` |
| | `P(s1=1, s3=1, s4=0) = 0.225*0.1 = 0.0225` |
| | `P(s1=1, s3=3, s4=0) = 0.225*0.9 = 0.2025` |
| | `P(s1=3, s3=1, s4=2) = 0.025` |
| | `P(s1=3, s3=0, s4=1) = 0.05*0.1 = 0.005` |
| | `P(s1=3, s3=2, s4=1) = 0.05*0.9 = 0.045` |
| | `P(s1=3, s3=1, s4=0) = 0.025*0.1 = 0.0025` |
| | `P(s1=3, s3=3, s4=0) = 0.025*0.9 = 0.0225` |
| *merge*, s1, s3 → s5 | `P(s4=2, s5=2) = 0.225` |
| | `P(s4=1, s5=1) = 0.045` |
| | `P(s4=1, s5=3) = 0.405 + 0.005 = 0.41` |
| | `P(s4=0, s5=2) = 0.0225` |
| | `P(s4=0, s5=4) = 0.2025 + 0.0025 = 0.205` |
| | `P(s4=2, s5=4) = 0.025` |
| | `P(s4=1, s5=5) = 0.045` |
| | `P(s4=0, s5=6) = 0.0225` |
| *disappear*, s5 | `P(s4=0) = 0.0769 (= 0.0225*0.5/((0.0225+0.045+0.225)*0.5))` |
| | `P(s4=1) = 0.1538` |
| | `P(s4=2) = 0.7692` |

The text also says that in the merge step the 3rd and 7th entries of the
previous step are combined, and so are the 5th and 9th. Exact values:
`P(s4=0,1,2) = 1/13, 2/13, 10/13`. Before the merge there are **10** entries
(ICRA 2010: "the original sequence has only 10 entries").
Fig. 13(a)–(f) plots the same six steps: balls sized by mass, light balls
before an event and dark balls after it.

Monte Carlo, from the paper, with 1000 successful trials:
`P(s4=0) = 0.079, P(s4=1) = 0.154, P(s4=2) = 0.767`.

Basic-truncation remark (T-RO VI-E). If the fourth entry after the merge,
`P(s4=0, s5=2) = 0.0225`, were truncated, the final `P(s4=0) = 0.0769` would be
lost. Concretely, keep only the 6 largest entries of the merge output and do
not truncate at any other step. This drops both `0.0225` entries, `(0,2)` and
`(0,6)`, and the result becomes `P(s4=0,1,2) = 0, 1/6, 5/6`. Keeping 7 entries
is a bad test, because two entries tie at `0.0225`. Applying TR-6 at *every*
step gives yet another result, because the 10-entry `y_e` step is truncated
too.

## 6. The "slightly more complicated" sequence (Fig. 12 light events, Fig. 14)

Paper: five targets each in `s1` and `s2`, with the same model as above
(binomial 0.5, `p = 0.9`, no null). "135 joint probability table entries are
obtained before the merge step" (T-RO and ICRA 2010). Fig. 14 caption: the axes
`s1, s3, s4` have ranges `[1, 9]`, `[0, 7]`, `[0, 5]`.

Our exact recomputation, in `tor_fig12_complex.json`:
- In the strict figure order (events 1–8, then the merge), there are **110**
  entries before the merge. The ranges are s1 `[1,9]`, s3 `[0,7]`, s4 `[0,7]`.
- If **both** `y_e` observations on s4 (events 8 and 10) are processed before
  the merge, there are exactly **135** entries. The ranges are s1 `[1,9]`,
  s3 `[0,7]`, s4 `[0,8]`. The s4 observation commutes with the merge, so the
  final joint is the same either way. This reproduces the paper's count, and
  it strongly confirms the model: no null events, binomial split, and summing
  of identical keys. We also searched other readings by brute force: initial
  counts 2–7, up to two observations each on s2, s3 and s4, and branches of
  ±1 or ±1/0. Nothing else gives `[1,9]` × `[0,7]` together with 135 entries.
- No reading gives the caption's s4 range `[0, 5]`. Entries with `s4 ≥ 6` carry
  very little mass and are probably not visible in the plot. **Do not test the
  s4 range.**
- After the merge and event 10 (`P(s4, s5)`, 63 entries): `E[s4] = 2973/800 = 3.71625`,
  `E[s5] = 25317/3200 = 7.9115625`. The fixture also gives the full marginals.
  The paper gives no `a5` pmf for this variant, so the fixture omits the final
  disappear.

## 7. Efficient propagation: Monte Carlo and truncation (T-RO VI-E, VII-B)

### Monte Carlo baseline (rejection)
Each trial:
1. sample the initial counts from the initial joint pmf;
2. propagate the targets through the observation queue in order:
   - split: route each target independently (for the binomial rule; in general,
     sample `(n_a, n_b)` from `rule(n)`);
   - merge: add the counts;
   - appear: sample `n_e ~ P(s_e)`;
   - FOV observation `y` on `s`: sample an event, then apply it. The event must
     be sampled from `P(e | y)` **restricted to the events feasible at the
     current count**, i.e. exit is excluded when the count is 0. This mirrors
     rule (A);
   - disappear: sample `n_v ~ P(s_v)`. If it differs from the simulated count,
     the outcome "contradicts an observation" and the trial is **discarded**.
     This is the same as accepting with probability `P(s_v = count)`;
3. continue until **1000 successful trials**, the number used for every Monte
   Carlo run in the paper. Report the empirical distribution and means of the
   final shadows (Tables IV/V report the means).

Verification on Table III:
- With feasible-event sampling (A), 10^6 trials give 0.0769 / 0.1543 / 0.7688,
  which converges to 1/13, 2/13, 10/13. Five independent 1000-trial runs gave
  `P(s4=0)` between 0.056 and 0.092, consistent with the paper's 0.079.
- With "reject when an exit hits an empty shadow" (B), the result converges to
  1/12, 1/6, 3/4.
- The paper's 0.079 / 0.154 / 0.767 fit (A) much better: their z-scores against
  (A) are about 0.3 / 0.0 / −0.2, against (B) about −0.5 / −1.1 / +1.2. The
  paper also says the MC result "matches closely" the exact one. **So
  implement (A).**

The paper says this Monte Carlo "does not depend on data", so its result is
probabilistically correct and serves as the baseline.

### Basic truncation, TR-X
"Retain the first X probability mass entries of largest value." This is done
after **each update**, i.e. after every processed observation, component or
FOV. If there are more than X entries, keep the X with the largest `p_j`. Each
step then costs constant time, and the total time is linear in the number of
observations. (In ICRA 2010 this same heuristic is called **RT-X**. In T-RO,
RT means *random* truncation, see below.)

### Random truncation, RT-X
"Randomly allow probability mass entries with low value to survive
truncation." Implementation choice (T-RO VII-B): "the entries are kept based on
their probability multiplied by a random number in (0, 1)". Draw an
independent `u_j ~ U(0,1)` for each entry, rank the entries by `p_j · u_j`, and
keep the top X.

### Random truncation with event lookahead, RT-LA-X
RT-X, but more entries are retained right before merge and disappear events,
since those events shrink the table. Implementation choice (T-RO VII-B): **do
not truncate** if a **disappear** event is within the **next four events**, or a
**merge** event is within the **next two events**. Otherwise truncate as in
RT-X.

Precise rule for the implementation. After processing the observation at queue
index `t`, skip truncation if any of `Q[t+1..t+4]` is a disappear or any of
`Q[t+1..t+2]` is a merge. The paper does not say whether FOV observations count
as "events" here. We count positions in the full observation queue (component
and FOV), because Algorithm 1 processes one queue Q. Make this a parameter,
e.g. `lookahead=(4, 2)` and `lookahead_counts_fov=True`.

Notes for all three heuristics:
- Every exact operation is linear in the table, and the ranking is
  scale-invariant. Renormalising after truncation therefore does not change
  the final normalised result. Renormalise anyway, so that intermediate
  expectations make sense.
- **failure**: a disappear leaves zero total mass, because every entry
  consistent with the revealed count was truncated, so no valid result exists.
  **frequent failure** means more than one third of the runs fail. These rows
  mark roughly the smallest X at which a method works.
- For the randomised methods (RT, RT-LA), the paper averages the time over 10
  runs and reports accuracy as mean (std) over those runs.

### Complexity (T-RO VII-B)
With `n` split events, `p` targets on average in a splitting shadow, `n_f` FOV
events, and merge/appear/disappear events of the same order as splits, the
exact algorithm takes `O(p^n 3^{n_f})` time. The resampling methods have a
large constant that depends on X, but are otherwise linear in the number of
critical events.

## 8. Fig. 16 sequence and Tables IV, V (T-RO VII-B; ICRA 2010 Fig. 5, Table V)

Fig. 16 (identical to ICRA 2010 Fig. 5, which is vector art and was used to
confirm the reading). There are 20 shadows and 14 component events. The
**32 FOV observations are "scattered along the sequence" and not marked**, so
the paper does not give them. Tables IV/V therefore cannot be reproduced
exactly. The split rule is not stated either; binomial(0.5) is our assumption.

The order below is top to bottom in the figure. Shadows 1, 2, 3 are present at
`t0`, 11 appears, 13, 16 and 20 disappear, and 17, 18 and 19 are final.

| # | event |
|---|---|
| 1 | split 2 → 4, 5 |
| 2 | split 1 → 6, 7 |
| 3 | split 7 → 8, 9 |
| 4 | merge 5, 3 → 10 |
| 5 | appear 11 |
| 6 | merge 9, 4 → 12 |
| 7 | split 8 → 13, 14 |
| 8 | split 10 → 15, 16 |
| 9 | merge 11, 6 → 17 |
| 10 | disappear 13 (`a13`) |
| 11 | merge 14, 15 → 18. The horizontal line from 14 to 15 crosses 12's vertical with a bridge, so 12 is **not** part of it. |
| 12 | split 12 → 19, 20 |
| 13 | disappear 16 (`a16`) |
| 14 | disappear 20 (`a20`) |

That is 6 splits, 4 merges, 1 appear and 3 disappears, which makes 14. Events
3 and 4, and events 13 and 14, are close vertically, but they act on disjoint
shadows and commute for the exact filter.

Counts, all with probability 1:

| setting | s1 | s2 | s3 | s11 (appear) | s13 (revealed) | s16 (revealed) | s20 (revealed) |
|---|---|---|---|---|---|---|---|
| Table IV (T-RO first test) | 10 | 7 | 8 | 9 | 6 | 9 | 4 |
| Table V (T-RO second test; ICRA 2010 Table V) | 25 | 22 | 23 | 9 | 8 | 15 | 9 |

T-RO says that for the second test "we change the number of targets in shadows
1, 2, 3, 13, 16, and 20 to 25, 22, 23, 8, 15, and 9", so s11 stays at 9. ICRA
2010 lists "1, 2, 3, 11, 13, 16, 20 … 25, 22, 23, 9, 8, 15, 9", which agrees.

### Table IV: results of the probabilistic methods (first test; E[·] at `tf`)

| heuristic | s17 | s18 | s19 | t(s) |
|---|---|---|---|---|
| none, precise | 11.12 | 5.84 | 5.60 | 329.5 |
| TR-10000 | failure | | | |
| TR-20000 | 10.67 | 4.85 | 5.33 | 6.1 |
| TR-50000 | 11.00 | 5.38 | 5.52 | 15.5 |
| TR-100000 | 10.96 | 5.59 | 5.53 | 29.4 |
| TR-200000 | 11.06 | 5.73 | 5.56 | 61.4 |
| RT-10000 | frequent failure | | | |
| RT-20000 | 11.38(0.20) | 5.31(0.23) | 5.67(0.20) | 6.1 |
| RT-50000 | 11.16(0.02) | 5.36(0.03) | 5.64(0.02) | 14.6 |
| RT-100000 | 11.03(0.01) | 5.62(0.01) | 5.56(0.01) | 28.3 |
| RT-LA-2000 | frequent failure | | | |
| RT-LA-5000 | 11.32(1.30) | 6.54(1.46) | 5.28(1.28) | 2.0 |
| RT-LA-10000 | 11.17(0.56) | 5.87(0.91) | 5.02(0.48) | 4.2 |
| RT-LA-20000 | 11.62(0.26) | 5.30(0.18) | 5.57(0.16) | 8.3 |
| RT-LA-50000 | 11.32(0.01) | 5.51(0.01) | 5.60(0.01) | 17.7 |
| Monte Carlo | 11.16 | 5.58 | 5.57 | 42.1 |

### Table V: results of selected methods on a larger instance (second test)

| heuristic | s17 | s18 | s19 | t(s) |
|---|---|---|---|---|
| none, exact | out of memory after 10 mins | | | |
| Monte Carlo | 18.26 | 22.25 | 11.85 | 43.2 |
| RT-100000 | 17.78(0.03) | 21.83(0.03) | 12.34(0.02) | 66.3 |
| RT-LA-50000 | 18.24(0.08) | 21.99(0.07) | 12.57(0.07) | 40.6 |

ICRA 2010 Table V has the same Monte Carlo and RT-100000 means, without the
std values or the RT-LA row. There, "RT-100000" means *basic* truncation. In
the failed exact run, "more than 2 × 10^7 probability mass entries" were in
memory at the peak; the run had not finished the third split (T-RO). Setup:
Java 1.6, 3.0 GHz Core 2 Quad, 1.5 GB JVM, single-threaded.

### Our derived values (no FOV events), for regression only
These are exact results, binomial(0.5) splits, with the unknown FOV events
omitted. They are in `tor_fig16_sequence.json` under `derived_no_fov_exact`.
- first test: `E[s17, s18, s19] = 10.773213, 2.970124, 1.256662` (peak 28 600 entries);
- second test: `18.560103, 20.489988, 7.949909` (peak 2 637 180 entries, about 3 s in pure Python).

Without FOV events the totals are fixed. In the first test 10+7+8+9−6−9−4 = 15
targets remain at the end; the paper's means sum to 22.56, so its FOV events
add about 7.5 targets net. In the second test 47 remain, against the paper's
52.36. Use these values only to sanity-check the truncation heuristics against
the exact filter, e.g. that RT-LA-X is close to exact for moderate X, not to
compare with the paper's numbers.

## 9. Further material from ICRA 2010 (not in T-RO)

- FOV observations are "attached to the components as they evolve". A null
  observation stands for "nothing happens at the boundary for a period of
  time". The paper suggests recording how long a null event lasts and adding
  a time parameter to the sensor model; it does not elaborate.
- With beam sensors, the FOV is a line segment, so one crossing produces **two
  consecutive FOV events**: the target enters the FOV, then leaves it into the
  neighbouring shadow (ICRA 2010 Fig. 2).
- An exact count of the targets inside the FOV would improve the observation
  statistics `P(e|y)`. It is treated only as a factor in that model.
- Regimes: (1) a few agents and few events: exact; (2) a few agents and many
  events: the probabilistic answer may drift, so the nondeterministic filter
  may be better; (3) many agents: approximations are fine. The Java exact
  implementation handled "tens of agents and events combined". The heuristics
  handled "up to a thousand agents and events".
- ICRA 2010 Fig. 3 is a further "typical event observation sequence" with no
  numbers. Top to bottom: split 2 → 3, 4; `y_x` on s1; split 1 → 5, 6; `y_e`
  on s3; split 6 → 7, 8; `y_x` on s4; `y_e` on s5; merge 8, 3 → 10; `y_x` on
  s4; appear 9; merge 9, 5 → 11; disappear 10 (`a10`). The final shadows are
  11, 7 and 4. It can be used as a smoke test, e.g. exact against Monte Carlo
  for any chosen numbers.
- ICRA 2010 describes the conditional as "the conditional probability of an
  observation given an event, P(e = e | y = y)". Both papers *use* it as the
  weight of the event given the observation, which is how the fixtures encode
  it (`obs_model[y][e]`).

## 10. Extensions (T-RO VI-F), not needed for the regression tests

- Imperfect component events: keep a distribution over all shadow sequences
  consistent with the observations, and take expectations over them.
  Resampling helps.
- Teams with a single attribute: one computation per team. Joint multi-attribute
  initial conditions: one computation per joint initial entry, with resampling
  if needed.

## 11. Checklist for the implementer

1. FOV observations use **per-entry** renormalisation at `x = 0` (rule A).
   Never renormalise globally after an FOV observation.
2. Sum identical keys after every operation. Check: 10 entries before the merge
   in Table III, and 135 in the complex sequence with the s4 observations moved
   before the merge.
3. Disappear: weight by `P(s_v = x)`, drop zero entries, renormalise. Zero
   total mass means inconsistent or a failure.
4. Table III final result: exactly `1/13, 2/13, 10/13`. A Monte Carlo test
   using feasible-event sampling should land within statistical tolerance.
5. Truncate after every processed observation. TR ranks by `p`; RT ranks by
   `p·U(0,1)`; RT-LA skips truncation when a disappear is within the next 4
   queue items or a merge within the next 2.
