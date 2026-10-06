# Fig. 11 max-flow example (T-RO 2012, Sec. V-C, pp. 448–449)

Fixture: `tests/fixtures/tor_fig11.json`.

## Reconstruction

Fig. 11 is a 400×439 JPEG embedded in the PDF at 150 dpi, so it was read at its
native resolution (`pdfimages -j`, upscaled). Panels (b), (c) and (d) agree.

Shadow sequence I-state, Fig. 11(b). Time runs downwards.

| kind | shadows |
|---|---|
| at `t0` | 1, 2, 3, 4, 5 |
| appear | 10, 12, 18 |
| split | 1 → 6, 7; 7 → 9, 8; 3 → 13, 14 |
| merge | 2 + 6 → 15; 4 + 10 → 11; 11 + 5 → 16; 8 + 16 → 17; 17 + 12 → 19 |
| disappear | 9, 14 |
| at `tf` | 13, 15, 18, 19 (18 both appears and survives) |

The events are ordered by their vertical position in the figure. The one
exception is `appear(18)`, which is moved before the merge that creates s19 so
that labels increase. Moving it does not change the bipartite graph.

Bipartite I-state, Fig. 11(c), produced by this sequence:

- left: 1, 2, 3, 4, 5, 10, 12, 18
- right: 9, 14 (disappeared) and 13, 15, 18, 19 (final)
- edges: 1–9, 1–15, 1–19, 2–15, 3–13, 3–14, 4–19, 5–19, 10–19, 12–19, 18–18

Input bounds from the text: 1:(2,4) 2:(0,3) 3:(5,5) 4:(2,6) 5:(4,5) 9:(2,3)
10:(1,3) 12:(3,8) 14:(2,4) 18:(5,7). The two disappear events keep these ranges
rather than exact counts.

## Verification (scipy `linprog`, plus networkx max-flow for the recipe)

| | lower s19 | upper s19 |
|---|---|---|
| paper | 10 | 24 |
| LP with the plain meaning (`l ≤ supply ≤ u` on every left vertex, `l ≤ revealed ≤ u` on 9 and 14) | 10 | 24 |
| paper's literal recipe (Eqs. 7, 8; `c(i,T)=l_i` for disappearing shadows) | 10 | 24 |

So for this example (a) and (b) give the same answer, and 24 and 10 are the
true optimum. By hand:

- **Upper bound.** 1 must send at least 2 to 9, which leaves 2 for 19. Then
  2 + 6 + 5 + 3 + 8 = 24.
- **Lower bound.** 4, 5, 10 and 12 can only reach 19, so 19 gets at least
  2 + 4 + 1 + 3 = 10. 1 can send everything to 9 and 15.

LP bounds for the other final shadows: 13:[1,3], 15:[0,5], 18:[5,7].

The quoted witness flows are max-flows of the recipe networks. They are not
feasible target distributions:

- **Lower-bound witness.** It has f(9,T)=f(14,T)=1 < l=2 and total 12. It is a
  different max-flow of the same value, 12. Eq. 8 uses only the total, so
  10 = 22 − 12 either way.
- **Upper-bound witness.** It has f(18,T)=7 even though the text sets
  c(18,T)=0. This is a small slip: the witness only fits if c(18,T)=∞. It does
  not change the result because 18 is not a disappearing shadow.

The fixture includes feasible witnesses for both bounds.

## Caveat for `bounds_paper`

On random bipartite instances the literal lower-bound recipe disagrees with the
LP in about 54% of feasible cases (752 of 1381). Here is why:

- It caps disappearing shadows at `l_i`, not `u_i`.
- Eq. 8 then counts any lower-bound supply that cannot be routed as if it went
  into the query.

With `c(i,T)=u_i` for disappearing shadows in the lower-bound run, the recipe
matched the LP on all tested instances. The upper-bound recipe matched the LP
on all feasible instances tested.

The flow-with-lower-bounds formulation in `DESIGN.md` §3.2 should be treated as
the reference. `bounds_paper` should only be required to match the LP on the
Fig. 11 example.
