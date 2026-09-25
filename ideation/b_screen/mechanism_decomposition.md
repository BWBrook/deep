# Mechanism diagnostic on the four E-class extremisers (Fable)

**Interpretation superseded by `mechanism_audit.md` (Astra, 2026-09-25).** The arithmetic below stands (independently
reconstructed). The interpretive section overstated: an all-four hit is exposure, not death (survival from full occupancy
is .3439 under C and .6647 under R); colonisation can refill vacancies in the same step; quiet-step colonisation shares
are 46-72%, not "almost entirely"; per-category nets are order-dependent; q differs between the two pairs and the q=.06
minimiser has tau3 = 2/75 > 0; cumulative sums are accounting prefixes, not shorter-horizon Delta. The sentence "buffering
insures against wipeouts, connectivity repairs partial losses" is withdrawn. The retained account is in the audit and in
`synthesis/results_identification.md`.

Per Astra #516. Files: `mechanism_decomposition.py`, `mechanism_decomposition.json`, `mechanism_decomposition.log`.
Run from the repo root: `uv run --python 3.12 --with numpy --with scipy python ideation/b_screen/mechanism_decomposition.py`
(about 40 s). n=4, e0=1/50, e1=9/10, c=1/10, T=50; exact rationals on the common labelled 16-state chain.

## Checks

- One-step identity Delta_1 = (s/2)[e0^3 p + 3 e0^2 s q + 3 e0 s^2 tau3 + s^3 tau4], s = e1 - e0: holds exactly at all
  four laws (2233/156250000 and 1503601/156250000 at q=.03 as Astra stated; 0.00057341 and 0.01924574 at q=.06).
- Telescoping totals: both orderings sum exactly to the certified Delta at all four laws (assertion).

## Decomposition (forward ordering: R past, C future; reverse in brackets)

| law (hit-count probabilities pi_0..pi_4) | Delta_50 | buffering gain | colonisation gain | net by hit count k=0..4 (forward) |
|---|---:|---:|---:|---|
| q=.03 E-min: (0.78, 0.04, 0.18, 0, 0) | -0.1793 | +0.163 [+0.153] | +0.342 [+0.332] | [-0.212, -0.007, +0.039, 0, 0] |
| q=.03 E-max: (0.69, 0.28, 0, 0, 0.03) | +0.0643 | +0.182 [+0.184] | +0.118 [+0.120] | [-0.067, -0.009, 0, 0, +0.141] |
| q=.06 E-min: (0.853, 0, 0.04, 0.107, 0) | -0.0534 | +0.184 [+0.182] | +0.238 [+0.235] | [-0.172, 0, +0.007, +0.111, 0] |
| q=.06 E-max: (0.78, 0.16, 0, 0, 0.06) | +0.1125 | +0.157 [+0.160] | +0.045 [+0.048] | [-0.031, +0.002, 0, 0, +0.142] |

"Buffering gain" is sum_t mu K_R^t (K_R - K_0) K_C^(T-1-t) h, the accumulated one-step contrast of buffering against the
baseline; "colonisation gain" the same with K_C - K_0; Delta = buffering - colonisation. Order dependence is small
(allocations move by about 0.01) and never changes the qualitative picture.

## What distinguishes the rankings

1. **Where each action earns its gain.** Buffering's gain accrues almost entirely at steps whose hit count is the
   law's largest (k=2 for the q=.03 minimiser, k=3 for the q=.06 minimiser, k=4 for both maximisers). Colonisation's
   gain accrues almost entirely at k=0 steps (0.21 of 0.34, 0.17 of 0.24, 0.07 of 0.12, 0.03 of 0.04): doubled
   colonisation refills patches emptied by earlier partial losses during quiet steps. Colonisation cannot act on a
   step's own mortality (the one-step identity) and cannot rescue an extinct region.
2. **Bulk versus dispersed co-hits.** All four laws share p and q, hence the same expected number of co-hit pairs per
   step, E[K(K-1)]/2 = 6q. The minimisers deliver those pairs dispersed over frequent partial-loss steps (k=2 or 3),
   which leave survivors and empties for colonisation to repair; the maximisers deliver them in rare all-four steps
   (pi_4 = tau4 = q: 0.03 and 0.06 here) that leave nothing to repair. With all four hit, the region
   survives the step with probability 1 - 0.9^4 = 0.344 unbuffered and 1 - 0.9^3 x 0.46 = 0.665 buffered; that single
   contrast is the k=4 buffering entry (+0.145, +0.145) and it is the whole R advantage.
3. **The occupancy trajectories confirm it.** Under the q=.03 minimiser, C keeps the region fuller at every time
   (P(4 occupied) at t=25: 0.49 under C versus 0.26 under R) and has less extinction (0.165 versus 0.264). Under the
   maximiser, C is extinct more often at every time (0.43 versus 0.36 at t=25): partial states are rare, so there is
   little to refill, while wipeouts are more frequent and only buffering mitigates them.
4. **Why pairwise information cannot decide.** q fixes the mean co-hit rate, not whether co-hits arrive in bulk. The
   scalar-hazard class S ties tau3 and tau4 to (p, q) through the Hausdorff lens, which excludes both extremes
   (tau3 = 0 and tau3 = q); E admits them. This is the mechanism behind the S-versus-E result: the conditional-iid
   representation constrains how concentrated the co-hits can be, and that concentration is what the two action
   packages respond to in opposite directions.

## Limits

Four laws at one ecological setting; the decomposition is order-dependent and neither allocation is a unique causal
mediation; "buffering insures against wipeouts, connectivity repairs partial losses" is a description of these
contrasts, not a theorem about other parameters or horizons. No new laws, grids or policies.
