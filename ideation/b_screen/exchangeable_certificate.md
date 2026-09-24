# Exact whole-class certificates over the finite-exchangeable class E at n=4, 16-case design (Fable)

Per Astra's contract (`synthesis/exchangeable_certificate_contract.md`). Files: `exchangeable_certificate.py`,
`exchangeable_certificate.json` (all proof data: root-box coefficient bounds, leaves, edge enclosures with attained
witnesses and outer bounds, validation points, source hashes, timings), `exchangeable_certificate.log`. Run from the
repo root: `uv run --python 3.12 --with gmpy2 --with numpy python ideation/b_screen/exchangeable_certificate.py`
(8 processes, 199 s wall on the Mac mini; 62-115 s per case).

## Method

E(p,q) in (x,y) = (tau3, tau4) is the polygon pi_k(x,y) >= 0, k = 0..4 (Astra's count-probability formulas). Vertices
enumerated exactly by pairwise line intersection and feasibility: q=3/100: (0,0), (3/200,0), (3/100,3/100);
q=3/50: (2/75,0), (3/100,0), (3/50,3/50), (1/25,1/25). The aggregated kernels K = base + x A3 + y A4 apply on all of E
(`exchangeable_kernel_audit.md`). Per case:

1. Tensor Bernstein of Delta on E's own bounding box (degree (50,50)); derivative in y; sign certificate by branch and
   bound with polygon pruning (a box is discarded only if some affine pi_k is strictly negative on it). The scalar-lens
   pruning of `n4_vertical.certify_dw` is not used.
2. With dDelta/dy > 0 on E, a_E is the minimum of Delta on the lower vertex chain and b_E the maximum on the upper
   chain. Each edge is an exact degree-100 Bernstein polynomial (affine kernel along the edge, middle control point the
   average of the endpoints), validated against direct exact propagation at u = 1/4 and 3/4 on every edge (assertion),
   and enclosed to 1e-9 with attained witnesses. Monotonicity along edges is recorded, not relied upon.
3. Classification: identified C if the outer b_E < 0; identified R if the outer a_E > 0; unidentified if attained values
   have opposite signs. Regret min(-a_E, b_E) from attained extrema; lottery regret -a_E b_E/(b_E - a_E).
4. An informational direct sign branch-and-bound of Delta itself over E (not used for classification).

## Result

Vertical certificate on E: uniformly positive derivative coefficients on the root box in 16/16 cases (no subdivision
needed; minimum coefficient recorded per case). All 56 edge enclosures met tolerance; all 112 edge validations passed.

| q | e0 | e1 | c | S class range (certified) | S | E class range (outer) | E | E regret |
|---|---|---|---|---|---|---|---|---|
| 3/100 | 1/50 | 9/10 | 1/10 | [-0.12052, +0.03209] | unidentified (0.03209) | [-0.17925, +0.06434] | unidentified | 0.06434 |
| 3/50 | 1/50 | 9/10 | 1/10 | [+0.03101, +0.11044] | R | [-0.05342, +0.11246] | **unidentified** | 0.05342 |
| 3/100 | 1/50 | 9/10 | 1/50 | [+0.00382, +0.01035] | R | [+0.00079, +0.01250] | R | 0 |
| 3/50 | 1/50 | 9/10 | 1/50 | [+0.02685, +0.03498] | R | [+0.02268, +0.03553] | R | 0 |
| other 12 cases | | | | C | C | (ranges widen 1-40% at each end) | C | 0 |

Full table with exact rationals in the JSON and the log. Lottery regrets at the two unidentified cases: 0.04735 and
0.03622 (attained extrema).

Every E extremum is attained at a polygon vertex, and the four extreme values coincide exactly with the labelled-class
probe witnesses (`labelled_class_probe.md`): -0.1792518955, +0.0643391587, -0.0534233678, +0.1124620892. Those
witnesses are therefore the E-class extrema, not merely lower bounds.

## What this establishes

- Within E (all exchangeable four-patch hit laws with the given pairwise tables), complete pairwise information
  identifies the ranking in 14 of the 16 design cases and fails in 2; within S it identified 15 and failed in 1. The
  scalar conditional-iid representation supplies identifying information in exactly one design case
  (q=3/50, e0=1/50, e1=9/10, c=1/10), where it certified R but E admits both signs with regret 0.053.
- At the reference case the E-class regret is 0.06434, twice the S-class value; at q=3/100, e1=9/10, c=1/50 the R
  margin shrinks from 0.00382 (S) to 0.00079 (E) but survives.
- Where S identifies, E widens the range but changes no sign except in the one case above. This is a certified
  statement for these 16 cases, not for other parameters, other n, or the labelled class L (an E-identified case is
  unresolved in L; an E-reversal is a reversal in L).

## Costs

Tensor construction 54-100 s per case (the dominant cost, as Astra's pilot audit predicted); derivative + sign
< 0.1 s (uniform at root); edges 6-16 s; direct sign B&B 0-26 s. Sixteen cases in 199 s wall on 8 processes.
