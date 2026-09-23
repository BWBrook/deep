# Shape restriction: the n=4 ranking is identified within the unimodal class (Fable)

Files: `unimodal_check.py`, `unimodal_check.json`, `unimodal_check.log`. Run from the repo root:
`uv run --python 3.12 --with gmpy2 --with numpy --with sympy python ideation/b_screen/unimodal_check.py` (about 2 s).

## Question

In the one unidentified design case (p = 1/10, q = 3/100, e0 = 1/50, e1 = 9/10, c = 1/10, T = 50) the sign of
Delta = P_R - P_C is not identified over the full conditional-iid class: certified Delta range [-0.12052, +0.03209].
Does restricting the yearly hazard law G to unimodal laws identify the sign?

## Shape class (answering Astra #461)

Khinchin's representation: theta on [0, 1] is unimodal about mode m iff theta = m + U Y with U ~ Uniform(0, 1)
independent of Y, and Y supported on [-m, 1-m]. The class is the closed convex hull of uniform laws sharing the
endpoint m. It **allows**: a boundary mode (m = 0 gives all non-increasing densities; m = 1 all non-decreasing), an
atom at the mode (Y = 0 with positive mass), and arbitrary mixtures of uniforms sharing the mode. It **excludes**
atoms away from the mode and any law with two local maxima. In particular the full-class extremizers Glo = {0, 3/10}
and Ghi = {a, 1} (two atoms each) are bimodal and lie outside the class.

## Argument

Two certified inputs, both Astra-audited: (i) dDelta/dm4 > 0 over the whole (m3, m4) lens for this case
(`n4_vertical_all.json`); (ii) along the upper curve m4 = U(m3), Delta is strictly increasing in m3 with a single
root, bracketed exactly at m3 in [4177210314157/206158430208000, 50126523769919/2473901162496000] ~ 0.0202621368
(`n4_boundary_results.json`). Hence any G in the class with m3(G) <= x* := 4177210314157/206158430208000 has
Delta(G) <= Delta(m3, U(m3)) < 0. (Symmetrically, m3 above the lower-curve root ~0.0215364 forces Delta > 0; the
whole-class ambiguity is exactly the m3 band between the two roots, as already used in `n4_reversal_map`.)

So the shape question reduces to: what is sup m3 over unimodal G with (p, q) = (1/10, 3/100)?

Khinchin gives E theta^k = sum_j C(k, j) m^(k-j) E[Y^j]/(j+1), so with y_j = E Y^j:
y1 = 2(p - m), y2 = 3(q - m^2 - m y1), m3 = m^3 + (3/2) m^2 y1 + m y2 + y3/4. Given (y1, y2) on [A, B] = [-m, 1-m],
y3 is maximised by the two-point law on {a, B}, a = y1 - Var(Y)/(B - y1) (upper principal representation), feasible
iff a >= A. This gives an exact rational function of the mode:

    m3max_uni(m) = (-7500 m^2 + 1780 m + 599) / (8000 (5 m + 4)),   feasible modes m in [0, 11/80].

(No unimodal law with these (p, q) can have its mode above 0.1375.)

## Result

- sup over unimodal laws of m3 = **0.0197404736** at mode m = sqrt(3)/2 - 4/5 ~ 0.0660 (the only stationary point).
- Exact certificate: x* - m3max_uni(m) = (193273528320000 m^2 - 24984199150495 m + 1272728794804) /
  (1030792151040000 m + 824633720832000); Sturm root count of numerator and denominator on [0, 11/80] is 0 for both,
  both positive at the midpoint. So m3(G) < x* for **every** unimodal G, hence **Delta(G) < 0: C is strictly preferred
  for every unimodal hazard law.** The ranking is identified within the unimodal class.
- Rigorous margin: Delta(G) <= Delta(79/4000, U(79/4000)) = **-0.00427** for every unimodal G (79/4000 = 0.01975 lies
  between sup m3_uni and x*; Delta increases in m4 and along the upper curve). Exploratory float search over
  (mode, y3) with y4 at its Hausdorff maximum found the largest unimodal Delta at about -0.0160 (m ~ 0.066, at the
  m3-maximising law). So the true unimodal supremum of Delta lies in [-0.0160, -0.0043]; the full-class supremum is
  +0.0321.
- The m3-maximising unimodal law is 0.9272 * Uniform[0.0660, 0.0661] + 0.0728 * Uniform[0.0660, 1]: essentially an
  atom at the mode plus a flat tail to 1. Direct exact evaluation: Delta = -0.015965 (P_R 0.42796, P_C 0.44392).
- Unimodal m3 range is [0.010125, 0.019740] versus the full-class [0.009, 0.024556]; the ambiguity band
  (0.020262, 0.021536) lies entirely above it, with a gap of 0.00052 in m3.
- For reference, at q = 6/100 the feasible modes are [0, 1/40] and the unimodal m3 range is [0.0405, 0.0449] versus
  the class range [0.036, 0.0582]; no ranking claim is made there (all q = 6/100 design cases were already identified).

## Sharpness

The restriction cannot be weakened to "at most two modes": Ghi is a mixture of two point masses, each unimodal, so the
class of two-component unimodal mixtures reaches the full-class m3 maximum and the whole ambiguity band. One mode
identifies the ranking; two modes restore the full ambiguity. This is consistent with D2 (the two-regime family spans
the certified ambiguity) and sharpens it: within this setting, a second regime is necessary for the ranking to be
ambiguous, not just sufficient.

## Checks performed

- Khinchin moment identity verified exactly on 200 random rational mixtures of uniforms sharing a mode (0 failures).
- Two-point y3 bounds verified on 20000 random discrete Y laws with 1-6 atoms (0 violations).
- The extremal law's moments (p, q, m3, m4) recomputed by direct integration of the uniform pieces; equal.
- The extremal law's m4 lies inside the lens [L(m3), U(m3)].
- Delta at the extremal law evaluated exactly with the affine moment kernels from `n4_vertical.py`.

## Scope

n = 4, this one parameter setting, T = 50. The argument is reusable for any setting with certified (i) and (ii): the
unimodal m3 range depends only on (p, q), so any ambiguity band lying above 0.019740 (at p = 1/10, q = 3/100) is
identified within the unimodal class without further computation. Not claimed: anything about n = 3, 5, 6, other
(e1, c), or laws with atoms away from the mode.
