# Failure-region pilot: eight exact qualification points (Fable)

Per Astra's contract (`synthesis/failure_region_contract.md`). Files: `failure_region_pilot.py`,
`failure_region_pilot_points.json` (points saved before evaluation), `failure_region_pilot.json`, `failure_region_pilot.log`.
Run from the repo root: `uv run --python 3.12 --with gmpy2 --with numpy python ideation/b_screen/failure_region_pilot.py`
(about 12 minutes on the Mac mini). Setting: n=4, p=1/10, q=3/100, e0=1/50, T=50, scalar class S.

## What was done at each point

(A) exact Delta at Glo, Ghi and the lens midpoint; (B) a direct two-dimensional sign certificate of Delta over the
lens (tensor Bernstein of degree (50,50), branch and bound with lens pruning, no monotonicity assumed); (C) the
vertical certificate dDelta/dm4 > 0; (D) degree-100 boundary-curve polynomials with exact min/max enclosures
(tolerance 1e-9), sign pieces, and a strict-monotonicity test (all derivative Bernstein coefficients of one sign);
(E) the shortcut classification from (C)+(D), compared with (B).

## Results (exact; decimals rounded)

| point | e1 | c | Delta(Glo) | Delta(Ghi) | direct 2-D certificate | vertical | both curves increasing | roots (lower, upper) | classification | class Delta range |
|---|---|---|---:|---:|---|---|---|---|---|---|
| P1 | 1/2 | 9/50 | -0.03586 | -0.00558 | identified C | pos | yes | 0, 0 | C | [-0.03586, -0.00558] |
| P2 | 1/2 | 11/50 | -0.02404 | +0.00529 | unidentified | pos | yes | 1, 1 | unidentified, regret 0.00529 | [-0.02404, +0.00529] |
| P3 | 7/10 | 1/10 | -0.12277 | -0.03092 | identified C | pos | yes | 0, 0 | C | [-0.12277, -0.03092] |
| P4 | 7/10 | 3/20 | -0.08510 | +0.01794 | unidentified | pos | yes | 1, 1 | unidentified, regret 0.01794 | [-0.08510, +0.01794] |
| P5 | 9/10 | 3/50 | -0.09418 | -0.01373 | identified C | pos | yes | 0, 0 | C | [-0.09418, -0.01373] |
| P6 | 17/20 | 1/50 | -0.00117 | +0.00591 | unidentified | pos | yes | 1, 1 | unidentified, regret 0.00117 | [-0.00117, +0.00591] |
| P7 | 1/2 | 2/5 | -0.00382 | +0.02427 | unidentified | pos | yes | 1, 1 | unidentified, regret 0.00382 | [-0.00382, +0.02427] |
| P8 | 19/20 | 9/20 | +0.00659 | +0.21652 | identified R | pos | yes | 0, 0 | R | [+0.00659, +0.21652] |

"Classification" is certified: identified cases by the direct 2-D sign certificate (and, independently, by the
shortcut); unidentified cases by exact opposite-sign attained witnesses. "Class Delta range" is the exact
whole-class range, valid because at every point the vertical certificate held and both boundary curves were strictly
increasing, so the class extrema are the two lens endpoints, i.e. Delta(Glo) and Delta(Ghi). Regret for unidentified
points is min(-a, b) of that range.

## Where the shortcuts held

- Vertical certificate: held at 8/8 points, each at branch-and-bound depth 0 (uniformly positive coefficients), cost
  under 0.1 s. It is not uniform over the domain by proof, only at these points.
- Strict monotonicity of Delta along both boundary curves: held at 8/8 points. Consequently the class extrema were the
  endpoint-law values at every point, and "both roots exist" agreed with "unidentified" at every point. Astra's
  logical counterexample to the root shortcut stands; it simply did not arise here.
- The direct 2-D certificate agreed with the shortcut at 8/8 points. It is the expensive stage (68-104 s per point);
  witnesses + vertical + curves cost 5-11 s per point.

## What the points say about the region

- The endpoint-pair scan's diagonal boundary is a whole-class boundary at the tested resolution: at e1 = 1/2 the
  class switches from identified C (c = 0.18) to unidentified (c = 0.22); at e1 = 7/10 between c = 0.10 and 0.15;
  at e1 = 9/10, c = 0.06 is still identified C (reference c = 0.10 unidentified).
- P6 (e1 = 0.85, c = 0.02) is a genuine second unidentified component near the low-c edge, not a cancellation
  artefact: Delta(Glo) = -0.00117 exactly negative. Its regret is small (0.00117).
- P8 (e1 = 0.95, c = 0.45) is identified R over the whole class with a large range: at strong colonisation and
  severe shocks, buffering wins throughout S.
- Failure is not confined to a corner in this slice: unidentified points occur at (0.5, 0.22), (0.7, 0.15), (0.85, 0.02),
  (0.5, 0.40) and the reference (0.9, 0.10). Regret within the unidentified set ranged 0.001-0.032 at the tested
  points; the reference is the largest so far.

## Method and resource recommendation

1. For a trace at declared resolution: per point, run (A) + (C) + (D) only (about 10 s). Accept the endpoint
   classification when both certificates hold; fall back to (B) (about 90 s) when either fails. On this evidence the
   fallback rate is 0/8, so a 30 x 30 point trace of the contract domain costs about 2.5 hours on this machine, less
   in parallel on deep-server. Point classifications do not certify cells (contract section 4); the plotted boundary
   would be a bracket at grid resolution.
2. For cell certificates: the rectangle machinery (`n4_rectangle_certificate.py`) bounds the endpoint pair over
   (e1, c) cells; upgrading a cell to a whole-class certificate would additionally need cell-wise vertical and
   curve-monotonicity certificates, which do not yet exist. Recommend cell certificates only along the located
   boundary, not over the domain.
3. Before committing to the trace: the pilot already shows the region is large (roughly the half-plane above the
   diagonal, plus a thin low-c band), so the scientific value of a full trace is the boundary's shape and the regret
   surface within the region, not the fact of failure.

## Not shown

Anything at other q, e0, T, n; certified cells; the E or L classes (the S kernels apply to E with tau3, tau4 in place of
m3, m4, so the same machinery could certify E-class extrema over its polygon, but that was not run here).
