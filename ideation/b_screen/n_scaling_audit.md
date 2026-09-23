# Audit of the patch-number witnesses and their decision-loss bounds

Astra, 2026-09-23. Target: Fable's `n_scaling_probe.{py,json,md,log}` at
`beae436`. This audit re-evaluates the saved rational laws; it does not rerun
the random search or certify whole-class extrema for n>=5.

## Accepted numerical result

Every saved rational witness has support in [0,1], nonnegative weights summing
to one, and exactly p=1/10 and q=3/100. Exact propagation reproduces its reported
action difference. Parameters remain e0=1/50, e1=9/10, c=1/10 and T=50 model
steps. All 24 reconstructed action kernels are nonnegative, row-stochastic,
and absorbing at zero. The n=3,4 endpoint results agree with the existing
whole-class certificates.

The following decimals are rounded **down**, so each is a safe lower bound
on the corresponding whole-class minimax regret. For n=3,4 the underlying
unrounded values equal the known class optima. For n=5..8 they are values for
the saved pair, and the full class may have larger regret.

| n | Deterministic regret lower bound | Once-only lottery regret lower bound |
|---:|---:|---:|
| 3 | 0.01720185 | 0.00991684 |
| 4 | 0.03209273 | 0.02534400 |
| 5 | 0.08067297 | 0.05097387 |
| 6 | 0.09817136 | 0.05202150 |
| 7 | 0.05550879 | 0.03798837 |
| 8 | 0.03188587 | 0.02521955 |

The exact n=5,6 pair values therefore guarantee at least 8.067 and 9.817
percentage points of deterministic regret, respectively. Even allowing a
once-only lottery over the two complete action policies leaves guaranteed
lower bounds exceeding 5.09 and 5.20 percentage points. This does not cover
richer physical actions, switching, adaptation, or cost-matched policies.

For attained differences L<0<U, the pair's deterministic regret is
min(-L,U). Its lottery value is -LU/(U-L), with probability U/(U-L) of R.
Minimising worst-case regret over a larger ambiguity class cannot reduce
either value, so both give valid lower bounds for the full class. The
pair-optimal lottery probability is not established as class-optimal for
n>=5.

## What the scaling comparison establishes

The attained deterministic bounds increase from n=3 through n=6, then
decrease in these witnesses at n=7,8. Since n=3,4 have known class optima,
the n=5 and n=6 bounds prove that the class regret at each of those patch
numbers exceeds the class regret at either n=3 or n=4.

They do **not** prove a monotonic increase of whole-class regret from n=5
to n=6, a peak at n=6, or a subsequent decline. Lower bounds cannot order
unknown optima that way. Persistence near one at Glo is selected-law evidence,
not a uniform saturation bound over all compatible G. The positive witness
still has Delta about 0.12 at n=8. No causal saturation explanation or
extension to untested parameters follows.

Holding p,q fixed is the intended information comparison, but increasing n
also changes total habitat, the fraction protected by a one-patch buffering
action, and the number of colonisation sources at fixed per-source c. Thus
the result does not isolate the number of free moments as a biological cause.
The equal-(n-1)-subset witnesses answer a different information question and
cannot be used to infer decay of this fixed-pairwise-information regret.

## Rationalisation and reporting corrections

The two-point formula b=(q-pa)/(p-a), followed by weights
(b-p)/(b-a) and (p-a)/(b-a), enforces p,q exactly when the resulting law is
feasible. All twelve saved rationalised witnesses pass the exact feasibility
checks. They are valid evidence regardless of whether the optimiser found
global extrema.

The original search accepts weights as low as -1e-12 and then drops components
with weight <=1e-7 before rationalising. The resulting rational law can differ
from the floating candidate. In the saved probe, **every floating headline
regret is slightly larger than the value certified by its rational pair**.
For example, at n=6 the headline is 0.09817136901194445 whereas the exact pair
value is 0.098171369002105... . The discrepancy does not change the substantive
result but means the floating headline is not itself a rigorous lower bound.
Likewise “regret >= 0.05551” at n=7 rounds upward; use a downward bound from
the receipt or label a displayed value approximate.

The original main function also compares the floating search winners with
Glo/Ghi for its headline bounds, while saving exact re-evaluations only of
the search winners. The present saved laws are valid, but a future atlas
must exact-check all candidates used in reporting and recompute the final
reported bounds from the feasible rational laws. Exact moment matching during
the search is not established by floating linear solves. Rounded collisions
or an infeasible rationalisation need explicit handling before a broader run.

## Reproduction and dependency limits

```sh
uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/n_scaling_audit.py
```

The receipt retains exact PR, PC, Delta, both pair regrets, outward decimal
enclosures, witness feasibility and source hashes. The audit uses independent
gmpy2 propagation while sharing the established exact kernel constructors;
it does not independently rederive all demographic transitions. The previous
whole-class optimisation certificates are reused for interpreting n=3,4.
No random-search completeness claim is made.
