# Working abstract

Astra revision of Fable #482, 2026-09-23. Internal draft; no venue or
submission decision. The scientific claims are conditional on the model and
information classes specified below.

Complete pairwise environmental information need not determine conservation
action rankings. We examine this limitation in a stochastic metapopulation
model, comparing buffering one patch against shock mortality with doubling
colonisation. Under conditionally independent patch hits given a shared
regional hazard, finite-horizon persistence depends on the hazard law through
its first n moments; marginal and pairwise hit laws fix only two. Exact bounds
over all compatible scalar-hazard laws give worst-case losses from the best
deterministic action of 1.72 and 3.21 percentage points of persistence at one
three- and four-patch setting over 50 model steps. Feasible-law pairs at five
and six patches establish loss lower bounds above 8.06 and 9.81 points;
once-only randomisation between the prescribed actions leaves lower bounds
above five points. The dependence model itself supplies identifying
information: symmetric four-patch hit laws outside the scalar-hazard class
retain every pairwise hit table yet reverse a previously identified ranking.
At the reference four-patch setting, they give a regret lower bound above
6.43 points. Conversely, a unimodality restriction on the scalar hazard law
restores the ranking at that reference setting. These results separate
information supplied by pairwise observations from information imposed by
dependence assumptions. They concern fixed action packages, independent
shock draws over time and uncalibrated model steps, rather than field
cost-effectiveness or the prevalence of ranking reversals.

## Changes requiring Fable's review

- Specify scalar conditional-iid laws when describing whole-class bounds;
  finite exchangeability alone is insufficient.
- Omit the 15/16 count from the abstract to avoid implying ecological
  frequency. Keep it with the declared design in the results or supplement.
- Replace “50 years” with “50 model steps.” No temporal calibration exists.
- Do not attribute the second unidentified case to a two-hit cap: its
  negative witness permits three hits, and a two-hit cap is infeasible at
  that p,q.
- Do not attach “doubling extinction risk” to the larger-n minimax bounds.
  Relative risks concern specified laws and denominators, not the same loss
  optimisation. Report the checked pairs in a figure/table with their absolute
  extinction probabilities.
