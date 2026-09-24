# Working abstract

Astra revision after Fable #499, 2026-09-24. Internal draft; no venue or
submission decision. The scientific claims are conditional on the model and
information classes specified below.

Two environmental hit laws can agree at every patch and every pair yet differ
in how often three or more patches are hit together, enough to reverse
conservation action rankings. We separate information supplied by pairwise
observations from information imposed by dependence assumptions. In a
stochastic metapopulation model, we compare buffering one patch against shock
mortality with doubling colonisation. For conditionally independent hits given
a shared scalar hazard, finite-horizon persistence depends on the first n
hazard moments; pairwise hit laws fix only two. Exact bounds over compatible
hazard laws give deterministic minimax regrets of 1.72 and 3.21 percentage
points of persistence at one three- and four-patch setting over 50 model steps.
Feasible-law pairs establish regret lower bounds above 8.06 and 9.81 points
at five and six patches; once-only randomisation between actions leaves bounds
above five points. Symmetric four-patch hit laws outside the scalar-hazard
class retain every pairwise table yet reverse a ranking identified within
that class. At the reference four-patch setting, they raise the regret lower
bound above 6.43 points. Conversely, restricting the scalar hazard law to
Khinchin unimodality restores identification at that setting. Results assume
fixed action packages, independent shocks across time and uncalibrated model
steps; they establish neither field cost-effectiveness nor the ecological
prevalence of reversals.

## Agreed scope corrections

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
