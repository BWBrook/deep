# Hazard-family restriction audit (D2)

Astra, 2026-09-21. This tests whether the existing reversal requires arbitrary
mixing distributions, boundary atoms, or uncertainty between named families.
It does not calibrate a hazard generator or assess ecological plausibility.

The setting is the existing whole-class reversal point: p = 1/10, q = 3/100,
e0 = 1/50, e1 = 9/10, c = 1/10, T = 50 model steps, n = 3 or 4, with the
unchanged prescribed R and C actions. Delta = P_R - P_C. Interior witnesses
were selected constructively after inspecting an initial numerical pair;
these are existence checks, not a prespecified sample of hazard families.

## A single two-regime family retains the ambiguity

Let theta take two values a < b with probability w on b. Fixing p and q gives

```
b = (q - p*a)/(p - a)
w = (p - a)/(b - a)
0 <= a <= (p - q)/(1 - p).
```

For interior moment pairs p² < q < p, this is a one-dimensional family of
genuine two-regime laws after matching the two moments. The endpoints are
exactly Glo = {0,q/p} and Ghi = {(p-q)/(1-p),1}. This family is an identifying
restriction only if further information fixes a, b, or w.

Writing s = p-a and v = q-p² gives b = p+v/s and

`m3 = p³ + 3*p*v + v*(v/s - s)`.

Thus dm3/da = v*(1+v/s²) > 0. Its range is the complete feasible m3 interval
`[q²/p, q-(p-q)²/(1-p)]`. At n = 4 these two-point laws run along the lower
moment boundary

`m4 = (m3² - 2*p*q*m3 + q³)/(q-p²)`.

They do not fill the two-dimensional (m3,m4) body. Nevertheless, at the
existing reversal setting both **global action-difference extrema** are
already certified to be attained at Glo and Ghi. Since both laws lie in this
restricted family, its extrema and deterministic minimax regret equal those
whole-class values at this setting: approximately 0.01720186 at n = 3 and
0.03209273 at n = 4. This conclusion uses the prior certificates; it is not
a general theorem that two-point laws optimise the finite-horizon objective.

## Interior supports and continuous laws

Boundary hazard probabilities are not necessary. Two exact interior witnesses
are (a,b,w) = (1/100,29/90,81/281) and (3/40,9/10,1/33).

We also replace each atom by a uniform interval with half-width h = 1/1000.
Let A take values a,b and let U be independent Uniform[-h,h]. For theta=A+U,
set `q_centres=q-h²/3`, use q_centres in the two-point formulas above, and
retain E(A)=p. Then E(theta)=p and E(theta²)=q exactly. Both intervals lie
strictly inside [0,1]. Higher moments are rational:

`E[(mu+U)^k] = sum over even j of binom(k,j)*mu^(k-j)*h^j/(j+1)`.

These distributions have no point masses. They have piecewise constant
densities on two disjoint intervals, not globally smooth or unimodal densities.
Both witnesses belong to the same two-uniform family with fixed h; the
remaining freedom is the lower regime centre a.

| Law | m3 | m4 | Delta, n=3 | Delta, n=4 |
|---|---:|---:|---:|---:|
| Glo | 0.00900000 | 0.00270000 | -0.02341633 | -0.12052031 |
| Interior two-regime, low | 0.00964444 | 0.00310743 | -0.02157174 | -0.11543327 |
| Interior two-regime, high | 0.02250000 | 0.01991250 | +0.01228764 | +0.01025984 |
| Two-uniform, low | 0.00964433 | 0.00310739 | -0.02157208 | -0.11543403 |
| Two-uniform, high | 0.02249948 | 0.01991170 | +0.01228637 | +0.01025518 |
| Ghi | 0.02455556 | 0.02413210 | +0.01720186 | +0.03209273 |
| Beta | 0.01281818 | 0.00660629 | -0.01269735 | -0.08336312 |
| Probit-normal (numerical) | 0.01346590 | 0.00742750 | -0.01092859 | -0.07657427 |
| Logit-normal (numerical) | 0.01410184 | 0.00818603 | -0.00920586 | -0.07027327 |

The interior two-regime pair has deterministic minimax regrets 0.01228764
(n=3), 0.01025984 (n=4). The continuous pair has 0.01228637 and 0.01025518.
These are selected-pair values and lower bounds for the corresponding
restricted-class regrets, not sharp bounds for those interior classes.

## The simpler single-component comparators agree

The moment-matched beta law has alpha=7/20, beta=63/20. Its moments, kernels,
persistence values and rankings are exact rational calculations. For
probit-normal theta=Phi((Phi^-1(p)+sqrt(rho)*Z)/sqrt(1-rho)), deterministic
quadrature gives rho=0.4598486442003875. For logit-normal
theta=logistic(mu+sigma*Z), the fitted parameters are mu=-3.102015655494 and
sigma=1.6678692665272061. In both cases Z is standard normal.

All three comparators favour C at both system sizes. Their union therefore
does not reverse the ranking at this setting, despite nonzero moment spread.
The logit-normal fit converged to the same solution from three starting
points; this is not a proof of global parameter uniqueness. Probit/logit
results are numerical fitted examples, not certified ranges over all roots.
We do not interpret this small comparison as evidence that the regime laws
are more credible, or that the single-component families are adequate.

## Evidence and limits

- [Code](hazard_family_audit.py) and [results](hazard_family_audit.json) contain
  exact rational moments and persistence for atomic, uniform-mixture and beta
  laws, numerical results for probit/logit laws, parameters, and source hashes.
- Exact transition kernels use the existing degree-n moment representation;
  for every atomic law they agree entry-for-entry with direct atom averaging.
  All exact rows are nonnegative and sum to one. Persistence is the direct
  surviving mass after 50 steps, with rational arithmetic for exact cases.
- Numerical moment fitting uses deterministic adaptive normal integration;
  1024-node Gauss-Hermite integration agrees within 1e-15 on the moments.
  Adaptive quadrature error estimates are about 2e-13; these are numerical
  diagnostics, not rigorous error enclosures.
- An independent audit integrated the separate full-labelled-state kernel
  (`b_screen_v2.kernel_labelled`) with 64-point Gauss-Legendre quadrature on
  each uniform interval, then propagated transient mass. All four uniform-law
  Delta values agreed within 1.30e-14; p and q agreed within 2e-18.
  Its [standalone check](check_hazard_uniform_labelled.py) can be rerun with
  `uv run --python 3.12 --with numpy --with scipy python ideation/b_screen/check_hazard_uniform_labelled.py`.
- The exact machinery shares n3_certificate.py and n4_vertical.py with earlier
  certificates. The separate labelled implementation checks the moment-based
  calculation numerically; it is not a fully independent exact proof system.

The unrestricted-class objection now has a more precise target: why exclude
even uncertain two-regime hazards, or why prefer one fitted single-component
family? That requires ecological evidence or substantive assumptions. The
mathematics rules out boundary atoms and arbitrary distributional complexity
as necessary causes, while leaving mechanism credibility and field decision
value unresolved. D2 is not a blanket pass for every interpretation of the
claim.
