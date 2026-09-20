# Three-patch moment and polynomial certificate

Scope: the fixed, iid regional-hazard model in correspondence thread 49,
message 355; the 16 three-patch parameter cases reported in thread 48,
message 352. This is a bounded ideation audit, not a characterization of all
ecological parameter values or a claim of novelty.

## Reduction to one moment

Write `p = E(theta)`, `q = E(theta^2)`, and `r = E(theta^3)`, with
`0 < p < 1` and `p^2 <= q <= p`. Nonnegativity of two expectations gives

```
E[theta (theta - q/p)^2] >= 0
    => r >= q^2/p = L.

E[(1-theta) (theta - (p-q)/(1-p))^2] >= 0
    => r <= q - (p-q)^2/(1-p) = U.
```

The lower endpoint is attained by the law on `{0, q/p}` with mass `p^2/q`
on its nonzero atom. The upper endpoint is attained by the law on `{a, 1}`,
where `a=(p-q)/(1-p)` and the mass at 1 is
`w=(q-p^2)/(1-2p+q)`. These formulas apply to the nondegenerate cases used
here. Mixtures `G_lambda=(1-lambda)G_lower+lambda G_upper` attain every
`r=L+lambda(U-L)`, `0 <= lambda <= 1`.

Each conditional transition entry is a polynomial of degree at most three
in theta: mortality involves at most three occupied patches, and
recolonisation conditional on the survivor state has no theta dependence.
Consequently every admissible G with a given r has the same integrated
transition matrix as the matching G_lambda. For each action,

```
K(lambda) = (1-lambda) K_lower + lambda K_upper.
P(lambda) = v_initial K(lambda)^50 s,
```

where `s` indicates nonempty states. Since extinction is absorbing,
occupancy at step 50 is equivalent to survival through step 50.
This parameterizes the entire attainable set of outputs for three patches.
It does not establish where the extrema of P or P_R-P_C occur.

## Exact Bernstein calculation

Use Bernstein basis `B[k,t](lambda) = choose(t,k) lambda^k (1-lambda)^(t-k)`.
If the state-distribution polynomial after t steps has vector coefficients
`v[t,k]`, its recurrence is

```
v[t,k] = ((t-k)/t) v[t-1,k] K_lower
       + (k/t)     v[t-1,k-1] K_upper.
```

Absent terms are zero. The initial distribution occupies all three patches.
All parameters, endpoint matrices, and coefficient operations use exact
rational arithmetic. Summing coefficients over nonempty states gives
coefficients `b[k]` for persistence; subtracting the two action coefficient
vectors gives coefficients for their difference. For degree T,

```
P'(lambda) = T sum_k (b[k+1]-b[k]) B[k,T-1](lambda).
```

The Bernstein basis is nonnegative and sums to one on [0,1]. Uniformly
positive or negative derivative coefficients therefore certify monotonicity
and place extrema at the endpoints. Uniformly signed objective coefficients
certify the sign of the objective. Exact subdivision can tighten these
sufficient tests if coefficient signs are mixed. A finite sweep alone is
not such a certificate, even when the polynomial degree is known.

If the action difference is strictly monotone and has opposite endpoint
signs, there is exactly one root; exact rational bisection gives a guaranteed
bracket. Decimal summaries are rounded displays of rational results.

## Reproduction and interpretation

Run `python3 ideation/b_screen/n3_certificate.py` from the repository root.
The companion JSON records the finite-case certificates and any unresolved
tests. A sign-identified ranking does not imply point-identified persistence
or a point-identified action difference. Equality at a threshold is a tie.
These results assume a single fixed G across years, independent yearly draws,
the stated homogeneous model, and the predeclared action packages.

For three patches the remaining moment has a direct interpretation:
`r = P(all three patches hit in a year)`. More generally, if H is the number
of patches hit among n, conditional independence gives
`E[choose(H,k)] / choose(n,k) = E(theta^k)`. This identifies a relevant
environmental-hit statistic; it does not establish its estimation precision
from short records, or make it interchangeable with occupancy statistics.

For n >= 4, the endpoint-mixture family generally does not parameterize the
whole feasible moment body. A finite-atom representation of each moment
vector does not imply that extrema of the nonlinear persistence objective
occur at extreme measures or at principal representations alone.

## Result on the 16 specified cases

All 50 derivative Bernstein coefficients are strictly negative for each of
P0, PR, and PC in every case. All 50 are strictly positive for PR-PC in
every case. Thus persistence under each action decreases with m3, while the
advantage of buffering over connectivity increases. Global extrema over G
occur at the endpoint laws in these cases; this conclusion follows from
the derivative certificates, not from the laws' extremality for m3 alone.

Buffering wins throughout the moment class in seven cases; connectivity
wins throughout in eight. The remaining case is
`n=3, p=.1, q=.03, e0=.02, e1=.9, c=.1, T=50`. Its action difference runs
from approximately `-0.0234163336923` to `+0.0172018595903`, with exactly
one zero. The certified lambda bracket is

```
148218197/268435456 < lambda_star < 592872789/1073741824.
```

Mapping through `m3 = 9/1000 + lambda*(221/9000 - 9/1000)` gives a
threshold approximately `0.017589090406`. Connectivity wins below the
threshold and buffering above it, with a tie at the threshold. Here m3 is
the annual probability of three simultaneous environmental hits.

These are whole-class results in G for 16 fixed ecological parameter
combinations. They do not characterize all combinations of mortality,
colonisation, horizon, or patch number. Exact rational arithmetic removes
floating-point sign ambiguity; validity still depends on the model and
correct implementation of its transition probabilities.
