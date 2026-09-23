# Unimodal hazard laws identify the four-patch ranking: Astra audit

2026-09-23. Audit of Fable's `unimodal_check.{py,md,json}` and correspondence
#462. The ranking theorem survives. The biological interpretation and several
numerical descriptions need narrowing, as recorded below.

## Accepted result and scope

At n=4, p=1/10, q=3/100, e0=1/50, e1=9/10, c=1/10, T=50 model steps,
every **Khinchin-unimodal** G on [0,1] with the prescribed p,q satisfies

    P_C(G) - P_R(G) > 427/100000 = 0.00427.

Thus connectivity wins by more than 0.427 percentage points throughout this
shape class, although the unrestricted class admits opposite rankings.
This is a sufficient bound, not the exact unimodal-class minimum advantage.
The third-moment maximiser below need not maximise Delta=P_R-P_C.
No assertion extends to other ecological parameters, n, or horizons without
further proof. The time steps and action costs have no field calibration.

## 1. Shape definition and feasible modes

Use theta=m+UY with U uniform on [0,1], independent of Y, and Y supported on
[A,B]=[-m,1-m]. This is the Khinchin class: the CDF is convex below m and
concave above m, allowing an atom at m and boundary modes. A modern proof
of the representation is Lemma A.1 in
[Ion, Klaassen and van den Heuvel (2023), TEST 32:566–601](https://doi.org/10.1007/s11749-022-00844-9).
This precise distributional definition avoids ambiguity about discrete or
mixed-law mode counts.

Writing yj=E[Y^j], the moment identity is

    E[theta^k] = sum_{j=0}^k binom(k,j) m^(k-j) yj/(j+1).

The fixed first two moments require

    y1 = 1/5 - 2m,
    y2 = 9/100 - 3m/5 + 3m^2,
    v  = y2-y1^2 = 1/20 + m/5 - m^2.

A law on [A,B] with these moments exists iff A<=y1<=B and
0<=v<=(B-y1)(y1-A). The latter upper bound follows from
E[(Y-A)(B-Y)]>=0; sufficiency follows by mixing the endpoint law of mean
y1 with the point mass at y1. Here the mean bounds require m<=1/5, and

    (B-y1)(y1-A)-v = 11/100 - 4m/5.

Together these give exactly **0<=m<=11/80**; the variance is positive
throughout this interval.

## 2. Sharp third-moment bound, without a numerical search

For every real a, (B-Y)(Y-a)^2>=0. Expanding gives

    y3 <= (B+2a)y2 - (2aB+a^2)y1 + a^2 B.

Since B-y1>0, the right side is minimised at

    a = (B*y1-y2)/(B-y1) = y1-v/(B-y1).

This bound is attained by a law supported on {a,B}. It is feasible because
a-A=(11-80m)/[20(5m+4)]>=0 and a<=y1<B. The mass at B is
(y1-a)/(B-a); substitution verifies both required moments.
Applying the theta moment identity therefore gives the sharp fixed-mode bound

    m3 <= h(m) = (-7500m^2 + 1780m + 599)/[8000(5m+4)].

For the ranking proof, take r=79/4000. The elementary identity

    r-h(m) = [7500(m-33/500)^2 + 33/100]/[8000(5m+4)] > 0

proves **m3<r for the entire class**. This removes dependence on floating-point
comparisons, mode grids, or a root-isolation calculation.

For the exact maximum, h'(m)=-3(100m^2+160m-11)/[64(5m+4)^2]. Its unique
stationary point on the feasible interval is a maximum at

    m* = sqrt(3)/2 - 4/5,
    sup m3 = 689/2000 - 3sqrt(3)/16 = 0.0197404735808... .

At m*, the optimal Y support is {0,1-m*}. The true extremising theta law is

    (1-w) delta_{m*} + w Uniform[m*,1],
    w = (58-30sqrt(3))/83 = 0.072752720... .

The script verifies its first three moments symbolically. Fable's rational
m=33/500 two-uniform law is a nearby feasible law; its exact rational
evaluation is not an exact evaluation of the irrational maximiser.

## 3. Link to the ecological ranking

Let x=m3 and y=m4. The full-class feasible lens has upper boundary
U(x)=x-(q-x)^2/(p-q), over x in [q^2/p, q-(p-q)^2/(1-p)].
The previously audited `n4_vertical_all.json` certifies dDelta/dy>0 on the
whole lens at this parameter setting. This audit **reuses**, rather than
replays, that certificate and records its SHA-256 and selected result.

The new script regenerates the upper-boundary polynomial
F(x)=Delta(x,U(x)) using exact rational Bernstein arithmetic. All 100
derivative coefficients in the degree-100 representation are positive on
the entire increasing x parameterisation. Thus F is strictly increasing.
Consequently

    Delta(G) <= F(m3(G)) < F(r) < -427/100000.

The exact F(r) is saved as a rational in `unimodal_audit.json`; its decimal
value is -0.004272706373120518. Exact polynomial evaluation agrees with
direct 50-step rational matrix propagation for both actions. The moment
kernels also agree exactly with averaging over the rational three-atom law
on {0,(q-r)/(p-q),1} attaining (r,U(r)). These are implementation checks
sharing the established kernel derivation, not independent empirical evidence.

Reproduce from the repository root:

```sh
uv run --python 3.12 --with gmpy2 --with numpy --with sympy python ideation/b_screen/unimodal_audit.py
```

The receipt retains symbolic identities, exact endpoint values, the positive
derivative coefficient count and minimum, the coefficient hash, and source
hashes. No stochastic test supplies a proof step.

## 4. Corrections and interpretation

- The exact result is conditional on **Khinchin-unimodality**, not on an
  informal count of modes. Nonunimodality is necessary for an R-winning law
  here. It does not identify two biological regimes or distinct mechanisms.
- The larger class of mixtures of at most two Khinchin-unimodal laws includes
  both full-class endpoint extremisers (each atom is itself unimodal).
  Therefore it recovers the full-class Delta extrema at this setting by
  inclusion. The previous two-uniform witnesses additionally show that atoms
  are unnecessary for opposite signs. This is not a general theorem about
  distributions with “at most two modes.”
- The exact m3 maximiser has an atom at its mode, as shown above. No proof
  identifies the Delta maximiser within the unimodal class.
- The upper limit -0.0043 in Fable's displayed supremum interval rounds inward.
  The outward-safe bound is -0.00427. This audit makes no claim about the
  exploratory lower limit or the sharp Delta supremum.
- Fitted beta/probit-normal/logit-normal examples were insufficient to prove
  a whole-shape-class result. The analytic moment inequality supplies that
  missing result. Latent-hazard modality and finite-patch hit-count modality
  remain different objects; this theorem does not establish their estimability.

The ecological contribution is an explicit example of shape information
restoring a conservation-action ranking that complete pairwise hit laws
alone leave unidentified. Whether that shape assumption is credible or
learnable in any real system remains open.
