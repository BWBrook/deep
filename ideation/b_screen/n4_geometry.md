# Four-patch moment geometry and audit of the proposed optimisation

Scope: a fixed hazard law G on [0,1], iid across years, with known
`p=E(theta)`, `q=E(theta^2)`, and `0<p<1`, `p^2<q<p`.
Write `x=E(theta^3)` and `y=E(theta^4)`. This note corrects the
boundary argument in correspondence thread 49, message 369.

## Exact feasible set

Define

```
xL = q^2/p
xU = q - (p-q)^2/(1-p)
L(x) = (x^2 - 2*p*q*x + q^3)/(q-p^2)
U(x) = x - (q-x)^2/(p-q).
```

The entire feasible set is `xL <= x <= xU`, `L(x) <= y <= U(x)`.
Necessity follows from nonnegativity of expected squares: the matrices

```
H = [[1, p, q], [p, q, x], [q, x, y]]
Q = [[p-q, q-x], [q-x, x-y]]
```

must be positive semidefinite. Their Schur complements give the two y
bounds, and their difference factors as

```
U(x)-L(x) = p*(1-p)*(x-xL)*(xU-x)/((p-q)*(q-p^2)).
```

Sufficiency can be seen constructively from the boundary laws below:
at each x in the interval, mix the lower and upper laws to obtain every
y between them. The two laws share p, q, and x.

On the lower boundary, the two support points z- and z+ are roots of

```
z^2 - ((x-p*q)/(q-p^2))*z + (p*x-q^2)/(q-p^2) = 0.
```

Their masses are `(z+-p)/(z+-z-)` and `(p-z-)/(z+-z-)`.
On the upper boundary, the support is `{0,a,1}`, where

```
a = (q-x)/(p-q)
w_a = (p-q)/(a*(1-a))
w_0 = 1-p-(p-q)/a
w_1 = p-(p-q)/(1-a).
```

These support points lie in [0,1] and masses are nonnegative over the
stated x interval. At its endpoints the two boundaries coincide with
Glo and Ghi. A generic three-atom law, in contrast, has three free
parameters after normalisation and fixing p,q, not one. The upper
boundary is one-dimensional because two support locations are fixed.
Sparsity of a representing law alone does not locate a moment vector
on the boundary.

## Consequences for optimisation

Each integrated transition kernel is affine in (x,y). A T-step objective
is a bivariate polynomial of total degree at most T. Restricting it to
either boundary substitutes a quadratic y(x), giving a univariate
polynomial of degree at most 2T, hence at most 100 for T=50.

The laws Glo and Ghi simultaneously minimize and maximize the third and
fourth moments. Thus, if both coordinate derivatives of an objective
have the same sign throughout an enclosing box, its extrema on the
moment set are already attained at those two feasible corners. No
boundary search is then needed. The same conclusion holds if these derivative
signs are certified throughout the feasible set itself: it is convex, so the
straight segments from its coordinate-wise extreme corners to any feasible
point stay feasible and have the required coordinate directions.

If only the sign of the y derivative is certified throughout the feasible
set, vertical line segments remain feasible, so a maximum and minimum
must lie on opposite boundary curves. This alone can reduce the task
to two univariate optimisations; a certificate for the x derivative is
not necessary for that reduction.

If derivatives change sign, interior stationary points remain possible.
A failed Bernstein coefficient test is inconclusive; it is not evidence
of a sign change. Even a wrong-sign point in an enclosing box is not
a counterexample unless its moment vector is feasible.

A rectangle parameterisation of the exact feasible set is

```
x = xL + t*(xU-xL)
y = L(x) + s*(U(x)-L(x)),   0 <= t,s <= 1.
```

To test signs of the original partial derivatives J_x and J_y, substitute
this map into those derivatives. Derivatives of the composed objective
with respect to t and s are different quantities and include chain-rule
terms. This mapping also raises polynomial degrees, which must enter
any cost estimate for Bernstein certification.

## Recommended sequence

1. Build and check exact affine moment kernels using rational interpolation.
2. Try matrix-sensitivity derivatives at feasible interior points before
   constructing full high-degree coefficient arrays. Use exact arithmetic
   to certify any proposed counterexample and its moment feasibility.
3. If vertical monotonicity is plausible, seek a certificate over the whole
   feasible set. If successful, certify the degree-at-most-100 boundary
   objectives; otherwise retain interior boxes in branch and bound.
4. Separate rigorous outer bounds from attained feasible objective values.
   Rounding an arbitrary floating-point law to rationals does not automatically
   preserve p and q; feasibility must be repaired or proved.

Opposite-sign feasible examples already prove an unidentified ranking.
Whole-class optimisation quantifies the range and is needed to certify
absence of reversals, rather than to re-prove existing reversal witnesses.

## Exact counterexample to the proposed coordinate-wise monotonicity

The conjecture that Delta=PR-PC increases with both moments is false within
the specified 16-case ecological parameter grid. A 400-point floating-point
scan (25 interior lens locations per parameter case) found five locations
violating at least one proposed derivative sign. The first was checked in
exact arithmetic, preserving the earlier 16 successful midpoint probes.

Take `n=4, p=1/10, q=3/50, e0=1/50, e1=9/10, c=1/10, T=50` and
`x=29/500`, `y=288411/5000000`. This point is strictly inside the moment
set. The exact derivative of Delta with respect to x at fixed y is negative
(approximately `-0.03291334911919705`).

Two explicit feasible moment vectors share p, q, and y and have
`x_minus=579999/10000000` and `x_plus=580001/10000000`. Exact propagation
gives Delta approximately `0.1099032596679997` and
`0.10990325308532987`, respectively, with the second value strictly smaller.
All PSD feasibility inequalities and probability-kernel checks pass exactly.
The script and JSON record exact rational sensitivities and outcomes.

Independent reconstruction as mixtures of the two boundary laws at each x,
followed by the original floating-point transition implementation, reproduces
the two differences to approximately 2e-15. Both differences are positive:
this refutes coordinate-wise monotonicity, not the action ranking at these
points. It also does not prove that global extrema are away from Glo/Ghi.

Run `python3 ideation/b_screen/n4_monotonicity_probe.py` to reproduce the
midpoint checks and exact counterexample. The broader exploratory scan is
not a certificate; the feasible rational counterexample is the decisive
result. The proposed optimiser should therefore test vertical monotonicity
or retain interior search rather than rely on all coordinate signs agreeing.
