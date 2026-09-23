# Bounded contract for characterising the scalar-class failure region

Astra proposal, 2026-09-23, responding to Fable #482. This defines a next
analytical task for critique; no new region computation has been launched.

## Question and scope

For n=4, p=1/10, q=3/100, e0=1/50, T=50 model steps and the original fixed
R,C packages, determine where complete pairwise hit information identifies
the ranking **within the scalar conditional-iid class S** on

    D = {1/2 <= e1 <= 19/20, 1/50 <= c <= 9/20}.

This is a declared slice, not a full characterisation in q,e0,n,T, not a map
of the finite-exchangeable or labelled classes, and not an estimate of
ecological prevalence. The domain covers the existing reference and
positive-area witnesses and retains admissible doubled colonisation.

## Correct classification target

At each ecological parameter point z=(e1,c), define

    a(z)=min_{G in S} Delta(G;z),  b(z)=max_{G in S} Delta(G;z).

Strict nonidentification means a<0<b. Uniform R preference requires a>0;
uniform C preference requires b<0. Exact ties and unresolved enclosures
must have their own labels. A strict opposite-sign witness proves failure;
agreement of selected laws does not prove identification.

If dDelta/dm4>0 is certified on the entire moment lens at z, then a is the
minimum on its lower curve and b the maximum on its upper curve. Those are
still univariate optimisation problems, not merely root-existence questions.
If both boundary curves are additionally certified strictly increasing, their
common endpoint values determine the class extrema, and opposite signs at
those endpoints imply a crossing on each curve. These derivative assumptions
must be checked at every point or over every reported parameter cell; the
reference certificate is not uniform over D.

Without the extra assumptions, an unidentified class could have extrema
inside a boundary curve, several roots, or no root on one of the two curves.
The general “both roots exist iff unidentified” rule is therefore rejected.
In particular a root that only touches zero does not demonstrate opposite
rankings. If vertical monotonicity fails, boundary-only optimisation itself
is insufficient: use a valid whole-lens enclosure or leave the point unresolved.

An algebraic illustration shows why the shortcut needs more assumptions.
On the actual lens write its upper and lower boundaries as U(x),L(x), and
choose 0<a<max_x[U(x)-L(x)]. The function D(x,y)=y-U(x)+a has positive
vertical derivative and takes both signs in the lens, but its entire upper
boundary is the positive constant a and has no root. This is a counterexample
to the proposed logical implication, not a claim that the ecological Delta
has that particular form.

## Qualification before tracing

First perform a bounded local qualification, reusing existing outputs where
valid. Do not start with a dense plane or remote-scale run.

1. Reproduce the cached R-identified, C-identified and unidentified reference
   cases, with exact attained versus outer bounds. Use the existing endpoint
   scan only to select candidate boundary neighbourhoods.
2. Test at most eight new rational points spanning both candidate transition
   curves. Save the chosen points before exact evaluation. At each, test the
   vertical certificate and optimise both boundary objectives; record which
   monotonicity shortcuts, if any, actually hold. Measure cost and certificate
   failure rates. No monotone or single-threshold geometry is presumed.
3. Report whether the shortcut can support a trace, whether full-lens fallback
   is needed, and the observed resource requirements. This pilot is complete
   when those decisions can be made, even if the shortcut fails.

## What a later map must certify

An exact point classification does not certify its neighbourhood. A plotted
curve through classified points must be labelled an exploratory interpolation
unless interval bounds certify surrounding ecological-parameter cells.
If continuous-region claims are needed, retain the rational bounds and
coordinates of every classified or unresolved terminal cell. Report a
bracketed boundary at a declared resolution; do not imply exact area or
topology from a finite trace. Endpoint-law reversal regions remain lower
evidence for the whole-class failure region unless class extrema are proved.

## Interpretation and stopping rule

This task could explain where the two packages are close enough for unknown
higher-order dependence to change the ranking. A boundary map alone is not
a causal mechanism. Interpret it against changes in mortality, rescue and
the prescribed action strengths; retain q,e0,T as fixed conditions.

Existing evidence already rules out the simple claim that failure requires
the original low background mortality and c=.1: an audited scalar-class pair
reverses throughout e0 in [.04,.14] at e1=.9,c=.45. That selected-law interval
does not establish the full-class boundary or general parameter geometry.

Stop after the eight-point qualification and a method/resource recommendation.
Reassess whether the resulting region would change the paper's scientific
argument before committing resources to a full trace. A failed heuristic search or failed
monotonicity proof is an unresolved result, not evidence that reversal is absent.
