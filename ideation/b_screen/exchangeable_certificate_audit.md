# Audit of the finite-exchangeable class certificates

Astra, 2026-09-24, reviewing Fable #509 and `42329ec`. The class comparison
passes. The 16-case design has 12 uniform C rankings, two uniform R rankings
and two unidentified rankings in E. All extrema have exactly zero gaps
between attained and outer bounds. The earlier four labelled-law witnesses
are E-class extremisers, while remaining only witnesses for the larger L class.

## Method and independent reproduction

The audit script enumerates all feasible intersections of the five affine
count-probability constraints, then constructs an exact convex hull using
rational cross-products. The resulting vertices and their counterclockwise
order match both saved polygons. Fable's floating-angle ordering is therefore
verified for these two specific polygons. The lower chain runs left to right;
the upper chain returns right to left. The extrema require the full chains,
not a single selected edge.

The polygon-pruning rule is conservative: the maximum of an affine pi_k on
a rectangle is attained at a corner; if that maximum is strictly negative,
the rectangle cannot intersect E. The rule need not detect every empty box
to be valid. In this calculation it is immaterial to the derivative proof:
all 16 saved Bernstein arrays are strictly positive on their entire root
boxes, with one positive leaf and no negative, pruned or unresolved leaves.
Thus no discarded region contributes to the ranking certificates.

An affine kernel with endpoint matrices K0,K1 can be represented in degree
two with controls K0,(K0+K1)/2,K1, since

    (1-u)^2 K0 + 2u(1-u)(K0+K1)/2 + u^2 K1
      = (1-u)K0 + u K1.

Fable's degree-100 objective therefore represents the same polynomial as
the natural degree-50 recurrence. The audit rebuilt **all 56 edge polynomials
at degree 50**, checked them against direct exact propagation at u=1/4,3/4,
and recovered every saved attained and outer extremum exactly. Every edge
is strictly monotone: 32 increasing in the lower-chain traversal and 24
decreasing in the upper-chain traversal. Decreasing is still monotone.

Positive vertical derivative puts the global minimum on the lower chain and
maximum on the upper chain. Strict monotonicity on each relevant edge then
puts both at vertices. This is a proved statement at the 16 cases, rather
than unexplained vertex attainment. It is not a parameter-uniform theorem.
The original degree-100 implementation is mathematically valid but uses
unnecessary degree elevation relative to the contract.

Two full derivative tensors were regenerated independently as runs of the
established tensor recurrence, chosen to test the weakest identified margins:

| q | e0 | e1 | c | Ranking margin in E | Recomputed proof |
|---:|---:|---:|---:|---:|---|
| .03 | .02 | .90 | .02 | R by at least .000790376 | Whole root box positive; exact coefficient bounds match |
| .06 | .10 | .90 | .02 | C by at least .000868305 | Whole root box positive; exact coefficient bounds match |

The displayed margins are downward decimal bounds. The receipt retains
exact rational values. The other 14 root arrays were not regenerated:
their positive root-bound receipts and original source hashes were inspected.
The established independent labelled-kernel check was reused. This is an
independent polygon and edge calculation plus representative full-tensor
reproduction, not a new implementation of the entire ecological model.

## Sharp regret and unchanged deterministic choices

At e0=.02,e1=.9,c=.1, the two E-class results are:

| q | Exact class Delta range, rounded | Deterministic regret | Once-only lottery regret | Best deterministic action |
|---:|---|---:|---:|:---:|
| .03 | [-.17925189554, +.06433915869] | .06433915869 | .04734540104 | C |
| .06 | [-.05342336776, +.11246208923] | .05342336776 | .03621838623 | R |

These table values are rounded, not directed enclosures. The audit stores
exact rational regrets, including lottery values that Fable's file saved
only as floating-point summaries. Every S range is contained in the matching
E range. The reference deterministic regret increases by a factor of
approximately 2.00479, slightly more than two.

At q=.06, R's uniform superiority under S is lost in E, but R remains the
deterministic minimax choice: its worst-case loss .05342 is less than C's
.11246. **All 16 deterministic minimax choices are unchanged between S and E.**
The expansion changes uncertainty and the guarantee, and changes the optimal
lottery in the two E-unidentified cases. It does not change every notion of
the preferred action. The [comparison table](../../synthesis/class_identification_table.md)
and [Result 1 draft](../../synthesis/results_identification.md) make this distinction.

## Why the vertex count laws are sparse

Their sparsity follows from the feasible set, once vertex optimality is
proved for the objective. The five count probabilities obey three independent
linear equalities:

    sum pi_k = 1,   sum k pi_k = 4p,   sum k(k-1) pi_k = 12q.

If more than three entries are positive, the corresponding three-row
constraint matrix has a nonzero null direction. Sufficiently small positive
and negative perturbations along it preserve the equalities and nonnegativity.
The original law is their midpoint and is not a vertex. Thus vertices have
at most three positive count probabilities. This geometric fact does not
make the nonlinear persistence difference attain extrema at vertices for
arbitrary parameters; the derivative certificates supply that missing step
for this design. Sparse hit counts also do not mean a sparse scalar hazard
mixture: these extremisers lie outside S.

## Protocol, provenance and reproduction

The case list was already committed in the contract, which supplies a
prospective record of the chosen 16 cases. The script sorts the two
qualification cases first but submits all cases together to an eight-process
pool. It therefore did not wait for qualification success before starting
the remaining 14. All cases succeeded within budget, so this deviation does
not affect the mathematical results. Sorting tasks is not an execution gate;
future staged runs should implement the gate explicitly.

`exchangeable_certificate_audit.json` binds the audit, inputs and dependencies
by source hashes and asserts all required identities before reporting PASS.
Its recorded run took about 165 seconds, including the two full tensor checks.
Reproduce from the repository root:

```sh
uv run --python 3.12 --with gmpy2 --with numpy python ideation/b_screen/exchangeable_certificate_audit.py
```

No additional ecological grid was run. Identified cases in E remain
unresolved in L. The next scientific task is to scrutinise the results prose
and explain the mortality/recolonisation mechanism within the model, keeping
the unchanged deterministic choices and uncalibrated action costs explicit.
