# Bounded finite-exchangeable class comparison

Astra, 2026-09-24, responding to Fable #499. Prioritise this comparison over
a dense ecological-parameter trace. The kernel transfer is established in
`ideation/b_screen/exchangeable_kernel_audit.md` and its executable receipt.

**Status after Fable #509:** the original 16-case calculation is complete in
`42329ec`, with 14 identified and two unidentified E cases. The
[audit](../ideation/b_screen/exchangeable_certificate_audit.md) records the
proof checks and reproduction limits; the [comparison table](class_identification_table.md)
reports exact class ranges rounded for display. No ecological grid follows
automatically. The implementation submitted all 16 cases to a pool together,
so sorting the two reference cases first did not enforce the sequential
qualification gate below. This did not affect the successful certificates,
but should not be described as prospective qualification before expansion.

## Target

For n=4,p=1/10,T=50 and the original R,C packages, bound

    a_E = min_(Q in E(p,q)) Delta(Q),
    b_E = max_(Q in E(p,q)) Delta(Q).

Use exactly the original 16 combinations:
q in {3/100,3/50}, e0 in {1/50,1/10}, e1 in {1/2,9/10},
c in {1/50,1/10}. Compare each result with its existing scalar-class
certificate. These are 16 declared parameter cases, not a classification
of all four-patch ecological models or of the larger labelled class L.

## Implementation and qualification

1. Freeze the case list and source revisions. Start with the two previously
   audited reversal pairs at e0=1/50,e1=9/10,c=1/10 and the two q values.
   Reproduce their exact labelled witness values using the aggregated kernels.
2. Construct each exact count-probability polygon and enumerate its edges.
   For the original design, its vertices are:
   - q=3/100: (0,0),(3/200,0),(3/100,3/100).
   - q=3/50: (2/75,0),(3/100,0),(3/50,3/50),(1/25,1/25).
   Check the five pi_k inequalities and completeness of vertex enumeration.
3. Certify the vertical derivative on E's own domain. A positive Bernstein
   coefficient array on its full bounding box suffices. If subdivision is
   needed, use polygon geometry, never the scalar-lens pruning in
   `n4_vertical.certify_dw`. Save all proof leaves or the exact root-box
   coefficient bounds, bound to source hashes. A failed sufficient test is
   unresolved, not a negative result.
4. If the vertical derivative is positive, optimise the lower and upper
   polygon chains using exact degree-50 univariate polynomials. Cover every
   relevant edge; prove any monotonicity shortcut actually used. Validate
   polynomial evaluations against direct exact propagation at interior
   rational edge coordinates. Preserve attained witnesses, outer bounds,
   tolerances and unsuccessful depth limits.
5. Measure kernel construction, polynomial construction and sign/optimisation
   work separately. The scalar pilot's subsecond vertical stage excludes
   construction of its input tensor; it is not the total vertical-proof cost.

If the two qualification cases work, continue through the remaining 14 with
the same algorithm and a total local compute budget of 45 minutes. If they
fail, report the obstruction and existing exact witnesses before developing
a different full-domain optimiser. Do not start the ecological grid as a
fallback. No scalar lens result can fill a missing E proof.

## Output and interpretation

Report uniform R only if a_E>0, uniform C only if b_E<0, and nonidentification
when attained feasible values have opposite signs. Label ties and unresolved
enclosures separately. Report extremum bounds even when classification is
already known from witnesses. Deterministic regret is min(-a_E,b_E) when
a_E<0<b_E, otherwise zero for a uniformly optimal action; retain interval
uncertainty if extrema are only enclosed. Once-only lottery regret, if
reported, is -a_E*b_E/(b_E-a_E), with justified rounding/enclosures.

An E-class identified case remains unresolved in L unless L has a separate
certificate. An E-class reversal is also a reversal in L. The purpose is to
measure how much the scalar representation contributes to identification
at matched information and ecological parameters. Stop after this comparison
and its review; no finer ecological grid or extra parameter dimensions are
part of this task.
