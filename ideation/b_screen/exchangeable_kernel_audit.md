# The moment kernels extend to finite exchangeability

Astra, 2026-09-24, checking Fable #499. The extension is valid. It does not
transfer any certificate whose domain was only the scalar-mixture lens.

## Why the same coefficient matrices apply

Let X be the four binary hit indicators. Conditional on X, a labelled
one-step transition probability is a multilinear polynomial in X: each
occupied patch contributes one survival or death factor, and subsequent
recolonisation depends on the survivor set. Write the coefficient matrices
of this polynomial as B_A, indexed by subsets of patches. For an exchangeable
hit law, every specified k-subset has the same joint-hit probability tau_k.
Therefore

    E[K(X)] = sum_(A subset {1,2,3,4}) B_A tau_|A|
            = sum_(k=0)^4 A_k tau_k,   A_k = sum_(|A|=k) B_A.

Under iid Bernoulli(theta) hits this is the same polynomial with
tau_k=theta^k. Thus extracting A_k by the existing exact polynomial
interpolation gives the same coefficient matrices for every exchangeable
law. A scalar mixture is not needed for this identity.

For C, permutation invariance makes the total occupied count a valid state.
For R, permutations of the three unbuffered patches preserve the model, so
the state is (buffered occupancy, number of other occupied patches). The
buffered status must be retained; a total-count-only R chain would not be
justified. Independent draws of the hit vector across model steps are still
assumed.

## Executable check on a spanning basis

`exchangeable_kernel_audit.py` uses the separately derived labelled kernel
from `labelled_class_audit.py`. For each fixed hit count k=0,...,4, it averages
uniformly over all patterns having k hits, verifies exact strong lumpability,
and compares the reduced kernel with the existing moment formula using

    tau_j = choose(k,j) / choose(4,j).

The five fixed-count laws span the entire exchangeable-law simplex. Equality
at each vertex and linear averaging imply equality for every mixture of
these laws, including the polygons at fixed p,q outside the scalar class.
The check covers all eight original combinations of e0,e1,c and both actions:
80 reduced kernels, 3,560 matrix entries, and another 4,760 comparisons of
rows within aggregate states. All comparisons agree exactly. The receipt
records source hashes and elapsed time. This is a check of the kernel
extension, not an optimisation over E.

At e0=.02,e1=.9,c=.1, the script also constructs degree-50 Delta polynomials
on all three edges at q=.03 and all four edges at q=.06. Direct propagation
at edge coordinates 1/4 and 3/4 matches the polynomial exactly for both
actions. Every edge is strictly monotone, with derivative signs certified
by its exact Bernstein coefficients. The previously audited negative and
positive witnesses are consequently the extrema **on each polygon boundary**.
The receipt stores endpoint values, derivative bounds and coefficient hashes.
An interior-domain proof is still required to establish full E-class extrema.

Reproduce from the repository root:

```sh
uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/exchangeable_kernel_audit.py
```

## Domain correction required for E certificates

At fixed p,q, let x=tau_3 and y=tau_4. The count probabilities are

    pi_0 = 1-4p+6q-4x+y
    pi_1 = 4p-12q+12x-4y
    pi_2 = 6q-12x+6y
    pi_3 = 4x-4y
    pi_4 = y.

Their nonnegativity defines E. Its bounding box is larger than the scalar
moment lens. A uniformly positive derivative Bernstein array on the whole
new box would certify positive vertical derivative on E, even if some box
points do not correspond to probability laws. Mixed coefficients are
inconclusive, not proof of a sign change.

The existing `n4_vertical.certify_dw` prunes rectangles disjoint from the
**scalar lens**. Calling it unchanged with E's bounding box would discard
the part of E central to the scientific comparison. Use a whole-box test,
or replace its geometry test with valid polygon pruning. A safe sufficient
pruning test is: some affine pi_k is strictly negative even at its maximum
over the candidate rectangle. Such a rectangle cannot intersect E.

With a valid positive vertical derivative certificate, minimise Delta over
the lower polygon chain and maximise it over the upper chain. Each edge
has an affine kernel, so persistence and Delta have degree at most T=50
in its interpolation coordinate. All required edges and vertices must be
covered; neither vertex-only optimisation nor edge monotonicity is automatic.
Retain attained and outer bounds separately until any exact coincidence is
proved. The proposed bounded computation is in
`synthesis/exchangeable_certificate_contract.md`.
