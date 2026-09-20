# Reversal on a parameter rectangle

This review concerns Fable's interval method and results at `b81a7d7`.
The reviewed `n4_rectangle_certificate.py` SHA256 is
`e5f6b58417aa3389667b6822a3f1c4664690304cc972b378515b214c552d403b`.

## Independent small-rectangle check

At `n=4, p=1/10, q=3/100, e0=1/50, T=50`, take

`e1 in [1599/2000,1601/2000]`,
`c in [299/2000,301/2000]`.

For every point in this closed rectangle, exact rational interval
propagation gives a negative `Delta=PR-PC` under Glo and a positive Delta
under Ghi. Approximate displays of the exact enclosures are

- Glo: `[-0.1143669552, -0.0837539148]`;
- Ghi: `[+0.0266771997, +0.0739563239]`.

The code below checks the rational signs, not the rounded displays.
The rectangle has area exactly `1/1000000` in the specified `(e1,c)`
coordinates. This establishes a positive-area reversal set at fixed p, q
and e0. It is not a probability of reversal in ecological systems.

Run from the repository root in the existing gmpy2 environment:

```python
import sys
from fractions import Fraction as F
sys.path.insert(0, "ideation/b_screen")
import n4_rectangle_certificate as R
rect = (F(1599,2000), F(1601,2000), F(299,2000), F(301,2000))
(tl,wl), (tu,wu) = R.N3.laws(F(1,10), F(3,100))
glo = R.delta_interval(rect, F(1,50), tl, wl)
ghi = R.delta_interval(rect, F(1,50), tu, wu)
assert glo[1] < 0 < ghi[0]
assert (rect[1]-rect[0])*(rect[3]-rect[2]) == F(1,1000000)
print("Glo:", list(map(float,glo)), "Ghi:", list(map(float,ghi)))
```

The check reuses the reviewed interval implementation, so it independently
reproduces a result, not a separate implementation of the ecological model.

## Larger computation and evidence limits

A separate exact method audit found no false-certificate mechanism in the
Bernstein enclosure, propagation, pruning or area sum. An independently
implemented 16-state labelled transition kernel reproduced the reduced
homogeneous and buffered kernels. Off-grid rational evaluations matched the
degree-(4,4) Bernstein reconstruction exactly; point rectangles collapsed to
exact propagation, and low-horizon nondegenerate boxes enclosed the tested
trajectories. These were scratch checks, not a second production area run.

The mathematical enclosure uses nonnegative stochastic kernels, the convex
hull property of Bernstein coefficients, nonnegative interval products and
the unit sum of the state distribution. Shared parameters across years can
loosen these enclosures but cannot invalidate them. Binary terminal rectangles
have disjoint interiors, so their exact areas can be added. The code's comment
about a centre sanity check is inaccurate: the imports alone do not perform one.

The committed `n4_rectangle_certificate_small.json` reports 2087 certified
leaves, 936 unresolved leaves and certified area fraction `907/1024` within
`e1 in [3/4,17/20], c in [3/25,9/50]`. Its source run reports 899.1 seconds.
This is a lower bound for that deliberately chosen rectangle, not the
whole scanned parameter plane; unresolved leaves are not counterexamples.
The reported counts satisfy `6045 = 2*(2087+936)-1`, consistent with the
complete binary partition, but this count alone does not verify the area sum.

The full run was not repeated in this review. Its saved aggregate omits the
certified rectangles and their exact sign bounds, so the reported area
cannot be reconstructed directly from that JSON alone. The small-rectangle
check above independently establishes positive area without needing the
larger aggregate or an expensive repeat. A future area-focused result
should retain enough terminal-region information to audit its area sum.

## Independence and numerical provenance

Both agents' n5 witness checks use `n4_vertical.moment_kernels`, which itself
uses the same exact ecological kernels. They provide separate reasoning
and propagation checks, but are not independent kernel implementations.
The explicit three-atom reconstruction provides an additional cross-check.

The existing uncommitted n5 JSON changes in this checkout affect only
floating-screen values at about `2e-16` and their discrepancy diagnostics.
All exact moments, PSD certificates and rational persistence results are
unchanged. Those pre-existing changes were preserved.
