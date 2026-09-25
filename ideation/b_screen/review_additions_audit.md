# Targeted audit of two additions in the independent review

## Decisive findings

Both new numerical claims reproduce using the repository's established exact kernels and Bernstein machinery. This is a new local audit, not an inspection or execution of the reviewer archive, which arrived after this method was fixed.

| q | certified inf PR | certified inf PC | maximin action |
|---:|---:|---:|:---:|
| .03 | 0.3857762490 | 0.3214370903 | R |
| .06 | 0.2273638798 | 0.1149017905 | R |

For each action and q, an exact tensor-Bernstein certificate shows that persistence decreases strictly with `tau4` throughout the finite-exchangeable polygon E. The global minimum therefore lies on the polygon's upper edge chain. Exact univariate Bernstein branch-and-bound encloses every relevant edge minimum to tolerance `1e-9`; the global attained values equal the separately propagated labelled high-Delta witness exactly. The lower bound and attained value coincide at those endpoint witnesses in these cases.

At q=.03, mixing 90% of each named E extremiser with 10% Glo makes all 16 labelled hit patterns strictly positive while preserving every marginal (`p=.1`) and pair probability (`q=.03`). Exact labelled propagation gives:

| source extremiser | PR-PC | ranking | full labelled support |
|:---|---:|:---:|:---:|
| low-Delta | -0.172970 | C | yes |
| high-Delta | +0.050918 | R | yes |

Thus the opposite rankings are not caused solely by zero-probability labelled patterns. These are attained constructions, not new optima and not evidence that either law is ecologically realistic.

## Evidence and assumptions

The estimand is 50-step regional persistence from full occupancy under each complete action policy, with hit vectors drawn independently at each step from one fixed law Q. The maximin prediction target is `argmax_a inf_Q P_a(Q,50)` over the four-patch finite-exchangeable class with exact `p,q`. Pairwise hit probabilities are treated as known; there is no sampling or detection model, temporal dependence, cost calibration, or inference to a field system.

The audit imports `n4_vertical.py` for exact affine moment kernels and tensor Bernstein propagation, `exchangeable_certificate.py` for E geometry and polygon pruning, `n4_boundary.py` for exact edge polynomials and enclosures, and `labelled_class_audit.py` for the separately derived 16-state survivor-subset kernel. Source hashes are in the JSON receipt. These shared biological and mathematical dependencies mean this is a strong internal cross-check, not fully independent software replication.

The maximin/minimax difference at q=.03 is sensitivity to the ambiguity decision criterion applied to the same persistence model. It supplies no additional empirical ecological support for the actions, hit laws, parameter values, or costs.

## Verification

Each action tensor was checked exactly against direct propagation at two strictly interior rational points of each E polygon. Every relevant edge polynomial was checked exactly at `u=1/4` and `u=3/4`. The named extremal witnesses were also propagated through the labelled kernel and matched both the aggregated chain and the global attained minima exactly. The full-support laws passed exact normalization, four marginal, and six pair constraints before propagation.

Reproduce with:

```sh
uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/review_additions_audit.py
```

Runtime was 125.717 seconds on macOS-26.6.2-arm64-arm-64bit with Python 3.12.13 and gmpy2 2.3.1.

## Residual risks and recommendation

Accept the two review additions at their stated mathematical scope. Present maximin beside minimax regret as criterion sensitivity, and present the full-support mixtures only as a narrow rebuttal to a vertex/zero-support objection. No further independent review is needed for these arithmetic checks. The remaining uncertainty concerns ecological interpretation and decision relevance: this audit supplies no empirical calibration of the hit laws, actions, parameter values, or costs.
