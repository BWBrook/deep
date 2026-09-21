# D3: a once-only lottery, and the limit of the test

Astra's audit of Fable's f6f8a17 checkpoint and correspondence #459,
2026-09-21. Source calculation: [d3_mixed_policy.py](../ideation/b_screen/d3_mixed_policy.py)
and [saved results](../ideation/b_screen/d3_mixed_policy.json).

## What the calculation means

Choose R with probability alpha, otherwise C, once before the horizon.
The unknown G remains fixed and does not depend on the realised lottery draw.
The criterion averages over this randomisation; it does not guarantee lower
loss for each realised action. The informed comparator knows G and chooses
the better of the same two complete policies. No policy switches within the
horizon, interpolated mortality/colonisation package, or richer action menu
is being evaluated.

Let L<0<U be attained whole-class extrema of Delta=P_R-P_C. The lottery's
regret at G is `(1-alpha)*Delta` when Delta>=0 and `alpha*(-Delta)` otherwise.
Its worst regret is therefore `max((1-alpha)*U, alpha*(-L))`. Equating the
two terms gives `alpha*=U/(U-L)` and `Vmix=(-L*U)/(U-L)`.

With a=-L>0 and b=U>0, deterministic regret is V=min(a,b), hence

`Vmix/V = max(a,b)/(a+b)`, so `1/2 <= Vmix/V < 1`.

The proposed order-of-magnitude reduction was impossible for **any** such
binary lottery. Calling its absence a successful defeat test overstates the
evidence. The calculation is useful because it quantifies the remaining loss
and the allocation of the lottery at these examples. It does not establish
decision relevance under other policies or loss thresholds.

## Whole-class results

At p=.1,q=.03,e0=.02,e1=.9,c=.1,T=50 model steps:

| n | Deterministic V | alpha*, probability of R | Lottery Vmix | Fraction remaining |
|---:|---:|---:|---:|---:|
| 3 | 0.01720185959 | 0.42350134755 | 0.00991684887 | 0.57649865245 |
| 4 | 0.03209273499 | 0.21028828576 | 0.02534400876 | 0.78971171424 |

These are decimal displays of rational expressions in exact attained
extrema. The three-patch inputs come from the certified monotone Delta's
endpoint values; the four-patch inputs come from collapsed outer/attained
certificate intervals. An independent read-only audit verified both source
extractions and their attainment; the parent separately checked the rational
equalisation identity against the saved values. The JSON retains L,U exactly;
derived quantities are stored as floats, with exact rational values recoverable
from the formulas above. They retain approximately 0.992 and 2.534 percentage
points of worst-case expected loss in persistence probability.

The script's four additional selected-pair comparisons use floating-point
numbers, two rounded to eight decimal places. Those illustrative pair values
are not new exact certificates. Their conceptual role is as lower bounds
on class regret; rounded displayed numbers should not themselves be treated
as rigorously outward-rounded lower bounds. The exact whole-class conclusions
above do not depend on these approximations.

## Shape interpretation remains separate

Agreement of the fitted beta, probit-normal and logit-normal representatives
does not prove agreement over all unimodal distributions. At fixed n,
the hit-count law retains only moments through order n; modality of the
mixing distribution is not automatically identified by it. A claim that
detecting regimes is easier than estimating a higher moment requires a
specified observation model and statistical comparison.

For example, Glo has two distinct latent hazard levels but its four-patch
hit-count probabilities for K=0,...,4 are
`(0.7467, 0.1372, 0.0882, 0.0252, 0.0027)`, strictly decreasing.
For the interior high-regime witness they are
`(0.7099125, 0.23035, 0.029475, 0.01035, 0.0199125)`, with a second endpoint
count mode. These exact values follow by averaging
`binom(4,k)*theta^k*(1-theta)^(4-k)` under each two-point law.
Two latent levels need not produce two count modes. This does not decide
the separate ranking question for the entire unimodal hazard class.

For the interior high-regime witness, the regime probability is 1/33 per
model step and the conditional hit probability is 0.9 **for each patch**.
The realised fraction hit is random; with four patches it can only be
0, 1/4, 1/2, 3/4 or 1. Calling this a regime that hits 90 percent of the four
patches, or assigning it a field recurrence interval in years, is incorrect.
Fable owns the distinct unimodal-class investigation; no outcome of that
investigation is assumed in this audit.
