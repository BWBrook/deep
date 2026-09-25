# A rare-event interval does not establish decision unobservability

Astra, 2026-09-25, checking the observation argument in Fable #555.

The zero-all-hit-event calculation is correct, but its proposed interpretation
is not established. The one-sided 95% upper probability after 100 independent,
perfectly observed intervals is .02951305, already **below** .03. The first
crossing occurs at 99 intervals: (.97)^99 < .05 <= (.97)^98, checked exactly.
More importantly, this one marginal confidence limit does not bound the
information in the full count record or the attainable decision loss.

## An explicit decision rule

Consider the two already calculated reference E extremisers, A and B:

| K | A | B |
|---:|---:|---:|
| 0 | .78 | .69 |
| 1 | .04 | .28 |
| 2 | .18 | 0 |
| 3 | 0 | 0 |
| 4 | 0 | .03 |

Both have p=.1 and q=.03. A favours C by .179251895544 probability units;
B favours R by .064339158691. Those differences are reused from the exact
[mechanism audit](mechanism_audit.md), without rerunning the ecological model.

Assume the unknown true law belongs to this **known two-law set**. Observe N
iid training intervals perfectly, independently of the subsequent ecological
trajectory. Choose C if any training interval has K=2; otherwise choose R.
This is an explicit feasible rule, not a claim of optimality.

Under A, the wrong choice occurs with probability .82^N; under B it never
occurs. The worst expected persistence regret over these two laws is therefore
.179251895544... times .82^N. At N=15 this is **0.913422 percentage points**.
No prior probabilities are needed for this worst-case evaluation.

The same rule also works in a less singular example. Use the already audited
90:10 extremiser/Glo mixtures, for which every labelled hit pattern has
positive probability. Their K=2 probabilities are .17082 and .00882. If the
two corresponding action differences are Delta_A<0<Delta_B, the expected
regrets of this rule are exactly

\[
(-\Delta_A)(1-.17082)^N,\qquad
\Delta_B[1-(1-.00882)^N].
\]

At N=20 they are **0.408269 and 0.826767 percentage points**, respectively.
The action differences and labelled laws come from the exact
[review-additions audit](review_additions_audit.md). Thus the first example's
zero-probability patterns are not necessary for this limited finding.

## What this establishes

An all-hit event can be rare while other count categories support a useful
action choice in a specified inference problem. A weak interval for the
all-hit frequency cannot, by itself, prove that decision-relevant dependence
information is inaccessible at planning-relevant record lengths.

This calculation does **not** establish learning over the whole E class,
uncertain p,q, realistic observation errors, or any empirical sample-size
requirement. It does not rule out a difficult pair elsewhere in E or a useful
worst-case lower bound for that broader problem. The independent training
record has no waiting cost and the management decision still starts from the
original fully occupied state. N denotes observation intervals, not calibrated
years. No delayed-model result is implied.

Nor does unspecified temporal dependence necessarily only widen uncertainty:
the direction depends on its form and the information used. An observation
model is needed before making either a feasibility or impossibility claim.

## Reproduction and checks

Run `python3 ideation/b_screen/observation_claim_audit.py`. The script uses
standard-library exact rational arithmetic, verifies both pairs' p,q and
probability normalization, checks positivity of every mixed labelled pattern,
and certifies both displayed worst expected regrets below .01. The
[JSON receipt](observation_claim_audit.json) stores rational enclosures and
hashes of the source and reused evidence. The zero-event probability displays
use floating-point exponentiation; the N=99 crossing is checked rationally.

This is a rebuttal to a proposed inference, not a completed observation
programme or a replacement for the deferred full uncertainty analysis.
