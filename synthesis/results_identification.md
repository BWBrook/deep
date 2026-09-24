# Result 1 draft: pairwise information and dependence assumptions

Astra, 2026-09-24. Internal manuscript text for Fable's substantive review.
Numbers derive from the [certified comparison table](class_identification_table.md)
and [audit](../ideation/b_screen/exchangeable_certificate_audit.md).

Complete pairwise environmental-hit laws did not always determine which of
the prescribed actions had higher persistence. With four patches, the
pairwise tables fix the marginal hit probability p and the joint-hit
probability q for every pair. They leave two further coordinates under
exchangeability: the probability that a specified triple is hit, tau_3, and
the probability that all four patches are hit, tau_4. The finite-exchangeable
class E contains every symmetric four-patch hit law satisfying those tables.
The scalar-hazard class S additionally requires that patch hits be independent
conditional on a shared hit probability. This extra representation restricts
the feasible (tau_3,tau_4) combinations to a proper subset of E, despite the
identical pairwise information supplied to both classes.

We compared the two classes at 16 prescribed combinations of q, background
mortality, shock mortality and colonisation, holding n=4, p=.1 and the horizon
at 50 model steps. R halves excess shock mortality at one patch; C doubles
colonisation. In S, C had higher persistence throughout the class in 12 cases,
R in three, and neither action ranked uniformly first in one. Expanding to E
retained the 12 uniform C rankings and two uniform R rankings, and left two
rankings unidentified. These are counts in a selected parameter design, not
estimates of how frequently each outcome occurs in ecological systems.

The additional unidentified case had q=.06, e0=.02, e1=.9 and c=.1. Under S,
the exact range of Delta=P_R-P_C was approximately [.03101,.11044], uniformly
favouring R. Under E it widened to [-.05342,.11246], admitting opposite
rankings while retaining every pairwise hit table. Thus the uniform-ranking
conclusion at this case used information supplied by the scalar conditional-iid
representation. It did not follow from pairwise environmental-hit information
and symmetry alone. The smallest surviving R advantage across the original
design was at q=.03,e0=.02,e1=.9,c=.02: the class minimum decreased from
approximately .00382 in S to .00079038 in E but remained strictly positive.

Expanding the class also increased the decision loss possible at the reference
case q=.03,e0=.02,e1=.9,c=.1. Delta ranged from approximately [-.12052,.03209]
in S and [-.17925,.06434] in E. The exact deterministic minimax regret rose
from approximately 3.20927 to 6.43392 percentage points of persistence, a
factor of 2.00479. At the additional unidentified case, deterministic minimax
regret increased from zero to 5.34234 points. These are sharp class values
for the two prescribed actions, rather than lower bounds from selected laws.

The best deterministic minimax action nevertheless remained unchanged between
S and E in every design case. It was C at the reference case and R at the
additional unidentified case. Loss of a uniform ranking therefore reduced the
guarantee without changing the best deterministic decision under this loss
criterion. Once-only randomisation responded to the larger uncertainty class:
at q=.06 the optimal probability of selecting R fell from one in S to
approximately .67795 in E, giving E-class regret of 3.62184 points. At the
reference q=.03 case, the corresponding E-class lottery regret was 4.73454
points. Randomisation is over complete action packages chosen before the
shock sequence, not adaptive switching during it.

The E extrema were attained at polygon vertices in all 16 cases. Positive
vertical derivatives throughout E and strict monotonicity on each of its
edges certify this statement case by case. Each extremising hit-count law
has at most three nonzero count probabilities, consistent with the three
linear constraints fixing total probability, mean count and second factorial
moment. This does not imply that vertex laws optimise persistence differences
for arbitrary ecological parameters. Nor do the E certificates establish
uniform rankings for the larger class of asymmetric labelled hit laws.

## Mechanistic anchor for the next review

One exact short-horizon identity separates direct mortality protection from
the delayed benefit of colonisation. Starting with all four patches occupied,
one-step regional persistence is the probability that at least one patch
survives mortality. Colonisation cannot restart an extinct region, and is
unnecessary to meet the nonextinction criterion when a survivor remains.
Thus changing c has no effect on this first-step persistence probability.
Write s=e1-e0. Buffering halves the hit-related mortality increment at one
patch, giving

    Delta_1 = (s/2) E[X_1 product_(j=2)^4 (e0+s X_j)]
            = (s/2) [e0^3 p + 3 e0^2 s q + 3 e0 s^2 tau_3 + s^3 tau_4].

This identity holds for every exchangeable hit law. Its terms are nonnegative,
and at the positive-background settings studied here Delta_1 is strictly
positive. The identity was also checked exactly from the transition matrices
at each of the five fixed-count laws and all eight original ecological
parameter combinations (40 checks). At q=.03's two E extremisers,
Delta_1 is respectively 2233/156250000 and 1503601/156250000, both positive;
their 50-step differences are approximately -.17925 and +.06434.

The observed longer-horizon C preference therefore requires benefits that
arise through subsequent occupancy dynamics. The formula also shows exactly
how higher-order hits affect R's immediate protection advantage at fixed
p,q. It does not decompose the full 50-step difference or prove that the same
monotonic dependence holds at other horizons. Explaining that accumulated
trade-off is a substantive next step, beyond describing where the polynomial
changes sign.

## Questions for substantive review

The mathematical distinction and loss scale are now supported. The remaining
scientific argument should explain why the prescribed action packages respond
differently to higher-order hit structure, using the model's mortality and
recolonisation mechanisms. A causal explanation requires more than naming
triples or all-patch hits: moving between the extreme laws changes the entire
hit-count distribution. The explanation should also acknowledge the unchanged
deterministic minimax choices and specify what additional observation or
dependence assumption would alter a defensible decision. No additional
parameter grid is needed to review this text.

A bounded route to that explanation uses the four existing E extremisers,
without changing ecological parameters or searching for further laws. Put
R, C and the unmodified baseline into the same labelled state space, with
transition matrices K_R,K_C,K_0. Let mu start with all patches occupied and
h indicate nonempty states. The telescoping identity

    Delta_T = sum_(t=0)^(T-1) mu K_R^t (K_R-K_C) K_C^(T-1-t) h

splits into a buffering gain and a colonisation gain by writing
K_R-K_C=(K_R-K_0)-(K_C-K_0). Each term is a model-defined one-step intervention
contrast, with R before that step and C after it. It can be further summed
by the current step's hit count using the count-conditioned kernels. The
reverse-order identity, with mu K_C^t and future K_R^(T-1-t), gives the same
total but a different allocation. Report both to expose order dependence;
neither allocation is a unique causal mediation decomposition.

Fable's next bounded task is to challenge this identity and interpretation,
then evaluate the two allocations for these four laws at T=50 if they pass.
Verify that each sums exactly to the already certified Delta, and explain
which mortality and recolonisation contributions distinguish the opposing
rankings. Keep this diagnostic to those laws; no extra grid or new policy
optimisation is part of it.
