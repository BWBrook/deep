# Result 1 draft: pairwise information and dependence assumptions

Astra, 2026-09-25. Internal manuscript text, revised after Fable #526 and
the independent four-law mechanism audit.
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
the per-source colonisation probability c. In S, C had higher persistence
throughout the class in 12 cases, R in three, and neither action ranked
uniformly first in one. Expanding to E
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
criterion. This conclusion concerns minimax regret; maximin persistence and
decisions averaged over a prior on hit laws were not evaluated. Once-only
randomisation responded to the larger uncertainty class:
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

## How the same pairwise information admits opposite rankings

One exact short-horizon identity separates direct mortality protection from
colonisation's effects on future persistence. Starting with all four patches
occupied, one-step regional persistence is the probability that at least one patch
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

We examined the four already certified E extremisers at e0=.02, e1=.9 and
c=.1. Within each fixed-q pair, both laws have the same expected hit count
4p=.4 and expected number of co-hit pairs 6q. The Delta-minimising laws put
co-hits into two-patch events at q=.03, or two- and three-patch events at q=.06.
The maximising laws put all co-hits into four-patch events with probability q,
alongside single-patch hits and steps without environmental hits. This
redistribution changes the full count law, while preserving every pairwise
hit table. The scalar representation excludes these extreme combinations of
triple and quadruple hits.

An all-patch hit does not necessarily cause regional extinction. From full
occupancy, conditional on all four patches being hit, the probability of
surviving mortality is .3439 under C and .66466 under R. If any patch survives,
colonisation can refill vacancies created in that same step. It cannot rescue
the region after its last occupied patch dies. Its absence from Delta_1
therefore concerns the regional persistence criterion, not the timing of
recolonisation.

Exact trajectories show R's initial advantage changing to a C advantage at
horizon 4 under the q=.03 minimiser and horizon 11 under the q=.06 minimiser.
C remains ahead through horizon 50 in those two cases. R stays ahead at
every horizon 1 through 50 under the two maximisers. These are four fixed-law
trajectories, not a sequence of class optima reselected at each horizon
([mechanism figure](figures/mechanism_extremisers.pdf)).

The unmodified baseline helps distinguish the actions' accumulated effects.
At q=.03, moving from the Delta-minimising to the Delta-maximising law reduces
C's 50-step persistence gain above baseline from .32167 to .09172; R's gain
changes from .14242 to .15606. At q=.06, C's gain falls from .19668 to .02416,
while R's changes from .14325 to .13662. Thus the reversal accompanies a much
smaller benefit from enhanced recolonisation; it need not accompany a larger
benefit from buffering. Both actions improve persistence above baseline at
each law. These direct comparisons include the occupancy histories induced
by each complete action package.

To locate the contributions by environmental hit count, we also used a
telescoping decomposition on the common labelled 16-state chain. Put R, C
and the baseline into transition matrices K_R,K_C,K_0; let mu start with all
patches occupied and h indicate nonempty states. Then

    Delta_T = sum_(t=0)^(T-1) mu K_R^t (K_R-K_C) K_C^(T-1-t) h

splits into a buffering allocation B and a colonisation allocation Gamma by
writing K_R-K_C=(K_R-K_0)-(K_C-K_0). Each term is a model-defined one-step intervention
contrast, with R before that step and C after it. It can be further summed
by the current step's hit count using the count-conditioned kernels. The
reverse-order identity, with mu K_C^t and future K_R^(T-1-t), gives the same
total but a different allocation. We report both to expose order dependence;
neither allocation is a unique causal mediation decomposition. Nor do these
hybrid totals equal the complete-action gains P_R-P_0 and P_C-P_0 above.
A prefix of the fixed-T sum is an accounting subtotal, not Delta at a shorter
horizon.

Both allocations reproduce each certified Delta_50 exactly. Across the four
laws and two orderings, approximately 80–94% of B is attributed to the largest
hit count admitted by that law. Steps without environmental hits account for
approximately 46–72% of Gamma; hit steps make substantial contributions too.
Under both orderings, the change from each minimiser to its maximiser is
associated mainly with a smaller Gamma. Detailed count attribution is less
stable: the single-hit net contribution at q=.06's maximiser is positive in
one order and negative in the other. The [audit](../ideation/b_screen/mechanism_audit.md)
reports both orders, exact checks and the distinction between hits and deaths.

These four contrasts explain how dependence left open by pairwise tables can
change the balance between retaining survivors and rebuilding occupancy.
They do not establish a universal division in which buffering acts only
against all-patch hits or colonisation benefits only quiet steps. The model
and the prescribed action costs still require ecological justification before
these conditional comparisons can inform an applied recommendation.
