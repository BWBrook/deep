# Model and decision problem

Repository specification for reading the draft results without the private
correspondence board. Prepared 2026-09-25 from the existing implementations
and derivations; it introduces no new model or analysis. The research question
is whether specified environmental-hit information determines the persistence
ranking of two prescribed interventions in this model.

## State, environment and within-step events

There are n labelled patches. A_t is the set occupied after model step t,
and A_0 contains every patch. Draw a binary environmental-hit vector
X_t=(X_1t,...,X_nt) independently at each step from one distribution Q.
The unknown Q is fixed over the whole horizon, and the interventions do not
change it. Optimisation ranges over a single fixed Q for each trajectory;
Q cannot change with time, occupancy state, preceding events or the realised
lottery draw.

Let 0 <= e0 <= e1 <= 1 and s=e1-e0. Conditional on the hit vector, mortality
events at occupied patches are independent. Under the baseline and action C,
the mortality probability at occupied patch i is

    d_i(X_t) = e0 + s X_it.

Under R, a single named patch (label 0) instead has mortality probability
e0 + (s/2) X_0t; all other patches retain the baseline rule. R thus halves
the excess mortality caused by a hit at that patch. It does not halve
background mortality or the patch's total mortality in a hit step.

Let U be the surviving occupied set and m=|U|. Each patch outside U is
recolonised independently with probability

    r_a(m) = 1 - (1 - c_a)^m,
    c_0 = c_R = c,    c_C = 2c.

The studied c values lie in [0,1/2]. Donors are the m survivors, so newly
colonised patches do not become additional donors during the same step.
Vacant patches include those emptied in that step as well as those empty
before it. If m=0, recolonisation cannot begin: there is no external immigration
and regional extinction is absorbing. A hit is an exposure, not a certain
death. Patches may die on a step with no hits because e0 can be positive.

Both actions apply for the complete horizon. C doubles the per-source
parameter c, not the complete recolonisation probability r_a(m). The two
packages have not been calibrated to equal cost or effort. Steps are not
calibrated to years or to a particular field system.

For old occupied set A, final set B and a possible survivor set U contained
in both A and B, the conditional transition contribution is

    product_(i in U) (1-d_i)
      * product_(i in A\U) d_i
      * r_a(|U|)^(|B\U|)
      * (1-r_a(|U|))^(n-|B|).

Sum over U and average over Q to obtain the labelled transition matrix K_a(Q),
where a denotes the action. The implementation
in [labelled_class_audit.py](../ideation/b_screen/labelled_class_audit.py)
uses this survivor-subset formula. The separate nested enumeration in
[labelled_class_probe.py](../ideation/b_screen/labelled_class_probe.py)
provides a comparator at n=4. The latter also contains exploratory searches;
their outputs are not all class certificates.

## What the information fixes

For each labelled patch and each distinct pair, the supplied information is

    Pr(X_i=1)=p,    Pr(X_i=X_j=1)=q.

It fixes the complete pairwise hit table:

| X_i, X_j | Probability |
|---|---:|
| 0,0 | 1-2p+q |
| 1,0 | p-q |
| 0,1 | p-q |
| 1,1 | q |

These are environmental-hit probabilities, not joint occupancy distributions
after mortality and recolonisation. p and q are treated as exactly known;
the main calculations do not estimate them from finite observations.

Three dependence classes have distinct meanings:

- **S_n, scalar mixtures:** draw theta independently each step from a single
  fixed distribution G on [0,1], then draw patch hits conditionally iid
  Bernoulli(theta). Here p=E(theta), q=E(theta^2), and the probability that
  any specified k patches are all hit is m_k=E(theta^k).
- **E_n, finite exchangeable laws:** any symmetric distribution of the n hit
  indicators satisfying the same pairwise tables. It is specified by hit-count
  probabilities pi_k=Pr(sum_i X_i=k); every labelled pattern with k hits has
  probability pi_k/choose(n,k). A scalar-mixture representation is not required.
- **L_n, labelled laws:** any distribution on the 2^n labelled hit patterns
  satisfying those marginal and pairwise constraints, without exchangeability.

Thus S_n is contained in E_n, which is contained in L_n. This inclusion is
within a fixed n and action menu; changing n also changes the ecological
system. The main small-n class calculations concern S_3, S_4 and E_4.
The broader labelled constructions and the n>=5 calculations have narrower
witness scopes, as specified in their receipts.

For n=4 in E, write x=tau_3 and y=tau_4 for specified triple- and quadruple-hit
probabilities. Then

    pi = (1-4p+6q-4x+y,
          4p-12q+12x-4y,
          6q-12x+6y,
          4x-4y,
          y).

Nonnegative pi defines the E polygon. The scalar-mixture constraints further
restrict this polygon to a truncated Hausdorff moment region; see
[n4_boundary.py](../ideation/b_screen/n4_boundary.py) and
[the class comparison](class_identification_table.md). S identifies transition
matrices through moments up to order n, not the full distribution G.

Two frequently used scalar laws are called Glo and Ghi. In the nondegenerate
cases studied here they are

    Glo = (1-p^2/q) delta_0 + (p^2/q) delta_(q/p),
    Ghi = (1-w) delta_a + w delta_1,
    a = (p-q)/(1-p),    w = (q-p^2)/(1-2p+q).

They attain the lower and upper feasible third moments. Their status as
extremisers of a persistence objective requires a separate proof for that
objective and parameter setting. Their convex mixtures cover the relevant
third-moment interval at n=3; they do not parameterise the full moment region
at larger n.

The shape restriction in Result 3 is a subclass of S: G is Khinchin-unimodal
if theta=m+UV with U uniform on [0,1], independent of V supported on
[-m,1-m], for some mode m. This definition includes boundary modes and an
atom at the mode. It is a restriction on a mixing distribution, not on
mortality severity e1 or directly on the finite hit-count probabilities.

## Persistence, ranking and regret

For action a, let h indicate nonempty states and mu put all mass on full
occupancy. The outcome is finite-horizon regional persistence,

    P_a(Q,T) = mu K_a(Q)^T h,
    Delta(Q,T) = P_R(Q,T) - P_C(Q,T).

With absorbing extinction this also equals survival through T. Summing
nonempty-state mass directly avoids floating-point subtraction from extinction
probabilities close to one; the exact rational forms are algebraically equal.
Symmetry
permits a count chain for C and the baseline and a buffered-patch/ordinary-count
chain for R; common labelled chains are used where needed to compare actions
step by step.

For a specified ambiguity class H, a uniform action ranking means that Delta
has one sign throughout H. It need not make the persistence probability or
Delta itself point identified. If the extrema are a=min_H Delta and
b=max_H Delta, deterministic minimax regret over the R/C menu is

    min(max(0,-a), max(0,b)).

Regret compares a chosen package with the better of the same two packages
under the true fixed law. For a<0<b, choosing R with probability alpha once
before the horizon has worst expected regret max((1-alpha)b,-alpha a), with
optimum alpha=b/(b-a). The unknown law is fixed independently of the realised
lottery draw. There is no switching, adaptation, spatial allocation search,
or combined physical intervention in this policy class. The baseline P_0
is used for explanatory contrasts and is not part of this optimised menu.

## Parameter and evidence scopes

The original small-n design fixes p=.1 and T=50 and crosses q in {.03,.06},
e0 in {.02,.10}, e1 in {.50,.90}, and c in {.02,.10}: 16 settings for each
analysed small-n class. The reference case is q=.03,e0=.02,e1=.9,c=.1.
Other points, witnesses and horizons have their own stated scopes. These
selected settings are not an ecological sampling design.

[Result 1](results_identification.md), [Result 2](results_decision_loss.md)
and [Result 3](results_shape_information.md) state the current claims.
[The evidence ledger](claim_evidence_antecedents.md) links them to exact
receipts, proofs, exploratory calculations and antecedent notes. Its labels
and the authors' use of “certificate” or “audit” describe intended evidence
types; a reviewer must assess the domain, logic and implementation supporting
each claim.
