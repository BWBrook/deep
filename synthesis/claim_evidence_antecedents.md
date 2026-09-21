# DEEP: claim, evidence and antecedents

Astra's draft of sections 2 and 3 for Fable's critique, 2026-09-21.
Sections 1 (claim) and 4 (defeat tests and exclusions) remain with Fable;
this is an internal synthesis, not an agreed manuscript or an established
novelty claim. The proposed division is recorded in thread 49, #450–452.

## 2. Evidence and what it establishes

### Common specification and decision loss

Each time step draws theta independently from one fixed G on [0,1]. Patch
hits are conditionally iid Bernoulli(theta); p=E(theta), q=E(theta^2).
Thus p and q fix each complete bivariate **environmental-hit** distribution.
They do not fix bivariate occupancy histories after mortality and recolonisation.
The model starts with all patches occupied and has absorbing regional extinction.
Time horizons below count model steps, without field calibration.

R halves excess shock mortality at one patch; C doubles all colonisation
probabilities. The packages are prescribed, not demonstrated to have equal
costs. Delta(G)=P_R(G)-P_C(G) concerns persistence at the stated finite T.
All laws in a comparison retain the same ecological parameters and action
packages. These are not comparisons of separately fitted models.

For a deterministic, once-only choice from {R,C}, if the *whole-class*
attained extrema of Delta are L and U, the exact minimax regret is

`V = min(max(0,-L), max(0,U))`.

This follows because the loss of choosing R is max(0,-Delta), and the loss
of choosing C is max(0,Delta). If L,U instead describe only a selected
pair of feasible laws, the same expression is the pair's exact value and a
**lower bound** for the full class. Zero pair regret does not establish
zero full-class regret. Randomisation, adaptive management and other actions
are excluded. This distinction makes the previously certified n3/n4 extrema
immediately useful for decision loss without another parameter search.

### Evidence table

Unless stated otherwise, p=.1, q=.03, e0=.02, e1=.9, c=.1, T=50.
Numbers displayed below are rounded; the linked records retain exact arithmetic
or explicitly identify a numerical/aggregate result.

| Result | Parameters and information fixed | Certified conclusion and decision-loss meaning | Evidence and limit |
| --- | --- | --- | --- |
| Moment reduction | Any finite n, fixed ecological parameters, iid fixed G | Kernels depend on G only through m0,...,mn. At n=2, p,q identify these model outputs. Complete hit-count distribution is sufficient for the kernels; it does not identify the full mixing law G. | [n3 derivation](../ideation/b_screen/n3_certificate.md). Mathematical reduction, not a novelty claim or an inference theorem for occupancy records. |
| Three-patch whole-class bounds | n=3; 16 combinations q in {.03,.06}, e0 in {.02,.1}, e1 in {.5,.9}, c in {.02,.1}; p=.1,T=50 | 8 cases favour C throughout, 7 favour R throughout, 1 has both signs. In the default reversal case Delta ranges from -.02341633369 to +.01720185959. Hence whole-class deterministic minimax regret is .01720185959. | [Exact certificates](../ideation/b_screen/n3_certificate_results.json), [derivation](../ideation/b_screen/n3_certificate.md). Whole-class only at these parameter cases; one exact action-switch threshold in m3 in the reversal case. |
| Four-patch whole-class bounds | n=4; same 16-case parameter design | 12 cases favour C throughout, 3 favour R throughout, 1 has both signs. In the default reversal case Delta ranges from -.12052030700 to +.03209273499. Whole-class deterministic minimax regret is .03209273499. | [Exact attained/outer bounds](../ideation/b_screen/n4_boundary_results.json), [audit](../ideation/b_screen/n4_audit.md). All 256 extremum certificates have zero gaps. This is not a theorem that the endpoint laws optimise persistence at arbitrary parameters. |
| Equal triple-hit distributions | n=4, default parameters; additionally m3=.021; m4=.0171 or 1389/70000 | Complete labelled triple-hit laws agree, but Delta is -.00570317978 versus +.00598197568. Pair regret .00570317978. | [Constructive laws and exact reproduction](../ideation/b_screen/n4_audit.md). Feasible witnesses, not whole-class optimisation conditional on m3. |
| Equal four-patch hit distributions | n=5, default parameters; common m3=337/18000,m4=773993/51840000 | Feasible fifth-moment endpoints give Delta=-.00283857955,+.00235757569; pair regret .00235757569. | [Witness](../ideation/b_screen/n5_equal_quadruples.json), [method](../ideation/b_screen/n5_equal_quadruples.md). Same complete four-patch environmental-hit laws. No n5 whole-class bounds. |
| Equal five-patch hit distributions | n=6; default parameters at c=.1 and a second case c=.2; lower moments fixed within each pair | Certified endpoint signs reverse. Pair regret .00012115843 at c=.1 and .00007086475 at c=.2. | [Exact witnesses](../ideation/b_screen/n6_witness_check.json). Existence only. No general-n theorem, universal decay rate or ecological-frequency conclusion. |
| Positive-area ecological-parameter witness | n=4,p=.1,q=.03,e0=.02,T=50; e1 in [.7995,.8005],c in [.1495,.1505] | Exact enclosures have opposite signs everywhere on a rectangle of coordinate area 1e-6. | [Audit and reproducible interval check](../ideation/b_screen/n4_rectangle_audit.md). Pair witness. The separate 907/1024 large-rectangle area report lacks saved leaves from that original run and is not needed for this conclusion. |
| Fixed-e1 background-mortality intervals | n=4,p=.1,q=.03,T=50: (e1,c)=(.8,.15) or (.9,.45) | Reversal throughout e0 in [0,.03] or [.04,.14], respectively. At e0=.1 in the second setting, pair regret .07583193349. | [48 exact leaf certificates](../ideation/b_screen/n4_background_mortality.json), [interpretation](../ideation/b_screen/n4_background_mortality.md). Fixed e1 changes shock amplitude and buffering benefit as e0 changes. Pair reversal sets are not nested; no whole-class geometric conclusion. |
| Constant-amplitude background intervals | n=4,p=.1,q=.03,T=50,c=.3,e1=e0+a | Reversal throughout e0 in [0,.04] for a=.5, [0,.08] for a=.7, [0,.12] for a=.88. The last is the whole admissible e0 interval at that a. | [Leaf certificates](../ideation/b_screen/n4_constant_amplitude_cert.json). Pair existence over these intervals; no certified absence beyond them, amplitude threshold, or universal shrinking-set theorem. |
| Finite-horizon change in ranking and loss | n=4,p=.3,q=41/300,e0=.02,e1=.8,c=.15 | Opposite signs at T20, both C at T50. At the separate e1=.9,c=.3 witness, pair regret falls from .07455511 at T20 to .00002570 at T200 while signs remain opposite. | [Exact anchors](../ideation/b_screen/n4_horizon_audit.json), [audit](../ideation/b_screen/n4_horizon_audit.md). Stable aggregate sign counts do not establish stable decisions or benefits. |
| Uniform long-horizon limit | n=4,e0>0,e1>=e0, same closed model/actions | P_action(T)<=(1-e0^4)^T, so persistence, Delta and full-class two-action minimax regret tend to zero. | [Proof](../ideation/b_screen/n4_horizon_audit.md). Crude structural bound; no rapid-convergence claim, no e0=0 assertion. |

### Constant-amplitude audit and an exact interpretation counterexample

All 12, 47 and 68 saved leaves for amplitudes .5, .7 and .88 have the
required exact signs and contiguously cover their stated targets, with no
unresolved leaves. A separate Fraction implementation reproduced the first
and last leaf bounds in each target. Degree four follows from at most four
affine mortality factors; exact fifth finite differences also vanish.

However, constant amplitude does not imply nested shrinking reversal sets.
At a=.7, p=.1, q=.03, n=4, T=50, exact propagation gives:

| e0 | c | Delta(Glo) | Delta(Ghi) | Endpoint pair reverses? |
| --- | --- | --- | --- | --- |
| 0 | .03 | -.02301329022 | +.01245409236 | Yes |
| 0 | .45 | +.00410151475 | +.10709598505 | No |
| .03 | .03 | -.06592752247 | -.05195853939 | No |
| .03 | .45 | -.00480621616 | +.10616537234 | Yes |

Raising e0 removes one c value from the pair's reversal set and adds the
other, so neither set contains the other. These four cases came from the
saved grid and were then checked by exact rational propagation; no new
parameter screen was added. As elsewhere, agreement of the pair at the
other two points does not establish whole-class agreement. This refutes
monotone pair-set shrinkage, not a claim about the full-class region.
The panel counts also rise and fall at a=.88. The three tested amplitudes
do not locate a critical amplitude near .88, and neither the grid nor inner
existence certificates characterize how an outer boundary disappears.

Reproduce the four exact signs from the repository root:

```python
from fractions import Fraction as F
import sys
sys.path.insert(0, 'ideation/b_screen')
import n3_certificate as N
cases = [(F(0), F(3,100), (-1,1)), (F(0), F(9,20), (1,1)),
         (F(3,100), F(3,100), (-1,-1)), (F(3,100), F(9,20), (-1,1))]
for e0, c, signs in cases:
    ds = []
    for law in N.laws(F(1,10), F(3,100)):
        P = {}
        for action in ('R', 'C'):
            K = N.averaged_kernel(action, 4, e0, e0+F(7,10), c, *law)
            v = [F(0)]*len(K)
            v[7 if action == 'R' else 4] = F(1)
            for _ in range(50):
                v = N.matvec(v, K)
            P[action] = sum(v[1:])
        ds.append(P['R']-P['C'])
    assert all(s*d > 0 for s,d in zip(signs,ds))
    print(e0, c, list(map(float,ds)))
```

### Evidence that supports design but does not certify identification

The [finite-record benchmark](../ideation/b_screen/n4_observability.md)
compares complete hit counts with all-hit indicators for two specified laws.
It illustrates information loss from retaining only all-hit years; it is not
an empirical sampling recommendation or an impossibility claim for occupancy
observations. The n4 p/T and constant-amplitude maps are exploratory grid counts
under one endpoint pair: **not areas, not prevalence, not whole-class identification**.

The corrected [p/T scan](../ideation/b_screen/n4_p_T_scan.json) contains
560/1750 marked cells at p=.3,T200; its minimum marked c is .05. The earlier
561 count and c=.01 endpoint included a cancellation artifact. The stable
routine and regression check are recorded at f8ca3b1. Historical unregenerated
float tables cannot substitute for the exact records above.

## 3. Antecedents and the remaining novelty question

This is a bounded comparison, not a systematic review. A difference from
these papers does not prove that no other paper has already solved the
specific problem. The generic claim that an identification perspective is
new to conservation or ecological risk must not appear in the manuscript.

| Antecedent and reading status | Established contribution relevant here | Defensible distinction for DEEP |
| --- | --- | --- |
| [Ferson & Burgman (1995), *Correlations, dependency bounds and extinction risks*](https://doi.org/10.1016/0006-3207(95)90031-4). Publisher PDF read; [audit](../ideation/b_screen/antecedent_audit.md), pp.103–104. | Ecological risk bounds under unspecified dependence; matching marginal summaries and correlations need not fix the joint distribution or risk. | Our candidate extension fixes complete binary bivariate shock laws and certifies finite-horizon action-ranking ranges under one stationary iid mixing law. Neither dependence bounds nor the broad identification viewpoint is new. |
| [Etienne & Heesterbeek (2001), *Rules of thumb for conservation of metapopulations based on a stochastic winking-patch model*](https://doi.org/10.1086/321986). Utrecht PDF rechecked, equations 3–5, Figure 1C, Appendix A and discussion pp.400–402. | Mortality then survivor-driven recolonisation; a common-Gaussian conditional-independence construction for correlated extinction; extinction-versus-colonisation interventions. The paper explicitly discusses limitations of comparing equal probability changes without effort/cost functions. | Our question holds pairwise shock information fixed while varying compatible higher moments; actions are finite packages and the objective is finite-horizon persistence. Their Gaussian threshold acts directly on extinction, whereas our general mortality is an affine transform of the latent hit probability. We do not overturn their specific derivative-based, eigenvalue/extinction-time ranking. |
| [Ghosh, Sheppard & Reuman (2020), *Tail associations in ecological variables and their impact on extinction risk*](https://doi.org/10.1002/ecs2.3132). Publisher full-text methods/results inspected, especially Figure 1 and equations/model descriptions around Figures 5–7. | Equal marginal distributions and rank correlations can coexist with different tail associations and extinction risks; density dependence changes the effect. The paper examines extinction risk against dispersal. | We fix the entire bivariate binary hit law, eliminating unmatched bivariate tail structure within that observable, and study action rankings and class bounds. Calling this predecessor simply “no actions” obscures its dispersal analyses; the narrower distinction is our prescribed R/C ranking under fixed complete pairwise laws. |
| [Heino et al. (1997), *Synchronous dynamics and rates of extinction in spatially structured populations*](https://doi.org/10.1098/rspb.1997.0069). Publisher PDF read; [audit](../ideation/b_screen/antecedent_audit.md). | Synchrony and extinction studied with nonlinear abundance dynamics, migration and common disturbances. | Their abundance synchrony and repeated-event extinction rate are different observables from our hit laws and absorbing finite-horizon survival. This is context, not evidence that our simplification is generally valid. |
| [McConnachie et al. (2016), *Improving credibility and transparency of conservation impact evaluations through the partial identification approach*](https://doi.org/10.1111/cobi.12610). Publisher PDF abstract/introduction inspected in this synthesis pass; detailed application not audited. | Explicit partial-identification framing in conservation impact evaluation, using bounds and transparent assumptions for missing counterfactual outcomes. | Our unidentified object is a specified environmental dependence law in a mechanistic decision problem, rather than a counterfactual treatment effect. This directly rules out “first identification framing in conservation” rhetoric without establishing an exact mathematical overlap. |

The hit-count/moment reduction, Hausdorff feasibility conditions and Bernstein
bounds are mathematical tools, not inventions claimed here. The
[n5 method note](../ideation/b_screen/n5_equal_quadruples.md) cites the
truncated moment criterion used for its feasibility certificate.

A defensible candidate contribution is therefore **certified partial
identification of a specified finite-horizon conservation-action ranking
under fixed complete bivariate environmental-hit laws**, with sharp small-n
class bounds, stronger matched-subset counterexamples, and an explicit
decision-loss scale. It remains to establish whether that combination is
sufficiently new and scientifically useful; neither the number of examples
nor model-to-model agreement answers that question.

Before a novelty statement is accepted, the remaining bounded search should
cover environmental-dependence ambiguity in stochastic population decision
models, finite-state robust/partial-identification methods with one common
law across time, and decision-value/learning results for higher-order event
information. Methods that choose a new worst-case law at each time or state
solve a different ambiguity problem and must not be conflated with ours.
This synthesis does not need another parameter atlas to expose those gaps.
