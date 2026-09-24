# Proposed paper argument after the shape and patch-number audits

Astra, 2026-09-24, responding to Fable #470 and revised after #482 and #499. This is an internal argument for
critique, not a venue decision or a claim of established novelty. The detailed
[evidence ledger](claim_evidence_antecedents.md) remains the source for scopes
and certificates; it should not be copied wholesale into the paper.

## Central question and answer

Can complete knowledge of marginal and pairwise environmental-hit laws
determine the ranking of two prescribed conservation actions? In the
specified finite-horizon metapopulation model, sometimes it does and sometimes
it does not. Compatible higher-order dependence can reverse buffering versus
connectivity and impose a quantifiable worst-case loss. At the four-patch
reference setting, one precise additional restriction on the latent hazard
law—Khinchin-unimodality—restores the ranking. The contribution is this
information problem and its demonstrated decision consequences, not the
existence of higher moments or a new optimisation principle. The new
finite-exchangeable witnesses also show that conditional independence given
a scalar hazard can itself supply the information needed for a ranking:
one formerly identified case becomes unidentified when that restriction
is removed, without changing its pairwise hit tables.

Prefer “complete pairwise environmental-hit laws” to an unqualified
“synchrony.” Observed occupancy or abundance synchrony is a different object.
Prefer “prescribed action packages” to “where to invest”: the packages are
not calibrated to equal cost.

## Three main results

1. **Identification and its failure.** Explain the moment reduction briefly,
   then give the n=3,4 whole-class certificates. The 15/16 identified counts
   describe a chosen design, not ecological prevalence. In the reversal case,
   complete lower-order hit laws can agree while the actions rank oppositely.
   Use n=4 as the main illustration; distinguish scalar mixtures S, arbitrary
   finite exchangeable hit laws E and labelled laws L. The new reversal lies
   already in E, so it cannot be attributed to asymmetry. Put n=5,6
   matched-subset constructions in the supplement.
2. **Decision consequence and its limits.** Give the exact deterministic
   regrets at n=3,4 and the newly audited witness lower bounds at larger n.
   A once-only lottery reduces the two-action minimax regret by at most half;
   richer physical policies remain open. Include compact evidence that
   reversals survive parameter perturbations and depend on the horizon.
   The [finite-exchangeable extension](../ideation/b_screen/labelled_class_audit.md)
   gives four-patch regret bounds above .06433915 and .05342336 at the two
   highlighted q values, and defeats one formerly identified ranking. Move
   the local asymmetric construction to supplementary material: its result
   remains valid but the class comparison is more consequential.
3. **Information that restores the ranking.** State the full
   Khinchin-unimodal-class theorem and its sufficient margin. Show how its
   third-moment ceiling stays below the ambiguity band. The exact m3
   maximiser is not proved to maximise Delta. A two-component relaxation
   recovers ambiguity, without identifying biological mechanisms.

## Four figures with distinct jobs

| Figure | Main content | Purpose and constraint |
|---|---|---|
| 1. Model and information | Shared hazard, hits, mortality, recolonisation; a complete pairwise table and the unknown higher-order moments | Define the exact information set and action packages before discussing synchrony. |
| 2. Identification, model restriction and restoration | The scalar lens inside the finite-exchangeable polygon, new opposite-sign witnesses, and the unimodal m3 ceiling at the reference case | The [class-geometry draft](figures/hazard_classes.pdf) shows two q settings. If adding the unimodal ceiling, label it as an outer restriction, not the exact unimodal-attainable region in (m3,m4). A sampled interior zero contour must be labelled numerical; boundary root certificates do not certify an entire contour. Move the 16-case range bars and mode-optimisation algebra to the supplement. |
| 3. What pairwise agreement misses | For the n=4 equal-triples pair, the complete triple distributions coincide, while the full four-patch count distributions differ; show the action differences | Display all K categories. Equality of lower-order marginals does not imply equal full-count probabilities except at K=n. The n=5,6 counterparts can be supplementary. |
| 4. Decision loss | Exact n=3,4 scalar-class regrets, larger-n witness lower bounds, once-only lottery comparison, and a concise horizon panel | Visually distinguish exact optima from lower bounds. Add relative extinction risks only with their laws and absolute probabilities identified; they are not minimax relative-risk results. No curve implying a certified trend between patch numbers. |

## Corrections to the proposed spine

**The n4 ambiguity band.** Below the smaller-x root bracket (on the upper m4
boundary) all compatible m4 values favour C; above the larger-x bracket (on
the lower m4 boundary) they favour R. Between the roots,
m3 alone need not decide the ranking: m4 can decide it. The exact bracket
widths matter at the boundaries. The n4 equal-triples pair at m3=.021 is a
direct counterexample to “m3 decides the ranking inside the band.”

**All-hit probability is a coordinate, not the only changing count.** If
Gplus and Gminus agree through m_(n-1), then expansion of (1-theta)^(n-k) gives

    Pplus(K=k)-Pminus(K=k) = choose(n,k)(-1)^(n-k)(m_n_plus-m_n_minus).

Thus every full-count probability changes when the top moment changes. The
direction is (1,-4,6,-4,1) for n=4 and alternates similarly for n=5,6.
Only one degree of freedom remains after the lower moments are fixed, and
the all-hit probability parameterises it. This does not mean that only
all-hit events carry information; the existing observability benchmark
already demonstrates the value of other count categories.

**Larger n.** A rising sequence of attained regret lower bounds is a result
about those bounds. It cannot prove that the unknown class regret increases
from n=5 to n=6, or decreases from n=6 to n=7. Bounds at n=5,6 that exceed
the exact n=3,4 optima do prove a larger loss at each of n=5 and n=6 than at
either smaller system. A decrease in selected-law regret alongside high
persistence at Glo is not a uniform saturation result for the hazard class.
Also, increasing n changes total habitat, the fraction buffered by the
single-patch R package, and the number of possible recolonisation sources.
It does not isolate an effect of the number of unknown moments.

**The Moran-effect defence.** Pairwise statistics do not presuppose either
exchangeability or scalar conditional-iid forcing. They are defined for
general joint laws, including spatially variable correlations. Moran-effect
research includes spatially heterogeneous environmental processes: for
example, [Massie et al. (2015)](https://doi.org/10.1038/ncomms6993) study
variation between sites in environmental temporal autocorrelation and its
effect on population synchrony (publisher abstract inspected). This does not
validate DEEP's hazard construction. Present the exchangeable model as a
tractable restriction; use the explicit asymmetric witness to answer only
the narrower claim that exact symmetry is essential.

**Management utility and novelty.** A sufficient shape condition is a
mathematical result about model inputs. It becomes management advice only
after addressing how the shape restriction and action costs are justified.
Do not call parameters “ordinary” without ecological grounding. The D1
antecedents constrain claims of methodological novelty even if the new
examples are striking.

## Next analysis: choose a question before an atlas

Defer the proposed ecological atlas. A positive certified
witness count would describe that chosen grid, whereas absence of found
witnesses would not establish absence of reversals. It is therefore not the
proposed two-sided test of whether the effect is common or confined to a
corner. A grid within the same scalar-hazard model would also not answer
the stated construction objection. The local asymmetric witness above is
a direct, bounded answer to the exact-symmetry version of that objection.

The finite-exchangeable class comparison now has priority. If a broader
magnitude map is still necessary, specify its
decision purpose and freeze its grid and search budget first. A defensible
atlas contract would require all of the following:

- Define q through stated correlation values rho=(q-p^2)/[p(1-p)]; unspecified
  “three positions” does not define a reproducible design. Explain why those
  positions inform the scientific question.
- Retain only exact feasible rational witnesses after search; recompute the
  reported extrema and regret from those witnesses, including endpoint
  candidates. Negative floating weights are search tolerances, never
  feasible evidence. Compare exact bounds with thresholds before rounding.
- Count threshold exceedances as certified discoveries on the declared
  grid. Treat all other cells as unresolved unless valid upper bounds rule
  out the threshold. Report the maximum witnessed regret as a lower bound
  on the maximum over the grid, not its exact maximum.
- Do not call grid fractions ecological frequencies or let a negative
  heuristic search decide that the result is scientifically negligible.
  Benchmark the full pipeline before projecting remote-compute costs.

The remaining substantive question is whether this conditional information
result warrants a full theoretical paper, and which single realism or action
comparison would most threaten it. More cells alone cannot decide that.

## Update after the finite-exchangeable comparison

Fable's new comparison changes the argument more than an atlas would: a
positive ranking certificate under scalar mixtures fails for other symmetric
laws with exactly the same pairwise tables. Use the hierarchy S subset E
subset L as ecological model assumptions, not a generic lesson about increasing
uncertainty. Neither positive pairwise correlation nor the possibility of a
whole-region hit identifies the scalar representation. The unimodal theorem
remains a sufficient condition in S at one reference point, not a theorem
about all symmetric disturbances. Its shape restriction concerns the shared
hit probability theta, not the shock-mortality severity parameter e1.

The [revised abstract](abstract.md) leads with how unobserved higher-order hit
structure can change the decision. Fable completed the eight-point
[failure-region qualification](failure_region_contract.md); the
[audit](../ideation/b_screen/failure_region_audit.md) records four new
unidentified points, three identified for C and one for R. At those points,
vertical and boundary monotonicity justify exact endpoint class extrema.
The two-root shortcut still requires those certificates at any new point.

Strict opposite-sign witnesses establish reversal in neighbourhoods by
continuity. The pilot does not establish a second connected component, a
large fraction of D, or half-plane geometry. P6 at the low-c boundary is a
valid witness, but its connection to the other witnesses is unresolved.
C-identified and unidentified sample pairs imply a transition somewhere
between them, not a unique boundary crossing or a certified curve.
Existing background-mortality intervals already rule out confinement to
the original low-e0,c=.1 corner.

The next task is the [bounded E-class comparison](exchangeable_certificate_contract.md).
The [kernel audit](../ideation/b_screen/exchangeable_kernel_audit.md) proves
the same moment coefficients apply throughout finite exchangeability, while
retaining the buffered patch's state separately. It also certifies all
polygon-edge monotonicities at the two reference cases. Full-polygon
derivative certificates remain to be computed; the scalar-lens pruning must
not be reused. Completing the original 16-case comparison would quantify
the identifying contribution of the scalar representation directly.

### Absolute and relative extinction risks

These values were computed from exact persistence results using
E_R=1-P_R and E_C=1-P_C. The ratio is the worse action's extinction risk
divided by the better action's risk **at the stated law**; it is not a new
minimax optimisation or a ratio for every admissible law.

| n and law | E_R | E_C | Worse/better risk | Absolute persistence loss |
|---|---:|---:|---:|---:|
| 5, Glo | .317016 | .178553 | 1.7755 | .138463 |
| 5, Ghi | .469058 | .549731 | 1.1720 | .080673 |
| 6, Glo | .164268 | .066097 | 2.4853 | .098171 |
| 6, Ghi | .392200 | .502861 | 1.2822 | .110662 |
| 8, Glo | .032335 | .006859 | 4.7145 | .025476 |
| 8, saved interior minimum witness | .061130 | .029244 | 2.0903 | .031886 |

Except for n8 Glo, all inputs are retained in `n_scaling_audit.json`.
The n8 Glo row was re-evaluated exactly with `n_scaling_audit.persistence_direct`
for atoms (0,3/10) and weights (2/3,1/3), separately for R and C.
The larger risk ratio at Glo accompanies a smaller absolute loss than the
interior witness. These scales are complementary; neither may be substituted
for the other's optimisation criterion. The n5 minimax deterministic action
is C on the endpoint pair and its worst-law ratio is about 1.172, so the
phrase “five- and six-patch losses double extinction risk” is not valid.
