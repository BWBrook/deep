# Independent scientific review of DEEP

**Repository:** `BWBrook/deep` (private).  
**Revision actually accessed:** `6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3`.  
**Review date:** 25 September 2026.  
**Reviewer:** GPT-6 Astra Pro.  
**Commission:** Scientific review following `synthesis/external_review_prompt.md`, with particular attention to ecological interest and usefulness.  
**Execution:** Yes. Fresh, independently written implementations were executed locally. They do not import repository modules. No repository changes were made.

## Overall assessment

**There is a sound, potentially publishable theoretical contribution here, but the current evidence is substantially stronger on mathematical validity than on ecological importance. I recommend a bounded, question-led next stage—not another parameter atlas and not an immediate broad conservation claim.**

I found no consequential error in the central four-patch calculations. More importantly, I independently reconstructed their global certification, rather than merely obtaining similar simulation outputs. The scalar and finite-exchangeable action-ranking ranges at all 16 reference parameter combinations reproduce. The equal-triples reversal reproduces, and the mathematical argument for the stated unimodal-class sufficient condition is sound within its declared domain.

My principal concern is not that ecologists will be unable to understand the work. It is that, after understanding it, many may reasonably ask: **“Which conservation inference should I now distrust, which decision should I change, and what should I measure instead?”** The current package answers the first question in a precise but idealised setting. It does not yet answer the latter two convincingly enough for a broad ecological audience.

A useful result from this review sharpens that concern. At the main four-patch reference setting, **minimax regret chooses connectivity, whereas maximising worst-case persistence chooses buffering**, over the same entire exchangeable uncertainty class. Both conclusions have exact supporting calculations. The objective function—not just the dependence assumption—therefore needs to be part of the ecological argument.

The strongest contribution is not “correlation is insufficient,” a new moment problem, a new robust-decision principle, or a universal recommendation about refuges and connectivity. It is a concrete, certified demonstration that **complete pairwise disturbance information can fail to determine a conservation-action ranking, with a measurable cost to that ambiguity, and that apparently innocuous dependence restrictions can supply the missing identifying information.**

## 1. What was reviewed and independently verified

I read the pinned review brief, README, formal model specification, three current Results drafts, their draft audit, the class-comparison table, relevant exact-method source, and supporting audits. I formed the initial model and claim assessment before reading the proposed paper argument and antecedent discussion. I subsequently read the abstract, paper argument, antecedent ledger, unimodality audit and observation benchmark. The source map at the end identifies the relevant paths.

The independent implementation starts from the labelled, 16-state occupancy process. It factors each transition into local extirpation and survivor-driven recolonisation. It checks exact stochasticity, the absorbing state, the prescribed marginal and pairwise hit probabilities, and the buffered-patch state aggregation. It does not reuse the repository's moment-kernel construction.

For certification, I used integer weighted-path polynomial recurrences and exact binomial normalisation to obtain Bernstein coefficients. This shares the mathematical bounding principle with the repository, but not its rational recurrence or source modules. Exact checks at interior points connect the polynomial calculations to direct propagation.

| Claim or component | Independent check in this review | Outcome and scope |
|---|---|---|
| Labelled mortality/recolonisation model | Fresh exact transition construction and 50-step propagation | Consistent with the formal specification; stochastic rows and absorption verified |
| Four-patch finite-exchangeable class E | Enumerated feasible count-law vertices; rebuilt 56 edge polynomials and all 16 bivariate derivative tensors | All 16 global action-ranking ranges certified; not merely vertex sampling |
| Four-patch scalar class S | Rebuilt both curved-boundary polynomials in each of 16 cases; used the independently established positive vertical derivative | All 32 boundary envelopes attain their extrema at endpoints; all 16 reported S ranges reproduced |
| Two-action deterministic and lottery regret | Independent algebra and evaluation from certified extrema | Reported four-patch values reproduced |
| Equal-triples four-patch witness | Exact labelled propagation from its rational moments | Opposite signs reproduced despite identical complete triple-hit distributions |
| Khinchin-unimodal sufficient condition | Checked the analytic moment argument, its polynomial cutoff identity, upper-boundary monotonicity and exact cutoff evaluation | The >0.427-percentage-point C advantage is supported; not a sharp unimodal optimum |
| Mechanism at the reference extremal laws | Independent action and baseline survival trajectories | Supports the main account of changing connectivity benefit |
| New reviewer comparison of decision criteria | Certified separate R and C persistence extrema over E at the two highlighted q values | Maximin persistence selects R at both; at q=.03 this differs from deterministic minimax regret |
| New full-support examples | Exact propagation after mixing each q=.03 E extremiser with 10% of the scalar lower-endpoint law | Both action rankings persist with every labelled hit pattern having positive probability |

### Explicit limits of the review

I did **not** rerun the repository's complete test suite, every historical audit, the three-patch global certificates, all larger-n constructions, or the full asymmetric-law derivation. Their status here is inspected or conditional, not independently replicated. I did not independently recalculate the finite-record discrimination table. The repository's rendered private binary figures could not be reliably inspected through the available retrieval/rendering route; comments about figure hierarchy concern the stated content and numerical results, not a completed visual-layout audit.

The literature check was targeted, not systematic. Publisher full text was accessible for Ghosh et al., Junges et al., Walter et al. and the Khinchin-representation source; the Etienne–Heesterbeek publisher abstract was accessible. Ferson–Burgman's full paper was not independently accessible in this review: the detailed comparison with it therefore remains dependent on the repository's documented reading, rather than a new full-text audit. No claim of exhaustive novelty clearance is warranted.

## 2. Major finding: the result can be made immediately intelligible

The cleanest explanation is already latent in the finite-exchangeable results. Begin with two possible disturbance processes affecting four patches. In both, an individual patch is hit with probability 0.10 per model step, and any specified pair is jointly hit with probability 0.03. Conditional on the number of patches hit, the affected subset is chosen uniformly.

| Number of patches hit in one step | Process A | Process B |
|---:|---:|---:|
| 0 | 78% | 69% |
| 1 | 4% | 28% |
| 2 | 18% | 0% |
| 3 | 0% | 0% |
| 4 | 0% | 3% |

These are event probabilities, not a prescribed sequence of 100 years. Both processes have the same complete pairwise table: neither hit = .83, first only = .07, second only = .07, both = .03. This is stronger than matching an overall correlation coefficient while leaving the rest of a continuous bivariate distribution free.

At the prescribed reference model parameters (`e0=.02, e1=.90, c=.10, T=50`), the independently reproduced persistence probabilities are:

| Outcome | Process A | Process B |
|---|---:|---:|
| No intervention | 35.817% | 22.971% |
| Buffer one patch, R | 50.059% | 38.578% |
| Increase colonisation, C | 67.984% | 32.144% |
| Better prescribed action | C by 17.925 percentage points | R by 6.434 percentage points |

This is the explanation I would put near the beginning of the paper. It lets an ecologist see the missing information before encountering a moment space or an optimisation certificate. These laws are the q=.03 E-class extremisers, not examples chosen independently of the class calculation. [R1], [TABLE], [EA]

There is a real ecological intuition: opportunities for survivors to rebuild the patch network depend on the configuration of disturbance, not merely pairwise coincidence. But the quantitative mechanism should remain faithful to the results. Across these two laws, the gain from C above its own baseline falls from about 32.167 to 9.172 percentage points; the corresponding R gain changes from about 14.242 to 15.606 points. **The action reversal is associated much more strongly with a collapse in connectivity's benefit than with a dramatic increase in buffering's benefit.** An all-patch hit is not necessarily an all-patch extirpation. [R1]

Nor does this phenomenon require exactly zero probabilities in the table. In a reviewer check, replacing each law by 90% of itself plus 10% of the scalar lower-endpoint law preserves p and q and makes every hit pattern possible. The resulting Delta values remain approximately −0.172970 and +0.050918. These are exact feasible examples, not newly optimised bounds. This answers a narrow “vertex pathology” objection; it does not establish that either law is a plausible model of a particular landscape.

**Recommendation:** lead with this event-level comparison. Use the geometry to establish what the example proves, rather than requiring readers to understand the geometry before discovering why the problem matters.

## 3. Major finding: decision relevance needs an explicit conservation objective

### 3.1 Ranking identification and decision choice are different

Let Delta(Q)=P_R(Q)−P_C(Q), with sharp class bounds a and b. A uniformly positive or negative sign identifies the better action under every law. If a<0<b, the best action depends on which compatible law is true.

That does not prevent a decision under a specified criterion. Deterministic minimax regret chooses between

- R's worst-case regret, max(0,−a), and
- C's worst-case regret, max(0,b).

The revised drafts correctly acknowledge that expanding S to E leaves the deterministic minimax action unchanged in all 16 cases, even though it destroys one uniform ranking and worsens two regret guarantees. This should remain prominent, not be relegated to a caveat. A headline suggesting that this expansion necessarily changes the recommended deterministic investment would be misleading. [R1], [R2], [TABLE]

At q=.06, for example, the result is not “a robust manager must abandon buffering.” It is “the scalar restriction made buffering appear uniformly superior; without it, buffering remains the minimax-regret choice but has positive worst-case regret.” That distinction is scientifically valuable, provided the loss of assurance itself is connected to a management tolerance or decision threshold.

### 3.2 A new exact comparison from this review

The repository explicitly says maximin persistence was not evaluated. I evaluated it at the two highlighted four-patch parameter settings, over the entire E class—not merely the four showcased laws. Separate exact vertical and edge certificates establish the extrema of P_R and P_C. Results are saved in `additional_checks.json`.

| q | inf_Q P_R(Q) | inf_Q P_C(Q) | Maximin-persistence action | Deterministic minimax-regret action |
|---:|---:|---:|:---:|:---:|
| .03 | .3857762490 | .3214370903 | R | C |
| .06 | .2273638798 | .1149017905 | R | R |

At q=.03, a manager seeking the highest guaranteed persistence chooses R. A manager seeking to minimise the largest avoidable loss relative to the law-informed choice chooses C. Neither criterion is mathematically wrong. They answer different questions.

This is a material interpretive issue because the ecological motivation naturally evokes protecting a population from extinction, whereas minimax regret protects against having made a comparatively poor choice. The manuscript should explain why the latter criterion is appropriate, and preferably show the former as a sensitivity analysis. A declared prior would yield yet another, prior-weighted decision problem; no prior is supplied by p and q alone.

### 3.3 Randomisation is a benchmark, not yet a management strategy

The once-only lottery calculation is correct under its stated information structure. At q=.03 in E, Pr(R) is about .26413 and worst-case expected regret falls from 6.434 to 4.735 percentage points. The unknown law cannot respond to the realised draw. [R2], [TABLE]

However, a random draw between two complete interventions on one endangered metapopulation is not the same as splitting a budget, protecting some patches while improving links, or adapting after observations. Those physical policies have different transition kernels and may have different costs. For a broad ecological paper, I would retain the lottery as a concise decision-theoretic benchmark, but not let it become the main practical remedy.

## 4. Major finding: the prescribed action menu is adequate for a theorem, insufficient for investment advice

The comparison halves excess hit-related local extirpation probability at one patch and doubles the per-source colonisation parameter throughout the network. These are coherent model interventions. They are not demonstrably equal-cost alternatives, and neither response magnitude is biologically calibrated. [MODEL]

That does **not** invalidate the information counterexample. It does limit how much ecological significance can be inferred from its size. A comparison between two arbitrary points in action space may reveal uncertainty close to a decision boundary without establishing a consequential planning dilemma.

The next analysis should therefore expose a small action-response trade-off: buffering effectiveness against the multiplier or increment in colonisation. Where possible, relate the alternatives through a defensible budget or effort curve. Without data, show a transparent range of cost ratios rather than assert equal cost. The key question is whether a material uncertainty penalty persists over a credible set of feasible choices, not whether a sufficiently fine grid finds more sign changes.

Be precise about the biological meanings. The “mortality” event in this occupancy model is **local population extirpation**, not individual mortality. At e1=.9, an unbuffered occupied patch has a 90% chance of losing its local population when hit. Likewise, doubling c does not double the probability of recolonisation: the latter is `1−(1−c)^m`, with m surviving source patches. “Connectivity” is a label for this particular improvement in recolonisation, not a general representation of corridors, movement, gene flow, or demographic rescue. [MODEL]

A useful decision contract would name the management unit, define a hit, define the interval represented by a model step, specify the terminal conservation objective, and state how the two interventions change that system. A compact hypothetical system with justified assumptions could suffice for a theoretical paper. A practical recommendation requires substantially more.

## 5. Major finding: the observation-to-decision bridge is still missing

### 5.1 The information restriction is meaningful, but narrower than “the data cannot decide”

The fixed information is the within-step environmental-hit distribution for each patch and pair. It is not an entire occupancy time series, abundance synchrony, all joint information in repeated observations, or a fitted dynamic model with every observed variable retained. Coupled occupancy histories can contain information not present in a collection of within-step pairwise environmental summaries. [MODEL], [R3]

Consequently, the result should not be presented as a general impossibility of learning the correct intervention from ecological records. It is a demonstration of information loss or unresolved dependence under a specified observation summary and model class. This is potentially a strong criticism of an actual modelling workflow—but the paper should identify that workflow rather than imply that all ecological models discard higher-order information in the same way.

There is also an important positive message: in most of the reference design, the action ranking is identified despite incomplete knowledge of dependence. Complete reconstruction of the entire disturbance law is not always necessary for a useful conservation decision. The 16-point design does not estimate how often that favourable situation occurs in nature. [TABLE]

### 5.2 Unimodality is assumed information, not automatically measured information

The unimodal result is mathematically convincing. At the reference scalar setting, the Khinchin representation yields the sharp third-moment ceiling

`sup m3 = 689/2000 − 3 sqrt(3)/16 ≈ .01974047358`,

and the convenient rational cutoff `79/4000` gives an upper-boundary Delta below −.00427270637. Thus C wins by more than .427 percentage points throughout that shape class. [R3], [UA]

But the shape restriction is on the latent shared hit probability theta. It is not unimodality of observed rainfall, event severity, local mortality, or the distribution of the number of patches hit. Monotone transformations do not generally preserve density modality, so a unimodal histogram for an environmental covariate is not by itself a justification for this assumption.

At fixed n, the process depends only on finitely many moments. Multiple latent distributions can share those moments. This prevents unique recovery of arbitrary latent shape from those observations, although the moment implications of a shape class can still be falsified. The important distinction is between **learning an observable decision-relevant bound** and **choosing a convenient unobserved shape restriction**.

I would make the operational third-moment bound more prominent than the Khinchin label. In E, an estimable higher-order summary is

`tau3 = E[binom(K,3)] / binom(n,3)`,

where K is the number of hit patches. The corresponding fourth-order quantity is `Pr(K=4)` when n=4. Bounds on such observables can be supplied directly; they do not require claiming that an unobserved distribution is unimodal. The scalar theorem does not automatically extend to E, so a new E-class restriction must be evaluated on its own domain.

### 5.3 The existing observation benchmark is a beginning, not an answer

The repository already compares full hit-count records with records reduced to all-four-hit events for two specified equal-triples laws. It correctly concludes that the full record is more informative, and explicitly distinguishes misclassification probability from persistence regret. [OBS]

Build on that rather than start another unrelated theoretical branch. A decision-focused extension would ask how a feasible observation programme reduces **expected or worst-case persistence loss**, allowing p and q to be estimated rather than known exactly. A confidence region for observable count probabilities could be propagated through the same decision calculation. Under a Bayesian approach, specify the prior and the loss explicitly.

Because the full-count probabilities all change when the top moment changes at fixed lower moments, information is not confined to waiting for all-patch disasters. Conversely, correlated pairs within a single regional event should not be treated as independent replicates when estimating uncertainty. The observation programme needs a sampling unit and an observation model, not simply a longer record length.

## 6. Major finding: test the ecological mechanism, not every imaginable realism extension

The model's simplicity is not a defect in an existence theorem. Four patches are sufficient to expose an information problem, and demanding a continent-scale simulation would obscure the result. The question is which simplification could most alter the ecological lesson.

My priority is **the temporal relation between disturbance and recovery**. The process uses independent environmental draws across steps, surviving donors recolonise within the same step, recently emptied patches are immediately available, and no external immigration is possible. [MODEL]

Those assumptions bear directly on the explanation that connectivity's benefit depends on rebuilding opportunities. For a disturbance that damages habitat or occurs in clusters, a connected patch network may not rebuild on the assumed timetable. Alternatively, rapid recovery in a naturally recolonising system could make the current representation reasonable. The paper has not yet established which interpretation its example represents.

Choose one stress test from a named ecological scenario: temporally persistent disturbance, or a recovery delay that makes newly damaged patches temporarily unavailable. Preserve the relevant marginal and pairwise within-step information when comparing alternative dependence laws. Report the action effects and regret scale, not only whether a reversal still exists. Full global optimisation of an enlarged model need not be a prerequisite: clearly labelled feasible witnesses and bounded sensitivity results can answer a targeted mechanism question.

A local asymmetric perturbation already answers the narrow question of whether perfect symmetry is logically necessary. It does not establish robustness to substantial habitat heterogeneity or spatially realistic disturbance footprints. [LEDGER] Likewise, increasing n changes habitat amount, the number of donors and the fraction of patches buffered; the larger-n examples do not isolate a pure “more unknown moments” effect. [R2]

The long-horizon limit is also correctly handled in the repository. In the closed absorbing model with positive background extirpation probability, persistence and regret eventually tend to zero. That does not mean uncertainty becomes benign: extinction becomes nearly certain. Choose a biologically meaningful horizon rather than use either very long or very short horizons to inflate or suppress the message. [R2]

## 7. Mathematical assessment and important scope boundaries

### 7.1 Central formal structure

The distinction `S ⊂ E ⊂ L` is essential. The scalar conditional-iid representation is more restrictive than symmetry of a finite hit vector, and symmetry is more restrictive than arbitrary labelled laws with the same pairwise tables. The E extension is therefore a substantive improvement over comparing a few scalar families. In particular, the q=.06 example demonstrates that a positive S-class ranking certificate used more than pairwise information. [MODEL], [R1]

The transition matrix is affine in the finite hit law, but finite-horizon persistence is nonlinear because the same law is used throughout repeated matrix multiplication. Optimising each time step or each matrix entry independently would solve a different problem. Similarly, extrema of persistence or Delta over a convex law class need not occur at vertices merely because the one-step kernel is affine. The domain-specific derivative and boundary certificates are doing real work here.

The independent calculation confirms that work at the reference grid: all 16 E derivative arrays have a strict positive vertical sign on a containing rectangle; all 56 E edges are strictly monotone for Delta; scalar curved-boundary Bernstein envelopes attain endpoint extrema. I found no missing interior extremum in those cases.

### 7.2 Information hierarchy and scope

The equal-triples witness is especially clean logically. At `m3=.021`, the two scalar laws with `m4=.0171` and `m4=1389/70000` reproduce Delta values about −.00570318 and +.00598198. Every labelled triple has the same complete hit distribution. This rules out an explanation based merely on unmatched triple-hit information. [N4A]

Its effect size is smaller than the pairwise-only E example. Use the matched-triples construction for logical strength and the E example for decision-loss magnitude; do not conflate their guarantees.

The n≥5 results have a different evidential status. Feasible opposite-sign laws prove ambiguity and lower-bound class regret. They do not identify the unknown extrema, establish a monotone n trend, or measure ecological prevalence. The current Results draft mostly maintains those distinctions. [R2]

Strict opposite signs imply a neighbourhood of reversals by continuity. Consequently, positive-area certificates are useful quantitative robustness statements, but a very small certified rectangle should not be sold as evidence that the effect occupies a broad ecological parameter region. Conversely, failure to find a reversal in a numerical search does not establish its absence. [LEDGER]

### 7.3 Reproducibility and presentation of certainty

The corrected files distinguish attained witnesses, outer bounds, rounded displays, finite samples and unproven extrapolations. That discipline is a genuine strength. The reviewed source also correctly treats the degree-100 E-edge representation as a valid degree elevation, not a different objective. [EA], [N4A]

For publication, give readers a short reproduction route for the current claims, separate from the extensive development history. Freeze dependencies, identify the authoritative scripts and outputs, include a compact regression suite, and fail clearly when certification is unresolved. The goal is not to erase the audit history but to prevent old exploratory outputs from being mistaken for current evidence.

Exact arithmetic eliminates numerical ambiguity within the formal model; it does not justify extra ecological decimal places. Main-text probabilities can be rounded sensibly, with exact fractions and proof receipts in the supplement.

## 8. Novelty and likely audience

The closest literature constrains how the contribution should be pitched.

**Ferson and Burgman (1995)** already addressed extinction-risk bounds under unknown dependence. The repository's reading therefore correctly rules out a claim that dependence ambiguity or its ecological relevance is new. I could not independently inspect their full paper, so I would retain that qualification when assessing exact overlap. [ANT]

**Etienne and Heesterbeek (2001)** directly compared reducing local extinction with increasing colonisation in a stochastic metapopulation model. DEEP does not invent that management contrast. Its distinction is the fixed-information action-ranking question, with different finite packages and objective; it should not be presented as overturning their original rules under their assumptions. [E1]

**Ghosh, Sheppard and Reuman (2020)** demonstrated consequences of tail associations for extinction risk, including spatial models with dispersal. DEEP's sharper distinction is fixing complete binary bivariate hit laws, rather than just overall correlations, and bounding the resulting action-ranking uncertainty. “Earlier studies had no dispersal or management-relevant mechanisms” would be an unfair contrast. [E2]

**Junges et al. (2024)** place fixed-parameter Markov reachability and parameter-region verification within established methodology. **Xu et al. (2021)** already apply minimax regret and policy mixtures to conservation decisions. These relationships narrow the novelty to the ecological information problem and its particular consequences, not generic optimisation or fixed-law uncertainty. [E3], [E4]

As an additional ecological connection, **Walter et al. (2022)** combine theory with evidence from kelp to relate environmental forcing and nonlinear responses to tail-dependent synchrony. This is not an exact antecedent of DEEP's matched-binary-law result. It is an example of how a dependence argument becomes ecologically recognisable: a named driver, a response mechanism and observable predictions. [E5]

I found no exact duplicate in the material inspected. That is not evidence of absence. More importantly, absence of an exact duplicate would not by itself establish ecological importance.

### Will ecologists care?

My judgement differs by audience.

**Theoretical ecologists and quantitative conservation scientists:** there is a credible, interesting paper here. The proof that full pairwise agreement can support opposite interventions, and the separation of information from model restrictions, are worth attention. The certified effect sizes prevent the result from being merely a symbolic counterexample.

**Empirical ecologists:** interest is conditional on identifying a real modelling habit or observation summary at risk, and showing a feasible alternative. Without that, the response may be “we already know covariance is not the whole distribution,” even though that understates the formal result.

**Conservation practitioners:** the present package is not yet an actionable decision method. The action costs, time scale, observational uncertainty and decision criterion remain unspecified. That is an evidence gap, not a failure of communication alone.

In short, readers should not be asked to care because the certificates are exact. They should care because the certificates establish the reliability—or unreliability—of a decision they recognise.

## 9. Recommended next steps, in priority order

### Priority 1 — Fix the decision contract and expose criterion sensitivity

Use the existing four-patch case to state exactly what is being optimised, which action packages are feasible and how horizon and costs are interpreted. Incorporate the maximin/minimax contrast from this review. Add a compact action-efficacy or cost-ratio comparison, rather than an arbitrary expansion of ecological parameter combinations.

**Success criterion:** a reader can identify a defensible conservation objective, understand why the two alternatives belong in the same choice set, and judge whether the unresolved loss is material.

### Priority 2 — Turn “missing higher-order information” into an observation or model-checking task

Choose one illustrative ecological setting and define the hit indicator operationally. Translate higher-order dependence into observable count statistics or a justified disturbance-footprint restriction. Extend the existing observation benchmark to persistence loss and uncertainty in the lower-order quantities. An empirically informed worked example would strengthen broad ecological relevance; a transparent hypothetical example may suffice for a narrower theoretical contribution.

**Success criterion:** the analysis identifies either an attainable observation that improves the decision guarantee or a defensible bound showing that further measurement would not be worth its cost.

### Priority 3 — Test one mechanism-threatening assumption

Select temporal clustering or delayed patch recovery according to the chosen ecological setting. Keep the information comparison controlled and distinguish full-class proofs from attained examples. Do not simultaneously add abundance dynamics, spatial heterogeneity, adaptation, observation error, multiple species and climate trends.

**Success criterion:** the ecological reason for the action reversal survives a meaningful change—or its failure yields a clear limitation that tells readers when the simplified lesson should not be used.

### What not to prioritise

Do not prioritise more n values, many more named hazard families, a broad grid with no sampling measure, a sharper decimal on the unimodal optimum, or a full solution over L merely for completeness. None is necessary to establish the current existential result. Each should earn its place by changing the scientific decision or interpretation.

No new field experiment is essential to validate the conditional theorem. Stronger real-world recommendations, however, require externally justified ecological and intervention assumptions, not another round of internal agreement.

## 10. Suggested manuscript structure and stopping rule

A more effective structure would start with the event table and the two recognisable conservation choices; establish the exact information held constant; show the certified action-ranking interval and its decision consequence; then show how an extra assumption or observation narrows that interval. Put detailed moment geometry, proof machinery, larger-n witnesses and historical repairs in supporting material unless the target audience is explicitly mathematical.

The current abstract has a strong opening but becomes an inventory of results: n=3,4,5,6, scalar versus exchangeable classes, regret, lotteries and unimodality. That makes the reader work to reconstruct the ecological point. I would retain one leading four-patch effect size, the model-assumption reversal, and a statement about what information would resolve the choice. The remaining numbers belong later. [ABSTRACT], [ARGUMENT]

A defensible headline is:

> **Identical pairwise disturbance statistics can support opposite conservation choices.**

The accompanying scientific claim must remain conditional: this has been established for a specified metapopulation process and action menu, with sharp small-system bounds and explicit examples—not as a frequency estimate for real landscapes.

The next stage should have a stopping rule. If the action comparison becomes uninteresting under plausible budgets, if realistic observation uncertainty overwhelms the proposed distinction without suggesting a feasible remedy, or if the chosen recovery mechanism eliminates material decision loss, narrow the paper to a theoretical warning or methods contribution. Do not rescue a broad ecological claim by searching indefinitely for a more spectacular corner case.

Conversely, one well-chosen demonstration that a familiar modelling restriction produces an unjustified recommendation or unacceptable guarantee, together with an observable way to diagnose it, would improve the paper more than dozens of additional exact examples.

## Final recommendation

**Continue, but redirect the effort from proving that ambiguity can exist to explaining when a conservation scientist should care about it.** The central mathematics is not presently the weak link. The strongest next contribution would be a disciplined connection between a conservation objective, a feasible intervention comparison and measurable disturbance information.

I would support development toward a focused theoretical ecology or quantitative conservation paper. I would not yet support broad claims about real-world investment priorities, the prevalence of decision reversals, or a generally superior dependence-aware management strategy.

The risk is not that ecologists “will not get it.” The risk is that they will get the conditional theorem and remain unconvinced that it changes a scientific or management task. That is addressable, but requires some new ecological substance—not only clearer prose.

---

## Source map and reading status

All repository links below are pinned to the reviewed revision. The independent computations and new results are in this review package, not in the repository commit.

[MODEL]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/model_specification.md
[R1]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/results_identification.md
[R2]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/results_decision_loss.md
[R3]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/results_shape_information.md
[TABLE]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/class_identification_table.md
[EA]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/ideation/b_screen/exchangeable_certificate_audit.md
[N4A]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/ideation/b_screen/n4_audit.md
[UA]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/ideation/b_screen/unimodal_audit.md
[OBS]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/ideation/b_screen/n4_observability.md
[LEDGER]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/claim_evidence_antecedents.md
[ANT]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/d1_antecedent_search.md
[ABSTRACT]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/abstract.md
[ARGUMENT]: https://github.com/BWBrook/deep/blob/6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3/synthesis/paper_argument.md
[E1]: https://www.journals.uchicago.edu/doi/10.1086/321986
[E2]: https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.3132
[E3]: https://link.springer.com/article/10.1007/s10703-023-00442-x
[E4]: https://proceedings.mlr.press/v161/xu21a.html
[E5]: https://onlinelibrary.wiley.com/doi/10.1111/ele.13991

| Key | Source and extent of use |
|---|---|
| MODEL | Formal model and information-set specification; read |
| R1–R3 | Current Results drafts; read and central n=4 claims independently checked as described |
| TABLE | All 16 reported scalar/exchangeable ranges; compared with independent calculations |
| EA, N4A, UA | Exchangeable, scalar and unimodal audits; read; independent verification scope above |
| OBS | Finite-record observation benchmark; read, not independently re-executed |
| LEDGER, ANT | Authors' claim/evidence and antecedent analyses; read after forming initial assessment |
| ABSTRACT, ARGUMENT | Current proposed narrative; read |
| E1 | Etienne & Heesterbeek (2001), *The American Naturalist* 158:389–407; publisher abstract and metadata checked, full paper not independently inspected |
| E2 | Ghosh, Sheppard & Reuman (2020), *Ecosphere* 11:e03132; publisher full-text sections on model, dependence and dispersal inspected |
| E3 | Junges et al. (2024), *Parameter synthesis for Markov models: covering the parameter space*; publisher definitions and parameter/reachability discussion inspected |
| E4 | Xu et al. (2021), *Robust reinforcement learning under minimax regret for green security*, PMLR 161; primary proceedings abstract inspected, not a theorem-by-theorem audit |
| E5 | Walter et al. (2022), *Ecology Letters* 25:1189–1201; publisher introduction, theory and empirical methods inspected |

The Khinchin representation was also checked against Ion, Klaassen & van den Heuvel (2023), *TEST* 32:566–601, DOI [10.1007/s11749-022-00844-9](https://link.springer.com/article/10.1007/s11749-022-00844-9), which identifies the representation and its appendix proof. Ferson & Burgman (1995), *Biological Conservation* 73:101–105, DOI `10.1016/0006-3207(95)90031-4`, is retained as an important antecedent with the access limitation stated above.
