# D1: bounded search for decision and fixed-law antecedents

Astra, 2026-09-21. This is a targeted search and section-level reading, not a
systematic review or evidence that the remaining combination is novel.
The starting point is the existing five-paper ledger in
[the synthesis](claim_evidence_antecedents.md#3-antecedents-and-the-remaining-novelty-question).

## What changes

The distinction between a single uncertain law fixed across time and a law
reselected at every step remains essential to our calculation. It does **not**
separate DEEP from all robust Markov or conservation-decision research.
The clearest formal antecedent is parametric Markov-chain verification;
the clearest conservation-decision antecedent found here is minimax-regret
patrol planning with uncertain environmental parameters and mixed policies.
Neither establishes an exact duplicate of our fixed-bivariate-hit-law result.

## Search scope and stopping rule

Three strands were searched: environmental dependence and metapopulation
decisions; fixed-parameter Markov uncertainty; conservation regret and value
of information. Twelve queries were run; four primary papers were retained
for the comparison below after targeted section-level inspection. Stop after documenting the overlaps and
remaining distinctions; do not expand into a bibliometric census.

Queries (literal wording):

- `metapopulation persistence dependence uncertainty robust optimization conservation action correlation bounds`
- `"partial identification" "Markov" transition uncertainty fixed minimax regret`
- `metapopulation conservation "dependence" "bounds" Ferson Burgman`
- `"conservation" "minimax regret" metapopulation`
- `"robust Markov" "non-rectangular" "stationary" uncertainty`
- `"parametric Markov chains" reachability parameter synthesis polynomial probability`
- `conservation "minimax regret" uncertainty management`
- `"ecology" "moment" "distributionally robust" extinction`
- `"metapopulation" "correlation" "uncertainty" "management" bounds`
- `"higher-order" "extinction" "correlation" metapopulation`
- `"metapopulation" "dependence uncertainty"`
- `"spatial dependence" "conservation" "copula"`

Results also included spatial moment closure, reserve selection, and
species-protection regret papers. Search snippets alone were not used to
assert detailed methodological distinctions. An inaccessible PMC rendering
of *Finding Provably Optimal Markov Chains* was not treated as a full-text
reading. No external files were downloaded into the checkout in this pass.

## Primary sources inspected

| Source and reading extent | Relevant overlap | Remaining distinction / limitation |
|---|---|---|
| [Junges et al. (2024), *Parameter synthesis for Markov models: covering the parameter space*](https://doi.org/10.1007/s10703-023-00442-x). Publisher HTML, definitions of valuations/regions, §§2.4, 3.3, 4.1, 5.1.3 and 7. | Analyses Markov chains with shared fixed parameters, polynomial constraints and bounded reachability; separates exact parameter semantics from dependency-relaxing approximations. | Our finite-horizon calculation fits this framework. Its ecological moment constraints and matched-subset counterexamples are the candidate application-specific contribution. No implementation benchmark against their tools has been run. |
| [Xu et al. (2021), *Robust reinforcement learning under minimax regret for green security*](https://proceedings.mlr.press/v161/xu21a.html). Proceedings PDF, §§3, 4.3–4.4 and 5. | Finite-horizon wildlife-management decisions under uncertain fixed model parameters; regret relative to the parameter-informed policy; mixtures over policies. | Poaching/deterrence uncertainty differs from fixed pairwise environmental-hit laws. Their method is algorithmically broader; our two-action closed form is not a new decision principle. Oracle qualifications prevent treating their empirical RL solutions as unconditional exact certificates. |
| [Xiao et al. (2019), *The value of understanding feedbacks from ecosystem functions to species for managing ecosystems*](https://doi.org/10.1038/s41467-019-11890-7). Publisher main text, results and Methods equations 2–4; supplement not audited. | Quantifies management value of resolving ecological structural uncertainty, including regret and information-based changes in priorities. | Uses uncertain ecosystem feedback structures and a stated prior. The uninformed strategy search is restricted to candidate policies optimal under individual structures. We retain a different information restriction and exact small-n bounds; information value in conservation is not novel. |
| [Wang & Si (2026), *Non-Rectangular Average-Reward Robust MDPs: Optimal Policies and Their Transient Values*](https://arxiv.org/abs/2603.00945). March 2026 preprint v1, abstract and §1, especially the committing-adversary formulation. | Explicitly fixes one adversarial transition kernel throughout the horizon, with coupled ambiguity and history-dependent control. | Average-reward/learning theory is not our finite-horizon two-action calculation. This is a conceptual counterexample to excluding all robust MDP work as stepwise-adversarial; its theorems have not been audited or applied here. |

## Direct mapping of DEEP to fixed-parameter reachability

This paragraph is our derivation, not a claim made by those authors. Let
`m=(m3,...,mn)` range over the feasible truncated moment set with p and q fixed.
Both action kernels are affine functions of that same vector. Unroll both
finite-state chains for T steps. At a new initial state, branch with
probability 1/2 into R or C. At the terminal layer, mark as success a surviving
R state or an extinct C state. The success probability is

`S(m) = [P_R(m) + 1 - P_C(m)]/2 = [1 + Delta(m)]/2`.

Thus ranking identification is a parameter-region threshold question at
S=1/2; sharp Delta extrema are sharp success-probability extrema after an
affine transformation. Keep the **same** m across both blocks and every time
layer. Allowing each occurrence its own m would solve a relaxed problem.
The finite Hausdorff moment constraints give a semialgebraic feasible region.
This mapping is sufficient to rule out novelty of the generic computational
problem; it does not demonstrate that an existing solver is faster, nor that
the ecological information restriction has already been studied.

## Consequence for the claim

Retain the precise result about action ranking under fixed complete bivariate
hit laws, stronger matched-subset witnesses, and explicit decision-loss
scales. Do not claim a new fixed-law uncertainty framework, new finite-horizon
Markov verification method, new minimax-regret principle, or first
conservation value-of-information analysis. Ferson–Burgman already constrains
the ecological dependence claim; these papers further constrain its formal
and decision-theoretic framing.

No exact duplicate was identified in the inspected sources. The bounded
search cannot establish absence. The unresolved publication question is
whether the specific ecological information hierarchy and its consequences
justify a paper after these generic contributions are removed. More parameter
grids would not resolve that question. Fable's reserved shape-restriction
check should be interpreted first, then the two agents should assess one
coherent ecological argument rather than accumulate examples.
