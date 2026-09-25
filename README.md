# DEEP

Dynamics of Eco-Evolutionary Patterns: an experiment in collaborative scientific research led by two
frontier models, Fable (Anthropic) and Astra (OpenAI Codex), with light-touch human input from
Barry W. Brook (University of Tasmania).

Cross-agent correspondence, decisions and the collaboration record are on the shared correspondence
board, project `deep`. This repository holds code, results and manuscript material.

Barry authorised both agents to commit and push when appropriate on 2026-09-20,
with a preference for bundled commits. Coordinate Git mutations on the board
and checkpoint completed, reviewed work together.

## Research topic

Barry approved candidate B on 2026-09-20: "Proceed with B" (correspondence
thread 50, message 373). The project studies what fixed marginal and pairwise
environmental-shock laws identify about finite-horizon metapopulation
persistence and the ranking of buffering versus connectivity, within the
specified conditional-iid regional-hazard model. Scientific discussion now
uses thread 49. Novelty and ecological applicability remain questions to test.

The three- and four-patch certificates under `ideation/b_screen/` cover the
full scalar conditional-iid hazard-law class at 16 fixed ecological parameter combinations
for each system size. For four patches, connectivity wins throughout the
class in 12 cases, buffering in three, and the ranking reverses in one.
This is narrower than arbitrary symmetric hit laws. Full certificates for
that larger finite-exchangeable class identify 14 of the same 16 cases and
leave two unidentified, with sharp regrets of approximately 5.342 and 6.434
percentage points (`exchangeable_certificate_audit.md`). The earlier labelled
witnesses attain those extrema. The [range comparison](synthesis/figures/class_identification.pdf)
and [class geometry](synthesis/figures/hazard_classes.pdf) show how the scalar
assumption contributes information. The deterministic minimax choice remains
unchanged in all 16 cases; its guarantee and optimal once-only lottery change.
These are full finite-exchangeable bounds, not bounds over all labelled laws.
The four-patch audit also gives two laws with identical complete triple-shock
distributions and opposite intervention rankings (`n4_audit.md`). These are
conditional model results, not claims about field systems. Sharp class bounds
for five or more patches remain unresolved. Verified feasible-law pairs give
deterministic regret lower bounds exceeding 8.067 and 9.817 percentage points
at five and six patches; once-only lotteries still have lower bounds exceeding
5.09 and 5.20 points (`n_scaling_audit.md`). These do not establish a monotonic
trend in class regret with patch number. A five-patch existence witness fixes all
complete four-patch shock distributions and still reverses the ranking
(`n5_equal_quadruples.md`). A separate finite-record benchmark compares
information in complete hit counts with all-hit event counts
(`n4_observability.md`); it does not establish field-data requirements.
An exact interval check also establishes a positive-area four-patch reversal
rectangle in the specified mortality/connectivity coordinates
(`n4_rectangle_audit.md`). The larger saved area report remains an aggregate
without leaf-level certificates.
The background-mortality study certifies reversal across an interval from
0.04 to 0.14 at one interior parameter setting, and shows that reversal sets
for the two fixed endpoint laws need not shrink monotonically with background
mortality (`n4_background_mortality.md`). Raising background mortality at fixed
shock mortality also weakens buffering's absolute benefit; this path does
not isolate one biological mechanism.
Exact horizon checks show that similar reversal-grid counts can conceal
changes in both ranking and decision-loss magnitude at fixed parameters
(`n4_horizon_audit.md`). In this closed model with positive background
mortality, persistence and absolute action differences tend to zero at long
horizons; stable signs do not establish stable benefits.
The hazard-family audit finds exact reversals within simple two-regime laws,
including interior supports and continuous two-uniform mixtures, while the
tested beta, probit-normal and logit-normal comparators agree on connectivity
(`hazard_family_audit.md`). This tests identifying restrictions, not ecological
credibility of those families.
An exact shape-class audit proves that connectivity wins for every compatible
Khinchin-unimodal hazard law at the four-patch reversal setting, with an
advantage greater than 0.427 percentage points (`unimodal_audit.md`). This
illustrates additional shape information restoring identification; it does
not establish the empirical credibility of unimodality or biological regimes.
A small explicit asymmetric perturbation also preserves the original complete
pairwise hit laws and opposite action rankings (`asymmetric_pairwise_witness.md`).
This establishes local robustness beyond exchangeability, not a calibrated
spatial model. The proposed paper narrative is in
[`synthesis/paper_argument.md`](synthesis/paper_argument.md).
The [working abstract](synthesis/abstract.md) records the current thesis.
The [failure-region pilot audit](ideation/b_screen/failure_region_audit.md)
checks eight new scalar-class points: four unidentified, three favouring C,
and one favouring R. These points do not establish area or connected
components. The [finite-exchangeable class comparison](synthesis/exchangeable_certificate_contract.md)
is now complete. The [Result 1 draft](synthesis/results_identification.md)
includes an [independent mechanism audit](ideation/b_screen/mechanism_audit.md)
of its four extremising laws and a [comparison figure](synthesis/figures/mechanism_extremisers.pdf).
Buffering wins initially at all four laws; connectivity overtakes it at
horizons 4 and 11 under the two minimisers. Complete-action baseline contrasts
show that the reversal accompanies a large reduction in connectivity's
benefit, with much smaller changes in buffering's. An all-patch environmental
hit need not kill every occupied patch, and recolonisation can occur in the
same step if survivors remain.
An [exact kernel-transfer check](ideation/b_screen/exchangeable_kernel_audit.md)
supports using the same moment coefficients on that larger class; the domain
of its certificates must expand with the class.

## Current reading set

Start with the [model specification](synthesis/model_specification.md), then
[Result 1](synthesis/results_identification.md),
[Result 2](synthesis/results_decision_loss.md) and
[Result 3](synthesis/results_shape_information.md).
The [draft review](synthesis/results_draft_audit.md) records the latest
source checks and corrections. It also distinguishes the added finite
boundary samples from continuous-domain certificates.
Astra and Fable agree this corrected Results 1–3 package is the point for the
first independent GPT-6 Pro review, before further exploration or journal
drafting. The [neutral review brief](synthesis/external_review_prompt.md)
incorporates their agreed edits and uses only this repository for context.
Use it with the immutable commit supplied in the commissioning message.
The review is ready to commission; no external assessment has yet occurred.

## Layout

- `synthesis/` — internal claim, evidence and antecedent draft; the evidence table separates whole-class results from selected-law witnesses.
- `ideation/` — record of the topic-selection phase, including the candidate-B feasibility screen
  (`ideation/b_screen/`: code, patch, log and results, with hashes posted on the board).
- `resources/` — downloaded external files, untracked; every file is recorded in `resources/LEDGER.md`.

## Round-1 constraints (Barry, 2026-09-20)

Work must be theoretical, simulation-based, or built on public-domain online datasets. No new field
data or experiments. Downloads go under `resources/` with a ledger entry. Compute: Mac mini for small
jobs; `deep-server` or `nectar-vm` for larger; `kunanyi` via PBS for memory- or CPU-bound work.
