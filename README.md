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
full admissible hazard-law class at 16 fixed ecological parameter combinations
for each system size. For four patches, connectivity wins throughout the
class in 12 cases, buffering in three, and the ranking reverses in one.
The four-patch audit also gives two laws with identical complete triple-shock
distributions and opposite intervention rankings (`n4_audit.md`). These are
conditional model results, not claims about field systems. Bounds for five
or more patches remain unresolved. A five-patch existence witness fixes all
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

## Layout

- `synthesis/` — internal claim, evidence and antecedent draft; the evidence table separates whole-class results from selected-law witnesses.
- `ideation/` — record of the topic-selection phase, including the candidate-B feasibility screen
  (`ideation/b_screen/`: code, patch, log and results, with hashes posted on the board).
- `resources/` — downloaded external files, untracked; every file is recorded in `resources/LEDGER.md`.

## Round-1 constraints (Barry, 2026-09-20)

Work must be theoretical, simulation-based, or built on public-domain online datasets. No new field
data or experiments. Downloads go under `resources/` with a ledger entry. Compute: Mac mini for small
jobs; `deep-server` or `nectar-vm` for larger; `kunanyi` via PBS for memory- or CPU-bound work.
