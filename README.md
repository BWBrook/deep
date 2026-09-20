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

The three-patch certificate under `ideation/b_screen/` covers the full
admissible hazard-law class at 16 fixed ecological parameter combinations.
The four-patch exact counterexample disproves coordinate-wise monotonicity of
the intervention difference; it does not establish full-class bounds.
Larger-system bounds are under methodological review.

## Layout

- `ideation/` — record of the topic-selection phase, including the candidate-B feasibility screen
  (`ideation/b_screen/`: code, patch, log and results, with hashes posted on the board).
- `resources/` — downloaded external files, untracked; every file is recorded in `resources/LEDGER.md`.

## Round-1 constraints (Barry, 2026-09-20)

Work must be theoretical, simulation-based, or built on public-domain online datasets. No new field
data or experiments. Downloads go under `resources/` with a ledger entry. Compute: Mac mini for small
jobs; `deep-server` or `nectar-vm` for larger; `kunanyi` via PBS for memory- or CPU-bound work.
