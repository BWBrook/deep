# Independent review: proposed timing and priming prompt

Prepared by Astra, 2026-09-25, following Barry's request. **Draft for Fable's
critique; no external review has been commissioned.**

## Proposed review point

Astra proposes commissioning the first independent review after the current
Results 1–3 source-checking pass and Fable's response to the corrected drafts,
before further exploratory calculations or a journal-specific rewrite.
The model specification below, draft results and source-linked evidence are
enough to expose the scientific proposition to criticism. Review readiness
does not require the authors to have resolved novelty or ecological value;
those are central matters for the reviewer.

Before asking Barry to commission it, agree the timing with Fable, confirm
the repository contains the complete reading set, and provide the full
immutable Git commit in place of the revision placeholder below. The review
should cover that commit even if main subsequently advances. This file is a
draft brief, not evidence that those steps or a joint decision have occurred.

The reviewer has GitHub access to the private repository, but no access to
the correspondence board or the authors' chat histories. The scientific
specification, claims, evidence types and material limitations must therefore
be understandable from repository files alone. Historical filenames such as
“audit” and “certificate” do not give their conclusions privileged status.

## Prompt to provide to GPT-6 Pro

Please conduct an independent scientific review of the private GitHub
repository **BWBrook/deep**, at **[FULL REVIEW COMMIT TO BE SUPPLIED]**.
Record the revision you can actually access. If your connector cannot
retrieve that revision, say which version or files it supplied and qualify
comparisons that could depend on later changes. You do not have access to
the private message board or the authors' conversation histories.

The scientific question is whether specified marginal and pairwise
environmental-hit information determines the finite-horizon persistence
ranking of two prescribed conservation interventions in a stochastic
metapopulation model. The interventions protect one patch from part of the
hit-related mortality or increase recolonisation. The work examines different
dependence assumptions, the loss associated with an unresolved ranking, and
additional restrictions that might determine it. The motivation is to
understand what information a conservation-action comparison needs, rather
than to calibrate a recommendation for a particular landscape at this stage.

The work was developed through collaboration between two models, Astra and
Fable, with limited input from Barry W. Brook, an ecologist at the University
of Tasmania. Barry is separately interested in how useful this mode of
research collaboration can be. Your primary task is to judge the scientific
work on its evidence. Neither model agreement nor the experiment's ambition
establishes correctness, novelty or significance. A positive, negative or
mixed assessment is equally useful when supported by reasons.

Begin with [synthesis/model_specification.md](model_specification.md) and the
three current draft results:

1. [synthesis/results_identification.md](results_identification.md)
2. [synthesis/results_decision_loss.md](results_decision_loss.md)
3. [synthesis/results_shape_information.md](results_shape_information.md)

Form your own account of the question, assumptions and strongest claims.
Then use [synthesis/claim_evidence_antecedents.md](claim_evidence_antecedents.md)
as an index to the underlying code, derivations, exact receipts and exploratory
outputs in `ideation/b_screen/`. The evidence index is the authors' account;
you should assess its labels and claimed scopes. Consult other repository
files and Git history as your connector permits. Historical probes sometimes have later
corrections, so resolve discrepancies against their calculations and stated
domains rather than choosing the more favourable prose.

[synthesis/d1_antecedent_search.md](d1_antecedent_search.md) records a bounded
literature search, and [synthesis/paper_argument.md](paper_argument.md) gives
the authors' proposed narrative. Read these after forming your initial
assessment if that helps preserve an independent perspective. Verify the
important literature relationships yourself where access permits, and
consider antecedents the authors may have missed. The existing search is
not a claim that no prior work overlaps.

Assess the following, adding other issues you consider consequential:

- Whether the model, information sets, dependence classes and decision
  criteria are defined consistently, and whether the equations, code and
  claimed certificate domains support the stated conclusions.
- Which claims are mathematical results, which are demonstrated examples or
  finite samples, and which are interpretation. Identify both sound central
  conclusions and any consequential gap or counterexample.
- The contribution relative to prior work: what is established methodology,
  what is an application or synthesis, and what, if anything, is a substantive
  new scientific insight.
- The ecological meaning and limits of the model and interventions, and
  whether the conditional findings matter for the motivating information
  problem. Judge the effect of simplifying assumptions on the actual claims.
- Whether further development is worthwhile, and which small number of
  additional analyses, comparisons or conceptual changes would most change
  that judgement. A case for narrowing or stopping is as useful as a case
  for extension.

Support findings with repository paths, equations, numerical values or
external sources as appropriate. Distinguish what you independently derived
or recomputed from what you inspected or accepted conditionally. If you lack
execution or literature access, identify the resulting limits rather than
implying those checks occurred. Do not assume that exact arithmetic by itself
validates a model or a claimed global optimisation.

Please return a review in the conversation, without editing the repository
or contacting the authors elsewhere. Include a short overall assessment,
major findings ordered by consequence, the strongest defensible contribution
as you see it, and a prioritised recommendation for what should happen next.
Use no quota of positive or negative findings. Preserve uncertainty where the
available evidence does not decide the issue.
