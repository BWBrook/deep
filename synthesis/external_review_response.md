# Response to the independent scientific review

Astra, 2026-09-25. The reviewer is GPT-6 Astra Pro; the reviewed revision is
`6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3`. This records review intake,
subsequent checks, and a proposed next stage for discussion with Fable. It is
not a claim that the ecological extensions below have been agreed or run.

**Subsequent discussion:** Fable's message 551 proposes fire-exposed refuge
networks, a fixed efficacy grid and delayed recovery. The
[proposal audit](fire_refuge_proposal_audit.md) retains fire as a candidate,
checks the observation and software claims, and specifies an alternative
establishment-delay rule for challenge. The subsequent
[contract](delayed_recovery_contract.md) fixes two hypothetical donor variants
and a bounded pilot; no annual calibration or delayed-kernel results are
available. The [observation check](../ideation/b_screen/observation_claim_audit.md)
separates an explicit known-law decision calculation from the still-open
problem of learning with full-class uncertainty.

## Assessment

The [review](../ideation/DEEP_independent_review_6ec79d8.md) supports continuing
with a focused theoretical ecology or quantitative conservation paper. Its
central criticism is persuasive: an exact information counterexample does
not establish that a realistic observation programme or conservation decision
would change. Clearer exposition helps readers understand the result; costs,
observations and recovery assumptions determine its ecological usefulness.

The existing conditional theorem remains valuable if a more realistic
comparison removes material decision loss. We should allow that narrowing
outcome. We should not respond by searching additional patch numbers or
parameter combinations until another large reversal appears.

## Evidence received and checked

Barry supplied the report and its accompanying ZIP. The
[import record](../ideation/reviews/6ec79d8/IMPORT.md) preserves the unchanged
archive, extracted source, environment, logs and receipts. The two supplied
copies of the report are identical, all 16 manifest hashes match, and a
replay of the package-validation script reproduces its saved result exactly:
32 joins between certificate and attained ranges, and 64 rounded endpoint
comparisons. This replay checks the supplied evidence joins, not all polynomial
calculations afresh.

Pro reports a separate implementation and full-domain checks of all 16
four-patch E cases and both S boundaries in all 16 cases. Its source imports
its own factorised transition implementation, not our analysis modules.
This is useful external computational evidence about the formal model. It
is neither empirical validation nor a guarantee that every historical result
or literature relationship has been independently checked. The review states
its exclusions, including larger-n results, the observation benchmark and
reliable inspection of rendered figures.

A targeted [local audit](../ideation/b_screen/review_additions_audit.md)
reconstructs the two new findings using our established exact kernels and
Bernstein machinery. It certifies separate whole-E action minima at q=.03
and .06, with zero attained/outer gaps, and checks the two full-support
mixtures by exact labelled propagation. Four action minima and all six
R/C/Delta values at the two mixtures equal Pro's supplied exact fractions.
This internal cross-check shares repository methods; it is distinct from
Pro's independently written implementation. We have not rerun all of Pro's
polynomial calculations.

| q | Whole-E minimum P_R | Whole-E minimum P_C | Maximin persistence | Minimax regret |
|---:|---:|---:|:---:|:---:|
| .03 | .3857762490 | .3214370903 | R | C |
| .06 | .2273638798 | .1149017905 | R | R |

At q=.03 the 90:10 extremiser/Glo mixtures give Delta approximately -.172970
and +.050918 while retaining identical pairwise information and positive
probability for every labelled hit pattern. These are feasible examples,
not additional optima or evidence of ecological plausibility. The certified
criterion comparison is now included in [Result 2](results_decision_loss.md).

A further algebraic implication is recorded there: both action minima share
the same law within each q, so every once-only lottery attains its worst
persistence at that law too. Its guarantee is the weighted average of the
two action minima. Pure R therefore also maximises worst-case persistence
over these lotteries, whereas a mixture improves minimax regret. This is a
deduction from the checked minima, not an additional numerical search.

## What we accept now

- Lead the eventual manuscript with the two q=.03 hit-count laws and their
  action outcomes. They expose identical complete pairwise tables and opposite
  rankings directly. They are the calculated E extremisers, not an independently
  chosen ecological case study.
- Explain ranking identification, the size of an unresolved loss, and a choice
  under a declared decision criterion separately. The S-to-E comparison leaves
  the deterministic minimax-regret choice unchanged in every reference case.
- Treat the objective-criterion comparison as essential context for regret.
  A manager can protect absolute worst-case persistence or the worst avoidable
  loss relative to the better package under the true law. Both use persistence
  as the outcome; the ambiguity criterion differs.
- Define the model event correctly: an occupied patch contains a local
  population, so the code's mortality is local extirpation. Connectivity denotes
  the prescribed change in per-donor recolonisation. The model specification
  now states these meanings explicitly without changing any kernel.
- Keep latent-shape assumptions separate from measured information. A
  restriction on theta is not established by a unimodal histogram of a driver
  or of hit counts. The observable moment implications can be tested even
  where arbitrary latent shape is not identifiable.
- Keep once-only lotteries secondary. They are not fractional investment or
  adaptive management. Keep the larger-n and named-family calculations as
  supporting results, with their existing scope limits.

## Where we should resist overinterpreting the review

The different maximin and minimax choices are not a defect in either rule,
nor another example of changing the dependence class. They use the same
uncertainty class and outcome. We should show both, rather than choose a
manager's tolerance for absolute loss or missed opportunity ourselves.
A rule cannot be justified merely by its yielding the preferred action.

The full-support mixtures answer a sparse-law objection. Positivity of every
hit pattern is not an ecological disturbance model. Similarly, a literal cost
ratio added to two unchanged kernels does not establish an equal-effort
comparison. Either costs constrain attainable response or their conversion
into the decision objective must be explicit.

The next observational calculation must go beyond reweighting the existing
two-known-law classification errors. That would change the loss function but
would still assume the alternatives and their p,q values are already known.
Likewise, adding a temporally correlated class that contains the original iid
laws automatically retains the old ambiguity. A useful stress test compares
specified nonzero clustering or recovery delay at controlled information.

## Proposed next stage: one decision problem, three linked steps

This is a proposal for Fable's challenge, not an adopted ecological model.
Use the existing four-patch example, deterministic R/C packages and absorbing
regional persistence as the starting point. Select one ecological interpretation
before choosing new numerical settings. An explicitly hypothetical worked
example remains acceptable for the narrower theoretical paper; a species name
and uncalibrated rates would not turn it into empirical evidence.

### 1. Specify the decision and action response

State the management unit, operational hit definition, time step, horizon,
feasible actions and ambiguity criterion. Retain both maximin persistence and
minimax regret as a sensitivity comparison unless a particular decision-maker's
objective has external justification.

For a compact response comparison, write R(r) as reducing the excess local
extirpation probability at the named patch by a fraction r, and C(g) as changing
c to gc, with 0<=r<=1 and 0<=gc<=1. The present packages have r=1/2 and g=2.
Specify a short set or interval of r,g values before computing their results,
using the chosen interpretation. Where an effort relation is available, compare
responses attainable at the same effort. Without it, describe an efficacy
sensitivity, not an equal-budget experiment. Do not choose the response range
because it preserves a reversal.

An additive cost sensitivity is possible only after defining costs in the
objective's units. If k_R and k_C are persistence-equivalent penalties, then
Delta shifts to Delta+d, where d=k_C-k_R. For fixed packages and sharp Delta
bounds [a,b], the deterministic regret rule switches at d=-(a+b)/2; maximin
switches at d=inf P_C-inf P_R. These are conditional thresholds, not estimated
prices or a substitute for an intervention response curve.

**Completion condition:** a reader can tell which outcome is protected, why
both packages are feasible competitors, and what magnitude of unresolved loss
would matter. If that contract cannot be defended, retain the conditional
information result and narrow the ecological claim.

### 2. Define an observable decision problem

A natural starting record is the number K of hit patches in each regional
observation interval, with all four patches observed. Under exchangeability,
its probability vector pi determines the hit law. It also determines
p=E[K]/4 and q=E[K(K-1)]/12, so lower moments must vary with pi in any honest
finite-record analysis. Pairs within one event are not independent replicates.

Use a specified observation model and uncertainty set or prior on pi. Explain
what coverage or repeated-sample performance is guaranteed, what remains
numerical, and how parameter uncertainty enters the action comparison. Report
persistence loss and residual decision guarantees against record length,
including the no-new-record baseline. A decision may improve without the
entire disturbance law becoming precisely known.

Do not claim a feasible sample size or a field monitoring recommendation until
the observation interval, dependence between records and measurement error are
justified. A full four-dimensional confidence-region optimisation also needs
an explicit computational plan; passing a two-law benchmark does not supply it.
This is the previously deferred D6 extension, now proposed for discussion.

**Completion condition:** identify either an attainable observation that
materially improves the chosen decision, or a defensible limit showing that
the proposed monitoring would not resolve useful uncertainty. Define materiality
before selecting favourable realised records or parameter values.

### 3. Challenge recovery timing once

Choose temporal clustering or delayed availability of disturbed patches from
the ecological interpretation in step 1. The choice precedes results. State
whether a delay follows the environmental hit or local extirpation, how existing
occupants and donors are treated, and when recovery occurs. These alternatives
have different kernels and should not be conflated.

Hold each compared law's within-step marginal and pairwise hit information
fixed, apply the same recovery rule and action response to both, and retain
the original iid/immediate-recovery calculation as a control. Report P_0,
P_R-P_0, P_C-P_0 and persistence regret. Feasible witnesses can answer a limited
mechanism question; they cannot be relabelled whole-class bounds.

**Completion condition:** either the proposed mechanism survives a meaningful
challenge at a material scale, or its failure supplies a clear boundary on the
model's use. Do not search across several recovery formulations simply to
retain the original ranking reversal.

## Coordination and stopping point

Astra has preserved the external review and prepared the targeted evidence
check and this response. The next joint scientific task is to settle the
compact decision contract and ecological interpretation with Fable. Fable's
most useful challenge is whether the proposed observation/recovery comparison
answers a recognisable ecological question at all, and what externally justified
assumption would make it do so.

No new field collection, broad parameter atlas, larger-n search or full labelled
class optimisation is proposed. No empirical calibration, general management
recommendation or publication decision has been made. If a substantive gap
requires Barry's ecological judgement, bring that specific gap with a concrete
proposal; do not ask him to approve routine implementation details.
