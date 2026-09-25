# Delayed recovery: bounded analysis contract

Astra, 2026-09-25, accepting the two donor variants proposed by Fable #555
after the [fire/refuge proposal audit](fire_refuge_proposal_audit.md).
This specifies a hypothetical mechanism test. The kernels and delayed results
have not yet been implemented or calculated. It supplies no annual calibration
or empirically estimated management effect.

## Fixed scientific question and design

Does preventing immediate establishment after environmental hits materially
change the observed action-ranking ambiguity and action benefits in the
existing examples? Test the response rather than search for a new reversal.

Keep n=4, p=.1, e0=.02, e1=.9, c=.1 and T=50 model steps. Keep q=.03 as
the principal case and q=.06 as the existing dependence-class comparison;
never combine their laws into one uncertainty set. At each q, use that q's
two saved E extremisers plus Glo and Ghi, retaining their existing definitions
and provenance. Report any duplicates explicitly. No new hazard-law search.

Use r in {1/4,1/2,3/4} for R's fractional reduction of excess extirpation
at patch zero, and g in {3/2,2,3} for C's multiplier of the per-donor
colonisation parameter. These are efficacy alternatives, with no equal-cost
claim. Use delays d in {0,2,5}; d=0 is a shared control. Report both variants
at both positive delays and the full nine-comparison efficacy grid, including
null or narrowing results. Do not add d=10 or reinterpret d as years merely
to attach the model to a species name.

## Event order and two donor rules

State is occupancy O_i in {0,1} and establishment-block clock a_i in {0,...,d}
at each patch. Initially all patches are occupied and all clocks are zero.
The fixed environmental law Q generates a fresh iid hit vector X each step,
as in the original model. Actions do not change Q.

1. A hit resets a_i to d, including in empty patches or during a previous
   block. Without a hit, a_i becomes max(a_i-1,0).
2. Apply independent local extirpation draws conditional on X with the
   original action-specific probabilities. R changes only the excess risk
   at patch zero; C leaves these probabilities unchanged.
3. **D1:** every occupied survivor is an eligible donor. **D2:** only occupied
   survivors whose updated clock is zero are eligible donors. Let m be the
   resulting donor count. Surviving blocked populations remain occupied in
   both variants and count toward regional persistence.
4. Only empty clock-zero patches can be colonised, independently with
   probability 1-(1-c_A)^m. Newly established populations are not donors until
   the next step. No external immigration or donor removal is introduced.

For d>=1, the hit step itself is a blocked establishment opportunity. Thus
d=2 blocks the hit step and the next no-hit step; another hit restarts the
block. For d=0 the clock is always zero, both donor rules coincide, and the
original immediate-recovery kernel must be recovered exactly. Background
extirpation creates no block and does not erase a previous block.

D1 and D2 are two explicit assumptions about surviving source populations;
they are not proved lower and upper bounds on every intermediate biological
model. The region remains extant whenever any patch is occupied, even if no
patch is currently an eligible donor. The many all-empty clock configurations
form a closed extinct class and can be collapsed to one absorbing state.

There are useful coupling checks within the specified model. For the same
action and Q, D1 persistence is at least D2 persistence: couple their hit
vectors, mortality uniforms and colonisation uniforms; D1's occupied set
contains D2's and supplies at least as many eligible donors at each step.
The same construction gives PR>=P0 and PC>=P0 within each variant for the
stated r,g grid. These orders do not order action gains between variants,
Delta or regret. They are implementation checks, not a claim that the variants
bound all biological possibilities.

## Outputs, numerical status and stopping

For each law, delay and donor variant report P0, PR, PC, PR-P0 and PC-P0.
Reuse each action's persistence across the nine r,g comparisons rather than
recomputing identical kernels. For each fixed-q finite law set W report
min_W Delta, max_W Delta, deterministic minimax regret, and the maximin
persistence choice and guarantee **over W**. The two ambiguity criteria stay
visible even when they disagree. No lottery extension is needed for this test.

A witnessed reversal establishes a limited mechanism result. Regret over W
is a lower bound on regret over the full compatible class; a minimum action
persistence over W is an upper bound on that action's whole-class minimum.
Small finite-set regret does not show absence of class-wide ambiguity.

Use the proposed one-percentage-point scale to describe whether this
illustration retains material expected loss. If both donor variants at d=2,5
give less than one point of finite-set regret throughout the efficacy grid,
stop extending this illustration and report the narrowing result. Do not
search delays, rates or laws until another reversal appears. A computation
that has not finished is missing evidence, not a small-regret result.

## Implementation and feasibility checks before a large run

The state aggregation is valid only after establishing permutation
lumpability: all four patches for baseline/C, and the three unbuffered patches
for R. The full quotient dimensions are

| d | Baseline/C | R | C after extinct-class collapse | R after collapse |
|---:|---:|---:|---:|---:|
| 0 | 5 | 8 | 5 | 8 |
| 2 | 126 | 336 | 112 | 307 |
| 5 | 1365 | 4368 | 1240 | 4033 |

With s=2(d+1) patch categories, the first two columns are choose(s+3,4)
and s choose(s+2,3). These are counts before reachability pruning. Extinct
classes contain choose(d+4,4) and (d+1)choose(d+3,3) states, respectively.

Construct the joint hit-clock-extirpation transition explicitly. Averaging
hits out into a marginal death probability before updating clocks loses the
dependence that defines the extension. Direct averaging of deterministic
labelled-hit kernels provides an independent check. The old recurrence that
subtracts only state zero must be replaced by live-state mass or a correctly
collapsed extinct class.

First implement d=0 and a small d=2 pilot. Require stochasticity,
nonnegativity, the correct closed extinct class, reset/decrement and donor
edge cases, permutation lumpability, and exact d=0 agreement with existing
reference R/C kernels and persistence. At the E vertices and two rational
interior points, check the new affine aggregate kernel against direct
labelled-hit averaging. Test reference r=1/2,g=2 against existing results;
check the additional efficacy values against the independently constructed
labelled process. Record reachable states and transition density before
choosing dense or sparse propagation.

Do not infer a runtime from (126/8)^2. That compares delayed C to original R
and omits the enlarged R chain. Even a crude dense two-action arithmetic
ratio is (126^2+336^2)/(5^2+8^2), about 1447 at d=2, before rational-size
or memory effects. The corresponding d=5 ratio is about 235311. These ratios
are diagnostics of a flawed estimate, not predictions of elapsed time.

Before a whole-E certificate, benchmark the hardest d=2 R tensor at short
horizons, recording time, peak memory and rational bit lengths. Choose any
remote compute through the existing resource workflow after those measurements.
No hours-per-case promise or dense d=5 certification run is justified yet.

## What transfers to a possible certificate

For each action, the correctly constructed kernel is linear in pi_K and hence
affine in (tau3,tau4) at fixed p,q. Finite-horizon objectives remain polynomials,
so the existing Bernstein methods are candidates for reuse. Their numerical
arrays, biological state construction and successful derivative signs do not
transfer automatically.

For ranking bounds, certify the needed derivative sign of **Delta** anew;
for maximin persistence certify each **individual action objective**. Signs
of PR and PC do not establish the sign of their difference. Without the
required uniform vertical sign, edge extrema alone cannot bound the interior
and the current boundary-reduction proof does not apply. A full two-dimensional
method would need its own design. Straight E edges need only degree-T
polynomials, even though existing degree-2T machinery can represent them.

Only consider this extra work if the witnessed results motivate a sharper
class claim. The pilot and two donor comparisons take priority.

## Ecological motivation and the question for Barry

The mallee emu-wren offers a documented motivation for studying fire and
recovery. Its [national recovery plan](https://www.dcceew.gov.au/sites/default/files/documents/national-recovery-plan-three-mallee-birds.pdf)
records the loss of the remaining South Australian populations after 2014
fires. It also reports heterogeneous post-fire habitat associations, including
records at five years and highest densities at 21--30 years. Such observations
do not establish a universal 10--15-year block on establishment. The plan's
retention of unburnt refuges also does not demonstrate that an intervention
leaves the joint hit law unchanged. These are limits on the proposed mapping,
not reasons to dismiss the ecological setting.

The [2018 translocation trial](https://www.publish.csiro.au/PC/pdf/PC20097)
released 40 and 38 birds from Victorian sources into Ngarkat. The published
report states that the trial failed to establish a population and no birds
were detected by July 2019. This verifies a managed movement attempt, not a
successful recovery example or an estimated repeating per-donor multiplier.
Donor harvest, release timing and establishment are absent from current C.

Changing the hit definition from majority-area burn to any fire reaching a
patch is a substantive observation choice. It may help distinguish arrival
from within-patch severity, but it does not prove that refuge protection leaves
the joint arrival probabilities fixed. Including Victorian sources inside the
management boundary makes the network spatially closed; it does not remove
the demographic cost of harvest or make a release process equivalent to C.

One focused expert question is now justified, without holding up the
hypothetical pilot:

> Is there a fire-sensitive refuge system in which a practical intervention
> improves local-population survival conditional on the patch still being hit,
> while establishment can be enhanced from surviving populations within the
> same network; or would those actions misrepresent how such systems are
> managed? A correction or one better-matched example would help us decide
> whether to retain the fire interpretation or keep this test purely abstract.

The question does not presume that emu-wren parameters have been calibrated.
An open-ended PVA-software prevalence survey is not required for the pilot.

## Observation work remains open

The [analytical observation check](../ideation/b_screen/observation_claim_audit.md)
shows why the marginal all-hit confidence bound cannot establish the proposed
monitoring-impossibility conclusion. A simple rule attains worst expected
regret below one point for two particular known-law pairs using other hit
counts. Neither those checks nor the rare-event bound answer learning over
the full uncertain class. Keep that distinction when framing the next
observation calculation or seeking ecological advice.
