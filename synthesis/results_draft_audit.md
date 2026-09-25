# Source review of Results 2–3 and the boundary sample

Astra, 2026-09-25, reviewing Fable #530 and commit `d894a24`. The main
small-n regret and reference unimodal-ranking results stand. The current
drafts correct the numerical and scope issues below. This note separates
new checks from previously established dependencies; it is not an external
scientific review.

## Decision-loss draft

The sharp n=3,4 values and lottery probabilities agree with the existing
D3 and class-comparison records. The following corrections matter:

- At the scalar pilot point P7, e1=.5 and c=.4, the exact attained endpoint
  differences are approximately -.00381685780474484 and +.02426974831980456.
  The deterministic minimax regret is their smaller opposing loss,
  **.00381685780474484**, or **.381685780474484 percentage points**. The original
  2.427 pp used C's worst loss rather than the minimax value; R is the
  minimax action at this point. The source is
  [failure_region_audit.json](../ideation/b_screen/failure_region_audit.json),
  `results[id=P7].direct_endpoint_delta`. All eight point regrets were
  recalculated from their retained exact endpoint fractions.
- Values labelled downward-rounded lower bounds must remain below the exact
  values after conversion to percentage points. The corrected n=7 row is
  **5.550 and 3.798 pp** (deterministic and lottery), and the n=8 row is
  **3.188 and 2.521 pp**. The n=5,6 displays already rounded down correctly.
  These checks use the exact pair-decision fields in
  [n_scaling_audit.json](../ideation/b_screen/n_scaling_audit.json).
- The interior lottery formula requires a<0<b. Its expected regret averages
  over a once-only draw, with the hit law fixed independently of the realised
  action. It is not a guarantee for each realised package.
- The uniform long-horizon extinction bound requires e0>0 and e1>=e0.
  The draft now states these assumptions and identifies the two horizon
  examples by parameters and laws. The bound gives no useful rate at small
  e0 and no extinction conclusion at e0=0.

Fable asked whether the n>=5 bounds also apply beyond S. For each fixed n,
the same feasible witnesses belong to E_n. With the same fixed-patch R/C
menu, enlarging the ambiguity class cannot reduce minimax regret, whether
deterministic or averaged over a once-only lottery. The lower bounds
therefore carry over by inclusion. The draft says so while keeping the
comparison with the known n=3,4 optima explicitly within S. No new E_n
optimisation has been performed.

## Shape-information draft

The feasible modes, fixed-mode third-moment ceiling, its exact maximiser and
the reference guarantee P_C-P_R>.00427 agree with
[unimodal_audit.md](../ideation/b_screen/unimodal_audit.md). The seven hashes
in its receipt match the current input files. The full vertical certificate
was not regenerated in this drafting pass.

The displayed supremum interval `[-.0160,-.0043]` was not justified:
rounding the negative upper bound down to -.0043 excludes values the
certificate still permits. An existing rational feasible law supplies a
rigorous lower comparison without another search. In the representation
theta=m+UY, with U uniform and independent, take

    m = 33/500,
    Y support = {11/216500, 467/500},
    Y weights = {187489/202200, 14711/202200}.

Its first four moments are

    m1=1/10, m2=3/100,
    m3=68381/3464000, m4=14621110637/937445000000.

Exact propagation at the reference setting gives
Delta=-.015965488433553693... and checks Delta>-.016 exactly. The existing
certified upper-envelope evaluation at r=79/4000 is
-.004272706373120518..., strictly below -.00427. Together they justify

    -.016 < sup_(compatible unimodal G) Delta < -.00427.

The witness was re-evaluated using the established exact transition
construction; this is not an independent rederivation of every transition
or a new optimisation of the unimodal class. Its reproducible parameters
are given above, and the earlier envelope proof remains a named dependency.

Reproduce that existing witness from the repository root:

```sh
uv run --python 3.12 --with gmpy2 --with numpy --with sympy python - <<'PY'
from fractions import Fraction as F
import sys
sys.path.insert(0, 'ideation/b_screen')
import unimodal_check as U
import unimodal_audit as A
import n4_vertical as V

m = F(33, 500)
_, _, _, law, feasible = U.m3_bounds(m, U.P, U.Qm)
mom = U.uniform_mixture_moments(m, law)
assert feasible and mom[:3] == [F(1), U.P, U.Qm]
values = {}
for action in ('R', 'C'):
    parts = V.affine_parts(action, 4, U.E0, U.E1, U.C, U.P, U.Qm)
    kernel = V.kernel_at(*parts, mom[3], mom[4])
    assert all(sum(row) == 1 and min(row) >= 0 for row in kernel)
    values[action] = A.propagate(kernel, V.start_state(action, 4))
delta = values['R'] - values['C']
assert delta > F(-16, 1000)
print('Y support and weights:', law)
print('moments:', mom)
print('Delta:', float(delta), 'exactly greater than -0.016')
PY
```

Further conceptual repairs are in the draft:

- The (p,q)-only ceiling can be reused at another ecological setting only
  with the needed sign and domain bounds for that setting's persistence
  envelope. A threshold above rounded .01974 is insufficient, since the
  true ceiling is .019740473580835506....
- Nonunimodality is necessary for R to win at this reference setting, but
  it is not sufficient. The excluded two-atom Glo favours C. The restriction
  does not identify two biological regimes as the cause of a reversal.
- Khinchin unimodality induces a subclass of S within E. E alone supplies
  no such latent-shape restriction; it is too strong to say that it has
  “no analogue” in E. Unimodality of the finite hit-count law is a different
  property and was not tested here.
- “Annual” and “severe years” have been replaced by model-step language.
  The precise mixture/CDF definition replaces an informal count of local
  maxima. The beta comparison is exact; probit/logit matches are numerical,
  and those integrations were not repeated in this pass.

## What the boundary computation establishes

The new [sample audit](../ideation/b_screen/baseline_gain_edges_audit.json)
and its [script](../ideation/b_screen/baseline_gain_edges_audit.py) reconstruct
every saved point using the independent labelled survivor-subset kernels
and compare them exactly with the source's reduced-kernel construction.
All **70 E records and 44 S records** match; the saved floats match their
exact evaluations. All laws are feasible, and the four previous extremisers
agree exactly with the mechanism receipt. The calculation takes about
12 seconds on the Mac mini and retains exact ranges, sample locations and
source hashes. It adds no sample points or parameter search.

There are 30 unique E-boundary laws at q=.03 and 40 at q=.06. Each edge has
its start vertex and nine strict interior points; its end vertex is the next
edge's start. The 22 S-boundary records at each q represent **20 unique laws**,
because the upper and lower curves coincide at both endpoints.

| q | E sampled buffering-gain span | E sampled colonisation-gain span | Best constant-minus-colonisation maximum residual, on these samples |
|---:|---:|---:|---:|
| .03 | .0136425200 | .2299485342 | .0068212600 |
| .06 | .0076454075 | .1725170582 | .0038227037 |

The last column is half the first: since Delta=(PR-P0)-(PC-P0), choosing
the midpoint of the observed buffering-gain range gives the best uniform
constant approximation on the finite sample. It is a descriptive identity,
not an additional causal explanation or a bound between sample points.
“Nearly invariant” also requires an explicit tolerance if intended as more
than a relative description of the displayed variation.

The source evaluates finitely many points **on** each boundary. It does not
evaluate or certify the whole continuous boundary. At q=.06 the largest
observed buffering gain occurs at an edge-interior sample,
(tau3,tau4)=(1/30,1/50); even vertex attainment for this gain is not established.
No continuous-boundary range, monotonicity, interior bound or class-wide
near-invariance follows from this sample. The source's “whole boundary”
language is superseded by this audit. Its original JSON retained only
floating summaries of gains and most coordinates; the new receipt retains
their reconstructed rational values.

The four-law figure already labels its trajectories as fixed-law evaluations.
It is retained as the main illustration. The additional sample is supporting
exploration; it is not promoted to a continuous curve certificate or a new
main result. No more boundary sampling is required before the proposed
independent review.

## Review readiness

The [model specification](model_specification.md) gives a repository-only
entry point to the process, information sets and decision criteria. Fable
accepted these corrections and agreed that this Results 1–3 package is ready
for the first independent GPT-6 Pro review, before further exploration or
journal drafting. The specification now includes the shared within-step
theta draw, the higher-moment feasibility conditions and a five-term evidence
glossary. The [external review prompt](external_review_prompt.md) incorporates
the agreed neutrality edits and lets the reviewer choose the subsequent
reading order after the README and model specification.

Fable's documentation-only commit `ee00417` added supersession notes to
the original mechanism, boundary-sample and unimodality notes. This changed
two inputs recorded in exact-audit receipts. To bind the receipts to the
corrected notes, `mechanism_audit.py` and then `baseline_gain_edges_audit.py`
were rerun. Both passed: 48 deterministic kernels and 12,288 entries for the
four mechanism laws, followed by all 70 E and 44 S boundary records. Comparing
the complete JSON objects before and after regeneration found changes only
in `source_hashes_sha256` and `seconds_wall`. All scientific payloads are
identical. The boundary receipt also records the regenerated mechanism
receipt's hash. The unchanged numerical inputs leave the existing four-law
figure unchanged; the unimodal receipt's seven source hashes still match.

The complete package is to be committed and pushed before Barry receives
the commissioning message. That message supplies the full immutable commit
for use with the brief, avoiding a self-referential revision inside this file.
No external model had assessed the package at that commissioning point.
The subsequent [GPT-6 Astra Pro review](../ideation/DEEP_independent_review_6ec79d8.md)
and [response](external_review_response.md) are recorded separately; they do
not change the historical scope of the source checks above.
