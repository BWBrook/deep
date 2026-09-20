# Background mortality: a bounded four-patch test

At fixed shock moments `p=0.1, q=0.03` and horizon `T=50`, reversal between
the same two admissible mixing laws survives an entire background-mortality
interval from **0.04 to 0.14** at `e1=0.9, c=0.45`. It therefore does not require
the original background mortality of 0.02. This is a conditional model result,
not an empirically calibrated mortality range.

## Exact evidence

`n4_background_mortality.py` uses the existing mortality-then-recolonisation
model, starts all four patches occupied, and compares buffering the focal
patch's excess mortality by half against doubling connectivity. The two
fixed endpoint laws `Glo,Ghi` have equal first two hazard moments. Define
`Delta = P_R - P_C` for 50-step regional persistence.

| e1 | c | e0 | Delta(Glo) | Delta(Ghi) | Result for these two laws |
| --- | --- | --- | --- | --- | --- |
| 0.9 | 0.45 | 0.02 | +0.00404733 | +0.19754190 | Both favour buffering |
| 0.9 | 0.45 | 0.10 | -0.07583193 | +0.11023902 | Opposite rankings |
| 0.9 | 0.45 | 0.16 | -0.17649726 | -0.00059318 | Both favour connectivity |

These displayed values are rounded from exact rational calculations. The
interval `e0 in [1/25,7/50]` is covered by 24 adjacent, certified leaves, with
`Delta(Glo)<0<Delta(Ghi)` throughout every leaf. The smallest certified sign
margin across the retained enclosures is about 0.0001663; this conservative
bound is not the minimum true effect size.

At the earlier slice `e1=0.8,c=0.15`, another 24 leaves certify reversal for
all `e0 in [0,0.03]`. At `e0=0.04`, exact differences are -0.15718508 and
-0.01348489, so these two laws both favour connectivity there.

An additional exact point at `e0=0.20,e1=1,c=0.5` also reverses ranking
(-0.16968590 versus +0.01269195). It lies on the maximal-shock and
maximal-doubled-connectivity boundaries and does not establish an interior
region at that mortality.

The second slice was selected after an exploratory screen. Its exact
certificates verify the selected result; they are not a prospectively chosen
independent ecological validation.

## Interpretation and boundaries

For this fixed pair of laws, increasing background mortality can introduce
and subsequently remove a ranking reversal. Thus this pair's reversal sets
in `(e1,c)` are not nested as e0 increases. Agreement of these two laws does
**not** prove agreement throughout the admissible mixing-law class, so this
argument does not establish non-nesting of the full class's unidentified
region. It also does not locate a unique boundary or establish monotonicity
between checked intervals and points.

Persistence is a polynomial in the ecological parameters at finite T.
Strict opposite signs therefore persist locally around an interior witness;
there is no discontinuous jump in either action's persistence along these
slices. That does not settle the global geometry or area of the reversal set.

Changing e0 at fixed e1 also decreases the excess `e1-e0` that buffering
acts upon. Mean ordinary mortality increases, while the environmental
amplitude and induced mortality covariance `(e1-e0)^2*(q-p^2)` decrease.
Consequently this is sensitivity to the existing model's e0
parameter, not an isolated experiment in extra mortality with constant
buffering benefit. Connectivity here is substantial (`c=0.45`, doubled to
0.9). Neither intervention has a financial cost model, so rankings are of
persistence under the specified actions, not cost effectiveness.

At the interior point e0=0.10, choosing buffering loses 0.07583 persistence
probability under Glo, while choosing connectivity loses 0.11024 under Ghi.
For a deterministic choice between these two actions, the minimum worst-case
regret over these two laws is therefore 0.07583. This is also a lower bound
on deterministic minimax regret over the full admissible class, with these
parameters and this two-action menu. It is not the regret bound for a
randomised policy or an expanded action menu.

## Screen and certificate construction

The saved float screen covers e0=0,0.02,...,0.30; e1=0.4,0.5,...,1;
and c=0.01,0.02,...,0.50: 350 cells per mortality value. A cell is counted
only if both required signs exceed a margin of `1e-10`. The counts for
successive e0 values are 229,229,197,152,110,78,51,32,17,7,2,0,0,0,0,0.
Those counts are not area estimates. Zero counts at e0>=0.22 are not an
absence certificate, either for this pair or for the whole law class.
The exploratory narrower c<=0.30 screen stopped finding reversals at e0=0.14;
extending the admissible connectivity range found more. This directly
illustrates why a screen boundary must not be called a biological limit.

At fixed e1,c each transition-kernel entry is a polynomial of degree at most
four in e0: each occupied patch's mortality is affine in e0 and at most four
patches enter a survival event. Recolonisation conditional on survivors has
no e0 dependence. The code interpolates degree-four Bernstein coefficients
using exact rationals, bounds entries by their coefficient extrema, and
uses the previously reviewed nonnegative interval propagation and simplex
bounds. Adaptive bisection stops only with certified signs or a depth cap;
all leaves in both reported intervals certified reversal before that cap.

The JSON saves all leaf endpoints and exact sign enclosures, plus source
hashes, screen cells and exact point results. Checks cover off-node exact
polynomial reconstruction, stochastic kernels, point-interval collapse and
agreement with floating calculations. Existing kernel and interval routines
are reused: this is a new parameter study, not an independent implementation
of those routines.

Reproduce from the repository root:

```sh
uv run --python 3.12 --with gmpy2 python ideation/b_screen/n4_background_mortality.py
```

The bounded run takes about 13 seconds on this Mac mini. The strongest next
comparison is Fable's p/T sensitivity, followed by a clearly defined
persistence-regret contrast, rather than treating grid counts as prevalence
or highest-moment widths as a proxy for decision importance. A companion
path holding `e1-e0` constant would preserve buffering's absolute mortality
benefit and help separate the mechanisms in the present fixed-e1 path.

A separate method review found no certificate defect. An independently
implemented 16-state labelled float kernel reproduced the point comparisons
at `(e0,e1,c)=(.02,.9,.45),(.16,.9,.45),(.04,.8,.15),(.2,1,.5)` within
`6e-15`. That is independent numerical kernel evidence; exact interval
certification still uses the shared reviewed rational implementation.
