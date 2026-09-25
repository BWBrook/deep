# Audit of the four finite-exchangeable mechanism contrasts

Astra, 2026-09-25, responding to Fable #526 and commit `917ab75`.
The four-law arithmetic passes. Several interpretations in the original
`mechanism_decomposition.md` require correction, as detailed below. The
original files are preserved as the audited record; the revised scientific
account is in [Result 1](../../synthesis/results_identification.md).

## Scope and independent reproduction

Only the four already certified E extremisers were evaluated: n=4, p=.1,
e0=.02, e1=.9, c=.1, T=50 model steps, all patches initially occupied.
No law search, new parameter grid, policy optimisation or class-certificate
generation was undertaken. The unmodified baseline is a diagnostic comparison,
not an addition to the optimised action menu.

The audit rebuilt the 48 deterministic labelled kernels (16 hit patterns
times baseline, R and C) using the established independent survivor-subset
formula. All **12,288 entries** agree exactly with the original nested
survivor/colonisation construction. The exchangeable pattern probabilities
Q(x)=pi_|x|/choose(4,|x|) normalise exactly and reproduce all marginal, pair,
triple and quadruple hit moments. Mortality precedes same-step recolonisation
from survivors; regional extinction is absorbing.

Both T=50 telescoping allocations, their time and hit-count components, and
actual horizon-specific trajectories were recomputed with exact rationals.
The one-step identity holds at all four laws. Every B_t and Gamma_t is
nonnegative, and the final-step Gamma is exactly zero. Each total Delta_50
matches three saved values exactly: the mechanism output, the earlier
independent labelled-law witness, and the attained E-class certificate.
The input receipts' ten declared source hashes were checked. The new JSON
also binds this audit's source and inputs by SHA-256.

This external binding closes a small provenance gap: the original source's
`CERT` dictionary is unused, and its telescoping assertion compares with
Delta propagated from its own kernels. That assertion is a valid internal
consistency check, but is not a comparison with an earlier certificate.
No numerical discrepancy was found on making that comparison independently.
The earlier independent kernel implementation is reused; this audit is not
a third new implementation of the ecological model.

## Actual dynamics and complete-action gains

Values below are rounded probabilities, not directed numerical enclosures.
P0 uses ordinary mortality and the original c; R buffers patch 0 and C uses
2c throughout the trajectory. The receipt retains exact values.

| q | E extremiser of Delta_50 | P0 | PR | PC | PR-P0 | PC-P0 | Delta_50 |
|---:|:---|---:|---:|---:|---:|---:|---:|
| .03 | Minimum | .358170 | .500589 | .679841 | .142420 | .321672 | -.179252 |
| .03 | Maximum | .229714 | .385776 | .321437 | .156062 | .091723 | +.064339 |
| .06 | Minimum | .194005 | .337258 | .390681 | .143252 | .196676 | -.053423 |
| .06 | Maximum | .090743 | .227364 | .114902 | .136621 | .024159 | +.112462 |

Within each fixed-q pair, moving from the Delta-minimising to maximising law
reduces C's gain above baseline much more than it changes R's. At q=.06,
even R's gain decreases. The switch to an R ranking therefore need not mean
that buffering becomes more beneficial: connectivity becomes much less
beneficial relative to the same law's baseline. These are unambiguous
complete-action contrasts, with each action's induced occupancy history.
They complement the order-dependent one-step allocations below.

All four laws favour R at horizon 1. The q=.03 minimiser changes from
Delta_3=+.000329145 to Delta_4=-.000489807; the q=.06 minimiser changes from
Delta_10=+.001392220 to Delta_11=-.001057959. Both then favour C through
horizon 50. The maximisers favour R at every horizon 1 through 50. Equality
holds at horizon zero. These exact sign checks concern four fixed laws,
not class extrema optimised again at each horizon.

Full occupancy and persistence must also be distinguished. Under the
q=.03 maximiser at t=25, C has more probability on four occupied patches
(.3682 versus .2354), yet also more extinction (.4303 versus .3615).
Those occupancy distributions do not establish stochastic dominance.

## What the telescoping allocation does and does not show

Let f_t=mu K_R^t and v_t=K_C^(T-1-t)h for the forward ordering. Then

    B_t = f_t (K_R-K_0) v_t,
    Gamma_t = f_t (K_C-K_0) v_t,
    Delta_T = sum_t (B_t-Gamma_t).

The reverse ordering uses C in the past and R in the future. Both are
exact decompositions of the same total. They are hybrid trajectory contrasts,
not unique causal mediation estimates, and B and Gamma are not PR-P0 and
PC-P0. For example, the q=.03 minimiser has forward B=.163209 while its
complete-action PR-P0=.142420.

| q | Law | B forward / reverse | Gamma forward / reverse | B at largest admitted k, forward / reverse | Gamma at k=0, forward / reverse |
|---:|:---|---:|---:|---:|---:|
| .03 | Minimum | .163209 / .153111 | .342461 / .332363 | 93.95% / 94.08% | 61.76% / 51.54% |
| .03 | Maximum | .182144 / .184313 | .117805 / .119974 | 79.78% / 81.29% | 57.20% / 46.03% |
| .06 | Minimum | .184330 / .181813 | .237754 / .235236 | 87.32% / 87.63% | 72.32% / 61.35% |
| .06 | Maximum | .157009 / .160088 | .044547 / .047626 | 92.43% / 92.84% | 69.42% / 55.43% |

Buffering is concentrated at the largest hit count. Colonisation has a
substantial quiet-step contribution, but it does not accrue "almost entirely"
there. For the q=.03 maximiser in reverse order, singles supply 47.74% of
Gamma, exceeding the quiet-step share of 46.03%. Quiet means no environmental
hits; positive background mortality can still create vacancies in that step.

For each fixed-q pair, the move from minimum to maximum reduces Gamma under
both orderings. At q=.03 the decrease is .22466 forward and .21239 reverse,
against B increases of .01894 and .03120. At q=.06 Gamma falls by .19321 or
.18761, against B decreases of .02732 or .02173. This pattern agrees with
the direct baseline comparison. Detailed count attribution is less stable:
the q=.06 maximiser's k=1 net is +.001764 forward and -.003243 reverse.
For the q=.03 minimiser, the k=0 net moves by .04020 between orders even
though each total allocation moves by only .01010. Thus "order never changes
the qualitative picture" is too broad at the category level.

The saved `cum_Delta_by_t` in the original output is a prefix of a sum
defined for **T=50**, not Delta_t. The q=.03 minimiser's first forward term
is approximately -.00090, whereas its actual Delta_1 is +.0000142912.
The [new figure](../../synthesis/figures/mechanism_extremisers.pdf) uses direct
horizon calculations and complete-action baseline gains; neither is a
telescoping prefix.

## Hits, deaths and recolonisation

The laws share p=.1 throughout, but share q only within each of the two
pairs. Their expected co-hit pairs are .18 at q=.03 and .36 at q=.06,
not a common value across all four laws. Within each pair the minimiser
puts these co-hits into partial-hit events, while the maximiser puts them
into four-patch hits with probability q. Singles and no-hit probabilities
change too. No claim isolates a change in the all-hit category alone.

An all-four hit is an exposure pattern, not an all-four mortality event.
Conditional on such a hit from full occupancy, survival is
1-.9^4=.3439 under the baseline or C and 1-.9^3*.46=.66466 under R.
The conditional contrast is **.32076**, not the approximately .145
weighted k=4 buffering allocation. That allocation sums over hit probability,
occupancy at each step and future persistence values. Nor is it the entire
R advantage: for the q=.03 maximiser, the k=4 net is about .14069, offset
by other categories to give Delta_50=.06434.

If survivors remain, enhanced colonisation can refill vacancies created in
that same step, including after a four-patch hit. Colonisation cannot change
that step's regional nonextinction indicator: either mortality leaves no
survivor, in which case colonisation cannot start, or the region already
persists. This explains both Delta_1's lack of c and the zero terminal Gamma
without claiming that recolonisation waits for a later or quiet step.

The defensible mechanism is that retaining survivors and rebuilding
occupancy respond differently to the higher-order hit patterns left open
by pairwise tables. The four contrasts give model-specific evidence for
that trade-off; they do not establish universal separation of each action's
benefit into all-hit versus no-hit events.

## Reproduction and figure

```sh
uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/mechanism_audit.py
uv run --python 3.12 --with matplotlib python synthesis/plot_mechanism_extremisers.py
```

The exact audit takes about four seconds on the Mac mini. The
[JSON receipt](mechanism_audit.json) retains exact horizon contrasts,
per-time and per-count allocations, complete-action values, exact occupancy
anchors and source hashes. No class-wide optimisation is repeated.

Figure caption: **Dependence changes the relative value of rebuilding
occupancy.** Rows show q=.03 and q=.06 at the same ecological parameters.
Within each row, the two symmetric hit laws have identical complete pairwise
tables and attain the previously certified minimum and maximum Delta_50.
(A) The full hit-count distributions differ. (B) Exact action differences at
each integer horizon; lines join these discrete values. The extremising laws
are held fixed throughout. (C) Each complete action's 50-step improvement
above the unmodified baseline under that same law. These gains are distinct
from the hybrid telescoping allocations, and the baseline is not included
in the two-action optimisation. Values in the figure are rounded from the
exact receipt.
