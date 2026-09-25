# Result 2 draft: decision consequence and its limits

Fable, 2026-09-25; corrected in Astra's same-day source review. Companion to
[Result 1](results_identification.md). Sources: [D3 lottery audit](d3_lottery_audit.md),
[patch-number audit](../ideation/b_screen/n_scaling_audit.md),
[class table](class_identification_table.md), [E certificate audit](../ideation/b_screen/exchangeable_certificate_audit.md),
[failure-region audit](../ideation/b_screen/failure_region_audit.md),
[horizon audit](../ideation/b_screen/n4_horizon_audit.md),
[background-mortality note](../ideation/b_screen/n4_background_mortality.md),
[rectangle audit](../ideation/b_screen/n4_rectangle_audit.md). All values are exact
rationals in those records; decimals below are rounded displays.

## Loss criterion

A manager who must choose once between the prescribed packages R and C, and who
knows the pairwise hit tables and adopts a specified dependence class, faces
a set of compatible laws. Let a and b be the smallest and largest values of
Delta = P_R - P_C over that set.
Choosing R loses at most max(0, -a); choosing C loses at most max(0, b). The
deterministic minimax regret is V = min(max(0, -a), max(0, b)), the worst-case
loss of persistence probability at the horizon for the better of the two
deterministic choices. When a and b are certified class extrema, V is the exact
class value; when they come from a selected pair of feasible laws, V is a lower
bound for the class. Zero pair regret does not imply zero class regret.

## Exact class values at small n

At the reference setting (p = .1, q = .03, e0 = .02, e1 = .9, c = .1, T = 50):

| class | n | Delta range | V (pp) | once-only lottery regret (pp) | Pr(R) in the lottery |
|---|---:|---|---:|---:|---:|
| S (scalar hazard) | 3 | [-.02342, +.01720] | 1.720 | 0.992 | .4235 |
| S | 4 | [-.12052, +.03209] | 3.209 | 2.534 | .2103 |
| E (finite exchangeable) | 4 | [-.17925, +.06434] | 6.434 | 4.735 | .2641 |
| E, q = .06 | 4 | [-.05342, +.11246] | 5.342 | 3.622 | .6780 |

The last row is the design case that S certifies for R (V = 0 in S) and E does
not. Between S and E the best deterministic choice under minimax regret is
unchanged in every one of the 16 design cases; what changes is the guarantee,
the loss scale (a factor 2.005 at the reference case) and the optimal lottery.
This statement is specific to minimax regret; maximin persistence and
prior-weighted criteria were not evaluated.

## The limit of hedging by lottery

When a < 0 < b, drawing R with probability alpha once, before the shock
sequence, has worst regret max((1 - alpha) b, -alpha a). The unknown hit law
is fixed independently of the realised draw, and the criterion averages
over that draw. The optimum alpha* = b/(b - a) gives
V_mix = -ab/(b - a). With -a and b both positive, V_mix / V = max(-a, b)/(b - a),
which lies in [1/2, 1) for every such problem. A once-only lottery over two
complete packages can therefore halve the worst-case loss at best, and does so
only when the two extreme losses are equal. This is a limit of that policy
class, not of adaptive or within-horizon switching policies, which are outside
the analysis.

## Larger metapopulations: lower bounds

Certified extrema exist only for n = 3 and 4. For n = 5 to 8 at the same setting
in S, exact evaluation of saved feasible-law pairs gives regret lower bounds
(rounded down; the class may have larger regret):

| n | V lower bound (pp) | lottery lower bound (pp) |
|---:|---:|---:|
| 5 | 8.067 | 5.097 |
| 6 | 9.817 | 5.202 |
| 7 | 5.550 | 3.798 |
| 8 | 3.188 | 2.521 |

Within S, the n = 5 and n = 6 bounds each exceed the exact n = 3 and n = 4
class values, so the loss at five and six patches is larger than at three or
four. These witnesses also belong to the corresponding finite-exchangeable
class E_n, so the displayed regret lower bounds carry over to E_n for the
same fixed-patch R/C menu. No whole-class E_n optimum at n >= 5 is calculated.
The bounds do not order the unknown n = 5 and n = 6 optima, do not establish a peak, and
do not isolate the number of free moments as the cause, because n also changes
habitat amount, the fraction buffered and the number of colonisation sources.
As a companion scale, at the n = 6 endpoint law Glo the worse action's
extinction risk is 2.49 times the better action's (.164 versus .066); at the
n = 6 law Ghi the ratio is 1.28. These are law-specific ratios, not minimax
quantities.

## Where and how far the loss extends

Within S at n = 4, p = .1, q = .03, e0 = .02, T = 50, exact point certificates
locate unidentified rankings at (e1, c) = (.5, .22), (.7, .15), (.85, .02),
(.5, .40) and the reference (.9, .10), with class regrets .529, 1.794, .117,
.382 and 3.209 pp, while (.5, .18), (.7, .10) and (.9, .06) are identified for
C and (.95, .45) for R. These are point classifications; they force a
transition between neighbouring identified and unidentified points but do not
certify a boundary, an area or a fraction of the plane. For the two endpoint
laws, opposite rankings are certified on a positive-area rectangle of
(e1, c), throughout e0 in [0, .03] at (e1, c) = (.8, .15) and e0 in [.04, .14] at
(.9, .45), and along constant-amplitude background-mortality paths. Endpoint-pair
reversals lower-bound class ambiguity; their sets are not nested in e0.

The loss depends on the horizon. For e0 > 0 and e1 >= e0, extinction is
absorbing and every nonempty state has extinction probability at least e0^n
per step. Thus P_action(T) <= (1 - e0^n)^T, so persistence, Delta and the class regret all tend
to zero as T grows. This bound proves a limit, not a useful convergence rate
at small e0; it does not give extinction at e0 = 0. At the witness pair with
p=.3, q=41/300, e0=.02, e1=.9 and c=.3, regret falls from 7.456 pp at
T = 20 to .00257 pp at T = 200 with the signs unchanged. At p=.3, q=41/300,
e0=.02, e1=.8 and c=.15, Ghi changes from favouring R at T=20 to favouring C
at T=50, while Glo favours C at both. Stable sign counts across horizons
therefore do not imply stable benefits.

## What Result 2 does not claim

The packages are prescribed and uncosted; "doubling colonisation" and "halving
one patch's excess shock mortality" are not calibrated to equal effort, so
regret in persistence units is not a cost-effectiveness statement. Model steps
are uncalibrated. Grid counts and point classifications are not ecological
prevalence. Nothing here evaluates adaptive management, partial buffering of
several patches, or mixtures of the two packages within a horizon.
