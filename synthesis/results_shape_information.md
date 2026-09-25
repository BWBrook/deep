# Result 3 draft: information that restores the ranking

Fable, 2026-09-25. Internal manuscript text for Astra's review. Companion to
[Result 1](results_identification.md) and [Result 2](results_decision_loss.md).
Sources: [unimodal audit](../ideation/b_screen/unimodal_audit.md),
[unimodal check](../ideation/b_screen/unimodal_check.md),
[hazard-family audit](../ideation/b_screen/hazard_family_audit.md),
[E certificate audit](../ideation/b_screen/exchangeable_certificate_audit.md).

## Two kinds of restoring information

Result 1 showed that, at matched pairwise information, the scalar-hazard
representation itself supplies identifying information: it excludes the
extreme concentrations of co-hits that E admits, and in one design case that
exclusion is what certifies the ranking. This section asks what further
information, within the scalar class S, identifies the ranking at the
reference case where S itself does not (n = 4, p = .1, q = .03, e0 = .02,
e1 = .9, c = .1, T = 50; Delta ranges over [-.12052, +.03209] in S).

## A shape restriction that identifies the ranking

Restrict the annual hazard law G to be unimodal in Khinchin's sense: theta =
m + U Y with U uniform on [0, 1] independent of Y on [-m, 1 - m]. Equivalently,
the distribution function is convex below the mode m and concave above it. The
class allows a mode at either boundary (all monotone densities), an atom at the
mode, and mixtures of uniforms sharing the mode; it excludes atoms away from
the mode and any law with two local maxima. The full-class extremisers Glo and
Ghi (two atoms each) lie outside it.

The moment identity E[theta^k] = sum_j C(k, j) m^(k - j) E[Y^j]/(j + 1) turns
the fixed p and q into fixed first and second moments of Y, and the feasible
modes are exactly 0 <= m <= 11/80. The third moment of Y on an interval is
bounded above by the two-point law at the interior point a = y1 - Var(Y)/(B - y1)
and the endpoint B (from E[(B - Y)(Y - a)^2] >= 0), giving the exact fixed-mode
ceiling

    m3 <= h(m) = (-7500 m^2 + 1780 m + 599) / [8000 (5m + 4)],

whose supremum over feasible modes is .0197404736 at m = sqrt(3)/2 - 4/5. The
certified upper-boundary root of Delta lies at m3 ~ .0202621; with r = 79/4000,
the identity r - h(m) = [7500 (m - 33/500)^2 + 33/100] / [8000 (5m + 4)] > 0
proves m3 < r for every compatible unimodal law. Since dDelta/dm4 > 0 on the
lens and Delta is increasing along the upper boundary, both certified,

    P_C - P_R > 427/100000

for every Khinchin-unimodal G compatible with (p, q). Connectivity is
identified throughout the unimodal class at this setting, by more than .427 pp.
The margin is sufficient, not sharp: an exploratory search over unimodal laws
found Delta as large as about -.0160, so the true unimodal supremum lies in
[-.0160, -.0043]. The m3-maximising unimodal law is an atom at the mode plus a
uniform tail to 1, (1 - w) delta_m* + w Uniform[m*, 1] with w = (58 - 30 sqrt(3))/83;
it is not shown to maximise Delta.

## Why the restriction works, and what it costs

The unimodal ceiling on m3 depends only on (p, q), not on the ecological
parameters. The restriction therefore identifies the ranking at any setting
whose certified ambiguity band in m3 lies above .01974 at these (p, q); it is
a general-purpose device evaluated here at one point. Its content, in the
mechanism of Result 1, is that a unimodal annual hazard cannot concentrate
co-hits into rare severe years to the extent that reverses the ranking: the
two-regime and footprint-type laws that favour buffering are exactly the laws
it excludes. The cost is that the exclusion is by assumption. Mixtures of two
Khinchin-unimodal components already contain Glo and Ghi (each atom is a
degenerate unimodal law) and so recover the full ambiguity; nonunimodality is
necessary for buffering to win here, but a second biological regime is not
thereby established as the ecological cause.

Three fitted parametric comparators (beta, probit-normal, logit-normal), each
matched to (p, q), all favour C at this setting. They are consistent with the
theorem but are not class results, and agreement of fitted families is not
evidence about the shape assumption's truth.

## Relation to the exchangeable class

Unimodality is a statement about a scalar hazard law and has no analogue in
E, where there is no latent theta. In E the ranking at this setting is
unidentified with regret 6.43 pp (Result 1); the unimodal theorem operates
inside S, whose own contribution to identification is the conditional-iid
representation. A reader who does not accept a regional scalar hazard cannot
use the shape restriction; a reader who accepts it gains identification from
a single qualitative assumption about the hazard's annual distribution.

## What Result 3 does not claim

Not a theorem at other (e1, c), n or T. Not a statement that unimodality of
theta is testable from hit records at n = 4, where the count law retains only
four moments; no estimability or observation model is analysed. Not an
empirical claim that regional hazards are unimodal. The identified action is
determined by the certified sign of Delta on the upper boundary below the
root; the analysis provides the sign and a sufficient margin, not the exact
unimodal-class optimum.
