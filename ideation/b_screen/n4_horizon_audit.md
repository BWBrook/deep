# Horizon changes the size, and sometimes the existence, of ranking ambiguity

Fable's p/T scan records similar fractions of grid cells with opposite action
rankings at T=20,...,200. That observation does not imply that persistence
differences have converged by T=20 or that the decision is insensitive to T.
This audit follows four fixed parameter settings and computes exact rational
anchors at T=20,50,200. The same annual iid model and two endpoint laws are
used throughout, with n=4 and q=p^2+(2/9)p(1-p).

## Exact finite-horizon results

| Setting | p | e0 | e1 | c | Pair minimax regret at T=20 | T=50 | T=200 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Reference | 0.1 | 0.02 | 0.8 | 0.15 | 0.03623049 | 0.05292696 | 0.03556145 |
| Higher p, same actions | 0.3 | 0.02 | 0.8 | 0.15 | 0.00144696 | 0 | 0 |
| Higher p, reversal | 0.3 | 0.02 | 0.9 | 0.30 | 0.07455511 | 0.03807481 | 0.00002570 |
| Higher background | 0.1 | 0.10 | 0.9 | 0.45 | 0.03291202 | 0.07583193 | 0.05113253 |

The table gives rounded exact values for a **deterministic, once-only choice**
between R and C and uncertainty restricted to Glo/Ghi. For differences
`dl=PR(Glo)-PC(Glo), dh=PR(Ghi)-PC(Ghi)`, the worst regrets are
`rR=max(0,-dl,-dh)` and `rC=max(0,dl,dh)`, and the minimax regret is
`min(rR,rC)`. Positive values lower-bound the corresponding full-class
minimax regret. A zero for this pair is not a full-class zero. Randomised
policies, adaptive decisions, new observations and other actions are excluded.

The higher-p/same-actions case supplies an exact change of ranking with T:

- T=20: Delta(Glo)=-0.06772142, Delta(Ghi)=+0.00144696;
- T=50: Delta(Glo)=-0.03878432, Delta(Ghi)=-0.00540348.

The pair reverses at T=20, while both laws favour connectivity at T=50.
Similar aggregate grid counts can therefore conceal consequential changes
at individual settings.

There are also exact existence anchors for the broader scan: the higher-p
reversal case at p=0.3,T=50 has differences -0.07102939 and +0.03807481;
the reference case at p=0.1,T=200 has -0.18094228 and +0.03556145.
Strict signs persist in a neighbourhood by continuity, but this audit does
not measure those neighbourhoods or certify any grid fraction.

## The long-horizon limit is extinction, not a nonzero plateau

For e0>0 and e1>=e0, every occupied patch has mortality probability at least
e0 under either action, conditional on the annual hazard. From any nonempty
state with k<=4 occupied patches, the probability of losing every occupant
before recolonisation is at least e0^k>=e0^4. Without survivors there are no
donors, and extinction is absorbing. Thus, uniformly over admissible G,

`0 <= P_action(T) <= (1-e0^4)^T -> 0`.

Consequently `abs(Delta(T)) <= (1-e0^4)^T -> 0`, and both the pair and
full-class deterministic minimax regrets tend to zero for this two-action
menu. The bound is intentionally crude; it proves the limit rather than a
useful extinction time at small e0. It does not imply rapid convergence or
monotonic Delta, and the stated argument does not apply at e0=0.

The reference float curve's Ghi difference peaks at T=80 (~0.05687), while
the magnitude of its Glo difference peaks at T=180 (~0.18207). At T=20 the
reference Delta values are -0.03623,+0.03735; at T=200 they are
-0.18094,+0.03556. A claim that the differences had already settled by T=20
would be false here. The higher-p reversal case retains opposite signs at
T=200 but its pair minimax regret has fallen by roughly three orders of
magnitude from T=20. Persistence-ranking signs and decision-loss magnitude
are distinct outcomes.

The reported peak horizons are numerical maxima over integer T=0,...,2000,
not exact global maximisers over unbounded T. See the accompanying
`n4_horizon_audit.svg` for the first 1000 steps.

## Reproduction and checks

An independent numerical audit of all 28,000 p/T grid cells found one false
reversal caused by `b_screen_v2.persist` computing survival as `1-v[0]`.
At p=0.3,q=41/300,e0=0.02,e1=0.95,c=0.01,T=200, that expression produced
opposite differences around 1e-15. Exact rational propagation instead gives
positive differences under both laws: about 1.08118091e-15 and
1.08039520e-15. Both favour buffering, with negligible absolute differences.
The corrected routine sums the nonabsorbing states directly. The p=0.3,
T=200 panel has 560/1750 reversals (0.320), rather than 561/1750; the other
27,999 saved classifications agree with the stable calculation. Its smallest
marked c changes from 0.01 to 0.05, illustrating how one numerical artifact
can distort a reported range even when the count changes little. A small
regression check preserves this failure case. No arbitrary absolute sign
tolerance is needed to repair the cancellation.

This audit does not make the surviving grid signs exact certificates. In the
stable calculation the smallest retained reversal margin is about 6.14e-12,
and 77 cells have at least one action/law survival probability below 1e-12.
Magnitude should accompany sign when those cells are interpreted. Earlier
artifacts computed with the old routine remain historical and have not all
been regenerated; substantive claims from them need their own verification.

`n4_horizon_audit.py` reuses the exact N3 kernels, verifies exact stochasticity,
and propagates with gmpy2 rational arithmetic for all finite-T anchors. The
floating curves propagate only transient-state mass and sum it directly,
avoiding subtraction from an extinction probability close to one. Exact and
float metrics agree to absolute tolerance 2e-12 at all 12 anchors. This is
an arithmetic cross-check using shared kernels, not an independent ecological
model implementation. Source hashes and all anchor fractions are retained
in `n4_horizon_audit.json`.

```sh
uv run --python 3.12 --with numpy --with gmpy2 --with matplotlib python ideation/b_screen/n4_horizon_audit.py
```

The selected settings include earlier witnesses and exploratory high-p checks;
these are targeted mechanism tests, not a representative ecological sample.
Changing p while retaining hazard correlation also changes q and the two
endpoint laws. The higher-background case retains the fixed-e1 confounding
documented in `n4_background_mortality.md`. None of these time horizons or
interventions has been calibrated to a field system.
