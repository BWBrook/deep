# How the (p,q)-ambiguity of the ranking scales with patch number: witness probe (Fable)

Files: `n_scaling_probe.py`, `n_scaling_probe.json`, `n_scaling_probe.log`. Run from the repo root:
`uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/n_scaling_probe.py` (about 10 s).

## Question

The certified whole-class results stop at n = 4. The matched-subset witnesses (n = 4, 5, 6) fix all but the top
moment and their regret shrinks roughly threefold per patch, which invites the reading that the phenomenon is
confined to tiny systems. But that is the wrong comparison for the paper's information set: with only (p, q) fixed,
the number of free moments grows with n (n - 2). Does the (p, q)-ambiguity of the R-vs-C ranking grow or shrink with n?

## Method

Setting p = 1/10, q = 3/100, e0 = 1/50, e1 = 9/10, c = 1/10, T = 50 for every n in 3..8. Float search (40000 random
3- and 4-atom laws with (p, q) matched exactly by solving the moment equations, nonnegative weights; Nelder-Mead
refinement of the 12 best in each direction), plus the two endpoint laws Glo = {0, 3/10} and Ghi = {7/90, 1}. The
best witnesses are then rationalised and re-evaluated exactly (direct averaging of the exact kernels,
`n4_vertical.avg_direct`, 50 exact steps). Every reported number is an attained Delta of a feasible law: a rigorous
lower bound on the class maximum or upper bound on the class minimum, and the implied regret
min(max(0, -L), max(0, U)) is a rigorous **lower bound** on whole-class deterministic minimax regret. Calibration:
at n = 3 and n = 4 the search recovers the certified extrema to all displayed digits.

## Result (exact re-evaluations; witnesses in the JSON)

| n | Delta min found | Delta max found | regret lower bound | extremes at Glo/Ghi? | P_R, P_C at Glo |
|---:|---:|---:|---:|:---:|---:|
| 3 | -0.023416 (Glo) | +0.017202 (Ghi) | 0.01720 (= certified class value) | yes | 0.257, 0.280 |
| 4 | -0.120520 (Glo) | +0.032093 (Ghi) | 0.03209 (= certified class value) | yes | 0.472, 0.593 |
| 5 | -0.138463 (Glo) | +0.080673 (Ghi) | >= 0.08067 | yes (among searched) | 0.683, 0.821 |
| 6 | -0.098171 (Glo) | +0.110662 (Ghi) | >= 0.09817 | yes (among searched) | 0.836, 0.934 |
| 7 | -0.055509 (interior two-point law {0.0288, 0.3810}) | +0.120356 (Ghi) | >= 0.05551 | no | 0.924, 0.978 |
| 8 | -0.031886 (interior two-point law {0.0484, 0.4878}) | +0.120628 (Ghi) | >= 0.03189 | no | 0.968, 0.993 |

- The (p, q)-ambiguity does not vanish with n. The regret lower bound rises from 0.017 (n = 3) to about 0.10 at
  n = 6, three times the n = 4 value, then declines at n = 7, 8 as persistence under both actions saturates near 1
  at this hazard level (the absolute scale of any Delta shrinks with 1 - P).
- The maximiser found is Ghi at every n; the minimiser is Glo up to n = 6 and an interior two-point law at n = 7, 8.
  So the endpoint-law pair is not the extremal pair in general, which also means the n = 5, 6 rows are lower bounds
  from a search, not attained class extrema.
- The decline at n >= 7 is a saturation effect at these parameters, not evidence that large metapopulations are
  immune: at higher hazard (larger p or e1) or shorter T the saturation point moves to larger n. Not tested here.

## What this does and does not show

Shows: at the reference setting, the ranking ambiguity under complete pairwise information is larger at n = 5 and
n = 6 than at n = 3 and n = 4, rigorously (attained witnesses, exact arithmetic). The paper's phenomenon is not a
small-n curiosity in this sense.

Does not show: whole-class extrema for n >= 5 (no certificate; the class may be wider than the witnesses), any
statement at other parameters, or anything about the matched-subset decay, which is a different question
(top-moment-only information) and stands as recorded.
