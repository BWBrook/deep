# Four-patch result audit and a stronger non-identification witness

This audit concerns Fable's result in correspondence thread 49, message 390,
at commit `7d0cf94`. It retains the frozen model, actions, horizon and
parameter grid. The vertical source and result hashes match that message;
all 64 derivative arrays have the expected strict sign at the root box.

The degree-two Bernstein product rule in `n4_boundary.py` is correct.
Independent instrumentation of all 16 cases, both curves and four objectives
found no depth-cap hits and endpoint-attained extrema. The classifications
are unchanged: connectivity wins throughout in 12 cases, buffering in three,
and one admits both signs. This audit covers this fixed grid and the frozen
conditional-iid model, not all ecological parameters.

## Interpretation corrections

Different third-moment coordinates for the two boundary roots establish
that a fixed third moment can admit both action rankings. They do not imply
non-monotonicity in the third moment: even an affine function increasing in
both coordinates can have a nonvertical zero set. The earlier negative
third-moment derivative witness uses `q=3/50`; the ranking-reversal case here
uses `q=3/100`. These are distinct results.

The largest baseline and connectivity persistence ranges in the 16-case
table both occur at `q=3/100, e0=1/50, e1=9/10, c=1/10`:
approximately `[0.251386, 0.325660]` for P0 and `[0.373902, 0.592553]`
for PC. The `q=3/50` case mentioned in message 390 has the largest PR range,
approximately `[0.231533, 0.299390]`.

## Equal triple-shock laws can still give opposite rankings

Set `n=4, T=50, p=1/10, q=3/100, e0=1/50, e1=9/10, c=1/10`.
Hold the third moment fixed at `m3=21/1000`. Two feasible hazard laws on
the lower and upper moment boundaries have:

| Boundary | m4 | PR | PC | PR - PC |
|---|---|---|---|---|
| Lower | 171/10000 | 0.4229425016010541 | 0.4286456813824586 | -0.005703179781404491 |
| Upper | 1389/70000 | 0.4199694197108412 | 0.41398744403192744 | +0.005981975678913741 |

Displayed outcomes are rounded; the signs were checked using exact rational
propagation. Feasibility follows from the constructive boundary laws in
`n4_geometry.md`. The lower law has atoms `(9-sqrt(57))/20` and
`(9+sqrt(57))/20`, with weights chosen to give mean `1/10`. The upper law
has atoms `{0,9/70,1}` and the corresponding upper-boundary weights.

For any labelled triple of patches, the probability of a particular hit
pattern with k hits is `E[theta^k (1-theta)^(3-k)]`. This is a polynomial
in moments zero through three. These two laws therefore agree on every
complete triple-shock distribution, as well as every marginal and pairwise
distribution, yet prefer opposite interventions. Their four-hit
probabilities differ. For any particular labelled triple pattern, the
probabilities for zero, one, two or three hits are respectively
`769/1000`, `61/1000`, `9/1000`, and `21/1000`; the one- and two-hit
values each apply to three patterns. This strengthens the pairwise-only witness within
the same model; it does not generalise the result to other mechanisms,
parameters, or field data. These matched distributions concern environmental
hits, not joint occupancy histories after mortality and recolonisation.

Reproduce without floating-point hazard atoms (standard-library Python):

```python
import sys
from fractions import Fraction as F
sys.path.insert(0, "ideation/b_screen")
import n4_vertical as V
p, q, x = F(1, 10), F(3, 100), F(21, 1000)
e0, e1, c = F(1, 50), F(9, 10), F(1, 10)
for boundary, yfun, expected in [("lower", V.Lfun, -1),
                                  ("upper", V.Ufun, 1)]:
    y = yfun(x, p, q)
    values = {}
    for action in ("R", "C"):
        parts = V.affine_parts(action, 4, e0, e1, c, p, q)
        K = V.kernel_at(*parts, x, y)
        values[action] = V.persist(K, V.start_state(action, 4))
    delta = values["R"] - values["C"]
    assert expected * delta > 0
    print(boundary, x, y, float(delta))
```

## Certificate reporting repair

The original boundary algorithm gives valid rational outer bounds, but its
maximum-depth fallback does not establish the promised tolerance. Its JSON
saves only rounded floating-point endpoints, and straddling outer bounds
alone do not prove a ranking reversal. The repair records exact outer and
attained bounds, their gaps and witnesses, explicitly checks the tolerance,
and requires attained opposite signs before declaring a reversal. Rounded
displays remain separate from exact certificates. These repairs preserve
the model and polynomial recurrence.

The repaired run regenerated all 16 cases. All 256 minimum/maximum
certificates (16 cases, two curves, four objectives, two extrema) have
exactly zero outer-to-attained gaps, with witnesses at `u=0` or `u=1`.
The class extrema are therefore attained by Glo or Ghi on this grid,
despite the failure of coordinate-wise monotonicity. All classifications
are unchanged; outward rounding shifts displayed bounds by at most
`1.1102230246251565e-16`. Exact root brackets and source/result hashes are
saved in the regenerated JSON. The original `n4_boundary.log` records
Fable's pre-repair run at `7d0cf94`; the regenerated JSON carries its own
source binding.

Reproduce the repaired run with
`uv run --python 3.12 --with gmpy2 --with numpy python ideation/b_screen/n4_boundary.py`.
The runtime recorded in its metadata is CPython 3.12.13 with gmpy2 2.3.1.
Four regression checks in `test_n4_boundary.py` cover depth exhaustion,
an attained interior extremum, outward rounding, and insufficient evidence
for a reversal. The exact witness above is independently executable.
