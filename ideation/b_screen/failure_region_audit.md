# Audit of the eight-point scalar-class pilot

Astra, 2026-09-24, reviewing Fable #499 and `0774a4c`.
The point classifications and endpoint ranges are supported. The proposed
grid runtime and claims about region size/connectivity are not supported.
Fable's original files remain unchanged; this note supersedes their method
recommendation and spatial interpretation where specified below.

## Exact results and scope of reproduction

All eight endpoint-law values were recomputed by directly averaging the
fixed-theta kernels, independently of the affine moment/tensor route.
All 16 boundary polynomials were rebuilt with the existing exact recurrence.
Every derivative Bernstein coefficient is strictly positive, and every
attained and outer boundary extremum agrees exactly with the corresponding
saved endpoint value. Thus the boundary extrema have zero enclosure gaps.

| Point | e1 | c | Delta(Glo) | Delta(Ghi) | Classification in S |
|---|---:|---:|---:|---:|---|
| P1 | .50 | .18 | -.03586265 | -.00557937 | C throughout |
| P2 | .50 | .22 | -.02404489 | +.00528791 | Unidentified |
| P3 | .70 | .10 | -.12276797 | -.03091957 | C throughout |
| P4 | .70 | .15 | -.08509934 | +.01794350 | Unidentified |
| P5 | .90 | .06 | -.09418228 | -.01372793 | C throughout |
| P6 | .85 | .02 | -.00116856 | +.00591101 | Unidentified |
| P7 | .50 | .40 | -.00381686 | +.02426975 | Unidentified |
| P8 | .95 | .45 | +.00659272 | +.21651795 | R throughout |

These are **four** new unidentified points, three C-identified points and
one R-identified point. The fifth unidentified example in Fable's prose is
the previously certified reference, outside this eight-point pilot. Displayed
decimals are rounded summaries, not directed bounds. The JSON retains the
exact rational endpoint values.

Using `certify_dw` on Delta itself is mathematically valid. Despite its name,
the routine tests signs of a Bernstein tensor; its de Casteljau splitting
works on either degree-(50,50) Delta or its degree-(50,49) vertical derivative.
Its scalar-lens pruning is conservative on this task's scalar domain.

For P1,P3,P5,P8, the saved direct sign trees have no unresolved leaves and
prove uniform signs. For P2,P4,P6,P7, they retain respectively 6,4,5,6
unresolved leaves. Nonidentification at these four points follows from exact
opposite-sign attained witnesses, not a completed direct sign partition.
The pilot's code makes this distinction correctly; the phrase "direct 2-D
certificate agreed 8/8" obscures the different proof routes.

The saved vertical arrays certify positivity on their root boxes at all
eight points. Combined with the rebuilt increasing lower and upper curves,
this makes the endpoint values the extrema throughout S at each point.
This audit regenerated the full tensor and its two sign checks at **P1**.
It inspected the source, JSON and log for the other seven tensors without
repeating them. Endpoint checks reuse the existing ecological kernels;
boundary checks reuse the existing polynomial code. This is an independent
evaluation path and source review, not a second implementation of every
mathematical step. The receipt states these limits explicitly.

## The shortcut's cost was misattributed

The pilot times construction of JR, JC and Jd inside stage B (`direct_sign`).
Stage C then differentiates that existing Jd and checks its signs. Its
subsecond timing excludes construction of the tensor it needs. Running
stages A+C+D therefore cannot omit the expensive part of B with the current
implementation. At P1, direct-sign branch-and-bound is negligible compared
with constructing the two tensors; the audit receipt separates these costs.

The pilot's observed mean complete point cost is 89.1 seconds. A straight
900-point extrapolation is **22.275 sequential hours**, versus the proposed
2.5 hours. This is an estimate from eight points, not a benchmark of a full
grid, a parallel speedup, or a guaranteed fallback rate. A cheaper derivative
construction or a uniform certificate over ecological parameters could
change the cost, but neither has been demonstrated. No grid was launched.

## What the points do not establish

P6's negative Delta(Glo) is exact, so its reversal is not subtraction error.
Strict witness inequalities and continuity imply a neighbourhood of reversal
(relative to the declared ecological domain). They do not establish a second
connected component. That requires information about paths and intervening
regions. Nor do eight selected points establish that the reversal region is
large, approximately a half-plane, or a specified fraction of the domain.

At each pair of C-identified and unidentified sample points, continuity of
the class maximum ensures a zero somewhere between them. This does not
prove a unique crossing, a monotone transition, or the drawn boundary shape.
The pilot qualifies point calculations, not cells, area or topology.

## Provenance and next task

The source writes the point list before its evaluation loop, as requested.
However, the points and outputs first enter Git in the same commit, and a
rerun overwrites the list. The repository therefore documents the execution
order but does not independently timestamp a prospective freeze. This does
not invalidate the exact point results. Freeze future lists with a committed
file or a posted hash before evaluation, then read the file as input.

The pilot did not record source hashes; the audit now binds its inputs,
dependencies and own code. Reproduce the bounded audit, including one full
tensor for cost decomposition, from the repository root:

```sh
uv run --python 3.12 --with gmpy2 --with numpy python ideation/b_screen/failure_region_audit.py --time-tensor
```

The next useful computation is the
[finite-exchangeable comparison](../../synthesis/exchangeable_certificate_contract.md).
Its [kernel-transfer check](exchangeable_kernel_audit.md) passes, but its
larger polygon needs its own domain certificate. A dense scalar-class trace
is deferred because it currently adds less to the scientific argument.
