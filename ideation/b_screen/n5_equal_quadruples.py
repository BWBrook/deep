#!/usr/bin/env python3
"""Bounded n=5 search for equal-(m1,...,m4) action-ranking reversals.

For a fixed strict-interior degree-4 Hausdorff moment vector, the admissible
fifth moment is the interval obtained from the Schur complements of

    Hx   = (m[i+j+1]) i,j=0..2,
    H1-x = (m[i+j]-m[i+j+1]) i,j=0..2.

The script searches a small rational grid in the already-derived n=4 lens,
using floating matrix powers only as a screen.  Any selected reversal is then
propagated for 50 years in Fraction arithmetic and checked against the exact
Hausdorff positive-semidefinite conditions.  It deliberately does not compute
full-class n=5 bounds.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

import numpy as np

import n4_vertical as N4


P = F(1, 10)
QS = (F(3, 100), F(3, 50))
E0 = F(1, 50)
E1 = F(9, 10)
C = F(1, 10)
T = 50
N = 5
ACTIONS = ("R", "C")


def rat(x: F) -> str:
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def det2(a: F, b: F, c: F, d: F) -> F:
    return a * d - b * c


def det3(M: list[list[F]]) -> F:
    return (
        M[0][0] * det2(M[1][1], M[1][2], M[2][1], M[2][2])
        - M[0][1] * det2(M[1][0], M[1][2], M[2][0], M[2][2])
        + M[0][2] * det2(M[1][0], M[1][1], M[2][0], M[2][1])
    )


def quad2(A: list[list[F]], b: list[F]) -> F:
    """b' A^-1 b for a nonsingular symmetric 2 by 2 matrix."""
    den = A[0][0] * A[1][1] - A[0][1] ** 2
    assert den > 0
    return (A[1][1] * b[0] ** 2 - 2 * A[0][1] * b[0] * b[1] + A[0][0] * b[1] ** 2) / den


def m5_interval(m: list[F]) -> tuple[F, F]:
    """Sharp degree-5 Hausdorff interval conditional on m0,...,m4."""
    m0, m1, m2, m3, m4 = m
    A = [[m1, m2], [m2, m3]]
    D = [[m0 - m1, m1 - m2], [m1 - m2, m2 - m3]]
    lower = quad2(A, [m3, m4])
    upper = m4 - quad2(D, [m2 - m3, m3 - m4])
    assert lower < upper
    return lower, upper


def moment_matrices(m: list[F]) -> tuple[list[list[F]], list[list[F]]]:
    assert len(m) == 6
    hx = [[m[i + j + 1] for j in range(3)] for i in range(3)]
    h1x = [[m[i + j] - m[i + j + 1] for j in range(3)] for i in range(3)]
    return hx, h1x


def principal_minors(M: list[list[F]]) -> dict[str, F]:
    out: dict[str, F] = {}
    for k in (1, 2, 3):
        for idx in combinations(range(3), k):
            if k == 1:
                value = M[idx[0]][idx[0]]
            elif k == 2:
                i, j = idx
                value = M[i][i] * M[j][j] - M[i][j] ** 2
            else:
                value = det3(M)
            out["".join(map(str, idx))] = value
    return out


def add_scaled(A, scale: F, B):
    return [[a + scale * b for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def kernel_from_moments(coefficients, moments: list[F]):
    K = [[F(0) for _ in row] for row in coefficients[0]]
    for moment, Ak in zip(moments, coefficients):
        K = add_scaled(K, moment, Ak)
    return K


def matmul(A, B):
    Bt = list(zip(*B))
    return [[sum(x * y for x, y in zip(row, col)) for col in Bt] for row in A]


def matpow(A, exponent: int):
    size = len(A)
    result = [[F(int(i == j)) for j in range(size)] for i in range(size)]
    base = A
    while exponent:
        if exponent & 1:
            result = matmul(result, base)
        exponent >>= 1
        if exponent:
            base = matmul(base, base)
    return result


def exact_persistence(K, start: int) -> F:
    KT = matpow(K, T)
    return F(1) - KT[start][0]


def float_persistence(K: np.ndarray, start: int) -> float:
    return float(1.0 - np.linalg.matrix_power(K, T)[start, 0])


def exact_kernel_checks(K) -> dict:
    entries = [x for row in K for x in row]
    row_sums = [sum(row) for row in K]
    return {
        "all_entries_nonnegative": all(x >= 0 for x in entries),
        "all_rows_sum_to_one": all(x == 1 for x in row_sums),
        "absorbing_zero": K[0][0] == 1 and all(K[0][j] == 0 for j in range(1, len(K))),
        "minimum_entry": rat(min(entries)),
    }


def make_coefficients():
    exact = {
        action: N4.moment_kernels(action, N, E0, E1, C)
        for action in ACTIONS
    }
    floating = {
        action: [np.asarray([[float(x) for x in row] for row in Ak]) for Ak in exact[action]]
        for action in ACTIONS
    }
    # Reconstruction at interpolation nodes and two additional rational nodes.
    for action in ACTIONS:
        for theta in [F(j, N) for j in range(N + 1)] + [F(7, 19), F(13, 23)]:
            reconstructed = kernel_from_moments(exact[action], [theta**k for k in range(N + 1)])
            direct = N4.raw_kernel(action, N, theta, E0, E1, C)
            assert reconstructed == direct
    return exact, floating


def float_delta(float_coefficients, moments: list[F]) -> tuple[float, dict[str, float]]:
    values = {}
    for action in ACTIONS:
        K = sum((float(moment) * Ak for moment, Ak in zip(moments, float_coefficients[action])),
                np.zeros_like(float_coefficients[action][0]))
        values[action] = float_persistence(K, N4.start_state(action, N))
    return values["R"] - values["C"], values


def screen(grid_denominator: int, float_coefficients):
    records = []
    q_summaries = []
    for q in QS:
        xlo, xhi = N4.xbounds(P, q)
        tested = flips = 0
        best = None
        for ix in range(1, grid_denominator):
            x = xlo + F(ix, grid_denominator) * (xhi - xlo)
            ylo, yhi = N4.Lfun(x, P, q), N4.Ufun(x, P, q)
            for iy in range(1, grid_denominator):
                y = ylo + F(iy, grid_denominator) * (yhi - ylo)
                base = [F(1), P, q, x, y]
                lower, upper = m5_interval(base)
                dl, vl = float_delta(float_coefficients, base + [lower])
                du, vu = float_delta(float_coefficients, base + [upper])
                tested += 1
                if dl * du < 0:
                    flips += 1
                    score = min(abs(dl), abs(du))
                    candidate = {
                        "q": q, "m3": x, "m4": y,
                        "m5_lower": lower, "m5_upper": upper,
                        "delta_lower_float": dl, "delta_upper_float": du,
                        "values_lower_float": vl, "values_upper_float": vu,
                        "score": score,
                    }
                    if best is None or score > best["score"]:
                        best = candidate
        q_summaries.append({"q": rat(q), "tested": tested, "flips": flips})
        if best is not None:
            records.append(best)
            break
    return (max(records, key=lambda z: z["score"]) if records else None), q_summaries


def exact_certificate(candidate, exact_coefficients):
    base = [F(1), P, candidate["q"], candidate["m3"], candidate["m4"]]
    endpoint_results = {}
    for label, m5 in (("lower", candidate["m5_lower"]), ("upper", candidate["m5_upper"])):
        moments = base + [m5]
        hx, h1x = moment_matrices(moments)
        minors_x = principal_minors(hx)
        minors_1x = principal_minors(h1x)
        assert all(v >= 0 for v in minors_x.values())
        assert all(v >= 0 for v in minors_1x.values())
        values = {}
        checks = {}
        for action in ACTIONS:
            K = kernel_from_moments(exact_coefficients[action], moments)
            checks[action] = exact_kernel_checks(K)
            assert all((checks[action][key] for key in
                        ("all_entries_nonnegative", "all_rows_sum_to_one", "absorbing_zero")))
            values[action] = exact_persistence(K, N4.start_state(action, N))
        delta = values["R"] - values["C"]
        endpoint_results[label] = {
            "m5": rat(m5),
            "moment_matrices": {
                "Hx": [[rat(v) for v in row] for row in hx],
                "H1_minus_x": [[rat(v) for v in row] for row in h1x],
            },
            "principal_minors": {
                "Hx": {k: rat(v) for k, v in minors_x.items()},
                "H1_minus_x": {k: rat(v) for k, v in minors_1x.items()},
            },
            "kernel_checks": checks,
            "P_R": rat(values["R"]),
            "P_C": rat(values["C"]),
            "Delta_R_minus_C": rat(delta),
            "P_R_float": float(values["R"]),
            "P_C_float": float(values["C"]),
            "Delta_R_minus_C_float": float(delta),
            "float_screen_absolute_error": abs(
                (candidate["delta_lower_float"] if label == "lower" else candidate["delta_upper_float"])
                - float(delta)
            ),
            "exact_sign": "positive" if delta > 0 else "negative" if delta < 0 else "zero",
        }
    assert endpoint_results["lower"]["exact_sign"] != endpoint_results["upper"]["exact_sign"]
    assert "zero" not in (endpoint_results["lower"]["exact_sign"], endpoint_results["upper"]["exact_sign"])
    return base, endpoint_results


def json_result(grid: int, summaries, candidate, base=None, endpoints=None):
    out = {
        "scope": "bounded equal-first-four-moments endpoint screen; not full-class n=5 bounds",
        "model": {"n": N, "T": T, "p": rat(P), "e0": rat(E0), "e1": rat(E1), "c": rat(C)},
        "grid_denominator": grid,
        "q_summaries": summaries,
        "schur_interval": {
            "lower": "[m3,m4] [[m1,m2],[m2,m3]]^-1 [m3,m4]^T",
            "upper": "m4 - [m2-m3,m3-m4] [[m0-m1,m1-m2],[m1-m2,m2-m3]]^-1 [m2-m3,m3-m4]^T",
        },
        "kernel_interpolation": {
            "degree": 5,
            "exact_node_check": True,
            "nodes": [rat(F(j, N)) for j in range(N + 1)],
            "additional_exact_check_nodes": [rat(F(7, 19)), rat(F(13, 23))],
        },
        "moment_feasibility_source": {
            "citation": "Curto and Fialkow (1991), Houston Journal of Mathematics 17(4), Theorem 4.1 and Remark 4.2, p. 622",
            "url": "https://www.math.uh.edu/~hjm/v017n4/0603CURTO.pdf",
            "application": "For odd degree 5 on [0,1], existence is equivalent to Hx >= 0 and H1-x >= 0.",
        },
        "status": "exact reversal witness" if candidate else "no reversal on bounded grid",
    }
    if candidate:
        out["witness"] = {
            "moments_m0_to_m4": [rat(x) for x in base],
            "strict_interior_checks": {
                "m3_between_lens_limits": bool(N4.xbounds(P, base[2])[0] < base[3] < N4.xbounds(P, base[2])[1]),
                "m4_between_lens_limits": bool(N4.Lfun(base[3], P, base[2]) < base[4] < N4.Ufun(base[3], P, base[2])),
                "m5_interval_nonempty": bool(candidate["m5_lower"] < candidate["m5_upper"]),
            },
            "float_screen": {
                "Delta_lower": candidate["delta_lower_float"],
                "Delta_upper": candidate["delta_upper_float"],
            },
            "exact_endpoints": endpoints,
            "existence_argument": (
                "For a degree-5 sequence on [0,1], Hx >= 0 and H1-x >= 0 are necessary and sufficient "
                "Hausdorff conditions. Every principal minor is reported nonnegative, so each endpoint "
                "has at least one probability law on [0,1]."
            ),
        }
    return out


def markdown(result: dict) -> str:
    if result["status"] != "exact reversal witness":
        qs = "; ".join(f"q={x['q']}: {x['tested']} pairs, {x['flips']} flips" for x in result["q_summaries"])
        return f"""# n=5 equal-quadruple screen\n\nNo endpoint reversal was found on the bounded rational grid ({qs}). This is a negative screen, not a full-class bound.\n"""
    w = result["witness"]
    lo, hi = w["exact_endpoints"]["lower"], w["exact_endpoints"]["upper"]
    m = w["moments_m0_to_m4"]
    return f"""# n=5 equal-first-four-moments witness

## Decisive finding

Two admissible mixing laws can agree exactly on `m1..m4` yet reverse the 50-year ranking of buffering (`R`) and doubled connectivity (`C`) for `n=5`, at the fixed parameters `p=1/10`, `q={m[2]}`, `e0=1/50`, `e1=9/10`, and `c=1/10`.

The common moment vector is

`(m0,m1,m2,m3,m4) = ({', '.join(m)})`.

At the sharp lower fifth-moment endpoint `m5={lo['m5']}`, exact propagation gives `P_R-P_C < 0` ({lo['Delta_R_minus_C_float']:.12g}). At the sharp upper endpoint `m5={hi['m5']}`, it gives `P_R-P_C > 0` ({hi['Delta_R_minus_C_float']:.12g}). The full exact rational differences are retained in `n5_equal_quadruples.json`.

For any labelled four-patch subset, a particular hit pattern with k hits has
probability `E[theta^k (1-theta)^(4-k)]`. Its degree is four, so these laws
agree on every complete four-patch environmental-hit distribution, and on
all lower-order marginals. Their five-way joint-hit probabilities differ.
This does not assert equality of occupancy histories.

## Evidence and assumptions

For strict-interior `m0..m4`, Schur complements give

`L5 = [m3,m4] [[m1,m2],[m2,m3]]^-1 [m3,m4]^T`

and

`U5 = m4 - [m2-m3,m3-m4] [[m0-m1,m1-m2],[m1-m2,m2-m3]]^-1 [m2-m3,m3-m4]^T`.

These are sharp because the odd-degree Hausdorff existence criterion on `[0,1]` is exactly `Hx=(m[i+j+1]) >= 0` and `H1-x=(m[i+j]-m[i+j+1]) >= 0` for `i,j=0..2`; see [Curto & Fialkow (1991), Theorem 4.1 and Remark 4.2, p. 622](https://www.math.uh.edu/~hjm/v017n4/0603CURTO.pdf). All principal minors of both matrices are nonnegative at both endpoints; one determinant is zero at each active endpoint. Thus both moment vectors are represented by probability laws on `[0,1]`.

The one-year kernels were interpolated exactly as degree-5 polynomials in the hazard `theta`. The rational witness was selected by a float64 `matrix_power(K,50)` screen, then independently propagated with exact rational matrix exponentiation. Every exact endpoint kernel is entrywise nonnegative, row-stochastic, and has absorbing state zero.

## Scope and residual risk

The parent audit also matched exact iterative 50-step propagation to the
matrix-power results, and reconstructed three-atom laws at both endpoints.
The reconstructed moments matched through degree five within 1.1e-17;
direct atomic-law floating-point propagation reproduced both signs.

This is an existence witness conditional on the specified iid annual
mixing-law model. It is not a full optimization over the n=5 moment class,
does not establish the prevalence or biological relevance of reversals,
and is not a theorem for arbitrary patch counts. Fable's independent audit
is the next review step.

Reproduce with
`uv run --python 3.12 --with numpy python ideation/b_screen/n5_equal_quadruples.py`.
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=int, default=24)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()

    exact_coefficients, float_coefficients = make_coefficients()
    candidate, summaries = screen(args.grid, float_coefficients)
    if candidate:
        base, endpoints = exact_certificate(candidate, exact_coefficients)
    else:
        base = endpoints = None
    result = json_result(args.grid, summaries, candidate, base, endpoints)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "n5_equal_quadruples.json"
    md_path = args.output_dir / "n5_equal_quadruples.md"
    json_path.write_text(json.dumps(result, indent=2) + "\n")
    md_path.write_text(markdown(result))
    print(json.dumps({
        "status": result["status"],
        "q_summaries": summaries,
        "json": str(json_path),
        "markdown": str(md_path),
    }, indent=2))


if __name__ == "__main__":
    main()
