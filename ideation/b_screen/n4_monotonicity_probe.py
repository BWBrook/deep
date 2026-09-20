#!/usr/bin/env python3
"""Bounded n=4 moment-sensitivity probe for candidate B.

For fixed p=E(theta) and q=E(theta**2), a one-year transition entry is a
quartic in theta.  Interpolation at five rational nodes therefore gives the
exact moment decomposition

    K = A0 + p*A1 + q*A2 + x*A3 + y*A4,

where x=E(theta**3) and y=E(theta**4).  This script evaluates the two
coordinate sensitivities at the midpoint of the lower and upper two-point
moment vectors.  The initial pass is floating point so that any contrary
sign can be reported quickly; the reported sensitivities are then recomputed
with Fraction arithmetic.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction as F
from pathlib import Path


# Allow reuse of the exact transition kernels and endpoint moment laws in the
# existing n=3 certificate without modifying that file.
sys.path.insert(0, str(Path(__file__).parent))
from n3_certificate import (  # noqa: E402
    CS,
    E0S,
    E1S,
    P,
    QS,
    T,
    add_scaled,
    averaged_kernel,
    check_kernel,
    kernel_hom,
    kernel_buffered,
    laws,
    matvec,
    rat,
)


N = 4
THETA_NODES = (F(0), F(1, 4), F(1, 2), F(3, 4), F(1))
EXPECTED = {"P0": (-1, -1), "PR": (-1, -1), "PC": (-1, -1), "delta_PR_minus_PC": (1, 1)}


def solve_square(A, b):
    """Exact Gauss-Jordan solve for a small square Fraction system."""
    n = len(A)
    aug = [[F(x) for x in row] + [F(y)] for row, y in zip(A, b)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if aug[i][col])
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for i in range(n):
            if i == col:
                continue
            scale = aug[i][col]
            if scale:
                aug[i] = [x - scale * y for x, y in zip(aug[i], aug[col])]
    return [aug[i][-1] for i in range(n)]


def interpolation_inverse(nodes):
    V = [[th**j for j in range(len(nodes))] for th in nodes]
    return [solve_square(V, [F(int(i == row)) for i in range(len(nodes))]) for row in range(len(nodes))]


INTERP_INV = interpolation_inverse(THETA_NODES)


def interpolate_matrix(values):
    """Return [A0,...,A4] from exact matrices at THETA_NODES."""
    nrow, ncol = len(values[0]), len(values[0][0])
    coeff = [[[F(0) for _ in range(ncol)] for _ in range(nrow)] for _ in THETA_NODES]
    for i in range(nrow):
        for j in range(ncol):
            ys = [values[k][i][j] for k in range(len(values))]
            for degree in range(len(THETA_NODES)):
                # INTERP_INV[node][degree] is the degree coefficient of the
                # Lagrange basis function that is one at ``node``.
                coeff[degree][i][j] = sum(INTERP_INV[k][degree] * ys[k] for k in range(len(ys)))
    return coeff


def combine_coefficients(coeff, moments):
    out = [[F(0) for _ in row] for row in coeff[0]]
    for A, scale in zip(coeff, moments):
        out = add_scaled(out, scale, A)
    return out


def moment_vector(theta, weights):
    return tuple(sum(w * th**k for th, w in zip(theta, weights)) for k in range(5))


def moment_midpoint(p, q):
    (tl, wl), (tu, wu) = laws(p, q)
    glo = moment_vector(tl, wl)
    ghi = moment_vector(tu, wu)
    mid = tuple((a + b) / 2 for a, b in zip(glo, ghi))
    return glo, ghi, mid


def determinant(M):
    """Exact determinant by elimination (M is at most 3 by 3 here)."""
    A = [[F(x) for x in row] for row in M]
    out = F(1)
    for col in range(len(A)):
        pivot = next((i for i in range(col, len(A)) if A[i][col]), None)
        if pivot is None:
            return F(0)
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            out = -out
        out *= A[col][col]
        scale = A[col][col]
        for i in range(col + 1, len(A)):
            factor = A[i][col] / scale
            for j in range(col + 1, len(A)):
                A[i][j] -= factor * A[col][j]
    return out


def interior_checks(mom):
    _, p, q, x, y = mom
    H = ((F(1), p, q), (p, q, x), (q, x, y))
    L = ((p - q, q - x), (q - x, x - y))
    H_principal = (H[0][0], H[1][1], H[2][2], determinant(((H[0][0], H[0][1]), (H[1][0], H[1][1]))), determinant(((H[0][0], H[0][2]), (H[2][0], H[2][2]))), determinant(((H[1][1], H[1][2]), (H[2][1], H[2][2]))), determinant(H))
    L_principal = (L[0][0], L[1][1], determinant(L))
    return {
        "H_principal_minors_strictly_positive": all(v > 0 for v in H_principal),
        "L_principal_minors_strictly_positive": all(v > 0 for v in L_principal),
        "H_principal_minors": [rat(v) for v in H_principal],
        "L_principal_minors": [rat(v) for v in L_principal],
    }


def action_kernel_at_theta(action, e0, e1, c, theta):
    d = e0 + (e1 - e0) * theta
    if action == "R":
        return kernel_buffered(N, d, e0 + (d - e0) * F(1, 2), c)
    return kernel_hom(N, d, c * (2 if action == "C" else 1))


def matrix_equal(A, B):
    return A == B


def reconstruct_coefficients(coeff, theta):
    out = [[F(0) for _ in row] for row in coeff[0]]
    for degree, A in enumerate(coeff):
        out = add_scaled(out, theta**degree, A)
    return out


def float_matvec(v, K):
    return [sum(x * y for x, y in zip(v, col)) for col in zip(*K)]


def float_sensitivity(K, A, start, horizon=T):
    v = [0.0] * len(K)
    v[start] = 1.0
    s = [0.0] * len(K)
    for _ in range(horizon):
        s = [a + b for a, b in zip(float_matvec(s, K), float_matvec(v, A))]
        v = float_matvec(v, K)
    return sum(s) - s[0]


def exact_sensitivity(K, A, start, horizon=T):
    v = [F(0)] * len(K)
    v[start] = F(1)
    s = [F(0)] * len(K)
    for _ in range(horizon):
        s = [a + b for a, b in zip(matvec(s, K), matvec(v, A))]
        v = matvec(v, K)
    return sum(s) - s[0]


def sign(value):
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def action_matrices(action, e0, e1, c, moments, verify_endpoint_mix=True):
    samples = [action_kernel_at_theta(action, e0, e1, c, th) for th in THETA_NODES]
    for K in samples:
        check_kernel(K)
    coeff = interpolate_matrix(samples)
    interpolation_ok = all(matrix_equal(reconstruct_coefficients(coeff, th), K) for th, K in zip(THETA_NODES, samples))
    Kmid = combine_coefficients(coeff, moments)
    check_kernel(Kmid)
    if not verify_endpoint_mix:
        return coeff, Kmid, interpolation_ok
    # The midpoint is exactly the equal mixture of the two endpoint laws.
    (tl, wl), (tu, wu) = laws(P, moments[2])
    Klo = averaged_kernel(action, N, e0, e1, c, tl, wl)
    Khi = averaged_kernel(action, N, e0, e1, c, tu, wu)
    direct_mix = [[(a + b) / 2 for a, b in zip(ra, rb)] for ra, rb in zip(Klo, Khi)]
    return coeff, Kmid, interpolation_ok and matrix_equal(Kmid, direct_mix)


def counterexample():
    """Exact check of the supplied off-midpoint lens counterexample."""
    p, q = F(1, 10), F(3, 50)
    e0, e1, c = F(1, 50), F(9, 10), F(1, 10)
    x, y = F(29, 500), F(288411, 5000000)
    step = F(1, 10000000)
    moments = (F(1), p, q, x, y)
    actions = {}
    for action, start in (("R", 2 * N - 1), ("C", N)):
        coeff, _, interpolation_ok = action_matrices(action, e0, e1, c, moments, verify_endpoint_mix=False)
        assert interpolation_ok
        actions[action] = {"coeff": coeff, "start": start, "interpolation_exact": interpolation_ok}
    delta_derivative = exact_sensitivity(
        combine_coefficients(actions["R"]["coeff"], moments), actions["R"]["coeff"][3], actions["R"]["start"]
    ) - exact_sensitivity(
        combine_coefficients(actions["C"]["coeff"], moments), actions["C"]["coeff"][3], actions["C"]["start"]
    )
    vectors = []
    assert delta_derivative < 0
    for label, xx in (("x_minus_step", x - step), ("x_plus_step", x + step)):
        vec = (F(1), p, q, xx, y)
        checks = interior_checks(vec)
        assert checks["H_principal_minors_strictly_positive"]
        assert checks["L_principal_minors_strictly_positive"]
        outcomes = {}
        for action in ("R", "C"):
            K = combine_coefficients(actions[action]["coeff"], vec)
            check_kernel(K)
            outcomes[action] = persistence(K, actions[action]["start"])
        vectors.append({"label": label, "moments": [rat(z) for z in vec], "interior_checks": interior_checks(vec), "outcomes": {"PR": rat(outcomes["R"]), "PC": rat(outcomes["C"]), "delta_PR_minus_PC": rat(outcomes["R"] - outcomes["C"])}, "outcomes_float": {"PR": float(outcomes["R"]), "PC": float(outcomes["C"]), "delta_PR_minus_PC": float(outcomes["R"] - outcomes["C"])}})
    lens_lo = (x * x - 2 * p * q * x + q**3) / (q - p * p)
    lens_hi = x - (q - x) ** 2 / (p - q)
    assert lens_lo < y < lens_hi
    assert F(vectors[1]["outcomes"]["delta_PR_minus_PC"]) < F(vectors[0]["outcomes"]["delta_PR_minus_PC"])
    return {"parameters": {"n": N, "T": T, "p": rat(p), "q": rat(q), "e0": rat(e0), "e1": rat(e1), "c": rat(c)}, "base_moments": [rat(z) for z in moments], "base_interior_checks": interior_checks(moments), "lens_y_bounds": {"lower": rat(lens_lo), "upper": rat(lens_hi), "strictly_between": lens_lo < y < lens_hi}, "interpolation_exact_for_actions": all(v["interpolation_exact"] for v in actions.values()), "delta_m3_derivative": rat(delta_derivative), "delta_m3_derivative_float": float(delta_derivative), "step": rat(step), "vectors": vectors, "delta_decreases_with_x": F(vectors[1]["outcomes"]["delta_PR_minus_PC"]) < F(vectors[0]["outcomes"]["delta_PR_minus_PC"])}


def feasible_perturbation(mom, j, step=F(1, 100000)):
    """Return two feasible vectors differing only in x (j=3) or y (j=4)."""
    for scale in (step, step / 10, step / 100, step / 1000):
        lo = list(mom)
        hi = list(mom)
        lo[j] -= scale
        hi[j] += scale
        if all(interior_checks(tuple(z))["H_principal_minors_strictly_positive"] and interior_checks(tuple(z))["L_principal_minors_strictly_positive"] for z in (lo, hi)):
            return tuple(lo), tuple(hi), scale
    return None


def persistence(K, start):
    v = [F(0)] * len(K)
    v[start] = F(1)
    for _ in range(T):
        v = matvec(v, K)
    return sum(v) - v[0]


def main():
    cases = []
    failures = []
    decomposition_ok = []
    all_interiors = True
    for q in QS:
        glo, ghi, mom = moment_midpoint(P, q)
        checks = interior_checks(mom)
        all_interiors &= checks["H_principal_minors_strictly_positive"] and checks["L_principal_minors_strictly_positive"]
        for e0 in E0S:
            for e1 in E1S:
                for c in CS:
                    metrics_float = {}
                    metrics_exact = {}
                    for action, key, start in (("0", "P0", N), ("R", "PR", 2 * N - 1), ("C", "PC", N)):
                        coeff, Kmid, decomp = action_matrices(action, e0, e1, c, (F(1), P, q, mom[3], mom[4]))
                        decomposition_ok.append(decomp)
                        vals_float = [float_sensitivity([[float(x) for x in row] for row in Kmid], [[float(x) for x in row] for row in coeff[j]], start) for j in (3, 4)]
                        vals_exact = [exact_sensitivity(Kmid, coeff[j], start) for j in (3, 4)]
                        metrics_float[key] = {"dm3": vals_float[0], "dm4": vals_float[1]}
                        metrics_exact[key] = {"dm3": rat(vals_exact[0]), "dm4": rat(vals_exact[1]), "dm3_sign": sign(vals_exact[0]), "dm4_sign": sign(vals_exact[1])}
                        expected = EXPECTED[key]
                        for j, (vf, ve, exp) in enumerate(zip(vals_float, vals_exact, expected), start=3):
                            if (ve > 0) != (exp > 0) or ve == 0:
                                failures.append({"case": f"n4_q{rat(q)}_e0{rat(e0)}_e1{rat(e1)}_c{rat(c)}", "metric": key, "moment": f"m{j}", "float_derivative": vf, "exact_derivative": rat(ve), "expected_sign": "negative" if exp < 0 else "positive"})
                    delta_float = {m: metrics_float["PR"][m] - metrics_float["PC"][m] for m in ("dm3", "dm4")}
                    dpr = [F(metrics_exact["PR"][m]) for m in ("dm3", "dm4")]
                    dpc = [F(metrics_exact["PC"][m]) for m in ("dm3", "dm4")]
                    de = [a - b for a, b in zip(dpr, dpc)]
                    for j, (vf, ve, exp) in enumerate(zip((delta_float["dm3"], delta_float["dm4"]), de, (1, 1)), start=3):
                        if (ve > 0) != (exp > 0) or ve == 0:
                            failures.append({"case": f"n4_q{rat(q)}_e0{rat(e0)}_e1{rat(e1)}_c{rat(c)}", "metric": "delta_PR_minus_PC", "moment": f"m{j}", "float_derivative": vf, "exact_derivative": rat(ve), "expected_sign": "positive"})
                    case_id = f"n4_q{rat(q)}_e0{rat(e0)}_e1{rat(e1)}_c{rat(c)}"
                    cases.append({"id": case_id, "parameters": {"n": N, "T": T, "p": rat(P), "q": rat(q), "e0": rat(e0), "e1": rat(e1), "c": rat(c)}, "float_derivatives": {**metrics_float, "delta_PR_minus_PC": delta_float}, "exact_derivatives": {**metrics_exact, "delta_PR_minus_PC": {"dm3": rat(de[0]), "dm4": rat(de[1]), "dm3_sign": sign(de[0]), "dm4_sign": sign(de[1])}}})
    assert all_interiors and all(decomposition_ok)
    result = {
        "status": "complete",
        "conjecture_status": "refuted_by_exact_feasible_off_midpoint_counterexample",
        "method": "exact quartic interpolation and Fraction sensitivities at 16 midpoints, plus an exact off-midpoint counterexample found by a separate feasible-lens numerical scan",
        "parameters": {"n": N, "T": T, "p": rat(P), "q": [rat(q) for q in QS], "e0": [rat(x) for x in E0S], "e1": [rat(x) for x in E1S], "c": [rat(x) for x in CS]},
        "moment_midpoints": {},
        "validation": {"all_moment_midpoints_strictly_interior": all_interiors, "all_interpolations_and_direct_endpoint_mix_exact": all(decomposition_ok), "interpolation_nodes": [rat(x) for x in THETA_NODES]},
        "midpoint_failures": failures,
        "cases": cases,
        "off_midpoint_counterexample": counterexample(),
    }
    for q in QS:
        glo, ghi, mom = moment_midpoint(P, q)
        result["moment_midpoints"][rat(q)] = {"Glo": [rat(x) for x in glo], "Ghi": [rat(x) for x in ghi], "midpoint": [rat(x) for x in mom], "interior_checks": interior_checks(mom)}

    # If a contrary sign exists, attach a small exact feasible witness for the
    # clearest failure.  No witness is fabricated when the conjecture survives.
    if failures:
        first = failures[0]
        q = F(first["case"].split("_q")[1].split("_e0")[0])
        mom = moment_midpoint(P, q)[2]
        j = int(first["moment"][1:])
        witness = feasible_perturbation(mom, j)
        if witness:
            lo, hi, step = witness
            first["witness"] = {"step": rat(step), "lower": [rat(x) for x in lo], "upper": [rat(x) for x in hi], "lower_interior_checks": interior_checks(lo), "upper_interior_checks": interior_checks(hi)}
    out = Path(__file__).with_name("n4_monotonicity_probe_results.json")
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"wrote {out} ({len(cases)} midpoint cases; midpoint failures={len(failures)})")
    if failures:
        print("first failure:", failures[0])
    else:
        print("Midpoints satisfy the conjectured signs; the exact off-midpoint counterexample refutes the global conjecture.")


if __name__ == "__main__":
    main()
