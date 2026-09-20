#!/usr/bin/env python3
"""Exact n=3 Bernstein certificates for the DEEP candidate-B screen.

The two endpoint mixing laws have identical first two moments and differ only
in their third moment.  For n=3 every one-year kernel is cubic in theta, so
the identified family is represented exactly by K(lambda)=(1-lambda)Klo+
lambda Khi.  All arithmetic below is Fraction arithmetic; floating point is
used only for the compact summaries written to the JSON result.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb
import json
from pathlib import Path


T = 50
P = F(1, 10)
QS = (F(3, 100), F(6, 100))
E0S = (F(2, 100), F(1, 10))
E1S = (F(1, 2), F(9, 10))
CS = (F(2, 100), F(1, 10))
DEPTH = 12


def rat(x: F) -> str:
    """Stable rational serialization (including denominator 1)."""
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def pmf(k: int, m: int, pr: F) -> F:
    return F(comb(m, k)) * pr**k * (1 - pr) ** (m - k)


def laws(p: F, q: F):
    """Return the lower and upper two-point theta laws."""
    tl = (F(0), q / p)
    wl = (1 - p * p / q, p * p / q)
    w = (q - p * p) / (1 - 2 * p + q)
    a = (p - q) / (1 - p)
    tu = (a, F(1))
    wu = (1 - w, w)
    return (tl, wl), (tu, wu)


def d_of(theta: F, e0: F, e1: F) -> F:
    return e0 + (e1 - e0) * theta


def kernel_hom(n: int, d: F, c: F):
    """Exact occupied-count transition kernel for one fixed theta."""
    K = [[F(0) for _ in range(n + 1)] for _ in range(n + 1)]
    for k in range(n + 1):
        for m in range(k + 1):
            ps = pmf(m, k, 1 - d)
            r = 1 - (1 - c) ** m
            for j in range(n - m + 1):
                K[k][m + j] += ps * pmf(j, n - m, r)
    return K


def kernel_buffered(n: int, d: F, d_buffered: F, c: F):
    """Exact 2*n state (b,j) chain, b buffered occupancy and j ordinary."""
    K = [[F(0) for _ in range(2 * n)] for _ in range(2 * n)]
    for b in (0, 1):
        for j in range(n):
            i = b * n + j
            for s in range(j + 1):
                ps = pmf(s, j, 1 - d)
                bnext = (0,) if b == 0 else (0, 1)
                for bp in bnext:
                    pb = F(1) if b == 0 else (d_buffered if bp == 0 else 1 - d_buffered)
                    m = s + bp
                    r = 1 - (1 - c) ** m
                    for jj in range(n - 1 - s + 1):
                        pj = pmf(jj, n - 1 - s, r)
                        if bp:
                            K[i][n + s + jj] += ps * pb * pj
                        else:
                            K[i][s + jj] += ps * pb * pj * (1 - r)
                            K[i][n + s + jj] += ps * pb * pj * r
    return K


def add_scaled(A, scale, B):
    return [[a + scale * b for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def averaged_kernel(action: str, n: int, e0: F, e1: F, c: F, theta, weights):
    """Average fixed-theta kernels under one endpoint theta law."""
    d_buffered_e1 = e0 + F(1, 2) * (e1 - e0)
    K = None
    for th, wt in zip(theta, weights):
        d = d_of(th, e0, e1)
        if action == "R":
            Ki = kernel_buffered(n, d, d_of(th, e0, d_buffered_e1), c)
        else:
            Ki = kernel_hom(n, d, c * (2 if action == "C" else 1))
        K = [[wt * x for x in row] for row in Ki] if K is None else add_scaled(K, wt, Ki)
    return K


def matvec(v, K):
    return [sum(x * y for x, y in zip(v, col)) for col in zip(*K)]


def persistence_coefficients(Klo, Khi, start: int, absorbing: int = 0):
    """Bernstein coefficients of survival through T under K(lambda)."""
    # v[t][k] is a row-vector Bernstein coefficient at degree t.
    coeff = [[F(0) for _ in range(len(Klo))] for _ in range(1)]
    coeff[0][start] = F(1)
    for t in range(1, T + 1):
        prev = coeff
        cur = [[F(0) for _ in range(len(Klo))] for _ in range(t + 1)]
        for k in range(t + 1):
            if k < t:
                a = F(t - k, t)
                left = matvec(prev[k], Klo)
            else:
                a, left = F(0), None
            if k:
                b = F(k, t)
                right = matvec(prev[k - 1], Khi)
            else:
                b, right = F(0), None
            for j in range(len(Klo)):
                cur[k][j] = (a * left[j] if left is not None else F(0)) + (b * right[j] if right is not None else F(0))
        coeff = cur
    return [sum(row) - row[absorbing] for row in coeff]


def sign_of_coeffs(values):
    nonneg = all(x >= 0 for x in values)
    nonpos = all(x <= 0 for x in values)
    if nonneg and nonpos:
        return "zero"
    if nonneg and any(x > 0 for x in values):
        return "nonnegative"
    if nonpos and any(x < 0 for x in values):
        return "nonpositive"
    return "mixed"


def subdivide(values):
    """One exact de Casteljau subdivision at x=1/2."""
    levels = [list(values)]
    while len(levels[-1]) > 1:
        old = levels[-1]
        levels.append([(old[i] + old[i + 1]) / 2 for i in range(len(old) - 1)])
    left = [levels[i][0] for i in range(len(values))]
    right = [levels[-1 - i][-1] for i in range(len(values))]
    return left, right


def derivative_certificate(values):
    derivative = [T * (values[i + 1] - values[i]) for i in range(T)]
    direct = sign_of_coeffs(derivative)
    if direct != "mixed":
        return {"sign": direct, "method": "bernstein_coefficients", "depth": 0, "derivative": derivative}

    # Resolve mixed control points with exact dyadic de Casteljau boxes.  A
    # certificate is emitted only if every terminal box has the same sign.
    boxes = [derivative]
    terminal_signs = []
    for depth in range(1, DEPTH + 1):
        nxt = []
        for box in boxes:
            sg = sign_of_coeffs(box)
            if sg != "mixed":
                terminal_signs.append(sg)
            else:
                left, right = subdivide(box)
                nxt.extend((left, right))
        boxes = nxt
        if not boxes:
            break
        if not terminal_signs or len(set(terminal_signs)) > 1:
            # A sign change is already certified, or unresolved boxes remain.
            if len(set(terminal_signs)) > 1:
                return {"sign": "mixed", "method": "de_casteljau", "depth": depth, "derivative": derivative}
    terminal_signs.extend(sign_of_coeffs(box) for box in boxes)
    if terminal_signs and len(set(terminal_signs)) == 1 and terminal_signs[0] != "mixed":
        return {"sign": terminal_signs[0], "method": "de_casteljau", "depth": DEPTH, "derivative": derivative}
    return {"sign": "mixed", "method": "de_casteljau", "depth": DEPTH, "derivative": derivative}


def evaluate_bernstein(values, x: F) -> F:
    work = list(values)
    while len(work) > 1:
        work = [(1 - x) * work[i] + x * work[i + 1] for i in range(len(work) - 1)]
    return work[0]


def root_certificate(values, derivative_sign):
    """Exact monotone root interval, with width at most 2^-30."""
    lo, hi = F(0), F(1)
    flo, fhi = values[0], values[-1]
    if flo == 0:
        return {"status": "endpoint", "lambda": rat(lo), "lambda_float": 0.0, "value": rat(flo)}
    if fhi == 0:
        return {"status": "endpoint", "lambda": rat(hi), "lambda_float": 1.0, "value": rat(fhi)}
    if flo * fhi > 0:
        return {"status": "none"}
    # The derivative certificate establishes monotonicity; bisection therefore
    # retains a single guaranteed root in the sign-changing bracket.
    while hi - lo > F(1, 2**30):
        mid = (lo + hi) / 2
        fm = evaluate_bernstein(values, mid)
        if fm == 0:
            lo, hi = mid, mid
            break
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return {
        "status": "interval",
        "lambda_lower": rat(lo),
        "lambda_upper": rat(hi),
        "lambda_lower_float": float(lo),
        "lambda_upper_float": float(hi),
        "width": rat(hi - lo),
        "delta_lower": rat(flo),
        "delta_upper": rat(fhi),
        "delta_lower_float": float(flo),
        "delta_upper_float": float(fhi),
    }


def endpoint_values(values):
    return {
        "lambda_0": rat(values[0]),
        "lambda_1": rat(values[-1]),
        "lambda_0_float": float(values[0]),
        "lambda_1_float": float(values[-1]),
    }


def certified_bounds(values):
    lo, hi = min(values[0], values[-1]), max(values[0], values[-1])
    assert lo <= hi
    return {
        "min": rat(lo),
        "max": rat(hi),
        "min_float": float(lo),
        "max_float": float(hi),
    }


def metric_summary(values, is_delta=False):
    cert = derivative_certificate(values)
    sg = cert["sign"]
    out = {
        "degree": T,
        "derivative_sign": sg,
        "derivative_method": cert["method"],
        "subdivision_depth": cert["depth"],
        "endpoint_values": endpoint_values(values),
        "coefficient_sign": sign_of_coeffs(values),
        "certified": sg != "mixed",
        "derivative_control_min": rat(min(cert["derivative"])),
        "derivative_control_max": rat(max(cert["derivative"])),
    }
    if sg != "mixed":
        out["monotonicity"] = {"nonnegative": "increasing", "nonpositive": "decreasing", "zero": "constant"}[sg]
        out["certified_bounds"] = certified_bounds(values)
    else:
        out["monotonicity"] = "unresolved"
        out["certified_bounds"] = None
    if is_delta:
        out["root_lambda"] = root_certificate(values, sg) if sg != "mixed" else {"status": "unresolved"}
    return out


def check_kernel(K):
    assert all(x >= 0 for row in K for x in row)
    assert all(sum(row) == 1 for row in K)
    assert K[0][0] == 1


def m3_interval(p: F, q: F):
    (tl, wl), (tu, wu) = laws(p, q)
    lo = sum(w * th**3 for th, w in zip(tl, wl))
    hi = sum(w * th**3 for th, w in zip(tu, wu))
    return {"lower": rat(lo), "upper": rat(hi), "lower_float": float(lo), "upper_float": float(hi)}


def m3_bounds(p: F, q: F):
    (tl, wl), (tu, wu) = laws(p, q)
    return (sum(w * th**3 for th, w in zip(tl, wl)),
            sum(w * th**3 for th, w in zip(tu, wu)))


def validate_n2(qs):
    """n=2 endpoint kernels are identical because only m1,m2 enter."""
    for q in qs:
        (tl, wl), (tu, wu) = laws(P, q)
        for e0 in E0S:
            for e1 in E1S:
                for c in CS:
                    for action in ("0", "R", "C"):
                        lo = averaged_kernel(action, 2, e0, e1, c, tl, wl)
                        hi = averaged_kernel(action, 2, e0, e1, c, tu, wu)
                        assert lo == hi, (q, e0, e1, c, action)


def make_results():
    validate_n2(QS)
    cases = []
    for q in QS:
        (tl, wl), (tu, wu) = laws(P, q)
        for e0 in E0S:
            for e1 in E1S:
                for c in CS:
                    kernels = {}
                    for action in ("0", "R", "C"):
                        klo = averaged_kernel(action, 3, e0, e1, c, tl, wl)
                        khi = averaged_kernel(action, 3, e0, e1, c, tu, wu)
                        check_kernel(klo)
                        check_kernel(khi)
                        start = 5 if action == "R" else 3
                        kernels[action] = persistence_coefficients(klo, khi, start)
                    delta = [a - b for a, b in zip(kernels["R"], kernels["C"])]
                    metrics = {
                        "P0": metric_summary(kernels["0"]),
                        "PR": metric_summary(kernels["R"]),
                        "PC": metric_summary(kernels["C"]),
                        "delta_PR_minus_PC": metric_summary(delta, is_delta=True),
                    }
                    root = metrics["delta_PR_minus_PC"]["root_lambda"]
                    if root["status"] == "interval":
                        m3lo, m3hi = m3_bounds(P, q)
                        rl, rh = F(root["lambda_lower"]), F(root["lambda_upper"])
                        ml, mh = m3lo + rl * (m3hi - m3lo), m3lo + rh * (m3hi - m3lo)
                        root["m3_interval"] = {
                            "lower": rat(ml), "upper": rat(mh),
                            "lower_float": float(ml), "upper_float": float(mh),
                        }
                    elif root["status"] == "endpoint":
                        m3lo, m3hi = m3_bounds(P, q)
                        rv = F(root["lambda"])
                        mv = m3lo + rv * (m3hi - m3lo)
                        root["m3"] = {"value": rat(mv), "value_float": float(mv)}
                    cases.append({
                        "id": f"n3_q{rat(q)}_e0{rat(e0)}_e1{rat(e1)}_c{rat(c)}",
                        "parameters": {"n": 3, "T": T, "p": rat(P), "q": rat(q), "e0": rat(e0), "e1": rat(e1), "c": rat(c)},
                        "metrics": metrics,
                    })
    certified = all(m["certified"] for case in cases for m in case["metrics"].values())
    return {
        "method": "exact rational endpoint kernels and Bernstein coefficients; dyadic exact de Casteljau when needed",
        "status": "complete" if certified else "complete_with_unresolved_certificates",
        "parameters": {"n": 3, "T": T, "p": rat(P), "q": [rat(x) for x in QS], "e0": [rat(x) for x in E0S], "e1": [rat(x) for x in E1S], "c": [rat(x) for x in CS]},
        "m3intervals": {rat(q): m3_interval(P, q) for q in QS},
        "validation": {"row_stochastic_and_nonnegative": True, "n2_endpoint_kernel_invariance": True},
        "cases": cases,
    }


def main():
    path = Path(__file__).with_name("n3_certificate_results.json")
    result = make_results()
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"wrote {path} ({len(result['cases'])} cases; status={result['status']})")
    for case in result["cases"]:
        d = case["metrics"]["delta_PR_minus_PC"]
        print(case["id"], d["derivative_sign"], d["coefficient_sign"], d["root_lambda"]["status"])


if __name__ == "__main__":
    main()
