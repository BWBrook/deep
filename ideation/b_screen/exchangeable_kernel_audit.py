#!/usr/bin/env python3
"""Check the scalar moment kernels on the entire finite-exchangeable simplex.

For each ecological parameter combination in the original design, construct
the independent labelled kernel for each of the five fixed hit counts, verify
strong lumpability, and compare with the moment kernel using factorial hit
moments. These five laws span every exchangeable four-patch hit law, including
laws outside scalar mixtures. This checks kernel transfer, not class extrema.

uv run --python 3.12 --with gmpy2 --with numpy --with scipy \
  python ideation/b_screen/exchangeable_kernel_audit.py
"""
from __future__ import annotations

import hashlib
import itertools
import json
import platform
import time
from fractions import Fraction as F
from math import comb
from pathlib import Path

import gmpy2
from gmpy2 import mpq as Q

import labelled_class_audit as labelled
import n4_vertical as moment
import n3_certificate as scalar
import n4_boundary as boundary

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "exchangeable_kernel_audit.json"
N = 4
S = 1 << N


def state_group(state, action):
    # R retains the buffered patch (bit 0) separately from the other patches.
    return (state & 1) * N + (state >> 1).bit_count() if action == "R" else state.bit_count()


def collapse(kernel, action):
    size = 2 * N if action == "R" else N + 1
    rows = [None] * size
    checks = 0
    for old in range(S):
        row = [Q(0)] * size
        for dest in range(S):
            row[state_group(dest, action)] += kernel[old][dest]
        group = state_group(old, action)
        if rows[group] is None:
            rows[group] = row
        else:
            assert rows[group] == row, ("strong lumpability failed", action, old)
            checks += size
    assert all(row is not None for row in rows)
    return rows, checks


def kernel_for_count(action, count, e0, e1, c):
    out = [[Q(0)] * S for _ in range(S)]
    patterns = [x for x in labelled.PATTERNS if sum(x) == count]
    assert len(patterns) == comb(N, count)
    for hit in patterns:
        death = tuple(e0 + (e1-e0) * hit[i] * (Q(1, 2) if action == "R" and i == 0 else 1)
                      for i in range(N))
        kernel = labelled.independent_kernel(death, c * (2 if action == "C" else 1))
        for old in range(S):
            for dest in range(S):
                out[old][dest] += kernel[old][dest] / len(patterns)
    labelled.check_kernel(out)
    return out


def edge_checks():
    """Qualify degree-T boundary polynomials at the two reference E polygons.

    These establish extrema on the boundary only. An interior/domain proof
    is still required before treating them as full-class extrema.
    """
    polygons = [
        (F(3, 100), [(F(0), F(0)), (F(3, 200), F(0)), (F(3, 100), F(3, 100))]),
        (F(3, 50), [(F(2, 75), F(0)), (F(3, 100), F(0)), (F(3, 50), F(3, 50)), (F(1, 25), F(1, 25))]),
    ]
    output = []
    for q, vertices in polygons:
        parts = {a: moment.affine_parts(a, N, F(1,50), F(9,10), F(1,10), moment.P, q) for a in ("R", "C")}
        edges = []
        for i, pt0 in enumerate(vertices):
            pt1 = vertices[(i+1) % len(vertices)]
            cs = {}
            for action in ("R", "C"):
                kernels = [[[Q(v) for v in row] for row in moment.kernel_at(*parts[action], *pt)]
                           for pt in (pt0, pt1)]
                cs[action] = scalar.persistence_coefficients(*kernels, moment.start_state(action, N))
                assert len(cs[action]) == moment.T+1
                for u in (F(1,4), F(3,4)):
                    pt = [(1-u)*a+u*b for a,b in zip(pt0,pt1)]
                    kernel = [[Q(v) for v in row] for row in moment.kernel_at(*parts[action], *pt)]
                    assert boundary.evaluate(cs[action],u) == moment.persist(kernel,moment.start_state(action,N))
            delta = [r-c for r,c in zip(cs["R"],cs["C"])]
            derivative = [moment.T*(delta[k+1]-delta[k]) for k in range(moment.T)]
            sign = "increasing" if min(derivative)>0 else "decreasing" if max(derivative)<0 else "unresolved"
            assert sign != "unresolved"
            edges.append({
                "from": list(map(str,pt0)), "to": list(map(str,pt1)), "degree": moment.T,
                "Delta_endpoints": [labelled.exact_value(delta[0]), labelled.exact_value(delta[-1])],
                "derivative_coefficient_bounds": [labelled.exact_value(min(derivative)), labelled.exact_value(max(derivative))],
                "strict_monotonicity": sign,
                "interior_exact_cross_checks": {"u": ["1/4","3/4"], "both_actions_match": True},
                "coefficient_sha256": hashlib.sha256("\n".join(map(str,delta)).encode()).hexdigest(),
            })
        output.append({"q": str(q), "e0": "1/50", "e1": "9/10", "c": "1/10", "edges": edges,
                       "scope": "Boundary extrema only; full-polygon vertical derivative not evaluated."})
    return output


def main():
    start = time.perf_counter()
    assert moment.T == scalar.T == 50
    cases = []
    entries = 0
    lumpability_checks = 0
    for e0, e1, c in itertools.product(moment.E0S, moment.E1S, moment.CS):
        case = {"e0": str(e0), "e1": str(e1), "c": str(c), "actions": {}}
        for action in ("R", "C"):
            matrices = moment.moment_kernels(action, N, e0, e1, c)
            size = len(matrices[0])
            for count in range(N + 1):
                direct = kernel_for_count(action, count, Q(e0), Q(e1), Q(c))
                reduced, checks = collapse(direct, action)
                lumpability_checks += checks
                # E[prod of j specified hits] under uniform count = C(count,j)/C(n,j).
                tau = [F(comb(count, j), comb(N, j)) if j <= count else F(0) for j in range(N+1)]
                expected = [[sum(tau[j] * matrices[j][r][s] for j in range(N+1))
                             for s in range(size)] for r in range(size)]
                assert reduced == expected, (e0, e1, c, action, count)
                entries += size * size
            case["actions"][action] = {"fixed_count_laws_checked": list(range(N+1)),
                                       "strong_lumpability": True, "exact_kernel_match": True}
        cases.append(case)
    edge_start = time.perf_counter()
    edges = edge_checks()
    edge_seconds = time.perf_counter()-edge_start
    sources = [Path(__file__).resolve(), HERE/"labelled_class_audit.py", HERE/"n4_vertical.py", HERE/"n3_certificate.py", HERE/"n4_boundary.py"]
    result = {
        "status": "PASS",
        "scope": "All finite-exchangeable n=4 hit laws; exact basis checks at all eight original ecological parameter combinations. Two reference polygon boundaries qualified; no interior E-class optimisation performed.",
        "n": N,
        "cases": cases,
        "summary": {"ecological_parameter_combinations": len(cases), "actions": 2,
                    "fixed_count_kernels_checked": len(cases)*2*(N+1),
                    "moment_kernel_entries_compared": entries,
                    "extra_group_row_entries_compared_for_lumpability": lumpability_checks},
        "mathematical_basis": "The five count laws form the vertices of the exchangeable hit-law simplex. Averaging is linear, and tau_j=E[C(K,j)]/C(4,j). Therefore exact equality at every count vertex implies equality throughout the simplex, including every fixed-p,q polygon.",
        "edge_qualification": edges,
        "edge_qualification_s": round(edge_seconds,3),
        "elapsed_s": round(time.perf_counter()-start, 3),
        "runtime": {"python": platform.python_version(), "gmpy2": gmpy2.version()},
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    }
    OUTPUT.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], **result["summary"], "elapsed_s": result["elapsed_s"]}))


if __name__ == "__main__":
    main()
