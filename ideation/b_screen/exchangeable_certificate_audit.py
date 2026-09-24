"""Bounded independent audit of the 16 finite-exchangeable certificates.

All 56 polygon edges are rebuilt at their natural degree T=50 using an
affine-kernel recurrence.  This is independent of the pilot's redundant
degree-100 edge representation.  Two full bivariate derivative tensors are
recomputed: the identified R and C cases with the weakest certified margins.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
import hashlib
import itertools
import json
from math import comb, isclose
from pathlib import Path
import time

from gmpy2 import mpq as Q

import n3_certificate as N3
import n4_boundary as Bd
import n4_vertical as V


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
INPUT = HERE / "exchangeable_certificate.json"
OUTPUT = HERE / "exchangeable_certificate_audit.json"
S_INPUT = HERE / "n4_boundary_results.json"
KERNEL_AUDIT = HERE / "exchangeable_kernel_audit.json"
SOURCE_PATHS = (
    HERE / "exchangeable_certificate_audit.py",
    HERE / "exchangeable_certificate.py",
    INPUT,
    HERE / "exchangeable_certificate.md",
    HERE / "exchangeable_certificate.log",
    HERE / "exchangeable_kernel_audit.py",
    KERNEL_AUDIT,
    HERE / "n4_vertical.py",
    HERE / "n4_boundary.py",
    HERE / "n3_certificate.py",
    S_INPUT,
)
EXPECTED_VERTICES = {
    F(3, 100): ((F(0), F(0)), (F(3, 200), F(0)), (F(3, 100), F(3, 100))),
    F(3, 50): (
        (F(2, 75), F(0)),
        (F(3, 100), F(0)),
        (F(3, 50), F(3, 50)),
        (F(1, 25), F(1, 25)),
    ),
}
REPRESENTATIVE_KEYS = (
    ("3/100", "1/50", "9/10", "1/50"),  # weakest identified R margin
    ("3/50", "1/10", "9/10", "1/50"),   # weakest identified C margin
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def constraints(p: F, q: F):
    return (
        (1 - 4 * p + 6 * q, F(-4), F(1)),
        (4 * p - 12 * q, F(12), F(-4)),
        (6 * q, F(-12), F(6)),
        (F(0), F(4), F(-4)),
        (F(0), F(0), F(1)),
    )


def cross(origin, left, right):
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


def polygon(p: F, q: F):
    """Enumerate feasible line intersections and form an exact CCW hull."""
    cons = constraints(p, q)
    points = set()
    for first, second in itertools.combinations(cons, 2):
        a1, b1, c1 = first
        a2, b2, c2 = second
        determinant = b1 * c2 - b2 * c1
        if determinant == 0:
            continue
        x = (-a1 * c2 + a2 * c1) / determinant
        y = (-b1 * a2 + b2 * a1) / determinant
        if all(a + b * x + c * y >= 0 for a, b, c in cons):
            points.add((x, y))
    ordered = sorted(points)
    lower = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    hull = tuple(lower[:-1] + upper[:-1])
    assert all(cross(hull[i - 1], hull[i], hull[(i + 1) % len(hull)]) > 0
               for i in range(len(hull)))
    assert all(all(a + b * x + c * y >= 0 for a, b, c in cons)
               for x, y in hull)
    return hull


def chains(vertices):
    right = max(range(len(vertices)), key=lambda index: vertices[index][0])
    lower = tuple((vertices[i], vertices[i + 1]) for i in range(right))
    upper = tuple(
        (vertices[i], vertices[(i + 1) % len(vertices)])
        for i in range(right, len(vertices))
    )
    return lower, upper


def exact_kernel(parts, point):
    return [[Q(value) for value in row]
            for row in V.kernel_at(*parts, point[0], point[1])]


def edge_delta_degree_t(parts, start, finish):
    coefficients = {}
    for action in ("R", "C"):
        first = exact_kernel(parts[action], start)
        second = exact_kernel(parts[action], finish)
        coefficients[action] = N3.persistence_coefficients(
            first, second, V.start_state(action, 4)
        )
    return [r - c for r, c in zip(coefficients["R"], coefficients["C"])]


def direct_delta(parts, point):
    values = {}
    for action in ("R", "C"):
        values[action] = V.persist(
            V.kernel_at(*parts[action], point[0], point[1]),
            V.start_state(action, 4),
        )
    return Q(values["R"] - values["C"])


def edge_key(chain, start, finish):
    return chain, tuple(map(str, start)), tuple(map(str, finish))


def saved_edges(case):
    return {
        edge_key(edge["chain"], tuple(map(F, edge["from"])),
                 tuple(map(F, edge["to"]))): edge
        for edge in case["edges"].values()
    }


def classify(lower, upper):
    if upper < 0:
        return "identified: C"
    if lower > 0:
        return "identified: R"
    if lower < 0 < upper:
        return "unidentified"
    return "unresolved"


def deterministic_minimax(lower, upper):
    regret_r = max(Q(0), -lower)
    regret_c = max(Q(0), upper)
    choice = "R" if regret_r < regret_c else "C" if regret_c < regret_r else "tie"
    return choice, regret_r, regret_c


def natural_edge_audit(case, scalar_case):
    q, e0, e1, c = (F(case[name]) for name in ("q", "e0", "e1", "c"))
    vertices = polygon(V.P, q)
    assert vertices == EXPECTED_VERTICES[q]
    assert tuple(tuple(map(F, point)) for point in case["polygon_vertices"]) == vertices
    lower, upper = chains(vertices)
    parts = {
        action: V.affine_parts(action, 4, e0, e1, c, V.P, q)
        for action in ("R", "C")
    }
    saved = saved_edges(case)
    expected_keys = {
        edge_key(name, first, second)
        for name, chain in (("lower", lower), ("upper", upper))
        for first, second in chain
    }
    assert set(saved) == expected_keys

    edges = []
    extrema = {"lower": [], "upper": []}
    for name, chain in (("lower", lower), ("upper", upper)):
        for first, second in chain:
            coefficients = edge_delta_degree_t(parts, first, second)
            assert len(coefficients) == V.T + 1
            derivative = [
                V.T * (coefficients[i + 1] - coefficients[i])
                for i in range(V.T)
            ]
            monotonicity = (
                "increasing" if min(derivative) > 0
                else "decreasing" if max(derivative) < 0
                else "not certified monotone"
            )
            assert monotonicity != "not certified monotone"
            minimum = coefficients[0] if monotonicity == "increasing" else coefficients[-1]
            maximum = coefficients[-1] if monotonicity == "increasing" else coefficients[0]
            extrema[name].append((minimum, maximum))
            saved_edge = saved[edge_key(name, first, second)]
            assert saved_edge["degree"] == 2 * V.T
            assert saved_edge["monotone"] == monotonicity
            assert Q(saved_edge["min_attained"]) == minimum
            assert Q(saved_edge["min_outer"]) == minimum
            assert Q(saved_edge["max_attained"]) == maximum
            assert Q(saved_edge["max_outer"]) == maximum
            assert saved_edge["tolerance_met"]
            for u in (F(1, 4), F(3, 4)):
                point = (
                    first[0] + u * (second[0] - first[0]),
                    first[1] + u * (second[1] - first[1]),
                )
                assert Bd.evaluate(coefficients, u) == direct_delta(parts, point)
            edges.append({
                "chain": name,
                "from": list(map(str, first)),
                "to": list(map(str, second)),
                "natural_degree": V.T,
                "saved_degree": saved_edge["degree"],
                "monotonicity": monotonicity,
                "minimum": str(minimum),
                "maximum": str(maximum),
                "minimum_derivative_coefficient": str(min(derivative)),
                "maximum_derivative_coefficient": str(max(derivative)),
                "coefficient_sha256": hashlib.sha256(
                    "\n".join(map(str, coefficients)).encode()
                ).hexdigest(),
                "interior_direct_checks": ["1/4", "3/4"],
            })

    lower_value = min(value[0] for value in extrema["lower"])
    upper_value = max(value[1] for value in extrema["upper"])
    E = case["E_class"]
    assert Q(E["a_E"]["outer"]) == Q(E["a_E"]["attained"]) == lower_value
    assert Q(E["b_E"]["outer"]) == Q(E["b_E"]["attained"]) == upper_value
    classification = classify(lower_value, upper_value)
    assert classification == E["classification"]

    s_min = Q(scalar_case["whole_class"]["Delta"]["min_certificate"]["outer"])
    s_max = Q(scalar_case["whole_class"]["Delta"]["max_certificate"]["outer"])
    assert lower_value <= s_min <= s_max <= upper_value
    assert case["S_class"]["classification"] == (
        "C" if s_max < 0 else "R" if s_min > 0 else "unidentified"
    )

    regret = min(-lower_value, upper_value) if classification == "unidentified" else Q(0)
    lottery = (-lower_value * upper_value / (upper_value - lower_value)
               if classification == "unidentified" else Q(0))
    assert Q(E["regret_attained"]) == regret
    assert isclose(E["regret_float"], float(regret), rel_tol=0, abs_tol=1e-15)
    assert isclose(E["lottery_regret_float"], float(lottery), rel_tol=0, abs_tol=1e-15)
    e_choice, e_regret_r, e_regret_c = deterministic_minimax(lower_value, upper_value)
    s_choice, s_regret_r, s_regret_c = deterministic_minimax(s_min, s_max)

    vertices_direct = {
        f"({point[0]},{point[1]})": direct_delta(parts, point)
        for point in vertices
    }
    for source_key, value in vertices_direct.items():
        assert Q(case["vertex_Delta"][source_key]["exact"]) == value

    return {
        "key": [case[name] for name in ("q", "e0", "e1", "c")],
        "vertices": [list(map(str, point)) for point in vertices],
        "edges": edges,
        "E_range": {"minimum": str(lower_value), "maximum": str(upper_value)},
        "classification": classification,
        "all_extrema_are_vertices": True,
        "S_range": {"minimum": str(s_min), "maximum": str(s_max)},
        "S_range_contained": True,
        "regret_exact": str(regret),
        "lottery_regret_exact": str(lottery),
        "saved_lottery_is_float_only": True,
        "E_deterministic_minimax": {
            "choice": e_choice,
            "worst_regret_if_R": str(e_regret_r),
            "worst_regret_if_C": str(e_regret_c),
        },
        "S_deterministic_minimax": {
            "choice": s_choice,
            "worst_regret_if_R": str(s_regret_r),
            "worst_regret_if_C": str(s_regret_c),
        },
        "deterministic_minimax_choice_unchanged_S_to_E": e_choice == s_choice,
    }


def recompute_root(case):
    q, e0, e1, c = (F(case[name]) for name in ("q", "e0", "e1", "c"))
    vertices = polygon(V.P, q)
    x_low, x_high = min(x for x, _ in vertices), max(x for x, _ in vertices)
    y_low, y_high = min(y for _, y in vertices), max(y for _, y in vertices)
    parts = {
        action: V.affine_parts(action, 4, e0, e1, c, V.P, q)
        for action in ("R", "C")
    }
    started = time.perf_counter()
    tensors = {
        action: V.bernstein_J(
            *parts[action], V.start_state(action, 4),
            x_low, x_high, y_low, y_high,
        )
        for action in ("R", "C")
    }
    delta = [[r - c for r, c in zip(row_r, row_c)]
             for row_r, row_c in zip(tensors["R"], tensors["C"])]
    derivative = V.dw_coeffs(delta)
    minimum = min(min(row) for row in derivative)
    maximum = max(max(row) for row in derivative)
    assert minimum > 0
    saved = case["vertical_on_E"]["rootbox"]
    assert saved["uniform_sign"] == "pos"
    assert Q(saved["min_coeff"]) == minimum
    assert Q(saved["max_coeff"]) == maximum
    assert case["vertical_on_E"]["leaves"] == {
        "pos": 1, "neg": 0, "pruned": 0, "unresolved": 0
    }
    return {
        "key": [case[name] for name in ("q", "e0", "e1", "c")],
        "selection": (
            "weakest identified R margin" if case["E_class"]["classification"] == "identified: R"
            else "weakest identified C margin"
        ),
        "minimum_derivative_coefficient": str(minimum),
        "maximum_derivative_coefficient": str(maximum),
        "uniformly_positive_on_full_bounding_box": True,
        "saved_rootbox_matches": True,
        "elapsed_s": round(time.perf_counter() - started, 3),
    }


def main():
    started = time.perf_counter()
    source = json.loads(INPUT.read_text())
    scalar = json.loads(S_INPUT.read_text())
    kernel = json.loads(KERNEL_AUDIT.read_text())
    assert source["n"] == 4 and source["T"] == V.T == N3.T == Bd.T == 50
    assert source["p"] == str(V.P) == "1/10"
    assert kernel["status"] == "PASS"
    assert kernel["summary"] == {
        "ecological_parameter_combinations": 8,
        "actions": 2,
        "fixed_count_kernels_checked": 80,
        "moment_kernel_entries_compared": 3560,
        "extra_group_row_entries_compared_for_lumpability": 4760,
    }
    for name, expected in source["metadata"]["source_sha256"].items():
        assert sha256(HERE / name) == expected

    expected_cases = {
        tuple(map(str, values))
        for values in itertools.product(V.QS, V.E0S, V.E1S, V.CS)
    }
    by_key = {
        tuple(case[name] for name in ("q", "e0", "e1", "c")): case
        for case in source["cases"]
    }
    scalar_by_key = {
        tuple(case[name] for name in ("q", "e0", "e1", "c")): case
        for case in scalar["cases"]
    }
    assert set(by_key) == set(scalar_by_key) == expected_cases

    edge_results = [natural_edge_audit(by_key[key], scalar_by_key[key])
                    for key in sorted(expected_cases)]
    classifications = Counter(row["classification"] for row in edge_results)
    assert classifications == Counter({
        "identified: C": 12, "identified: R": 2, "unidentified": 2
    })
    assert sum(len(row["edges"]) for row in edge_results) == 56
    assert all(row["all_extrema_are_vertices"] and row["S_range_contained"]
               for row in edge_results)
    assert all(row["deterministic_minimax_choice_unchanged_S_to_E"]
               for row in edge_results)

    margins = {
        key: (-Q(by_key[key]["E_class"]["b_E"]["outer"])
              if by_key[key]["E_class"]["classification"] == "identified: C"
              else Q(by_key[key]["E_class"]["a_E"]["outer"])
              if by_key[key]["E_class"]["classification"] == "identified: R"
              else None)
        for key in expected_cases
    }
    weakest_r = min((margin, key) for key, margin in margins.items()
                    if by_key[key]["E_class"]["classification"] == "identified: R")[1]
    weakest_c = min((margin, key) for key, margin in margins.items()
                    if by_key[key]["E_class"]["classification"] == "identified: C")[1]
    assert (weakest_r, weakest_c) == REPRESENTATIVE_KEYS
    root_recomputations = [recompute_root(by_key[key]) for key in REPRESENTATIVE_KEYS]

    vertical_receipts = []
    for key in sorted(expected_cases):
        case = by_key[key]
        proof = case["vertical_on_E"]
        root = proof["rootbox"]
        assert proof["status"] == "certified pos"
        assert proof["leaves"] == {"pos": 1, "neg": 0, "pruned": 0, "unresolved": 0}
        assert root["uniform_sign"] == "pos" and Q(root["min_coeff"]) > 0
        vertical_receipts.append({
            "key": list(key),
            "minimum_root_coefficient": root["min_coeff"],
            "maximum_root_coefficient": root["max_coeff"],
            "proof_scope": "full bounding box; polygon pruning not invoked",
            "recomputed": key in REPRESENTATIVE_KEYS,
        })

    claims = {
        "all_16_cases_present": len(by_key) == 16,
        "all_56_edges_naturally_rebuilt": sum(len(row["edges"]) for row in edge_results) == 56,
        "all_edge_extrema_match_saved_exactly": True,
        "all_E_extrema_are_attained_vertices": True,
        "all_S_ranges_contained_in_E_ranges": True,
        "all_deterministic_minimax_choices_unchanged_S_to_E": all(
            row["deterministic_minimax_choice_unchanged_S_to_E"]
            for row in edge_results
        ),
        "classification_counts_match": classifications == Counter({
            "identified: C": 12, "identified: R": 2, "unidentified": 2
        }),
        "two_representative_rootboxes_match": len(root_recomputations) == 2,
        "kernel_extension_receipt_passes": kernel["status"] == "PASS",
    }
    assert all(claims.values()), claims
    output = {
        "status": "PASS",
        "scope": (
            "All 56 edge polynomials rebuilt at natural degree 50; all saved root "
            "receipts inspected; full derivative tensors independently rebuilt only "
            "for the weakest identified R and C margins."
        ),
        "source_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in SOURCE_PATHS
        },
        "claims": claims,
        "polygon_method": {
            "enumeration": "all feasible pairwise affine-constraint intersections",
            "ordering": "exact monotone convex hull with strictly positive CCW cross-products",
            "pruning_check": (
                "Source box pruning is conservative: a box is excluded only when "
                "the exact affine maximum of some pi_k over all four corners is negative. "
                "All saved derivative proofs terminate positive at the root, so pruning "
                "is not invoked in any of the 16 certificates."
            ),
        },
        "classification_counts": dict(classifications),
        "case_audits": edge_results,
        "vertical_receipts": vertical_receipts,
        "root_tensor_recomputations": root_recomputations,
        "rounding_scope": (
            "E extrema and deterministic regrets are retained here as exact rationals. "
            "The source JSON stores lottery regret only as a binary float; this audit "
            "recomputes and stores both lottery regrets exactly."
        ),
        "protocol_deviation": (
            "The contract required the two qualification cases to succeed before the "
            "remaining 14. The source sorts those cases first but submits all 16 at once "
            "to Pool(8). This weakens the prospective compute gate but does not alter "
            "the exact numerical certificates."
        ),
        "recomputation_limits": (
            "The established exchangeable-kernel audit is reused rather than repeated. "
            "The edge recurrence and polygon hull are independent of the certificate's "
            "degree-100 edge route. Bivariate tensors for the other 14 cases are inspected "
            "through exact root-box receipts and source bindings, not regenerated."
        ),
        "elapsed_s": round(time.perf_counter() - started, 3),
    }
    OUTPUT.write_text(json.dumps(output, indent=1) + "\n")
    print(json.dumps({
        "status": output["status"],
        "cases": len(edge_results),
        "edges": sum(len(row["edges"]) for row in edge_results),
        "classification_counts": output["classification_counts"],
        "root_recomputations": len(root_recomputations),
        "elapsed_s": output["elapsed_s"],
    }))


if __name__ == "__main__":
    main()
