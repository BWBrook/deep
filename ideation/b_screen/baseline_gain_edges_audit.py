#!/usr/bin/env python3
"""Exact audit of Fable's finite E/S boundary gain samples.

This reconstructs only the points already sampled in baseline_gain_edges.py.
Every saved point is evaluated with the independent labelled survivor-subset
kernel and compared exactly with the reduced kernels used by the source.  It
does not add points or claim a continuous-boundary certificate.

Run from the repository root:
  uv run --python 3.12 --with gmpy2 --with numpy --with scipy \
    python ideation/b_screen/baseline_gain_edges_audit.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from fractions import Fraction
from math import comb
from pathlib import Path

from gmpy2 import mpq


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "baseline_gain_edges_audit.json"
SOURCE_JSON = HERE / "baseline_gain_edges.json"
N, S, T = 4, 16, 50
P, E0, E1, C = mpq(1, 10), mpq(1, 50), mpq(9, 10), mpq(1, 10)
ZERO, ONE = mpq(0), mpq(1)

sys.path.insert(0, str(HERE))
import exchangeable_certificate as exchangeable  # noqa: E402
import labelled_class_audit as independent  # noqa: E402
import n4_vertical as reduced  # noqa: E402

PATTERNS = independent.PATTERNS
HASHED_INPUTS = (
    "ideation/b_screen/baseline_gain_edges_audit.py",
    "ideation/b_screen/baseline_gain_edges.py",
    "ideation/b_screen/baseline_gain_edges.json",
    "ideation/b_screen/baseline_gain_edges.md",
    "ideation/b_screen/baseline_gain_edges.log",
    "ideation/b_screen/mechanism_audit.py",
    "ideation/b_screen/mechanism_audit.json",
    "ideation/b_screen/labelled_class_audit.py",
    "ideation/b_screen/n4_vertical.py",
    "ideation/b_screen/exchangeable_certificate.py",
)
EXPECTED_VERTICES = {
    mpq(3, 100): ((mpq(0), mpq(0)), (mpq(3, 200), mpq(0)), (mpq(3, 100), mpq(3, 100))),
    mpq(3, 50): ((mpq(2, 75), mpq(0)), (mpq(3, 100), mpq(0)), (mpq(3, 50), mpq(3, 50)), (mpq(1, 25), mpq(1, 25))),
}


def q(value) -> mpq:
    return mpq(str(value))


def rat(value) -> str:
    value = mpq(value)
    return f"{value.numerator}/{value.denominator}"


def val(value) -> dict[str, object]:
    return {"exact": rat(value), "float": float(value)}


def frac(value) -> Fraction:
    value = mpq(value)
    return Fraction(int(value.numerator), int(value.denominator))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def vm(vector, matrix):
    return [sum((vector[i] * matrix[i][j] for i in range(S)), ZERO) for j in range(S)]


def persist(matrix):
    vector = [ZERO] * S
    vector[-1] = ONE
    for _ in range(T):
        vector = vm(vector, matrix)
    return sum(vector[1:], ZERO)


def weighted(law, kernels):
    return [
        [sum((law[x] * kernels[x][i][j] for x in PATTERNS), ZERO) for j in range(S)]
        for i in range(S)
    ]


def pattern_law(q2, tau3, tau4):
    pi = (
        ONE - 4 * P + 6 * q2 - 4 * tau3 + tau4,
        4 * P - 12 * q2 + 12 * tau3 - 4 * tau4,
        6 * q2 - 12 * tau3 + 6 * tau4,
        4 * tau3 - 4 * tau4,
        tau4,
    )
    assert all(x >= 0 for x in pi) and sum(pi, ZERO) == ONE
    law = {x: pi[sum(x)] / comb(N, sum(x)) for x in PATTERNS}
    assert sum(law.values(), ZERO) == ONE
    assert all(sum((law[x] * x[i] for x in PATTERNS), ZERO) == P for i in range(N))
    assert all(
        sum((law[x] * x[i] * x[j] for x in PATTERNS), ZERO) == q2
        for i in range(N) for j in range(i + 1, N)
    )
    return pi, law


def action_parameters(action, hit):
    excess = E1 - E0
    if action == "C":
        return tuple(E0 + excess * x for x in hit), 2 * C
    if action == "R":
        return tuple(E0 + excess * x * (mpq(1, 2) if i == 0 else ONE) for i, x in enumerate(hit)), C
    assert action == "0"
    return tuple(E0 + excess * x for x in hit), C


def deterministic_kernels():
    kernels = {action: {} for action in ("0", "R", "C")}
    for action in kernels:
        for hit in PATTERNS:
            death, colonisation = action_parameters(action, hit)
            kernel = independent.independent_kernel(death, colonisation)
            assert all(sum(row, ZERO) == ONE for row in kernel)
            kernels[action][hit] = kernel
    return kernels


def reduced_parts(q2):
    return {
        action: reduced.affine_parts(action, N, frac(E0), frac(E1), frac(C), frac(P), frac(q2))
        for action in ("0", "R", "C")
    }


def evaluate(q2, tau3, tau4, kernels, parts):
    pi, law = pattern_law(q2, tau3, tau4)
    labelled = {}
    reduced_values = {}
    for action in ("0", "R", "C"):
        labelled[action] = persist(weighted(law, kernels[action]))
        kernel = reduced.kernel_at(*parts[action], frac(tau3), frac(tau4))
        reduced_values[action] = q(reduced.persist(kernel, reduced.start_state(action, N)))
        assert labelled[action] == reduced_values[action]
    return {
        "tau3": rat(tau3),
        "tau4": rat(tau4),
        "pi_by_hit_count": [rat(x) for x in pi],
        "P0": val(labelled["0"]),
        "PR_minus_P0": val(labelled["R"] - labelled["0"]),
        "PC_minus_P0": val(labelled["C"] - labelled["0"]),
        "Delta_PR_minus_PC": val(labelled["R"] - labelled["C"]),
        "exact_match_to_source_reduced_kernel": True,
    }


def ranges(records):
    out = {}
    for key in ("P0", "PR_minus_P0", "PC_minus_P0", "Delta_PR_minus_PC"):
        vals = [q(row[key]["exact"]) for row in records]
        lo, hi = min(vals), max(vals)
        out[key] = {
            "min": val(lo),
            "max": val(hi),
            "span": val(hi - lo),
            "min_attainers": [{"tau3": row["tau3"], "tau4": row["tau4"]} for row in records if q(row[key]["exact"]) == lo],
            "max_attainers": [{"tau3": row["tau3"], "tau4": row["tau4"]} for row in records if q(row[key]["exact"]) == hi],
        }
    gain = out["PR_minus_P0"]
    lo, hi = q(gain["min"]["exact"]), q(gain["max"]["exact"])
    constant = (lo + hi) / 2
    out["finite_sample_best_constant_minus_connectivity_fit"] = {
        "form": "Delta_hat = b - (PC-P0)",
        "Chebyshev_best_b": val(constant),
        "maximum_absolute_residual_on_saved_sample": val((hi - lo) / 2),
        "identity": "residual = (PR-P0)-b; this is finite-sample algebra, not a fitted continuous-boundary model",
    }
    return out


def compare_saved_float(exact_row, saved_row, boundary):
    assert abs(float(q(exact_row["tau3"])) - saved_row["tau3"]) < 5e-17
    assert abs(float(q(exact_row["tau4"])) - saved_row["tau4"]) < 5e-17
    pairs = (
        (exact_row["P0"]["float"], saved_row["P0"]),
        (exact_row["PR_minus_P0"]["float"], saved_row["PR-P0"]),
        (exact_row["PC_minus_P0"]["float"], saved_row["PC-P0"]),
    )
    error = max(abs(a - b) for a, b in pairs)
    assert error < 5e-15
    if boundary == "E":
        assert exact_row["edge"] == saved_row["edge"] and exact_row["u"] == saved_row["u"]
    else:
        assert exact_row["curve"] == saved_row["curve"]
    return error


def tie_to_four_extremisers(q2, records, mechanism):
    specs = (
        (mpq(0), mpq(0), "q=3/100 E-min (0,0)"),
        (mpq(3, 100), mpq(3, 100), "q=3/100 E-max (3/100,3/100)"),
    ) if q2 == mpq(3, 100) else (
        (mpq(2, 75), mpq(0), "q=3/50 E-min (2/75,0)"),
        (mpq(3, 50), mpq(3, 50), "q=3/50 E-max (3/50,3/50)"),
    )
    out = []
    mechanism_by_name = {case["name"]: case for case in mechanism["cases"]}
    for tau3, tau4, name in specs:
        row = next(x for x in records if q(x["tau3"]) == tau3 and q(x["tau4"]) == tau4)
        old = mechanism_by_name[name]["direct_persistence"]
        assert q(row["P0"]["exact"]) == q(old["P0"]["exact"])
        assert q(row["PR_minus_P0"]["exact"]) == q(old["PR_minus_P0"]["exact"])
        assert q(row["PC_minus_P0"]["exact"]) == q(old["PC_minus_P0"]["exact"])
        out.append({"name": name, "tau3": rat(tau3), "tau4": rat(tau4), "exact_match": True})
    return out


def audit_q(q2, saved, kernels, mechanism):
    key = rat(q2)
    saved_q = saved[key]
    vertices = [(q(x), q(y)) for x, y in saved_q["vertices"]]
    independently_enumerated = [(q(x), q(y)) for x, y in exchangeable.polygon(frac(P), frac(q2))]
    assert tuple(vertices) == EXPECTED_VERTICES[q2]
    assert vertices == independently_enumerated
    parts = reduced_parts(q2)

    polygon_records = []
    for edge, (start, end) in enumerate(zip(vertices, vertices[1:] + vertices[:1])):
        for j in range(10):
            u = mpq(j, 10)
            tau3 = start[0] + u * (end[0] - start[0])
            tau4 = start[1] + u * (end[1] - start[1])
            row = {"edge": edge, "u": str(frac(u)), **evaluate(q2, tau3, tau4, kernels, parts)}
            polygon_records.append(row)
    assert len(polygon_records) == len(vertices) * 10 == len(saved_q["polygon_boundary"])
    assert len({(x["tau3"], x["tau4"]) for x in polygon_records}) == len(polygon_records)

    lens_records = []
    xlo_f, xhi_f = reduced.xbounds(frac(P), frac(q2))
    xlo, xhi = q(xlo_f), q(xhi_f)
    for j in range(11):
        tau3 = xlo + mpq(j, 10) * (xhi - xlo)
        for curve, function in (("L", reduced.Lfun), ("U", reduced.Ufun)):
            tau4 = q(function(frac(tau3), frac(P), frac(q2)))
            lens_records.append({"curve": curve, "x_index": j, **evaluate(q2, tau3, tau4, kernels, parts)})
    assert len(lens_records) == len(saved_q["lens_boundary"]) == 22
    unique_lens = {(x["tau3"], x["tau4"]) for x in lens_records}
    assert len(unique_lens) == 20
    duplicate_lens = sorted(
        ({"tau3": x, "tau4": y, "multiplicity": sum(r["tau3"] == x and r["tau4"] == y for r in lens_records)}
         for x, y in unique_lens if sum(r["tau3"] == x and r["tau4"] == y for r in lens_records) > 1),
        key=lambda row: (q(row["tau3"]), q(row["tau4"])),
    )
    assert len(duplicate_lens) == 2 and all(x["multiplicity"] == 2 for x in duplicate_lens)

    max_error = ZERO
    for exact_row, saved_row in zip(polygon_records, saved_q["polygon_boundary"]):
        max_error = max(max_error, q(compare_saved_float(exact_row, saved_row, "E")))
    for exact_row, saved_row in zip(lens_records, saved_q["lens_boundary"]):
        max_error = max(max_error, q(compare_saved_float(exact_row, saved_row, "S")))

    return {
        "q": rat(q2),
        "E_polygon_sample": {
            "vertices": [{"tau3": rat(x), "tau4": rat(y)} for x, y in vertices],
            "edges": len(vertices),
            "records": len(polygon_records),
            "unique_points": len(polygon_records),
            "sampling": "u=0,1/10,...,9/10 on each edge: each start vertex once plus nine strict interior points; edge endpoints at u=1 are represented as the next edge's u=0",
            "records_exact": polygon_records,
            "sample_ranges": ranges(polygon_records),
        },
        "S_lens_boundary_sample": {
            "records": len(lens_records),
            "unique_points": len(unique_lens),
            "duplicate_endpoint_records": duplicate_lens,
            "sampling": "11 shared tau3 values on each of L and U; L=U at both endpoints, yielding 22 records but 20 unique points",
            "records_exact": lens_records,
            "sample_ranges": ranges(lens_records),
        },
        "four_extremiser_exact_ties": tie_to_four_extremisers(q2, polygon_records, mechanism),
        "maximum_absolute_error_against_saved_float_fields": float(max_error),
        "all_sampled_laws_feasible_for_E": True,
        "all_labelled_persistences_match_source_reduced_persistences_exactly": True,
    }


def main():
    started = time.perf_counter()
    saved = json.loads(SOURCE_JSON.read_text())
    mechanism = json.loads((HERE / "mechanism_audit.json").read_text())
    assert set(saved) == {"3/100", "3/50"} and mechanism["status"] == "PASS"
    kernels = deterministic_kernels()
    cases = [audit_q(q2, saved, kernels, mechanism) for q2 in (mpq(3, 100), mpq(3, 50))]
    assert sum(x["E_polygon_sample"]["records"] for x in cases) == 70
    assert sum(x["S_lens_boundary_sample"]["records"] for x in cases) == 44
    assert sum(x["S_lens_boundary_sample"]["unique_points"] for x in cases) == 40

    receipt = {
        "status": "PASS",
        "scientific_target": "verify the finite saved E-polygon and S-lens boundary samples of complete-action persistence gains",
        "estimands": ["P_R-P_0", "P_C-P_0", "P_0", "Delta=(P_R-P_0)-(P_C-P_0)"],
        "scope": "only the 70 saved E records and 44 saved S records at q=3/100 and 3/50; 40 unique S points because both endpoints are duplicated",
        "process": "iid exchangeable hazard law at each step; mortality then recolonisation from survivors; all occupied start; T=50; absorbing zero state",
        "actions": "baseline uses ordinary mortality and c; R halves hazard excess mortality at patch 0 and uses c; C uses ordinary mortality and 2c",
        "kernel_reconstruction": "independent exact labelled 16-state closed survivor-subset kernels; exact equality to source reduced kernels at every sampled point",
        "cases": cases,
        "claims_boundary": {
            "supported": "the reported numerical ranges describe the finite saved samples and reproduce exactly before float conversion",
            "unsupported": [
                "the samples do not certify ranges over every point of the continuous E boundary",
                "the samples do not certify the E interior or the whole E class",
                "the samples do not certify monotonicity between sample points",
                "terms such as nearly invariant require an explicitly chosen scientific tolerance; a narrow sampled range is descriptive evidence, not a theorem",
            ],
            "source_storage": "vertices and E interpolation u are retained exactly, but coordinates and all gain/persistence values are stored only as binary floats; S interpolation indices and exact coordinates are not retained in source JSON",
        },
        "recomputation_scope": {
            "new_sample_points_added": 0,
            "continuous_certificate_or_optimisation_run": False,
            "all_saved_records_recomputed": True,
            "exact_rational_sample_ranges_retained": True,
        },
        "source_hashes_sha256": {path: sha256(ROOT / path) for path in HASHED_INPUTS},
        "seconds_wall": time.perf_counter() - started,
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: 70 E records and 44 S records (40 unique) reconstructed exactly in {receipt['seconds_wall']:.2f} s")
    for case in cases:
        er = case["E_polygon_sample"]["sample_ranges"]
        sr = case["S_lens_boundary_sample"]["sample_ranges"]
        print(
            f"  q={case['q']}: E buffering span={er['PR_minus_P0']['span']['float']:.9f}, "
            f"connectivity span={er['PC_minus_P0']['span']['float']:.9f}, "
            f"best-constant residual={er['finite_sample_best_constant_minus_connectivity_fit']['maximum_absolute_residual_on_saved_sample']['float']:.9f}; "
            f"S buffering span={sr['PR_minus_P0']['span']['float']:.9f}"
        )


if __name__ == "__main__":
    main()
