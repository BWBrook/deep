"""Independent, bounded audit of the eight-point scalar failure-region pilot.

The audit deliberately avoids repeating all eight expensive two-dimensional
Bernstein tensors.  It independently reconstructs the endpoint laws with the
direct fixed-theta kernels, reconstructs both exact boundary polynomials, and
checks every saved classification.  With ``--time-tensor`` it additionally
decomposes the tensor/vertical cost at P1 to audit the pilot's trace estimate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import time

import n3_certificate as N3
import n4_boundary as Bd
import n4_vertical as V


HERE = Path(__file__).resolve().parent
POINTS_PATH = HERE / "failure_region_pilot_points.json"
PILOT_PATH = HERE / "failure_region_pilot.json"
OUTPUT_PATH = HERE / "failure_region_audit.json"
DEPENDENCIES = (
    "failure_region_audit.py",
    "failure_region_pilot.py",
    "failure_region_pilot_points.json",
    "failure_region_pilot.json",
    "failure_region_pilot.md",
    "failure_region_pilot.log",
    "n3_certificate.py",
    "n4_vertical.py",
    "n4_boundary.py",
)
EXPECTED_POINTS = (
    ("P1", "1/2", "9/50"),
    ("P2", "1/2", "11/50"),
    ("P3", "7/10", "1/10"),
    ("P4", "7/10", "3/20"),
    ("P5", "9/10", "3/50"),
    ("P6", "17/20", "1/50"),
    ("P7", "1/2", "2/5"),
    ("P8", "19/20", "9/20"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_persistence(kernel, start: int, steps: int) -> F:
    vector = [F(0)] * len(kernel)
    vector[start] = F(1)
    for _ in range(steps):
        vector = N3.matvec(vector, kernel)
    return sum(vector[1:])


def direct_delta(e0: F, e1: F, c: F, law, steps: int) -> F:
    atoms, weights = law
    values = {}
    for action in ("R", "C"):
        kernel = N3.averaged_kernel(action, 4, e0, e1, c, atoms, weights)
        values[action] = exact_persistence(
            kernel, V.start_state(action, 4), steps
        )
    return values["R"] - values["C"]


def curve_delta(parts, x_lower: F, x_upper: F, boundary, p: F, q: F):
    r = Bd.objective_coeffs(
        *parts["R"], V.start_state("R", 4),
        x_lower, x_upper, boundary, p, q,
    )
    c = Bd.objective_coeffs(
        *parts["C"], V.start_state("C", 4),
        x_lower, x_upper, boundary, p, q,
    )
    return [a - b for a, b in zip(r, c)]


def classify(lower: F, upper: F) -> str:
    if upper < 0:
        return "identified: C"
    if lower > 0:
        return "identified: R"
    if lower < 0 < upper:
        return "unidentified"
    return "unresolved"


def audit_point(point, saved, p: F, q: F, e0: F, steps: int,
                laws, x_lower: F, x_upper: F):
    e1, c = F(point["e1"]), F(point["c"])
    direct = {
        name: direct_delta(e0, e1, c, law, steps)
        for name, law in zip(("Glo", "Ghi"), laws)
    }
    endpoint_matches = {
        name: direct[name] == F(saved["witnesses"][name]["exact"])
        for name in direct
    }

    parts = {
        action: V.affine_parts(action, 4, e0, e1, c, p, q)
        for action in ("R", "C")
    }
    curves = {}
    for name, boundary in (("lower", V.Lfun), ("upper", V.Ufun)):
        coefficients = curve_delta(
            parts, x_lower, x_upper, boundary, p, q
        )
        degree = len(coefficients) - 1
        derivatives = [
            degree * (coefficients[index + 1] - coefficients[index])
            for index in range(degree)
        ]
        saved_curve = saved["curves"][name]
        curves[name] = {
            "degree": degree,
            "all_derivative_coefficients_strictly_positive": all(
                value > 0 for value in derivatives
            ),
            "minimum_derivative_coefficient": str(min(derivatives)),
            "minimum_derivative_coefficient_float": float(min(derivatives)),
            "left_endpoint": str(coefficients[0]),
            "right_endpoint": str(coefficients[-1]),
            "saved_attained_endpoints_match": (
                coefficients[0] == F(saved_curve["min_attained"])
                and coefficients[-1] == F(saved_curve["max_attained"])
            ),
            "saved_outer_endpoints_match": (
                coefficients[0] == F(saved_curve["min_outer"])
                and coefficients[-1] == F(saved_curve["max_outer"])
            ),
        }

    classification = classify(direct["Glo"], direct["Ghi"])
    leaves = saved["direct_sign_certificate"]["leaves"]
    identified_by_complete_tensor_sign_proof = (
        leaves["unresolved"] == 0
        and ((leaves["pos"] > 0) != (leaves["neg"] > 0))
    )
    classification_basis = (
        "complete two-dimensional tensor sign proof"
        if classification.startswith("identified")
        else "opposite-sign exact attained endpoint witnesses"
    )
    return {
        "id": point["id"],
        "e1": point["e1"],
        "c": point["c"],
        "direct_endpoint_delta": {
            name: {"exact": str(value), "float": float(value)}
            for name, value in direct.items()
        },
        "direct_endpoint_matches_saved": endpoint_matches,
        "curves": curves,
        "saved_vertical_certificate_complete_positive": (
            saved["vertical"]["status"] == "certified pos"
            and saved["vertical"]["leaves"]["pos"] > 0
            and saved["vertical"]["leaves"]["neg"] == 0
            and saved["vertical"]["leaves"]["unresolved"] == 0
        ),
        "saved_direct_tensor_leaves": leaves,
        "identified_by_complete_tensor_sign_proof": (
            identified_by_complete_tensor_sign_proof
        ),
        "classification": classification,
        "classification_basis": classification_basis,
        "matches_saved_shortcut_classification": (
            classification == saved["shortcut_classification"]
        ),
    }


def time_p1_tensor(p: F, q: F, e0: F, steps: int, point):
    e1, c = F(point["e1"]), F(point["c"])
    x_lower, x_upper = V.xbounds(p, q)
    y_lower, y_upper = V.lens_range(x_lower, x_upper, p, q)
    timings = {}

    start = time.perf_counter()
    parts = {
        action: V.affine_parts(action, 4, e0, e1, c, p, q)
        for action in ("R", "C")
    }
    timings["affine_parts_s"] = time.perf_counter() - start

    tensors = {}
    for action in ("R", "C"):
        start = time.perf_counter()
        tensors[action] = V.bernstein_J(
            *parts[action], V.start_state(action, 4),
            x_lower, x_upper, y_lower, y_upper, T=steps,
        )
        timings[f"tensor_{action}_s"] = time.perf_counter() - start

    start = time.perf_counter()
    delta = [
        [a - b for a, b in zip(row_r, row_c)]
        for row_r, row_c in zip(tensors["R"], tensors["C"])
    ]
    timings["tensor_difference_s"] = time.perf_counter() - start

    start = time.perf_counter()
    direct_leaves, _ = V.certify_dw(
        delta, x_lower, x_upper, y_lower, y_upper, p, q, max_depth=10
    )
    timings["direct_sign_branch_and_bound_s"] = time.perf_counter() - start

    start = time.perf_counter()
    derivative = V.dw_coeffs(delta)
    vertical_leaves, _ = V.certify_dw(
        derivative, x_lower, x_upper, y_lower, y_upper,
        p, q, max_depth=10,
    )
    timings["derivative_and_vertical_certificate_s"] = (
        time.perf_counter() - start
    )
    timings["tensor_required_for_vertical_s"] = (
        timings["tensor_R_s"]
        + timings["tensor_C_s"]
        + timings["tensor_difference_s"]
        + timings["derivative_and_vertical_certificate_s"]
    )
    return {
        "point": point["id"],
        "timings": {name: round(value, 6) for name, value in timings.items()},
        "direct_leaves": direct_leaves,
        "vertical_leaves": vertical_leaves,
    }


def main(time_tensor: bool):
    points_document = json.loads(POINTS_PATH.read_text())
    pilot = json.loads(PILOT_PATH.read_text())
    expected_metadata = {"n": 4, "p": "1/10", "q": "3/100", "e0": "1/50", "T": 50}
    assert all(points_document.get(key) == value for key, value in expected_metadata.items())
    assert all(pilot.get(key) == value for key, value in expected_metadata.items())
    point_spec = tuple(
        (row["id"], row["e1"], row["c"])
        for row in points_document["points"]
    )
    pilot_spec = tuple(
        (row["id"], row["e1"], row["c"])
        for row in pilot["results"]
    )
    assert point_spec == EXPECTED_POINTS
    assert pilot_spec == EXPECTED_POINTS
    assert len({row[0] for row in point_spec}) == 8
    saved_by_id = {row["id"]: row for row in pilot["results"]}
    p, q = F(points_document["p"]), F(points_document["q"])
    e0, steps = F(points_document["e0"]), points_document["T"]
    assert V.P == p and V.T == steps and Bd.P == p and Bd.T == steps
    laws = N3.laws(p, q)
    x_lower, x_upper = V.xbounds(p, q)
    y_lower, y_upper = V.lens_range(x_lower, x_upper, p, q)
    assert pilot["lens"] == [str(x_lower), str(x_upper), str(y_lower), str(y_upper)]

    started = time.perf_counter()
    results = [
        audit_point(
            point, saved_by_id[point["id"]], p, q, e0, steps,
            laws, x_lower, x_upper,
        )
        for point in points_document["points"]
    ]
    classifications = {}
    for row in results:
        classifications[row["classification"]] = (
            classifications.get(row["classification"], 0) + 1
        )

    reported = [row["cost_s"] for row in pilot["results"]]
    mean_total = sum(row["total"] for row in reported) / len(reported)
    mean_direct = sum(row["direct_sign"] for row in reported) / len(reported)
    mean_curves = sum(row["curves"] for row in reported) / len(reported)
    pass_flags = {
        "all_endpoint_values_match": all(
            all(row["direct_endpoint_matches_saved"].values()) for row in results
        ),
        "all_boundary_derivative_coefficients_strictly_positive": all(
            all(curve["all_derivative_coefficients_strictly_positive"]
                for curve in row["curves"].values())
            for row in results
        ),
        "all_saved_boundary_endpoints_match": all(
            all(curve["saved_attained_endpoints_match"]
                and curve["saved_outer_endpoints_match"]
                for curve in row["curves"].values())
            for row in results
        ),
        "all_saved_shortcut_classifications_match": all(
            row["matches_saved_shortcut_classification"] for row in results
        ),
        "classification_counts_match": classifications == {
            "identified: C": 3,
            "unidentified": 4,
            "identified: R": 1,
        },
    }
    assert all(pass_flags.values()), pass_flags

    audit = {
        "status": "PASS",
        "scope": "eight scalar-class points; no E/L transfer and no full tensor rerun",
        "source_sha256": {
            name: sha256(HERE / name) for name in DEPENDENCIES
        },
        "results": results,
        "classification_counts": classifications,
        **pass_flags,
        "resource_audit": {
            "pilot_mean_direct_sign_stage_s": mean_direct,
            "pilot_mean_curves_stage_s": mean_curves,
            "pilot_mean_total_s": mean_total,
            "pilot_projected_900_point_sequential_hours": mean_total * 900 / 3600,
            "reason_claimed_shortcut_is_not_cheap": (
                "The direct_sign timer includes construction of JR, JC and Jd. "
                "Stage C reuses Jd, so omitting the direct sign branch-and-bound "
                "does not omit the expensive tensor construction required by the "
                "current vertical-certificate implementation."
            ),
        },
        "provenance_limits": [
            "The pilot output stores no source hashes or Git revision.",
            "The points file and all evaluated outputs first appear in one commit; "
            "the script writes the points file before its evaluation loop, but Git "
            "history alone cannot establish prospective freezing before results were known.",
            "Exact endpoint reconstruction here uses the shared n3_certificate kernel "
            "implementation; it is independent of the affine moment/tensor route, not "
            "an independent implementation of the ecological transition model.",
            "P1's full Delta tensor and both sign checks were recomputed once for "
            "resource decomposition when --time-tensor was used. P2-P8 tensors were "
            "not repeated; their saved leaf counts remain bound only by the committed "
            "source, JSON and log."
        ],
        "tensor_recomputation": {
            "P1": (
                "full JR, JC, Jd, direct sign and vertical sign recomputed"
                if time_tensor else "not recomputed in this invocation"
            ),
            "P2-P8": "saved tensor leaf summaries inspected, not recomputed",
        },
        "cheap_audit_seconds": round(time.perf_counter() - started, 3),
    }
    if time_tensor:
        audit["resource_audit"]["measured_P1_decomposition"] = time_p1_tensor(
            p, q, e0, steps, points_document["points"][0]
        )
    OUTPUT_PATH.write_text(json.dumps(audit, indent=1) + "\n")
    print(
        "PASS: 8/8 endpoint and boundary audits; "
        "classifications 4 unidentified, 3 C, 1 R; "
        f"cheap stage {audit['cheap_audit_seconds']} s"
    )
    if time_tensor:
        required = audit["resource_audit"]["measured_P1_decomposition"][
            "timings"
        ]["tensor_required_for_vertical_s"]
        print(f"P1 tensor required for vertical certificate: {required:.3f} s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--time-tensor", action="store_true")
    arguments = parser.parse_args()
    main(arguments.time_tensor)
