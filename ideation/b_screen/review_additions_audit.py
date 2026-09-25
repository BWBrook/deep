#!/usr/bin/env python3
"""Targeted exact audit of two new claims in the independent DEEP review.

The method was fixed before the reviewer archive arrived. This script does not
inspect, import, or execute that implementation. It reuses the repository's exact affine moment kernels,
tensor Bernstein recurrence, polygon geometry, one-dimensional Bernstein
enclosure, and independently derived labelled survivor-subset kernel.

Run from the repository root:
  uv run --python 3.12 --with gmpy2 --with numpy --with scipy \
    python ideation/b_screen/review_additions_audit.py
"""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from fractions import Fraction as F
from math import comb
from pathlib import Path

import gmpy2
from gmpy2 import mpq as Q


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT_JSON = HERE / "review_additions_audit.json"
OUTPUT_MD = HERE / "review_additions_audit.md"
LABELLED_JSON = HERE / "labelled_class_audit.json"
REVIEW = ROOT / "ideation/DEEP_independent_review_6ec79d8.md"

sys.path.insert(0, str(HERE))
import exchangeable_certificate as exchangeable  # noqa: E402
import labelled_class_audit as labelled  # noqa: E402
import n4_boundary as boundary  # noqa: E402
import n4_vertical as vertical  # noqa: E402


N = 4
T = 50
P = F(1, 10)
E0 = F(1, 50)
E1 = F(9, 10)
C = F(1, 10)
QS = (F(3, 100), F(3, 50))
TOL = F(1, 10**9)
SOURCE_PATHS = (
    "synthesis/model_specification.md",
    "ideation/DEEP_independent_review_6ec79d8.md",
    "ideation/b_screen/review_additions_audit.py",
    "ideation/b_screen/exchangeable_certificate.py",
    "ideation/b_screen/n3_certificate.py",
    "ideation/b_screen/n4_vertical.py",
    "ideation/b_screen/n4_boundary.py",
    "ideation/b_screen/labelled_class_audit.py",
    "ideation/b_screen/labelled_class_audit.json",
)

def q(value) -> Q:
    if isinstance(value, F):
        return Q(value.numerator, value.denominator)
    return Q(value)


assert (labelled.N, labelled.T) == (N, T)
assert (labelled.P, labelled.E0, labelled.E1, labelled.C) == tuple(map(q, (P, E0, E1, C)))
assert (vertical.T, vertical.P) == (T, P)


def rat(value) -> str:
    value = Q(value)
    return f"{value.numerator}/{value.denominator}"


def exact_value(value) -> dict[str, object]:
    value = Q(value)
    return {
        "exact": rat(value),
        "float": float(value),
        "decimal_enclosure_15_places": labelled.decimal_bracket(value),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tensor_evaluate(coeffs, u: F, w: F) -> Q:
    """Evaluate a square tensor-product Bernstein polynomial exactly."""
    degree_u = len(coeffs) - 1
    degree_w = len(coeffs[0]) - 1
    uq, wq = q(u), q(w)
    bu = [Q(comb(degree_u, i)) * uq**i * (1 - uq) ** (degree_u - i) for i in range(degree_u + 1)]
    bw = [Q(comb(degree_w, j)) * wq**j * (1 - wq) ** (degree_w - j) for j in range(degree_w + 1)]
    return sum((coeffs[i][j] * bu[i] * bw[j] for i in range(degree_u + 1) for j in range(degree_w + 1)), Q(0))


def edge_objective(parts, action: str, p0: tuple[F, F], p1: tuple[F, F]) -> list[Q]:
    """Degree-elevated Bernstein coefficients for persistence along one E edge."""
    def kernel(point):
        matrix = vertical.kernel_at(*parts[action], point[0], point[1])
        return [[q(value) for value in row] for row in matrix]

    k0, k2 = kernel(p0), kernel(p1)
    k1 = [[(a + b) / 2 for a, b in zip(row0, row2)] for row0, row2 in zip(k0, k2)]
    states = len(k0)
    vector = [[Q(0)] * states]
    vector[0][vertical.start_state(action, N)] = Q(1)
    for _ in range(T):
        vector = boundary.bern_mul_deg2(vector, k0, k1, k2)
    return [sum(row, Q(0)) - row[0] for row in vector]


def direct_persistence(parts, action: str, x: F, y: F) -> Q:
    matrix = vertical.kernel_at(*parts[action], x, y)
    return q(vertical.persist(matrix, vertical.start_state(action, N)))


def parse_saved_law(saved: dict[str, str]) -> tuple[Q, ...]:
    law = [Q(0)] * len(labelled.PATTERNS)
    for pattern, value in saved.items():
        law[labelled.PATTERN_INDEX[pattern]] = Q(value)
    return tuple(law)


def law_moments(law: tuple[Q, ...]) -> tuple[Q, Q, Q, Q, Q]:
    total = sum(law, Q(0))
    marginal = sum((weight * pattern[0] for weight, pattern in zip(law, labelled.PATTERNS)), Q(0))
    pair = sum((weight * pattern[0] * pattern[1] for weight, pattern in zip(law, labelled.PATTERNS)), Q(0))
    triple = sum((weight * pattern[0] * pattern[1] * pattern[2] for weight, pattern in zip(law, labelled.PATTERNS)), Q(0))
    fourth = law[-1]
    return total, marginal, pair, triple, fourth


DETERMINISTIC = {
    action: tuple(labelled.deterministic_action_kernel(action, pattern) for pattern in labelled.PATTERNS)
    for action in ("R0", "C")
}


def labelled_persistence(law: tuple[Q, ...], action: str) -> Q:
    zero = Q(0)
    kernel = [[zero for _ in range(labelled.S)] for _ in range(labelled.S)]
    for weight, deterministic in zip(law, DETERMINISTIC[action]):
        if not weight:
            continue
        for old in range(labelled.S):
            for dest in range(labelled.S):
                kernel[old][dest] += weight * deterministic[old][dest]
    labelled.check_kernel(kernel)
    return labelled.persistence(kernel)


def check_law(law: tuple[Q, ...], target_q: Q) -> dict[str, object]:
    assert all(weight >= 0 for weight in law)
    assert sum(law, Q(0)) == 1
    marginals = [
        sum((weight * pattern[i] for weight, pattern in zip(law, labelled.PATTERNS)), Q(0))
        for i in range(N)
    ]
    pairs = [
        sum((weight * pattern[i] * pattern[j] for weight, pattern in zip(law, labelled.PATTERNS)), Q(0))
        for i in range(N) for j in range(i + 1, N)
    ]
    assert all(value == q(P) for value in marginals)
    assert all(value == target_q for value in pairs)
    positive = [weight for weight in law if weight > 0]
    _, _, _, m3, m4 = law_moments(law)
    return {
        "support_size": len(positive),
        "all_16_patterns_positive": len(positive) == 16,
        "minimum_pattern_probability": exact_value(min(positive)),
        "marginals": [rat(value) for value in marginals],
        "pair_probabilities": [rat(value) for value in pairs],
        "specific_triple_probability": rat(m3),
        "all_four_probability": rat(m4),
    }


def vertical_record(tensor, constraints, xlo, xhi, ylo, yhi) -> dict[str, object]:
    derivative = vertical.dw_coeffs(tensor)
    leaves, unresolved, root = exchangeable.certify_dy_polygon(
        derivative, constraints, xlo, xhi, ylo, yhi
    )
    certified_negative = leaves["neg"] > 0 and leaves["pos"] == 0 and leaves["unresolved"] == 0
    flat = [value for row in derivative for value in row]
    return {
        "status": "CERTIFIED_NEGATIVE" if certified_negative else "UNRESOLVED",
        "meaning": "persistence strictly decreases with tau4 over E" if certified_negative else "no global vertical reduction established",
        "leaves": leaves,
        "unresolved_boxes": unresolved[:10],
        "root_box_uniform_sign": root["uniform_sign"],
        "root_coefficient_range_float": [float(min(flat)), float(max(flat))],
        "coefficient_sha256": hashlib.sha256("\n".join(map(str, flat)).encode()).hexdigest(),
    }


def edge_record(parts, action: str, p0: tuple[F, F], p1: tuple[F, F]) -> dict[str, object]:
    coeffs = edge_objective(parts, action, p0, p1)
    validations = []
    for u in (F(1, 4), F(3, 4)):
        x = p0[0] + u * (p1[0] - p0[0])
        y = p0[1] + u * (p1[1] - p0[1])
        bernstein = boundary.evaluate(coeffs, u)
        direct = direct_persistence(parts, action, x, y)
        assert bernstein == direct
        validations.append({"u": str(u), "x": str(x), "y": str(y), "exact_match": True})
    enclosure = boundary.enclose(coeffs, tol=TOL)
    minimum = enclosure["min"]
    u = F(int(minimum["witness_u"].numerator), int(minimum["witness_u"].denominator))
    x = p0[0] + u * (p1[0] - p0[0])
    y = p0[1] + u * (p1[1] - p0[1])
    assert boundary.evaluate(coeffs, u) == minimum["attained"] == direct_persistence(parts, action, x, y)
    return {
        "from": [str(p0[0]), str(p0[1])],
        "to": [str(p1[0]), str(p1[1])],
        "degree": len(coeffs) - 1,
        "minimum": {
            "outer": rat(minimum["outer"]),
            "attained": rat(minimum["attained"]),
            "gap": rat(minimum["gap"]),
            "tolerance": rat(minimum["tolerance"]),
            "tolerance_met": minimum["tolerance_met"],
            "depth_limit_reached": minimum["depth_limit_reached"],
            "witness": {"u": str(u), "tau3": str(x), "tau4": str(y)},
        },
        "direct_validations": validations,
        "coefficient_sha256": hashlib.sha256("\n".join(map(str, coeffs)).encode()).hexdigest(),
    }


def audit_exchangeable_case(target_q: F, saved_case: dict[str, object]) -> dict[str, object]:
    started = time.time()
    parts = {
        action: vertical.affine_parts(action, N, E0, E1, C, P, target_q)
        for action in ("R", "C")
    }
    constraints = exchangeable.pis(P, target_q)
    vertices = exchangeable.polygon(P, target_q)
    _, upper = exchangeable.chains(vertices)
    xlo, xhi = min(x for x, _ in vertices), max(x for x, _ in vertices)
    ylo, yhi = min(y for _, y in vertices), max(y for _, y in vertices)

    tensors = {
        action: vertical.bernstein_J(
            *parts[action], vertical.start_state(action, N), xlo, xhi, ylo, yhi
        )
        for action in ("R", "C")
    }

    # Two strictly interior convex combinations of all vertices.
    weight_sets = (
        [F(1, len(vertices))] * len(vertices),
        [F(i + 1, sum(range(1, len(vertices) + 1))) for i in range(len(vertices))],
    )
    tensor_validations = []
    for weights in weight_sets:
        x = sum((weight * point[0] for weight, point in zip(weights, vertices)), F(0))
        y = sum((weight * point[1] for weight, point in zip(weights, vertices)), F(0))
        assert all(a + b * x + c * y > 0 for a, b, c in constraints)
        u = (x - xlo) / (xhi - xlo)
        w = (y - ylo) / (yhi - ylo)
        checks = {}
        for action in ("R", "C"):
            tensor_value = tensor_evaluate(tensors[action], u, w)
            direct_value = direct_persistence(parts, action, x, y)
            assert tensor_value == direct_value
            checks[action] = {"exact_match": True, "persistence": rat(direct_value)}
        tensor_validations.append({"tau3": str(x), "tau4": str(y), "actions": checks})

    action_records = {}
    global_minima = {}
    for action in ("R", "C"):
        vertical_check = vertical_record(tensors[action], constraints, xlo, xhi, ylo, yhi)
        edge_checks = [edge_record(parts, action, p0, p1) for p0, p1 in upper]
        if vertical_check["status"] != "CERTIFIED_NEGATIVE":
            global_minima[action] = {
                "status": "UNRESOLVED",
                "reason": "vertical derivative was not certified negative; upper-edge minima are attained checks only",
            }
        else:
            outer = min(Q(edge["minimum"]["outer"]) for edge in edge_checks)
            attained = min(Q(edge["minimum"]["attained"]) for edge in edge_checks)
            witness_edge = min(edge_checks, key=lambda edge: Q(edge["minimum"]["attained"]))
            tolerance_met = all(edge["minimum"]["tolerance_met"] for edge in edge_checks)
            global_minima[action] = {
                "status": "CERTIFIED" if tolerance_met else "UNRESOLVED",
                "outer": exact_value(outer),
                "attained": exact_value(attained),
                "global_gap": exact_value(attained - outer),
                "witness": witness_edge["minimum"]["witness"],
                "all_upper_edge_tolerances_met": tolerance_met,
            }
        action_records[action] = {"vertical_certificate": vertical_check, "upper_edge_certificates": edge_checks}

    # Bind each certified minimum to the named high-Delta exchangeable law through
    # the separately derived labelled 16-state kernel.
    high_law = parse_saved_law(saved_case["witnesses"]["high"]["law"])
    law_check = check_law(high_law, q(target_q))
    _, _, _, m3, m4 = law_moments(high_law)
    named = {
        "name": f"q={target_q} high-Delta exchangeable witness",
        "law_constraints": law_check,
        "moments": {"tau3": rat(m3), "tau4": rat(m4)},
        "persistence": {},
    }
    for action, labelled_action in (("R", "R0"), ("C", "C")):
        value = labelled_persistence(high_law, labelled_action)
        direct = direct_persistence(parts, action, F(int(m3.numerator), int(m3.denominator)), F(int(m4.numerator), int(m4.denominator)))
        assert value == direct
        assert global_minima[action]["status"] == "CERTIFIED"
        assert value == Q(global_minima[action]["attained"]["exact"])
        named["persistence"][action] = {**exact_value(value), "exact_match_to_aggregated_chain_and_global_attained_minimum": True}

    expected = {
        F(3, 100): {"R": "0.3857762490", "C": "0.3214370903"},
        F(3, 50): {"R": "0.2273638798", "C": "0.1149017905"},
    }[target_q]
    local_10dp = {
        action: f"{float(Q(global_minima[action]['attained']['exact'])):.10f}"
        for action in ("R", "C")
    }
    assert local_10dp == expected

    return {
        "q": str(target_q),
        "parameters": {"n": N, "T": T, "p": str(P), "e0": str(E0), "e1": str(E1), "c": str(C)},
        "E_polygon_vertices": [[str(x), str(y)] for x, y in vertices],
        "tensor_direct_validations": tensor_validations,
        "actions": action_records,
        "whole_E_minima": global_minima,
        "named_witness_labelled_validation": named,
        "reviewer_display_10dp": expected,
        "local_display_10dp": local_10dp,
        "maximin_action": "R" if Q(global_minima["R"]["outer"]["exact"]) > Q(global_minima["C"]["attained"]["exact"]) else "UNRESOLVED",
        "seconds": round(time.time() - started, 3),
    }


def audit_full_support(labelled_saved: dict[str, object]) -> dict[str, object]:
    target_q = Q(3, 100)
    case = next(row for row in labelled_saved["cases"] if row["q"] == "3/100")
    glo = labelled.iid_mixture_pattern_law((Q(0), Q(3, 10)), (Q(2, 3), Q(1, 3)))
    glo_constraints = check_law(glo, target_q)
    results = {}
    expected = {"low": "-0.172970", "high": "+0.050918"}
    for name in ("low", "high"):
        extremiser = parse_saved_law(case["witnesses"][name]["law"])
        extremiser_constraints = check_law(extremiser, target_q)
        mixed = tuple(Q(9, 10) * a + Q(1, 10) * b for a, b in zip(extremiser, glo))
        mixed_constraints = check_law(mixed, target_q)
        assert mixed_constraints["all_16_patterns_positive"]
        pr = labelled_persistence(mixed, "R0")
        pc = labelled_persistence(mixed, "C")
        delta = pr - pc
        observed = f"{float(delta):+.6f}"
        assert observed == expected[name]
        results[name] = {
            "construction": "9/10 named q=.03 E extremiser + 1/10 Glo",
            "extremiser_constraints": extremiser_constraints,
            "mixed_law": {"".join(map(str, pattern)): rat(weight) for pattern, weight in zip(labelled.PATTERNS, mixed)},
            "mixed_constraints": mixed_constraints,
            "persistence": {"PR": exact_value(pr), "PC": exact_value(pc), "Delta_PR_minus_PC": exact_value(delta)},
            "reviewer_display_6dp": expected[name],
            "local_display_6dp": observed,
            "ranking": "C" if delta < 0 else "R",
        }
    assert results["low"]["ranking"] == "C" and results["high"]["ranking"] == "R"
    return {
        "q": "3/100",
        "Glo": {
            "theta_atoms": ["0", "3/10"],
            "weights": ["2/3", "1/3"],
            "constraints": glo_constraints,
        },
        "mixtures": results,
        "scope": "two attained feasible full-support laws; no class-extremum, realism, or empirical-frequency claim",
    }


def render_markdown(result: dict[str, object]) -> str:
    cases = {row["q"]: row for row in result["whole_E_maximin"]}
    q03, q06 = cases["3/100"], cases["3/50"]
    full = result["full_support_q03"]["mixtures"]
    assert all(
        Q(action["global_gap"]["exact"]) == 0
        for case in cases.values()
        for action in case["whole_E_minima"].values()
    )
    return f"""# Targeted audit of two additions in the independent review

## Decisive findings

Both new numerical claims reproduce using the repository's established exact kernels and Bernstein machinery. This is a new local audit, not an inspection or execution of the reviewer archive, which arrived after this method was fixed.

| q | certified inf PR | certified inf PC | maximin action |
|---:|---:|---:|:---:|
| .03 | {q03['local_display_10dp']['R']} | {q03['local_display_10dp']['C']} | {q03['maximin_action']} |
| .06 | {q06['local_display_10dp']['R']} | {q06['local_display_10dp']['C']} | {q06['maximin_action']} |

For each action and q, an exact tensor-Bernstein certificate shows that persistence decreases strictly with `tau4` throughout the finite-exchangeable polygon E. The global minimum therefore lies on the polygon's upper edge chain. Exact univariate Bernstein branch-and-bound encloses every relevant edge minimum to tolerance `1e-9`; the global attained values equal the separately propagated labelled high-Delta witness exactly. The lower bound and attained value coincide at those endpoint witnesses in these cases.

At q=.03, mixing 90% of each named E extremiser with 10% Glo makes all 16 labelled hit patterns strictly positive while preserving every marginal (`p=.1`) and pair probability (`q=.03`). Exact labelled propagation gives:

| source extremiser | PR-PC | ranking | full labelled support |
|:---|---:|:---:|:---:|
| low-Delta | {full['low']['local_display_6dp']} | {full['low']['ranking']} | yes |
| high-Delta | {full['high']['local_display_6dp']} | {full['high']['ranking']} | yes |

Thus the opposite rankings are not caused solely by zero-probability labelled patterns. These are attained constructions, not new optima and not evidence that either law is ecologically realistic.

## Evidence and assumptions

The estimand is 50-step regional persistence from full occupancy under each complete action policy, with hit vectors drawn independently at each step from one fixed law Q. The maximin prediction target is `argmax_a inf_Q P_a(Q,50)` over the four-patch finite-exchangeable class with exact `p,q`. Pairwise hit probabilities are treated as known; there is no sampling or detection model, temporal dependence, cost calibration, or inference to a field system.

The audit imports `n4_vertical.py` for exact affine moment kernels and tensor Bernstein propagation, `exchangeable_certificate.py` for E geometry and polygon pruning, `n4_boundary.py` for exact edge polynomials and enclosures, and `labelled_class_audit.py` for the separately derived 16-state survivor-subset kernel. Source hashes are in the JSON receipt. These shared biological and mathematical dependencies mean this is a strong internal cross-check, not fully independent software replication.

The maximin/minimax difference at q=.03 is sensitivity to the ambiguity decision criterion applied to the same persistence model. It supplies no additional empirical ecological support for the actions, hit laws, parameter values, or costs.

## Verification

Each action tensor was checked exactly against direct propagation at two strictly interior rational points of each E polygon. Every relevant edge polynomial was checked exactly at `u=1/4` and `u=3/4`. The named extremal witnesses were also propagated through the labelled kernel and matched both the aggregated chain and the global attained minima exactly. The full-support laws passed exact normalization, four marginal, and six pair constraints before propagation.

Reproduce with:

```sh
{result['verification']['command']}
```

Runtime was {result['verification']['seconds_wall']:.3f} seconds on {result['verification']['platform']} with Python {result['verification']['python']} and gmpy2 {result['verification']['gmpy2']}.

## Residual risks and recommendation

Accept the two review additions at their stated mathematical scope. Present maximin beside minimax regret as criterion sensitivity, and present the full-support mixtures only as a narrow rebuttal to a vertex/zero-support objection. No further independent review is needed for these arithmetic checks. The remaining uncertainty concerns ecological interpretation and decision relevance: this audit supplies no empirical calibration of the hit laws, actions, parameter values, or costs.
"""


def main() -> dict[str, object]:
    started = time.time()
    labelled_saved = json.loads(LABELLED_JSON.read_text())
    cases_by_q = {row["q"]: row for row in labelled_saved["cases"]}
    whole_e = [
        audit_exchangeable_case(target_q, cases_by_q[str(target_q)])
        for target_q in QS
    ]
    full_support = audit_full_support(labelled_saved)
    elapsed = time.time() - started
    result = {
        "audit_status": "PASS",
        "scientific_target": {
            "estimand": "50-step regional persistence from full occupancy under complete action R or C",
            "maximin_prediction_target": "argmax_a inf_Q P_a(Q,50) over the n=4 finite-exchangeable class E with exact p and q",
            "full_support_target": "ranking at two specified q=.03 mixtures, not an extremum",
            "data_generating_process": "one fixed environmental-hit law Q; hit vectors iid across steps; conditional patch deaths followed by survivor-driven recolonisation; extinction absorbing",
            "observation_and_sampling": "p and q are supplied exactly; no finite-sample observation, detection, censoring, or missingness model",
            "scope_of_inference": "n=4, p=.1, e0=.02, e1=.9, c=.1, T=50, q in {.03,.06}; prescribed actions only",
        },
        "whole_E_maximin": whole_e,
        "full_support_q03": full_support,
        "interpretation": {
            "criterion_sensitivity": "maximin persistence and minimax regret apply different ambiguity decision rules to the same persistence outcome",
            "empirical_support": "none added by this audit",
            "full_support_meaning": "zero labelled probabilities are unnecessary for opposite rankings at these two attained q=.03 constructions",
        },
        "limitations": [
            "The reviewer archive arrived after this audit method was fixed; its implementation was not inspected, imported, or executed here.",
            "The proof reuses repository kernels and Bernstein machinery; labelled propagation changes state representation and transition construction but shares the model specification.",
            "No temporal dependence, nonexchangeability, parameter uncertainty, cost equivalence, observation error, or empirical plausibility was evaluated.",
            "The full-support laws are feasible attained witnesses, not extrema or frequency statements.",
        ],
        "verification": {
            "command": "uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/review_additions_audit.py",
            "seconds_wall": round(elapsed, 3),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "gmpy2": gmpy2.version(),
            "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True, capture_output=True).stdout.strip(),
            "source_sha256": {path: sha256(ROOT / path) for path in SOURCE_PATHS},
        },
    }
    OUTPUT_JSON.write_text(json.dumps(result, indent=2) + "\n")
    OUTPUT_MD.write_text(render_markdown(result))
    return result


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "audit_status": result["audit_status"],
        "seconds_wall": result["verification"]["seconds_wall"],
        "maximin": [
            {"q": row["q"], "min_PR": row["local_display_10dp"]["R"], "min_PC": row["local_display_10dp"]["C"], "action": row["maximin_action"]}
            for row in result["whole_E_maximin"]
        ],
        "full_support_deltas": {
            name: row["local_display_6dp"]
            for name, row in result["full_support_q03"]["mixtures"].items()
        },
        "outputs": [str(OUTPUT_JSON.relative_to(ROOT)), str(OUTPUT_MD.relative_to(ROOT))],
    }, indent=2))
