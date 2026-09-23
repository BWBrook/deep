#!/usr/bin/env python3
"""Independent exact audit of the saved n-scaling rational witnesses.

This script deliberately does not rerun the floating-point search.  It reads
the rational laws saved by n_scaling_probe.json, verifies their probability
and moment constraints exactly, constructs each averaged transition kernel via
n4_vertical.avg_direct, then propagates 50 model steps with gmpy2.mpq arithmetic.

Run from the repository root:
  uv run --python 3.12 --with gmpy2 --with numpy --with scipy \
    python ideation/b_screen/n_scaling_audit.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

from gmpy2 import mpq


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROBE_JSON = HERE / "n_scaling_probe.json"
OUTPUT = HERE / "n_scaling_audit.json"

# Imported only for the audited exact kernel constructor and state indexing.
import sys

sys.path.insert(0, str(HERE))
import n4_vertical as V  # noqa: E402


P = mpq(1, 10)
Q2 = mpq(3, 100)
E0 = mpq(1, 50)
E1 = mpq(9, 10)
C = mpq(1, 10)
T = 50
DECIMAL_PLACES = 15

HASHED_SOURCES = (
    "ideation/b_screen/n_scaling_audit.py",
    "ideation/b_screen/n_scaling_probe.py",
    "ideation/b_screen/n_scaling_probe.json",
    "ideation/b_screen/n_scaling_probe.md",
    "ideation/b_screen/n_scaling_probe.log",
    "ideation/b_screen/n3_certificate.py",
    "ideation/b_screen/n3_certificate_results.json",
    "ideation/b_screen/n4_vertical.py",
    "ideation/b_screen/n4_boundary.py",
    "ideation/b_screen/n4_boundary_results.json",
)


def q(value) -> mpq:
    """Parse a saved integer/fraction/decimal as an exact mpq."""
    return mpq(str(value))


def rat(value: mpq) -> str:
    value = mpq(value)
    return f"{value.numerator}/{value.denominator}"


def fixed(integer: int, places: int = DECIMAL_PLACES) -> str:
    sign = "-" if integer < 0 else ""
    digits = str(abs(integer)).rjust(places + 1, "0")
    return f"{sign}{digits[:-places]}.{digits[-places:]}"


def decimal_bracket(value: mpq, places: int = DECIMAL_PLACES) -> dict[str, str]:
    """Outward decimal enclosure: floor <= exact value <= ceiling."""
    value = mpq(value)
    scale = 10**places
    numerator = int(value.numerator) * scale
    denominator = int(value.denominator)
    lower = numerator // denominator
    upper = -((-numerator) // denominator)
    assert mpq(lower, scale) <= value <= mpq(upper, scale)
    return {"places": places, "floor": fixed(lower, places), "ceiling": fixed(upper, places)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def matrix_mpq(matrix):
    return [[mpq(x.numerator, x.denominator) for x in row] for row in matrix]


def check_stochastic(kernel) -> None:
    zero, one = mpq(0), mpq(1)
    assert all(x >= zero for row in kernel for x in row)
    assert all(sum(row, zero) == one for row in kernel)
    assert kernel[0][0] == one and all(x == zero for x in kernel[0][1:])


def persistence_direct(action: str, n: int, atoms, weights) -> mpq:
    """Exact direct average plus exact propagation, summing live-state mass."""
    atoms_f = [V.F(int(x.numerator), int(x.denominator)) for x in atoms]
    weights_f = [V.F(int(x.numerator), int(x.denominator)) for x in weights]
    kernel = matrix_mpq(V.avg_direct(action, n, V.F(1, 50), V.F(9, 10), V.F(1, 10), atoms_f, weights_f))
    check_stochastic(kernel)
    state_count = len(kernel)
    vector = [mpq(0)] * state_count
    vector[V.start_state(action, n)] = mpq(1)
    columns = list(zip(*kernel))
    for _ in range(T):
        vector = [sum((vector[i] * column[i] for i in range(state_count)), mpq(0)) for column in columns]
    assert sum(vector, mpq(0)) == 1
    # Direct live-mass summation avoids subtractive cancellation near persistence 1.
    return sum(vector[1:], mpq(0))


def feasible(atoms, weights) -> dict:
    checks = {
        "same_length": len(atoms) == len(weights),
        "atom_count_at_least_two": len(atoms) >= 2,
        "atoms_sorted": atoms == sorted(atoms),
        "atoms_in_unit_interval": all(0 <= x <= 1 for x in atoms),
        "weights_nonnegative": all(x >= 0 for x in weights),
        "weights_sum_one": sum(weights, mpq(0)) == 1,
        "mean_equals_p": sum((w * x for x, w in zip(atoms, weights)), mpq(0)) == P,
        "second_moment_equals_q": sum((w * x * x for x, w in zip(atoms, weights)), mpq(0)) == Q2,
    }
    checks["all"] = all(checks.values())
    return checks


def close_to_saved(exact: mpq, saved: float, tolerance: float = 5e-12) -> dict:
    discrepancy = abs(float(exact) - float(saved))
    return {"saved_float": saved, "absolute_discrepancy": discrepancy, "within_5e-12": discrepancy <= tolerance}


def evaluate_law(n: int, tag: str, saved: dict) -> dict:
    atoms = [q(x) for x in saved["atoms"]]
    weights = [q(x) for x in saved["weights"]]
    checks = feasible(atoms, weights)
    assert checks["all"], (n, tag, checks)
    pr = persistence_direct("R", n, atoms, weights)
    pc = persistence_direct("C", n, atoms, weights)
    delta = pr - pc
    return {
        "tag": tag,
        "atoms": [rat(x) for x in atoms],
        "weights": [rat(x) for x in weights],
        "feasibility": checks,
        "exact": {"PR": rat(pr), "PC": rat(pc), "Delta": rat(delta)},
        "decimal_enclosures": {
            "PR": decimal_bracket(pr),
            "PC": decimal_bracket(pc),
            "Delta": decimal_bracket(delta),
        },
        "agreement_with_probe_exact_reevaluation_floats": {
            "PR": close_to_saved(pr, saved["PR"]),
            "PC": close_to_saved(pc, saved["PC"]),
            "Delta": close_to_saved(delta, saved["Delta_exact_float"]),
        },
    }


def main() -> dict:
    probe = json.loads(PROBE_JSON.read_text())
    n3_certificate = json.loads((HERE / "n3_certificate_results.json").read_text())
    n3_case = next(
        case for case in n3_certificate["cases"]
        if case["id"] == "n3_q3/100_e01/50_e19/10_c1/10"
    )
    n3_delta_certificate = n3_case["metrics"]["delta_PR_minus_PC"]
    assert n3_delta_certificate["certified"]
    n4_certificate = json.loads((HERE / "n4_boundary_results.json").read_text())
    n4_case = next(
        case for case in n4_certificate["cases"]
        if (case["q"], case["e0"], case["e1"], case["c"])
        == ("3/100", "1/50", "9/10", "1/10")
    )
    prior_exact_extrema = {
        3: (
            q(n3_delta_certificate["certified_bounds"]["min"]),
            q(n3_delta_certificate["certified_bounds"]["max"]),
        ),
        4: (
            q(n4_case["curves"]["lower"]["Delta"]["at_u0_exact"]),
            q(n4_case["curves"]["upper"]["Delta"]["at_u1_exact"]),
        ),
    }
    rows = []
    for source_row in probe["rows"]:
        n = int(source_row["n"])
        minimum = evaluate_law(n, "min", source_row["exact_reevaluation"]["min"])
        maximum = evaluate_law(n, "max", source_row["exact_reevaluation"]["max"])
        L = q(minimum["exact"]["Delta"])
        U = q(maximum["exact"]["Delta"])
        assert L < 0 < U
        regret_r = -L
        regret_c = U
        pair_regret = min(regret_r, regret_c)
        lottery_r_probability = U / (U - L)
        pair_lottery_regret = (-L * U) / (U - L)
        floor_bound = decimal_bracket(pair_regret)["floor"]
        exact_not_above_probe = pair_regret <= q(source_row["regret_lower_bound"])
        # This inequality is expected because the probe headline is derived from
        # a float optimizer; the exact rational witness is the certifying object.
        calibration = None
        if n in (3, 4):
            prior_L, prior_U = prior_exact_extrema[n]
            assert (L, U) == (prior_L, prior_U)
            summary_L, summary_U = map(q, source_row["certified_class_range"])
            calibration = {
                "prior_certificate_scope": "whole compatible class",
                "endpoint_witnesses_are_prior_certified_class_extrema": True,
                "exact_match_to_prior_certificate_artifact": True,
                "exact_class_min": rat(L),
                "exact_class_max": rat(U),
                "exact_deterministic_minimax_regret": rat(pair_regret),
                "probe_rounded_certificate_summary": [str(summary_L), str(summary_U)],
                "summary_is_not_an_outward_exact_enclosure": not (summary_L <= L and summary_U >= U),
            }
        rows.append(
            {
                "n": n,
                "witnesses": {"min": minimum, "max": maximum},
                "exact_pair_decision": {
                    "regret_if_choose_R": rat(regret_r),
                    "regret_if_choose_C": rat(regret_c),
                    "deterministic_minimax_regret": rat(pair_regret),
                    "lower_bound_on_whole_class_regret": True,
                    "only_a_lower_bound_not_a_class_value": n >= 5,
                    "equals_prior_certified_whole_class_regret": n in (3, 4),
                    "safe_downward_decimal_floor": floor_bound,
                    "decimal_enclosure": decimal_bracket(pair_regret),
                    "once_only_lottery": {
                        "probability_choose_R": rat(lottery_r_probability),
                        "probability_choose_C": rat(1 - lottery_r_probability),
                        "pair_minimax_regret": rat(pair_lottery_regret),
                        "lower_bound_on_whole_class_lottery_regret": True,
                        "only_a_lower_bound_not_a_class_value": n >= 5,
                        "equals_prior_certified_whole_class_lottery_regret": n in (3, 4),
                        "safe_downward_decimal_floor": decimal_bracket(pair_lottery_regret)["floor"],
                        "decimal_enclosure": decimal_bracket(pair_lottery_regret),
                    },
                },
                "probe_float_headline": {
                    "Delta_min_found": source_row["Delta_min_found"],
                    "Delta_max_found": source_row["Delta_max_found"],
                    "regret_lower_bound": source_row["regret_lower_bound"],
                    "exact_pair_regret_not_above_float_headline": exact_not_above_probe,
                    "consequence": (
                        "The float headline is not itself a certified bound. The exact rational pair and its "
                        "downward decimal floor are the auditable bound."
                    ),
                },
                "n3_n4_calibration": calibration,
            }
        )

    exact_regrets = [q(row["exact_pair_decision"]["deterministic_minimax_regret"]) for row in rows]
    attained_increases_through_n6 = all(exact_regrets[i] < exact_regrets[i + 1] for i in range(3))
    declines_after_n6 = exact_regrets[3] > exact_regrets[4] > exact_regrets[5]
    assert attained_increases_through_n6 and declines_after_n6
    assert all(row["witnesses"][tag]["feasibility"]["all"] for row in rows for tag in ("min", "max"))
    assert all(
        value["within_5e-12"]
        for row in rows
        for tag in ("min", "max")
        for value in row["witnesses"][tag]["agreement_with_probe_exact_reevaluation_floats"].values()
    )

    source_hashes = {name: sha256(ROOT / name) for name in HASHED_SOURCES}
    git_head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout.strip()
    return {
        "audit": {
            "target_commit": git_head,
            "expected_target_commit": "beae4360156717e5e311eeffe33c39efa938554b",
            "target_commit_matches": git_head == "beae4360156717e5e311eeffe33c39efa938554b",
            "method": (
                "Read only the saved rational laws; verify support, nonnegative weights, total mass, and moments "
                "m1=p and m2=q exactly; form kernels with n4_vertical.avg_direct; convert to gmpy2.mpq; propagate "
                "50 model steps; sum nonabsorbing mass directly. No optimization or random search."
            ),
            "parameters": {"p": rat(P), "q": rat(Q2), "e0": rat(E0), "e1": rat(E1), "c": rat(C), "T": T},
            "source_sha256": source_hashes,
        },
        "rows": rows,
        "cross_n_findings": {
            "six_exact_rational_witness_pairs_feasible": True,
            "attained_pair_regret_strictly_increases_n3_through_n6": attained_increases_through_n6,
            "attained_pair_regret_strictly_declines_n6_through_n8": declines_after_n6,
            "whole_class_monotonicity_proved": False,
            "reason": (
                "For n>=5 these are two feasible attained laws, not certified class extrema. Increasing attained "
                "lower bounds through n=6 cannot establish monotonicity of the unknown whole-class regret."
            ),
            "saturation_interpretation": (
                "The saved witnesses show high absolute persistence and a narrowing negative witness after n=6 "
                "at this one parameter setting. They do not identify saturation as the unique cause, locate a "
                "general saturation threshold, or imply that the positive side or whole compatible class narrows."
            ),
        },
        "search_pipeline_audit": {
            "float_feasibility_tolerance": "weights3/weights4 accept weights down to -1e-12",
            "rationalisation_rule": (
                "rationalise drops every searched component with weight <=1e-7, rounds retained atoms, and solves "
                "new exact weights; its result can therefore be a nearby different law rather than an exact "
                "certificate for the floating candidate."
            ),
            "selector_mismatch_risk": (
                "main chooses L and U from each float search winner plus Glo and Ghi, but exact_reevaluation "
                "serializes only rationalisations of the float search winners. It does not explicitly exact-check "
                "which candidate supplied each headline L or U."
            ),
            "concrete_consequence_here": (
                "All twelve saved rationalisations are feasible and reproduce the displayed exact-reevaluation "
                "floats. Their exact pair regrets are slightly below the float regret headlines, so those headlines "
                "must not be used as rigorous lower bounds. The exact fractions or recorded downward floors are safe."
            ),
        },
        "verification": {
            "all_12_laws_feasible": True,
            "all_24_kernels_exactly_stochastic_with_absorbing_zero_state": True,
            "all_exact_results_match_saved_exact_reevaluation_floats_within_5e-12": True,
            "n3_n4_exact_values_match_prior_certified_endpoint_results": True,
            "n5_n8_values_are_certified_attained_bounds_only": True,
        },
    }


if __name__ == "__main__":
    result = main()
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {OUTPUT}")
    for row in result["rows"]:
        pair = row["exact_pair_decision"]
        print(
            f"n={row['n']}: exact pair regret={float(q(pair['deterministic_minimax_regret'])):.15f}; "
            f"safe floor={pair['safe_downward_decimal_floor']}"
        )
