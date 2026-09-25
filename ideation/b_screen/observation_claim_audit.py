"""Analytical check of the monitoring-impossibility inference in Fable #555.

This evaluates one explicit rule for two already specified pairs of known
laws. It neither fits an observation model nor solves learning over all E.
Run: python3 ideation/b_screen/observation_claim_audit.py
"""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUTS = (HERE / "mechanism_audit.json", HERE / "review_additions_audit.json")


def enclosure(value):
    scale = 10**15
    lower = value.numerator * scale // value.denominator
    upper = -((-value.numerator * scale) // value.denominator)
    return {"lower": str(F(lower, scale)), "upper": str(F(upper, scale)),
            "display": float(value)}


def rule_record(name, pi_a, pi_b, delta_a, delta_b, intervals):
    """Choose C if any observed K=2, otherwise R; exact expected regrets."""
    assert delta_a < 0 < delta_b
    for pi in (pi_a, pi_b):
        assert len(pi) == 5 and min(pi) >= 0 and sum(pi) == 1
        assert sum(k * v for k, v in enumerate(pi)) / 4 == F(1, 10)
        assert sum(k * (k - 1) * v for k, v in enumerate(pi)) / 12 == F(3, 100)
    error_a = (1 - pi_a[2])**intervals
    error_b = 1 - (1 - pi_b[2])**intervals
    regret_a, regret_b = -delta_a * error_a, delta_b * error_b
    worst = max(regret_a, regret_b)
    assert worst < F(1, 100)
    return {"name": name, "intervals": intervals,
            "pi_A": list(map(str, pi_a)), "pi_B": list(map(str, pi_b)),
            "all_count_categories_positive": min(pi_a + pi_b) > 0,
            "error_A": enclosure(error_a), "error_B": enclosure(error_b),
            "expected_regret_A": enclosure(regret_a),
            "expected_regret_B": enclosure(regret_b),
            "worst_expected_regret": enclosure(worst),
            "exact_worst_regret_below_one_percentage_point": True}


def main():
    mechanism = json.loads(INPUTS[0].read_text())["cases"]
    selected = [r for r in mechanism if F(r["parameters"]["q"]) == F(3, 100)]
    selected.sort(key=lambda r: F(r["parameters"]["tau4"]))
    assert len(selected) == 2
    pi = [list(map(F, r["parameters"]["pi_by_hit_count"])) for r in selected]
    delta = [F(r["direct_persistence"]["Delta_PR_minus_PC"]["exact"])
             for r in selected]
    rows = [rule_record("Original E extremisers", *pi, *delta, 15)]

    mixtures = json.loads(INPUTS[1].read_text())["full_support_q03"]["mixtures"]
    mixed_pi, mixed_delta = [], []
    for key in ("low", "high"):
        row = mixtures[key]
        labelled = {x: F(v) for x, v in row["mixed_law"].items()}
        assert len(labelled) == 16 and min(labelled.values()) > 0
        mixed_pi.append([sum(v for x, v in labelled.items() if x.count("1") == k)
                         for k in range(5)])
        mixed_delta.append(F(row["persistence"]["Delta_PR_minus_PC"]["exact"]))
    rows.append(rule_record("90:10 extremiser/Glo mixtures", *mixed_pi, *mixed_delta, 20))

    # Zero K=4 at N intervals excludes probability .03 at one-sided 95%
    # exactly when (.97)**N < .05. Certify the first crossing without logs.
    n_first = next(n for n in range(1, 1000) if F(97, 100)**n < F(1, 20))
    assert n_first == 99 and F(97, 100)**(n_first - 1) >= F(1, 20)
    result = {
        "scope": "Two known-law pairs; original immediate-recovery kernel only",
        "assumptions": ["fixed true law belongs to the declared two-law set",
                        "iid, perfectly observed training intervals",
                        "training record independent of future ecological draws",
                        "original fully occupied initial state at the decision",
                        "no waiting cost or nuisance-parameter estimation"],
        "rule": "Choose C if any training interval has K=2; otherwise choose R",
        "loss": "Expected 50-step persistence regret over training records, at each true law",
        "method": "Exact rational arithmetic using previously certified action differences",
        "zero_all_hit_95_upper": {str(n): 1 - .05**(1/n) for n in (50, 99, 100)},
        "first_zero_event_upper_below_003_exact": n_first,
        "cases": rows,
        "limitation": "No optimality claim, whole-E learning bound, field sample requirement or delayed-model result",
        "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (*INPUTS, Path(__file__))},
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    for row in rows:
        print(row["name"], row["intervals"],
              "worst expected regret (pp):", 100 * row["worst_expected_regret"]["display"])


if __name__ == "__main__":
    main()
