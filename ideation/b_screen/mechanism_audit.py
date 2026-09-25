#!/usr/bin/env python3
"""Independent exact audit of the four saved mechanism decompositions.

The audit rebuilds every deterministic labelled kernel from the closed
survivor-subset expression in labelled_class_audit.py, rather than reusing the
nested constructor used by mechanism_decomposition.py.  It then reconstructs
the four fixed exchangeable laws, both T=50 telescoping allocations, the true
horizon-specific persistence contrasts, and occupancy-count trajectories.

Run from the repository root:
  uv run --python 3.12 --with gmpy2 --with numpy --with scipy \
    python ideation/b_screen/mechanism_audit.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from math import comb
from pathlib import Path

from gmpy2 import mpq


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "mechanism_audit.json"
sys.path.insert(0, str(HERE))

import labelled_class_audit as independent  # noqa: E402
import mechanism_decomposition as source  # noqa: E402


N, S, T = 4, 16, 50
P, E0, E1, C = mpq(1, 10), mpq(1, 50), mpq(9, 10), mpq(1, 10)
PATTERNS = independent.PATTERNS
ZERO, ONE = mpq(0), mpq(1)
CASES = (
    ("q=3/100 E-min (0,0)", mpq(3, 100), mpq(0), mpq(0), "low", "a_E"),
    ("q=3/100 E-max (3/100,3/100)", mpq(3, 100), mpq(3, 100), mpq(3, 100), "high", "b_E"),
    ("q=3/50 E-min (2/75,0)", mpq(3, 50), mpq(2, 75), mpq(0), "low", "a_E"),
    ("q=3/50 E-max (3/50,3/50)", mpq(3, 50), mpq(3, 50), mpq(3, 50), "high", "b_E"),
)
HASHED_INPUTS = (
    "ideation/b_screen/mechanism_audit.py",
    "ideation/b_screen/mechanism_decomposition.py",
    "ideation/b_screen/mechanism_decomposition.json",
    "ideation/b_screen/mechanism_decomposition.md",
    "ideation/b_screen/mechanism_decomposition.log",
    "ideation/b_screen/labelled_class_audit.py",
    "ideation/b_screen/labelled_class_audit.json",
    "ideation/b_screen/exchangeable_certificate.py",
    "ideation/b_screen/exchangeable_certificate.json",
    "ideation/b_screen/exchangeable_certificate_audit.py",
    "ideation/b_screen/exchangeable_certificate_audit.json",
)


def q(value) -> mpq:
    return mpq(str(value))


def rat(value) -> str:
    value = mpq(value)
    return f"{value.numerator}/{value.denominator}"


def value(value) -> dict[str, object]:
    value = mpq(value)
    return {"exact": rat(value), "float": float(value)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def vm(vector, matrix):
    return [sum((vector[i] * matrix[i][j] for i in range(S)), ZERO) for j in range(S)]


def mv(matrix, vector):
    return [sum((matrix[i][j] * vector[j] for j in range(S)), ZERO) for i in range(S)]


def dot(left, right):
    return sum((a * b for a, b in zip(left, right)), ZERO)


def madd(left, right, sign=1):
    return [[left[i][j] + sign * right[i][j] for j in range(S)] for i in range(S)]


def weighted(law, kernels):
    return [
        [sum((law[x] * kernels[x][i][j] for x in PATTERNS), ZERO) for j in range(S)]
        for i in range(S)
    ]


def law_from_moments(q2, tau3, tau4):
    pi = (
        1 - 4 * P + 6 * q2 - 4 * tau3 + tau4,
        4 * P - 12 * q2 + 12 * tau3 - 4 * tau4,
        6 * q2 - 12 * tau3 + 6 * tau4,
        4 * tau3 - 4 * tau4,
        tau4,
    )
    assert all(x >= 0 for x in pi) and sum(pi, ZERO) == ONE
    law = {pattern: pi[sum(pattern)] / comb(N, sum(pattern)) for pattern in PATTERNS}
    assert sum(law.values(), ZERO) == ONE
    assert all(sum((law[x] * x[i] for x in PATTERNS), ZERO) == P for i in range(N))
    assert all(
        sum((law[x] * x[i] * x[j] for x in PATTERNS), ZERO) == q2
        for i in range(N) for j in range(i + 1, N)
    )
    assert all(
        sum((law[x] * x[i] * x[j] * x[k] for x in PATTERNS), ZERO) == tau3
        for i in range(N) for j in range(i + 1, N) for k in range(j + 1, N)
    )
    assert law[(1, 1, 1, 1)] == tau4
    return pi, law


def action_parameters(action, pattern):
    excess = E1 - E0
    if action == "C":
        return tuple(E0 + excess * x for x in pattern), 2 * C
    if action == "R":
        return tuple(
            E0 + excess * x * (mpq(1, 2) if i == 0 else ONE)
            for i, x in enumerate(pattern)
        ), C
    assert action == "0"
    return tuple(E0 + excess * x for x in pattern), C


def deterministic_kernels():
    ours = {action: {} for action in ("0", "R", "C")}
    for action in ours:
        for pattern in PATTERNS:
            death, colonisation = action_parameters(action, pattern)
            ours[action][pattern] = independent.independent_kernel(death, colonisation)
    source_kernels = source.det_kernels()
    entries = 0
    max_difference = ZERO
    for action in ours:
        for pattern in PATTERNS:
            for i in range(S):
                assert sum(ours[action][pattern][i], ZERO) == ONE
                for j in range(S):
                    difference = abs(ours[action][pattern][i][j] - q(source_kernels[action][pattern][i][j]))
                    max_difference = max(max_difference, difference)
                    entries += 1
    assert entries == 3 * 16 * 256 and max_difference == 0
    return ours, {
        "deterministic_kernels": 48,
        "entries_compared": entries,
        "exact_max_absolute_difference": rat(max_difference),
        "exact_match": True,
        "independent_derivation": "closed sum over exact survivor subsets versus source nested survivor/colonisation enumeration",
    }


def sign(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


def sign_runs(sequence, start_index=0):
    runs = []
    for index, x in enumerate(sequence, start_index):
        sx = sign(x)
        if not runs or runs[-1]["sign"] != sx:
            runs.append({"start": index, "end": index, "sign": sx})
        else:
            runs[-1]["end"] = index
    return runs


def exact_series(sequence):
    return [{"t": i, **value(x)} for i, x in enumerate(sequence)]


def saved_external_delta(labelled, certificates, q2, witness, extremum):
    labelled_case = next(case for case in labelled["cases"] if q(case["q"]) == q2)
    labelled_delta = q(labelled_case["witnesses"][witness]["persistence"]["Delta_PR_minus_PC"]["exact"])
    certificate_case = next(
        case for case in certificates["cases"]
        if q(case["q"]) == q2 and q(case["e0"]) == E0 and q(case["e1"]) == E1 and q(case["c"]) == C
    )
    certificate_delta = q(certificate_case["E_class"][extremum]["attained"])
    return labelled_delta, certificate_delta


def validate_upstream_hash_receipts(labelled, certificates):
    checked = {}
    for declared_path, declared_hash in labelled["verification"]["source_sha256"].items():
        observed = sha256(ROOT / declared_path)
        assert observed == declared_hash
        checked[declared_path] = True
    for basename, declared_hash in certificates["metadata"]["source_sha256"].items():
        path = HERE / basename
        observed = sha256(path)
        assert observed == declared_hash
        checked[str(path.relative_to(ROOT))] = True
    return {"files_checked": len(checked), "all_declared_hashes_match_current_files": True}


def audit_case(name, q2, tau3, tau4, witness, extremum, kernels, saved, labelled, certificates):
    pi, law = law_from_moments(q2, tau3, tau4)
    K = {action: weighted(law, kernels[action]) for action in ("0", "R", "C")}
    mu = [ZERO] * S
    mu[-1] = ONE
    h = [ZERO] + [ONE] * (S - 1)

    trajectories = {action: [mu] for action in ("0", "R", "C")}
    for action in trajectories:
        for _ in range(T):
            trajectories[action].append(vm(trajectories[action][-1], K[action]))
    persistence = {
        action: [sum(v[1:], ZERO) for v in trajectories[action]]
        for action in trajectories
    }
    delta_horizon = [r - c for r, c in zip(persistence["R"], persistence["C"])]
    delta = delta_horizon[T]
    d1_formula = ((E1 - E0) / 2) * (
        E0**3 * P + 3 * E0**2 * (E1 - E0) * q2
        + 3 * E0 * (E1 - E0) ** 2 * tau3 + (E1 - E0) ** 3 * tau4
    )
    assert delta_horizon[1] == d1_formula > 0

    saved_case = saved["results"][name]
    assert [q(x) for x in saved_case["pi"]] == list(pi)
    assert q(saved_case["Delta_1_exact"]) == d1_formula
    assert q(saved_case["Delta_T_exact"]) == delta
    labelled_delta, certificate_delta = saved_external_delta(labelled, certificates, q2, witness, extremum)
    assert delta == labelled_delta == certificate_delta

    # Backward vectors for the two exact telescoping orders.
    backwards = {}
    for action in ("R", "C"):
        backwards[action] = [h]
        for _ in range(T):
            backwards[action].append(mv(K[action], backwards[action][-1]))
    dR0 = madd(K["R"], K["0"], -1)
    dC0 = madd(K["C"], K["0"], -1)
    dR0x, dC0x = {}, {}
    for pattern in PATTERNS:
        dR0x[pattern] = [
            [law[pattern] * (kernels["R"][pattern][i][j] - kernels["0"][pattern][i][j]) for j in range(S)]
            for i in range(S)
        ]
        dC0x[pattern] = [
            [law[pattern] * (kernels["C"][pattern][i][j] - kernels["0"][pattern][i][j]) for j in range(S)]
            for i in range(S)
        ]

    decompositions = {}
    order_specs = (
        ("forward_R_past_C_future", trajectories["R"], backwards["C"]),
        ("reverse_C_past_R_future", trajectories["C"], backwards["R"]),
    )
    source_rows = {row["ordering"].split()[0]: row for row in saved_case["decompositions"]}
    raw_by_order = {}
    for key, left, right in order_specs:
        B, G, Bk_t, Gk_t = [], [], [], []
        for t in range(T):
            a, b = left[t], right[T - 1 - t]
            B.append(dot(vm(a, dR0), b))
            G.append(dot(vm(a, dC0), b))
            bk, gk = [ZERO] * 5, [ZERO] * 5
            for pattern in PATTERNS:
                k = sum(pattern)
                bk[k] += dot(vm(a, dR0x[pattern]), b)
                gk[k] += dot(vm(a, dC0x[pattern]), b)
            assert sum(bk, ZERO) == B[-1] and sum(gk, ZERO) == G[-1]
            Bk_t.append(bk)
            Gk_t.append(gk)
        assert all(x >= 0 for x in B) and all(x >= 0 for x in G)
        assert G[-1] == 0
        Bsum, Gsum = sum(B, ZERO), sum(G, ZERO)
        Bk = [sum((row[k] for row in Bk_t), ZERO) for k in range(5)]
        Gk = [sum((row[k] for row in Gk_t), ZERO) for k in range(5)]
        assert sum(Bk, ZERO) == Bsum and sum(Gk, ZERO) == Gsum
        assert Bsum - Gsum == delta

        source_key = "forward" if key.startswith("forward") else "reverse"
        source_row = source_rows[source_key]
        floats_to_check = (
            [(Bsum, source_row["buffering_total"]), (Gsum, source_row["colonisation_total"])]
            + list(zip(B, source_row["buffering_by_t"]))
            + list(zip(G, source_row["colonisation_by_t"]))
            + list(zip(Bk, source_row["buffering_by_hitcount"]))
            + list(zip(Gk, source_row["colonisation_by_hitcount"]))
        )
        decomposition_float_error = max(abs(float(exact) - observed) for exact, observed in floats_to_check)
        assert decomposition_float_error < 5e-15
        raw_by_order[key] = (Bsum, Gsum, Bk, Gk)
        decompositions[key] = {
            "definition": "fixed-T=50 hybrid telescoping allocation; prefixes are not shorter-horizon Delta",
            "buffering_total": value(Bsum),
            "colonisation_Gamma_total": value(Gsum),
            "net_total": value(Bsum - Gsum),
            "components_nonnegative": True,
            "final_step_Gamma_exactly_zero": True,
            "by_hit_count": [
                {
                    "k": k,
                    "buffering": value(Bk[k]),
                    "buffering_fraction": value(Bk[k] / Bsum),
                    "colonisation_Gamma": value(Gk[k]),
                    "colonisation_Gamma_fraction": value(Gk[k] / Gsum),
                    "net": value(Bk[k] - Gk[k]),
                }
                for k in range(5)
            ],
            "by_time": [
                {
                    "t": t,
                    "buffering": value(B[t]),
                    "colonisation_Gamma": value(G[t]),
                    "net_float": float(B[t] - G[t]),
                }
                for t in range(T)
            ],
            "time_sign_runs": {
                "buffering": sign_runs(B),
                "colonisation_Gamma": sign_runs(G),
                "net": sign_runs([b - g for b, g in zip(B, G)]),
            },
            "max_abs_error_against_saved_float_fields": decomposition_float_error,
        }

    fB, fG, fBk, fGk = raw_by_order["forward_R_past_C_future"]
    rB, rG, rBk, rGk = raw_by_order["reverse_C_past_R_future"]
    order_difference = {
        "buffering_forward_minus_reverse": value(fB - rB),
        "colonisation_forward_minus_reverse": value(fG - rG),
        "net_difference_exactly_zero": fB - fG == rB - rG,
        "by_hit_count_forward_minus_reverse": [
            {"k": k, "buffering": value(fBk[k] - rBk[k]), "colonisation_Gamma": value(fGk[k] - rGk[k]),
             "net": value((fBk[k] - fGk[k]) - (rBk[k] - rGk[k]))}
            for k in range(5)
        ],
    }

    occupancy = []
    max_saved_occupancy_error = 0.0
    for t in range(T + 1):
        row = {"t": t}
        for action in ("R", "C"):
            counts = [
                sum((trajectories[action][t][state] for state in range(S) if state.bit_count() == k), ZERO)
                for k in range(5)
            ]
            assert sum(counts, ZERO) == ONE
            # Full float trajectories support plotting; exact anchors reproduce
            # every time point the source elected to save.
            row[action] = [float(x) for x in counts]
            if t in (0, 1, 5, 10, 25, 49, 50):
                row[f"{action}_exact"] = [rat(x) for x in counts]
            if str(t) in saved_case["occupancy_count_distribution"]:
                observed = saved_case["occupancy_count_distribution"][str(t)][action]
                max_saved_occupancy_error = max(max_saved_occupancy_error, *(abs(float(x) - y) for x, y in zip(counts, observed)))
        occupancy.append(row)
    assert max_saved_occupancy_error < 5e-15

    direct = {
        "P0": value(persistence["0"][T]),
        "PR": value(persistence["R"][T]),
        "PC": value(persistence["C"][T]),
        "PR_minus_P0": value(persistence["R"][T] - persistence["0"][T]),
        "PC_minus_P0": value(persistence["C"][T] - persistence["0"][T]),
        "Delta_PR_minus_PC": value(delta),
        "interpretation": "direct pure-chain contrasts; distinct from the order-dependent hybrid telescoping totals",
    }

    # Conditional on a full starting occupancy and an all-four hit, extinction
    # means all current occupants die.  The event label itself is not mortality.
    unbuffered_survival = ONE - E1**4
    buffered_survival = ONE - E1**3 * (E0 + (E1 - E0) / 2)
    conditional_contrast = buffered_survival - unbuffered_survival

    return {
        "name": name,
        "parameters": {"q": rat(q2), "tau3": rat(tau3), "tau4": rat(tau4), "pi_by_hit_count": [rat(x) for x in pi]},
        "exchangeable_pattern_law_checks": {
            "labelled_pattern_probability": "Q(x)=pi_|x|/choose(4,|x|)",
            "normalises_exactly": True,
            "all_four_marginals_equal_p_exactly": True,
            "all_six_pair_moments_equal_q_exactly": True,
            "all_four_triple_moments_equal_tau3_exactly": True,
            "four_way_moment_equals_tau4_exactly": True,
        },
        "external_exact_value_matches": {
            "mechanism_decomposition_json": True,
            "labelled_class_audit_json": True,
            "exchangeable_certificate_attained_value": True,
            "Delta": value(delta),
        },
        "one_step": {"Delta": value(delta_horizon[1]), "formula_exact_match": True, "R_advantage": True},
        "true_horizon_Delta_PR_minus_PC": exact_series(delta_horizon),
        "true_horizon_sign_runs": sign_runs(delta_horizon),
        "sign_transition_horizons": [
            t for t in range(1, T + 1) if sign(delta_horizon[t]) != sign(delta_horizon[t - 1])
        ],
        "strict_R_advantage_to_C_advantage_crossings": [
            t for t in range(2, T + 1) if delta_horizon[t - 1] > 0 > delta_horizon[t]
        ],
        "direct_persistence": direct,
        "decompositions": decompositions,
        "order_difference": order_difference,
        "occupancy_count_trajectories": occupancy,
        "max_abs_error_against_saved_occupancy_float_fields": max_saved_occupancy_error,
        "all_four_hit_conditional_from_full_occupancy": {
            "unbuffered_or_C_survival": value(unbuffered_survival),
            "buffered_R_survival": value(buffered_survival),
            "conditional_survival_contrast": value(conditional_contrast),
            "event_semantics": "all four patches receive the hazard hit; extinction additionally requires death of every currently occupied patch; survivors may remain because e1=9/10 and the buffered hit death probability is 23/50",
            "not_equal_to_T50_k4_telescope_entry": conditional_contrast != fBk[4],
        },
    }


def main():
    started = time.perf_counter()
    saved = json.loads((HERE / "mechanism_decomposition.json").read_text())
    labelled = json.loads((HERE / "labelled_class_audit.json").read_text())
    certificates = json.loads((HERE / "exchangeable_certificate.json").read_text())
    assert (saved["n"], saved["T"], q(saved["e0"]), q(saved["e1"]), q(saved["c"]), q(saved["p"])) == (N, T, E0, E1, C, P)
    assert set(saved["results"]) == {case[0] for case in CASES}
    upstream_hashes = validate_upstream_hash_receipts(labelled, certificates)

    kernels, kernel_receipt = deterministic_kernels()
    results = [audit_case(*case, kernels, saved, labelled, certificates) for case in CASES]
    assert all(case["external_exact_value_matches"]["mechanism_decomposition_json"] for case in results)
    elapsed = time.perf_counter() - started
    receipt = {
        "status": "PASS",
        "scientific_scope": "exact reconstruction of only the four pre-existing E extremisers; no search, optimisation, parameter grid, or new law",
        "process_and_estimand": {
            "estimand": "Delta_t = Pr(nonextinct at t | buffer patch 0) - Pr(nonextinct at t | doubled colonisation), from all four patches occupied",
            "data_generating_process": "iid exchangeable hazard pattern at each model step; mortality precedes same-step recolonisation from surviving occupied patches; zero occupancy is absorbing",
            "scope_of_inference": "n=4, t=0..50, p=1/10, e0=1/50, e1=9/10, c=1/10 and the four named laws only",
            "actions_and_baseline": "K0 uses ordinary mortality and c; KR halves hazard-attributable excess mortality at labelled patch 0 and uses c; KC uses ordinary mortality and 2c",
        },
        "kernel_validation": kernel_receipt,
        "upstream_hash_receipt_validation": upstream_hashes,
        "cases": results,
        "interpretive_findings": {
            "fixed_T_prefix_warning": "cum_Delta_by_t in the source is a prefix of a T=50 hybrid allocation and is not Delta at a shorter horizon",
            "hit_is_not_death": "a k=4 hazard hit is not an all-four mortality event; extinction depends on occupancy and realised deaths",
            "conditional_vs_weighted_warning": "the conditional one-step survival contrast under an all-four hit is not the weighted T=50 k=4 contribution",
            "source_provenance_gap": "mechanism_decomposition.json and .log contain no source hashes; exact values are externally bound here to two independently saved certificate/audit artefacts",
            "colonisation_quiet_step_fraction": "forward k=0 Gamma shares are 57.2% to 72.3%; reverse shares are 46.0% to 61.3%, so 'almost entirely at k=0' is unsupported and q=.03 E-max reverse has a larger k=1 share",
            "order_dependence": "overall Delta is invariant, but q=.06 E-max k=1 net changes sign from positive forward to negative reverse, so ordering can change a category-level qualitative statement",
            "extinction_timing": "R has the exact one-step advantage at all four laws. The minimisers cross to C at horizons 4 and 11; maximisers retain R's advantage for every t=1..50, with equality only at t=0",
            "last_death_semantics": "extinction occurs when every currently occupied patch dies in the mortality phase; same-step colonisation needs at least one survivor, so zero is absorbing. A k=4 hazard hit need not be extinction",
        },
        "recomputation_scope": {
            "deterministic_kernel_entries_exactly_compared": kernel_receipt["entries_compared"],
            "weighted_kernels_rebuilt": 12,
            "horizon_steps_per_action_and_law": T,
            "both_telescoping_orders_rebuilt_exactly": True,
            "all_time_and_hit_count_components_rebuilt_exactly": True,
            "all_occupancy_count_trajectories_t0_to_t50_rebuilt_exactly": True,
            "external_search_or_certificate_generation_repeated": False,
        },
        "source_hashes_sha256": {path: sha256(ROOT / path) for path in HASHED_INPUTS},
        "seconds_wall": elapsed,
    }
    OUTPUT.write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"PASS: 48 deterministic kernels / {kernel_receipt['entries_compared']} entries; 4 laws; exact T=50 decompositions and t=0..50 trajectories; {elapsed:.2f} s")
    for case in results:
        d = case["external_exact_value_matches"]["Delta"]["float"]
        quiet = case["decompositions"]["forward_R_past_C_future"]["by_hit_count"][0]["colonisation_Gamma_fraction"]["float"]
        runs = case["true_horizon_sign_runs"]
        print(f"  {case['name']}: Delta50={d:+.9f}; forward Gamma k0 share={quiet:.3%}; horizon sign runs={runs}")


if __name__ == "__main__":
    main()
