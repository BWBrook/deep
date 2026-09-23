#!/usr/bin/env python3
"""Independent exact audit of the four highlighted labelled-class witnesses.

This audit deliberately does not rerun the LP/random/Nelder-Mead search.  It
reads the four saved rational witnesses, checks the labelled moment polytope
with independent exact linear algebra, and evaluates them with a separately
derived 16-state kernel.  The kernel sums a closed transition probability over
all survivor subsets; it does not reuse the nested survivor/colonisation loops
in labelled_class_probe.kernel_det.

Run from the repository root:
  uv run --python 3.12 --with gmpy2 --with numpy --with scipy \
    python ideation/b_screen/labelled_class_audit.py
"""

from __future__ import annotations

import hashlib
import itertools
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from gmpy2 import mpq


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROBE_JSON = HERE / "labelled_class_probe.json"
BOUNDARY_JSON = HERE / "n4_boundary_results.json"
OUTPUT = HERE / "labelled_class_audit.json"

sys.path.insert(0, str(HERE))
import labelled_class_probe as probe  # noqa: E402  (audited comparator only)
import n4_vertical as vertical  # noqa: E402  (independent count-chain comparator)


N = 4
S = 1 << N
T = 50
P = mpq(1, 10)
E0 = mpq(1, 50)
E1 = mpq(9, 10)
C = mpq(1, 10)
PATTERNS = tuple(itertools.product((0, 1), repeat=N))
PATTERN_INDEX = {"".join(map(str, x)): i for i, x in enumerate(PATTERNS)}
DECIMAL_PLACES = 15

HASHED_SOURCES = (
    "ideation/b_screen/labelled_class_audit.py",
    "ideation/b_screen/labelled_class_probe.py",
    "ideation/b_screen/labelled_class_probe.json",
    "ideation/b_screen/labelled_class_probe.md",
    "ideation/b_screen/labelled_class_probe.log",
    "ideation/b_screen/n3_certificate.py",
    "ideation/b_screen/n4_vertical.py",
    "ideation/b_screen/n4_boundary.py",
    "ideation/b_screen/n4_boundary_results.json",
)


def q(value) -> mpq:
    """Parse integers, fractions and saved decimals without binary float."""
    return mpq(str(value))


def rat(value) -> str:
    value = mpq(value)
    return f"{value.numerator}/{value.denominator}"


def frac(value) -> Fraction:
    value = mpq(value)
    return Fraction(int(value.numerator), int(value.denominator))


def fixed(integer: int, places: int = DECIMAL_PLACES) -> str:
    sign = "-" if integer < 0 else ""
    digits = str(abs(integer)).rjust(places + 1, "0")
    return f"{sign}{digits[:-places]}.{digits[-places:]}"


def decimal_bracket(value, places: int = DECIMAL_PLACES) -> dict[str, object]:
    """Outward decimal enclosure, including a safe downward floor."""
    value = mpq(value)
    scale = 10**places
    numerator = int(value.numerator) * scale
    denominator = int(value.denominator)
    lower = numerator // denominator
    upper = -((-numerator) // denominator)
    assert mpq(lower, scale) <= value <= mpq(upper, scale)
    return {"places": places, "floor": fixed(lower, places), "ceiling": fixed(upper, places)}


def exact_value(value) -> dict[str, object]:
    return {"exact": rat(value), "decimal_enclosure": decimal_bracket(value)}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def subsets(mask: int):
    subset = mask
    while True:
        yield subset
        if subset == 0:
            break
        subset = (subset - 1) & mask


def independent_kernel(death: tuple[mpq, ...], colonisation: mpq) -> list[list[mpq]]:
    """Exact labelled kernel from a closed sum over survivor subsets.

    For old occupied set A, destination B and exact survivor set U subset A,
    the term is
      prod_U(1-d_i) prod_(A\\U)d_i * r_U^|B\\U| * (1-r_U)^(N-|B|),
    provided U subset B, with r_U = 1-(1-c)^|U|.
    """
    zero, one = mpq(0), mpq(1)
    assert len(death) == N and all(zero <= d <= one for d in death)
    assert zero <= colonisation <= one
    kernel = [[zero for _ in range(S)] for _ in range(S)]
    for old in range(S):
        for dest in range(S):
            probability = zero
            for survivors in subsets(old):
                if survivors & ~dest:
                    continue
                survival_term = one
                for i, d in enumerate(death):
                    if old & (1 << i):
                        survival_term *= one - d if survivors & (1 << i) else d
                count_survivors = survivors.bit_count()
                recolonisation = one - (one - colonisation) ** count_survivors
                count_colonised = dest.bit_count() - count_survivors
                count_uncolonised = N - dest.bit_count()
                probability += (
                    survival_term
                    * recolonisation**count_colonised
                    * (one - recolonisation) ** count_uncolonised
                )
            kernel[old][dest] = probability
    return kernel


def action_parameters(action: str, hit: tuple[int, ...]) -> tuple[tuple[mpq, ...], mpq]:
    excess = E1 - E0
    if action == "C":
        return tuple(E0 + excess * x for x in hit), 2 * C
    site = int(action[1:])
    death = tuple(
        E0 + excess * x * (mpq(1, 2) if i == site else mpq(1))
        for i, x in enumerate(hit)
    )
    return death, C


def deterministic_action_kernel(action: str, hit: tuple[int, ...]) -> list[list[mpq]]:
    death, colonisation = action_parameters(action, hit)
    return independent_kernel(death, colonisation)


def weighted_kernel(law: tuple[mpq, ...], action: str) -> list[list[mpq]]:
    zero = mpq(0)
    out = [[zero for _ in range(S)] for _ in range(S)]
    for weight, hit in zip(law, PATTERNS):
        if weight == 0:
            continue
        kernel = deterministic_action_kernel(action, hit)
        for old in range(S):
            for dest in range(S):
                out[old][dest] += weight * kernel[old][dest]
    return out


def check_kernel(kernel: list[list[mpq]]) -> None:
    zero, one = mpq(0), mpq(1)
    assert all(value >= zero for row in kernel for value in row)
    assert all(sum(row, zero) == one for row in kernel)
    assert kernel[0][0] == one and all(kernel[0][j] == zero for j in range(1, S))


def persistence(kernel: list[list[mpq]], horizon: int = T) -> mpq:
    zero, one = mpq(0), mpq(1)
    vector = [zero for _ in range(S)]
    vector[-1] = one
    for _ in range(horizon):
        vector = [sum((vector[i] * kernel[i][j] for i in range(S)), zero) for j in range(S)]
    # Direct live-mass summation avoids a subtraction from absorbing mass.
    return sum(vector[1:], zero)


def average_and_persist(law: tuple[mpq, ...]) -> dict[str, mpq]:
    values = {}
    for action in ("C", "R0", "R1", "R2", "R3"):
        kernel = weighted_kernel(law, action)
        check_kernel(kernel)
        values[action] = persistence(kernel)
    return values


def constraint_matrix() -> list[list[mpq]]:
    rows = [[mpq(1) for _ in PATTERNS]]
    rows.extend([[mpq(x[i]) for x in PATTERNS] for i in range(N)])
    rows.extend(
        [[mpq(x[i] * x[j]) for x in PATTERNS] for i, j in itertools.combinations(range(N), 2)]
    )
    return rows


def exact_rank(matrix: list[list[mpq]]) -> int:
    """Independent exact Gaussian rank; no call to the probe vertex routine."""
    if not matrix:
        return 0
    work = [list(map(mpq, row)) for row in matrix]
    nrow, ncol = len(work), len(work[0])
    pivot_row = 0
    for col in range(ncol):
        pivot = next((r for r in range(pivot_row, nrow) if work[r][col] != 0), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][col]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for r in range(nrow):
            if r != pivot_row and work[r][col] != 0:
                factor = work[r][col]
                work[r] = [a - factor * b for a, b in zip(work[r], work[pivot_row])]
        pivot_row += 1
        if pivot_row == nrow:
            break
    return pivot_row


def parse_law(saved: dict[str, str]) -> tuple[mpq, ...]:
    law = [mpq(0) for _ in PATTERNS]
    for label, value in saved.items():
        law[PATTERN_INDEX[label]] = q(value)
    return tuple(law)


def iid_mixture_pattern_law(atoms: tuple[mpq, ...], weights: tuple[mpq, ...]) -> tuple[mpq, ...]:
    law = []
    for hit in PATTERNS:
        k = sum(hit)
        law.append(sum((w * theta**k * (1 - theta) ** (N - k) for theta, w in zip(atoms, weights)), mpq(0)))
    assert sum(law, mpq(0)) == 1
    return tuple(law)


def law_constraints(law: tuple[mpq, ...], target_q: mpq) -> dict[str, object]:
    total = sum(law, mpq(0))
    marginals = [sum((law[k] * PATTERNS[k][i] for k in range(S)), mpq(0)) for i in range(N)]
    pairs = {
        f"{i}{j}": sum((law[k] * PATTERNS[k][i] * PATTERNS[k][j] for k in range(S)), mpq(0))
        for i, j in itertools.combinations(range(N), 2)
    }
    triples = {
        "".join(map(str, ijk)): sum(
            (law[k] * PATTERNS[k][ijk[0]] * PATTERNS[k][ijk[1]] * PATTERNS[k][ijk[2]] for k in range(S)),
            mpq(0),
        )
        for ijk in itertools.combinations(range(N), 3)
    }
    all_four = law[-1]
    by_count = [sum((law[k] for k, x in enumerate(PATTERNS) if sum(x) == count), mpq(0)) for count in range(N + 1)]
    exchangeable = all(
        len({law[k] for k, x in enumerate(PATTERNS) if sum(x) == count}) == 1
        for count in range(N + 1)
    )
    support = [k for k, weight in enumerate(law) if weight != 0]
    full = constraint_matrix()
    support_matrix = [[row[k] for k in support] for row in full]
    full_rank = exact_rank(full)
    support_rank = exact_rank(support_matrix)
    vertex = support_rank == len(support)
    assert total == 1 and all(value >= 0 for value in law)
    assert all(value == P for value in marginals)
    assert all(value == target_q for value in pairs.values())
    assert full_rank == 11 and support_rank == len(support)
    assert exchangeable
    triple_values = set(triples.values())
    assert len(triple_values) == 1
    m3 = next(iter(triple_values))
    count_formula = (
        1 - 4 * P + 6 * target_q - 4 * m3 + all_four,
        4 * P - 12 * target_q + 12 * m3 - 4 * all_four,
        6 * target_q - 12 * m3 + 6 * all_four,
        4 * m3 - 4 * all_four,
        all_four,
    )
    assert tuple(by_count) == count_formula
    return {
        "support": {"size": len(support), "patterns": ["".join(map(str, PATTERNS[k])) for k in support]},
        "probability": {"nonnegative": True, "total": rat(total)},
        "marginals": [rat(value) for value in marginals],
        "pairs": {key: rat(value) for key, value in pairs.items()},
        "specific_triples": {key: rat(value) for key, value in triples.items()},
        "all_four": rat(all_four),
        "hit_count_probabilities": [rat(value) for value in by_count],
        "exchangeable": exchangeable,
        "vertex_check": {
            "full_equality_matrix_rank": full_rank,
            "affine_nullity": S - full_rank,
            "support_column_rank": support_rank,
            "support_size": len(support),
            "is_vertex": vertex,
            "criterion": "feasible support columns of the equality matrix are linearly independent",
        },
    }


def scalar_moment_lens(target_q: mpq) -> dict[str, object]:
    x_lo = target_q**2 / P
    x_hi = target_q - (P - target_q) ** 2 / (1 - P)
    lower_y = lambda x: (x**2 - 2 * P * target_q * x + target_q**3) / (target_q - P**2)
    upper_y = lambda x: x - (target_q - x) ** 2 / (P - target_q)
    y_lo, y_hi = lower_y(x_lo), upper_y(x_hi)
    return {
        "m3_specific_triple_range": [rat(x_lo), rat(x_hi)],
        "m4_all_four_global_range": [rat(y_lo), rat(y_hi)],
        "basis": "Hausdorff moment feasibility for E(theta)=p, E(theta^2)=q, theta in [0,1]",
    }


def witness_audit(name: str, law: tuple[mpq, ...], target_q: mpq, saved: dict) -> dict[str, object]:
    constraints = law_constraints(law, target_q)
    lens = scalar_moment_lens(target_q)
    m3 = q(next(iter(constraints["specific_triples"].values())))
    m4 = q(constraints["all_four"])
    x_lo, x_hi = map(q, lens["m3_specific_triple_range"])
    y_lo, y_hi = map(q, lens["m4_all_four_global_range"])
    assert not (x_lo <= m3 <= x_hi)
    assert not (y_lo <= m4 <= y_hi)

    values = average_and_persist(law)
    prs = [values[f"R{i}"] for i in range(N)]
    assert len(set(prs)) == 1
    pc, pr = values["C"], prs[0]
    delta = pr - pc
    saved_pr = [q(value) for value in saved["PR"]]
    assert all(abs(float(observed) - expected) < 5e-15 for observed, expected in zip(prs, saved_pr))
    assert abs(float(pc) - saved["PC"]) < 5e-15
    assert abs(float(delta) - saved["Delta_1"]) < 5e-15
    return {
        "name": name,
        "law": {"".join(map(str, PATTERNS[k])): rat(weight) for k, weight in enumerate(law) if weight},
        "constraints": constraints,
        "scalar_hazard_lens": lens,
        "outside_scalar_hazard_lens": {
            "specific_triple": not (x_lo <= m3 <= x_hi),
            "all_four": not (y_lo <= m4 <= y_hi),
            "therefore_no_conditionally_iid_scalar_hazard_representation": True,
        },
        "persistence": {
            "PC": exact_value(pc),
            "PR_by_buffer_site": [exact_value(value) for value in prs],
            "all_buffer_sites_equal": len(set(prs)) == 1,
            "Delta_PR_minus_PC": exact_value(delta),
            "matches_saved_float_summary_within_5e-15": True,
        },
    }


def pair_decision_audit(low: dict, high: dict) -> dict[str, object]:
    low_delta = q(low["persistence"]["Delta_PR_minus_PC"]["exact"])
    high_delta = q(high["persistence"]["Delta_PR_minus_PC"]["exact"])
    assert low_delta < 0 < high_delta
    regret_if_r = -low_delta
    regret_if_c = high_delta
    deterministic = min(regret_if_r, regret_if_c)
    deterministic_action = "R (any buffer site)" if regret_if_r <= regret_if_c else "C"
    alpha_r = high_delta / (high_delta - low_delta)
    lottery = (-low_delta * high_delta) / (high_delta - low_delta)
    assert alpha_r * (-low_delta) == (1 - alpha_r) * high_delta == lottery
    return {
        "deterministic_two_law_minimax": {
            "regret_if_R": exact_value(regret_if_r),
            "regret_if_C": exact_value(regret_if_c),
            "best_action_for_pair": deterministic_action,
            "regret": exact_value(deterministic),
            "safe_decimal_floor": decimal_bracket(deterministic)["floor"],
        },
        "once_only_lottery_two_law_minimax": {
            "probability_R": exact_value(alpha_r),
            "probability_C": exact_value(1 - alpha_r),
            "regret": exact_value(lottery),
            "safe_decimal_floor": decimal_bracket(lottery)["floor"],
        },
        "scope": {
            "pair_value_is_lower_bound_for_full_labelled_class": True,
            "buffer_site_randomisation_cannot_improve_pair_value": True,
            "reason": "both witness laws are exchangeable and give identical persistence at every buffer site",
        },
    }


def source_kernel_validation() -> dict[str, object]:
    """Compare every deterministic pattern/action kernel to the saved constructor exactly."""
    comparisons = 0
    max_difference = mpq(0)
    for action in ("C", "R0", "R1", "R2", "R3"):
        for hit in PATTERNS:
            death, colonisation = action_parameters(action, hit)
            ours = independent_kernel(death, colonisation)
            source = probe.kernel_det([frac(value) for value in death], frac(colonisation))
            for old in range(S):
                for dest in range(S):
                    difference = abs(ours[old][dest] - q(source[old][dest]))
                    max_difference = max(max_difference, difference)
            comparisons += 1
    assert comparisons == 80 and max_difference == 0
    return {
        "deterministic_kernels_compared": comparisons,
        "entries_per_kernel": S * S,
        "exact_max_absolute_difference": rat(max_difference),
        "exact_match": True,
        "independence": "closed survivor-subset transition sum versus probe nested survivor/colonisation enumeration",
        "probability_type": "gmpy2.mpq from the first zero/one and throughout",
    }


def control_audit(boundary: dict) -> dict[str, object]:
    """Reproduce Glo/Ghi in labelled form and prove exact count-chain lumping."""
    controls = {
        "Glo": ((mpq(0), mpq(3, 10)), (mpq(2, 3), mpq(1, 3)), "min_certificate"),
        "Ghi": ((mpq(7, 90), mpq(1)), (mpq(81, 83), mpq(2, 83)), "max_certificate"),
    }
    out = {}
    delta_certificate = boundary["whole_class"]["Delta"]
    for name, (atoms, weights, certificate_name) in controls.items():
        law = iid_mixture_pattern_law(atoms, weights)
        values = average_and_persist(law)
        pc, pr = values["C"], values["R0"]
        assert len({values[f"R{i}"] for i in range(N)}) == 1
        delta = pr - pc

        atoms_f = tuple(frac(value) for value in atoms)
        weights_f = tuple(frac(value) for value in weights)
        count_pc = vertical.persist(vertical.avg_direct("C", N, frac(E0), frac(E1), frac(C), atoms_f, weights_f), vertical.start_state("C", N))
        count_pr = vertical.persist(vertical.avg_direct("R", N, frac(E0), frac(E1), frac(C), atoms_f, weights_f), vertical.start_state("R", N))
        assert pc == q(count_pc) and pr == q(count_pr)

        certified = q(delta_certificate[certificate_name]["attained"])
        assert delta == certified
        moments = [sum((weight * atom**power for atom, weight in zip(atoms, weights)), mpq(0)) for power in range(5)]
        out[name] = {
            "theta_atoms": [rat(value) for value in atoms],
            "weights": [rat(value) for value in weights],
            "theta_moments_m0_to_m4": [rat(value) for value in moments],
            "PC": exact_value(pc),
            "PR": exact_value(pr),
            "Delta_PR_minus_PC": exact_value(delta),
            "exact_match_to_n4_vertical_count_chain": True,
            "exact_match_to_boundary_certificate_attained_value": True,
        }
    return {
        "controls": out,
        "shared_dependencies": [
            "same rational Glo/Ghi theta laws and annual iid law",
            "same e0=1/50, e1=9/10, c=1/10 action definitions",
            "same death-then-recolonisation observation/process order",
            "same all-occupied start, absorbing zero state and T=50 horizon",
        ],
        "algorithmic_comparison": "independent 16-state labelled kernel versus n4_vertical.avg_direct using n3_certificate count/flag kernels",
    }


def main() -> dict[str, object]:
    saved = json.loads(PROBE_JSON.read_text())
    boundary_results = json.loads(BOUNDARY_JSON.read_text())
    selected = {
        row["q"]: row
        for row in saved["rows"]
        if row["q"] in {"3/100", "3/50"}
        and (row["e0"], row["e1"], row["c"]) == ("1/50", "9/10", "1/10")
    }
    assert set(selected) == {"3/100", "3/50"}
    reference_boundary = next(
        row
        for row in boundary_results["cases"]
        if (row["q"], row["e0"], row["e1"], row["c"]) == ("3/100", "1/50", "9/10", "1/10")
    )

    cases = []
    for q_label in ("3/100", "3/50"):
        row = selected[q_label]
        exact = row["exact_vertex_reevaluation"]
        assert exact["min_D1"] == exact["min_Dstar"]
        assert exact["max_D1"] == exact["max_Dstar"]
        target_q = q(q_label)
        low = witness_audit(f"q={q_label} low-Delta witness", parse_law(exact["min_D1"]["Q"]), target_q, exact["min_D1"])
        high = witness_audit(f"q={q_label} high-Delta witness", parse_law(exact["max_D1"]["Q"]), target_q, exact["max_D1"])
        cases.append(
            {
                "q": q_label,
                "parameters": {"p": rat(P), "e0": rat(E0), "e1": rat(E1), "c": rat(C), "n": N, "T": T},
                "witnesses": {"low": low, "high": high},
                "decision_regret": pair_decision_audit(low, high),
            }
        )

    kernel_validation = source_kernel_validation()
    controls = control_audit(reference_boundary)
    hashes = {relative: sha256(ROOT / relative) for relative in HASHED_SOURCES}
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, text=True, capture_output=True
    ).stdout.strip()
    result = {
        "audit_status": "PASS_WITH_REPORT_CORRECTIONS",
        "scientific_target": {
            "estimand": "50-step persistence under each complete action policy and Delta=PR-PC",
            "prediction_target": "ranking/regret for annual iid labelled hit laws satisfying all four marginals p and all six pair probabilities q",
            "scope_of_inference": "the four saved rational witnesses at n=4, T=50, e0=1/50, e1=9/10, c=1/10; no class-extremum claim",
            "observation_and_process": "one hit vector is drawn each year; deaths depend on that hit; survivors recolonise independently; extinction is absorbing",
        },
        "kernel_validation": kernel_validation,
        "original_scalar_hazard_controls": controls,
        "cases": cases,
        "policy_interpretation": {
            "Delta_star_definition": "max_i PR_i - PC is a law-specific best-buffer-site oracle comparison",
            "not_general_unknown_law_policy": True,
            "reason": "the maximising site can depend on the unknown law; implementable robust policies must fix or randomise sites before that law is known",
            "symmetric_witness_implication": "for these four witnesses PR_i is identical for every i, so neither choosing nor randomising among buffer sites removes either two-law regret bound",
        },
        "source_report_inaccuracies": [
            {
                "issue": "labelled_class_probe.md claims Glo/Ghi validation, but labelled_class_probe.py contains no such executed check",
                "resolution": "this audit independently constructs both pattern laws and matches the exact n4 boundary certificate and aggregated count chain",
            },
            {
                "issue": "labelled_class_probe.md says the full-labelled fixed-action objective has degree 2T",
                "resolution": "an averaged one-step kernel is linear in Q, so fixed-action persistence after T steps has degree at most T; Delta_star is a maximum of four such polynomials",
            },
            {
                "issue": "the q=3/50 low witness has hits capped at three, not at two",
                "resolution": "its support contains pair and triple patterns; its hit-count probabilities are recorded exactly above",
            },
        ],
        "limitations": [
            "The four witnesses prove attained sign reversals and regret lower bounds, not extrema over the five-dimensional labelled polytope.",
            "The audit retains the model assumption of independent, identically distributed annual hit laws; temporal dependence is outside scope.",
            "The exact controls share biological transition assumptions with the source count chain even though the state representation and transition calculation differ.",
        ],
        "verification": {
            "command": "uv run --python 3.12 --with gmpy2 --with numpy --with scipy python ideation/b_screen/labelled_class_audit.py",
            "assertions": [
                "80 deterministic kernels exactly equal the source kernel_det",
                "all averaged kernels are nonnegative, row-stochastic and absorb at zero",
                "all four laws exactly meet normalization, marginal and pair constraints",
                "all four support-column ranks establish vertices independently",
                "all four laws are exchangeable and lie outside both scalar-hazard moment ranges",
                "all four buffering sites have exactly equal persistence",
                "Glo/Ghi exactly match the count-chain persistence and saved boundary certificate",
            ],
            "git_head": head,
            "source_sha256": hashes,
        },
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    result = main()
    print(json.dumps({
        "audit_status": result["audit_status"],
        "git_head": result["verification"]["git_head"],
        "output": str(OUTPUT.relative_to(ROOT)),
        "kernel_exact_match": result["kernel_validation"]["exact_match"],
        "cases": [case["q"] for case in result["cases"]],
    }, indent=2))
