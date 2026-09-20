"""Two-known-law discrimination for the equal-triples witness; no fitted prior.

Equal prior probabilities and zero-one classification loss. Enumerate every
multinomial count vector for a complete iid yearly K=0..4 record, then compare
with retaining only the number of K=4 years. Finite sums use floating point,
not Monte Carlo; the underlying one-year probabilities are exact rationals.
Run: uv run --python 3.12 --with numpy python ideation/b_screen/n4_observability.py
"""
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np


def hit_count_law(moments):
    n = len(moments) - 1
    probabilities = [math.comb(n, k) * sum(
        (-1)**j * math.comb(n-k, j) * moments[k+j]
        for j in range(n-k+1)) for k in range(n+1)]
    assert sum(probabilities) == 1 and min(probabilities) > 0
    return probabilities


def full_record_error(years, p, q):
    """Equal-prior Bayes error = one half the sum of the smaller likelihoods."""
    lp, lq = np.log(p), np.log(q)
    logfact = np.array([math.lgamma(i+1) for i in range(years+1)])
    triangles = {r: np.tril_indices(r+1) for r in range(years+1)}
    error = 0.0
    masses = np.zeros(2)
    for a in range(years+1):
        for b in range(years-a+1):
            r = years-a-b
            row, column = triangles[r]
            c, d, e = column, row-column, r-row
            multiplier = (logfact[years]-logfact[a]-logfact[b]
                          -logfact[c]-logfact[d]-logfact[e])
            l = np.exp(multiplier+a*lp[0]+b*lp[1]+c*lp[2]+d*lp[3]+e*lp[4])
            u = np.exp(multiplier+a*lq[0]+b*lq[1]+c*lq[2]+d*lq[3]+e*lq[4])
            error += 0.5 * np.minimum(l, u).sum()
            masses += [l.sum(), u.sum()]
    assert np.max(np.abs(masses-1)) < 2e-11
    assert 0 <= error <= 0.5+2e-11
    return float(error), masses.tolist()


def all_four_error(years, p4, q4):
    return 0.5 * sum(min(
        math.comb(years, k)*p4**k*(1-p4)**(years-k),
        math.comb(years, k)*q4**k*(1-q4)**(years-k))
        for k in range(years+1))


def main():
    common = [F(1), F(1, 10), F(3, 100), F(21, 1000)]
    lower = hit_count_law(common+[F(171, 10000)])
    upper = hit_count_law(common+[F(1389, 70000)])
    p, q = np.array(lower, dtype=float), np.array(upper, dtype=float)
    # Independent closed-form one-year check and a no-information control.
    one_year = 0.5 * sum(min(a, b) for a, b in zip(lower, upper))
    assert abs(full_record_error(1, p, q)[0]-float(one_year)) < 1e-14
    assert abs(full_record_error(10, p, p)[0]-0.5) < 1e-13
    assert abs(full_record_error(10, p[::-1], q[::-1])[0]
               -full_record_error(10, p, q)[0]) < 1e-13
    rows = []
    for years in (1, 25, 50, 100, 200):
        full, masses = full_record_error(years, p, q)
        reduced = all_four_error(years, p[-1], q[-1])
        assert full <= reduced+2e-11  # discarding observations cannot improve optimal error
        rows.append({"years": years, "full_hit_count_error": full,
                     "all_four_only_error": reduced,
                     "enumerated_probability_masses": masses})
    result = {
        "question": "Discriminate two known laws with equal prior odds and zero-one loss",
        "assumptions": ["iid years", "perfect observation of all four patch hits",
                        "candidate laws known", "no nuisance-parameter estimation"],
        "method": "multinomial likelihood enumeration in floating point; no Monte Carlo",
        "lower_K_probabilities_exact": list(map(str, lower)),
        "upper_K_probabilities_exact": list(map(str, upper)),
        "one_year_total_variation_exact": str(sum(abs(a-b) for a, b in zip(lower, upper))/2),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "numpy_version": np.__version__,
        "results": rows,
    }
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(result, indent=2)+"\n")
    for row in rows:
        print(row["years"], row["full_hit_count_error"], row["all_four_only_error"])


if __name__ == "__main__":
    main()
