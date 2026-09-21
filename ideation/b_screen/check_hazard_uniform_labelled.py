"""Independent labelled-state check of DEEP's two-uniform witnesses."""
import json
from pathlib import Path
from fractions import Fraction as F

import numpy as np

ROOT = Path(__file__).resolve().parent
import b_screen_v2 as B

with (ROOT / "hazard_family_audit.json").open() as src:
    rows = {r["name"]: r for r in json.load(src)["families"]}

p, q, h = F(1, 10), F(3, 100), F(1, 1000)
e0, e1, c = 0.02, 0.9, 0.1
xg, wg = np.polynomial.legendre.leggauss(64)
worst = 0.0

for name, a in (("uniform_low", F(1, 100)),
                ("uniform_high", F(3, 40))):
    q_centres = q - h * h / 3
    b = (q_centres - p * a) / (p - a)
    weight_b = (p - a) / (b - a)
    theta = np.concatenate((float(a) + float(h) * xg,
                            float(b) + float(h) * xg))
    weights = np.concatenate((float(1 - weight_b) * wg / 2,
                              float(weight_b) * wg / 2))
    assert abs(weights.sum() - 1) < 3e-16
    assert abs(weights @ theta - 0.1) < 3e-16
    assert abs(weights @ theta**2 - 0.03) < 3e-16

    for n in (3, 4):
        full = (1 << n) - 1
        e1_buffered = e0 + B.BUFFER * (e1 - e0)
        kr = B.avg_kernel(
            lambda th: B.kernel_labelled(
                n,
                [B.d_of(th, e0, e1_buffered)]
                + [B.d_of(th, e0, e1)] * (n - 1),
                c,
            ),
            theta,
            weights,
        )
        kc = B.avg_kernel(
            lambda th: B.kernel_labelled(
                n, [B.d_of(th, e0, e1)] * n, B.KAPPA * c
            ),
            theta,
            weights,
        )
        B.check_kernel(kr)
        B.check_kernel(kc)
        pr = B.persist(kr, full, 0, T=50)
        pc = B.persist(kc, full, 0, T=50)
        got = {"R": pr, "C": pc, "delta": pr - pc}
        ref = rows[name]["results"][str(n)]["display"]
        discrepancy = max(abs(got[k] - ref[k]) for k in got)
        worst = max(worst, discrepancy)
        print(name, f"n={n}", f"delta={got['delta']:+.16g}",
              f"max_abs_discrepancy={discrepancy:.3e}")

assert worst < 2e-12
print(f"PASS worst_abs_discrepancy={worst:.3e}")
