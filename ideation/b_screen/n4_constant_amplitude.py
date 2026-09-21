"""n=4 constant-amplitude path: e1 = e0 + a, isolating background mortality from shock amplitude
(Astra, thread 49 #418: the fixed-e1 path compresses e1 - e0 as e0 rises).

Part 1 (float scan, exploratory): for a in {0.5, 0.7, 0.88}, grid e0 in [0, 1 - a] and c in [0.01, 0.50];
mark a cell when Delta = PR - PC has strictly opposite signs under the two endpoint laws Glo and Ghi, with
both |Delta| > 1e-9 (numerical margin; Astra #424 found a cancellation artefact in the unmargined scan).
Grid counts, not areas, not prevalence, not whole-class identification.

Part 2 (exact certificate): at fixed (a, c), certify e0 intervals on which Delta(Glo) < 0 < Delta(Ghi) by exact
rational interval branch and bound in e0. With a fixed, d(theta) = e0 + a*theta is affine in e0, so every
averaged kernel entry is a polynomial of degree <= 4 in e0; entry ranges over an e0 interval come from exact
degree-4 Bernstein coefficients (5-node interpolation), then interval propagation with simplex tightening as in
n4_rectangle_certificate.py. Leaf-level exact evidence is saved.
p = 1/10, q = 3/100, T = 50 throughout.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
from math import comb
import numpy as np
import b_screen_v2 as B
import n3_certificate as N3
import n4_rectangle_certificate as RC   # Minv, propagate, qq
from gmpy2 import mpq as Q

P, QQ = 0.1, 0.03
AS = (0.5, 0.7, 0.88)
CC = np.round(np.arange(0.01, 0.501, 0.01), 3)
MARGIN = 1e-9

def scan():
    Glo, Ghi = B.G_mix(0.0, P, QQ), B.G_mix(1.0, P, QQ)
    panels = []
    for a in AS:
        E0 = np.round(np.arange(0.0, 1 - a + 1e-9, 0.01), 3)
        grid = np.zeros((len(E0), len(CC)), dtype=int); dvals = np.zeros((len(E0), len(CC), 2))
        for i, e0 in enumerate(E0):
            e1 = e0 + a
            for j, c in enumerate(CC):
                _, PRl, PCl = B.run_case(4, e0, e1, c, *Glo); _, PRh, PCh = B.run_case(4, e0, e1, c, *Ghi)
                dl, dh = PRl - PCl, PRh - PCh; dvals[i, j] = (dl, dh)
                if abs(dl) > MARGIN and abs(dh) > MARGIN:
                    grid[i, j] = 1 if (dl < 0 < dh) else 2 if (dh < 0 < dl) else 0
        idx = np.argwhere(grid > 0)
        per_e0 = [(float(E0[i]), int((grid[i] > 0).sum())) for i in range(len(E0))]
        panel = {"a": a, "e0_grid": E0.tolist(), "n_reversal": int((grid > 0).sum()), "n_cells": int(grid.size),
                 "opposite_direction": int((grid == 2).sum()),
                 "e0_range": [float(E0[idx[:, 0]].min()), float(E0[idx[:, 0]].max())] if len(idx) else None,
                 "c_range": [float(CC[idx[:, 1]].min()), float(CC[idx[:, 1]].max())] if len(idx) else None,
                 "reversals_per_e0": per_e0, "grid": grid.tolist(),
                 "min_abs_delta_marked": float(np.min(np.abs(dvals[grid > 0]))) if len(idx) else None}
        panels.append(panel)
        print(f"a={a}: reversal cells {panel['n_reversal']}/{panel['n_cells']}; e0 range {panel['e0_range']}; c range {panel['c_range']}; opposite-direction {panel['opposite_direction']}")
        print("   e0 :", " ".join(f"{e:.2f}" for e, _ in per_e0[::max(1, len(per_e0)//12)]))
        print("   rev:", " ".join(f"{k:>4d}" for _, k in per_e0[::max(1, len(per_e0)//12)]))
    return panels

# ---------- exact 1-D certificate in e0 at fixed (a, c) ----------
def kernel_bounds_e0(action, a, c, lo, hi, atoms, weights):
    samples = [RC.avg_kernel(action, lo + F(i, 4) * (hi - lo) + a, c, lo + F(i, 4) * (hi - lo), atoms, weights) for i in range(5)]
    # avg_kernel(action, e1, c, e0, atoms, weights): e1 = e0 + a
    S = len(samples[0]); Klo = [[None] * S for _ in range(S)]; Khi = [[None] * S for _ in range(S)]
    for r in range(S):
        for s in range(S):
            vals = [RC.qq(samples[i][r][s]) for i in range(5)]
            beta = [sum(RC.Minv[k][i] * vals[i] for i in range(5)) for k in range(5)]
            Klo[r][s] = max(Q(0), min(beta)); Khi[r][s] = min(Q(1), max(beta))
    return Klo, Khi

def delta_interval_e0(a, c, lo, hi, atoms, weights):
    KloR, KhiR = kernel_bounds_e0("R", a, c, lo, hi, atoms, weights); KloC, KhiC = kernel_bounds_e0("C", a, c, lo, hi, atoms, weights)
    PRlo, PRhi = RC.propagate(KloR, KhiR, 2 * RC.n - 1); PClo, PChi = RC.propagate(KloC, KhiC, RC.n)
    return PRlo - PChi, PRhi - PClo

def certify_e0(a, c, lo, hi, q=F(3, 100), maxdepth=14):
    (tl, wl), (tu, wu) = N3.laws(F(1, 10), q)
    stack = [((lo, hi), 0)]; certified = []; unresolved = []; t0 = time.time(); evals = 0
    while stack:
        (u0, u1), depth = stack.pop(); evals += 1
        glo = delta_interval_e0(a, c, u0, u1, tl, wl); ghi = delta_interval_e0(a, c, u0, u1, tu, wu)
        if glo[1] < 0 and ghi[0] > 0:
            certified.append(((u0, u1), glo, ghi)); continue
        if glo[0] > 0 or ghi[1] < 0:
            unresolved.append(((u0, u1), "whole-interval-wrong-sign", glo, ghi)); continue
        if depth >= maxdepth:
            unresolved.append(((u0, u1), "depth", glo, ghi)); continue
        um = (u0 + u1) / 2; stack += [((u0, um), depth + 1), ((um, u1), depth + 1)]
    length_c = sum(u1 - u0 for (u0, u1), _, _ in certified)
    # merge certified intervals into maximal contiguous runs
    runs = []
    for (u0, u1), _, _ in sorted(certified):
        if runs and runs[-1][1] == u0: runs[-1][1] = u1
        else: runs.append([u0, u1])
    return dict(a=str(a), c=str(c), e0_interval=[str(lo), str(hi)], maxdepth=maxdepth, evaluations=evals,
                certified_length_fraction=str(F(length_c / (hi - lo))), certified_runs=[[str(x), str(y)] for x, y in runs],
                certified_runs_float=[[float(x), float(y)] for x, y in runs],
                certified_leaves=[{"e0": [str(u0), str(u1)], "delta_glo": [str(g) for g in glo], "delta_ghi": [str(g) for g in ghi]} for (u0, u1), glo, ghi in certified],
                unresolved=[{"e0": [str(u0), str(u1)], "why": why, "delta_glo": [float(g) for g in glo], "delta_ghi": [float(g) for g in ghi]} for (u0, u1), why, glo, ghi in unresolved],
                seconds=round(time.time() - t0, 1))

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "scan"
    out = {"p": "1/10", "q": "3/100", "T": B.T, "margin": MARGIN}
    if mode in ("scan", "all"):
        out["scan"] = scan()
        json.dump(out, open(os.path.join(os.path.dirname(__file__), "n4_constant_amplitude_scan.json"), "w"))
    if mode in ("cert", "all"):
        # certification targets: chosen from the scan (exploratory selection, then exact verification)
        targets = json.loads(sys.argv[2]) if len(sys.argv) > 2 else [["7/10", "3/10", "0", "3/10"], ["1/2", "3/10", "0", "1/2"], ["22/25", "3/10", "0", "3/25"]]
        certs = []
        for a, c, lo, hi in targets:
            r = certify_e0(F(a), F(c), F(lo), F(hi))
            certs.append(r)
            print(f"a={a} c={c} e0 in [{lo},{hi}]: certified length fraction {float(F(r['certified_length_fraction'])):.4f}; runs {r['certified_runs_float']}; unresolved {len(r['unresolved'])}; {r['seconds']}s", flush=True)
        json.dump(certs, open(os.path.join(os.path.dirname(__file__), "n4_constant_amplitude_cert.json"), "w"), indent=1)
