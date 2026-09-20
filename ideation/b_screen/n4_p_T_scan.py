"""n=4 witness scan over the two untouched axes: marginal hazard p and horizon T.
For each (p, T), pairwise q is set to keep the implied pairwise correlation rho = (q - p^2)/(p - p^2) at 2/9,
the value of the reference case (p=0.1, q=0.03). e0 = 0.02. Grid e1 in [0.10, 0.95] (35), c in [0.01, 0.50] (50).
A point is marked when Delta = PR - PC has opposite signs under Glo and Ghi (sufficient witness of an unidentified
ranking over the class; grid counts are not areas). Floating point; exploratory."""
import sys, json, itertools; sys.path.insert(0, ".")
import numpy as np
import b_screen_v2 as B
RHO = 2/9; E0 = 0.02
E1 = np.linspace(0.10, 0.95, 35); CC = np.linspace(0.01, 0.50, 50)
out = {"rho": RHO, "e0": E0, "e1_grid": E1.tolist(), "c_grid": CC.tolist(), "panels": []}
print(f"{'p':>5} {'q':>7} {'T':>4} | reversal points / valid | fraction | e1 range | c range | opposite-direction")
for p in (0.05, 0.10, 0.20, 0.30):
    q = p*p + RHO*(p - p*p)
    Glo, Ghi = B.G_mix(0.0, p, q), B.G_mix(1.0, p, q)
    for T in (20, 50, 100, 200):
        grid = np.full((len(E1), len(CC)), -9, dtype=int)
        for i, e1 in enumerate(E1):
            if e1 <= E0: continue
            for j, c in enumerate(CC):
                _, PRl, PCl = B.run_case(4, E0, e1, c, *Glo, T=T); _, PRh, PCh = B.run_case(4, E0, e1, c, *Ghi, T=T)
                dl, dh = PRl - PCl, PRh - PCh
                grid[i, j] = 1 if (dl < 0 < dh) else 2 if (dh < 0 < dl) else 0
        valid = grid >= 0; nrev = int((grid > 0).sum()); idx = np.argwhere(grid > 0)
        e1r = [float(E1[idx[:,0]].min()), float(E1[idx[:,0]].max())] if len(idx) else None
        cr = [float(CC[idx[:,1]].min()), float(CC[idx[:,1]].max())] if len(idx) else None
        out["panels"].append({"p": p, "q": q, "T": T, "n_reversal": nrev, "n_valid": int(valid.sum()), "grid_fraction": nrev/valid.sum(), "e1_range": e1r, "c_range": cr, "opposite_direction": int((grid == 2).sum()), "grid": grid.tolist()})
        print(f"{p:>5} {q:>7.4f} {T:>4} | {nrev:>4} / {int(valid.sum()):>4} | {nrev/valid.sum():.3f} | {e1r} | {cr} | {int((grid==2).sum())}")
json.dump(out, open("n4_p_T_scan.json", "w"))
