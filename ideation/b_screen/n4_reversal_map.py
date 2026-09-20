"""n=4 witness-based reversal map over (e1, c), floating point, exploratory.

For each (q, e0) and a grid of (e1, c), compute Delta = PR - PC under the two endpoint laws Glo and Ghi
(the class extremizers of m3 and m4). Opposite signs at Glo and Ghi are a sufficient witness that the
action ranking is unidentified over the conditional-iid class at fixed (p, q); they are not necessary in
general, although on the 16 audited grid cases all class extrema were attained at Glo or Ghi (n4_audit.md).
Also computes, from n4_boundary_results.json, the m3 band over which the sign is ambiguous in the one
unidentified grid case. Uses the floating-point implementation b_screen_v2.run_case; this is a scan, not
a certificate.
"""
import sys, json, os
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import b_screen_v2 as B

P = 0.1
E1 = np.linspace(0.10, 0.95, 35)
CC = np.linspace(0.01, 0.30, 30)

def main():
    out = {"n": 4, "T": B.T, "p": P, "e1_grid": E1.tolist(), "c_grid": CC.tolist(), "panels": []}
    for q in (0.03, 0.06):
        for e0 in (0.02, 0.10):
            Glo, Ghi = B.G_mix(0.0, P, q), B.G_mix(1.0, P, q)
            grid = np.full((len(E1), len(CC)), -9, dtype=int)
            for i, e1 in enumerate(E1):
                if e1 <= e0:
                    continue
                for j, c in enumerate(CC):
                    _, PRl, PCl = B.run_case(4, e0, e1, c, *Glo)
                    _, PRh, PCh = B.run_case(4, e0, e1, c, *Ghi)
                    dl, dh = PRl - PCl, PRh - PCh
                    grid[i, j] = 1 if (dl < 0 < dh) else 2 if (dh < 0 < dl) else 0
            valid = grid >= 0
            frac = float((grid[valid] > 0).mean())
            idx = np.argwhere(grid > 0)
            panel = {"q": q, "e0": e0, "reversal_fraction": frac, "n_reversal": int((grid > 0).sum()), "n_valid": int(valid.sum()),
                     "opposite_direction_points": int((grid == 2).sum()),
                     "e1_range": [float(E1[idx[:, 0]].min()), float(E1[idx[:, 0]].max())] if len(idx) else None,
                     "c_range": [float(CC[idx[:, 1]].min()), float(CC[idx[:, 1]].max())] if len(idx) else None,
                     "grid": grid.tolist()}
            out["panels"].append(panel)
            print(f"q={q} e0={e0}: reversal fraction {frac:.3f} ({panel['n_reversal']} of {panel['n_valid']}); e1 range {panel['e1_range']}; c range {panel['c_range']}; opposite-direction {panel['opposite_direction_points']}")
            for i in range(0, len(E1), 5):
                print(f"   e1={E1[i]:.2f} " + "".join("X" if grid[i, j] == 1 else "o" if grid[i, j] == 2 else "." if grid[i, j] == 0 else " " for j in range(0, len(CC), 2)))
            print("   c: " + " ".join(f"{CC[j]:.2f}" for j in range(0, len(CC), 6)))
    # m3 ambiguity band in the unidentified grid case, from the certified boundary roots
    try:
        J = json.load(open(os.path.join(os.path.dirname(__file__), "n4_boundary_results.json")))
        cases = J.get("cases", [])
        for r in cases:
            if str(r.get("q")) in ("3/100", "0.03") and str(r.get("e0")) in ("1/50", "0.02") and str(r.get("e1")) in ("9/10", "0.9") and str(r.get("c")) in ("1/10", "0.1"):
                out["unidentified_case_roots"] = {k: r["curves"][k].get("Delta_roots_x") or r["curves"][k].get("delta_roots_x") for k in r["curves"]}
                print("unidentified case boundary roots in m3:", out["unidentified_case_roots"])
    except Exception as e:  # results schema may have changed in the audited version
        print("could not read boundary roots:", e)
    with open(os.path.join(os.path.dirname(__file__), "n4_reversal_map.json"), "w") as f:
        json.dump(out, f)

if __name__ == "__main__":
    main()
