"""Do the complete-action gains P_R - P_0 and P_C - P_0 vary along the whole E boundary the way the four extremisers
suggest (connectivity's gain falling with co-hit concentration, buffering's roughly flat)? Exact evaluation at the
polygon vertices and at 9 interior points of every edge, for the two reference cases (e0=1/50, e1=9/10, c=1/10;
q=3/100 and 3/50), plus along the S lens boundary at the same x for comparison. Fable, 2026-09-25."""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
import n4_vertical as V
import exchangeable_certificate as EC

n = 4; P = V.P; e0, e1, c = F(1, 50), F(9, 10), F(1, 10)
out = {}
for q in (F(3, 100), F(3, 50)):
    parts = {a: V.affine_parts(a, n, e0, e1, c, P, q) for a in ("0", "R", "C")}
    def gains(x, y):
        Pv = {a: V.persist(V.kernel_at(*parts[a], x, y), V.start_state(a, n)) for a in ("0", "R", "C")}
        return float(Pv["R"] - Pv["0"]), float(Pv["C"] - Pv["0"]), float(Pv["0"])
    verts = EC.polygon(P, q); rows = []
    for i in range(len(verts)):
        p0, p1 = verts[i], verts[(i + 1) % len(verts)]
        for j in range(10):
            u = F(j, 10); x = p0[0] + u * (p1[0] - p0[0]); y = p0[1] + u * (p1[1] - p0[1])
            gR, gC, P0 = gains(x, y); rows.append({"edge": i, "u": str(u), "tau3": float(x), "tau4": float(y), "PR-P0": gR, "PC-P0": gC, "P0": P0})
    # S lens boundary at 11 x values
    xL, xU = V.xbounds(P, q); lens = []
    for j in range(11):
        x = xL + F(j, 10) * (xU - xL)
        for nm, yf in (("L", V.Lfun), ("U", V.Ufun)):
            gR, gC, P0 = gains(x, yf(x, P, q)); lens.append({"curve": nm, "tau3": float(x), "tau4": float(yf(x, P, q)), "PR-P0": gR, "PC-P0": gC, "P0": P0})
    out[str(q)] = {"vertices": [[str(a), str(b)] for a, b in verts], "polygon_boundary": rows, "lens_boundary": lens}
    print(f"\nq={q}: vertices {[(float(a), float(b)) for a, b in verts]}")
    print("  polygon boundary (edge, u, tau3, tau4): PR-P0, PC-P0, P0")
    for r in rows:
        if r["u"] in ("0", "1/2"): print(f"   e{r['edge']} u={r['u']:<3} ({r['tau3']:.4f},{r['tau4']:.4f}): {r['PR-P0']:+.4f} {r['PC-P0']:+.4f}  P0 {r['P0']:.4f}")
    rng = lambda key, rs: (min(r[key] for r in rs), max(r[key] for r in rs))
    print(f"  polygon boundary range: PR-P0 {rng('PR-P0', rows)}  PC-P0 {rng('PC-P0', rows)}")
    print(f"  S lens boundary range:  PR-P0 {rng('PR-P0', lens)}  PC-P0 {rng('PC-P0', lens)}")
    # monotonicity of PC-P0 in tau4 at fixed tau3 across the polygon: compare lower vs upper boundary at matching tau3
    for r in lens[::2][:0]: pass
json.dump(out, open(os.path.join(os.path.dirname(__file__), "baseline_gain_edges.json"), "w"), indent=1)
