"""n=6: endpoint swing of Delta across the admissible m6 range at fixed (m1..m5), for nine (e1, c) pairs, e0=1/50,
p=1/10, q=3/100, T=50; 441 interior (m3,m4,m5) points. Float scan (numpy); reports flips, max endpoint swing,
minimum |Delta| at endpoints. Also the like-for-like n=5 (m5 swing at fixed m3,m4) and n=4 (m4 swing at fixed m3)
endpoint swings at the witness parameters e1=9/10, c=1/10, for a matched-lower-moments comparison across n.
Endpoint spans on sampled fibres are not global bounds (Astra, thread 49 #418)."""
import sys, json, itertools; sys.path.insert(0, ".")
from fractions import Fraction as F
import numpy as np
import n4_vertical as V
exec(open("n6_screen_fable.py").read().split("# interior (m3, m4) grid")[0].split("parts = {a: V.moment_kernels")[0])  # constants + schur/feasible5/m6_bounds defs come below
exec("\n".join(l for l in open("n6_screen_fable.py").read().split("# interior (m3, m4) grid")[0].splitlines() if l.startswith("def ") or l.startswith("    ") or l.startswith("        ") or l.startswith("            ")))
n = 6; p, q = F(1,10), F(3,100); e0 = F(1,50)
def kernels_f(nn, e1_, c_):
    P = {a: V.moment_kernels(a, nn, e0, e1_, c_) for a in ("R","C")}
    return {a: [np.array([[float(x) for x in row] for row in Ak]) for Ak in P[a]] for a in ("R","C")}
def delta_f(nn, pf, mom):
    out = {}
    for a in ("R","C"):
        K = sum(float(mom[k]) * pf[a][k] for k in range(nn+1)); v = np.zeros(len(K)); v[V.start_state(a, nn)] = 1.0
        for _ in range(V.T): v = v @ K
        out[a] = 1 - v[0]
    return out["R"] - out["C"]
xL, xU = V.xbounds(p, q); pts = []; pts5 = []; pts4 = []
for tx in (F(i,10) for i in range(1,10)):
    m3 = xL + tx*(xU-xL); yL, yU = V.Lfun(m3,p,q), V.Ufun(m3,p,q); pts4.append(([F(1),p,q,m3], yL, yU))
    for ty in (F(i,8) for i in range(1,8)):
        m4 = yL + ty*(yU-yL); m = [F(1), p, q, m3, m4]
        A = [[m[1],m[2]],[m[2],m[3]]]; L5 = -schur(A, [m[3],m[4]], F(0))
        Bm = [[m[0]-m[1], m[1]-m[2]],[m[1]-m[2], m[2]-m[3]]]; U5 = m[4] + schur(Bm, [m[2]-m[3], m[3]-m[4]], F(0))
        if not (L5 < U5): continue
        pts5.append((m, L5, U5))
        for tz in (F(i,8) for i in range(1,8)):
            m5 = L5 + tz*(U5-L5); mm = m + [m5]
            if not feasible5(mm): continue
            L6, U6 = m6_bounds(mm)
            if L6 < U6: pts.append((mm, L6, U6))
out = {"n6_points": len(pts), "rows": []}
print(f"n=6 interior points: {len(pts)}")
for e1_, c_ in itertools.product((F(9,10), F(7,10), F(1,2)), (F(1,20), F(1,10), F(1,5))):
    pf = kernels_f(6, e1_, c_); fl = 0; sw = 0.0; mn = 9.0
    for mm, L6, U6 in pts:
        a, b = delta_f(6, pf, mm+[L6]), delta_f(6, pf, mm+[U6]); fl += int(a*b < 0); sw = max(sw, abs(b-a)); mn = min(mn, abs(a), abs(b))
    out["rows"].append({"e1": str(e1_), "c": str(c_), "flips": fl, "max_endpoint_swing": sw, "min_abs_delta_endpoints": mn})
    print(f"e1={float(e1_):.1f} c={float(c_):.2f}: flips {fl:>2}; max endpoint swing {sw:.6f}; min |Delta| {mn:.6f}")
# like-for-like at the witness parameters
pf5 = kernels_f(5, F(9,10), F(1,10)); sw5 = max(abs(delta_f(5, pf5, m+[U5]) - delta_f(5, pf5, m+[L5])) for m, L5, U5 in pts5)
pf4 = kernels_f(4, F(9,10), F(1,10)); sw4 = max(abs(delta_f(4, pf4, m+[yU]) - delta_f(4, pf4, m+[yL])) for m, yL, yU in pts4)
sw6 = [r["max_endpoint_swing"] for r in out["rows"] if r["e1"] == "9/10" and r["c"] == "1/10"][0]
out["matched_lower_moment_swings_e1_0.9_c_0.1"] = {"n4_m4_swing_max_over_9_m3": sw4, "n5_m5_swing_max_over_63_m3m4": sw5, "n6_m6_swing_max_over_441": sw6}
print(f"matched-lower-moment endpoint swings at e1=0.9, c=0.1: n=4 {sw4:.6f} (over 9 m3 values), n=5 {sw5:.6f} (63 points), n=6 {sw6:.6f} (441 points)")
json.dump(out, open("n6_swing_scan.json", "w"), indent=1)
