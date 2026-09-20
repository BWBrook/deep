"""Exact check of the best n=6 sixth-moment-only flip found by the float screen (e0=1/50, e1=9/10, c in {1/10, 1/5})."""
import sys; sys.path.insert(0, ".")
exec(open("n6_screen_fable.py").read().split("# interior (m3, m4) grid")[0])
import numpy as np, json
from fractions import Fraction as F
def build_exact(c_): return {a: V.moment_kernels(a, n, F(1,50), F(9,10), c_) for a in ("R","C")}
def delta_exact(P, mom):
    out = {}
    for a in ("R","C"):
        Ak = P[a]; S = len(Ak[0]); K = [[sum(mom[k]*Ak[k][r][s] for k in range(n+1)) for s in range(S)] for r in range(S)]
        assert all(sum(row) == 1 for row in K) and all(x >= 0 for row in K for x in row), "kernel"
        out[a] = V.persist(K, V.start_state(a, n))
    return out["R"] - out["C"]
def minors_psd(M):
    import itertools as it
    def det(A):
        if len(A) == 1: return A[0][0]
        return sum((-1)**j * A[0][j] * det([row[:j]+row[j+1:] for row in A[1:]]) for j in range(len(A)))
    k = len(M)
    return all(det([[M[i][j] for j in idx] for i in idx]) >= 0 for r in range(1,k+1) for idx in it.combinations(range(k), r))
xL, xU = V.xbounds(p, q); results = []
for c_ in (F(1,10), F(1,5)):
    P = build_exact(c_); Pf = {a: [np.array([[float(x) for x in row] for row in Ak]) for Ak in P[a]] for a in ("R","C")}
    def dfl(mom):
        out = {}
        for a in ("R","C"):
            K = sum(float(mom[k]) * Pf[a][k] for k in range(n+1)); v = np.zeros(len(K)); v[V.start_state(a, n)] = 1.0
            for _ in range(V.T): v = v @ K
            out[a] = 1 - v[0]
        return out["R"] - out["C"]
    best = None
    for tx in (F(i,10) for i in range(1,10)):
        m3 = xL + tx*(xU-xL); yL, yU = V.Lfun(m3,p,q), V.Ufun(m3,p,q)
        for ty in (F(i,8) for i in range(1,8)):
            m4 = yL + ty*(yU-yL); m = [F(1), p, q, m3, m4]
            A = [[m[1],m[2]],[m[2],m[3]]]; L5 = -schur(A, [m[3],m[4]], F(0))
            Bm = [[m[0]-m[1], m[1]-m[2]],[m[1]-m[2], m[2]-m[3]]]; U5 = m[4] + schur(Bm, [m[2]-m[3], m[3]-m[4]], F(0))
            if not (L5 < U5): continue
            for tz in (F(i,8) for i in range(1,8)):
                m5 = L5 + tz*(U5-L5); mm = m + [m5]
                if not feasible5(mm): continue
                L6, U6 = m6_bounds(mm)
                if not (L6 < U6): continue
                a, b = dfl(mm+[L6]), dfl(mm+[U6])
                if a*b < 0 and (best is None or min(abs(a),abs(b)) > best[0]): best = (min(abs(a),abs(b)), mm, L6, U6, a, b)
    if best is None: print(f"c={c_}: no float flip"); continue
    _, mm, L6, U6, a, b = best
    eL, eU = delta_exact(P, mm+[L6]), delta_exact(P, mm+[U6])
    full_L, full_U = mm+[L6], mm+[U6]
    H = lambda full: [[full[i+j] for j in range(4)] for i in range(4)]; Hp = lambda full: [[full[i+j+1]-full[i+j+2] for j in range(3)] for i in range(3)]
    rec = {"c": str(c_), "m3_m4_m5": [str(x) for x in mm[3:]], "L6": str(L6), "U6": str(U6), "Delta_L6_exact": str(eL), "Delta_U6_exact": str(eU),
           "Delta_L6_float": float(eL), "Delta_U6_float": float(eU), "signs_opposite": (eL < 0 < eU) or (eU < 0 < eL),
           "L6_feasible_minors": minors_psd(H(full_L)) and minors_psd(Hp(full_L)), "U6_feasible_minors": minors_psd(H(full_U)) and minors_psd(Hp(full_U))}
    results.append(rec)
    print(f"c={c_}: m3..m5 = {rec['m3_m4_m5']}; m6 in [{float(L6):.8f}, {float(U6):.8f}]; EXACT Delta(L6) = {float(eL):+.10f}, Delta(U6) = {float(eU):+.10f}; opposite signs: {rec['signs_opposite']}; endpoints feasible (all principal minors >= 0): {rec['L6_feasible_minors']}, {rec['U6_feasible_minors']}")
json.dump(results, open("n6_witness_check.json", "w"), indent=1)
