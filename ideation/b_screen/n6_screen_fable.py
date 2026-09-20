"""n=6 screen: at fixed (m1..m5), does the sixth moment alone flip PR-PC? Float screen over interior points, then
exact check of the best witness. Even-degree Hausdorff conditions on [0,1] for moments to degree 6:
  H  = (m_{i+j})_{i,j=0..3} >= 0  (m6 at (3,3): gives LOWER bound on m6 via Schur complement)
  H' = (m_{i+j+1} - m_{i+j+2})_{i,j=0..2} >= 0  (m5 - m6 at (2,2): gives UPPER bound on m6)
Parameters as the n=4/5 witnesses: p=1/10, q=3/100, e0=1/50, e1=9/10, c=1/10, T=50."""
import sys, itertools; sys.path.insert(0, ".")
from fractions import Fraction as F
import numpy as np
import n4_vertical as V
n = 6; p, q = F(1,10), F(3,100); e0, e1, c = F(1,50), F(9,10), F(1,10)
parts = {a: V.moment_kernels(a, n, e0, e1, c) for a in ("R","C")}
partsf = {a: [np.array([[float(x) for x in row] for row in Ak]) for Ak in parts[a]] for a in ("R","C")}
def delta_float(mom):
    out = {}
    for a in ("R","C"):
        K = sum(float(mom[k]) * partsf[a][k] for k in range(n+1))
        v = np.zeros(len(K)); v[V.start_state(a, n)] = 1.0
        for _ in range(V.T): v = v @ K
        out[a] = 1 - v[0]
    return out["R"] - out["C"]
def schur(Mat, vec, corner):
    """corner - vec^T Mat^{-1} vec, exact"""
    M = np.array(Mat, dtype=object); k = len(vec)
    # exact inverse via Fraction Gaussian elimination
    A = [list(row) + [F(int(i==j)) for j in range(k)] for i, row in enumerate(Mat)]
    for col in range(k):
        piv = next(r for r in range(col, k) if A[r][col] != 0); A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]; A[col] = [x/pv for x in A[col]]
        for r in range(k):
            if r != col and A[r][col] != 0:
                f = A[r][col]; A[r] = [x - f*y for x, y in zip(A[r], A[col])]
    Inv = [row[k:] for row in A]
    return corner - sum(vec[i]*Inv[i][j]*vec[j] for i in range(k) for j in range(k))
def m6_bounds(m):  # m = [m0..m5]
    # H: rows/cols 0..3; m6 at (3,3); block = H[0:3,0:3], vec = (m3,m4,m5)
    Hb = [[m[i+j] for j in range(3)] for i in range(3)]
    L6 = -schur(Hb, [m[3], m[4], m[5]], F(0))   # m6 >= vec^T Hb^-1 vec  (schur returns 0 - vec^T Hb^-1 vec)
    assert L6 >= 0
    # H': entries m_{i+j+1}-m_{i+j+2}, i,j=0..2; (2,2) = m5 - m6; block i,j=0..1, vec = (m3-m4, m4-m5)
    Hp = [[m[i+j+1]-m[i+j+2] for j in range(2)] for i in range(2)]
    U6 = m[5] - (-schur(Hp, [m[3]-m[4], m[4]-m[5]], F(0)))  # m5 - m6 >= vec^T Hp^-1 vec  ->  m6 <= m5 - vec^T Hp^-1 vec
    return L6, U6
def feasible5(m):  # degree-5 odd conditions for m0..m5
    Hx = [[m[i+j+1] for j in range(3)] for i in range(3)]; H1 = [[m[i+j]-m[i+j+1] for j in range(3)] for i in range(3)]
    def psd3(M):
        d = [M[i][i] >= 0 for i in range(3)]
        d2 = [M[i][i]*M[j][j]-M[i][j]*M[j][i] >= 0 for i in range(3) for j in range(i+1,3)]
        det = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
        return all(d) and all(d2) and det >= 0
    return psd3(Hx) and psd3(H1)
# interior (m3, m4) grid in the n=4 lens, m5 grid in its [L5, U5] range, then m6 endpoints
xL, xU = V.xbounds(p, q)
flips = []; tested = 0
for tx in (F(i,8) for i in range(1,8)):
    m3 = xL + tx*(xU-xL)
    yL, yU = V.Lfun(m3,p,q), V.Ufun(m3,p,q)
    for ty in (F(i,6) for i in range(1,6)):
        m4 = yL + ty*(yU-yL)
        m = [F(1), p, q, m3, m4]
        # m5 range via Astra's Schur formulas
        A = [[m[1],m[2]],[m[2],m[3]]]; L5 = -schur(A, [m[3],m[4]], F(0))
        Bm = [[m[0]-m[1], m[1]-m[2]],[m[1]-m[2], m[2]-m[3]]]; U5 = m[4] - (-schur(Bm, [m[2]-m[3], m[3]-m[4]], F(0)))
        if not (L5 < U5): continue
        for tz in (F(i,6) for i in range(1,6)):
            m5 = L5 + tz*(U5-L5); mm = m + [m5]
            if not feasible5(mm): continue
            L6, U6 = m6_bounds(mm)
            if not (L6 < U6): continue
            tested += 1
            dL, dU = delta_float(mm + [L6]), delta_float(mm + [U6])
            if dL * dU < 0:
                flips.append((mm, L6, U6, dL, dU))
print(f"n=6: tested {tested} interior (m3,m4,m5) points; sign flips across the m6 range: {len(flips)}")
if flips:
    best = max(flips, key=lambda f: min(abs(f[3]), abs(f[4])))
    mm, L6, U6, dL, dU = best
    print("best witness (float screen): m3..m5 =", [str(x) for x in mm[3:]], "; m6 range [", L6, ",", U6, "]")
    print(f"  Delta at L6 = {dL:+.10f}, at U6 = {dU:+.10f}")
    # exact check of signs
    def delta_exact(mom):
        out = {}
        for a in ("R","C"):
            Ak = parts[a]; S = len(Ak[0])
            K = [[sum(mom[k]*Ak[k][r][s] for k in range(n+1)) for s in range(S)] for r in range(S)]
            assert all(sum(row) == 1 for row in K) and all(x >= 0 for row in K for x in row)
            out[a] = V.persist(K, V.start_state(a, n))
        return out["R"] - out["C"]
    eL, eU = delta_exact(mm + [L6]), delta_exact(mm + [U6])
    print(f"  EXACT: Delta at L6 = {float(eL):+.16f} (sign {'+' if eL>0 else '-'}), at U6 = {float(eU):+.16f} (sign {'+' if eU>0 else '-'})")
    # degree-6 feasibility of both endpoints: H (4x4) and H' (3x3) PSD via principal minors
    def minors_psd(M):
        import itertools as it
        k = len(M)
        for r in range(1, k+1):
            for idx in it.combinations(range(k), r):
                sub = [[M[i][j] for j in idx] for i in idx]
                # determinant exact (Bareiss-free simple expansion for k<=4)
                def det(A):
                    if len(A) == 1: return A[0][0]
                    return sum((-1)**j * A[0][j] * det([row[:j]+row[j+1:] for row in A[1:]]) for j in range(len(A)))
                if det(sub) < 0: return False
        return True
    for name, m6 in (("L6", L6), ("U6", U6)):
        full = mm + [m6]
        H = [[full[i+j] for j in range(4)] for i in range(4)]; Hp = [[full[i+j+1]-full[i+j+2] for j in range(3)] for i in range(3)]
        print(f"  {name}: all principal minors of H (4x4) >= 0: {minors_psd(H)}; of H' (3x3) >= 0: {minors_psd(Hp)}")
    print("  5-patch pattern probabilities E[th^k(1-th)^(5-k)] depend on m0..m5 only -> identical for both laws by construction.")
