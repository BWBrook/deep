"""Independent audit of Astra's n=5 equal-quadruples witness (thread 49, #406) via Fable's own exact
affine-in-moments kernel decomposition (n4_vertical.moment_kernels generalised to n=5)."""
import sys; sys.path.insert(0, ".")
from fractions import Fraction as F
import n4_vertical as V
n = 5; p, q = F(1,10), F(3,100); e0, e1, c = F(1,50), F(9,10), F(1,10)
m = [F(1), p, q, F(337,18000), F(773993,51840000)]
# Schur bounds for m5 (odd-degree Hausdorff on [0,1]): Hx = (m[i+j+1]) >= 0, H1x = (m[i+j]-m[i+j+1]) >= 0, i,j=0..2
def inv2(a,b,c_,d): det = a*d - b*c_; return [[d/det, -b/det], [-c_/det, a/det]]
A = inv2(m[1], m[2], m[2], m[3]); v = [m[3], m[4]]
L5 = sum(v[i]*A[i][j]*v[j] for i in range(2) for j in range(2))
Bm = inv2(m[0]-m[1], m[1]-m[2], m[1]-m[2], m[2]-m[3]); w = [m[2]-m[3], m[3]-m[4]]
U5 = m[4] - sum(w[i]*Bm[i][j]*w[j] for i in range(2) for j in range(2))
print("L5 =", L5, "==", F(1855440127,149299200000), L5 == F(1855440127,149299200000))
print("U5 =", U5, "==", F(1990575127,149299200000), U5 == F(1990575127,149299200000))
def minors_ok(M):
    d1 = all(M[i][i] >= 0 for i in range(3))
    d2 = all(M[i][i]*M[j][j] - M[i][j]*M[j][i] >= 0 for i in range(3) for j in range(i+1,3))
    det = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    return d1 and d2 and det >= 0, det
for name, m5 in (("L5", L5), ("U5", U5)):
    mm = m + [m5]
    Hx = [[mm[i+j+1] for j in range(3)] for i in range(3)]
    H1x = [[mm[i+j]-mm[i+j+1] for j in range(3)] for i in range(3)]
    ok1, d1 = minors_ok(Hx); ok2, d2 = minors_ok(H1x)
    print(f"{name}: Hx minors ok {ok1} (det {float(d1):.3e}); H1-x minors ok {ok2} (det {float(d2):.3e})")
# sanity: (m3, m4) inside the n=4 lens
xL, xU = V.xbounds(p, q); print("m3 in [xL,xU]:", xL <= m[3] <= xU, "; m4 in [L(m3),U(m3)]:", V.Lfun(m[3],p,q) <= m[4] <= V.Ufun(m[3],p,q))
# exact propagation with Fable's decomposition
parts = {a: V.moment_kernels(a, n, e0, e1, c) for a in ("R","C")}
def persist_at(a, mom):
    Ak = parts[a]; S = len(Ak[0])
    K = [[sum(mom[k]*Ak[k][r][s] for k in range(n+1)) for s in range(S)] for r in range(S)]
    assert all(abs(sum(row)-1) == 0 for row in K) and all(x >= 0 for row in K for x in row), "kernel check"
    return V.persist(K, V.start_state(a, n))
for name, m5 in (("L5", L5), ("U5", U5)):
    mom = m + [m5]
    d = persist_at("R", mom) - persist_at("C", mom)
    print(f"{name}: PR - PC = {float(d):+.16f}  sign {'+' if d > 0 else '-'}")
# equal-quadruples identity: all labelled 4-patch pattern probabilities depend only on m0..m4 -> identical by construction; verify numerically via the two reconstructed 3-atom laws Astra reported (rounded), degree-4 pattern probabilities agree to ~1e-16
import numpy as np
lo_t = np.array([0, .11946456, .83852155]); lo_w = np.array([.34315425, .62690347, .02994228])
up_t = np.array([.05987322, .43977956, 1.0]); up_w = np.array([.91220039, .07571303, .01208658])
for k in range(5):
    a = np.dot(lo_w, lo_t**k*(1-lo_t)**(4-k)); b = np.dot(up_w, up_t**k*(1-up_t)**(4-k))
    print(f"4-patch pattern prob, k={k} hits: lower {a:.9f} upper {b:.9f} diff {a-b:+.1e}")
print("5-way hit prob: lower", float(np.dot(lo_w, lo_t**5)), "upper", float(np.dot(up_w, up_t**5)), " (exact L5, U5:", float(L5), float(U5), ")")
