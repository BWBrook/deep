"""Certified positive-area reversal rectangle at n=4 (thread 49: Astra's #406 asks for a rectangle first).

Claim to certify: for all (e1, c) in the rectangle R0, at fixed q, e0, p=1/10, n=4, T=50:
  Delta(Glo) = PR - PC under the lower endpoint law is < 0, and Delta(Ghi) under the upper endpoint law is > 0.
Both laws are fixed atomic laws, so the averaged kernels are polynomials in (e1, c) of degree <= 4 in each variable
(d is affine in e1 and appears to power <= n; recolonisation factors have degree m(n-m) <= 4 in c).

Method: exact rational interval arithmetic with branch and bound over sub-rectangles.
 - Kernel entry ranges over a sub-rectangle from exact tensor Bernstein coefficients of degree (4,4), obtained by
   interpolation of exact kernels at the 5x5 grid of nodes (i/4, j/4) in rectangle coordinates.
 - State-vector intervals propagated 50 steps: lo = lo_v . Klo, hi = hi_v . Khi (all quantities nonnegative), with the
   simplex tightening hi_i <= 1 - sum_{j != i} lo_j and clipping to [0,1].
 - Persistence P in [1 - hi_0, 1 - lo_0]; Delta interval from the two actions. A sub-rectangle is certified when the
   Delta(Glo) interval is strictly negative and the Delta(Ghi) interval strictly positive. Otherwise it is bisected on its
   longer side (relative), to a depth cap; leaves at the cap are reported as unresolved (not as counterexamples).
Output: certified area (exact rational) as a fraction of R0, leaf counts, runtime. Exact throughout (gmpy2.mpq).
"""
import sys, time, json, itertools
sys.path.insert(0, ".")
from fractions import Fraction as F
from math import comb
from gmpy2 import mpq as Q
import n3_certificate as N3

T = 50; n = 4; p = F(1,10)
def qq(f): return Q(f.numerator, f.denominator)

# exact 5x5 Bernstein-Vandermonde inverse for degree 4 on nodes i/4
Mb = [[Q(comb(4,k)) * Q(i,4)**k * (1 - Q(i,4))**(4-k) for k in range(5)] for i in range(5)]
def inv(M):
    nn = len(M); A = [row[:] + [Q(int(i==j)) for j in range(nn)] for i, row in enumerate(M)]
    for col in range(nn):
        piv = next(r for r in range(col, nn) if A[r][col] != 0); A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]; A[col] = [v / pv for v in A[col]]
        for r in range(nn):
            if r != col and A[r][col] != 0:
                f = A[r][col]; A[r] = [a - f*b for a, b in zip(A[r], A[col])]
    return [row[nn:] for row in A]
Minv = inv(Mb)

def avg_kernel(action, e1, c, e0, atoms, weights):
    K = None
    for th, wt in zip(atoms, weights):
        d = N3.d_of(th, e0, e1)
        Ki = N3.kernel_buffered(n, d, N3.d_of(th, e0, e0 + F(1,2)*(e1-e0)), c) if action == "R" else N3.kernel_hom(n, d, c * (2 if action == "C" else 1))
        K = [[wt*x for x in row] for row in Ki] if K is None else [[a + wt*b for a, b in zip(ra, rb)] for ra, rb in zip(K, Ki)]
    return K

def kernel_bounds(action, rect, e0, atoms, weights):
    """Entrywise [lo, hi] over rect = (a0, a1, c0, c1) via exact degree-(4,4) Bernstein coefficients."""
    a0, a1, c0, c1 = rect
    samples = [[avg_kernel(action, a0 + F(i,4)*(a1-a0), c0 + F(j,4)*(c1-c0), e0, atoms, weights) for j in range(5)] for i in range(5)]
    S = len(samples[0][0])
    lo = [[None]*S for _ in range(S)]; hi = [[None]*S for _ in range(S)]
    for r in range(S):
        for s in range(S):
            Vals = [[qq(samples[i][j][r][s]) for j in range(5)] for i in range(5)]
            # beta = Minv * Vals * Minv^T
            tmp = [[sum(Minv[i][k]*Vals[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
            beta = [sum(tmp[i][k]*Minv[j][k] for k in range(5)) for i in range(5) for j in range(5)]
            lo[r][s] = max(Q(0), min(beta)); hi[r][s] = min(Q(1), max(beta))
    return lo, hi

def propagate(Klo, Khi, start):
    S = len(Klo); vlo = [Q(0)]*S; vhi = [Q(0)]*S; vlo[start] = vhi[start] = Q(1)
    KloT = list(zip(*Klo)); KhiT = list(zip(*Khi))
    for _ in range(T):
        nlo = [sum(vlo[r]*KloT[s][r] for r in range(S)) for s in range(S)]
        nhi = [min(Q(1), sum(vhi[r]*KhiT[s][r] for r in range(S))) for s in range(S)]
        tot_lo = sum(nlo)
        nhi = [min(nhi[s], 1 - (tot_lo - nlo[s])) for s in range(S)]
        vlo, vhi = nlo, nhi
    return 1 - vhi[0], 1 - vlo[0]   # persistence [lo, hi]

def delta_interval(rect, e0, atoms, weights):
    KloR, KhiR = kernel_bounds("R", rect, e0, atoms, weights); KloC, KhiC = kernel_bounds("C", rect, e0, atoms, weights)
    PRlo, PRhi = propagate(KloR, KhiR, 2*n - 1); PClo, PChi = propagate(KloC, KhiC, n)
    return PRlo - PChi, PRhi - PClo

def certify(rect0, q, e0, maxdepth=12):
    (tl, wl), (tu, wu) = N3.laws(p, q)
    stack = [(rect0, 0)]; certified = []; unresolved = []; t0 = time.time(); evals = 0
    while stack:
        rect, depth = stack.pop(); evals += 1
        glo = delta_interval(rect, e0, tl, wl); ghi = delta_interval(rect, e0, tu, wu)
        if glo[1] < 0 and ghi[0] > 0:
            certified.append((rect, glo, ghi)); continue
        if glo[0] > 0 or ghi[1] < 0:   # certified NON-reversal somewhere? only if whole interval has wrong sign
            unresolved.append((rect, "wrong-sign-whole-rect", glo, ghi)); continue
        if depth >= maxdepth:
            unresolved.append((rect, "depth", glo, ghi)); continue
        a0, a1, c0, c1 = rect
        if (a1 - a0) / (rect0[1] - rect0[0]) >= (c1 - c0) / (rect0[3] - rect0[2]):
            am = (a0 + a1)/2; stack += [((a0, am, c0, c1), depth+1), ((am, a1, c0, c1), depth+1)]
        else:
            cm = (c0 + c1)/2; stack += [((a0, a1, c0, cm), depth+1), ((a0, a1, cm, c1), depth+1)]
    area0 = (rect0[1]-rect0[0])*(rect0[3]-rect0[2])
    area_c = sum((r[1]-r[0])*(r[3]-r[2]) for r, _, _ in certified)
    return dict(rect0=[str(v) for v in rect0], q=str(q), e0=str(e0), n=n, T=T, maxdepth=maxdepth, evaluations=evals,
                certified_leaves=len(certified), unresolved_leaves=len(unresolved), certified_area_fraction=str(F(area_c/area0)),
                certified_leaf_evidence=[{"rect": [str(v) for v in r], "delta_glo_interval": [str(g) for g in glo], "delta_ghi_interval": [str(g) for g in ghi]} for r, glo, ghi in certified],
                certified_area_fraction_float=float(area_c/area0), unresolved=[(tuple(str(v) for v in r), why, [float(g) for g in glo], [float(g) for g in ghi]) for r, why, glo, ghi in unresolved[:10]],
                seconds=round(time.time()-t0, 1))

if __name__ == "__main__":
    q, e0 = F(3,100), F(1,50)
    import os
    rect0 = (F(7,10), F(17,20), F(1,10), F(1,5))   # e1 in [0.70, 0.85], c in [0.10, 0.20]
    if os.environ.get("RECT") == "small":
        rect0 = (F(3,4), F(17,20), F(3,25), F(9,50))   # e1 in [0.75, 0.85], c in [0.12, 0.18]
    res = certify(rect0, q, e0, maxdepth=int(sys.argv[1]) if len(sys.argv) > 1 else 10)
    print(json.dumps(res, indent=1))
    out = "n4_rectangle_certificate_small.json" if os.environ.get("RECT") == "small" else "n4_rectangle_certificate.json"
    with open(out, "w") as f: json.dump(res, f, indent=1)
