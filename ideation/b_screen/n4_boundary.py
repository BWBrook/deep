"""n=4, second stage: exact univariate Bernstein certificates on the two boundary curves of the (m3, m4) feasible set.

Given (n4_vertical.py): dJ/dy < 0 for J in {P0, PR, PC} and dDelta/dy > 0 over the whole feasible set, for all
16 cases. Hence over the class: max P_a lies on the lower curve y = L(x), min P_a on the upper curve y = U(x);
max Delta on the upper curve, min Delta on the lower curve. Each curve gives a kernel quadratic in x, so the
T-step objective is a univariate polynomial of degree <= 2T = 100 in u, x = xL + u (xU - xL), u in [0, 1].

Outputs per case: exact Bernstein coefficients (degree 100) of P0, PR, PC, Delta on both curves; sign
certificates for Delta (uniform sign, or isolated root brackets via monotone pieces after exact subdivision);
certified enclosures for min/max of each objective on each curve to width <= 1e-9; the resulting whole-class
bounds. Exact rationals throughout (gmpy2.mpq). Cross-check: Bernstein evaluation at u = 1/4, 3/4 equals
direct exact propagation with the affine kernel at (x, y(x)).
"""
import sys, json, time, itertools
from fractions import Fraction as F
from math import comb
sys.path.insert(0, __import__("os").path.dirname(__file__))
import n4_vertical as V
from gmpy2 import mpq as Q

T = V.T; P = V.P
def q_(f): return Q(f.numerator, f.denominator)

def bern_mul_deg2(v, K0, K1, K2):
    """v: list (len m+1) of state row-vectors (Bernstein degree m). Returns degree m+2 coefficients of v(u) K(u),
    K(u) = B02 K0 + B12 K1 + B22 K2. Uses B_{i,m} B_{k,2} = C(m,i)C(2,k)/C(m+2,i+k) B_{i+k,m+2}."""
    m = len(v) - 1; S = len(K0)
    Ks = (K0, K1, K2); Kt = [list(zip(*K)) for K in Ks]
    out = [[Q(0)] * S for _ in range(m + 3)]
    for i in range(m + 1):
        for k in range(3):
            coef = Q(comb(m, i) * comb(2, k), comb(m + 2, i + k))
            row = v[i]; col = Kt[k]
            tgt = out[i + k]
            for s in range(S):
                tgt[s] += coef * sum(row[r] * col[s][r] for r in range(S))
    return out

def curve_kernels(base, Bx, Cy, xL, xU, yfun, p, q):
    """Bernstein-degree-2 control matrices of K(u) = base + x(u) B + y(x(u)) C on u in [0,1]."""
    def Kat(u):
        x = xL + u * (xU - xL); y = yfun(x, p, q)
        return [[q_(base[r][s] + x * Bx[r][s] + y * Cy[r][s]) for s in range(len(base))] for r in range(len(base))]
    K0, Kh, K2 = Kat(F(0)), Kat(F(1, 2)), Kat(F(1))
    K1 = [[2 * Kh[r][s] - (K0[r][s] + K2[r][s]) / 2 for s in range(len(K0))] for r in range(len(K0))]
    return K0, K1, K2

def objective_coeffs(base, Bx, Cy, start, xL, xU, yfun, p, q):
    K0, K1, K2 = curve_kernels(base, Bx, Cy, xL, xU, yfun, p, q)
    S = len(K0)
    v = [[Q(0)] * S]; v[0][start] = Q(1)
    for _ in range(T):
        v = bern_mul_deg2(v, K0, K1, K2)
    return [sum(row) - row[0] for row in v]   # degree 2T Bernstein coefficients of persistence

def evaluate(coeffs, u):
    n = len(coeffs) - 1; u = Q(u.numerator, u.denominator) if isinstance(u, F) else Q(u)
    return sum(c * comb(n, k) * u ** k * (1 - u) ** (n - k) for k, c in enumerate(coeffs))

def split(coeffs):
    n = len(coeffs) - 1; tri = [list(coeffs)]
    for k in range(1, n + 1):
        prev = tri[-1]; tri.append([(prev[m] + prev[m + 1]) / 2 for m in range(n - k + 1)])
    return [tri[k][0] for k in range(n + 1)], [tri[n - k][k] for k in range(n + 1)]

def enclose(coeffs, tol=F(1, 10**9), maxd=60):
    """Certified [min, max] of the Bernstein polynomial on [0,1] by branch and bound.
    A sub-interval is discarded for the max search when its coefficient maximum is below the best exact
    value found (endpoint values are exact evaluations); symmetrically for the min. Returns (lo, hi) such that
    lo <= min f <= lo + tol and hi - tol <= max f <= hi."""
    tolq = q_(tol)
    def bound(kind):
        best = max(coeffs[0], coeffs[-1]) if kind == "max" else min(coeffs[0], coeffs[-1])  # exact attained values
        stack = [(coeffs, 0)]; outer = best
        while stack:
            cf, depth = stack.pop()
            cmax, cmin = max(cf), min(cf)
            if kind == "max":
                if cmax <= best + tolq: outer = max(outer, cmax); continue
                best = max(best, cf[0], cf[-1])
            else:
                if cmin >= best - tolq: outer = min(outer, cmin); continue
                best = min(best, cf[0], cf[-1])
            if depth >= maxd:
                outer = max(outer, cmax) if kind == "max" else min(outer, cmin); continue
            a, b = split(cf); stack += [(a, depth + 1), (b, depth + 1)]
        return best, outer
    lo_att, lo_out = bound("min"); hi_att, hi_out = bound("max")
    return lo_out, hi_out

def sign_certificate(coeffs, u0=F(0), u1=F(1), depth=0, maxd=60):
    """Return list of (u0, u1, status) pieces: status in {'pos','neg','root'}; 'root' pieces are certified
    to contain exactly one simple root (polynomial strictly monotone with opposite endpoint signs) and are
    narrowed to width <= 2^-40."""
    if all(c > 0 for c in coeffs): return [(u0, u1, "pos")]
    if all(c < 0 for c in coeffs): return [(u0, u1, "neg")]
    n = len(coeffs) - 1
    d = [n * (coeffs[k + 1] - coeffs[k]) for k in range(n)]
    f0, f1 = coeffs[0], coeffs[-1]
    if (all(x > 0 for x in d) or all(x < 0 for x in d)) and f0 * f1 < 0:
        # monotone with a sign change: bisect exactly to isolate the root
        lo, hi, cl, ch = u0, u1, coeffs, None
        while hi - lo > F(1, 2**40):
            a, b = split(cl); mid = (lo + hi) / 2
            if a[0] * a[-1] <= 0: cl, hi = a, mid
            else: cl, lo = b, mid
        return [(lo, hi, "root")]
    if depth >= maxd: return [(u0, u1, "unresolved")]
    a, b = split(coeffs); um = (u0 + u1) / 2
    return sign_certificate(a, u0, um, depth + 1, maxd) + sign_certificate(b, um, u1, depth + 1, maxd)

def main():
    n = 4; t0 = time.time()
    out = {"T": T, "p": str(P), "cases": []}
    for q in V.QS:
        xL, xU = V.xbounds(P, q)
        for e0, e1, c in itertools.product(V.E0S, V.E1S, V.CS):
            t1 = time.time()
            parts = {a: V.affine_parts(a, n, e0, e1, c, P, q) for a in ("0", "R", "C")}
            rec = {"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), "curves": {}}
            co = {}
            for cname, yfun in (("lower", V.Lfun), ("upper", V.Ufun)):
                cc = {}
                for a in ("0", "R", "C"):
                    base, Bx, Cy = parts[a]
                    cc[a] = objective_coeffs(base, Bx, Cy, V.start_state(a, n), xL, xU, yfun, P, q)
                    # cross-check against direct exact propagation at u = 1/4, 3/4
                    for u in (F(1, 4), F(3, 4)):
                        x = xL + u * (xU - xL); y = yfun(x, P, q)
                        direct = V.persist(V.kernel_at(base, Bx, Cy, x, y), V.start_state(a, n))
                        assert evaluate(cc[a], u) == q_(direct), ("cross-check failed", cname, a, u)
                cc["D"] = [r - s for r, s in zip(cc["R"], cc["C"])]
                co[cname] = cc
                crv = {}
                for name, key in (("P0", "0"), ("PR", "R"), ("PC", "C"), ("Delta", "D")):
                    lo, hi = enclose(cc[key])
                    crv[name] = {"degree": len(cc[key]) - 1, "min": float(lo), "max": float(hi), "at_u0": float(cc[key][0]), "at_u1": float(cc[key][-1])}
                pieces = sign_certificate(cc["D"])
                crv["Delta_sign"] = [(float(a), float(b), s) for a, b, s in pieces]
                crv["Delta_roots_x"] = [(float(xL + a * (xU - xL)), float(xL + b * (xU - xL))) for a, b, s in pieces if s == "root"]
                rec["curves"][cname] = crv
            # whole-class bounds from vertical monotonicity: P decreasing in y -> max on lower, min on upper; Delta increasing in y -> max upper, min lower
            wc = {}
            for name in ("P0", "PR", "PC"):
                wc[name] = {"min": rec["curves"]["upper"][name]["min"], "max": rec["curves"]["lower"][name]["max"]}
            wc["Delta"] = {"min": rec["curves"]["lower"]["Delta"]["min"], "max": rec["curves"]["upper"]["Delta"]["max"]}
            dmin, dmax = wc["Delta"]["min"], wc["Delta"]["max"]
            wc["Delta_sign_over_class"] = "R wins throughout" if dmin > 0 else "C wins throughout" if dmax < 0 else "UNIDENTIFIED"
            rec["whole_class"] = wc; rec["seconds"] = round(time.time() - t1, 1)
            out["cases"].append(rec)
            print(f"q={str(q):>5} e0={str(e0):>4} e1={str(e1):>4} c={str(c):>4} | P0 [{wc['P0']['min']:.4f},{wc['P0']['max']:.4f}] PR [{wc['PR']['min']:.4f},{wc['PR']['max']:.4f}] PC [{wc['PC']['min']:.4f},{wc['PC']['max']:.4f}] | Delta [{dmin:+.5f},{dmax:+.5f}] {wc['Delta_sign_over_class']:<18} | lower-curve sign {[s for _,_,s in rec['curves']['lower']['Delta_sign']]} upper-curve sign {[s for _,_,s in rec['curves']['upper']['Delta_sign']]} [{rec['seconds']}s]", flush=True)
    with open("n4_boundary_results.json", "w") as f:
        json.dump(out, f, indent=1)
    print(f"total {time.time()-t0:.1f}s; wrote n4_boundary_results.json")

if __name__ == "__main__":
    main()
