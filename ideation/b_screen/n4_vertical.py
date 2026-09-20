"""n=4: exact affine moment kernels, boundary laws of the (m3, m4) feasible set, independent witness check,
float scan of dJ/dy sign over the feasible lens, and an exact Bernstein branch-and-bound certificate of
sign(dJ/dy) over the feasible set for J in {P0, PR, PC, Delta}.

Geometry (Astra, n4_geometry.md): x=m3 in [xL, xU]; y=m4 in [L(x), U(x)].
Kernel construction reuses kernel_hom / kernel_buffered / d_of from n3_certificate.py (Astra) via import.
"""
import sys, json, time, itertools
from fractions import Fraction as F
from math import comb
sys.path.insert(0, __import__("os").path.dirname(__file__))
import n3_certificate as N3  # Astra's exact kernels: kernel_hom(n, d, c), kernel_buffered(n, d, dR, c), d_of(theta, e0, e1)

T = 50
P = F(1, 10)
QS = (F(3, 100), F(6, 100)); E0S = (F(1, 50), F(1, 10)); E1S = (F(1, 2), F(9, 10)); CS = (F(1, 50), F(1, 10))

# ---------- geometry ----------
def xbounds(p, q):
    return q * q / p, q - (p - q) ** 2 / (1 - p)
def Lfun(x, p, q): return (x * x - 2 * p * q * x + q ** 3) / (q - p * p)
def Ufun(x, p, q): return x - (q - x) ** 2 / (p - q)
def lower_law(x, p, q):
    # two-point law with moments p, q, x: atoms roots of z^2 - B z + C
    B = (x - p * q) / (q - p * p); C = (p * x - q * q) / (q - p * p)
    disc = B * B - 4 * C
    from math import isqrt
    # exact sqrt if rational square, else float fallback (only used for float checks)
    num, den = disc.numerator, disc.denominator
    rn, rd = isqrt(num), isqrt(den)
    if rn * rn == num and rd * rd == den:
        r = F(rn, rd); zm, zp = (B - r) / 2, (B + r) / 2
        return [zm, zp], [(zp - p) / (zp - zm), (p - zm) / (zp - zm)]
    r = disc ** 0.5; zm, zp = (float(B) - r) / 2, (float(B) + r) / 2
    return [zm, zp], [(zp - float(p)) / (zp - zm), (float(p) - zm) / (zp - zm)]
def upper_law(x, p, q):
    a = (q - x) / (p - q)
    wa = (p - q) / (a * (1 - a)); w0 = 1 - p - (p - q) / a; w1 = p - (p - q) / (1 - a)
    return [F(0), a, F(1)], [w0, wa, w1]

# ---------- exact affine kernels ----------
def raw_kernel(action, n, th, e0, e1, c):
    d = N3.d_of(th, e0, e1)
    if action == "R":
        return N3.kernel_buffered(n, d, N3.d_of(th, e0, e0 + F(1, 2) * (e1 - e0)), c)
    return N3.kernel_hom(n, d, c * (2 if action == "C" else 1))

def moment_kernels(action, n, e0, e1, c):
    """A_k, k=0..n, with K_G = sum_k m_k A_k. Exact via Vandermonde at nodes theta_j = j/n."""
    nodes = [F(j, n) for j in range(n + 1)]
    Ks = [raw_kernel(action, n, th, e0, e1, c) for th in nodes]
    S = len(Ks[0])
    # solve V a = k for each entry, V[j][k] = nodes[j]^k
    V = [[th ** k for k in range(n + 1)] for th in nodes]
    # invert V by Gaussian elimination once
    m = [row[:] + [F(int(i == j)) for j in range(n + 1)] for i, row in enumerate(V)]
    for col in range(n + 1):
        piv = next(r for r in range(col, n + 1) if m[r][col] != 0); m[col], m[piv] = m[piv], m[col]
        pv = m[col][col]; m[col] = [v / pv for v in m[col]]
        for r in range(n + 1):
            if r != col and m[r][col] != 0:
                f = m[r][col]; m[r] = [a - f * b for a, b in zip(m[r], m[col])]
    Vinv = [row[n + 1:] for row in m]
    A = [[[sum(Vinv[k][j] * Ks[j][r][s] for j in range(n + 1)) for s in range(S)] for r in range(S)] for k in range(n + 1)]
    return A

def affine_parts(action, n, e0, e1, c, p, q):
    A = moment_kernels(action, n, e0, e1, c)
    S = len(A[0])
    base = [[A[0][r][s] + p * A[1][r][s] + q * A[2][r][s] for s in range(S)] for r in range(S)]
    return base, A[3], A[4]

def avg_direct(action, n, e0, e1, c, atoms, weights):
    K = None
    for th, wt in zip(atoms, weights):
        Ki = raw_kernel(action, n, th, e0, e1, c)
        K = [[wt * v for v in row] for row in Ki] if K is None else [[a + wt * b for a, b in zip(ra, rb)] for ra, rb in zip(K, Ki)]
    return K

def kernel_at(base, Bx, Cy, x, y):
    return [[base[r][s] + x * Bx[r][s] + y * Cy[r][s] for s in range(len(base))] for r in range(len(base))]

def persist(K, start, T=T):
    v = [F(0)] * len(K); v[start] = F(1)
    for _ in range(T):
        v = [sum(v[r] * K[r][s] for r in range(len(K))) for s in range(len(K))]
    return 1 - v[0]

def start_state(action, n):
    return (2 * n - 1) if action == "R" else n   # buffered chain index b*n + j with b=1, j=n-1

# ---------- Bernstein machinery on a box, scalar-valued after summing live states ----------
try:
    from gmpy2 import mpq as Q  # exact rationals, much faster than fractions.Fraction
except Exception:  # pragma: no cover
    Q = F

def bernstein_J(base, Bx, Cy, start, xlo, xhi, ylo, yhi, T=T):
    """Tensor Bernstein coefficients (degree T in u, T in w) of J(u,w) = persistence, box-normalised coords.
    Exact rational arithmetic (gmpy2.mpq when available, else Fraction)."""
    S = len(base)
    K = {(a, b): [[Q(v.numerator, v.denominator) for v in row] for row in kernel_at(base, Bx, Cy, xlo if a == 0 else xhi, ylo if b == 0 else yhi)] for a in (0, 1) for b in (0, 1)}
    # transpose kernels for fast row-vector multiplication: Kt[(a,b)][s] = column s
    Kt = {k: list(zip(*v)) for k, v in K.items()}
    v = {(0, 0): [Q(0)] * S}; v[(0, 0)][start] = Q(1)
    for t in range(1, T + 1):
        nv = {}
        for i in range(t + 1):
            for j in range(t + 1):
                acc = [Q(0)] * S
                for a in (0, 1):
                    ca = Q(t - i, t) if a == 0 else Q(i, t)
                    if ca == 0: continue
                    for b in (0, 1):
                        cb = Q(t - j, t) if b == 0 else Q(j, t)
                        if cb == 0: continue
                        prev = v.get((i - a, j - b))
                        if prev is None: continue
                        col = Kt[(a, b)]; coef = ca * cb
                        for s in range(S):
                            acc[s] += coef * sum(prev[r] * col[s][r] for r in range(S))
                nv[(i, j)] = acc
        v = nv
    return [[sum(v[(i, j)]) - v[(i, j)][0] for j in range(T + 1)] for i in range(T + 1)]  # J[i][j]

def dw_coeffs(J):
    T_ = len(J) - 1
    return [[T_ * (J[i][j + 1] - J[i][j]) for j in range(T_)] for i in range(T_ + 1)]

def sign_of(arr):
    flat = [v for row in arr for v in row]
    if all(v > 0 for v in flat): return "pos"
    if all(v < 0 for v in flat): return "neg"
    return "mixed"

def de_casteljau_split(arr, axis):
    """Split a 2-D Bernstein coefficient array at the midpoint along axis (0: rows/u, 1: cols/w). Exact."""
    if axis == 1:
        arr = [list(r) for r in arr]
        left, right = [], []
        for row in arr:
            n = len(row) - 1; tri = [row[:]]
            for k in range(1, n + 1):
                prev = tri[-1]; tri.append([(prev[m] + prev[m + 1]) / 2 for m in range(n - k + 1)])
            left.append([tri[k][0] for k in range(n + 1)]); right.append([tri[n - k][k] for k in range(n + 1)])
        return left, right
    cols = list(map(list, zip(*arr)))
    l, r = de_casteljau_split(cols, 1)
    return list(map(list, zip(*l))), list(map(list, zip(*r)))

def lens_range(x0, x1, p, q):
    """min of L and max of U over [x0, x1], exact."""
    Lv = [Lfun(x0, p, q), Lfun(x1, p, q)]; xv = p * q
    if x0 <= xv <= x1: Lv.append(Lfun(xv, p, q))
    Uv = [Ufun(x0, p, q), Ufun(x1, p, q)]; xu = q + (p - q) / 2
    if x0 <= xu <= x1: Uv.append(Ufun(xu, p, q))
    return min(Lv), max(Uv)

def certify_dw(D, xlo, xhi, ylo, yhi, p, q, max_depth=10):
    """Branch and bound over the box in (u,w); prune boxes disjoint from the lens. Returns dict.
    Box coordinates and lens bounds stay in Fraction; coefficient arrays are Q."""
    stack = [(D, F(0), F(1), F(0), F(1), 0)]
    leaves = {"pos": 0, "neg": 0, "pruned": 0, "unresolved": 0}; unresolved_boxes = []
    while stack:
        arr, u0, u1, w0, w1, depth = stack.pop()
        x0, x1 = xlo + u0 * (xhi - xlo), xlo + u1 * (xhi - xlo)
        y0, y1 = ylo + w0 * (yhi - ylo), ylo + w1 * (yhi - ylo)
        lmin, umax = lens_range(x0, x1, p, q)
        if y1 < lmin or y0 > umax:
            leaves["pruned"] += 1; continue
        s = sign_of(arr)
        if s in ("pos", "neg"):
            leaves[s] += 1; continue
        if depth >= max_depth:
            leaves["unresolved"] += 1; unresolved_boxes.append((float(x0), float(x1), float(y0), float(y1))); continue
        axis = 0 if (u1 - u0) * (xhi - xlo) / (xhi) >= (w1 - w0) * (yhi - ylo) / max(yhi, F(1, 10**9)) else 1
        a, b = de_casteljau_split(arr, axis)
        if axis == 0:
            um = (u0 + u1) / 2; stack += [(a, u0, um, w0, w1, depth + 1), (b, um, u1, w0, w1, depth + 1)]
        else:
            wm = (w0 + w1) / 2; stack += [(a, u0, u1, w0, wm, depth + 1), (b, u0, u1, wm, w1, depth + 1)]
    return leaves, unresolved_boxes

# ---------- main ----------
def main(mode):
    n = 4
    out = {"T": T, "p": str(P), "stage0": [], "witness": {}, "float_scan": [], "certificates": []}
    t0 = time.time()
    # Stage 0: exact affine kernels and checks
    parts = {}
    for q in QS:
        xL, xU = xbounds(P, q)
        for e0, e1, c in itertools.product(E0S, E1S, CS):
            for action in ("0", "R", "C"):
                base, Bx, Cy = affine_parts(action, n, e0, e1, c, P, q)
                parts[(q, e0, e1, c, action)] = (base, Bx, Cy)
                # check against direct averaging for Glo, Ghi and an interior mixture law
                (tl, wl), (tu, wu) = N3.laws(P, q)
                for atoms, weights in ((tl, wl), (tu, wu)):
                    x = sum(w * a ** 3 for a, w in zip(atoms, weights)); y = sum(w * a ** 4 for a, w in zip(atoms, weights))
                    Kd = avg_direct(action, n, e0, e1, c, atoms, weights); Ka = kernel_at(base, Bx, Cy, x, y)
                    assert all(Kd[r][s] == Ka[r][s] for r in range(len(Kd)) for s in range(len(Kd))), "affine kernel mismatch"
            out["stage0"].append({"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), "affine_exact_match": True})
    print(f"Stage 0: exact affine kernels for 16 cases x 3 actions verified against direct averaging at Glo, Ghi ({time.time()-t0:.1f}s)")
    # Witness check (Astra #380): q=3/50, e0=1/50, e1=9/10, c=1/10; x=29/500, y=288411/5000000; x+-
    q, e0, e1, c = F(3, 50), F(1, 50), F(9, 10), F(1, 10)
    y = F(288411, 5000000)
    def delta_at(x, y):
        vals = {}
        for action in ("R", "C"):
            base, Bx, Cy = parts[(q, e0, e1, c, action)]
            vals[action] = persist(kernel_at(base, Bx, Cy, x, y), start_state(action, n))
        return vals["R"] - vals["C"]
    xm, xp = F(579999, 10**7), F(580001, 10**7)
    xL, xU = xbounds(P, q)
    feas = lambda x, y: xL <= x <= xU and Lfun(x, P, q) <= y <= Ufun(x, P, q)
    dm, dp = delta_at(xm, y), delta_at(xp, y)
    out["witness"] = {"feasible_minus": feas(xm, y), "feasible_plus": feas(xp, y), "Delta_minus": float(dm), "Delta_plus": float(dp), "Delta_plus_lt_minus": dp < dm,
                      "dDelta_dx_finite_diff": float((dp - dm) / (xp - xm))}
    print("Witness:", out["witness"])
    if mode == "stage0": return out
    # Float scan of dJ/dy over the lens (25 x 25 grid in (t,s)) for all 16 cases and 4 objectives
    import numpy as np
    for q in QS:
        xL, xU = xbounds(P, q)
        for e0, e1, c in itertools.product(E0S, E1S, CS):
            bp = {a: tuple(np.array([[float(v) for v in row] for row in M]) for M in parts[(q, e0, e1, c, a)]) for a in ("0", "R", "C")}
            def Jf(a, x, y):
                base, Bx, Cy = bp[a]; K = base + x * Bx + y * Cy
                v = np.zeros(len(K)); v[start_state(a, n)] = 1.0
                for _ in range(T): v = v @ K
                return 1 - v[0]
            signs = {k: set() for k in ("P0", "PR", "PC", "Delta")}; minabs = {k: 1.0 for k in signs}
            for t in np.linspace(0.02, 0.98, 25):
                x = float(xL) + t * float(xU - xL); L_, U_ = float(Lfun(F(x), P, q)), float(Ufun(F(x), P, q))
                if U_ - L_ <= 0: continue
                h = (U_ - L_) * 1e-4
                for s in np.linspace(0.05, 0.95, 25):
                    yv = L_ + s * (U_ - L_)
                    d = {a: (Jf(a, x, yv + h) - Jf(a, x, yv - h)) / (2 * h) for a in ("0", "R", "C")}
                    for k, val in (("P0", d["0"]), ("PR", d["R"]), ("PC", d["C"]), ("Delta", d["R"] - d["C"])):
                        signs[k].add(1 if val > 0 else -1 if val < 0 else 0); minabs[k] = min(minabs[k], abs(val))
            rec = {"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), **{k: ("uniform " + ("pos" if 1 in v else "neg") if len(v) == 1 else "MIXED") for k, v in signs.items()},
                   "min_abs": {k: float(v) for k, v in minabs.items()}}
            out["float_scan"].append(rec)
            print(f"scan q={str(q):>5} e0={str(e0):>4} e1={str(e1):>4} c={str(c):>4}: " + "  ".join(f"{k}:{rec[k]}" for k in ("P0", "PR", "PC", "Delta")))
    print(f"float scan done ({time.time()-t0:.1f}s)")
    if mode == "scan": return out
    # Exact certificate for cases listed in argv (default: all), objective Delta and P0/PR/PC
    for q in QS:
        xL, xU = xbounds(P, q); ylo, yhi = lens_range(xL, xU, P, q)
        for e0, e1, c in itertools.product(E0S, E1S, CS):
            t1 = time.time()
            Js = {}
            for action in ("0", "R", "C"):
                base, Bx, Cy = parts[(q, e0, e1, c, action)]
                Js[action] = bernstein_J(base, Bx, Cy, start_state(action, n), xL, xU, ylo, yhi)
            res = {"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), "box": [float(xL), float(xU), float(ylo), float(yhi)], "objectives": {}}
            for name, J in (("P0", Js["0"]), ("PR", Js["R"]), ("PC", Js["C"]), ("Delta", [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(Js["R"], Js["C"])])):
                D = dw_coeffs(J)
                leaves, unres = certify_dw(D, xL, xU, ylo, yhi, P, q, max_depth=MAXD)
                status = ("certified " + ("pos" if leaves["pos"] and not leaves["neg"] else "neg" if leaves["neg"] and not leaves["pos"] else "MIXED-SIGN-LEAVES")) if leaves["unresolved"] == 0 else "unresolved"
                res["objectives"][name] = {"status": status, "leaves": leaves, "unresolved_boxes": unres[:5]}
            res["seconds"] = round(time.time() - t1, 1)
            out["certificates"].append(res)
            print(f"cert q={str(q):>5} e0={str(e0):>4} e1={str(e1):>4} c={str(c):>4}: " + "  ".join(f"{k}:{v['status']}({v['leaves']['pos']}+,{v['leaves']['neg']}-,{v['leaves']['pruned']}x,{v['leaves']['unresolved']}?)" for k, v in res["objectives"].items()) + f"  [{res['seconds']}s]", flush=True)
    return out

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    MAXD = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    out = main(mode)
    with open(f"n4_vertical_{mode}.json", "w") as f:
        json.dump(out, f, indent=1)
    print("wrote", f"n4_vertical_{mode}.json")
