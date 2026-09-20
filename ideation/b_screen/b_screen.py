"""DEEP candidate B: exact screening of action ranking under an unidentified regional-hazard mixing law.

Spec: correspondence thread 48, message 338 (Astra), accepted by Fable.
Model: n patches; state = occupancy before the annual hazard; all occupied at t=0.
Each year theta ~ G iid. Given theta, occupied patch i dies w.p. d_i(theta) = e0 + (e1_i - e0) theta,
independently. m = number of survivors. Each empty patch (incl. just vacated) is recolonised
independently w.p. r_m = 1 - (1-c)^m; only survivors supply colonists; r_0 = 0 so m = 0 is absorbing.
Actions: 0 baseline; R buffers one designated patch, e1_R = e0 + 0.5 (e1 - e0); C sets c_C = 2c.
Output: P_a(G) = P(no regional extinction through T = 50 steps).
"""
import itertools, sys, json
import numpy as np
from scipy.stats import norm, multivariate_normal
from math import comb

def bpmf(k, m, pr):
    return comb(m, k) * pr ** k * (1 - pr) ** (m - k)

from scipy.optimize import brentq

T = 50
BUFFER = 0.5
KAPPA = 2.0

# ---------------- mixing laws ----------------
def G_lower(p, q):
    return np.array([0.0, q / p]), np.array([1 - p * p / q, p * p / q])

def G_upper(p, q):
    w = (q - p * p) / (1 - 2 * p + q)
    a = (p - q) / (1 - p)
    return np.array([a, 1.0]), np.array([1 - w, w])

def G_mix(lam, p, q):
    tl, wl = G_lower(p, q); tu, wu = G_upper(p, q)
    return np.concatenate([tl, tu]), np.concatenate([(1 - lam) * wl, lam * wu])

def latent_rho(p, q):
    z = norm.ppf(p)
    f = lambda rho: multivariate_normal(mean=[0, 0], cov=[[1, rho], [rho, 1]]).cdf([z, z]) - q
    return brentq(f, 1e-9, 1 - 1e-9, xtol=1e-14)

def G_gauss(p, q, npts=80):
    """Moment-matched Gaussian-threshold mixing law, Gauss-Hermite (probabilists') quadrature on Z.
    Kernel entries are polynomials of degree <= n in theta, and theta(Z) is smooth, so this converges fast."""
    rho = latent_rho(p, q)
    z0 = norm.ppf(p)
    Z, w = np.polynomial.hermite_e.hermegauss(npts)
    w = w / w.sum()
    theta = norm.cdf((z0 - np.sqrt(rho) * Z) / np.sqrt(1 - rho))
    return theta, w, rho

def moments(theta, w):
    return float(np.dot(w, theta)), float(np.dot(w, theta ** 2))

# ---------------- kernels ----------------
def kernel_hom(n, d, c):
    """Occupied-count chain, homogeneous patches, given theta (death prob d)."""
    K = np.zeros((n + 1, n + 1))
    for k in range(n + 1):
        for m in range(k + 1):
            ps = bpmf(m, k, 1 - d)
            r = 1 - (1 - c) ** m
            for j in range(n - m + 1):
                K[k, m + j] += ps * bpmf(j, n - m, r)
    return K

def kernel_R(n, d, dR, c):
    """(b, j) chain: b = buffered patch occupied, j = ordinary occupied (0..n-1). index = b*n + j."""
    S = 2 * n
    K = np.zeros((S, S))
    for b in (0, 1):
        for j in range(n):
            i = b * n + j
            for s in range(j + 1):
                ps = bpmf(s, j, 1 - d)
                for bp in ((0,) if b == 0 else (0, 1)):
                    pb = 1.0 if b == 0 else (dR if bp == 0 else 1 - dR)
                    m = s + bp
                    r = 1 - (1 - c) ** m
                    for jj in range(n - 1 - s + 1):
                        pj = bpmf(jj, n - 1 - s, r)
                        if bp == 1:
                            K[i, n + s + jj] += ps * pb * pj
                        else:
                            K[i, 0 * n + s + jj] += ps * pb * pj * (1 - r)
                            K[i, 1 * n + s + jj] += ps * pb * pj * r
    return K

def kernel_labelled(n, dvec, c):
    """Full 2^n labelled-patch chain for validation. dvec[i] = death prob of patch i given theta."""
    S = 1 << n
    K = np.zeros((S, S))
    for state in range(S):
        occ = [i for i in range(n) if state >> i & 1]
        for surv_mask in range(1 << len(occ)):
            ps = 1.0; survivors = []
            for idx, i in enumerate(occ):
                if surv_mask >> idx & 1:
                    ps *= 1 - dvec[i]; survivors.append(i)
                else:
                    ps *= dvec[i]
            m = len(survivors)
            r = 1 - (1 - c) ** m
            empties = [i for i in range(n) if i not in survivors]
            for col_mask in range(1 << len(empties)):
                pc = 1.0; new = 0
                for i in survivors:
                    new |= 1 << i
                for idx, i in enumerate(empties):
                    if col_mask >> idx & 1:
                        pc *= r; new |= 1 << i
                    else:
                        pc *= 1 - r
                K[state, new] += ps * pc
    return K

def avg_kernel(build, theta, w):
    K = None
    for th, wt in zip(theta, w):
        Ki = build(th)
        K = wt * Ki if K is None else K + wt * Ki
    return K

def persist(K, start, absorb, T=T):
    v = np.zeros(K.shape[0]); v[start] = 1.0
    for _ in range(T):
        v = v @ K
    return 1.0 - v[absorb]

def d_of(theta, e0, e1):
    return e0 + (e1 - e0) * theta

def run_case(n, e0, e1, c, theta, w, T=T):
    dR_e1 = e0 + BUFFER * (e1 - e0)
    K0 = avg_kernel(lambda th: kernel_hom(n, d_of(th, e0, e1), c), theta, w)
    KC = avg_kernel(lambda th: kernel_hom(n, d_of(th, e0, e1), KAPPA * c), theta, w)
    KR = avg_kernel(lambda th: kernel_R(n, d_of(th, e0, e1), d_of(th, e0, dR_e1), c), theta, w)
    for K in (K0, KC, KR):
        assert np.allclose(K.sum(1), 1.0, atol=1e-12), "non-stochastic row"
    assert abs(K0[0, 0] - 1) < 1e-10 and abs(KC[0, 0] - 1) < 1e-10 and abs(KR[0, 0] - 1) < 1e-10, "absorption"
    P0 = persist(K0, n, 0, T)
    PC = persist(KC, n, 0, T)
    PR = persist(KR, 1 * n + (n - 1), 0, T)
    return P0, PR, PC

def labelled_case(n, e0, e1, c, theta, w, T=T):
    dR_e1 = e0 + BUFFER * (e1 - e0)
    full = (1 << n) - 1
    K0 = avg_kernel(lambda th: kernel_labelled(n, [d_of(th, e0, e1)] * n, c), theta, w)
    KC = avg_kernel(lambda th: kernel_labelled(n, [d_of(th, e0, e1)] * n, KAPPA * c), theta, w)
    KR = avg_kernel(lambda th: kernel_labelled(n, [d_of(th, e0, dR_e1)] + [d_of(th, e0, e1)] * (n - 1), c), theta, w)
    return persist(K0, full, 0, T), persist(KR, full, 0, T), persist(KC, full, 0, T)

# ---------------- validation ----------------
def laws(p, q, npts=80):
    out = [(f"lam={lam}", *G_mix(lam, p, q)) for lam in (0, 0.25, 0.5, 0.75, 1)]
    th, w, rho = G_gauss(p, q, npts)
    out.append(("gauss", th, w))
    return out, rho

def validate():
    p = 0.10
    print("== validation ==")
    for q in (0.03, 0.06):
        L, rho = laws(p, q)
        for name, th, w in L:
            m1, m2 = moments(th, w)
            assert abs(m1 - p) < 1e-9 and abs(m2 - q) < 1e-9, (name, m1, m2)
        gm1, gm2 = moments(L[-1][1], L[-1][2])
        print(f"moments matched for all six laws at q={q}; latent rho={rho:.10f}; gauss moment errors {gm1-p:.1e}, {gm2-q:.1e}")
        # Gaussian quadrature convergence
        thA, wA, _ = G_gauss(p, q, 80); thB, wB, _ = G_gauss(p, q, 160)
        a = run_case(5, 0.10, 0.90, 0.10, thA, wA); b = run_case(5, 0.10, 0.90, 0.10, thB, wB)
        print(f"  gauss convergence (n=5 case): max |dP| = {max(abs(x-y) for x,y in zip(a,b)):.2e}")
        # n=2 negative control: all laws identical
        vals = np.array([run_case(2, 0.10, 0.90, 0.10, th, w) for _, th, w in L])
        print(f"  n=2 control: spread across laws = {np.ptp(vals, axis=0).max():.2e} (must be ~0)")
        # e1 = e0: laws irrelevant
        vals = np.array([run_case(5, 0.10, 0.10, 0.10, th, w) for _, th, w in L])
        print(f"  e1=e0 control: spread across laws = {np.ptp(vals, axis=0).max():.2e} (must be ~0)")
        # c = 1 analytic check
        for name, th, w in L[:2] + L[-1:]:
            n, e0, e1 = 5, 0.10, 0.90
            P0, PR, PC = run_case(n, e0, e1, 1.0, th, w)
            d = d_of(th, e0, e1); dR = d_of(th, e0, e0 + BUFFER * (e1 - e0))
            a0 = (1 - np.dot(w, d ** n)) ** T; aR = (1 - np.dot(w, d ** (n - 1) * dR)) ** T
            print(f"  c=1 analytic ({name}): |P0-analytic|={abs(P0-a0):.2e} |PR-analytic|={abs(PR-aR):.2e} |PC-P0|={abs(PC-P0):.2e}")
        # reduced vs labelled kernels
        for n, e0, e1, c in ((3, 0.10, 0.90, 0.10), (3, 0.02, 0.50, 0.02), (5, 0.10, 0.90, 0.10)):
            for name, th, w in (L[0], L[4]):
                r = run_case(n, e0, e1, c, th, w); l = labelled_case(n, e0, e1, c, th, w)
                print(f"  reduced vs labelled n={n} e0={e0} e1={e1} c={c} {name}: max diff {max(abs(x-y) for x,y in zip(r,l)):.2e}")
        # T=1: C equals baseline; R >= baseline
        for name, th, w in (L[0], L[4]):
            P0, PR, PC = run_case(5, 0.10, 0.90, 0.10, th, w, T=1)
            print(f"  T=1 ({name}): PC-P0={PC-P0:.2e} (must be 0), PR-P0={PR-P0:.3e} (>=0)")

# ---------------- main grid ----------------
def main():
    validate()
    p = 0.10
    rows = []
    print("\n== screen ==")
    for n, q, e0, e1, c in itertools.product((5, 10), (0.03, 0.06), (0.02, 0.10), (0.50, 0.90), (0.02, 0.10)):
        L, rho = laws(p, q)
        res = {}
        for name, th, w in L:
            P0, PR, PC = run_case(n, e0, e1, c, th, w)
            assert PR >= P0 - 1e-12 and PC >= P0 - 1e-12, "action worse than baseline"
            res[name] = (P0, PR, PC)
        gP0, gPR, gPC = res["gauss"]
        if abs(gPR - gPC) < 1e-12:
            sel = "tie"
        else:
            sel = "R" if gPR > gPC else "C"
        for name, (P0, PR, PC) in res.items():
            best = max(PR, PC)
            Psel = {"R": PR, "C": PC, "tie": best}[sel]
            regret = best - Psel
            bestgain = best - P0
            rows.append(dict(n=n, q=q, e0=e0, e1=e1, c=c, law=name, P0=P0, PR=PR, PC=PC,
                             dR=PR - P0, dC=PC - P0, R_minus_C=PR - PC, gauss_sel=sel,
                             regret=regret, rel_regret=(regret / bestgain if bestgain > 0 else float("nan"))))
    import csv
    with open("b_screen_results.csv", "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wri.writeheader(); wri.writerows(rows)
    # summary per case
    print(f"{'n':>2} {'q':>4} {'e0':>4} {'e1':>4} {'c':>4} | {'P0 lam0':>8} {'P0 lam1':>8} {'P0 gauss':>8} | {'R-C lam0':>9} {'R-C lam1':>9} {'R-C gauss':>9} | sel  max regret  max relreg  sign-flip")
    flips = 0
    for key in sorted(set((r['n'], r['q'], r['e0'], r['e1'], r['c']) for r in rows)):
        sub = {r['law']: r for r in rows if (r['n'], r['q'], r['e0'], r['e1'], r['c']) == key}
        signs = set(np.sign(round(r['R_minus_C'], 12)) for r in sub.values() if abs(r['R_minus_C']) > 1e-12)
        flip = len(signs) > 1
        flips += flip
        mr = max(r['regret'] for r in sub.values()); mrr = np.nanmax([r['rel_regret'] for r in sub.values()])
        print(f"{key[0]:>2} {key[1]:>4} {key[2]:>4} {key[3]:>4} {key[4]:>4} | {sub['lam=0']['P0']:8.4f} {sub['lam=1']['P0']:8.4f} {sub['gauss']['P0']:8.4f} | {sub['lam=0']['R_minus_C']:+9.5f} {sub['lam=1']['R_minus_C']:+9.5f} {sub['gauss']['R_minus_C']:+9.5f} | {sub['gauss']['gauss_sel']:>3}  {mr:10.5f}  {mrr:10.3f}  {'YES' if flip else 'no'}")
    print(f"\ncases with sign reversal of P_R - P_C across the six laws: {flips} / 32")

if __name__ == "__main__":
    main()
