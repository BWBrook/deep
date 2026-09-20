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
from scipy.special import roots_hermitenorm
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

def G_gauss(p, q, npts=1280):
    """Moment-matched Gaussian-threshold mixing law, Gauss-Hermite (probabilists') quadrature on Z.
    Kernel entries are polynomials of degree <= n in theta, and theta(Z) is smooth, so this converges fast."""
    rho = latent_rho(p, q)
    z0 = norm.ppf(p)
    Z, w = roots_hermitenorm(npts)  # stable for large npts, unlike numpy's hermegauss
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
    # Sum transient mass directly.  Subtracting the absorbing mass from one
    # loses the signal when persistence is close to machine precision.
    return v[:absorb].sum() + v[absorb + 1:].sum()

def d_of(theta, e0, e1):
    return e0 + (e1 - e0) * theta

def check_prob(x, name):
    if not (0.0 <= x <= 1.0):
        raise ValueError(f"{name}={x} is not a probability")

def check_kernel(K):
    assert np.all(np.isfinite(K)), "non-finite kernel entry"
    assert np.all(K >= -1e-15), "negative kernel entry"
    assert np.allclose(K.sum(1), 1.0, rtol=0, atol=1e-12), "non-stochastic row"
    assert abs(K[0, 0] - 1) < 1e-10, "absorption at zero survivors"

def run_case(n, e0, e1, c, theta, w, T=T, with_C=True):
    for x, nm in ((e0, "e0"), (e1, "e1"), (c, "c")) + (((KAPPA * c, "kappa*c"),) if with_C else ()):
        check_prob(x, nm)
    assert np.all((theta >= 0) & (theta <= 1)) and abs(w.sum() - 1) < 1e-12, "invalid mixing law"
    dR_e1 = e0 + BUFFER * (e1 - e0)
    K0 = avg_kernel(lambda th: kernel_hom(n, d_of(th, e0, e1), c), theta, w)
    KR = avg_kernel(lambda th: kernel_R(n, d_of(th, e0, e1), d_of(th, e0, dR_e1), c), theta, w)
    check_kernel(K0); check_kernel(KR)
    P0 = persist(K0, n, 0, T)
    PR = persist(KR, 1 * n + (n - 1), 0, T)
    if not with_C:
        return P0, PR, None
    KC = avg_kernel(lambda th: kernel_hom(n, d_of(th, e0, e1), KAPPA * c), theta, w)
    check_kernel(KC)
    PC = persist(KC, n, 0, T)
    return P0, PR, PC

def labelled_case(n, e0, e1, c, theta, w, T=T):
    dR_e1 = e0 + BUFFER * (e1 - e0)
    full = (1 << n) - 1
    K0 = avg_kernel(lambda th: kernel_labelled(n, [d_of(th, e0, e1)] * n, c), theta, w)
    KC = avg_kernel(lambda th: kernel_labelled(n, [d_of(th, e0, e1)] * n, KAPPA * c), theta, w)
    KR = avg_kernel(lambda th: kernel_labelled(n, [d_of(th, e0, dR_e1)] + [d_of(th, e0, e1)] * (n - 1), c), theta, w)
    return persist(K0, full, 0, T), persist(KR, full, 0, T), persist(KC, full, 0, T)

# ---------------- validation ----------------
def laws(p, q, npts=1280):
    out = [(f"lam={lam}", *G_mix(lam, p, q)) for lam in (0, 0.25, 0.5, 0.75, 1)]
    th, w, rho = G_gauss(p, q, npts)
    out.append(("gauss", th, w))
    return out, rho

def validate():
    p = 0.10
    print("== validation (every control is an assertion with a declared tolerance) ==")
    for q in (0.03, 0.06):
        L, rho = laws(p, q)
        for name, th, w in L:
            m1, m2 = moments(th, w)
            assert abs(m1 - p) < 1e-9 and abs(m2 - q) < 1e-9, (name, m1, m2)
        gm1, gm2 = moments(L[-1][1], L[-1][2])
        print(f"q={q}: moments matched (tol 1e-9); latent rho={rho:.10f}; gauss moment errors {gm1-p:.1e}, {gm2-q:.1e}")
        # independent check of the Gaussian law's first 10 moments (all that enter any kernel here) by trapezoid on Z
        Zt = np.linspace(-12, 12, 480001); wt = norm.pdf(Zt); wt[0] *= .5; wt[-1] *= .5; wt /= wt.sum()
        tht = norm.cdf((norm.ppf(p) - np.sqrt(rho) * Zt) / np.sqrt(1 - rho))
        mq = np.array([np.dot(L[-1][2], L[-1][1] ** k) for k in range(1, 11)])
        mt = np.array([np.dot(wt, tht ** k) for k in range(1, 11)])
        assert np.max(np.abs(mq - mt)) < 1e-12, np.max(np.abs(mq - mt))
        print(f"  gauss law moments 1..10, quadrature vs 480001-point trapezoid: max diff {np.max(np.abs(mq - mt)):.2e} (tol 1e-12)")
        vals = np.array([run_case(2, 0.10, 0.90, 0.10, th, w) for _, th, w in L])
        sp = np.ptp(vals, axis=0).max(); assert sp < 1e-9, sp
        print(f"  n=2 control: spread across laws {sp:.2e} (tol 1e-9)")
        vals = np.array([run_case(5, 0.10, 0.10, 0.10, th, w) for _, th, w in L])
        sp = np.ptp(vals, axis=0).max(); assert sp < 1e-12, sp
        print(f"  e1=e0 control: spread across laws {sp:.2e} (tol 1e-12)")
        for name, th, w in L:
            for n in (5, 10):
                e0, e1 = 0.10, 0.90
                d = d_of(th, e0, e1); dR = d_of(th, e0, e0 + BUFFER * (e1 - e0))
                P0, PR, _ = run_case(n, e0, e1, 1.0, th, w, with_C=False)
                a0 = (1 - np.dot(w, d ** n)) ** T; aR = (1 - np.dot(w, d ** (n - 1) * dR)) ** T
                assert abs(P0 - a0) < 1e-12 and abs(PR - aR) < 1e-12, (name, n, P0 - a0, PR - aR)
                _, _, PC = run_case(n, e0, e1, 0.5, th, w)
                assert abs(PC - a0) < 1e-12, (name, n, PC - a0)
        print(f"  c=1 analytic controls (baseline and R at c=1; C at base c=0.5 so kappa*c=1): within 1e-12, six laws, n=5,10")
        for n, e0, e1, c in ((3, 0.10, 0.90, 0.10), (3, 0.02, 0.50, 0.02), (5, 0.10, 0.90, 0.10)):
            for name, th, w in (L[0], L[4], L[5]):
                r = run_case(n, e0, e1, c, th, w); l = labelled_case(n, e0, e1, c, th, w)
                md = max(abs(x - y) for x, y in zip(r, l)); assert md < 1e-12, md
        print(f"  reduced vs labelled persistence (n=3 x2, n=5; lower, upper, gauss laws): within 1e-12")
        for n in (3, 5):
            for th_ in (0.0, 0.37, 1.0):
                for e0, e1, c in ((0.02, 0.50, 0.02), (0.10, 0.90, 0.10)):
                    d = d_of(th_, e0, e1)
                    Kl = kernel_labelled(n, [d] * n, c); Kh = kernel_hom(n, d, c)
                    cnt = np.array([bin(s_).count("1") for s_ in range(1 << n)])
                    agg = np.zeros((n + 1, n + 1))
                    for s_ in range(1 << n):
                        for t_ in range(1 << n):
                            agg[cnt[s_], cnt[t_]] += Kl[s_, t_]
                    mult = np.array([comb(n, k) for k in range(n + 1)])
                    assert np.max(np.abs(agg / mult[:, None] - Kh)) < 1e-13
        print(f"  reduced vs labelled whole-kernel (hom chain; n=3,5; theta 0, .37, 1; two parameter sets): within 1e-13")
        for name, th, w in L:
            P0, PR, PC = run_case(5, 0.10, 0.90, 0.10, th, w, T=1)
            assert PC == P0 and PR >= P0, (name, PC - P0, PR - P0)
        print(f"  T=1: PC == P0 exactly and PR >= P0, six laws")
    worst = 0.0; where = None
    for n, q, e0, e1, c in itertools.product((5, 10), (0.03, 0.06), (0.02, 0.10), (0.50, 0.90), (0.02, 0.10)):
        thA, wA, _ = G_gauss(p, q, 640); thB, wB, _ = G_gauss(p, q, 1280)
        a_ = run_case(n, e0, e1, c, thA, wA); b_ = run_case(n, e0, e1, c, thB, wB)
        dmax = max(abs(x - y) for x, y in zip(a_, b_))
        if dmax > worst: worst, where = dmax, (n, q, e0, e1, c)
    assert worst < 1e-9, (worst, where)
    print(f"  gauss convergence over all 32 cases, 640 vs 1280 nodes: max |dP| {worst:.2e} at {where} (tol 1e-9)")

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
        TOL = 1e-9
        sel = "tie" if abs(gPR - gPC) < TOL else ("R" if gPR > gPC else "C")
        diffs = [PR - PC for (_, PR, PC) in res.values()]
        reversal = (min(diffs) < -TOL) and (max(diffs) > TOL)
        for name, (P0, PR, PC) in res.items():
            best = max(PR, PC)
            regret_if_R = best - PR; regret_if_C = best - PC
            regret = {"R": regret_if_R, "C": regret_if_C, "tie": max(regret_if_R, regret_if_C)}[sel]
            bestgain = best - P0
            rows.append(dict(n=n, q=q, e0=e0, e1=e1, c=c, law=name, P0=P0, PR=PR, PC=PC,
                             dR=PR - P0, dC=PC - P0, R_minus_C=PR - PC, gauss_sel=sel,
                             regret=regret, regret_if_R=regret_if_R, regret_if_C=regret_if_C,
                             rel_regret=(regret / bestgain if bestgain > 0 else float("nan")), reversal=reversal))
    import csv
    with open("b_screen_v2_results.csv", "w", newline="") as f:
        wri = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wri.writeheader(); wri.writerows(rows)
    # summary per case
    print(f"{'n':>2} {'q':>4} {'e0':>4} {'e1':>4} {'c':>4} | {'P0 lam0':>8} {'P0 lam1':>8} {'P0 gauss':>8} | {'R-C lam0':>9} {'R-C lam1':>9} {'R-C gauss':>9} | sel  max regret  max relreg  sign-flip")
    flips = 0
    for key in sorted(set((r['n'], r['q'], r['e0'], r['e1'], r['c']) for r in rows)):
        sub = {r['law']: r for r in rows if (r['n'], r['q'], r['e0'], r['e1'], r['c']) == key}
        flip = sub['lam=0']['reversal']
        flips += flip
        mr = max(r['regret'] for r in sub.values()); mrr = np.nanmax([r['rel_regret'] for r in sub.values()])
        print(f"{key[0]:>2} {key[1]:>4} {key[2]:>4} {key[3]:>4} {key[4]:>4} | {sub['lam=0']['P0']:8.4f} {sub['lam=1']['P0']:8.4f} {sub['gauss']['P0']:8.4f} | {sub['lam=0']['R_minus_C']:+9.5f} {sub['lam=1']['R_minus_C']:+9.5f} {sub['gauss']['R_minus_C']:+9.5f} | {sub['gauss']['gauss_sel']:>3}  {mr:10.5f}  {mrr:10.3f}  {'YES' if flip else 'no'}")
    print(f"\ncases with sign reversal of P_R - P_C across the six laws: {flips} / 32")

def n3_exhaustive():
    """At n=3 every kernel entry is a polynomial of degree <= 3 in theta, so P_a depends on G only through
    (m1, m2, m3). With m1=p, m2=q fixed, m3 ranges over [m3(G_lower), m3(G_upper)] (theta^3 is 3-convex),
    and G_lambda has m3 linear in lambda. Hence the lambda sweep covers the whole identified set at n=3."""
    p = 0.10
    print("\n== n=3: exhaustive sweep of the identified set (401 values of lambda) ==")
    lams = np.linspace(0, 1, 401)
    print(f"{'q':>4} {'e0':>4} {'e1':>4} {'c':>4} | {'m3 range':>20} | {'P0 min':>7} {'P0 max':>7} | {'R-C min':>9} {'R-C max':>9} | flip | {'gauss lam*':>10} {'|P(gauss)-P(G_lam*)|':>20}")
    for q, e0, e1, c in itertools.product((0.03, 0.06), (0.02, 0.10), (0.50, 0.90), (0.02, 0.10)):
        tl, wl = G_lower(p, q); tu, wu = G_upper(p, q)
        m3lo = float(np.dot(wl, tl ** 3)); m3hi = float(np.dot(wu, tu ** 3))
        P = np.array([run_case(3, e0, e1, c, *G_mix(l, p, q)) for l in lams])
        diff = P[:, 1] - P[:, 2]
        flip = (diff.min() < -1e-9) and (diff.max() > 1e-9)
        thg, wg, _ = G_gauss(p, q)
        m3g = float(np.dot(wg, thg ** 3)); lam_star = (m3g - m3lo) / (m3hi - m3lo)
        Pg = np.array(run_case(3, e0, e1, c, thg, wg)); Pl = np.array(run_case(3, e0, e1, c, *G_mix(lam_star, p, q)))
        md = np.max(np.abs(Pg - Pl)); assert md < 1e-9, md
        print(f"{q:>4} {e0:>4} {e1:>4} {c:>4} | [{m3lo:.6f}, {m3hi:.6f}] | {P[:,0].min():7.4f} {P[:,0].max():7.4f} | {diff.min():+9.5f} {diff.max():+9.5f} | {'YES' if flip else 'no ':>4} | {lam_star:10.5f} {md:20.2e}")

if __name__ == "__main__":
    main()
    n3_exhaustive()
