"""How does the whole-class (p,q)-ambiguity of the R-vs-C ranking scale with patch number n?

Setting: p=1/10, q=3/100, e0=1/50, e1=9/10, c=1/10, T=50 (the n=3/n=4 reversal setting). For each n = 3..8, search
over discrete hazard laws with exactly matched (p, q) (3- and 4-atom laws; weights solved from the moment equations,
nonnegativity enforced) for the largest and smallest Delta = P_R - P_C. The (p,q)-compatible class is the same for
every n; only the number of free moments (n-2) and the chains change.

Every value reported is a WITNESS: an attained Delta of a feasible law, hence a lower bound on the class maximum
(or upper bound on the class minimum), and the implied deterministic minimax regret min(max(0,-L), max(0,U)) is a
LOWER BOUND on the whole-class regret. Certified whole-class extrema exist only for n=3 and n=4 (n3_certificate,
n4_boundary); the probe compares its n=3,4 witnesses against those as a calibration of search quality. Float search
(numpy, scipy Nelder-Mead), then exact rational re-evaluation of the best witnesses with n3_certificate kernels.
This is not a theorem in n. Fable, 2026-09-23.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
import n4_vertical as V
import n3_certificate as N3

P, Qm = F(1, 10), F(3, 100); E0, E1, C = F(1, 50), F(9, 10), F(1, 10); T = V.T
pf, qf = float(P), float(Qm)
CERTIFIED = {3: (-0.02341633369, 0.01720185959), 4: (-0.12052030700, 0.03209273499)}

def float_parts(n):
    return {a: [np.array([[float(x) for x in row] for row in Ak]) for Ak in V.moment_kernels(a, n, E0, E1, C)] for a in ("R", "C")}

def delta_from_moments(n, parts, mom):
    out = {}
    for a in ("R", "C"):
        K = sum(mom[k] * parts[a][k] for k in range(n + 1)); v = np.zeros(len(K)); v[V.start_state(a, n)] = 1.0
        for _ in range(T): v = v @ K
        out[a] = v[1:].sum()
    return out["R"] - out["C"], out

def weights3(atoms):
    a = np.asarray(atoms, float)
    M = np.vstack([np.ones(3), a, a * a]); 
    try: w = np.linalg.solve(M, np.array([1.0, pf, qf]))
    except np.linalg.LinAlgError: return None
    return w if np.all(w >= -1e-12) else None

def weights4(atoms, w4):
    a = np.asarray(atoms, float); rhs = np.array([1.0, pf, qf]) - w4 * np.array([1.0, a[3], a[3] ** 2])
    M = np.vstack([np.ones(3), a[:3], a[:3] ** 2])
    try: w = np.linalg.solve(M, rhs)
    except np.linalg.LinAlgError: return None
    w = np.append(w, w4); return w if np.all(w >= -1e-12) else None

def law_delta(n, parts, atoms, w):
    mom = np.array([np.sum(w * np.asarray(atoms) ** k) for k in range(n + 1)])
    return delta_from_moments(n, parts, mom)[0]

def search(n, parts, rng, nrand=40000, nrefine=12):
    cands = []
    for _ in range(nrand):
        if rng.random() < 0.6:
            atoms = np.sort(rng.random(3)); w = weights3(atoms)
        else:
            atoms = np.sort(rng.random(4)); w = weights4(atoms, rng.random() * 0.5)
        if w is None or np.ptp(atoms) < 1e-9: continue
        cands.append((law_delta(n, parts, atoms, w), atoms, w))
    cands.sort(key=lambda t: t[0])
    def refine(seed, sign):
        atoms0, w0 = seed[1], seed[2]
        if len(atoms0) == 3:
            def obj(x):
                a = np.clip(np.sort(x), 0, 1); w = weights3(a)
                return 10.0 if w is None else sign * law_delta(n, parts, a, w)
            r = minimize(obj, atoms0, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=4000))
            a = np.clip(np.sort(r.x), 0, 1); w = weights3(a)
        else:
            def obj(x):
                a = np.clip(np.sort(x[:4]), 0, 1); w = weights4(a, min(max(x[4], 0), 1))
                return 10.0 if w is None else sign * law_delta(n, parts, a, w)
            r = minimize(obj, np.append(atoms0, w0[3]), method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-13, maxiter=6000))
            a = np.clip(np.sort(r.x[:4]), 0, 1); w = weights4(a, min(max(r.x[4], 0), 1))
        return (law_delta(n, parts, a, w), a, w) if w is not None else seed
    lo = min((refine(s, +1) for s in cands[:nrefine]), key=lambda t: t[0])
    hi = max((refine(s, -1) for s in cands[-nrefine:]), key=lambda t: t[0])
    return lo, hi, len(cands)

def exact_delta(n, atoms, weights):
    """Exact Delta for a rational law via direct averaging of the exact kernels (n3_certificate)."""
    out = {}
    for a in ("R", "C"):
        K = V.avg_direct(a, n, E0, E1, C, atoms, weights); out[a] = V.persist(K, V.start_state(a, n))
    return out["R"] - out["C"], out

def rationalise(atoms, w, den=10**6):
    """Round atoms to multiples of 1/den, then re-solve the weights EXACTLY so that (p, q) hold exactly; the 4-atom
    case keeps the fourth weight rounded and re-solves the other three. Returns None if any weight is negative."""
    atoms = list(atoms); w = list(w)
    keep = [i for i in range(len(w)) if w[i] > 1e-7]
    if len(keep) == 2:
        # effectively a two-point law {a, b}: round a, then b = (q - p a)/(p - a) makes (p, q) exact
        lo_, hi_ = float(atoms[keep[0]]), float(atoms[keep[1]])
        if hi_ > 1 - 1e-6: b = F(1); a = (P - Qm) / (1 - P)                 # Ghi
        elif lo_ < 1e-6: a = F(0); b = Qm / P                                # Glo
        else: a = F(int(round(lo_ * den)), den); b = (Qm - P * a) / (P - a)
        if not (0 <= a < b <= 1): return None
        return ([a, b], [(b - P) / (b - a), (P - a) / (b - a)])
    if len(keep) == 1: return None
    atoms = [atoms[i] for i in keep]; w = [w[i] for i in keep]
    A = [F(int(round(float(a) * den)), den) for a in atoms]
    if len(A) == 3:
        rhs = [F(1), P, Qm]
    else:
        w4 = F(int(round(float(w[3]) * den)), den); rhs = [F(1) - w4, P - w4 * A[3], Qm - w4 * A[3] ** 2]
    a, b, c = A[:3]
    # solve Vandermonde 3x3 exactly (Lagrange-type formulas)
    def wt(x, y, z, r0, r1, r2): return (r2 - (y + z) * r1 + y * z * r0) / ((x - y) * (x - z))
    W = [wt(a, b, c, *rhs), wt(b, a, c, *rhs), wt(c, a, b, *rhs)] + ([w4] if len(A) == 4 else [])
    return (A, W) if all(x >= 0 for x in W) and all(0 <= x <= 1 for x in A) else None

def main():
    t0 = time.time(); rng = np.random.default_rng(11); rows = []
    for n in range(3, 9):
        parts = float_parts(n)
        (tl, wl), (tu, wu) = N3.laws(P, Qm)
        dlo_end = law_delta(n, parts, [float(x) for x in tl], np.array([float(x) for x in wl]))
        dhi_end = law_delta(n, parts, [float(x) for x in tu], np.array([float(x) for x in wu]))
        lo, hi, ncand = search(n, parts, rng)
        L, U = min(lo[0], dlo_end, dhi_end), max(hi[0], dlo_end, dhi_end)
        Vreg = min(max(0.0, -L), max(0.0, U))
        # exact re-evaluation of the searched witnesses
        ex = {}
        for tag, wit in (("min", lo), ("max", hi)):
            rat = rationalise(wit[1], wit[2])
            if rat is None: ex[tag] = None; continue
            d, PP = exact_delta(n, *rat)
            ex[tag] = {"atoms": [str(x) for x in rat[0]], "weights": [str(x) for x in rat[1]], "Delta_exact_float": float(d), "PR": float(PP["R"]), "PC": float(PP["C"])}
        # persistence level under the endpoint laws for context
        _, outs = delta_from_moments(n, parts, np.array([np.sum(np.array([float(x) for x in wl]) * np.array([float(x) for x in tl]) ** k) for k in range(n + 1)]))
        row = {"n": n, "candidates": ncand, "Delta_min_found": L, "Delta_max_found": U, "regret_lower_bound": Vreg,
               "Delta_Glo": dlo_end, "Delta_Ghi": dhi_end, "endpoints_are_extremes": bool(abs(L - min(dlo_end, dhi_end)) < 1e-9 and abs(U - max(dlo_end, dhi_end)) < 1e-9),
               "min_witness": {"atoms": lo[1].tolist(), "weights": lo[2].tolist()}, "max_witness": {"atoms": hi[1].tolist(), "weights": hi[2].tolist()},
               "exact_reevaluation": ex, "PR_PC_at_Glo": {k: float(v) for k, v in outs.items()},
               "certified_class_range": CERTIFIED.get(n)}
        rows.append(row)
        cert = f"  certified [{CERTIFIED[n][0]:+.5f}, {CERTIFIED[n][1]:+.5f}]" if n in CERTIFIED else ""
        print(f"n={n}: Delta found in [{L:+.5f}, {U:+.5f}]  regret>= {Vreg:.5f}  endpoints Glo {dlo_end:+.5f} Ghi {dhi_end:+.5f}  endpoints-extreme={row['endpoints_are_extremes']}{cert}  exact min {ex['min'] and round(ex['min']['Delta_exact_float'],6)} max {ex['max'] and round(ex['max']['Delta_exact_float'],6)}  P at Glo R {outs['R']:.4f} C {outs['C']:.4f}  ({time.time()-t0:.0f}s)", flush=True)
        print(f"     min witness atoms {np.round(lo[1],5).tolist()} w {np.round(lo[2],5).tolist()} | max witness atoms {np.round(hi[1],5).tolist()} w {np.round(hi[2],5).tolist()}")
    json.dump({"p": str(P), "q": str(Qm), "e0": str(E0), "e1": str(E1), "c": str(C), "T": T, "rows": rows, "seconds": round(time.time() - t0, 1)},
              open(os.path.join(os.path.dirname(__file__), "n_scaling_probe.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
