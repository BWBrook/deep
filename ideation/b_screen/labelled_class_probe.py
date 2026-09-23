"""Does dropping exchangeability change what complete pairwise hit laws identify? Witness probe, n=4.

Class L(p,q): ALL laws Q on labelled hit vectors x in {0,1}^4 with P(x_i=1)=p for every patch and P(x_i=x_j=1)=q for
every pair, i.e. every complete pairwise hit table equal to the exchangeable one. Dimension 16-11=5 (the
exchangeable (p,q)-class is a 2-parameter subset: m3, m4). The exact machinery here is the labelled 16-state chain
(no occupancy-count aggregation): each occupied patch i dies with probability e0+(e1-e0)x_i (R_i: excess halved at
patch i), survivors recolonise empties independently w.p. 1-(1-c)^m (C: 2c), all-occupied start, absorbing 0, T=50.

Two ranking objects: Delta_1 = P_{R_1} - P_C (buffer a fixed patch), and Delta* = max_i P_{R_i} - P_C (buffer the best
patch; equals Delta_1 on exchangeable laws). Search: polytope vertices from random-objective LPs (exact rational
reconstruction from the support), random convex combinations, Nelder-Mead refinement in null-space coordinates.
Every reported extreme is an attained value (float, or exact for vertex witnesses): lower bounds on the class range,
never class extrema. Compared with the certified exchangeable-class Delta ranges of the 16-case design
(n4_boundary_results.json). Fable, 2026-09-23.
"""
import sys, os, json, time, itertools
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog, minimize

n = 4; T = 50; PATTERNS = list(itertools.product((0, 1), repeat=n)); S = 1 << n

def kernel_det(dvec, c):
    """Labelled kernel for a deterministic hit vector, i.e. per-patch death probabilities dvec; works for float or Fraction."""
    one = 1 if isinstance(c, F) else 1.0
    K = [[0 * one for _ in range(S)] for _ in range(S)]
    for state in range(S):
        occ = [i for i in range(n) if state >> i & 1]
        for surv_mask in range(1 << len(occ)):
            ps = one; survivors = []
            for idx, i in enumerate(occ):
                if surv_mask >> idx & 1: ps = ps * (1 - dvec[i]); survivors.append(i)
                else: ps = ps * dvec[i]
            m = len(survivors); r = 1 - (1 - c) ** m
            empties = [i for i in range(n) if i not in survivors]
            for col_mask in range(1 << len(empties)):
                pc = one; new = 0
                for i in survivors: new |= 1 << i
                for idx, i in enumerate(empties):
                    if col_mask >> idx & 1: pc = pc * r; new |= 1 << i
                    else: pc = pc * (1 - r)
                K[state][new] += ps * pc
    return K

def action_kernels(e0, e1, c):
    """For each action a in {C, R_0..R_3}: list over patterns of deterministic kernels."""
    out = {}
    for a in ("C", 0, 1, 2, 3):
        ks = []
        for x in PATTERNS:
            if a == "C": dvec = [e0 + (e1 - e0) * xi for xi in x]; cc = 2 * c
            else: dvec = [e0 + (e1 - e0) * xi * (F(1, 2) if isinstance(e0, F) else 0.5) if i == a else e0 + (e1 - e0) * xi for i, xi in enumerate(x)]; cc = c
            ks.append(kernel_det(dvec, cc))
        out[a] = ks
    return out

def persist(K, T=T):
    v = [0 * K[0][0]] * S; v[S - 1] = 1 if isinstance(K[0][0], F) else 1.0
    for _ in range(T): v = [sum(v[r] * K[r][s] for r in range(S)) for s in range(S)]
    return 1 - v[0]

def persist_np(K):
    v = np.zeros(S); v[S - 1] = 1.0
    for _ in range(T): v = v @ K
    return 1 - v[0]

def deltas_float(Qvec, AK):
    """AK: {action: 3-D array [pattern, S, S]}; returns (Delta_1, Delta*, PC, [PR_i])."""
    Ks = {a: np.tensordot(Qvec, AK[a], axes=(0, 0)) for a in AK}
    PC = persist_np(Ks["C"]); PR = [persist_np(Ks[i]) for i in range(n)]
    return PR[0] - PC, max(PR) - PC, PC, PR

def constraints(p, q):
    A = [[1.0] * S]; b = [1.0]
    for i in range(n): A.append([float(x[i]) for x in PATTERNS]); b.append(p)
    for i, j in itertools.combinations(range(n), 2): A.append([float(x[i] * x[j]) for x in PATTERNS]); b.append(q)
    return np.array(A), np.array(b)

def vertices(A, b, rng, ntry=600):
    seen = {}; 
    for _ in range(ntry):
        cobj = rng.normal(size=S)
        r = linprog(cobj, A_eq=A, b_eq=b, bounds=[(0, 1)] * S, method="highs")
        if r.status == 0:
            key = tuple(np.round(r.x, 9)); seen[key] = r.x
    return list(seen.values())

def exact_vertex(Qf, p, q):
    """Reconstruct a vertex exactly from its support by solving the equality system with Fractions."""
    supp = [k for k in range(S) if Qf[k] > 1e-9]
    rows = [[F(1)] * len(supp)]; rhs = [F(1)]
    for i in range(n): rows.append([F(PATTERNS[k][i]) for k in supp]); rhs.append(p)
    for i, j in itertools.combinations(range(n), 2): rows.append([F(PATTERNS[k][i] * PATTERNS[k][j]) for k in supp]); rhs.append(q)
    M = [row[:] + [r] for row, r in zip(rows, rhs)]; m_, k_ = len(M), len(supp); piv_cols = []; r_ = 0
    for col in range(k_):
        pr = next((i for i in range(r_, m_) if M[i][col] != 0), None)
        if pr is None: continue
        M[r_], M[pr] = M[pr], M[r_]; pv = M[r_][col]; M[r_] = [v / pv for v in M[r_]]
        for i in range(m_):
            if i != r_ and M[i][col] != 0:
                f = M[i][col]; M[i] = [a - f * bb for a, bb in zip(M[i], M[r_])]
        piv_cols.append(col); r_ += 1
    if r_ < k_: return None   # not a vertex (free variables remain)
    if any(all(v == 0 for v in M[i][:k_]) and M[i][k_] != 0 for i in range(r_, m_)): return None
    sol = [F(0)] * k_
    for i, col in enumerate(piv_cols): sol[col] = M[i][k_]
    if any(v < 0 for v in sol): return None
    Q = [F(0)] * S
    for k, v in zip(supp, sol): Q[k] = v
    if any(abs(float(Q[k]) - Qf[k]) > 1e-7 for k in range(S)): return None
    return Q

def exact_deltas(Q, e0, e1, c):
    AK = action_kernels(e0, e1, c)
    Ks = {a: [[sum(Q[k] * AK[a][k][r][s] for k in range(S)) for s in range(S)] for r in range(S)] for a in AK}
    PC = persist(Ks["C"]); PR = [persist(Ks[i]) for i in range(n)]
    return PR[0] - PC, max(PR) - PC, PC, PR

def main():
    t0 = time.time(); rng = np.random.default_rng(3)
    J = json.load(open(os.path.join(os.path.dirname(__file__), "n4_boundary_results.json")))
    cert = {(cs["q"], cs["e0"], cs["e1"], cs["c"]): (cs["whole_class"]["Delta"]["min"], cs["whole_class"]["Delta"]["max"]) for cs in J["cases"]}
    p = F(1, 10); rows = []
    for (q, e0, e1, c) in itertools.product((F(3, 100), F(6, 100)), (F(1, 50), F(1, 10)), (F(1, 2), F(9, 10)), (F(1, 50), F(1, 10))):
        pf, qf = float(p), float(q); e0f, e1f, cf = float(e0), float(e1), float(c)
        AKf = {a: np.array([[[float(v) for v in row] for row in K] for K in ks]) for a, ks in action_kernels(e0f, e1f, cf).items()}
        A, b = constraints(pf, qf); verts = vertices(A, b, rng)
        evals = []
        for V_ in verts:
            d1, ds, PC, PR = deltas_float(V_, AKf); evals.append((d1, ds, V_, "vertex"))
        # random convex combinations
        Vm = np.array(verts)
        for _ in range(1500):
            w = rng.dirichlet(np.full(len(verts), 0.3)); Qv = w @ Vm
            d1, ds, PC, PR = deltas_float(Qv, AKf); evals.append((d1, ds, Qv, "mixture"))
        # Nelder-Mead refinement in null-space coordinates from the best/worst few
        N_ = np.linalg.svd(A)[2][np.linalg.matrix_rank(A):].T   # 16 x 5 null-space basis
        def refine(Q0, sign, which):
            def obj(z):
                Qv = Q0 + N_ @ z
                if Qv.min() < 0: return 10.0 + 100 * (-Qv.min())
                d1, ds, _, _ = deltas_float(Qv, AKf); return sign * (d1 if which == 0 else ds)
            r = minimize(obj, np.zeros(5), method="Nelder-Mead", options=dict(xatol=1e-9, fatol=1e-12, maxiter=3000))
            Qv = Q0 + N_ @ r.x
            if Qv.min() < 0: return None
            d1, ds, _, _ = deltas_float(Qv, AKf); return (d1, ds, Qv, "refined")
        for which in (0, 1):
            key = (lambda t: t[which])
            for seed in sorted(evals, key=key)[:4]: 
                rr = refine(seed[2], +1, which); evals.append(rr) if rr else None
            for seed in sorted(evals, key=key)[-4:]:
                rr = refine(seed[2], -1, which); evals.append(rr) if rr else None
        lo1 = min(evals, key=lambda t: t[0]); hi1 = max(evals, key=lambda t: t[0]); los = min(evals, key=lambda t: t[1]); his = max(evals, key=lambda t: t[1])
        # exact re-evaluation for vertex-type witnesses of Delta_1 extremes
        ex = {}
        for tag, wit in (("min_D1", lo1), ("max_D1", hi1), ("min_Dstar", los), ("max_Dstar", his)):
            Qx = exact_vertex(wit[2], p, q) if wit[3] == "vertex" else None
            if Qx is not None:
                d1, ds, PC, PR = exact_deltas(Qx, e0, e1, c)
                ex[tag] = {"Q": {"".join(map(str, PATTERNS[k])): str(Qx[k]) for k in range(S) if Qx[k] != 0}, "Delta_1": float(d1), "Delta_star": float(ds), "PC": float(PC), "PR": [float(v) for v in PR]}
            else: ex[tag] = None
        cmin, cmax = cert[(str(q), str(e0), str(e1), str(c))]
        exch_ident = "C" if cmax < 0 else "R" if cmin > 0 else "unidentified"
        lab_ident1 = "C" if hi1[0] < 0 else "R" if lo1[0] > 0 else "unidentified(witnessed)"
        lab_idents = "C" if his[1] < 0 else "R" if los[1] > 0 else "unidentified(witnessed)"
        row = {"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), "n_vertices": len(verts),
               "exchangeable_certified_Delta_range": [cmin, cmax], "exchangeable_identification": exch_ident,
               "labelled_Delta1_found": [lo1[0], hi1[0]], "labelled_Dstar_found": [los[1], his[1]],
               "labelled_identification_Delta1_from_witnesses": lab_ident1, "labelled_identification_Dstar_from_witnesses": lab_idents,
               "witness_types": {"min_D1": lo1[3], "max_D1": hi1[3], "min_Dstar": los[3], "max_Dstar": his[3]},
               "exact_vertex_reevaluation": ex}
        rows.append(row)
        flag = "  <-- CHANGES" if (exch_ident != "unidentified" and (lab_ident1.startswith("unid") or lab_idents.startswith("unid"))) else ""
        print(f"q={str(q):>5} e0={str(e0):>4} e1={str(e1):>4} c={str(c):>4} | exch [{cmin:+.5f},{cmax:+.5f}] {exch_ident:<12} | labelled D1 [{lo1[0]:+.5f},{hi1[0]:+.5f}] D* [{los[1]:+.5f},{his[1]:+.5f}] -> D1 {lab_ident1}, D* {lab_idents}{flag}  ({len(verts)} vertices, {time.time()-t0:.0f}s)", flush=True)
    json.dump({"n": n, "T": T, "p": "1/10", "rows": rows, "seconds": round(time.time() - t0, 1)}, open(os.path.join(os.path.dirname(__file__), "labelled_class_probe.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
