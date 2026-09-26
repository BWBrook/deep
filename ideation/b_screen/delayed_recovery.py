"""Delayed-recovery pilot (synthesis/delayed_recovery_contract.md). n=4, p=1/10, e0=1/50, e1=9/10, c=1/10, T=50.

State per patch: occupancy O_i in {0,1}, establishment-block clock a_i in {0..d}. Step order given the hit vector X:
 1. clocks: a_i <- d if X_i = 1 else max(a_i - 1, 0);
 2. extirpation, independent, occupied patches: baseline/C: e0 + s X_i; R(r): patch 0 has e0 + (1 - r) s X_0;
 3. donors m: D1 = all occupied survivors; D2 = occupied survivors with updated clock 0;
 4. colonisation of empty clock-zero patches, independently w.p. 1 - (1 - c_A)^m, c_C = g c, c_0 = c_R = c.
Regional extinction (all empty) is absorbing; persistence = mass on states with some O_i = 1 after T steps.

Implementation: explicit labelled kernels per deterministic hit pattern (exact gmpy2 rationals, sparse); the aggregate
kernel for an exchangeable law is sum_k pi_k K_k where K_k averages over the C(4,k) patterns with k hits. Aggregation
classes: C/baseline = multiset of (O,a) over the four patches; R = (O,a) of patch 0 times the multiset over the other
three. Lumpability is checked over EVERY labelled state (class-level transition distributions identical within a class)
for each fixed-count law. All-empty configurations are collapsed to one absorbing state. d = 0 must reproduce the
existing kernels and persistence values exactly. Fable, 2026-09-26.
"""
import sys, os, json, time, itertools, resource, argparse
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
from math import comb
from gmpy2 import mpq as Q
import n4_vertical as V
import exchangeable_certificate as EC

n = 4; T = 50; P = F(1, 10); E0, E1, C = Q(1, 50), Q(9, 10), Q(1, 10); S_ = E1 - E0
PATTERNS = list(itertools.product((0, 1), repeat=n))

def clocks_after(state, x, d):
    return tuple(d if x[i] else max(state[i][1] - 1, 0) for i in range(n))

def step_dist(state, x, d, action, r, g, variant):
    """Distribution over next labelled states after one step with hit vector x."""
    a2 = clocks_after(state, x, d)
    occ = [i for i in range(n) if state[i][0] == 1]
    cA = C * g if action == "C" else C
    out = {}
    for surv_mask in range(1 << len(occ)):
        ps = Q(1); surv = []
        for idx, i in enumerate(occ):
            di = E0 + (1 - r) * S_ * x[i] if (action == "R" and i == 0) else E0 + S_ * x[i]
            if surv_mask >> idx & 1: ps *= (1 - di); surv.append(i)
            else: ps *= di
        if ps == 0: continue
        m = len(surv) if variant == "D1" else sum(1 for i in surv if a2[i] == 0)
        rcol = 1 - (1 - cA) ** m
        empties = [i for i in range(n) if i not in surv and a2[i] == 0]
        for col_mask in range(1 << len(empties)):
            pc = Q(1); occ2 = set(surv)
            for idx, i in enumerate(empties):
                if col_mask >> idx & 1: pc *= rcol; occ2.add(i)
                else: pc *= (1 - rcol)
            if pc == 0: continue
            nxt = tuple((1 if i in occ2 else 0, a2[i]) for i in range(n))
            out[nxt] = out.get(nxt, Q(0)) + ps * pc
    return out

def cls_of(state, action):
    if action == "R": return (state[0], tuple(sorted(state[1:])))
    return tuple(sorted(state))

def collapse(state):
    """All-empty configurations -> the single absorbing class 'EXT'."""
    return "EXT" if all(o == 0 for o, _ in state) else state

def build(d, action, r, g, variant):
    """Labelled reachable states, per-count aggregate kernels K_k (k=0..4) on classes, lumpability check, stats."""
    t0 = time.time()
    start = tuple((1, 0) for _ in range(n))
    # reachable labelled states under the union of all hit patterns
    seen = {start}; frontier = [start]; trans = {}   # trans[(state, x)] = dist
    while frontier:
        s = frontier.pop()
        if collapse(s) == "EXT": continue
        for x in PATTERNS:
            dist = step_dist(s, x, d, action, r, g, variant); trans[(s, x)] = dist
            for nxt in dist:
                if nxt not in seen: seen.add(nxt); frontier.append(nxt)
    live = [s for s in seen if collapse(s) != "EXT"]
    classes = {}; 
    for s in live: classes.setdefault(cls_of(s, action), []).append(s)
    idx = {c: i for i, c in enumerate(sorted(classes, key=repr))}; ext = len(idx)   # extinct class last
    S = ext + 1
    Ks = []
    for k in range(n + 1):
        pats = [x for x in PATTERNS if sum(x) == k]; w = Q(1, len(pats))
        K = [[Q(0)] * S for _ in range(S)]; K[ext][ext] = Q(1)
        for c, members in classes.items():
            rows = []
            for s in members:
                row = {}
                for x in pats:
                    for nxt, pr in trans[(s, x)].items():
                        j = ext if collapse(nxt) == "EXT" else idx[cls_of(nxt, action)]
                        row[j] = row.get(j, Q(0)) + w * pr
                rows.append(row)
            assert all(rr == rows[0] for rr in rows[1:]), ("lumpability failure", d, action, k, c)
            for j, pr in rows[0].items(): K[idx[c]][j] = pr
        for i in range(S): assert sum(K[i]) == 1 and min(K[i]) >= 0
        Ks.append(K)
    nnz = sum(1 for K in Ks for row in K for v in row if v != 0)
    stats = {"labelled_reachable": len(seen), "classes": ext, "states_with_ext": S, "nonzeros_per_count_kernel": nnz // (n + 1),
             "density": nnz / ((n + 1) * S * S), "build_s": round(time.time() - t0, 1)}
    return Ks, idx, ext, start, stats

def pi_of(p, q, x, y):
    return [1 - 4 * p + 6 * q - 4 * x + y, 4 * p - 12 * q + 12 * x - 4 * y, 6 * q - 12 * x + 6 * y, 4 * x - 4 * y, y]

def persist(Ks, idx, ext, start, pi, action):
    S = ext + 1; K = [[sum(Q(pi[k]) * Ks[k][i][j] for k in range(n + 1)) for j in range(S)] for i in range(S)]
    v = [Q(0)] * S; v[idx[cls_of(start, action)]] = Q(1)
    for _ in range(T): v = [sum(v[i] * K[i][j] for i in range(S)) for j in range(S)]
    return 1 - v[ext]

def labelled_persist(d, action, r, g, variant, Qlaw):
    """Independent check: propagate the labelled chain directly with a pattern law Qlaw (dict pattern -> prob)."""
    start = tuple((1, 0) for _ in range(n)); v = {start: Q(1)}
    for _ in range(T):
        nv = {}
        for s, pr in v.items():
            if pr == 0 or collapse(s) == "EXT": nv[s] = nv.get(s, Q(0)) + pr; continue
            for x, qx in Qlaw.items():
                if qx == 0: continue
                for nxt, p2 in step_dist(s, x, d, action, r, g, variant).items():
                    nv[nxt] = nv.get(nxt, Q(0)) + pr * qx * p2
        v = nv
    return sum(pr for s, pr in v.items() if collapse(s) != "EXT")

def laws_for_q(q):
    """Glo, Ghi (scalar, as (tau3,tau4)) and the two E extremisers (polygon vertices) for this q; duplicates reported."""
    xL, xU = V.xbounds(P, q); L = {}
    for nm, x in (("Glo", xL), ("Ghi", xU)):
        y = V.Lfun(x, P, q); L[nm] = (F(x), F(y))
    verts = EC.polygon(P, q)
    L["E-min"] = verts[0]; L["E-max"] = max(verts, key=lambda t: (t[0], t[1]))
    return L

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--delays", default="0,2"); ap.add_argument("--grid", action="store_true"); ap.add_argument("--qs", default="3/100,3/50")
    args = ap.parse_args(); delays = [int(v) for v in args.delays.split(",")]; qs = [F(v) for v in args.qs.split(",")]
    out = {"n": n, "T": T, "p": str(P), "e0": str(E0), "e1": str(E1), "c": str(C), "checks": {}, "results": [], "footprint": {}}
    RG = [(F(1, 4), F(2)), (F(1, 2), F(2)), (F(3, 4), F(2)), (F(1, 2), F(3, 2)), (F(1, 2), F(3))] if args.grid else [(F(1, 2), F(2))]
    rs = sorted({r for r, _ in RG}); gs = sorted({g for _, g in RG})
    # ---- d = 0 control: exact agreement with existing reference values ----
    t = time.time()
    Ks0 = {a: build(0, a, F(1, 2) if a == "R" else F(0), F(2) if a == "C" else F(1), "D1") for a in ("0", "R", "C")}
    ref = {}
    for q in qs:
        laws = laws_for_q(q)
        for nm, (x, y) in laws.items():
            pi = pi_of(F(P), q, x, y)
            vals = {a: persist(*Ks0[a][:4], pi, a) for a in ("0", "R", "C")}
            partsR = V.affine_parts("R", n, F(E0), F(E1), F(C), P, q); partsC = V.affine_parts("C", n, F(E0), F(E1), F(C), P, q)
            refR = V.persist(V.kernel_at(*partsR, x, y), V.start_state("R", n)); refC = V.persist(V.kernel_at(*partsC, x, y), V.start_state("C", n))
            assert Q(refR.numerator, refR.denominator) == vals["R"] and Q(refC.numerator, refC.denominator) == vals["C"], ("d=0 mismatch", q, nm)
            ref[(str(q), nm)] = {a: float(v) for a, v in vals.items()}
    out["checks"]["d0_exact_agreement_with_reference"] = True; out["checks"]["d0_seconds"] = round(time.time() - t, 1)
    print(f"d=0: exact agreement with existing R/C kernels and persistence at all {len(ref)} laws ({time.time()-t:.1f}s)", flush=True)
    # ---- delayed runs ----
    for d in delays:
        if d == 0: continue
        for variant in ("D1", "D2"):
            t = time.time(); built = {}
            built["0"] = build(d, "0", F(0), F(1), variant)
            for r in rs: built[("R", r)] = build(d, "R", r, F(1), variant)
            for g in gs: built[("C", g)] = build(d, "C", F(0), g, variant)
            out["footprint"][f"d={d} {variant}"] = {"C_or_0": built["0"][4], "R": built[("R", rs[0])][4], "build_all_kernels_s": round(time.time() - t, 1),
                                                    "peak_rss_MB": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, 1)}
            print(f"d={d} {variant}: classes C/0 {built['0'][4]['classes']} (+ext), R {built[('R', rs[0])][4]['classes']} (+ext); labelled reachable C {built['0'][4]['labelled_reachable']} R {built[('R', rs[0])][4]['labelled_reachable']}; density C {built['0'][4]['density']:.3f} R {built[('R', rs[0])][4]['density']:.3f}; build {time.time()-t:.1f}s; peak RSS {out['footprint'][f'd={d} {variant}']['peak_rss_MB']} MB", flush=True)
            # affine/lumpability cross-check against direct labelled propagation at two interior rational points (q=3/100)
            if variant == "D1" or True:
                q = qs[0]; t2 = time.time()
                for (x, y) in ((F(1, 100), F(1, 200)), (F(2, 100), F(3, 200))):
                    pi = pi_of(F(P), q, x, y); assert min(pi) >= 0
                    Qlaw = {xx: Q(pi[sum(xx)].numerator, pi[sum(xx)].denominator) / comb(n, sum(xx)) for xx in PATTERNS}
                    for key, (a, r, g) in (("0", ("0", F(0), F(1))), (("R", rs[0]), ("R", rs[0], F(1))), (("C", gs[0]), ("C", F(0), gs[0]))):
                        agg = persist(*built[key][:4], pi, a); lab = labelled_persist(d, a, r, g, variant, Qlaw)
                        assert agg == lab, ("aggregate vs labelled mismatch", d, variant, a, x, y)
                out["checks"][f"d={d} {variant} aggregate_equals_labelled_at_2_interior_points"] = True
                print(f"   aggregate == direct labelled propagation at 2 interior points, actions 0/R/C ({time.time()-t2:.1f}s)", flush=True)
            t3 = time.time()
            for q in qs:
                laws = laws_for_q(q); W = {}
                for nm, (x, y) in laws.items():
                    pi = pi_of(F(P), q, x, y)
                    P0 = persist(*built["0"][:4], pi, "0")
                    PR = {str(r): persist(*built[("R", r)][:4], pi, "R") for r in rs}; PC = {str(g): persist(*built[("C", g)][:4], pi, "C") for g in gs}
                    assert all(v >= P0 for v in PR.values()) and all(v >= P0 for v in PC.values())
                    W[nm] = (P0, PR, PC)
                    rec = {"d": d, "variant": variant, "q": str(q), "law": nm, "tau3": str(x), "tau4": str(y), "P0": float(P0),
                           "PR": {k: float(v) for k, v in PR.items()}, "PC": {k: float(v) for k, v in PC.items()},
                           "P0_d0": ref[(str(q), nm)]["0"], "PR_d0": ref[(str(q), nm)]["R"], "PC_d0": ref[(str(q), nm)]["C"]}
                    out["results"].append(rec)
                # finite-set summaries per (r, g)
                for r in rs:
                    for g in gs:
                        if (r, g) not in RG: continue
                        deltas = {nm: W[nm][1][str(r)] - W[nm][2][str(g)] for nm in W}
                        a_, b_ = min(deltas.values()), max(deltas.values())
                        regret = min(max(Q(0), -a_), max(Q(0), b_)); minR = min(W[nm][1][str(r)] for nm in W); minC = min(W[nm][2][str(g)] for nm in W)
                        summ = {"d": d, "variant": variant, "q": str(q), "r": str(r), "g": str(g), "Delta_by_law": {nm: float(v) for nm, v in deltas.items()},
                                "min_Delta": float(a_), "max_Delta": float(b_), "regret_W_pp": float(100 * regret),
                                "maximin_W": "R" if minR > minC else "C" if minC > minR else "tie", "min_PR_W": float(minR), "min_PC_W": float(minC),
                                "minimax_regret_W": ("R" if -a_ < b_ else "C") if a_ < 0 < b_ else ("R" if a_ >= 0 else "C")}
                        out["results"].append({"summary": summ})
                        print(f"   q={q} r={r} g={g}: Delta by law " + ", ".join(f"{nm} {float(v):+.4f}" for nm, v in deltas.items()) + f" | regret over W {float(100*regret):.3f} pp, minimax-regret {summ['minimax_regret_W']}, maximin {summ['maximin_W']} (min PR {float(minR):.4f}, min PC {float(minC):.4f})", flush=True)
                for nm in W:
                    print(f"      {nm} (tau3={laws[nm][0]}, tau4={laws[nm][1]}): P0 {float(W[nm][0]):.4f} [d0 {ref[(str(q), nm)]['0']:.4f}]  PR(1/2) {float(W[nm][1]['1/2']):.4f} [d0 {ref[(str(q), nm)]['R']:.4f}]  PC(2) {float(W[nm][2]['2']):.4f} [d0 {ref[(str(q), nm)]['C']:.4f}]")
            print(f"   propagation for all laws/grid: {time.time()-t3:.1f}s", flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(__file__), f"delayed_recovery_pilot.json"), "w"), indent=1, default=str)

if __name__ == "__main__":
    main()
