"""Bounded mechanism diagnostic on the four E-class extremisers (Astra #516). n=4, e0=1/50, e1=9/10, c=1/10, T=50.

Common labelled 16-state occupancy chain (labelled_class_probe.kernel_det), exact rationals. For each extremiser Q:
  K0 = baseline (no action), KR = buffer patch 0, KC = doubled colonisation; mu = all occupied; h = indicator of nonempty.
  (1) One-step identity check: Delta_1 = mu (KR - KC) h  versus  (s/2)[e0^3 p + 3 e0^2 s q + 3 e0 s^2 tau3 + s^3 tau4], s = e1 - e0.
  (2) Telescoping decomposition, forward ordering:  Delta_T = sum_t mu KR^t (KR - KC) KC^(T-1-t) h,
      split as buffering B_t = mu KR^t (KR - K0) KC^(T-1-t) h and colonisation C_t = mu KR^t (KC - K0) KC^(T-1-t) h, so Delta_T = sum_t (B_t - C_t);
      and by the hit count k at step t (K_a - K_b = sum_x Q(x) (K_a(x) - K_b(x)), grouped by |x|).
  (3) Reverse ordering: Delta_T = sum_t mu KC^t (KR - KC) KR^(T-1-t) h, same splits. Totals must equal the certified Delta exactly.
  (4) Also reported: cumulative sums over t, and the occupancy-count distribution of mu KR^t and mu KC^t at t = 0, 5, 10, 25, 49.
Fable, 2026-09-25.
"""
import sys, os, json, itertools
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
from math import comb
import labelled_class_probe as L

n = 4; T = 50; S = 1 << n
e0, e1, c = F(1, 50), F(9, 10), F(1, 10); p = F(1, 10)
CASES = {"q=3/100 E-min (0,0)": (F(3, 100), F(0), F(0)), "q=3/100 E-max (3/100,3/100)": (F(3, 100), F(3, 100), F(3, 100)),
         "q=3/50 E-min (2/75,0)": (F(3, 50), F(2, 75), F(0)), "q=3/50 E-max (3/50,3/50)": (F(3, 50), F(3, 50), F(3, 50))}
CERT = {"q=3/100 E-min (0,0)": None, "q=3/100 E-max (3/100,3/100)": None, "q=3/50 E-min (2/75,0)": None, "q=3/50 E-max (3/50,3/50)": None}

def pattern_law(p, q, x, y):
    pi = [1 - 4 * p + 6 * q - 4 * x + y, 4 * p - 12 * q + 12 * x - 4 * y, 6 * q - 12 * x + 6 * y, 4 * x - 4 * y, y]
    assert all(v >= 0 for v in pi) and sum(pi) == 1
    return {xx: pi[sum(xx)] / comb(4, sum(xx)) for xx in L.PATTERNS}

def det_kernels():
    """Deterministic-pattern kernels for actions 0 (baseline), R (patch 0 buffered), C (2c)."""
    out = {}
    for a in ("0", "R", "C"):
        ks = {}
        for xx in L.PATTERNS:
            if a == "R": dvec = [e0 + (e1 - e0) * xi * (F(1, 2) if i == 0 else 1) for i, xi in enumerate(xx)]; cc = c
            else: dvec = [e0 + (e1 - e0) * xi for xi in xx]; cc = 2 * c if a == "C" else c
            ks[xx] = L.kernel_det(dvec, cc)
        out[a] = ks
    return out

def mat(Q, ks): return [[sum(Q[xx] * ks[xx][r][s] for xx in L.PATTERNS) for s in range(S)] for r in range(S)]
def madd(A, B, sgn=1): return [[a + sgn * b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]
def vm(v, M): return [sum(v[r] * M[r][s] for r in range(S)) for s in range(S)]
def mv(M, v): return [sum(M[r][s] * v[s] for s in range(S)) for r in range(S)]
def dot(u, v): return sum(a * b for a, b in zip(u, v))

def main():
    DK = det_kernels(); results = {}
    for name, (q, x, y) in CASES.items():
        Q = pattern_law(p, q, x, y)
        K = {a: mat(Q, DK[a]) for a in ("0", "R", "C")}
        mu = [F(0)] * S; mu[S - 1] = F(1); h = [F(1)] * S; h[0] = F(0)
        # one-step identity
        d1 = dot(vm(mu, madd(K["R"], K["C"], -1)), h); s = e1 - e0
        d1_formula = (s / 2) * (e0 ** 3 * p + 3 * e0 ** 2 * s * q + 3 * e0 * s ** 2 * x + s ** 3 * y)
        # certified T-step Delta
        vR, vC = mu[:], mu[:]
        for _ in range(T): vR = vm(vR, K["R"]); vC = vm(vC, K["C"])
        Delta = dot(vR, h) - dot(vC, h)
        # forward and backward vectors
        fwdR = [mu]; fwdC = [mu]
        for _ in range(T): fwdR.append(vm(fwdR[-1], K["R"])); fwdC.append(vm(fwdC[-1], K["C"]))
        bwdC = [h]; bwdR = [h]
        for _ in range(T): bwdC.append(mv(K["C"], bwdC[-1])); bwdR.append(mv(K["R"], bwdR[-1]))   # bwdC[j] = KC^j h
        dR0 = madd(K["R"], K["0"], -1); dC0 = madd(K["C"], K["0"], -1)
        # per-pattern differences for hit-count split
        dR0x = {xx: Q[xx] and [[Q[xx] * (DK["R"][xx][r][s_] - DK["0"][xx][r][s_]) for s_ in range(S)] for r in range(S)] for xx in L.PATTERNS}
        dC0x = {xx: Q[xx] and [[Q[xx] * (DK["C"][xx][r][s_] - DK["0"][xx][r][s_]) for s_ in range(S)] for r in range(S)] for xx in L.PATTERNS}
        rows = []
        for ordering, left, right in (("forward (R past, C future)", fwdR, bwdC), ("reverse (C past, R future)", fwdC, bwdR)):
            Bt, Ct, byk = [], [], []
            for t in range(T):
                a = left[t]; b = right[T - 1 - t]
                Bt.append(dot(vm(a, dR0), b)); Ct.append(dot(vm(a, dC0), b))
                bk = [F(0)] * (n + 1); ck = [F(0)] * (n + 1)
                for xx in L.PATTERNS:
                    if Q[xx] == 0: continue
                    k = sum(xx); bk[k] += dot(vm(a, dR0x[xx]), b); ck[k] += dot(vm(a, dC0x[xx]), b)
                byk.append((bk, ck))
            total = sum(Bt) - sum(Ct); assert total == Delta, (ordering, float(total), float(Delta))
            Bsum, Csum = sum(Bt), sum(Ct)
            byk_B = [sum(bk[k] for bk, _ in byk) for k in range(n + 1)]; byk_C = [sum(ck[k] for _, ck in byk) for k in range(n + 1)]
            rows.append({"ordering": ordering, "buffering_total": float(Bsum), "colonisation_total": float(Csum), "Delta": float(total),
                         "buffering_by_t": [float(v) for v in Bt], "colonisation_by_t": [float(v) for v in Ct],
                         "cum_Delta_by_t": [float(sum(Bt[:i + 1]) - sum(Ct[:i + 1])) for i in range(T)],
                         "buffering_by_hitcount": [float(v) for v in byk_B], "colonisation_by_hitcount": [float(v) for v in byk_C],
                         "net_by_hitcount": [float(bb - cc) for bb, cc in zip(byk_B, byk_C)]})
        occ = {}
        for t in (0, 1, 5, 10, 25, 49):
            occ[t] = {"R": [float(sum(fwdR[t][st] for st in range(S) if bin(st).count("1") == k)) for k in range(n + 1)],
                      "C": [float(sum(fwdC[t][st] for st in range(S) if bin(st).count("1") == k)) for k in range(n + 1)]}
        results[name] = {"q": str(q), "tau3": str(x), "tau4": str(y), "pi": [str(v) for v in [Q[xx] * comb(4, sum(xx)) for xx in [(0,0,0,0),(1,0,0,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)]]],
                         "Delta_1_exact": str(d1), "Delta_1_formula": str(d1_formula), "identity_holds": d1 == d1_formula,
                         "Delta_T_exact": str(Delta), "Delta_T_float": float(Delta), "PR": float(dot(vR, h)), "PC": float(dot(vC, h)),
                         "decompositions": rows, "occupancy_count_distribution": occ}
        print(f"\n== {name}: pi = {[str(v) for v in results[name]['pi']]}")
        print(f"   Delta_1 = {float(d1):+.8f} (formula {float(d1_formula):+.8f}, identity {d1 == d1_formula}); Delta_50 = {float(Delta):+.8f}  P_R {float(dot(vR,h)):.5f} P_C {float(dot(vC,h)):.5f}")
        for r in rows:
            print(f"   {r['ordering']}: buffering {r['buffering_total']:+.6f}  colonisation {r['colonisation_total']:+.6f}  -> Delta {r['Delta']:+.6f}")
            print(f"      by hit count k=0..4: buffering {[round(v,5) for v in r['buffering_by_hitcount']]}  colonisation {[round(v,5) for v in r['colonisation_by_hitcount']]}  net {[round(v,5) for v in r['net_by_hitcount']]}")
            cum = r['cum_Delta_by_t']; print(f"      cumulative Delta after t=1,5,10,20,30,50: {[round(cum[i],5) for i in (0,4,9,19,29,49)]}")
        for t, o in occ.items(): print(f"      occupancy count dist at t={t}: R {[round(v,4) for v in o['R']]}  C {[round(v,4) for v in o['C']]}")
    json.dump({"n": n, "T": T, "e0": str(e0), "e1": str(e1), "c": str(c), "p": str(p), "results": results}, open(os.path.join(os.path.dirname(__file__), "mechanism_decomposition.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
