"""Shape restriction: is the n=4 action ranking identified within the class of UNIMODAL hit laws?

Setting: the one unidentified design case, p = 1/10, q = 3/100, e0 = 1/50, e1 = 9/10, c = 1/10, T = 50.

Class (Khinchin): theta on [0,1] is unimodal about a mode m in [0,1] iff theta = m + U*Y with U ~ Uniform(0,1)
independent of Y, Y supported on [-m, 1-m] (Khinchin 1938; Feller vol. II, XV.9; Dharmadhikari & Joag-dev 1988,
Thm 1.3). This class is the closed convex hull of uniforms sharing endpoint m: it admits a boundary mode (m = 0 or
m = 1, i.e. monotone densities), an atom at the mode (Y = 0 with positive mass), and mixtures of such. It excludes
atoms away from the mode and any density with two local maxima. All results below are exact rationals.

Argument. Certified inputs (n4_vertical_all.json, n4_boundary_results.json, both Astra-audited):
  (i)  dDelta/dm4 > 0 over the whole (m3, m4) lens for this case;
  (ii) on the upper curve m4 = U(m3), Delta is strictly increasing in m3 with a single root, bracketed exactly in
       m3 in [4177210314157/206158430208000, 50126523769919/2473901162496000] ~ 0.0202621368.
Hence for any hit law G in the conditional-iid class with (p, q) fixed and m3(G) <= x* := 4177210314157/206158430208000,
Delta(G) <= Delta(m3, U(m3)) < 0, i.e. C is strictly preferred. So it suffices to bound sup m3 over unimodal G.

Moments under Khinchin: E theta^k = sum_j C(k,j) m^{k-j} E[Y^j]/(j+1). With y_j = E Y^j:
  p  = m + y1/2                                  -> y1 = 2(p - m)
  q  = m^2 + m y1 + y2/3                         -> y2 = 3(q - m^2 - m y1)
  m3 = m^3 + (3/2) m^2 y1 + m y2 + y3/4.
Given (y1, y2) on [A, B] = [-m, 1-m], y3 is maximised by the two-point law on {a, B}, a = y1 - var/(B - y1)
(upper principal representation; Karlin & Studden 1966, ch. II), feasible iff a >= A. It is the shifted analogue
of the construction giving Ghi = {a, 1} on the full class; Ghi itself (two atoms) is bimodal, not unimodal, so the
full-class extreme is not attainable here. m3max_uni(m) is a rational function of m; its supremum over the feasible
mode set is certified by exact root counting (Sturm sequences) of polynomial numerators.

Also reported: a rigorous upper bound on Delta over the unimodal class, Delta(G) <= Delta(x, U(x)) at any rational
x >= sup m3_uni (Delta on the upper curve is certified increasing in m3), i.e. the smallest margin by which C beats R.

Outputs: feasible mode interval; m3max_uni(m) and m3min_uni(m); exact certificate that m3max_uni(m) < x* for all
feasible m; the extremal unimodal law (mixture of two uniforms) at the maximising mode with its exact m4 and exact
Delta; random-search sanity checks of the y3 bound and the Khinchin moment identities; the same unimodal m3 range
for q = 6/100 (p = 1/10) for reference.
"""
import sys, os, json, time, random
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
from math import comb
import sympy as sp
import n4_vertical as V

P, Qm = F(1, 10), F(3, 100)
E0, E1, C = F(1, 50), F(9, 10), F(1, 10)
XSTAR = F(4177210314157, 206158430208000)        # lower end of the certified upper-curve root bracket in m3
XSTAR_HI = F(50126523769919, 2473901162496000)   # upper end
LOWER_ROOT = (F(8879837253229, 412316860416000), F(53279023519409, 2473901162496000))

# ---------- Khinchin reduction ----------
def y12(m, p, q):
    y1 = 2 * (p - m); y2 = 3 * (q - m * m - m * y1)
    return y1, y2

def y3_extremes(m, p, q):
    """(y3min, y3max, lower law, upper law) for Y on [A,B]=[-m,1-m] with moments y1, y2. Laws as (atoms, weights)."""
    A, B = -m, 1 - m
    y1, y2 = y12(m, p, q); var = y2 - y1 * y1
    # upper: atoms {a, B}
    a = y1 - var / (B - y1); wB = (y1 - a) / (B - a)
    up = ([a, B], [1 - wB, wB]); y3max = (1 - wB) * a ** 3 + wB * B ** 3
    # lower: atoms {A, b}
    b = y1 + var / (y1 - A); wA = (b - y1) / (b - A)
    lo = ([A, b], [wA, 1 - wA]); y3min = wA * A ** 3 + (1 - wA) * b ** 3
    feasible = (var >= 0) and (a >= A) and (b <= B)
    return y3min, y3max, lo, up, feasible

def theta_moment_from_Y(k, m, ymoms):
    """E theta^k with theta = m + U Y; ymoms[j] = E Y^j, j = 0..k."""
    return sum(comb(k, j) * m ** (k - j) * ymoms[j] / (j + 1) for j in range(k + 1))

def m3_bounds(m, p, q):
    y1, y2 = y12(m, p, q)
    y3min, y3max, lo, up, feas = y3_extremes(m, p, q)
    m3max = theta_moment_from_Y(3, m, [F(1), y1, y2, y3max]); m3min = theta_moment_from_Y(3, m, [F(1), y1, y2, y3min])
    return m3min, m3max, lo, up, feas

def mode_interval(p, q):
    """Feasible modes: y1 in [A,B] and 0 <= var <= (B-y1)(y1-A). Returns exact [m_lo, m_hi] (solved by sympy)."""
    m = sp.symbols('m')
    y1 = 2 * (sp.Rational(p) - m); y2 = 3 * (sp.Rational(q) - m ** 2 - m * y1); var = y2 - y1 ** 2
    A, B = -m, 1 - m
    conds = [y1 - A >= 0, B - y1 >= 0, var >= 0, (B - y1) * (y1 - A) - var >= 0, m >= 0, m <= 1]
    sol = sp.solve_univariate_inequality(sp.And(*conds), m, relational=False) if False else None
    # do it by hand: each condition is polynomial of degree <= 2; intersect exactly
    S = sp.Interval(0, 1)
    for cnd in conds:
        S = S.intersect(sp.solve_univariate_inequality(cnd, m, relational=False))
    return S

# ---------- exact Delta from (m3, m4) via the affine moment kernels ----------
def delta_exact(m3, m4, p=P, q=Qm, e0=E0, e1=E1, c=C):
    out = {}
    for act in ("R", "C"):
        base, Bx, Cy = V.affine_parts(act, 4, e0, e1, c, p, q)
        K = V.kernel_at(base, Bx, Cy, m3, m4)
        out[act] = V.persist(K, V.start_state(act, 4))
    return out["R"] - out["C"], out["R"], out["C"]

def uniform_mixture_moments(m, law, kmax=4):
    """theta = m + U Y with Y two-point law (atoms, weights): exact E theta^k, k<=kmax, by direct integration of the
    uniform pieces [m, m+a] (independent check of the Khinchin moment identity)."""
    atoms, weights = law; moms = []
    for k in range(kmax + 1):
        tot = F(0)
        for a, w in zip(atoms, weights):
            if a == 0: tot += w * m ** k
            else:
                lo, hi = sorted((m, m + a)); tot += w * (hi ** (k + 1) - lo ** (k + 1)) / ((k + 1) * (hi - lo))
        moms.append(tot)
    return moms

# ---------- exact certificate: m3max_uni(m) < x* on the feasible mode set ----------
def certify_max(p, q, xstar, modes):
    m = sp.symbols('m')
    p_, q_ = sp.Rational(p), sp.Rational(q)
    y1 = 2 * (p_ - m); y2 = 3 * (q_ - m ** 2 - m * y1); var = sp.expand(y2 - y1 ** 2)
    A, B = -m, 1 - m
    a = y1 - var / (B - y1); wB = sp.simplify((y1 - a) / (B - a))
    y3max = sp.simplify((1 - wB) * a ** 3 + wB * B ** 3)
    m3max = sp.simplify(m ** 3 + sp.Rational(3, 2) * m ** 2 * y1 + m * y2 + y3max / 4)
    gap = sp.together(sp.Rational(xstar) - m3max)
    num, den = sp.fraction(gap)
    num, den = sp.Poly(sp.expand(num), m), sp.Poly(sp.expand(den), m)
    lo, hi = modes.inf, modes.sup
    # sign of denominator and numerator on [lo, hi]: no roots inside and positive at an interior point
    mid = (lo + hi) / 2
    res = {"m3max_uni_expr": str(m3max), "gap_numerator": str(num.as_expr()), "gap_denominator": str(den.as_expr()),
           "modes": [str(lo), str(hi)]}
    for name, poly in (("numerator", num), ("denominator", den)):
        nroots = poly.count_roots(lo, hi)   # Sturm, exact, closed interval
        val = poly.eval(mid)
        res[name + "_roots_in_mode_interval"] = int(nroots); res[name + "_sign_at_midpoint"] = int(sp.sign(val))
        res[name + "_at_endpoints"] = [str(poly.eval(lo)), str(poly.eval(hi))]
    res["certified_gap_positive"] = (res["numerator_roots_in_mode_interval"] == 0 and res["denominator_roots_in_mode_interval"] == 0
                                     and res["numerator_sign_at_midpoint"] == res["denominator_sign_at_midpoint"])
    # location of the unimodal maximum: stationary points of m3max on [lo, hi]
    d = sp.together(sp.diff(m3max, m)); dn, dd = sp.fraction(d)
    dn = sp.Poly(sp.expand(dn), m)
    crit = [r for r in sp.real_roots(dn) if lo <= r <= hi]
    cands = [lo, hi] + crit
    vals = [(sp.nsimplify(cnd) if cnd.is_Rational else cnd, m3max.subs(m, cnd)) for cnd in cands]
    best = max(vals, key=lambda t: float(t[1]))
    res["stationary_modes_in_interval"] = [str(sp.nsimplify(r)) + f" ~ {float(r):.10f}" for r in crit]
    res["argmax_mode"] = str(best[0]) + f" ~ {float(best[0]):.10f}"; res["sup_m3_uni"] = f"{float(best[1]):.12f}"
    res["sup_m3_uni_exact_if_rational"] = str(best[1]) if best[1].is_Rational else "algebraic (see stationary_modes)"
    return res, m3max, best

def main():
    t0 = time.time(); out = {"p": str(P), "q": str(Qm), "e0": str(E0), "e1": str(E1), "c": str(C), "T": V.T,
                             "xstar_upper_curve_root_bracket_m3": [str(XSTAR), str(XSTAR_HI)],
                             "lower_curve_root_bracket_m3": [str(LOWER_ROOT[0]), str(LOWER_ROOT[1])],
                             "class_m3_range": [str(v) for v in V.xbounds(P, Qm)]}
    modes = mode_interval(P, Qm); out["feasible_modes"] = [str(modes.inf), str(modes.sup)]
    print(f"feasible modes for unimodal laws with (p,q)=({P},{Qm}): m in [{modes.inf}, {modes.sup}] = [{float(modes.inf):.4f}, {float(modes.sup):.4f}]")
    print(f"full conditional-iid class m3 range: [{float(V.xbounds(P,Qm)[0]):.6f}, {float(V.xbounds(P,Qm)[1]):.6f}]; ambiguity band (certified roots) m3 in ({float(XSTAR):.10f}, {float(LOWER_ROOT[1]):.10f})")
    # table over modes
    rows = []
    lo, hi = F(modes.inf.p, modes.inf.q), F(modes.sup.p, modes.sup.q)
    for i in range(0, 12):
        m = lo + (hi - lo) * F(i, 11)
        m3min, m3max, lawlo, lawup, feas = m3_bounds(m, P, Qm)
        rows.append({"mode": str(m), "mode_float": float(m), "m3min": float(m3min), "m3max": float(m3max), "feasible": feas,
                     "upper_law_Y": {"atoms": [str(v) for v in lawup[0]], "weights": [str(v) for v in lawup[1]]}})
        print(f"  m={float(m):.4f}  m3_uni in [{float(m3min):.6f}, {float(m3max):.6f}]  feasible={feas}   Y-law for max: atoms {[float(v) for v in lawup[0]]} weights {[float(v) for v in lawup[1]]}")
    out["mode_table"] = rows
    # exact certificate
    cert, m3max_expr, best = certify_max(P, Qm, XSTAR, modes); out["certificate"] = cert
    print("certificate:", json.dumps({k: v for k, v in cert.items() if k not in ("m3max_uni_expr",)}, indent=1))
    print("m3max_uni(m) =", cert["m3max_uni_expr"])
    # minimum of m3 over unimodal laws (for the full unimodal m3 range)
    m = sp.symbols('m')
    m3min_vals = [float(m3_bounds(lo + (hi - lo) * F(i, 400), P, Qm)[0]) for i in range(401)]
    out["inf_m3_uni_grid"] = min(m3min_vals); print(f"inf m3 over unimodal laws (401-point mode grid, float): {min(m3min_vals):.8f}  [full class min {float(V.xbounds(P,Qm)[0]):.8f}]")
    # extremal unimodal law at the maximising mode; exact Delta there and at the boundary-mode laws
    checks = []
    for label, mm in (("argmax_mode", best[0]), ("mode_0", sp.Integer(0)), ("mode_sup", modes.sup)):
        mmF = F(int(mm.p), int(mm.q)) if mm.is_Rational else None
        if mmF is None:
            continue
        m3min, m3max, lawlo, lawup, feas = m3_bounds(mmF, P, Qm)
        moms = uniform_mixture_moments(mmF, lawup)
        assert moms[1] == P and moms[2] == Qm and moms[3] == m3max, (moms, m3max)
        y4 = sum(w * a ** 4 for a, w in zip(*lawup)); m4 = theta_moment_from_Y(4, mmF, [F(1), *y12(mmF, P, Qm), lawup and (lawup[1][0] * lawup[0][0] ** 3 + lawup[1][1] * lawup[0][1] ** 3), y4])
        assert m4 == moms[4]
        d, pr, pc = delta_exact(m3max, m4)
        chk = {"label": label, "mode": str(mmF), "law_theta": f"mixture of Uniform[{float(min(mmF, mmF + lawup[0][0])):.6f},{float(max(mmF, mmF + lawup[0][0])):.6f}] (w={float(lawup[1][0]):.6f}) and Uniform[{float(mmF):.6f},1] (w={float(lawup[1][1]):.6f})",
               "m3": str(m3max), "m3_float": float(m3max), "m4": str(m4), "m4_float": float(m4), "m3_below_xstar": m3max < XSTAR,
               "Delta_float": float(d), "PR_float": float(pr), "PC_float": float(pc), "Delta_negative": d < 0,
               "m4_within_lens": V.Lfun(m3max, P, Qm) <= m4 <= V.Ufun(m3max, P, Qm)}
        checks.append(chk); print("extremal unimodal law:", json.dumps(chk, indent=1))
    # law at a rational mode near the irrational argmax (exact Delta there)
    for mmF in (F(33, 500), F(66, 1000)):
        m3min, m3max, lawlo, lawup, feas = m3_bounds(mmF, P, Qm); moms = uniform_mixture_moments(mmF, lawup)
        assert moms[3] == m3max; d, pr, pc = delta_exact(moms[3], moms[4])
        checks.append({"label": "near_argmax_mode", "mode": str(mmF), "m3": str(moms[3]), "m3_float": float(moms[3]), "m4": str(moms[4]), "m4_float": float(moms[4]),
                       "law_theta": f"mixture of Uniform[{float(min(mmF, mmF + lawup[0][0])):.6f},{float(max(mmF, mmF + lawup[0][0])):.6f}] (w={float(lawup[1][0]):.6f}) and Uniform[{float(mmF):.6f},1] (w={float(lawup[1][1]):.6f})",
                       "Delta_float": float(d), "PR_float": float(pr), "PC_float": float(pc), "Delta_negative": d < 0, "m3_below_xstar": moms[3] < XSTAR})
        print("near-argmax unimodal law:", json.dumps(checks[-1], indent=1)); break
    out["extremal_law_checks"] = checks
    # rigorous upper bound on Delta over the unimodal class: Delta(x, U(x)) at rational x >= sup m3_uni
    xb = F(79, 4000); assert float(xb) > float(best[1]) and xb < XSTAR
    dU, prU, pcU = delta_exact(xb, V.Ufun(xb, P, Qm))
    out["unimodal_class_Delta_upper_bound"] = {"x": str(xb), "x_float": float(xb), "m4_on_upper_curve": str(V.Ufun(xb, P, Qm)), "Delta_float": float(dU), "PR_float": float(prU), "PC_float": float(pcU),
                                              "argument": "Delta increasing in m4 over lens (certified) and increasing in m3 along the upper curve (certified); sup m3_uni < x < x*"}
    print(f"unimodal-class upper bound on Delta: Delta(G) <= Delta({xb}, U({xb})) = {float(dU):+.8f} for every unimodal G  (PR {float(prU):.6f}, PC {float(pcU):.6f})")
    # exploratory float search for sup Delta over unimodal laws: over (mode m, y3), with y4 at its Hausdorff maximum
    # given (y1,y2,y3) on [A,B] (upper principal representation, atoms {A, a, B}: a = E[(Y-A)(Y-B)Y]/E[(Y-A)(Y-B)],
    # y4 = (A+B+a) y3 - (AB + a(A+B)) y2 + AB a y1). Delta is increasing in m4, so the sup lies on this boundary.
    import numpy as np
    parts = {act: [np.array([[float(v) for v in row] for row in M]) for M in V.affine_parts(act, 4, E0, E1, C, P, Qm)] for act in ("R", "C")}
    def delta_float(x, y):
        d = {}
        for act in ("R", "C"):
            base, Bx, Cy = parts[act]; K = base + x * Bx + y * Cy
            v = np.zeros(len(K)); v[V.start_state(act, 4)] = 1.0
            for _ in range(V.T): v = v @ K
            d[act] = v[1:].sum()
        return d["R"] - d["C"]
    bestD = (-1.0, None); pf, qf = float(P), float(Qm)
    for mm in np.linspace(0.0, float(hi), 111):
        A, B = -mm, 1 - mm; y1 = 2 * (pf - mm); y2 = 3 * (qf - mm * mm - mm * y1); var = y2 - y1 * y1
        a = y1 - var / (B - y1); wB = (y1 - a) / (B - a); y3max = (1 - wB) * a ** 3 + wB * B ** 3
        b = y1 + var / (y1 - A); wA = (b - y1) / (b - A); y3min = wA * A ** 3 + (1 - wA) * b ** 3
        for y3 in np.linspace(y3min, y3max, 81):
            den = y2 - (A + B) * y1 + A * B
            aa = (y3 - (A + B) * y2 + A * B * y1) / den
            y4 = (A + B + aa) * y3 - (A * B + aa * (A + B)) * y2 + A * B * aa * y1
            ym = [1.0, y1, y2, y3, y4]
            x3 = sum(comb(3, j) * mm ** (3 - j) * ym[j] / (j + 1) for j in range(4)); x4 = sum(comb(4, j) * mm ** (4 - j) * ym[j] / (j + 1) for j in range(5))
            dd = delta_float(x3, x4)
            if dd > bestD[0]: bestD = (dd, {"mode": float(mm), "y3": float(y3), "m3": float(x3), "m4": float(x4), "interior_atom_a": float(aa)})
    out["unimodal_sup_Delta_float_search"] = {"grid": "111 modes x 81 y3 values, y4 at Hausdorff max", "sup_Delta_found": bestD[0], "at": bestD[1],
                                              "note": "exploratory float grid; a lower estimate of the true unimodal sup, which the certified bound above exceeds"}
    print(f"float search: largest Delta over unimodal laws found {bestD[0]:+.8f} at {bestD[1]}")
    # full-class comparison at the same setting
    out["full_class_Delta_range"] = {"min": -0.1205203070043943, "max": 0.03209273498500964, "source": "n4_boundary_results.json (certified, width <= 1e-9)"}
    # random sanity check of the y3 upper/lower bound on an interval (float)
    rng = random.Random(7); viol = 0; ntest = 20000; worst = 0.0
    for _ in range(ntest):
        mm = rng.uniform(0, 0.99); A, B = -mm, 1 - mm
        k = rng.randint(1, 6); atoms = [rng.uniform(A, B) for _ in range(k)]; w = [rng.random() for _ in range(k)]; s = sum(w); w = [x / s for x in w]
        y1 = sum(a * x for a, x in zip(atoms, w)); y2 = sum(a * a * x for a, x in zip(atoms, w)); y3 = sum(a ** 3 * x for a, x in zip(atoms, w)); var = y2 - y1 * y1
        if var <= 1e-12: continue
        a = y1 - var / (B - y1); wB = (y1 - a) / (B - a); y3max = (1 - wB) * a ** 3 + wB * B ** 3
        b = y1 + var / (y1 - A); wA = (b - y1) / (b - A); y3min = wA * A ** 3 + (1 - wA) * b ** 3
        if y3 > y3max + 1e-12 or y3 < y3min - 1e-12: viol += 1
        worst = max(worst, y3 - y3max, y3min - y3)
    out["y3_bound_random_check"] = {"tests": ntest, "violations": viol, "worst_excess": worst}
    print(f"y3 bound random check: {viol} violations in {ntest} random discrete Y laws (worst excess {worst:.2e})")
    # Khinchin identity check: random mixture of uniforms sharing endpoint m vs formula, exact
    bad = 0
    for _ in range(200):
        mm = F(rng.randint(0, 100), 100); k = rng.randint(1, 4)
        atoms = [F(rng.randint(-int(mm * 100), int((1 - mm) * 100)), 100) for _ in range(k)]; w = [F(rng.randint(1, 9)) for _ in range(k)]; s = sum(w); w = [x / s for x in w]
        moms = uniform_mixture_moments(mm, (atoms, w)); ym = [sum(x * a ** j for a, x in zip(atoms, w)) for j in range(5)]
        if any(moms[kk] != theta_moment_from_Y(kk, mm, ym) for kk in range(5)): bad += 1
    out["khinchin_identity_exact_check_failures"] = bad; print(f"Khinchin moment identity: {bad} failures in 200 exact random uniform mixtures")
    # reference: unimodal m3 range for q = 6/100
    modes6 = mode_interval(P, F(6, 100))
    if modes6.measure > 0 or not modes6.is_EmptySet:
        lo6, hi6 = F(modes6.inf.p, modes6.inf.q), F(modes6.sup.p, modes6.sup.q)
        vals = [m3_bounds(lo6 + (hi6 - lo6) * F(i, 400), P, F(6, 100))[:2] for i in range(401)]
        out["q6_reference"] = {"feasible_modes": [str(lo6), str(hi6)], "inf_m3_uni_grid": float(min(v[0] for v in vals)), "sup_m3_uni_grid": float(max(v[1] for v in vals)),
                               "class_m3_range": [float(v) for v in V.xbounds(P, F(6, 100))]}
        print("q=6/100 reference:", out["q6_reference"])
    else:
        out["q6_reference"] = "no unimodal law has (p,q) = (1/10, 6/100)"; print(out["q6_reference"])
    out["seconds"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(os.path.dirname(__file__), "unimodal_check.json"), "w"), indent=1)
    print(f"done in {out['seconds']}s")

if __name__ == "__main__":
    main()
