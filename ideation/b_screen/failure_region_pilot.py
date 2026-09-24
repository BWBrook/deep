"""Failure-region pilot (Astra's contract, synthesis/failure_region_contract.md): eight new rational (e1, c) points
at n=4, p=1/10, q=3/100, e0=1/50, T=50, scalar class S. Points were chosen from the endpoint-pair scan
(n4_reversal_map.json, exploratory) to straddle the two candidate transition curves and to probe the interior, and
are saved to failure_region_pilot_points.json BEFORE any exact evaluation.

At each point, exactly (gmpy2 rationals):
  (A) attained witnesses: Delta at Glo, Ghi and the lens midpoint;
  (B) direct 2-D sign certificate of Delta over the lens (Bernstein tensor of degree (T,T), branch and bound with
      lens pruning; needs NO monotonicity): uniform sign -> identified; mixed leaves + opposite-sign witnesses -> unidentified;
  (C) vertical certificate dDelta/dm4 > 0 over the lens (n4_vertical.certify_dw);
  (D) boundary curves: degree-100 Bernstein coefficients of Delta on the lower and upper curves; exact min/max
      enclosures (n4_boundary.enclose, tol 1e-9); sign pieces; and whether Delta is strictly monotone along each curve
      (all derivative coefficients of one sign);
  (E) classification from (C)+(D) when (C) holds: a = min on lower curve, b = max on upper curve; and the check
      whether the "both roots exist" shortcut agrees with the classification.
Costs are recorded per stage. Fable, 2026-09-24.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
import n4_vertical as V
import n4_boundary as Bd
import n3_certificate as N3
from gmpy2 import mpq as Q

n = 4; P = V.P; q = F(3, 100); e0 = F(1, 50); T = V.T
POINTS = [("P1", F(1, 2), F(9, 50)), ("P2", F(1, 2), F(11, 50)), ("P3", F(7, 10), F(1, 10)), ("P4", F(7, 10), F(3, 20)),
          ("P5", F(9, 10), F(3, 50)), ("P6", F(17, 20), F(1, 50)), ("P7", F(1, 2), F(2, 5)), ("P8", F(19, 20), F(9, 20))]
HERE = os.path.dirname(__file__)

def sign_bb(Jd, xlo, xhi, ylo, yhi, max_depth=10):
    """Direct sign certificate of the tensor Bernstein polynomial Jd over the lens (no derivative)."""
    return V.certify_dw(Jd, xlo, xhi, ylo, yhi, P, q, max_depth=max_depth)

def curve_delta(parts, xL, xU, yfun):
    cR = Bd.objective_coeffs(*parts["R"], V.start_state("R", n), xL, xU, yfun, P, q)
    cC = Bd.objective_coeffs(*parts["C"], V.start_state("C", n), xL, xU, yfun, P, q)
    return [r - s for r, s in zip(cR, cC)]

def monotone(coeffs):
    m = len(coeffs) - 1; d = [m * (coeffs[k + 1] - coeffs[k]) for k in range(m)]
    return "increasing" if all(v > 0 for v in d) else "decreasing" if all(v < 0 for v in d) else "not certified monotone"

def main():
    json.dump({"n": n, "p": str(P), "q": str(q), "e0": str(e0), "T": T, "points": [{"id": i, "e1": str(a), "c": str(b)} for i, a, b in POINTS],
               "selection": "endpoint-pair scan n4_reversal_map.json (q=.03,e0=.02): P1/P2 straddle the diagonal boundary at e1=.5 (X from c=.20); P3/P4 at e1=.7 (X from .13); P5 '.' side at e1=.9 (X from .08); P6 isolated marked cell at c=.02 (suspected cancellation artefact); P7 high-c interior; P8 top corner of the contract domain"},
              open(os.path.join(HERE, "failure_region_pilot_points.json"), "w"), indent=1)
    xL, xU = V.xbounds(P, q); ylo, yhi = V.lens_range(xL, xU, P, q)
    (tl, wl), (tu, wu) = N3.laws(P, q)
    results = []
    for pid, e1, c in POINTS:
        t0 = time.time(); rec = {"id": pid, "e1": str(e1), "c": str(c), "cost_s": {}}
        parts = {a: V.affine_parts(a, n, e0, e1, c, P, q) for a in ("R", "C")}
        def delta_at(x, y):
            return V.persist(V.kernel_at(*parts["R"], x, y), V.start_state("R", n)) - V.persist(V.kernel_at(*parts["C"], x, y), V.start_state("C", n))
        # (A) witnesses
        t = time.time()
        wit = {}
        for nm, (atoms, wts) in (("Glo", (tl, wl)), ("Ghi", (tu, wu))):
            x = sum(w * a ** 3 for a, w in zip(atoms, wts)); y = sum(w * a ** 4 for a, w in zip(atoms, wts)); wit[nm] = delta_at(x, y)
        xm = (xL + xU) / 2; wit["mid"] = delta_at(xm, (V.Lfun(xm, P, q) + V.Ufun(xm, P, q)) / 2)
        rec["witnesses"] = {k: {"exact": str(v), "float": float(v)} for k, v in wit.items()}; rec["cost_s"]["witnesses"] = round(time.time() - t, 1)
        signs = {1 if v > 0 else -1 if v < 0 else 0 for v in wit.values()}
        # (B) direct 2-D sign certificate
        t = time.time()
        JR = V.bernstein_J(*parts["R"], V.start_state("R", n), xL, xU, ylo, yhi); JC = V.bernstein_J(*parts["C"], V.start_state("C", n), xL, xU, ylo, yhi)
        Jd = [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(JR, JC)]
        leaves, unres = sign_bb(Jd, xL, xU, ylo, yhi)
        if leaves["unresolved"] == 0 and (leaves["pos"] == 0 or leaves["neg"] == 0):
            direct = "identified: " + ("R" if leaves["pos"] else "C")
        elif len(signs - {0}) == 2:
            direct = "unidentified (opposite-sign exact witnesses)"
        else:
            direct = "unresolved"
        rec["direct_sign_certificate"] = {"leaves": leaves, "classification": direct}; rec["cost_s"]["direct_sign"] = round(time.time() - t, 1)
        # (C) vertical certificate
        t = time.time()
        D = V.dw_coeffs(Jd); vleaves, vunres = V.certify_dw(D, xL, xU, ylo, yhi, P, q, max_depth=10)
        vert = ("certified " + ("pos" if vleaves["pos"] and not vleaves["neg"] else "neg" if vleaves["neg"] and not vleaves["pos"] else "MIXED")) if vleaves["unresolved"] == 0 else "unresolved"
        rec["vertical"] = {"status": vert, "leaves": vleaves}; rec["cost_s"]["vertical"] = round(time.time() - t, 1)
        # (D) boundary curves
        t = time.time(); curves = {}
        for cname, yfun in (("lower", V.Lfun), ("upper", V.Ufun)):
            cd = curve_delta(parts, xL, xU, yfun)
            enc = Bd.enclose(cd); pieces = Bd.sign_certificate(cd)
            curves[cname] = {"min_attained": str(enc["min"]["attained"]), "min_outer": str(enc["min"]["outer"]), "max_attained": str(enc["max"]["attained"]), "max_outer": str(enc["max"]["outer"]),
                             "min_float": float(enc["min"]["outer"]), "max_float": float(enc["max"]["outer"]), "tolerance_met": enc["min"]["tolerance_met"] and enc["max"]["tolerance_met"],
                             "monotone": monotone(cd), "sign_pieces": [[str(a), str(b), s] for a, b, s in pieces], "roots_x": [[float(xL + a * (xU - xL)), float(xL + b * (xU - xL))] for a, b, s in pieces if s == "root"],
                             "n_roots": sum(1 for _, _, s in pieces if s == "root"), "unresolved_pieces": sum(1 for _, _, s in pieces if s == "unresolved")}
        rec["curves"] = curves; rec["cost_s"]["curves"] = round(time.time() - t, 1)
        # (E) classification via shortcut
        if vert == "certified pos":
            a_out, b_out = Q(curves["lower"]["min_outer"]), Q(curves["upper"]["max_outer"]); a_att, b_att = Q(curves["lower"]["min_attained"]), Q(curves["upper"]["max_attained"])
            if b_out < 0: sc = "identified: C"
            elif a_out > 0: sc = "identified: R"
            elif a_att < 0 < b_att: sc = "unidentified"
            else: sc = "unresolved (enclosure straddles 0 without attained opposite signs)"
            rec["shortcut_classification"] = sc
            rec["class_Delta_range"] = [float(a_out), float(b_out)]
            rec["both_roots_exist"] = curves["lower"]["n_roots"] >= 1 and curves["upper"]["n_roots"] >= 1
            rec["both_curves_increasing"] = curves["lower"]["monotone"] == "increasing" and curves["upper"]["monotone"] == "increasing"
        else:
            rec["shortcut_classification"] = "not applicable (vertical certificate failed)"; rec["both_roots_exist"] = None; rec["both_curves_increasing"] = None
        rec["cost_s"]["total"] = round(time.time() - t0, 1)
        results.append(rec)
        print(f"{pid} e1={e1} c={c}: Glo {float(wit['Glo']):+.5f} Ghi {float(wit['Ghi']):+.5f} | direct: {direct} | vertical: {vert} | lower Δ[{curves['lower']['min_float']:+.5f},{curves['lower']['max_float']:+.5f}] {curves['lower']['monotone']}, roots {curves['lower']['n_roots']} | upper Δ[{curves['upper']['min_float']:+.5f},{curves['upper']['max_float']:+.5f}] {curves['upper']['monotone']}, roots {curves['upper']['n_roots']} | shortcut: {rec['shortcut_classification']} | cost {rec['cost_s']}", flush=True)
        json.dump({"n": n, "p": str(P), "q": str(q), "e0": str(e0), "T": T, "lens": [str(xL), str(xU), str(ylo), str(yhi)], "results": results}, open(os.path.join(HERE, "failure_region_pilot.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
