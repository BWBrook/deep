"""Exact whole-class certificates over the finite-exchangeable class E(p,q) at n=4, for the original 16-case design
(Astra's contract, synthesis/exchangeable_certificate_contract.md).

E(p,q) in coordinates (x, y) = (tau3, tau4) is the polygon {pi_k(x,y) >= 0, k=0..4} with
  pi0 = 1-4p+6q-4x+y, pi1 = 4p-12q+12x-4y, pi2 = 6q-12x+6y, pi3 = 4x-4y, pi4 = y.
The aggregated kernels are K = base + x*A3 + y*A4 on all of E (exchangeable_kernel_audit.md), the same A3, A4 as on the
scalar lens S, so Delta(x, y) is the same bivariate function extended to the larger domain.

Per case: (1) exact polygon vertices and edges; (2) tensor Bernstein of Delta on E's own bounding box, derivative in y,
sign certificate by branch and bound with POLYGON pruning (a box is discarded only if some affine pi_k is strictly
negative on all of it); (3) if dDelta/dy > 0 on E, a_E = min of Delta on the lower vertex chain and b_E = max on the
upper chain, each edge an exact degree-2T Bernstein polynomial (affine kernel along the edge; validated against direct
exact propagation at u = 1/4, 3/4), enclosed to 1e-9 with attained witnesses; monotonicity along each edge recorded but
not relied upon; (4) classification, regret, and comparison with the scalar-class certificate (n4_boundary_results.json).
Construction and proof timings are separated. gmpy2 rationals throughout. Fable, 2026-09-24.
"""
import sys, os, json, time, itertools, hashlib, platform
sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
import gmpy2
from gmpy2 import mpq as Q
import n4_vertical as V
import n4_boundary as Bd
import n3_certificate as N3

n = 4; T = V.T; P = V.P; HERE = os.path.dirname(__file__)
CASES = [(q, e0, e1, c) for q in V.QS for e0, e1, c in itertools.product(V.E0S, V.E1S, V.CS)]
def q_(f): return Q(f.numerator, f.denominator) if isinstance(f, F) else Q(f)

# ---------- polygon ----------
def pis(p, q):
    """Affine constraints pi_k(x,y) = a + b x + c y as (a, b, c)."""
    return [(1 - 4 * p + 6 * q, -4, 1), (4 * p - 12 * q, 12, -4), (6 * q, -12, 6), (0, 4, -4), (0, 0, 1)]

def polygon(p, q):
    cons = pis(p, q); pts = set()
    for (a1, b1, c1), (a2, b2, c2) in itertools.combinations(cons, 2):
        det = b1 * c2 - b2 * c1
        if det == 0: continue
        x = (-a1 * c2 + a2 * c1) / det; y = (-b1 * a2 + b2 * a1) / det
        if all(a + b * x + c * y >= 0 for a, b, c in cons): pts.add((F(x), F(y)))
    pts = sorted(pts)
    # CCW order around centroid
    cx = sum(x for x, _ in pts) / len(pts); cy = sum(y for _, y in pts) / len(pts)
    import math
    pts.sort(key=lambda t: math.atan2(float(t[1] - cy), float(t[0] - cx)))
    # rotate so that the leftmost (then lowest) vertex is first
    i0 = min(range(len(pts)), key=lambda i: (pts[i][0], pts[i][1])); pts = pts[i0:] + pts[:i0]
    # drop collinear middle vertices
    keep = []
    for i in range(len(pts)):
        a, b, c = pts[i - 1], pts[i], pts[(i + 1) % len(pts)]
        cross = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        if cross != 0: keep.append(b)
    return keep

def chains(verts):
    """Lower chain: from leftmost vertex CCW to the rightmost; upper chain: the rest (rightmost back to leftmost)."""
    ir = max(range(len(verts)), key=lambda i: (verts[i][0], -verts[i][1]))
    lower = [(verts[i], verts[i + 1]) for i in range(ir)]
    upper = [(verts[i], verts[(i + 1) % len(verts)]) for i in range(ir, len(verts))]
    return lower, upper

def box_excluded(cons, x0, x1, y0, y1):
    for a, b, c in cons:
        mx = a + max(b * x0, b * x1) + max(c * y0, c * y1)
        if mx < 0: return True
    return False

# ---------- vertical certificate over the polygon ----------
def certify_dy_polygon(D, cons, xlo, xhi, ylo, yhi, max_depth=10):
    stack = [(D, F(0), F(1), F(0), F(1), 0)]; leaves = {"pos": 0, "neg": 0, "pruned": 0, "unresolved": 0}; unres = []; rootbox = None
    while stack:
        arr, u0, u1, w0, w1, depth = stack.pop()
        x0, x1 = xlo + u0 * (xhi - xlo), xlo + u1 * (xhi - xlo); y0, y1 = ylo + w0 * (yhi - ylo), ylo + w1 * (yhi - ylo)
        if box_excluded(cons, x0, x1, y0, y1): leaves["pruned"] += 1; continue
        s = V.sign_of(arr)
        if depth == 0: rootbox = {"min_coeff": str(min(min(r) for r in arr)), "max_coeff": str(max(max(r) for r in arr)), "uniform_sign": s}
        if s in ("pos", "neg"): leaves[s] += 1; continue
        if depth >= max_depth: leaves["unresolved"] += 1; unres.append([float(x0), float(x1), float(y0), float(y1)]); continue
        axis = 0 if (u1 - u0) >= (w1 - w0) else 1
        a, b = V.de_casteljau_split(arr, axis)
        if axis == 0: um = (u0 + u1) / 2; stack += [(a, u0, um, w0, w1, depth + 1), (b, um, u1, w0, w1, depth + 1)]
        else: wm = (w0 + w1) / 2; stack += [(a, u0, u1, w0, wm, depth + 1), (b, u0, u1, wm, w1, depth + 1)]
    return leaves, unres, rootbox

# ---------- edge polynomials ----------
def edge_delta(parts, p0, p1):
    def Kq(action, pt):
        base, Bx, Cy = parts[action]; K = V.kernel_at(base, Bx, Cy, pt[0], pt[1])
        return [[q_(v) for v in row] for row in K]
    out = {}
    for action in ("R", "C"):
        K0, K2 = Kq(action, p0), Kq(action, p1); K1 = [[(a + b) / 2 for a, b in zip(ra, rb)] for ra, rb in zip(K0, K2)]
        S = len(K0); v = [[Q(0)] * S]; v[0][V.start_state(action, n)] = Q(1)
        for _ in range(T): v = Bd.bern_mul_deg2(v, K0, K1, K2)
        out[action] = [sum(row) - row[0] for row in v]
    return [r - s for r, s in zip(out["R"], out["C"])]

def monotone(coeffs):
    m = len(coeffs) - 1; d = [m * (coeffs[k + 1] - coeffs[k]) for k in range(m)]
    return "increasing" if all(v > 0 for v in d) else "decreasing" if all(v < 0 for v in d) else "not certified monotone"

def delta_direct(parts, x, y):
    return V.persist(V.kernel_at(*parts["R"], x, y), V.start_state("R", n)) - V.persist(V.kernel_at(*parts["C"], x, y), V.start_state("C", n))

def run_case(case):
    q, e0, e1, c = case; t0 = time.time(); rec = {"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), "cost_s": {}}
    t = time.time(); parts = {a: V.affine_parts(a, n, e0, e1, c, P, q) for a in ("R", "C")}; rec["cost_s"]["kernels"] = round(time.time() - t, 2)
    cons = pis(P, q); verts = polygon(P, q); lower, upper = chains(verts)
    rec["polygon_vertices"] = [[str(x), str(y)] for x, y in verts]
    rec["vertex_Delta"] = {f"({x},{y})": {"exact": str(delta_direct(parts, x, y)), "float": float(delta_direct(parts, x, y))} for x, y in verts}
    xlo, xhi = min(x for x, _ in verts), max(x for x, _ in verts); ylo, yhi = min(y for _, y in verts), max(y for _, y in verts)
    rec["bounding_box"] = [str(xlo), str(xhi), str(ylo), str(yhi)]
    # tensor and vertical certificate
    t = time.time(); JR = V.bernstein_J(*parts["R"], V.start_state("R", n), xlo, xhi, ylo, yhi); JC = V.bernstein_J(*parts["C"], V.start_state("C", n), xlo, xhi, ylo, yhi)
    Jd = [[a - b for a, b in zip(ra, rb)] for ra, rb in zip(JR, JC)]; rec["cost_s"]["tensor_construction"] = round(time.time() - t, 1)
    t = time.time(); D = V.dw_coeffs(Jd); leaves, unres, rootbox = certify_dy_polygon(D, cons, xlo, xhi, ylo, yhi); rec["cost_s"]["derivative_and_sign"] = round(time.time() - t, 1)
    vert = ("certified " + ("pos" if leaves["pos"] and not leaves["neg"] else "neg" if leaves["neg"] and not leaves["pos"] else "MIXED")) if leaves["unresolved"] == 0 else "unresolved"
    rec["vertical_on_E"] = {"status": vert, "leaves": leaves, "unresolved_boxes": unres[:10], "rootbox": rootbox}
    # direct sign certificate of Delta on E (informational; classification below uses extrema)
    t = time.time(); sleaves, sunres, _ = certify_dy_polygon(Jd, cons, xlo, xhi, ylo, yhi); rec["cost_s"]["direct_sign"] = round(time.time() - t, 1)
    rec["direct_sign_on_E"] = {"leaves": sleaves, "unresolved_boxes": sunres[:10]}
    # edges
    t = time.time(); edges = {}
    for cname, chain in (("lower", lower), ("upper", upper)):
        for (p0, p1) in chain:
            cd = edge_delta(parts, p0, p1)
            for u in (F(1, 4), F(3, 4)):
                x = p0[0] + u * (p1[0] - p0[0]); y = p0[1] + u * (p1[1] - p0[1])
                assert Bd.evaluate(cd, u) == q_(delta_direct(parts, x, y)), ("edge validation failed", cname, p0, p1, u)
            enc = Bd.enclose(cd); pieces = Bd.sign_certificate(cd)
            edges[f"{cname} ({p0[0]},{p0[1]})->({p1[0]},{p1[1]})"] = {
                "chain": cname, "from": [str(p0[0]), str(p0[1])], "to": [str(p1[0]), str(p1[1])], "degree": len(cd) - 1,
                "min_attained": str(enc["min"]["attained"]), "min_outer": str(enc["min"]["outer"]), "max_attained": str(enc["max"]["attained"]), "max_outer": str(enc["max"]["outer"]),
                "min_witness_u": str(enc["min"]["witness_u"]), "max_witness_u": str(enc["max"]["witness_u"]),
                "tolerance_met": enc["min"]["tolerance_met"] and enc["max"]["tolerance_met"], "monotone": monotone(cd),
                "sign_pieces": [[str(a), str(b), s] for a, b, s in pieces], "validated_u": ["1/4", "3/4"]}
    rec["cost_s"]["edges"] = round(time.time() - t, 1); rec["edges"] = edges
    if vert == "certified pos":
        lo_edges = [e for e in edges.values() if e["chain"] == "lower"]; up_edges = [e for e in edges.values() if e["chain"] == "upper"]
        a_out = min(Q(e["min_outer"]) for e in lo_edges); a_att = min(Q(e["min_attained"]) for e in lo_edges)
        b_out = max(Q(e["max_outer"]) for e in up_edges); b_att = max(Q(e["max_attained"]) for e in up_edges)
        tol_ok = all(e["tolerance_met"] for e in edges.values())
        if b_out < 0: cls = "identified: C"
        elif a_out > 0: cls = "identified: R"
        elif a_att < 0 < b_att: cls = "unidentified"
        else: cls = "unresolved"
        regret = min(-a_att, b_att) if cls == "unidentified" else Q(0)
        lottery = (-a_att * b_att / (b_att - a_att)) if cls == "unidentified" else Q(0)
        rec["E_class"] = {"a_E": {"outer": str(a_out), "attained": str(a_att), "float_outer": Bd._float_outward(a_out, "down")}, "b_E": {"outer": str(b_out), "attained": str(b_att), "float_outer": Bd._float_outward(b_out, "up")},
                          "classification": cls, "regret_attained": str(regret), "regret_float": float(regret), "lottery_regret_float": float(lottery), "all_tolerances_met": tol_ok}
    else:
        rec["E_class"] = {"classification": "unresolved (vertical certificate on E not obtained)"}
    rec["cost_s"]["total"] = round(time.time() - t0, 1)
    return rec

def main():
    t0 = time.time()
    S = json.load(open(os.path.join(HERE, "n4_boundary_results.json")))
    scal = {(cs["q"], cs["e0"], cs["e1"], cs["c"]): cs for cs in S["cases"]}
    order = sorted(CASES, key=lambda cs: (cs[1], cs[2], cs[3]) != (F(1, 50), F(9, 10), F(1, 10)))   # the two reference cases first
    import multiprocessing as mp
    with mp.get_context("fork").Pool(8) as pool:
        recs = pool.map(run_case, order)
    out = {"n": n, "T": T, "p": str(P), "metadata": {"python": platform.python_version(), "gmpy2": gmpy2.version(),
           "source_sha256": {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in ("exchangeable_certificate.py", "n4_vertical.py", "n4_boundary.py", "n3_certificate.py", "n4_boundary_results.json")}}, "cases": []}
    print(f"{'q':>5} {'e0':>4} {'e1':>4} {'c':>4} | S class Delta range, S class | E vertical | E class Delta range (outer), E class, regret | costs")
    for rec in recs:
        cs = scal[(rec["q"], rec["e0"], rec["e1"], rec["c"])]; smin, smax = cs["whole_class"]["Delta"]["min"], cs["whole_class"]["Delta"]["max"]
        scls = "C" if smax < 0 else "R" if smin > 0 else "unidentified"
        rec["S_class"] = {"Delta_min": smin, "Delta_max": smax, "classification": scls, "regret": (min(-smin, smax) if scls == "unidentified" else 0.0)}
        E = rec["E_class"]; out["cases"].append(rec)
        es = f"[{E['a_E']['float_outer']:+.5f},{E['b_E']['float_outer']:+.5f}] {E['classification']:<16} regret {E['regret_float']:.5f}" if "a_E" in E else E["classification"]
        print(f"{rec['q']:>5} {rec['e0']:>4} {rec['e1']:>4} {rec['c']:>4} | [{smin:+.5f},{smax:+.5f}] {scls:<12} | {rec['vertical_on_E']['status']:<14} | {es} | {rec['cost_s']}", flush=True)
    out["seconds_wall"] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(HERE, "exchangeable_certificate.json"), "w"), indent=1)
    print(f"wall {out['seconds_wall']}s")

if __name__ == "__main__":
    main()
