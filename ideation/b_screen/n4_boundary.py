"""n=4, second stage: exact univariate Bernstein certificates on the two boundary curves of the (m3, m4) feasible set.

Given (n4_vertical.py): dJ/dy < 0 for J in {P0, PR, PC} and dDelta/dy > 0 over the whole feasible set, for all
16 cases. Hence over the class: max P_a lies on the lower curve y = L(x), min P_a on the upper curve y = U(x);
max Delta on the upper curve, min Delta on the lower curve. Each curve gives a kernel quadratic in x, so the
T-step objective is a univariate polynomial of degree <= 2T = 100 in u, x = xL + u (xU - xL), u in [0, 1].

Outputs per case: sign certificates for Delta (uniform sign, or isolated root brackets via monotone pieces
after exact subdivision); exact attained/outer enclosures for min/max of each objective on each curve to width
<= 1e-9; the resulting whole-class bounds. Bernstein coefficients are retained internally rather than
exported. Exact rationals throughout (gmpy2.mpq). Cross-check: Bernstein evaluation at u = 1/4, 3/4 equals
direct exact propagation with the affine kernel at (x, y(x)).
"""
import sys, json, time, itertools, hashlib, platform
from fractions import Fraction as F
from math import comb
from math import nextafter
sys.path.insert(0, __import__("os").path.dirname(__file__))
import n4_vertical as V
import gmpy2
from gmpy2 import mpq as Q

T = V.T; P = V.P
def q_(f): return Q(f.numerator, f.denominator)

def bern_mul_deg2(v, K0, K1, K2):
    """v: list (len m+1) of state row-vectors (Bernstein degree m). Returns degree m+2 coefficients of v(u) K(u),
    K(u) = B02 K0 + B12 K1 + B22 K2. Uses B_{i,m} B_{k,2} = C(m,i)C(2,k)/C(m+2,i+k) B_{i+k,m+2}."""
    m = len(v) - 1; S = len(K0)
    Ks = (K0, K1, K2); Kt = [list(zip(*K)) for K in Ks]
    out = [[Q(0)] * S for _ in range(m + 3)]
    for i in range(m + 1):
        for k in range(3):
            coef = Q(comb(m, i) * comb(2, k), comb(m + 2, i + k))
            row = v[i]; col = Kt[k]
            tgt = out[i + k]
            for s in range(S):
                tgt[s] += coef * sum(row[r] * col[s][r] for r in range(S))
    return out

def curve_kernels(base, Bx, Cy, xL, xU, yfun, p, q):
    """Bernstein-degree-2 control matrices of K(u) = base + x(u) B + y(x(u)) C on u in [0,1]."""
    def Kat(u):
        x = xL + u * (xU - xL); y = yfun(x, p, q)
        return [[q_(base[r][s] + x * Bx[r][s] + y * Cy[r][s]) for s in range(len(base))] for r in range(len(base))]
    K0, Kh, K2 = Kat(F(0)), Kat(F(1, 2)), Kat(F(1))
    K1 = [[2 * Kh[r][s] - (K0[r][s] + K2[r][s]) / 2 for s in range(len(K0))] for r in range(len(K0))]
    return K0, K1, K2

def objective_coeffs(base, Bx, Cy, start, xL, xU, yfun, p, q):
    K0, K1, K2 = curve_kernels(base, Bx, Cy, xL, xU, yfun, p, q)
    S = len(K0)
    v = [[Q(0)] * S]; v[0][start] = Q(1)
    for _ in range(T):
        v = bern_mul_deg2(v, K0, K1, K2)
    return [sum(row) - row[0] for row in v]   # degree 2T Bernstein coefficients of persistence

def evaluate(coeffs, u):
    n = len(coeffs) - 1; u = Q(u.numerator, u.denominator) if isinstance(u, F) else Q(u)
    return sum(c * comb(n, k) * u ** k * (1 - u) ** (n - k) for k, c in enumerate(coeffs))

def split(coeffs):
    n = len(coeffs) - 1; tri = [list(coeffs)]
    for k in range(1, n + 1):
        prev = tri[-1]; tri.append([(prev[m] + prev[m + 1]) / 2 for m in range(n - k + 1)])
    return [tri[k][0] for k in range(n + 1)], [tri[n - k][k] for k in range(n + 1)]

def _qstr(value):
    """Serialize an exact rational without passing through binary floating point."""
    return str(value)


def _float_outward(value, direction):
    """Convert an exact rational to a float, moving one ULP only when rounding is inward."""
    x = float(value)
    if direction == "down" and Q(x) > value:
        return nextafter(x, float("-inf"))
    if direction == "up" and Q(x) < value:
        return nextafter(x, float("inf"))
    return x


def enclose(coeffs, tol=F(1, 10**9), maxd=60):
    """Return exact attained/outer intervals for a Bernstein polynomial on ``[0, 1]``.

    Branch-and-bound coefficient ranges are mathematically valid outer bounds. Endpoint
    evaluations supply feasible attained values and their global ``u`` witnesses. A
    depth-limited search remains an explicit non-converged certificate rather than a
    silently successful result. The result is a dictionary of exact ``gmpy2.mpq``
    values; JSON serialization is handled by :func:`certificate_json`.
    """
    tolq = q_(tol)
    def bound(kind):
        if kind == "max":
            best = max(coeffs[0], coeffs[-1]); best_u = F(0) if coeffs[0] >= coeffs[-1] else F(1)
        else:
            best = min(coeffs[0], coeffs[-1]); best_u = F(0) if coeffs[0] <= coeffs[-1] else F(1)
        stack = [(coeffs, F(0), F(1), 0)]
        outer = None
        depth_limit_nodes = 0
        max_depth_seen = 0
        while stack:
            cf, u0, u1, depth = stack.pop()
            max_depth_seen = max(max_depth_seen, depth)
            cmax, cmin = max(cf), min(cf)
            # Every interval endpoint is feasible and carries a global witness.
            for value, u in ((cf[0], u0), (cf[-1], u1)):
                if (kind == "max" and value > best) or (kind == "min" and value < best):
                    best, best_u = value, u
            if kind == "max":
                if cmax <= best + tolq:
                    outer = cmax if outer is None else max(outer, cmax)
                    continue
            else:
                if cmin >= best - tolq:
                    outer = cmin if outer is None else min(outer, cmin)
                    continue
            if depth >= maxd:
                depth_limit_nodes += 1
                bound_value = cmax if kind == "max" else cmin
                outer = bound_value if outer is None else (max(outer, bound_value) if kind == "max" else min(outer, bound_value))
                continue
            a, b = split(cf)
            um = (u0 + u1) / 2
            stack += [(a, u0, um, depth + 1), (b, um, u1, depth + 1)]
        gap = best - outer if kind == "min" else outer - best
        return {
            "outer": outer,
            "attained": best,
            "witness_u": best_u,
            "gap": gap,
            "tolerance": tolq,
            "tolerance_met": gap <= tolq,
            "max_depth": maxd,
            "max_depth_seen": max_depth_seen,
            "depth_limit_reached": depth_limit_nodes > 0,
            "depth_limit_nodes": depth_limit_nodes,
        }
    return {"min": bound("min"), "max": bound("max")}


def _exact_interval_json(cert, kind):
    """JSON-safe exact certificate, retaining the interval orientation explicitly."""
    if kind == "min":
        interval = [cert["outer"], cert["attained"]]
    else:
        interval = [cert["attained"], cert["outer"]]
    return {
        "interval": [_qstr(v) for v in interval],
        "outer": _qstr(cert["outer"]),
        "attained": _qstr(cert["attained"]),
        "witness_u": _qstr(cert["witness_u"]),
        "gap": _qstr(cert["gap"]),
        "tolerance": _qstr(cert["tolerance"]),
        "tolerance_met": cert["tolerance_met"],
        "max_depth": cert["max_depth"],
        "max_depth_seen": cert["max_depth_seen"],
        "depth_limit_reached": cert["depth_limit_reached"],
        "depth_limit_nodes": cert["depth_limit_nodes"],
    }


def certificate_json(cert):
    """Serialize an exact min/max result while keeping float compatibility fields."""
    lo, hi = cert["min"], cert["max"]
    return {
        "min": _float_outward(lo["outer"], "down"),
        "max": _float_outward(hi["outer"], "up"),
        "exact": {"min": _exact_interval_json(lo, "min"), "max": _exact_interval_json(hi, "max")},
        "tolerance_met": lo["tolerance_met"] and hi["tolerance_met"],
        "status": "CERTIFIED" if lo["tolerance_met"] and hi["tolerance_met"] else "UNRESOLVED",
    }


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _case_key(row):
    return tuple(row[k] for k in ("q", "e0", "e1", "c"))


def validate_vertical_result(path="n4_vertical_all.json"):
    """Validate the exact 16-case vertical certificate before boundary reduction."""
    with open(path) as stream:
        vertical = json.load(stream)
    if vertical.get("T") != T or vertical.get("p") != str(P):
        raise RuntimeError("n4_vertical_all.json T/p do not match the current boundary recurrence")
    expected = {(str(q), str(e0), str(e1), str(c)) for q in V.QS for e0, e1, c in itertools.product(V.E0S, V.E1S, V.CS)}
    if len(vertical.get("stage0", [])) != 16 or len(vertical.get("certificates", [])) != 16:
        raise RuntimeError("n4_vertical_all.json must contain exactly 16 stage0 and certificate records")
    observed_stage0 = {_case_key(row) for row in vertical.get("stage0", [])}
    observed_certs = {_case_key(row) for row in vertical.get("certificates", [])}
    if observed_stage0 != expected or observed_certs != expected:
        raise RuntimeError("n4_vertical_all.json does not contain the matching 16 parameter keys")
    required = {"P0": "certified neg", "PR": "certified neg", "PC": "certified neg", "Delta": "certified pos"}
    failures = []
    for row in vertical["certificates"]:
        for objective, status in required.items():
            observed = row.get("objectives", {}).get(objective, {})
            if observed.get("status") != status:
                failures.append((_case_key(row), objective, observed.get("status")))
            if observed.get("leaves", {}).get("unresolved", 0) != 0:
                failures.append((_case_key(row), objective, "unresolved leaves"))
            if status == "certified neg" and (observed.get("leaves", {}).get("pos", 0) != 0 or observed.get("leaves", {}).get("neg", 0) == 0):
                failures.append((_case_key(row), objective, "inconsistent negative leaves"))
            if status == "certified pos" and (observed.get("leaves", {}).get("neg", 0) != 0 or observed.get("leaves", {}).get("pos", 0) == 0):
                failures.append((_case_key(row), objective, "inconsistent positive leaves"))
    if failures:
        raise RuntimeError(f"vertical root-box statuses do not support boundary reduction: {failures}")
    return vertical


def classify_whole_class(delta_min, delta_max):
    """Classify ranking using exact outer bounds first, then exact attained signs."""
    if delta_min["outer"] > 0:
        return "R wins throughout"
    if delta_max["outer"] < 0:
        return "C wins throughout"
    if delta_min["attained"] < 0 and delta_max["attained"] > 0:
        return "UNIDENTIFIED"
    return "UNRESOLVED"


def _whole_certificate(curve_name, objective, kind, curve_cert):
    """Bind a whole-class endpoint to its exact boundary-curve certificate."""
    return {
        "curve": curve_name,
        "objective": objective,
        "bound": kind,
        **curve_cert["exact"][kind],
    }

def sign_certificate(coeffs, u0=F(0), u1=F(1), depth=0, maxd=60):
    """Return exact sign tests and isolated root brackets on ``[u0, u1]``.

    Root pieces are certified to contain exactly one simple root (polynomial strictly
    monotone with opposite endpoint signs) and are narrowed to width ``2^-40``.
    At ``maxd`` an ``unresolved`` piece is returned explicitly; root pieces alone
    are not a claim that a complete covering partition was returned.
    """
    if all(c > 0 for c in coeffs): return [(u0, u1, "pos")]
    if all(c < 0 for c in coeffs): return [(u0, u1, "neg")]
    n = len(coeffs) - 1
    d = [n * (coeffs[k + 1] - coeffs[k]) for k in range(n)]
    f0, f1 = coeffs[0], coeffs[-1]
    if (all(x > 0 for x in d) or all(x < 0 for x in d)) and f0 * f1 < 0:
        # monotone with a sign change: bisect exactly to isolate the root
        lo, hi, cl, ch = u0, u1, coeffs, None
        while hi - lo > F(1, 2**40):
            a, b = split(cl); mid = (lo + hi) / 2
            if a[0] * a[-1] <= 0: cl, hi = a, mid
            else: cl, lo = b, mid
        return [(lo, hi, "root")]
    if depth >= maxd: return [(u0, u1, "unresolved")]
    a, b = split(coeffs); um = (u0 + u1) / 2
    return sign_certificate(a, u0, um, depth + 1, maxd) + sign_certificate(b, um, u1, depth + 1, maxd)

def main():
    n = 4; t0 = time.time()
    here = __import__("os").path.dirname(__file__)
    vertical_path = __import__("os").path.join(here, "n4_vertical_all.json")
    output_path = __import__("os").path.join(here, "n4_boundary_results.json")
    validate_vertical_result(vertical_path)
    out = {
        "T": T,
        "p": str(P),
        "metadata": {
            "vertical_source_sha256": _sha256(__import__("os").path.join(here, "n4_vertical.py")),
            "vertical_result_sha256": _sha256(vertical_path),
            "boundary_source_sha256": _sha256(__file__),
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "gmpy2_version": gmpy2.version(),
            "vertical_rootbox_statuses": {"P0": "certified neg", "PR": "certified neg", "PC": "certified neg", "Delta": "certified pos"},
            "tolerance": "1/1000000000",
        },
        "cases": [],
    }
    all_tolerances_met = True
    for q in V.QS:
        xL, xU = V.xbounds(P, q)
        for e0, e1, c in itertools.product(V.E0S, V.E1S, V.CS):
            t1 = time.time()
            parts = {a: V.affine_parts(a, n, e0, e1, c, P, q) for a in ("0", "R", "C")}
            rec = {"q": str(q), "e0": str(e0), "e1": str(e1), "c": str(c), "curves": {}}
            co = {}
            for cname, yfun in (("lower", V.Lfun), ("upper", V.Ufun)):
                cc = {}
                for a in ("0", "R", "C"):
                    base, Bx, Cy = parts[a]
                    cc[a] = objective_coeffs(base, Bx, Cy, V.start_state(a, n), xL, xU, yfun, P, q)
                    # cross-check against direct exact propagation at u = 1/4, 3/4
                    for u in (F(1, 4), F(3, 4)):
                        x = xL + u * (xU - xL); y = yfun(x, P, q)
                        direct = V.persist(V.kernel_at(base, Bx, Cy, x, y), V.start_state(a, n))
                        assert evaluate(cc[a], u) == q_(direct), ("cross-check failed", cname, a, u)
                cc["D"] = [r - s for r, s in zip(cc["R"], cc["C"])]
                co[cname] = cc
                crv = {}
                for name, key in (("P0", "0"), ("PR", "R"), ("PC", "C"), ("Delta", "D")):
                    cert = enclose(cc[key])
                    summary = certificate_json(cert)
                    summary.update({
                        "degree": len(cc[key]) - 1,
                        "at_u0": float(cc[key][0]),
                        "at_u1": float(cc[key][-1]),
                        "at_u0_approx": float(cc[key][0]),
                        "at_u1_approx": float(cc[key][-1]),
                        "decimal_display": "approximate",
                        "at_u0_exact": _qstr(cc[key][0]),
                        "at_u1_exact": _qstr(cc[key][-1]),
                    })
                    all_tolerances_met = all_tolerances_met and summary["tolerance_met"]
                    crv[name] = summary
                pieces = sign_certificate(cc["D"])
                crv["Delta_sign"] = [(float(a), float(b), s) for a, b, s in pieces]
                crv["Delta_roots_x"] = [(float(xL + a * (xU - xL)), float(xL + b * (xU - xL))) for a, b, s in pieces if s == "root"]
                crv["Delta_sign_exact"] = [[_qstr(a), _qstr(b), s] for a, b, s in pieces]
                crv["Delta_roots_u_exact"] = [[_qstr(a), _qstr(b)] for a, b, s in pieces if s == "root"]
                crv["Delta_roots_x_exact"] = [[_qstr(xL + a * (xU - xL)), _qstr(xL + b * (xU - xL))] for a, b, s in pieces if s == "root"]
                crv["Delta_sign_unresolved"] = any(s == "unresolved" for _, _, s in pieces)
                rec["curves"][cname] = crv
            # whole-class bounds from vertical monotonicity: P decreasing in y -> max on lower, min on upper; Delta increasing in y -> max upper, min lower
            wc = {}
            for name in ("P0", "PR", "PC"):
                lower = rec["curves"]["lower"][name]; upper = rec["curves"]["upper"][name]
                wc[name] = {
                    "min": upper["min"], "max": lower["max"],
                    "min_certificate": _whole_certificate("upper", name, "min", upper),
                    "max_certificate": _whole_certificate("lower", name, "max", lower),
                    "tolerance_met": upper["tolerance_met"] and lower["tolerance_met"],
                    "status": "CERTIFIED" if upper["tolerance_met"] and lower["tolerance_met"] else "UNRESOLVED",
                }
            lower = rec["curves"]["lower"]["Delta"]; upper = rec["curves"]["upper"]["Delta"]
            wc["Delta"] = {
                "min": lower["min"], "max": upper["max"],
                "min_certificate": _whole_certificate("lower", "Delta", "min", lower),
                "max_certificate": _whole_certificate("upper", "Delta", "max", upper),
                "tolerance_met": lower["tolerance_met"] and upper["tolerance_met"],
                "status": "CERTIFIED" if lower["tolerance_met"] and upper["tolerance_met"] else "UNRESOLVED",
            }
            dmin = {"outer": Q(lower["exact"]["min"]["outer"]), "attained": Q(lower["exact"]["min"]["attained"])}
            dmax = {"outer": Q(upper["exact"]["max"]["outer"]), "attained": Q(upper["exact"]["max"]["attained"])}
            wc["Delta_sign_over_class"] = classify_whole_class(dmin, dmax)
            rec["whole_class"] = wc; rec["seconds"] = round(time.time() - t1, 1)
            out["cases"].append(rec)
            print(f"q={str(q):>5} e0={str(e0):>4} e1={str(e1):>4} c={str(c):>4} | P0 [{wc['P0']['min']:.4f},{wc['P0']['max']:.4f}] PR [{wc['PR']['min']:.4f},{wc['PR']['max']:.4f}] PC [{wc['PC']['min']:.4f},{wc['PC']['max']:.4f}] | Delta [{float(dmin['outer']):+.5f},{float(dmax['outer']):+.5f}] {wc['Delta_sign_over_class']:<18} | lower-curve sign {[s for _,_,s in rec['curves']['lower']['Delta_sign']]} upper-curve sign {[s for _,_,s in rec['curves']['upper']['Delta_sign']]} [{rec['seconds']}s]", flush=True)
    out["metadata"]["all_tolerances_met"] = all_tolerances_met
    with open(output_path, "w") as f:
        json.dump(out, f, indent=1)
    if not all_tolerances_met:
        raise RuntimeError("boundary certificate tolerance was not met; statuses were persisted in n4_boundary_results.json")
    print(f"total {time.time()-t0:.1f}s; wrote {output_path}")

if __name__ == "__main__":
    main()
