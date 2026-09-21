"""D3: does a once-only randomised choice between R and C hedge away the unidentified ranking?
With certified whole-class extrema L < 0 < U of Delta = PR - PC, choosing R with probability alpha has worst-case
regret max((1 - alpha) U, alpha (-L)); the minimax mixture is alpha* = U/(U - L) with value Vmix = -L U/(U - L),
against the deterministic value V = min(-L, U) (Astra, thread 49 #457). Exact rationals from the certified records.
Also cross-checks two D2 numbers (Astra #458) with Fable's floating-point implementation."""
import sys, json, os; sys.path.insert(0, os.path.dirname(__file__))
from fractions import Fraction as F
import numpy as np
import b_screen_v2 as B

def frac(s): return F(str(s))
out = {}
# ---- n=3 default reversal case: Delta monotone in lambda, extrema at endpoints (certified) ----
J3 = json.load(open("n3_certificate_results.json"))
c3 = next(c for c in J3["cases"] if c["id"] == "n3_q3/100_e01/50_e19/10_c1/10")
ev = c3["metrics"]["delta_PR_minus_PC"]["endpoint_values"]
L3, U3 = sorted([frac(ev["lambda_0"]), frac(ev["lambda_1"])])
# ---- n=4 default reversal case: exact certified min/max (audited results carry exact outer/attained intervals) ----
J4 = json.load(open("n4_boundary_results.json"))
c4 = next(c for c in J4["cases"] if c["q"] == "3/100" and c["e0"] == "1/50" and c["e1"] == "9/10" and c["c"] == "1/10")
D = c4["whole_class"]["Delta"]
def pick(cert):
    iv = cert.get("interval") or cert.get("outer") or cert.get("bound_interval")
    return [frac(x) for x in iv]
mn = pick(D["min_certificate"]); mx = pick(D["max_certificate"]) if "max_certificate" in D else None
L4 = mn[0]  # outer lower bound (valid for regret bounds; zero gap certified in audit)
U4 = mx[1] if mx else frac(D["max"])
def analyse(name, L, U):
    assert L < 0 < U
    V = min(-L, U); alpha = U / (U - L); Vmix = (-L) * U / (U - L)
    rec = {"L": str(L), "U": str(U), "L_float": float(L), "U_float": float(U), "V_deterministic": float(V), "alpha_star_prob_R": float(alpha),
           "V_mixed": float(Vmix), "ratio_mixed_over_deterministic": float(Vmix / V)}
    out[name] = rec
    print(f"{name}: L={float(L):+.8f} U={float(U):+.8f} | deterministic minimax regret V={float(V):.6f} | lottery alpha*(R)={float(alpha):.4f}, Vmix={float(Vmix):.6f} | Vmix/V={float(Vmix/V):.3f}")
analyse("n3_default_reversal_whole_class", L3, U3)
analyse("n4_default_reversal_whole_class", L4, U4)
# ---- selected-pair witnesses (lower bounds on class values); values as reported in Astra's audited notes ----
pairs = {"n4_equal_triples_pair": (-0.005703179781404491, 0.005981975678913741),
         "n5_equal_quadruples_pair": (-0.0028385795533267424, 0.0023575756860951293),
         "n4_background_e0_0.1_e1_0.9_c_0.45_pair": (-0.07583193, 0.11023902),
         "n4_p0.3_e1_0.9_c_0.3_T50_pair": (-0.07102939, 0.03807481)}
for k, (l, u) in pairs.items():
    V = min(-l, u); Vmix = (-l) * u / (u - l)
    out[k] = {"L": l, "U": u, "V_pair": V, "V_mixed_pair": Vmix, "ratio": Vmix / V}
    print(f"{k}: pair V={V:.6f}, lottery Vmix={Vmix:.6f}, ratio {Vmix/V:.3f}  (pair values are lower bounds for the class)")
# ---- D2 cross-checks with Fable's float implementation ----
p, q = 0.1, 0.03
# interior two-regime high witness (a, b, w) = (3/40, 9/10, 1/33): Delta at n=4 reported +0.01025984
a, b, w = 3/40, 9/10, 1/33
th = np.array([a, b]); wt = np.array([1 - w, w])
assert abs(np.dot(wt, th) - p) < 1e-15 and abs(np.dot(wt, th**2) - q) < 1e-15
P0, PR, PC = B.run_case(4, 0.02, 0.9, 0.1, th, wt); out["check_two_regime_high_n4_Delta"] = PR - PC
print(f"cross-check two-regime (3/40, 9/10, 1/33) n=4 Delta = {PR-PC:+.8f}  (Astra: +0.01025984)")
# beta law with mean p, second moment q: alpha = p (p(1-p)/var - 1), beta = (1-p) (p(1-p)/var - 1)
var = q - p*p; k = p*(1-p)/var - 1; al, be = p*k, (1-p)*k
from scipy.special import roots_jacobi
# Gauss-Jacobi on [-1,1] with weight (1-x)^(be-1) (1+x)^(al-1); theta = (1+x)/2. Exact for polynomial integrands of
# degree <= 2N-1, and every kernel entry has degree <= 4 in theta.
xj, wj = roots_jacobi(64, be - 1, al - 1); thb = (1 + xj) / 2; wb = wj / wj.sum()
mom = [float(np.dot(wb, thb**j)) for j in (1,2,3,4)]
P0, PR, PC = B.run_case(4, 0.02, 0.9, 0.1, thb, wb); out["check_beta_n4_Delta"] = PR - PC; out["beta_params"] = [al, be]; out["beta_moments_gauss_jacobi"] = mom
print(f"cross-check beta(alpha={al:.4f}, beta={be:.4f}), 64-node Gauss-Jacobi: moments {mom[0]:.8f} {mom[1]:.8f} {mom[2]:.8f} {mom[3]:.8f} (exact m3 0.01281818); n=4 Delta = {PR-PC:+.8f}  (Astra: -0.08336312)")
json.dump(out, open("d3_mixed_policy.json", "w"), indent=1)
