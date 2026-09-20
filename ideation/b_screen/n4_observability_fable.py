"""Independent recomputation of Astra's two-law discrimination benchmark (thread 49, #406; n4_observability.md).
Two known laws (the n=4 equal-triples witness at m3=21/1000), equal priors, iid perfectly observed years, 0-1 loss.
Bayes error = 0.5 * sum over count vectors of min(P, Q). Exact enumeration of multinomial count vectors, float64."""
import numpy as np
from math import lgamma, comb
p, q, x = 0.1, 0.03, 0.021
Bq = (x - p*q)/(q - p*p); Cq = (p*x - q*q)/(q - p*p); r = (Bq*Bq - 4*Cq)**0.5
zm, zp = (Bq - r)/2, (Bq + r)/2; lo_t = np.array([zm, zp]); lo_w = np.array([(zp - p)/(zp - zm), (p - zm)/(zp - zm)])
a = (q - x)/(p - q); up_t = np.array([0, a, 1.0]); up_w = np.array([1 - p - (p - q)/a, (p - q)/(a*(1 - a)), p - (p - q)/(1 - a)])
def Klaw(t, w): return np.array([comb(4, k) * np.dot(w, t**k * (1 - t)**(4 - k)) for k in range(5)])
PL, PU = Klaw(lo_t, lo_w), Klaw(up_t, up_w)
print("K-law lower:", np.round(PL, 6), " upper:", np.round(PU, 6), " sums", PL.sum(), PU.sum())
print("one-year Bayes error = 0.5*(1 - TV) =", 0.5 * (1 - 0.5 * np.abs(PL - PU).sum()))
def bayes_error_full(Y, P, Q):
    logP, logQ = np.log(P), np.log(Q); lgY = lgamma(Y + 1)
    lg = np.array([lgamma(k + 1) for k in range(Y + 1)])
    total = 0.0
    for k0 in range(Y + 1):
        rest = Y - k0
        # enumerate k1,k2,k3 with k1+k2+k3 <= rest; k4 = rest - k1 - k2 - k3
        k1 = np.arange(rest + 1)
        k1g, k2g = np.meshgrid(k1, k1, indexing="ij"); mask = k1g + k2g <= rest
        k1v, k2v = k1g[mask], k2g[mask]
        # for each (k1,k2), k3 from 0..rest-k1-k2
        maxk3 = rest - k1v - k2v
        # expand
        reps = maxk3 + 1
        k1e = np.repeat(k1v, reps); k2e = np.repeat(k2v, reps)
        k3e = np.concatenate([np.arange(m + 1) for m in maxk3]) if len(maxk3) else np.array([], dtype=int)
        k4e = rest - k1e - k2e - k3e
        ks = np.stack([np.full_like(k1e, k0), k1e, k2e, k3e, k4e])
        logmult = lgY - lg[ks].sum(0)
        lp = logmult + (ks * logP[:, None]).sum(0); lq = logmult + (ks * logQ[:, None]).sum(0)
        total += np.exp(np.minimum(lp, lq)).sum()
    return 0.5 * total
def bayes_error_k4only(Y, P, Q):
    pl, pu = P[4], Q[4]; ks = np.arange(Y + 1)
    lgb = np.array([lgamma(Y + 1) - lgamma(k + 1) - lgamma(Y - k + 1) for k in ks])
    lp = lgb + ks * np.log(pl) + (Y - ks) * np.log(1 - pl); lq = lgb + ks * np.log(pu) + (Y - ks) * np.log(1 - pu)
    return 0.5 * np.exp(np.minimum(lp, lq)).sum()
# controls: identical laws -> 0.5
print("identical-law control (Y=30):", bayes_error_full(30, PL, PL))
for Y in (25, 50, 100, 200):
    ef = bayes_error_full(Y, PL, PU); ek = bayes_error_k4only(Y, PL, PU)
    print(f"Y={Y:>3}: Bayes error, full K record {ef:.4f}; K=4-only record {ek:.4f}   (Astra: see #406)")
