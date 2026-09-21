"""D2: restricted hazard families at the existing n=3,4 reversal setting.

uv run --python 3.12 --with numpy --with scipy --with gmpy2 python ideation/b_screen/hazard_family_audit.py

Atomic, uniform-mixture and beta results use exact rational moments/kernels.
Probit/logit normal results are numerical, with independent quadrature checks.
This is a constructed sensitivity audit, not calibration or family selection.
"""
from fractions import Fraction as F
from math import comb, pi, exp, sqrt
from pathlib import Path
import hashlib
import json

from gmpy2 import mpq as Q
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, root
from scipy.special import expit, ndtr, ndtri, roots_hermitenorm

import n3_certificate as N
import n4_vertical as V

ROOT = Path(__file__).resolve().parent
P, Q2 = F(1, 10), F(3, 100)
E0, E1, C, T = F(1, 50), F(9, 10), F(1, 10), 50


def two_regime(a, h=F(0)):
    """Theta = A + U, U independent Uniform[-h,h]; h=0 gives two atoms."""
    q_centres = Q2 - h*h/3
    b = (q_centres - P*a)/(P-a)
    w = (P-a)/(b-a)
    assert 0 <= a-h <= a+h < b-h <= b+h <= 1
    assert 0 < w < 1
    moments = [sum(wt * sum(F(comb(k,j))*mu**(k-j)*h**j/(j+1)
                           for j in range(0,k+1,2))
                   for mu,wt in ((a,1-w),(b,w))) for k in range(5)]
    assert moments[:3] == [F(1), P, Q2]
    return moments, dict(a=str(a), b=str(b), weight_b=str(w), half_width=str(h))


def beta_moments():
    concentration = P*(1-P)/(Q2-P*P)-1
    alpha, beta = P*concentration, (1-P)*concentration
    m = [F(1)]
    for k in range(4):
        m.append(m[-1]*(alpha+k)/(alpha+beta+k))
    assert m[:3] == [F(1), P, Q2]
    return m, dict(alpha=str(alpha), beta=str(beta))


def normal_moments(transform):
    pairs = [quad(lambda z: float(transform(z))**k*exp(-z*z/2)/sqrt(2*pi),
                  -np.inf, np.inf, epsabs=2e-13, epsrel=2e-13, limit=200)
             for k in range(1,5)]
    return np.array([1]+[x[0] for x in pairs]), max(x[1] for x in pairs)


def numerical_families():
    z0 = ndtri(float(P))
    def probit(rho):
        return lambda z: ndtr((z0+sqrt(rho)*z)/sqrt(1-rho))
    rho = brentq(lambda r: normal_moments(probit(r))[0][2]-float(Q2),
                 1e-8, 1-1e-8, xtol=1e-14)
    def logit(params):
        mu, log_sigma = params
        return lambda z: expit(mu+exp(log_sigma)*z)
    def residual(params):
        return normal_moments(logit(params))[0][1:3]-[float(P),float(Q2)]
    solutions = [root(residual, guess, tol=1e-10)
                 for guess in ((-3.,0.),(-5.,1.),(-2.,-.5))]
    assert all(s.success and np.max(np.abs(residual(s.x)))<1e-11 for s in solutions)
    assert max(np.max(np.abs(s.x-solutions[0].x)) for s in solutions)<1e-8
    params=solutions[0].x
    laws = [('probit_normal',probit(rho),dict(latent_rho=rho)),
            ('logit_normal',logit(params),dict(mu=params[0],sigma=exp(params[1])))]
    out=[]
    for name,transform,parameters in laws:
        m,err=normal_moments(transform)
        z,w=roots_hermitenorm(1024)
        w=w/w.sum()
        gh=np.array([np.dot(w,transform(z)**k) for k in range(5)])
        discrepancy=float(np.max(np.abs(m-gh)))
        assert discrepancy<2e-11
        assert np.max(np.abs(m[1:3]-[float(P),float(Q2)]))<1e-11
        out.append((name,m,dict(**parameters,quad_error_estimate=err,
                              independent_hermite_moment_discrepancy=discrepancy)))
    return out


def evaluate(moments,n,coeffs,exact):
    values={}
    for action in ('R','C'):
        A=coeffs[n,action]
        size=len(A[0])
        K=[[sum(moments[k]*A[k][i][j] for k in range(n+1))
            for j in range(size)] for i in range(size)]
        start=2*n-1 if action=='R' else n
        if exact:
            assert all(sum(row)==1 and min(row)>=0 for row in K)
            K=[[Q(x.numerator,x.denominator) for x in row] for row in K]
            v=[Q(int(i==start)) for i in range(size)]
            for _ in range(T):
                v=N.matvec(v,K)
            values[action]=sum(v[1:])
        else:
            K=np.array(K,dtype=float)
            assert K.min()>-1e-12 and np.max(abs(K.sum(1)-1))<1e-12
            values[action]=np.linalg.matrix_power(K[1:,1:],T)[start-1].sum()
    values['delta']=values['R']-values['C']
    return dict(display={k:float(v) for k,v in values.items()},
                exact={k:str(v) for k,v in values.items()} if exact else None)


def main():
    coeffs={(n,action):V.moment_kernels(action,n,E0,E1,C)
            for n in (3,4) for action in ('R','C')}
    families=[]
    for name,a in (('endpoint_low',F(0)),('interior_low',F(1,100)),
                   ('interior_high',F(3,40)),('endpoint_high',F(7,90))):
        m,params=two_regime(a)
        # Exact agreement with direct atom averaging protects moment assembly.
        for n in (3,4):
            for action in ('R','C'):
                direct=N.averaged_kernel(action,n,E0,E1,C,
                    (a,F(params['b'])),(1-F(params['weight_b']),F(params['weight_b'])))
                A=coeffs[n,action]
                assert direct==[[sum(m[k]*A[k][i][j] for k in range(n+1))
                                 for j in range(len(direct))] for i in range(len(direct))]
        families.append((name,m,params,True))
    for name,a in (('uniform_low',F(1,100)),('uniform_high',F(3,40))):
        m,params=two_regime(a,F(1,1000))
        families.append((name,m,params,True))
    m,params=beta_moments()
    families.append(('beta',m,params,True))
    families.extend((name,m,params,False) for name,m,params in numerical_families())
    rows=[]
    for name,m,params,exact in families:
        results={str(n):evaluate(m,n,coeffs,exact) for n in (3,4)}
        rows.append(dict(name=name,parameters=params,moments=[str(v) for v in m],
                         exact_moments=exact,results=results))
        print(name,'m3,m4',float(m[3]),float(m[4]),
              'Delta3,4',*[results[str(n)]['display']['delta'] for n in (3,4)])
    for lower,upper in (('interior_low','interior_high'),('uniform_low','uniform_high')):
        for n in ('3','4'):
            lo=next(r for r in rows if r['name']==lower)['results'][n]['exact']['delta']
            hi=next(r for r in rows if r['name']==upper)['results'][n]['exact']['delta']
            assert Q(lo)<0<Q(hi)
    result=dict(p=str(P),q=str(Q2),e0=str(E0),e1=str(E1),c=str(C),T=T,
                families=rows,source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                for name in (Path(__file__).name,'n3_certificate.py','n4_vertical.py')})
    (ROOT/'hazard_family_audit.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
