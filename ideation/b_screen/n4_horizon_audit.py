"""Horizon audit: exact finite-T anchors and stable transient-state curves.

uv run --python 3.12 --with numpy --with gmpy2 --with matplotlib python ideation/b_screen/n4_horizon_audit.py
Float curves are exploratory. Exact anchors use the unchanged rational kernels.
Regret is for a deterministic, once-only choice of R or C under two fixed laws.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np
from gmpy2 import mpq as Q
import n3_certificate as N

ROOT = Path(__file__).resolve().parent
CASES = [
    ('reference', F(1,10), F(1,50), F(4,5), F(3,20)),
    ('higher_p_same_actions', F(3,10), F(1,50), F(4,5), F(3,20)),
    ('higher_p_reversal', F(3,10), F(1,50), F(9,10), F(3,10)),
    ('higher_background', F(1,10), F(1,10), F(9,10), F(9,20)),
]
HORIZONS = [0,1,5,10,20,50,100,200,500,1000,2000,10000]


def metrics(P):
    dl, dh = P['Glo_R']-P['Glo_C'], P['Ghi_R']-P['Ghi_C']
    regret_R, regret_C = max(0,-dl,-dh), max(0,dl,dh)
    return dict(delta_glo=dl, delta_ghi=dh, regret_R=regret_R,
                regret_C=regret_C, deterministic_minimax_regret=min(regret_R,regret_C))


def kernels(p,e0,e1,c):
    q=p*p+F(2,9)*p*(1-p)
    result={}
    for name,law in zip(('Glo','Ghi'),N.laws(p,q)):
        for action in ('R','C'):
            K=N.averaged_kernel(action,4,e0,e1,c,*law)
            assert all(sum(row)==1 and all(x>=0 for x in row) for row in K)
            assert K[0][0]==1
            result[name+'_'+action]=([[Q(v.numerator,v.denominator) for v in row] for row in K],7 if action=='R' else 4)
    return q,result


def exact_anchors(Ks):
    vectors={k:[Q(int(i==start)) for i in range(len(K))] for k,(K,start) in Ks.items()}
    anchors=[]
    for T in range(1,201):
        vectors={k:N.matvec(v,Ks[k][0]) for k,v in vectors.items()}
        if T in (20,50,200):
            P={k:sum(v[1:]) for k,v in vectors.items()}
            assert all(Q(0)<=x<=Q(1) for x in P.values())
            m=metrics(P)
            anchors.append(dict(T=T,persistence={k:str(v) for k,v in P.items()},
                                metrics={k:str(v) for k,v in m.items()},
                                display={k:float(v) for k,v in m.items()}))
    return anchors


def run_case(case):
    name,p,e0,e1,c=case
    q,Ks=kernels(p,e0,e1,c)
    transient={k:np.array(K,dtype=float)[1:,1:] for k,(K,_) in Ks.items()}
    vectors={k:np.eye(len(K))[start-1] for k,K in transient.items() for start in [Ks[k][1]]}
    rows=[]
    peaks={}
    for T in range(2001):
        P={k:float(v.sum()) for k,v in vectors.items()}
        m=metrics(P)
        if T in HORIZONS or (T<=1000 and T%5==0):
            rows.append(dict(T=T,persistence=P,**m))
        for k in ('delta_glo','delta_ghi','deterministic_minimax_regret'):
            if k not in peaks or abs(m[k])>abs(peaks[k]['value']):
                peaks[k]=dict(T=T,value=m[k])
        vectors={k:v@transient[k] for k,v in vectors.items()}
    # Long horizon uses transient powers, avoiding subtraction from an absorbing mass near one.
    P={k:float(np.linalg.matrix_power(K,10000)[Ks[k][1]-1].sum()) for k,K in transient.items()}
    rows.append(dict(T=10000,persistence=P,**metrics(P)))
    exact=exact_anchors(Ks)
    for a in exact:
        r=next(r for r in rows if r['T']==a['T'])
        assert all(abs(r[k]-v)<2e-12 for k,v in a['display'].items())
    return dict(name=name,p=str(p),q=str(q),e0=str(e0),e1=str(e1),c=str(c),
                curve=rows,peaks_over_T_0_to_2000=peaks,exact_anchors=exact)


def plot(cases):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for case in cases:
        rows=[r for r in case['curve'] if r['T']<=1000]
        t=[r['T'] for r in rows]
        label=case['name'].replace('_',' ')
        axes[0].plot(t,[r['deterministic_minimax_regret'] for r in rows],label=label)
    ref=cases[0]
    rows=[r for r in ref['curve'] if r['T']<=1000]
    for key,label in [('delta_glo','Glo'),('delta_ghi','Ghi')]:
        axes[1].plot([r['T'] for r in rows],[r[key] for r in rows],label=label)
    axes[0].set_ylabel('Deterministic minimax regret\n(two fixed laws; persistence probability)')
    axes[1].set_ylabel('Buffering minus connectivity persistence')
    axes[1].set_title('Reference setting')
    axes[1].axhline(0,color='grey',lw=.7)
    for ax in axes:
        ax.set_xlabel('Horizon (annual steps)')
        ax.legend(frameon=False,fontsize=8)
    fig.savefig(ROOT/'n4_horizon_audit.svg')
    plt.close(fig)


def main():
    cases=[run_case(case) for case in CASES]
    result=dict(n=4,rho='2/9',cases=cases,source_sha256={
        name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        for name in (Path(__file__).name,'n3_certificate.py')})
    (ROOT/'n4_horizon_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(cases)
    for case in cases:
        print(case['name'],case['peaks_over_T_0_to_2000'])
        for anchor in case['exact_anchors']:
            print(anchor['T'],anchor['display'])


if __name__=='__main__':
    main()
