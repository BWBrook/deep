#!/usr/bin/env python3
"""Independent scalar-boundary, shape, and decision-criterion checks.
No repository modules are imported. NumPy is required. All reported signs
and class extrema use exact arithmetic; floating values are displays only.
"""
from fractions import Fraction as F
from math import comb,lcm
from pathlib import Path
from itertools import product
import json,time,hashlib
import numpy as np
import independent_checks as I
import independent_certificates as C

T=50

def matrices_integer(controls):
    D=lcm(*(v.denominator for d in controls for K in d.values() for row in K for v in row))
    mats={a:[np.array([[v.numerator*(D//v.denominator) for v in row] for row in d[a]],dtype=object) for d in controls] for a in ['R','C']}
    return mats,D

def raw_by_action(mats,degree=1,bivariate=False):
    out={}
    for a,Ms in mats.items():
        size=len(Ms[0])
        v=np.zeros((1,1,size) if bivariate else (1,size),dtype=object)
        if bivariate:v[0,0,-1]=1
        else:v[0,-1]=1
        for t in range(1,T+1):
            if bivariate:
                w=np.zeros((t+1,t+1,size),dtype=object)
                for (i,j),M in zip([(0,0),(1,0),(0,1),(1,1)],Ms):
                    w[i:i+t,j:j+t]+=(v.reshape(-1,size)@M).reshape(t,t,size)
            else:
                w=np.zeros((degree*t+1,size),dtype=object)
                for j,M in enumerate(Ms):w[j:j+len(v)]+=comb(degree,j)*(v@M)
            v=w
        out[a]=np.sum(v[...,1:],axis=-1)
    out['Delta']=out['R']-out['C']
    return out

def direct(parts,p,q,pt):
    Ks=C.at_point(parts,p,q,pt)
    # Feasible points only: row sums and signs checked by trajectory + here.
    assert all(sum(row)==1 and min(row)>=0 for K in Ks.values() for row in K)
    v={a:I.trajectory(K)[-1] for a,K in Ks.items()}
    v['Delta']=v['R']-v['C']
    return v

def unieval(raw,D,u):
    n=len(raw)-1
    return sum(F(int(raw[j]))*u**j*(1-u)**(n-j) for j in range(n+1))/D**T

def bieval(raw,D,u,v):
    return sum(F(int(raw[i,j]))*u**i*(1-u)**(T-i)*v**j*(1-v)**(T-j) for i in range(T+1) for j in range(T+1))/D**T

def dsig(raw):
    n=len(raw)-1
    vals=[int(raw[j+1])*comb(n,j)-int(raw[j])*comb(n,j+1) for j in range(n)]
    return 1 if min(vals)>0 else -1 if max(vals)<0 else 0

def xlimits(p,q):return q*q/p,q-(p-q)**2/(1-p)
def lower(x,p,q):return (x*x-2*p*q*x+q**3)/(q-p*p)
def upper(x,p,q):return x-(q-x)**2/(p-q)

def main():
    start=time.time();p=F(1,10)
    out={'commit':I.COMMIT,'independent':True,'scalar_cases':[],'reference_maximin':[]}
    parts_cache={}
    for e0,e1,c in product([F(1,50),F(1,10)],[F(1,2),F(9,10)],[F(1,50),F(1,10)]):
        parts=C.conditioned_kernels(e0,e1,c);parts_cache[e0,e1,c]=parts
        for q in [F(3,100),F(6,100)]:
            xl,xh=xlimits(p,q);rec={'q':str(q),'e0':str(e0),'e1':str(e1),'c':str(c),'boundaries':[]}
            for name,fun in [('lower',lower),('upper',upper)]:
                pts=[(xl+u*(xh-xl),fun(xl+u*(xh-xl),p,q)) for u in [F(0),F(1,2),F(1)]]
                K0,Kmid,K2=[C.at_point(parts,p,q,pt) for pt in pts]
                K1={a:[[2*b-(s+t)/2 for s,b,t in zip(r0,rm,r2)] for r0,rm,r2 in zip(K0[a],Kmid[a],K2[a])] for a in ['R','C']}
                mats,D=matrices_integer([K0,K1,K2]);raw=raw_by_action(mats,degree=2)['Delta']
                sign=dsig(raw)
                coeff=[F(int(v),comb(2*T,j)*D**T) for j,v in enumerate(raw)]
                endmin,endmax=min(coeff[0],coeff[-1]),max(coeff[0],coeff[-1])
                # Even if not monotone, an endpoint envelope suffices for extrema.
                assert min(coeff)==endmin and max(coeff)==endmax
                for u in [F(1,4),F(3,4)]:
                    x=xl+u*(xh-xl);expected=direct(parts,p,q,(x,fun(x,p,q)))['Delta']
                    assert unieval(raw,D,u)==expected
                rec['boundaries'].append({'name':name,'degree':2*T,'derivative_sign':sign,
                    'all_coefficients_between_endpoint_values':True,'min':str(endmin),'max':str(endmax),'interior_checks':['1/4','3/4']})
            out['scalar_cases'].append(rec)
    print('Scalar boundaries certified:',len(out['scalar_cases'])*2,flush=True)
    parts=parts_cache[F(1,50),F(9,10),F(1,10)]
    # New reviewer analysis: exact reference maximin persistence on E.
    for q in [F(3,100),F(6,100)]:
        verts=C.polygon(p,q);xl,xh=min(x for x,y in verts),max(x for x,y in verts);yl,yh=min(y for x,y in verts),max(y for x,y in verts)
        mats,D=C.controls_integer([(xl,yl),(xh,yl),(xl,yh),(xh,yh)],parts,p,q)
        raws=raw_by_action(mats,bivariate=True)
        rec={'q':str(q),'vertical_signs':{},'edges':[]}
        for a,raw in raws.items():
            dif=[int(raw[i,j+1])*comb(T,j)-int(raw[i,j])*comb(T,j+1) for i in range(T+1) for j in range(T)]
            sign=1 if min(dif)>0 else -1 if max(dif)<0 else 0
            assert sign!=0
            rec['vertical_signs'][a]=sign
        # Verify every bivariate objective against direct propagation at two feasible points.
        for pt in [(sum(x for x,y in verts)/len(verts),sum(y for x,y in verts)/len(verts)),verts[0]]:
            expected=direct(parts,p,q,pt);u,v=(pt[0]-xl)/(xh-xl),(pt[1]-yl)/(yh-yl)
            for a,raw in raws.items():assert bieval(raw,D,u,v)==expected[a]
        for j,v0 in enumerate(verts):
            v1=verts[(j+1)%len(verts)];mats1,D1=C.controls_integer([v0,v1],parts,p,q);raw1=raw_by_action(mats1)
            signs={a:dsig(raw) for a,raw in raw1.items()};assert all(signs.values())
            for u in [F(1,4),F(3,4)]:
                pt=(v0[0]+u*(v1[0]-v0[0]),v0[1]+u*(v1[1]-v0[1]));expected=direct(parts,p,q,pt)
                for a,raw in raw1.items():assert unieval(raw,D1,u)==expected[a]
            rec['edges'].append(signs)
        values=[direct(parts,p,q,pt) for pt in verts]
        rec['minima']={a:I.pack(min(d[a] for d in values)) for a in ['R','C']}
        rec['maximin_action']=max(['R','C'],key=lambda a:min(d[a] for d in values))
        rec['minimax_regret']=I.regret(min(d['Delta'] for d in values),max(d['Delta'] for d in values))
        rec['global_proof']='Nonzero uniform vertical derivative and monotone polygon edges for each objective; vertex extrema.'
        out['reference_maximin'].append(rec)
        print('Reference maximin q',q,rec['maximin_action'],{a:d['value'] for a,d in rec['minima'].items()},flush=True)
    # Equal-triples witness, full exact labelled propagation.
    q=F(3,100);x=F(21,1000)
    out['equal_triples']=[]
    for fun in [lower,upper]:
        y=fun(x,p,q);law=I.moment_to_law([F(1),p,q,x,y]);info=I.information(law)
        assert all(v==x for v in info[2])
        vals=direct(parts,p,q,(x,y))
        out['equal_triples'].append({'m3':str(x),'m4':str(y),'values':{k:I.pack(v) for k,v in vals.items()}})
    # Reference Khinchin envelope: rational witness m and upper-boundary cutoff.
    r=F(79,4000);vals=direct(parts,p,q,(r,upper(r,p,q)));assert vals['Delta'] < -F(427,100000)
    out['unimodal_envelope']={k:I.pack(v) for k,v in vals.items()}
    # Strictly interior E witnesses: give all 16 patterns positive probability.
    # These are attained examples, not newly optimised extrema.
    base=I.scalar_law(p,F(3,100))
    out['full_support_witnesses']=[]
    for pi in [[F(78,100),F(4,100),F(18,100),F(0),F(0)],
               [F(69,100),F(28,100),F(0),F(0),F(3,100)]]:
        extreme=I.count_to_law(pi)
        law=tuple(F(9,10)*a+F(1,10)*b for a,b in zip(extreme,base))
        assert min(law)>0
        inf=I.information(law)
        assert set(inf[0])=={p} and set(inf[1])=={F(3,100)}
        values,_=I.evaluate(law,F(1,50),F(9,10),F(1,10),baseline=True)
        out['full_support_witnesses'].append({
            'hit_count_probabilities':[str(sum(law[h] for h in range(16) if h.bit_count()==k)) for k in range(5)],
            'values':{k:I.pack(v) for k,v in values.items()}})
    assert F(out['full_support_witnesses'][0]['values']['Delta']['exact'])<0
    assert F(out['full_support_witnesses'][1]['values']['Delta']['exact'])>0
    # Coefficient identity for the global Khinchin third-moment cutoff.
    r=F(79,4000);m0=F(33,500)
    lhs=[32000*r-599,40000*r-1780,F(7500)]
    rhs=[7500*m0*m0+F(33,100),-15000*m0,F(7500)]
    assert lhs==rhs
    out['shape_cutoff_identity']={'coefficient_identity_verified':True,
        'polynomial_coefficients_low_to_high':[str(v) for v in lhs],
        'denominator':'8000*(5*m+4)>0 for 0<=m<=11/80',
        'consequence':'m3<h_cutoff=79/4000 under the derived Khinchin moment bound'}
    out['seconds']=time.time()-start
    out['source_sha256']={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['additional_checks.py','independent_checks.py','independent_certificates.py']}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('All additional checks passed',round(time.time()-start,2),'seconds',flush=True)

if __name__=='__main__':main()
