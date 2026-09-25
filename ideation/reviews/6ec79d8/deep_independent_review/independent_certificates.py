#!/usr/bin/env python3
"""Independent exact certificates for the four-patch exchangeable class.

Uses the independently factorised labelled model in independent_checks.py.
Reconstructs count-conditioned kernels, verifies exact lumpability, and
forms polynomial extensions from count probabilities. Bernstein coefficients
are built as INTEGER weighted-path sums, not with the repository's rational
Bernstein recurrence. Each coefficient is the weighted sum divided by a
binomial normalisation and a common kernel denominator to the horizon power.

Default: all 56 polygon edges, both signs of their derivatives.
--vertical: also certify the vertical derivative over each entire bounding
rectangle. Combined with monotone edges, this certifies global extrema at
vertices, without assuming vertex optimality for a nonlinear objective.

Requires NumPy, in addition to the Python standard library.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from itertools import product
from math import comb,lcm
from pathlib import Path
import hashlib,json,time
import numpy as np
import independent_checks as I


def lump(K,action,n=4):
    ids=[((h&1)*n+(h>>1).bit_count()) if action=='R' else h.bit_count() for h in range(1<<n)]
    groups=[[h for h,i in enumerate(ids) if i==g] for g in range(max(ids)+1)]
    out=[]
    for members in groups:
        rows=[tuple(sum(K[h][j] for j in dest) for dest in groups) for h in members]
        assert len(set(rows))==1, 'Exact lumpability failed'
        out.append(list(rows[0]))
    return out


def cross(o,a,b):
    return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])


def hull(pts):
    pts=sorted(set(pts))
    lo=[]
    for p in pts:
        while len(lo)>=2 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    hi=[]
    for p in reversed(pts):
        while len(hi)>=2 and cross(hi[-2],hi[-1],p)<=0: hi.pop()
        hi.append(p)
    return lo[:-1]+hi[:-1]


def polygon(p,q):
    return hull([(pi[3]/4+pi[4],pi[4]) for pi in I.e_vertices(p,q)])


def pi_xy(p,q,x,y):
    return [1-4*p+6*q-4*x+y,4*p-12*q+12*x-4*y,6*q-12*x+6*y,4*x-4*y,y]


def conditioned_kernels(e0,e1,c):
    out={}
    for action in ['R','C']:
        out[action]=[]
        for count in range(5):
            pi=[F(int(k==count)) for k in range(5)]
            out[action].append(lump(I.kernel(I.count_to_law(pi),e0,e1,c,action),action))
    return out


def at_point(parts,p,q,point):
    pi=pi_xy(p,q,*point)
    return {a:[[sum(pi[k]*parts[a][k][r][s] for k in range(5))
                for s in range(len(parts[a][0]))] for r in range(len(parts[a][0]))] for a in ['R','C']}


def controls_integer(points,parts,p,q):
    controls=[at_point(parts,p,q,pt) for pt in points]
    D=lcm(*(v.denominator for d in controls for K in d.values() for row in K for v in row))
    matrices={a:[np.array([[v.numerator*(D//v.denominator) for v in row] for row in d[a]],dtype=object)
                  for d in controls] for a in ['R','C']}
    return matrices,D


def univariate_raw(matrices,T=50):
    out={}
    for a,Ms in matrices.items():
        size=len(Ms[0]); v=np.zeros((1,size),dtype=object); v[0,-1]=1
        for t in range(1,T+1):
            w=np.zeros((t+1,size),dtype=object)
            w[:-1]+=v@Ms[0]; w[1:]+=v@Ms[1]; v=w
        out[a]=np.sum(v[:,1:],axis=1)
    return out['R']-out['C']


def bivariate_raw(matrices,T=50):
    out={}
    for a,Ms in matrices.items():
        size=len(Ms[0]); v=np.zeros((1,1,size),dtype=object); v[0,0,-1]=1
        for t in range(1,T+1):
            w=np.zeros((t+1,t+1,size),dtype=object)
            for (i,j),M in zip([(0,0),(1,0),(0,1),(1,1)],Ms):
                w[i:i+t,j:j+t]+=(v.reshape(-1,size)@M).reshape(t,t,size)
            v=w
        out[a]=np.sum(v[:,:,1:],axis=2)
    return out['R']-out['C']


def edge_certificate(matrices,D,T=50):
    raw=univariate_raw(matrices,T)
    d=[raw[i+1]*comb(T,i)-raw[i]*comb(T,i+1) for i in range(T)]
    sign='positive' if min(d)>0 else 'negative' if max(d)<0 else 'not_uniform'
    assert sign!='not_uniform'
    # Store exact end values to connect certificate to attained extrema.
    return {'derivative_sign':sign,'degree':T,'endpoint_deltas':[str(F(raw[0],D**T)),str(F(raw[-1],D**T))],
            'integer_coefficient_hash':hashlib.sha256(str(raw.tolist()).encode()).hexdigest()}


def vertical_certificate(matrices,D,ywidth,T=50):
    raw=bivariate_raw(matrices,T)
    d=[[raw[i,j+1]*comb(T,j)-raw[i,j]*comb(T,j+1) for j in range(T)] for i in range(T+1)]
    assert all(v>0 for row in d for v in row),'Nonpositive vertical Bernstein coefficient'
    # Exact minimum coefficient of physical derivative dDelta/dy.
    coeff=[F(T*d[i][j],comb(T,i)*comb(T,j)*comb(T,j+1)*D**T)/ywidth
           for i in range(T+1) for j in range(T)]
    bound=min(coeff)
    return {'all_coefficients_strictly_positive':True,'coefficient_count':len(coeff),
            'min_physical_derivative_coefficient':str(bound),'min_physical_derivative_float':float(bound),
            'integer_coefficient_hash':hashlib.sha256(str(raw.tolist()).encode()).hexdigest()}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--vertical',action='store_true');ap.add_argument('--reference-only',action='store_true');args=ap.parse_args()
    start=time.time();p=F(1,10);out={'commit':I.COMMIT,'independent_implementation':True,'cases':[],
        'proof':'positive vertical Bernstein coefficients on containing rectangle plus monotone edges',
        'vertical_run':args.vertical,'reference_only':args.reference_only}
    for e0,e1,c in product([F(1,50),F(1,10)],[F(1,2),F(9,10)],[F(1,50),F(1,10)]):
        if args.reference_only and (e0,e1,c)!=(F(1,50),F(9,10),F(1,10)):continue
        parts=conditioned_kernels(e0,e1,c)
        for q in [F(3,100),F(6,100)]:
            verts=polygon(p,q)
            rec={'q':str(q),'e0':str(e0),'e1':str(e1),'c':str(c),'vertices':[[str(x),str(y)] for x,y in verts],'edges':[]}
            for k,a in enumerate(verts):
                b=verts[(k+1)%len(verts)]
                matrices,D=controls_integer([a,b],parts,p,q)
                rec['edges'].append(edge_certificate(matrices,D))
            if args.vertical:
                xl,xh=min(x for x,y in verts),max(x for x,y in verts)
                yl,yh=min(y for x,y in verts),max(y for x,y in verts)
                pts=[(xl,yl),(xh,yl),(xl,yh),(xh,yh)]
                matrices,D=controls_integer(pts,parts,p,q)
                rec['vertical']=vertical_certificate(matrices,D,yh-yl)
            out['cases'].append(rec)
            print('certified',q,e0,e1,c,'edges',len(verts),'vertical',args.vertical,'seconds',round(time.time()-start,2),flush=True)
    out['runtime_seconds']=time.time()-start
    out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    name='independent_certificates'+('_vertical' if args.vertical else '_edges')+('_reference' if args.reference_only else '')+'.json'
    Path(__file__).with_name(name).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
