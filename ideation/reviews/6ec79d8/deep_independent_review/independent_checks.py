#!/usr/bin/env python3
"""Independent, dependency-light checks for BWBrook/deep, review commit 6ec79d8.

This is a fresh implementation from synthesis/model_specification.md. It does
not import or execute repository modules. Exact arithmetic uses Fraction for
kernel construction and common-denominator integer propagation. The two
within-step processes are factored into mortality and recolonisation matrices.

Run: python independent_checks.py
Writes independent_checks.json beside this script. Optional exploratory
checks in other scripts depend on NumPy/SciPy and are not class certificates.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, lcm
from pathlib import Path
import hashlib
import json
import time

COMMIT = '6ec79d8c1878d1afdd4044b21af0d8b8d8de2cc3'
N = 4


def subsets(mask):
    s = mask
    while True:
        yield s
        if s == 0:
            return
        s = (s-1) & mask


def zeros(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


@lru_cache(None)
def mortality_parts(n, e0, e1, buffer_fraction=F(1)):
    """Conditional matrices old occupancy -> survivors for each hit mask."""
    size = 1 << n
    parts = []
    for hits in range(size):
        d = [e0 + (e1-e0) * bool(hits & (1 << i)) *
             (buffer_fraction if i == 0 else 1) for i in range(n)]
        M = zeros(size)
        for old in range(size):
            for surviving in subsets(old):
                prob = F(1)
                for i in range(n):
                    if old & (1 << i):
                        prob *= 1-d[i] if surviving & (1 << i) else d[i]
                M[old][surviving] = prob
        assert all(sum(row) == 1 for row in M)
        parts.append(M)
    return parts


@lru_cache(None)
def recolonisation(n, c):
    """Survivors -> final occupancy, using only original survivors as donors."""
    if not 0 <= c <= 1:
        raise ValueError('Per-donor colonisation must lie in [0,1].')
    size = 1 << n
    R = zeros(size)
    full = size - 1
    for survivors in range(size):
        vacant = full ^ survivors
        r = 1-(1-c)**survivors.bit_count()
        for added in subsets(vacant):
            k = added.bit_count()
            R[survivors][survivors | added] = r**k * (1-r)**(vacant.bit_count()-k)
    assert all(sum(row) == 1 for row in R)
    return R


def kernel(law, e0, e1, c, action='R', buffer_fraction=F(1,2), c_multiplier=F(2)):
    """Form E_Q[mortality matrix] @ recolonisation matrix, exactly."""
    n = (len(law)-1).bit_length()
    size = 1 << n
    assert len(law) == size and sum(law) == 1 and min(law) >= 0
    parts = mortality_parts(n, e0, e1, buffer_fraction if action == 'R' else F(1))
    M = zeros(size)
    for h, weight in enumerate(law):
        if weight:
            for old in range(size):
                for surviving in subsets(old):
                    M[old][surviving] += weight * parts[h][old][surviving]
    R = recolonisation(n, c * (c_multiplier if action == 'C' else 1))
    K = zeros(size)
    for old in range(size):
        for surviving in subsets(old):
            if M[old][surviving]:
                for dest, prob in enumerate(R[surviving]):
                    if prob:
                        K[old][dest] += M[old][surviving] * prob
    assert K[0][0] == 1 and all(v == 0 for v in K[0][1:])
    assert all(sum(row) == 1 and min(row) >= 0 for row in K)
    return K


def integer_matrix(K):
    D = lcm(*(v.denominator for row in K for v in row))
    return [[v.numerator * (D//v.denominator) for v in row] for row in K], D


def trajectory(K, horizon=50, start=None):
    """Exact probabilities, without recurrent rational reductions."""
    M, D = integer_matrix(K)
    size = len(K)
    start = size-1 if start is None else start
    v = [0]*size
    v[start] = 1
    out = [F(int(start != 0))]
    scale = 1
    cols = list(zip(*M))
    for t in range(horizon):
        v = [sum(a*b for a,b in zip(v,col)) for col in cols]
        scale *= D
        assert sum(v) == scale
        out.append(F(sum(v[1:]),scale))
    return out


def count_to_law(pi):
    n = len(pi)-1
    return tuple(F(pi[h.bit_count()],comb(n,h.bit_count())) for h in range(1 << n))


def moment_to_law(m):
    n = len(m)-1
    law = tuple(sum((-1)**j * comb(n-h.bit_count(), j) * m[h.bit_count()+j]
                    for j in range(n-h.bit_count()+1)) for h in range(1 << n))
    assert min(law) >= 0 and sum(law) == 1
    return law


def scalar_law(p,q,upper=False,n=4):
    if upper:
        a = (p-q)/(1-p)
        w = (q-p*p)/(1-2*p+q)
        atoms, weights = (a,F(1)), (1-w,w)
    else:
        atoms, weights = (F(0),q/p),(1-p*p/q,p*p/q)
    moments = [sum(w*a**k for a,w in zip(atoms,weights)) for k in range(n+1)]
    return moment_to_law(moments)


def solve(A,b):
    k = len(b)
    rows = [[F(x) for x in row] + [F(y)] for row,y in zip(A,b)]
    for col in range(k):
        pivot = next((r for r in range(col,k) if rows[r][col]), None)
        if pivot is None:
            return None
        rows[col], rows[pivot] = rows[pivot],rows[col]
        z = rows[col][col]
        rows[col] = [a/z for a in rows[col]]
        for r in range(k):
            if r != col:
                z = rows[r][col]
                rows[r] = [a-z*b for a,b in zip(rows[r],rows[col])]
    return [row[-1] for row in rows]


def e_vertices(p,q,n=4):
    out = set()
    for support in combinations(range(n+1),3):
        A = [[1]*3,list(support),[k*(k-1) for k in support]]
        weights = solve(A,[1,n*p,n*(n-1)*q])
        if weights is not None and min(weights) >= 0:
            pi = [F(0)]*(n+1)
            for k,w in zip(support,weights):
                pi[k] = w
            out.add(tuple(pi))
    return sorted(out)


def information(law):
    n = (len(law)-1).bit_length()
    def cohit(inds):
        mask = sum(1 << i for i in inds)
        return sum(w for h,w in enumerate(law) if h & mask == mask)
    return [[cohit(inds) for inds in combinations(range(n),k)] for k in range(1,n+1)]


def pack(x):
    return {'exact':str(x),'value':float(x)}


def evaluate(law,e0,e1,c,baseline=False):
    values = {}
    traces = {}
    for action in (('0','R','C') if baseline else ('R','C')):
        traces[action] = trajectory(kernel(law,e0,e1,c,action))
        values[action] = traces[action][-1]
    values['Delta'] = values['R']-values['C']
    s = e1-e0
    one_step = sum(w*s/F(2)*bool(h&1)*
                   product_values(e0+s*bool(h & (1 << j)) for j in range(1,4))
                   for h,w in enumerate(law))
    assert traces['R'][1]-traces['C'][1] == one_step
    return values,traces


def product_values(values):
    out = F(1)
    for v in values:
        out *= v
    return out


def regret(a,b):
    rR,rC = max(F(0),-a),max(F(0),b)
    out = {'deterministic':pack(min(rR,rC)), 'choice':'R' if rR<rC else 'C' if rC<rR else 'tie'}
    if a < 0 < b:
        out.update(lottery=pack(-a*b/(b-a)),pR=pack(b/(b-a)))
    return out


def main():
    started = time.time()
    p = F(1,10)
    output = {'commit':COMMIT,'implementation':'independent; no repository imports',
              'arithmetic':'exact rational kernels; integer propagation',
              'scope':'attained values only; does not by itself certify class extrema',
              'design':[],'reference_laws':[]}
    for q in [F(3,100),F(6,100)]:
        vertices = e_vertices(p,q)
        for e0,e1,c in product([F(1,50),F(1,10)],[F(1,2),F(9,10)],[F(1,50),F(1,10)]):
            reference = (e0,e1,c) == (F(1,50),F(9,10),F(1,10))
            E,S = [],[]
            for pi in vertices:
                law = count_to_law(pi)
                inf = information(law)
                assert set(inf[0]) == {p} and set(inf[1]) == {q}
                values,traces = evaluate(law,e0,e1,c,baseline=reference)
                E.append(values['Delta'])
                if reference:
                    output['reference_laws'].append({'q':str(q),'class':'E','pi':list(map(str,pi)),
                        'information':[[str(x) for x in v] for v in inf],
                        'outcomes':{k:pack(v) for k,v in values.items()},
                        'trajectory':{k:[float(x) for x in v] for k,v in traces.items()},
                        'one_step_delta':pack(traces['R'][1]-traces['C'][1])})
            for upper in [False,True]:
                law = scalar_law(p,q,upper)
                values,traces = evaluate(law,e0,e1,c,baseline=reference)
                S.append(values['Delta'])
                if reference:
                    pi = [sum(law[h] for h in range(16) if h.bit_count()==k) for k in range(5)]
                    output['reference_laws'].append({'q':str(q),'class':'S','name':'Ghi' if upper else 'Glo',
                        'pi':list(map(str,pi)), 'outcomes':{k:pack(v) for k,v in values.items()},
                        'trajectory':{k:[float(x) for x in v] for k,v in traces.items()}})
            output['design'].append({'q':str(q),'e0':str(e0),'e1':str(e1),'c':str(c),
                'E_vertex_range':list(map(pack,[min(E),max(E)])),
                'S_endpoint_range':list(map(pack,[min(S),max(S)])),
                'E_pair_regret':regret(min(E),max(E)),'S_pair_regret':regret(min(S),max(S))})
            print('checked',q,e0,e1,c, 'E',float(min(E)),float(max(E)),flush=True)
    # The published feasible unimodal witness and the containing-envelope point.
    moments=[F(1),p,F(3,100),F(68381,3464000),F(14621110637,937445000000)]
    vals,_=evaluate(moment_to_law(moments),F(1,50),F(9,10),F(1,10))
    output['unimodal_witness']={k:pack(v) for k,v in vals.items()}
    r=F(79,4000)
    ymax=r-(F(3,100)-r)**2/(p-F(3,100))
    vals,_=evaluate(moment_to_law([F(1),p,F(3,100),r,ymax]),F(1,50),F(9,10),F(1,10))
    assert vals['Delta'] < -F(427,100000)
    output['unimodal_envelope_point']={k:pack(v) for k,v in vals.items()}
    output['runtime_seconds']=time.time()-started
    output['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dest=Path(__file__).with_suffix('.json')
    dest.write_text(json.dumps(output,indent=2)+'\n')
    print('Wrote',dest,'seconds',output['runtime_seconds'])


if __name__=='__main__':
    main()
