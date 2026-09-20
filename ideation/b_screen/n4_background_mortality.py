"""Bounded n4 e0 screen and exact univariate interval certificates.

Run: uv run --python 3.12 --with gmpy2 python ideation/b_screen/n4_background_mortality.py
Fixed iid annual mixing laws: p=.1, q=.03. Both actions, initial occupancy,
and T=50 are inherited unchanged. Screen counts are neither area nor absence proofs.
The exact kernel has degree <=4 in e0 at fixed e1,c. Bernstein entry bounds
and probability-simplex propagation enclose every e0 in each retained leaf.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import time
from gmpy2 import mpq as Q
import n4_rectangle_certificate as R

N = R.N3
LAWS = N.laws(F(1, 10), F(3, 100))
ROOT = Path(__file__).resolve().parent


def delta_float(e0, e1, c, law):
    atoms, weights = law
    result = []
    for action in ('R', 'C'):
        K = N.averaged_kernel(action, 4, e0, e1, c,
                              list(map(float, atoms)), list(map(float, weights)))
        v = [0.] * len(K)
        v[7 if action == 'R' else 4] = 1.
        for _ in range(50):
            v = N.matvec(v, K)
        result.append(sum(v[1:]))  # avoid cancellation near extinction
    return result[0] - result[1]


def bounds(action, interval, e1, c, law):
    lo, hi = interval
    samples = [R.avg_kernel(action, e1, c, lo + F(j, 4)*(hi-lo), *law)
               for j in range(5)]
    size = len(samples[0])
    lower, upper = [], []
    for r in range(size):
        lr, ur = [], []
        for s in range(size):
            coeff = [sum(R.Minv[k][j] * R.qq(samples[j][r][s])
                         for j in range(5)) for k in range(5)]
            lr.append(max(Q(0), min(coeff)))
            ur.append(min(Q(1), max(coeff)))
        lower.append(lr)
        upper.append(ur)
    return lower, upper


def delta_interval(interval, e1, c, law):
    rlo, rhi = R.propagate(*bounds('R', interval, e1, c, law), 7)
    clo, chi = R.propagate(*bounds('C', interval, e1, c, law), 4)
    return rlo-chi, rhi-clo


def classify(ds):
    if ds[0][1] < 0 < ds[1][0]:
        return 'reversal'
    if all(x[1] < 0 for x in ds):
        return 'C_for_both_endpoint_laws'
    if all(x[0] > 0 for x in ds):
        return 'R_for_both_endpoint_laws'
    return 'unresolved'


def record(interval, ds):
    return dict(e0=[str(x) for x in interval], classification=classify(ds),
                delta_intervals=[[str(x) for x in d] for d in ds],
                delta_display=[[float(x) for x in d] for d in ds])


def certify(interval, e1, c, maxdepth=14):
    stack = [(interval, 0)]
    leaves = []
    evaluations = 0
    while stack:
        box, depth = stack.pop()
        ds = [delta_interval(box, e1, c, law) for law in LAWS]
        evaluations += 1
        if classify(ds) != 'unresolved' or depth == maxdepth:
            leaves.append(record(box, ds))
        else:
            mid = sum(box)/2
            stack.extend([((mid, box[1]), depth+1), ((box[0], mid), depth+1)])
    leaves.sort(key=lambda x: F(x['e0'][0]))
    assert F(leaves[0]['e0'][0]) == interval[0]
    assert F(leaves[-1]['e0'][1]) == interval[1]
    assert all(a['e0'][1] == b['e0'][0] for a, b in zip(leaves, leaves[1:]))
    return dict(e1=str(e1), c=str(c), e0=[str(x) for x in interval],
                maxdepth=maxdepth, evaluations=evaluations, leaves=leaves)


def screen():
    rows = []
    for i in range(16):
        e0 = i/50
        cells = []
        for j in range(4, 11):
            for k in range(1, 51):
                e1, c = j/10, k/100
                ds = [delta_float(e0, e1, c, law) for law in LAWS]
                cells.append(dict(e1=e1, c=c, delta=ds,
                                  margin=min(-ds[0], ds[1])))
        rows.append(dict(e0=e0, tested=len(cells),
                         reversal_count=sum(x['margin'] > 1e-10 for x in cells),
                         best=max(cells, key=lambda x: x['margin']), cells=cells))
    return dict(e0_grid='i/50, i=0..15', e1_grid='j/10, j=4..10',
                c_grid='k/100, k=1..50', strict_margin_threshold=1e-10, rows=rows)


def checks():
    # Degree-four reconstruction at a node not used for fitting, all matrix entries.
    interval, e1, c = (F(1, 25), F(7, 50)), F(9, 10), F(9, 20)
    u = F(2, 7)
    from math import comb
    for law in LAWS:
        for action in ('R', 'C'):
            samples = [R.avg_kernel(action, e1, c, interval[0]+F(j,4)*(interval[1]-interval[0]), *law) for j in range(5)]
            direct = R.avg_kernel(action, e1, c, interval[0]+u*(interval[1]-interval[0]), *law)
            for r, row in enumerate(direct):
                assert sum(row) == 1 and all(x >= 0 for x in row)
                for s, value in enumerate(row):
                    coeff = [sum(R.Minv[k][j]*R.qq(samples[j][r][s]) for j in range(5)) for k in range(5)]
                    assert sum(coeff[k]*comb(4,k)*u**k*(1-u)**(4-k) for k in range(5)) == value
        exact = delta_interval((F(1,10),)*2, e1, c, law)
        assert exact[0] == exact[1]
        assert abs(float(exact[0])-delta_float(.1, float(e1), float(c), law)) < 1e-12
    return 'Exact off-node polynomial reconstruction, stochastic kernels, point collapse and float comparison passed.'


def main():
    t0 = time.monotonic()
    result = dict(n=4, p='1/10', q='3/100', T=50, checks=checks(), screen=screen())
    result['certified_slices'] = []
    # Original slice and a second slice selected after the exploratory screen.
    for interval, e1, c in [((F(0),F(3,100)),F(4,5),F(3,20)),
                             ((F(1,25),F(7,50)),F(9,10),F(9,20))]:
        report = certify(interval,e1,c)
        assert all(x['classification']=='reversal' for x in report['leaves'])
        result['certified_slices'].append(report)
        print('Certified', report['e0'], 'at', report['e1'], report['c'], 'leaves', len(report['leaves']), flush=True)
    result['exact_points'] = []
    for e0,e1,c in [(F(1,50),F(9,10),F(9,20)), (F(1,10),F(9,10),F(9,20)),
                     (F(4,25),F(9,10),F(9,20)),
                     (F(1,25),F(4,5),F(3,20)), (F(1,5),F(1),F(1,2))]:
        entry=record((e0,e0),[delta_interval((e0,e0),e1,c,law) for law in LAWS])
        entry.update(e1=str(e1),c=str(c))
        result['exact_points'].append(entry)
    result['source_sha256']={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                            for name in (Path(__file__).name,'n4_rectangle_certificate.py','n3_certificate.py')}
    result['seconds']=time.monotonic()-t0
    out=ROOT/'n4_background_mortality.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(out, 'seconds',result['seconds'],flush=True)


if __name__ == '__main__':
    main()
