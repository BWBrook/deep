#!/usr/bin/env python3
"""Join independent exact certificates to attained values and reported tables.
Run after the other three scripts. Standard library only.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
B=Path(__file__).resolve().parent

def load(name):return json.loads((B/name).read_text())
def key(r):return tuple(r[k] for k in ['q','e0','e1','c'])
def classify(bounds):
    a,b=map(F,bounds)
    return 'R' if a>0 else 'C' if b<0 else 'unidentified'

def main():
    direct=load('independent_checks.json');edge=load('independent_certificates_vertical.json');extra=load('additional_checks.json')
    D={key(r):r for r in direct['design']};E={key(r):r for r in edge['cases']};S={key(r):r for r in extra['scalar_cases']}
    assert len(D)==len(E)==len(S)==16 and D.keys()==E.keys()==S.keys()
    counts={c:{'R':0,'C':0,'unidentified':0} for c in ['S','E']}
    for k,d in D.items():
        e=E[k];s=S[k]
        ends=[F(x) for r in e['edges'] for x in r['endpoint_deltas']]
        expected=[F(x['exact']) for x in d['E_vertex_range']]
        assert [min(ends),max(ends)]==expected
        assert e['vertical']['all_coefficients_strictly_positive']
        assert all(r['derivative_sign'] in ['positive','negative'] for r in e['edges'])
        vals=[F(r[x]) for r in s['boundaries'] for x in ['min','max']]
        assert [min(vals),max(vals)]==[F(x['exact']) for x in d['S_endpoint_range']]
        assert all(r['all_coefficients_between_endpoint_values'] for r in s['boundaries'])
        counts['E'][classify(expected)]+=1
        counts['S'][classify([min(vals),max(vals)])]+=1
    assert counts=={'S':{'R':3,'C':12,'unidentified':1},'E':{'R':2,'C':12,'unidentified':2}}
    # Published table in synthesis/class_identification_table.md at the commit.
    # Entries are rounded percentage points, NOT exact enclosures.
    rows='''
.03 .02 .5 .02 -6.95960 -6.27160 -7.29636 -6.04025
.03 .02 .5 .1 -8.76850 -5.56242 -9.73029 -4.64800
.03 .02 .9 .02 .38225 1.03521 .07904 1.25002
.03 .02 .9 .1 -12.05203 3.20927 -17.92519 6.43392
.03 .1 .5 .02 -2.04517 -2.03348 -2.05156 -2.02946
.03 .1 .5 .1 -31.04345 -29.58085 -31.72573 -29.10718
.03 .1 .9 .02 -.19751 -.18842 -.20231 -.18534
.03 .1 .9 .1 -15.02100 -10.25591 -17.27362 -8.93188
.06 .02 .5 .02 -3.50174 -2.57338 -3.93329 -2.50590
.06 .02 .5 .1 -5.98947 -1.74084 -8.49499 -1.50722
.06 .02 .9 .02 2.68465 3.49793 2.26760 3.55289
.06 .02 .9 .1 3.10145 11.04449 -5.34234 11.24621
.06 .1 .5 .02 -1.83207 -1.81506 -1.83948 -1.81375
.06 .1 .5 .1 -26.00483 -24.11382 -26.92652 -23.98269
.06 .1 .9 .02 -.10174 -.08786 -.10810 -.08683
.06 .1 .9 .1 -5.25686 -1.96467 -7.40774 -1.78926
'''
    checked=0
    for line in rows.strip().splitlines():
        vals=list(map(F,line.split()));d=D[tuple(str(v) for v in vals[:4])]
        computed=[F(r['exact'])*100 for name in ['S_endpoint_range','E_vertex_range'] for r in d[name]]
        assert all(abs(a-b)<=F(5,10**6) for a,b in zip(computed,vals[4:])),line
        checked+=4
    assert extra['reference_maximin'][0]['maximin_action']=='R'
    assert extra['reference_maximin'][0]['minimax_regret']['choice']=='C'
    for name,expected in extra['source_sha256'].items():assert hashlib.sha256((B/name).read_bytes()).hexdigest()==expected
    report={'all_checks_passed':True,'counts':counts,'rounded_table_endpoints_matched':checked,
        'exact_certificate_to_attained_range_joins':32,'new_maximin_contrast_confirmed':True,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (B/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
