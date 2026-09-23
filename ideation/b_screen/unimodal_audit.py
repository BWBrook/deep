"""Astra's exact audit of the n=4 unimodal-class ranking certificate.

Run: uv run --python 3.12 --with gmpy2 --with numpy --with sympy python ideation/b_screen/unimodal_audit.py
Reuses the audited vertical certificate; regenerates upper-curve monotonicity.
"""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

import sympy as S
from gmpy2 import mpq as Q
import n4_boundary as B
import n4_vertical as V

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def propagate(kernel, start):
    kernel = [[Q(v.numerator, v.denominator) for v in row] for row in kernel]
    v = [Q(0)] * len(kernel)
    v[start] = Q(1)
    for _ in range(V.T):
        v = [sum(v[r] * kernel[r][s] for r in range(len(v))) for s in range(len(v))]
    assert sum(v) == 1
    return sum(v[1:])


def main():
    m = S.symbols('m', real=True)
    p, q = S.Rational(1, 10), S.Rational(3, 100)
    A, top = -m, 1-m
    y1 = 2*(p-m)
    y2 = 3*(q-m*m-m*y1)
    variance = S.expand(y2-y1*y1)
    slack = S.expand((top-y1)*(y1-A)-variance)
    assert variance == S.Rational(1, 20)+m/5-m*m
    assert slack == S.Rational(11, 100)-4*m/5
    a = S.factor((top*y1-y2)/(top-y1))
    y3_upper = (top+2*a)*y2-(2*a*top+a*a)*y1+a*a*top
    h = S.factor(m**3+S.Rational(3, 2)*m*m*y1+m*y2+y3_upper/4)
    h_expected = (-7500*m*m+1780*m+599)/(8000*(5*m+4))
    assert S.cancel(h-h_expected) == 0
    rational_cap = S.Rational(79, 4000)
    gap = (7500*(m-S.Rational(33, 500))**2+S.Rational(33, 100))/(8000*(5*m+4))
    assert S.cancel(rational_cap-h-gap) == 0
    mstar = S.sqrt(3)/2-S.Rational(4, 5)
    assert S.simplify(S.diff(h, m).subs(m, mstar)) == 0
    assert S.simplify(a.subs(m, mstar)) == 0
    # Derivative numerator strictly decreases for m >= 0: the positive root is unique.
    assert S.cancel(S.diff(h, m)+3*(100*m*m+160*m-11)/(64*(5*m+4)**2)) == 0
    tail_weight = S.simplify(2*(p-mstar)/(1-mstar))
    max_m3 = S.simplify(h.subs(m, mstar))
    for k, target in [(1, p), (2, q), (3, max_m3)]:
        uniform_moment = sum(mstar**j for j in range(k+1))/(k+1)
        assert S.simplify((1-tail_weight)*mstar**k+tail_weight*uniform_moment-target) == 0

    params = dict(q='3/100', e0='1/50', e1='9/10', c='1/10')
    prior = json.loads((HERE/'n4_vertical_all.json').read_text())
    cases = [v for v in prior['certificates'] if all(v[k] == val for k, val in params.items())]
    assert len(cases) == 1
    vertical = cases[0]['objectives']['Delta']
    assert vertical['status'] == 'certified pos'
    assert vertical['leaves']['neg'] == vertical['leaves']['unresolved'] == 0
    assert vertical['unresolved_boxes'] == []

    p, q, e0, e1, c = F(1, 10), F(3, 100), F(1, 50), F(9, 10), F(1, 10)
    xlo, xhi = V.xbounds(p, q)
    r = F(79, 4000)
    y = V.Ufun(r, p, q)
    assert xlo < r < xhi
    coeffs, values = {}, {}
    for action in ('R', 'C'):
        base, bx, cy = V.affine_parts(action, 4, e0, e1, c, p, q)
        start = V.start_state(action, 4)
        coeffs[action] = B.objective_coeffs(base, bx, cy, start, xlo, xhi, V.Ufun, p, q)
        kernel = V.kernel_at(base, bx, cy, r, y)
        atoms, weights = V.upper_law(r, p, q)
        assert all(w >= 0 for w in weights) and sum(weights) == 1
        assert kernel == V.avg_direct(action, 4, e0, e1, c, atoms, weights)
        values[action] = propagate(kernel, start)
        assert values[action] == B.evaluate(coeffs[action], (r-xlo)/(xhi-xlo))
    delta_coeffs = [r-c for r, c in zip(coeffs['R'], coeffs['C'])]
    degree = len(delta_coeffs)-1
    derivative = [degree*(b-a) for a, b in zip(delta_coeffs, delta_coeffs[1:])]
    assert all(v > 0 for v in derivative)
    delta = values['R']-values['C']
    safe_upper = Q(-427, 100000)
    assert delta < safe_upper
    sources = ['unimodal_audit.py', 'n4_boundary.py', 'n4_vertical.py',
               'n3_certificate.py', 'n4_vertical_all.json', 'unimodal_check.py', 'unimodal_check.json']
    output = {
        'scope': {'n': 4, 'p': str(p), **params, 'T': V.T},
        'mode_geometry': {'feasible_modes': ['0', '11/80'], 'variance': str(variance),
                          'bounded_variance_slack': str(slack), 'm3_upper_at_mode': str(h),
                          'rational_cap': str(r), 'positive_gap_identity': str(gap),
                          'maximising_mode': str(mstar), 'max_m3': str(max_m3),
                          'max_m3_float': float(max_m3), 'maximiser': 'atom at mode plus Uniform[mode,1]',
                          'uniform_weight': str(tail_weight)},
        'vertical_dependency': {'replayed': False, 'certificate': 'n4_vertical_all.json', **vertical},
        'upper_curve': {'x_interval': [str(xlo), str(xhi)], 'degree': degree,
                        'positive_derivative_coefficients': len(derivative),
                        'minimum_derivative_coefficient': str(min(derivative)),
                        'coefficient_sha256': hashlib.sha256('\n'.join(map(str, delta_coeffs)).encode()).hexdigest()},
        'evaluation': {'m3': str(r), 'm4': str(y), 'PR_exact': str(values['R']),
                       'PC_exact': str(values['C']), 'Delta_exact': str(delta), 'Delta_float': float(delta),
                       'safe_Delta_upper_bound': str(safe_upper),
                       'polynomial_equals_matrix': True, 'moment_kernel_equals_atomic_average': True},
        'source_sha256': {name: digest(HERE/name) for name in sources},
    }
    (HERE/'unimodal_audit.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'checks': 'passed', 'upper_curve_derivative_coefficients_positive': len(derivative),
                      'Delta_at_cap': float(delta), 'guaranteed_C_minus_R': '427/100000',
                      'vertical_certificate_replayed': False}, indent=2))


if __name__ == '__main__':
    main()
