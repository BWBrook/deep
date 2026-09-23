"""Exact local robustness witness with nonexchangeable shocks and unchanged p,q.

Uses a path-coupling bound and the existing n4 endpoint certificates; does not
apply the exchangeable occupancy-count kernels to nonexchangeable shocks.
Run with standard-library Python from any working directory.
"""
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    source = HERE / 'n4_boundary_results.json'
    previous = json.loads(source.read_text())
    case = next(c for c in previous['cases'] if
                (c['q'], c['e0'], c['e1'], c['c']) == ('3/100', '1/50', '9/10', '1/10'))
    dminus = F(case['curves']['lower']['Delta']['at_u0_exact'])
    dplus = F(case['curves']['upper']['Delta']['at_u1_exact'])
    assert dminus < 0 < dplus
    original_regret = min(-dminus, dplus)
    assert original_regret > F(3209, 100000)

    patterns = list(product((0, 1), repeat=4))
    laws = {'Qminus': {}, 'Qplus': {}, 'H': {}}
    eps = F(1, 4000)
    for b in patterns:
        k = sum(b)
        laws['Qminus'][b] = F(2, 3)*int(k == 0) + F(1, 3)*F(3, 10)**k*F(7, 10)**(4-k)
        laws['Qplus'][b] = F(81, 83)*F(7, 90)**k*F(83, 90)**(4-k) + F(2, 83)*int(k == 4)
        laws['H'][b] = laws['Qminus'][b] + F(1, 1600)*(-1)**sum(b[:3])
    for tag in ('minus', 'plus'):
        laws['mixed_'+tag] = {b: (1-eps)*laws['Q'+tag][b]+eps*laws['H'][b] for b in patterns}

    checks = 0
    for name, law in laws.items():
        assert sum(law.values()) == 1 and min(law.values()) >= 0
        for i, j in combinations(range(4), 2):
            for hit_i, hit_j in product((0, 1), repeat=2):
                value = sum(w for b, w in law.items() if (b[i], b[j]) == (hit_i, hit_j))
                expected = {(0, 0): F(83, 100), (0, 1): F(7, 100),
                            (1, 0): F(7, 100), (1, 1): F(3, 100)}[(hit_i, hit_j)]
                assert value == expected
                checks += 1
        if name in ('H', 'mixed_minus', 'mixed_plus'):
            assert law[(1, 0, 0, 0)] != law[(0, 0, 0, 1)]
    eta = 1-(1-eps)**50
    delta_error_bound = 2*eta
    assert delta_error_bound < F(1, 40)
    lower_regret = min(-dminus-delta_error_bound, dplus-delta_error_bound)
    safe_lower = F(709, 100000)
    assert lower_regret > safe_lower
    out = {
        'scope': 'n=4, p=1/10, q=3/100, e0=1/50, e1=9/10, c=1/10, T=50; fixed labelled actions',
        'method': 'exact pairwise-marginal checks and path-coupling inequality; no new kernel propagation',
        'epsilon': str(eps), 'H_perturbation': '(-1)^(x1+x2+x3)/1600',
        'H_minimum_pattern_probability': str(min(laws['H'].values())),
        'pairwise_pattern_equalities_checked': checks,
        'laws': {name: {''.join(map(str, b)): str(w) for b, w in law.items()} for name, law in laws.items()},
        'original_Delta_minus': str(dminus), 'original_Delta_plus': str(dplus),
        'Delta_error_bound': str(delta_error_bound), 'Delta_error_bound_float': float(delta_error_bound),
        'pair_regret_lower_bound': str(lower_regret), 'pair_regret_lower_bound_float': float(lower_regret),
        'safe_pair_regret_lower_bound': str(safe_lower),
        'dependency_replayed': False,
        'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (source, Path(__file__))},
    }
    (HERE / 'asymmetric_pairwise_witness.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'checks': checks, 'safe_pair_regret_lower_bound': str(safe_lower),
                      'actual_coupling_lower_bound_float': float(lower_regret)}, indent=2))


if __name__ == '__main__':
    main()
