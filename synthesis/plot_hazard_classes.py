"""Plot exact n4 class geometry; signs refer to the audited labelled witnesses.

uv run --python 3.12 --with numpy --with matplotlib python synthesis/plot_hazard_classes.py
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'figures'
P = F(1, 10)
CASES = [
    (F(3, 100), [(F(0), F(0)), (F(3, 200), F(0)), (F(3, 100), F(3, 100))], 'Both rankings occur in S'),
    (F(3, 50), [(F(2, 75), F(0)), (F(3, 100), F(0)), (F(3, 50), F(3, 50)), (F(1, 25), F(1, 25))], 'R wins throughout S'),
]


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 5.2), layout='constrained')
    for ax, (q, vertices, claim) in zip(axs, CASES):
        # Enumerate all feasible intersections of count-probability zero lines.
        lines = [(-4, 1, 1-4*P+6*q), (12, -4, 4*P-12*q),
                 (-12, 6, 6*q), (4, -4, F(0)), (0, 1, F(0))]
        found = set()
        for (a, b, c), (d, e, f) in combinations(lines, 2):
            det = a*e-b*d
            if det:
                x, y = (-c*e+b*f)/det, (-a*f+c*d)/det
                if all(ai*x+bi*y+ci >= 0 for ai, bi, ci in lines):
                    found.add((x, y))
        assert found == set(vertices)
        for x, y in vertices:
            pi = [1-4*P+6*q-4*x+y, 4*P-12*q+12*x-4*y, 6*q-12*x+6*y, 4*x-4*y, y]
            assert all(z >= 0 for z in pi) and sum(pi) == 1
        poly = Polygon(np.array(vertices, dtype=float), facecolor='#eeeeee', edgecolor='#555555', lw=1.3)
        ax.add_patch(poly)
        lo, hi = q*q/P, q-(P-q)**2/(1-P)
        x = np.linspace(float(lo), float(hi), 1001)
        lower = (x*x-2*float(P*q)*x+float(q**3))/float(q-P*P)
        upper = x-(float(q)-x)**2/float(P-q)
        ax.fill_between(x, lower, upper, color='#a9cedf', alpha=1)
        ax.plot(x, lower, color='#356b83', lw=1.5)
        ax.plot(x, upper, color='#356b83', lw=1.5)
        low = vertices[0]
        high = (q, q)
        ax.scatter(*map(float, low), s=85, marker='v', color='#0072B2', zorder=4, clip_on=False)
        ax.scatter(*map(float, high), s=85, marker='^', color='#D55E00', zorder=4)
        ax.annotate('C witness', tuple(map(float, low)), xytext=(6, 9), textcoords='offset points', color='#006396')
        ax.annotate('R witness', tuple(map(float, high)), xytext=(-6, -20), ha='right', textcoords='offset points', color='#af4d00')
        ax.text(.03, .95, claim, transform=ax.transAxes, va='top', weight='bold', fontsize=10)
        span = float(q)-float(low[0])
        ax.set_xlim(float(low[0])-span*.07, float(q)+span*.07)
        ax.set_ylim(-float(q)*.05, float(q)*1.15)
        ax.set_title(f'p = 0.10, q = {float(q):.2f}', pad=12)
        ax.set_xlabel('Triple-hit probability, τ₃')
        ax.set_ylabel('Four-hit probability, τ₄')
        ax.grid(color='white', lw=.6)
        ax.set_axisbelow(True)
    handles = [Patch(facecolor='#eeeeee', edgecolor='#555555', label='E: all symmetric four-patch hit laws'),
               Patch(facecolor='#a9cedf', edgecolor='#356b83', label='S: scalar-hazard mixtures')]
    fig.legend(handles=handles, loc='outside lower center', ncols=2, frameon=False)
    fig.suptitle('Identical pairwise hit laws permit different higher-order structures\n'
                 'n = 4; e₀ = 0.02, e₁ = 0.90, c = 0.10; T = 50 model steps', fontsize=12)
    for suffix in ('png', 'pdf'):
        fig.savefig(OUT/f'hazard_classes.{suffix}', dpi=220)
    plt.close(fig)


if __name__ == '__main__':
    main()
