"""Plot the four certified E extremisers and their independently audited dynamics.

uv run --python 3.12 --with matplotlib python synthesis/plot_mechanism_extremisers.py
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RECEIPT = ROOT / "ideation/b_screen/mechanism_audit.json"
OUT = HERE / "figures"
MIN_COLOR, MAX_COLOR = "#187d8d", "#b75d35"
R_COLOR, C_COLOR = "#8056a0", "#609047"


def exact(entry):
    return F(entry["exact"])


def read_cases():
    receipt = json.loads(RECEIPT.read_text())
    assert receipt["status"] == "PASS"
    for path, digest in receipt["source_hashes_sha256"].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == digest, path
    assert len(receipt["cases"]) == 4
    grouped = {}
    for row in receipt["cases"]:
        q = F(row["parameters"]["q"])
        pi = [F(x) for x in row["parameters"]["pi_by_hit_count"]]
        assert sum(pi) == 1
        assert sum(k * v for k, v in enumerate(pi)) == F(2, 5)
        assert sum(k * (k - 1) * v / 2 for k, v in enumerate(pi)) == 6 * q
        trajectory = row["true_horizon_Delta_PR_minus_PC"]
        assert [x["t"] for x in trajectory] == list(range(51))
        d = row["direct_persistence"]
        assert exact(trajectory[-1]) == exact(d["PR"]) - exact(d["PC"])
        assert exact(d["PR_minus_P0"]) == exact(d["PR"]) - exact(d["P0"])
        assert exact(d["PC_minus_P0"]) == exact(d["PC"]) - exact(d["P0"])
        for action in ("R", "C"):
            assert exact(d[f"P{action}_minus_P0"]) > 0
        grouped.setdefault(q, []).append(row)
    assert set(grouped) == {F(3, 100), F(3, 50)}
    for q, rows in grouped.items():
        rows.sort(key=lambda r: exact(r["direct_persistence"]["Delta_PR_minus_PC"]))
        assert len(rows) == 2
        assert exact(rows[0]["direct_persistence"]["Delta_PR_minus_PC"]) < 0
        assert exact(rows[1]["direct_persistence"]["Delta_PR_minus_PC"]) > 0
    return grouped


def plot(groups):
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({
        "font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "axes.labelsize": 10,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.3))
    for i, (q, rows) in enumerate(sorted(groups.items())):
        counts, horizons, gains = axes[i]
        for law, color, offset, ls in zip(rows, (MIN_COLOR, MAX_COLOR), (-.19, .19), ("-", "--")):
            pi = [float(F(x)) for x in law["parameters"]["pi_by_hit_count"]]
            counts.bar([k + offset for k in range(5)], pi, width=.36, color=color)
            ts = law["true_horizon_Delta_PR_minus_PC"]
            delta = [float(100 * exact(x)) for x in ts]
            horizons.plot(range(51), delta, color=color, lw=2.3, ls=ls)
            horizons.scatter([50], [delta[-1]], s=23, color=color, zorder=3)
            horizons.annotate(f"{delta[-1]:+.2f}", (50, delta[-1]), xytext=(5, 0),
                              textcoords="offset points", color=color, va="center", fontsize=9)
        counts.set_ylim(0, .93)
        counts.set_yticks([0, .2, .4, .6, .8])
        counts.set_xticks(range(5))
        counts.set_xlabel("Patches hit in one model step, k")
        counts.set_ylabel("Hit-count probability")
        counts.text(.40, .72, f"q = {float(q):.2f}\nExpected co-hit pairs = {float(6*q):.2f}",
                    transform=counts.transAxes, va="top", fontsize=9)

        horizons.axhline(0, color="#777777", lw=.8)
        horizons.set_xlim(0, 61)
        horizons.set_xticks([0, 10, 20, 30, 40, 50])
        horizons.set_ylim(-20, 16)
        horizons.set_xlabel("Actual horizon, t (model steps)")
        horizons.set_ylabel(r"$P_R(t)-P_C(t)$ (percentage points)")
        crossing, = rows[0]["strict_R_advantage_to_C_advantage_crossings"]
        horizons.annotate(f"C first ahead at t = {crossing}", (crossing, 0),
                          xytext=(12, -4.5 if i == 0 else -10), textcoords="data",
                          arrowprops={"arrowstyle": "-", "color": MIN_COLOR, "lw": .9},
                          color=MIN_COLOR, fontsize=9)

        for j, row in enumerate(rows):
            d = row["direct_persistence"]
            for action, offset, color in (("R", -.19, R_COLOR), ("C", .19, C_COLOR)):
                value = float(100 * exact(d[f"P{action}_minus_P0"]))
                gains.bar(j + offset, value, width=.36, color=color)
                gains.text(j + offset, value + .7, f"{value:.2f}", ha="center", fontsize=9)
        gains.set_xlim(-.6, 1.6)
        gains.set_ylim(0, 37)
        gains.set_yticks([0, 10, 20, 30])
        gains.set_xticks([0, 1], ["Minimising law", "Maximising law"])
        gains.set_ylabel("Gain above baseline at t = 50 (pp)")
        gains.set_xlabel(r"Law minimising or maximising $\Delta_{50}$")
        for ax in (counts, horizons, gains):
            ax.grid(axis="y", color="#ededed", lw=.7)
            ax.set_axisbelow(True)

    for ax, title in zip(axes[0], ("A  Pairwise agreement permits different hit counts",
                                  "B  The ranking changes along actual trajectories",
                                  "C  Connectivity's benefit falls most")):
        ax.set_title(title, fontsize=10.5, pad=14)
    axes[0, 0].legend(handles=[Patch(color=MIN_COLOR, label=r"Minimises $\Delta_{50}$"),
                              Patch(color=MAX_COLOR, label=r"Maximises $\Delta_{50}$")],
                       loc="upper right", frameon=False, fontsize=9)
    axes[0, 1].legend(handles=[Line2D([], [], color=MIN_COLOR, lw=2, label="Minimising law"),
                              Line2D([], [], color=MAX_COLOR, lw=2, ls="--", label="Maximising law")],
                       loc="lower left", frameon=False, fontsize=9)
    axes[0, 2].legend(handles=[Patch(color=R_COLOR, label=r"Buffering: $P_R-P_0$"),
                              Patch(color=C_COLOR, label=r"Connectivity: $P_C-P_0$")],
                       loc="upper right", frameon=False, fontsize=9)
    fig.suptitle("Preserving survivors and rebuilding occupancy respond differently to higher-order hits",
                 fontsize=14, y=.985)
    fig.text(.5, .943, "Four certified E extremisers · n = 4, p = .10, e₀ = .02, e₁ = .90, c = .10 · initially all occupied",
             ha="center", fontsize=10)
    fig.text(.5, .045, "Within each row, both laws preserve all pairwise hit tables. A hit is an exposure, not necessarily a death.",
             ha="center", fontsize=10)
    fig.text(.5, .020, "B uses direct horizon calculations; C uses complete-action gains. Neither shows a prefix of the fixed-T telescoping allocation.",
             ha="center", fontsize=9)
    fig.subplots_adjust(left=.065, right=.985, top=.855, bottom=.14, hspace=.48, wspace=.37)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"mechanism_extremisers.{suffix}", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    plot(read_cases())
    print("Rendered four audited laws, actual horizon trajectories and direct baseline contrasts.")
