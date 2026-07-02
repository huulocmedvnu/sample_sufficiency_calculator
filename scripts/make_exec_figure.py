#!/usr/bin/env python
"""Render the 'arrow and jitter' schematic for the plain-language executive summary.

Two panels: (left) few cells -> the measured direction wobbles widely, with one estimate decomposed
into harmless along-arrow jitter (blue) and harmful sideways jitter (vermillion) that rotates it;
(right) many cells -> jitter cancels and the arrow locks onto the true direction. Colorblind-safe
Okabe-Ito palette; labels placed in clear zones with thin leader lines. Writes outputs/exec_arrow_jitter.png.
"""
import os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, FancyArrowPatch

OUT = os.path.join(os.path.dirname(__file__), "..", "outputs", "exec_arrow_jitter.png")
INK = "#222222"; TRUE = "#000000"; MEAS = "#C4C4C4"
BLUE = "#0072B2"   # harmless: along the arrow (just length)
VERM = "#D55E00"   # harmful: sideways (rotates the direction)
CONE = "#0072B2"
LEAD = dict(arrowstyle="-", color="#888888", lw=0.8, shrinkA=2, shrinkB=2)

ang = np.deg2rad(30.0); L = 4.6
tip = np.array([L*np.cos(ang), L*np.sin(ang)])
u = tip/np.linalg.norm(tip); perp = np.array([-u[1], u[0]])   # perp points up-left; -perp down-right


def arrow(ax, p0, p1, color, lw=1.6, alpha=1.0, z=3, mut=14):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=mut, lw=lw,
                                 color=color, alpha=alpha, zorder=z, shrinkA=0, shrinkB=0))


def panel(ax, sigma, title, sub, decomp):
    rng = np.random.default_rng(1)
    n = 12
    para = rng.normal(0, sigma*0.4, n); ortho = rng.normal(0, sigma, n)
    tips = tip[None] + para[:, None]*u[None] + ortho[:, None]*perp[None]
    half = np.degrees(np.arctan2(2*sigma, L))
    ax.add_patch(Wedge((0, 0), L*1.13, np.degrees(ang)-half, np.degrees(ang)+half,
                       color=CONE, alpha=0.11, zorder=0))
    for t in tips:
        arrow(ax, (0, 0), t, MEAS, lw=1.0, alpha=0.75, z=1, mut=9)
    arrow(ax, (0, 0), tip, TRUE, lw=2.6, z=4, mut=16)
    ax.plot(0, 0, "o", color=INK, ms=6, zorder=5)
    ax.annotate("untreated\ncell state", (0, 0), (-1.35, -0.6), color=INK, fontsize=9,
                ha="center", va="top")
    # true-effect label: clear zone ABOVE the tip, leader to the tip
    ax.annotate("true drug effect\n(its DIRECTION = the mechanism)", tip, (1.35, 4.15),
                color=TRUE, fontsize=9.5, fontweight="bold", ha="center", va="bottom",
                arrowprops=LEAD)

    if decomp:
        foot = tip + 0.95*u                       # along-arrow endpoint
        est = foot - 1.65*perp                    # sideways offset -> DOWN-RIGHT, clear zone
        arrow(ax, tip, foot, BLUE, lw=3.0, z=6, mut=13)      # along (harmless)
        arrow(ax, foot, est, VERM, lw=3.0, z=6, mut=13)      # sideways (harmful)
        ax.plot(*est, "o", color="#6a6a6a", ms=4, zorder=6)
        ax.annotate("one noisy measurement", est, (est[0]+0.15, est[1]-0.55), color="#5a5a5a",
                    fontsize=8, ha="left", va="top", arrowprops=LEAD)
        ax.annotate("along-arrow jitter\nHARMLESS (just length)", (tip+foot)/2, (5.0, 3.35),
                    color=BLUE, fontsize=9, fontweight="bold", ha="left", va="center",
                    arrowprops=LEAD)
        ax.annotate("sideways jitter\nROTATES the arrow\n(changes the answer)", (foot+est)/2,
                    (5.05, 1.55), color=VERM, fontsize=9, fontweight="bold", ha="left",
                    va="center", arrowprops=LEAD)

    ax.set_title(title, fontsize=12.5, fontweight="bold", color=INK, pad=6)
    ax.text(0.5, -0.14, sub, transform=ax.transAxes, ha="center", va="top", fontsize=10, color=INK)
    ax.set_xlim(-1.9, 8.4); ax.set_ylim(-1.9, 4.6); ax.set_aspect("equal"); ax.axis("off")


fig, (axl, axr) = plt.subplots(1, 2, figsize=(12.5, 5.2))
panel(axl, 1.05, "A few cells", "the measured arrow could point many ways —\nseveral drugs look alike", True)
panel(axr, 0.30, "Many more cells", "the jitter averages out —\nthe arrow locks onto the true mechanism", False)
fig.suptitle("What a drug does is an ARROW;  each cell you measure adds a little JITTER",
             fontsize=14, fontweight="bold", color=INK, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.94])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white")
print("wrote", os.path.abspath(OUT))
