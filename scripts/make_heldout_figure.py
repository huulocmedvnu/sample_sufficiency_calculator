"""
Per-dataset held-out falsification figure: one panel per screen, each dataset INDEPENDENTLY, on ALL its
N>=400 groups. Plots the realized/predicted slope RATIO against 1/rho^2 (log x). The parameter-free law
predicts ratio -> 1 at high signal-to-noise (small 1/rho^2) and a monotone deficit as 1/rho^2 grows
(the second-order low-SNR breakdown). A dashed line marks ratio = 1. No shared-cell-line showcase.

Reads the committed per-group held-out fixtures (see make_heldout_summary.py). Orion panels render only
after scripts/orion_recompute/pass4_heldout_stream.py has written orion_<LINE>_heldout.json.

Output -> figures/fig5_heldout_per_dataset.{pdf,png}
"""
import os, json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

FX = os.path.join(os.path.dirname(__file__), "..", "fixtures")
FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
NMIN = 400
GREEN, VERM, GREY, BLACK = "#0F9D58", "#D55E00", "#9AA0A6", "#202124"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "figure.facecolor": "white",
                     "savefig.dpi": 300, "savefig.bbox": "tight", "axes.edgecolor": BLACK})


def load(name):
    p = os.path.join(FX, name)
    return json.load(open(p)) if os.path.exists(p) else None


def rows_of(kind):
    """Return list of (rho2, ratio) for N>=400 groups, or None if the fixture is absent."""
    if kind == "tahoe":
        d = load("tahoe_direct_curves_full.json"); rows = d["per_group"] if d else None
    elif kind == "emb":
        d = load("emeraldbay_falsification_full.json"); rows = d["per_group_slope"] if d else None
    elif kind in ("HCT116", "HEK293T"):
        d = load(f"orion_{kind}_heldout.json"); rows = d["per_group"] if d else None
    else:
        d = load(f"trade_{kind}_falsification.json"); rows = d["rows"] if d else None
    if rows is None:
        return None
    return [r for r in rows if r.get("N", 0) >= NMIN]


PANELS = [("tahoe", "Tahoe-100M (chemical)"), ("emb", "EmeraldBay (chemical)"),
          ("HCT116", "Orion HCT116 (genetic)"), ("HEK293T", "Orion HEK293T (genetic)"),
          ("jurkat", "TRADE Jurkat (genetic)"), ("hepg2", "TRADE HepG2 (genetic)")]

fig, axes = plt.subplots(2, 3, figsize=(11, 6.4))
for ax, (kind, title) in zip(axes.ravel(), PANELS):
    rows = rows_of(kind)
    if not rows:
        ax.text(0.5, 0.5, "re-stream running" if kind in ("HCT116", "HEK293T") else "n/a",
                ha="center", va="center", color=GREY, fontsize=11, transform=ax.transAxes)
        ax.set_title(title); ax.set_xticks([]); ax.set_yticks([]); continue
    rho2 = np.array([r["rho2"] for r in rows], float)
    ratio = np.array([r["ratio"] for r in rows], float)
    x = 1.0 / np.maximum(rho2, 1e-3)
    inr = rho2 >= 3
    ax.axhline(1.0, color=BLACK, lw=1, ls="--", zorder=1)
    ax.axvline(1 / 3, color=GREY, lw=0.8, ls=":", zorder=1)   # rho2 = 3 boundary
    ax.scatter(x[inr], ratio[inr], s=14, c=GREEN, alpha=0.6, edgecolor="none", label="ρ²≥3 (in-regime)", zorder=3)
    ax.scatter(x[~inr], ratio[~inr], s=10, c=VERM, alpha=0.35, edgecolor="none", label="ρ²<3", zorder=2)
    sp = spearmanr(x, 1 - ratio).correlation if len(rows) > 2 else np.nan
    if inr.any():
        reg = f"ρ²≥3 median {np.median(ratio[inr]):.2f} ({int(inr.sum())} in regime)"
    else:
        reg = f"all ρ²<3: max ρ²={rho2.max():.1f}"    # no group reaches the first-order regime
    ax.set_xscale("log"); ax.set_ylim(-0.05, 1.5)
    ax.set_title(f"{title}   ({len(rows):,} groups N≥400)\n{reg},  Spearman={sp:.2f}", fontsize=7.4, pad=3)
    ax.set_xlabel("1 / ρ²  (low SNR →)", fontsize=8)
    ax.set_ylabel("realized / predicted slope", fontsize=8)
    if kind == "tahoe":
        ax.legend(fontsize=6.5, loc="lower left", framealpha=0.9)

fig.suptitle("Held-out falsification, each dataset independently on all its N≥400 groups: the parameter-free "
             "slope ratio → 1\nat high signal-to-noise, and the genome-wide Orion screens are entirely below "
             "the first-order regime (all ρ²<3)", fontsize=9.5, y=1.03)
fig.tight_layout(rect=(0, 0, 1, 0.99))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig5_heldout_per_dataset.{ext}"))
print("wrote figures/fig5_heldout_per_dataset.{pdf,png}")
