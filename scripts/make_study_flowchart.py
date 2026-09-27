#!/usr/bin/env python3
"""
Figure 1: study overview flowchart (matplotlib, classic monochrome, serif, real math typesetting).

    question -> closed-form two-arm quota -> verified four ways -> one from-raw embedding
    -> two parallel validation columns (chemical: Tahoe-100M -> EmeraldBay ; genetic: Orion -> TRADE)

Black on white, thin rules, Liberation Serif (metric clone of Times New Roman), italic variables via
mathtext (STIX). Numbers mirror docs/SUPPLEMENT.md. No Graphviz dependency.

Run:  python scripts/make_study_flowchart.py  ->  figures/fig1_study_flowchart.{pdf,png}
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

plt.rcParams.update({
    "font.family": "Liberation Serif",
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "text.color": "black", "axes.edgecolor": "black",
})

W, H = 7.2, 8.5                      # inches; axes use inch coordinates directly
LW = 0.8                              # rule weight
GAP = 0.06                            # arrowhead standoff


class Canvas:
    def __init__(self):
        self.fig = plt.figure(figsize=(W, H))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, W); self.ax.set_ylim(0, H); self.ax.set_aspect("equal"); self.ax.axis("off")

    # ---- boxes -----------------------------------------------------------------
    def box(self, cx, cy, w, h, title, lines=(), title_pt=10.5, body_pt=8.6, lw=LW, title_gap=0.155):
        """Rectangular node: bold title, then centred body lines."""
        self.ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h, boxstyle="square,pad=0",
                                         fc="white", ec="black", lw=lw, zorder=2))
        n = len(lines)
        step = 0.155 if body_pt <= 9 else 0.19
        block = title_gap + step * max(n - 1, 0)      # title baseline to last body baseline
        y = cy + block / 2
        self.ax.text(cx, y, title, ha="center", va="center", fontsize=title_pt, fontweight="bold", zorder=3)
        for i, ln in enumerate(lines):
            self.ax.text(cx, y - title_gap - step * i, ln, ha="center", va="center", fontsize=body_pt, zorder=3)
        return dict(cx=cx, cy=cy, w=w, h=h, top=cy + h / 2, bot=cy - h / 2, l=cx - w / 2, r=cx + w / 2)

    def group(self, x0, y0, x1, y1, title, pt=8.8, align="left"):
        """Thin grouping rectangle with a small bold title set into its top rule."""
        self.ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="white", ec="black", lw=LW, zorder=1))
        x = x0 + 0.12 if align == "left" else x1 - 0.12
        self.ax.text(x, y1, "  " + title + "  ", ha=align, va="center", fontsize=pt, fontweight="bold",
                     bbox=dict(fc="white", ec="none", pad=1.0), zorder=3)

    # ---- connectors ---------------------------------------------------------------
    def arrow(self, x0, y0, x1, y1, label=None, lx=0.08, ly=0.0, ha="left"):
        self.ax.annotate("", xy=(x1, y1 + GAP if y1 < y0 else y1 - GAP), xytext=(x0, y0),
                         arrowprops=dict(arrowstyle="-|>", lw=LW, color="black", mutation_scale=9,
                                         shrinkA=0, shrinkB=0), zorder=1)
        if label:
            self.ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, fontsize=8, ha=ha, va="center", zorder=3)

    def elbow(self, x0, y0, x1, y1, ymid, label=None, side="right"):
        """Orthogonal route: down from (x0,y0) to ymid, across to x1, down to y1 (arrowhead).
        The label sits just above the horizontal segment, next to the corner at x1."""
        self.ax.plot([x0, x0, x1], [y0, ymid, ymid], color="black", lw=LW, zorder=1, solid_capstyle="round")
        self.arrow(x1, ymid, x1, y1)
        if label:
            if side == "right":
                self.ax.text(x1 + 0.10, ymid + 0.09, label, fontsize=8, ha="left", va="center", zorder=3)
            else:
                self.ax.text(x1 - 0.10, ymid + 0.09, label, fontsize=8, ha="right", va="center", zorder=3)

    def save(self, stem):
        for fmt in ("pdf", "png"):
            self.fig.savefig(f"{stem}.{fmt}", dpi=600 if fmt == "png" else None, facecolor="white")
            print(f"wrote {stem}.{fmt}")


def build():
    c = Canvas()
    cx = W / 2

    # 1. question ---------------------------------------------------------------------
    q = c.box(cx, 8.10, 5.9, 0.62, "Question",
              [r"How many cells resolve a perturbation's $\it{direction}$ to an angular tolerance $\theta_\star$?"],
              body_pt=9.4, title_gap=0.19)

    # 2. closed-form quota ------------------------------------------------------------
    t = c.box(cx, 6.94, 5.9, 1.22, "Closed-form two-arm cell quota", body_pt=9.0, title_gap=0.34,
              lines=[r"$n_t^{\star} \;=\; \dfrac{1}{\,m^{2}\theta_\star^{2}\,/\,\mathrm{tr}(P\Sigma P)\;-\;1/n_c\,}$",
                     "",
                     r"only noise perpendicular to the effect rotates it $\;\cdot\;$ control-pool floor "
                     r"$m_{\min}=\sqrt{\mathrm{tr}(P\Sigma P)/(n_c\theta_\star^{2})}$, below it $n_t^{\star}=\infty$"])
    # nudge: the formula line is taller than a text line; the "" spacer above absorbs it.
    c.arrow(cx, q["bot"], cx, t["top"])

    # 3. verification strip -----------------------------------------------------------
    gy0, gy1 = 5.24, 6.06
    c.group(0.45, gy0, W - 0.45, gy1, "Verified four independent ways")
    chips = [("Symbolic Jacobian", "SymPy, residual 0"),
             ("Monte-Carlo", r"$K=50{,}000$, error 0.13%"),
             ("Independent", "implementation"),
             ("Lean 4 / Mathlib", "deterministic core")]
    cw, ch, xs = 1.45, 0.46, [1.25, 2.83, 4.37, 5.95]
    for (ttl, sub), x in zip(chips, xs):
        c.box(x, gy0 + 0.34, cw, ch, ttl, [sub], title_pt=9.2, body_pt=8.0, title_gap=0.16)
    c.arrow(cx, t["bot"], cx, gy1)

    # 4. shared embedding -------------------------------------------------------------
    e = c.box(cx, 4.64, 5.9, 0.68, "One from-raw embedding for every atlas", body_pt=8.8, title_gap=0.19,
              lines=[r"normalize to $10^{4}$ $\rightarrow$ $\log(1+x)$ $\rightarrow$ 2,000 HVG $\rightarrow$ PCA(50); "
                     r"within-condition $\Sigma$ from single cells; real control pool $n_c$"])
    c.arrow(cx, gy0, cx, e["top"])

    # 5. two validation columns -------------------------------------------------------
    oy0, oy1 = 0.28, 3.80
    c.group(0.30, oy0, W - 0.30, oy1, "Four public atlases, six screens, two modalities, three platforms")
    lx0, lx1 = 0.52, 3.50
    rx0, rx1 = 3.70, W - 0.52
    iy1 = oy1 - 0.34
    c.group(lx0, oy0 + 0.16, lx1, iy1, "Chemical (Vevo Mosaic)", align="right")
    c.group(rx0, oy0 + 0.16, rx1, iy1, "Genetic (CRISPRi)", align="right")
    lcx, rcx = (lx0 + lx1) / 2, (rx0 + rx1) / 2
    bw, bh = 2.62, 0.86
    ytop, ybot = 2.74, 1.08

    ta = c.box(lcx, ytop, bw, bh, "Tahoe-100M",
               ["56,827 conditions, 95.6 M cells", r"$\sigma^{2}=0.957$, held-out slope 0.94",
                "shared DMSO vehicle: 55% control-pool-limited"])
    eb = c.box(lcx, ybot, bw, bh, "EmeraldBay",
               ["4,912 conditions, 1.83 M cells", r"$\sigma^{2}=0.917$, held-out slope 0.98",
                "large per-line pool: 83% treated-depth-limited"])
    orn = c.box(rcx, ytop, bw, bh, "X-Atlas/Orion (HCT116, HEK293T)",
                ["35,441 knockdowns, about 8 M cells", r"$\sigma^{2}=0.91\,/\,1.03$",
                 "genome-wide: 15 to 35% detectable, none deep enough"])
    tr = c.box(rcx, ybot, bw, bh, "TRADE (Jurkat, HepG2)",
               ["4,070 essential-gene knockdowns", r"$\sigma^{2}=1.50\,/\,1.87$, slope 1.02 / 0.96",
                "strong effects at 48 to 85 cells per gene"])

    c.arrow(lcx, ta["bot"], lcx, eb["top"], label=r"$\sigma^{2}$ transfers across timepoints", lx=0.08)
    c.arrow(rcx, orn["bot"], rcx, tr["top"], label=r"direct test at $m>4$", lx=0.08)

    ymid = (e["bot"] + oy1) / 2
    c.elbow(cx, e["bot"], lcx, ta["top"], ymid, label="calibrate", side="right")
    c.elbow(cx, e["bot"], rcx, orn["top"], ymid, label="transfer", side="left")

    return c


if __name__ == "__main__":
    import sys
    if "--legacy" not in sys.argv:
        print("Figure 1 is produced as native Word shapes by scripts/make_fig1_word.py since 2026-09-27 "
              "(it also writes figures/fig1_study_flowchart.{docx,pdf,png}). Pass --legacy to draw the old "
              "matplotlib version.")
        sys.exit(0)
    os.makedirs("figures", exist_ok=True)
    build().save("figures/fig1_study_flowchart")
