#!/usr/bin/env python3
"""Manuscript Figures 2-6 in plotly (Figure 1 is native Word shapes: scripts/make_fig1_word.py).

Style: plotly's own "plotly_white" template and its default qualitative colorway, no marker or bar
outlines, no black frames. One colour per sufficiency regime, shared by Figures 4 and 5b so a single
key serves both; the three modalities of Figure 5a take later colours of the same colorway so they never
coincide with a regime colour. Figure 3 (simulation study) uses the colorway in order within each panel.

Inputs (all committed): fixtures/simulation_summary.json, fixtures/fig3_hist.json, fixtures/fig4_kde.json
(the two histogram/density fixtures keep their historical names), fixtures/unified_spectrum.json,
fixtures/tahoe_direct_curves_full.json, fixtures/emeraldbay_falsification_full.json,
fixtures/orion_{HCT116,HEK293T}_heldout.json, fixtures/trade_{jurkat,hepg2}_falsification.json.
Figure 2 is a synthetic schematic (fixed seed). Numbers printed on the figures come from the fixtures.

    python3 scripts/make_figures_plotly.py        # writes figures/fig{2,3,4,5,6}_*.{pdf,png}

Needs plotly >= 5.20, kaleido 0.2.x, numpy, scipy, pillow.
"""
import json
import math
import os

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PIL import Image
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FX = os.path.join(ROOT, "fixtures")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

C = px.colors.qualitative.Plotly          # plotly's default colorway
TEMPLATE = "plotly_white"
FONT = dict(family="Arial", size=11)
PX_PER_IN = 100                           # layout pixels per inch; PNG export at scale 6 -> 600 dpi
INK = "#2a3f5f"                           # plotly_white's own text/line colour, used for reference lines

REGIMES = ["over-sampled", "treated-depth-limited", "ghost", "control-pool-limited", "not detectable"]
RCOL = dict(zip(REGIMES, C[:5]))
RKEY = {
    "over-sampled": "over-sampled<br><span style='font-size:9px'><i>n</i>* within the acquired depth</span>",
    "treated-depth-limited": "treated-depth-limited<br><span style='font-size:9px'><i>n</i>* above the acquired depth</span>",
    "ghost": "ghost<br><span style='font-size:9px'><i>n</i>* > 50,000 treated cells</span>",
    "control-pool-limited": "control-pool-limited<br><span style='font-size:9px'><i>n</i>* = ∞, control pool too small</span>",
    "not detectable": "not detectable<br><span style='font-size:9px'>effect below the sampling floor</span>",
}
SCREENS = ["Tahoe-100M", "EmeraldBay", "Orion HCT116", "Orion HEK293T", "TRADE Jurkat", "TRADE HepG2"]


def load(name):
    with open(os.path.join(FX, name)) as fh:
        return json.load(fh)


def rgba(hex_colour, a):
    h = hex_colour.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{a})"


def sup(k):
    return "10" + "".join("⁰¹²³⁴⁵⁶⁷⁸⁹"[int(c)] for c in str(k))


def save(fig, name, w_in, h_in):
    fig.update_layout(template=TEMPLATE, font=FONT, width=int(w_in * PX_PER_IN), height=int(h_in * PX_PER_IN),
                      paper_bgcolor="white", plot_bgcolor="white")
    pdf = os.path.join(FIG, name + ".pdf")
    png = os.path.join(FIG, name + ".png")
    fig.write_image(pdf)
    fig.write_image(png, scale=6)
    im = Image.open(png)
    im.save(png, dpi=(600, 600))          # journals check the pHYs tag
    print(f"  wrote figures/{name}.{{pdf,png}}")


# ------------------------------------------------------------------------------------------------
# Figure 2: geometry of the angular error (schematic)
# ------------------------------------------------------------------------------------------------
def fig2():
    rng = np.random.default_rng(7)
    v = np.array([3.0, 0.0])
    noise = 2.6
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.04,
                        subplot_titles=("(a)  few cells per arm,  <i>n</i> = 12", "(b)  many cells per arm,  <i>n</i> = 200"))
    for col, n in ((1, 12), (2, 200)):
        sd = noise / math.sqrt(n)
        xs = v[0] + rng.normal(0, sd, 60)
        ys = v[1] + rng.normal(0, sd, 60)
        xr, yr = f"x{'' if col == 1 else col}", f"y{'' if col == 1 else col}"
        fig.add_shape(type="circle", x0=v[0] - sd, x1=v[0] + sd, y0=-sd, y1=sd, xref=xr, yref=yr,
                      fillcolor=rgba(C[0], 0.13), line=dict(width=0), layer="below")
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="markers", showlegend=False, hoverinfo="skip",
                                 marker=dict(size=4, color=C[0], opacity=0.55, line=dict(width=0))), row=1, col=col)
        # true effect
        fig.add_annotation(x=v[0], y=v[1], ax=0, ay=0, xref=xr, yref=yr, axref=xr, ayref=yr, showarrow=True,
                           arrowhead=2, arrowsize=1, arrowwidth=2.4, arrowcolor=INK, text="")
        fig.add_annotation(x=1.2, y=-0.3, xref=xr, yref=yr, text="<b><i>v</i></b>  true effect", showarrow=False,
                           font=dict(size=11), yanchor="top")
        if n == 12:
            vh = np.array([3.5, 0.7])
            foot = np.array([vh[0], 0.0])
            fig.add_annotation(x=vh[0], y=vh[1], ax=0, ay=0, xref=xr, yref=yr, axref=xr, ayref=yr, showarrow=True,
                               arrowhead=2, arrowwidth=1.6, arrowcolor=C[1], text="")
            fig.add_trace(go.Scatter(x=[v[0], foot[0]], y=[0, 0], mode="lines", line=dict(color=C[2], width=4),
                                     name="along-signal: changes length only", showlegend=False), row=1, col=1)
            fig.add_trace(go.Scatter(x=[foot[0], vh[0]], y=[0, vh[1]], mode="lines", line=dict(color=C[3], width=4),
                                     name="across-signal: rotates the direction", showlegend=False), row=1, col=1)
            th = math.atan2(vh[1], vh[0])
            t = np.linspace(0, th, 40)
            arc_r = 1.25                      # drawing radius of the theta arc (schematic only)
            fig.add_trace(go.Scatter(x=arc_r * np.cos(t), y=arc_r * np.sin(t), mode="lines", showlegend=False,
                                     line=dict(color=INK, width=1)), row=1, col=1)
            fig.add_annotation(x=1.5, y=0.2, xref=xr, yref=yr, text="<i>θ</i>", showarrow=False, font=dict(size=14))
            fig.add_annotation(x=vh[0] - 0.1, y=vh[1] + 0.1, xref=xr, yref=yr, text="<b><i>v̂</i></b>  estimate",
                               showarrow=False, xanchor="right", yanchor="bottom", font=dict(color=C[1]))
            fig.add_annotation(x=vh[0] + 0.03, y=vh[1] / 2, xref=xr, yref=yr, showarrow=True, ax=5.0, ay=1.2,
                               axref=xr, ayref=yr, arrowhead=0, arrowwidth=0.8, arrowcolor=C[3], yanchor="bottom",
                               text="across-signal<br>rotates the direction", font=dict(size=10, color=C[3]))
            fig.add_annotation(x=(v[0] + foot[0]) / 2, y=-0.12, xref=xr, yref=yr, showarrow=True, ax=5.0, ay=-1.0,
                               axref=xr, ayref=yr, arrowhead=0, arrowwidth=0.8, arrowcolor=C[2], yanchor="top",
                               text="along-signal<br>changes length only", font=dict(size=10, color=C[2]))
            fig.add_annotation(x=v[0] - 1.0, y=-sd - 0.12, xref=xr, yref=yr, showarrow=False, yanchor="top",
                               text="RMS scatter ∝ 1/√<i>n</i>", font=dict(size=10))
        else:
            fig.add_annotation(x=v[0], y=sd + 0.3, xref=xr, yref=yr, showarrow=False, yanchor="bottom",
                               text="RMS scatter ∝ 1/√<i>n</i>", font=dict(size=10))
            fig.add_annotation(x=v[0] + 0.1, y=-0.55, xref=xr, yref=yr, showarrow=False, yanchor="top",
                               text="<b><i>v̂</i></b> locks onto <b><i>v</i></b>,  <i>θ</i> → 0")
        fig.add_trace(go.Scatter(x=[0], y=[0], mode="markers", marker=dict(size=6, color=INK), showlegend=False),
                      row=1, col=col)
    for col in (1, 2):
        fig.update_xaxes(range=[-0.3, 6.2], visible=False, row=1, col=col)
        fig.update_yaxes(range=[-1.75, 1.6], visible=False, scaleanchor=f"x{'' if col == 1 else col}",
                         scaleratio=1, row=1, col=col)
    fig.update_annotations(selector=dict(text="(a)  few cells per arm,  <i>n</i> = 12"), x=0, xanchor="left", font_size=12)
    fig.update_annotations(selector=dict(text="(b)  many cells per arm,  <i>n</i> = 200"), x=0.52, xanchor="left", font_size=12)
    fig.update_layout(margin=dict(l=8, r=8, t=30, b=8), showlegend=False)
    save(fig, "fig2_geometry", 7.4, 2.7)


# ------------------------------------------------------------------------------------------------
# Figure 3: simulation study with known truth (scripts/simulation_study.py)
# ------------------------------------------------------------------------------------------------
def fig3():
    S = load("simulation_summary.json")
    kb = S["meta"]["k_base"]
    fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.11, vertical_spacing=0.2,
                        subplot_titles=("(a)  Two-arm quota: magnitude and control pool",
                                        "(b)  Noise law (quota from the true covariance)",
                                        f"(c)  Plug-in quota from a pilot, <i>m</i> = {kb:g} <i>m</i><sub>min</sub>",
                                        "(d)  Baseline rules, <i>n</i><sub>c</sub> = 3,100"))
    ref = dict(color=INK, width=1, dash="dash")
    ebar = lambda se: dict(type="data", array=se, visible=True, thickness=1, width=2)
    # (a) ratio vs k for each control-pool size
    for i, (nc, rows) in enumerate(S["block_a"]["panel"].items()):
        rows = [r for r in rows if r["ratio"] is not None]
        lab = "∞" if nc == "inf" else f"{float(nc):,.0f}"
        fig.add_trace(go.Scatter(x=[r["k_nominal"] for r in rows], y=[r["ratio"] for r in rows],
                                 error_y=ebar([r["ratio_se"] for r in rows]), mode="lines+markers",
                                 name=f"<i>n</i><sub>c</sub> = {lab}", legend="legend",
                                 marker=dict(size=6, color=C[i], line=dict(width=0)), line=dict(color=C[i], width=1.4),
                                 hoverinfo="skip"), row=1, col=1)
    fig.add_hline(y=1, line=ref, row=1, col=1)
    # (b) ratio by noise law, one marker series per magnitude
    laws = ["gaussian", "t5", "t3", "mixture (full covariance)", "mixture (within-component covariance)",
            "negative binomial counts"]
    lab_law = {"gaussian": "Gaussian", "t5": "<i>t</i><sub>5</sub>", "t3": "<i>t</i><sub>3</sub>",
               "mixture (full covariance)": "mixture,<br>full Σ", "mixture (within-component covariance)": "mixture,<br>within Σ",
               "negative binomial counts": "NB counts,<br>gene space"}
    B = [r for r in S["block_b"]["rows"] if r["ratio"] is not None]
    ks = sorted({r["k_nominal"] for r in B})
    for j, k in enumerate(ks):
        rows = [r for r in B if r["k_nominal"] == k]
        fig.add_trace(go.Scatter(x=[lab_law[r["law"]] for r in rows], y=[r["ratio"] for r in rows],
                                 error_y=ebar([r["ratio_se"] for r in rows]), mode="markers",
                                 name=f"<i>m</i> = {k:g} <i>m</i><sub>min</sub>", legend="legend2",
                                 marker=dict(size=7, color=C[j], line=dict(width=0), symbol=["circle", "diamond", "square", "triangle-up"][j % 4]),
                                 hoverinfo="skip"), row=1, col=2)
    fig.add_hline(y=1, line=ref, row=1, col=2)
    fig.update_xaxes(categoryorder="array", categoryarray=[lab_law[l] for l in laws], tickfont=dict(size=9), row=1, col=2)
    # (c) plug-in quota / true quota from pilots of 50, 200, 1,000 cells per arm (5-25-50-75-95 percentiles)
    Cr = [r for r in S["block_c"]["rows"] if r["k_nominal"] == kb and r["n_hat_over_true_q"]]
    for j, est in enumerate(["sample", "ledoit-wolf"]):
        rows = sorted([r for r in Cr if r["estimator"] == est], key=lambda r: r["n_pilot"])
        q = np.array([r["n_hat_over_true_q"] for r in rows])
        fig.add_trace(go.Box(x=[f"{r['n_pilot']:,}" for r in rows], lowerfence=q[:, 0], q1=q[:, 1], median=q[:, 2],
                             q3=q[:, 3], upperfence=q[:, 4], name={"sample": "sample covariance", "ledoit-wolf": "Ledoit–Wolf"}[est],
                             legend="legend3", marker=dict(color=C[j]), line=dict(color=C[j], width=1.2),
                             fillcolor=rgba(C[j], 0.35), hoverinfo="skip", offsetgroup=est), row=2, col=1)
    fig.add_hline(y=1, line=ref, row=2, col=1)
    fig.update_layout(boxmode="group", boxgap=0.35, boxgroupgap=0.15)
    # (d) baselines at the base pool
    A3 = [r for r in S["block_a"]["panel"]["3100.0"] if r["ratio"] is not None]
    series = [("two-arm (this paper)", [(r["k_nominal"], r["ratio"], r["ratio_se"]) for r in A3])]
    for rule in ("equal-arm", "large-pool", "isotropic"):
        rows = [r for r in S["block_d"]["base_pool"] if r["rule"] == rule and r["ratio"] is not None]
        series.append((rule, [(r["k_nominal"], r["ratio"], r["ratio_se"]) for r in rows]))
    for j, (name, pts) in enumerate(series):
        pts.sort()
        fig.add_trace(go.Scatter(x=[p[0] for p in pts], y=[p[1] for p in pts], error_y=ebar([p[2] for p in pts]),
                                 mode="lines+markers", name=name, legend="legend4",
                                 marker=dict(size=6, color=C[j], line=dict(width=0),
                                             symbol=["circle", "diamond", "square", "triangle-up"][j]),
                                 line=dict(color=C[j], width=1.4, dash="solid" if j == 0 else "dot"),
                                 hoverinfo="skip"), row=2, col=2)
    fig.add_hline(y=1, line=ref, row=2, col=2)
    kt = [k for k in S["meta"]["k_grid"] if k != 1.25]
    for (r, c) in ((1, 1), (2, 2)):
        fig.update_xaxes(type="log", tickvals=kt, ticktext=[f"{k:g}" for k in kt], tickangle=0,
                         title_text="effect magnitude <i>m</i> / floor <i>m</i><sub>min</sub>", row=r, col=c)
    fig.update_xaxes(title_text="pilot cells per arm", row=2, col=1)
    fig.update_yaxes(title_text="realized RMS angle / tolerance", range=[0.905, 1.06], row=1, col=1)
    fig.update_yaxes(title_text="realized RMS angle / tolerance", range=[0.9, 1.15], row=1, col=2)
    fig.update_yaxes(title_text="plug-in quota / true quota", type="log", row=2, col=1)
    fig.update_yaxes(title_text="realized RMS angle / tolerance", range=[0.6, 1.45], row=2, col=2)
    for a in fig.layout.annotations[:4]:
        a.update(xanchor="left", font=dict(size=11))
    fig.layout.annotations[0].x = 0.0; fig.layout.annotations[1].x = 0.555
    fig.layout.annotations[2].x = 0.0; fig.layout.annotations[3].x = 0.555
    lg = dict(font=dict(size=9), bgcolor="rgba(255,255,255,0.6)", itemsizing="constant", tracegroupgap=0)
    fig.update_layout(legend=dict(lg, x=0.44, y=0.615, xanchor="right", yanchor="bottom", font=dict(size=8.5)),
                      legend2=dict(lg, x=0.555, y=0.99, xanchor="left", yanchor="top"),
                      legend3=dict(lg, x=0.44, y=0.40, xanchor="right", yanchor="top"),
                      legend4=dict(lg, x=0.99, y=0.40, xanchor="right", yanchor="top"),
                      margin=dict(l=60, r=10, t=30, b=50))
    save(fig, "fig3_simulation_study", 7.4, 5.8)


# ------------------------------------------------------------------------------------------------
# Figure 4: Tahoe-100M sufficiency spectrum
# ------------------------------------------------------------------------------------------------
def fig4():
    H = load("fig3_hist.json")
    pct = H["pct"]
    vals = [pct["OVER"], pct["UNDER"], pct["GHOST"], pct["POOL-LIMITED"], pct["NOT-DETECTABLE"]]
    fig = make_subplots(rows=2, cols=1, row_heights=[0.16, 0.84], vertical_spacing=0.17,
                        subplot_titles=(f"(a)  All {H['n_all']:,} conditions by regime",
                                        f"(b)  Finite quotas of the {H['pct_detectable']:.1f}% detectable conditions"))
    left = 0.0
    for r, p in zip(REGIMES, vals):
        fig.add_trace(go.Bar(x=[p], y=[""], orientation="h", width=0.6, marker=dict(color=RCOL[r], line=dict(width=0)),
                             name=RKEY[r], legendgroup=r, showlegend=True, hoverinfo="skip"), row=1, col=1)
        mid = left + p / 2
        if p >= 5:
            fig.add_annotation(x=mid, y=0, text=f"<b>{p:.1f}%</b>", showarrow=False, font=dict(color="white", size=11),
                               xref="x", yref="y")
        else:
            fig.add_annotation(x=mid, y=0.32, ax=mid, ay=0.55, xref="x", yref="y", axref="x", ayref="y",
                               text=f"{p:.1f}%", showarrow=True, arrowhead=0, arrowwidth=0.7, arrowcolor="#8a96a8",
                               yanchor="bottom", font=dict(size=10))
        left += p
    fig.update_xaxes(range=[0, 100], visible=False, row=1, col=1)
    fig.update_yaxes(visible=False, range=[-0.32, 0.85], row=1, col=1)

    edges = np.log10(np.array(H["bin_edges"]))
    mids, widths = (edges[:-1] + edges[1:]) / 2, np.diff(edges)
    keep = edges[:-1] < math.log10(2e6)
    for k, r in (("OVER", "over-sampled"), ("UNDER", "treated-depth-limited"), ("GHOST", "ghost")):
        fig.add_trace(go.Bar(x=mids[keep], y=np.array(H["counts"][k])[keep], width=widths[keep] * 0.96,
                             marker=dict(color=RCOL[r], line=dict(width=0)), legendgroup=r, showlegend=False,
                             hoverinfo="skip"), row=2, col=1)
    tot = np.sum([np.array(H["counts"][k])[keep] for k in ("OVER", "UNDER", "GHOST")], axis=0)
    ytop = tot.max()
    refs = [(H["N0"], "dash", f"median depth <i>N</i><sub>0</sub> = {H['N0']:,}", "right", 1.06),
            (H["median_n_star"], "solid", f"median <i>n</i>* = {H['median_n_star']:,}", "left", 1.17),
            (H["ghost"], "dot", f"ghost <i>n</i>* ≥ {H['ghost']:,}", "left", 1.06)]
    for x, dash, lab, anchor, yf in refs:
        lx = math.log10(x)
        fig.add_shape(type="line", x0=lx, x1=lx, y0=0, y1=ytop * 1.12, xref="x2", yref="y2",
                      line=dict(color=INK, width=1.3, dash=dash))
        fig.add_annotation(x=lx, y=ytop * yf, xref="x2", yref="y2", text=lab, showarrow=False,
                           xanchor="right" if anchor == "right" else "left", xshift=-3 if anchor == "right" else 3,
                           font=dict(size=10))
    fig.update_layout(barmode="stack", bargap=0)
    fig.update_xaxes(range=[math.log10(20), math.log10(2e6)], tickvals=list(range(2, 7)),
                     ticktext=[sup(k) for k in range(2, 7)], showgrid=False, showline=True, linecolor="#c8d0dc",
                     ticks="outside", tickcolor="#c8d0dc",
                     title_text="required treated cells per condition, two-arm quota <i>n</i>*", row=2, col=1)
    fig.update_yaxes(range=[0, ytop * 1.25], title_text="number of conditions", separatethousands=True,
                     tickformat=",", row=2, col=1)
    for a in fig.layout.annotations[:2]:
        a.update(x=0, xanchor="left", font=dict(size=12))
    fig.update_layout(legend=dict(x=1.02, y=1.0, xanchor="left", yanchor="top", tracegroupgap=6,
                                  font=dict(size=10), itemsizing="constant"),
                      margin=dict(l=60, r=10, t=30, b=50))
    save(fig, "fig4_tahoe_spectrum", 7.4, 4.5)


# ------------------------------------------------------------------------------------------------
# Figure 5: cross-modality summary
# ------------------------------------------------------------------------------------------------
def fig5():
    K = load("fig4_kde.json")
    S = {s["name"]: s for s in load("unified_spectrum.json")["headline"]}
    grid = np.array(K["grid_log10m"])
    modlab = {"chemical": "chemical, small-molecule (Vevo Mosaic)",
              "genome-wide": "genome-wide CRISPRi (X-Atlas/Orion)",
              "essential-gene": "essential-gene CRISPRi (TRADE)"}
    mcol = {"chemical": C[5], "genome-wide": C[6], "essential-gene": C[9]}
    fig = make_subplots(rows=1, cols=2, column_widths=[0.47, 0.53], horizontal_spacing=0.09,
                        subplot_titles=("(a)  Effect magnitude per screen", "(b)  Sufficiency regime per screen"))
    order = list(reversed(SCREENS))
    shown = set()
    sel = (grid >= math.log10(0.03)) & (grid <= math.log10(30))
    for sc in K["screens"]:
        i = order.index(sc["name"])
        d = np.array(sc["density"])
        h = d / d.max() * 0.85
        gx, gh = grid[sel], h[sel]
        col = mcol[sc["modality"]]
        fig.add_trace(go.Scatter(x=np.r_[gx, gx[::-1]], y=np.r_[i + gh, np.full(gx.size, i)], fill="toself",
                                 fillcolor=rgba(col, 0.38), line=dict(width=0), hoverinfo="skip",
                                 name=modlab[sc["modality"]], legendgroup=sc["modality"], legend="legend2",
                                 showlegend=sc["modality"] not in shown), row=1, col=1)
        fig.add_trace(go.Scatter(x=gx, y=i + gh, mode="lines", line=dict(color=col, width=1.6), showlegend=False,
                                 hoverinfo="skip"), row=1, col=1)
        shown.add(sc["modality"])
        lm = math.log10(sc["median_m"])
        fig.add_shape(type="line", x0=lm, x1=lm, y0=i, y1=i + 0.85, xref="x", yref="y", line=dict(color=INK, width=1.4))
        fig.add_annotation(x=lm, y=i + 0.85, xref="x", yref="y", text=f"{sc['median_m']:.2f}", showarrow=False,
                           xanchor="left", yanchor="top", xshift=3, font=dict(size=9.5))
    fig.update_xaxes(range=[math.log10(0.03), math.log10(30)], tickvals=[-1, 0, 1], ticktext=["0.1", "1", "10"],
                     title_text="bias-corrected effect magnitude <i>m</i> (log scale)", showline=True,
                     linecolor="#c8d0dc", ticks="outside", tickcolor="#c8d0dc", row=1, col=1)
    fig.update_yaxes(tickvals=list(range(len(order))), ticktext=order, range=[-0.15, len(order) - 0.1],
                     showgrid=False, row=1, col=1)

    keys = ["pct_over", "pct_under", "pct_ghost", "pct_pool_limited", "not_detectable_pct"]
    xl = ["Tahoe-<br>100M", "Emerald-<br>Bay", "Orion<br>HCT116", "Orion<br>HEK293T", "TRADE<br>Jurkat", "TRADE<br>HepG2"]
    for r, k in reversed(list(zip(REGIMES, keys))):
        ys = [S[nm][k] for nm in SCREENS]
        txt = [f"{y:.0f}%" if (y >= 7 and r != "over-sampled") else "" for y in ys]
        fig.add_trace(go.Bar(x=xl, y=ys, marker=dict(color=RCOL[r], line=dict(width=0)), name=RKEY[r],
                             legendgroup=r, text=txt, textposition="inside", insidetextanchor="middle",
                             textfont=dict(color="white", size=10), hoverinfo="skip", width=0.62), row=1, col=2)
    for x, nm in zip(xl, SCREENS):
        fig.add_annotation(x=x, y=101, xref="x2", yref="y2", text=f"<b>{S[nm]['pct_over']:.1f}%</b>",
                           showarrow=False, yanchor="bottom", font=dict(size=10, color=RCOL["over-sampled"]))
    fig.add_annotation(x=-0.5, y=101, xref="x2", yref="y2", text="over-<br>sampled", showarrow=False,
                       xanchor="right", yanchor="bottom", font=dict(size=9, color=RCOL["over-sampled"]))
    for x0, x1, lab in ((-0.35, 1.35, "chemical"), (1.65, 5.35, "genetic (CRISPRi)")):
        fig.add_shape(type="line", x0=x0, x1=x1, y0=-0.2, y1=-0.2, xref="x2", yref="paper",
                      line=dict(color="#8a96a8", width=1))
        fig.add_annotation(x=(x0 + x1) / 2, y=-0.22, xref="x2", yref="paper", text=lab, showarrow=False,
                           yanchor="top", font=dict(size=10, color="#5b6b82"))
    fig.update_xaxes(range=[-1.1, 5.5], tickfont=dict(size=10), row=1, col=2)
    fig.update_yaxes(range=[0, 112], tickvals=[0, 25, 50, 75, 100], title_text="% of conditions", row=1, col=2)
    for a in fig.layout.annotations[:2]:
        a.update(x=0 if a.text.startswith("(a)") else 0.555, xanchor="left", font=dict(size=12))
    fig.update_layout(barmode="stack", bargap=0.3,
                      legend=dict(x=1.01, y=1.0, xanchor="left", yanchor="top", tracegroupgap=6, font=dict(size=10),
                                  traceorder="reversed"),
                      legend2=dict(x=0, y=-0.2, xanchor="left", yanchor="top", font=dict(size=10),
                                   orientation="v", tracegroupgap=0),
                      margin=dict(l=95, r=10, t=30, b=95))
    save(fig, "fig5_crossmodality_summary", 8.6, 4.2)


# ------------------------------------------------------------------------------------------------
# Figure 6: held-out downsample-and-measure test, one panel per screen
# ------------------------------------------------------------------------------------------------
def fig6():
    src = [("Tahoe-100M (chemical)", "tahoe_direct_curves_full.json", "per_group"),
           ("EmeraldBay (chemical)", "emeraldbay_falsification_full.json", "per_group_slope"),
           ("Orion HCT116 (genetic)", "orion_HCT116_heldout.json", "per_group"),
           ("Orion HEK293T (genetic)", "orion_HEK293T_heldout.json", "per_group"),
           ("TRADE Jurkat (genetic)", "trade_jurkat_falsification.json", "rows"),
           ("TRADE HepG2 (genetic)", "trade_hepg2_falsification.json", "rows")]
    panels = []
    for name, f, key in src:
        rows = [r for r in load(f)[key] if r["N"] >= 400]
        rho2 = np.array([r["rho2"] for r in rows])
        ratio = np.array([r["ratio"] for r in rows])
        x = 1 / np.maximum(rho2, 1e-3)
        inr = rho2 >= 3
        sp = spearmanr(x, 1 - ratio).statistic
        if inr.any():
            mid = f"ρ² ≥ 3: median ratio {np.median(ratio[inr]):.2f} (n = {inr.sum():,})"
        else:
            mid = f"all ρ² < 3 (max {rho2.max():.1f})"
        title = f"<b>{name}</b>, {len(rows):,} groups<br>{mid}, Spearman {sp:.2f}"
        panels.append((title, x, ratio, inr))
    fig = make_subplots(rows=2, cols=3, subplot_titles=[p[0] for p in panels], horizontal_spacing=0.06,
                        vertical_spacing=0.2, shared_yaxes=True)
    cls = [("first-order regime (ρ² ≥ 3)", C[0], True), ("below it (ρ² < 3)", C[1], False)]
    for i, (_, x, ratio, inr) in enumerate(panels):
        r, c = i // 3 + 1, i % 3 + 1
        for lab, col, flag in cls:
            m = inr if flag else ~inr
            fig.add_trace(go.Scatter(x=x[m], y=ratio[m], mode="markers", name=lab, legendgroup=lab,
                                     showlegend=(i == 0), hoverinfo="skip",
                                     marker=dict(size=3.2, color=col, opacity=0.5, line=dict(width=0))), row=r, col=c)
        fig.add_hline(y=1, line=dict(color=INK, width=1, dash="dash"), row=r, col=c)
        fig.add_vline(x=1 / 3, line=dict(color="#8a96a8", width=1, dash="dot"), row=r, col=c)
        fig.update_xaxes(type="log", range=[-4.3, 3.3], tickvals=[1e-4, 1e-2, 1, 1e2], ticktext=["10⁻⁴", "10⁻²", "1", "10²"],
                         tickangle=0, row=r, col=c)
        fig.update_yaxes(range=[-0.05, 1.5], row=r, col=c)
    for a in fig.layout.annotations:
        if a.text.startswith("<b>"):
            a.update(font=dict(size=10), xanchor="left", x=a.x - 0.145, align="left")
    fig.update_xaxes(title_text="1 / ρ²  (signal-to-noise decreases to the right)", row=2, col=2)
    fig.update_yaxes(title_text="realized / predicted slope", row=1, col=1)
    fig.update_yaxes(title_text="realized / predicted slope", row=2, col=1)
    fig.update_layout(legend=dict(orientation="h", x=0.5, xanchor="center", y=-0.14, yanchor="top",
                                  itemsizing="constant", font=dict(size=11)),
                      margin=dict(l=60, r=10, t=50, b=70))
    save(fig, "fig6_heldout_per_dataset", 8.6, 5.6)


if __name__ == "__main__":
    print("Generating plotly figures ->", FIG)
    fig2(); fig3(); fig4(); fig5(); fig6()
    print("done.")
