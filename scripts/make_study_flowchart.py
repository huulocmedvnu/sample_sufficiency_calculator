#!/usr/bin/env python3
"""
Publication-ready study flowchart (Graphviz Python API, HTML-table node labels).

    question -> closed-form quota -> airtight 4-way proof -> shared embedding engine
    -> cross-modality validation (chemical | genetic columns) -> one threshold, two ledgers
    -> wet-lab / dry-lab downstream utility.

Numbers mirror the manuscript headline constants of record. The equivalent Mermaid.js source
(for web / GitHub rendering) lives at docs/study_flowchart.mmd.

Requires:  pip install graphviz   +   the Graphviz `dot` binary on PATH.
Run:       python scripts/make_study_flowchart.py  ->  figures/fig1_study_flowchart.{pdf,png}
"""
import os
import graphviz

# ---- modern academic palette: muted slates, corporate grays, soft green/blue ----
INK, MUTE = "#1f2933", "#5b6875"          # primary title / metadata text
PAL = {                                    # (fill, border)
    "question": ("#eef1f4", "#8a97a6"),
    "theory":   ("#e8eef5", "#5b7a9d"),
    "proof":    ("#eef0f2", "#98a3b0"),
    "engine":   ("#e7efe9", "#6f9a80"),
    "chem":     ("#e9f0f7", "#6f8caa"),
    "gen":      ("#e9f1ec", "#6f9a84"),
    "hub":      ("#e8eef5", "#5b7a9d"),
    "wet":      ("#e7eef3", "#5f8aa0"),
    "dry":      ("#eef0f2", "#8a97a6"),
}


def html(title, subtitle=None, rows=(), title_pt=12, meta_pt=9):
    """Build an HTML-like label: bold title, optional muted subtitle, muted metadata rows."""
    tr = [f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{title_pt}" COLOR="{INK}">'
          f'<B>{title}</B></FONT></TD></TR>']
    if subtitle:
        tr.append(f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{meta_pt + 0.5}" COLOR="{MUTE}">'
                  f'<B>{subtitle}</B></FONT></TD></TR>')
    for r in rows:
        tr.append(f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{meta_pt}" COLOR="{MUTE}">'
                  f'{r}</FONT></TD></TR>')
    body = "".join(tr)
    return (f'<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2">'
            f'{body}</TABLE>>')


def add(target, nid, kind, title, subtitle=None, rows=(), sharp=False):
    """Add a styled node. sharp=True -> data-table rectangle; else rounded step."""
    fill, border = PAL[kind]
    target.node(nid, label=html(title, subtitle, rows),
                shape="box", style="filled" if sharp else "rounded,filled",
                fillcolor=fill, color=border)


def build():
    g = graphviz.Digraph("study")
    g.attr(rankdir="TB", splines="true", compound="true", newrank="true",
           bgcolor="white", nodesep="0.4", ranksep="0.55", pad="0.3")
    g.attr("node", fontname="Helvetica", penwidth="1.1")
    g.attr("edge", fontname="Helvetica", fontsize="9", fontcolor="#4a545e",
           color="#69747f", arrowsize="0.7", penwidth="1.1")

    # ---- top anchor + theory ------------------------------------------------
    add(g, "Q", "question", "Core research question",
        rows=["How many cells per arm resolve a perturbation's <I>direction</I>",
              "(mechanism) to an angular tolerance &#952;?"])
    add(g, "TH", "theory", "Closed-form cell quota",
        rows=["n* = 2&#183;tr(P&#931;P) / (m&#178; &#952;&#178;),   P = I &#8722; uu&#7488;",
              "only noise &#8869; to the effect rotates the direction",
              "calibrated &#8594; n* &#8776; 9,376 / m&#178;   (&#952; = 0.1 rad)"])

    # ---- airtight proof foundation (4 aligned boxes) ------------------------
    with g.subgraph(name="cluster_proof") as c:
        c.attr(label='<<B>Airtight foundation</B>  &#183;  proven four independent ways>',
               labeljust="l", style="rounded,filled", color="#c2cad3",
               fillcolor="#f4f6f8", fontname="Helvetica", fontsize="10", margin="12")
        c.attr(rank="same")
        add(c, "P1", "proof", "Independent impl.", rows=["cross-check"])
        add(c, "P2", "proof", "Monte-Carlo", rows=["&lt; 0.13% error"])
        add(c, "P3", "proof", "SymPy", rows=["symbolic Jacobian"])
        add(c, "P4", "proof", "Lean 4 / Mathlib", rows=["deterministic core, no <I>sorry</I>"])
        for a, b in (("P1", "P2"), ("P2", "P3"), ("P3", "P4")):  # lock left-to-right order
            c.edge(a, b, style="invis")

    # ---- shared embedding engine (horizontal infrastructure bridge) ---------
    add(g, "ENG", "engine", "Shared embedding engine",
        rows=["Normalize &#8594; log1p &#8594; HVG 2000 &#8594; PCA 50",
              "within-condition &#931; &#8594; &#963;&#178;   &#183;   centroids &#8594; magnitude m"])

    # ---- cross-modality validation: dashed container, two strict columns ----
    with g.subgraph(name="cluster_emp") as emp:
        emp.attr(label='<<B>Cross-Modality Empirical Validation</B>>', labeljust="l",
                 style="rounded,dashed,filled", color="#b3bcc6", fillcolor="#fbfcfd",
                 fontname="Helvetica", fontsize="11", margin="16")
        with emp.subgraph(name="cluster_chem") as ch:
            ch.attr(label='<<B>Chemical modality</B>  &#183;  Vevo Mosaic>', labeljust="l",
                    style="rounded,filled", color="#9cb2c6", fillcolor="#eef3f8",
                    fontname="Helvetica", fontsize="10", margin="12")
            add(ch, "TA", "chem", "Tahoe-100M", "Calibration anchor",
                ["100.6M cells &#183; 56,827 conditions",
                 "&#963;&#178; = 0.96 &#8594; n* = 9,376 / m&#178;",
                 "10.6% over &#183; 89.0% under &#183; 0.4% ghost"], sharp=True)
            add(ch, "EB", "chem", "EmeraldBay", "Out-of-distribution validation",
                ["1.83M cells &#183; 52 lines",
                 "&#963;&#178; &#8776; 0.94 &#183; held-out slope R&#178; &#8776; 0.999",
                 "gating 347 / 347 met"], sharp=True)
            ch.edge("TA", "EB", label="  &#963;&#178; transfers &#8776; 0.94&#8211;0.96")
        with emp.subgraph(name="cluster_gen") as ge:
            ge.attr(label='<<B>Genetic modality</B>  &#183;  CRISPRi Perturb-seq>', labeljust="l",
                    style="rounded,filled", color="#8fb0a0", fillcolor="#eef4f1",
                    fontname="Helvetica", fontsize="10", margin="12")
            add(ge, "OR", "gen", "X-Atlas/Orion", "Cross-modality transfer",
                ["genome-wide CRISPRi &#183; ~8M cells",
                 "18,903 knockdowns &#215; 2 lines",
                 "&#963;&#178; = 0.91 / 1.03 &#183; median m 0.14&#8211;0.35",
                 "0% over-sampled (magnitude-limited)"], sharp=True)
            add(ge, "TR", "gen", "TRADE", "Validated genetic falsification",
                ["essential-gene CRISPRi &#183; Jurkat + HepG2",
                 "2,393 genes/line &#183; slope 1.02 / 0.96 &#183; R&#178; = 0.996",
                 "median 45&#8211;85 cells/gene &#8594; depth-limited"], sharp=True)
            ge.edge("OR", "TR", label="  direct falsification")

    # ---- transition node + downstream utility (two aligned outcomes) --------
    add(g, "HUB", "hub", "One threshold n*, two ledgers")
    with g.subgraph(name="cluster_out") as out:
        out.attr(label='<<B>Downstream utility</B>>', labeljust="l",
                 style="rounded,filled", color="#c2cad3", fillcolor="#f4f6f8",
                 fontname="Helvetica", fontsize="11", margin="14")
        out.attr(rank="same")
        add(out, "WET", "wet", "Wet-lab",
            rows=["magnitude-adaptive budgeting", "~3.1&#215; fewer cells &#183; triage / multiplex"])
        add(out, "DRY", "dry", "Dry-lab",
            rows=["provably-safe downsampling", "triage of low-SNR similarity edges"])
        out.edge("WET", "DRY", style="invis")  # lock left-to-right order

    # ---- horizontal alignment of the two parallel columns -------------------
    for a, b in (("TA", "OR"), ("EB", "TR")):
        with g.subgraph() as s:
            s.attr(rank="same")
            s.node(a)
            s.node(b)

    # ---- clean, cluster-aware edges ----------------------------------------
    g.edge("Q", "TH")
    g.edge("TH", "P3", label="  proven", lhead="cluster_proof")
    g.edge("P3", "ENG", ltail="cluster_proof")
    g.edge("ENG", "TA", label="  calibrate", lhead="cluster_chem")
    g.edge("ENG", "OR", lhead="cluster_gen")
    g.edge("EB", "HUB", ltail="cluster_chem")
    g.edge("TR", "HUB", ltail="cluster_gen")
    g.edge("HUB", "WET", lhead="cluster_out")
    g.edge("HUB", "DRY", lhead="cluster_out")
    return g


if __name__ == "__main__":
    os.makedirs("figures", exist_ok=True)
    g = build()
    stem = "figures/fig1_study_flowchart"
    with open(stem + ".dot", "w") as f:
        f.write(g.source)
    for fmt in ("pdf", "png"):
        g.render(stem, format=fmt, cleanup=True, engine="dot")
        print(f"wrote {stem}.{fmt}")
