#!/usr/bin/env python3
"""
Publication-ready study flowchart (Graphviz Python API, HTML-table node labels).

    question -> closed-form quota -> airtight 4-way proof -> shared embedding engine
    -> cross-modality validation (chemical | genetic columns) -> one threshold, two ledgers
    -> wet-lab / dry-lab downstream utility.

High-contrast profile: white node fills, thick semantic-colored borders, dark text.
Numbers mirror the manuscript headline constants of record. The equivalent Mermaid.js source
(for web / GitHub rendering) lives at docs/study_flowchart.mmd.

Requires:  pip install graphviz   +   the Graphviz `dot` binary on PATH.
Run:       python scripts/make_study_flowchart.py  ->  figures/fig1_study_flowchart.{pdf,png}
"""
import os
import graphviz

# ---- high-contrast text + border palette -------------------------------------
TITLE = "#0F172A"     # titles/headers: pure dark charcoal (bold)
META  = "#334155"     # secondary text / statistical metrics: solid dark slate
FILL  = "#FFFFFF"     # all node fills: clean white so text pops
SLATE, BLUE, GREEN = "#334155", "#2563EB", "#16A34A"
BORDER = {            # node border colour by role
    "question": SLATE, "theory": SLATE, "proof": SLATE,
    "engine":   GREEN,                       # infrastructure
    "chem":     BLUE,                         # chemical modality
    "gen":      GREEN,                        # genetic modality
    "hub":      SLATE, "wet": SLATE, "dry": SLATE,
}


def html(title, subtitle=None, rows=(), title_pt=13, meta_pt=10):
    """HTML-like label: bold dark title, optional bold subtitle, dark metadata rows."""
    tr = [f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{title_pt}" COLOR="{TITLE}">'
          f'<B>{title}</B></FONT></TD></TR>']
    if subtitle:
        tr.append(f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{meta_pt}" COLOR="{META}">'
                  f'<B>{subtitle}</B></FONT></TD></TR>')
    for r in rows:
        tr.append(f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{meta_pt}" COLOR="{META}">'
                  f'{r}</FONT></TD></TR>')
    body = "".join(tr)
    return (f'<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2">'
            f'{body}</TABLE>>')


def add(target, nid, kind, title, subtitle=None, rows=(), sharp=False):
    """White-filled node with a thick coloured border. sharp=True -> data-table rectangle."""
    target.node(nid, label=html(title, subtitle, rows),
                shape="box", style="filled" if sharp else "rounded,filled",
                fillcolor=FILL, color=BORDER[kind], penwidth="2")


def build():
    g = graphviz.Digraph("study")
    g.attr(rankdir="TB", splines="true", compound="true", newrank="true",
           bgcolor="#FFFFFF", nodesep="0.42", ranksep="0.58", pad="0.3")
    g.attr("node", fontname="Helvetica")
    g.attr("edge", fontname="Helvetica", fontsize="10", fontcolor=META,
           color="#64748B", arrowsize="0.8", penwidth="1.5")

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
               labeljust="l", style="rounded", color="#94A3B8", penwidth="1.4",
               fontname="Helvetica", fontsize="11", fontcolor=TITLE, margin="12")
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

    # ---- cross-modality validation: transparent dashed container ------------
    with g.subgraph(name="cluster_emp") as emp:
        emp.attr(label='<<B>Cross-Modality Empirical Validation</B>>', labeljust="l",
                 style="rounded,dashed", color="#94A3B8", penwidth="1.5",
                 fontname="Helvetica", fontsize="11", fontcolor=TITLE, margin="16")
        with emp.subgraph(name="cluster_chem") as ch:
            ch.attr(label='<<B>Chemical modality</B>  &#183;  Vevo Mosaic>', labeljust="l",
                    style="rounded,filled", color=BLUE, fillcolor="#EFF6FF", penwidth="1.5",
                    fontname="Helvetica", fontsize="11", fontcolor=TITLE, margin="12")
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
                    style="rounded,filled", color=GREEN, fillcolor="#F0FDF4", penwidth="1.5",
                    fontname="Helvetica", fontsize="11", fontcolor=TITLE, margin="12")
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
                 style="rounded", color="#94A3B8", penwidth="1.4",
                 fontname="Helvetica", fontsize="11", fontcolor=TITLE, margin="14")
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
