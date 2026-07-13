#!/usr/bin/env python3
"""
Publication-ready study flowchart (Graphviz Python API, HTML-table node labels).

Minimalist, high-contrast overview that ends at the experimental validations:

    question -> closed-form quota -> airtight 4-way proof -> shared embedding engine
    -> two clean parallel columns (chemical: Tahoe-100M -> EmeraldBay ; genetic: Orion -> TRADE)

White node fills, thin semantic-coloured borders, dark text, borderless (label-only) containers.
Numbers mirror the manuscript headline constants of record.

Requires:  pip install graphviz   +   the Graphviz `dot` binary on PATH.
Run:       python scripts/make_study_flowchart.py  ->  figures/fig1_study_flowchart.{pdf,png}
"""
import os
import graphviz

# ---- high-contrast text + 3-colour border palette ----------------------------
TITLE = "#0F172A"     # main titles: pure dark charcoal (bold)
META  = "#475569"     # key statistical parameters: dark slate gray
FILL  = "#FFFFFF"     # all node fills: crisp white
SLATE, BLUE, GREEN = "#334155", "#1E40AF", "#065F46"
BORDER = {            # node outline by role
    "question": SLATE, "theory": SLATE, "proof": SLATE,   # theory / foundation
    "engine":   GREEN,                                    # infrastructure
    "chem":     BLUE,                                     # chemical modality
    "gen":      GREEN,                                    # genetic modality
}


def html(title, rows=(), title_pt=13, meta_pt=10):
    """HTML-like label: bold dark title over dark-slate metric rows."""
    tr = [f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{title_pt}" COLOR="{TITLE}">'
          f'<B>{title}</B></FONT></TD></TR>']
    for r in rows:
        tr.append(f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{meta_pt}" COLOR="{META}">'
                  f'{r}</FONT></TD></TR>')
    return (f'<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="2">'
            f'{"".join(tr)}</TABLE>>')


def add(target, nid, kind, title, rows=(), sharp=False):
    """White node with a thin coloured border. sharp=True -> data rectangle, else rounded step."""
    target.node(nid, label=html(title, rows), shape="box",
                style="filled" if sharp else "rounded,filled",
                fillcolor=FILL, color=BORDER[kind], penwidth="1.5")


def header(cluster, text, color="#CBD5E1"):
    """Clean minimalist cluster: a thin rounded border, no fill, bold header label."""
    cluster.attr(label=f'<<B>{text}</B>>', labeljust="l",
                 style="rounded", color=color, penwidth="1.5",
                 fontname="Helvetica", fontsize="11", fontcolor=TITLE, margin="14")


def build():
    g = graphviz.Digraph("study")
    g.attr(rankdir="TB", splines="true", compound="true", newrank="true",
           bgcolor="#FFFFFF", nodesep="0.45", ranksep="0.6", pad="0.3")
    g.attr("node", fontname="Helvetica")
    g.attr("edge", fontname="Helvetica", fontsize="10", fontcolor=SLATE,
           color="#64748B", arrowsize="0.8", penwidth="1.5")

    # ---- linear spine: question -> quota ------------------------------------
    add(g, "Q", "question", "Core research question",
        ["How many cells resolve a perturbation's <I>direction</I> to tolerance &#952;?"])
    add(g, "TH", "theory", "Closed-form cell quota",
        ["n* = 2&#183;tr(P&#931;P) / (m&#178; &#952;&#178;)",
         "calibrated &#8594; n* &#8776; 9,376 / m&#178;"])

    # ---- airtight proof foundation (borderless header + aligned name chips) --
    with g.subgraph(name="cluster_proof") as c:
        header(c, "Airtight foundation  &#183;  proven four independent ways")
        c.attr(rank="same")
        add(c, "P1", "proof", "Independent impl.")
        add(c, "P2", "proof", "Monte-Carlo")
        add(c, "P3", "proof", "SymPy")
        add(c, "P4", "proof", "Lean 4 / Mathlib")
        for a, b in (("P1", "P2"), ("P2", "P3"), ("P3", "P4")):  # lock left-to-right order
            c.edge(a, b, style="invis")

    # ---- shared embedding engine (infrastructure) ---------------------------
    add(g, "ENG", "engine", "Shared embedding engine", ["Estimates within-condition &#963;&#178;"])

    # ---- two clean parallel columns in a borderless space -------------------
    with g.subgraph(name="cluster_emp") as emp:
        header(emp, "Cross-Modality Empirical Validation")
        with emp.subgraph(name="cluster_chem") as ch:
            header(ch, "Chemical modality", BLUE)
            add(ch, "TA", "chem", "Tahoe-100M",
                ["&#963;&#178; = 0.96", "in-regime slope 0.94 &#183; R&#178; &#8776; 0.999"], sharp=True)
            add(ch, "EB", "chem", "EmeraldBay",
                ["&#963;&#178; &#8776; 0.94", "held-out slope 0.98 &#183; R&#178; &#8776; 0.999"], sharp=True)
            ch.edge("TA", "EB", label="  &#963;&#178; transfers &#8776; 0.94&#8211;0.96")
        with emp.subgraph(name="cluster_gen") as ge:
            header(ge, "Genetic modality", GREEN)
            add(ge, "OR", "gen", "X-Atlas/Orion",
                ["&#963;&#178; = 0.91 / 1.03", "n* &#8776; 9,000 / m&#178; &#183; &#963;&#178; transfers"], sharp=True)
            add(ge, "TR", "gen", "TRADE",
                ["&#963;&#178; = 1.5 / 1.9", "slope 1.02 / 0.96 &#183; R&#178; = 0.996"], sharp=True)
            ge.edge("OR", "TR", label="  direct falsification")

    # ---- strict horizontal alignment of the two columns ---------------------
    for a, b in (("TA", "OR"), ("EB", "TR")):
        with g.subgraph() as s:
            s.attr(rank="same")
            s.node(a)
            s.node(b)

    # ---- clean, cluster-aware edges (spine ends at the validations) ---------
    g.edge("Q", "TH")
    g.edge("TH", "P3", label="  proven", lhead="cluster_proof")
    g.edge("P3", "ENG", ltail="cluster_proof")
    g.edge("ENG", "TA", label="  calibrate", lhead="cluster_chem")
    g.edge("ENG", "OR", label="  transfer", lhead="cluster_gen")
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
