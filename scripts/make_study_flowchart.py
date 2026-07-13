#!/usr/bin/env python3
"""
Publication-ready study flowchart (Graphviz Python API, HTML-table node labels).

Minimalist, high-contrast overview that ends at the experimental validations:

    question -> closed-form quota -> airtight 4-way proof -> shared embedding engine
    -> two clean parallel columns (chemical: Tahoe-100M -> EmeraldBay ; genetic: Orion -> TRADE)

Typography/colour tuned for print: Helvetica, italic scalar variables, white node fills,
thin semantic borders in the Okabe-Ito palette shared with the paper's other figures,
light-gray rounded grouping boxes. Numbers mirror the manuscript constants of record.

Requires:  pip install graphviz   +   the Graphviz `dot` binary on PATH.
Run:       python scripts/make_study_flowchart.py  ->  figures/fig1_study_flowchart.{pdf,png}
"""
import os
import graphviz

# ---- palette (Okabe-Ito accents, matching scripts/make_manuscript_figures.py) ----
INK, META, FILL = "#0F172A", "#475569", "#FFFFFF"   # title / metric / node fill
SLATE = "#334155"                                    # theory & foundation (neutral)
BLUE  = "#0072B2"                                    # chemical modality
GREEN = "#009E73"                                    # infrastructure + genetic modality
LGRAY = "#CBD5E1"                                    # grouping-box borders
BORDER = {"question": SLATE, "theory": SLATE, "proof": SLATE,
          "engine": GREEN, "chem": BLUE, "gen": GREEN}

# italic scalar variables (proper math typography); Greek/entities render in labels
N, M2, S2, R2, THETA = ("<I>n</I>*", "<I>m</I>&#178;", "<I>&#963;</I>&#178;",
                        "<I>R</I>&#178;", "<I>&#952;</I>")


def html(title, rows=(), title_pt=13, meta_pt=10):
    """HTML-like label: bold dark title, a hairline of space, then dark-slate metric rows."""
    cells = [f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{title_pt}" COLOR="{INK}">'
             f'<B>{title}</B></FONT></TD></TR>']
    if rows:
        cells.append('<TR><TD HEIGHT="5"></TD></TR>')          # title/metric separation
    for r in rows:
        cells.append(f'<TR><TD ALIGN="CENTER"><FONT POINT-SIZE="{meta_pt}" COLOR="{META}">'
                     f'{r}</FONT></TD></TR>')
    return (f'<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="3">'
            f'{"".join(cells)}</TABLE>>')


def add(target, nid, kind, title, rows=(), sharp=False):
    """White node with a thin coloured border. sharp=True -> data rectangle, else rounded step."""
    target.node(nid, label=html(title, rows), shape="box",
                style="filled" if sharp else "rounded,filled",
                fillcolor=FILL, color=BORDER[kind], penwidth="1.7")


def header(cluster, text, color=LGRAY, pen="1.3"):
    """Clean grouping box: thin rounded border, no fill, bold header label."""
    cluster.attr(label=f'<<B>{text}</B>>', labeljust="l", style="rounded",
                 color=color, penwidth=pen, fontname="Helvetica",
                 fontsize="10.5", fontcolor=INK, margin="14")


def build():
    g = graphviz.Digraph("study")
    g.attr(rankdir="TB", splines="true", compound="true", newrank="true",
           bgcolor="#FFFFFF", nodesep="0.5", ranksep="0.62", pad="0.35")
    g.attr("node", fontname="Helvetica", margin="0.20,0.13")
    g.attr("edge", fontname="Helvetica-Oblique", fontsize="9.5", fontcolor=META,
           color="#64748B", arrowsize="0.85", penwidth="1.4")

    # ---- linear spine: question -> quota ------------------------------------
    add(g, "Q", "question", "Core research question",
        [f"How many cells resolve a perturbation's <I>direction</I> to tolerance {THETA}?"])
    add(g, "TH", "theory", "Closed-form cell quota",
        [f"{N} = 2&#183;tr(P&#931;P) / ({M2} {THETA}&#178;)",
         f"calibrated &#8594; {N} &#8776; 9,376 / {M2}"])

    # ---- airtight proof foundation (grouping box + aligned name chips) ------
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
    add(g, "ENG", "engine", "Shared embedding engine", [f"Estimates within-condition {S2}"])

    # ---- two parallel columns inside a light-gray grouping box --------------
    with g.subgraph(name="cluster_emp") as emp:
        header(emp, "Cross-Modality Empirical Validation")
        with emp.subgraph(name="cluster_chem") as ch:
            header(ch, "Chemical modality", BLUE, pen="1.5")
            add(ch, "TA", "chem", "Tahoe-100M",
                [f"{S2} = 0.96", f"in-regime slope 0.94 &#183; {R2} &#8776; 0.999"], sharp=True)
            add(ch, "EB", "chem", "EmeraldBay",
                [f"{S2} &#8776; 0.94", f"held-out slope 0.98 &#183; {R2} &#8776; 0.999"], sharp=True)
            ch.edge("TA", "EB", label="  &#963;&#178; transfers &#8776; 0.94&#8211;0.96")
        with emp.subgraph(name="cluster_gen") as ge:
            header(ge, "Genetic modality", GREEN, pen="1.5")
            add(ge, "OR", "gen", "X-Atlas/Orion",
                [f"{S2} = 0.91 / 1.03", f"{N} &#8776; 9,000 / {M2} &#183; &#963;&#178; transfers"], sharp=True)
            add(ge, "TR", "gen", "TRADE",
                [f"{S2} = 1.5 / 1.9", f"slope 1.02 / 0.96 &#183; {R2} = 0.996"], sharp=True)
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
