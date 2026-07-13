#!/usr/bin/env python3
"""
Study flowchart (Graphviz Python API, HTML-table node labels) — classic monochrome style.

    question -> closed-form quota -> airtight 4-way proof -> shared embedding engine
    -> two parallel columns (chemical: Tahoe-100M -> EmeraldBay ; genetic: Orion -> TRADE)

Black-and-white, print-classic: white background, black text and borders, Times New Roman,
italic scalar variables. Numbers mirror the manuscript constants of record.

Requires:  pip install graphviz   +   the Graphviz `dot` binary on PATH.
Run:       python scripts/make_study_flowchart.py  ->  figures/fig1_study_flowchart.{pdf,png}
"""
import os
import graphviz

# Real Microsoft Times New Roman IS installed (msttcorefonts, ~/.local/share/fonts/), but this
# build box's graphviz/pango cannot embed it -- it substitutes DejaVu Sans for "Times New Roman"
# (and even for the system "Nimbus Roman"), a local pango font-matching defect. "Liberation Serif"
# is Red Hat's metric- and shape-identical open clone of Times New Roman and DOES embed here, so we
# render with it (visually indistinguishable). On a machine with working pango, set "Times New Roman".
FONT = "Liberation Serif"
BLACK, WHITE = "#000000", "#FFFFFF"

# italic scalar variables (classic math typography); Greek/entities render in labels
N, M2, S2, R2, THETA = ("<I>n</I>*", "<I>m</I>&#178;", "<I>&#963;</I>&#178;",
                        "<I>R</I>&#178;", "<I>&#952;</I>")


def html(title, rows=(), title_pt=14, meta_pt=11):
    """HTML-like label: bold title, a hairline of space, then metric rows. All black text."""
    cells = [f'<TR><TD ALIGN="CENTER"><FONT FACE="{FONT}" POINT-SIZE="{title_pt}" '
             f'COLOR="{BLACK}"><B>{title}</B></FONT></TD></TR>']
    if rows:
        cells.append('<TR><TD HEIGHT="5"></TD></TR>')
    for r in rows:
        cells.append(f'<TR><TD ALIGN="CENTER"><FONT FACE="{FONT}" POINT-SIZE="{meta_pt}" '
                     f'COLOR="{BLACK}">{r}</FONT></TD></TR>')
    return (f'<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="3">'
            f'{"".join(cells)}</TABLE>>')


def add(target, nid, title, rows=()):
    """Plain rectangular box: white fill, thin black border, black text."""
    target.node(nid, label=html(title, rows), shape="box", style="filled",
                fillcolor=WHITE, color=BLACK, penwidth="1")


def header(cluster, text):
    """Grouping box: thin black rectangular border, bold black header label."""
    cluster.attr(label=f'<<FONT FACE="{FONT}"><B>{text}</B></FONT>>', labeljust="l",
                 style="solid", color=BLACK, penwidth="1", fontname=FONT, fontsize="11",
                 fontcolor=BLACK, margin="14")


def build():
    g = graphviz.Digraph("study")
    g.attr(rankdir="TB", splines="true", compound="true", newrank="true",
           bgcolor=WHITE, nodesep="0.5", ranksep="0.62", pad="0.35")
    g.attr("node", fontname=FONT, margin="0.20,0.13")
    g.attr("edge", fontname=FONT, fontsize="10", fontcolor=BLACK,
           color=BLACK, arrowsize="0.85", penwidth="1")

    # ---- linear spine: question -> quota ------------------------------------
    add(g, "Q", "Core research question",
        [f"How many cells resolve a perturbation's <I>direction</I> to tolerance {THETA}?"])
    add(g, "TH", "Closed-form cell quota",
        [f"{N} = 2&#183;tr(P&#931;P) / ({M2} {THETA}&#178;)",
         f"calibrated &#8594; {N} &#8776; 9,376 / {M2}"])

    # ---- airtight proof foundation (grouping box + aligned name chips) ------
    with g.subgraph(name="cluster_proof") as c:
        header(c, "Airtight foundation  &#183;  proven four independent ways")
        c.attr(rank="same")
        add(c, "P1", "Independent impl.")
        add(c, "P2", "Monte-Carlo")
        add(c, "P3", "SymPy")
        add(c, "P4", "Lean 4 / Mathlib")
        for a, b in (("P1", "P2"), ("P2", "P3"), ("P3", "P4")):  # lock left-to-right order
            c.edge(a, b, style="invis")

    # ---- shared embedding engine --------------------------------------------
    add(g, "ENG", "Shared embedding engine", [f"Estimates within-condition {S2}"])

    # ---- two parallel columns inside a grouping box -------------------------
    with g.subgraph(name="cluster_emp") as emp:
        header(emp, "Cross-Modality Empirical Validation")
        with emp.subgraph(name="cluster_chem") as ch:
            header(ch, "Chemical modality")
            add(ch, "TA", "Tahoe-100M",
                [f"{S2} = 0.96", f"in-regime slope 0.94 &#183; {R2} &#8776; 0.999"])
            add(ch, "EB", "EmeraldBay",
                [f"{S2} &#8776; 0.94", f"held-out slope 0.98 &#183; {R2} &#8776; 0.999"])
            ch.edge("TA", "EB", label="  &#963;&#178; transfers &#8776; 0.94&#8211;0.96")
        with emp.subgraph(name="cluster_gen") as ge:
            header(ge, "Genetic modality")
            add(ge, "OR", "X-Atlas/Orion",
                [f"{S2} = 0.91 / 1.03", f"{N} &#8776; 9,000 / {M2} &#183; &#963;&#178; transfers"])
            add(ge, "TR", "TRADE",
                [f"{S2} = 1.5 / 1.9", f"slope 1.02 / 0.96 &#183; {R2} = 0.996"])
            ge.edge("OR", "TR", label="  direct falsification")

    # ---- strict horizontal alignment of the two columns ---------------------
    for a, b in (("TA", "OR"), ("EB", "TR")):
        with g.subgraph() as s:
            s.attr(rank="same")
            s.node(a)
            s.node(b)

    # ---- edges (spine ends at the validations) ------------------------------
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
