"""
Generate the study flowchart (Graphviz) — the redesigned, container-structured overview:
  research question -> closed-form quota -> AIRTIGHT PROOF box (Lean 4, SymPy, Monte-Carlo,
  independent impl.) -> shared embedding ENGINE -> CROSS-MODALITY VALIDATION container with two
  parallel modality columns (chemical: Tahoe-100M -> EmeraldBay ; genetic: X-Atlas/Orion -> TRADE)
  -> DOWNSTREAM VALUE split into wet-lab and dry-lab ledgers.

Muted, professional palette (cool grays, muted blues, pale teals). Numbers mirror the manuscript
headline constants of record. Writes figures/fig1_study_flowchart.{pdf,png} and figures/study_flowchart.dot.
The equivalent Mermaid.js source (for web / GitHub rendering) lives at docs/study_flowchart.mmd.

Run:  python scripts/make_study_flowchart.py     (needs the `dot` CLI on PATH)
"""
import os, subprocess, shutil, sys

FIG = "figures"
os.makedirs(FIG, exist_ok=True)

# muted palette: (fill, border, text)
C = {
    "q":      ("#eef1f4", "#7d8b9a", "#26313d"),   # cool gray
    "theory": ("#e6eef6", "#5b7a9d", "#1f2d3d"),   # muted blue
    "proof":  ("#eceef1", "#95a0ad", "#2b333c"),   # cool gray
    "engine": ("#e0edea", "#5f9a90", "#22322f"),   # pale teal
    "chem":   ("#eaf1f7", "#6f8caa", "#22303d"),   # muted blue (data table)
    "gen":    ("#e4efec", "#6f9a90", "#22322f"),   # pale teal (data table)
    "hub":    ("#e6eef6", "#5b7a9d", "#1f2d3d"),
    "wet":    ("#e4eef2", "#5f8aa0", "#1f2f38"),
    "dry":    ("#eceef1", "#8a97a6", "#2b333c"),
}


def node(nid, kind, title, lines, shape="box", rounded=True, pen=1.3):
    fill, border, txt = C[kind]
    detail = "".join(
        f'<BR/><FONT POINT-SIZE="9" COLOR="#4a545e">{ln}</FONT>' for ln in lines)
    label = f'<<B>{title}</B>{detail}>'
    style = "rounded,filled" if rounded else "filled"
    return (f'  {nid} [shape={shape}, style="{style}", fillcolor="{fill}", '
            f'color="{border}", fontcolor="{txt}", penwidth={pen}, label={label}];\n')


d = []
d.append('digraph study {\n')
d.append('  compound=true; rankdir=TB; bgcolor="white"; splines=true; nodesep=0.4; ranksep=0.5;\n')
d.append('  node [fontname="Helvetica", fontsize=11, margin="0.16,0.10"];\n')
d.append('  edge [fontname="Helvetica", fontsize=9, color="#5f6b78", penwidth=1.2, arrowsize=0.75];\n')

# top anchor + theory
d.append(node("Q", "q", "Core research question",
    ["How many cells per arm resolve a perturbation's <I>direction</I>",
     "(mechanism) to an angular tolerance θ?"], shape="box"))
d.append(node("TH", "theory", "Closed-form cell quota",
    ["n★ = 2·tr(PΣP) / (m² θ²),   P = I − uuᵀ",
     "only noise ⊥ to the effect rotates the direction",
     "calibrated → n★ ≈ 9,376 / m²   (θ = 0.1 rad)"]))

# section 1 — airtight proof foundation
d.append('  subgraph cluster_proof {\n')
d.append('    label=<<B>Airtight foundation</B>  ·  proven four independent ways>;\n')
d.append('    fontname="Helvetica"; fontsize=10; labeljust="l"; style="rounded,filled";\n')
d.append('    fillcolor="#f3f5f7"; pencolor="#aeb9c4"; penwidth=1.4; margin=12;\n')
d.append(node("P1", "proof", "Lean 4 / Mathlib", ["deterministic core, no <I>sorry</I>"]))
d.append(node("P2", "proof", "SymPy", ["symbolic Jacobian"]))
d.append(node("P3", "proof", "Monte-Carlo", ["&lt; 0.13% error"]))
d.append(node("P4", "proof", "Independent impl.", ["cross-check"]))
d.append('    { rank=same; P1; P2; P3; P4; }\n')
d.append('  }\n')

# section 2 — embedding engine (sleek horizontal bridge)
d.append(node("ENG", "engine", "Shared embedding engine",
    ["Normalize → log1p → HVG 2000 → PCA 50      ·      "
     "within-condition Σ → σ²      ·      centroids → magnitude m"]))

# section 3 — cross-modality validation container, two parallel columns
d.append('  subgraph cluster_emp {\n')
d.append('    label=<<B>Cross-Modality Empirical Validation</B>>;\n')
d.append('    fontname="Helvetica"; fontsize=11; labeljust="l"; style="rounded,filled";\n')
d.append('    fillcolor="#f5f7f9"; pencolor="#aeb9c4"; penwidth=1.4; margin=16;\n')
d.append('    subgraph cluster_chem {\n')
d.append('      label=<<B>Chemical modality</B>  ·  Vevo Mosaic>;\n')
d.append('      fontname="Helvetica"; fontsize=10; labeljust="l"; style="rounded,filled";\n')
d.append('      fillcolor="#eef3f8"; pencolor="#8ea6bd"; penwidth=1.3; margin=12;\n')
d.append(node("TA", "chem", "Tahoe-100M  ·  calibration anchor",
    ["100.6M cells · 56,827 conditions",
     "σ² = 0.96 → n★ = 9,376 / m²",
     "spectrum: 10.6% over / 89.0% under / 0.4% ghost"], rounded=False))
d.append(node("EB", "chem", "EmeraldBay  ·  out-of-distribution validation",
    ["1.83M cells · 52 lines",
     "σ² ≈ 0.94 · held-out slope R² ≈ 0.999",
     "gating 347 / 347 met"], rounded=False))
d.append('      TA -> EB [label="σ² transfers ≈ 0.94–0.96"];\n')
d.append('    }\n')
d.append('    subgraph cluster_gen {\n')
d.append('      label=<<B>Genetic modality</B>  ·  CRISPRi Perturb-seq>;\n')
d.append('      fontname="Helvetica"; fontsize=10; labeljust="l"; style="rounded,filled";\n')
d.append('      fillcolor="#eef4f2"; pencolor="#8bb0a6"; penwidth=1.3; margin=12;\n')
d.append(node("OR", "gen", "X-Atlas/Orion  ·  cross-modality transfer",
    ["genome-wide CRISPRi · ~8M cells",
     "18,903 knockdowns × 2 lines",
     "σ² = 0.91 / 1.03 · median m 0.14–0.35",
     "0% over-sampled (magnitude-limited)"], rounded=False))
d.append(node("TR", "gen", "TRADE  ·  validated genetic falsification",
    ["essential-gene CRISPRi · Jurkat + HepG2",
     "2,393 genes/line · slope 1.02 / 0.96 · R² = 0.996",
     "median 45–85 cells/gene → depth-limited"], rounded=False))
d.append('      OR -> TR [label="direct falsification"];\n')
d.append('    }\n')
d.append('  }\n')

# bottom anchor — downstream value split
d.append(node("HUB", "hub", "One threshold n★, two ledgers", []))
d.append('  subgraph cluster_out {\n')
d.append('    label=<<B>Downstream Value</B>>;\n')
d.append('    fontname="Helvetica"; fontsize=11; labeljust="l"; style="rounded,filled";\n')
d.append('    fillcolor="#f3f5f7"; pencolor="#aeb9c4"; penwidth=1.4; margin=14;\n')
d.append(node("WET", "wet", "Wet-lab",
    ["magnitude-adaptive budgeting", "~3.1× fewer cells · triage / multiplex"]))
d.append(node("DRY", "dry", "Dry-lab",
    ["provably-safe downsampling", "triage of low-SNR similarity edges"]))
d.append('    { rank=same; WET; DRY; }\n')
d.append('  }\n')

# edges
d.append('  Q  -> TH;\n')
d.append('  TH -> P2 [lhead=cluster_proof, label="  proven"];\n')
d.append('  P2 -> ENG [ltail=cluster_proof];\n')
d.append('  ENG -> TA [lhead=cluster_chem, label="  calibrate"];\n')
d.append('  ENG -> OR [lhead=cluster_gen];\n')
d.append('  EB -> HUB [ltail=cluster_chem];\n')
d.append('  TR -> HUB [ltail=cluster_gen];\n')
d.append('  HUB -> WET [lhead=cluster_out];\n')
d.append('  HUB -> DRY [lhead=cluster_out];\n')
d.append('}\n')

dot_src = "".join(d)
dot_path = os.path.join(FIG, "study_flowchart.dot")
with open(dot_path, "w") as f:
    f.write(dot_src)

if not shutil.which("dot"):
    sys.exit("graphviz `dot` not found on PATH")
for ext, args in (("pdf", ["-Tpdf"]), ("png", ["-Tpng", "-Gdpi=200"])):
    out = os.path.join(FIG, f"fig1_study_flowchart.{ext}")
    subprocess.run(["dot", *args, dot_path, "-o", out], check=True)
    print(f"  wrote {out}")
print(f"  wrote {dot_path}")
