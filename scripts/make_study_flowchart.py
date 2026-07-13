"""
Generate the study flowchart (Graphviz) summarizing the whole pipeline:
  question -> closed-form quota -> theorem verification -> from-raw recompute
  -> chemical calibration/validation (Tahoe-100M, EmeraldBay)
  -> genetic transfer/falsification (X-Atlas/Orion, TRADE) -> two-ledger applications.

Numbers mirror the manuscript headline constants of record (within-condition calibration).
Writes figures/fig0_study_flowchart.{pdf,png} and figures/study_flowchart.dot.

Run:  python scripts/make_study_flowchart.py     (needs the `dot` CLI on PATH)
"""
import os, subprocess, shutil, sys

FIG = "figures"
os.makedirs(FIG, exist_ok=True)

# Okabe-Ito-derived tints (fill, border) — consistent with make_manuscript_figures.py
C = {
    "theory": ("#D6E6F5", "#0072B2"),   # blue   = core method
    "verify": ("#ECECEC", "#555555"),   # gray   = proof/verification
    "pipe":   ("#F6F6F6", "#8A8A8A"),   # neutral= shared processing
    "tahoe":  ("#D6E6F5", "#0072B2"),   # blue   = calibration anchor
    "emb":    ("#D6EFE2", "#009E73"),   # green  = out-of-distribution validation
    "orion":  ("#FBE7CF", "#E69F00"),   # orange = cross-modality transfer
    "trade":  ("#F1DEEA", "#CC79A7"),   # purple = genetic falsification
    "apps":   ("#E7F3EC", "#2E7D5B"),   # teal   = applications
    "q":      ("#FFFFFF", "#444444"),
}


def node(nid, kind, title, lines, shape="box"):
    fill, border = C[kind]
    detail = "".join(
        f'<BR/><FONT POINT-SIZE="9" COLOR="#333333">{ln}</FONT>' for ln in lines)
    label = f'<<B>{title}</B>{detail}>'
    return (f'  {nid} [shape={shape}, style="rounded,filled", '
            f'fillcolor="{fill}", color="{border}", penwidth=1.4, label={label}];\n')


dot = []
dot.append('digraph study {\n')
dot.append('  rankdir=TB; bgcolor="white"; splines=true; nodesep=0.45; ranksep=0.55;\n')
dot.append('  node [fontname="Helvetica", fontsize=11, margin="0.14,0.09"];\n')
dot.append('  edge [fontname="Helvetica", fontsize=9, color="#333333", penwidth=1.4, arrowsize=0.8];\n')

# --- backbone entry: the motivating question + method ---
dot.append(node("Q", "q", "Question",
    ["How many cells per arm resolve a perturbation's",
     "<I>direction</I> (mechanism) to an angular tolerance θ?"]))
dot.append(node("TH", "theory", "Closed-form cell quota  (Δ-method geometry)",
    ["n★ = 2·tr(PΣP) / (m² θ²),   P = I − uuᵀ",
     "only noise ⊥ to the effect rotates its direction",
     "+ tail bound gives a confidence quota (level 1−δ)"]))
dot.append(node("VE", "verify", "Established as a theorem: four independent ways",
    ["symbolic (SymPy) · Monte-Carlo (&lt;0.13%)",
     "independent implementation",
     "Lean 4 / Mathlib deterministic core (no <I>sorry</I>)"]))
dot.append(node("PI", "pipe", "Shared from-raw recompute (every dataset)",
    ["normalize → log1p → HVG(2000) → PCA(50)",
     "within-condition Σ → σ²   ·   centroids → magnitude m"]))

# --- chemical modality cluster ---
dot.append('  subgraph cluster_chem {\n')
dot.append('    label=<<B>Chemical modality (Vevo Mosaic platform)</B>>; fontname="Helvetica"; fontsize=10;\n')
dot.append('    labeljust="l"; style="rounded,filled"; fillcolor="#F1F8FF"; pencolor="#8FB9DD"; penwidth=1.6; margin=14;\n')
dot.append(node("TA", "tahoe", "Tahoe-100M:  calibration anchor",
    ["100.6M cells · 379 × 3 × 50 = 56,827 conditions",
     "σ² = 0.9567  →  n★ = 9,376 / m²   (θ = 0.1 rad)",
     "spectrum: 10.6% over · 89.0% under · 0.4% ghost"]))
dot.append(node("EB", "emb", "EmeraldBay:  out-of-distribution validation",
    ["1.83M cells · 52 lines · frozen HVG basis",
     "within σ² ≈ 0.938  (matches Tahoe to ~2%)",
     "held-out angular-error slope R² ≈ 0.999 · gating 347/347"]))
dot.append('    { rank=same; TA; EB; }\n')
dot.append('  }\n')

# --- genetic modality cluster ---
dot.append('  subgraph cluster_gen {\n')
dot.append('    label=<<B>Genetic modality (CRISPRi Perturb-seq)</B>>; fontname="Helvetica"; fontsize=10;\n')
dot.append('    labeljust="l"; style="rounded,filled"; fillcolor="#FFFBF2"; pencolor="#E0C48F"; penwidth=1.6; margin=14;\n')
dot.append(node("OR", "orion", "X-Atlas/Orion:  cross-modality transfer",
    ["genome-wide CRISPRi · ~8M cells · 18,903 KD × 2 lines",
     "σ² = 0.91 / 1.03 transfers  →  n★ ≈ 9,000 / m²",
     "0% over-sampled (knockdowns move ~4–10× less than drugs)"]))
dot.append(node("TR", "trade", "TRADE:  validated genetic falsification",
    ["essential-gene CRISPRi · Jurkat + HepG2 · 2,393 genes/line",
     "realized vs a-priori slope 1.02 / 0.96 · R² = 0.996",
     "strong but shallowly sampled → depth-limited regime"]))
dot.append('    { rank=same; OR; TR; }\n')
dot.append('  }\n')

# --- applications ---
dot.append(node("AP", "apps", "Two ledgers, one threshold n★",
    ["wet-lab: magnitude-aware budget (~3.1× fewer cells) · triage / multiplex",
     "dry-lab: provably-safe downsampling (polynomial compute savings)",
     "reliability audit of MoA / drug-similarity graph"]))

# --- backbone edges (solid) ---
dot.append('  Q  -> TH;\n')
dot.append('  TH -> VE [label="  proven"];\n')
dot.append('  VE -> PI;\n')
dot.append('  PI -> TA [label="  calibrate"];\n')
dot.append('  TA -> EB [label="out-of-distribution test"];\n')
dot.append('  TA -> OR [label="  transfer σ² across modality"];\n')
dot.append('  OR -> TR [label="direct falsification"];\n')
dot.append('  TA -> AP [label="  sufficiency spectrum"];\n')

# --- shared-pipeline provenance (dashed, non-structural) ---
dot.append('  edge [style=dashed, color="#9AA0A6", penwidth=1.0, arrowsize=0.7];\n')
dot.append('  PI -> EB [constraint=false];\n')
dot.append('  PI -> OR [constraint=false];\n')
dot.append('  PI -> TR [constraint=false];\n')
dot.append('  TR -> AP [constraint=false];\n')
dot.append('}\n')

dot_src = "".join(dot)
dot_path = os.path.join(FIG, "study_flowchart.dot")
with open(dot_path, "w") as f:
    f.write(dot_src)

if not shutil.which("dot"):
    sys.exit("graphviz `dot` not found on PATH")

for ext, args in (("pdf", ["-Tpdf"]), ("png", ["-Tpng", "-Gdpi=200"])):
    out = os.path.join(FIG, f"fig0_study_flowchart.{ext}")
    subprocess.run(["dot", *args, dot_path, "-o", out], check=True)
    print(f"  wrote {out}")
print(f"  wrote {dot_path}")
