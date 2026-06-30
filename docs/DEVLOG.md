# Manuscript development log

Chronological record of how the manuscript and its tooling were built. Newest phase last. Every
quantitative claim is gated against the Constants of Record ([`SUPPLEMENT.md`](SUPPLEMENT.md)); every
reference is Crossref-verified; no language model was trusted as a source of fact.

## Pipeline (current)

```
docs/MANUSCRIPT_DEEPSEEK.md      (source: prose + \(..\)/\[..\] math + [@cite] keys, BiB section order)
        │  scripts/normalize_manuscript.py   (unify math to $..$, ASCII prose, fix pandoc edge cases)
        ▼
docs/manuscript.md               (compilable; YAML carries bibliography + Vancouver CSL)
        │  pandoc --citeproc      (csl/vancouver.csl + references.bib)
        ├─► manuscript.pdf        (pandoc + tectonic; 14 pp)
        └─► manuscript.docx       (native OMML equations + Word tables; editable)
```
Rebuild from repo root: `pandoc docs/manuscript.md -o manuscript.pdf --citeproc --pdf-engine=xelatex`
and `pandoc docs/manuscript.md -o manuscript.docx --citeproc`. Reference check: `python
scripts/verify_references.py` (32/32 DOIs).

## Phases

1. **Agent team** (`a23603e`). Built `agents/` — a fact-gated drafting pipeline: a deterministic
   SourcePacker injects the Constants of Record, DeepSeek `deepseek-v4-pro` writes each section, an
   adversarial Fact-Checker labels every claim SUPPORTED/UNSUPPORTED/DRIFTED against the pack, a Reviser
   repairs flagged claims, an Editor stitches. The run auto-caught real drifts (a wrong MSE formula, a
   per-well-vs-per-condition median); one section whose JSON came back empty was escalated to manual
   (Claude) audit. Output: `docs/MANUSCRIPT_DEEPSEEK.md`.

2. **Verified bibliography** (`11e4757`, later `d102493`). Every reference checked against the Crossref
   REST API (DOI + first author + year) by `scripts/verify_references.py` — never generated from model
   recall. Tahoe-100M (bioRxiv) and EmeraldBay (Hugging Face) verified against their primary records.

3. **Compilable + typeset** (`bd9ee41`, `d0adb91`, `9b07356`, `a7a7edc`). `normalize_manuscript.py`
   strips the editor's outer code fence, unifies three competing math notations to `$…$`/`$$…$$`, and
   removes breaking Unicode. Real compiles with pandoc + tectonic; rendered pages inspected and math-
   rendering bugs fixed (`$…$` adjacent to a digit, trailing space before a closing `$`, accented author
   names, `{,}` leaking into prose).

4. **Briefings in Bioinformatics format + 34 refs** (`d102493`). Reordered to Abstract → Introduction →
   Materials and methods → Results → Discussion → Conflicts/Funding/Data availability → References.
   Switched to pandoc `[@key]` + citeproc + a Vancouver CSL so citations auto-number in order of
   appearance and the reference list is generated; expanded 17 → **34** verified references (Connectivity
   Map/L1000, Perturb-seq variants, scPerturb, demuxlet/DoubletFinder, sctransform/Seurat/best-practices,
   scVI/CPA/GEARS, Johnstone, Vershynin, Mardia–Jupp).

5. **De-anonymization** (`f16a3a2`, `8afa36f`). Reversed the earlier generalization at the source: removed
   the `ALIAS` map from `agents/team.py` and the `REF_BY_*` relabelling in `scripts/dataset_inventory.py`
   (the actual cause the agent team mirrored), then propagated real names repo-wide — **Tahoe-100M**,
   **EmeraldBay**, and the cell lines **HS-578T, AN3-CA, HEC-1-A, BT-474, C-33 A**.

6. **Prose rewrites** (`f16a3a2`, `aa7c335`, `1d647b2`). Abstract + Introduction, then Results/Methods/
   Discussion, rewritten in neutral scientific register and made more application-oriented (a Methods
   "how to apply" recipe; a Discussion three-stage workflow; Results read as planning rules). The
   Introduction was then trimmed to motivation only — context → gap → high-level aim — with the
   derivation, the closed-form equation, and all result numbers moved out to Results/Methods.

7. **Editable Word output** (`8c375e3`). `manuscript.docx` built via pandoc; all math is **native Office
   Math (OMML)** equation objects (editable in Word, not images), plus native Word tables and the
   numbered reference list. Fixed a normalizer bug that double-escaped `\arccos` (LaTeX tolerated it; the
   docx math parser did not).

8. **Key Points + Limitations** (`646f1c1`). Added a Briefings in Bioinformatics-style **Key Points** box
   (5 high-level bullets) between the Abstract and Introduction. Rewrote the Limitations in the neutral
   register and removed the redundant inline "Limitations." paragraph the Discussion still carried,
   folding its unique points (Gaussian-centroid assumption; additive vs multiplicative/rotational batch)
   into the list — now 10 numbered items.

## Honesty ledger (carried in the manuscript banner + `SUPPLEMENT.md`)

- The two high-magnitude validation groups are DMSO time-zero populations, not drug effects; only one
  group (Encorafenib × HEC-1-A) is a drug — the held-out test validates geometry, not biology.
- EmeraldBay σ² ≈ 2.12 (within-condition) and Tahoe-100M σ² = 7.66 (marginal) are not a clean platform
  head-to-head; both only establish σ² as a plug-in.
- Gating tolerance in the EmeraldBay experiment is 0.20 rad (not the standard 0.1).
- Resource figures are Tahoe-100M-derived; compute gains are polynomial, never exponential.
- Author lists beyond the third name are abbreviated and should be expanded before submission.
