"""Guard against the root cause of audit findings #1-#9: the sufficiency math must live in ONE place
(src/engine.py). Any OTHER module -- in src/ OR scripts/ -- that recomputes the quota n*, the
bias-correction m^2 - tr(S), the snr>1.5 detection floor, or hardcodes a dead per-modality constant is
a duplication that can drift the way pass6_invariance.py did (it fed the equal-arm n*=184 straight into
the manuscript while the src/-only guard never looked at scripts/).

This guard scans src/ AND scripts/ (every *.py except src/engine.py and this test). It matches by
MEANING, not by the exact spelling of engine.py, on comment/string-stripped code:

  * QUOTA formula   -- a perpendicular-noise trace (tr(P.Sigma.P)) divided by a squared magnitude and
                       theta^2 (any spelling): this is solving for n*, which only the engine may do.
                       The held-out MEASUREMENT law tr(P.Sigma.P)/m^2 * (1/n - 1/N) is NOT this: it
                       multiplies by (1/n - 1/N) and never divides by theta^2, so it is allowed.
  * BIAS-correction -- max(0, m^2 - tr(S)) in any spelling.
  * DETECTION floor -- snr > <const>, or m/sqrt(tr(S)) compared to a constant.
  * DEAD constants  -- 9376 / 9192 (old isotropic quotas) and 4688 / 4495 / 4475 (large-pool C-bounds
                       that must be COMPUTED by the engine, never hardcoded), and the baseless
                       detection/threshold constants 1.25 / 2.25 and the N>=300 cell-count filter.

A single legitimate exception is whitelisted by exact (file, line-substring): the held-out subsample
grid in the falsification uses 300 as a SUBSAMPLE-DEPTH floor, not a quota.
"""
import glob, io, os, re, tokenize

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
SCAN_DIRS = ("src", "scripts")
ALLOWED_FILES = {"src/engine.py"}                    # the sole owner of the quota math

# exact (relpath, needle-in-line) pairs that are legitimate despite matching a DEAD-constant pattern
LINE_WHITELIST = {
    # the held-out falsification's subsample-depth floor (>=300 cells to fit a slope), not a quota
    ("scripts/tahoe_recompute/pass4_falsification.py", "300"),
    ("scripts/tahoe_recompute/pass4b_curves.py", "300"),
    ("scripts/tahoe_recompute/pass4c_gating.py", "300"),
}

# NOTE: all case-insensitive -- THETA/Theta/theta and trPSP spellings all count (a THETA-only quota,
# 2.0*trPSP/(m**2*THETA**2), slipped past an earlier lowercase-only version of this guard).
PERP = r"(?i:trPSP|tr_?PSP|np\.trace\s*\(\s*P|tr\s*\(\s*P|P\s*@\s*Sig|Sig\w*\s*@\s*P|sum\w*\(\s*\(?\s*1\s*-\s*u)"
THETA = r"(?i:theta|th)\b"
MSQ = r"(?i:m\w*\s*\*\*\s*2)"
# solving-for-n: divide by theta^2 (quota), as opposed to multiplying by (1/n-1/N) (measurement law)
DIV_THETA2 = r"(?i:/\s*\(?[^)\n]{0,60}?theta\s*\*\*\s*2|theta\s*\*\*\s*2\s*\)|/\s*\(?[^)\n]{0,60}?th\s*\*\*\s*2)"

RE_BIAS = re.compile(r"m\w*\s*\*\*\s*2\s*-\s*(?:np\.)?tr|maximum\s*\(\s*0[.,][^)]*-\s*(?:np\.)?tr|\bm2_corr\b")
RE_DETECT = re.compile(r"\bsnr\b\s*>|\bsnr_detect\b|/\s*np\.sqrt\s*\(\s*trS")
RE_DEAD_NUM = re.compile(r"\b(?:9376|9192|4688|4495|4475)\b")
# dead thresholds/factors: the baseless snr factor 1.25 / 2.25 and the N>=300 cell-count filter,
# in COMPARISON or MULTIPLIER position (not as an incidental plot coordinate literal).
RE_DEAD_THR = re.compile(r"[<>]=?\s*(?:300|1\.25|2\.25)\b|\b(?:300|1\.25|2\.25)\s*[<>*]|"
                         r"min[_-]?cells?\s*[=:]\s*300\b")


def code_only(path):
    """Source with comments and string literals (docstrings) removed."""
    out = []
    try:
        with open(path, "rb") as fh:
            for tok in tokenize.tokenize(fh.readline):
                if tok.type in (tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.NEWLINE,
                                tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING):
                    continue
                out.append(tok.string)
    except (tokenize.TokenError, SyntaxError):
        return open(path, encoding="utf-8", errors="ignore").read()
    return " ".join(out)


def _scan_files():
    files = []
    for d in SCAN_DIRS:
        files += glob.glob(os.path.join(ROOT, d, "**", "*.py"), recursive=True)
    return sorted(files)


def find_violations():
    viol = []
    for f in _scan_files():
        rel = os.path.relpath(f, ROOT).replace(os.sep, "/")
        if rel in ALLOWED_FILES:
            continue
        code = code_only(f)
        raw = open(f, encoding="utf-8", errors="ignore").read()

        if re.search(PERP, code) and re.search(THETA, code) and re.search(MSQ, code) \
                and re.search(DIV_THETA2, code):
            viol.append((rel, "QUOTA formula (tr(PSP) / (m^2 theta^2)) -- migrate to engine.compute"))
        if RE_BIAS.search(code):
            viol.append((rel, "BIAS-correction (m^2 - tr(S)) -- engine only"))
        if RE_DETECT.search(code):
            viol.append((rel, "DETECTION floor (snr > const) -- engine only"))
        if RE_DEAD_NUM.search(code):
            viol.append((rel, f"DEAD constant {RE_DEAD_NUM.search(code).group(0)} -- compute in engine"))
        for m in RE_DEAD_THR.finditer(code):
            tok = m.group(0)
            # whitelist exact legitimate lines
            wl = any(rel == wf and wn in tok for wf, wn in LINE_WHITELIST) or \
                any(rel == wf and wn in raw and ("300" in tok) for wf, wn in LINE_WHITELIST)
            if not wl:
                viol.append((rel, f"DEAD/threshold constant '{tok.strip()}' -- baseless filter/const"))
    # de-dup (file, reason)
    seen, out = set(), []
    for v in viol:
        if v not in seen:
            seen.add(v); out.append(v)
    return out


def test_no_duplicate_sufficiency_math():
    viol = find_violations()
    msg = "\n".join(f"  {f}: {r}" for f, r in viol)
    assert not viol, f"legacy sufficiency math outside src/engine.py:\n{msg}"


# ---- retired-thesis constants must not reappear in prose/code (docs + scripts + src) ----
# Scanned in .md and .py only (data .json/.csv hold engine floats whose digits coincide with these).
# Exempt: the audit/architecture/readme docs (which document the retired numbers) and this guard.
DEAD_DOC_EXEMPT = ("docs/audit_log.md", "docs/audit_denominators.md", "docs/architecture.md",
                   "readme.md", "tests/test_no_duplicate_math.py")
# comma-formatted old constants + old spectrum phrases + retired exemplars/thresholds
RE_DEAD_TEXT = re.compile(
    r"9,376|9,192|23,577|18,146|14,570|5,794|"                       # retired equal-arm/marginal constants
    r"\b9376\b|\b9192\b|\b23577\b|\b14570\b|"                        # ... as bare code literals
    r"10\.6\s*%?\s*over|89\.0\s*%?\s*under|0\.4\s*%?\s*ghost|"       # retired within-condition spectrum
    r"89\.3\s*%?\s*under|2\.5\s*%?\s*over|8\.2\s*%?\s*ghost|"        # retired marginal spectrum
    r"n_d\s*=\s*130|n\*\s*=\s*184|n\*\s*=\s*130")                    # retired invariance exemplar quotas


def test_no_retired_thesis_constants_in_prose_or_code():
    files = glob.glob(os.path.join(ROOT, "docs", "*.md")) + _scan_files()
    viol = []
    for f in files:
        rel = os.path.relpath(f, ROOT).replace(os.sep, "/")
        if rel.lower() in DEAD_DOC_EXEMPT:
            continue
        for i, line in enumerate(open(f, encoding="utf-8", errors="ignore"), 1):
            m = RE_DEAD_TEXT.search(line)
            if m:
                viol.append(f"  {rel}:{i}: retired-thesis token '{m.group(0)}'")
    assert not viol, ("retired-thesis numbers reappeared (only docs/AUDIT_LOG.md, "
                      "AUDIT_DENOMINATORS.md, ARCHITECTURE.md, README.md may contain them):\n"
                      + "\n".join(viol))


if __name__ == "__main__":
    v = find_violations()
    print(f"{len(v)} violation(s):")
    for f, r in v:
        print(f"  {f}: {r}")
