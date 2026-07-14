"""
Guard against the root cause of audit findings #1-#9: the sufficiency math must live in ONE place
(src/engine.py). If any other module recomputes m^2 - tr(S), the snr>1.5 detection threshold, or the
quota n*, or hardcodes a per-modality constant, the two branches can drift apart again.

The checks inspect CODE only (docstrings and comments are stripped), so prose that merely mentions a
formula does not trip them.

Fails if:
  1. A data pipeline (src/pipelines/*.py) contains engine math -- it must only assemble the standard
     structure (mu_t, mu_c, n_t, n_c, Sigma) and hand it to the engine.
  2. The two-arm quota formula appears in any src/ module other than engine.py.
  3. A live src/ module reintroduces the dead isotropic constants 9,376 / 9,192.
"""
import glob, io, os, tokenize

ROOT = os.path.join(os.path.dirname(__file__), "..")


def code_only(path):
    """Source with comments and string literals (docstrings) removed."""
    out = []
    with open(path, "rb") as fh:
        for tok in tokenize.tokenize(fh.readline):
            if tok.type in (tokenize.COMMENT, tokenize.STRING, tokenize.NL, tokenize.NEWLINE,
                            tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING):
                continue
            out.append(tok.string)
    return " ".join(out)


def _src_py():
    return glob.glob(os.path.join(ROOT, "src", "**", "*.py"), recursive=True)


def test_pipelines_assemble_only_no_math():
    math_tokens = ["m_raw", "n_star", "trPSP", "snr", "m2_corr", "detectable"]
    for f in glob.glob(os.path.join(ROOT, "src", "pipelines", "*.py")):
        code = code_only(f)
        hits = [t for t in math_tokens if t in code]
        assert not hits, f"{os.path.basename(f)} computes {hits}; pipelines only assemble data, math is in engine.py"


def test_engine_is_sole_owner_of_quota_formula():
    owners = [os.path.relpath(f, ROOT) for f in _src_py()
              if "theta ** 2 / trPSP" in open(f).read() or "m2_corr * theta" in open(f).read()]
    assert owners == ["src/engine.py"], f"quota formula found outside engine.py: {owners}"


def test_no_dead_isotropic_constants_in_src_code():
    for f in _src_py():
        code = code_only(f)
        for c in ["9376", "9192"]:
            assert c not in code, f"{os.path.relpath(f, ROOT)} reintroduced dead constant {c} in code"
