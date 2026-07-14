"""Structural (AST) guard — a second, regex-independent line of defence for the single-source-of-truth
rule. Three legacy files (pass6_invariance, then pass4c_gating, then pass5_gating_full / pass3_validate /
moa_recovery) each fed equal-arm quota numbers into the paper and each slipped past a *string* guard
(case, spelling). This guard reasons over the parse tree instead:

A file under scripts/ VIOLATES if it either
  (1) DEFINES a function whose name looks like a sufficiency quantity (quota / nstar / n_star / detect /
      bias / mmin / m_min / regime) AND that function body contains a quota formula, or
  (2) ASSIGNS to a variable named like a sufficiency quantity (nstar / n_star / quota / regime / m_min /
      mmin / detectable / m2_corr / snr) FROM A QUOTA FORMULA
where "quota formula" = an arithmetic expression (contains `** 2`) that references an angular tolerance
(theta/th/theta_gate/theta_star) together with a perpendicular-noise/variance term (trPSP/trSig/sigma2/
sig/ell) or a magnitude (m/m2/m_raw). Reading a value back from the engine -- `nstar = res["n_star"][0]`,
`reg = res["regime"]` -- is a Subscript, not a formula, so it is allowed. src/engine.py is the sole owner
and is not scanned; src/calculator.py is the sanctioned standalone library.
"""
import ast, glob, os

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))
SUSPECT_TARGETS = {"nstar", "n_star", "quota", "regime", "m_min", "mmin", "detectable", "m2_corr", "snr"}
FUNC_SUBSTR = ("quota", "nstar", "n_star", "detect", "bias", "mmin", "m_min", "regime")
THETA_NAMES = {"theta", "theta_star", "th", "theta_gate", "tol", "tolerance"}
PERP_NAMES = {"trpsp", "tr_psp", "trsig", "trsigma", "sigma2", "sig", "ell", "psp"}
MAG_NAMES = {"m", "m2", "m_raw", "mag", "m2_corr", "m_corr"}


def _name_ids(node):
    return {n.id.lower() for n in ast.walk(node) if isinstance(n, ast.Name)}


def _has_pow2(node):
    for n in ast.walk(node):
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Pow):
            r = n.right
            if isinstance(r, ast.Constant) and r.value == 2:
                return True
    return False


def _is_quota_formula(value):
    """An arithmetic expression solving for cells: has **2 AND a tolerance AND a perp/variance-or-magnitude."""
    if not _has_pow2(value):
        return False
    names = _name_ids(value)
    has_theta = bool(names & THETA_NAMES)
    has_perp_or_mag = bool(names & PERP_NAMES) or bool(names & MAG_NAMES)
    return has_theta and has_perp_or_mag


ALLOWED = {"src/engine.py"}                                # the sole definition site of the quota formula


def find_ast_violations():
    viol = []
    files = (glob.glob(os.path.join(ROOT, "scripts", "**", "*.py"), recursive=True)
             + glob.glob(os.path.join(ROOT, "src", "**", "*.py"), recursive=True))
    for f in sorted(files):
        rel = os.path.relpath(f, ROOT).replace(os.sep, "/")
        if rel in ALLOWED:
            continue
        try:
            tree = ast.parse(open(f, encoding="utf-8", errors="ignore").read())
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                if any(s in node.name.lower() for s in FUNC_SUBSTR) and \
                        any(_is_quota_formula(b) for b in ast.walk(node) if isinstance(b, (ast.BinOp,))):
                    viol.append(f"  {rel}:{node.lineno}: function '{node.name}' computes a quota formula locally")
            if isinstance(node, ast.Assign):
                tgts = {t.id.lower() for t in node.targets if isinstance(t, ast.Name)}
                if (tgts & SUSPECT_TARGETS) and _is_quota_formula(node.value):
                    viol.append(f"  {rel}:{node.lineno}: assigns {sorted(tgts & SUSPECT_TARGETS)} from a quota formula")
    return viol


def test_scripts_do_not_recompute_the_quota():
    viol = find_ast_violations()
    assert not viol, ("scripts/ must delegate the quota/regime to src/engine.py (read res['n_star'] / "
                      "res['regime']), not recompute it:\n" + "\n".join(viol))


if __name__ == "__main__":
    v = find_ast_violations()
    print(f"{len(v)} AST violation(s):")
    print("\n".join(v))
