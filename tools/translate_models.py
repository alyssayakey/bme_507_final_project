"""Reproducibly translate the two upstream straight-line C equation kernels.

This deliberately supports only the C syntax present in the pinned sources.
Generated modules contain ordinary Python, with optional Numba acceleration.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def translate(source, target, size, outputs):
    text = source.read_text(encoding="utf-8-sig")
    params = re.findall(r"#define\s+(\w+)\s+parms\[(\d+)\]", text)
    body = text.split("void derivs", 1)[1].split("{", 1)[1]
    body = re.sub(r"//[^\n]*", "", body)
    body = re.sub(r'if\s*\(ip\[0\].*?;', '', body)
    body = re.sub(r'time_t s = time\(NULL\);', '', body)
    body = re.sub(r'if\(\(int\) s.*?;', '', body)
    body = body.rsplit("}", 1)[0]
    body = re.sub(r"\bdouble\s+", "", body)
    body = re.sub(r"\*t\b", "t", body)
    # D=0 has the exact continuous limit exp(n*log(D))=0 for n>0.
    body = body.replace("exp(n*log(D))", "D**n")
    body = body.replace("exp(n*log(y[9]))", "y[9]**n")
    lines = ['"""Generated Python equations derived from FDA/CiPA; see licenses/."""',
             "import numpy as np", "from numpy import exp, log, sqrt, power as pow",
             "from .acceleration import njit", "", f"PARAMETER_NAMES = {tuple(p[0] for p in params)!r}",
             "", '@njit(cache=True, error_model="numpy")', "def evaluate(t, y, p):",
             f"    ydot = np.zeros({size})", f"    yout = np.zeros({outputs})"]
    lines += [f"    {name} = p[{idx}]" for name, idx in params]
    level = 1
    single = False
    for token in re.split(r"([{};])", body):
        token = token.strip()
        if not token:
            continue
        if token == "}":
            level -= 1
        elif token == "{":
            level += 1
        elif token == ";":
            if single:
                level -= 1
                single = False
        elif token.startswith(("if", "else")):
            unbraced = re.match(r"if\s*\(([^\n]+)\)\s*\n(.*)", token, re.S)
            if unbraced:
                lines += ["    "*level + "if " + unbraced[1] + ":",
                          "    "*(level+1) + unbraced[2].strip()]
                continue
            match = re.fullmatch(r"(else\s+)?if\s*\((.*)\)", token, re.S)
            if match:
                line = ("elif" if match[1] else "if") + " " + match[2] + ":"
            elif token == "else":
                line = "else:"
            else:
                raise ValueError(f"Unsupported condition: {token}")
            lines.append("    " * level + line)
            # The only unbraced conditional in these sources.
            if token == "if (celltype==1)\n    GNaL=GNaL*0.6":
                raise ValueError("Unexpected token")
        else:
            # Separate an unbraced if from its assignment.
            match = re.match(r"if\s*\(([^\n]+)\)\s*\n(.*)", token, re.S)
            if match:
                lines += ["    "*level + "if " + match[1] + ":",
                          "    "*(level+1) + match[2].strip()]
            else:
                lines.append("    "*level + token)
    lines.append("    return ydot, yout")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    for folder, model, size, outputs, target in [
        ("hERG_fitting", "hergmod", 11, 0, "_herg_equations.py"),
        ("AP_simulation", "newordherg_qNet", 50, 8, "_ord_equations.py"),
    ]:
        translate(ROOT / "reference/CiPA" / folder / "models" / (model + ".c"),
                  ROOT / "cipa" / target, size, outputs)
