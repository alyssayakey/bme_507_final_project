"""Compile upstream C without R and compare derivatives against the Python port.

Run: python tools/validate_c_models.py
Development-only compiler dependency: pip install ziglang
The Python package itself does not call C or need this compiler.
"""
from pathlib import Path
import ctypes
import json
import re
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cipa.models import herg, ord, initial_ap, ap_parameters, parameter_vector
from cipa.common import DATA, read_pairs


def compile_reference(folder, model, size, nout):
    source = (ROOT / "reference/CiPA" / folder / "models" / (model + ".c")).read_text(encoding="utf-8-sig")
    source = re.sub(r"#include[^\n]*", "", source)
    source = re.sub(r'if\s*\(ip\[0\].*?;', '', source)
    source = re.sub(r'time_t s = time\(NULL\);', '', source)
    source = re.sub(r'if\(\(int\) s.*?;', '', source)
    source = '#include <math.h>\n' + source
    source += '''
#ifdef _WIN32
__declspec(dllexport)
#endif
void reference_evaluate(double t, double *y, double *p, double *dy, double *out) {
  for (int i=0; i<sizeof(parms)/sizeof(double); i++) parms[i]=p[i];
  int neq=NSTATES; int ip=NOUT;
  derivs(&neq, &t, y, dy, out, &ip);
}
'''.replace("NSTATES", str(size)).replace("NOUT", str(nout))
    build = ROOT / "build/validation"
    build.mkdir(parents=True, exist_ok=True)
    c = build / (model + ".c")
    dll = build / (model + (".dll" if sys.platform == "win32" else ".so"))
    c.write_text(source, encoding="utf-8")
    subprocess.run([sys.executable, "-m", "ziglang", "cc", "-shared", "-O2", str(c), "-o", str(dll)], check=True)
    lib = ctypes.CDLL(str(dll))
    func = lib.reference_evaluate
    ptr = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
    func.argtypes = [ctypes.c_double, ptr, ptr, ptr, ptr]
    return lib, func


def main():
    rng = np.random.default_rng(10)
    report = {}
    for folder, model, module, size, nout in [("hERG_fitting", "hergmod", herg, 11, 0),
                                            ("AP_simulation", "newordherg_qNet", ord, 50, 8)]:
        lib, func = compile_reference(folder, model, size, nout)
        errors = []
        for i in range(60):
            if size == 11:
                p = read_pairs(DATA / "herg/hergmod_pars.txt")
                p.update(Kmax=5., Ku=.001, n=1.2, halfmax=100., Vhalf=-8.)
                y = np.r_[rng.dirichlet(np.ones(9)), 0. if i % 3 == 0 else 100., rng.uniform(-90, 40)]
            else:
                p = ap_parameters()
                p.update(celltype=i % 3, Kmax=5., Ku=.001, n=1.2, halfmax=100., Vhalf=-8.)
                y = initial_ap()
                y[0] = rng.uniform(-90, 40)
                y[1:9] *= rng.uniform(.8, 1.2, 8)
                y[9:36] = rng.uniform(.01, .99, 27)
                y[39:48] = rng.dirichlet(np.ones(9))
                y[48] = 0. if i % 3 == 0 else 100.
            vector = parameter_vector(p, "herg" if size == 11 else "ap")
            time = 0. if i % 2 else 1.
            dy, out = np.zeros(size), np.zeros(nout)
            func(time, y, vector, dy, out)
            py_dy, py_out = module.evaluate(time, y, vector)
            np.testing.assert_allclose(py_dy, dy, rtol=1e-11, atol=1e-10)
            np.testing.assert_allclose(py_out, out, rtol=1e-11, atol=1e-10)
            errors.append(float(np.max(np.abs(dy-py_dy))))
        report[model] = dict(cases=len(errors), max_absolute_error=max(errors))
    path = ROOT / "validation_c_models.json"
    path.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
