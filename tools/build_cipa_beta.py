"""Build the hybrid from reviewed local CiPA/Gong blocks, without editing either.

Run from repository root: python tools/build_cipa_beta.py
The generated kernel uses CiPA geometry, baseline scaling, and Markov IKr.
"""
from pathlib import Path
import ast
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "cipa_beta"
EXTRA = ('hPf', 'hPs', 'jP', 'hBPf', 'hBPs', 'jBP', 'dP', 'fPf', 'fPs',
         'fcaPf', 'fcaPs', 'fBPf', 'fcaBPf', 'xs1P', 'xs2P', 'JrelP', 'JrelBP')


def between(text, start, end):
    return text[text.index(start):text.index(end, text.index(start))]


def replace_block(text, start, end, new):
    old = between(text, start, end)
    return text.replace(old, new, 1)


def main():
    DEST.mkdir(exist_ok=True)
    cipa = (ROOT / "cipa/_ord_equations.py").read_text()
    gong = (ROOT / "beta-signaling/gong_beta.py").read_text()
    tree = ast.parse(gong)
    functions = {n.name: ast.get_source_segment(gong, n) for n in tree.body if isinstance(n, ast.FunctionDef)}
    signal_names = re.findall(r'^    (\w+) = y\[\((\d+)\)-1\]', functions['signaling'], re.M)
    signal_names = tuple(name for name, _ in sorted(signal_names, key=lambda item: int(item[1])))
    assert len(signal_names) == 57 and len(set(signal_names)) == 57
    init = next(n.value for n in ast.walk(ast.parse(functions['initial_state']))
                if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'initial_Signaling' for t in n.targets))
    initial = ast.unparse(init)
    effective = re.sub(r'y\[\((\d+)\+58\)-1\]', lambda m: f'y[{int(m[1])-1}]', functions['effective'])
    signaling = '\n\n'.join(functions[n] for n in ('ifthenelse', 'constants', 'signaling'))
    header = '''"""Gong 2020 signaling, extracted from the partner-supplied Python translation.
See NOTICE and provenance.json. Derivatives already use milliseconds.
The source translation has not established MATLAB equation parity.
"""
import numpy as np
from numpy import exp, log, sqrt, sin, cos, arctan as atan
from types import SimpleNamespace
'''
    (DEST / 'signaling.py').write_text(header + f'\nSTATE_NAMES = {signal_names!r}\n' + signaling + '\n\n' + effective +
        '\n\ndef initial_state():\n    return ' + initial + '\n', encoding='utf-8')

    # Transplant only regulation blocks into the CiPA kernel, retaining its
    # common ion transport, geometry, physical constants, and Markov hERG.
    electrical = functions['ephys']
    names = dict(V='v', Nai='nai', Nass='nass', Ki='ki', Kss='kss', Cai='cai',
                 Cass='cass', Cansr='cansr', Cajsr='cajsr', Nao='nao', Cao='cao', Ko='ko')
    derivatives = {'dJrelnp':36, 'dJrelp':37, **{'d'+n:50+i for i,n in enumerate(EXTRA)}}
    def mapped(block):
        for name, value in derivatives.items():
            block = re.sub(r'\b'+name+r'\b', f'ydot[{value}]', block)
        for name, value in names.items():
            block = re.sub(r'\b'+name+r'\b', value, block)
        return block
    kernel = cipa.replace('from .acceleration', 'from cipa.acceleration')
    kernel = kernel.replace('def evaluate(t, y, p):', 'def evaluate(t, y, p, control, pp1, enabled):')
    kernel = kernel.replace('ydot = np.zeros(50)', 'ydot = np.zeros(67)')
    assignments = ''.join(f'    {name} = y[{50+i}]\n' for i,name in enumerate(EXTRA))
    assignments += '''    ICaL_P, IKs_P, SERCA_P, TnI_P, INa_P, INaK_P, RyR_P, IKb_P = control
    GKS_mul = 8.
    GCal_mul = 1.9
    VP_d = 12.
    VP_f = 8.
    RyR_total = 1.
    SERCA_total = 1.
    Leak_total = 1.
'''
    kernel = kernel.replace('    celltype = p[0]', assignments + '    celltype = p[0]')
    line = '    ydot[38]=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*CaMKt'
    assert line in kernel
    kernel = kernel.replace(line, line + '\n    if enabled:\n        ydot[38]=aCaMK*CaMKb*(CaMKb+CaMKt)-bCaMK*(0.1+0.9*pp1/0.1371)*CaMKt')
    kernel = replace_block(kernel, '    fINap=', '    mLss=', mapped(between(electrical, '    GNaP =', '    mLss=')))
    ca = between(electrical, '    PCaP =', '    xrss=')
    # Preserve CiPA GHK fluxes. Gong's extra 0.025 mM Cass ceiling is not
    # needed for coupling, and would change the disabled baseline.
    ca = replace_block(ca, '    Cass_ceiling =', '    ICaL_NP=', '')
    kernel = replace_block(kernel, '    fICaLp=', '    GKr=', mapped(ca))
    kernel = replace_block(kernel, '    IKs=GKs*', '    xk1ss=', mapped(between(electrical, '    GKsP =', '    xk1ss=')))
    kernel = kernel.replace('    Knai0=9.073', '    Knai0=9.073*(1.-0.3*INaK_P)')
    kernel = replace_block(kernel, '    IKb=GKb*', '    PNab=', mapped(between(electrical, '    IKb_NP=', '    INab=')))
    kernel = replace_block(kernel, '    bt=4.75', '    Jtr=', mapped(between(electrical, '    bt=4.75', '    Jtr=')))
    kernel = kernel.replace('    kmtrpn=0.0005', '    kmtrpn=0.0005*(1.+0.6*TnI_P)')
    # Validate all 17 new state derivatives were assigned.
    for i in range(50, 67):
        assert re.search(rf'ydot\[{i}\]\s*=', kernel), i
    kernel = kernel.replace('Generated Python equations derived from FDA/CiPA; see licenses/.',
        'Generated CiPA + Gong hybrid. See NOTICE, README.md and tools/build_cipa_beta.py.')
    ast.parse(kernel)
    (DEST / '_equations.py').write_text(kernel, encoding='utf-8')
    provenance = dict(extra_states=EXTRA, signaling_states=57, total_states=124,
        sources={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in
                 ['cipa/_ord_equations.py', 'beta-signaling/gong_beta.py']},
        modifications=['CiPA baseline parameters and Markov IKr preserved',
                       'Gong phosphorylation and PP1-CaMK coupling',
                       'Gong Cass ceiling excluded to preserve CiPA baseline',
                       'Effective mapping uses signaling-local indices'])
    (DEST/'provenance.json').write_text(json.dumps(provenance, indent=2))


if __name__ == '__main__':
    main()
