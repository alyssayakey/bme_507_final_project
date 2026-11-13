"""Short four-condition integration checks, not steady-state validation."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cipa_beta.cli import main


if __name__ == '__main__':
    out = ROOT/'results/cipa_beta_validation'
    rows = []
    for iso in (0, 1):
        for drug in ('control', 'quinidine'):
            main(['ap', '-d', drug, '--iso', str(iso), '--beats', '5', '--nsave', '2',
                  '--trace-beats', '5', '--output', str(out)])
            folder = out/f'iso_{iso}_uM/fixed'/drug
            metric = pd.read_csv(folder/'metrics.csv').iloc[0].to_dict()
            trace = pd.read_csv(folder/f'dose_{0 if drug == "control" else 1}/trace_cycles.csv')
            assert np.isfinite(trace.to_numpy()).all()
            assert (trace[['nai','nass','ki','kss','cai','cass','cansr','cajsr']] > 0).all().all()
            assert trace.beat.nunique() == 5
            assert np.max(np.abs(trace[['IC1','IC2','C1','C2','O','IO','IObound','Obound','Cbound']].sum(axis=1)-.9999999779145)) < 1e-5
            rows.append(metric)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, layout='constrained')
    for row in rows:
        drug, iso = row['drug'], row['iso_uM']
        trace = pd.read_csv(out/f'iso_{iso:g}_uM/fixed'/drug/f'dose_{0 if drug == "control" else 1}/trace.csv')
        label = f'{drug}, ISO={iso:g} uM'
        axes[0].plot(trace.time, trace.v, label=label)
        axes[1].plot(trace.time, trace.cai*1000, label=label)
    axes[0].set(ylabel='Voltage (mV)', title='Five-beat operational checks; selected traces')
    axes[1].set(xlabel='Time (ms)', ylabel='Calcium (microM)', xlim=(0, 1000))
    axes[0].legend(fontsize=8)
    fig.savefig(out/'comparison.png')
    plt.close(fig)
    report = dict(status='passed', conditions=rows, beats_per_condition=5,
                  limitation='Short operational checks, not MATLAB parity or steady-state validation')
    (ROOT/'cipa_beta/validation.json').write_text(json.dumps(report, indent=2))
    print('Four-condition validation passed:', out)
