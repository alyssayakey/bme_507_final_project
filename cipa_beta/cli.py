"""AP-only entry point; fitting remains in the unchanged cipa package."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from cipa.common import parse_range
from cipa.study import DRUG_TABLE
from cipa.analysis import plot_trace
from .cell import CiPABetaCell
from .ap import simulate


def main(argv=None):
    p = argparse.ArgumentParser(prog='python -m cipa_beta')
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('ap', help='Run CiPA + Gong beta-adrenergic single-cell APs')
    a.add_argument('-d', '--drug', choices=['control', 'quinidine'], default='quinidine')
    a.add_argument('-x', '--dose', default='1', help='Multiples of Cmax; quinidine Cmax=3237 nM')
    a.add_argument('-i', '--sample', default='0')
    a.add_argument('--iso', type=float, default=0., help='Constant isoproterenol in micromolar')
    a.add_argument('--no-regulation', action='store_true', help='Disable beta regulation for CiPA parity checks')
    a.add_argument('--precondition-seconds', type=float, default=10000., help='Signaling-only equilibration before pacing')
    a.add_argument('--beats', type=int, default=1000)
    a.add_argument('--nsave', type=int, default=250)
    a.add_argument('--trace-beats', type=int, default=0)
    a.add_argument('--cycle-length', type=float, default=2000.)
    a.add_argument('--dt', type=float, default=1.)
    a.add_argument('--hergpath', type=Path, default=Path('results/herg'))
    a.add_argument('--hillpath', type=Path, default=Path('results/hill'))
    a.add_argument('--cmaxfile', type=Path, default=DRUG_TABLE)
    a.add_argument('--no-opt-scaling', action='store_true')
    a.add_argument('--dropcurrent', choices=['hERG', 'ICaL', 'INaL', 'INa', 'IKs', 'Ito', 'IK1'])
    a.add_argument('--output', type=Path, default=Path('results/cipa_beta'))
    args = p.parse_args(argv)
    if args.no_regulation and args.iso:
        p.error('--iso must be zero when --no-regulation is specified')
    if args.beats < 1 or not 1 <= args.nsave <= args.beats or not 0 <= args.trace_beats <= args.beats:
        p.error('Require beats >= nsave >= 1 and 0 <= trace-beats <= beats')
    model = CiPABetaCell(args.iso, not args.no_regulation, args.precondition_seconds)
    cmax = 0.
    if args.drug != 'control':
        table = pd.read_csv(args.cmaxfile)
        values = table.loc[table.drug == args.drug, 'therapeutic']
        if len(values) != 1 or not np.isfinite(values.iloc[0]) or values.iloc[0] <= 0:
            p.error('Cmax table must contain exactly one positive finite quinidine concentration')
        cmax = float(values.iloc[0])
    condition = 'regulation_off' if args.no_regulation else f'iso_{args.iso:g}_uM'
    for sample in parse_range(args.sample, integer=True):
        out = args.output / condition / ('uncertainty' if sample else 'fixed') / args.drug
        if sample:
            out /= f'{sample:05d}'
        rows = []
        doses = [0.] if args.drug == 'control' else parse_range(args.dose)
        for dose in doses:
            dest = out / f'dose_{dose:g}'
            print(f'{args.drug}: ISO={args.iso:g} uM, dose={dose:g}x, sample={sample}', flush=True)
            metric, trace = simulate(args.drug, cmax*dose, sample, args.beats, args.nsave,
                args.cycle_length, args.dt, args.hergpath, args.hillpath, args.dropcurrent,
                not args.no_opt_scaling, dest, args.trace_beats, model)
            row = dict(drug=args.drug, dose=dose, sample=sample, concentration=cmax*dose,
                       iso_uM=args.iso, regulation_enabled=not args.no_regulation,
                       CL=args.cycle_length, **metric)
            rows.append(row)
            pd.DataFrame(rows).to_csv(out/'metrics.csv', index=False)
            plot_trace(trace, dest/'trace.png')
            if args.trace_beats:
                plot_trace(pd.read_csv(dest/'trace_cycles.csv'), dest/'trace_cycles.png')
            print(json.dumps(row), flush=True)
