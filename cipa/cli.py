"""Cross-platform replacements for the R scripts and shell launchers."""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .common import DATA, DRUGS, parse_range, read_pairs


def parser():
    p = argparse.ArgumentParser(prog="cipa", description="Python port of FDA/CiPA")
    sub = p.add_subparsers(dest="command", required=True)
    b = sub.add_parser("bootstrap", help="Generate stratified hERG cell bootstrap indices")
    b.add_argument("-d", "--drug", default=DRUGS)
    b.add_argument("-n", "--nboots", type=int, default=2000)
    b.add_argument("-s", "--seed", type=int, default=100)
    b.add_argument("--data", type=Path, default=DATA / "herg")
    b.add_argument("--output", type=Path, default=Path("results/herg"))
    h = sub.add_parser("herg-fit", help="Fit hERG drug binding with CMA-ES")
    h.add_argument("-d", "--drug", required=True)
    h.add_argument("-i", "--sample", default="0")
    h.add_argument("-s", "--seed", type=int, default=100)
    h.add_argument("-l", "--population", type=int, default=80)
    h.add_argument("-m", "--maxiter", type=int, default=2500)
    h.add_argument("-t", "--tol", type=float, default=.001)
    h.add_argument("--data", type=Path, default=DATA / "herg")
    h.add_argument("--output", type=Path, default=Path("results/herg"))
    for name, helptext, default in [("herg-summary", "Combine hERG bootstrap fits and compute intervals", "herg"),
                                    ("herg-sensitivity", "Propagate hERG bootstrap variation through Milnes protocol", "herg"),
                                    ("hill-summary", "Compute Hill MCMC credible intervals", "hill")]:
        s = sub.add_parser(name, help=helptext)
        s.add_argument("-d", "--drug", default=DRUGS)
        s.add_argument("-a", "--alpha", type=float, default=.05)
        s.add_argument("--output", type=Path, default=Path("results") / default)
        if name == "herg-sensitivity":
            s.add_argument("--data", type=Path, default=DATA / "herg")
    h = sub.add_parser("hill-fit", help="Fit Hill dose-response curves and sample uncertainty")
    h.add_argument("-d", "--drug", required=True)
    h.add_argument("-f", "--filepath", type=Path, default=DATA / "hill/drug_block.csv")
    h.add_argument("-n", "--nsamples", type=int, default=2000)
    h.add_argument("-b", "--burnin", type=int, default=10000)
    h.add_argument("-t", "--thin", type=int, default=10)
    h.add_argument("-s", "--seed", type=int, default=100)
    h.add_argument("--output", type=Path, default=Path("results/hill"))
    s = sub.add_parser("process-patchclamp", help="Prepare patch clamp input CSV (concentrations in nM)")
    s.add_argument("--data", type=Path, default=DATA / "hill")
    s.add_argument("--output", type=Path, default=Path("results/hill/drug_block.csv"))
    s = sub.add_parser("filter-li2017", help="Mask Hill values absent from Li 2017; creates backups")
    s.add_argument("--output", type=Path, default=Path("results/hill"))
    s.add_argument("--table", type=Path, default=DATA / "hill/Li2017_IC50.csv")
    a = sub.add_parser("ap", help="Run action potentials with fixed or uncertainty inputs")
    a.add_argument("-d", "--drug", default="control")
    a.add_argument("-i", "--sample", default="0")
    a.add_argument("-x", "--dose", default="0")
    a.add_argument("--beats", type=int, default=1000)
    a.add_argument("--nsave", type=int, default=250)
    a.add_argument("--trace-beats", type=int, default=0,
                   help="Additionally save the final N consecutive beats as trace_cycles.csv/png [default 0]")
    a.add_argument("--cycle-length", type=float, default=2000)
    a.add_argument("--dt", type=float, default=1)
    a.add_argument("-c", "--cmaxfile", type=Path, default=DATA / "ap/CiPA_training_drugs.csv")
    a.add_argument("-b", "--hergpath", type=Path, default=Path("results/herg"))
    a.add_argument("-m", "--hillpath", type=Path, default=Path("results/hill"))
    a.add_argument("-r", "--dropcurrent", choices=["hERG", "ICaL", "INaL", "INa", "IKs", "Ito", "IK1"])
    a.add_argument("-o", "--omit", action="store_true")
    a.add_argument("--no-opt-scaling", action="store_true")
    a.add_argument("--output", type=Path, default=Path("results/ap"))
    s = sub.add_parser("combine", help="Combine AP metrics; optionally verify all sample indices")
    s.add_argument("directory", type=Path)
    s.add_argument("-n", "--nsamples", type=int)
    s = sub.add_parser("qnet-ci", help="Compute and plot qNet credible intervals")
    s.add_argument("--results", type=Path, default=Path("results/ap"))
    s.add_argument("-a", "--alpha", type=float, default=.05)
    s = sub.add_parser("tdp-error", help="Ordinal logistic classification and leave-one-drug-out validation")
    s.add_argument("directory", type=Path)
    s.add_argument("-t", "--tdpfile", type=Path, default=DATA / "ap/CiPA_training_drugs.csv")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    from . import herg, hill, analysis
    command = args.command
    if command == "bootstrap":
        for index, drug in enumerate(args.drug.replace(" ", "").split(",")):
            data = pd.read_csv(args.data / f"{drug}.csv")
            out = args.output / drug
            out.mkdir(parents=True, exist_ok=True)
            np.save(out / "bootstrap_indices.npy", herg.bootstrap(data, args.nboots, args.seed+index))
            herg.experiments(data).to_csv(out / "experiments.csv", index=False)
            print(f"{drug}: saved {args.nboots} bootstrap samples")
    elif command == "herg-fit":
        data = pd.read_csv(args.data / f"{args.drug}.csv")
        out = args.output / args.drug
        for sample in parse_range(args.sample, integer=True):
            indices, start, dest = None, None, out
            if sample:
                indices = np.load(out / "bootstrap_indices.npy")[sample-1]
                cells = pd.read_csv(out / "experiments.csv")
                if not cells.equals(herg.experiments(data)):
                    raise ValueError("Bootstrap data changed; regenerate indices")
                start = read_pairs(out / "pars.txt")
                dest = out / f"boot/{sample:05d}"
            values = herg.fit(data, dest, indices, args.seed, args.population, args.maxiter, args.tol, start)
            print(args.drug, sample, values)
    elif command in ("herg-summary", "hill-summary", "herg-sensitivity"):
        for drug in args.drug.replace(" ", "").split(","):
            out = args.output / drug
            if command == "herg-summary":
                herg.summarize(out, args.alpha)
            elif command == "hill-summary":
                hill.summarize(out, args.alpha)
            else:
                df = herg.sensitivity(pd.read_csv(args.data / f"{drug}.csv"), pd.read_csv(out / "boot_pars.csv"), args.alpha)
                df.to_csv(out / "sensRange_summary.csv", index=False)
            print(out)
    elif command == "hill-fit":
        hill.fit_drug(args.drug, args.filepath, args.output, args.nsamples, args.burnin, args.thin, args.seed)
        print(args.output / args.drug)
    elif command == "process-patchclamp":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        hill.process_patchclamp(args.data).to_csv(args.output, index=False)
        print(args.output)
    elif command == "filter-li2017":
        hill.filter_li2017(args.output, args.table)
    elif command == "ap":
        from .ap import simulate
        doses = [0.] if args.drug == "control" else parse_range(args.dose)
        samples = parse_range(args.sample, integer=True)
        cmax = 1.
        if args.drug != "control":
            rows = pd.read_csv(args.cmaxfile) if args.cmaxfile.exists() else pd.DataFrame(columns=["drug", "therapeutic"])
            values = rows.loc[rows.drug == args.drug, "therapeutic"]
            if len(values) > 1:
                raise ValueError(f"Duplicate therapeutic concentrations for {args.drug}")
            if len(values):
                cmax = float(values.iloc[0])
                if not np.isfinite(cmax) or cmax <= 0:
                    raise ValueError("Cmax must be finite and positive")
            else:
                print("Cmax unavailable; doses interpreted as nM")
        drops = ["hERG", "ICaL", "INaL", "INa", "IKs", "Ito", "IK1"] if args.omit else [args.dropcurrent]
        for sample in samples:
            for drop in drops:
                out = args.output / ("uncertainty" if sample else "fixed")
                if drop:
                    out /= f"drop_{drop}"
                out /= args.drug
                if sample:
                    out /= f"{sample:05d}"
                out.mkdir(parents=True, exist_ok=True)
                rows = []
                for dose in doses:
                    metric, trace = simulate(args.drug, dose*cmax, sample, args.beats, args.nsave,
                        args.cycle_length, args.dt, args.hergpath, args.hillpath, drop,
                        not args.no_opt_scaling, out / f"dose_{dose:g}", trace_beats=args.trace_beats)
                    row = dict(drug=args.drug, dose=dose, sample=sample, CL=args.cycle_length,
                               concentration=dose*cmax, **metric)
                    rows.append(row)
                    pd.DataFrame(rows).to_csv(out / "metrics.csv", index=False)
                    analysis.plot_trace(trace, out / f"dose_{dose:g}/trace.png")
                    if args.trace_beats:
                        cycles = pd.read_csv(out / f"dose_{dose:g}/trace_cycles.csv")
                        analysis.plot_trace(cycles, out / f"dose_{dose:g}/trace_cycles.png")
                    print(json.dumps(row))
    elif command == "combine":
        print(analysis.combine(args.directory, args.nsamples))
    elif command == "qnet-ci":
        df = analysis.qnet_intervals(pd.read_csv(args.results / "fixed/metrics.csv"),
                                    pd.read_csv(args.results / "uncertainty/metrics.csv"), args.alpha)
        name = args.results / f"qNet_CI{100*(1-args.alpha):g}"
        df.to_csv(name.with_suffix(".csv"), index=False)
        analysis.plot_intervals(df, name.with_suffix(".png"))
        print(name)
    elif command == "tdp-error":
        errors, predictions = analysis.tdp_errors(pd.read_csv(args.directory / "metrics.csv"), args.tdpfile)
        for mode in ("training", "LOOCV"):
            errors[errors["mode"] == mode].to_csv(args.directory / f"{mode}_errors.csv", index=False)
            predictions[predictions["mode"] == mode].to_csv(args.directory / f"{mode}_probs.csv", index=False)
        print(errors)
