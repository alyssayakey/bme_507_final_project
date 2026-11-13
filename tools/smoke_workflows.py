"""Small end-to-end checks. These deliberately short fits are not research results."""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cipa.cli import main
from cipa.analysis import tdp_errors


def run(*args):
    main(list(args))


if __name__ == "__main__":
    out = ROOT / "results/smoke"
    h, hill, ap = (str(out / n) for n in ("herg", "hill", "ap"))
    run("bootstrap", "-d", "quinidine", "-n", "2", "--output", h)
    run("herg-fit", "-d", "quinidine", "--population", "6", "--maxiter", "1", "--output", h)
    run("herg-fit", "-d", "quinidine", "-i", "1-2", "--population", "6", "--maxiter", "1", "--output", h)
    run("herg-summary", "-d", "quinidine", "--output", h)
    run("herg-sensitivity", "-d", "quinidine", "--output", h)
    run("hill-fit", "-d", "quinidine", "-n", "50", "-b", "100", "-t", "2", "--output", hill)
    run("hill-summary", "-d", "quinidine", "--output", hill)
    run("filter-li2017", "--output", hill)
    for samples in ("0", "1-2"):
        run("ap", "-d", "quinidine", "-i", samples, "-x", "1", "--beats", "2", "--nsave", "1",
            "--hergpath", h, "--hillpath", hill, "--output", ap)
    run("ap", "-d", "control", "--beats", "2", "--nsave", "1", "--output", ap)
    run("combine", str(out / "ap/fixed"))
    run("combine", str(out / "ap/uncertainty"), "-n", "2")
    run("qnet-ci", "--results", ap)
    # Classification smoke uses labelled synthetic observations in all categories.
    rng = np.random.default_rng(5)
    risk = pd.read_csv(ROOT / "cipa/data/ap/CiPA_training_drugs.csv")
    rows = []
    for row in risk.itertuples():
        for sample in range(1, 9):
            rows.append(dict(drug=row.drug, dose=1, sample=sample, max_dv=-.1,
                             qNet=float(100-row.CiPA*15+rng.normal(0, 20))))
    errors, predictions = tdp_errors(pd.DataFrame(rows))
    assert {"training", "LOOCV"} == set(errors["mode"])
    assert predictions[["prob_0", "prob_1", "prob_2"]].sum(axis=1).between(.999999, 1.000001).all()
    report = dict(workflows="passed", drug="quinidine", note="Tiny smoke fits only; not converged research estimates",
                  risk_fits=len(errors), classification_observations=len(predictions))
    (ROOT / "validation_workflows.json").write_text(json.dumps(report, indent=2))
    print(report)
