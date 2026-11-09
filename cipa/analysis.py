"""Result aggregation, qNet intervals, ordinal risk classification, and plots."""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from .common import DATA


def combine(directory, nsamples=None):
    directory = Path(directory)
    files = sorted(p for p in directory.rglob("metrics.csv") if p != directory / "metrics.csv")
    # Avoid folding current-omission runs into the main analysis.
    files = [p for p in files if not any(s.startswith("drop_") for s in p.relative_to(directory).parts)]
    if not files:
        raise ValueError(f"No simulation metrics under {directory}")
    result = pd.concat([pd.read_csv(p) for p in files], ignore_index=True)
    keys = ["drug", "dose", "sample"]
    if result.duplicated(keys).any():
        raise ValueError("Duplicate drug/dose/sample simulations")
    if nsamples is not None:
        expected = set(range(1, nsamples+1))
        for (drug, dose), group in result.groupby(["drug", "dose"]):
            if drug != "control" and set(group["sample"]) != expected:
                raise ValueError(f"Missing or unexpected samples for {drug} dose {dose}")
    result = result.sort_values(keys)
    result.to_csv(directory / "metrics.csv", index=False)
    return result


def qnet_intervals(fixed, uncertainty, alpha=.05):
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    rows = []
    for row in fixed[fixed.drug != "control"].itertuples():
        group = uncertainty[(uncertainty.drug == row.drug) & (uncertainty.dose == row.dose)].sort_values("sample")
        if group.empty or group["sample"].tolist() != list(range(1, len(group)+1)):
            raise ValueError(f"Missing or duplicate uncertainty samples for {row.drug}, dose {row.dose}")
        valid = group.loc[group.max_dv.notna(), "qNet"]
        lo, hi = valid.quantile([alpha/2, 1-alpha/2]) if len(valid) else (np.nan, np.nan)
        rows.append(dict(drug=row.drug, dose=row.dose, value=row.qNet if np.isfinite(row.max_dv) else np.nan,
                         samples=len(group), conf=1-alpha, lower=lo, upper=hi))
    return pd.DataFrame(rows)


def tdp_errors(metrics, riskfile=DATA / "ap/CiPA_training_drugs.csv"):
    from statsmodels.miscmodels.ordinal_model import OrderedModel
    risk = pd.read_csv(riskfile)
    df = metrics[(metrics.drug != "control") & (metrics.dose != 0) & metrics.max_dv.notna()].merge(risk[["drug", "CiPA"]], on="drug", validate="many_to_one")
    errors, predictions = [], []
    for dose, group in df.groupby("dose", sort=True):
        for heldout in [None, *sorted(group.drug.unique())]:
            train = group if heldout is None else group[group.drug != heldout]
            test = group if heldout is None else group[group.drug == heldout]
            if train.CiPA.nunique() < 3:
                warnings.warn(f"Dose {dose}, heldout {heldout}: training requires all three risk categories", stacklevel=2)
                continue
            # Scaling improves conditioning; it does not alter predictions.
            center, scale = train.qNet.mean(), train.qNet.std()
            if not np.isfinite(scale) or scale == 0:
                raise ValueError("qNet has no variation")
            x = pd.DataFrame({"qNet": (train.qNet-center)/scale}, index=train.index)
            model = OrderedModel(train.CiPA.astype(int), x, distr="logit")
            fit = model.fit(method="bfgs", maxiter=10000, disp=False)
            if not fit.mle_retvals.get("converged", False):
                warnings.warn(f"Ordinal fit did not converge: dose {dose}, heldout {heldout}", stacklevel=2)
            probabilities = np.asarray(model.predict(fit.params, exog=pd.DataFrame({"qNet": (test.qNet-center)/scale})))
            predicted = np.argmax(probabilities, axis=1)
            cutpoints = model.transform_threshold_params(fit.params)[1:-1]
            slope = float(fit.params.iloc[0])/scale
            intercepts = -cutpoints - float(fit.params.iloc[0])*center/scale
            errors.append(dict(dose=dose, drug=heldout, mode="training" if heldout is None else "LOOCV",
                               intercept1=intercepts[0], intercept2=intercepts[1], slope=slope,
                               error=float(np.mean(np.abs(predicted-test.CiPA.to_numpy())))))
            for (_, row), pred, probs in zip(test.iterrows(), predicted, probabilities):
                predictions.append(dict(drug=row.drug, dose=dose, sample=row["sample"],
                    mode="training" if heldout is None else "LOOCV", CiPA=int(row.CiPA), predicted=int(pred),
                    **{f"prob_{i}": value for i, value in enumerate(probs)}))
    if not errors:
        raise ValueError("No risk models could be fitted; check categories, doses, and valid AP metrics")
    return pd.DataFrame(errors), pd.DataFrame(predictions)


def plot_intervals(df, path, x="dose", group="drug", ylabel="qNet (microC/microF)", factor=.001):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for name, data in df.groupby(group, sort=True):
        data = data.sort_values(x)
        ax.plot(data[x], data.value*factor, label=str(name))
        ax.fill_between(data[x], data.lower*factor, data.upper*factor, alpha=.2)
    ax.set(xlabel=x, ylabel=ylabel)
    ax.legend(fontsize=7, ncol=2)
    fig.savefig(path)
    plt.close(fig)


def plot_trace(trace, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # print(trace.columns.tolist())
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, layout="constrained")
    axes[0].plot(trace.time, trace.v)
    axes[0].set_ylabel("Voltage (mV)")
    axes[1].plot(trace.time, trace.cai*1000)
    axes[1].set(xlabel="Time (ms)", ylabel="Calcium (microM)")
    fig.savefig(path)
    plt.close(fig)
