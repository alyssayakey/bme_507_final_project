"""Bounded Hill fitting and adaptive Metropolis uncertainty sampling."""
from pathlib import Path
import json
import shutil
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import expit
from .common import DATA, intervals

LOW = np.array([-10., 0.])
HIGH = np.array([10., 10.])


def hill_block(concentration, log10_ic50, h):
    concentration = np.asarray(concentration, dtype=float)
    if np.any(concentration < 0) or np.any(np.asarray(h) < 0):
        raise ValueError("Concentrations and Hill coefficient must be nonnegative")
    with np.errstate(divide="ignore", invalid="ignore"):
        z = h * (np.log(concentration)-np.log(10.)*log10_ic50)
    # R's 0^0 is 1, giving 50% block at h=0.
    z = np.where(h == 0, 0., z)
    return 100*expit(z)


def fit_hill(concentration, block):
    x, y = np.asarray(concentration, float), np.asarray(block, float)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0)
    x, y = x[valid], y[valid]
    if len(np.unique(x)) < 4:
        raise ValueError("At least four distinct positive concentrations are required")
    starts = np.r_[np.mean(x)*np.array([1, .5, .3, .75]),
                   [1e-5, 1e-4, .001, .01, .1, 1, 100, 300, 500, 1000, 3000, 5000,
                    10000, 30000, 50000, 100000, 300000, 500000, 1e6, 3e6, 5e6]]
    fits = [least_squares(lambda p: hill_block(x, *p)-y,
                         np.clip([np.log10(s), .9], LOW+1e-8, HIGH-1e-8),
                         bounds=(LOW, HIGH), max_nfev=2000) for s in starts]
    best = min((f for f in fits if f.success), key=lambda f: np.sum(f.fun**2), default=None)
    if best is None:
        raise RuntimeError("Hill least-squares fit failed")
    variance = max(float(np.sum(best.fun**2)/max(len(x)-2, 1)), 1e-8)
    covariance = np.linalg.pinv(best.jac.T @ best.jac)*variance
    return best.x, covariance, variance, x, y


def geweke(chain):
    """Geweke-style early/late z-score using batch-mean spectral variances."""
    def stats(values):
        n = len(values)
        batch = max(2, int(np.sqrt(n)))
        count = n//batch
        if count < 2:
            return values.mean(axis=0), np.full(2, np.nan)
        means = values[:count*batch].reshape(count, batch, 2).mean(axis=1)
        return values.mean(axis=0), means.var(axis=0, ddof=1)/count
    if len(chain) < 40:
        return np.full(2, np.nan)
    m1, s1 = stats(chain[:max(4, len(chain)//10)])
    m2, s2 = stats(chain[len(chain)//2:])
    return (m1-m2)/np.sqrt(s1+s2+1e-30)


def sample_hill(x, y, initial, covariance, variance, nsamples=2000, burnin=10000, thin=10, seed=100):
    if nsamples < 1 or burnin < 0 or thin < 1:
        raise ValueError("Require nsamples >= 1, burnin >= 0, thin >= 1")
    rng = np.random.default_rng(seed)
    current = initial.copy()
    residual = y-hill_block(x, *current)
    sse = residual @ residual
    sigma2 = variance
    proposal_cov = covariance*(2.4**2/2) + np.eye(2)*1e-6
    total = burnin + nsamples*thin
    history = np.empty((total, 2))
    accepted = 0
    # Inverse-gamma variance prior: n0=N, scale=n0*estimated variance/2.
    n0 = len(x)
    for i in range(total):
        candidate = rng.multivariate_normal(current, proposal_cov)
        if np.all(candidate >= LOW) and np.all(candidate <= HIGH):
            r = y-hill_block(x, *candidate)
            new_sse = r @ r
            if np.log(rng.uniform()) < -(new_sse-sse)/(2*sigma2):
                current, sse = candidate, new_sse
                accepted += 1
        sigma2 = (n0*variance+sse)/2/rng.gamma((n0+len(x))/2)
        history[i] = current
        # Adapt only during burn-in, leaving a fixed sampling kernel afterward.
        if i >= 100 and i < burnin and (i+1) % 10 == 0:
            proposal_cov = np.cov(history[:i+1].T)*(2.4**2/2) + np.eye(2)*1e-8
    chain = history[burnin+thin-1::thin]
    z = geweke(chain)
    return chain, dict(acceptance_rate=accepted/total, geweke_z=z.tolist(),
                       convergence_passed=bool(np.isfinite(z).all() and (np.abs(z) < 1.96).all()))


def fit_drug(drug, datafile=DATA / "hill/drug_block.csv", output="results/hill", nsamples=2000,
             burnin=10000, thin=10, seed=100):
    data = pd.read_csv(datafile)
    data = data[data.drug == drug]
    if data.empty:
        raise ValueError(f"No data for {drug}")
    output = Path(output) / drug
    output.mkdir(parents=True, exist_ok=True)
    optimal, samples, diagnostics = {}, {}, {}
    for channel, group in data.groupby("channel", sort=True):
        names = [channel+"_IC50", channel+"_h"]
        try:
            par, cov, variance, x, y = fit_hill(group.conc, group.block)
            chain, diag = sample_hill(x, y, par, cov, variance, nsamples, burnin, thin, seed)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            warnings.warn(f"{drug}/{channel}: {exc}", stacklevel=2)
            for name in names:
                optimal[name], samples[name] = np.nan, np.full(nsamples, np.nan)
            diagnostics[channel] = dict(error=str(exc))
            continue
        optimal[names[0]], optimal[names[1]] = 10**par[0], par[1]
        samples[names[0]], samples[names[1]] = 10**chain[:, 0], chain[:, 1]
        diagnostics[channel] = diag
        pd.DataFrame(chain, columns=["log10IC50", "h"]).to_csv(output / f"{channel}_chain.csv", index=False)
        if not diag["convergence_passed"]:
            warnings.warn(f"{drug}/{channel}: convergence diagnostic did not pass; inspect chain before using samples", stacklevel=2)
        grid = np.geomspace(x.min(), x.max(), 100)
        curves = hill_block(grid[None, :], chain[:, 0, None], chain[:, 1, None])
        lo, hi = np.quantile(curves, [.025, .975], axis=0)
        pd.DataFrame(dict(conc=grid, block=hill_block(grid, *par), lower=lo, upper=hi)).to_csv(output / f"{channel}_sensitivity.csv", index=False)
    pd.DataFrame([optimal]).to_csv(output / "IC50_optimal.csv", index=False)
    pd.DataFrame(samples).to_csv(output / "IC50_samples.csv", index=False)
    (output / "mcmc_diagnostics.json").write_text(json.dumps(diagnostics, indent=2))
    return pd.DataFrame([optimal]), pd.DataFrame(samples)


def summarize(output, alpha=.05):
    output = Path(output)
    optimal = pd.read_csv(output / "IC50_optimal.csv").iloc[0].to_dict()
    samples = pd.read_csv(output / "IC50_samples.csv")
    opt = {}
    for name in list(samples):
        label = name.replace("_IC50", "_log10IC50")
        if name.endswith("_IC50"):
            samples[name] = np.log10(samples[name])
            optimal[name] = np.log10(optimal[name])
        samples = samples.rename(columns={name: label})
        opt[label] = optimal[name]
    ci = intervals(samples, opt, alpha)
    ci["channel"] = ci.param.str.rsplit("_", n=1).str[0]
    ci["param"] = ci.param.str.rsplit("_", n=1).str[-1]
    ci.to_csv(output / f"MCMC_CI{100*(1-alpha):g}.csv", index=False)
    return ci


def filter_li2017(output="results/hill", table=DATA / "hill/Li2017_IC50.csv", drugs=None):
    for row in pd.read_csv(table).to_dict("records"):
        if drugs is not None and row["drug"] not in drugs:
            continue
        for filename in ("IC50_optimal.csv", "IC50_samples.csv"):
            path = Path(output) / row["drug"] / filename
            if not path.exists():
                continue
            df = pd.read_csv(path)
            backup = Path(str(path) + ".bkup.0")
            i = 0
            while backup.exists():
                i += 1
                backup = Path(str(path) + f".bkup.{i}")
            shutil.copy2(path, backup)
            for name in df:
                if name in row and pd.isna(row[name]):
                    df[name] = np.nan
            df.to_csv(path, index=False)


def process_patchclamp(directory=DATA / "hill"):
    directory = Path(directory)
    data = pd.read_csv(directory / "mergedpatchclampdata-20160514.csv")
    whiteoak = pd.read_csv(directory / "WO_hERG_block.csv")
    for df in (data, whiteoak):
        df.columns = df.columns.str.lower()
        df["drug"] = df.drug.str.lower()
    data = pd.concat([data[data.channel != "hERG"], whiteoak], ignore_index=True)
    # Upstream mistakenly used $Units after lowercasing names. Correct to units.
    micro = data.units.eq("uM")
    data.loc[micro, "conc"] *= 1000
    data.loc[micro, "units"] = "nM"
    data["channel"] = data.channel.str.replace(" ", "_").replace(
        {"Calcium": "ICaL", "Kv4.3": "Ito", "Peak_sodium": "INa", "Late_sodium": "INaL"})
    return data
