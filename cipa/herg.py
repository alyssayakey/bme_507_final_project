"""Milnes voltage-clamp model, stratified bootstrap, and CMA-ES fitting."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.linalg import expm
from .common import DATA, read_pairs, write_pairs, intervals
from .models import herg, parameter_vector

PARAMETERS = ("Kmax", "Ku", "n", "halfmax", "Vhalf")
BOUNDS = pd.read_csv(DATA / "herg/hergmod_drug_param_bounds.txt", sep=r"\s+")
LOW, HIGH = (BOUNDS[n].to_numpy() for n in ("Low", "High"))


def encode(values):
    return 10*np.log(np.asarray(values)/LOW)/np.log(HIGH/LOW)


def decode(values):
    return LOW*(HIGH/LOW)**(np.asarray(values)/10)


def rate_matrix(parameters, concentration, voltage):
    p = parameter_vector(parameters, "herg")
    matrix = np.empty((9, 9))
    for i in range(9):
        state = np.zeros(11)
        state[i], state[9], state[10] = 1, concentration, voltage
        matrix[:, i] = herg.evaluate(0., state, p)[0][:9]
    return matrix


class MilnesModel:
    """Ten sweeps: -80 mV to 940 ms, 0 mV for 10 s, then -80 to 25 s.

    The first 900 ms and next 40 ms are both -80 mV, as in stepprotocol.R.
    Matrix exponentials solve the linear nine-state subsystem exactly.
    """
    def __init__(self, nbeats=10):
        self.nbeats = nbeats
        self.parameters = read_pairs(DATA / "herg/hergmod_pars.txt")
        state = np.zeros(9)
        state[2] = 1
        self.control_state = expm(rate_matrix(self.parameters, 0, -85)*180000) @ state

    def sweeps(self, binding, concentration, times):
        times = np.asarray(times, dtype=float)
        if not np.isfinite(times).all() or np.any(np.diff(times) <= 0) or times[0] < 0 or times[-1] > 25000:
            raise ValueError("Milnes times must be unique, increasing, and within 0..25000 ms")
        p = {**self.parameters, **binding}
        matrices = {v: rate_matrix(p, concentration, v) for v in (-80, 0)}
        state = self.control_state.copy()
        if concentration > 0:
            state = expm(matrices[-80]*180000) @ state
        knots = np.unique(np.r_[0, times, 940, 10940, 25000])
        transitions = []
        cache = {}
        for a, b in zip(knots[:-1], knots[1:]):
            key = (0 if 940 <= a < 10940 else -80, round(b-a, 9))
            if key not in cache:
                cache[key] = expm(matrices[key[0]]*key[1])
            transitions.append(cache[key])
        out = []
        indices = np.searchsorted(knots, times)
        for _ in range(self.nbeats):
            states = [state.copy()]
            for transition in transitions:
                state = transition @ state
                states.append(state.copy())
            out.append(np.asarray(states)[indices])
        return np.asarray(out)

    def fractional_current(self, binding, concentration, times):
        control = self.sweeps({}, 0, times)[:, :, 4]
        drug = self.sweeps(binding, concentration, times)[:, :, 4]
        if np.any(control <= 0):
            raise ValueError("Fractional current requested where control open probability is zero")
        return drug/control


def experiments(data):
    return data[["conc", "exp"]].drop_duplicates().sort_values(["conc", "exp"]).reset_index(drop=True)


def bootstrap(data, nboots=2000, seed=100):
    if nboots < 1:
        raise ValueError("nboots must be positive")
    cells = experiments(data)
    rng = np.random.default_rng(seed)
    indices = np.empty((nboots, len(cells)), dtype=int)
    for _, group in cells.groupby("conc", sort=True):
        idx = group.index.to_numpy()
        indices[:, idx] = rng.choice(idx, size=(nboots, len(idx)), replace=True)
    return indices


def mean_data(data, indices=None):
    cells = experiments(data)
    if indices is not None:
        cells = cells.iloc[np.asarray(indices)]
    sampled = pd.concat([data[(data.conc == row.conc) & (data.exp == row.exp)]
                         for row in cells.itertuples()], ignore_index=True)
    mean = sampled.groupby(["conc", "sweep", "time"], as_index=False).frac.mean()
    selected = []
    for _, group in mean.groupby("conc", sort=True):
        first = group[group.sweep == 1].sort_values("time")
        crossings = np.flatnonzero(first.frac.to_numpy() <= .5)
        start = min(10, int(crossings[0])) if len(crossings) else 10
        for _, sweep in group.groupby("sweep", sort=True):
            selected.append(sweep.sort_values("time").iloc[start:])
    return pd.concat(selected, ignore_index=True)


class HergObjective:
    def __init__(self, data):
        self.data = data
        self.model = MilnesModel()
        self.groups = []
        for conc, group in data.groupby("conc", sort=True):
            times = np.sort(group.time.unique())
            control = self.model.sweeps({}, 0, times)[:, :, 4]
            self.groups.append((conc, group, times, control))

    def __call__(self, encoded):
        try:
            binding = dict(zip(PARAMETERS, decode(encoded)))
            error, special, negative = 0., 0., 0.
            for conc, group, times, control in self.groups:
                states = self.model.sweeps(binding, conc, times)
                prediction = states[:, :, 4] / control
                first_obs, first_pred = {}, {}
                significant = (group.frac <= .5).any()
                for sweep, obs in group.groupby("sweep", sort=True):
                    idx = np.searchsorted(times, obs.time)
                    pred = prediction[int(sweep)-1, idx]
                    error += np.sum((pred-obs.frac.to_numpy())**2)
                    first_obs[sweep], first_pred[sweep] = obs.frac.iloc[0], pred[0]
                negative += np.mean(np.minimum(states[-1], 0)**2, axis=0).sum()
                if significant:
                    for last in (2, 10):
                        if last in first_obs and first_obs[1] != 0 and first_pred[1] != 0:
                            delta_obs = max(0., (first_obs[1]-first_obs[last])/first_obs[1])
                            delta_pred = (first_pred[1]-first_pred[last])/first_pred[1]
                            special += .2*len(self.data)/len(self.groups)*(delta_obs-delta_pred)**2
            value = error + negative + special
            return float(value) if np.isfinite(value) else 1e50
        except (ValueError, FloatingPointError, OverflowError):
            return 1e50


def fit(data, output, indices=None, seed=100, population=80, maxiter=2500, tolerance=.001, start=None):
    import cma
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    objective = HergObjective(mean_data(data, indices))
    initial = encode([start[n] for n in PARAMETERS]) if start is not None else np.random.default_rng(seed).uniform(0, 10, 5)
    es = cma.CMAEvolutionStrategy(initial, 2, dict(bounds=[0, 10], popsize=population,
                                  maxiter=maxiter, tolx=tolerance, seed=seed, verbose=-9))
    history = []
    while not es.stop():
        candidates = es.ask()
        errors = [objective(x) for x in candidates]
        es.tell(candidates, errors)
        history.append(dict(generation=es.countiter, error=min(errors),
                            **dict(zip(PARAMETERS, decode(candidates[np.argmin(errors)])))))
    if es.result.fbest >= 1e49:
        raise RuntimeError("All hERG simulations failed during fitting")
    values = dict(zip(PARAMETERS, decode(es.result.xbest)))
    write_pairs(values, output / "pars.txt")
    pd.DataFrame(history).to_csv(output / "CMAES_history.csv", index=False)
    (output / "fit_diagnostics.json").write_text(json.dumps(dict(objective=es.result.fbest,
        evaluations=es.result.evaluations, stop={k: str(v) for k, v in es.stop().items()}), indent=2))
    return values


def summarize(output, alpha=.05):
    output = Path(output)
    nboots = len(np.load(output / "bootstrap_indices.npy"))
    optimal = read_pairs(output / "pars.txt")
    rows = []
    for i in range(1, nboots+1):
        path = output / f"boot/{i:05d}/pars.txt"
        rows.append(read_pairs(path) if path.exists() else {n: np.nan for n in PARAMETERS})
    samples = pd.DataFrame(rows)
    samples["slope"] = samples.Kmax/samples.halfmax
    samples.to_csv(output / "boot_pars.csv", index=False)
    optimal["slope"] = optimal["Kmax"]/optimal["halfmax"]
    transformed, opt = samples.copy(), {}
    for name in samples:
        log = name not in ("n", "Vhalf")
        label = "log10" + name if log else name
        transformed = transformed.rename(columns={name: label})
        if log:
            transformed[label] = np.log10(transformed[label])
        opt[label] = np.log10(optimal[name]) if log else optimal[name]
    ci = intervals(transformed, opt, alpha)
    ci.to_csv(output / f"boot_CI{100*(1-alpha):g}.csv", index=False)
    return samples, ci


def sensitivity(data, samples, alpha=.05):
    model = MilnesModel()
    rows = []
    for conc, group in mean_data(data).groupby("conc", sort=True):
        times = np.sort(group[group.sweep == 1].time.unique())
        traces = np.array([model.fractional_current(row, conc, times)
                           for row in samples[list(PARAMETERS)].dropna().to_dict("records")])
        if not len(traces):
            raise ValueError("No complete bootstrap parameter samples")
        lo, hi = np.quantile(traces, [alpha/2, 1-alpha/2], axis=0)
        mean = traces.mean(axis=0)
        for sweep in range(model.nbeats):
            rows.extend(dict(conc=conc, sweep=sweep+1, time=t, mean=m, lower=l, upper=h)
                        for t, m, l, h in zip(times, mean[sweep], lo[sweep], hi[sweep]))
    return pd.DataFrame(rows)
