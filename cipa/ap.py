"""Action potential propagation and original CiPA beat-selection metrics."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from .models import AP_STATES, CURRENTS, initial_ap, parameter_vector, ap_parameters, ord


def simulate_beat(state, parameters, cycle_length=2000, dt=1, rtol=1e-6, atol=1e-6):
    if cycle_length <= 0.5 or dt <= 0:
        raise ValueError("cycle_length must exceed 0.5 ms and dt must be positive")
    p = parameter_vector(parameters) if isinstance(parameters, dict) else np.asarray(parameters)
    grid = np.unique(np.r_[np.arange(0, cycle_length, dt), cycle_length])
    result = np.empty((len(grid), 50))
    state = np.asarray(state, dtype=float).copy()
    # Split the stimulus discontinuity explicitly so the solver cannot skip it.
    for start, stop, stimulating in [(0, 0.5, True), (0.5, cycle_length, False)]:
        def rhs(t, y):
            return ord.evaluate(0.0 if stimulating else 1.0, y, p)[0]
        sol = solve_ivp(rhs, (start, stop), state, method="LSODA", rtol=rtol, atol=atol, dense_output=True)
        if not sol.success or not np.isfinite(sol.y).all():
            raise RuntimeError(f"AP integration failed: {sol.message}")
        mask = (grid >= start) & (grid <= stop)
        result[mask] = sol.sol(grid[mask]).T
        state = sol.y[:, -1]
    currents = np.array([ord.evaluate(t, y, p)[1] for t, y in zip(grid, result)])
    return pd.DataFrame(np.column_stack([grid, result, currents]), columns=("time",) + AP_STATES + CURRENTS)


def beat_metrics(trace):
    t, v, ca = (trace[n].to_numpy() for n in ("time", "v", "cai"))
    def rep(x, rest, peak, pct):
        idx = np.flatnonzero((x <= rest + (x[peak]-rest)*(1-pct)) & (t > t[peak]))
        return t[idx[0]] if len(idx) else np.nan
    up = np.argmax(np.diff(v)/np.diff(t))
    peak = np.argmax(v)
    rest = np.min(v[:peak]) if peak else v[0]
    out = dict(qNet=trace.qNet.iloc[-1], dvdtmax=(np.diff(v)/np.diff(t))[up],
               vmax=v[peak], vrest=rest, APD50=np.nan, APD90=np.nan, max_dv=np.nan)
    start = end = t[-1]
    if v[peak] >= 0:
        start = t[peak]
    if v[peak] >= 0 and rest < -40:
        r30 = rep(v, rest, peak, .3)
        if np.isfinite(r30):
            start = r30
            r50 = rep(v, rest, peak, .5)
            if np.isfinite(r50):
                out["APD50"] = r50-t[up]
                r90 = rep(v, rest, peak, .9)
                if np.isfinite(r90):
                    out["APD90"] = r90-t[up]
                    end = r90
    mask = (t > start) & (t < end)
    if end-start > 1 and mask.any():
        out["max_dv"] = trace.dv.to_numpy()[mask].max()
    peak = np.argmax(ca)
    rest = ca[:peak].min() if peak else ca[0]
    idx = np.flatnonzero(ca > rest + .9*(ca[peak]-rest))
    t_up = t[idx[0]] if len(idx) else np.nan
    out.update(camax=ca[peak], carest=rest,
               CaTD50=rep(ca, rest, peak, .5)-t_up, CaTD90=rep(ca, rest, peak, .9)-t_up)
    return out


def simulate(drug="control", concentration=0, sample=0, beats=1000, nsave=250,
             cycle_length=2000, dt=1, hergpath="results/herg", hillpath="results/hill",
             drop=None, optimized=True, output=None, trace_beats=0):
    if beats < 1 or not 1 <= nsave <= beats:
        raise ValueError("Require beats >= nsave >= 1")
    if not isinstance(trace_beats, (int, np.integer)) or not 0 <= trace_beats <= beats:
        raise ValueError("trace_beats must be an integer between 0 and beats")
    if drug == "control":
        concentration = 0
    p = ap_parameters(drug, concentration, sample, hergpath, hillpath, drop, optimized)
    state = initial_ap(concentration)
    metrics, starts, checkpoints = [], [], []
    cycle_traces = []
    for beat in range(1, beats+1):
        if beat > 1 and beat % 250 == 1:
            checkpoints.append(dict(beat=beat, **dict(zip(AP_STATES, state))))
        trace = simulate_beat(state, p, cycle_length, dt)
        if output is not None and trace_beats and beat > beats-trace_beats:
            # Keep the actual consecutive trajectories, independent of metric selection.
            cycle = trace.copy()
            cycle.insert(0, "beat", beat)
            cycle["beat_time"] = cycle.time
            cycle["simulation_time"] = cycle.time + (beat-1)*cycle_length
            cycle["time"] = cycle.time + (beat-(beats-trace_beats+1))*cycle_length
            # At shared boundaries keep the next beat's t=0 sample, which
            # includes the new stimulus and reset qNet, without duplicate times.
            cycle_traces.append(cycle.iloc[:-1] if beat < beats else cycle)
        if beat > beats-nsave:
            metrics.append(dict(beat=beat, **beat_metrics(trace)))
            starts.append(state.copy())
        state = trace[list(AP_STATES)].iloc[-1].to_numpy().copy()
        state[-1] = 0
    df = pd.DataFrame(metrics)
    keep = df.max_dv.idxmax() if df.max_dv.notna().any() else len(df)-1
    selected = dict(df.loc[keep])
    selected.update(EADs=int((df.max_dv > 0).sum()), noDepols=int(df.max_dv.isna().sum()),
                    noRepols=int(df.APD90.isna().sum()))
    # Save the trace from the selected beat, which need not be the last beat.
    selected_trace = simulate_beat(starts[keep], p, cycle_length, dt)
    if output is not None:
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        df.to_csv(output / "beat_metrics.csv", index=False)
        selected_trace.to_csv(output / "trace.csv", index=False)
        if cycle_traces:
            pd.concat(cycle_traces, ignore_index=True).to_csv(output / "trace_cycles.csv", index=False)
        pd.DataFrame(checkpoints).to_csv(output / "checkpoints.csv", index=False)
        pd.DataFrame([dict(beat=selected["beat"], **dict(zip(AP_STATES, starts[keep]))) ]).to_csv(output / "states_max_dv.csv", index=False)
    return selected, selected_trace
