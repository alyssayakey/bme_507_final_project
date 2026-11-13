"""Action potential propagation and original CiPA beat-selection metrics."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from cipa.models import parameter_vector, ap_parameters
from .cell import CiPABetaCell
from cipa.ap import beat_metrics


def simulate_beat(state, parameters, cycle_length=2000, dt=1, rtol=1e-6, atol=1e-9,
                  model=None, start_time=0):
    model = CiPABetaCell() if model is None else model
    if cycle_length <= model.stimulus_duration or dt <= 0:
        raise ValueError("cycle_length must exceed stimulus duration and dt must be positive")
    p = parameter_vector(parameters) if isinstance(parameters, dict) else np.asarray(parameters)
    grid = np.unique(np.r_[np.arange(0, cycle_length, dt), cycle_length])
    result = np.empty((len(grid), len(model.state_names)))
    state = np.asarray(state, dtype=float).copy()
    if state.shape != (len(model.state_names),):
        raise ValueError("State vector does not match cell model state_names")
    # Split the stimulus discontinuity explicitly so the solver cannot skip it.
    for start, stop, stimulating in [(0, model.stimulus_duration, True), (model.stimulus_duration, cycle_length, False)]:
        def rhs(t, y):
            return model.evaluate(start_time+t, y, p, stimulating)[0]
        sol = solve_ivp(rhs, (start, stop), state, method="BDF", rtol=rtol, atol=atol, dense_output=True)
        if not sol.success or not np.isfinite(sol.y).all():
            raise RuntimeError(f"AP integration failed: {sol.message}")
        mask = (grid >= start) & (grid <= stop)
        result[mask] = sol.sol(grid[mask]).T
        state = sol.y[:, -1]
    currents = np.array([model.evaluate(start_time+t, y, p, t <= model.stimulus_duration)[1]
                         for t, y in zip(grid, result)])
    return pd.DataFrame(np.column_stack([grid, result, currents]),
                        columns=("time",) + model.state_names + model.output_names)



def simulate(drug="control", concentration=0, sample=0, beats=1000, nsave=250,
             cycle_length=2000, dt=1, hergpath="results/herg", hillpath="results/hill",
             drop=None, optimized=True, output=None, trace_beats=0, model=None):
    if beats < 1 or not 1 <= nsave <= beats:
        raise ValueError("Require beats >= nsave >= 1")
    if not isinstance(trace_beats, (int, np.integer)) or not 0 <= trace_beats <= beats:
        raise ValueError("trace_beats must be an integer between 0 and beats")
    if drug == "control":
        concentration = 0
    model = CiPABetaCell() if model is None else model
    p = ap_parameters(drug, concentration, sample, hergpath, hillpath, drop, optimized)
    state = model.initial_state(concentration)
    metrics, starts, checkpoints = [], [], []
    cycle_traces = []
    for beat in range(1, beats+1):
        if beat > 1 and beat % 250 == 1:
            print("concentration: {} beat: {}".format(concentration, beat))
            checkpoints.append(dict(beat=beat, **dict(zip(model.state_names, state))))
        trace = simulate_beat(state, p, cycle_length, dt, model=model, start_time=(beat-1)*cycle_length)
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
        state = model.reset_beat(trace[list(model.state_names)].iloc[-1].to_numpy())
    df = pd.DataFrame(metrics)
    keep = df.max_dv.idxmax() if df.max_dv.notna().any() else len(df)-1
    selected = dict(df.loc[keep])
    selected.update(EADs=int((df.max_dv > 0).sum()), noDepols=int(df.max_dv.isna().sum()),
                    noRepols=int(df.APD90.isna().sum()))
    # Save the trace from the selected beat, which need not be the last beat.
    selected_trace = simulate_beat(starts[keep], p, cycle_length, dt, model=model,
                                    start_time=(selected["beat"]-1)*cycle_length)
    if output is not None:
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        df.to_csv(output / "beat_metrics.csv", index=False)
        selected_trace.to_csv(output / "trace.csv", index=False)
        if cycle_traces:
            pd.concat(cycle_traces, ignore_index=True).to_csv(output / "trace_cycles.csv", index=False)
        pd.DataFrame(checkpoints).to_csv(output / "checkpoints.csv", index=False)
        pd.DataFrame([dict(beat=selected["beat"], **dict(zip(model.state_names, starts[keep]))) ]).to_csv(output / "states_max_dv.csv", index=False)
        (output / "run_config.json").write_text(json.dumps(dict(
            model=model.name, drug=drug, concentration_nM=concentration, sample=sample,
            beats=beats, nsave=nsave, trace_beats=trace_beats, cycle_length_ms=cycle_length,
            dt_ms=dt, optimized=optimized, drop_current=drop, 
            hergpath=str(hergpath), hillpath=str(hillpath), **model.metadata()), indent=2))
    return selected, selected_trace
