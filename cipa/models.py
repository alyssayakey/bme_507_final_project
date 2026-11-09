"""Named initial states and parameter vectors for the translated equations."""
import numpy as np
import pandas as pd
from .common import DATA, read_pairs
from . import _herg_equations as herg, _ord_equations as ord

AP_STATES = tuple(read_pairs(DATA / "ap/newordherg_states_CL2000.txt")) + ("qNet",)
HERG_STATES = tuple(read_pairs(DATA / "herg/hergmod_states.txt"))
CURRENTS = ("INa", "INaL", "Ito", "ICaL", "IKr", "IKs", "IK1", "dv")
SCALE_NAMES = dict(INa="GNafc", INaL="GNaLfc", ICaL="PCafc", IK1="GK1fc", IKs="GKsfc", Ito="Gtofc")


def initial_ap(concentration=0):
    values = read_pairs(DATA / "ap/newordherg_states_CL2000.txt")
    values.update(D=concentration, qNet=0)
    return np.array([values[n] for n in AP_STATES])


def parameter_vector(values, model="ap"):
    names = ord.PARAMETER_NAMES if model == "ap" else herg.PARAMETER_NAMES
    return np.array([values["T" if n == "Temp" else n] for n in names], dtype=float)


def ap_parameters(drug="control", concentration=0, sample=0,
                  hergpath="results/herg", hillpath="results/hill", drop=None, optimized=True):
    from pathlib import Path
    if concentration < 0 or sample < 0:
        raise ValueError("concentration and sample must be nonnegative")
    if drop is not None and drop not in {"hERG", *SCALE_NAMES}:
        raise ValueError(f"Unknown current: {drop}")
    p = read_pairs(DATA / "ap/newordherg_pars.txt")
    p.update(GKrfc=1.1, T=37)
    if drug == "control":
        p.update(Kmax=0, Kt=0)
    else:
        hp, ip = Path(hergpath) / drug, Path(hillpath) / drug
        if sample == 0:
            binding = read_pairs(hp / "pars.txt")
            hill = pd.read_csv(ip / "IC50_optimal.csv").iloc[0]
        else:
            binding = pd.read_csv(hp / "boot_pars.csv").iloc[sample-1].to_dict()
            hill = pd.read_csv(ip / "IC50_samples.csv").iloc[sample-1]
        for name, value in binding.items():
            if name in p:
                if not np.isfinite(value):
                    raise ValueError(f"Missing hERG parameter {name}, sample {sample}")
                p[name] = value
        p["Kt"] = 3.5e-5
        if drop == "hERG":
            p.update(Kmax=0, Kt=0)
        for channel, scale in SCALE_NAMES.items():
            ic50, h = hill.get(channel + "_IC50", np.nan), hill.get(channel + "_h", np.nan)
            if channel != drop and concentration > 0 and np.isfinite(ic50) and np.isfinite(h):
                if ic50 <= 0 or h < 0:
                    raise ValueError(f"Invalid Hill parameters for {channel}")
                p[scale] = 1 + (concentration / ic50)**h
    if optimized:
        for name, value in dict(GKrfc=1.114, GNaLfc=2.661, GKsfc=1.870, GK1fc=1.698, PCafc=1.007).items():
            p[name] /= value
    return p
