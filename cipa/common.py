"""Shared paths, input validation, and portable text formats."""
from pathlib import Path
import re
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parent / "data"
DRUGS = "dofetilide,cisapride,bepridil,verapamil,terfenadine,ranolazine,sotalol,mexiletine,quinidine,ondansetron,diltiazem,chlorpromazine"


def read_pairs(path):
    return {name: float(value) for name, value in
            (line.split() for line in Path(path).read_text().splitlines() if line.strip())}


def write_pairs(values, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("".join(f"{k} {v:.12g}\n" for k, v in values.items()))


def parse_range(spec, integer=False):
    values = []
    for item in spec.split(","):
        match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)(?:-(\d+(?:\.\d+)?))?\s*", item)
        if not match:
            raise ValueError(f"Invalid nonnegative range: {item!r}")
        start = float(match[1])
        stop = float(match[2]) if match[2] else start
        if stop < start or (integer and (not start.is_integer() or not stop.is_integer())):
            raise ValueError(f"Invalid range: {item!r}")
        values.extend(np.arange(start, stop + 1e-10, 1).tolist())
    return [int(v) for v in values] if integer else values


def intervals(samples, optimal=None, alpha=0.05):
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    rows = []
    for name in samples:
        valid = samples[name].dropna()
        limits = valid.quantile([alpha / 2, 1 - alpha / 2]).to_numpy() if len(valid) else [np.nan]*2
        rows.append(dict(param=name, value=optimal.get(name, np.nan) if optimal is not None else np.nan,
                         lower=limits[0], upper=limits[1]))
    return pd.DataFrame(rows)
