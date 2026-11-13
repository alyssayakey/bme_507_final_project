"""Regression checks for future signaling states and experiment-wide time."""
import numpy as np
import pandas as pd
import pytest
from cipa.cell import ORdCell
from cipa.ap import simulate
from cipa.cli import parser, main


class ClockCell(ORdCell):
    # A mathematical test state, not a biological signaling approximation.
    name = "test-clock"
    state_names = ORdCell.state_names + ("clock",)

    def initial_state(self, concentration=0):
        return np.r_[super().initial_state(concentration), 0.]

    def evaluate(self, time, state, parameters, stimulating):
        # LSODA perturbs states to estimate its Jacobian.
        assert abs(state[-1]-time) < .01
        dy, out = super().evaluate(time, state[:-1], parameters, stimulating)
        return np.r_[dy, 1.], out


def test_extra_state_persists_and_replay_uses_absolute_time(tmp_path):
    model = ClockCell()
    _, trace = simulate(beats=3, nsave=2, trace_beats=3, output=tmp_path, model=model)
    cycles = pd.read_csv(tmp_path / "trace_cycles.csv")
    np.testing.assert_allclose(cycles.clock, cycles.simulation_time, atol=1e-5)
    assert cycles.clock.iloc[-1] == pytest.approx(6000)
    assert "clock" in trace
    starts = cycles.groupby("beat").first()
    np.testing.assert_allclose(starts.qNet, 0.)


def test_quinidine_cli_defaults_and_control():
    p = parser()
    for command in ("ap", "herg-fit", "hill-fit", "bootstrap", "herg-summary"):
        assert p.parse_args([command]).drug == "quinidine"
    assert p.parse_args(["ap"]).dose == "1"
    assert p.parse_args(["ap", "-d", "control"]).drug == "control"
    with pytest.raises(SystemExit):
        p.parse_args(["ap", "-d", "bepridil"])
    with pytest.raises(SystemExit):
        main(["tdp-error", "results/ap/fixed"])
