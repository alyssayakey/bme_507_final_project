"""Cell-model boundary shared by pacing, output, and future signaling models.

All times are ms. evaluate receives absolute experiment time plus an explicit
stimulus flag; signaling must not use a beat-local clock. Extra signaling
states belong in state_names and persist between beats.
"""
from typing import Protocol
import numpy as np
from .models import AP_STATES, CURRENTS, initial_ap, parameter_vector, ord


class CellModel(Protocol):
    name: str
    state_names: tuple[str, ...]
    output_names: tuple[str, ...]
    stimulus_duration: float

    def initial_state(self, concentration: float) -> np.ndarray: ...
    def evaluate(self, time: float, state: np.ndarray, parameters,
                 stimulating: bool) -> tuple[np.ndarray, np.ndarray]: ...
    def reset_beat(self, state: np.ndarray) -> np.ndarray: ...


class ORdCell:
    """Unmodified CiPA ORd equations behind the extensible cell interface."""
    name = "cipa-ord"
    state_names = AP_STATES
    output_names = CURRENTS
    stimulus_duration = .5

    def initial_state(self, concentration=0):
        return initial_ap(concentration)

    def evaluate(self, time, state, parameters, stimulating):
        # Only the original stimulus calculation uses the kernel's time input.
        vector = parameter_vector(parameters) if isinstance(parameters, dict) else parameters
        return ord.evaluate(0. if stimulating else 1., state, vector)

    def reset_beat(self, state):
        state = state.copy()
        state[self.state_names.index("qNet")] = 0.
        return state
