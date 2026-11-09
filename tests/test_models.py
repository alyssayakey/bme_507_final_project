import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from cipa.models import ord, herg, initial_ap, parameter_vector, ap_parameters
from cipa.common import DATA, read_pairs
from cipa.herg import rate_matrix, MilnesModel
from cipa.ap import simulate_beat, beat_metrics


def test_herg_mass_and_matrix_solution():
    from scipy.linalg import expm
    p = read_pairs(DATA / "herg/hergmod_pars.txt")
    p.update(Kmax=5., Ku=.001, n=1., halfmax=100., Vhalf=-8.)
    y = np.r_[np.arange(1, 10)/45, 100., 0.]
    dy, _ = herg.evaluate(0., y, parameter_vector(p, "herg"))
    assert abs(dy[:9].sum()) < 1e-12
    m = rate_matrix(p, 100., 0.)
    np.testing.assert_allclose(m @ y[:9], dy[:9], atol=1e-14)
    sol = solve_ivp(lambda t, x: m @ x, (0, 1000), y[:9], method="LSODA", rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(expm(m*1000) @ y[:9], sol.y[:, -1], atol=1e-9)


def test_milnes_control_fraction_is_one():
    model = MilnesModel()
    result = model.fractional_current({}, 0, np.array([945., 955., 1000., 10940.]))
    np.testing.assert_allclose(result, 1.)
    states = model.sweeps({}, 0, np.array([0., 940., 10940., 25000.]))
    np.testing.assert_allclose(states.sum(axis=2), 1., atol=1e-9)


def test_zero_drug_and_qnet_current_identity():
    state = initial_ap()
    p = parameter_vector(ap_parameters())
    dy, currents = ord.evaluate(1., state, p)
    assert np.isfinite(dy).all()
    np.testing.assert_allclose(dy[49], currents[[1, 2, 3, 4, 5, 6]].sum())
    assert abs(dy[39:48].sum()) < 1e-12


def test_control_action_potential_and_qnet_integral():
    trace = simulate_beat(initial_ap(), ap_parameters(), dt=.25)
    metric = beat_metrics(trace)
    assert 20 < metric["vmax"] < 60
    assert 150 < metric["APD90"] < 600
    assert metric["APD90"] > metric["APD50"] > 0
    assert trace.cai.min() > 0
    # Integral from output currents agrees with the extra ODE state.
    net = trace[["INaL", "ICaL", "Ito", "IKr", "IKs", "IK1"]].sum(axis=1)
    np.testing.assert_allclose(trapezoid(net, trace.time), metric["qNet"], rtol=.005)
