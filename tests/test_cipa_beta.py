"""Hybrid equation, conservation, pharmacology, and trajectory checks."""
from pathlib import Path
import ast
import importlib.util
import numpy as np
import pandas as pd
import pytest
from cipa.models import initial_ap, ap_parameters, parameter_vector, ord, AP_STATES
from cipa_beta.cell import CiPABetaCell, EXTRA_STATES, SIGNAL_OFFSET
from cipa_beta.ap import simulate, simulate_beat
from cipa_beta import signaling


@pytest.mark.parametrize('celltype', [0, 1, 2])
def test_disabled_hybrid_matches_cipa_with_drug_and_scaling(celltype):
    rng = np.random.default_rng(12)
    model = CiPABetaCell(regulation=False)
    p = ap_parameters()
    p.update(celltype=celltype, Kmax=400., Ku=.006, n=.85, halfmax=78000., Kt=3.5e-5, Vhalf=-65.,
             GNafc=1.6, GNaLfc=.8, PCafc=1.2, GKsfc=1.4, Gtofc=1.8)
    for time in (0., 1.):
        y = model.initial_state(3237.)
        y[0] = rng.uniform(-85, 30)
        a, b = model.evaluate(time, y, p, time == 0)
        expected, currents = ord.evaluate(time, y[:50], parameter_vector(p))
        np.testing.assert_allclose(a[:50], expected, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(b[:8], currents, rtol=1e-12, atol=1e-12)
        assert abs(a[39:48].sum()) < 1e-12
        assert a[49] == pytest.approx(b[[1,2,3,4,5,6]].sum())


def test_signaling_extraction_matches_supplied_translation():
    path = Path(__file__).resolve().parents[1]/'beta-signaling/gong_beta.py'
    spec = importlib.util.spec_from_file_location('source_gong', path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    for iso in (0., 1.):
        a, b = signaling.constants(iso), source.constants(iso)
        y = signaling.initial_state()
        np.testing.assert_array_equal(y, source.initial_state()[58:])
        np.testing.assert_allclose(signaling.signaling(0., y, a), source.signaling(0., y, b), atol=0, rtol=0)
        sa = signaling.effective(y, signaling.SimpleNamespace(const_signaling=a))
        sb = source.effective(np.r_[np.zeros(58), y], signaling.SimpleNamespace(const_signaling=b))
        for name in ('fICaLP','fIKsP','fPLBP','fTnIP','fINaP','fINaKP','fRyRP','fIKurP'):
            assert getattr(sa, name) == getattr(sb, name)


def test_active_equations_match_gong_regulation_with_cipa_baseline():
    # Independent evaluation of the supplied Gong ephys function. Match CiPA
    # constants/scaling, then substitute the independently computed CiPA IKr.
    path = Path(__file__).resolve().parents[1]/'beta-signaling/gong_beta.py'
    spec = importlib.util.spec_from_file_location('reference_gong', path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    tree = ast.parse(path.read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'ephys')
    for n in ast.walk(function):
        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'IKr' for t in n.targets):
            n.value = ast.Name(id='reference_ikr', ctx=ast.Load())
    namespace = source.__dict__.copy()
    namespace.update(source.initialize_constants(0, 0))
    vc = 1000*3.14*.0011**2*.01
    namespace.update(R=8314., T=310., F=96485., Acap=2*(2*3.14*.0011**2+2*3.14*.0011*.01),
                     vmyo=.68*vc, vss=.02*vc, vnsr=.0552*vc, vjsr=.0048*vc)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])), '<reference>', 'exec'), namespace)
    model = CiPABetaCell(iso=1., precondition_seconds=0)
    p = ap_parameters()
    scales = dict(GNa=('GNafc',75.), GNaL=('GNaLfc',.0075), Gto=('Gtofc',.02),
                  GKs=('GKsfc',.0034), GK1=('GK1fc',.1908), PCa=('PCafc',.0001))
    for name,(key,value) in scales.items():
        namespace[name] = value/p[key]
    gong_names = ast.parse(path.read_text())
    ef = next(n for n in gong_names.body if isinstance(n, ast.FunctionDef) and n.name == 'ephys')
    state_names = tuple(n.id for n in ef.body[0].targets[0].elts)
    mapping = dict(V='v', Nai='nai', Nass='nass', Ki='ki', Kss='kss', Cai='cai',
                   Cass='cass', Cansr='cansr', Cajsr='cajsr')
    for voltage in (-80., -20., 25.):
        y = model.initial_state()
        y[0] = voltage
        # Use nonzero phosphorylation to exercise all transferred populations.
        controls = np.array([.7,.6,.5,.4,.3,.8,.9,.2])
        from cipa_beta import _equations
        dy, out = _equations.evaluate(1., y[:SIGNAL_OFFSET], parameter_vector(p), controls, .12, True)
        namespace['reference_ikr'] = out[4]
        settings = signaling.SimpleNamespace(**dict(zip(
            ('fICaLP','fIKsP','fPLBP','fTnIP','fINaP','fINaKP','fRyRP','fIKurP'), controls)), Whole_cell_PP1=.12)
        values = dict(zip(model.state_names, y))
        refstate = np.array([values.get(mapping.get(n,n), .2) for n in state_names])
        expected = namespace['ephys'](1., refstate, 0., settings)
        for i, name in enumerate(state_names):
            if name in ('xrf', 'xrs'):
                continue
            index = model.state_names.index(mapping.get(name,name))
            np.testing.assert_allclose(dy[index], expected[i], rtol=1e-10, atol=1e-10, err_msg=name)


def test_active_run_saves_signals_cycles_and_resets_only_qnet(tmp_path):
    model = CiPABetaCell(iso=1.)
    metric, trace = simulate(beats=2, nsave=1, trace_beats=2, model=model, output=tmp_path)
    assert np.isfinite(trace.to_numpy()).all()
    assert trace.cai.min() > 0
    assert np.isfinite(metric['APD90'])
    cycles = pd.read_csv(tmp_path/'trace_cycles.csv')
    assert cycles.beat.nunique() == 2
    assert cycles.time.iloc[-1] == 4000
    assert 'beta_cAMP_cav' in cycles and 'fICaLP' in cycles
    assert cycles.loc[cycles.beat == 2, 'qNet'].iloc[0] == 0
    y = model.initial_state()
    reset = model.reset_beat(y)
    np.testing.assert_array_equal(reset[SIGNAL_OFFSET:], y[SIGNAL_OFFSET:])


def test_active_drug_block_scales_all_phosphorylated_populations():
    m = CiPABetaCell(iso=1.)
    y = m.initial_state(3237.)
    y[0] = -20.
    p = ap_parameters()
    p.update(Kmax=400., Ku=.006, n=.85, halfmax=78000., Kt=3.5e-5, Vhalf=-65.)
    dy, currents = m.evaluate(1., y, p, False)
    reference, _ = ord.evaluate(1., y[:50], parameter_vector(p))
    np.testing.assert_allclose(dy[39:49], reference[39:49], rtol=1e-12, atol=1e-12)
    for scale, index in [('GNafc',0), ('GNaLfc',1), ('Gtofc',2), ('PCafc',3), ('GKsfc',5), ('GK1fc',6)]:
        altered = dict(p)
        altered[scale] *= 2
        _, blocked = m.evaluate(1., y, altered, False)
        assert blocked[index] == pytest.approx(currents[index]/2, abs=1e-12)
