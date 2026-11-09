import numpy as np
import pandas as pd
import pytest
from cipa.hill import hill_block, fit_hill, sample_hill, process_patchclamp
from cipa.herg import bootstrap, experiments, mean_data, encode, decode, PARAMETERS
from cipa.common import DATA, parse_range
from cipa.models import ap_parameters
from cipa.analysis import qnet_intervals, combine


def test_hill_known_parameters_and_reproducible_chain():
    x = np.geomspace(1, 10000, 20)
    y = hill_block(x, 2., 1.3)
    par, cov, var, x, y = fit_hill(x, y)
    np.testing.assert_allclose(par, [2, 1.3], atol=1e-6)
    a, _ = sample_hill(x, y, par, cov, var, 50, 100, 2, 20)
    b, _ = sample_hill(x, y, par, cov, var, 50, 100, 2, 20)
    np.testing.assert_array_equal(a, b)
    assert a.shape == (50, 2)
    np.testing.assert_allclose(hill_block([0., 100.], 2., 1.), [0., 50.])


def test_bootstrap_stratification_and_duplicate_cells():
    data = pd.read_csv(DATA / "herg/bepridil.csv")
    a, b = bootstrap(data, 5, 100), bootstrap(data, 5, 100)
    np.testing.assert_array_equal(a, b)
    cells = experiments(data)
    for indices in a:
        np.testing.assert_array_equal(cells.conc.to_numpy()[indices], cells.conc.to_numpy())
    original = mean_data(data)
    repeated = mean_data(data, np.repeat(np.arange(len(cells)), 2))
    pd.testing.assert_frame_equal(original, repeated)
    values = [5, .001, 1, .15, -8]
    np.testing.assert_allclose(decode(encode(values)), values)


def test_range_errors():
    assert parse_range("1-3,5", integer=True) == [1, 2, 3, 5]
    assert parse_range("0.5,1.5") == [.5, 1.5]
    for spec in ("-1", "3-1", "abc", "1,,2", "nan"):
        with pytest.raises(ValueError):
            parse_range(spec)
    with pytest.raises(ValueError):
        parse_range("1.5", integer=True)


def test_parameter_scaling_and_drop(tmp_path):
    from cipa.common import write_pairs
    folder = tmp_path / "drug"
    folder.mkdir()
    write_pairs(dict(Kmax=5, Ku=.001, n=1, halfmax=100, Vhalf=-8), folder / "pars.txt")
    pd.DataFrame([dict(INa_IC50=100., INa_h=2.)]).to_csv(folder / "IC50_optimal.csv", index=False)
    p = ap_parameters("drug", 100, hergpath=tmp_path, hillpath=tmp_path)
    assert p["GNafc"] == 2
    assert p["GKrfc"] == pytest.approx(1.1/1.114)
    p = ap_parameters("drug", 100, hergpath=tmp_path, hillpath=tmp_path, drop="INa")
    assert p["GNafc"] == 1
    p = ap_parameters("drug", 100, hergpath=tmp_path, hillpath=tmp_path, drop="hERG")
    assert p["Kmax"] == p["Kt"] == 0


def test_qnet_excludes_depolarization_failure_and_detects_gaps():
    fixed = pd.DataFrame([dict(drug="drug", dose=1, qNet=100, max_dv=-1)])
    samples = pd.DataFrame(dict(drug=["drug"]*3, dose=[1]*3, sample=[1, 2, 3],
                                qNet=[90., 100., 1e6], max_dv=[-1, -1, np.nan]))
    result = qnet_intervals(fixed, samples)
    assert result.upper.iloc[0] < 100
    with pytest.raises(ValueError, match="Missing"):
        qnet_intervals(fixed, samples.iloc[[0, 2]])


def test_combine_detects_missing_samples(tmp_path):
    folder = tmp_path / "drug/00001"
    folder.mkdir(parents=True)
    pd.DataFrame([dict(drug="drug", sample=1, dose=1, qNet=10)]).to_csv(folder / "metrics.csv", index=False)
    assert len(combine(tmp_path, 1)) == 1
    with pytest.raises(ValueError, match="Missing"):
        combine(tmp_path, 2)


def test_patchclamp_concentrations_are_converted():
    result = process_patchclamp()
    assert not result.units.eq("uM").any()
    assert {"ICaL", "INa", "INaL", "Ito"} <= set(result.channel)
