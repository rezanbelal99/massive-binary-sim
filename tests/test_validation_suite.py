"""The Section-5 validation suite as pytest tests.

Run:  pytest -q
Each test prints its quantitative detail; failures are NOT suppressed.
"""
import numpy as np
import pytest

from massive_binary_sim.analysis import validation as V


def _report(d):
    print("\n  " + d["detail"])
    return d["passed"]


def test_1_kepler_closure():
    d = V.test_kepler_closure(n_orbits=300)
    assert _report(d), "energy/AM/a must close to < 1e-6 over 300 dissipationless orbits"


def test_2_lane_emden_analytic():
    d = V.test_lane_emden_analytic()
    assert _report(d), "Lane-Emden n=0,1,5 must match closed forms to < 1e-4"


def test_3_solar_structure_eddington():
    d = V.test_solar_structure()
    assert _report(d), "Eddington solar model P_c,T_c,rho_c within factor 2; Teff within 15%"


def test_4_periastron_advance_1pn():
    d = V.test_periastron_advance()
    assert _report(d), "1PN periastron advance must match 6 pi GM/(c^2 a (1-e^2)) to < 2%"


def test_5_hulse_taylor_decay():
    d = V.test_hulse_taylor()
    assert _report(d), "Peters dPb/dt must match the PSR B1913+16 GR prediction to < 1%"


def test_6_energy_conservation_2pn():
    d = V.test_energy_conservation_2pn(n_orbits=200)
    assert _report(d), "2PN conservative integration: no secular energy drift above 1e-7"


def test_7_eddington_limit_flag():
    d = V.test_eddington_flag()
    assert _report(d), "code must flag Gamma_Edd > 1"


@pytest.mark.slow
def test_8_convergence_order():
    from massive_binary_sim.analysis.convergence import orbit_convergence
    c = orbit_convergence()
    print("\n  " + c["detail"])
    # observed order should be positive and the solution should converge
    assert c["rel_errors"][-1] < 1e-6
    assert np.nanmean(c["observed_order_per_decade"]) > 0.5
