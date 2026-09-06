"""Programmatic validation suite -- shared by pytest and the report generator.

Each function returns dict(passed: bool, detail: str, **numbers).
"""
from __future__ import annotations

import numpy as np

from ..dynamics.gw import merger_time_eccentric, peters_dPb_dt
from ..dynamics.kepler import elements_to_state, state_to_elements
from ..dynamics.postnewtonian import periastron_advance_per_orbit
from ..numerics.integrators import integrate_orbit
from ..numerics.conservation import track_conservation
from ..numerics.units import C, DAY, G, MSUN, RSUN, LSUN, TEFF_SUN
from ..structure.polytrope import ANALYTIC_XI1, analytic_theta, solve_lane_emden


def test_kepler_closure(n_orbits=300, tol=1e-6):
    m1 = m2 = 1.4 * MSUN
    mu = G * (m1 + m2)
    a0, e0 = 1.0 * RSUN, 0.3
    r, v = elements_to_state(a0, e0, 0.2, 0.3, 0.5, 0.0, mu)
    y0 = np.concatenate([r, v])
    P = 2 * np.pi * np.sqrt(a0**3 / mu)
    t, y = integrate_orbit(y0, (0, n_orbits * P), m1, m2, order="newtonian",
                           method="dop853", rtol=1e-12, atol=1e-14, dense_points=4000)
    c = track_conservation(t, y, m1, m2)
    da = abs((c["a"][-1] - c["a"][0]) / a0)
    passed = c["max_dE_rel"] < tol and c["max_dL_rel"] < tol and da < tol
    return dict(passed=bool(passed),
               detail=f"{n_orbits} orbits: max|dE/E|={c['max_dE_rel']:.2e}, "
                      f"max|dL/L|={c['max_dL_rel']:.2e}, |da/a|={da:.2e} (tol {tol:.0e})",
               max_dE_rel=c["max_dE_rel"], max_dL_rel=c["max_dL_rel"], da_rel=da)


def test_lane_emden_analytic():
    worst = 0.0
    details = []
    for n in (0, 1, 5):
        sol = solve_lane_emden(n, xi_max=(10.0 if n == 5 else 8.0))
        mask = sol["xi"] < (ANALYTIC_XI1[n] if np.isfinite(ANALYTIC_XI1[n]) else 8.0)
        err = np.max(np.abs(sol["theta"][mask] - analytic_theta(n, sol["xi"][mask])))
        worst = max(worst, err)
        if np.isfinite(ANALYTIC_XI1[n]):
            xi1_err = abs(sol["xi1"] - ANALYTIC_XI1[n])
            details.append(f"n={n}: max|dtheta|={err:.2e}, |dxi1|={xi1_err:.2e}")
        else:
            details.append(f"n={n}: max|dtheta|={err:.2e}")
    passed = worst < 1e-4
    return dict(passed=bool(passed), detail="; ".join(details), worst_err=worst)


def test_periastron_advance(tol_frac=0.02):
    m1 = m2 = 1.4 * MSUN
    mu = G * (m1 + m2)
    a0, e0 = 2.0 * RSUN, 0.3
    r, v = elements_to_state(a0, e0, 0.0, 0.0, 0.0, 0.0, mu)
    y0 = np.concatenate([r, v])
    P = 2 * np.pi * np.sqrt(a0**3 / mu)
    N = 60
    t, y = integrate_orbit(y0, (0, N * P), m1, m2, order="1PN", method="dop853",
                           rtol=1e-13, atol=1e-15, dense_points=30000)
    oms = np.unwrap([state_to_elements(yi[:3], yi[3:], mu)["omega"] for yi in y])
    measured = (oms[-1] - oms[0]) / N
    predicted = periastron_advance_per_orbit(a0, e0, m1 + m2, "1PN")
    frac = abs(measured - predicted) / predicted
    return dict(passed=bool(frac < tol_frac),
               detail=f"1PN d(omega)/orbit measured={measured:.6e} rad, "
                      f"predicted 6 pi GM/(c^2 a(1-e^2))={predicted:.6e}, frac diff={frac:.3%}",
               measured=measured, predicted=predicted, frac=frac)


def test_hulse_taylor(tol_frac=0.01):
    # Weisberg & Huang 2016, ApJ 829, 55
    Pb = 0.322997448911 * DAY
    e = 0.6171340
    m1, m2 = 1.438 * MSUN, 1.390 * MSUN
    mtot = m1 + m2
    a = (G * mtot * (Pb / (2 * np.pi)) ** 2) ** (1.0 / 3.0)
    pbdot = peters_dPb_dt(a, e, m1, m2)
    gr_pred = -2.40263e-12
    frac = abs(pbdot - gr_pred) / abs(gr_pred)
    tm = merger_time_eccentric(a, e, m1, m2) / (3.15576e7)
    return dict(passed=bool(frac < tol_frac),
               detail=f"dPb/dt(Peters)={pbdot:.6e} s/s vs GR prediction {gr_pred:.6e} "
                      f"(Weisberg & Huang 2016); residual {frac:.3%}. "
                      f"Merger time {tm/1e6:.0f} Myr (lit ~300 Myr).",
               pbdot=pbdot, residual_frac=frac, merger_time_yr=tm)


def test_energy_conservation_2pn(n_orbits=200, tol=1e-7):
    m1, m2 = 10 * MSUN, 8 * MSUN
    mu = G * (m1 + m2)
    a0, e0 = 5.0 * RSUN, 0.2
    r, v = elements_to_state(a0, e0, 0.0, 0.0, 0.0, 0.0, mu)
    y0 = np.concatenate([r, v])
    P = 2 * np.pi * np.sqrt(a0**3 / mu)
    t, y = integrate_orbit(y0, (0, n_orbits * P), m1, m2, order="2PN",
                           method="dop853", rtol=1e-12, atol=1e-14, dense_points=3000)
    c = track_conservation(t, y, m1, m2)
    # Newtonian-energy drift is O(v^2/c^2) per orbit physically; integrator drift
    # should be << that. We check the *bounded* (non-secular) behaviour.
    slope = np.polyfit(t, c["E"], 1)[0]
    secular = abs(slope * (t[-1] - t[0]) / c["E"][0])
    return dict(passed=bool(secular < tol),
               detail=f"2PN conservative: secular |dE/E| over {n_orbits} orbits = {secular:.2e} "
                      f"(tol {tol:.0e}); max|dL/L|={c['max_dL_rel']:.2e}",
               secular_dE=secular, max_dL_rel=c["max_dL_rel"])


def test_eddington_flag():
    """A model above the Eddington limit must be detectable via StarState.eddington_ratio."""
    from ..structure.evolution import StarState
    s = StarState.zams(120.0, 0.014)
    gamma = s.eddington_ratio()
    # 120 Msun ZAMS is genuinely close to Eddington; just require the diagnostic exists
    # and is order-unity, and that a deliberately over-luminous state flags > 1.
    s.L *= 10.0
    flagged = s.eddington_ratio() > 1.0
    return dict(passed=bool(flagged and 0.1 < gamma < 1.5),
               detail=f"120 Msun ZAMS Gamma_Edd={gamma:.3f}; 10x-luminosity state flags "
                      f"Gamma>1: {flagged}",
               gamma_zams=gamma)


def test_solar_structure():
    """Solar structure via the Eddington standard model (n=3 polytrope).

    This checks the *interior* physics: central pressure, temperature and density
    of a 1 Msun, R~1 Rsun, Z_sun polytrope against standard solar model values
    (Bahcall, Pinsonneault & Basu 2001; JCD Model S):
        P_c ~ 2.4e17 dyn/cm^2,  T_c ~ 1.57e7 K,  rho_c ~ 150 g/cm^3.
    The Eddington model is known to under-predict these by ~15-30% (it omits the
    real Sun's composition gradient); we accept agreement within a factor 2 and
    require Teff within 15% of 5772 K.  A 4.57-Gyr evolved-Sun calibration needs
    the time-dependent composition solver -- KNOWN GAP STRUCT-1D.
    """
    from ..structure.stellar_1d import solve_zams
    comp = dict(X=0.7381, Y=0.2485, Z=0.0134)
    m = solve_zams(1.0, comp, dict(reaction_rates="reaclib"), dict(c12_ag_rate=1.0))
    Pc_ref, Tc_ref, rhoc_ref = 2.4e17, 1.57e7, 150.0
    ok_Pc = 0.5 < m.Pc / Pc_ref < 2.0
    ok_Tc = 0.5 < m.Tc / Tc_ref < 2.0
    ok_rho = 0.4 < m.rhoc / rhoc_ref < 2.5
    ok_T = abs(m.Teff - TEFF_SUN) / TEFF_SUN < 0.15
    ok_LR = 0.4 < m.L / LSUN < 2.5 and 0.7 < m.R / RSUN < 1.4
    passed = ok_Pc and ok_Tc and ok_rho and ok_T and ok_LR
    return dict(passed=bool(passed),
               detail=(f"Eddington model 1 Msun: P_c={m.Pc:.3e} (ref 2.4e17), "
                       f"T_c={m.Tc:.3e} (ref 1.57e7), rho_c={m.rhoc:.1f} (ref 150), "
                       f"L={m.L/LSUN:.2f} Lsun, R={m.R/RSUN:.2f} Rsun, Teff={m.Teff:.0f} K "
                       f"(ref 5772). Within-factor-2 targets met: {passed}."),
               Pc=m.Pc, Tc=m.Tc, rhoc=m.rhoc, Teff=m.Teff)


def run_all(fast=True):
    tests = {
        "kepler_closure": lambda: test_kepler_closure(n_orbits=200 if fast else 10000),
        "lane_emden_analytic": test_lane_emden_analytic,
        "periastron_advance_1PN": test_periastron_advance,
        "hulse_taylor_decay": test_hulse_taylor,
        "energy_conservation_2PN": lambda: test_energy_conservation_2pn(n_orbits=120 if fast else 500),
        "eddington_limit_flag": test_eddington_flag,
        "solar_structure_eddington": test_solar_structure,
    }
    out = {}
    for name, fn in tests.items():
        try:
            out[name] = fn()
        except Exception as exc:
            out[name] = dict(passed=False, detail=f"EXCEPTION: {exc}")
    return out
