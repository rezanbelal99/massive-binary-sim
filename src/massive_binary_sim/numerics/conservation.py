"""Conservation monitors for the direct-orbit integrator.

For conservative PN orders (newtonian, 1PN, 2PN) the Newtonian mechanical energy
is not exactly conserved -- only the PN-accurate energy is -- but its *relative
drift* over an orbit is a sensitive integrator-quality diagnostic and is what the
Kepler-closure test checks.  For 2.5PN the energy must *decrease*; we then monitor
the agreement between the numerical dE/dt and the Peters (1964) luminosity.
"""
from __future__ import annotations

import numpy as np

from ..dynamics.gw import peters_dadt
from ..dynamics.kepler import orbital_energy_angular_momentum, state_to_elements
from .units import C, G


def track_conservation(t, y, m1, m2):
    """Return dict of arrays: E(t), L(t), relative drifts, and (a,e)(t)."""
    mu = G * (m1 + m2)
    E = np.empty_like(t)
    L = np.empty_like(t)
    a = np.empty_like(t)
    e = np.empty_like(t)
    for i, (ti, yi) in enumerate(zip(t, y)):
        Ei, Li = orbital_energy_angular_momentum(yi[:3], yi[3:], m1, m2)
        E[i], L[i] = Ei, Li
        el = state_to_elements(yi[:3], yi[3:], mu)
        a[i], e[i] = el["a"], el["e"]
    E0 = E[0] if E[0] != 0 else 1.0
    L0 = L[0] if L[0] != 0 else 1.0
    return dict(t=t, E=E, L=L, a=a, e=e,
                dE_rel=(E - E[0]) / abs(E0),
                dL_rel=(L - L[0]) / abs(L0),
                max_dE_rel=float(np.max(np.abs((E - E[0]) / abs(E0)))),
                max_dL_rel=float(np.max(np.abs((L - L[0]) / abs(L0)))))


def gw_luminosity_numeric_vs_peters(cons, m1, m2):
    """Compare the numerically-measured secular energy-loss rate to Peters (1964).

    Returns (dEdt_numeric, dEdt_peters, fractional_residual) using an orbit-averaged
    finite difference of E(t) versus  dE/dt = (dE/da) * (da/dt)_Peters.
    """
    t, E, a, e = cons["t"], cons["E"], cons["a"], cons["e"]
    if len(t) < 10:
        return np.nan, np.nan, np.nan
    # linear fit to E(t) over the span (secular slope)
    slope = np.polyfit(t, E, 1)[0]
    a_mean = float(np.mean(a))
    e_mean = float(np.mean(e))
    # E = -G m1 m2 / (2 a)  =>  dE/da = G m1 m2 / (2 a^2)
    dEda = G * m1 * m2 / (2.0 * a_mean**2)
    dadt = peters_dadt(a_mean, e_mean, m1, m2)
    dEdt_peters = dEda * dadt
    resid = (slope - dEdt_peters) / dEdt_peters if dEdt_peters != 0 else np.nan
    return slope, dEdt_peters, resid
