"""Tidal evolution: circularisation and synchronisation.

Hut (1981), A&A 99, 126 -- equilibrium-tide (weak-friction) secular equations.
Zahn (1977, 1989) -- calibration of the tidal dissipation:
  * radiative envelopes: dynamical tide, E_2 coupling constant (Zahn 1975);
  * convective envelopes: turbulent viscosity, tau_conv.
Hurley, Tout & Pols (2002), MNRAS 329, 897, sec. 2.3 -- practical (k/T) formulae
  used here for the two dissipation regimes.
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import G, RSUN, YEAR


def _k_over_T_convective(star, m_other, a):
    """(k/T) for a convective envelope (Hurley+ 2002 eq. 30, Rasio+ 1996)."""
    M, R, L = star.m, star.R, star.L
    Menv = max(0.1 * M, M - star.m_he_core) if star.m_he_core else 0.3 * M
    Renv = 0.3 * R
    tau_conv = (Menv * Renv * (R - Renv) / (3.0 * L)) ** (1.0 / 3.0)
    P_orb = 2 * np.pi * np.sqrt(a**3 / (G * (M + m_other)))
    P_tid = 1.0 / abs(1.0 / P_orb - star.spin * star.omega_crit() / (2 * np.pi) + 1e-30)
    f_conv = min(1.0, (P_tid / (2.0 * tau_conv)) ** 2)
    return (2.0 / 21.0) * (f_conv / tau_conv) * (Menv / M)


def _k_over_T_radiative(star, m_other, a):
    """(k/T) for a radiative envelope, per second (Hurley+ 2002 eqs. 42-43).

    (k/T) = 1.9782e4 * sqrt(M R^2 / a^5) * (1 + q2)^(5/6) * E2   [yr^-1]
    with M, R, a in solar units;  E2 = 1.592e-9 M^2.84 (M in Msun).
    """
    from ..numerics.units import MSUN
    M_s, R_s, a_s = star.m / MSUN, star.R / RSUN, a / RSUN
    q2 = m_other / star.m
    E2 = 1.592e-9 * M_s ** 2.84
    kT_per_yr = 1.9782e4 * np.sqrt(M_s * R_s**2 / a_s**5) * (1.0 + q2) ** (5.0 / 6.0) * E2
    return kT_per_yr / YEAR


def tidal_derivatives(star, m_other, a, e, envelope="radiative"):
    """Return (da/dt, de/dt, dOmega/dt) from equilibrium tides (Hut 1981 eqs. 9-11).

    Uses the standard f1..f5(e) polynomials (Hut 1981 eq. 2.9-2.13).
    """
    M, R = star.m, star.R
    kT = (_k_over_T_convective(star, m_other, a) if envelope == "convective"
          else _k_over_T_radiative(star, m_other, a))
    n_orb = np.sqrt(G * (M + m_other) / a**3)
    Omega = star.spin * star.omega_crit()

    e2 = e * e
    oe = 1.0 - e2
    f1 = 1 + 31/2*e2 + 255/8*e2**2 + 185/16*e2**3 + 25/64*e2**4
    f2 = 1 + 15/2*e2 + 45/8*e2**2 + 5/16*e2**3
    f3 = 1 + 15/4*e2 + 15/8*e2**2 + 5/64*e2**3
    f4 = 1 + 3/2*e2 + 1/8*e2**2
    f5 = 1 + 3*e2 + 3/8*e2**2

    q_other = m_other / M
    R_a = R / a
    common = kT * q_other * (1.0 + q_other) * R_a ** 8

    dadt = -6.0 * common * a / oe ** (15.0 / 2.0) * (
        f1 - oe ** 1.5 * f2 * Omega / n_orb)
    dedt = -27.0 * common * e / oe ** (13.0 / 2.0) * (
        f3 - (11.0 / 18.0) * oe ** 1.5 * f4 * Omega / n_orb)
    # spin: dOmega/dt (Hut 1981 eq. 11), I ~ 0.1 M R^2 (rough moment of inertia)
    I = 0.08 * M * R**2
    dOmegadt = 3.0 * common * (M * R**2) / I * n_orb / oe ** 6 * (
        f2 - oe ** 1.5 * f5 * Omega / n_orb)
    return dadt, dedt, dOmegadt


def circularisation_timescale(star, m_other, a, e, envelope="radiative"):
    _, dedt, _ = tidal_derivatives(star, m_other, a, max(e, 1e-3), envelope)
    if dedt == 0:
        return np.inf
    return abs(max(e, 1e-3) / dedt)


def synchronisation_timescale(star, m_other, a, e, envelope="radiative"):
    Omega = star.spin * star.omega_crit()
    n_orb = np.sqrt(G * (star.m + m_other) / a**3)
    _, _, dOmegadt = tidal_derivatives(star, m_other, a, e, envelope)
    if dOmegadt == 0:
        return np.inf
    return abs((n_orb - Omega) / dOmegadt)
