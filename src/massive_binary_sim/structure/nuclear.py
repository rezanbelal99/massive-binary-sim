"""Nuclear energy generation rates.

Analytic rate formulae from Kippenhahn, Weigert & Weiss (2012), "Stellar
Structure and Evolution" (2nd ed.), sec. 18.4-18.5, which in turn follow
Angulo et al. (1999, NACRE) and Caughlan & Fowler (1988).  The
`theories.reaction_rates` switch selects the screening/temperature-fit set:

  nacre   : NACRE (Angulo+ 1999) coefficients
  reaclib : JINA REACLIB v2 -- here approximated by the same analytic fits with
            a +5% CNO / -3% 3-alpha offset (representative of NACRE->REACLIB
            differences; the true REACLIB rates are tabulated, see KNOWN_GAPS NUC-REACLIB)

The 12C(alpha,gamma)16O rate carries a user multiplier `c12_ag_rate` because its
uncertainty dominates the final C/O ratio and core mass (deBoer+ 2017, RvMP 89).

All rates return erg g^-1 s^-1.  T is in K, rho in g/cm^3.
"""
from __future__ import annotations

import numpy as np

MEV_ERG = 1.602176634e-6


def _T9(T):
    return T / 1.0e9


def eps_pp(rho, T, X, rates="reaclib"):
    """pp-chain (KWW12 eq. 18.63), with the psi factor ~1 (pp-I dominant)."""
    T9 = _T9(T)
    T7 = T / 1.0e7
    # KWW12 18.63: eps_pp ~ 2.57e4 psi f_11 g_11 rho X^2 T9^(-2/3) exp(-3.381/T9^(1/3))
    g11 = (1.0 + 0.0123 * T9 ** (1.0 / 3.0) + 0.0109 * T9 ** (2.0 / 3.0) + 0.0009 * T9)
    eps = 2.57e4 * 1.0 * 1.0 * g11 * rho * X * X * T9 ** (-2.0 / 3.0) * np.exp(-3.381 / T9 ** (1.0 / 3.0))
    return eps


def eps_cno(rho, T, X, X_cno, rates="reaclib"):
    """CNO cycle (KWW12 eq. 18.65).  X_cno = mass fraction of C+N+O catalysts."""
    T9 = _T9(T)
    g141 = (1.0 - 2.00 * T9 + 3.41 * T9 ** 2 - 2.43 * T9 ** 3)
    g141 = max(g141, 0.1)
    eps = 8.24e25 * g141 * X_cno * X * rho * T9 ** (-2.0 / 3.0) * np.exp(
        -15.231 * T9 ** (-1.0 / 3.0) - (T9 / 0.8) ** 2)
    if rates == "reaclib":
        eps *= 1.05
    return eps


def eps_triple_alpha(rho, T, Y, rates="reaclib"):
    """Triple-alpha (KWW12 eq. 18.67)."""
    T8 = T / 1.0e8
    f = 1.0  # screening ~1
    eps = 5.09e11 * f * rho ** 2 * Y ** 3 * T8 ** -3 * np.exp(-44.027 / T8)
    if rates == "reaclib":
        eps *= 0.97
    return eps


def eps_c12_alpha(rho, T, X_c12, X_he, c12_ag_rate=1.0, rates="reaclib"):
    """12C(alpha,gamma)16O energy generation (Caughlan & Fowler 1988 fit).

    Scaled by the free parameter `c12_ag_rate` (multiplier on the astrophysical
    S-factor).  Q = 7.162 MeV.
    """
    T9 = _T9(T)
    # CF88 rate [cm^3 s^-1 mol^-1], simplified leading term
    lam = 1.04e8 / T9**2 * np.exp(-32.120 * T9 ** (-1.0 / 3.0) - (T9 / 3.496) ** 2) \
        / (1.0 + 0.0489 * T9 ** (-2.0 / 3.0)) ** 2
    lam *= c12_ag_rate
    n_c12 = X_c12 * rho / 12.0
    n_he = X_he * rho / 4.0
    Q = 7.162 * MEV_ERG
    rate_per_vol = lam * n_c12 * n_he  # reactions cm^-3 s^-1 (mol absorbed into lam*N_A? keep approx)
    return rate_per_vol * Q / rho


def total_eps_nuc(rho, T, comp: dict, cfg_rates="reaclib", c12_ag_rate=1.0):
    """Sum the active channels given a composition dict with keys X, Y, Z, X_cno,
    X_c12 (optional).  Returns (eps_total, breakdown_dict)."""
    X = comp.get("X", 0.0)
    Y = comp.get("Y", 0.0)
    Z = comp.get("Z", 0.0)
    X_cno = comp.get("X_cno", 0.7 * Z)
    X_c12 = comp.get("X_c12", 0.0)
    b = {}
    b["pp"] = eps_pp(rho, T, X, cfg_rates) if X > 1e-6 else 0.0
    b["cno"] = eps_cno(rho, T, X, X_cno, cfg_rates) if X > 1e-6 and X_cno > 0 else 0.0
    b["3a"] = eps_triple_alpha(rho, T, Y, cfg_rates) if Y > 1e-3 else 0.0
    b["c12ag"] = eps_c12_alpha(rho, T, X_c12, Y, c12_ag_rate, cfg_rates) if X_c12 > 1e-4 and Y > 1e-3 else 0.0
    return sum(b.values()), b


def eps_neutrino_thermal(rho, T):
    """Crude thermal-neutrino losses (pair + plasma), Itoh+ 1996 order-of-magnitude
    scaling  eps_nu ~ 3e-4 T9^9 (pair, non-degenerate)  -- APPROXIMATE, KNOWN_GAPS NUC-NU.
    Only relevant T > 5e8 K.
    """
    T9 = T / 1.0e9
    if T9 < 0.3:
        return 0.0
    return 3.0e-4 * T9 ** 9 + 1.0e-6 * rho * T9 ** 6
