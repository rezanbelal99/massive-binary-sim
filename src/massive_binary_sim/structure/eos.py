"""Equation of state.

Levels (selectable via `theories.opacity` is separate; EOS level follows
structure fidelity and the `ideal`/`ion`/`deg` flag passed in):

  ideal : ideal gas + radiation pressure
          P = rho k T / (mu m_H) + (1/3) a T^4          (KWW12, eqs. 13.2, 13.19)
  ion   : + Saha partial ionisation of H and He (KWW12, sec. 14.1) -> effective mu(T,rho)
  deg   : + partial non-relativistic electron degeneracy via the
          Chandrasekhar interpolation using Fermi-Dirac integrals
          (approximate; flagged in KNOWN_GAPS EOS-DEG)

Only H, He and a mean "metal" are tracked.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ..numerics.units import A_RAD, K_B, M_H, M_E, HBAR, C

CHI_H = 13.598 * 1.602176634e-12        # H ionisation energy [erg]
CHI_HE1 = 24.587 * 1.602176634e-12      # He I
CHI_HE2 = 54.418 * 1.602176634e-12      # He II


def pressure_ideal(rho, T, mu):
    """Gas + radiation pressure."""
    return rho * K_B * T / (mu * M_H) + A_RAD * T**4 / 3.0


def beta_gas(rho, T, mu):
    """P_gas / P_total."""
    pg = rho * K_B * T / (mu * M_H)
    pr = A_RAD * T**4 / 3.0
    return pg / (pg + pr)


def internal_energy_ideal(rho, T, mu):
    """Specific internal energy: (3/2) kT/(mu m_H) + a T^4 / rho  (KWW12 eq. 13.21)."""
    return 1.5 * K_B * T / (mu * M_H) + A_RAD * T**4 / rho


def adiabatic_gamma(rho, T, mu):
    """Generalised adiabatic exponent Gamma_1 for a gas+radiation mixture.

    Chandrasekhar (1939), sec. 56;  KWW12 eq. 13.20:
        Gamma_1 = beta + (4 - 3 beta)^2 (gamma_g - 1) / (beta + 12(gamma_g - 1)(1-beta))
    with gamma_g = 5/3 for a monatomic ideal gas.
    """
    b = beta_gas(rho, T, mu)
    g = 5.0 / 3.0
    return b + (4.0 - 3.0 * b) ** 2 * (g - 1.0) / (b + 12.0 * (g - 1.0) * (1.0 - b))


def grad_ad(rho, T, mu):
    """Adiabatic temperature gradient nabla_ad = (dlnT/dlnP)_s.

    KWW12 eq. 13.24:  nabla_ad = (4 - 3 beta) / (beta * (32 - 24 beta - 3 beta^2) / (4-3beta) + ... )
    We use the compact standard form
        nabla_ad = (1 + (1-beta)(4+beta)/beta^2 * (5/2 + ...)) ^ -1  ... use KWW12 13.24 directly.
    """
    b = beta_gas(rho, T, mu)
    # KWW12 (13.24): nabla_ad = (4 - 3 beta) / (32 - 24 beta - 3 beta^2)
    return (4.0 - 3.0 * b) / (32.0 - 24.0 * b - 3.0 * b * b)


# --------------------------------------------------------------------------------------
# Partial ionisation (Saha) -> effective mean molecular weight
# --------------------------------------------------------------------------------------
def _saha_ratio(chi, T, n_e):
    """n_{i+1} n_e / n_i  from the Saha equation (KWW12 eq. 14.7), g-factors ~1."""
    lam3 = (2.0 * np.pi * M_E * K_B * T / (2.0 * np.pi * HBAR) ** 2) ** 1.5
    return 2.0 * lam3 / n_e * np.exp(-chi / (K_B * T))


def mu_ionised(rho, T, X, Y, Z, max_iter=60):
    """Effective mean molecular weight including partial H/He ionisation.

    Fixed-point iteration on the electron density.  Metals assumed fully ionised
    contributing <Z/A> ~ 0.5 electrons per nucleon (approx; KNOWN_GAPS EOS-ION).
    """
    n_H = X * rho / M_H
    n_He = Y * rho / (4.0 * M_H)
    n_e = 0.5 * (n_H + 2 * n_He) + 0.5 * Z * rho / M_H + 1e-30
    for _ in range(max_iter):
        rH = _saha_ratio(CHI_H, T, n_e)
        yH = rH / (1.0 + rH)                       # H ionisation fraction
        r1 = _saha_ratio(CHI_HE1, T, n_e)
        r2 = _saha_ratio(CHI_HE2, T, n_e)
        # He: fractions neutral:once:twice
        denom = 1.0 + r1 + r1 * r2
        f1 = r1 / denom
        f2 = r1 * r2 / denom
        n_e_new = yH * n_H + (f1 + 2.0 * f2) * n_He + 0.5 * Z * rho / M_H
        n_e_new = max(n_e_new, 1e-30)
        if abs(n_e_new - n_e) < 1e-6 * n_e_new:
            n_e = n_e_new
            break
        n_e = 0.5 * (n_e + n_e_new)
    n_ion = n_H + n_He + Z * rho / M_H
    n_tot = n_ion + n_e
    return rho / (M_H * n_tot)


# --------------------------------------------------------------------------------------
# Partial electron degeneracy (approximate)
# --------------------------------------------------------------------------------------
def electron_degeneracy_pressure(rho, T, mu_e):
    """Non-relativistic partially-degenerate electron pressure via a smooth bridge
    between the ideal and the fully-degenerate limits (Paczynski 1983 style).

    P_e = sqrt( P_ideal^2 + P_deg^2 )    -- APPROXIMATE, see KNOWN_GAPS EOS-DEG.
    """
    n_e = rho / (mu_e * M_H)
    P_ideal = n_e * K_B * T
    K_nr = (3.0 / np.pi) ** (2.0 / 3.0) * (2.0 * np.pi * HBAR) ** 2 / (20.0 * M_E)
    P_deg = K_nr * n_e ** (5.0 / 3.0)
    return np.sqrt(P_ideal**2 + P_deg**2)
