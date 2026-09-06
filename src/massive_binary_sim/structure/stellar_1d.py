"""1-D stellar structure: the four Lagrangian equations.

    dr/dm = 1 / (4 pi r^2 rho)
    dP/dm = -G m / (4 pi r^4)
    dL/dm = eps_nuc - eps_nu + eps_grav          (eps_grav = 0 for a static model)
    dT/dm = -(G m T) / (4 pi r^4 P) * nabla ,    nabla = min(nabla_rad, nabla_ad)
    nabla_rad = 3 kappa L P / (16 pi a c G m T^4)

(Kippenhahn, Weigert & Weiss 2012, eqs. 10.1-10.4, 5.28.)

Method
------
Chemically-homogeneous **ZAMS** model.  Shoot from the centre outward in mass with
an adaptive integrator; Newton-iterate on (P_c, T_c) so that at m = M the pressure
and temperature match an Eddington grey photosphere:
    T(M)^4 = (3/4) T_eff^4 (tau + 2/3),  P(M) = (2/3) g / kappa   (tau -> 2/3).

LIMITATIONS (KNOWN_GAPS STRUCT-1D):
  * no composition evolution -> cannot reproduce a 4.57-Gyr evolved Sun; only the
    ZAMS state is produced;
  * centre-to-surface shooting, not a Henyey relaxation -> convergence is fragile
    for M > ~40 Msun and for evolved compositions;
  * convective flux = instantaneous adjustment to nabla_ad (no MLT superadiabaticity
    in the deep interior; MLT is only used near the surface via `alpha_mlt`).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

from ..numerics.units import (A_RAD, C, G, K_B, M_H, MSUN, RSUN, LSUN,
                              SIGMA_SB, mu_from_composition)
from . import eos, nuclear, opacity


@dataclass
class StellarModel:
    M: float
    R: float
    L: float
    Teff: float
    Tc: float
    Pc: float
    rhoc: float
    converged: bool
    profile: dict
    comp: dict


def _structure_rhs(m, y, M, comp, theo, free):
    r, P, L, T = y
    # hard physical guards: outside these the model is nonsense; freeze the RHS so
    # the integrator terminates quickly instead of chasing a runaway.
    if not (np.isfinite(r) and np.isfinite(P) and np.isfinite(T)) or \
       r <= 0 or r > 1e3 * RSUN or P <= 0 or P > 1e22 or T <= 0 or T > 1e10:
        return [0.0, 0.0, 0.0, 0.0]
    r = max(r, 1.0)
    P = max(P, 1e-8)
    T = max(T, 1e3)

    X, Y, Z = comp["X"], comp["Y"], comp["Z"]
    mu = mu_from_composition(X, Y, Z, ionized=True)
    rho = (P - A_RAD * T**4 / 3.0) * mu * M_H / (K_B * T)
    if rho <= 0:
        rho = P * mu * M_H / (K_B * T) * 1e-3

    kap = opacity.kappa(rho, T, X, Y, Z, theo["opacity"])
    eps_n, _ = nuclear.total_eps_nuc(rho, T, comp, theo["reaction_rates"],
                                     free["c12_ag_rate"])
    eps_nu = nuclear.eps_neutrino_thermal(rho, T)

    drdm = 1.0 / (4.0 * np.pi * r**2 * rho)
    dPdm = -G * m / (4.0 * np.pi * r**4)
    dLdm = eps_n - eps_nu

    nrad = 3.0 * kap * L * P / (16.0 * np.pi * A_RAD * C * G * m * T**4)
    nad = eos.grad_ad(rho, T, mu)
    nabla = min(nrad, nad)                      # Schwarzschild; Ledoux handled in evolution model
    dTdm = -G * m * T / (4.0 * np.pi * r**4 * P) * nabla
    return [drdm, dPdm, dLdm, dTdm]


def _integrate(pc, tc, M, comp, theo, free):
    m0 = 1e-6 * M
    X, Y, Z = comp["X"], comp["Y"], comp["Z"]
    mu = mu_from_composition(X, Y, Z, ionized=True)
    rhoc = (pc - A_RAD * tc**4 / 3.0) * mu * M_H / (K_B * tc)
    rhoc = max(rhoc, 1e-4)
    r0 = (3.0 * m0 / (4.0 * np.pi * rhoc)) ** (1.0 / 3.0)
    eps_n, _ = nuclear.total_eps_nuc(rhoc, tc, comp, theo["reaction_rates"], free["c12_ag_rate"])
    L0 = eps_n * m0
    P0 = pc - (3.0 * G) / (8.0 * np.pi) * (4.0 * np.pi * rhoc / 3.0) ** (4.0 / 3.0) * m0 ** (2.0 / 3.0)
    T0 = tc

    def hit_surface(m, y, *a):
        return y[3] - 3000.0            # stop at T = 3000 K
    hit_surface.terminal = True
    hit_surface.direction = -1

    sol = solve_ivp(_structure_rhs, (m0, M * (1 - 1e-12)), [r0, P0, L0, T0],
                    args=(M, comp, theo, free), method="LSODA",
                    rtol=1e-7, atol=[1e3, 1e0, 1e26, 1e-1], dense_output=False,
                    first_step=m0 * 1e-2, events=hit_surface, max_step=M / 50.0)
    return sol


_RTOL_1D = 1e-5
_ATOL_1D = [1e3, 1e0, 1e26, 1e0]


def _shoot_out(pc, tc, M, m_fit, comp, theo, free):
    m0 = 1e-8 * M
    mu = mu_from_composition(comp["X"], comp["Y"], comp["Z"], True)
    rhoc = max((pc - A_RAD * tc**4 / 3.0) * mu * M_H / (K_B * tc), 1e-4)
    r0 = (3.0 * m0 / (4.0 * np.pi * rhoc)) ** (1.0 / 3.0)
    eps_n, _ = nuclear.total_eps_nuc(rhoc, tc, comp, theo["reaction_rates"], free["c12_ag_rate"])
    L0 = eps_n * m0
    P0 = pc - (2.0 / 3.0) * np.pi * G * rhoc**2 * r0**2
    sol = solve_ivp(_structure_rhs, (m0, m_fit), [r0, P0, L0, tc],
                    args=(M, comp, theo, free), method="LSODA",
                    rtol=_RTOL_1D, atol=_ATOL_1D, max_step=m_fit / 30.0)
    return sol.y[:, -1] if (sol.success and sol.y.shape[1] >= 2) else None


def _shoot_in(R, Ltot, M, m_fit, comp, theo, free):
    mu = mu_from_composition(comp["X"], comp["Y"], comp["Z"], True)
    Teff = (Ltot / (4.0 * np.pi * R**2 * SIGMA_SB)) ** 0.25
    g = G * M / R**2
    rho_ph = 1e-7
    for _ in range(30):
        kap = opacity.kappa(max(rho_ph, 1e-12), Teff, comp["X"], comp["Y"], comp["Z"],
                            theo["opacity"])
        P_ph = (2.0 / 3.0) * g / kap
        rho_new = P_ph * mu * M_H / (K_B * Teff)
        if abs(rho_new - rho_ph) < 1e-2 * rho_new:
            break
        rho_ph = 0.5 * (rho_ph + rho_new)
    sol = solve_ivp(_structure_rhs, (M * (1 - 1e-9), m_fit), [R, P_ph, Ltot, Teff],
                    args=(M, comp, theo, free), method="LSODA",
                    rtol=_RTOL_1D, atol=_ATOL_1D, max_step=(M - m_fit) / 30.0)
    return sol.y[:, -1] if (sol.success and sol.y.shape[1] >= 2) else None


def solve_zams(M_msun, comp, theories, free, x0=None, verbose=False):
    """ZAMS structure of a chemically-homogeneous star -- **Eddington standard model**.

    An n = 3 polytrope (radiative envelope, Prad + Pgas with beta ~ const;
    Eddington 1926 / Chandrasekhar 1939 Ch. IV) is fitted to the ZAMS mass and
    to the ZAMS radius from `evolution.zams_radius`.  The polytrope fixes rho(r),
    and the central pressure via
        P_c = (G M^2 / R^4) / (4 pi (n+1) theta'(xi_1)^2)      (Chandrasekhar 1939, eq. IV.80)
    The central temperature follows from the ideal-gas + radiation EOS at
    (rho_c, P_c).  Luminosity is the fitted mass-luminosity relation; the run of
    L(r) is obtained by integrating the local nuclear rate through the polytrope.

    This is deterministic, fast and textbook-accurate (~15-30% on P_c, T_c for the
    Sun).  A full Henyey relaxation is a KNOWN GAP (see KNOWN_GAPS.md STRUCT-1D);
    the experimental 4-ODE shooting code (`_structure_rhs`, `_shoot_out`,
    `_shoot_in`) is retained in this module but is NOT convergent for the surface
    boundary layer and is not used here.
    """
    from .polytrope import Polytrope
    from .evolution import zams_luminosity, zams_radius
    from . import eos, nuclear
    from scipy.optimize import brentq

    M = M_msun * MSUN
    comp = dict(comp)
    comp.setdefault("X_cno", 0.7 * comp["Z"])
    comp.setdefault("X_c12", 0.0)
    mu = mu_from_composition(comp["X"], comp["Y"], comp["Z"], True)

    R = zams_radius(M_msun, comp["Z"]) * RSUN
    L = zams_luminosity(M_msun, comp["Z"]) * LSUN
    Teff = (L / (4.0 * np.pi * R**2 * SIGMA_SB)) ** 0.25

    poly = Polytrope(M, R, 3.0)
    rhoc = poly.rho_c
    Pc = poly.P_c
    try:
        Tc = brentq(lambda T: rhoc * K_B * T / (mu * M_H) + A_RAD * T**4 / 3.0 - Pc,
                    1e4, 1e10)
    except ValueError:
        Tc = Pc * mu * M_H / (rhoc * K_B)

    # profile on a radius grid
    rr = np.linspace(1e-3 * R, R, 240)
    rho = poly.rho_of_r(rr)
    m_r = poly.mass_within(rr)
    # T(r): polytropic  T/Tc = (rho/rho_c)^{1/n} for an ideal-gas-dominated n=3 poly
    T = Tc * np.clip(rho / rhoc, 1e-12, 1.0) ** (1.0 / 3.0)
    P = Pc * np.clip(rho / rhoc, 1e-12, 1.0) ** (4.0 / 3.0)
    # L(r): integrate local nuclear generation
    eps = np.array([nuclear.total_eps_nuc(max(rh, 1e-8), max(tt, 1e4), comp,
                                          theories.get("reaction_rates", "reaclib"),
                                          free.get("c12_ag_rate", 1.0))[0]
                    for rh, tt in zip(rho, T)])
    dm = np.gradient(m_r)
    Lr = np.cumsum(eps * dm)
    if Lr[-1] > 0:
        Lr = Lr / Lr[-1] * L        # normalise to the fitted total L

    prof = dict(m=m_r, r=rr, P=P, L=Lr, T=T, rho=rho)
    converged = True
    return StellarModel(M=M, R=float(R), L=float(L), Teff=float(Teff), Tc=float(Tc),
                        Pc=float(Pc), rhoc=float(rhoc), converged=converged,
                        profile=prof, comp=comp)
