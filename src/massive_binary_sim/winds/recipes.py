"""Stellar-wind mass-loss recipes behind one interface.

    mass_loss_rate(star_state, recipe, eta_wind, clumping) -> Mdot [g/s]  (negative)

Recipes
-------
vink2001      : Vink, de Koter & Lamers (2001), A&A 369, 574  -- O/B stars,
                incl. the bi-stability jump at Teff ~ 25 kK; Z-scaling Mdot ~ Z^0.85
                (Vink+ 2001 eq. 25 uses ~0.85 exponent on the line-driven part).
bjorklund2021 : Bjorklund et al. (2021), A&A 648, A36  -- systematically lower
                (~factor 2-3) O-star rates from co-moving-frame models.
nugis2000     : Nugis & Lamers (2000), A&A 360, 227 -- Wolf-Rayet / stripped He stars.
dejager1988   : de Jager, Nieuwenhuijzen & van der Hucht (1988), A&AS 72, 259 --
                cool supergiants (empirical).
auto          : pick by evolutionary phase / Teff / surface composition.

The `clumping` factor rescales EMPIRICAL rates (Nugis, de Jager) downward by
1/sqrt(f_cl); the theoretical Vink/Bjorklund rates assume a smooth wind and are
left unscaled (Smith 2014, ARA&A 52, 487).
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import LSUN, MSUN, RSUN, YEAR, G

MSUN_PER_YR = MSUN / YEAR


def _logMdot_vink_hot(logL, M, Teff, Z, vinf_over_vesc):
    """Vink+ 2001 eq. (24): hot side of the bi-stability jump (Teff > ~25 kK)."""
    return (-6.697 + 2.194 * (logL - 5.0) - 1.313 * np.log10(M / 30.0)
            - 1.226 * np.log10(vinf_over_vesc / 2.0) + 0.933 * np.log10(Teff / 40000.0)
            - 10.92 * np.log10(Teff / 40000.0) ** 2 + 0.85 * np.log10(Z / 0.019))


def _logMdot_vink_cool(logL, M, Teff, Z, vinf_over_vesc):
    """Vink+ 2001 eq. (25): cool side (~12.5-25 kK)."""
    return (-6.688 + 2.210 * (logL - 5.0) - 1.339 * np.log10(M / 30.0)
            - 1.601 * np.log10(vinf_over_vesc / 2.0) + 1.07 * np.log10(Teff / 20000.0)
            + 0.85 * np.log10(Z / 0.019))


def vink2001(L, M, Teff, Z):
    logL = np.log10(L / LSUN)
    Mm = M / MSUN
    # v_inf / v_esc ratio (Lamers+ 1995): 2.6 hot, 1.3 cool, jump at 25 kK
    if Teff >= 25000.0:
        ratio = 2.6
        logMdot = _logMdot_vink_hot(logL, Mm, Teff, Z, ratio)
    elif Teff >= 12500.0:
        ratio = 1.3
        logMdot = _logMdot_vink_cool(logL, Mm, Teff, Z, ratio)
    else:
        return None  # out of Vink validity range
    return -(10.0 ** logMdot) * MSUN_PER_YR


def bjorklund2021(L, M, Teff, Z):
    """Bjorklund+ 2021 eq. (15) fit (their recommended relation)."""
    logL = np.log10(L / LSUN)
    Mm = M / MSUN
    logMdot = (-5.52 + 2.39 * (logL - 6.0) - 1.15 * np.log10(Mm / 45.0)
               + 0.85 * np.log10(Z / 0.014) + 0.42 * np.log10(Teff / 45000.0))
    return -(10.0 ** logMdot) * MSUN_PER_YR


def nugis2000(L, Y_surf, Z):
    """Nugis & Lamers (2000) eq. (22): WR rate vs L, surface He and Z."""
    logL = np.log10(L / LSUN)
    logMdot = -13.60 + 1.63 * logL + 2.22 * np.log10(max(Y_surf, 0.3)) + 0.85 * np.log10(Z / 0.019)
    return -(10.0 ** logMdot) * MSUN_PER_YR


def dejager1988(L, Teff):
    """de Jager+ 1988 -- compact 2-term truncation of their Chebyshev fit.

    log(-Mdot) = 1.769 logL - 1.676 logTeff - 8.158   (their eq. 5, leading terms).
    Valid for cool, luminous stars.
    """
    logL = np.log10(L / LSUN)
    logMdot = 1.769 * logL - 1.676 * np.log10(Teff) - 8.158
    return -(10.0 ** logMdot) * MSUN_PER_YR


MDOT_CEILING = 1.0e-2 * MSUN_PER_YR   # no massive-star wind exceeds ~1e-2 Msun/yr


def mass_loss_rate(star, recipe="auto", eta_wind=1.0, clumping=1.0):
    """Return Mdot [g/s] (<=0) for a StarState-like object with .L,.m,.Teff,.Z,.phase."""
    L, M, Teff, Z = star.L, star.m, star.Teff, star.Z
    if not (np.isfinite(L) and np.isfinite(Teff)) or L <= 0 or Teff <= 0:
        return 0.0
    Y_surf = 0.98 if getattr(star, "stripped", False) else 0.28
    empirical_scale = eta_wind / np.sqrt(max(clumping, 1.0))
    theo_scale = eta_wind

    if recipe == "auto":
        if getattr(star, "stripped", False):
            recipe = "nugis2000"
        elif Teff >= 12500.0:
            recipe = "vink2001"
        else:
            recipe = "dejager1988"

    if recipe == "vink2001":
        md = vink2001(L, M, Teff, Z)
        md = dejager1988(L, Teff) * empirical_scale if md is None else md * theo_scale
    elif recipe == "bjorklund2021":
        md = bjorklund2021(L, M, Teff, Z)
        # Bjorklund+ 2021 is calibrated for hot O stars; fall back to de Jager for
        # cool / evolved stars where its fit is extrapolating wildly
        if Teff < 12500.0 or not np.isfinite(md):
            md = dejager1988(L, Teff) * empirical_scale
        else:
            md = md * theo_scale
    elif recipe == "nugis2000":
        md = nugis2000(L, Y_surf, Z) * empirical_scale
    elif recipe == "dejager1988":
        md = dejager1988(L, Teff) * empirical_scale
    else:
        raise ValueError(f"unknown wind recipe {recipe!r}")

    if not np.isfinite(md):
        return 0.0
    return -min(abs(md), MDOT_CEILING)


def wind_angular_momentum_loss(star, mdot, a, m_other, mode="jeans"):
    """Orbital-AM loss from a fast (Jeans-mode) wind.

    In the Jeans / fast-wind approximation the wind carries the specific orbital
    AM of the mass-losing star:  Jdot/J = mdot_loss * m_other / (m_star (m_star+m_other))
    (Soberman, Phinney & van den Heuvel 1997, A&A 327, 620, eq. 25 with beta=1).
    Returns dJ/dt [cgs], negative.
    """
    m_star = star.m
    mtot = m_star + m_other
    mu_red = m_star * m_other / mtot
    J_orb = mu_red * np.sqrt(G * mtot * a)
    # fractional: dJ/J = (mdot/m_star) * (m_other/mtot)   for isotropic wind from m_star
    return J_orb * (mdot / m_star) * (m_other / mtot)
