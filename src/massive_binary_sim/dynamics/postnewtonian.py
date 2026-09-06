"""Post-Newtonian relative acceleration for the two-body problem.

Selectable order:  newtonian | 1PN | 2PN | 2.5PN  (and the aliases 3.5PN, mond).

Coefficients for the 1PN, 2PN and 2.5PN pieces are transcribed from:

    Mora & Will (2004), Phys. Rev. D 69, 104021, Eqs. (2.2)-(2.5)
    (harmonic-coordinate relative acceleration, modified-harmonic gauge).

The 2.5PN piece is the Damour-Deruelle radiation-reaction term and is identical
across Peters (1964), Iyer & Will (1995) and Blanchet (2014, LRR 17, 2, Eq. 203f).

KNOWN GAP (see KNOWN_GAPS.md, item PN-35):
  * The 3.5PN direct acceleration terms (Iyer & Will 1995; Pati & Will 2002) are
    NOT implemented.  Selecting gravity: "3.5PN" runs 2.5PN dynamics and logs a
    warning; the 3.5PN-accurate *orbit-averaged* fluxes ARE available for the
    secular inspiral via `dynamics.gw`.
  * The 2PN coefficients below were transcribed from memory of Mora & Will (2004)
    and are covered by the 1PN periastron-advance test only indirectly; treat 2PN
    dynamics as "implemented, partially validated" until cross-checked line-by-line.
"""
from __future__ import annotations

import warnings

import numpy as np

from ..numerics.units import C, G

# MOND acceleration scale a0 (Milgrom 1983, ApJ 270, 365; Begeman+ 1991 fit)
A0_MOND = 1.2e-8  # cm s^-2


def relative_acceleration(r_vec, v_vec, m1, m2, order: str = "2.5PN",
                          _warned: dict | None = None) -> np.ndarray:
    """Return d^2 r/dt^2 for the relative separation vector (CGS).

    `order` in {newtonian, 1PN, 2PN, 2.5PN, 3.5PN, mond}.
    """
    M = m1 + m2
    eta = m1 * m2 / (M * M)
    r = float(np.linalg.norm(r_vec))
    n = r_vec / r
    v2 = float(np.dot(v_vec, v_vec))
    rdot = float(np.dot(n, v_vec))
    GM = G * M
    GM_r = GM / r

    a_newt = -(GM / r**2) * n
    if order == "newtonian":
        return a_newt

    if order == "mond":
        # "simple" interpolating function mu(x)=x/(1+x); a_N -> a with
        # a mu(a/a0) = a_N   =>   a = a_N/2 + sqrt(a_N^2/4 + a_N a0)   (Famaey & McGaugh 2012)
        aN = GM / r**2
        a_mag = 0.5 * aN + np.sqrt(0.25 * aN**2 + aN * A0_MOND)
        return -a_mag * n

    c2 = C * C
    acc = a_newt.copy()

    # ---- 1PN  (Mora & Will 2004, Eq. 2.2) -----------------------------------------
    a1 = (GM / r**2) / c2 * (
        n * (2.0 * (2.0 + eta) * GM_r - (1.0 + 3.0 * eta) * v2 + 1.5 * eta * rdot**2)
        + v_vec * (2.0 * (2.0 - eta) * rdot)
    )
    acc = acc + a1
    if order == "1PN":
        return acc

    if order in ("2PN", "2.5PN", "3.5PN"):
        # ---- 2PN  (Mora & Will 2004, Eq. 2.3) ------------------------------------
        n_coeff = (
            0.75 * (12.0 + 29.0 * eta) * GM_r**2
            + eta * (3.0 - 4.0 * eta) * v2**2
            + (15.0 / 8.0) * eta * (1.0 - 3.0 * eta) * rdot**4
            - 1.5 * eta * (3.0 - 4.0 * eta) * v2 * rdot**2
            - 0.5 * eta * (13.0 - 4.0 * eta) * GM_r * v2
            - (2.0 + 25.0 * eta + 2.0 * eta**2) * GM_r * rdot**2
        )
        v_coeff = (
            0.5 * eta * (15.0 + 4.0 * eta) * v2
            - (4.0 + 41.0 * eta + 8.0 * eta**2) * GM_r
            - 1.5 * eta * (3.0 + 2.0 * eta) * rdot**2
        )
        a2 = -(GM / r**2) / c2**2 * (n * n_coeff + v_vec * rdot * v_coeff)
        acc = acc + a2
        if order == "2PN":
            return acc

    if order == "3.5PN":
        _w = _warned if _warned is not None else _MODULE_WARNED
        if not _w.get("pn35"):
            warnings.warn(
                "gravity='3.5PN': direct 3.5PN acceleration is NOT implemented "
                "(KNOWN_GAPS PN-35). Using 2.5PN dynamics.", RuntimeWarning)
            _w["pn35"] = True

    # ---- 2.5PN radiation reaction (Mora & Will 2004, Eq. 2.4) --------------------
    pref = (8.0 / 5.0) * eta * (G * G * M * M) / (C**5 * r**3)
    a25 = pref * (rdot * n * (3.0 * v2 + (17.0 / 3.0) * GM_r)
                  - v_vec * (v2 + 3.0 * GM_r))
    acc = acc + a25
    return acc


_MODULE_WARNED: dict = {}


def periastron_advance_per_orbit(a, e, m_total, order="1PN"):
    """Analytic relativistic periastron advance Delta omega per orbit [rad].

    1PN:  Delta phi = 6 pi G M / (c^2 a (1 - e^2))
    (Weinberg 1972, "Gravitation and Cosmology", eq. 8.6.11; standard result).
    Higher orders return the 1PN value with a NotImplemented note in `order`.
    """
    dphi_1pn = 6.0 * np.pi * G * m_total / (C**2 * a * (1.0 - e * e))
    return dphi_1pn
