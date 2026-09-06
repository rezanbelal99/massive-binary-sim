"""Kozai-Lidov secular evolution from a distant third body (optional, toggle).

Quadrupole-order, test-particle (outer orbit fixed) secular equations in the
eccentric Kozai-Lidov formalism.

Reference
---------
Naoz, S. (2016), ARA&A 54, 441  -- Eqs. (11)-(17), quadrupole test-particle limit.
Kozai (1962), AJ 67, 591 ; Lidov (1962), Planet. Space Sci. 9, 719.

We evolve (e_in, i_mut, omega_in) with the conserved quantity
    H = sqrt(1-e_in^2) cos(i_mut)   (Kozai integral, quadrupole test-particle).
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import G


def lk_timescale(a_in, a_out, e_out, m0, m1, m2):
    """Characteristic LK oscillation period (Naoz 2016, Eq. 27), seconds."""
    P_in = 2 * np.pi * np.sqrt(a_in**3 / (G * (m0 + m1)))
    return (8.0 / (15.0 * np.pi) * (m0 + m1 + m2) / m2
            * (a_out * np.sqrt(1 - e_out**2) / a_in) ** 3 * P_in)


def lk_derivs(t, y, a_in, a_out, e_out, m0, m1, m2):
    """d/dt of [e_in, i_mut(rad), omega_in(rad)] at quadrupole order.

    Naoz (2016) Eqs. (A)(quadrupole, test particle), with time normalised by the
    outer-orbit mean motion; here written dimensionally (per second).
    """
    e, inc, w = y
    e = min(max(e, 1e-8), 1 - 1e-8)
    n_out = np.sqrt(G * (m0 + m1 + m2) / a_out**3)
    eps = (m2 / (m0 + m1)) * (a_in / a_out) ** 3 / (1 - e_out**2) ** 1.5
    tau = n_out * eps  # inverse of ~ (3/4) LK time

    j = np.sqrt(1 - e**2)
    de = tau * (15.0 / 8.0) * e * j * np.sin(inc) ** 2 * np.sin(2 * w)
    di = -tau * (15.0 / 8.0) * (e**2 / j) * np.sin(inc) * np.cos(inc) * np.sin(2 * w)
    dw = tau * (3.0 / 4.0) * (1.0 / j) * (
        2.0 * (1 - e**2) + 5.0 * np.sin(w) ** 2 * (e**2 - np.sin(inc) ** 2))
    return [de, di, dw]
