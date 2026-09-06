"""Roche geometry: Eggleton radius, full dimensionless potential, Lagrange points.

References
---------
Eggleton, P. P. (1983), ApJ 268, 368 -- R_L/a fit.
Kopal (1959), "Close Binary Systems" -- Roche potential in the corotating frame.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from ..numerics.units import roche_lobe_eggleton  # noqa: F401


def roche_lobe_radius(a, m_lobe, m_other):
    """Volume-equivalent Roche-lobe radius [same units as a]. Eggleton (1983)."""
    q = m_lobe / m_other
    return a * roche_lobe_eggleton(q)


def roche_potential(x, y, z, q):
    """Dimensionless Roche potential in units where a=1, G(M1+M2)=1, primary at
    origin (mass fraction 1/(1+q)) and secondary at (1,0,0) (mass fraction q/(1+q)),
    corotating frame.  q = M2/M1.

    Phi = -mu1/r1 - mu2/r2 - 1/2 [(x-x_cm)^2 + y^2] * (omega^2, =1 here)
    with mu1 = 1/(1+q), mu2 = q/(1+q), x_cm = q/(1+q).
    """
    mu1 = 1.0 / (1.0 + q)
    mu2 = q / (1.0 + q)
    x_cm = mu2
    r1 = np.sqrt(x**2 + y**2 + z**2)
    r2 = np.sqrt((x - 1.0) ** 2 + y**2 + z**2)
    return -mu1 / r1 - mu2 / r2 - 0.5 * ((x - x_cm) ** 2 + y**2)


def _dPhi_dx_on_axis(x, q):
    mu1 = 1.0 / (1.0 + q)
    mu2 = q / (1.0 + q)
    x_cm = mu2
    t1 = mu1 * np.sign(x) / x**2 if x != 0 else np.inf
    t2 = mu2 * np.sign(x - 1.0) / (x - 1.0) ** 2 if x != 1.0 else np.inf
    return t1 + t2 - (x - x_cm)


def lagrange_points(q):
    """Return dict L1..L5 as (x, y) in the a=1 frame (primary at 0, secondary at 1)."""
    f = lambda x: _dPhi_dx_on_axis(x, q)
    L1 = brentq(f, 1e-6, 1.0 - 1e-6)
    L2 = brentq(f, 1.0 + 1e-6, 3.0)
    L3 = brentq(f, -3.0, -1e-6)
    x_cm = q / (1.0 + q)
    L4 = (0.5 - x_cm + 0.5, np.sqrt(3.0) / 2.0)  # equilateral: x = 0.5 (midpoint), offset
    # equilateral points are at x = 1/2 (relative to the two masses) ... in our frame
    L4 = (0.5, np.sqrt(3.0) / 2.0)
    L5 = (0.5, -np.sqrt(3.0) / 2.0)
    return dict(L1=(L1, 0.0), L2=(L2, 0.0), L3=(L3, 0.0), L4=L4, L5=L5)


def roche_lobe_fill_factor(star_radius, a, m_lobe, m_other, eccentricity=0.0,
                           use_periastron=True):
    """R_star / R_L.  For eccentric orbits the lobe at periastron uses a(1-e)
    (Sepinsky, Willems & Kalogera 2007, ApJ 660, 1624 -- instantaneous-separation
    approximation; full non-synchronous treatment not implemented, KNOWN_GAPS BIN-ECC).
    """
    sep = a * (1.0 - eccentricity) if use_periastron else a
    return star_radius / roche_lobe_radius(sep, m_lobe, m_other)


def equipotential_grid(q, nx=400, ny=400, extent=1.8):
    """Return (X, Y, Phi) meshgrid for contour plotting in the orbital plane."""
    xs = np.linspace(-extent + 0.5, extent + 0.5, nx)
    ys = np.linspace(-extent, extent, ny)
    X, Y = np.meshgrid(xs, ys)
    Phi = roche_potential(X, Y, 0.0 * X, q)
    return X, Y, Phi
