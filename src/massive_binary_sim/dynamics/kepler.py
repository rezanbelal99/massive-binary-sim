"""Newtonian two-body kinematics: orbital elements <-> Cartesian state.

Reduced / relative description: we integrate the *relative* separation vector
r = r1 - r2 with the standard equation  d^2 r / dt^2 = -G M r / |r|^3  (M = m1+m2),
then reconstruct individual positions from the mass ratio.

References
----------
Murray & Dermott (1999), "Solar System Dynamics", Ch. 2  -- element definitions,
Kepler-equation solution.
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import G


def elements_to_state(a, e, inc, Omega, omega, nu, mu):
    """(a,e,i,Omega,omega,true-anomaly) -> relative position & velocity vectors.

    Angles in radians, `mu = G(m1+m2)` in CGS, `a` in cm.  Returns (r_vec, v_vec).
    Murray & Dermott (1999) eq. 2.122 / 2.126.
    """
    p = a * (1.0 - e * e)
    r = p / (1.0 + e * np.cos(nu))

    # perifocal frame
    r_pf = np.array([r * np.cos(nu), r * np.sin(nu), 0.0])
    v_pf = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])

    cO, sO = np.cos(Omega), np.sin(Omega)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(omega), np.sin(omega)
    R = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si, cw * si, ci],
    ])
    return R @ r_pf, R @ v_pf


def state_to_elements(r_vec, v_vec, mu):
    """Relative state -> classical elements.  Returns dict.

    Standard algorithm (Vallado 2013, "Fundamentals of Astrodynamics", Alg. 9).
    """
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    h_vec = np.cross(r_vec, v_vec)
    h = np.linalg.norm(h_vec)
    n_vec = np.cross([0, 0, 1.0], h_vec)
    n = np.linalg.norm(n_vec)

    e_vec = (np.cross(v_vec, h_vec) / mu) - r_vec / r
    e = np.linalg.norm(e_vec)

    energy = v * v / 2.0 - mu / r
    if abs(energy) > 1e-30:
        a = -mu / (2.0 * energy)
    else:
        a = np.inf

    inc = np.arccos(np.clip(h_vec[2] / h, -1.0, 1.0))

    if n > 1e-30:
        Omega = np.arccos(np.clip(n_vec[0] / n, -1.0, 1.0))
        if n_vec[1] < 0:
            Omega = 2 * np.pi - Omega
    else:
        Omega = 0.0

    if n > 1e-12 * h and e > 1e-12:
        omega = np.arccos(np.clip(np.dot(n_vec, e_vec) / (n * e), -1.0, 1.0))
        if e_vec[2] < 0:
            omega = 2 * np.pi - omega
    elif e > 1e-12:
        # equatorial orbit: node undefined -> use longitude of periapsis in the
        # x-y plane (retrograde-corrected). Adequate for tracking precession.
        omega = np.arctan2(e_vec[1], e_vec[0]) % (2 * np.pi)
        if h_vec[2] < 0:
            omega = (2 * np.pi - omega) % (2 * np.pi)
    else:
        omega = 0.0

    if e > 1e-12:
        nu = np.arccos(np.clip(np.dot(e_vec, r_vec) / (e * r), -1.0, 1.0))
        if np.dot(r_vec, v_vec) < 0:
            nu = 2 * np.pi - nu
    else:
        nu = 0.0

    return dict(a=a, e=e, inc=inc, Omega=Omega, omega=omega, nu=nu, h=h, energy=energy,
               period=(2 * np.pi * np.sqrt(a**3 / mu) if np.isfinite(a) and a > 0 else np.inf))


def split_to_bodies(r_rel, v_rel, m1, m2):
    """Relative -> individual body states in the centre-of-mass frame."""
    mt = m1 + m2
    r1 = (m2 / mt) * r_rel
    r2 = -(m1 / mt) * r_rel
    v1 = (m2 / mt) * v_rel
    v2 = -(m1 / mt) * v_rel
    return (r1, v1), (r2, v2)


def orbital_energy_angular_momentum(r_rel, v_rel, m1, m2):
    """Total (Newtonian) orbital energy and |L| in the CM frame, CGS."""
    mu_red = m1 * m2 / (m1 + m2)
    r = np.linalg.norm(r_rel)
    v2 = float(np.dot(v_rel, v_rel))
    E = 0.5 * mu_red * v2 - G * m1 * m2 / r
    L = mu_red * np.linalg.norm(np.cross(r_rel, v_rel))
    return E, L
