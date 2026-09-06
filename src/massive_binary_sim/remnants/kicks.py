"""Supernova natal kicks and the Blaauw mass-loss recoil, and their effect on the orbit.

Blaauw (1961), BAN 15, 265        -- symmetric instantaneous mass loss.
Hobbs et al. (2005), MNRAS 360, 974 -- single Maxwellian, sigma = 265 km/s (NS).
Fallback-scaled kick  v_kick = v_maxwell * (1 - f_fb)  (Fryer+ 2012) -- BHs formed
with large fallback get little or no kick.

Orbit update: Brandt & Podsiadlowski (1995), MNRAS 274, 461; Kalogera (1996),
ApJ 471, 352 -- instantaneous change of a, e (and system CM velocity) given the
pre-SN circular/eccentric state and a kick vector applied to the collapsing star.
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import G, KMS, MSUN


def draw_kick(rng, sigma_kms=265.0, f_fallback=0.0, scale_by_fallback=True):
    """Draw a 3-D kick vector [cm/s].  Maxwellian speed, isotropic direction."""
    v = rng.normal(0.0, sigma_kms, size=3) * KMS
    speed = np.linalg.norm(v)
    if scale_by_fallback:
        speed *= (1.0 - f_fallback)
    direction = v / (np.linalg.norm(v) + 1e-30)
    return speed * direction


def apply_sn_to_orbit(a_i, e_i, m1_pre, m2, m1_post, v_kick, rng, mean_anomaly=None):
    """Instantaneous SN of star 1.  Returns dict with a_f, e_f, bound flag,
    v_system (CM recoil, cm/s), and the pre-SN orbital speed used.

    Uses the standard vis-viva + kick composition at a random orbital phase
    (Kalogera 1996 eqs. 1-5).  For an initially eccentric orbit we sample the
    separation r from the eccentric anomaly.
    """
    if mean_anomaly is None:
        M_anom = rng.uniform(0, 2 * np.pi)
    else:
        M_anom = mean_anomaly
    # solve Kepler for eccentric anomaly
    E = M_anom
    for _ in range(50):
        E = E - (E - e_i * np.sin(E) - M_anom) / (1 - e_i * np.cos(E))
    r = a_i * (1 - e_i * np.cos(E))
    mu_i = G * (m1_pre + m2)
    v_rel = np.sqrt(mu_i * (2.0 / r - 1.0 / a_i))

    # place relative velocity in-plane; decompose into radial/tangential
    nu = 2.0 * np.arctan2(np.sqrt(1 + e_i) * np.sin(E / 2), np.sqrt(1 - e_i) * np.cos(E / 2))
    # flight-path angle
    vr = np.sqrt(mu_i / (a_i * (1 - e_i**2))) * e_i * np.sin(nu)
    vt = np.sqrt(mu_i / (a_i * (1 - e_i**2))) * (1 + e_i * np.cos(nu))
    v_vec = np.array([vr, vt, 0.0])

    v_new = v_vec + v_kick
    mu_f = G * (m1_post + m2)
    v2 = float(np.dot(v_new, v_new))
    energy = 0.5 * v2 - mu_f / r
    if energy >= 0:
        return dict(a_f=np.inf, e_f=np.nan, bound=False, v_system=np.nan,
                    v_rel_pre=v_rel, note="orbit unbound by SN")
    a_f = -mu_f / (2.0 * energy)
    h_vec = np.cross(np.array([r, 0.0, 0.0]), v_new)
    h = np.linalg.norm(h_vec)
    e_f = np.sqrt(max(0.0, 1.0 - h**2 / (mu_f * a_f)))

    # CM recoil (Blaauw + kick): dP = dM * v_1_orb_pre + m1_post * v_kick, over Mtot_post
    m1_frac = m2 / (m1_pre + m2)
    v1_orb = m1_frac * v_vec           # star-1 velocity about CM pre-SN
    dM = m1_pre - m1_post
    v_system = np.linalg.norm(dM * v1_orb + m1_post * v_kick) / (m1_post + m2)

    return dict(a_f=float(a_f), e_f=float(e_f), bound=True, v_system=float(v_system),
                v_rel_pre=float(v_rel), kick_speed=float(np.linalg.norm(v_kick)),
                note="bound after SN")
