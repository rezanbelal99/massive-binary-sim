"""Common-envelope evolution: alpha-lambda energy formalism and gamma AM formalism.

alpha-lambda  (Webbink 1984, ApJ 277, 355; de Kool 1990, ApJ 358, 189):
    E_bind = -G M_donor M_env / (lambda R_donor)
           = alpha_CE [ G M_core M_2 / (2 a_f)  -  G M_donor M_2 / (2 a_i) ]
    -> solve for a_f.  Merger if a_f implies the stripped core + companion would
       already be in contact (R_core + R_2 > a_f) or a_f <= 0.

gamma  (Nelemans & Tout 2005, MNRAS 356, 753):
    (J_i - J_f) / J_i = gamma * dM / (M_donor + M_2)
    -> a_f / a_i from angular-momentum balance.

lambda:  either the config constant `lambda_ce`, or a fitted structural value from
Xu & Li (2010), ApJ 716, 114 / Klencki et al. (2021), A&A 645, A54 (envelope
binding energy including internal energy).  Fitted mode is APPROXIMATE
(KNOWN_GAPS CE-LAMBDA).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..numerics.units import G, MSUN, RSUN


@dataclass
class CEResult:
    survived: bool
    a_final: float           # cm (nan if merger)
    a_initial: float
    m_core_donor: float      # g
    m_env_removed: float     # g
    lambda_used: float
    verdict: str


def lambda_fit(donor, include_internal_energy=True):
    """Fitted envelope-structure parameter.

    Very rough: lambda rises from ~0.5 (compact BSG) to ~1.5 (extended RSG with
    recombination energy).  Klencki+ 2021 show lambda can reach ~2-20 for cool
    RSGs when internal energy is credited; we cap at 5 and flag it.
    """
    import numpy as np
    R = donor.R / RSUN
    base = 0.5 + 0.4 * np.log10(max(R, 1.0))
    if include_internal_energy and getattr(donor, "phase", "") == "CHeB" and R > 200:
        base *= 3.0
    return float(np.clip(base, 0.05, 5.0))


def _effective_donor_radius(donor):
    """Radius to use for the envelope binding energy.

    A star that initiates CE while crossing the Hertzsprung gap will, on the
    thermal timescale, be a giant by the time the spiral-in completes; the
    envelope binding energy is set by that (larger) radius, not the HG radius.
    We use the star's own CHeB radius estimate when it is post-MS.
    (Klencki et al. 2021 emphasise this radius sensitivity.)
    """
    from ..numerics.units import MSUN, RSUN
    R = donor.R
    if getattr(donor, "phase", "") in ("HG", "CHeB") and not getattr(donor, "stripped", False):
        M0 = donor.m0 / MSUN
        R_rsg = (700.0 * (M0 / 20.0) ** 0.5 if M0 < 40 else 60.0 * (M0 / 40.0) ** 0.4) * RSUN
        R = max(R, R_rsg)
    return R


def binding_energy(donor, lam):
    """E_bind = -G M_donor M_env / (lambda R_donor)  (negative)."""
    m_env = donor.m - (donor.m_he_core if donor.m_he_core > 0 else 0.5 * donor.m)
    m_env = max(m_env, 0.0)
    R_eff = _effective_donor_radius(donor)
    return -G * donor.m * m_env / (lam * R_eff), m_env


def alpha_lambda_outcome(donor, m_companion, a_i, alpha_ce, lambda_ce=None,
                         r_companion=None):
    m_core = donor.m_he_core if donor.m_he_core > 0 else 0.5 * donor.m
    lam = lambda_ce if lambda_ce is not None else lambda_fit(donor)
    Ebind, m_env = binding_energy(donor, lam)  # Ebind < 0

    # alpha_CE ( -G Mc M2/(2 af) + G Md M2/(2 ai) ) = Ebind  (Ebind<0)
    # => G Mc M2 / (2 af) = G Md M2/(2 ai) - Ebind/alpha
    rhs = G * m_core * m_companion / 2.0
    denom = G * donor.m * m_companion / (2.0 * a_i) - Ebind / alpha_ce
    if denom <= 0:
        return CEResult(False, np.nan, a_i, m_core, m_env, lam,
                        "merger: envelope binding energy exceeds available orbital energy")
    a_f = rhs / denom
    # companion radius: a stripped-core / MS star if given, else a compact-object
    # scale (~0.01 Rsun stand-in so the periastron check still has a floor)
    r_c = r_companion if (r_companion is not None and r_companion > 0) else 0.02 * RSUN
    r_core = 0.2 * (m_core / MSUN) ** 0.6 * RSUN
    # require the post-CE orbit to clear periastron contact with a factor-1.5 margin
    if a_f <= 1.5 * (r_core + r_c):
        return CEResult(False, a_f, a_i, m_core, m_env, lam,
                        f"merger: a_f={a_f/RSUN:.3f} Rsun < 1.5(R_core+R_comp)")
    return CEResult(True, a_f, a_i, m_core, m_env, lam,
                    f"survived: a_f={a_f/RSUN:.3f} Rsun (spiral-in x{a_i/a_f:.1f})")


def gamma_outcome(donor, m_companion, a_i, gamma):
    m_core = donor.m_he_core if donor.m_he_core > 0 else 0.5 * donor.m
    dM = donor.m - m_core
    Mi = donor.m + m_companion
    Mf = m_core + m_companion
    # (Ji - Jf)/Ji = gamma dM / Mi ; J ~ M1 M2 sqrt(a/Mtot) (circular)
    Ji_factor = donor.m * m_companion * np.sqrt(a_i / Mi)
    Jf_target = Ji_factor * (1.0 - gamma * dM / Mi)
    if Jf_target <= 0:
        return CEResult(False, np.nan, a_i, m_core, dM, np.nan,
                        "merger: gamma-formalism removes all orbital AM")
    a_f = Mf * (Jf_target / (m_core * m_companion)) ** 2
    r_core = 0.2 * (m_core / MSUN) ** 0.6 * RSUN
    r_c = 0.1 * m_companion / MSUN * RSUN
    if a_f <= (r_core + r_c):
        return CEResult(False, a_f, a_i, m_core, dM, np.nan,
                        f"merger: a_f={a_f/RSUN:.2f} Rsun too small")
    return CEResult(True, a_f, a_i, m_core, dM, np.nan,
                    f"survived (gamma): a_f={a_f/RSUN:.3f} Rsun (shrink x{a_i/a_f:.1f})")


def common_envelope(donor, m_companion, a_i, formalism="alpha_lambda",
                    alpha_ce=1.0, lambda_ce=0.5, gamma=1.3, r_companion=None):
    if formalism == "gamma":
        return gamma_outcome(donor, m_companion, a_i, gamma)
    return alpha_lambda_outcome(donor, m_companion, a_i, alpha_ce, lambda_ce, r_companion)
