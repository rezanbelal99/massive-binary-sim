"""Core-collapse remnant masses.

Fryer et al. (2012), ApJ 749, 91 -- "rapid" and "delayed" explosion engines,
Eqs. (10)-(20).  Inputs: pre-SN total mass M and CO-core mass M_CO (Msun).
Returns remnant baryonic mass, then gravitational mass via the Lattimer & Yahil
(1989) / Timmes+ (1996) binding-energy correction
    M_grav = (sqrt(1 + 0.336 M_rem_bary) - 1) / 0.168      (for neutron stars)
    M_grav ~ 0.9 M_rem_bary                                (for black holes, approx)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

M_NS_MAX = 2.5   # Msun, maximum (gravitational) NS mass adopted (Fryer+ 2012 use ~2.0-3.0)


@dataclass
class RemnantResult:
    kind: str            # "NS" | "BH" | "PISN_none" | "PPISN_BH"
    m_baryonic: float    # Msun
    m_gravitational: float
    m_fallback: float    # Msun that fell back
    f_fallback: float    # fallback fraction of the ejected-then-returned envelope
    engine: str
    note: str = ""


def _proto_ns(M_CO):
    """Proto-compact-object mass M_proto (Fryer+ 2012 Eq. 15/18)."""
    return np.where(M_CO < 3.5, 1.2,
           np.where(M_CO < 6.0, 1.3,
           np.where(M_CO < 11.0, 1.4, 1.6)))


def fryer2012(M, M_CO, engine="delayed"):
    """Return RemnantResult. M, M_CO in Msun."""
    M_CO = float(M_CO)
    M = float(M)
    if engine == "direct_collapse":
        m_bar = M
        return _finish("direct_collapse", m_bar, M, 1.0, "full direct collapse (no SN)")

    if engine == "rapid":
        # Fryer+ 2012 Eqs. (15)-(17): proto-mass fixed at 1.0 Msun
        M_proto = 1.0
        if M_CO < 2.5:
            f_fb = 0.2 / max(M - M_proto, 1e-3)
        elif M_CO < 6.0:
            f_fb = 0.286 * M_CO - 0.514
        else:
            f_fb = 1.0
        f_fb = float(np.clip(f_fb, 0.0, 1.0))
        m_bar = M_proto + f_fb * (M - M_proto)
        return _finish("rapid", m_bar, M, f_fb, f"rapid engine, M_CO={M_CO:.2f}")

    # delayed engine (Fryer+ 2012 Eqs. 18-20)
    if M_CO < 2.5:
        M_proto = 1.2
        m_fb = 0.2
    elif M_CO < 3.5:
        M_proto = 1.3
        m_fb = 0.5 * M_CO - 1.05
    elif M_CO < 11.0:
        M_proto = 1.4
        a2 = 0.133 - 0.093 / (M - M_proto if M - M_proto > 0.1 else 0.1)
        b2 = -11.0 * a2 + 1.0
        f_fb = a2 * M + b2
        m_fb = float(np.clip(f_fb, 0.0, 1.0)) * (M - M_proto)
    else:
        M_proto = 1.6
        m_fb = M - M_proto
    f_fb = float(np.clip(m_fb / max(M - M_proto, 1e-3), 0.0, 1.0))
    m_bar = M_proto + f_fb * (M - M_proto)
    return _finish("delayed", m_bar, M, f_fb, f"delayed engine, M_CO={M_CO:.2f}")


def _finish(engine, m_bar, M, f_fb, note):
    # NS vs BH by gravitational mass
    m_grav_ns = (np.sqrt(1.0 + 0.336 * m_bar) - 1.0) / 0.168
    if m_grav_ns <= M_NS_MAX and m_bar <= 2.8:
        return RemnantResult("NS", m_bar, float(m_grav_ns), (M - m_bar), f_fb, engine, note)
    m_grav_bh = 0.9 * m_bar   # ~10% neutrino losses at BH formation (approx)
    return RemnantResult("BH", m_bar, float(m_grav_bh), (M - m_bar), f_fb, engine, note)
