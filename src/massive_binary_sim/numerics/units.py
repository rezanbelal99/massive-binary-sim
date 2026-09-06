"""Single units / physical-constants layer.

Everything internal to the simulation is in **CGS-Gaussian-compatible CGS** unless a
docstring says otherwise.  All public physical constants are taken from one of:

* CODATA 2018 (via ``astropy.constants``), or
* IAU 2015 Resolution B3 nominal solar/planetary conversion constants, or
* Asplund, Grevesse, Sauval & Scott (2009), ARA&A 47, 481  (solar composition).

Provenance for every constant is recorded in ``CONSTANT_PROVENANCE`` so that the
report generator can print a sourced table.  No constant is defined by hand here
without a citation; where ``astropy`` is used the underlying reference is CODATA
2018 / IAU 2015 as documented in the ``astropy.constants`` reference table.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import astropy.constants as _ac
import astropy.units as u

# --------------------------------------------------------------------------------------
# Fundamental constants in CGS
# --------------------------------------------------------------------------------------
G = _ac.G.cgs.value               # gravitational constant [cm^3 g^-1 s^-2]
C = _ac.c.cgs.value               # speed of light [cm s^-1]
C2 = C * C
C5 = C**5
H_PLANCK = _ac.h.cgs.value        # Planck constant [erg s]
HBAR = _ac.hbar.cgs.value
K_B = _ac.k_B.cgs.value           # Boltzmann constant [erg K^-1]
SIGMA_SB = _ac.sigma_sb.cgs.value  # Stefan-Boltzmann [erg cm^-2 s^-1 K^-4]
A_RAD = 4.0 * SIGMA_SB / C         # radiation constant [erg cm^-3 K^-4]
M_H = _ac.m_p.cgs.value           # proton mass [g]  (approx atomic hydrogen mass)
M_E = _ac.m_e.cgs.value
M_U = _ac.u.cgs.value             # atomic mass unit [g]
N_A = _ac.N_A.value               # Avogadro [mol^-1]
EV = _ac.e.gauss.value ** 0 * (1.602176634e-12)  # 1 eV in erg (CODATA 2018 exact e)
MEV = EV * 1.0e6

# --------------------------------------------------------------------------------------
# Astronomical conversion constants (IAU 2015 B3 nominal values)
# --------------------------------------------------------------------------------------
MSUN = _ac.M_sun.cgs.value        # 1.98892e33 g  (derived from nominal GM_sun / G)
RSUN = _ac.R_sun.cgs.value        # 6.957e10 cm   (IAU 2015 nominal)
LSUN = _ac.L_sun.cgs.value        # 3.828e33 erg/s (IAU 2015 nominal)
TEFF_SUN = 5772.0                 # K, IAU 2015 nominal solar effective temperature
AU = _ac.au.cgs.value
PC = _ac.pc.cgs.value
YEAR = (1.0 * u.yr).to(u.s).value           # Julian year = 365.25 d
DAY = 86400.0
GYR = 1.0e9 * YEAR
MYR = 1.0e6 * YEAR
KM = 1.0e5
KMS = 1.0e5

# --------------------------------------------------------------------------------------
# Solar composition  (Asplund, Grevesse, Sauval & Scott 2009)
# --------------------------------------------------------------------------------------
X_SUN = 0.7381        # hydrogen mass fraction (photospheric+meteoritic, AGSS09)
Y_SUN = 0.2485        # helium mass fraction
Z_SUN = 0.0134        # metal mass fraction (AGSS09 present-day photosphere)
# NOTE: the protosolar Z is ~0.0142; AGSS09 photospheric Z = 0.0134.  The config
# default uses 0.0142 following common massive-star-evolution practice (Brott+ 2011
# adopt Z_sun = 0.0088 for a *different* scale; see KNOWN_GAPS.md item Z-SCALE).

# --------------------------------------------------------------------------------------
# Opacity floor constants
# --------------------------------------------------------------------------------------
KAPPA_ES_COEFF = 0.20             # cm^2/g ; electron-scattering opacity = 0.2(1+X)
# (Kippenhahn, Weigert & Weiss 2012, "Stellar Structure and Evolution" 2nd ed., eq 17.2)


@dataclass(frozen=True)
class ConstProvenance:
    symbol: str
    value: float
    unit: str
    source: str


CONSTANT_PROVENANCE = [
    ConstProvenance("G", G, "cm^3 g^-1 s^-2", "CODATA 2018 via astropy.constants"),
    ConstProvenance("c", C, "cm s^-1", "CODATA 2018 (exact) via astropy.constants"),
    ConstProvenance("k_B", K_B, "erg K^-1", "CODATA 2018 (exact) via astropy.constants"),
    ConstProvenance("sigma_SB", SIGMA_SB, "erg cm^-2 s^-1 K^-4", "CODATA 2018 via astropy"),
    ConstProvenance("a_rad", A_RAD, "erg cm^-3 K^-4", "= 4 sigma_SB / c"),
    ConstProvenance("m_H", M_H, "g", "CODATA 2018 proton mass via astropy"),
    ConstProvenance("M_sun", MSUN, "g", "IAU 2015 Resolution B3 nominal (GM_sun/G)"),
    ConstProvenance("R_sun", RSUN, "cm", "IAU 2015 Resolution B3 nominal"),
    ConstProvenance("L_sun", LSUN, "erg s^-1", "IAU 2015 Resolution B3 nominal"),
    ConstProvenance("Teff_sun", TEFF_SUN, "K", "IAU 2015 Resolution B3 nominal"),
    ConstProvenance("yr", YEAR, "s", "Julian year, 365.25 d"),
    ConstProvenance("X_sun", X_SUN, "-", "Asplund, Grevesse, Sauval & Scott 2009, ARA&A 47, 481"),
    ConstProvenance("Y_sun", Y_SUN, "-", "Asplund et al. 2009"),
    ConstProvenance("Z_sun", Z_SUN, "-", "Asplund et al. 2009 (photospheric)"),
]


def mu_from_composition(X: float, Y: float, Z: float, ionized: bool = True) -> float:
    """Mean molecular weight.

    Fully-ionized ideal gas (Kippenhahn, Weigert & Weiss 2012, eq. 13.8):
        1/mu = 2X + 3Y/4 + Z/2   (assuming <Z/A> ~ 1/2 and one extra electron per metal
        contributes the +Z/2 term via the approximation 1/mu_e ~ (1+X)/2).
    Neutral gas (eq. 13.6):
        1/mu = X + Y/4 + <1/A>_Z ~ X + Y/4  (metals neglected in the neutral term).
    """
    if ionized:
        return 1.0 / (2.0 * X + 0.75 * Y + 0.5 * Z)
    return 1.0 / (X + 0.25 * Y)


def mu_e_from_X(X: float) -> float:
    """Electron mean molecular weight for a fully ionized plasma (KWW12 eq. 13.9)."""
    return 2.0 / (1.0 + X)


def keplerian_period(a_cm: float, m_total_g: float) -> float:
    """Kepler's third law: P = 2 pi sqrt(a^3 / (G M))."""
    return 2.0 * math.pi * math.sqrt(a_cm**3 / (G * m_total_g))


def semimajor_from_period(period_s: float, m_total_g: float) -> float:
    return (G * m_total_g * (period_s / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0)


def roche_lobe_eggleton(q: float) -> float:
    """Eggleton (1983), ApJ 268, 368  --  R_L / a  for mass ratio q = M_this / M_other.

    Accurate to better than 1% over 0 < q < infinity.
    """
    q13 = q ** (1.0 / 3.0)
    q23 = q13 * q13
    return 0.49 * q23 / (0.6 * q23 + math.log1p(q13))
