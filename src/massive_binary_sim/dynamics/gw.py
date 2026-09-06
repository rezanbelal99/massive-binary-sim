"""Gravitational-wave emission: orbit-averaged inspiral and leading-order waveform.

References
----------
Peters, P. C. (1964), Phys. Rev. 136, B1224   -- orbit-averaged da/dt, de/dt, T_merge.
Peters & Mathews (1963), Phys. Rev. 131, 435  -- power radiated, enhancement factor.
Maggiore (2007), "Gravitational Waves, Vol. 1", Ch. 4  -- restricted PN waveform,
    chirp mass, characteristic strain.
Robson, Cornish & Liu (2019), CQG 36, 105011  -- analytic LISA sensitivity curve.
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import C, G, MSUN, PC  # noqa: F401

MPC = 1.0e6 * PC


# ---------------------------------------------------------------------------------------
# Orbit-averaged evolution -- Peters (1964)
# ---------------------------------------------------------------------------------------
def peters_dadt(a, e, m1, m2):
    """da/dt  (Peters 1964, Eq. 5.6).  CGS; returns cm/s (negative)."""
    M1, M2 = m1, m2
    fe = 1.0 + (73.0 / 24.0) * e**2 + (37.0 / 96.0) * e**4
    return -(64.0 / 5.0) * G**3 * M1 * M2 * (M1 + M2) / (
        C**5 * a**3 * (1.0 - e**2) ** 3.5) * fe


def peters_dedt(a, e, m1, m2):
    """de/dt  (Peters 1964, Eq. 5.7).  CGS; returns 1/s (negative)."""
    if e <= 0.0:
        return 0.0
    M1, M2 = m1, m2
    ge = 1.0 + (121.0 / 304.0) * e**2
    return -(304.0 / 15.0) * e * G**3 * M1 * M2 * (M1 + M2) / (
        C**5 * a**4 * (1.0 - e**2) ** 2.5) * ge


def merger_time_circular(a, m1, m2):
    """T_merge for a circular orbit  (Peters 1964, Eq. 5.10 with e=0).

    T = 5 c^5 a^4 / (256 G^3 m1 m2 (m1+m2))   [s]
    """
    return 5.0 * C**5 * a**4 / (256.0 * G**3 * m1 * m2 * (m1 + m2))


def merger_time_eccentric(a0, e0, m1, m2, n_steps=4000):
    """T_merge for an eccentric orbit by integrating da/de (Peters 1964, Eq. 5.14).

    Uses the closed relation a(e) = c0 * e^{12/19}/(1-e^2) * [1 + (121/304)e^2]^{870/2299}
    then  T = (12/19) c0^4 / beta * integral_0^{e0} de' * ...   (Peters 1964 Eq. 5.14).
    """
    if e0 <= 1e-6:
        return merger_time_circular(a0, m1, m2)
    beta = (64.0 / 5.0) * G**3 * m1 * m2 * (m1 + m2) / C**5

    def a_of_e(e):
        return (e ** (12.0 / 19.0) / (1.0 - e**2)
                * (1.0 + (121.0 / 304.0) * e**2) ** (870.0 / 2299.0))

    c0 = a0 / a_of_e(e0)
    e = np.linspace(1e-8, e0, n_steps)
    integrand = (e ** (29.0 / 19.0) * (1.0 + (121.0 / 304.0) * e**2) ** (1181.0 / 2299.0)
                 / (1.0 - e**2) ** 1.5)
    integral = np.trapezoid(integrand, e)
    return (12.0 / 19.0) * c0**4 / beta * integral


def peters_dPb_dt(a, e, m1, m2):
    """Orbit-period decay rate  dP_b/dt  [s/s]  from Peters (1964).

    P_b = 2 pi a^{3/2} (G M)^{-1/2}  ->  dP_b/dt = 3 pi sqrt(a / (G M)) * da/dt.
    This is the quantity measured for PSR B1913+16.
    """
    M = m1 + m2
    return 3.0 * np.pi * np.sqrt(a / (G * M)) * peters_dadt(a, e, m1, m2)


def chirp_mass(m1, m2):
    """M_c = (m1 m2)^{3/5} / (m1+m2)^{1/5}."""
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def symmetric_mass_ratio(m1, m2):
    return m1 * m2 / (m1 + m2) ** 2


def peak_gw_frequency_isco(m1, m2):
    """GW frequency at the Schwarzschild ISCO of the total mass (2 x orbital).

    f_GW = 2 f_orb(r=6GM/c^2) = c^3 / (6^{3/2} pi G M).
    """
    M = m1 + m2
    return C**3 / (6.0 ** 1.5 * np.pi * G * M)


def strain_amplitude(m1, m2, f_gw, distance):
    """Leading-order (Newtonian-quadrupole) sky/inclination-averaged strain amplitude.

    h_0 = (4/sqrt(5)) (G M_c / c^2 D) (pi f_gw G M_c / c^3)^{2/3}
    (Maggiore 2007 eq. 4.34, angle-averaged over pattern and inclination).
    """
    mc = chirp_mass(m1, m2)
    x = (np.pi * f_gw * G * mc / C**3) ** (2.0 / 3.0)
    return (4.0 / np.sqrt(5.0)) * (G * mc / (C**2 * distance)) * x


def characteristic_strain_inspiral(m1, m2, f_gw, distance):
    """Characteristic strain of a Newtonian inspiral (Maggiore 2007 eq. 4.44):

        h_c(f) = (2/pi^{2/3}) * sqrt(5/24) * (G M_c)^{5/6} / (c^{3/2} D) * f^{-1/6}

    (sky/orientation-averaged; valid f_start << f << f_ISCO).
    """
    mc = chirp_mass(m1, m2)
    return ((2.0 / np.pi ** (2.0 / 3.0)) * np.sqrt(5.0 / 24.0)
            * (G * mc) ** (5.0 / 6.0) / (C ** 1.5 * distance) * f_gw ** (-1.0 / 6.0))


def inspiral_waveform(m1, m2, f_start, distance, t_grid, e0=0.0):
    """Restricted 0PN (Newtonian chirp) time-domain waveform h_+(t), h_x(t).

    Uses the Newtonian-quadrupole frequency sweep
        f(t) = f_start * (1 - t/tau)^(-3/8)
    with tau the Newtonian coalescence time from f_start (Maggiore 2007 eq. 4.21).
    Eccentricity is ignored in the waveform (circularises fast); a warning-level
    caveat is recorded in the report.  This is a qualitative chirp, adequate for
    the characteristic-strain figure, NOT for matched filtering.
    """
    mc = chirp_mass(m1, m2)
    tau = (5.0 / 256.0) * (C**3 / (G * mc)) ** (5.0 / 3.0) * (np.pi * f_start) ** (-8.0 / 3.0)
    tt = np.clip(t_grid, 0.0, tau * (1 - 1e-6))
    f = f_start * (1.0 - tt / tau) ** (-3.0 / 8.0)
    # accumulated GW phase
    phase = -2.0 * np.pi * f_start * tau * (8.0 / 5.0) * ((1.0 - tt / tau) ** (5.0 / 8.0) - 1.0)
    x = (np.pi * f * G * mc / C**3) ** (2.0 / 3.0)
    amp = (4.0 / np.sqrt(5.0)) * (G * mc / (C**2 * distance)) * x
    h_plus = amp * np.cos(phase)
    h_cross = amp * np.sin(phase)
    return dict(t=tt, f=f, h_plus=h_plus, h_cross=h_cross, tau=tau)


# ---------------------------------------------------------------------------------------
# Detector sensitivity curves (analytic fits)
# ---------------------------------------------------------------------------------------
def lisa_characteristic_strain(f):
    """Sky-averaged LISA characteristic strain sqrt(f S_n(f)).

    Robson, Cornish & Liu (2019), CQG 36, 105011, Eqs. (1)-(2) & Table 1
    (L = 2.5e9 m, 4-year mission; galactic confusion noise included).
    """
    f = np.asarray(f, dtype=float)
    L = 2.5e9  # m
    f_star = 19.09e-3
    # single-link optical metrology + acceleration noise (m^2/Hz)
    P_oms = (1.5e-11) ** 2 * (1.0 + (2.0e-3 / f) ** 4)
    P_acc = (3.0e-15) ** 2 * (1.0 + (0.4e-3 / f) ** 2) * (1.0 + (f / 8.0e-3) ** 4)
    Pn = (P_oms + 2.0 * (1.0 + np.cos(f / f_star) ** 2) * P_acc / (2.0 * np.pi * f) ** 4) / L**2
    R = 3.0 / 10.0 / (1.0 + 0.6 * (f / f_star) ** 2)
    Sn = Pn / R
    # galactic binary confusion (Eq. 14, 4-yr coefficients)
    A = 9e-45
    Sc = (A * f ** (-7.0 / 3.0) * np.exp(-f ** 0.138 - 221.0 * f * np.sin(521.0 * f))
          * (1.0 + np.tanh(1680.0 * (0.00113 - f))))
    Sn = Sn + Sc
    return np.sqrt(f * Sn)


def ligo_aligo_characteristic_strain(f):
    """Approximate Advanced LIGO design characteristic strain sqrt(f S_n).

    Analytic broken-power-law fit to the aLIGO design ZERO_DET_high_P curve
    (Abbott+ 2018, LRR 21, 3, Fig. 1).  Accurate to a factor ~2; used only for
    the visual overlay.  Marked APPROXIMATE in the figure caption.
    """
    f = np.asarray(f, dtype=float)
    x = f / 245.4
    Sn = 1.0e-48 * (0.0152 * x ** -4.0 + 0.2935 * x ** (9.0 / 4.0)
                    + 2.7951 * x ** (3.0 / 2.0) - 6.5080 * x ** (3.0 / 4.0) + 17.7622)
    Sn = np.where(f < 9.0, np.inf, Sn)
    return np.sqrt(f * Sn)


def einstein_telescope_characteristic_strain(f):
    """Einstein Telescope ET-D sensitivity, crude analytic proxy.

    Rescales the aLIGO fit by ~10x amplitude and extends the low-f wall to 3 Hz
    (Hild+ 2011, CQG 28, 094013).  APPROXIMATE -- overlay only.
    """
    f = np.asarray(f, dtype=float)
    base = ligo_aligo_characteristic_strain(np.maximum(f, 3.0001)) / 10.0
    return np.where(f < 3.0, np.inf, base)
