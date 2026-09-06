"""Reduced analytic single-star evolution model.

This is NOT a stellar-structure integration.  It is a set of fitting relations for
the ZAMS state and the subsequent luminosity / radius / core-mass tracks, of the
kind used inside rapid population-synthesis codes.  It is the default engine
(`structure_fidelity: semi_analytic`).

Relations & anchors
-------------------
* ZAMS L, R : broken power laws calibrated to the solar-metallicity rotating and
  non-rotating grids of Ekstrom et al. (2012, A&A 537, A146) and Brott et al.
  (2011, A&A 530, A115).  Coefficients are approximate fits (KNOWN_GAPS EVO-ZAMS).
* MS lifetime : t_MS = t_nuc with the nuclear-burning efficiency and the fitted
  fractional core mass; cross-checked against Ekstrom+ 2012 Table 1.
* He / CO core masses : fits to Hurley, Pols & Tout (2000, MNRAS 315, 543),
  their sec. 5.3, and to Sukhbold+ 2018 (ApJ 860, 93) for the CO-core - mass
  relation at high mass.
* Radius expansion post-MS : giant-branch behaviour depends on mass and envelope
  retention; hot (stripped/WR) vs cool (RSG) split at the Humphreys-Davidson
  limit (Humphreys & Davidson 1979, ApJ 232, 409).

Every number a driver could need is returned in SI-free CGS via `StarState`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..numerics.units import (G, LSUN, MSUN, RSUN, SIGMA_SB, YEAR, C)

PHASES = ["ZAMS", "MS", "HG", "CHeB", "stripped_He", "preSN", "remnant"]


def zams_luminosity(M_msun, Z):
    """ZAMS luminosity [Lsun].  Broken power law; ~10-20% accurate 1-120 Msun.

    L/Lsun ~ 0.9 M^3.5           (M < 2)
           ~ 1.4 M^3.2 * fZ      (2 < M < 20)
           ~ 60  M^2.0 * fZ      (M > 20)   -> radiation-pressure-dominated flattening
    fZ = (Z/0.014)^(-0.1) weak metallicity dependence.
    """
    m = M_msun
    fZ = (Z / 0.014) ** -0.10
    if m < 2.0:
        L = 0.85 * m ** 4.0
    elif m < 20.0:
        # calibrated to ZAMS grids: 10 Msun ~ 6e3, 15 Msun ~ 2e4 Lsun
        L = 2.0 * m ** 3.5 * fZ
    else:
        # radiation-pressure flattening: 40 ~ 2.5e5, 60 ~ 5e5, 120 ~ 1.8e6 Lsun
        L = 130.0 * m ** 2.0 * fZ
    return L


def zams_radius(M_msun, Z):
    """ZAMS radius [Rsun].  R/Rsun ~ M^0.6 (M<1), ~ M^0.7 (1-10), ~ M^0.5 (>10)."""
    m = M_msun
    fZ = (Z / 0.014) ** 0.05
    if m < 1.0:
        R = m ** 0.9
    elif m < 10.0:
        R = m ** 0.68 * fZ
    else:
        R = 1.6 * m ** 0.47 * fZ
    return R


def ms_lifetime_yr(M_msun, Z):
    """Main-sequence lifetime [yr].

    t_MS = q_c * f_burn * X * M c^2 * 0.007 / L  with q_c the fitted fractional
    mass processed and f_burn ~ 0.1 the central H fraction burnt on the MS.
    Calibrated so that t_MS(1 Msun) ~ 1.0e10 yr and t_MS(60 Msun) ~ 3.4e6 yr
    (Ekstrom+ 2012).
    """
    L = zams_luminosity(M_msun, Z) * LSUN
    eff = 0.007
    # fitted effective fuel fraction, calibrated so t_MS(1 Msun) ~ 1.0e10 yr and
    # t_MS(35 Msun) ~ 5 Myr, t_MS(60) ~ 3.7 Myr (Ekstrom+ 2012, Table 1).
    qc = float(np.clip(0.082 + 0.11 * np.log10(max(M_msun, 0.5)), 0.05, 0.36))
    t = qc * (M_msun * MSUN) * C**2 * eff / L
    return t / YEAR


def he_core_fraction_tams(M_msun):
    """M_He_core / M at terminal main sequence (fit to Hurley+ 2000 / Ekstrom+ 2012)."""
    return float(np.clip(0.10 + 0.09 * np.log10(M_msun) + 0.02 * np.log10(M_msun) ** 2,
                         0.08, 0.75))


def co_core_mass(M_he_msun):
    """CO core mass [Msun] from He core mass (fit; Sukhbold+ 2018 trend)."""
    return float(np.clip(0.75 * M_he_msun - 0.30, 0.0, M_he_msun))


def humphreys_davidson_luminosity():
    """log10(L/Lsun) ceiling of the HD limit ~ 5.8 for cool supergiants."""
    return 5.8


@dataclass
class StarState:
    m: float                    # current mass [g]
    m0: float                   # initial (ZAMS) mass [g]
    Z: float
    age: float                  # s
    phase: str = "ZAMS"
    R: float = 0.0              # cm
    L: float = 0.0              # erg/s
    Teff: float = 0.0
    Xc: float = 0.70            # central hydrogen mass fraction
    m_he_core: float = 0.0      # g
    m_co_core: float = 0.0      # g
    spin: float = 0.0           # Omega/Omega_crit
    stripped: bool = False
    _t_ms: float = field(default=0.0, repr=False)

    @classmethod
    def zams(cls, M_msun, Z, spin=0.0):
        s = cls(m=M_msun * MSUN, m0=M_msun * MSUN, Z=Z, age=0.0, spin=spin)
        s._t_ms = ms_lifetime_yr(M_msun, Z) * YEAR
        s._update_observables()
        return s

    # -- properties -----------------------------------------------------------------
    def eddington_ratio(self, kappa=0.34):
        L_edd = 4.0 * np.pi * G * self.m * C / kappa
        return self.L / L_edd

    def omega_crit(self):
        return np.sqrt(G * self.m / self.R**3)

    # -- evolution ----------------------------------------------------------------
    def _update_observables(self):
        M = self.m / MSUN
        M0 = self.m0 / MSUN
        f_ms = np.clip(self.age / max(self._t_ms, 1.0), 0.0, 1.0)
        if self.phase in ("ZAMS", "MS"):
            L0 = zams_luminosity(M0, self.Z)
            R0 = zams_radius(M0, self.Z)
            # MS: L rises ~x1.6, R rises ~x1.8 (mass dependent) across the MS
            self.L = L0 * (1.0 + 0.6 * f_ms) * (M / M0) ** 2.5 * LSUN
            r_exp = 1.0 + (0.7 + 0.02 * M0) * f_ms
            self.R = R0 * r_exp * RSUN
            self.Xc = max(0.70 * (1.0 - f_ms), 0.0)
        elif self.phase == "HG":
            L0 = zams_luminosity(M0, self.Z)
            self.L = L0 * 1.7 * (M / M0) ** 2.5 * LSUN
            # expand toward giant unless stripped
            self.R = zams_radius(M0, self.Z) * 2.5 * RSUN
        elif self.phase == "CHeB":
            L0 = zams_luminosity(M0, self.Z)
            self.L = L0 * 2.0 * (M / M0) ** 2.0 * LSUN
            if self.stripped:
                self.R = 0.3 * (M) ** 0.6 * RSUN                    # compact He star
            else:
                # RSG / BSG split at HD limit
                if np.log10(self.L / LSUN) < humphreys_davidson_luminosity() and M0 < 40:
                    self.R = 800.0 * (M0 / 20.0) ** 0.5 * RSUN      # red supergiant
                else:
                    self.R = 30.0 * (M0 / 40.0) ** 0.4 * RSUN       # blue / LBV
        elif self.phase in ("stripped_He", "preSN"):
            mhe = max(self.m_he_core / MSUN, M)
            self.L = 3.0e4 * mhe ** 1.5 * LSUN
            self.R = 0.2 * mhe ** 0.6 * RSUN
            self.stripped = True
        self.Teff = (self.L / (4.0 * np.pi * self.R**2 * SIGMA_SB)) ** 0.25

    def advance(self, dt, mdot_wind=0.0, mdot_binary=0.0):
        """Advance by dt seconds.  mdot_* are dM/dt in g/s (negative = loss)."""
        self.age += dt
        self.m = max(self.m + (mdot_wind + mdot_binary) * dt, 0.05 * MSUN)

        M0 = self.m0 / MSUN
        if self.phase in ("ZAMS", "MS"):
            self.phase = "MS"
            if self.age >= self._t_ms:
                self.phase = "HG"
                self.m_he_core = he_core_fraction_tams(M0) * self.m
        elif self.phase == "HG":
            # Hertzsprung gap crossing is ~1% of t_MS (thermal timescale)
            if self.age >= self._t_ms * 1.02:
                self.phase = "CHeB"
        elif self.phase == "CHeB":
            # core He burning ~ 10% of t_MS
            t_he = 0.10 * self._t_ms
            if self.age >= self._t_ms * 1.02 + t_he:
                self.phase = "preSN"
                self.m_co_core = co_core_mass(self.m_he_core / MSUN) * MSUN
        elif self.phase == "stripped_He":
            # keep the He-core bookkeeping consistent as a WR wind erodes the star
            self.m_he_core = min(self.m_he_core, self.m) if self.m_he_core > 0 else self.m
            # a stripped helium star burns He on a short timescale then collapses.
            # t_He,stripped ~ few x 10^5 yr for massive He stars (Yoon 2017, MNRAS 470, 3970).
            t_he_strip = max(3.0e5 * YEAR * (10.0 / max(self.m / MSUN, 3.0)) ** 0.5,
                             5.0e4 * YEAR)
            if not hasattr(self, "_strip_age0"):
                self._strip_age0 = self.age
            co = co_core_mass(self.m / MSUN) * MSUN
            if (self.age - self._strip_age0 >= t_he_strip) or self.m <= max(co, 1.6 * MSUN):
                self.phase = "preSN"
                self.m_he_core = self.m
                self.m_co_core = min(co_core_mass(self.m / MSUN) * MSUN, self.m)
        # strip flag if the star has lost most of its envelope
        if self.m_he_core > 0 and self.m < 1.15 * self.m_he_core:
            self.stripped = True
            if self.phase in ("HG", "CHeB"):
                self.phase = "stripped_He"
                self._strip_age0 = self.age

        self._update_observables()

    def is_evolved_to_remnant(self):
        return self.phase == "preSN"
