"""Derived quantities with the formula shown alongside each value."""
from __future__ import annotations

import numpy as np

from ..dynamics import gw
from ..numerics import units as U


def summarise(result, cfg):
    h = result.history.as_arrays()
    fin = result.final
    out = []

    def add(name, value, unit, formula):
        out.append(dict(name=name, value=value, unit=unit, formula=formula))

    m1 = fin["m1_Msun"] * U.MSUN
    m2 = fin["m2_Msun"] * U.MSUN
    add("chirp mass M_c", gw.chirp_mass(m1, m2) / U.MSUN, "Msun",
        "M_c = (m1 m2)^{3/5} / (m1+m2)^{1/5}")
    add("mass ratio q", min(m1, m2) / max(m1, m2), "-", "q = m_2 / m_1 (<=1)")
    add("symmetric mass ratio eta", gw.symmetric_mass_ratio(m1, m2), "-",
        "eta = m1 m2 / (m1+m2)^2")
    wind1 = float(np.trapezoid(-h["mdot_wind1"], h["t"]))
    wind2 = float(np.trapezoid(-h["mdot_wind2"], h["t"]))
    rlof = float(np.trapezoid(np.abs(h["mdot_rlof"]), h["t"]))
    add("integrated wind loss (star1 / star2)", wind1, "Msun",
        f"int |Mdot_wind,1| dt ; star2 = {wind2:.3g} Msun")
    add("integrated RLOF transfer", rlof, "Msun", "int |Mdot_RLOF| dt (gross, donor side)")
    dm1 = h["m1"][0] - h["m1"][-1]
    dm2 = h["m2"][-1] - h["m2"][0]
    add("net Delta M_1", dm1, "Msun", "M_1(ZAMS) - M_1(final)")
    add("net Delta M_2", dm2, "Msun", "M_2(final) - M_2(ZAMS)")

    a = fin.get("a_Rsun", np.inf)
    if np.isfinite(a):
        a_cm = a * U.RSUN
        add("orbital separation (final)", a, "Rsun", "state variable")
        add("orbital period (final)", fin.get("period_day", np.inf), "day",
            "P = 2 pi sqrt(a^3 / G(m1+m2))")
        E_orb = -U.G * m1 * m2 / (2 * a_cm)
        add("orbital binding energy", E_orb, "erg", "E = -G m1 m2 / (2a)")
        mu_red = m1 * m2 / (m1 + m2)
        J_orb = mu_red * np.sqrt(U.G * (m1 + m2) * a_cm * (1 - fin["eccentricity"] ** 2))
        add("orbital angular momentum", J_orb, "g cm^2 / s",
            "J = mu sqrt(G M a (1-e^2))")

    if fin.get("remnant1") and fin.get("remnant2") and np.isfinite(a):
        tm = gw.merger_time_eccentric(a * U.RSUN, fin["eccentricity"], m1, m2)
        add("GW merger time", tm / U.YEAR, "yr",
            "Peters 1964 eq. 5.14 (eccentric)")
        add("GW merger time", tm / U.GYR, "Gyr", "as above")
        fpk = gw.peak_gw_frequency_isco(m1, m2)
        add("peak GW frequency (ISCO)", fpk, "Hz",
            "f = c^3 / (6^{3/2} pi G M)")
        D = 100e6 * U.PC
        add("char. strain at ISCO (D=100 Mpc)",
            gw.characteristic_strain_inspiral(m1, m2, fpk, D), "-",
            "Maggiore 2007 eq. 4.44")
    return out
