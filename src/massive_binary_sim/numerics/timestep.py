"""Adaptive timestep control for the evolution driver.

The step is the minimum of a set of physically-motivated limits.  The binding
criterion (the one that set the step) is logged at every step so the report can
show a timestep-history / binding-criterion figure.

Timescales
----------
  orbital     : f_orb * P_orb                         (resolve the orbit's secular drift)
  nuclear     : f_nuc * M_core c^2 X / L_nuc  ~  E_available / L         (burning)
  thermal     : f_th * tau_KH = f_th * G M^2 / (R L)   (Kelvin-Helmholtz)
  mass_transfer: f_mt * min(M_donor/|Mdot_RLOF|)
  wind        : f_w  * min(M / |Mdot_wind|)
  gw          : f_gw * a / |da/dt|                     (Peters inspiral)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .units import C, G, YEAR


@dataclass
class StepLimits:
    f_orbit: float = 0.02
    f_nuclear: float = 0.02
    f_thermal: float = 0.05
    f_mt: float = 0.1
    f_wind: float = 0.1
    f_gw: float = 0.01
    dt_min: float = 1.0e-3 * YEAR
    dt_max: float = 1.0e5 * YEAR


def kelvin_helmholtz_time(M, R, L):
    """tau_KH = G M^2 / (R L)   (Kippenhahn, Weigert & Weiss 2012, eq. 3.20)."""
    return G * M * M / (R * L)


def nuclear_time(M_core, L, efficiency=0.007, X_fuel=0.7):
    """Rough nuclear timescale: rest-mass energy of available fuel / luminosity.

    efficiency 0.007 for H->He (Kippenhahn+ 2012, sec. 18.1).
    """
    if L <= 0:
        return np.inf
    return efficiency * X_fuel * M_core * C**2 / L


def choose_timestep(state, limits: StepLimits):
    """`state` is a dict with keys (subset ok):
        P_orbit, a, dadt_gw, M1,R1,L1, M2,R2,L2, Mdot_rlof, Mdot_wind1, Mdot_wind2,
        M_core1, M_core2
    Returns (dt_seconds, binding_criterion_str).
    """
    cands = {}
    # Orbit criterion: fraction of the a-evolution timescale a/|da/dt| from ALL
    # contributions (winds + tides + MT + GW), NOT tied to the orbital period.
    # The bare orbital period only limits the step during active RLOF.
    a = state.get("a")
    if a and state.get("dadt_total"):
        da = abs(state["dadt_total"])
        if da > 0:
            cands["orbit(a/adot)"] = limits.f_orbit * a / da
    if state.get("rlof_active") and state.get("P_orbit"):
        # during RLOF also cap at a few thousand orbits so contact is not overshot
        cands["orbital_period"] = 4000.0 * state["P_orbit"]
    if a and state.get("dadt_gw"):
        da = abs(state["dadt_gw"])
        if da > 0:
            cands["gw"] = limits.f_gw * a / da
    for i in (1, 2):
        M, R, L = state.get(f"M{i}"), state.get(f"R{i}"), state.get(f"L{i}")
        if M and R and L and L > 0:
            # the Kelvin-Helmholtz limit only bites when the star is thermally
            # *relaxing* (HG crossing, post-CE readjustment); a star in thermal
            # equilibrium on the MS or in core He burning is not step-limited by it.
            if state.get(f"relaxing{i}"):
                cands[f"thermal{i}"] = limits.f_thermal * kelvin_helmholtz_time(M, R, L)
            mc = state.get(f"M_core{i}", 0.1 * M)
            cands[f"nuclear{i}"] = limits.f_nuclear * nuclear_time(max(mc, 0.05 * M), L)
    for i in (1, 2):
        tphase = state.get(f"t_to_phase_end{i}")
        tms = state.get(f"t_ms{i}")
        if tphase and tphase > 0:
            # pace toward the next phase boundary, but never let this drive the step
            # below ~1% of the star's MS lifetime (it is a guide, not a stiff limit)
            floor = max(0.01 * tms, 100.0 * YEAR) if tms else 1.0e3 * YEAR
            cands[f"phase_end{i}"] = max(0.1 * tphase, floor)
    if state.get("Mdot_rlof"):
        md = abs(state["Mdot_rlof"])
        if md > 0:
            cands["mass_transfer"] = limits.f_mt * state.get("M_donor", state.get("M1", 1)) / md
    for i in (1, 2):
        mw = state.get(f"Mdot_wind{i}")
        if mw and abs(mw) > 0:
            cands[f"wind{i}"] = limits.f_wind * state.get(f"M{i}", 1.0) / abs(mw)

    if not cands:
        return limits.dt_max, "dt_max(default)"
    crit = min(cands, key=cands.get)
    dt = cands[crit]
    dt = float(np.clip(dt, limits.dt_min, limits.dt_max))
    if dt == limits.dt_min:
        crit = "dt_min(floor)"
    elif dt == limits.dt_max:
        crit = "dt_max(ceiling)"
    return dt, crit
