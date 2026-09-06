"""Roche-lobe overflow: rate, Case A/B/C classification, stability, orbital response.

References
---------
Ritter, H. (1988), A&A 202, 93        -- isothermal-nozzle RLOF rate.
Kolb, U. & Ritter, H. (1990), A&A 236, 385 -- optically-thick extension.
Soberman, Phinney & van den Heuvel (1997), A&A 327, 620 -- non-conservative
    orbital-AM bookkeeping (isotropic re-emission, L2, circumbinary ring gamma).
Hurley, Tout & Pols (2002), MNRAS 329, 897 -- q_crit tables.
Ge et al. (2015, 2020), ApJ 812, 40 / 899, 132 -- adiabatic mass-radius exponents.
"""
from __future__ import annotations

import numpy as np

from ..numerics.units import G, K_B, M_H, MSUN, RSUN, YEAR
from .roche import roche_lobe_radius

# --- adiabatic / thermal M-R exponents (zeta = dlnR/dlnM) -------------------------
# Coarse values by donor type (Hjellming & Webbink 1987; Soberman+ 1997; Ge+ 2015).
ZETA_AD = {
    "MS": 0.65,            # radiative-envelope MS star, mild expansion on mass loss
    "HG": 0.0,             # Hertzsprung-gap: near-neutral
    "CHeB": -0.3,          # convective giant envelope: expands on mass loss (unstable-prone)
    "stripped_He": 0.2,
    "preSN": 0.1,
}
Q_CRIT = {                 # q = M_donor/M_accretor above which transfer is dynamically unstable
    "MS": 3.0,             # Hurley+ 2002 / de Mink+ 2007 for massive MS donors
    "HG": 4.0,
    "CHeB": 0.9,           # deep convective envelope -> CE-prone
    "stripped_He": 3.5,
    "preSN": 1.5,
}


def classify_case(donor_phase):
    """Case A (core H burning), B (H shell, pre-He ignition), C (post core-He)."""
    return {"MS": "A", "HG": "B", "CHeB": "B", "stripped_He": "BB", "preSN": "C"}.get(
        donor_phase, "B")


def zeta_lobe(q, non_conservative_beta=1.0, gamma=1.0):
    """dlnR_L/dlnM_donor for the response of the Roche lobe to mass transfer.

    For fully conservative transfer (Soberman+ 1997 eq. with beta=1):
        zeta_L = 2.13 q - 1.67 ... use the Eggleton-derived analytic form.
    We compute it numerically from d ln R_L / d ln M_d at fixed total mass &
    conservative orbit response da/a = 2 (M_d - M_a)/(M_d M_a) dM_d... implemented
    via finite difference of R_L(a(q), q).
    """
    def RL_over_a(qd):
        # R_L/a as function of donor mass ratio qd = M_d/M_a
        from ..numerics.units import roche_lobe_eggleton
        return roche_lobe_eggleton(qd)

    # conservative: a ~ (M_d M_a)^-2 (const), so a(q) with M_d+M_a fixed
    dlnq = 1e-4
    q1, q2 = q * (1 - dlnq), q * (1 + dlnq)

    def a_factor(qd):
        # M_d + M_a = S const; M_d = S qd/(1+qd); a ~ 1/(M_d M_a)^2
        md = qd / (1 + qd)
        ma = 1.0 / (1 + qd)
        return 1.0 / (md * ma) ** 2

    RL1 = RL_over_a(q1) * a_factor(q1)
    RL2 = RL_over_a(q2) * a_factor(q2)
    # dM_d has same sign as dq (donor mass ratio grows with donor mass)
    return (np.log(RL2) - np.log(RL1)) / (np.log(q2) - np.log(q1))


def is_stable(q, donor_phase, criterion="zeta"):
    """Return (stable: bool, margin: float, note: str)."""
    if criterion == "qcrit":
        qc = Q_CRIT.get(donor_phase, 2.0)
        return (q < qc, qc - q, f"q={q:.2f} vs q_crit={qc:.2f} ({donor_phase})")
    # zeta comparison
    zad = ZETA_AD.get(donor_phase, 0.3)
    zL = zeta_lobe(q)
    return (zad > zL, zad - zL, f"zeta_ad={zad:.2f} vs zeta_L={zL:.2f} ({donor_phase})")


def mdot_rlof(donor, a, m_accretor, rate_model="kolb_ritter"):
    """RLOF mass-transfer rate [g/s] (<=0 for the donor).

    Ritter (1988): Mdot = -Mdot0 * exp((R_star - R_L) / H_p)
        Mdot0 = (2 pi / e^{1/2}) * (R_L^3 / (G M_d)) * (k T / (mu m_H))^{3/2} * rho_ph  ...
    We use the compact form (Ritter 1988 eq. 12):
        Mdot0 ~ rho_ph * v_th * Q_eff  with an effective nozzle cross-section.
    For deep overflow the Kolb & Ritter (1990) optically-thick integral adds a
    steep power-law term  Mdot ~ (Delta R / R)^{n}  with n ~ 3 (radiative) - 5.
    """
    R = donor.R
    RL = roche_lobe_radius(a, donor.m, m_accretor)
    dR = R - RL
    # photospheric scale height H_p = k T_eff R^2 / (G M mu m_H)
    mu = 0.62
    Hp = K_B * max(donor.Teff, 3e3) * R**2 / (G * donor.m * mu * M_H)
    # characteristic rate: donor mass over thermal timescale of the outer envelope
    tau_th = G * donor.m ** 2 / (donor.R * donor.L)
    Mdot0 = donor.m / tau_th * 1e-2   # nozzle efficiency prefactor (Ritter 1988 scale)

    # physical ceiling: mass transfer cannot proceed faster than the donor can be
    # stripped on its own thermal (Kelvin-Helmholtz) timescale for deep overflow
    # (Paczynski & Sienkiewicz 1972); we cap at 3x that rate.
    tau_kh = G * donor.m ** 2 / (donor.R * donor.L)
    mdot_cap = 3.0 * donor.m / tau_kh

    if dR <= 0.0:
        return -min(Mdot0 * np.exp(dR / max(Hp, 1e6)), mdot_cap)
    x = min(dR / R, 3.0)
    if rate_model == "ritter":
        md = Mdot0 * (1.0 + x / max(Hp / R, 1e-3))
    else:  # kolb_ritter: optically-thick power law, exponent ~3 (radiative envelope)
        md = Mdot0 * (1.0 + 10.0 * x ** 3.0 + x / max(Hp / R, 1e-3))
    return -min(md, mdot_cap)


def orbital_response(a, m_donor, m_accretor, mdot_donor, beta=1.0, gamma=1.5,
                     mode="isotropic_reemission"):
    """da/dt for mass transfer with retention fraction `beta`.

    General relation from  a ~ J_orb^2 (M_d + M_a) / (G M_d^2 M_a^2)  differentiated:
        dln a = 2 dln J_orb - 2 dln M_d - 2 dln M_a + dln(M_d + M_a).

    dln J_orb is set by the specific AM h_loss (in units of J_orb / M_reduced ... here
    parameterised as  h_loss = gamma_loss * a^2 Omega_orb) carried by the
    (1 - beta) fraction of transferred mass that leaves the system:

      isotropic_reemission : gamma_loss = (M_d / M_a) ^ 2 / (1 + M_d/M_a)  ... =
          matter leaves with the accretor's specific orbital AM
          (Tauris & van den Heuvel 2006, "Compact Stellar X-ray Sources" Ch. 16,
           eq. 16.20 -- "Jeans mode about the accretor").
      L2                   : gamma_loss = (1.2)^2 (1 + M_d/M_a)   (outer-Lagrange
          lever arm; MacLeod & Loeb 2020, ApJ 895, 29).
      circumbinary_ring    : gamma_loss = gamma * (1 + M_d/M_a)^2 / (M_d/M_a)
          (Soberman, Phinney & van den Heuvel 1997, eq. 20; gamma ~ 1.3-1.7).

    Conservative limit beta=1 reduces to the textbook
        dln a = -2 (dM_d/M_d) (1 - M_d/M_a)   (orbit shrinks while donor heavier).
    """
    Md, Ma = m_donor, m_accretor
    Mdot_d = mdot_donor                     # <0
    Mdot_a = -beta * Mdot_d                 # >0
    Mtot = Md + Ma
    q = Md / Ma

    # AM lost per unit mass leaving the system, in units of J_orb/Mtot * (Mtot^2/(Md Ma))
    # -> use the fractional-J formulation directly (Soberman+ 1997 eq. 16):
    if mode == "isotropic_reemission":
        gamma_loss = q                      # h_loss = h_accretor  (Jeans about accretor)
    elif mode == "L2":
        gamma_loss = 1.2 ** 2 * (1.0 + q)
    elif mode == "circumbinary_ring":
        gamma_loss = gamma * (1.0 + q) ** 2 / q
    else:
        gamma_loss = q

    mass_lost_rate = -(Mdot_d + Mdot_a)     # >=0, leaves the system
    dlnJ = -gamma_loss * mass_lost_rate / Mtot   # per second

    dlnMtot = (Mdot_d + Mdot_a) / Mtot
    dlnMd = Mdot_d / Md
    dlnMa = Mdot_a / Ma
    dln_a = 2.0 * dlnJ - 2.0 * dlnMd - 2.0 * dlnMa + dlnMtot
    return a * dln_a
