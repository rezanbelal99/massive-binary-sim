# Simulation results

*massive-binary-sim v0.1.0  ·  git `UNKNOWN(no`  ·  seed 42  ·  wallclock 0.07 s*

## 1. Initial conditions

| quantity | value |
|---|---|
| M1, M2 (ZAMS) | 40, 30 Msun |
| separation / period | 9e+03 Rsun / 1.182e+04 d |
| eccentricity | 0.1 |
| metallicity Z | 0.0142 |
| spins (Omega/Omega_crit) | 0, 0 |

### Active prescriptions

| toggle | choice |
|---|---|
| gravity | `2.5PN` |
| wind | `vink2001` |
| convection | `ledoux` |
| overshoot | `exponential` |
| am_transport | `tayler_spruit` |
| mt_stability | `zeta` |
| mt_rate | `kolb_ritter` |
| ce_formalism | `alpha_lambda` |
| opacity | `kramers` |
| reaction_rates | `reaclib` |
| sn_engine | `delayed` |
| tides | `hut1981` |
| structure_fidelity | `semi_analytic` |
| third_body | `False` |

## 2. Evolutionary timeline

| t [Myr] | event | detail |
|---|---|---|
| 0.014 | star 1 -> MS | reduced-model phase |
| 5.158 | star 1 -> HG | reduced-model phase |
| 5.253 | star 1 -> CHeB | reduced-model phase |
| 5.768 | star 1 -> preSN | reduced-model phase |
| 5.781 | star 1 -> remnant | reduced-model phase |
| 0.014 | star 2 -> MS | reduced-model phase |
| 6.504 | star 2 -> HG | reduced-model phase |
| 6.631 | star 2 -> CHeB | reduced-model phase |
| 7.322 | star 2 -> preSN | reduced-model phase |
| 7.322 | star 2 -> remnant | reduced-model phase |
| (at SN) | **SN1** | star 1 -> BH 33.91 Msun (delayed); bound after SN |
| (at SN) | **SN2** | star 2 -> BH 18.74 Msun (delayed); bound after SN |
| 7.322 | **dco_formed** | wide double compact binary (BH+BH), GW merger time 49661796998.9 Gyr > Hubble |

**Fate:** wide double compact binary (BH+BH), GW merger time 49661796998.9 Gyr > Hubble

## 3. Final-state table

| quantity | star 1 | star 2 |
|---|---|---|
| mass [Msun] | 33.91 | 18.74 |
| phase | remnant | remnant |
| remnant | BH 33.91 Msun (delayed) | |
| remnant (2) | | BH 18.74 Msun |
| separation | 1.0363e+04 Rsun | |
| period | 1.6839e+04 d | |
| eccentricity | 0.1057 | |

## 4. Derived quantities

| quantity | value | unit | formula |
|---|---|---|---|
| chirp mass M_c | 21.76 | Msun | `M_c = (m1 m2)^{3/5} / (m1+m2)^{1/5}` |
| mass ratio q | 0.5527 | - | `q = m_2 / m_1 (<=1)` |
| symmetric mass ratio eta | 0.2293 | - | `eta = m1 m2 / (m1+m2)^2` |
| integrated wind loss (star1 / star2) | 2.315 | Msun | `int |Mdot_wind,1| dt ; star2 = 9.02 Msun` |
| integrated RLOF transfer | 0 | Msun | `int |Mdot_RLOF| dt (gross, donor side)` |
| net Delta M_1 | 6.08 | Msun | `M_1(ZAMS) - M_1(final)` |
| net Delta M_2 | -11.25 | Msun | `M_2(final) - M_2(ZAMS)` |
| orbital separation (final) | 1.0363e+04 | Rsun | `state variable` |
| orbital period (final) | 1.6839e+04 | day | `P = 2 pi sqrt(a^3 / G(m1+m2))` |
| orbital binding energy | -1.1634e+47 | erg | `E = -G m1 m2 / (2a)` |
| orbital angular momentum | 5.3576e+55 | g cm^2 / s | `J = mu sqrt(G M a (1-e^2))` |
| GW merger time | 4.9662e+19 | yr | `Peters 1964 eq. 5.14 (eccentric)` |
| GW merger time | 4.9662e+10 | Gyr | `as above` |
| peak GW frequency (ISCO) | 83.5 | Hz | `f = c^3 / (6^{3/2} pi G M)` |
| char. strain at ISCO (D=100 Mpc) | 9.7255e-21 | - | `Maggiore 2007 eq. 4.44` |

## 5. Error budget

**Numerical error** (from the convergence study):

- not run in this invocation (`--convergence` flag)

**Physical / model uncertainty** dominates and is quantified in §6-7. The reduced single-star evolution model (fitting formulae, not a structure integration) is itself a systematic; see KNOWN_GAPS.md items EVO-ZAMS, STRUCT-1D.

## 6. Free-parameter sensitivity (one-at-a-time)

_not run (`--sensitivity`)._

## 7. Theory comparison (competing prescriptions)

_not run (`--theory-compare`)._

## 8. Validation suite results

| test | status | detail |
|---|---|---|
| kepler_closure | PASS | 200 orbits: max|dE/E|=1.77e-11, max|dL/L|=8.76e-12, |da/a|=1.52e-11 (tol 1e-06) |
| lane_emden_analytic | PASS | n=0: max|dtheta|=5.72e-16, |dxi1|=8.88e-16; n=1: max|dtheta|=1.55e-14, |dxi1|=1.82e-14; n=5: max|dtheta|=5.44e-15 |
| periastron_advance_1PN | PASS | 1PN d(omega)/orbit measured=6.157175e-05 rad, predicted 6 pi GM/(c^2 a(1-e^2))=6.155112e-05, frac diff=0.034% |
| hulse_taylor_decay | PASS | dPb/dt(Peters)=-2.402102e-12 s/s vs GR prediction -2.402630e-12 (Weisberg & Huang 2016); residual 0.022%. Merger time 301 Myr (lit ~300 Myr). |
| energy_conservation_2PN | PASS | 2PN conservative: secular |dE/E| over 120 orbits = 1.68e-08 (tol 1e-07); max|dL/L|=1.12e-05 |
| eddington_limit_flag | PASS | 120 Msun ZAMS Gamma_Edd=0.406; 10x-luminosity state flags Gamma>1: True |
| solar_structure_eddington | PASS | Eddington model 1 Msun: P_c=1.256e+17 (ref 2.4e17), T_c=1.185e+07 (ref 1.57e7), rho_c=76.9 (ref 150), L=0.85 Lsun, R=1.00 Rsun, Teff=5548 K (ref 5772). Within-factor-2 targets met: True. |

## 9. Assumptions & limitations

- **1-D, orbit-averaged.** No 3-D hydrodynamics; the orbit is treated by
  secular equations except in the dedicated Kepler-closure test.
- **Reduced stellar evolution.** Single-star tracks are fitting formulae
  (Hurley+ 2000-style), not Henyey structure integrations. Core masses,
  radii and lifetimes carry ~10-30% systematic uncertainty.
- **No magnetic fields** beyond the Tayler-Spruit AM-transport toggle.
- **No rotational mixing feedback** on the tracks (rotation only sets
  Omega_crit and tidal rates here).
- **Mass-transfer rate** uses a Ritter/Kolb-Ritter-style analytic nozzle
  with an efficiency prefactor (KNOWN_GAPS BIN-MDOT).
- **GW waveform** is a 0PN Newtonian chirp -- qualitative, not for matched
  filtering.
- **3.5PN dynamics** not implemented (KNOWN_GAPS PN-35).
- Full list: KNOWN_GAPS.md.

## 11. Bibliography

See `REFERENCES.bib` for the full list with DOIs. Key prescriptions used in this run: Peters (1964); Mora & Will (2004); Eggleton (1983); Ritter (1988) / Kolb & Ritter (1990); Hut (1981); Webbink (1984); Vink et al. (2001) / Bjorklund et al. (2021); Nugis & Lamers (2000); Fryer et al. (2012); Hobbs et al. (2005); Farmer et al. (2019); Kippenhahn, Weigert & Weiss (2012).
