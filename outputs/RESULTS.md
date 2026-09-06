# Simulation results

*massive-binary-sim v0.1.0  ·  git `3b6e4a4840`  ·  seed 42  ·  wallclock 0.34 s*

## 1. Initial conditions

| quantity | value |
|---|---|
| M1, M2 (ZAMS) | 35, 25 Msun |
| separation / period | 50 Rsun / 5.287 d |
| eccentricity | 0.05 |
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
| 0.016 | star 1 -> MS | reduced-model phase |
| 2.600 | star 1 -> stripped_He | reduced-model phase |
| 2.852 | star 1 -> preSN | reduced-model phase |
| 2.863 | star 1 -> remnant | reduced-model phase |
| 0.016 | star 2 -> MS | reduced-model phase |
| 7.533 | star 2 -> HG | reduced-model phase |
| 7.674 | star 2 -> CHeB | reduced-model phase |
| 7.674 | star 2 -> stripped_He | reduced-model phase |
| 8.007 | star 2 -> preSN | reduced-model phase |
| 8.007 | star 2 -> remnant | reduced-model phase |
| 2.599 | **mass_transfer_end** | Case A transfer complete; donor -> 14.70 Msun stripped He star |
| (at SN) | **SN1** | star 1 -> BH 13.23 Msun (delayed); bound after SN |
| 7.674 | **common_envelope** | survived: a_f=9.034 Rsun (spiral-in x5.2) |
| (at SN) | **SN2** | star 2 -> BH 5.65 Msun (delayed); bound after SN |
| 8.007 | **dco_formed** | double compact binary (BH+BH) merging via GW in 326.0 Myr |

**Fate:** double compact binary (BH+BH) merging via GW in 326.0 Myr

## 3. Final-state table

| quantity | star 1 | star 2 |
|---|---|---|
| mass [Msun] | 13.23 | 5.653 |
| phase | remnant | remnant |
| remnant | BH 13.23 Msun (delayed) | |
| remnant (2) | | BH 5.653 Msun |
| separation | 8.064 Rsun | |
| period | 0.6105 d | |
| eccentricity | 0.2916 | |

## 4. Derived quantities

| quantity | value | unit | formula |
|---|---|---|---|
| chirp mass M_c | 7.397 | Msun | `M_c = (m1 m2)^{3/5} / (m1+m2)^{1/5}` |
| mass ratio q | 0.4273 | - | `q = m_2 / m_1 (<=1)` |
| symmetric mass ratio eta | 0.2097 | - | `eta = m1 m2 / (m1+m2)^2` |
| integrated wind loss (star1 / star2) | 0.3863 | Msun | `int |Mdot_wind,1| dt ; star2 = 3.15 Msun` |
| integrated RLOF transfer | 23.47 | Msun | `int |Mdot_RLOF| dt (gross, donor side)` |
| net Delta M_1 | 21.77 | Msun | `M_1(ZAMS) - M_1(final)` |
| net Delta M_2 | -19.35 | Msun | `M_2(final) - M_2(ZAMS)` |
| orbital separation (final) | 8.064 | Rsun | `state variable` |
| orbital period (final) | 0.6105 | day | `P = 2 pi sqrt(a^3 / G(m1+m2))` |
| orbital binding energy | -1.7587e+49 | erg | `E = -G m1 m2 / (2a)` |
| orbital angular momentum | 2.8244e+53 | g cm^2 / s | `J = mu sqrt(G M a (1-e^2))` |
| GW merger time | 3.2603e+08 | yr | `Peters 1964 eq. 5.14 (eccentric)` |
| GW merger time | 0.326 | Gyr | `as above` |
| peak GW frequency (ISCO) | 232.9 | Hz | `f = c^3 / (6^{3/2} pi G M)` |
| char. strain at ISCO (D=100 Mpc) | 3.3357e-21 | - | `Maggiore 2007 eq. 4.44` |

## 5. Error budget

**Numerical error** (from the convergence study):

- final a converges to 3.000069 Rsun; rel error at rtol=1e-8 is 1.75e-06, at 1e-11 is 1.07e-09; mean observed order ~ 1.07 per tol-decade
- relative error in final `a` vs the tightest tolerance: 1.7505e-06 at rtol=1e-8 -> 1.0713e-09 at rtol=1e-11
- mean observed order of convergence: 1.07 per tolerance decade (DOP853, 8th-order adaptive)
- **The numerical error (~1e-9 relative) is far below the physical/model uncertainty quantified in sections 6-7.**

**Physical / model uncertainty** dominates and is quantified in §6-7. The reduced single-star evolution model (fitting formulae, not a structure integration) is itself a systematic; see KNOWN_GAPS.md items EVO-ZAMS, STRUCT-1D.

## 6. Free-parameter sensitivity (one-at-a-time)

Outcome metric: **chirp_mass**. Parameters ranked by outcome spread across their literature range.

| rank | parameter | range | outcome spread | provenance |
|---|---|---|---|---|
| 1 | `beta_accretion` | 0-1 | 3.05 | fraction retained; poorly constrained, see Petrovic+ 2005 |
| 2 | `eta_wind` | 0.1-3 | 0.5548 | overall mass-loss multiplier (convention) |
| 3 | `clumping` | 1-10 | 8.3445e-04 | wind clumping factor; Vink 2001 rates assume smooth |
| 4 | `alpha_mlt` | 1.5-2.2 | 0 | solar-model calibrated; e.g. Magic+ 2015 A&A 573 A89 |
| 5 | `f_overshoot` | 0-0.04 | 0 | Brott+ 2011 A&A 530 A115; Claret & Torres 2018 ApJ 859 100 |
| 6 | `f0_overshoot` | 0-0.01 | 0 | Herwig 2000 A&A 360 952 (exponential-diffusive) |
| 7 | `alpha_sc` | 0.001-100 | 0 | Langer+ 1983 A&A 126 207 -- essentially unconstrained |
| 8 | `gamma_amloss` | 1-4 | 0 | specific AM of non-conservative outflow; Nelemans+ 2000 A&A 360 1011 |
| 9 | `alpha_ce` | 0.1-10 | 0 | CE efficiency; Ivanova+ 2013 A&ARv 21 59 |
| 10 | `lambda_ce` | 0.05-2 | 0 | envelope structure param; Xu & Li 2010 ApJ 716 114; Dewi & Tauris 2000 |
| 11 | `c12_ag_rate` | 0.6-1.6 | 0 | multiplier on 12C(a,g)16O S-factor; deBoer+ 2017 RvMP 89 035007 |
| 12 | `kick_sigma` | 0-400 | 0 | Hobbs+ 2005 MNRAS 360 974 (km/s, Maxwellian 1D sigma) |
| 13 | `overshoot_step` | 0-0.5 | 0 | step overshoot in H_p; Schaller+ 1992 / Ekstrom+ 2012 |

**Dominant uncertain input:** `beta_accretion` (outcome varies by 3.05 across its plausible range).

## 7. Theory comparison (competing prescriptions)

| toggle | option | default? | fate | M1 | M2 | Mchirp | a [Rsun] | log10 t_merge |
|---|---|---|---|---|---|---|---|---|
| wind | `vink2001` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| wind | `bjorklund2021` |  | common-envelope merger (Case A, q=12.12) | 2.86 | 34.68 | n/a | 58.13 | n/a |
| convection | `schwarzschild` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| convection | `ledoux` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| overshoot | `none` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| overshoot | `step` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| overshoot | `exponential` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| am_transport | `none` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| am_transport | `tayler_spruit` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| mt_stability | `qcrit` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| mt_stability | `zeta` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| mt_rate | `ritter` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| mt_rate | `kolb_ritter` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| ce_formalism | `alpha_lambda` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| ce_formalism | `gamma` |  | wide double compact binary (BH+BH), GW merger  | 13.23 | 5.653 | 7.397 | 27.97 | 10.21 |
| opacity | `electron_scattering` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| opacity | `kramers` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| reaction_rates | `nacre` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| reaction_rates | `reaclib` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| sn_engine | `rapid` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 7.596 | 8.661 | 9.414 | 8.748 |
| sn_engine | `delayed` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| sn_engine | `direct_collapse` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 7.596 | 8.661 | 9.414 | 8.748 |
| gravity | `newtonian` |  | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| gravity | `2.5PN` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |
| tides | `none` |  | common-envelope merger (Case B, q=1.37) | 33.52 | 24.53 | n/a | 51.48 | n/a |
| tides | `hut1981` | YES | double compact binary (BH+BH) merging via GW i | 13.23 | 5.653 | 7.397 | 8.064 | 8.513 |

*The scientific payload: the outcome is a mapped dependence on prescription choice, not a single number.*

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

## 10. Figures

- **hr_diagram** — `outputs/hr_diagram.png`
- **orbital_evolution** — `outputs/orbital_evolution.png`
- **mass_evolution** — `outputs/mass_evolution.png`
- **roche** — `outputs/roche.png`
- **kippenhahn1** — `outputs/kippenhahn1.png`
- **kippenhahn2** — `outputs/kippenhahn2.png`
- **gw** — `outputs/gw.png`
- **diagnostics** — `outputs/diagnostics.png`
- **orbit_frame** — `outputs/orbit_frame.png`

## 11. Bibliography

See `REFERENCES.bib` for the full list with DOIs. Key prescriptions used in this run: Peters (1964); Mora & Will (2004); Eggleton (1983); Ritter (1988) / Kolb & Ritter (1990); Hut (1981); Webbink (1984); Vink et al. (2001) / Bjorklund et al. (2021); Nugis & Lamers (2000); Fryer et al. (2012); Hobbs et al. (2005); Farmer et al. (2019); Kippenhahn, Weigert & Weiss (2012).
