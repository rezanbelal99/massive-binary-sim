"""Coupled binary-evolution driver -- the time-integration loop.

State machine over phases:  detached -> (RLOF: stable MT | common envelope) ->
detached/stripped -> SN1 -> post-SN binary -> SN2 -> double-compact -> GW inspiral.

Every sub-process is a pluggable module; this file only orchestrates them, chooses
the timestep (numerics.timestep), and records history + provenance.

Coupling is operator-split (first-order): within a step we (1) evolve stellar
structure, (2) apply winds + wind AM loss, (3) apply RLOF or CE, (4) apply tides,
(5) apply GW radiation, (6) apply third-body secular terms.  Split error is
controlled by the timestep limiter, not by iteration -- see KNOWN_GAPS DRV-SPLIT.
"""
from __future__ import annotations

import time as _time
from dataclasses import dataclass, field

import numpy as np

from ..numerics import units as U
from ..numerics.timestep import StepLimits, choose_timestep
from ..structure.evolution import StarState
from ..winds import recipes as winds
from . import common_envelope as ce
from . import mass_transfer as mt
from . import tides
from .roche import roche_lobe_radius
from ..dynamics.gw import peters_dadt, peters_dedt, merger_time_eccentric
from ..remnants.core_collapse import fryer2012
from ..remnants.pair_instability import classify as pi_classify
from ..remnants import kicks


@dataclass
class BinaryHistory:
    t: list = field(default_factory=list)          # yr
    a: list = field(default_factory=list)          # Rsun
    e: list = field(default_factory=list)
    m1: list = field(default_factory=list)         # Msun
    m2: list = field(default_factory=list)
    R1: list = field(default_factory=list)
    R2: list = field(default_factory=list)
    L1: list = field(default_factory=list)
    L2: list = field(default_factory=list)
    Teff1: list = field(default_factory=list)
    Teff2: list = field(default_factory=list)
    phase1: list = field(default_factory=list)
    phase2: list = field(default_factory=list)
    RL1: list = field(default_factory=list)
    RL2: list = field(default_factory=list)
    mdot_wind1: list = field(default_factory=list)
    mdot_wind2: list = field(default_factory=list)
    mdot_rlof: list = field(default_factory=list)
    mcore1: list = field(default_factory=list)
    mcore2: list = field(default_factory=list)
    dt: list = field(default_factory=list)
    binding_criterion: list = field(default_factory=list)
    state: list = field(default_factory=list)

    def as_arrays(self):
        return {k: np.asarray(v) for k, v in self.__dict__.items()}


@dataclass
class BinaryResult:
    history: BinaryHistory
    events: list
    fate: str
    final: dict
    wallclock_s: float
    warnings: list


def _wind_recipe_for(cfg, star):
    r = cfg.theories["wind"]
    return r


def evolve(cfg, verbose=False):
    """Run the coupled evolution for a Config. Returns BinaryResult."""
    t0 = _time.time()
    rng = np.random.default_rng(cfg.numerics.seed)
    warn_list = []
    events = []

    s = cfg.system
    star1 = StarState.zams(s.m1, s.metallicity, spin=s.spin1)
    star2 = StarState.zams(s.m2, s.metallicity, spin=s.spin2)
    # by convention star1 is the initially more massive (primary)
    if star2.m > star1.m:
        star1, star2 = star2, star1

    a = cfg.separation_cm()
    e = s.eccentricity

    hist = BinaryHistory()
    limits = StepLimits(dt_min=cfg.numerics.min_dt * U.YEAR,
                        dt_max=cfg.numerics.max_dt * U.YEAR)

    t = 0.0
    t_max = cfg.numerics.max_time * U.YEAR
    state = "detached"
    remnant1 = remnant2 = None
    fate = "unfinished"
    n_steps = 0
    n_ce = 0
    ce_cooldown_until = -1.0
    mt_episode_start = None
    mt_step_count = 0
    total_rlof_steps = 0
    prev_donor = None

    def record(dt, crit, mdw1, mdw2, mdrl):
        hist.t.append(t / U.YEAR)
        hist.a.append(a / U.RSUN)
        hist.e.append(e)
        hist.m1.append(star1.m / U.MSUN)
        hist.m2.append(star2.m / U.MSUN)
        hist.R1.append(star1.R / U.RSUN)
        hist.R2.append(star2.R / U.RSUN)
        hist.L1.append(star1.L / U.LSUN)
        hist.L2.append(star2.L / U.LSUN)
        hist.Teff1.append(star1.Teff)
        hist.Teff2.append(star2.Teff)
        hist.phase1.append(star1.phase)
        hist.phase2.append(star2.phase)
        hist.RL1.append(roche_lobe_radius(a, star1.m, star2.m) / U.RSUN)
        hist.RL2.append(roche_lobe_radius(a, star2.m, star1.m) / U.RSUN)
        hist.mdot_wind1.append(mdw1 * U.YEAR / U.MSUN)
        hist.mdot_wind2.append(mdw2 * U.YEAR / U.MSUN)
        hist.mdot_rlof.append(mdrl * U.YEAR / U.MSUN)
        hist.mcore1.append(star1.m_he_core / U.MSUN)
        hist.mcore2.append(star2.m_he_core / U.MSUN)
        hist.dt.append(dt / U.YEAR)
        hist.binding_criterion.append(crit)
        hist.state.append(state)

    MAX_STEPS = 200000
    while t < t_max and n_steps < MAX_STEPS:
        n_steps += 1
        P_orb = U.keplerian_period(a, star1.m + star2.m)

        # ---- (1) winds -----------------------------------------------------------
        def wind_mdot(st):
            if st.phase == "remnant":
                return 0.0
            try:
                return winds.mass_loss_rate(st, cfg.theories["wind"],
                                            cfg.p("eta_wind"), cfg.p("clumping"))
            except Exception as exc:  # keep going, log
                warn_list.append(f"wind calc failed ({exc}); set 0")
                return 0.0
        mdw1 = wind_mdot(star1) if remnant1 is None else 0.0
        mdw2 = wind_mdot(star2) if remnant2 is None else 0.0

        # ---- RLOF check --------------------------------------------------------
        RL1 = roche_lobe_radius(a * (1 - e), star1.m, star2.m)
        RL2 = roche_lobe_radius(a * (1 - e), star2.m, star1.m)
        donor = acc = None
        if remnant1 is None and star1.R >= RL1 and star1.R / RL1 >= star2.R / max(RL2, 1e-30):
            donor, acc, RLd = star1, star2, RL1
        elif remnant2 is None and star2.R >= RL2:
            donor, acc, RLd = star2, star1, RL2

        mdot_rlof = 0.0
        da_mt = 0.0
        if donor is not None and state not in ("merged", "disrupted"):
            if donor is not prev_donor:
                mt_episode_start = t
                mt_step_count = 0
                prev_donor = donor
            mt_step_count += 1
            q = donor.m / acc.m
            stable, margin, note = mt.is_stable(q, donor.phase, cfg.theories["mt_stability"])
            case = mt.classify_case(donor.phase)
            m_core_est = (donor.m_he_core if donor.m_he_core > 0
                          else 0.35 * donor.m0 if donor.phase in ("MS", "ZAMS")
                          else 0.5 * donor.m)
            total_rlof_steps += 1
            episode_too_long = ((t - (mt_episode_start or t)) > 0.15 * t_max
                                or mt_step_count > 20000 or total_rlof_steps > 40000)
            envelope_gone = donor.m <= 1.20 * m_core_est

            if donor.stripped and donor.phase == "stripped_He" and donor.R > RLd:
                # an already-stripped helium star that still overflows -> merges
                # with the companion (double-He / He-star + compact merger)
                state = "merged"
                fate = (f"merger of stripped He star with "
                        f"{'compact object' if (remnant1 or remnant2) else 'companion'} "
                        f"(Case {case})")
                events.append(dict(t_yr=t / U.YEAR, kind="merger", detail=fate))
                record(limits.dt_min, "event:He-merger", mdw1, mdw2, 0.0)
                break

            if envelope_gone or episode_too_long:
                # transfer has run to completion -> donor is a stripped helium star
                donor.m = max(donor.m, m_core_est)
                donor.stripped = True
                donor.phase = "stripped_He"
                donor.m_he_core = donor.m
                donor._strip_t = t
                donor._strip_age0 = donor.age
                donor._update_observables()
                events.append(dict(t_yr=t / U.YEAR, kind="mass_transfer_end", case=case,
                                   detail=f"Case {case} transfer complete; donor -> "
                                          f"{donor.m/U.MSUN:.2f} Msun stripped He star"))
                state = "post-MT"
                prev_donor = None
            elif (not stable) and n_ce < 2 and t >= ce_cooldown_until:
                res = ce.common_envelope(
                    donor, acc.m, a, cfg.theories["ce_formalism"],
                    alpha_ce=cfg.p("alpha_ce"), lambda_ce=cfg.p("lambda_ce"),
                    gamma=cfg.p("gamma_amloss"),
                    r_companion=(acc.R if remnant2 is None and remnant1 is None else 0.0))
                events.append(dict(t_yr=t / U.YEAR, kind="common_envelope", case=case,
                                   q=q, detail=res.verdict))
                n_ce += 1
                ce_cooldown_until = t + 0.02 * t_max
                if not res.survived:
                    state = "merged"
                    fate = f"common-envelope merger (Case {case}, q={q:.2f})"
                    record(limits.dt_min, "event:CE-merger", mdw1, mdw2, 0.0)
                    break
                a = res.a_final
                e = 0.0
                donor.m = res.m_core_donor
                donor.stripped = True
                donor.phase = "stripped_He"
                donor.m_he_core = res.m_core_donor
                donor._strip_t = t
                donor._update_observables()
                state = "post-CE"
                prev_donor = None
            else:
                state = f"RLOF-Case{case}"
                mdot_rlof = mt.mdot_rlof(donor, a, acc.m, cfg.theories["mt_rate"])
                beta = cfg.p("beta_accretion")
                da_mt = mt.orbital_response(a, donor.m, acc.m, mdot_rlof, beta=beta,
                                            gamma=cfg.p("gamma_amloss"),
                                            mode="isotropic_reemission")

        # ---- pick timestep ----------------------------------------------------
        dadt_gw = peters_dadt(a, e, star1.m, star2.m)
        # crude estimate of total |da/dt| for the step limiter
        dadt_wind_est = 0.0
        for st, md, other in ((star1, mdw1, star2), (star2, mdw2, star1)):
            if md < 0 and remnant1 is None and remnant2 is None:
                dadt_wind_est += abs(a * (md / st.m) * (other.m / (st.m + other.m)))
        dadt_total = abs(da_mt) + dadt_wind_est + (abs(dadt_gw) if (remnant1 or remnant2 or a / U.RSUN < 20) else 0.0)
        st_state = dict(P_orbit=P_orb, a=a, dadt_gw=(dadt_gw if (remnant1 and remnant2) else 0.0),
                        dadt_total=dadt_total, rlof_active=(donor is not None),
                        Mdot_wind1=mdw1, Mdot_wind2=mdw2,
                        Mdot_rlof=mdot_rlof, M_donor=(donor.m if donor else star1.m))
        # only non-degenerate stars set thermal / nuclear timestep limits
        for i, st in ((1, star1), (2, star2)):
            if (remnant1 if i == 1 else remnant2) is not None:
                continue
            st_state[f"M{i}"] = st.m
            st_state[f"R{i}"] = st.R
            st_state[f"L{i}"] = max(st.L, 1.0)
            st_state[f"M_core{i}"] = st.m_he_core
            st_state[f"t_ms{i}"] = st._t_ms
            st_state[f"relaxing{i}"] = st.phase in ("HG",) or (
                st.phase == "stripped_He" and (t - getattr(st, "_strip_t", 0.0)) < 0.02 * st._t_ms)
            # rough time to the end of the current phase (paces the step)
            if st.phase in ("ZAMS", "MS"):
                st_state[f"t_to_phase_end{i}"] = max(st._t_ms - st.age, limits.dt_min)
            elif st.phase == "HG":
                st_state[f"t_to_phase_end{i}"] = max(1.02 * st._t_ms - st.age, limits.dt_min)
            elif st.phase == "CHeB":
                st_state[f"t_to_phase_end{i}"] = max(
                    1.02 * st._t_ms + 0.10 * st._t_ms - st.age, limits.dt_min)
        dt, crit = choose_timestep(st_state, limits)
        dt = min(dt, t_max - t)

        # ---- (2..6) apply over dt -------------------------------------------------
        # never let a single step's wind remove more than 15% of a star's mass
        if dt > 0:
            mdw1 = -min(abs(mdw1), 0.15 * star1.m / dt)
            mdw2 = -min(abs(mdw2), 0.15 * star2.m / dt)

        # stellar structure
        if remnant1 is None:
            star1.advance(dt, mdot_wind=mdw1,
                          mdot_binary=(mdot_rlof if donor is star1 else
                                       (-cfg.p("beta_accretion") * mdot_rlof if acc is star1 else 0.0)))
        if remnant2 is None:
            star2.advance(dt, mdot_wind=mdw2,
                          mdot_binary=(mdot_rlof if donor is star2 else
                                       (-cfg.p("beta_accretion") * mdot_rlof if acc is star2 else 0.0)))

        # wind orbital AM loss (Jeans) -> a grows
        for st, md, other in ((star1, mdw1, star2), (star2, mdw2, star1)):
            if md < 0 and remnant1 is None and remnant2 is None:
                # fast wind: da/a = -(dM/M) (Jeans, beta=1); dM = md*dt < 0 -> a grows
                factor = 1.0 - (md * dt) / st.m * (other.m / (st.m + other.m))
                a *= float(np.clip(factor, 0.8, 1.25))
        # RLOF orbital response (cap the per-step change so a runaway MT episode
        # cannot blow the orbit up before the completion check fires)
        if da_mt != 0.0:
            da_step = float(np.clip(da_mt * dt, -0.25 * a, 0.25 * a))
            a = max(a + da_step, 1.1 * (star1.R + star2.R))

        # tides
        if cfg.theories["tides"] != "none" and remnant1 is None and remnant2 is None and e > 1e-4:
            for st, other in ((star1, star2), (star2, star1)):
                env = "convective" if st.phase in ("CHeB",) and not st.stripped else "radiative"
                try:
                    da_t, de_t, dOm = tides.tidal_derivatives(st, other.m, a, e, env)
                    a = a + np.clip(da_t * dt, -0.1 * a, 0.1 * a)
                    e = float(np.clip(e + np.clip(de_t * dt, -0.1, 0.1), 0.0, 0.99))
                except Exception:
                    pass

        # GW radiation (mainly post-second-SN)
        if a > 0 and (remnant1 or remnant2 or a / U.RSUN < 20.0):
            da_g = peters_dadt(a, e, star1.m, star2.m) * dt
            de_g = peters_dedt(a, e, star1.m, star2.m) * dt
            a = max(a + da_g, 1.0e6)
            e = float(np.clip(e + de_g, 0.0, 0.999))

        t += dt

        # ---- orbit sanity: a very wide orbit will neither interact nor merge -----
        if not np.isfinite(a) or a > 2.0e16:   # ~1400 AU
            state = "wide-detached"
            r1s = f"{remnant1['kind']}" if remnant1 else star1.phase
            r2s = f"{remnant2['kind']}" if remnant2 else star2.phase
            fate = (f"wide non-interacting binary ({r1s}+{r2s}), "
                    f"a = {a/U.AU:.0f} AU -- no further interaction or GW merger")
            record(dt, "event:wide", mdw1, mdw2, mdot_rlof)
            break

        record(dt, crit, mdw1, mdw2, mdot_rlof)

        # ---- contact / MS-merger check --------------------------------------------
        # True contact = BOTH stars overflow their Roche lobes simultaneously
        # (an over-contact system on a dynamical timescale).  A single deeply
        # overflowing donor is routed through the RLOF/CE logic above, not here.
        RL1_now = roche_lobe_radius(a * (1 - e), star1.m, star2.m)
        RL2_now = roche_lobe_radius(a * (1 - e), star2.m, star1.m)
        if (remnant1 is None and remnant2 is None
                and star1.R > RL1_now and star2.R > RL2_now):
            state = "merged"
            both_ms = star1.phase in ("MS", "ZAMS") and star2.phase in ("MS", "ZAMS")
            fate = ("main-sequence merger" if both_ms else "contact merger (post-MS)")
            events.append(dict(t_yr=t / U.YEAR, kind="merger", detail=fate))
            break

        # ---- supernova of star 1 --------------------------------------------------
        if remnant1 is None and star1.is_evolved_to_remnant():
            remnant1, a, e, ev = _do_supernova(star1, star2, a, e, cfg, rng, which=1)
            events.append(ev)
            if not np.isfinite(a):
                state = "disrupted"
                fate = "binary disrupted by first SN"
                record(limits.dt_min, "event:SN1-disrupt", 0, 0, 0)
                break
            star1.phase = "remnant"
            state = "post-SN1"

        if remnant2 is None and star2.is_evolved_to_remnant():
            remnant2, a, e, ev = _do_supernova(star2, star1, a, e, cfg, rng, which=2)
            events.append(ev)
            if not np.isfinite(a):
                state = "disrupted"
                fate = "binary disrupted by second SN"
                record(limits.dt_min, "event:SN2-disrupt", 0, 0, 0)
                break
            star2.phase = "remnant"
            state = "post-SN2"

        # ---- double compact object: check GW merger within Hubble time -----------
        if remnant1 is not None and remnant2 is not None:
            state = "double-compact"
            t_gw = merger_time_eccentric(a, e, star1.m, star2.m)
            HUBBLE = 13.8e9 * U.YEAR
            if t + t_gw < HUBBLE or a / U.RSUN < 0.05:
                fate = (f"double compact binary ({remnant1['kind']}+{remnant2['kind']}) "
                        f"merging via GW in {t_gw/U.YEAR/1e6:.1f} Myr")
                events.append(dict(t_yr=t / U.YEAR, kind="dco_formed",
                                   detail=fate, t_merge_yr=t_gw / U.YEAR))
            else:
                fate = (f"wide double compact binary ({remnant1['kind']}+{remnant2['kind']}), "
                        f"GW merger time {t_gw/U.YEAR/1e9:.1f} Gyr > Hubble")
                events.append(dict(t_yr=t / U.YEAR, kind="dco_formed", detail=fate,
                                   t_merge_yr=t_gw / U.YEAR))
            record(dt, "event:DCO", 0, 0, 0)
            break

    if fate == "unfinished":
        if n_steps >= MAX_STEPS:
            fate = "stopped at MAX_STEPS"
            warn_list.append("hit MAX_STEPS before a terminal state")
        else:
            fate = f"reached max_time ({cfg.numerics.max_time:.1e} yr) in state '{state}'"

    final = _final_dict(star1, star2, a, e, remnant1, remnant2, state)
    return BinaryResult(history=hist, events=events, fate=fate, final=final,
                        wallclock_s=_time.time() - t0, warnings=warn_list)


def _do_supernova(star, companion, a, e, cfg, rng, which):
    M = star.m / U.MSUN
    M_CO = (star.m_co_core / U.MSUN) if star.m_co_core > 0 else max(0.3 * M, 1.5)
    M_He = (star.m_he_core / U.MSUN) if star.m_he_core > 0 else max(0.5 * M, 2.0)

    pi = pi_classify(M_He, cfg.p("c12_ag_rate"))
    if pi.regime == "PISN":
        ev = dict(t_yr=None, kind=f"SN{which}", engine="PISN",
                  detail=f"star {which}: {pi.note}; no remnant, orbit disrupted")
        return None, np.inf, np.nan, ev
    if pi.regime in ("PPISN", "photodisintegration"):
        m_rem_grav = pi.remnant_bh_mass
        rem = dict(kind="BH", m_gravitational=m_rem_grav, m_baryonic=m_rem_grav / 0.9,
                   f_fallback=1.0, engine=pi.regime)
    else:
        fr = fryer2012(M, M_CO, cfg.theories["sn_engine"])
        rem = dict(kind=fr.kind, m_gravitational=fr.m_gravitational,
                   m_baryonic=fr.m_baryonic, f_fallback=fr.f_fallback, engine=fr.engine)

    m_post = rem["m_gravitational"] * U.MSUN
    sigma = cfg.p("kick_sigma")
    vk = kicks.draw_kick(rng, sigma_kms=sigma, f_fallback=rem["f_fallback"],
                         scale_by_fallback=True)
    orb = kicks.apply_sn_to_orbit(a, e, star.m, companion.m, m_post, vk, rng)
    star.m = m_post
    star.m0 = m_post
    ev = dict(t_yr=None, kind=f"SN{which}", engine=rem["engine"], remnant=rem["kind"],
              m_remnant=rem["m_gravitational"], kick_kms=orb.get("kick_speed", 0) / U.KMS,
              v_system_kms=orb.get("v_system", np.nan) / U.KMS if orb["bound"] else np.nan,
              bound=orb["bound"],
              detail=(f"star {which} -> {rem['kind']} {rem['m_gravitational']:.2f} Msun "
                      f"({rem['engine']}); " + orb["note"]))
    if not orb["bound"]:
        return None, np.inf, np.nan, ev
    return rem, orb["a_f"], orb["e_f"], ev


def _final_dict(star1, star2, a, e, r1, r2, state):
    from ..dynamics.gw import chirp_mass, merger_time_eccentric, peak_gw_frequency_isco
    d = dict(state=state,
             a_Rsun=a / U.RSUN if np.isfinite(a) else np.inf,
             period_day=(U.keplerian_period(a, star1.m + star2.m) / U.DAY
                         if np.isfinite(a) else np.inf),
             eccentricity=e,
             m1_Msun=star1.m / U.MSUN, m2_Msun=star2.m / U.MSUN,
             phase1=star1.phase, phase2=star2.phase,
             remnant1=r1, remnant2=r2)
    if r1 and r2 and np.isfinite(a):
        mc = chirp_mass(star1.m, star2.m)
        d["chirp_mass_Msun"] = mc / U.MSUN
        d["t_merge_yr"] = merger_time_eccentric(a, e, star1.m, star2.m) / U.YEAR
        d["f_gw_peak_Hz"] = peak_gw_frequency_isco(star1.m, star2.m)
    return d
