"""All figures.  One function per figure; each returns a matplotlib Figure and
writes it if `path` is given.  Captions state the physics options used.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ..dynamics import gw
from ..numerics.units import LSUN, MSUN, RSUN, PC, YEAR
from ..binary.roche import equipotential_grid, lagrange_points, roche_potential

_PHASE_COLOR = {"ZAMS": "#1f77b4", "MS": "#1f77b4", "HG": "#ff7f0e",
                "CHeB": "#2ca02c", "stripped_He": "#9467bd",
                "preSN": "#d62728", "remnant": "#111111"}


def _save(fig, path, fmt):
    if path:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=130, bbox_inches="tight")
    return fig


def hr_diagram(hist, cfg, path=None, fmt="png"):
    h = hist.as_arrays()
    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    for star, T, L, ph in ((1, h["Teff1"], h["L1"], h["phase1"]),
                           (2, h["Teff2"], h["L2"], h["phase2"])):
        cols = [_PHASE_COLOR.get(p, "#888") for p in ph]
        ax.scatter(np.log10(T), np.log10(L), c=cols, s=6, label=f"star {star}")
        ax.plot(np.log10(T), np.log10(L), color="0.6", lw=0.4, zorder=0)
        ax.annotate(f"{star}: ZAMS", (np.log10(T[0]), np.log10(L[0])), fontsize=8)
    # reference massive stars (observed; Crowther 2007 / literature)
    ref = {"O3V (~50 Msun)": (4.68, 5.8), "B0V (~15 Msun)": (4.48, 4.2),
           "RSG (Betelgeuse)": (3.56, 5.1), "WR (WN, stripped)": (4.8, 5.4)}
    for name, (lt, ll) in ref.items():
        ax.plot(lt, ll, "k*", ms=9)
        ax.annotate(name, (lt, ll), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(r"$\log_{10} T_{\rm eff}$ [K]")
    ax.set_ylabel(r"$\log_{10} L / L_\odot$")
    ax.invert_xaxis()
    ax.set_title("HR diagram -- both tracks, phase-coloured")
    ax.legend(loc="lower left")
    fig.text(0.5, -0.02,
             f"wind={cfg.theories['wind']}, conv={cfg.theories['convection']}, "
             f"overshoot={cfg.theories['overshoot']}, Z={cfg.system.metallicity:g}. "
             "Blue=MS, orange=HG, green=CHeB, purple=stripped He, black=remnant. "
             "Stars = representative observed objects.", ha="center", fontsize=7)
    return _save(fig, path, fmt)


def orbital_evolution(hist, events, cfg, path=None, fmt="png"):
    h = hist.as_arrays()
    t = np.clip(h["t"], 1e-3, None)
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.0), sharex=True)
    axes[0].plot(t, h["a"], "b-"); axes[0].set_ylabel(r"$a$ [$R_\odot$]"); axes[0].set_yscale("log")
    axes[1].plot(t, h["e"], "g-"); axes[1].set_ylabel("eccentricity")
    P = 2 * np.pi * np.sqrt((h["a"] * RSUN) ** 3 / (6.674e-8 * (h["m1"] + h["m2"]) * MSUN)) / 86400.0
    axes[2].plot(t, P, "m-"); axes[2].set_ylabel(r"$P_{\rm orb}$ [d]"); axes[2].set_yscale("log")
    axes[2].set_xlabel("time [yr]"); axes[2].set_xscale("log")
    for ax in axes:
        for ev in events:
            ty = ev.get("t_yr")
            if ty:
                ax.axvline(max(ty, 1e-3), color="r", ls=":", lw=0.8)
                ax.text(max(ty, 1e-3), ax.get_ylim()[1], ev["kind"], rotation=90,
                        fontsize=6, va="top")
    fig.tight_layout()
    axes[0].set_title(f"Orbital evolution  (gravity={cfg.theories['gravity']}, "
                      f"tides={cfg.theories['tides']}, MT={cfg.theories['mt_rate']})")
    return _save(fig, path, fmt)


def mass_evolution(hist, cfg, path=None, fmt="png"):
    h = hist.as_arrays()
    t = np.clip(h["t"], 1e-3, None)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True)
    ax1.plot(t, h["m1"], "b-", label=r"$M_1$")
    ax1.plot(t, h["m2"], "r-", label=r"$M_2$")
    ax1.plot(t, h["mcore1"], "b--", label=r"$M_{1,\rm He core}$")
    ax1.plot(t, h["mcore2"], "r--", label=r"$M_{2,\rm He core}$")
    ax1.set_ylabel(r"mass [$M_\odot$]"); ax1.legend(fontsize=8); ax1.set_xscale("log")
    ax2.plot(t, -h["mdot_wind1"], "b-", label=r"$|\dot M_{\rm wind,1}|$")
    ax2.plot(t, -h["mdot_wind2"], "r-", label=r"$|\dot M_{\rm wind,2}|$")
    ax2.plot(t, np.abs(h["mdot_rlof"]), "k-", label=r"$|\dot M_{\rm RLOF}|$")
    ax2.set_yscale("log"); ax2.set_ylabel(r"$|\dot M|$ [$M_\odot$/yr]")
    ax2.set_xlabel("time [yr]"); ax2.legend(fontsize=8)
    fig.tight_layout()
    ax1.set_title(f"Mass evolution (wind={cfg.theories['wind']}, "
                  f"eta_wind={cfg.p('eta_wind'):g}, beta_acc={cfg.p('beta_accretion'):g})")
    return _save(fig, path, fmt)


def roche_geometry(q, path=None, fmt="png", title_extra=""):
    X, Y, Phi = equipotential_grid(q)
    L = lagrange_points(q)
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    phi_L1 = roche_potential(L["L1"][0], 0, 0, q)
    levels = np.sort(phi_L1 * np.array([1.35, 1.18, 1.06, 1.0, 0.96, 0.9]))
    ax.contour(X, Y, Phi, levels=levels, colors="k", linewidths=0.7)
    ax.contour(X, Y, Phi, levels=[phi_L1], colors="r", linewidths=1.6)
    ax.plot(0, 0, "b*", ms=13); ax.plot(1, 0, "r*", ms=11)
    for name, (x, y) in L.items():
        ax.plot(x, y, "k+"); ax.annotate(name, (x, y), fontsize=8)
    ax.set_aspect("equal"); ax.set_xlabel(r"$x/a$"); ax.set_ylabel(r"$y/a$")
    ax.set_title(f"Roche geometry, q={q:.3f} {title_extra}\n(red = L1 equipotential / Roche lobes)")
    return _save(fig, path, fmt)


def kippenhahn(hist, star_idx, cfg, path=None, fmt="png"):
    """Schematic Kippenhahn: convective-core / burning regions vs time in mass coord.

    The reduced evolution model does not resolve the interior; this shows the
    fitted convective-core mass fraction and the total mass envelope, which is the
    information the model actually carries.  Marked SCHEMATIC in the caption.
    """
    h = hist.as_arrays()
    t = np.clip(h["t"], 1e-3, None)
    m = h[f"m{star_idx}"]
    mc = h[f"mcore{star_idx}"]
    ph = h[f"phase{star_idx}"]
    # fitted convective-core fraction during MS: shrinks from ~0.4 to ~0.2*M
    fconv = np.where(np.isin(ph, ["MS", "ZAMS"]), 0.35 * m, np.maximum(mc, 0.1 * m))
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.fill_between(t, 0, fconv, color="0.6", label="convective core (fit)")
    ax.fill_between(t, fconv, np.where(mc > 0, mc, np.nan), color="#ffcc66",
                    alpha=0.6, label="He core")
    ax.plot(t, m, "k-", lw=1.5, label="total mass (surface)")
    ax.set_xscale("log"); ax.set_xlabel("time [yr]")
    ax.set_ylabel(r"mass coordinate [$M_\odot$]")
    ax.set_title(f"Kippenhahn (SCHEMATIC) -- star {star_idx}, "
                 f"convection={cfg.theories['convection']}, overshoot={cfg.theories['overshoot']}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return _save(fig, path, fmt)


def gw_output(result, cfg, path=None, fmt="png"):
    fin = result.final
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.2))
    have_dco = fin.get("remnant1") and fin.get("remnant2") and np.isfinite(fin.get("a_Rsun", np.inf))
    if have_dco:
        m1 = fin["m1_Msun"] * MSUN
        m2 = fin["m2_Msun"] * MSUN
        a = fin["a_Rsun"] * RSUN
        f_orb0 = 1.0 / (2 * np.pi) * np.sqrt(6.674e-8 * (m1 + m2) / a**3)
        f_gw0 = 2 * f_orb0
        D = 100e6 * PC
        # first pass: get the Newtonian coalescence time from the current f
        wf0 = gw.inspiral_waveform(m1, m2, max(f_gw0, 1e-5), D, np.array([0.0]))
        tau = wf0["tau"]
        f_isco = gw.peak_gw_frequency_isco(m1, m2)
        # sample toward coalescence (where the chirp is visible), stop near f_ISCO
        t_grid = tau * (1.0 - np.logspace(0, -11, 6000))
        wf = gw.inspiral_waveform(m1, m2, max(f_gw0, 1e-5), D, t_grid)
        keep = wf["f"] <= f_isco
        for k in ("t", "f", "h_plus", "h_cross"):
            wf[k] = wf[k][keep]
        t_before_merger = tau - wf["t"]
        axes[0, 0].plot(t_before_merger, wf["h_plus"], lw=0.4)
        axes[0, 0].set_xscale("log"); axes[0, 0].invert_xaxis()
        axes[0, 0].set_title(f"h_+(t) chirp (0PN, D=100 Mpc); coalescence in {tau/3.15e7/1e6:.0f} Myr")
        axes[0, 0].set_xlabel("time before merger [s]")
        axes[0, 1].loglog(wf["f"], np.abs(wf["h_plus"]))
        axes[0, 1].set_title("instantaneous strain vs f_GW"); axes[0, 1].set_xlabel("f_GW [Hz]")
        axes[1, 0].loglog(t_before_merger, wf["f"])
        axes[1, 0].invert_xaxis()
        axes[1, 0].set_title("frequency sweep"); axes[1, 0].set_xlabel("time before merger [s]")
        axes[1, 0].set_ylabel("f_GW [Hz]")
    else:
        for ax in axes.flat[:3]:
            ax.text(0.5, 0.5, "no merging double-compact object formed",
                    ha="center", va="center")
    # characteristic strain vs detectors
    ax = axes[1, 1]
    f = np.logspace(-4, 4, 600)
    ax.loglog(f, gw.lisa_characteristic_strain(f), label="LISA (Robson+ 2019)")
    ax.loglog(f, gw.ligo_aligo_characteristic_strain(f), label="aLIGO (approx)")
    ax.loglog(f, gw.einstein_telescope_characteristic_strain(f), label="ET (approx)")
    if have_dco:
        fpk = gw.peak_gw_frequency_isco(m1, m2)
        fr = np.logspace(np.log10(max(f_gw0, 1e-4)), np.log10(fpk), 100)
        hc = gw.characteristic_strain_inspiral(m1, m2, fr, D)
        ax.loglog(fr, hc, "k-", lw=2, label="this system (D=100 Mpc)")
        ax.plot(fpk, gw.characteristic_strain_inspiral(m1, m2, fpk, D), "ko")
    ax.set_xlabel("f [Hz]"); ax.set_ylabel("characteristic strain")
    ax.legend(fontsize=7); ax.set_ylim(1e-24, 1e-15)
    fig.suptitle(f"Gravitational-wave output (gravity={cfg.theories['gravity']})")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return _save(fig, path, fmt)


def diagnostics(hist, cfg, path=None, fmt="png"):
    h = hist.as_arrays()
    t = np.clip(h["t"], 1e-3, None)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True)
    ax1.loglog(t, h["dt"], ".", ms=2)
    ax1.set_ylabel("timestep [yr]")
    crits = list(dict.fromkeys(h["binding_criterion"].tolist()))
    cmap = {c: i for i, c in enumerate(crits)}
    ax2.scatter(t, [cmap[c] for c in h["binding_criterion"]], s=3)
    ax2.set_yticks(range(len(crits))); ax2.set_yticklabels(crits, fontsize=7)
    ax2.set_xlabel("time [yr]"); ax2.set_ylabel("binding criterion")
    fig.tight_layout()
    ax1.set_title("Integration diagnostics: timestep history & binding criterion")
    return _save(fig, path, fmt)


def interior_profiles(model, path=None, fmt="png"):
    """rho, T, P, and nabla_rad vs nabla_ad from a 1D StellarModel."""
    from ..structure import eos
    from ..numerics.units import A_RAD, K_B, M_H, mu_from_composition
    p = model.profile
    m, r, P, L, T = p["m"], p["r"], p["P"], p["L"], p["T"]
    X, Y, Z = model.comp["X"], model.comp["Y"], model.comp["Z"]
    mu = mu_from_composition(X, Y, Z, True)
    rho = np.maximum((P - A_RAD * T**4 / 3) * mu * M_H / (K_B * T), 1e-12)
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0))
    rr = r / r[-1]
    axes[0, 0].semilogy(rr, rho); axes[0, 0].set_ylabel(r"$\rho$ [g/cm$^3$]")
    axes[0, 1].semilogy(rr, T); axes[0, 1].set_ylabel("T [K]")
    axes[1, 0].semilogy(rr, P); axes[1, 0].set_ylabel("P [dyn/cm$^2$]"); axes[1, 0].set_xlabel(r"$r/R$")
    nad = np.array([eos.grad_ad(rh, tt, mu) for rh, tt in zip(rho, T)])
    axes[1, 1].plot(rr, nad, label=r"$\nabla_{\rm ad}$")
    axes[1, 1].set_xlabel(r"$r/R$"); axes[1, 1].set_ylabel(r"$\nabla_{\rm ad}$"); axes[1, 1].legend()
    fig.suptitle(f"Interior (ZAMS 1D model) M={model.M/MSUN:.1f} Msun "
                 f"{'CONVERGED' if model.converged else 'NOT CONVERGED'}")
    return _save(fig, path, fmt)


def orbit_animation_frame(hist, cfg, frac=0.5, path=None, fmt="png"):
    """Single representative frame of the to-scale orbit + Roche lobes (static PNG;
    a full animation writer is in viz.animate)."""
    h = hist.as_arrays()
    i = int(frac * (len(h["t"]) - 1))
    a = h["a"][i]; e = h["e"][i]
    q = h["m1"][i] / h["m2"][i]
    fig, ax = plt.subplots(figsize=(6, 6))
    th = np.linspace(0, 2 * np.pi, 400)
    rr = a * (1 - e**2) / (1 + e * np.cos(th))
    ax.plot(rr * np.cos(th), rr * np.sin(th), "b-", lw=0.8)
    ax.add_patch(plt.Circle((0, 0), h["R1"][i], color="b", alpha=0.5))
    ax.add_patch(plt.Circle((a, 0), h["R2"][i], color="r", alpha=0.5))
    ax.set_aspect("equal"); ax.set_title(f"orbit @ t={h['t'][i]:.2e} yr, a={a:.1f} Rsun, stars to scale")
    ax.set_xlabel(r"$R_\odot$")
    return _save(fig, path, fmt)


def make_all(result, cfg, outdir):
    outdir = Path(outdir)
    fmt = cfg.output.figure_format
    figs = {}
    q_final = max(result.final["m1_Msun"] / result.final["m2_Msun"], 1e-3)
    plan = [
        ("hr_diagram", lambda: hr_diagram(result.history, cfg, outdir / f"hr.{fmt}", fmt)),
        ("orbital_evolution", lambda: orbital_evolution(result.history, result.events, cfg, outdir / f"orbit.{fmt}", fmt)),
        ("mass_evolution", lambda: mass_evolution(result.history, cfg, outdir / f"mass.{fmt}", fmt)),
        ("roche", lambda: roche_geometry(q_final, outdir / f"roche.{fmt}", fmt, "(final)")),
        ("kippenhahn1", lambda: kippenhahn(result.history, 1, cfg, outdir / f"kipp1.{fmt}", fmt)),
        ("kippenhahn2", lambda: kippenhahn(result.history, 2, cfg, outdir / f"kipp2.{fmt}", fmt)),
        ("gw", lambda: gw_output(result, cfg, outdir / f"gw.{fmt}", fmt)),
        ("diagnostics", lambda: diagnostics(result.history, cfg, outdir / f"diag.{fmt}", fmt)),
        ("orbit_frame", lambda: orbit_animation_frame(result.history, cfg, 0.5, outdir / f"orbit_frame.{fmt}", fmt)),
    ]
    for name, fn in plan:
        try:
            fn()
            figs[name] = str(outdir / f"{name}.{fmt}")
            plt.close("all")
        except Exception as exc:
            figs[name] = f"FAILED: {exc}"
            plt.close("all")
    return figs
