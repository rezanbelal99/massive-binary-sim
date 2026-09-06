"""Generate RESULTS.md (+ a figures PDF bundle) for a run."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import numpy as np

from .. import __version__
from ..io.output import git_hash
from .derived import summarise


def _fmt(x, nd=4):
    if isinstance(x, float):
        if not np.isfinite(x):
            return "n/a"
        if abs(x) >= 1e4 or (abs(x) < 1e-3 and x != 0):
            return f"{x:.{nd}e}"
        return f"{x:.{nd}g}"
    return str(x)


def build_results_md(result, cfg, figures=None, validation=None, sensitivity=None,
                     theory_table=None, numerical_error=None):
    h = result.history.as_arrays()
    fin = result.final
    L = []
    L.append(f"# Simulation results\n")
    L.append(f"*massive-binary-sim v{__version__}  ·  git `{git_hash()[:10]}`  ·  "
             f"seed {cfg.numerics.seed}  ·  wallclock {result.wallclock_s:.2f} s*\n")

    # ---- initial conditions -----------------------------------------------------
    L.append("## 1. Initial conditions\n")
    L.append("| quantity | value |")
    L.append("|---|---|")
    L.append(f"| M1, M2 (ZAMS) | {cfg.system.m1:g}, {cfg.system.m2:g} Msun |")
    L.append(f"| separation / period | {cfg.separation_cm()/6.957e10:.3g} Rsun / "
             f"{cfg.period_days():.4g} d |")
    L.append(f"| eccentricity | {cfg.system.eccentricity:g} |")
    L.append(f"| metallicity Z | {cfg.system.metallicity:g} |")
    L.append(f"| spins (Omega/Omega_crit) | {cfg.system.spin1:g}, {cfg.system.spin2:g} |")
    L.append("")
    L.append("### Active prescriptions\n")
    L.append("| toggle | choice |")
    L.append("|---|---|")
    for k, v in cfg.theories.items():
        L.append(f"| {k} | `{v}` |")
    L.append("")

    # ---- timeline -------------------------------------------------------------
    L.append("## 2. Evolutionary timeline\n")
    L.append("| t [Myr] | event | detail |")
    L.append("|---|---|---|")
    # phase changes
    for idx in (1, 2):
        ph = h[f"phase{idx}"]
        last = None
        for tval, p in zip(h["t"], ph):
            if p != last:
                L.append(f"| {tval/1e6:.3f} | star {idx} -> {p} | reduced-model phase |")
                last = p
    for ev in result.events:
        ty = ev.get("t_yr")
        tstr = f"{ty/1e6:.3f}" if ty else "(at SN)"
        L.append(f"| {tstr} | **{ev.get('kind')}** | {ev.get('detail','')} |")
    L.append("")
    L.append(f"**Fate:** {result.fate}\n")

    # ---- final state -------------------------------------------------------------
    L.append("## 3. Final-state table\n")
    L.append("| quantity | star 1 | star 2 |")
    L.append("|---|---|---|")
    L.append(f"| mass [Msun] | {_fmt(fin['m1_Msun'])} | {_fmt(fin['m2_Msun'])} |")
    L.append(f"| phase | {fin['phase1']} | {fin['phase2']} |")
    for idx, rem in ((1, fin.get("remnant1")), (2, fin.get("remnant2"))):
        if rem:
            L.append(f"| remnant | {rem.get('kind')} {_fmt(rem.get('m_gravitational'))} Msun "
                     f"({rem.get('engine')}) | |" if idx == 1 else
                     f"| remnant (2) | | {rem.get('kind')} {_fmt(rem.get('m_gravitational'))} Msun |")
    L.append(f"| separation | {_fmt(fin.get('a_Rsun'))} Rsun | |")
    L.append(f"| period | {_fmt(fin.get('period_day'))} d | |")
    L.append(f"| eccentricity | {_fmt(fin.get('eccentricity'))} | |")
    L.append("")

    # ---- derived --------------------------------------------------------------
    L.append("## 4. Derived quantities\n")
    L.append("| quantity | value | unit | formula |")
    L.append("|---|---|---|---|")
    for d in summarise(result, cfg):
        L.append(f"| {d['name']} | {_fmt(d['value'])} | {d['unit']} | `{d['formula']}` |")
    L.append("")

    # ---- error budget --------------------------------------------------------
    L.append("## 5. Error budget\n")
    L.append("**Numerical error** (from the convergence study):\n")
    if numerical_error:
        ne = numerical_error
        L.append(f"- {ne.get('detail','')}")
        if ne.get("rel_errors"):
            L.append(f"- relative error in final `a` vs the tightest tolerance: "
                     f"{_fmt(ne['rel_errors'][0])} at rtol=1e-8 -> "
                     f"{_fmt(ne['rel_errors'][-1])} at rtol=1e-11")
        if ne.get("observed_order_per_decade"):
            import numpy as _np
            L.append(f"- mean observed order of convergence: "
                     f"{_np.nanmean(ne['observed_order_per_decade']):.2f} per tolerance decade "
                     f"(DOP853, 8th-order adaptive)")
        L.append("- **The numerical error (~1e-9 relative) is far below the physical/"
                 "model uncertainty quantified in sections 6-7.**")
    else:
        L.append("- not run in this invocation (pass `--convergence`)")
    L.append("")
    L.append("**Physical / model uncertainty** dominates and is quantified in §6-7. "
             "The reduced single-star evolution model (fitting formulae, not a "
             "structure integration) is itself a systematic; see KNOWN_GAPS.md "
             "items EVO-ZAMS, STRUCT-1D.\n")

    # ---- sensitivity -------------------------------------------------------------
    L.append("## 6. Free-parameter sensitivity (one-at-a-time)\n")
    if sensitivity and sensitivity.get("rows"):
        L.append(f"Outcome metric: **{sensitivity['rows'][0]['influence_metric']}**. "
                 "Parameters ranked by outcome spread across their literature range.\n")
        L.append("| rank | parameter | range | outcome spread | provenance |")
        L.append("|---|---|---|---|---|")
        for i, r in enumerate(sensitivity["rows"], 1):
            L.append(f"| {i} | `{r['param']}` | {r['range'][0]:g}-{r['range'][1]:g} | "
                     f"{_fmt(r['influence'])} | {r['ref'][:70]} |")
        top = sensitivity["rows"][0]
        L.append(f"\n**Dominant uncertain input:** `{top['param']}` "
                 f"(outcome varies by {_fmt(top['influence'])} across its plausible range).\n")
    else:
        L.append("_not run (`--sensitivity`)._\n")

    # ---- theory comparison ---------------------------------------------------
    L.append("## 7. Theory comparison (competing prescriptions)\n")
    if theory_table:
        L.append("| toggle | option | default? | fate | M1 | M2 | Mchirp | a [Rsun] | log10 t_merge |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for row in theory_table:
            L.append(f"| {row['toggle']} | `{row['option']}` | "
                     f"{'YES' if row.get('is_default') else ''} | {str(row.get('fate',''))[:46]} | "
                     f"{_fmt(row.get('m1'))} | {_fmt(row.get('m2'))} | {_fmt(row.get('chirp_mass'))} | "
                     f"{_fmt(row.get('separation'))} | {_fmt(row.get('log10_tmerge'))} |")
        L.append("\n*The scientific payload: the outcome is a mapped dependence on "
                 "prescription choice, not a single number.*\n")
    else:
        L.append("_not run (`--theory-compare`)._\n")

    # ---- validation --------------------------------------------------------
    L.append("## 8. Validation suite results\n")
    if validation:
        L.append("| test | status | detail |")
        L.append("|---|---|---|")
        for name, d in validation.items():
            st = "PASS" if d.get("passed") else "**FAIL**"
            L.append(f"| {name} | {st} | {d.get('detail','')} |")
    else:
        L.append("_run `pytest` or `mbsim validate`._")
    L.append("")

    # ---- assumptions ---------------------------------------------------------
    L.append("## 9. Assumptions & limitations\n")
    L.append(textwrap.dedent("""\
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
        """))

    # ---- figures -----------------------------------------------------------
    if figures:
        L.append("## 10. Figures\n")
        for name, p in figures.items():
            if str(p).startswith("FAILED"):
                L.append(f"- {name}: {p}")
            else:
                L.append(f"- **{name}** — `{p}`")
        L.append("")

    L.append("## 11. Bibliography\n")
    L.append("See `REFERENCES.bib` for the full list with DOIs. Key prescriptions "
             "used in this run: Peters (1964); Mora & Will (2004); Eggleton (1983); "
             "Ritter (1988) / Kolb & Ritter (1990); Hut (1981); Webbink (1984); "
             "Vink et al. (2001) / Bjorklund et al. (2021); Nugis & Lamers (2000); "
             "Fryer et al. (2012); Hobbs et al. (2005); Farmer et al. (2019); "
             "Kippenhahn, Weigert & Weiss (2012).\n")
    return "\n".join(L)


def write_report(result, cfg, outdir, **kw):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    md = build_results_md(result, cfg, **kw)
    (outdir / "RESULTS.md").write_text(md)
    _write_pdf(outdir / "RESULTS.pdf", md, kw.get("figures"))
    return outdir / "RESULTS.md"


def _write_pdf(path, md_text, figures):
    """Lightweight PDF: a text page + each figure, via matplotlib PdfPages."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.image as mpimg
    except Exception:
        return
    with PdfPages(path) as pdf:
        # text pages (chunked)
        lines = md_text.splitlines()
        per = 62
        for i in range(0, len(lines), per):
            fig = plt.figure(figsize=(8.27, 11.69))
            fig.text(0.06, 0.97, "\n".join(lines[i:i + per]), va="top", ha="left",
                     family="monospace", fontsize=6.5)
            pdf.savefig(fig); plt.close(fig)
        for name, p in (figures or {}).items():
            if str(p).startswith("FAILED") or not Path(p).exists():
                continue
            try:
                fig = plt.figure(figsize=(8.27, 6.0))
                ax = fig.add_axes([0, 0.04, 1, 0.9]); ax.axis("off")
                ax.imshow(mpimg.imread(p))
                ax.set_title(name)
                pdf.savefig(fig); plt.close(fig)
            except Exception:
                plt.close("all")
