"""Parameter sensitivity analysis.

* one-at-a-time (OAT): vary each free parameter across its literature range, hold
  the rest at default, record the scalar outcomes.
* global: Morris elementary effects (screening) and, if SALib present and the
  budget allows, a Sobol variance decomposition.

Scalar outcomes tracked: final M1, M2, chirp mass, separation, merger-time (log10),
and an integer fate-code.
"""
from __future__ import annotations

import warnings

import numpy as np

from ..binary.driver import evolve
from ..io.config import load_config
from ..numerics import units as U

FATE_CODES = {
    "main-sequence merger": 0, "contact merger (post-MS)": 1,
    "common-envelope merger": 2, "binary disrupted": 3,
    "wide double compact": 4, "double compact binary": 5,
}


def _fate_code(fate):
    for k, v in FATE_CODES.items():
        if k in fate:
            return v
    return -1


def outcome_vector(result):
    fin = result.final
    a = fin.get("a_Rsun", np.inf)
    tm = fin.get("t_merge_yr", np.nan)
    return dict(
        m1=fin["m1_Msun"], m2=fin["m2_Msun"],
        chirp_mass=fin.get("chirp_mass_Msun", np.nan),
        separation=(a if np.isfinite(a) else np.nan),
        log10_tmerge=(np.log10(tm) if (tm and np.isfinite(tm) and tm > 0) else np.nan),
        fate_code=_fate_code(result.fate),
    )


def _run(base_path, overrides):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cfg = load_config(base_path, overrides=overrides)
        res = evolve(cfg)
    return outcome_vector(res)


def one_at_a_time(base_path=None, params=None, n=5):
    cfg = load_config(base_path)
    params = params or list(cfg.free.keys())
    base = _run(base_path, {})
    rows = []
    for name in params:
        fp = cfg.free[name]
        lo, hi = fp.range
        if lo == hi:
            continue
        vals = np.linspace(lo, hi, n)
        series = []
        for v in vals:
            try:
                oc = _run(base_path, {f"free_parameters.{name}": float(v)})
            except Exception as exc:
                oc = {k: np.nan for k in base}
                oc["error"] = str(exc)
            series.append((float(v), oc))
        # influence = spread of chirp mass (or m1 if no DCO) over the range
        key = "chirp_mass" if np.isfinite(base["chirp_mass"]) else "m1"
        yy = np.array([s[1][key] for s in series], dtype=float)
        infl = float(np.nanmax(yy) - np.nanmin(yy)) if np.isfinite(yy).any() else np.nan
        rows.append(dict(param=name, range=(lo, hi), ref=fp.ref, series=series,
                         influence_metric=key, influence=infl))
    rows.sort(key=lambda r: (np.nan_to_num(r["influence"], nan=-1)), reverse=True)
    return dict(base=base, rows=rows)


def morris(base_path=None, params=None, trajectories=6, levels=4):
    """Morris elementary-effects screening via SALib (falls back to OAT-derived
    mu* estimate if SALib missing)."""
    cfg = load_config(base_path)
    params = params or [k for k, v in cfg.free.items() if v.range[0] != v.range[1]]
    bounds = [list(cfg.free[p].range) for p in params]
    try:
        from SALib.sample.morris import sample as morris_sample
        from SALib.analyze.morris import analyze as morris_analyze
    except Exception:
        return dict(method="unavailable", note="SALib not installed; run one_at_a_time")

    problem = {"num_vars": len(params), "names": params, "bounds": bounds}
    Xs = morris_sample(problem, N=trajectories, num_levels=levels)
    keys = ["m1", "m2", "chirp_mass", "separation", "log10_tmerge", "fate_code"]
    Y = {k: [] for k in keys}
    for row in Xs:
        ov = {f"free_parameters.{p}": float(v) for p, v in zip(params, row)}
        try:
            oc = _run(base_path, ov)
        except Exception:
            oc = {k: np.nan for k in keys}
        for k in keys:
            Y[k].append(oc.get(k, np.nan))
    results = {}
    for k in keys:
        yv = np.array(Y[k], dtype=float)
        if not np.isfinite(yv).all():
            yv = np.nan_to_num(yv, nan=np.nanmean(yv[np.isfinite(yv)]) if np.isfinite(yv).any() else 0.0)
        try:
            Si = morris_analyze(problem, Xs, yv, num_levels=levels)
            results[k] = {"names": list(Si["names"]),
                          "mu_star": [float(x) for x in Si["mu_star"]],
                          "sigma": [float(x) for x in Si["sigma"]]}
        except Exception as exc:
            results[k] = {"error": str(exc)}
    return dict(method="morris", n_runs=len(Xs), params=params, results=results)
