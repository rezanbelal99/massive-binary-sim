"""Command-line interface:  mbsim <command> [options]

    mbsim run        [--config C] [--set k=v ...] [--sensitivity] [--theory-compare]
                     [--convergence] [--no-figures] [--outdir D]
    mbsim validate   [--full]
    mbsim compare    [--config C]
    mbsim sensitivity[--config C] [--morris]
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

from .io.config import load_config
from .io.output import write_hdf5


def _parse_sets(pairs):
    out = {}
    for p in pairs or []:
        if "=" not in p:
            raise SystemExit(f"--set expects key=value, got {p!r}")
        k, v = p.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def cmd_run(args):
    from .binary.driver import evolve
    from .viz import figures as figmod
    from .analysis import report as repmod
    from .analysis.validation import run_all
    from .analysis.sensitivity import one_at_a_time
    from .analysis.theory_compare import compare
    from .analysis.convergence import orbit_convergence

    overrides = _parse_sets(args.set)
    cfg = load_config(args.config, overrides=overrides)
    outdir = Path(args.outdir or cfg.output.dir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[run] M1={cfg.system.m1} M2={cfg.system.m2} a={cfg.separation_cm()/6.957e10:.2f} Rsun "
          f"P={cfg.period_days():.3f} d  Z={cfg.system.metallicity}")
    if cfg.unverified_params():
        print(f"[run] UNVERIFIED free params: {cfg.unverified_params()}")

    result = evolve(cfg, verbose=True)
    print(f"[run] fate: {result.fate}   ({result.wallclock_s:.2f} s, "
          f"{len(result.history.t)} steps)")
    for w in result.warnings[:10]:
        print(f"[warn] {w}")

    figures = None
    if not args.no_figures and cfg.output.make_figures:
        figures = figmod.make_all(result, cfg, outdir)
        print(f"[run] wrote {sum(1 for v in figures.values() if not str(v).startswith('FAILED'))} figures")

    validation = run_all(fast=True)
    npass = sum(1 for d in validation.values() if d["passed"])
    print(f"[run] validation: {npass}/{len(validation)} passed")

    sens = one_at_a_time(args.config) if args.sensitivity else None
    theory = compare(args.config) if args.theory_compare else None
    conv = orbit_convergence() if args.convergence else None

    hdf5_path = outdir / cfg.output.hdf5_name
    write_hdf5(hdf5_path, result, cfg, validation_results=validation)
    print(f"[run] wrote {hdf5_path}")

    repmod.write_report(result, cfg, outdir, figures=figures, validation=validation,
                        sensitivity=sens, theory_table=theory, numerical_error=conv)
    print(f"[run] wrote {outdir/'RESULTS.md'} and RESULTS.pdf")


def cmd_validate(args):
    from .analysis.validation import run_all
    res = run_all(fast=not args.full)
    ok = True
    for name, d in res.items():
        flag = "PASS" if d["passed"] else "FAIL"
        ok &= d["passed"]
        print(f"[{flag}] {name}: {d['detail']}")
    sys.exit(0 if ok else 1)


def cmd_compare(args):
    from .analysis.theory_compare import compare
    table = compare(args.config)
    print(json.dumps(table, indent=2, default=str))


def cmd_sensitivity(args):
    from .analysis.sensitivity import one_at_a_time, morris
    if args.morris:
        print(json.dumps(morris(args.config), indent=2, default=str))
    else:
        s = one_at_a_time(args.config)
        for i, r in enumerate(s["rows"], 1):
            print(f"{i:2d}. {r['param']:16s} influence({r['influence_metric']})={r['influence']:.4g}")


def main(argv=None):
    warnings.filterwarnings("ignore")
    p = argparse.ArgumentParser(prog="mbsim", description="massive-binary-sim")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run")
    r.add_argument("--config"); r.add_argument("--set", nargs="*")
    r.add_argument("--outdir")
    r.add_argument("--sensitivity", action="store_true")
    r.add_argument("--theory-compare", action="store_true")
    r.add_argument("--convergence", action="store_true")
    r.add_argument("--no-figures", action="store_true")
    r.set_defaults(func=cmd_run)

    v = sub.add_parser("validate")
    v.add_argument("--full", action="store_true")
    v.set_defaults(func=cmd_validate)

    c = sub.add_parser("compare")
    c.add_argument("--config")
    c.set_defaults(func=cmd_compare)

    s = sub.add_parser("sensitivity")
    s.add_argument("--config"); s.add_argument("--morris", action="store_true")
    s.set_defaults(func=cmd_sensitivity)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
