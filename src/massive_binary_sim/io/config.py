"""Configuration parsing, validation and provenance tracking.

The config file is the *only* place physical choices and free parameters live.
Nothing downstream is allowed to hardcode a number that appears here.
"""
from __future__ import annotations

import copy
import dataclasses
import importlib.resources as _res
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from ..numerics import units as U

# ---- allowed vocabulary for each theory toggle -------------------------------------
THEORY_OPTIONS: dict[str, set[str]] = {
    "gravity": {"newtonian", "1PN", "2PN", "2.5PN", "3.5PN", "mond"},
    "wind": {"vink2001", "bjorklund2021", "nugis2000", "dejager1988", "auto"},
    "convection": {"schwarzschild", "ledoux"},
    "overshoot": {"none", "step", "exponential"},
    "am_transport": {"none", "diffusive", "tayler_spruit"},
    "mt_stability": {"qcrit", "zeta"},
    "mt_rate": {"ritter", "kolb_ritter"},
    "ce_formalism": {"alpha_lambda", "gamma"},
    "opacity": {"electron_scattering", "kramers", "opal"},
    "reaction_rates": {"nacre", "reaclib"},
    "sn_engine": {"rapid", "delayed", "direct_collapse"},
    "tides": {"none", "hut1981", "zahn"},
    "structure_fidelity": {"semi_analytic", "polytrope", "onedim"},
}


@dataclass
class FreeParam:
    name: str
    value: float
    range: tuple[float, float]
    ref: str

    def clipped(self, x: float) -> float:
        lo, hi = self.range
        return min(max(x, lo), hi)

    @property
    def unverified(self) -> bool:
        return "UNVERIFIED" in self.ref.upper()


@dataclass
class SystemConfig:
    m1: float
    m2: float
    separation: float | None
    period: float | None
    eccentricity: float
    metallicity: float
    spin1: float
    spin2: float
    max_age: float | None


@dataclass
class NumericsConfig:
    orbit_integrator: str
    secular_integrator: str
    network_integrator: str
    rtol: float
    atol: float
    max_time: float
    energy_drift_tol: float
    am_drift_tol: float
    min_dt: float
    max_dt: float
    seed: int


@dataclass
class OutputConfig:
    dir: str
    hdf5_name: str
    make_figures: bool
    figure_format: str


@dataclass
class Config:
    system: SystemConfig
    theories: dict[str, str]
    free: dict[str, FreeParam]
    numerics: NumericsConfig
    output: OutputConfig
    raw: dict[str, Any] = field(repr=False, default_factory=dict)

    # -- convenience -------------------------------------------------------------
    def p(self, name: str) -> float:
        """Return the (validated, clipped-to-range) value of a free parameter."""
        fp = self.free[name]
        return fp.clipped(fp.value)

    def separation_cm(self) -> float:
        if self.system.period is not None:
            mtot = (self.system.m1 + self.system.m2) * U.MSUN
            return U.semimajor_from_period(self.system.period * U.DAY, mtot)
        return self.system.separation * U.RSUN

    def period_days(self) -> float:
        mtot = (self.system.m1 + self.system.m2) * U.MSUN
        return U.keplerian_period(self.separation_cm(), mtot) / U.DAY

    def unverified_params(self) -> list[str]:
        return [k for k, v in self.free.items() if v.unverified]

    def to_dict(self) -> dict[str, Any]:
        d = {
            "system": dataclasses.asdict(self.system),
            "theories": dict(self.theories),
            "free_parameters": {
                k: {"value": v.value, "range": list(v.range), "ref": v.ref}
                for k, v in self.free.items()
            },
            "numerics": dataclasses.asdict(self.numerics),
            "output": dataclasses.asdict(self.output),
        }
        return d


class ConfigError(ValueError):
    pass


def _validate(cfg: Config) -> None:
    s = cfg.system
    for label, m in (("m1", s.m1), ("m2", s.m2)):
        if not (0.5 <= m <= 300.0):
            raise ConfigError(f"{label}={m} M_sun outside supported 0.5-300 range")
        if not (8.0 <= m <= 150.0):
            # allowed but warned by the driver; massive-star regime is 8-150
            pass
    if s.separation is None and s.period is None:
        raise ConfigError("system: give either `separation` (R_sun) or `period` (days)")
    if not (0.0 <= s.eccentricity < 1.0):
        raise ConfigError("eccentricity must be in [0, 1)")
    if not (1e-4 <= s.metallicity <= 0.05):
        raise ConfigError("metallicity (mass fraction Z) must be in [1e-4, 0.05]")

    for key, allowed in THEORY_OPTIONS.items():
        if key not in cfg.theories:
            raise ConfigError(f"theories: missing key '{key}'")
        if cfg.theories[key] not in allowed:
            raise ConfigError(
                f"theories.{key}='{cfg.theories[key]}' not in {sorted(allowed)}"
            )

    for name, fp in cfg.free.items():
        lo, hi = fp.range
        if lo > hi:
            raise ConfigError(f"free_parameters.{name}: range {fp.range} inverted")
        if not (lo <= fp.value <= hi):
            # not fatal: driver clips and records, but flag loudly
            fp.ref = fp.ref + f"  [WARNING: default {fp.value} outside range {fp.range}]"

    n = cfg.numerics
    if n.rtol <= 0 or n.atol <= 0:
        raise ConfigError("numerics: rtol/atol must be positive")
    if n.min_dt >= n.max_dt:
        raise ConfigError("numerics: min_dt must be < max_dt")


def default_config_path() -> Path:
    # repo-root config.yaml, two levels up from this file's package dir
    here = Path(__file__).resolve()
    return here.parents[3] / "config.yaml"


def load_config(path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> Config:
    """Load a YAML config, apply dotted-key overrides, validate, return a Config."""
    path = Path(path) if path is not None else default_config_path()
    with open(path) as fh:
        raw = yaml.safe_load(fh)

    raw = copy.deepcopy(raw)
    if overrides:
        for dotted, val in overrides.items():
            _apply_override(raw, dotted, val)

    sys_raw = raw["system"]
    system = SystemConfig(
        m1=float(sys_raw["m1"]),
        m2=float(sys_raw["m2"]),
        separation=(None if sys_raw.get("separation") is None else float(sys_raw["separation"])),
        period=(None if sys_raw.get("period") is None else float(sys_raw["period"])),
        eccentricity=float(sys_raw["eccentricity"]),
        metallicity=float(sys_raw["metallicity"]),
        spin1=float(sys_raw.get("spin1", 0.0)),
        spin2=float(sys_raw.get("spin2", 0.0)),
        max_age=(None if sys_raw.get("max_age") is None else float(sys_raw["max_age"])),
    )

    free = {}
    for name, d in raw["free_parameters"].items():
        rng = tuple(float(x) for x in d["range"])
        free[name] = FreeParam(name=name, value=float(d["value"]), range=(rng[0], rng[1]),
                               ref=str(d.get("ref", d.get("note", "UNVERIFIED"))))

    numerics = NumericsConfig(
        orbit_integrator=str(raw["numerics"]["orbit_integrator"]),
        secular_integrator=str(raw["numerics"]["secular_integrator"]),
        network_integrator=str(raw["numerics"]["network_integrator"]),
        rtol=float(raw["numerics"]["rtol"]),
        atol=float(raw["numerics"]["atol"]),
        max_time=float(raw["numerics"]["max_time"]),
        energy_drift_tol=float(raw["numerics"]["energy_drift_tol"]),
        am_drift_tol=float(raw["numerics"]["am_drift_tol"]),
        min_dt=float(raw["numerics"]["min_dt"]),
        max_dt=float(raw["numerics"]["max_dt"]),
        seed=int(raw["numerics"]["seed"]),
    )
    out = OutputConfig(
        dir=str(raw["output"]["dir"]),
        hdf5_name=str(raw["output"]["hdf5_name"]),
        make_figures=bool(raw["output"]["make_figures"]),
        figure_format=str(raw["output"].get("figure_format", "png")),
    )

    cfg = Config(system=system, theories=dict(raw["theories"]), free=free,
                 numerics=numerics, output=out, raw=raw)
    _validate(cfg)
    return cfg


def _apply_override(d: dict, dotted: str, val: Any) -> None:
    keys = dotted.split(".")
    node = d
    for k in keys[:-1]:
        node = node.setdefault(k, {})
    leaf = keys[-1]
    # allow "free_parameters.alpha_ce=2.0" as a shorthand for its .value
    if keys[0] == "free_parameters" and isinstance(node.get(leaf), dict):
        node[leaf]["value"] = _coerce(val)
    else:
        node[leaf] = _coerce(val)


def _coerce(v: Any) -> Any:
    if not isinstance(v, str):
        return v
    low = v.lower()
    if low in ("true", "false"):
        return low == "true"
    if low in ("null", "~"):          # NB: "none" is a valid theory-option string
        return None
    try:
        if any(c in v for c in ".eE") and low not in ("inf",):
            return float(v)
        return int(v)
    except ValueError:
        return v
