"""Run identical initial conditions under competing prescriptions per toggle."""
from __future__ import annotations

import warnings

import numpy as np

from ..binary.driver import evolve
from ..io.config import load_config
from .sensitivity import outcome_vector

# competing options to sweep, one toggle at a time
TOGGLE_ALTERNATIVES = {
    "wind": ["vink2001", "bjorklund2021"],
    "convection": ["schwarzschild", "ledoux"],
    "overshoot": ["none", "step", "exponential"],
    "am_transport": ["none", "tayler_spruit"],
    "mt_stability": ["qcrit", "zeta"],
    "mt_rate": ["ritter", "kolb_ritter"],
    "ce_formalism": ["alpha_lambda", "gamma"],
    "opacity": ["electron_scattering", "kramers"],
    "reaction_rates": ["nacre", "reaclib"],
    "sn_engine": ["rapid", "delayed", "direct_collapse"],
    "gravity": ["newtonian", "2.5PN"],
    "tides": ["none", "hut1981"],
}


def compare(base_path=None, toggles=None):
    toggles = toggles or list(TOGGLE_ALTERNATIVES)
    base_cfg = load_config(base_path)
    table = []
    for toggle in toggles:
        for opt in TOGGLE_ALTERNATIVES[toggle]:
            ov = {f"theories.{toggle}": opt}
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    cfg = load_config(base_path, overrides=ov)
                    res = evolve(cfg)
                oc = outcome_vector(res)
                oc["fate"] = res.fate
                oc["wallclock_s"] = round(res.wallclock_s, 2)
            except Exception as exc:
                oc = {"fate": f"ERROR: {exc}"}
            table.append(dict(toggle=toggle, option=opt,
                              is_default=(base_cfg.theories[toggle] == opt), **oc))
    return table
