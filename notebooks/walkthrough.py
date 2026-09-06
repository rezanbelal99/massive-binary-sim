"""Guided walkthrough -- reproduces the headline result and the theory-comparison payload.

Run as a script:      python notebooks/walkthrough.py
or convert to .ipynb:  jupytext --to notebook notebooks/walkthrough.py
                       (a pre-built walkthrough.ipynb is also committed)

Each '# %%' marks a notebook cell.
"""

# %% [markdown]
# # A gravitationally bound pair of massive stars
#
# We evolve the default system (35 + 25 M_sun, 50 R_sun, Z = 0.0142) with the
# default prescriptions, validate the framework, then map how the outcome depends
# on the physics choices.

# %%
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib.pyplot as plt

from massive_binary_sim import load_config
from massive_binary_sim.binary.driver import evolve
from massive_binary_sim.analysis import validation, derived, sensitivity, theory_compare
from massive_binary_sim.viz import figures

cfg = load_config()
print("system:", cfg.system)
print("prescriptions:", cfg.theories)
print("UNVERIFIED free params:", cfg.unverified_params() or "none")

# %% [markdown]
# ## 1. Validation first -- do not trust the science until these pass

# %%
vres = validation.run_all(fast=True)
for name, d in vres.items():
    print(f"[{'PASS' if d['passed'] else 'FAIL'}] {name}: {d['detail']}")
assert all(d["passed"] for d in vres.values()), "validation failed -- stop here"

# %% [markdown]
# The Hulse-Taylor check reproduces the observed orbital decay of PSR B1913+16 to
# ~0.02 %, and the Lane-Emden solver matches the analytic n = 0, 1, 5 polytropes
# to machine precision.

# %% [markdown]
# ## 2. Evolve the default system

# %%
result = evolve(cfg)
print("FATE:", result.fate)
print(f"wall-clock {result.wallclock_s:.2f} s, {len(result.history.t)} steps")
for ev in result.events:
    print("  -", ev["kind"], ":", ev.get("detail", ""))

# %%
for d in derived.summarise(result, cfg):
    print(f"{d['name']:38s} = {d['value']:.4g} {d['unit']:12s}  [{d['formula']}]")

# %% [markdown]
# ## 3. Figures

# %%
figs = figures.make_all(result, cfg, "outputs")
for name, p in figs.items():
    print(name, "->", p)

# %% [markdown]
# ## 4. The scientific payload: outcome vs prescription
#
# The answer is not one number -- it is a mapped dependence.

# %%
table = theory_compare.compare()
print(f"{'toggle':14s} {'option':16s} {'fate':50s} {'Mchirp':>8s}")
for row in table:
    mc = row.get("chirp_mass")
    mc = f"{mc:.2f}" if isinstance(mc, float) and np.isfinite(mc) else "-"
    print(f"{row['toggle']:14s} {row['option']:16s} {str(row.get('fate',''))[:50]:50s} {mc:>8s}")

# %% [markdown]
# ## 5. Which uncertain input dominates?

# %%
sens = sensitivity.one_at_a_time()
print(f"outcome metric: {sens['rows'][0]['influence_metric']}")
for i, r in enumerate(sens["rows"][:6], 1):
    print(f"{i}. {r['param']:16s} spread={r['influence']:.4g}   ({r['ref'][:60]})")
print("\nDominant uncertain input:", sens["rows"][0]["param"])

# %% [markdown]
# ## 6. Vary the wind prescription by hand
#
# Vink 2001 vs Bjorklund 2021 is one of the largest current systematics in
# massive-star evolution.

# %%
for w in ["vink2001", "bjorklund2021", "nugis2000", "dejager1988"]:
    r = evolve(load_config(overrides={"theories.wind": w}))
    print(f"{w:14s} -> {r.fate}")

# %% [markdown]
# With the default (Vink) the system forms a BH+BH that merges via GWs within a
# Hubble time; the weaker Bjorklund rates change the mass-loss history and hence
# the interaction outcome; the strong WR (Nugis & Lamers) rates strip so much mass
# that the remnants become neutron stars. **The prescription choice changes the
# qualitative answer** -- this is the point of the framework.
