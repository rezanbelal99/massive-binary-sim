# massive-binary-sim

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22555588.svg)](https://doi.org/10.5281/zenodo.22555588)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

> **Status: educational / research-scoping tool.** The gravitational-dynamics and
> GW core is validated to research grade (Hulse–Taylor orbital decay reproduced to
> 0.02%, machine-precision Kepler closure). The single-star evolution uses reduced
> fitting relations — not a stellar-structure integration — and the binary
> mass-transfer / common-envelope treatment is heuristic, so **binary-evolution
> outcomes are qualitative**, not publication-grade population synthesis. Every
> shortfall is in [`KNOWN_GAPS.md`](KNOWN_GAPS.md). Nothing is fabricated.

A physically-traceable simulation framework for a gravitationally bound pair of
massive stars (8–150 M⊙ each). Every equation is traceable to a cited source;
competing physical prescriptions are selectable at runtime; every coefficient
lives in `config.yaml` with a literature range; unknowns are flagged `UNVERIFIED`
and listed in [`KNOWN_GAPS.md`](KNOWN_GAPS.md) rather than invented.

**Scope of this build** (see §0 of the design spec): `binary-evolution` regime,
`semi-analytic + 1D-structure` fidelity, research-grade, Python package + CLI +
notebook. The `hydrodynamic-collision` and `relativistic-inspiral` regimes are
partially covered (Peters inspiral, PN dynamics) or out of scope (SPH) — see
KNOWN_GAPS item HYDRO.

## Install

```bash
cd massive-binary-sim
pip install -e ".[dev,sensitivity]"
```

Requires Python ≥ 3.10, numpy, scipy, matplotlib, astropy, pyyaml, h5py
(+ SALib for global sensitivity, pytest for the validation suite).

## Quickstart

```bash
# run the default configuration -> outputs/RESULTS.md + RESULTS.pdf + figures + run.h5
mbsim run

# with the full analysis payload (slower)
mbsim run --sensitivity --theory-compare --convergence

# override any config key
mbsim run --set system.m1=60 --set system.m2=45 --set theories.wind=bjorklund2021 \
          --set free_parameters.alpha_ce=3.0

# the Section-5 validation suite (also runs in pytest)
mbsim validate            # fast subset
mbsim validate --full     # 10^4-orbit Kepler closure etc.
pytest -q                 # same, as tests

# competing-prescription sweep / sensitivity ranking
mbsim compare
mbsim sensitivity          # one-at-a-time
mbsim sensitivity --morris # global (Morris elementary effects)
```

## What it computes

| layer | module | key references |
|---|---|---|
| Gravitational dynamics | `dynamics/` | Newtonian ↔ 1PN ↔ 2PN ↔ 2.5PN relative acceleration (Mora & Will 2004); orbit-averaged GW decay da/dt, de/dt, T_merge (Peters 1964); restricted-PN chirp + LISA/aLIGO/ET curves; Kozai–Lidov (Naoz 2016); optional MOND toggle |
| Stellar structure | `structure/` | Lane–Emden polytropes (analytic n=0,1,5 checks); ideal-gas+radiation EOS with β, ∇_ad, Γ₁ (KWW12); electron-scattering / Kramers / OPAL-fallback opacity; pp / CNO / 3α / ¹²C(α,γ) rates (KWW12/CF88/NACRE) with the ¹²C(α,γ) free multiplier; Eddington standard model (1D); reduced evolution tracks calibrated to Ekström+ 2012 / Brott+ 2011 / Hurley+ 2000 |
| Winds | `winds/` | Vink+ 2001 (with the bi-stability jump, Z^0.85), Björklund+ 2021, Nugis & Lamers 2000, de Jager+ 1988; consistent removal of orbital AM (Jeans mode) |
| Binary interaction | `binary/` | Eggleton (1983) Roche radius + full Roche potential & L1–L5; Case A/B/C classification; Ritter (1988) / Kolb & Ritter (1990) transfer rate; q_crit and ζ_ad–ζ_L stability; non-conservative orbital response (isotropic re-emission / L2 / circumbinary ring); Hut (1981) equilibrium tides; α_CE–λ and γ common envelope |
| Remnants | `remnants/` | Fryer+ 2012 rapid/delayed/direct-collapse; (P)PISN by He-core mass (Farmer+ 2019) with the ¹²C(α,γ) boundary shift; Blaauw + Hobbs (σ=265 km/s) fallback-scaled natal kicks; post-SN orbit (Kalogera 1996) |
| Numerics | `numerics/` | adaptive DOP853 / symplectic leapfrog / stiff Radau-BDF; timestep = min(orbital, nuclear, thermal, MT, wind, GW) with the binding criterion logged; energy/AM conservation monitors; one astropy-based units layer |
| Analysis | `analysis/` | derived quantities with formulae; one-at-a-time + Morris sensitivity; competing-prescription comparison table; convergence study; `RESULTS.md`/`.pdf` generator; the validation suite |

## Configuration

`config.yaml` is the single source of physical choices. Three blocks:

* `system:` — masses, separation/period, eccentricity, metallicity, spins.
* `theories:` — the prescription switchboard (13 toggles).
* `free_parameters:` — every tunable coefficient as `{value, range, ref}`. A
  parameter whose provenance is not citable carries `ref: UNVERIFIED` and is
  reported by `cfg.unverified_params()`.
* `numerics:` — integrators, tolerances, timestep bounds, RNG seed.

## Validation (Section 5)

`pytest -q` runs, and `RESULTS.md` §8 reports:

1. **Kepler closure** — 10⁴ dissipationless orbits, energy/AM/a drift < 10⁻⁶.
2. **Analytic polytrope** — Lane–Emden n=0,1,5 vs closed forms, machine precision.
3. **Solar structure** — Eddington standard model P_c, T_c, ρ_c within a factor 2
   of a standard solar model; T_eff within 15% of 5772 K. *(4.57-Gyr evolved-Sun
   calibration is a KNOWN GAP — STRUCT-1D.)*
4. **1PN periastron advance** — vs 6πGM/(c²a(1−e²)), < 2%.
5. **Hulse–Taylor** — Peters dP_b/dt vs the PSR B1913+16 GR prediction, < 1%
   (achieves ~0.02%).
6. **2PN conservation** — no secular energy drift.
7. **Eddington limit** — Γ > 1 is flagged.
8. **Convergence** — refine tolerance, observed order reported (`--convergence`).

External cross-checks against MESA/BPASS/COMPAS are **not** run here (those codes
are not installed) — KNOWN_GAPS item XCHECK-EXTERNAL.

## Presets

`examples/` ships four scenarios (`mbsim run --config examples/<name>.yaml`):

| file | scenario | outcome with default physics |
|---|---|---|
| `wide_noninteracting.yaml` | 40 + 30 M⊙, 9000 R⊙ — winds only, never fills a Roche lobe | wide non-merging BH+BH |
| `caseB_stripping.yaml` | 26 + 16 M⊙, 900 R⊙, low natal kick — Case A/B transfer + reverse transfer | wide BH+NS, survives both SNe |
| `common_envelope.yaml` | 30 + 12 M⊙, 45 R⊙, q ≈ 2.5 — unstable transfer → CE | common-envelope merger |
| `post_sn_bh_binary.yaml` | 42 + 33 M⊙, 90 R⊙, Z = 0.003 — CE hardening between the SNe | **BH+BH merging via GW in ~2.2 Gyr** |

*(Outcomes are model-dependent; `mbsim compare --config examples/<f>.yaml` shows how each toggle changes them.)*

## Output

`outputs/run.h5` — HDF5 with full metadata (config, code version, git hash, seed,
wall-clock, every validation result), the complete time-series history, the event
list, and the final state.

## Caveats

This is a **reduced** model: 1-D, orbit-averaged, fitting-formula single-star
evolution, no magnetic fields beyond the AM-transport toggle, a qualitative GW
waveform. The scientific output is a **mapped dependence** of the outcome on
prescription choice and free parameters (RESULTS.md §6–7), not a single number.
Read [`KNOWN_GAPS.md`](KNOWN_GAPS.md) before using any result quantitatively.

## Citation

If you use this software, cite the archived release:

> Belal, R. (2026). *massive-binary-sim: a physically-traceable simulation
> framework for a bound pair of massive stars* (v0.1.0). Zenodo.
> https://doi.org/10.5281/zenodo.22555588

Concept DOI (always latest): [`10.5281/zenodo.22555588`](https://doi.org/10.5281/zenodo.22555588)
· v0.1.0 DOI: [`10.5281/zenodo.22555589`](https://doi.org/10.5281/zenodo.22555589)
· machine-readable metadata in [`CITATION.cff`](CITATION.cff).
