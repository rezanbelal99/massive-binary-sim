# KNOWN GAPS

Every item here is a place where the code does **less** than the specification, or
where a value could not be traced to a citable source. Nothing is silently
substituted: each gap is flagged in the code (search for the tag) and listed here
with what would close it. A flagged gap is acceptable; a fabricated citation is not.

Legend: **impact** = how much it affects the default-run conclusion.

---

## Dynamics

### PN-35 — 3.5PN direct acceleration not implemented
* **Where:** `dynamics/postnewtonian.py`, `order="3.5PN"`.
* **What:** Selecting `gravity: "3.5PN"` runs 2.5PN dynamics and emits a
  `RuntimeWarning`. The 1PN, 2PN and 2.5PN relative-acceleration terms
  (Mora & Will 2004, PRD 69, 104021, Eqs. 2.2–2.4) *are* implemented.
* **Close it:** transcribe and test the 3.5PN terms from Iyer & Will (1995,
  PRD 52, 6882) / Pati & Will (2002, PRD 65, 104008), or from Blanchet (2014,
  LRR 17, 2) Eq. (203).
* **Impact:** negligible for the default run (widest separations are Newtonian to
  1 part in 1e6); matters only for the last ~100 orbits before a compact merger.

### PN-2PN — 2PN coefficients partially validated
* **Where:** `dynamics/postnewtonian.py`, 2PN block.
* **What:** transcribed from Mora & Will (2004). Covered by the
  `energy_conservation_2PN` test (no secular drift) but not by an independent
  periastron-advance-at-2PN check.
* **Close it:** add a 2PN periastron-advance analytic comparison (e.g. Damour &
  Schäfer 1988).
* **Impact:** low.

### INT-IAS15 — no true Everhart/IAS15 integrator
* **Where:** `numerics/integrators.py`.
* **What:** `orbit_integrator: "ias15"` maps to SciPy `DOP853` (8th-order adaptive)
  at `rtol <= 1e-12`. `"wisdom_holman"` maps to a 2nd-order symplectic KDK
  leapfrog on the relative problem (the Hamiltonian is **not** operator-split into
  H_Kepler + H_PN).
* **Close it:** port `rebound`'s IAS15 (Rein & Spiegel 2015, MNRAS 446, 1424) or a
  Wisdom–Holman mapping with a universal-variable Kepler drift.
* **Impact:** low — the Kepler-closure and convergence tests pass to 1e-11; DOP853
  is not symplectic so very-long (>1e6 orbit) conservative integrations would
  drift, but the driver never runs those.

---

## Stellar structure & evolution

### STRUCT-1D — no Henyey relaxation; 1D model is the Eddington polytrope
* **Where:** `structure/stellar_1d.py::solve_zams`.
* **What:** `structure_fidelity: "onedim"` and the solar-structure test use an
  **n = 3 Eddington standard model** (polytrope fitted to M and the ZAMS radius),
  which gives central P, T, ρ within ~15–30 % of a standard solar model. The
  four-ODE shooting code (`_structure_rhs`, `_shoot_out`, `_shoot_in`,
  fitting-point residuals) is present and runs, but does **not** converge through
  the surface boundary layer with outward/inward shooting and is not used.
* **Close it:** implement a Henyey (Newton–Raphson block-tridiagonal) relaxation
  on an adaptive mass mesh with a proper atmosphere boundary condition
  (Kippenhahn, Weigert & Weiss 2012, Ch. 11).
* **Impact:** **high** for anything needing real interior profiles (Kippenhahn
  diagram, exact core masses). The default run uses the reduced evolution model
  (EVO-*), not this.

### EVO-ZAMS — ZAMS L, R and MS lifetime are calibrated fits, not integrations
* **Where:** `structure/evolution.py` (`zams_luminosity`, `zams_radius`,
  `ms_lifetime_yr`, `he_core_fraction_tams`, `co_core_mass`).
* **What:** broken power laws calibrated "by eye" to the Ekström et al. (2012) and
  Brott et al. (2011) Z ≈ 0.014 grids and to Hurley, Pols & Tout (2000). They are
  **not** the published Hurley+ 2000 fitting formulae (which are ~2000 lines of
  coefficient tables) and carry ~10–30 % systematic error.
* **Close it:** implement the full Hurley, Pols & Tout (2000, MNRAS 315, 543)
  SSE formulae, or interpolate published MESA/Geneva tracks.
* **Impact:** **high** — this is the dominant model systematic; §6 of RESULTS.md
  quantifies the free-parameter part but not this structural part.

### EVO-PHASE — post-MS radius evolution is schematic
* **Where:** `structure/evolution.py::StarState._update_observables`.
* **What:** HG/CHeB radii use single power laws with a blue/red split at the
  Humphreys–Davidson limit; there is no real giant branch, no thermal-pulse AGB
  analogue, no blue-loop.
* **Close it:** same as EVO-ZAMS.
* **Impact:** medium–high for the RLOF case classification (Case B vs C).

### NUC-REACLIB — REACLIB rates approximated
* **Where:** `structure/nuclear.py`.
* **What:** `reaction_rates: "reaclib"` uses the same analytic KWW12/CF88 fits as
  `"nacre"` with a representative +5 % CNO / −3 % 3α offset. The true JINA
  REACLIB v2 tabulated rates are not bundled.
* **Close it:** ship the REACLIB snapshot and a table interpolator.
* **Impact:** low for the default run (burning only sets timescales here).

### NUC-NU — thermal neutrino losses are order-of-magnitude
* **Where:** `structure/nuclear.py::eps_neutrino_thermal`.
* **Close it:** implement Itoh et al. (1996, ApJS 102, 411) pair/photo/plasma fits.
* **Impact:** low (only matters for the last ~1000 yr before collapse, which the
  reduced model does not resolve).

### OPAC-OPAL — OPAL/OP tables not bundled
* **Where:** `structure/opacity.py`.
* **What:** `opacity: "opal"` falls back to Kramers + electron scattering with a
  one-time warning (licensing prevents redistribution).
* **Close it:** user supplies OPAL Type-1/2 tables; add a bilinear interpolator.
* **Impact:** medium for envelope structure, low for the reduced default run.

### EOS-ION / EOS-DEG — partial ionisation and degeneracy are approximate
* **Where:** `structure/eos.py`.
* **What:** Saha treats only H and He with g-factors ≈ 1 and metals as a fixed
  0.5 e⁻/nucleon; degeneracy uses a `sqrt(P_ideal² + P_deg²)` bridge rather than
  Fermi–Dirac integrals.
* **Close it:** use the Timmes & Swesty (2000, ApJS 126, 501) `helmholtz` EOS.
* **Impact:** low for the default massive-star regime (radiation + fully-ionised
  ideal gas dominate).

---

## Binary interaction

### BIN-MDOT — RLOF rate has an efficiency prefactor
* **Where:** `binary/mass_transfer.py::mdot_rlof`.
* **What:** a Ritter (1988) / Kolb & Ritter (1990)-style analytic nozzle with an
  order-unity efficiency prefactor (`1e-2 * M_donor / tau_thermal`); the absolute
  normalisation is not calibrated against detailed models.
* **Close it:** calibrate against MESA `roche_lobe` test suite or Kolb & Ritter
  (1990) Table 1.
* **Impact:** medium — sets how fast Case A/B transfer proceeds, hence the orbit
  at contact.

### BIN-ECC — eccentric / asynchronous RLOF uses instantaneous separation
* **Where:** `binary/roche.py::roche_lobe_fill_factor`.
* **Close it:** Sepinsky, Willems, Kalogera & Rasio (2007, ApJ 660, 1624) full
  non-synchronous, eccentric Roche geometry.
* **Impact:** low for the near-circular default; higher for eccentric presets.

### CE-LAMBDA — fitted λ is crude
* **Where:** `binary/common_envelope.py::lambda_fit`.
* **What:** `log R`-based interpolation with a ×3 boost for extended RSGs; not
  from a structure table.
* **Close it:** Xu & Li (2010, ApJ 716, 114) / Klencki et al. (2021, A&A 645,
  A54) binding-energy tables.
* **Impact:** **high** when the system enters CE (changes survive/merge verdict);
  §7 theory-comparison shows the α_CE and λ_CE spread.

### DRV-STATEMACHINE — the RLOF/CE state machine is heuristic
* **Where:** `binary/driver.py`.
* **What:** stable mass transfer is integrated with the Ritter/Kolb-Ritter rate and
  the Soberman+ 1997 orbital response, but an episode is *declared complete* when
  the donor mass reaches ~1.2x its (fitted) core mass, or after a step/time budget
  (`mt_step_count`, `0.15 t_max`, `total_rlof_steps`). CE is limited to at most two
  episodes with a cooldown. A donor that keeps re-overflowing as it expands is
  eventually force-stripped. These are pragmatic guards, not physics: the detailed
  MT/CE bookkeeping (accretor spin-up, non-conservative fraction evolution,
  repeated RLOF) is not modelled.
* **Consequence:** the *sequence* of outcomes across parameter space is physically
  sensible (the theory-comparison table shows wind prescription flipping the
  outcome between BH+BH, BH+NS and NS+NS, and CE formalism setting whether the DCO
  merges within a Hubble time), but individual remnant masses carry large
  systematic error and the exact separation at DCO formation is only good to a
  factor ~2.
* **Close it:** a proper detailed-binary-evolution treatment (MESA `binary`, or the
  full COMPAS/binary_c prescription set with self-consistent orbital-AM accounting).

### DRV-SPLIT — first-order operator splitting in the driver
* **Where:** `binary/driver.py`.
* **What:** structure → winds → RLOF/CE → tides → GW → third-body applied
  sequentially within a step; split error is controlled only by the timestep
  limiter, not by iterating to convergence.
* **Close it:** Strang splitting or a coupled implicit step.
* **Impact:** low–medium; bounded by the `f_*` timestep fractions in
  `numerics/timestep.py`.

---

## Gravitational waves

### GW-WAVEFORM — 0PN Newtonian chirp only
* **Where:** `dynamics/gw.py::inspiral_waveform`.
* **What:** frequency sweep and phase are leading-order; eccentricity is dropped
  in the waveform (not in the da/dt, de/dt evolution). Fine for the
  characteristic-strain figure, **not** for matched filtering.
* **Close it:** TaylorF2 (3.5PN phase) or IMRPhenom.
* **Impact:** cosmetic for the stated science goal.

### GW-DETECTORS — aLIGO/ET curves are analytic approximations
* **Where:** `dynamics/gw.py`. LISA is the Robson, Cornish & Liu (2019) analytic
  curve (good); aLIGO/ET are broken-power-law fits accurate to ~×2 and labelled
  APPROXIMATE in the figure.
* **Close it:** ship the official ASD text files.
* **Impact:** cosmetic.

---

## Analysis

### SENS-SOBOL — global sensitivity is Morris, not Sobol, by default
* **Where:** `analysis/sensitivity.py`.
* **What:** Morris elementary-effects screening runs by default (cheap); a full
  Sobol variance decomposition is not wired to the CLI because the run count
  (`(2k+2)·N`) exceeds the "under 60 s" budget for k ≈ 12 parameters.
* **Close it:** `SALib.sample.sobol` + a longer budget or a surrogate.
* **Impact:** Morris `mu*` ranking is usually consistent with Sobol total-order
  indices; the dominant parameter identification is robust.

### XCHECK-EXTERNAL — no live MESA / BPASS / COMPAS cross-check
* **Where:** validation suite.
* **What:** the Hulse–Taylor, Lane–Emden, Kepler and Eddington-model checks are
  against analytic / observational benchmarks. A single-star track and a binary
  track are **not** compared against a run of MESA/BPASS/COMPAS (those codes are
  not installed here).
* **Close it:** add `tests/external/` that runs MESA `1M_pre_ms_to_wd` and a
  COMPAS grid and diffs the tracks.
* **Impact:** **high** for research use — the reduced evolution model (EVO-*) is
  currently unvalidated against a detailed code.

### Z-SCALE — solar metallicity convention
* **Where:** `config.yaml` (`metallicity: 0.0142`) vs `units.Z_SUN = 0.0134`.
* **What:** the config default uses the protosolar Z ≈ 0.0142; AGSS09 photospheric
  Z = 0.0134. Different grids (Brott+ 2011 use a bespoke scale) are not
  reconciled.
* **Impact:** low (a few percent on wind rates via the Z^0.85 scaling).

---

## Hydrodynamics

### HYDRO — not implemented
* The `hydrodynamic-collision` / `all-three-coupled` regimes from §0 of the spec
  are **out of scope** for this build (`binary-evolution` regime only, per the
  chosen defaults). No SPH/grid solver, no resolution-convergence study.
* **Close it:** wrap an external SPH code (e.g. `gadget`, `phantom`) or implement
  a self-gravitating SPH module; map the reduced 1D profile to 3D initial
  conditions.
