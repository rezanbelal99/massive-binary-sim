"""(Pulsational) pair-instability regime by helium-core mass.

Regime boundaries (He-core mass at core collapse), following
Heger & Woosley (2002), ApJ 567, 532 and Woosley (2017), ApJ 836, 244, as
compiled by Farmer et al. (2019), ApJ 887, 53:

    M_He <  ~32 Msun   : normal core collapse
    ~32 - ~64 Msun     : pulsational pair-instability (PPISN) -> mass-losing
                          pulses, remnant BH mass capped near ~35-45 Msun
    ~64 - ~135 Msun     : full pair-instability SN (PISN) -> NO remnant
    M_He >  ~135 Msun   : photodisintegration -> direct collapse to BH

These He-core boundaries translate to an initial-mass "black-hole mass gap"
roughly 45 - 120 Msun (Msun) at solar metallicity; the exact edges are uncertain
at the +-5 Msun level and depend on the 12C(a,g)16O rate (Farmer+ 2019, 2020).
"""
from __future__ import annotations

from dataclasses import dataclass

PPISN_LOW = 32.0
PISN_LOW = 64.0
PISN_HIGH = 135.0

# lower / upper edge of the resulting BH mass gap (Msun), Farmer+ 2019 fiducial
BH_GAP_LOW = 45.0
BH_GAP_HIGH = 120.0


@dataclass
class PairInstabilityVerdict:
    regime: str
    remnant_bh_mass: float | None   # Msun; None if PISN leaves nothing
    note: str


def classify(m_he_core_msun, c12_ag_rate=1.0):
    """Return a PairInstabilityVerdict.

    The C12(a,g) multiplier shifts the boundaries: a higher rate -> more O, less C,
    weaker pulses -> boundaries move up.  We apply the linear sensitivity from
    Farmer+ 2020 (ApJ 902, L36): d(M_gap_low)/d(ln rate) ~ +14 Msun.
    """
    shift = 14.0 * (c12_ag_rate - 1.0)
    lo = PPISN_LOW + 0.4 * shift
    pisn_lo = PISN_LOW + shift
    pisn_hi = PISN_HIGH + shift

    m = m_he_core_msun
    if m < lo:
        return PairInstabilityVerdict("core_collapse", None,
                                      f"M_He={m:.1f} < {lo:.1f}: no pair instability")
    if m < pisn_lo:
        # PPISN: remnant BH mass capped; Woosley 2017 gives ~ 0.9 * (pulse-stripped He core)
        cap = min(m * 0.90, BH_GAP_LOW + 0.4 * shift)
        return PairInstabilityVerdict("PPISN", cap,
                                      f"M_He={m:.1f}: pulsational PI, BH mass capped ~{cap:.1f} Msun")
    if m < pisn_hi:
        return PairInstabilityVerdict("PISN", None,
                                      f"M_He={m:.1f}: full pair-instability SN, NO remnant")
    return PairInstabilityVerdict("photodisintegration", m * 0.95,
                                  f"M_He={m:.1f}: photodisintegration -> direct collapse")
