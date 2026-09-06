"""massive-binary-sim: a physically-traceable simulation framework for a bound
pair of massive stars (8-150 Msun each).

Design rules (see README):
  * every equation is traceable to a cited source;
  * competing prescriptions are runtime-selectable, never hardcoded;
  * every coefficient lives in config.yaml with a literature range;
  * nothing is fabricated -- unknowns are marked UNVERIFIED and listed in
    KNOWN_GAPS.md.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .io.config import Config, load_config  # noqa: E402

__all__ = ["Config", "load_config", "__version__"]
