"""Rosseland-mean opacity prescriptions (selectable via theories.opacity).

  electron_scattering : kappa_es = 0.2 (1 + X)                     [cm^2/g]
       (Thomson scattering, KWW12 eq. 17.2)
  kramers             : kappa_es  +  bound-free + free-free Kramers law
       kappa_bf = 4.34e25 Z (1+X) rho T^-3.5    (KWW12 eq. 17.5, g_bf/t ~ 1)
       kappa_ff = 3.68e22 (1-Z)(1+X) rho T^-3.5 (KWW12 eq. 17.6)
       + H- term  kappa_H- = 2.5e-31 (Z/0.02) rho^0.5 T^9   (KWW12 eq. 17.8, cool)
  opal                : OPAL/OP table interpolation.  NOT bundled (licensing);
       falls back to `kramers` with a one-time warning.  See KNOWN_GAPS OPAC-OPAL.
"""
from __future__ import annotations

import warnings

import numpy as np

_OPAL_WARNED = {"done": False}


def kappa_electron_scattering(X):
    return 0.20 * (1.0 + X)


def kappa_kramers(rho, T, X, Y, Z):
    kes = 0.20 * (1.0 + X)
    kbf = 4.34e25 * Z * (1.0 + X) * rho * T ** -3.5
    kff = 3.68e22 * (1.0 - Z) * (1.0 + X) * rho * T ** -3.5
    kHm = 2.5e-31 * (Z / 0.02) * np.sqrt(rho) * T ** 9 if T < 1.0e4 else 0.0
    return kes + kbf + kff + kHm


def kappa(rho, T, X, Y, Z, prescription="kramers"):
    if prescription == "electron_scattering":
        return kappa_electron_scattering(X)
    if prescription == "kramers":
        return kappa_kramers(rho, T, X, Y, Z)
    if prescription == "opal":
        if not _OPAL_WARNED["done"]:
            warnings.warn("opacity='opal': OPAL tables not bundled (licensing); "
                          "using Kramers fallback (KNOWN_GAPS OPAC-OPAL).", RuntimeWarning)
            _OPAL_WARNED["done"] = True
        return kappa_kramers(rho, T, X, Y, Z)
    raise ValueError(f"unknown opacity prescription {prescription!r}")
