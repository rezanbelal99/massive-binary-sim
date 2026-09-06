"""Convergence study: halve tolerance & timestep, measure observed order."""
from __future__ import annotations

import numpy as np

from ..dynamics.kepler import elements_to_state, state_to_elements
from ..numerics.integrators import integrate_orbit
from ..numerics.units import G, MSUN, RSUN


def orbit_convergence():
    """Refine rtol for a 2.5PN inspiral; measure convergence of the final state."""
    m1, m2 = 20 * MSUN, 15 * MSUN
    mu = G * (m1 + m2)
    a0, e0 = 3.0 * RSUN, 0.4
    r, v = elements_to_state(a0, e0, 0.0, 0.0, 0.0, 0.0, mu)
    y0 = np.concatenate([r, v])
    P = 2 * np.pi * np.sqrt(a0**3 / mu)

    tols = [1e-8, 1e-9, 1e-10, 1e-11, 1e-12]
    finals = []
    for rt in tols:
        t, y = integrate_orbit(y0, (0, 40 * P), m1, m2, order="2.5PN", method="dop853",
                               rtol=rt, atol=rt * 1e-2, dense_points=200)
        el = state_to_elements(y[-1][:3], y[-1][3:], mu)
        finals.append(el["a"])
    finals = np.array(finals)
    ref = finals[-1]
    errs = np.abs(finals[:-1] - ref) / ref
    # observed order between successive halvings (tol drops by 10x here)
    orders = np.log(errs[:-1] / errs[1:]) / np.log(10.0)
    return dict(tols=tols, final_a=finals.tolist(),
               rel_errors=errs.tolist(),
               observed_order_per_decade=orders.tolist(),
               richardson_ref_a=float(ref),
               detail=f"final a converges to {ref/RSUN:.6f} Rsun; "
                      f"rel error at rtol=1e-8 is {errs[0]:.2e}, at 1e-11 is {errs[-1]:.2e}; "
                      f"mean observed order ~ {np.nanmean(orders):.2f} per tol-decade")
