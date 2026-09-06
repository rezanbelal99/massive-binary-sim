"""Time integrators.

Two families are exposed:

* **Direct orbit integration** of the relative separation vector under the chosen
  PN acceleration -- used for the Kepler-closure and periastron-advance tests and
  for short "zoom-in" orbit segments.

* **Secular / structural integration** -- stiff and non-stiff ODE systems handed
  to SciPy's battle-tested solvers.

Integrator map
--------------
  leapfrog        : 2nd-order symplectic kick-drift-kick (implemented here)
  wisdom_holman   : alias -> leapfrog KDK on the relative problem (H = H_Kep + H_PN
                    is NOT operator-split here; see KNOWN_GAPS INT-WH)
  ias15           : high-accuracy adaptive -> SciPy 'DOP853' (8th order).  A true
                    Everhart/Rein-Spiegel IAS15 is NOT implemented (KNOWN_GAPS
                    INT-IAS15); DOP853 at rtol<=1e-12 is the substitute.
  dop853          : SciPy 'DOP853'
  radau / bdf     : SciPy implicit stiff solvers (for nuclear network / thermal)

No fixed-step RK4 is used for any long-term evolution, per the spec.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..dynamics.postnewtonian import relative_acceleration


# ---------------------------------------------------------------------------------------
# Direct orbit integration
# ---------------------------------------------------------------------------------------
def _rhs_relative(t, y, m1, m2, order):
    r = y[:3]
    v = y[3:]
    a = relative_acceleration(r, v, m1, m2, order=order)
    return np.concatenate([v, a])


def leapfrog_kdk(y0, t_span, dt, m1, m2, order, n_out=2000):
    """2nd-order symplectic KDK leapfrog for the relative two-body problem.

    For pure Newtonian / conservative PN orders this conserves a shadow
    Hamiltonian: energy error oscillates but does not secularly grow.  The 2.5PN
    term is dissipative by construction, so with order='2.5PN' the energy *should*
    decrease monotonically -- that is physics, not integrator error.
    """
    t0, t1 = t_span
    n_steps = int(np.ceil((t1 - t0) / dt))
    dt = (t1 - t0) / n_steps
    y = np.array(y0, dtype=float)
    r, v = y[:3], y[3:]

    out_idx = np.unique(np.linspace(0, n_steps, min(n_out, n_steps + 1)).astype(int))
    ts, ys = [], []
    a = relative_acceleration(r, v, m1, m2, order=order)
    for i in range(n_steps + 1):
        if i in out_idx:
            ts.append(t0 + i * dt)
            ys.append(np.concatenate([r, v]))
        if i == n_steps:
            break
        v = v + 0.5 * dt * a
        r = r + dt * v
        a = relative_acceleration(r, v, m1, m2, order=order)
        v = v + 0.5 * dt * a
    return np.array(ts), np.array(ys)


def integrate_orbit(y0, t_span, m1, m2, order="2.5PN", method="ias15",
                    rtol=1e-12, atol=1e-14, dense_points=4000, max_step=np.inf):
    """Dispatch orbit integration. Returns (t, y[N,6])."""
    method = method.lower()
    if method in ("leapfrog", "wisdom_holman"):
        # choose dt ~ P/2000 from the initial state
        r0 = np.linalg.norm(y0[:3])
        v0 = np.linalg.norm(y0[3:])
        from ..numerics.units import G
        mu = G * (m1 + m2)
        energy = 0.5 * v0**2 - mu / r0
        a_sma = -mu / (2 * energy) if energy < 0 else r0
        P = 2 * np.pi * np.sqrt(abs(a_sma) ** 3 / mu)
        dt = P / 3000.0
        return leapfrog_kdk(y0, t_span, dt, m1, m2, order, n_out=dense_points)

    scipy_method = {"ias15": "DOP853", "dop853": "DOP853",
                    "radau": "Radau", "bdf": "BDF"}.get(method, "DOP853")
    t_eval = np.linspace(t_span[0], t_span[1], dense_points)
    sol = solve_ivp(_rhs_relative, t_span, y0, method=scipy_method, t_eval=t_eval,
                    args=(m1, m2, order), rtol=rtol, atol=atol, max_step=max_step,
                    dense_output=False)
    if not sol.success:
        raise RuntimeError(f"orbit integration failed: {sol.message}")
    return sol.t, sol.y.T


def integrate_secular(rhs, t_span, y0, method="DOP853", rtol=1e-10, atol=1e-12,
                      t_eval=None, events=None, args=()):
    """Thin wrapper around solve_ivp for secular / structural ODE systems."""
    sol = solve_ivp(rhs, t_span, y0, method=method, rtol=rtol, atol=atol,
                    t_eval=t_eval, events=events, args=args, dense_output=True)
    return sol
