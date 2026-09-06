"""Lane-Emden polytropes.

Solve   (1/xi^2) d/dxi ( xi^2 dtheta/dxi ) = -theta^n ,   theta(0)=1, theta'(0)=0.

Analytic solutions used for validation (Chandrasekhar 1939,
"An Introduction to the Study of Stellar Structure", Ch. IV):
    n = 0 :  theta = 1 - xi^2/6 ,             xi_1 = sqrt(6)
    n = 1 :  theta = sin(xi)/xi ,             xi_1 = pi
    n = 5 :  theta = (1 + xi^2/3)^(-1/2) ,    xi_1 = infinity

Physical use in this project:
    n = 3   -> radiative massive-star envelope (Eddington standard model)
    n = 1.5 -> fully convective star / convective core
The polytrope gives the run of density used for tidal and Roche-geometry
calculations at the `polytrope` structure-fidelity level.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from ..numerics.units import G


def _lane_emden_rhs(xi, y, n):
    theta, dtheta = y
    theta_pos = max(theta, 0.0)
    return [dtheta, -(theta_pos ** n) - 2.0 / xi * dtheta]


def solve_lane_emden(n: float, xi_max: float = 40.0, rtol: float = 1e-10):
    """Integrate the Lane-Emden equation for index `n`.

    Returns dict with xi grid, theta, dtheta, and the first zero xi_1 (np.inf if
    none in [0, xi_max]) plus the derived structure constants:
        theta'(xi_1)              -- needed for mass
        rho_c/<rho> = -xi_1 / (3 theta'(xi_1))
    """
    # series start to avoid the xi=0 singularity: theta ~ 1 - xi^2/6 + n xi^4/120
    xi0 = 1e-4
    th0 = 1.0 - xi0**2 / 6.0 + n * xi0**4 / 120.0
    dth0 = -xi0 / 3.0 + n * xi0**3 / 30.0

    def event_zero(xi, y, n):
        return y[0]
    event_zero.terminal = True
    event_zero.direction = -1

    sol = solve_ivp(_lane_emden_rhs, (xi0, xi_max), [th0, dth0], args=(n,),
                    method="DOP853", rtol=rtol, atol=1e-12, dense_output=True,
                    events=event_zero, max_step=xi_max / 500.0)
    xi = np.linspace(xi0, sol.t[-1], 4000)
    Y = sol.sol(xi)
    theta, dtheta = Y[0], Y[1]

    if sol.t_events[0].size:
        xi1 = float(sol.t_events[0][0])
        dtheta_1 = float(sol.sol(xi1)[1])
        rhoc_over_mean = -xi1 / (3.0 * dtheta_1)
    else:
        xi1 = np.inf
        dtheta_1 = np.nan
        rhoc_over_mean = np.inf

    return dict(n=n, xi=xi, theta=np.clip(theta, 0, None), dtheta=dtheta,
                xi1=xi1, dtheta1=dtheta_1, rhoc_over_mean=rhoc_over_mean)


def analytic_theta(n: float, xi: np.ndarray) -> np.ndarray:
    if n == 0:
        return 1.0 - xi**2 / 6.0
    if n == 1:
        return np.sinc(xi / np.pi)  # sin(xi)/xi
    if n == 5:
        return (1.0 + xi**2 / 3.0) ** -0.5
    raise ValueError("analytic Lane-Emden only for n in {0,1,5}")


ANALYTIC_XI1 = {0: np.sqrt(6.0), 1: np.pi, 5: np.inf}


class Polytrope:
    """A polytropic stellar model of given mass, radius and index."""

    def __init__(self, mass_g: float, radius_cm: float, n: float):
        self.M = mass_g
        self.R = radius_cm
        self.n = n
        self._le = solve_lane_emden(n)
        xi1 = self._le["xi1"]
        if not np.isfinite(xi1):
            raise ValueError(f"polytrope n={n} has infinite radius; choose n<5")
        self.xi1 = xi1
        self.dtheta1 = self._le["dtheta1"]
        # scale factors:  <rho> = M / (4/3 pi R^3);  rho_c/<rho> = -xi_1/(3 theta'_1)
        mean_rho = self.M / (4.0 / 3.0 * np.pi * self.R**3)
        self.rho_c = mean_rho * (-xi1 / (3.0 * self.dtheta1))
        # central pressure via  P_c = G M^2 / R^4 / (4 pi (n+1) theta'(xi_1)^2)
        self.P_c = G * self.M**2 / self.R**4 / (4.0 * np.pi * (n + 1.0) * self.dtheta1**2)

    def rho_of_r(self, r_cm):
        r_cm = np.atleast_1d(r_cm)
        xi = self.xi1 * r_cm / self.R
        th = np.interp(xi, self._le["xi"], self._le["theta"], right=0.0)
        return self.rho_c * np.clip(th, 0, None) ** self.n

    def mass_within(self, r_cm):
        r_cm = np.atleast_1d(r_cm)
        xi = np.clip(self.xi1 * r_cm / self.R, self._le["xi"][0], self.xi1)
        dth = np.interp(xi, self._le["xi"], self._le["dtheta"])
        return -4.0 * np.pi * self.R**3 / self.xi1**3 * self.rho_c * xi**2 * dth
