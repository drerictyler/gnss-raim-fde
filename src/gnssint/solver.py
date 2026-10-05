"""Measurement preparation and snapshot weighted least-squares (WLS) positioning."""

from dataclasses import dataclass

import numpy as np

from .constants import C, seconds_of_week
from .corrections import iono_free, klobuchar_delay, sigma_model, tropo_delay
from .orbits import (azimuth_elevation, ecef_to_geodetic, enu_rotation,
                     sagnac_rotate, sat_position_clock, travel_time)


@dataclass
class EpochMeasurements:
    """Corrected pseudoranges ready for the WLS model  y_i = |s_i - x| + c*dt_rx + e_i."""
    t: float
    sats: list
    sat_pos: np.ndarray    # (N, 3) ECEF at receive time (Sagnac-corrected)
    pr: np.ndarray         # (N,) corrected pseudoranges [m]
    sigma: np.ndarray      # (N,) 1-sigma errors [m]
    el: np.ndarray         # (N,) elevations [rad]


def prepare_epoch(epoch, nav, rx_approx, mode="L1", mask_deg=10.0, sigma_ab=(0.6, 0.9)):
    """Compute satellite positions and correct raw pseudoranges for one epoch.

    mode="L1":  C1 code with broadcast clock (minus TGD) and Klobuchar ionosphere,
                the classic single-frequency aviation receiver.
    mode="IF":  iono-free C1/P2 combination (no TGD, no iono model).
    """
    lat, lon, h = ecef_to_geodetic(rx_approx)
    sats, sat_pos, pr_out, els = [], [], [], []
    for sat, obs in sorted(epoch.obs.items()):
        c1 = obs.get("C1") or obs.get("P1")
        if not c1:
            continue
        if mode == "IF":
            p2 = obs.get("P2")
            if not p2:
                continue
            raw = iono_free(c1, p2)
        else:
            raw = c1
        eph = nav.select(sat, epoch.t)
        if eph is None:
            continue

        t_tx = epoch.t - raw / C                       # transmit time (receiver-clock frame)
        _, dts = sat_position_clock(eph, t_tx)
        t_tx -= dts
        pos, dts = sat_position_clock(eph, t_tx)
        if mode == "L1":
            dts -= eph.tgd
        pos = sagnac_rotate(pos, travel_time(rx_approx, pos))

        az, el = azimuth_elevation(rx_approx, pos)
        if np.degrees(el) < mask_deg:
            continue
        corrected = raw + C * dts - tropo_delay(lat, h, el)
        if mode == "L1":
            corrected -= klobuchar_delay(nav.ion_alpha, nav.ion_beta, lat, lon, az, el,
                                         seconds_of_week(epoch.t))
        sats.append(sat)
        sat_pos.append(pos)
        pr_out.append(corrected)
        els.append(el)

    els = np.array(els)
    return EpochMeasurements(t=epoch.t, sats=sats, sat_pos=np.array(sat_pos).reshape(-1, 3),
                             pr=np.array(pr_out), sigma=sigma_model(els, *sigma_ab), el=els)


@dataclass
class WlsSolution:
    x: np.ndarray          # [X, Y, Z, c*dt_rx] in metres
    G: np.ndarray          # (N, 4) geometry matrix at the solution
    residuals: np.ndarray  # (N,) post-fit residuals y - h(x)
    cov: np.ndarray        # (4, 4) (G^T W G)^-1
    converged: bool


def wls(sat_pos, pr, sigma, x0=None, max_iter=10, tol=1e-4):
    """Gauss-Newton weighted least squares for position and receiver clock."""
    x = np.zeros(4) if x0 is None else np.array(x0, dtype=float)
    w = 1.0 / sigma ** 2
    converged = False
    for _ in range(max_iter):
        d = sat_pos - x[:3]
        rho = np.linalg.norm(d, axis=1)
        G = np.hstack([-d / rho[:, None], np.ones((len(pr), 1))])
        y = pr - (rho + x[3])
        N = G.T @ (w[:, None] * G)
        dx = np.linalg.solve(N, G.T @ (w * y))
        x += dx
        if np.linalg.norm(dx) < tol:
            converged = True
            break
    d = sat_pos - x[:3]
    rho = np.linalg.norm(d, axis=1)
    G = np.hstack([-d / rho[:, None], np.ones((len(pr), 1))])
    residuals = pr - (rho + x[3])
    cov = np.linalg.inv(G.T @ (w[:, None] * G))
    return WlsSolution(x=x, G=G, residuals=residuals, cov=cov, converged=converged)


def enu_error(x_est, x_ref):
    lat, lon, _ = ecef_to_geodetic(x_ref)
    return enu_rotation(lat, lon) @ (x_est - x_ref)
