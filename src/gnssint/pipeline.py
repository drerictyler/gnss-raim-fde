"""End-to-end processing of an observation file, with optional fault injection."""

from dataclasses import dataclass

import numpy as np

from .orbits import ecef_to_geodetic, enu_rotation
from .raim import RaimConfig, raim_fde
from .solver import prepare_epoch, wls

# Applied to the fitted a and b, so it multiplies sigma (the variance by 1.25^2).
SIGMA_INFLATION = 1.25


@dataclass
class EpochRecord:
    t: float
    n_sat: int
    enu_error: np.ndarray       # error of the *reported* solution (after any exclusion)
    clock: float                # c * dt_rx [m] of the reported solution
    clock_var: float
    test_stat: float
    threshold: float
    alarm: bool
    excluded: str
    usable: bool                # False when an alarm could not be resolved by exclusion
    hpl: float
    vpl: float


def process(obs, nav, ref_xyz, sigma_ab, raim_cfg=RaimConfig(), mode="L1",
            mask_deg=10.0, inject=None, exclude=True):
    """Run WLS + RAIM/FDE on every epoch.

    ``inject(meas, R_enu) -> pr`` may replace the corrected pseudoranges to
    simulate a fault or a spoofing attack. With ``exclude=False`` detection still
    runs but the all-satellite solution is reported (useful as a baseline).
    """
    lat, lon, _ = ecef_to_geodetic(ref_xyz)
    R = enu_rotation(lat, lon)
    x0 = np.r_[ref_xyz, 0.0]
    out = []
    for ep in obs.epochs:
        meas = prepare_epoch(ep, nav, ref_xyz, mode=mode, mask_deg=mask_deg, sigma_ab=sigma_ab)
        if len(meas.sats) < 4:
            continue
        if inject is not None:
            meas.pr = inject(meas, R)
        fde = raim_fde(meas, R, x0, raim_cfg)
        det = fde.detection
        sol = fde.solution if exclude else wls(meas.sat_pos, meas.pr, meas.sigma, x0=x0)
        out.append(EpochRecord(
            t=ep.t, n_sat=len(meas.sats),
            enu_error=R @ (sol.x[:3] - ref_xyz),
            clock=sol.x[3], clock_var=sol.cov[3, 3],
            test_stat=det.test_stat, threshold=det.threshold, alarm=det.alarm,
            excluded=fde.excluded if exclude else None,
            usable=fde.usable if exclude else True, hpl=det.hpl, vpl=det.vpl,
        ))
        x0 = np.r_[ref_xyz, sol.x[3]]
    return out


def calibrate_sigma(obs, nav, ref_xyz, sigma_ab=(0.6, 0.9), mode="L1", mask_deg=10.0,
                    inflation=SIGMA_INFLATION, iterations=2):
    """Fit sigma^2(el) = a^2 + b^2 / sin^2(el) to fault-free residuals.

    Uses r_i / sqrt(1 - P_ii), whose variance is sigma_i^2 when the model holds,
    then multiplies a and b (and therefore sigma) by ``inflation`` so the model is
    deliberately conservative rather than a bare fit. The returned coefficients
    already include the inflation.
    """
    a, b = sigma_ab
    for _ in range(iterations):
        rows, targets = [], []
        for ep in obs.epochs:
            meas = prepare_epoch(ep, nav, ref_xyz, mode=mode, mask_deg=mask_deg, sigma_ab=(a, b))
            if len(meas.sats) < 6:
                continue
            sol = wls(meas.sat_pos, meas.pr, meas.sigma, x0=np.r_[ref_xyz, 0.0])
            w = 1.0 / meas.sigma ** 2
            S = np.linalg.solve(sol.G.T @ (w[:, None] * sol.G), sol.G.T * w)
            p_diag = np.einsum("ij,ji->i", sol.G, S)
            e_hat = sol.residuals / np.sqrt(np.clip(1.0 - p_diag, 1e-6, None))
            for el, e in zip(meas.el, e_hat):
                rows.append([1.0, 1.0 / np.sin(el) ** 2])
                targets.append(e ** 2)
        coef, *_ = np.linalg.lstsq(np.array(rows), np.array(targets), rcond=None)
        a2, b2 = np.clip(coef, 1e-4, None)
        a, b = np.sqrt(a2), np.sqrt(b2)
    return inflation * a, inflation * b
