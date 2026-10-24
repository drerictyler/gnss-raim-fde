"""Snapshot residual-based RAIM with fault detection and exclusion (FDE).

Theory (weighted least squares, N satellites, 4 unknowns):

    S = (G^T W G)^-1 G^T W            estimate sensitivity, x_hat = S y
    r = (I - G S) y                   post-fit residuals
    T = r^T W r                       test statistic

Fault-free and with correct sigmas, T ~ chi^2(N - 4). An alarm is raised when
T > T_D = chi2.isf(P_fa, N - 4).

A bias b on satellite i makes T non-central with lambda_i = b^2 (1 - P_ii) / sigma_i^2,
where P = G S, and shifts the estimate by S[:, i] * b. Hence the error per unit of
sqrt(lambda) ("slope") is

    slope_i = |S_k,i| * sigma_i / sqrt(1 - P_ii)      (k = vertical, or the
                                                       horizontal norm of E and N)

Choosing lambda_det so that the missed-detection probability at the threshold
equals P_md gives p_bias = sqrt(lambda_det) and the slope-based protection levels

    HPL = max_i Hslope_i * p_bias,   VPL = max_i Vslope_i * p_bias

(Brown & Chin 1998; Walter & Enge 1995). Some formulations add a fault-free noise
term k * sigma; this implementation keeps the classic slope form and says so.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2, ncx2

from .solver import wls


@dataclass
class RaimConfig:
    p_fa: float = 1e-5        # false-alert probability per test
    p_md: float = 1e-3        # missed-detection probability for protection levels
    n_states: int = 4         # x, y, z, clock (single constellation)


def detection_threshold(dof, p_fa):
    return chi2.isf(p_fa, dof)


def pbias(dof, threshold, p_md):
    """sqrt of the non-centrality at which P(T < threshold) = p_md."""
    if ncx2.cdf(threshold, dof, 1e-9) <= p_md:
        return 0.0
    lam = brentq(lambda nc: ncx2.cdf(threshold, dof, nc) - p_md, 1e-9, 1e5)
    return np.sqrt(lam)


@dataclass
class RaimResult:
    available: bool
    alarm: bool
    test_stat: float
    threshold: float
    dof: int
    hpl: float
    vpl: float
    hslope: np.ndarray = field(repr=False)
    vslope: np.ndarray = field(repr=False)


def raim_check(G, sigma, residuals, R_enu, cfg=RaimConfig()):
    """Fault detection and protection levels for one epoch.

    G: (N, 4) ECEF geometry matrix, sigma: (N,) 1-sigma errors,
    residuals: (N,) post-fit residuals, R_enu: (3, 3) ECEF->ENU rotation.
    """
    n = len(residuals)
    dof = n - cfg.n_states
    if dof < 1:
        nan = float("nan")
        return RaimResult(False, False, nan, nan, dof, nan, nan, np.array([]), np.array([]))

    w = 1.0 / sigma ** 2
    S = np.linalg.solve(G.T @ (w[:, None] * G), G.T * w)      # (4, N)
    P_diag = np.einsum("ij,ji->i", G, S)                        # diag(G S)
    S_enu = R_enu @ S[:3]                                        # rows: E, N, U
    denom = np.sqrt(np.clip(1.0 - P_diag, 1e-12, None))
    hslope = np.hypot(S_enu[0], S_enu[1]) * sigma / denom
    vslope = np.abs(S_enu[2]) * sigma / denom

    t_stat = float(residuals @ (w * residuals))
    thr = detection_threshold(dof, cfg.p_fa)
    pb = pbias(dof, thr, cfg.p_md)
    return RaimResult(True, t_stat > thr, t_stat, thr, dof,
                      float(hslope.max() * pb), float(vslope.max() * pb), hslope, vslope)


@dataclass
class FdeResult:
    detection: RaimResult
    excluded: str = None          # satellite removed, if any
    exclusion_ok: bool = False
    ambiguous: bool = False       # more than one subset passed: cannot tell which is faulty
    solution: object = None       # WLS solution after exclusion (or original)

    @property
    def usable(self):
        """True when the output position may be used: no alarm, or a unique exclusion."""
        return (not self.detection.alarm) or self.exclusion_ok


def raim_fde(meas, R_enu, x0, cfg=RaimConfig()):
    """Detect; if alarmed, test every single-satellite exclusion.

    The exclusion is accepted only if exactly one subset passes. With little
    redundancy several subsets can pass, and picking the "best" of them can
    remove a healthy satellite while keeping the faulty one (this happens on
    the sample data with 6 satellites). An ambiguous exclusion is reported as an
    alarm with no usable position, which is the safe outcome.
    """
    sol = wls(meas.sat_pos, meas.pr, meas.sigma, x0=x0)
    det = raim_check(sol.G, meas.sigma, sol.residuals, R_enu, cfg)
    if not det.available or not det.alarm:
        return FdeResult(det, solution=sol)
    if det.dof < 2:                      # after exclusion there would be no redundancy
        return FdeResult(det, solution=sol)

    passing = []
    for k in range(len(meas.sats)):
        keep = np.arange(len(meas.sats)) != k
        sub = wls(meas.sat_pos[keep], meas.pr[keep], meas.sigma[keep], x0=sol.x)
        chk = raim_check(sub.G, meas.sigma[keep], sub.residuals, R_enu, cfg)
        if chk.available and not chk.alarm:
            passing.append((meas.sats[k], sub))
    if len(passing) == 1:
        return FdeResult(det, excluded=passing[0][0], exclusion_ok=True, solution=passing[0][1])
    return FdeResult(det, ambiguous=len(passing) > 1, solution=sol)
