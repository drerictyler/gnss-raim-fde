import numpy as np
from scipy.stats import ncx2

from gnssint import faults
from gnssint.raim import RaimConfig, detection_threshold, pbias, raim_check, raim_fde
from gnssint.solver import wls


def test_pbias_inverts_the_missed_detection_probability():
    for dof in (1, 3, 6):
        thr = detection_threshold(dof, 1e-5)
        pb = pbias(dof, thr, 1e-3)
        assert abs(ncx2.cdf(thr, dof, pb ** 2) - 1e-3) < 1e-7


def test_time_push_is_invisible_to_raim(epoch, data):
    meas, sol = epoch
    shifted = wls(meas.sat_pos, meas.pr + 150.0, meas.sigma, x0=sol.x)
    base = raim_check(sol.G, meas.sigma, sol.residuals, data["R"])
    push = raim_check(shifted.G, meas.sigma, shifted.residuals, data["R"])
    assert abs(push.test_stat - base.test_stat) < 1e-6
    assert abs((shifted.x[3] - sol.x[3]) - 150.0) < 1e-6     # absorbed by the clock unknown


def test_position_pull_off_moves_the_fix_without_changing_the_test(epoch, data):
    meas, sol = epoch
    meas_t = meas.t
    pr = faults.position_pull_off(meas, data["ref"], meas_t - 100.0, (0.5, 0.0, 0.0), data["R"])
    pulled = wls(meas.sat_pos, pr, meas.sigma, x0=sol.x)
    shift_enu = data["R"] @ (pulled.x[:3] - sol.x[:3])
    np.testing.assert_allclose(shift_enu, [50.0, 0.0, 0.0], atol=0.05)
    base = raim_check(sol.G, meas.sigma, sol.residuals, data["R"])
    moved = raim_check(pulled.G, meas.sigma, pulled.residuals, data["R"])
    assert abs(moved.test_stat - base.test_stat) < 1e-3


def test_empirical_false_alert_rate(epoch):
    meas, sol = epoch
    G, sigma = sol.G, meas.sigma
    w = 1.0 / sigma ** 2
    M = np.eye(len(sigma)) - G @ np.linalg.solve(G.T @ (w[:, None] * G), G.T * w)
    rng = np.random.default_rng(1)
    r = (rng.standard_normal((100_000, len(sigma))) * sigma) @ M.T
    stats = (r ** 2 * w).sum(axis=1)
    p_fa = 0.01
    rate = (stats > detection_threshold(len(sigma) - 4, p_fa)).mean()
    assert abs(rate - p_fa) < 4 * np.sqrt(p_fa * (1 - p_fa) / 100_000)


def test_large_fault_is_excluded(epoch, data):
    meas, sol = epoch
    for k, sat in enumerate(meas.sats):
        faulty = meas.__class__(**{**meas.__dict__, "pr": meas.pr.copy()})
        faulty.pr[k] += 200.0
        result = raim_fde(faulty, data["R"], sol.x, RaimConfig())
        assert result.detection.alarm
        # With only 7 satellites exclusion can be ambiguous; it must never pick the wrong one.
        assert result.excluded in (sat, None)
        assert result.exclusion_ok or result.ambiguous
