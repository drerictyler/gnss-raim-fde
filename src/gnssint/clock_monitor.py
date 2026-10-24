"""Receiver-clock consistency monitor for spoofing detection.

Why it is needed: a spoofer that shifts *every* pseudorange by the same amount
(a "time push") is invisible to snapshot RAIM, because the common-mode error is
absorbed by the receiver-clock unknown and the residuals do not change. The
receiver's own oscillator, however, is predictable over short horizons. If the
clock estimate suddenly stops following its oscillator model, something other
than the oscillator is moving it.

Model: three-state clock [bias b (m), drift d (m/s), drift rate a (m/s^2)], each
driven by white noise with power spectral densities q_bias, q_drift, q_accel
(m^2/s, m^2/s^3, m^2/s^5). The first two are the classic two-state oscillator
model (q_bias = c^2 h0 / 2, q_drift = 2 pi^2 c^2 h_-2, Brown & Hwang). The third
state is there because a real receiver oscillator's frequency also wanders: on the
sample data the clock drift moves from 418.1 to 420.5 m/s (about 8 ppb) over the
hour, and a two-state model tuned on the first half raises false alarms in the second.

Measurements are the WLS clock estimates. Two tests:
  * per-epoch normalised innovation squared, NIS_k ~ chi^2(1)
  * a sliding-window sum over M epochs, ~ chi^2(M), more sensitive to slow ramps.
"""

from dataclasses import dataclass, replace

import numpy as np
from scipy.stats import chi2


@dataclass
class ClockMonitorConfig:
    q_bias: float = 1e-4      # m^2/s
    q_drift: float = 1e-5     # m^2/s^3
    q_accel: float = 1e-10    # m^2/s^5
    p_fa: float = 1e-4        # per-test false-alert probability
    window: int = 10          # epochs in the sliding-window test


def transition(dt):
    return np.array([[1.0, dt, dt * dt / 2.0],
                     [0.0, 1.0, dt],
                     [0.0, 0.0, 1.0]])


def process_noise(dt, cfg):
    qb, qd, qa = cfg.q_bias, cfg.q_drift, cfg.q_accel
    return np.array([
        [qb * dt + qd * dt ** 3 / 3 + qa * dt ** 5 / 20, qd * dt ** 2 / 2 + qa * dt ** 4 / 8, qa * dt ** 3 / 6],
        [qd * dt ** 2 / 2 + qa * dt ** 4 / 8, qd * dt + qa * dt ** 3 / 3, qa * dt ** 2 / 2],
        [qa * dt ** 3 / 6, qa * dt ** 2 / 2, qa * dt],
    ])


@dataclass
class ClockTrace:
    t: np.ndarray
    innovation: np.ndarray
    nis: np.ndarray
    window_stat: np.ndarray
    epoch_alarm: np.ndarray
    window_alarm: np.ndarray
    drift: np.ndarray

    @property
    def alarm(self):
        return self.epoch_alarm | self.window_alarm


def run_clock_monitor(t, bias_meas, meas_var, cfg=ClockMonitorConfig()):
    """Filter a clock-bias series [m] and return innovation statistics and alarms.

    ``meas_var`` is the clock-measurement variance [m^2], scalar or per epoch.
    """
    t = np.asarray(t, float)
    z = np.asarray(bias_meas, float)
    r = np.broadcast_to(np.asarray(meas_var, float), z.shape)
    n = len(z)

    dt0 = t[1] - t[0]                                   # initialise from 3 epochs
    x = np.array([z[2], (z[2] - z[1]) / dt0, (z[2] - 2 * z[1] + z[0]) / dt0 ** 2])
    P = np.diag([r[2], 2 * r[2] / dt0 ** 2, 6 * r[2] / dt0 ** 4])

    innov, nis, drift = (np.full(n, np.nan) for _ in range(3))
    for k in range(3, n):
        dt = t[k] - t[k - 1]
        F = transition(dt)
        x = F @ x
        P = F @ P @ F.T + process_noise(dt, cfg)
        s = P[0, 0] + r[k]                               # H = [1, 0, 0]
        nu = z[k] - x[0]
        innov[k], nis[k] = nu, nu ** 2 / s
        gain = P[:, 0] / s
        x = x + gain * nu
        P = P - np.outer(gain, P[0, :])
        drift[k] = x[1]

    window = np.full(n, np.nan)
    for k in range(3 + cfg.window - 1, n):
        window[k] = np.nansum(nis[k - cfg.window + 1:k + 1])
    return ClockTrace(t=t, innovation=innov, nis=nis, window_stat=window,
                      epoch_alarm=nis > chi2.isf(cfg.p_fa, 1),
                      window_alarm=window > chi2.isf(cfg.p_fa, cfg.window),
                      drift=drift)


def estimate_measurement_variance(bias_meas):
    """White-noise variance of a smooth clock-bias series from its second differences.

    For b_k = smooth(t_k) + v_k with white v, var(b_{k+1} - 2 b_k + b_{k-1}) = 6 var(v).
    The WLS covariance is not used on purpose: it comes from the *inflated*
    pseudorange sigmas that RAIM needs, which would make this monitor insensitive.
    """
    return float(np.var(np.diff(np.asarray(bias_meas, float), 2)) / 6.0)


def calibrate(t, bias_meas, meas_var, cfg=ClockMonitorConfig()):
    """Grid-search q_drift and q_accel so the fault-free NIS averages closest to 1.

    This is the clock-model analogue of calibrating sigma: tune on data you
    trust to be clean, then validate on data you did not tune on.
    """
    best, best_err = cfg, np.inf
    for qd in np.logspace(-8, -2, 13):
        for qa in np.logspace(-14, -6, 17):
            trial = replace(cfg, q_drift=qd, q_accel=qa)
            err = abs(np.log(np.nanmean(run_clock_monitor(t, bias_meas, meas_var, trial).nis)))
            if err < best_err:
                best, best_err = trial, err
    return best
