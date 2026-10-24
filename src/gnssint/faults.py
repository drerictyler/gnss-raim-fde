"""Fault and spoofing injection on corrected pseudoranges.

Three scenarios, chosen because each monitor sees a different subset of them:

  single_satellite_ramp   one satellite drifts (clock or signal fault)
                          -> RAIM detects and excludes it
  time_push               every pseudorange shifted by the same growing amount
                          -> invisible to RAIM, visible to the clock monitor
  position_pull_off       pseudoranges consistent with a false, drifting position
                          and a clean clock -> invisible to both; needs an
                          independent sensor (e.g. inertial cross-check)
"""

import numpy as np


def ramp(t, t_start, rate):
    """Linear ramp [m] that is zero before ``t_start`` and grows at ``rate`` m/s."""
    return max(0.0, t - t_start) * rate


def single_satellite_ramp(meas, sat, t_start, rate):
    pr = meas.pr.copy()
    if sat in meas.sats:
        pr[meas.sats.index(sat)] += ramp(meas.t, t_start, rate)
    return pr


def time_push(meas, t_start, rate):
    return meas.pr + ramp(meas.t, t_start, rate)


def position_pull_off(meas, rx_xyz, t_start, velocity_enu, R_enu):
    """Shift pseudoranges so the receiver solves for rx + offset(t).

    To first order a displacement d changes each range by -u_i . d, where u_i is
    the receiver-to-satellite unit vector.
    """
    offset_enu = np.array(velocity_enu, float) * max(0.0, meas.t - t_start)
    offset_ecef = R_enu.T @ offset_enu
    los = meas.sat_pos - rx_xyz
    u = los / np.linalg.norm(los, axis=1)[:, None]
    return meas.pr - u @ offset_ecef
