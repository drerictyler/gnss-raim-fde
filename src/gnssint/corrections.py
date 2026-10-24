"""Atmospheric corrections, iono-free combination and the measurement error model."""

import numpy as np

from .constants import C, F1, F2


def klobuchar_delay(alpha, beta, lat, lon, az, el, t_sow):
    """GPS broadcast (Klobuchar) L1 ionospheric delay [m], IS-GPS-200 20.3.3.5.2.5.

    Angles in radians; ``t_sow`` is GPS seconds of week.
    """
    phi_u, lam_u, e = lat / np.pi, lon / np.pi, el / np.pi        # semicircles
    psi = 0.0137 / (e + 0.11) - 0.022
    phi_i = np.clip(phi_u + psi * np.cos(az), -0.416, 0.416)
    lam_i = lam_u + psi * np.sin(az) / np.cos(phi_i * np.pi)
    phi_m = phi_i + 0.064 * np.cos((lam_i - 1.617) * np.pi)
    t_local = (43_200.0 * lam_i + t_sow) % 86_400.0
    slant = 1.0 + 16.0 * (0.53 - e) ** 3
    amp = max(sum(alpha[n] * phi_m ** n for n in range(4)), 0.0)
    per = max(sum(beta[n] * phi_m ** n for n in range(4)), 72_000.0)
    x = 2.0 * np.pi * (t_local - 50_400.0) / per
    if abs(x) < 1.57:
        delay = slant * (5e-9 + amp * (1.0 - x ** 2 / 2.0 + x ** 4 / 24.0))
    else:
        delay = slant * 5e-9
    return C * delay


def tropo_delay(lat, height, el):
    """Saastamoinen zenith delay (standard atmosphere, 50 % humidity) with the
    RTCA DO-229 (MOPS) mapping function. Adequate at the decimetre level."""
    if height < -100.0 or height > 10_000.0 or el <= 0.0:
        return 0.0
    pressure = 1013.25 * (1.0 - 2.2557e-5 * height) ** 5.2568
    temp = 15.0 - 6.5e-3 * height + 273.16
    e_vap = 0.5 * 6.108 * np.exp((17.15 * temp - 4684.0) / (temp - 38.45))
    zhd = 0.0022768 * pressure / (1.0 - 0.00266 * np.cos(2 * lat) - 0.00028 * height / 1e3)
    zwd = 0.002277 * (1255.0 / temp + 0.05) * e_vap
    mapping = 1.001 / np.sqrt(0.002001 + np.sin(el) ** 2)
    return (zhd + zwd) * mapping


def iono_free(p1, p2):
    """First-order ionosphere-free pseudorange combination."""
    g1, g2 = F1 ** 2, F2 ** 2
    return (g1 * p1 - g2 * p2) / (g1 - g2)


def sigma_model(el, a=0.6, b=0.9):
    """1-sigma pseudorange error [m] as a function of elevation [rad]:
    sigma^2 = a^2 + (b / sin el)^2.

    The defaults are placeholders: calibrate a and b from fault-free data
    (see experiments/run_real_data.py) so the model *overbounds* the real
    errors. RAIM thresholds and protection levels are only as honest as this model.
    """
    return np.sqrt(a ** 2 + (b / np.sin(el)) ** 2)
