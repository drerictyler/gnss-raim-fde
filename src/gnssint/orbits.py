"""GPS broadcast orbit/clock evaluation (IS-GPS-200 Table 20-IV) and geodesy helpers."""

import numpy as np

from .constants import C, F_REL, MU, OMEGA_E, WGS84_A, WGS84_F


def sat_position_clock(eph, t):
    """ECEF position [m] (at transmission time, Earth-fixed frame of that instant)
    and satellite clock offset [s] for ephemeris ``eph`` at GPS time ``t``.

    The clock offset includes the relativistic term but not TGD; subtract TGD
    for single-frequency L1 C/A users, leave it for the P1/P2 iono-free combination.
    """
    a = eph.sqrt_a ** 2
    tk = t - eph.toe
    n = np.sqrt(MU / a ** 3) + eph.dn
    mk = eph.m0 + n * tk

    ek = mk
    for _ in range(30):                        # Newton iteration on Kepler's equation
        step = (ek - eph.e * np.sin(ek) - mk) / (1.0 - eph.e * np.cos(ek))
        ek -= step
        if abs(step) < 1e-14:
            break

    vk = np.arctan2(np.sqrt(1.0 - eph.e ** 2) * np.sin(ek), np.cos(ek) - eph.e)
    phi = vk + eph.omega
    s2, c2 = np.sin(2 * phi), np.cos(2 * phi)
    uk = phi + eph.cus * s2 + eph.cuc * c2
    rk = a * (1.0 - eph.e * np.cos(ek)) + eph.crs * s2 + eph.crc * c2
    ik = eph.i0 + eph.cis * s2 + eph.cic * c2 + eph.idot * tk

    x_orb, y_orb = rk * np.cos(uk), rk * np.sin(uk)
    toe_sow = eph.toe - eph.week * 604_800.0
    om = eph.omega0 + (eph.omega_dot - OMEGA_E) * tk - OMEGA_E * toe_sow
    pos = np.array([
        x_orb * np.cos(om) - y_orb * np.cos(ik) * np.sin(om),
        x_orb * np.sin(om) + y_orb * np.cos(ik) * np.cos(om),
        y_orb * np.sin(ik),
    ])

    dt = t - eph.toc
    clock = eph.af0 + eph.af1 * dt + eph.af2 * dt ** 2 + F_REL * eph.e * eph.sqrt_a * np.sin(ek)
    return pos, clock


def sagnac_rotate(sat_pos, travel_time):
    """Rotate a transmit-time ECEF position into the receive-time ECEF frame."""
    theta = OMEGA_E * travel_time
    c, s = np.cos(theta), np.sin(theta)
    x, y, z = sat_pos
    return np.array([c * x + s * y, -s * x + c * y, z])


def ecef_to_geodetic(xyz):
    """WGS-84 latitude, longitude [rad] and ellipsoidal height [m]."""
    x, y, z = xyz
    e2 = WGS84_F * (2 - WGS84_F)
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(10):
        n = WGS84_A / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1 - e2 * n / (n + h)))
    n = WGS84_A / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - n
    return lat, lon, h


def enu_rotation(lat, lon):
    """Rows are the East, North, Up unit vectors expressed in ECEF."""
    sl, cl = np.sin(lat), np.cos(lat)
    so, co = np.sin(lon), np.cos(lon)
    return np.array([
        [-so, co, 0.0],
        [-sl * co, -sl * so, cl],
        [cl * co, cl * so, sl],
    ])


def azimuth_elevation(rx_xyz, sat_xyz):
    lat, lon, _ = ecef_to_geodetic(rx_xyz)
    e, n, u = enu_rotation(lat, lon) @ (sat_xyz - rx_xyz)
    return np.arctan2(e, n) % (2 * np.pi), np.arctan2(u, np.hypot(e, n))


def travel_time(rx_xyz, sat_xyz):
    return np.linalg.norm(sat_xyz - rx_xyz) / C
