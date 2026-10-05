"""Physical constants (IS-GPS-200, WGS-84) and GPS time helpers."""

from datetime import datetime

C = 299_792_458.0                 # speed of light [m/s]
MU = 3.986005e14                  # WGS-84 Earth gravitational constant [m^3/s^2] (IS-GPS-200 value)
OMEGA_E = 7.2921151467e-5         # WGS-84 Earth rotation rate [rad/s]
F_REL = -4.442807633e-10          # relativistic clock correction constant [s/m^0.5]
F1 = 1575.42e6                    # GPS L1 [Hz]
F2 = 1227.60e6                    # GPS L2 [Hz]
WGS84_A = 6_378_137.0             # semi-major axis [m]
WGS84_F = 1.0 / 298.257223563     # flattening
SECONDS_PER_WEEK = 604_800.0

_GPS_EPOCH = datetime(1980, 1, 6)


def gps_seconds(year, month, day, hour, minute, second):
    """Continuous GPS time in seconds since 1980-01-06 00:00:00 (no leap seconds).

    Working in continuous seconds avoids week-rollover bookkeeping: an ephemeris
    reference time becomes ``week * 604800 + toe`` and ``tk = t - toe_abs``.
    """
    whole = datetime(year, month, day, hour, minute) - _GPS_EPOCH
    return whole.days * 86_400.0 + whole.seconds + second


def seconds_of_week(t_abs):
    return t_abs % SECONDS_PER_WEEK
