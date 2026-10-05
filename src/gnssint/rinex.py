"""Minimal RINEX 2.x reader for GPS observation and navigation files.

Deliberately small and dependency-free so every step is inspectable. It handles
what integrity experiments need (GPS pseudoranges, broadcast ephemerides,
Klobuchar coefficients) and skips everything else. For RINEX 3 or
multi-constellation work, swap in a full parser such as ``georinex``.
"""

from dataclasses import dataclass, field

import numpy as np

from .constants import SECONDS_PER_WEEK, gps_seconds


def _year4(yy):
    yy = int(yy)
    return 1900 + yy if yy >= 80 else 2000 + yy


def _num(text):
    """Parse a Fortran-style float field ('D' exponents, blanks -> 0)."""
    text = text.strip().replace("D", "E").replace("d", "e")
    return float(text) if text else 0.0


# --------------------------------------------------------------------------- #
# Observation files
# --------------------------------------------------------------------------- #
@dataclass
class ObsEpoch:
    t: float                                   # continuous GPS seconds
    obs: dict = field(default_factory=dict)    # "G07" -> {"C1": ..., "P2": ...}


@dataclass
class ObsFile:
    obs_types: list
    approx_xyz: np.ndarray
    epochs: list


def read_obs(path):
    """Read a RINEX 2.x observation file (GPS satellites only)."""
    with open(path, "r", errors="replace") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]

    obs_types, approx = [], np.zeros(3)
    i = 0
    while i < len(lines):
        ln = lines[i]
        label = ln[60:].strip()
        if label == "# / TYPES OF OBSERV":
            if ln[:6].strip():
                n_types = int(ln[:6])
            obs_types += ln[6:60].split()
        elif label == "APPROX POSITION XYZ":
            approx = np.array([float(v) for v in ln[:42].split()])
        elif label == "END OF HEADER":
            i += 1
            break
        i += 1
    obs_types = obs_types[:n_types]
    lines_per_sat = (len(obs_types) + 4) // 5

    epochs = []
    while i < len(lines):
        ln = lines[i]
        if len(ln) < 32 or not ln[28].isdigit() or not ln[29:32].strip().isdigit():
            i += 1
            continue
        flag = int(ln[28])
        n_sat = int(ln[29:32])
        if flag > 1:                          # event record (e.g. inserted comments): skip it
            i += 1 + n_sat
            continue
        t = gps_seconds(_year4(ln[1:3]), int(ln[4:6]), int(ln[7:9]),
                        int(ln[10:12]), int(ln[13:15]), float(ln[15:26]))
        sat_text = ln[32:68]
        i += 1
        while len(sat_text) < 3 * n_sat:      # continuation lines for >12 satellites
            sat_text += lines[i][32:68]
            i += 1
        sats = []
        for k in range(n_sat):
            chunk = sat_text[3 * k:3 * k + 3]
            system = chunk[0] if chunk[0].strip() else "G"
            sats.append(f"{system}{int(chunk[1:3]):02d}")

        epoch = ObsEpoch(t=t)
        for sat in sats:
            record = "".join(lines[i + j].ljust(80) for j in range(lines_per_sat))
            i += lines_per_sat
            if not sat.startswith("G"):
                continue
            values = {}
            for k, name in enumerate(obs_types):
                field_text = record[16 * k:16 * k + 14]
                if field_text.strip():
                    values[name] = float(field_text)
            epoch.obs[sat] = values
        epochs.append(epoch)

    return ObsFile(obs_types=obs_types, approx_xyz=approx, epochs=epochs)


# --------------------------------------------------------------------------- #
# Navigation files
# --------------------------------------------------------------------------- #
@dataclass
class GpsEphemeris:
    prn: str
    toc: float          # clock reference time, continuous GPS seconds
    af0: float
    af1: float
    af2: float
    iode: float
    crs: float
    dn: float
    m0: float
    cuc: float
    e: float
    cus: float
    sqrt_a: float
    toe: float          # ephemeris reference time, continuous GPS seconds
    cic: float
    omega0: float
    cis: float
    i0: float
    crc: float
    omega: float
    omega_dot: float
    idot: float
    week: int
    ura: float
    health: float
    tgd: float


@dataclass
class NavFile:
    ion_alpha: np.ndarray
    ion_beta: np.ndarray
    ephemerides: dict   # "G07" -> [GpsEphemeris, ...]

    def select(self, sat, t, max_age=7200.0):
        """Healthy ephemeris with the closest toe, or None if none is valid."""
        best, best_dt = None, max_age
        for eph in self.ephemerides.get(sat, []):
            dt = abs(t - eph.toe)
            if eph.health == 0 and dt <= best_dt:
                best, best_dt = eph, dt
        return best


def read_nav(path):
    """Read a RINEX 2.x GPS navigation file."""
    with open(path, "r", errors="replace") as fh:
        lines = [ln.rstrip("\r\n") for ln in fh]

    alpha, beta = np.zeros(4), np.zeros(4)
    i = 0
    while i < len(lines):
        ln = lines[i]
        label = ln[60:].strip()
        if label == "ION ALPHA":
            alpha = np.array([_num(ln[2 + 12 * k:14 + 12 * k]) for k in range(4)])
        elif label == "ION BETA":
            beta = np.array([_num(ln[2 + 12 * k:14 + 12 * k]) for k in range(4)])
        elif label == "END OF HEADER":
            i += 1
            break
        i += 1

    ephs = {}
    while i + 7 < len(lines):
        first = lines[i].ljust(79)
        if not first[:2].strip():
            i += 1
            continue
        prn = f"G{int(first[:2]):02d}"
        toc_abs = gps_seconds(_year4(first[3:5]), int(first[6:8]), int(first[9:11]),
                              int(first[12:14]), int(first[15:17]), float(first[17:22]))
        af = [_num(first[22 + 19 * k:41 + 19 * k]) for k in range(3)]
        rows = []
        for j in range(1, 8):
            ln = lines[i + j].ljust(79)
            rows.append([_num(ln[3 + 19 * k:22 + 19 * k]) for k in range(4)])
        i += 8

        (iode, crs, dn, m0), (cuc, e, cus, sqrt_a), (toe_sow, cic, omega0, cis), \
            (i0, crc, omega, omega_dot), (idot, _, week, _), (ura, health, tgd, _), _ = rows
        week = int(week)
        eph = GpsEphemeris(
            prn=prn, toc=toc_abs, af0=af[0], af1=af[1], af2=af[2],
            iode=iode, crs=crs, dn=dn, m0=m0, cuc=cuc, e=e, cus=cus, sqrt_a=sqrt_a,
            toe=week * SECONDS_PER_WEEK + toe_sow, cic=cic, omega0=omega0, cis=cis,
            i0=i0, crc=crc, omega=omega, omega_dot=omega_dot, idot=idot,
            week=week, ura=ura, health=health, tgd=tgd,
        )
        ephs.setdefault(prn, []).append(eph)

    return NavFile(ion_alpha=alpha, ion_beta=beta, ephemerides=ephs)
