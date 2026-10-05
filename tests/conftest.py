import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssint.orbits import ecef_to_geodetic, enu_rotation  # noqa: E402
from gnssint.pipeline import calibrate_sigma  # noqa: E402
from gnssint.rinex import read_nav, read_obs  # noqa: E402
from gnssint.solver import prepare_epoch, wls  # noqa: E402

OBS = ROOT / "data" / "07590920.05o"
NAV = ROOT / "data" / "07590920.05n"


@pytest.fixture(scope="session")
def data():
    if not OBS.exists() or not NAV.exists():
        pytest.skip("sample data missing: run python scripts/get_sample_data.py")
    obs, nav = read_obs(OBS), read_nav(NAV)
    ref = obs.approx_xyz
    lat, lon, _ = ecef_to_geodetic(ref)
    sigma_ab = calibrate_sigma(obs, nav, ref)
    return dict(obs=obs, nav=nav, ref=ref, R=enu_rotation(lat, lon), sigma_ab=sigma_ab)


@pytest.fixture(scope="session")
def epoch(data):
    """Measurements and fault-free solution at t = 30 min."""
    ep = next(e for e in data["obs"].epochs if e.t - data["obs"].epochs[0].t >= 1800.0)
    meas = prepare_epoch(ep, data["nav"], data["ref"], sigma_ab=data["sigma_ab"])
    sol = wls(meas.sat_pos, meas.pr, meas.sigma, x0=np.r_[data["ref"], 0.0])
    return meas, sol
