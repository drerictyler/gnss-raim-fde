import numpy as np

from gnssint.orbits import sat_position_clock


def test_observation_file_parsed(data):
    obs = data["obs"]
    assert obs.obs_types == ["L1", "C1", "L2", "P2"]
    assert len(obs.epochs) == 120            # includes epochs after three file-splice event records
    assert all(len(e.obs) >= 7 for e in obs.epochs)
    np.testing.assert_allclose(obs.approx_xyz, [-3976219.5082, 3382372.5671, 3652512.9849])


def test_broadcast_orbits_are_physical(data):
    for ephs in data["nav"].ephemerides.values():
        for eph in ephs:
            pos, clock = sat_position_clock(eph, eph.toe)
            radius_km = np.linalg.norm(pos) / 1e3
            assert 26_000 < radius_km < 27_100          # GPS orbital radius ~26,560 km
            assert abs(clock) < 1e-3                     # satellite clock offsets are sub-millisecond
