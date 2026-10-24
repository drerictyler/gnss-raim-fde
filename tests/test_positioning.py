import numpy as np

from gnssint.pipeline import process


def test_sub_metre_positioning_on_real_data(data):
    recs = process(data["obs"], data["nav"], data["ref"], data["sigma_ab"])
    herr = np.array([np.hypot(*r.enu_error[:2]) for r in recs])
    verr = np.array([abs(r.enu_error[2]) for r in recs])
    assert herr.mean() < 1.0 and herr.max() < 2.0
    assert verr.mean() < 1.5
    assert not any(r.alarm for r in recs)                # no false alarms on clean data
    assert all(r.hpl > np.hypot(*r.enu_error[:2]) for r in recs)
    assert all(r.vpl > abs(r.enu_error[2]) for r in recs)
