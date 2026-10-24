import numpy as np

from gnssint.corrections import sigma_model
from gnssint.pipeline import SIGMA_INFLATION, calibrate_sigma


def test_sigma_inflation_applies_to_sigma(data):
    # The returned coefficients include the inflation, so sigma (not the variance) scales by 1.25.
    el = np.radians(30.0)
    s_with = sigma_model(el, *data["sigma_ab"])
    s_without = sigma_model(el, *calibrate_sigma(data["obs"], data["nav"], data["ref"], inflation=1.0))
    assert SIGMA_INFLATION == 1.25
    assert np.isclose(s_with / s_without, 1.25)
