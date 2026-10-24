import numpy as np

from gnssint.clock_monitor import calibrate, estimate_measurement_variance, run_clock_monitor


def synthetic_clock(n=240, dt=30.0, noise=0.4, seed=3):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * dt
    bias = -7.7e4 + 418.1 * t + 0.5 * 6e-4 * t ** 2       # like the sample receiver
    return t, bias + rng.normal(0.0, noise, n)


def test_noise_estimate_from_second_differences():
    _, z = synthetic_clock(noise=0.4)
    assert abs(np.sqrt(estimate_measurement_variance(z)) - 0.4) < 0.05


def test_calibrated_monitor_is_quiet_on_clean_data():
    t, z = synthetic_clock()
    r = estimate_measurement_variance(z[:120])
    cfg = calibrate(t[:120], z[:120], r)
    trace = run_clock_monitor(t, z, r, cfg)
    assert 0.5 < np.nanmean(trace.nis[120:]) < 2.0
    assert np.nansum(trace.alarm[120:]) <= 1


def test_time_push_is_caught_quickly():
    t, z = synthetic_clock()
    r = estimate_measurement_variance(z[:120])
    cfg = calibrate(t[:120], z[:120], r)
    start = t[150]
    pushed = z + np.maximum(0.0, t - start) * 0.5            # 0.5 m/s = 1.7 ns/s
    trace = run_clock_monitor(t, pushed, r, cfg)
    first = np.where(trace.alarm & (t > start))[0]
    assert len(first) and t[first[0]] - start <= 60.0
