# Real-data results (station 0759, GSI Japan, 2005-04-02, 1 h, 30 s)

Calibrated error model (sigma inflated 1.25x): sigma^2 = 0.44^2 + (0.41 / sin el)^2 m^2

## Nominal (no injected faults)
- RAIM alarms: 0 of 120 epochs (P_fa = 1e-05); max T/T_D = 0.19
- Horizontal error: mean 0.40 m, max 0.91 m; HPL 8.9-24.2 m
- Vertical error: mean 0.53 m, max 2.23 m; VPL 13.6-65.3 m

## Clock monitor calibration
- Clock-estimate noise from second differences: 0.41 m (the WLS covariance says 1.37 m, because RAIM's sigmas are deliberately inflated)
- Process noise tuned on the first 60 epochs (q_drift = 3.2e-07, q_accel = 3.2e-09); on the untouched second half: mean NIS 0.85, false alarms 0 per-epoch, 0 window

## Scenario A: G11 ramps at 0.2 m/s from t = 30 min
- Detected after 60 s, when the bias was 12.0 m
- Alarmed epochs: 58; G11 correctly excluded in 51, wrong satellite excluded in 0, ambiguous (alarm, no position output) in 7
- Max horizontal error of any position actually output after onset: 6.2 m (an unmonitored receiver reaches 357 m); epochs with undetected error above HPL: 0

## Scenario B: time push (every pseudorange shifted equally)
| Ramp rate (m/s) | Max change in RAIM T | Clock alarm after | Offset at alarm |
|---|---|---|---|
| 0.02 | 7.2e-10 | not detected | - |
| 0.05 | 6.3e-10 | not detected | - |
| 0.1 | 1.0e-09 | 30 s | 3.0 m (10 ns) |
| 0.2 | 6.8e-10 | 30 s | 6.0 m (20 ns) |
| 0.5 | 2.9e-10 | 30 s | 15.0 m (50 ns) |
| 2.0 | 7.7e-10 | 30 s | 60.0 m (200 ns) |

## Scenario C: position pull-off at 0.5 m/s east, clock left consistent
- RAIM alarms: 0; clock alarms: 0
- Horizontal error at end: 885 m vs HPL 11 m; epochs with undetected error > HPL (misleading information): 59
- Neither snapshot RAIM nor the clock monitor can see this; it needs an independent
  reference such as an inertial cross-check.

