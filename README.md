# GNSS Integrity Lab

Receiver autonomous integrity monitoring (RAIM) and a receiver-clock spoofing monitor,
written from first principles in Python and tested on recorded GPS observations with
injected faults and spoofing attacks. It is a research and learning prototype, not a
certified navigation or integrity-monitoring system.

The aim is to show, with numbers, what each monitor can and cannot see:

| Threat | Snapshot RAIM | Clock monitor |
|---|---|---|
| Single-satellite fault (drifting range) | Detects; excludes the satellite when only one exclusion makes the rest consistent | Partly: the fault leaks into the clock estimate |
| Time push (all ranges shifted equally) | **Blind**: absorbed by the clock estimate | Detects tested ramps ≥ 0.1 m/s |
| Position pull-off with a consistent clock | **Blind** | **Blind** |

The last row is why operational receivers add an independent reference, such as an
inertial cross-check. These results come from one recording and the modelled scenarios;
they are not general performance guarantees.

## The experiment in plain English

The data is recorded, not live, and the receiver never loses signal. The project replays
one continuous hour of GPS measurements recorded in 2005 and computes a position every
30 seconds, 120 times in all. Each attack is a separate replay of that same hour:

- **Minutes 0–30, normal.** The measurements are used as recorded. This is the baseline:
  small errors and no alarms.
- **Minute 30 to the end, under attack.** The software adds a growing error to the recorded
  measurements until the last epoch, at 59.5 minutes. The question is whether each monitor
  notices, and how quickly.

The two monitors look for different things:

- **RAIM** checks whether the satellites agree with each other. A position has four
  unknowns (three coordinates and the receiver's clock offset), so with more than four
  satellites the solution is over-determined. One bad measurement makes the satellites
  disagree; RAIM detects that and, when it can tell which satellite is at fault, removes it.
- **The clock monitor** follows the receiver's clock estimate over time. A real oscillator
  drifts smoothly, so a sudden change in how fast the clock appears to run is suspicious.

**A. One bad satellite.** From minute 30, satellite G11's range grows by 0.2 m every second,
as it would with a drifting satellite clock or a single spoofed signal.
- *Expected:* the satellites now disagree, so RAIM should raise an alarm once the error is
  large compared with the noise, then remove G11. Part of the error leaks into the clock
  estimate, so the clock monitor may react too.
- *Result:* RAIM raised an alarm after 60 s (12 m of error) and removed G11. The worst
  position it output was 6.2 m off, against 357 m with no monitoring.

**B. Time push.** From minute 30, every satellite's range grows by the same amount each
second, as if the receiver's clock were being dragged.
- *Expected:* the satellites still agree, because an equal shift on every range looks
  exactly like a clock error. RAIM therefore cannot see it, and the position stays correct;
  what the attacker controls is the receiver's time. The clock monitor should flag the
  sudden change in clock rate if it is larger than the oscillator's natural drift.
- *Result:* RAIM saw nothing. The clock monitor raised an alarm within one epoch (30 s) for
  pushes of 0.1 m/s and faster; slower pushes hid in the oscillator's natural drift.

**C. Position pull-off.** From minute 30, the ranges are altered together so that they agree
with a false position moving away at 0.5 m/s, while the clock behaves normally.
- *Expected:* the satellites agree and the clock looks normal, so neither monitor should
  notice. Only an independent reference, such as an inertial sensor, could.
- *Result:* no alarms of either kind. By the end of the hour the position was 885 m off
  while the horizontal protection level reported 11 m.

![RAIM test statistic during attack A, clock-monitor statistic during attack B, and position error against the protection level during attack C](results/real_data_scenarios.png)

*Top: RAIM test statistic against its threshold during attack A. Middle: clock-monitor
statistic during attack B (RAIM sees nothing). Bottom: position error against the
protection level during attack C (log scale). In each panel the vertical line at 30 minutes marks the
start of the attack.*

## Data

One hour of GPS observations at 30 s from GSI station 0759 (Japan), 2 April 2005, with the
matching broadcast ephemeris, as distributed with
[RTKLIB](https://github.com/tomojitakasu/RTKLIB). The recording is dual-frequency; this
project uses the L1 pseudoranges only. Faults are added to the recorded pseudoranges, so
geometry, noise, multipath and the receiver's own clock behaviour are all real. The injection
works at the measurement level: it does not simulate a transmitter, signal acquisition or
tracking loops.

Position errors are measured against the station coordinates in the RINEX header
(APPROX POSITION XYZ, written by GSI, with a zero antenna offset); they are not derived
from this project's solutions.
<!-- TODO(Eric): CHECK 2, section 2.4. Compare the header value with GSI's published coordinates for station 0759 and, if they agree, say "surveyed" here. -->

## Results

**Fault-free performance** (L1 only, broadcast ionosphere, 10° elevation mask)

- Mean horizontal error 0.40 m (max 0.91 m); mean vertical error 0.53 m (max 2.2 m).
- No RAIM alarms in 120 epochs at a per-epoch false-alarm probability of 1e-5. One hour
  cannot confirm a probability that small; the Monte Carlo below tests the detector itself.
- Horizontal protection level 9–24 m, vertical 14–65 m; the observed errors stayed below
  both. These slope-based levels bound the fault-induced error only, not the total error.
- Error model calibrated from the fault-free data, with σ inflated by 1.25 (the
  variance by 1.5625): σ² = 0.44² + (0.41 / sin el)² m².

**A. Single-satellite ramp (G11, 0.2 m/s from t = 30 min)**

- Detected after 60 s, at 12 m of bias.
- G11 correctly excluded in 51 of 58 alarmed epochs, and a healthy satellite never. The
  other 7 epochs were ambiguous and were reported as alarms with no position rather than
  guessed.
- Largest error in any position actually output: 6.2 m, and no output error exceeded the
  protection level without an alarm. An unmonitored receiver reaches 357 m.

**B. Time push (every pseudorange shifted by the same ramp, from t = 30 min)**

- The RAIM test statistic changes by about 1e-9, which is numerically zero: the attack is
  absorbed by the clock unknown.
- The clock monitor raised an alarm within one epoch (30 s) for ramps of 0.1 m/s
  (0.33 ns/s) and above.
- Ramps of 0.05 m/s and below were not detected: they are smaller than the oscillator's
  natural wander in this recording. Detection depends on the receiver's clock, the monitor's
  tuning and the sampling, so these thresholds apply to this data only.
- Ramp rates tested: 0.02, 0.05, 0.1, 0.2, 0.5 and 2.0 m/s, each in a separate run on the
  clean recording.

**C. Position pull-off at 0.5 m/s with a consistent clock (from t = 30 min)**

- No RAIM alarms and no clock alarms; the error reaches 885 m while the HPL reports 11 m.

**Monte Carlo check against theory** (real geometry, 7 satellites, simulated Gaussian noise)

![Probability of detection against bias on one satellite: simulation versus non-central chi-square theory, for the largest- and smallest-slope satellites](results/monte_carlo.png)

This experiment uses a false-alarm probability of 1e-3 rather than 1e-5. About 400 false
alarms are then expected in 400,000 trials, enough to measure the rate; at 1e-5 there would
be about 4. Agreement validates the implementation against the assumed noise model, not
field performance.

- 400,000 fault-free trials: empirical P_fa = 9.9e-4 against a design value of 1e-3.
- Detection probability against bias matches the non-central χ² prediction to within 0.003.
- At the bias the protection level is built from, the fault goes undetected in 1.0e-3 of
  trials, matching P_md.

Full numbers: [`results/real_data_summary.md`](results/real_data_summary.md),
[`results/monte_carlo_summary.md`](results/monte_carlo_summary.md).

## Lessons from building it

1. **Exclusion can remove the wrong satellite.** The first version excluded the satellite
   whose removal gave the smallest test statistic. With six satellites, several subsets can
   pass, and at one epoch it removed a healthy satellite and output a 184 m error. Exclusion
   is now accepted only when exactly one subset passes; otherwise the epoch is an alarm.
2. **One noise model does not fit two monitors.** RAIM needs inflated sigmas, which put the
   clock estimate's standard deviation at 1.37 m. Its actual scatter is 0.41 m, so a clock
   monitor using the least-squares covariance was three times too insensitive. It now
   estimates its own noise from second differences.
3. **Oscillators drift.** A two-state clock model tuned on the first half-hour raised 11
   false alarms in the second, because the receiver's frequency offset moved by about 8 ppb
   over the hour (about 2.4 m/s in range-rate terms). Adding a third state (drift rate)
   removed them. Because the second half-hour informed that change, this is development
   evidence rather than independent validation; testing on other recordings is still to do.

## Method

| Module | What it does |
|---|---|
| `rinex.py` | Minimal RINEX 2.x GPS observation and navigation reader |
| `orbits.py` | Broadcast orbit and clock (IS-GPS-200), Earth-rotation correction, geodesy |
| `corrections.py` | Klobuchar ionosphere, Saastamoinen troposphere with the DO-229 mapping function, error model |
| `solver.py` | Measurement preparation and weighted least squares |
| `raim.py` | χ² fault detection, exclusion, slope-based protection levels |
| `clock_monitor.py` | Three-state clock Kalman filter with per-epoch and windowed innovation tests |
| `faults.py` | Satellite ramp, time push and position pull-off injection |
| `pipeline.py` | End-to-end processing and error-model calibration |

**RAIM.** With WLS sensitivity S = (GᵀWG)⁻¹GᵀW and residuals r = (I − GS)y, the statistic
T = rᵀWr is χ² with N − 4 degrees of freedom when fault-free. A bias b on satellite i makes
it non-central with λ = b²(1 − P_ii)/σ_i², where P = GS, and moves the position by S_i·b.
The protection level is the largest error per unit √λ over all satellites ("slope") times the
√λ that is missed with probability P_md. This is the classic slope form, which bounds the
fault-induced error only. The χ² result assumes Gaussian errors with weights equal to the
inverse variances; the non-centrality assumes a single biased satellite.

**Clock monitor.** The receiver clock is modelled as bias, drift and drift rate driven by
white noise. Its process noise is tuned on clean data so that the normalised innovation
squared (NIS) averages 1, then frozen. Alarms come from a per-epoch χ²(1) test on the NIS and
a 10-epoch χ²(10) window test. The window test assumes white innovations, which has not yet
been checked.

## Run it

```bash
python -m venv .venv
source .venv/bin/activate              # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python scripts/get_sample_data.py      # downloads the two RINEX files into data/
python experiments/run_real_data.py    # scenarios A-C, writes results/
python experiments/run_monte_carlo.py  # validation against theory
python -m pytest                       # 12 tests
```

Tested with Python 3.12; exact package versions in `requirements-lock.txt`. The Monte Carlo
uses a fixed seed, so its numbers reproduce exactly.

## Limitations

- GPS only, single frequency, snapshot RAIM. Advanced RAIM (ARAIM), with multiple
  constellations and dual frequency, is the modern equivalent and is not implemented.
- One hour from one station. The clock monitor was tuned on 30 minutes; a fielded monitor
  would be tuned over many days and temperatures.
- Faults are injected in software at the measurement level; real spoofing also distorts
  signal tracking, which this model does not capture.
- The exclusion rule avoided wrong exclusions in these scenarios; it does not guarantee
  correct exclusion for arbitrary faults.
- The innovation whiteness and tail probabilities behind the clock monitor's thresholds have
  not been checked.

## Next steps

- Check the clock monitor's innovation whiteness and tail behaviour.
- An Allan-variance clock test.
- An inertial consistency check for the pull-off case.
- Longer, multi-constellation datasets (IGS and CORS stations, the Norwegian Jammertest data).

## References

- B. W. Parkinson and P. Axelrad, "Autonomous GPS integrity monitoring using the pseudorange residual," *Navigation*, 1988.
- R. G. Brown, "A baseline GPS RAIM scheme and a note on the equivalence of three RAIM methods," *Navigation*, 1992.
- R. G. Brown and G. Y. Chin, "GPS RAIM: calculation of thresholds and protection radius using chi-square methods," ION, 1998.
- T. Walter and P. Enge, "Weighted RAIM for precision approach," *Proc. ION GPS-95*, 1995.
- R. G. Brown and P. Y. C. Hwang, *Introduction to Random Signals and Applied Kalman Filtering*.
- IS-GPS-200, Navstar GPS Space Segment / Navigation User Interfaces.
- RTCA DO-229, MOPS for GPS/SBAS airborne equipment (tropospheric mapping function).
- A. Joseph et al., "GNSS Radio Frequency Interference Mitigation in Collins Commercial Airborne Receivers," *Engineering Proceedings* 54, 2023 (ENC 2023): clock-based and inertial-consistency spoofing monitors that motivated scenarios B and C.

## Licence

Code: MIT (see `LICENSE`). Sample data courtesy of GSI Japan via RTKLIB; it remains under
its original terms.
