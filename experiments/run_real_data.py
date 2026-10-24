"""Run the integrity monitors on real GPS data and three injected scenarios.

    python experiments/run_real_data.py            # uses data/07590920.05o/.05n

Outputs results/real_data_scenarios.png and results/real_data_summary.md.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssint import faults                                              # noqa: E402
from gnssint.clock_monitor import (calibrate, estimate_measurement_variance,  # noqa: E402
                                   run_clock_monitor)
from gnssint.orbits import ecef_to_geodetic, enu_rotation              # noqa: E402
from gnssint.pipeline import calibrate_sigma, process                  # noqa: E402
from gnssint.raim import RaimConfig, raim_check                        # noqa: E402
from gnssint.rinex import read_nav, read_obs                           # noqa: E402
from gnssint.solver import prepare_epoch, wls                          # noqa: E402

OBS = ROOT / "data" / "07590920.05o"
NAV = ROOT / "data" / "07590920.05n"
RESULTS = ROOT / "results"

RAIM = RaimConfig(p_fa=1e-5, p_md=1e-3)
FAULT_START_S = 1800.0          # inject 30 min into the hour
CALIBRATION_EPOCHS = 60         # clock model tuned on the first 30 min only


def arrays(recs, t0):
    return dict(
        t=np.array([r.t - t0 for r in recs]),
        ratio=np.array([r.test_stat / r.threshold for r in recs]),
        alarm=np.array([r.alarm for r in recs]),
        usable=np.array([r.usable for r in recs]),
        herr=np.array([np.hypot(*r.enu_error[:2]) for r in recs]),
        verr=np.array([abs(r.enu_error[2]) for r in recs]),
        hpl=np.array([r.hpl for r in recs]),
        vpl=np.array([r.vpl for r in recs]),
        clock=np.array([r.clock for r in recs]),
        clock_var=np.array([r.clock_var for r in recs]),
        excluded=[r.excluded for r in recs],
    )


def first_true(mask, t, after):
    idx = np.where(mask & (t >= after))[0]
    return None if len(idx) == 0 else idx[0]


def main():
    if not OBS.exists():
        sys.exit("Sample data missing: run  python scripts/get_sample_data.py")
    RESULTS.mkdir(exist_ok=True)
    obs, nav = read_obs(OBS), read_nav(NAV)
    ref = obs.approx_xyz
    t0 = obs.epochs[0].t
    lat, lon, _ = ecef_to_geodetic(ref)
    R = enu_rotation(lat, lon)
    lines = ["# Real-data results (station 0759, GSI Japan, 2005-04-02, 1 h, 30 s)", ""]

    # 1. Calibrate the pseudorange error model on the fault-free data -------------
    a, b = calibrate_sigma(obs, nav, ref)
    lines += [f"Calibrated error model (1.25x overbound): sigma^2 = {a:.2f}^2 + ({b:.2f} / sin el)^2 m^2", ""]

    # 2. Nominal run ----------------------------------------------------------
    nom = arrays(process(obs, nav, ref, (a, b), RAIM), t0)
    lines += ["## Nominal (no injected faults)",
              f"- RAIM alarms: {int(nom['alarm'].sum())} of {len(nom['t'])} epochs "
              f"(P_fa = {RAIM.p_fa:g}); max T/T_D = {nom['ratio'].max():.2f}",
              f"- Horizontal error: mean {nom['herr'].mean():.2f} m, max {nom['herr'].max():.2f} m; "
              f"HPL {nom['hpl'].min():.1f}-{nom['hpl'].max():.1f} m",
              f"- Vertical error: mean {nom['verr'].mean():.2f} m, max {nom['verr'].max():.2f} m; "
              f"VPL {nom['vpl'].min():.1f}-{nom['vpl'].max():.1f} m", ""]

    # Clock model: measurement noise and Q tuned on the first half only, then frozen.
    cal = slice(0, CALIBRATION_EPOCHS)
    r_clk = estimate_measurement_variance(nom["clock"][cal])
    ccfg = calibrate(nom["t"][cal], nom["clock"][cal], r_clk)
    clk_nom = run_clock_monitor(nom["t"], nom["clock"], r_clk, ccfg)
    fa_epoch = int(np.nansum(clk_nom.epoch_alarm[CALIBRATION_EPOCHS:]))
    fa_window = int(np.nansum(clk_nom.window_alarm[CALIBRATION_EPOCHS:]))
    lines += ["## Clock monitor calibration",
              f"- Clock-estimate noise from second differences: {np.sqrt(r_clk):.2f} m "
              f"(the WLS covariance says {np.sqrt(np.median(nom['clock_var'])):.2f} m, because RAIM's "
              "sigmas are deliberately overbounded)",
              f"- Process noise tuned on the first {CALIBRATION_EPOCHS} epochs (q_drift = {ccfg.q_drift:.1e}, "
              f"q_accel = {ccfg.q_accel:.1e}); on the untouched second half: mean NIS "
              f"{np.nanmean(clk_nom.nis[CALIBRATION_EPOCHS:]):.2f}, false alarms {fa_epoch} per-epoch, "
              f"{fa_window} window", ""]

    # 3. Scenario A: slow ramp on the most dangerous satellite (largest slope) ---
    ep_start = next(e for e in obs.epochs if e.t - t0 >= FAULT_START_S)
    meas = prepare_epoch(ep_start, nav, ref, sigma_ab=(a, b))
    sol = wls(meas.sat_pos, meas.pr, meas.sigma, x0=np.r_[ref, 0.0])
    chk = raim_check(sol.G, meas.sigma, sol.residuals, R, RAIM)
    always = [s for s in meas.sats if all(s in e.obs for e in obs.epochs)]
    sat = max(always, key=lambda s: chk.hslope[meas.sats.index(s)])
    rate_a = 0.2
    inj_a = lambda m, _R: faults.single_satellite_ramp(m, sat, t0 + FAULT_START_S, rate_a)
    fde_a = arrays(process(obs, nav, ref, (a, b), RAIM, inject=inj_a), t0)
    raw_a = arrays(process(obs, nav, ref, (a, b), RAIM, inject=inj_a, exclude=False), t0)
    k = first_true(fde_a["alarm"], fde_a["t"], FAULT_START_S)
    lines += [f"## Scenario A: {sat} ramps at {rate_a} m/s from t = {FAULT_START_S / 60:.0f} min"]
    if k is not None:
        dt = fde_a["t"][k] - FAULT_START_S
        post = fde_a["t"] >= FAULT_START_S
        alarmed = post & fde_a["alarm"]
        right = sum(1 for i in np.where(alarmed)[0] if fde_a["excluded"][i] == sat)
        wrong = sum(1 for i in np.where(alarmed)[0] if fde_a["excluded"][i] not in (None, sat))
        unresolved = int((alarmed & ~fde_a["usable"]).sum())
        usable_err = fde_a["herr"][post & fde_a["usable"]]
        mi_a = int((post & ~fde_a["alarm"] & (fde_a["herr"] > fde_a["hpl"])).sum())
        lines += [f"- Detected after {dt:.0f} s, when the bias was {rate_a * dt:.1f} m",
                  f"- Alarmed epochs: {int(alarmed.sum())}; {sat} correctly excluded in {right}, wrong "
                  f"satellite excluded in {wrong}, ambiguous (alarm, no position output) in {unresolved}",
                  f"- Max horizontal error of any position actually output after onset: {usable_err.max():.1f} m "
                  f"(an unmonitored receiver reaches {raw_a['herr'][post].max():.0f} m); "
                  f"epochs with undetected error above HPL: {mi_a}", ""]
    else:
        lines += ["- Not detected within the hour", ""]

    # 4. Scenario B: time push (common-mode), sweep of ramp rates -----------------
    lines += ["## Scenario B: time push (every pseudorange shifted equally)",
              "| Ramp rate (m/s) | Max change in RAIM T | Clock alarm after | Offset at alarm |",
              "|---|---|---|---|"]
    trace_b = None
    for rate in [0.02, 0.05, 0.1, 0.2, 0.5, 2.0]:
        inj = lambda m, _R, r=rate: faults.time_push(m, t0 + FAULT_START_S, r)
        spoof = arrays(process(obs, nav, ref, (a, b), RAIM, inject=inj), t0)
        d_stat = np.max(np.abs(spoof["ratio"] - nom["ratio"]))
        tr = run_clock_monitor(spoof["t"], spoof["clock"], r_clk, ccfg)
        k = first_true(tr.alarm, spoof["t"], FAULT_START_S)
        if rate == 0.2:
            trace_b = tr
        if k is None:
            lines.append(f"| {rate} | {d_stat:.1e} | not detected | - |")
        else:
            dt = spoof["t"][k] - FAULT_START_S
            lines.append(f"| {rate} | {d_stat:.1e} | {dt:.0f} s | {rate * dt:.1f} m ({rate * dt / 0.2998:.0f} ns) |")
    lines.append("")

    # 5. Scenario C: consistent position pull-off with a clean clock -------------
    vel = (0.5, 0.0, 0.0)
    inj_c = lambda m, R_: faults.position_pull_off(m, ref, t0 + FAULT_START_S, vel, R_)
    pull = arrays(process(obs, nav, ref, (a, b), RAIM, inject=inj_c), t0)
    tr_c = run_clock_monitor(pull["t"], pull["clock"], r_clk, ccfg)
    after = pull["t"] >= FAULT_START_S
    mi = after & (pull["herr"] > pull["hpl"]) & ~pull["alarm"]
    lines += ["## Scenario C: position pull-off at 0.5 m/s east, clock left consistent",
              f"- RAIM alarms: {int(pull['alarm'][after].sum())}; clock alarms: "
              f"{int(np.nansum(tr_c.alarm[after]))}",
              f"- Horizontal error at end: {pull['herr'][-1]:.0f} m vs HPL {pull['hpl'][-1]:.0f} m; "
              f"epochs with undetected error > HPL (misleading information): {int(mi.sum())}",
              "- Neither snapshot RAIM nor the clock monitor can see this; it needs an independent",
              "  reference such as an inertial cross-check.", ""]

    # 6. Figure -----------------------------------------------------------------
    fig, ax = plt.subplots(3, 1, figsize=(8, 10), constrained_layout=True)
    m = lambda s: s / 60.0
    ax[0].semilogy(m(nom["t"]), nom["ratio"], label="nominal")
    ax[0].semilogy(m(fde_a["t"]), fde_a["ratio"], label=f"{sat} ramp {rate_a} m/s")
    ax[0].axhline(1.0, color="k", ls="--", lw=1, label="threshold")
    ax[0].axvline(m(FAULT_START_S), color="grey", lw=0.8)
    ax[0].set(title="A. RAIM test statistic / threshold", ylabel="T / T_D")
    ax[0].legend(fontsize=8)

    ax[1].semilogy(m(clk_nom.t), clk_nom.nis, label="nominal")
    ax[1].semilogy(m(trace_b.t), trace_b.nis, label="time push 0.2 m/s")
    from scipy.stats import chi2
    ax[1].axhline(chi2.isf(ccfg.p_fa, 1), color="k", ls="--", lw=1, label="threshold")
    ax[1].axvline(m(FAULT_START_S), color="grey", lw=0.8)
    ax[1].set(title="B. Clock monitor: normalised innovation squared (RAIM sees nothing here)",
              ylabel="NIS")
    ax[1].legend(fontsize=8, loc="lower left")

    ax[2].plot(m(pull["t"]), pull["herr"], label="horizontal error")
    ax[2].plot(m(pull["t"]), pull["hpl"], label="HPL")
    ax[2].axvline(m(FAULT_START_S), color="grey", lw=0.8)
    ax[2].set(title="C. Position pull-off: error escapes HPL with no alarm",
              xlabel="minutes", ylabel="metres")
    ax[2].legend(fontsize=8)
    fig.savefig(RESULTS / "real_data_scenarios.png", dpi=130)

    (RESULTS / "real_data_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
