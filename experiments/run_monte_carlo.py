"""Monte Carlo check that the RAIM implementation matches its own theory.

Uses the real satellite geometry and calibrated error model from the sample data
at t = 30 min, then simulates Gaussian pseudorange errors:

  1. False-alert rate: fault-free trials, empirical P_fa vs the design value.
  2. Detection probability vs bias on the satellite with the largest horizontal
     slope (the most position error per unit of detectability) and the one with
     the smallest, against the non-central chi-square prediction.
  3. The protection-level definition: at the bias the HPL is built from, the
     missed-detection rate should equal P_md.

    python experiments/run_monte_carlo.py
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ncx2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnssint.orbits import ecef_to_geodetic, enu_rotation              # noqa: E402
from gnssint.pipeline import calibrate_sigma                           # noqa: E402
from gnssint.raim import RaimConfig, detection_threshold, pbias, raim_check  # noqa: E402
from gnssint.rinex import read_nav, read_obs                           # noqa: E402
from gnssint.solver import prepare_epoch, wls                          # noqa: E402

OBS = ROOT / "data" / "07590920.05o"
NAV = ROOT / "data" / "07590920.05n"
RESULTS = ROOT / "results"
CFG = RaimConfig(p_fa=1e-3, p_md=1e-3)     # loose enough to verify by simulation
RNG = np.random.default_rng(2026)


def simulate_stats(M, w, sigma, bias_vec, trials):
    """Test statistics for `trials` draws of y = e + bias, with r = M y."""
    e = RNG.standard_normal((trials, len(sigma))) * sigma
    r = (e + bias_vec) @ M.T
    return (r ** 2 * w).sum(axis=1)


def binomial_ci(k, n):
    p = k / n
    half = 1.96 * np.sqrt(max(p * (1 - p), 1e-12) / n)
    return p, max(p - half, 0.0), p + half


def main():
    if not OBS.exists():
        sys.exit("Sample data missing: run  python scripts/get_sample_data.py")
    RESULTS.mkdir(exist_ok=True)
    obs, nav = read_obs(OBS), read_nav(NAV)
    ref = obs.approx_xyz
    lat, lon, _ = ecef_to_geodetic(ref)
    R = enu_rotation(lat, lon)
    a, b = calibrate_sigma(obs, nav, ref)

    epoch = next(e for e in obs.epochs if e.t - obs.epochs[0].t >= 1800.0)
    meas = prepare_epoch(epoch, nav, ref, sigma_ab=(a, b))
    sol = wls(meas.sat_pos, meas.pr, meas.sigma, x0=np.r_[ref, 0.0])
    chk = raim_check(sol.G, meas.sigma, sol.residuals, R, CFG)

    G, sigma = sol.G, meas.sigma
    w = 1.0 / sigma ** 2
    S = np.linalg.solve(G.T @ (w[:, None] * G), G.T * w)
    M = np.eye(len(sigma)) - G @ S
    p_diag = np.einsum("ij,ji->i", G, S)
    dof = len(sigma) - 4
    thr = detection_threshold(dof, CFG.p_fa)
    pb = pbias(dof, thr, CFG.p_md)
    lines = [f"# Monte Carlo validation ({len(sigma)} satellites, {dof} degrees of freedom, "
             f"P_fa = {CFG.p_fa:g}, P_md = {CFG.p_md:g})", ""]

    # 1. False-alert rate ------------------------------------------------------
    n_fa = 400_000
    t_ff = simulate_stats(M, w, sigma, np.zeros(len(sigma)), n_fa)
    p, lo, hi = binomial_ci(int((t_ff > thr).sum()), n_fa)
    lines += ["## False alerts",
              f"- {n_fa:,} fault-free trials: empirical P_fa = {p:.2e} (95% CI {lo:.2e}-{hi:.2e}) "
              f"vs design {CFG.p_fa:.0e}", ""]

    # 2. Detection probability vs bias -----------------------------------------
    worst, mildest = int(np.argmax(chk.hslope)), int(np.argmin(chk.hslope))
    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    lines += ["## Detection probability vs bias (empirical vs non-central chi-square theory)",
              "| Satellite | Bias (m) | Empirical P_d | Theory P_d |", "|---|---|---|---|"]
    worst_gap = 0.0
    for idx, label in [(worst, "largest slope"), (mildest, "smallest slope")]:
        scale = sigma[idx] / np.sqrt(1.0 - p_diag[idx])     # bias per unit sqrt(lambda)
        biases = np.linspace(0.0, 2.2 * pb * scale, 12)
        emp, theory = [], []
        for bias in biases:
            vec = np.zeros(len(sigma))
            vec[idx] = bias
            stats = simulate_stats(M, w, sigma, vec, 20_000)
            emp.append((stats > thr).mean())
            theory.append(ncx2.sf(thr, dof, (bias / scale) ** 2) if bias > 0 else CFG.p_fa)
        emp, theory = np.array(emp), np.array(theory)
        worst_gap = max(worst_gap, np.abs(emp - theory).max())
        for bias, e_, t_ in list(zip(biases, emp, theory))[::3]:
            lines.append(f"| {meas.sats[idx]} ({label}) | {bias:.1f} | {e_:.3f} | {t_:.3f} |")
        line, = ax.plot(biases, theory, label=f"{meas.sats[idx]} ({label}) theory")
        ax.plot(biases, emp, "o", color=line.get_color(), ms=4, label=f"{meas.sats[idx]} simulated")
        ax.axvline(pb * scale, color=line.get_color(), ls=":", lw=1)
    lines += ["", f"- Largest |empirical - theory| over all points: {worst_gap:.3f} "
              "(20,000 trials per point, so ~0.01 is sampling noise)", ""]
    ax.set(xlabel="bias on one satellite (m)", ylabel="probability of detection",
           title="RAIM detection: simulation vs theory (dotted: bias used for protection levels)")
    ax.legend(fontsize=8)
    fig.savefig(RESULTS / "monte_carlo.png", dpi=130)

    # 3. What the HPL means ----------------------------------------------------
    scale = sigma[worst] / np.sqrt(1.0 - p_diag[worst])
    vec = np.zeros(len(sigma))
    vec[worst] = pb * scale
    n_md = 400_000
    missed = int((simulate_stats(M, w, sigma, vec, n_md) <= thr).sum())
    p, lo, hi = binomial_ci(missed, n_md)
    shift = np.hypot(*(R @ S[:3, worst] * vec[worst])[:2])
    lines += ["## Protection level",
              f"- A {vec[worst]:.1f} m bias on {meas.sats[worst]} shifts the horizontal position by "
              f"{shift:.1f} m, which equals the HPL ({chk.hpl:.1f} m at these settings)",
              f"- It goes undetected in {p:.2e} of trials (95% CI {lo:.2e}-{hi:.2e}) vs P_md = {CFG.p_md:.0e}",
              "- The slope HPL bounds the fault-induced shift only; noise adds on top, which is why",
              "  some formulations add a fault-free term.", ""]

    (RESULTS / "monte_carlo_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
