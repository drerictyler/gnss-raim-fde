# Monte Carlo validation (7 satellites, 3 degrees of freedom, P_fa = 0.001, P_md = 0.001)

## False alerts
- 400,000 fault-free trials: empirical P_fa = 9.92e-04 (95% CI 8.95e-04-1.09e-03) vs design 1e-03

## Detection probability vs bias (empirical vs non-central chi-square theory)
| Satellite | Bias (m) | Empirical P_d | Theory P_d |
|---|---|---|---|
| G11 (largest slope) | 0.0 | 0.001 | 0.001 |
| G11 (largest slope) | 5.5 | 0.646 | 0.646 |
| G11 (largest slope) | 11.0 | 1.000 | 1.000 |
| G11 (largest slope) | 16.4 | 1.000 | 1.000 |
| G08 (smallest slope) | 0.0 | 0.001 | 0.001 |
| G08 (smallest slope) | 10.7 | 0.649 | 0.646 |
| G08 (smallest slope) | 21.4 | 1.000 | 1.000 |
| G08 (smallest slope) | 32.2 | 1.000 | 1.000 |

- Largest |empirical - theory| over all points: 0.003 (20,000 trials per point, so ~0.01 is sampling noise)

## Protection level
- A 9.1 m bias on G11 shifts the horizontal position by 7.9 m, which equals the HPL (7.9 m at these settings)
- It goes undetected in 1.01e-03 of trials (95% CI 9.09e-04-1.11e-03) vs P_md = 1e-03
- The slope HPL bounds the fault-induced shift only; noise adds on top, which is why
  some formulations add a fault-free term.

