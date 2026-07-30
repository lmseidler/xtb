# O1NumHess Profiling Summary

Generated: 2026-06-26T10:57:35

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.202 | 1.128 | 0.069 | 0.972 | 39 | 0 |
| caffeine | 28 | 0.329 | 0.258 | 0.010 | 0.164 | 39 | 0 |
| taxol | 1 | 127.361 | 125.147 | 7.561 | 93.977 | 118 | 2 |
| taxol | 28 | 22.361 | 20.605 | 0.659 | 8.869 | 118 | 2 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.202 | 2.073 | -42.0% | 1.135 | +5.9% | 1.128 | 1.880 | 0.943 |
| caffeine | 28 | 0.329 | 1.442 | -77.2% | 0.501 | -34.3% | 0.258 | 1.224 | 0.313 |
| taxol | 1 | 127.361 | 87.180 | +46.1% | 64.550 | +97.3% | 125.147 | 72.708 | 50.259 |
| taxol | 28 | 22.361 | 26.923 | -16.9% | 23.360 | -4.3% | 20.605 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 127.361 | 60.000 | +112.3% | 55.000 | +131.6% |
| 28 | 22.361 | 22.000 | +1.6% | 10.000 | +123.6% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.000 | 0.029 | 0.022 | 0.017 | 0.001 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.000 | 0.021 | 0.023 | 0.036 | 0.001 |
| taxol | 1 | 7.561 | 15762 | 62.850 | 108 | 0.081 | 3.328 | 2.801 | 1.327 | 0.009 |
| taxol | 28 | 0.658 | 15762 | 62.850 | 108 | 0.203 | 5.321 | 5.726 | 2.571 | 0.022 |
