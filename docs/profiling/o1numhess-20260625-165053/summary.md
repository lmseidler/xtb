# O1NumHess Profiling Summary

Generated: 2026-06-25T16:53:17

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.415 | 1.339 | 0.068 | 1.037 | 40 | 1 |
| caffeine | 28 | 0.483 | 0.407 | 0.011 | 0.178 | 40 | 1 |
| taxol | 1 | 121.133 | 118.868 | 7.799 | 96.678 | 117 | 2 |
| taxol | 28 | 20.433 | 18.683 | 0.621 | 8.771 | 117 | 2 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.415 | 2.073 | -31.7% | 1.135 | +24.7% | 1.339 | 1.880 | 0.943 |
| caffeine | 28 | 0.483 | 1.442 | -66.5% | 0.501 | -3.7% | 0.407 | 1.224 | 0.313 |
| taxol | 1 | 121.133 | 87.180 | +38.9% | 64.550 | +87.7% | 118.868 | 72.708 | 50.259 |
| taxol | 28 | 20.433 | 26.923 | -24.1% | 23.360 | -12.5% | 18.683 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 121.133 | 60.000 | +101.9% | 55.000 | +120.2% |
| 28 | 20.433 | 22.000 | -7.1% | 10.000 | +104.3% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.068 | 1440 | 29.440 | 39 | 0.002 | 0.033 | 0.020 | 0.013 | 0.000 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.005 | 0.079 | 0.055 | 0.039 | 0.000 |
| taxol | 1 | 7.798 | 15762 | 62.850 | 108 | 0.068 | 3.573 | 2.791 | 1.337 | 0.008 |
| taxol | 28 | 0.620 | 15762 | 62.850 | 108 | 0.147 | 5.743 | 5.620 | 2.472 | 0.018 |
