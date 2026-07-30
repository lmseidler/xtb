# O1NumHess Profiling Summary

Generated: 2026-06-26T15:00:26

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.154 | 1.078 | 0.069 | 0.923 | 39 | 0 |
| caffeine | 28 | 0.310 | 0.229 | 0.011 | 0.137 | 39 | 0 |
| taxol | 1 | 115.441 | 113.279 | 6.932 | 89.790 | 114 | 0 |
| taxol | 28 | 20.719 | 18.894 | 0.579 | 8.056 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.154 | 2.073 | -44.4% | 1.135 | +1.6% | 1.078 | 1.880 | 0.943 |
| caffeine | 28 | 0.310 | 1.442 | -78.5% | 0.501 | -38.1% | 0.229 | 1.224 | 0.313 |
| taxol | 1 | 115.441 | 87.180 | +32.4% | 64.550 | +78.8% | 113.279 | 72.708 | 50.259 |
| taxol | 28 | 20.719 | 26.923 | -23.0% | 23.360 | -11.3% | 18.894 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 115.441 | 60.000 | +92.4% | 55.000 | +109.9% |
| 28 | 20.719 | 22.000 | -5.8% | 10.000 | +107.2% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.000 | 0.030 | 0.018 | 0.021 | 0.000 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.003 | 0.029 | 0.019 | 0.020 | 0.001 |
| taxol | 1 | 6.931 | 15348 | 61.370 | 108 | 0.073 | 3.008 | 2.603 | 1.231 | 0.007 |
| taxol | 28 | 0.578 | 15348 | 61.370 | 108 | 0.152 | 4.891 | 5.217 | 2.201 | 0.006 |
