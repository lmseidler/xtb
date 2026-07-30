# O1NumHess Profiling Summary

Generated: 2026-06-26T10:53:44

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.207 | 1.133 | 0.069 | 0.976 | 39 | 0 |
| caffeine | 28 | 0.342 | 0.262 | 0.010 | 0.164 | 39 | 0 |
| taxol | 1 | 106.055 | 103.922 | 7.556 | 93.427 | 108 | 5 |
| taxol | 28 | 14.595 | 12.803 | 0.678 | 9.276 | 108 | 5 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.207 | 2.073 | -41.8% | 1.135 | +6.4% | 1.133 | 1.880 | 0.943 |
| caffeine | 28 | 0.342 | 1.442 | -76.3% | 0.501 | -31.7% | 0.262 | 1.224 | 0.313 |
| taxol | 1 | 106.055 | 87.180 | +21.7% | 64.550 | +64.3% | 103.922 | 72.708 | 50.259 |
| taxol | 28 | 14.595 | 26.923 | -45.8% | 23.360 | -37.5% | 12.803 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 106.055 | 60.000 | +76.8% | 55.000 | +92.8% |
| 28 | 14.595 | 22.000 | -33.7% | 10.000 | +45.9% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.001 | 0.031 | 0.025 | 0.011 | 0.000 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.000 | 0.049 | 0.054 | 0.035 | 0.000 |
| taxol | 1 | 7.556 | 15762 | 62.850 | 108 | 0.068 | 3.251 | 2.912 | 1.295 | 0.012 |
| taxol | 28 | 0.677 | 15762 | 62.850 | 108 | 0.167 | 5.380 | 5.721 | 2.563 | 0.017 |
