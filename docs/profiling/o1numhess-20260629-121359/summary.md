# O1NumHess Profiling Summary

Generated: 2026-06-29T12:16:30

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.168 | 1.089 | 0.070 | 0.932 | 39 | 0 |
| caffeine | 28 | 0.320 | 0.242 | 0.012 | 0.146 | 39 | 0 |
| taxol | 1 | 116.127 | 113.967 | 6.939 | 90.403 | 114 | 0 |
| taxol | 28 | 20.992 | 19.165 | 0.615 | 8.241 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.168 | 2.073 | -43.6% | 1.135 | +2.9% | 1.089 | 1.880 | 0.943 |
| caffeine | 28 | 0.320 | 1.442 | -77.8% | 0.501 | -36.1% | 0.242 | 1.224 | 0.313 |
| taxol | 1 | 116.127 | 87.180 | +33.2% | 64.550 | +79.9% | 113.967 | 72.708 | 50.259 |
| taxol | 28 | 20.992 | 26.923 | -22.0% | 23.360 | -10.1% | 19.165 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 116.127 | 60.000 | +93.5% | 55.000 | +111.1% |
| 28 | 20.992 | 22.000 | -4.6% | 10.000 | +109.9% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.070 | 1440 | 29.440 | 39 | 0.001 | 0.030 | 0.019 | 0.019 | 0.000 |
| caffeine | 28 | 0.012 | 1440 | 29.440 | 39 | 0.004 | 0.069 | 0.065 | 0.037 | 0.000 |
| taxol | 1 | 6.939 | 15348 | 61.370 | 108 | 0.079 | 2.963 | 2.652 | 1.218 | 0.006 |
| taxol | 28 | 0.614 | 15348 | 61.370 | 108 | 0.150 | 4.890 | 5.236 | 2.328 | 0.017 |
