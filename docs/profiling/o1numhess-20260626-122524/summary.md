# O1NumHess Profiling Summary

Generated: 2026-06-26T12:28:28

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.282 | 1.207 | 0.069 | 0.967 | 39 | 0 |
| caffeine | 28 | 0.433 | 0.351 | 0.010 | 0.179 | 39 | 0 |
| taxol | 1 | 137.104 | 134.964 | 6.963 | 92.528 | 114 | 0 |
| taxol | 28 | 40.187 | 38.433 | 0.621 | 8.858 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.282 | 2.073 | -38.2% | 1.135 | +12.9% | 1.207 | 1.880 | 0.943 |
| caffeine | 28 | 0.433 | 1.442 | -70.0% | 0.501 | -13.6% | 0.351 | 1.224 | 0.313 |
| taxol | 1 | 137.104 | 87.180 | +57.3% | 64.550 | +112.4% | 134.964 | 72.708 | 50.259 |
| taxol | 28 | 40.187 | 26.923 | +49.3% | 23.360 | +72.0% | 38.433 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 137.104 | 60.000 | +128.5% | 55.000 | +149.3% |
| 28 | 40.187 | 22.000 | +82.7% | 10.000 | +301.9% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.002 | 0.031 | 0.023 | 0.012 | 0.000 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.007 | 0.064 | 0.076 | 0.026 | 0.000 |
| taxol | 1 | 6.963 | 15348 | 61.370 | 108 | 0.086 | 3.001 | 2.634 | 1.214 | 0.008 |
| taxol | 28 | 0.620 | 15348 | 61.370 | 108 | 0.129 | 4.876 | 5.138 | 2.365 | 0.018 |
