# O1NumHess Profiling Summary

Generated: 2026-06-29T12:41:13

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.162 | 1.083 | 0.069 | 0.928 | 39 | 0 |
| caffeine | 28 | 0.308 | 0.231 | 0.011 | 0.140 | 39 | 0 |
| taxol | 1 | 115.899 | 113.718 | 6.938 | 90.170 | 114 | 0 |
| taxol | 28 | 20.724 | 18.878 | 0.622 | 8.087 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.162 | 2.073 | -43.9% | 1.135 | +2.4% | 1.083 | 1.880 | 0.943 |
| caffeine | 28 | 0.308 | 1.442 | -78.6% | 0.501 | -38.5% | 0.231 | 1.224 | 0.313 |
| taxol | 1 | 115.899 | 87.180 | +32.9% | 64.550 | +79.5% | 113.718 | 72.708 | 50.259 |
| taxol | 28 | 20.724 | 26.923 | -23.0% | 23.360 | -11.3% | 18.878 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 115.899 | 60.000 | +93.2% | 55.000 | +110.7% |
| 28 | 20.724 | 22.000 | -5.8% | 10.000 | +107.2% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.001 | 0.026 | 0.021 | 0.020 | 0.001 |
| caffeine | 28 | 0.011 | 1440 | 29.440 | 39 | 0.001 | 0.029 | 0.046 | 0.045 | 0.001 |
| taxol | 1 | 6.937 | 15348 | 61.370 | 108 | 0.082 | 3.060 | 2.619 | 1.140 | 0.014 |
| taxol | 28 | 0.621 | 15348 | 61.370 | 108 | 0.203 | 4.982 | 5.187 | 2.372 | 0.009 |
