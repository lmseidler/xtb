# O1NumHess Profiling Summary

Generated: 2026-06-25T16:51:51

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.406 | 1.306 | 0.069 | 1.018 | 40 | 1 |
| caffeine | 28 | 0.578 | 0.480 | 0.016 | 0.223 | 40 | 1 |
| taxol | 1 | 124.526 | 122.376 | 8.799 | 99.198 | 117 | 2 |
| taxol | 28 | 20.531 | 18.751 | 0.661 | 8.970 | 117 | 2 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.406 | 2.073 | -32.2% | 1.135 | +23.8% | 1.306 | 1.880 | 0.943 |
| caffeine | 28 | 0.578 | 1.442 | -59.9% | 0.501 | +15.4% | 0.480 | 1.224 | 0.313 |
| taxol | 1 | 124.526 | 87.180 | +42.8% | 64.550 | +92.9% | 122.376 | 72.708 | 50.259 |
| taxol | 28 | 20.531 | 26.923 | -23.7% | 23.360 | -12.1% | 18.751 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 124.526 | 60.000 | +107.5% | 55.000 | +126.4% |
| 28 | 20.531 | 22.000 | -6.7% | 10.000 | +105.3% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.000 | 0.025 | 0.021 | 0.022 | 0.001 |
| caffeine | 28 | 0.016 | 1440 | 29.440 | 39 | 0.001 | 0.046 | 0.038 | 0.025 | 0.000 |
| taxol | 1 | 8.798 | 15762 | 62.850 | 108 | 0.115 | 4.067 | 3.192 | 1.388 | 0.012 |
| taxol | 28 | 0.660 | 15762 | 62.850 | 108 | 0.193 | 5.405 | 5.854 | 2.514 | 0.012 |
