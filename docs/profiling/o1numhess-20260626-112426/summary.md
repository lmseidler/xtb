# O1NumHess Profiling Summary

Generated: 2026-06-26T11:26:53

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.238 | 1.165 | 0.069 | 0.969 | 39 | 0 |
| caffeine | 28 | 0.375 | 0.298 | 0.010 | 0.163 | 39 | 0 |
| taxol | 1 | 120.797 | 118.677 | 6.972 | 92.585 | 114 | 0 |
| taxol | 28 | 24.336 | 22.510 | 0.629 | 9.133 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.238 | 2.073 | -40.3% | 1.135 | +9.1% | 1.165 | 1.880 | 0.943 |
| caffeine | 28 | 0.375 | 1.442 | -74.0% | 0.501 | -25.1% | 0.298 | 1.224 | 0.313 |
| taxol | 1 | 120.797 | 87.180 | +38.6% | 64.550 | +87.1% | 118.677 | 72.708 | 50.259 |
| taxol | 28 | 24.336 | 26.923 | -9.6% | 23.360 | +4.2% | 22.510 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 120.797 | 60.000 | +101.3% | 55.000 | +119.6% |
| 28 | 24.336 | 22.000 | +10.6% | 10.000 | +143.4% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.001 | 0.027 | 0.024 | 0.015 | 0.001 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.000 | 0.043 | 0.044 | 0.045 | 0.001 |
| taxol | 1 | 6.972 | 15348 | 61.370 | 108 | 0.067 | 3.018 | 2.616 | 1.241 | 0.007 |
| taxol | 28 | 0.628 | 15348 | 61.370 | 108 | 0.168 | 4.995 | 5.255 | 2.402 | 0.014 |
