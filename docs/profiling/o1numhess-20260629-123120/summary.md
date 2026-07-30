# O1NumHess Profiling Summary

Generated: 2026-06-29T12:33:43

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.160 | 1.083 | 0.069 | 0.927 | 39 | 0 |
| caffeine | 28 | 0.301 | 0.222 | 0.011 | 0.131 | 39 | 0 |
| taxol | 1 | 115.849 | 113.691 | 6.937 | 90.156 | 114 | 0 |
| taxol | 28 | 20.879 | 19.067 | 0.609 | 8.149 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.160 | 2.073 | -44.1% | 1.135 | +2.2% | 1.083 | 1.880 | 0.943 |
| caffeine | 28 | 0.301 | 1.442 | -79.1% | 0.501 | -39.8% | 0.222 | 1.224 | 0.313 |
| taxol | 1 | 115.849 | 87.180 | +32.9% | 64.550 | +79.5% | 113.691 | 72.708 | 50.259 |
| taxol | 28 | 20.879 | 26.923 | -22.4% | 23.360 | -10.6% | 19.067 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 115.849 | 60.000 | +93.1% | 55.000 | +110.6% |
| 28 | 20.879 | 22.000 | -5.1% | 10.000 | +108.8% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.003 | 0.028 | 0.020 | 0.018 | 0.000 |
| caffeine | 28 | 0.011 | 1440 | 29.440 | 39 | 0.002 | 0.020 | 0.032 | 0.011 | 0.001 |
| taxol | 1 | 6.936 | 15348 | 61.370 | 108 | 0.077 | 3.013 | 2.634 | 1.193 | 0.003 |
| taxol | 28 | 0.609 | 15348 | 61.370 | 108 | 0.161 | 4.787 | 5.358 | 2.333 | 0.015 |
