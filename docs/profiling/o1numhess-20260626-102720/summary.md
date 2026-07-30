# O1NumHess Profiling Summary

Generated: 2026-06-26T10:29:45

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.339 | 1.265 | 0.069 | 0.979 | 40 | 1 |
| caffeine | 28 | 0.467 | 0.396 | 0.013 | 0.169 | 40 | 1 |
| taxol | 1 | 117.666 | 115.525 | 7.571 | 93.594 | 117 | 2 |
| taxol | 28 | 20.585 | 18.809 | 0.658 | 8.870 | 117 | 2 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.339 | 2.073 | -35.4% | 1.135 | +18.0% | 1.265 | 1.880 | 0.943 |
| caffeine | 28 | 0.467 | 1.442 | -67.6% | 0.501 | -6.9% | 0.396 | 1.224 | 0.313 |
| taxol | 1 | 117.666 | 87.180 | +35.0% | 64.550 | +82.3% | 115.525 | 72.708 | 50.259 |
| taxol | 28 | 20.585 | 26.923 | -23.5% | 23.360 | -11.9% | 18.809 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 117.666 | 60.000 | +96.1% | 55.000 | +113.9% |
| 28 | 20.585 | 22.000 | -6.4% | 10.000 | +105.9% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.002 | 0.027 | 0.022 | 0.017 | 0.000 |
| caffeine | 28 | 0.013 | 1440 | 29.440 | 39 | 0.004 | 0.024 | 0.037 | 0.028 | 0.000 |
| taxol | 1 | 7.570 | 15762 | 62.850 | 108 | 0.082 | 3.301 | 2.849 | 1.296 | 0.016 |
| taxol | 28 | 0.657 | 15762 | 62.850 | 108 | 0.237 | 5.323 | 5.626 | 2.565 | 0.013 |
