# O1NumHess Profiling Summary

Generated: 2026-06-26T12:40:17

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.202 | 1.125 | 0.069 | 0.970 | 39 | 0 |
| caffeine | 28 | 0.352 | 0.269 | 0.011 | 0.174 | 39 | 0 |
| taxol | 1 | 118.448 | 116.288 | 6.950 | 92.653 | 114 | 0 |
| taxol | 28 | 21.691 | 19.894 | 0.634 | 8.946 | 114 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.202 | 2.073 | -42.0% | 1.135 | +5.9% | 1.125 | 1.880 | 0.943 |
| caffeine | 28 | 0.352 | 1.442 | -75.6% | 0.501 | -29.7% | 0.269 | 1.224 | 0.313 |
| taxol | 1 | 118.448 | 87.180 | +35.9% | 64.550 | +83.5% | 116.288 | 72.708 | 50.259 |
| taxol | 28 | 21.691 | 26.923 | -19.4% | 23.360 | -7.1% | 19.894 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 118.448 | 60.000 | +97.4% | 55.000 | +115.4% |
| 28 | 21.691 | 22.000 | -1.4% | 10.000 | +116.9% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.001 | 0.037 | 0.019 | 0.012 | 0.000 |
| caffeine | 28 | 0.011 | 1440 | 29.440 | 39 | 0.007 | 0.067 | 0.044 | 0.047 | 0.000 |
| taxol | 1 | 6.950 | 15348 | 61.370 | 108 | 0.080 | 3.013 | 2.612 | 1.223 | 0.007 |
| taxol | 28 | 0.633 | 15348 | 61.370 | 108 | 0.146 | 4.980 | 5.289 | 2.300 | 0.013 |
