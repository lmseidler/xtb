# O1NumHess Profiling Summary

Generated: 2026-07-30T15:42:56

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.902 | 0.822 | 0.053 | 0.665 | 40 | 1 |
| caffeine | 28 | 0.372 | 0.263 | 0.017 | 0.116 | 40 | 1 |
| taxol | 1 | 49.651 | 48.670 | 4.015 | 37.172 | 116 | 0 |
| taxol | 28 | 7.891 | 7.206 | 0.830 | 3.616 | 116 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.902 | 2.073 | -56.5% | 1.135 | -20.5% | 0.822 | 1.880 | 0.943 |
| caffeine | 28 | 0.372 | 1.442 | -74.2% | 0.501 | -25.7% | 0.263 | 1.224 | 0.313 |
| taxol | 1 | 49.651 | 87.180 | -43.0% | 64.550 | -23.1% | 48.670 | 72.708 | 50.259 |
| taxol | 28 | 7.891 | 26.923 | -70.7% | 23.360 | -66.2% | 7.206 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 49.651 | 60.000 | -17.2% | 55.000 | -9.7% |
| 28 | 7.891 | 22.000 | -64.1% | 10.000 | -21.1% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.053 | 1440 | 29.440 | 39 | 0.001 | 0.032 | 0.005 | 0.014 | 0.001 |
| caffeine | 28 | 0.017 | 1440 | 29.440 | 39 | 0.003 | 0.040 | 0.015 | 0.023 | 0.000 |
| taxol | 1 | 4.015 | 15348 | 61.370 | 108 | 0.068 | 2.951 | 0.463 | 0.511 | 0.005 |
| taxol | 28 | 0.830 | 15348 | 61.370 | 108 | 0.150 | 4.354 | 0.807 | 0.797 | 0.010 |
