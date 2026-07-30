# O1NumHess Profiling Summary

Generated: 2026-06-26T11:12:10

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.212 | 1.126 | 0.069 | 0.970 | 39 | 0 |
| caffeine | 28 | 0.342 | 0.267 | 0.010 | 0.168 | 39 | 0 |
| taxol | 1 | 114.826 | 112.700 | 6.953 | 92.738 | 113 | 0 |
| taxol | 28 | 17.859 | 16.082 | 0.614 | 8.886 | 113 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 1.212 | 2.073 | -41.5% | 1.135 | +6.8% | 1.126 | 1.880 | 0.943 |
| caffeine | 28 | 0.342 | 1.442 | -76.3% | 0.501 | -31.8% | 0.267 | 1.224 | 0.313 |
| taxol | 1 | 114.826 | 87.180 | +31.7% | 64.550 | +77.9% | 112.700 | 72.708 | 50.259 |
| taxol | 28 | 17.859 | 26.923 | -33.7% | 23.360 | -23.5% | 16.082 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 114.826 | 60.000 | +91.4% | 55.000 | +108.8% |
| 28 | 17.859 | 22.000 | -18.8% | 10.000 | +78.6% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.069 | 1440 | 29.440 | 39 | 0.001 | 0.030 | 0.021 | 0.016 | 0.000 |
| caffeine | 28 | 0.010 | 1440 | 29.440 | 39 | 0.002 | 0.071 | 0.080 | 0.023 | 0.000 |
| taxol | 1 | 6.952 | 15348 | 61.370 | 108 | 0.059 | 2.967 | 2.676 | 1.224 | 0.009 |
| taxol | 28 | 0.614 | 15348 | 61.370 | 108 | 0.137 | 4.838 | 5.255 | 2.283 | 0.019 |
