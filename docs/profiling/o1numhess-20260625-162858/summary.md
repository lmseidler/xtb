# O1NumHess Profiling Summary

Generated: 2026-06-25T16:32:49

## Current Runs

| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 2.092 | 1.943 | 0.298 | 1.509 | 62 | 0 |
| caffeine | 28 | 0.462 | 0.340 | 0.037 | 0.195 | 62 | 0 |
| taxol | 1 | 173.309 | 130.862 | 12.922 | 109.627 | 131 | 0 |
| taxol | 28 | 50.167 | 16.970 | 1.118 | 10.167 | 131 | 0 |

## Baseline Comparison

| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 2.092 | 2.073 | +0.9% | 1.135 | +84.3% | 1.943 | 1.880 | 0.943 |
| caffeine | 28 | 0.462 | 1.442 | -68.0% | 0.501 | -7.8% | 0.340 | 1.224 | 0.313 |
| taxol | 1 | 173.309 | 87.180 | +98.8% | 64.550 | +168.5% | 130.862 | 72.708 | 50.259 |
| taxol | 28 | 50.167 | 26.923 | +86.3% | 23.360 | +114.8% | 16.970 | 19.116 | 15.428 |

## Taxol Plan Targets

| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |
|---:|---:|---:|---:|---:|---:|
| 1 | 173.309 | 60.000 | +188.8% | 55.000 | +215.1% |
| 28 | 50.167 | 22.000 | +128.0% | 10.000 | +401.7% |

## Direction Generation Breakdown

| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.298 | 2214 | 41.600 | 60 | 0.002 | 0.118 | 0.103 | 0.073 | 0.000 |
| caffeine | 28 | 0.037 | 2214 | 41.600 | 60 | 0.013 | 0.195 | 0.170 | 0.110 | 0.000 |
| taxol | 1 | 12.921 | 18192 | 71.610 | 126 | 0.110 | 5.747 | 4.953 | 2.083 | 0.013 |
| taxol | 28 | 1.116 | 18192 | 71.610 | 126 | 0.213 | 9.366 | 10.071 | 4.058 | 0.012 |
