# vdW Neighbor O1NumHess Profiling Notes

Run date: 2026-06-25

Command:

```bash
python skills/o1numhess-profiling/scripts/profile_o1numhess.py \
  --repo /home/polt/xtb \
  --baseline /home/polt/compare/profiling/results.txt
```

The runner used `--hess --o1nh` and a clean rebuild was performed first to
avoid stale Fortran `.mod` ABI mismatches.

## Headline

The vdW-radius direction neighborhoods reduce direction count and direction
generation cost:

| Case | Threads | Previous current dirgen | vdW dirgen | Old xtb dirgen | wicz dirgen |
|---|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.298 | 0.069 | 0.695 | 0.067 |
| caffeine | 28 | 0.037 | 0.016 | 0.150 | 0.068 |
| taxol | 1 | 12.922 | 8.799 | 23.267 | 5.673 |
| taxol | 28 | 1.118 | 0.661 | 3.959 | 8.650 |

Taxol direction count dropped from 131 to 117, close to wicz's 114.
Caffeine direction count dropped from 62 to 40, slightly below wicz's 42.

## Caveat

The Hessian-only runs now report imaginary frequencies:

| Case | Threads | Imaginary frequencies |
|---|---:|---:|
| caffeine | 1 | 1 |
| caffeine | 28 | 1 |
| taxol | 1 | 2 |
| taxol | 28 | 2 |

So the vdW neighborhood change improves direction generation performance, but
the current combined implementation is not yet quality-equivalent to the prior
0-imaginary-frequency runs.

## Remaining Bottleneck

Taxol 1T remains slower than old xtb/wicz overall because `grad_evals` is still
high:

| Case | Threads | grad_evals |
|---|---:|---:|
| taxol | 1 | 99.198 |
| taxol | 28 | 8.970 |

Direction generation is no longer the dominant taxol cost.
