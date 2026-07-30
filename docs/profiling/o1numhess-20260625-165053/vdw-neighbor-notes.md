# vdW Neighbor Profiling Notes

Run date: 2026-06-25

Command:

```bash
python skills/o1numhess-profiling/scripts/profile_o1numhess.py \
  --repo /home/polt/xtb \
  --build-dir build-o1prof \
  --skip-build \
  --baseline /home/polt/compare/profiling/results.txt
```

This run uses `--hess --o1nh -1`, not `--ohess`.

## Direction Generation

The vdW-radius N^(1) neighborhood makes direction generation faster and
reduces the generated direction count.

| Case | Threads | Previous current dirgen | vdW dirgen | wicz dirgen | Previous ndispl | vdW ndispl | wicz ndirs |
|---|---:|---:|---:|---:|---:|---:|---:|
| caffeine | 1 | 0.298 | 0.068 | 0.067 | 62 | 40 | 42 |
| caffeine | 28 | 0.037 | 0.011 | 0.068 | 62 | 40 | 42 |
| taxol | 1 | 12.922 | 7.799 | 5.673 | 131 | 117 | 114 |
| taxol | 28 | 1.118 | 0.621 | 8.650 | 131 | 117 | 114 |

Compared with wicz, this implementation is now essentially tied for
caffeine 1T, still slower for taxol 1T, and faster at 28T because xtb's
direction generator parallelizes.

## Quality Caveat

The current vdW-neighborhood run reports imaginary frequencies:

| Case | Threads | Imaginary frequencies |
|---|---:|---:|
| caffeine | 1 | 1 |
| caffeine | 28 | 1 |
| taxol | 1 | 2 |
| taxol | 28 | 2 |

The earlier non-vdW run had zero imaginary frequencies for these cases.
So the pure direction-generation speedup is real, but this exact port is
not quality-equivalent yet. The likely missing pieces are wicz's matching
`near2` local-solve neighborhood/regularization and progressive negative-mode
handling, rather than direction generation alone.

## Hessian-Only Wall

Because this run uses `--hess`, wall time is no longer inflated by the
pre-Hessian optimizer. Taxol 28T wall is now 20.433 s, faster than the old
xtb `--ohess` baseline wall of 26.923 s and wicz baseline wall of 23.360 s,
but this wall comparison is not perfectly apples-to-apples because the
stored baseline was collected with `--ohess`.
