# O1NumHess Profiling Notes

Run date: 2026-06-25

Command:

```bash
python skills/o1numhess-profiling/scripts/profile_o1numhess.py \
  --repo /home/polt/xtb \
  --baseline /home/polt/compare/profiling/results.txt
```

## Headline

The current worktree is mixed:

- Caffeine at 28 threads is now slightly faster than xtb_wicz wall time
  (0.462 s vs 0.501 s).
- Taxol O1NumHess `TOTAL` at 28 threads improved versus old xtb
  (16.970 s vs 19.116 s), but full command wall time regressed
  (50.167 s vs 26.923 s) because the pre-Hessian ANC optimizer took
  31.444 s.
- Taxol at 1 thread regressed hard in O1NumHess itself:
  `TOTAL` 130.862 s vs old xtb 72.708 s and xtb_wicz 50.259 s.

## Taxol Phase Comparison

| Phase | Old xtb 1T | Current 1T | Delta | Old xtb 28T | Current 28T | Delta |
|---|---:|---:|---:|---:|---:|---:|
| dirgen | 23.267 | 12.922 | -44.5% | 3.959 | 1.118 | -71.8% |
| grad_evals | 46.515 | 109.627 | +135.7% | 14.478 | 10.167 | -29.8% |
| local_solve | 0.288 | 1.972 | +584.7% | 0.163 | 2.068 | +1168.7% |
| neg_mode | 2.322 | 5.778 | +148.8% | 0.400 | 3.243 | +710.8% |
| TOTAL | 72.708 | 130.862 | +80.0% | 19.116 | 16.970 | -11.2% |
| wall | 87.180 | 173.309 | +98.8% | 26.923 | 50.167 | +86.3% |

## Plan Comparison

Taxol wall targets in `PLAN.md`:

| Threads | Current wall | Phase 1+2 target | All-plan target |
|---:|---:|---:|---:|
| 1 | 173.309 | ~60 | ~55 |
| 28 | 50.167 | ~22 | ~10 |

Current taxol wall misses both targets. For 1 thread, the main blocker is
`grad_evals` increasing from 46.515 s to 109.627 s. For 28 threads, the
instrumented O1NumHess `TOTAL` is good relative to old xtb, but the driver
wall time is dominated by the optimizer before the Hessian.

## Likely Next Checks

- Investigate why current taxol 1-thread gradient evaluations are much
  slower despite fewer displacement directions (131 vs old 135).
- Re-run a Hessian-only path from already optimized geometry if the goal is
  isolating O1NumHess from the `--ohess` pre-optimization cost.
- Check whether current changes to `get_gradient_derivs` or nested OpenMP
  settings changed single-thread SCF behavior.
- Revisit `local_solve` and `neg_mode`: both are now much slower than old xtb
  in the taxol run, especially at 28 threads.
