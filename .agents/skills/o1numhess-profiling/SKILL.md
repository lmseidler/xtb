---
name: o1numhess-profiling
description: Reproducible profiling workflow for xtb O1NumHess changes. Use when comparing current `TCalculator%odlrhessian()` performance against the stored old xtb implementation, xtb_wicz baseline, or PLAN.md targets; rerunning caffeine/taxol timing matrices; collecting `PROF odlrhessian` / `PROF gen_displdir` output; or summarizing O1NumHess wall-time, phase-time, direction-count, and imaginary-frequency results.
---

# O1NumHess Profiling

## Workflow

Use `scripts/profile_o1numhess.py` from the repository root. It builds with
Meson, runs the instrumented `build/xtb` binary on caffeine and taxol at 1
and 28 OpenMP threads, stores full logs, and writes Markdown/JSON summaries.
The runner uses `--hess --o1nh` so the profiling isolates the O1NumHess
Hessian path. Do not use `--ohess` for this workflow because it includes the
pre-Hessian geometry optimization and makes wall-time comparisons misleading.
Do not carry over the `-1` optimization-level argument from `--ohess`; `--hess`
does not take that argument.

For taxol, pass `--taxol-geom` pointing to the optimized geometry
(`compare/profiling/xtbopt.xyz`, gnorm ≈ 0.0025) so the Hessian is evaluated
at a stationary point. The wicz baseline in `results.txt` was generated on
this same optimized geometry; using the unoptimized `taxol.xyz` (gnorm ≈ 0.14)
produces spurious imaginary frequencies from the non-zero gradient.

Default command:

```bash
python skills/o1numhess-profiling/scripts/profile_o1numhess.py \
  --repo /home/polt/xtb \
  --baseline /home/polt/compare/profiling/results.txt \
  --taxol-geom /home/polt/compare/profiling/xtbopt.xyz
```

The script follows the local repo instruction to compile with the `gnu-15`
conda environment. It uses `conda run -n gnu-15 meson compile -C build` for
the build step, then runs `build/xtb` directly with the build directory and
conda environment prepended to `PATH`.

## Outputs

Expect an output directory under `docs/profiling/o1numhess-<timestamp>/` with:

- `summary.md`: current timings plus old xtb and xtb_wicz comparisons when a
  baseline file is available.
- `summary.json`: machine-readable parsed timings.
- `*/xtb.log`: full stdout/stderr for each case/thread run.
- copied input molecules for auditability.

## Interpretation

Compare `TOTAL` and external `wall` separately. `TOTAL` comes from the
instrumented O1NumHess phases and excludes some process/output overhead;
`wall` is the full command runtime measured by the runner.

For taxol, the relevant PLAN targets are:

- Phase 1+2: about 60 s wall at 1 thread and about 22 s wall at 28 threads.
- All planned improvements: about 55 s wall at 1 thread and about 10 s wall
  at 28 threads.

The old baseline file uses labels like `taxol xtb 1T` and
`taxol wicz 28T`. Treat `xtb` there as the old implementation before the
current worktree changes, and `wicz` as the xtb_wicz comparison point.
