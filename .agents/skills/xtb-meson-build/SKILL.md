---
name: xtb-meson-build
description: Build, test, and profile xtb worktrees with Meson using the conda environment gnu-15. Use when Codex needs to configure a fresh xtb build directory, compile xtb, run Meson tests, recover from stale Fortran .mod/object issues after branch or interface changes, avoid dependency clone failures by reusing populated subprojects, or run O1NumHess profiling in this repository.
---

# xTB Meson Build

## Workflow

Use conda env `gnu-15` and Meson/Ninja for compilation and tests. Prefer a separate `build` directory per worktree/branch. For Fortran module interface changes, do a clean rebuild before trusting runtime or unit-test failures.

## Configure

If `build` is missing or Meson must be reset, run setup with explicit compiler and search-path overrides. This avoids inherited nested-conda variables that can point at `/opt/miniforge/bin/x86_64-conda-linux-gnu-gfortran`, which may not exist even though `gnu-15` contains a working compiler.

```bash
conda run -n gnu-15 env \
  FC=/home/polt/.conda/envs/gnu-15/bin/gfortran \
  F90=/home/polt/.conda/envs/gnu-15/bin/gfortran \
  F77=/home/polt/.conda/envs/gnu-15/bin/gfortran \
  CC=/home/polt/.conda/envs/gnu-15/bin/x86_64-conda-linux-gnu-cc \
  CXX=/home/polt/.conda/envs/gnu-15/bin/x86_64-conda-linux-gnu-c++ \
  PKG_CONFIG_PATH=/home/polt/.conda/envs/gnu-15/lib/pkgconfig \
  CMAKE_PREFIX_PATH=/home/polt/.conda/envs/gnu-15 \
  meson setup build --wipe
```

If setup tries to clone wrap subprojects and network or SSH config fails, check whether `/home/polt/xtb/subprojects` already has populated dependency directories. For temporary worktrees, copying those populated subprojects into the worktree is acceptable:

```bash
cp -a /home/polt/xtb/subprojects/. /path/to/worktree/subprojects/
```

## Build And Test

Compile:

```bash
conda run -n gnu-15 meson compile -C build
```

Run the focused Hessian test:

```bash
conda run -n gnu-15 meson test -C build hessian --print-errorlogs
```

## Stale Module Failures

After branch switches or public Fortran interface/type changes, stale `.mod` files and objects can cause false failures or runtime crashes. Symptoms include SIGSEGV/SIGABRT, `free(): invalid pointer`, corrupted arguments, or a unit failure that disappears after a clean rebuild.

Before debugging such failures, clean and rebuild:

```bash
conda run -n gnu-15 ninja -C build -t clean
conda run -n gnu-15 meson compile -C build
conda run -n gnu-15 meson test -C build hessian --print-errorlogs
```

If setup itself is suspect, wipe/reconfigure instead of only cleaning.

## O1NumHess Profiling

Use the repository profiling helper when available:

```bash
conda run -n gnu-15 python /home/polt/xtb/skills/o1numhess-profiling/scripts/profile_o1numhess.py \
  --repo /path/to/worktree \
  --baseline /home/polt/compare/profiling/results.txt
```

Read `summary.md` plus `xtb.log` files under the generated `docs/profiling/o1numhess-*` directory. Compare `dirgen`, `ndispl`, `imag`, `avg_nnb`, and `max_nnb`; speedups that increase `ndispl` or imaginary-mode reruns may regress total wall time.
