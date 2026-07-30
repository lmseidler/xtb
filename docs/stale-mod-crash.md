# Stale `.mod` crash after branch switch

## Symptom

`xtb` crashes with SIGSEGV or SIGABRT (`free(): invalid pointer`) at runtime,
typically in `singlepoint` or `hessian`. The crash only appears in optimized
builds; debug and ASan builds may run fine. Backtrace shows corrupted
arguments (e.g. `list=<error: Cannot access memory at 0x39>`).

## Cause

Switching branches with `git checkout` changes `.f90` source files but leaves
compiled `.mod` files in the build directory. If module layouts differ between
branches (added/removed fields, changed types), the stale `.mod` files create
an ABI mismatch: object files compiled against the old modules are linked with
new ones, corrupting argument passing at call sites.

## Fix

Wipe the build directory and reconfigure from scratch:

```bash
rm -rf build
meson setup build --wipe
ninja -C build
```

Or simply reconfigure (meson detects stale mods in some cases):

```bash
meson setup build --reconfigure
ninja -C build
```

## Prevention

After any `git checkout` that touches Fortran source, do a clean rebuild:

```bash
ninja -C build -t clean && ninja -C build
```

Or keep separate build directories per branch.
