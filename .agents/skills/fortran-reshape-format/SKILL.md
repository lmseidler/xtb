---
name: fortran-reshape-format
description: Format Fortran reshape() array literals with aligned columns matching xtb codebase style. Use when generating or updating parameter arrays (hessian_ref, dipgrad_ref, density) in test files or source.
user-invocable: false
---

# Fortran reshape() Array Formatter

## Quick start

```python
from fortran_reshape_format import cols_to_fortran_reshape

# rows[i] = i-th column (column-major dump)
block = cols_to_fortran_reshape(rows, "hessian_ref", "9, 9")
# -> "   real(wp), parameter :: hessian_ref(9, 9) = reshape([&\n      & ..."
```

## Why

Fortran `reshape()` takes a column-major flat array. Test reference arrays must be both readable and diff-friendly. This skill ensures consistent alignment: 3 values per line, signed-vs-unsigned padding, and trailing `&` alignment.

## Format rules

1. Header: `   real(wp), parameter :: VAR(SHAPE) = reshape([&` (no space before `&`)
2. 3 values per data line
3. Each first-on-line entry padded to 20 chars: leading space if positive, raw `-XX...` if negative
4. Subsequent entries keep natural width (19 chars positive, 20 negative)
5. Separator between entries: `, ` (2 chars) before negative, `,  ` (3 chars) before positive — keeps columns aligned
6. Line ending: `, &` for regular lines, `],&` for last data line
7. Footer: `      & shape(VAR))`
8. Value format: `ES22.15E2` — 15 significant digits, 2-digit exponent (e.g. `1.234567890123456E-05`)

## Workflow

1. Collect values as column-major flat list (column 0 rows first, then column 1, etc.)
2. Format each value to 15-digit mantissa + 2-digit exponent
3. Group into triplets, apply padding rules
4. Emit header, data lines, footer

## Anti-pattern

WRONG: one value per line, inconsistent indentation
```
   & 1.0E-05_wp,
   & 2.0E-05_wp,
```

RIGHT: 3 per line, aligned
```
      & 1.000000000000000E-05_wp,  2.000000000000000E-05_wp,  3.000000000000000E-05_wp, &
```

## Reference

See `extract_density.f90` in repo root for the canonical format generator.

## Checklist

- 3 values per line
- Positive first-entry has leading space (20-char width)
- Negative entries use `, ` separator; positive use `,  `
- Last data line ends with `],&` not `, &`
- Footer has `      & ` (6-space indent) before `shape(VAR))`