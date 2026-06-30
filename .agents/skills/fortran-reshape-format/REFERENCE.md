# Reference: Fortran reshape() Formatter

## Python implementation

```python
def fmt_val(v):
    """Format float as Fortran ES literal, e.g. 1.234567890123456E-05."""
    s = f"{v:.15E}"
    mant, exp = s.split('E')
    sign = exp[0]
    exp_num = int(exp[1:])
    return f"{mant}E{sign}{exp_num:02d}_wp"


def cols_to_fortran_reshape(rows, var_name, shape_str):
    """Convert column-major dump to Fortran reshape literal."""
    flat = []
    for col in rows:
        flat.extend(col)
    n = len(flat)
    header = f"   real(wp), parameter :: {var_name}({shape_str}) = reshape([&"
    footer = f"      & shape({var_name}))"
    lines_out = [header]
    for i in range(0, n, 3):
        is_last = (i == n - 3)
        entries = []
        for j in range(3):
            idx = i + j
            raw = fmt_val(flat[idx])
            if j == 0:
                if raw[0] != '-':
                    entries.append(' ' + raw)
                else:
                    entries.append(raw)
            else:
                entries.append(raw)
        line = "      &" + entries[0]
        for j in range(1, 3):
            if entries[j][0] == '-':
                line += ", " + entries[j]
            else:
                line += ",  " + entries[j]
        if is_last:
            line += "],&"
        else:
            line += ", &"
        lines_out.append(line)
    lines_out.append(footer)
    return "\n".join(lines_out)
```

## Column-major note

Fortran `reshape(flat, shape)` fills the output in column-major order.
If you dump a matrix column-by-column (each `print` statement = one column),
the flat list is already in the correct order for `reshape`.

## Value format

- Mantissa: 15 significant digits after the decimal point
- Exponent: always 2 digits with sign (E+02, E-05, E+00)
- Suffix: `_wp` (working precision, typically `real(8)`)