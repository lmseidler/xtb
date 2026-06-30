#!/usr/bin/env python3
"""Run and summarize xtb O1NumHess profiling cases."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


PROF_RE = re.compile(r"PROF (odlrhessian|gen_displdir): (.*)")
NATOMS_RE = re.compile(r"natoms=(\d+),\s*n(?:displ|dirs)=(\d+)")
PHASE_RE = re.compile(r"([A-Za-z_+]+)\s*=\s*([-+0-9.]+)\s*s")
GEN_WALL_RE = re.compile(
    r"wall=\s*([-+0-9.]+)\s*s,\s*evals=(\d+),\s*avg_nnb=\s*([-+0-9.]+),\s*max_nnb=(\d+)"
)
GEN_BREAKDOWN_RE = re.compile(r"([A-Za-z_]+)=\s*([-+0-9.]+)\s*s")
IMAG_RE = re.compile(r"# imaginary freq\.\s+(\d+)")
BASE_HEADER_RE = re.compile(r"^(\S+)\s+(xtb|wicz)\s+(\d+)T:\s+wall=([-+0-9.]+)s")

PLAN_TARGETS = {
    ("taxol", 1): {"phase12_wall": 60.0, "all_wall": 55.0},
    ("taxol", 28): {"phase12_wall": 22.0, "all_wall": 10.0},
}


def run_command(cmd: list[str], cwd: Path, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def prepend_path(value: str | None, prefix: Path) -> str:
    return f"{prefix}{os.pathsep}{value}" if value else str(prefix)


def build_repo(repo: Path, build_dir: str, skip_build: bool) -> str:
    if skip_build:
        return "build skipped"

    conda = shutil.which("conda")
    if conda:
        cmd = [conda, "run", "-n", "gnu-15", "meson", "compile", "-C", build_dir]
    else:
        cmd = ["meson", "compile", "-C", build_dir]

    proc = run_command(cmd, repo)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        raise SystemExit(f"Build failed with exit code {proc.returncode}: {' '.join(cmd)}")
    return proc.stdout


def parse_profile(text: str) -> dict:
    data: dict = {"phases": {}, "gen_displdir": {}}
    for line in text.splitlines():
        prof = PROF_RE.search(line)
        if prof:
            section, payload = prof.groups()
            if section == "odlrhessian":
                natoms = NATOMS_RE.search(payload)
                if natoms:
                    data["natoms"] = int(natoms.group(1))
                    data["ndispl"] = int(natoms.group(2))
                    continue
                phase = PHASE_RE.search(payload)
                if phase:
                    data["phases"][phase.group(1)] = float(phase.group(2))
                    continue
            elif section == "gen_displdir":
                gen = GEN_WALL_RE.search(payload)
                if gen:
                    data["gen_displdir"].update(
                        {
                            "wall": float(gen.group(1)),
                            "evals": int(gen.group(2)),
                            "avg_nnb": float(gen.group(3)),
                            "max_nnb": int(gen.group(4)),
                        }
                    )
                    continue
                for key, value in GEN_BREAKDOWN_RE.findall(payload):
                    data["gen_displdir"][key] = float(value)

        imag = IMAG_RE.search(line)
        if imag:
            data["imag_freqs"] = int(imag.group(1))

    return data


def parse_baseline(path: Path | None) -> dict:
    if not path or not path.exists():
        return {}

    baseline: dict = {}
    current_key: tuple[str, str, int] | None = None
    for line in path.read_text(errors="replace").splitlines():
        header = BASE_HEADER_RE.match(line)
        if header:
            case, impl, threads, wall = header.groups()
            current_key = (case, impl, int(threads))
            baseline[current_key] = {"wall": float(wall), "phases": {}}
            continue

        if current_key is None:
            continue

        natoms = NATOMS_RE.search(line)
        if natoms:
            baseline[current_key]["natoms"] = int(natoms.group(1))
            baseline[current_key]["ndispl"] = int(natoms.group(2))
            continue

        phase = PHASE_RE.search(line)
        if phase:
            baseline[current_key]["phases"][phase.group(1)] = float(phase.group(2))

    return baseline


def fmt(value: float | int | None, suffix: str = "") -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return f"{value}{suffix}"
    return f"{value:.3f}{suffix}"


def pct_change(new: float | None, old: float | None) -> str:
    if new is None or old in (None, 0):
        return "-"
    return f"{(new - old) / old * 100.0:+.1f}%"


def make_markdown(results: list[dict], baseline: dict) -> str:
    lines = [
        "# O1NumHess Profiling Summary",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Current Runs",
        "",
        "| Case | Threads | Wall s | TOTAL s | dirgen s | grad_evals s | ndispl | imag |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]

    for result in results:
        phases = result.get("phases", {})
        lines.append(
            "| {case} | {threads} | {wall} | {total} | {dirgen} | {grad} | {ndispl} | {imag} |".format(
                case=result["case"],
                threads=result["threads"],
                wall=fmt(result.get("wall")),
                total=fmt(phases.get("TOTAL")),
                dirgen=fmt(phases.get("dirgen")),
                grad=fmt(phases.get("grad_evals")),
                ndispl=fmt(result.get("ndispl")),
                imag=fmt(result.get("imag_freqs")),
            )
        )

    if baseline:
        lines += [
            "",
            "## Baseline Comparison",
            "",
            "| Case | Threads | New wall | Old xtb wall | vs old | wicz wall | vs wicz | New TOTAL | Old xtb TOTAL | wicz TOTAL |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for result in results:
            case = result["case"]
            threads = result["threads"]
            old = baseline.get((case, "xtb", threads), {})
            wicz = baseline.get((case, "wicz", threads), {})
            new_total = result.get("phases", {}).get("TOTAL")
            old_total = old.get("phases", {}).get("TOTAL")
            wicz_total = wicz.get("phases", {}).get("TOTAL")
            lines.append(
                "| {case} | {threads} | {new_wall} | {old_wall} | {old_delta} | {wicz_wall} | {wicz_delta} | {new_total} | {old_total} | {wicz_total} |".format(
                    case=case,
                    threads=threads,
                    new_wall=fmt(result.get("wall")),
                    old_wall=fmt(old.get("wall")),
                    old_delta=pct_change(result.get("wall"), old.get("wall")),
                    wicz_wall=fmt(wicz.get("wall")),
                    wicz_delta=pct_change(result.get("wall"), wicz.get("wall")),
                    new_total=fmt(new_total),
                    old_total=fmt(old_total),
                    wicz_total=fmt(wicz_total),
                )
            )

    lines += [
        "",
        "## Taxol Plan Targets",
        "",
        "| Threads | New wall | Phase 1+2 target | vs Phase 1+2 | All-plan target | vs all-plan |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for result in results:
        key = (result["case"], result["threads"])
        if key not in PLAN_TARGETS:
            continue
        targets = PLAN_TARGETS[key]
        wall = result.get("wall")
        lines.append(
            "| {threads} | {wall} | {phase12} | {phase12_delta} | {all_target} | {all_delta} |".format(
                threads=result["threads"],
                wall=fmt(wall),
                phase12=fmt(targets["phase12_wall"]),
                phase12_delta=pct_change(wall, targets["phase12_wall"]),
                all_target=fmt(targets["all_wall"]),
                all_delta=pct_change(wall, targets["all_wall"]),
            )
        )

    gen_rows = [r for r in results if r.get("gen_displdir")]
    if gen_rows:
        lines += [
            "",
            "## Direction Generation Breakdown",
            "",
            "| Case | Threads | wall | evals | avg_nnb | max_nnb | extract | orth | proj | diag | sign |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for result in gen_rows:
            gen = result["gen_displdir"]
            lines.append(
                "| {case} | {threads} | {wall} | {evals} | {avg_nnb} | {max_nnb} | {extract} | {orth} | {proj} | {diag} | {sign} |".format(
                    case=result["case"],
                    threads=result["threads"],
                    wall=fmt(gen.get("wall")),
                    evals=fmt(gen.get("evals")),
                    avg_nnb=fmt(gen.get("avg_nnb")),
                    max_nnb=fmt(gen.get("max_nnb")),
                    extract=fmt(gen.get("extract")),
                    orth=fmt(gen.get("orth")),
                    proj=fmt(gen.get("proj")),
                    diag=fmt(gen.get("diag")),
                    sign=fmt(gen.get("sign")),
                )
            )

    lines.append("")
    return "\n".join(lines)


def default_cases(repo: Path, taxol_geom: Path | None = None) -> list[tuple[str, Path]]:
    cases = [
        ("caffeine", repo / "subprojects/dftd4/app/04-caffeine.xyz"),
        ("taxol", taxol_geom if taxol_geom else repo / "assets/inputs/xyz/taxol.xyz"),
    ]
    return [(name, path) for name, path in cases if path is not None]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--build-dir", default="build")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--threads", default="1,28")
    parser.add_argument("--taxol-geom", type=Path, default=None,
                        help="Path to optimized taxol geometry (e.g. compare/profiling/xtbopt.xyz). "
                             "Default: assets/inputs/xyz/taxol.xyz")
    parser.add_argument("--skip-build", action="store_true")
    args = parser.parse_args()

    repo = args.repo.resolve()
    threads = [int(item) for item in args.threads.split(",") if item.strip()]
    out_dir = args.out_dir or repo / "docs/profiling" / f"o1numhess-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    out_dir.mkdir(parents=True, exist_ok=True)

    build_log = build_repo(repo, args.build_dir, args.skip_build)
    (out_dir / "build.log").write_text(build_log)

    xtb = repo / args.build_dir / "xtb"
    if not xtb.exists():
        raise SystemExit(f"Missing built xtb binary: {xtb}")

    conda_prefix = Path.home() / ".conda/envs/gnu-15"
    base_env = os.environ.copy()
    base_env["PATH"] = prepend_path(prepend_path(base_env.get("PATH"), conda_prefix / "bin"), xtb.parent)

    results = []
    for case, molecule in default_cases(repo, args.taxol_geom):
        if not molecule.exists():
            raise SystemExit(f"Missing molecule for {case}: {molecule}")

        for nth in threads:
            run_dir = out_dir / f"{case}-{nth}T"
            run_dir.mkdir(parents=True, exist_ok=True)
            local_mol = run_dir / molecule.name
            shutil.copyfile(molecule, local_mol)

            env = base_env.copy()
            env["OMP_NUM_THREADS"] = str(nth)
            env.setdefault("OMP_PROC_BIND", "spread")
            env.setdefault("OMP_PLACES", "sockets")

            cmd = [str(xtb), local_mol.name, "--hess", "--o1nh"]
            start = time.perf_counter()
            proc = run_command(cmd, run_dir, env)
            wall = time.perf_counter() - start

            log = proc.stdout
            (run_dir / "xtb.log").write_text(log)
            parsed = parse_profile(log)
            parsed.update({"case": case, "threads": nth, "wall": wall, "returncode": proc.returncode})
            results.append(parsed)

            if proc.returncode != 0:
                raise SystemExit(f"xtb failed for {case} {nth}T with exit code {proc.returncode}; see {run_dir / 'xtb.log'}")

    baseline = parse_baseline(args.baseline)
    summary = {
        "repo": str(repo),
        "build_dir": args.build_dir,
        "out_dir": str(out_dir),
        "threads": threads,
        "results": results,
        "baseline": {"/".join(map(str, key)): value for key, value in baseline.items()},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))
    markdown = make_markdown(results, baseline)
    (out_dir / "summary.md").write_text(markdown)
    print(markdown)
    print(f"Output directory: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
