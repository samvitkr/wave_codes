#!/usr/bin/env python3
"""Run VKFFT ddx benchmark sweep and export results as CSV.

This script runs:
  pixi run -e cuda build-local-cuda-single/benchmark/fft/bench_fft \
    --nx <size> --ny <batch> --nz 1 --ops ddx --backend vkfft

for all configured transform sizes and batch sizes, parses the benchmark table,
and writes a CSV with:
  size, batch_size, real_type, time_us, err_percent

Time values are normalized to microseconds in the CSV.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import re
import subprocess
import sys
from typing import Iterable


TRANSFORM_SIZES: list[int] = [
    128,
    160,
    144,
    192,
    224,
    256,
    288,
    320,
    384,
    420,
    512,
    576,
    640,
    768,
    896,
    1024,
    1152,
    1280,
    1536,
    1792,
    2048,
    4096,
    5120,
]

BATCH_SIZES: list[int] = [
    1,
    2,
    3,
    6,
    7,
    9,
    12,
    128,
    160,
    144,
    192,
    224,
    256,
    288,
    320,
    384,
    420,
    512,
    576,
    640,
    768,
    896,
    1024,
]

NUMBER_RE = re.compile(r"[-+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?")


def _parse_number(cell: str) -> float:
    match = NUMBER_RE.search(cell)
    if match is None:
        raise ValueError(f"Could not parse numeric value from cell: {cell!r}")
    return float(match.group(0).replace(",", ""))


def _extract_real_type(name_cell: str) -> str:
    # Name cell is usually markdown code, e.g. `double Cuda VkFFT ddx`
    match = re.search(r"`([^`]+)`", name_cell)
    if match is None:
        return ""

    name = match.group(1).strip().lower()
    if name.startswith("float ") or name == "float":
        return "float"
    if name.startswith("double ") or name == "double":
        return "double"
    return ""


def _time_unit_to_us_multiplier(unit_cell: str) -> float:
    unit = unit_cell.strip().lower()
    if unit.endswith("/op"):
        unit = unit[:-3]
    unit = unit.strip()

    # Support ASCII and unicode spellings.
    if unit in {"s", "sec", "second", "seconds"}:
        return 1_000_000.0
    if unit in {"ms"}:
        return 1_000.0
    if unit in {"us", "µs", "μs"}:
        return 1.0
    if unit in {"ns"}:
        return 0.001

    raise ValueError(f"Unsupported time unit in benchmark header: {unit_cell!r}")


def _parse_benchmark_output(text: str) -> list[tuple[str, float, float]]:
    """Parse bench markdown table rows.

    Returns list of tuples: (real_type, time_us, err_percent).
    """
    idx_time: int | None = None
    idx_err: int | None = None
    time_to_us: float | None = None
    parsed: list[tuple[str, float, float]] = []

    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        if not line.startswith("|"):
            continue

        parts = line.split("|")
        if len(parts) < 3:
            continue
        # Bench markdown rows usually do not end with a trailing '|'.
        # Keep all cells after the first delimiter, but drop a final empty token
        # when a trailing delimiter is present.
        raw_cells = parts[1:]
        if raw_cells and raw_cells[-1].strip() == "":
            raw_cells = raw_cells[:-1]
        cells = [cell.strip() for cell in raw_cells]
        if not cells:
            continue

        # Header line example: | us/op | op/s | err% | ... |
        if any("err%" == cell for cell in cells):
            idx_err = next((i for i, cell in enumerate(cells) if cell == "err%"), None)
            idx_time = next(
                (i for i, cell in enumerate(cells) if cell.endswith("/op")), None
            )
            if idx_time is not None:
                time_to_us = _time_unit_to_us_multiplier(cells[idx_time])
            continue

        if idx_time is None or idx_err is None or time_to_us is None:
            continue

        if len(cells) <= max(idx_time, idx_err):
            continue

        if "`" not in line:
            continue

        try:
            time_per_op = _parse_number(cells[idx_time])
            err_percent = _parse_number(cells[idx_err])
        except ValueError:
            continue

        real_type = _extract_real_type(cells[-1])
        if real_type not in {"float", "double"}:
            continue

        parsed.append((real_type, time_per_op * time_to_us, err_percent))

    return parsed


def _extract_suggested_min_epoch_iterations(text: str) -> int | None:
    # Example unstable hint in benchmark output:
    # "Increase `minEpochIterations` to e.g. 12345"
    match = re.search(r"Increase\s+`minEpochIterations`\s+to\s+e\.g\.\s+([\d,]+)", text)
    if match is None:
        return None
    return int(match.group(1).replace(",", ""))


def _parse_csv_values(raw: str) -> list[str]:
    return [token.strip() for token in raw.split(",") if token.strip()]


def _run_single_case(
    exe: str,
    nx: int,
    ny: int,
    nz: int,
    ops: str,
    backend: str,
    real_types: str,
    err_threshold_percent: float,
    max_reruns: int,
) -> list[tuple[str, float, float]]:
    if "," in ops or "," in backend or "," in real_types:
        raise ValueError(
            "_run_single_case() requires single values for ops/backend/real_types"
        )

    cmd = [
        "pixi",
        "run",
        "-e",
        "cuda",
        exe,
        "--nx",
        str(nx),
        "--ny",
        str(ny),
        "--nz",
        str(nz),
        "--ops",
        ops,
        "--backend",
        backend,
        "--real-type",
        real_types,
    ]

    attempts = 1 + max_reruns
    last_parsed: list[tuple[str, float, float]] = []
    min_epoch_iterations: int | None = None
    for attempt in range(1, attempts + 1):
        cmd_attempt = list(cmd)
        if min_epoch_iterations is not None:
            cmd_attempt.extend(
                [
                    "--min-epoch-iterations",
                    str(min_epoch_iterations),
                ]
            )

        completed = subprocess.run(
            cmd_attempt,
            check=False,
            capture_output=True,
            text=True,
        )

        combined_output = f"{completed.stdout}\n{completed.stderr}"

        if completed.returncode != 0:
            raise RuntimeError(
                "Benchmark command failed for "
                f"nx={nx}, ny={ny}, returncode={completed.returncode}\n"
                f"Command: {' '.join(cmd_attempt)}\n"
                f"Output:\n{combined_output}"
            )

        parsed = _parse_benchmark_output(combined_output)
        if not parsed:
            raise RuntimeError(
                "Could not parse benchmark output for "
                f"nx={nx}, ny={ny}.\nCommand: {' '.join(cmd_attempt)}\n"
                f"Output:\n{combined_output}"
            )

        last_parsed = parsed
        max_err = max(err_percent for _, _, err_percent in parsed)
        if max_err <= err_threshold_percent:
            return parsed

        if attempt < attempts:
            suggested_iters = _extract_suggested_min_epoch_iterations(combined_output)
            if suggested_iters is not None:
                if min_epoch_iterations is None:
                    min_epoch_iterations = suggested_iters
                else:
                    min_epoch_iterations = max(min_epoch_iterations, suggested_iters)

            print(
                "  err% exceeded threshold "
                f"({max_err:.3f} > {err_threshold_percent:.3f}), "
                f"re-running attempt {attempt + 1}/{attempts}"
                + (
                    f" with --min-epoch-iterations {min_epoch_iterations}"
                    if min_epoch_iterations is not None
                    else ""
                ),
                file=sys.stderr,
            )

    print(
        "  warning: err% remained above threshold after retries; "
        f"using last run for nx={nx}, ny={ny}",
        file=sys.stderr,
    )
    return last_parsed


def _iter_cases() -> Iterable[tuple[int, int]]:
    for nx in TRANSFORM_SIZES:
        for ny in BATCH_SIZES:
            yield nx, ny


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="benchmark/fft/bench_fft_vkfft_ddx_sweep.csv",
        help="Output CSV file path",
    )
    parser.add_argument(
        "--exe",
        default="build-local-cuda-single/benchmark/fft/bench_fft",
        help="Benchmark executable path",
    )
    parser.add_argument("--ops", default="ddx", help="Operation(s) for --ops")
    parser.add_argument(
        "--backend",
        default="vkfft",
        help="Backend value for --backend",
    )
    parser.add_argument(
        "--real-type",
        default="float,double",
        help="Value passed to --real-type (e.g., float,double or double)",
    )
    parser.add_argument("--nz", type=int, default=1, help="Value passed to --nz")
    parser.add_argument(
        "--err-threshold",
        type=float,
        default=10.0,
        help="Maximum allowed err%% before re-running a case",
    )
    parser.add_argument(
        "--max-reruns",
        type=int,
        default=3,
        help="Maximum number of re-runs when err%% exceeds threshold",
    )
    args = parser.parse_args()

    ops_values = _parse_csv_values(args.ops)
    backend_values = _parse_csv_values(args.backend)
    real_type_values = _parse_csv_values(args.real_type)
    if not ops_values or not backend_values or not real_type_values:
        raise ValueError("--ops, --backend, and --real-type must not be empty")

    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["size", "batch_size", "real_type", "time_us", "err_percent"])

        total = len(TRANSFORM_SIZES) * len(BATCH_SIZES)
        index = 0
        for nx, ny in _iter_cases():
            index += 1
            print(f"[{index:4d}/{total}] Running nx={nx}, ny={ny}...", file=sys.stderr)
            for op in ops_values:
                for backend in backend_values:
                    for real_type in real_type_values:
                        rows = _run_single_case(
                            exe=args.exe,
                            nx=nx,
                            ny=ny,
                            nz=args.nz,
                            ops=op,
                            backend=backend,
                            real_types=real_type,
                            err_threshold_percent=args.err_threshold,
                            max_reruns=args.max_reruns,
                        )
                        for row_real_type, time_us, err_percent in rows:
                            writer.writerow(
                                [
                                    nx,
                                    ny,
                                    row_real_type,
                                    f"{time_us:.6f}",
                                    f"{err_percent:.6f}",
                                ]
                            )
            f.flush()

    print(f"Wrote CSV results to: {out_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
