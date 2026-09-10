#!/usr/bin/env python3
"""Compare two FFT benchmark CSV files and sort by relative performance.

The script expects both CSVs to contain these columns:
  size,batch_size,real_type,time_us,err_percent

It matches rows by (size, batch_size, real_type), computes:
  B_over_A_time = B_time_us / A_time_us
  delta_pct = (B_time_us - A_time_us) / A_time_us * 100

Then it sorts by descending B_over_A_time so worse B performance appears first.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


KEY_COLUMNS = ("size", "batch_size", "real_type")
REQUIRED_COLUMNS = (*KEY_COLUMNS, "time_us")


@dataclass(frozen=True)
class ComparisonRow:
    size: int
    batch_size: int
    real_type: str
    a_time_us: float
    b_time_us: float
    b_over_a_time: float
    delta_pct: float


def _read_csv(path: Path) -> tuple[list[str], dict[tuple[str, str, str], dict[str, str]]]:
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")

        fieldnames = list(reader.fieldnames)
        missing = [name for name in REQUIRED_COLUMNS if name not in fieldnames]
        if missing:
            raise ValueError(
                f"CSV missing required columns {missing} in {path}. "
                f"Found: {fieldnames}"
            )

        rows: dict[tuple[str, str, str], dict[str, str]] = {}
        for row in reader:
            key = (row["size"], row["batch_size"], row["real_type"])
            rows[key] = row

        return fieldnames, rows


def _to_float(value: str, *, field_name: str, path: Path, key: tuple[str, str, str]) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(
            f"Cannot parse {field_name}={value!r} in {path} for key={key}"
        ) from exc


def _format_table(rows: list[ComparisonRow], limit: int | None) -> str:
    head = [
        "| size | batch_size | real_type | A time (us) | B time (us) | B/A time | delta % |",
        "|---:|---:|:---|---:|---:|---:|---:|",
    ]

    body: list[str] = []
    subset = rows if limit is None else rows[:limit]
    for row in subset:
        body.append(
            f"| {row.size} | {row.batch_size} | {row.real_type} | "
            f"{row.a_time_us:.6f} | {row.b_time_us:.6f} | "
            f"{row.b_over_a_time:.6f} | {row.delta_pct:.2f} |"
        )

    return "\n".join(head + body)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", required=True, help="Path to baseline CSV (A)")
    parser.add_argument("--b", required=True, help="Path to comparison CSV (B)")
    parser.add_argument(
        "--output-csv",
        default="benchmark/fft/bench_fft_B_vs_A_sorted.csv",
        help="Output CSV path",
    )
    parser.add_argument(
        "--output-md",
        default="benchmark/fft/bench_fft_B_vs_A_sorted.md",
        help="Output Markdown table path",
    )
    parser.add_argument(
        "--print-top",
        type=int,
        default=20,
        help="Print top N rows to stdout (set <=0 to print all)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if keys do not match exactly between A and B",
    )
    args = parser.parse_args()

    a_path = Path(args.a)
    b_path = Path(args.b)
    out_csv = Path(args.output_csv)
    out_md = Path(args.output_md)

    _, a_rows = _read_csv(a_path)
    _, b_rows = _read_csv(b_path)

    keys_a = set(a_rows.keys())
    keys_b = set(b_rows.keys())

    only_a = keys_a - keys_b
    only_b = keys_b - keys_a

    if args.strict and (only_a or only_b):
        raise ValueError(
            "Key mismatch between A and B. "
            f"only in A: {len(only_a)}, only in B: {len(only_b)}"
        )

    shared_keys = sorted(keys_a & keys_b, key=lambda k: (int(k[0]), int(k[1]), k[2]))

    compared: list[ComparisonRow] = []
    for key in shared_keys:
        a_time_us = _to_float(a_rows[key]["time_us"], field_name="time_us", path=a_path, key=key)
        b_time_us = _to_float(b_rows[key]["time_us"], field_name="time_us", path=b_path, key=key)

        if a_time_us == 0.0:
            b_over_a = float("inf")
            delta_pct = float("inf")
        else:
            b_over_a = b_time_us / a_time_us
            delta_pct = (b_time_us - a_time_us) / a_time_us * 100.0

        compared.append(
            ComparisonRow(
                size=int(key[0]),
                batch_size=int(key[1]),
                real_type=key[2],
                a_time_us=a_time_us,
                b_time_us=b_time_us,
                b_over_a_time=b_over_a,
                delta_pct=delta_pct,
            )
        )

    # Worse B performance at top (larger B/A means slower B)
    compared.sort(key=lambda row: row.b_over_a_time, reverse=True)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "size",
                "batch_size",
                "real_type",
                "A_time_us",
                "B_time_us",
                "B_over_A_time",
                "delta_pct",
            ]
        )
        for row in compared:
            writer.writerow(
                [
                    row.size,
                    row.batch_size,
                    row.real_type,
                    f"{row.a_time_us:.6f}",
                    f"{row.b_time_us:.6f}",
                    f"{row.b_over_a_time:.6f}",
                    f"{row.delta_pct:.2f}",
                ]
            )

    out_md.parent.mkdir(parents=True, exist_ok=True)
    limit = None if args.print_top <= 0 else args.print_top
    md = _format_table(compared, limit=None)
    out_md.write_text(md + "\n")

    print(f"A rows: {len(a_rows)}")
    print(f"B rows: {len(b_rows)}")
    print(f"Shared rows compared: {len(compared)}")
    print(f"Rows only in A: {len(only_a)}")
    print(f"Rows only in B: {len(only_b)}")
    print(f"Wrote CSV: {out_csv}")
    print(f"Wrote Markdown table: {out_md}")

    print()
    print("Top rows (B slower than A at top):")
    print(_format_table(compared, limit=limit))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
