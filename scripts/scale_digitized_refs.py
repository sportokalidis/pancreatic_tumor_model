#!/usr/bin/env python3
"""
Scale digitized paper reference files to different ABM scales.

Original files (data-export/*-Cells.csv) are at S=1e5 (paper scale).
Generate scaled versions for S=1e4, 1e3, etc. in data-export/digitized/

Usage:
  python3 scripts/scale_digitized_refs.py \\
    --scales 1e4 1e3 \\
    --source-dir data-export \\
    --out-dir data-export/digitized
"""

import argparse
from pathlib import Path


def scale_reference(input_csv: Path, scale_s: float, output_csv: Path):
    """
    Convert digitized cell counts to agent counts for a given ABM scale.

    Digitized references contain CELL COUNTS from paper figures.
    ABM uses AGENT COUNTS where 1 agent = scale_s cells.

    To convert: agent_count = cell_count / scale_s

    Args:
        input_csv: Original *-Cells.csv (cell counts from paper)
        scale_s: Target ABM scale (e.g., 1e4, 1e3)
        output_csv: Output path (agent counts for ABM)
    """
    # Read original (no header, col0=day, col1=cell_count)
    rows = []
    with open(input_csv) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            day_str, count_str = line.split(',')
            day = float(day_str)
            cell_count = float(count_str)
            rows.append((day, cell_count))

    # Convert cell counts to agent counts: agent_count = cell_count / scale
    # (because 1 agent represents scale_s cells)

    # Write with formatting
    with open(output_csv, "w") as f:
        for day, cell_count in rows:
            agent_count = int(cell_count / scale_s)
            f.write(f"{day:.1f},{agent_count}\n")

    return len(rows)


def main():
    parser = argparse.ArgumentParser(description="Scale digitized paper references")
    parser.add_argument("--scales", nargs="+", type=float, default=[1e4, 1e3],
                        help="Scale factors to generate (default: 1e4 1e3)")
    parser.add_argument("--source-dir", default="data-export",
                        help="Directory with original *-Cells.csv files")
    parser.add_argument("--out-dir", default="data-export/digitized",
                        help="Output directory for scaled files")
    args = parser.parse_args()

    source_dir = Path(args.source_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Population names
    pops = ["C", "P", "E", "N", "H", "R"]

    print(f"\n{'='*60}")
    print(f"Scale digitized references")
    print(f"Source: {source_dir}")
    print(f"Output: {out_dir}")
    print(f"{'='*60}\n")

    for scale in args.scales:
        print(f"Scale: S={scale:.0e}")
        scale_str = f"S{scale:.0e}".replace("+0", "").replace("e", "e")

        for pop in pops:
            input_csv = source_dir / f"{pop}-Cells.csv"
            if not input_csv.exists():
                print(f"  {pop}: ✗ NOT FOUND ({input_csv})")
                continue

            output_csv = out_dir / f"{pop}-Cells-{scale_str}.csv"
            n_rows = scale_reference(input_csv, scale, output_csv)
            print(f"  {pop}: ✓ ({n_rows} rows) → {output_csv.name}")

        print()

    print(f"{'='*60}")
    print(f"✓ All scaled references generated")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
