#!/usr/bin/env python3
"""Generate the gap-free mirror-tiling representative file.

For R rows the heads are C = [R, 2R-1].  With giant stride 2R-1,
the signed residues +/-C cover every residue class, so the same table
covers almost twice as many consecutive giant positions as the default
stride-R layout.  The engine must be run with --glv-orbit and
--giant-stride 2R-1.
"""

import argparse
import hashlib
import json
import time
from datetime import date
from pathlib import Path


N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def canonical_representative(value):
    best = None
    current = value % N
    for _ in range(3):
        best = min(best, current, N - current) if best is not None else min(current, N - current)
        current = current * LAMBDA % N
    return best


def build_file(rows, output):
    if rows < 1:
        raise SystemExit("rows must be positive")
    first = rows
    last = 2 * rows - 1
    t0 = time.perf_counter()
    heads = list(range(first, last + 1))
    if any(canonical_representative(value) != value for value in heads):
        raise SystemExit("mirror heads are not canonical representatives")

    header = [
        "# Canonical GLV mirror-tiling representatives",
        f"# construction date: {date.isoformat(date.today())}",
        f"# rows: {rows}; heads: [{first}, {last}]",
        f"# giant stride: {2 * rows - 1} = 2*rows-1",
        "# signed heads cover every residue class modulo the giant stride",
        "# usage: --glv-orbit --phase-signature --partner-table",
        f"# n = {N:064x}",
        f"# lambda = {LAMBDA:064x}",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(header + [hex(value) for value in heads]) + "\n"
    output.write_text(body)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    return {
        "kind": "canonical mirror tiling",
        "rows": rows,
        "first head": first,
        "last head": last,
        "giant stride": 2 * rows - 1,
        "file": output.name,
        "bytes": output.stat().st_size,
        "sha256": digest,
        "generation seconds": round(time.perf_counter() - t0, 3),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=int, required=True,
                        help="number of table rows R")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    stats = build_file(args.rows, args.output)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
