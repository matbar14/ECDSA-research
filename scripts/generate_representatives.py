#!/usr/bin/env python3
"""Файл канонических представителей орбит для длинных прогонов.

Головы j = 1..H (j=0 даёт нулевую орбиту, движок её отвергает),
представитель c(j) = min_r min(lambda^r*j mod n, n - lambda^r*j mod n) -
тот же canonical_orbit_scalar, что в движке и в sweep-скрипте.
Представители разбросаны по всему [1, n): таблица покрывает окно
W = stride * giant_count, stride = числу строк файла (= H), поэтому
одна и та же таблица годится для любого --size <= H * giant_max.

Использование:
  generate_representatives.py --heads 1000000 --output canonical_representatives_H1000000.txt
"""
import argparse
from pathlib import Path

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def canonical_representative(j):
    best = None
    v = j % N
    for _ in range(3):
        candidate = min(v, N - v) if v else 0
        if best is None or candidate < best:
            best = candidate
        v = v * LAMBDA % N
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--heads", type=int, default=1_000_000)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    out = Path(args.output)
    chunk = []
    with out.open("w") as f:
        f.write(f"# canonical orbit representatives, heads j=1..{args.heads}, "
                f"c(j)=min_r min(lambda^r*j mod n, n-lambda^r*j mod n), "
                f"rows={args.heads}, stride=rows\n")
        for j in range(1, args.heads + 1):
            # формат как в table_representative_sweep.write_scalars:
            # hex(s) -> "0x..." (parse_int движка берёт ветку SetBase16)
            chunk.append(f"0x{canonical_representative(j):x}")
            if len(chunk) == 100_000:
                f.write("\n".join(chunk) + "\n")
                chunk.clear()
        if chunk:
            f.write("\n".join(chunk) + "\n")
    size_mb = out.stat().st_size / 2**20
    print(f"[done] {out}: {args.heads} строк, {size_mb:.1f} MB")


if __name__ == "__main__":
    main()
