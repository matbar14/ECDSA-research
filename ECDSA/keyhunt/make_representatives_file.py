#!/usr/bin/env python3
"""Генератор файла канонических представителей полного диапазона.

PROGRESS_MEMO.md («Прямой canonical-проход», следующий шаг): перенос
канонического списка представителей в файл таблицы полного диапазона для
длинных прогонов. Скрипт строит persistent-файл для `--table-scalars`:
один канонический представитель на GLV-орбиту, H строк,

    c(j) = min over r in {0,1,2} of min(lambda^r*j mod n, n - lambda^r*j mod n),

головы j = 1..H. Представители разных голов различны (орбиты не пересекаются
для j < lambda-порядка), что проверяется генератором.

Файл предназначен для длинных прогонов движка:

    ./eng_bsgs --target PUB --start ... --size ... --glv-orbit \
      --phase-signature --partner-table --table-scalars FILE

Гигант-стрид равен числу строк файла (H), каждая строка покрывает все
6 скаляров своей орбиты: K3 = x^3 инвариантен к phi-вращению и отрицанию,
partner-expansion восстанавливает цель из маршрутизованной пары.

Головы j = 1..H детерминированы, поэтому файлы вложены: список голов
1..H - префикс списка 1..H' при H < H'. Отдельные файлы пишутся для
удобства запуска.
"""

import argparse
import hashlib
import json
import time
from datetime import date
from pathlib import Path

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72


def canonical_representative(j):
    """c = min_r min(lambda^r*j mod n, n - lambda^r*j mod n) - зеркало
    canonical_orbit_scalar из eng_bsgs.cpp."""
    best = None
    v = j % N
    for _ in range(3):
        candidate = min(v, N - v) if v else 0
        if best is None or candidate < best:
            best = candidate
        v = v * LAMBDA % N
    return best


def bitlen_stats(values):
    bits = [v.bit_length() for v in values]
    bits_sorted = sorted(bits)
    mid = len(bits_sorted) // 2
    median = bits_sorted[mid]
    if len(bits_sorted) % 2 == 0:
        median = (bits_sorted[mid - 1] + bits_sorted[mid]) / 2
    return {
        "min bit length": min(bits),
        "max bit length": max(bits),
        "mean bit length": round(sum(bits) / len(bits), 2),
        "median bit length": median,
        "below 2^254 share": round(
            sum(1 for v in values if v < (1 << 254)) / len(values), 4),
    }


def build_file(H, out_path):
    t0 = time.perf_counter()
    reps = [canonical_representative(j) for j in range(1, H + 1)]
    gen_seconds = round(time.perf_counter() - t0, 3)

    distinct = len(set(reps))
    if distinct != H:
        raise SystemExit(
            f"[E] canonical collision: {H} heads -> {distinct} distinct "
            f"representatives (orbits must be disjoint)")

    # Для голов j <= H канонический представитель - сама голова: партнёры
    # lambda^r*j mod n порядка n и минимум по орбите выбирает j. Ядро
    # покрытия - непрерывный тайлинг i*H + [1, H]; полный диапазон добавляют
    # орбитальные партнёры +-lambda^a*j на каждой гигант-позиции.
    equal_head = sum(1 for j, c in enumerate(reps, start=1) if c == j)

    header = [
        "# Canonical GLV-orbit representatives: full-range table file",
        f"# construction date: {date.isoformat(date.today())}",
        f"# heads: j = 1..{H}  (deterministic; files for smaller H are prefixes)",
        "# c(j) = min over r in {0,1,2} of min(lambda^r*j mod n, n - lambda^r*j mod n)",
        f"# rows: {H}  (one per orbit; covers 6 scalars per row via K3 = x^3)",
        f"# core structure: c(j) = j for {equal_head}/{H} heads -> contiguous"
        f" tiling i*stride + [1, {H}]",
        "# partner coverage: offsets +-lambda^a*j (a in {0,1,2}) recoverable"
        " at every giant position via --partner-table",
        f"# curve: secp256k1",
        f"# n = {N:064x}",
        f"# lambda = {LAMBDA:064x}",
        "# usage: eng_bsgs --glv-orbit --phase-signature --partner-table",
        "#        --table-scalars THIS_FILE  (giant stride = row count = H)",
        "# engine skips '#' comments; one hex scalar per line; zero rejected",
    ]
    body = "\n".join(header + [hex(c) for c in reps]) + "\n"
    out_path.write_text(body)
    digest = hashlib.sha256(out_path.read_bytes()).hexdigest()
    file_bytes = out_path.stat().st_size

    stats = {
        "heads": H,
        "rows": H,
        "distinct representatives": distinct,
        "representatives equal heads": equal_head,
        "file": out_path.name,
        "file bytes": file_bytes,
        "sha256": digest,
        "generation seconds": gen_seconds,
        **bitlen_stats(reps),
        "largest representative share of n/2": round(max(reps) / (N / 2), 12),
    }
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--heads", type=int, nargs="+",
                        default=[50_000, 100_000, 200_000, 400_000],
                        help="число голов орбит H для каждого файла")
    parser.add_argument("--output-dir", default="../tables",
                        help="каталог для representatives_*.txt и summary JSON")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    all_stats = {"curve": "secp256k1", "files": []}
    for H in args.heads:
        path = out_dir / f"representatives_canonical_H{H}.txt"
        stats = build_file(H, path)
        all_stats["files"].append(stats)
        print(f"[file] {path.name}: rows={stats['rows']}, "
              f"{stats['file bytes']} bytes, "
              f"gen {stats['generation seconds']}s, "
              f"bits mean={stats['mean bit length']}, "
              f"sha256={stats['sha256'][:16]}...")

    summary = out_dir / "representatives_files.json"
    summary.write_text(json.dumps(all_stats, indent=2) + "\n")
    print(f"[summary] {summary}")


if __name__ == "__main__":
    main()
