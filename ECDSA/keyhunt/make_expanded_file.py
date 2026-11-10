#!/usr/bin/env python3
"""Генератор expanded-файла представителей (вход merge-режима движка).

Прототип CBT1 замеряется против текстового пути: expanded-файл - это
pre-canonicalization сырьё merge-pipeline, 3 строки на голову
(j, lambda*j mod n, lambda^2*j mod n), H голов -> 3H строк. Конструкция
зеркалит существующий representatives_expanded_H400000.txt (заголовок,
формат "0x%x", порядок строк). Файлы детерминированы и вложены: головы
1..H - префикс голов 1..H'.

Использование движком (merge-режим, результат - те же H канонических строк):

    ./eng_bsgs --glv-orbit --canonical-k3 --table-scalars FILE

Канонический минимум орбиты c(j) = j для всех голов j = 1..H проверен
make_representatives_file.py (equal-heads 100%); здесь каноничность
входа дополнительно контролируется на sample каждые SAMPLE_EVERY голов.
"""

import argparse
import hashlib
import json
import time
from datetime import date
from pathlib import Path

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
SAMPLE_EVERY = 997          # простое число: несинхронно с периодом орбиты
CHUNK_HEADS = 100_000       # строк в буфере записи


def canonical_min(v_first: int) -> int:
    """c = min over r of min(lambda^r*v, n - lambda^r*v) - контроль входа."""
    best = None
    v = v_first % N
    for _ in range(3):
        for cand in (v, (-v) % N):
            if best is None or cand < best:
                best = cand
        v = v * LAMBDA % N
    return best


def build_file(H: int, out_path: Path):
    t0 = time.perf_counter()
    digest = hashlib.sha256()
    rows_total = 3 * H
    sample_checked = 0
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")

    header = (
        f"# expanded representatives, heads 1..{H}, 3 rows per head\n"
        f"# rows: {rows_total} (stride = rows = {rows_total})\n"
    ).encode("ascii")
    digest.update(header)

    mul = LAMBDA
    with open(tmp_path, "wb") as f:
        f.write(header)
        buf = []
        append = buf.append
        for j in range(1, H + 1):
            j2 = j * mul % N          # lambda*j mod n
            j3 = j2 * mul % N         # lambda^2*j mod n
            if j % SAMPLE_EVERY == 0:
                c = canonical_min(j)
                if c != j:
                    tmp_path.unlink(missing_ok=True)
                    raise SystemExit(
                        f"[E] head {j}: canonical min {c} != head j "
                        f"(expanded input must canonicalize to heads)")
                sample_checked += 1
            append(f"0x{j:x}\n")
            append(f"0x{j2:x}\n")
            append(f"0x{j3:x}\n")
            if len(buf) >= 3 * CHUNK_HEADS:
                chunk = "".join(buf).encode("ascii")
                f.write(chunk)
                digest.update(chunk)
                buf.clear()
        if buf:
            chunk = "".join(buf).encode("ascii")
            f.write(chunk)
            digest.update(chunk)
            buf.clear()
    tmp_path.rename(out_path)

    gen_seconds = round(time.perf_counter() - t0, 3)
    stats = {
        "kind": "expanded (3 rows per head, merge-mode input)",
        "heads": H,
        "rows": rows_total,
        "file": out_path.name,
        "file bytes": out_path.stat().st_size,
        "sha256": digest.hexdigest(),
        "generation seconds": gen_seconds,
        "sample canonical checks": sample_checked,
        "construction date": date.isoformat(date.today()),
    }
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--heads", type=int, default=4_000_000,
                        help="число голов орбит H (3H строк в файле)")
    parser.add_argument("--output-dir", default="../tables",
                        help="каталог для representatives_expanded_H*.txt")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    H = args.heads
    path = out_dir / f"representatives_expanded_H{H}.txt"
    if path.exists():
        raise SystemExit(f"[E] refuse to overwrite existing file: {path}")
    stats = build_file(H, path)

    summary = out_dir / "representatives_files.json"
    summary_data = json.loads(summary.read_text()) if summary.exists() else {}
    summary_data.setdefault("curve", "secp256k1")
    summary_data.setdefault("expanded_files", []).append(stats)
    summary.write_text(json.dumps(summary_data, indent=2) + "\n")

    print(f"[file] {path.name}: rows={stats['rows']}, "
          f"{stats['file bytes']} bytes, gen {stats['generation seconds']}s, "
          f"sample checks {stats['sample canonical checks']}, "
          f"sha256={stats['sha256'][:16]}...")
    print(f"[summary] {summary}")


if __name__ == "__main__":
    main()
