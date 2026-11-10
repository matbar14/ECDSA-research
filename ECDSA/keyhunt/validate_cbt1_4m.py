#!/usr/bin/env python3
"""Независительная валидация CBT1 на H=4M (полный файл, sample по ключам).

Проверки (зеркалят test_binary_table.py на H=100k, но на полном 4M-файле):
  [A] header: все поля big-endian по смещениям dump_table_binary +
      header FNV-1a (248 байт, чистый python);
  [B] payload FNV-1a на ВСЕХ 260 MB - пересчёт независимым C-инструментом
      (fnv1a64_tool.c, чистый цикл без SIMD);
  [C] сортировка по K3-ключу на всех 4M строках + отсутствие дубликатов
      ключа (строгий <);
  [D] скалярная колонка: multiset {скаляры} == {1..H} - полное сравнение
      (канонический минимум каждой орбиты = голове, проверено
      make_representatives_file.py: c(j) = j на 4M/4M голов);
  [E] K3-ключи: случайный sample 512 строк - key == (x^3 mod p) от
      x(c*G), пересчёт в python (cryptography), 33 байта бит-в-бит.

Артефакт пишется атомарно. Exit 0 iff все обязательные проверки зелёные.
"""

import json
import random
import struct
import subprocess
import sys
import time
from pathlib import Path

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
ROOT = Path(__file__).resolve().parent.parent
TABLE = ROOT / "tables" / "table_cbt1_H4000000.cbt1"
TOOL = ROOT / "keyhunt" / "fnv1a64_tool"
ARTIFACT = ROOT / "results" / "cbt1_4m_validation.json"
SAMPLE = 512

failures = []


def check(name, ok, detail=""):
    print(f"[{name}] {'PASS' if ok else 'FAIL'}{(' ' + detail) if detail else ''}")
    if not ok:
        failures.append(name)


def fnv1a64_bytes(data: bytes) -> int:
    h = FNV_OFFSET
    for b in data:
        h = ((h ^ b) * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return h


def main():
    t0 = time.monotonic()
    raw = TABLE.read_bytes()
    size = len(raw)
    hdr = raw[:256]
    payload = raw[256:]

    row_count = struct.unpack_from(">Q", hdr, 32)[0]
    check("size", size == 256 + 65 * row_count,
          f"{size} B == 256 + 65*{row_count}")

    # ---- [A] header fields ----
    fields = {
        "magic": hdr[:4] == b"CBT1",
        "canonical_version": struct.unpack_from(">I", hdr, 4)[0] == 1,
        "payload_kind(giant_grid)": struct.unpack_from(">I", hdr, 8)[0] == 1,
        "header_size": struct.unpack_from(">I", hdr, 12)[0] == 256,
        "row_size": struct.unpack_from(">I", hdr, 16)[0] == 65,
        "curve_id": struct.unpack_from(">I", hdr, 20)[0] == 1,
        "orbit_size": struct.unpack_from(">I", hdr, 24)[0] == 6,
        "flags(sorted)": struct.unpack_from(">I", hdr, 28)[0] == 1,
        "row_count": row_count == 4_000_000,
        "stride": struct.unpack_from(">Q", hdr, 40)[0] == 4_000_000,
        "m": struct.unpack_from(">Q", hdr, 48)[0] == row_count,
        "n": int.from_bytes(hdr[80:112], "big") == N,
        "p": int.from_bytes(hdr[112:144], "big") == P,
        "lambda": int.from_bytes(hdr[144:176], "big") == LAMBDA,
        "rows_per_head_source": struct.unpack_from(">I", hdr, 72)[0] == 1,
        "created_unix": 1_760_000_000 < struct.unpack_from(">Q", hdr, 176)[0]
                        < 2_000_000_000,
    }
    bad = [k for k, v in fields.items() if not v]
    hdr_fnv_ok = struct.unpack_from(">Q", hdr, 248)[0] == fnv1a64_bytes(hdr[:248])
    check("header", not bad and hdr_fnv_ok,
          f"bad={bad}, header_fnv={hdr_fnv_ok}, "
          f"created={struct.unpack_from('>Q', hdr, 176)[0]}")

    # ---- [B] payload FNV-1a via independent C tool ----
    tool_out = subprocess.run(
        [str(TOOL), str(TABLE), "256", str(len(payload))],
        capture_output=True, text=True, timeout=120)
    tool_digest = tool_out.stdout.strip()
    stored_payload_fnv = f"{struct.unpack_from('>Q', hdr, 240)[0]:016x}"
    check("payload_fnv", tool_out.returncode == 0
          and tool_digest == stored_payload_fnv,
          f"{tool_digest} vs stored {stored_payload_fnv}")

    # ---- [C] sorted ascending, no duplicate keys (all 4M rows) ----
    rows = memoryview(payload)
    sorted_ok = True
    dup_keys = 0
    prev_key = None
    for i in range(row_count):
        off = i * 65
        key = payload[off:off + 33]
        if prev_key is not None:
            if key < prev_key:
                sorted_ok = False
                print(f"  order violation at row {i}")
                break
            if key == prev_key:
                dup_keys += 1
        prev_key = key
    check("sorted+unique", sorted_ok and dup_keys == 0,
          f"duplicate keys {dup_keys}")

    # ---- [D] scalar column: multiset == heads 1..H (full 4M) ----
    t1 = time.monotonic()
    scalars = [int.from_bytes(payload[i * 65 + 33:i * 65 + 65], "big")
               for i in range(row_count)]
    scalars.sort()
    expected = list(range(1, row_count + 1))
    scalars_ok = scalars == expected
    check("scalar_column", scalars_ok,
          f"4M scalars == heads 1..{row_count} bit-for-bit "
          f"({time.monotonic() - t1:.1f} s)")

    # ---- [E] K3 keys on a random sample ----
    from cryptography.hazmat.primitives.asymmetric import ec
    rng = random.Random(20261031)
    sample_idx = rng.sample(range(row_count), SAMPLE)
    keys_ok = True
    first_bad = None
    for i in sample_idx:
        off = i * 65
        c = int.from_bytes(payload[off + 33:off + 65], "big")  # scalar OF row i
        x = ec.derive_private_key(c, ec.SECP256K1()) \
            .public_key().public_numbers().x
        key = b"\x00" + pow(x, 3, P).to_bytes(32, "big")
        if key != payload[off:off + 33]:
            keys_ok = False
            first_bad = (i, c)
            break
    check("k3_keys_sample", keys_ok,
          f"{SAMPLE}/{SAMPLE} rows: key == x^3 mod p of x(c*G)"
          + ("" if keys_ok else f", first bad {first_bad}"))

    wall = time.monotonic() - t0
    print(f"\nCBT1 H=4M validation: "
          f"{'ALL PASS' if not failures else 'FAILED: ' + str(failures)} "
          f"({wall:.1f} s)")

    artifact = {
        "file": TABLE.name,
        "bytes": size,
        "sha256": json.loads((ROOT / "results" / "cbt1_4m_build.json").read_text())
        ["dump file"]["sha256"],
        "row_count": row_count,
        "checks": {k: v if isinstance(v, bool) else "see stdout"
                   for k, v in fields.items()},
        "header_fnv ok": hdr_fnv_ok,
        "payload_fnv tool": tool_digest,
        "payload_fnv stored": stored_payload_fnv,
        "sorted ascending": sorted_ok,
        "duplicate keys": dup_keys,
        "scalar column == heads 1..H": scalars_ok,
        "k3 keys sample": {"rows": SAMPLE, "ok": keys_ok,
                           "seed": 20261031},
        "validation seconds": round(wall, 3),
        "failures": failures,
    }
    tmp = ARTIFACT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(artifact, indent=1))
    tmp.rename(ARTIFACT)
    print(f"[artifact] {ARTIFACT}")
    sys.exit(0 if not failures else 1)


if __name__ == "__main__":
    main()
