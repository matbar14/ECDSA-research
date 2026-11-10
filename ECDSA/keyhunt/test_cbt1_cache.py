#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CBT1 validation-cache + plan-fingerprint: пять тестов новой механики.

  [1] cache write + hit: первый старт делает full pass и пишет sidecar
      "<table>.cache"; второй старт с неизменным size+mtime делает cache
      hit и восстанавливает ту же цель (скаляр бит-в-бит);
  [2] mtime-инвалидация: utime по файлу -> cache miss -> полный проход
      заново, sidecar перезаписан с новым mtime;
  [3] нет протухшего пропуска: битый payload при свежем mtime отвергается
      (exit 2, payload checksum), кеш не может «провалить» битый файл;
  [4] fingerprint-привязка плана к таблице: resume с ДРУГОЙ валидной
      CBT1 (иные байты -> иной payload FNV) отвергается как чужой план;
      resume с той же таблицей принимается (already-recovered путь);
  [5] контракт sidecar: все поля на месте, sha256 совпадает с реальным
      дайджестом файла, row count/file bytes/payload fnv согласованы с
      заголовком.

Exit 0 iff все обязательные тесты зелёные.
"""

import os
import re
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
START = 0x4000000000000000000000000000000000
CANON_100K = Path("../tables/representatives_canonical_H100000.txt")

failures = []


def check(name, ok, detail=""):
    print(f"[{name}] {'PASS' if ok else 'FAIL'}{(' ' + detail) if detail else ''}")
    if not ok:
        failures.append(name)


def fnv1a64(data: bytes) -> int:
    h = FNV_OFFSET
    for b in data:
        h = ((h ^ b) * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return h


def get_u32(buf, off):
    return struct.unpack_from(">I", buf, off)[0]


def get_u64(buf, off):
    return struct.unpack_from(">Q", buf, off)[0]


def pubkey_hex(scalar: int) -> str:
    from cryptography.hazmat.primitives.asymmetric import ec
    numbers = ec.derive_private_key(scalar % N, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(args_list, timeout=180):
    proc = subprocess.run(args_list, text=True, capture_output=True,
                          timeout=timeout)
    return proc


def parse_kv(path: Path):
    out = {}
    for line in path.read_text().splitlines():
        if ": " in line:
            k, v = line.split(": ", 1)
            out[k] = v
    return out


def counter(out, name):
    m = re.search(rf"{re.escape(name)}: (\d+)", out)
    return int(m.group(1)) if m else None


def validation_mode(out):
    m = re.search(r"table binary validation: (.+)", out)
    return m.group(1).strip() if m else None


def main():
    binary = Path("./eng_bsgs").resolve()
    tmp = Path(tempfile.mkdtemp(prefix="cbt1_cache_"))

    # ---- fixture: fresh CBT1 from the canonical 100k file ----
    heads = 100_000
    target_scalar = START + 2 * heads + 777
    base = [str(binary), "--target", pubkey_hex(target_scalar),
            "--start", hex(START), "--size", "500000000000",
            "--segment", "500000000000", "--threads", "2",
            "--progress-interval", "0"]
    dump_path = tmp / "fixture.cbt1"
    proc = run(base + ["--table-scalars", str(CANON_100K.resolve()),
                       "--glv-orbit", "--canonical-k3", "--partner-table",
                       "--dump-table-binary", str(dump_path)])
    assert proc.returncode == 0 and dump_path.exists(), proc.stderr[-500:]
    raw = dump_path.read_bytes()
    print(f"[fixture] {dump_path} {len(raw)} B")

    # ---- [1] cache write + hit ----
    t1 = tmp / "t1.cbt1"
    t1.write_bytes(raw)
    cache1 = Path(str(t1) + ".cache")
    p_a = run(base + ["--table-binary", str(t1), "--glv-orbit",
                      "--partner-table"])
    mode_a = validation_mode(p_a.stdout)
    side_ok = cache1.exists() and (side_before := parse_kv(cache1)).get(
        "cbt1 validation cache version") == "1"
    p_b = run(base + ["--table-binary", str(t1), "--glv-orbit",
                      "--partner-table"])
    mode_b = validation_mode(p_b.stdout)
    side_unchanged = parse_kv(cache1) == side_before
    rec_a = re.search(r"recovered private key: 0x([0-9a-f]+)", p_a.stdout)
    rec_b = re.search(r"recovered private key: 0x([0-9a-f]+)", p_b.stdout)
    ok = (p_a.returncode == 0 and p_b.returncode == 0
          and mode_a is not None and mode_a.startswith("full pass")
          and mode_b is not None and mode_b.startswith("cache hit")
          and side_ok and side_unchanged
          and rec_a and rec_b and rec_a.group(1) == rec_b.group(1)
          and int(rec_a.group(1), 16) == target_scalar)
    check("1", ok, f"first='{mode_a}', second='{mode_b}', "
                   f"sidecar={'yes' if side_ok else 'NO'}, "
                   f"recovery equal={bool(rec_a and rec_b and rec_a.group(1) == rec_b.group(1))}")

    # ---- [2] mtime invalidation -> full revalidation ----
    t2 = tmp / "t2.cbt1"
    t2.write_bytes(raw)
    cache2 = Path(str(t2) + ".cache")
    run(base + ["--table-binary", str(t2), "--glv-orbit", "--partner-table"])
    old_mtime = parse_kv(cache2)["mtime ns"]
    st = t2.stat()
    new_ns = (st.st_atime_ns, st.st_mtime_ns + 10**9)
    os.utime(t2, ns=new_ns)
    p_c = run(base + ["--table-binary", str(t2), "--glv-orbit",
                      "--partner-table"])
    mode_c = validation_mode(p_c.stdout)
    new_mtime = parse_kv(cache2)["mtime ns"] if cache2.exists() else "?"
    ok = (p_c.returncode == 0 and mode_c is not None
          and mode_c.startswith("full pass")
          and new_mtime == str(new_ns[1]))
    check("2", ok, f"mode='{mode_c}', sidecar mtime updated={new_mtime}")

    # ---- [3] corrupted payload with fresh mtime -> reject, no stale pass ----
    t3 = tmp / "t3.cbt1"
    buf = bytearray(raw)
    buf[300] ^= 0xFF
    t3.write_bytes(bytes(buf))
    p_d = run(base + ["--table-binary", str(t3), "--glv-orbit"])
    ok = (p_d.returncode == 2
          and "payload checksum" in p_d.stderr
          and validation_mode(p_d.stdout) is None)
    check("3", ok, f"exit={p_d.returncode}, "
                   f"msg={'yes' if 'payload checksum' in p_d.stderr else 'MISSING'}")

    # ---- [4] fingerprint binding: resume rejects a different valid table ----
    t4 = tmp / "t4.cbt1"
    t4.write_bytes(raw)
    chunk = 10_000_000_000
    target4 = START + chunk + 7_777_777          # chunk 1, offset < H
    cover_extra = ["--start", hex(START), "--size", str(chunk),
                   "--segment", str(chunk), "--threads", "2",
                   "--progress-interval", "0",
                   "--cover-start", hex(START),
                   "--cover-size", hex(4 * chunk)]
    cp4 = tmp / "cp4.txt"
    p_e = run([str(binary), "--target", pubkey_hex(target4)]
              + cover_extra
              + ["--table-binary", str(t4), "--glv-orbit", "--partner-table",
                 "--checkpoint-file", str(cp4), "--checkpoint-steps", "1"])
    run_ok = (p_e.returncode == 0
              and re.search(r"recovered private key: 0x([0-9a-f]+)", p_e.stdout)
              and int(re.search(r"recovered private key: 0x([0-9a-f]+)",
                                p_e.stdout).group(1), 16) == target4)
    # Craft a VALID second table: flip one scalar byte, recompute both FNVs.
    t4b = tmp / "t4b.cbt1"
    buf = bytearray(raw)
    buf[256 + 5000 * 65 + 40] ^= 0x01
    payload_fnv = fnv1a64(bytes(buf[256:]))
    struct.pack_into(">Q", buf, 240, payload_fnv)
    struct.pack_into(">Q", buf, 248, fnv1a64(bytes(buf[:248])))
    t4b.write_bytes(bytes(buf))
    p_f = run([str(binary), "--target", pubkey_hex(target4)] + cover_extra
              + ["--table-binary", str(t4b), "--glv-orbit", "--partner-table",
                 "--checkpoint-file", str(cp4), "--resume"])
    reject_ok = (p_f.returncode == 2
                 and "different plan" in p_f.stderr)
    p_g = run([str(binary), "--target", pubkey_hex(target4)] + cover_extra
              + ["--table-binary", str(t4), "--glv-orbit", "--partner-table",
                 "--checkpoint-file", str(cp4), "--resume"])
    accept_ok = (p_g.returncode == 0
                 and "already recovered" in p_g.stdout)
    check("4", run_ok and reject_ok and accept_ok,
          f"cover run={run_ok}, other table rejected={reject_ok}, "
          f"same table accepted={accept_ok}")

    # ---- [5] sidecar contract ----
    fields = parse_kv(cache1)
    import hashlib
    sha_real = hashlib.sha256(t1.read_bytes()).hexdigest()
    hdr_fnv = get_u64(raw, 240)
    ok = (fields.get("file bytes") == str(len(raw))
          and fields.get("row count") == str(heads)
          and fields.get("file sha256") == sha_real
          and fields.get("payload fnv1a64") == f"{hdr_fnv:016x}"
          and int(fields.get("mtime ns", "-1")) == t1.stat().st_mtime_ns)
    check("5", ok, f"sha256 match={fields.get('file sha256') == sha_real}, "
                   f"rows={fields.get('row count')}, "
                   f"fnv match={fields.get('payload fnv1a64') == f'{hdr_fnv:016x}'}")

    print(f"\nCBT1 cache/fingerprint tests: "
          f"{'ALL PASS' if not failures else 'FAILED: ' + str(failures)}"
          f"  (files in {tmp})")
    sys.exit(0 if not failures else 1)


if __name__ == "__main__":
    main()
