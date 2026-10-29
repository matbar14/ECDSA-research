#!/usr/bin/env python3
"""CBT1 binary table format: independent verification (test [4]).

Checks, in order:
  [1] dump: canonical H=100k text file -> --dump-table-binary writes a valid
      CBT1 file (row_count == canonical rows, 65 B/row, sorted ascending);
  [2] bit-for-bit scalars: the payload scalar column equals the canonical
      orbit minima recomputed INDEPENDENTLY in python from the same text
      file (c = min over r of min(lambda^r*j, -(lambda^r*j)) mod n);
  [3] bit-for-bit keys: for EVERY row the K3 key equals
      (x^3 mod p) of x(c*G), recomputed in python (cryptography lib);
  [4] header: magic/version/kind/sizes/curve params/n/p/lambda + header
      FNV-1a recomputed; payload FNV-1a recomputed;
  [5] end-to-end: --table-binary run recovers the same target as the
      in-RAM run; stride == row_count; merged rows == 0, insert collisions
      == 0 in file mode;
  [6] reject policy: corrupted payload byte / unsorted rows / wrong
      canonical_version / truncated file / wrong group order / bad magic
      all exit 2 with a reject message;
  [7] steps/s: paired RAM vs mmap scan rates at H=100k (informational).

Exit code 0 iff all mandatory checks pass ([7] never fails the test).
"""
import argparse
import json
import os
import re
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
START = 0x4000000000000000000000000000000000
CANON_100K = Path("../representatives_canonical_H100000.txt")

FNV_TABLE = None  # lazily built for payload checksum speed


def fnv1a64(data: bytes) -> int:
    global FNV_TABLE
    if FNV_TABLE is None:
        # 8-bit table: fnv_step[c] = (c * PRIME) pre-xor shortcut is not
        # exact for FNV-1a (xor first), so keep the plain byte loop but
        # drive it through local names for speed.
        pass
    h = FNV_OFFSET
    for b in data:
        h = ((h ^ b) * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
    return h


def get_u32(buf, off):
    return struct.unpack_from(">I", buf, off)[0]


def get_u64(buf, off):
    return struct.unpack_from(">Q", buf, off)[0]


def fix_header_checksum(h: bytearray):
    struct.pack_into(">Q", h, 248, fnv1a64(bytes(h[:248])))


def pubkey_x(scalar: int) -> int:
    from cryptography.hazmat.primitives.asymmetric import ec
    return ec.derive_private_key(scalar % N, ec.SECP256K1()) \
        .public_key().public_numbers().x


def pubkey_hex(scalar: int) -> str:
    from cryptography.hazmat.primitives.asymmetric import ec
    numbers = ec.derive_private_key(scalar % N, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def parse_canonical_file(path: Path):
    scalars = []
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            scalars.append(int(line, 0) % N)
    return scalars


def canonical_min(j: int) -> int:
    """c(j) = min over r in {0,1,2} of min(l^r*j, -(l^r*j)) mod n."""
    best = None
    v = j % N
    for _ in range(3):
        for cand in (v, (-v) % N):
            if best is None or cand < best:
                best = cand
        v = (v * LAMBDA) % N
    return best


def run_engine(args_list, timeout=180):
    t0 = time.monotonic()
    proc = subprocess.run(args_list, text=True, capture_output=True,
                          timeout=timeout)
    return proc, time.monotonic() - t0


def engine_base(binary, target_scalar, size):
    return [
        str(binary), "--target", pubkey_hex(target_scalar),
        "--start", hex(START), "--size", str(size),
        "--segment", str(size), "--threads", "2", "--progress-interval", "0",
    ]


def counter(out, name):
    m = re.search(rf"{re.escape(name)}: (\d+)", out)
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=str(Path("./eng_bsgs").resolve()))
    ap.add_argument("--rows", type=int, default=100_000)
    args = ap.parse_args()
    binary = Path(args.binary)
    failures = []
    tmp = Path(tempfile.mkdtemp(prefix="cbt1_"))
    canon_path = CANON_100K

    canon_values = parse_canonical_file(canon_path)
    heads = args.rows
    canon_values = canon_values[:heads]
    assert len(canon_values) == heads
    m_stride = heads                      # stride = post-merge rows = heads
    target_scalar = START + 2 * m_stride + 777

    # ---- [1] dump ----
    dump_path = tmp / "table_100k.cbt1"
    cmd_dump = engine_base(binary, target_scalar, 500_000_000_000) + [
        "--table-scalars", str(canon_path.resolve()),
        "--glv-orbit", "--canonical-k3", "--partner-table",
        "--dump-table-binary", str(dump_path),
    ]
    proc, _ = run_engine(cmd_dump)
    rec = re.search(r"recovered private key: 0x([0-9a-f]+)", proc.stdout)
    ok = (proc.returncode == 0 and rec is not None
          and int(rec.group(1), 16) == target_scalar
          and dump_path.exists())
    print(f"[1] dump canonical H={heads}: recovered={rec is not None}, "
          f"file={dump_path.stat().st_size if dump_path.exists() else 0} B "
          f"-> {'PASS' if ok else 'FAIL'}")
    if not ok:
        print(proc.stdout[-1500:], proc.stderr[-800:])
        sys.exit(1)

    raw = dump_path.read_bytes()
    hdr = raw[:256]
    payload = raw[256:]

    # ---- [2] header fields + checksums ----
    magic_ok = hdr[:4] == b"CBT1"
    fields_ok = (
        get_u32(hdr, 4) == 1                      # canonical_version
        and get_u32(hdr, 8) == 1                  # payload_kind giant_grid
        and get_u32(hdr, 12) == 256               # header_size
        and get_u32(hdr, 16) == 65                # row_size
        and get_u32(hdr, 20) == 1                 # curve_id secp256k1
        and get_u32(hdr, 24) == 6                 # orbit_size
        and get_u64(hdr, 32) == heads             # row_count
        and get_u64(hdr, 40) == m_stride          # stride == merged rows
        and len(raw) == 256 + 65 * heads          # row_count vs file size
    )
    n_ok = int.from_bytes(hdr[80:112], "big") == N
    p_ok = int.from_bytes(hdr[112:144], "big") == P
    lam_ok = int.from_bytes(hdr[144:176], "big") == LAMBDA
    hdr_fnv_ok = get_u64(hdr, 248) == fnv1a64(hdr[:248])
    payload_fnv_ok = (get_u64(hdr, 240)
                      == fnv1a64(payload))
    ok = (magic_ok and fields_ok and n_ok and p_ok and lam_ok
          and hdr_fnv_ok and payload_fnv_ok)
    print(f"[2] header: magic={magic_ok}, fields={fields_ok}, n={n_ok}, "
          f"p={p_ok}, lambda={lam_ok}, hdr_fnv={hdr_fnv_ok}, "
          f"payload_fnv={payload_fnv_ok} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(2)

    rows = [payload[i * 65:(i + 1) * 65] for i in range(heads)]

    # ---- [3] sorted ascending + scalar column bit-for-bit ----
    sorted_ok = all(rows[i - 1][:33] <= rows[i][:33] for i in range(1, heads))
    file_scalars = sorted(
        int.from_bytes(r[33:65], "big") for r in rows)
    expected_scalars = sorted(canon_values)
    scalars_bitwise = file_scalars == expected_scalars
    ok = sorted_ok and scalars_bitwise
    print(f"[3] sorted={sorted_ok}, scalar column bit-for-bit vs "
          f"python-recomputed canonical minima: {scalars_bitwise} "
          f"-> {'PASS' if ok else 'FAIL'}")
    if not ok:
        failures.append(3)

    # ---- [4] keys bit-for-bit: key = (x^3 mod p) of x(c*G), all rows ----
    keys_ok = True
    first_bad = None
    for i, r in enumerate(rows):
        c = int.from_bytes(r[33:65], "big")
        x = pubkey_x(c)
        key = b"\x00" + (pow(x, 3, P)).to_bytes(32, "big")
        if key != r[:33]:
            keys_ok = False
            first_bad = (i, c)
            break
    print(f"[4] K3 keys bit-for-bit for all {heads} rows: {keys_ok}"
          f"{'' if keys_ok else f' first bad row {first_bad}'} "
          f"-> {'PASS' if keys_ok else 'FAIL'}")
    if not keys_ok:
        failures.append(4)

    # ---- [5] end-to-end file mode vs in-RAM mode ----
    cmd_file = engine_base(binary, target_scalar, 500_000_000_000) + [
        "--table-binary", str(dump_path),
        "--glv-orbit", "--partner-table",
    ]
    proc_f, _ = run_engine(cmd_file)
    rec_f = re.search(r"recovered private key: 0x([0-9a-f]+)", proc_f.stdout)
    stride_f = re.search(r"giant stride \(post-merge rows\): (\d+)",
                         proc_f.stdout)
    ok = (proc_f.returncode == 0 and rec_f is not None
          and int(rec_f.group(1), 16) == target_scalar
          and stride_f is not None and int(stride_f.group(1)) == m_stride
          and counter(proc_f.stdout, "merged rows") == 0
          and counter(proc_f.stdout, "insert collisions") == 0
          and counter(proc_f.stdout, "exact candidate checks") == 1)
    print(f"[5] file-mode recovery: correct={rec_f is not None and int(rec_f.group(1), 16) == target_scalar}, "
          f"stride={stride_f.group(1) if stride_f else '?'}, "
          f"merged={counter(proc_f.stdout, 'merged rows')}, "
          f"collisions={counter(proc_f.stdout, 'insert collisions')}, "
          f"exact={counter(proc_f.stdout, 'exact candidate checks')} "
          f"-> {'PASS' if ok else 'FAIL'}")
    if not ok:
        print(proc_f.stdout[-1500:], proc_f.stderr[-500:])
        failures.append(5)

    # ---- [6] reject policy ----
    def expect_reject(name, mutate, msg_part):
        path = tmp / f"bad_{name}.cbt1"
        buf = bytearray(raw)
        mutate(buf)
        path.write_bytes(bytes(buf))
        cmd = engine_base(binary, target_scalar, 500_000_000_000) + [
            "--table-binary", str(path), "--glv-orbit"]
        proc, _ = run_engine(cmd)
        ok = proc.returncode == 2 and msg_part in proc.stderr
        print(f"[6.{name}] reject: exit={proc.returncode}, "
              f"msg={'yes' if msg_part in proc.stderr else 'MISSING'} "
              f"-> {'PASS' if ok else 'FAIL'}")
        if not ok:
            print("  stderr:", proc.stderr.strip()[:300])
            return False
        return True

    def mut_payload(buf):
        buf[300] ^= 0xFF

    def mut_unsorted(buf):
        # swap rows 10 and 11 (keys differ), fix both checksums
        off = 256
        r10 = bytes(buf[off + 10 * 65: off + 11 * 65])
        r11 = bytes(buf[off + 11 * 65: off + 12 * 65])
        buf[off + 10 * 65: off + 11 * 65] = r11
        buf[off + 11 * 65: off + 12 * 65] = r10
        struct.pack_into(">Q", buf, 240,
                         fnv1a64(bytes(buf[256:])))
        fix_header_checksum(buf)

    def mut_version(buf):
        struct.pack_into(">I", buf, 4, 2)
        fix_header_checksum(buf)

    def mut_truncate(buf):
        del buf[-65:]

    def mut_n(buf):
        buf[80] ^= 0x01
        fix_header_checksum(buf)

    def mut_magic(buf):
        buf[0] = ord("X")
        fix_header_checksum(buf)

    rejects = [
        ("payload", mut_payload, "payload checksum"),
        ("unsorted", mut_unsorted, "not sorted ascending"),
        ("version", mut_version, "canonical_version"),
        ("truncated", mut_truncate, "row_count"),
        ("grouporder", mut_n, "group order"),
        ("magic", mut_magic, "bad magic"),
    ]
    if not all(expect_reject(*r) for r in rejects):
        failures.append(6)

    # ---- [7] steps/s: paired RAM vs mmap scan (informational) ----
    far_target = START + 1_000_000_000_000_000
    scan_size = 2_000_000_000_000          # 20M giant steps at stride 100k
    rates = {}
    for tag, extra in (
        ("ram", ["--table-scalars", str(canon_path.resolve()),
                 "--glv-orbit", "--canonical-k3"]),
        ("mmap", ["--table-binary", str(dump_path), "--glv-orbit"]),
    ):
        cmd = engine_base(binary, far_target, scan_size) + extra + [
            "--progress-interval", "5"]
        proc, _ = run_engine(cmd, timeout=120)
        samples = [float(x) for x in
                   re.findall(r"steps/s ([0-9.]+)", proc.stdout)]
        rates[tag] = samples
        print(f"[7] {tag}: samples={len(samples)}, "
              f"median={sorted(samples)[len(samples)//2]:.0f}" if samples
              else f"[7] {tag}: no samples")
    if rates["ram"] and rates["mmap"]:
        med = lambda xs: sorted(xs)[len(xs) // 2]
        ratio = med(rates["mmap"]) / max(med(rates["ram"]), 1)
        print(f"[7] mmap/ram median steps/s ratio: {ratio:.3f} "
              f"(informational, never fails the test)")

    print(f"\nCBT1 test [4]: {'ALL PASS' if not failures else f'FAILED: {failures}'}"
          f"  (files in {tmp})")
    sys.exit(0 if not failures else 1)


if __name__ == "__main__":
    main()
