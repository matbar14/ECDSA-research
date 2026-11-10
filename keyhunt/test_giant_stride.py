#!/usr/bin/env python3
"""Black-box test for --giant-stride: mirror-tiling of the giant grid.

Coverage model (mod-stride counting): with candidate rows C and giant
positions i*stride, a relative offset t is coverable iff t mod stride lies
in +/-C mod stride. Gapless tiling therefore needs |±C mod stride| = stride,
bounded by 2*rows-1 (the even stride 2*rows is impossible: residues 0 and
rows are both self-paired and need two dedicated rows). The construction
C = [R, 2R-1] (consecutive heads, canonical scalar == head) reaches the
ceiling: stride 2R-1, gapless, half the giant steps of the default
stride == rows on the same table.

Cases:
  [1] regression: continuous heads [1, R], default stride, target recovered
  [2] mirror-tiling heads [R, 2R-1], default stride (== R): recovered
  [3] mirror-tiling + --giant-stride 2R-1: targets across residue classes
      (boundaries, +/-1, tail, first offset) recovered; giant steps halved
  [4] --giant-stride 4R-2 > 2R-1: gap target NOT recovered (negative control),
      in-cover target recovered (positive control)
  [5] --giant-stride without --glv-orbit rejected with exit code 2
"""

import argparse
import random
import re
import subprocess
import tempfile
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

START = 1 << 134
SIZE = 1 << 20
R = 1000
STRIDE = 2 * R - 1
BAD_STRIDE = 4 * R - 2


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(binary, target, size, extra):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(size), "--segment", str(size), "--threads", "1"] + extra
    result = subprocess.run(cmd, text=True, capture_output=True)
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]+)",
                      result.stdout)
    recovered = int(match.group(1), 16) if match else None
    counters = {}
    for name in ("giant steps", "giant steps planned"):
        m = re.search(rf"{re.escape(name)}: (\d+)", result.stdout)
        if m:
            counters[name] = int(m.group(1))
    return recovered, counters, result


def write_heads_file(directory, first, count, name):
    path = Path(directory) / name
    lines = [f"# consecutive heads [{first}, {first + count - 1}]"]
    lines += [hex(s) for s in range(first, first + count)]
    path.write_text("\n".join(lines) + "\n")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    args = parser.parse_args()

    failures = []
    with tempfile.TemporaryDirectory() as tmpdir:
        core_rows = write_heads_file(tmpdir, 1, R, "heads_core.txt")
        mirror_rows = write_heads_file(tmpdir, R, R, "heads_mirror.txt")

        table_common = ["--glv-orbit", "--canonical-k3"]

        def check(name, cond, detail=""):
            status = "PASS" if cond else "FAIL"
            print(f"[{name}] {status} {detail}")
            if not cond:
                failures.append(name)

        # [1] regression: continuous heads, default stride
        t = 777
        rec, cnt, res = run(args.binary, pubkey_hex(START + t), SIZE,
                            ["--table-scalars", str(core_rows)] + table_common)
        check(1, rec == START + t and cnt.get("giant steps planned") == (SIZE - 1) // R + 1,
              f"t={t} recovered={rec == START + t} planned={cnt.get('giant steps planned')}")

        # [2] mirror-tiling heads, default stride (== R): still gapless
        t = 7 * R + 500
        rec, cnt, res = run(args.binary, pubkey_hex(START + t), SIZE,
                            ["--table-scalars", str(mirror_rows)] + table_common)
        giants_def = cnt.get("giant steps planned")
        check(2, rec == START + t,
              f"t={t} recovered={rec == START + t} planned={giants_def}")

        # [3] mirror-tiling + --giant-stride 2R-1: half the giants, all
        # residue classes covered
        rng = random.Random(20261110)
        targets = {1, STRIDE, 2 * STRIDE, 3 * STRIDE, STRIDE - 1,
                   STRIDE + 1, SIZE - 1, 7 * R + 500}
        targets |= {rng.randrange(1, SIZE) for _ in range(6)}
        ok = True
        giants_mirror = None
        for t in sorted(targets):
            rec, cnt, res = run(args.binary, pubkey_hex(START + t), SIZE,
                                ["--table-scalars", str(mirror_rows),
                                 "--giant-stride", str(STRIDE)] + table_common)
            if rec != START + t:
                ok = False
                print(f"    target t={t}: recovered={rec}")
            if t == 7 * R + 500:
                giants_mirror = cnt.get("giant steps planned")
        check(3, ok and giants_mirror is not None and giants_def is not None
              and giants_mirror == (SIZE - 1) // STRIDE + 2
              and 2 * giants_mirror <= giants_def + 4,
              f"{len(targets)} targets, planned default={giants_def} "
              f"mirror={giants_mirror} (expect ~half incl. tail giant)")

        # [4] stride beyond the ceiling: gaps guaranteed
        gap_t = 2 * BAD_STRIDE + 50          # t mod BAD_STRIDE = 50 -> gap
        hit_t = BAD_STRIDE + 1500            # t mod BAD_STRIDE in [R, 2R-1]
        rec_gap, cnt_gap, _ = run(args.binary, pubkey_hex(START + gap_t), SIZE,
                                  ["--table-scalars", str(mirror_rows),
                                   "--giant-stride", str(BAD_STRIDE)] + table_common)
        rec_hit, _, _ = run(args.binary, pubkey_hex(START + hit_t), SIZE,
                            ["--table-scalars", str(mirror_rows),
                             "--giant-stride", str(BAD_STRIDE)] + table_common)
        check(4, rec_gap is None and rec_hit == START + hit_t
              and cnt_gap.get("giant steps planned") == (SIZE - 1) // BAD_STRIDE + 2,
              f"gap t={gap_t} recovered={rec_gap is not None} (expect False, "
              f"full scan planned={cnt_gap.get('giant steps planned')}), "
              f"hit t={hit_t} recovered={rec_hit == START + hit_t}")

        # [5] --giant-stride requires --glv-orbit
        result = subprocess.run(
            [args.binary, "--target", pubkey_hex(START + 1), "--start",
             hex(START), "--size", str(SIZE), "--giant-stride", str(STRIDE)],
            text=True, capture_output=True)
        check(5, result.returncode == 2 and
              "--giant-stride requires --glv-orbit" in result.stderr,
              f"exit={result.returncode}")

    print("ALL PASS" if not failures else f"FAILURES: {failures}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
