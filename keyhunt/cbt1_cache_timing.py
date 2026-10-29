#!/usr/bin/env python3
"""Холодный (full validation) против кешированного старта на CBT1 4M."""
import re
import subprocess
import sys
import time
from pathlib import Path

ENGINE = Path("./eng_bsgs").resolve()
TABLE = Path("../table_cbt1_H4000000.cbt1")
CACHE = Path(str(TABLE) + ".cache")
SCALAR = 0x4000000000000000000000000000000000 + 7 * 4_000_000 + 4242

from cryptography.hazmat.primitives.asymmetric import ec
numbers = ec.derive_private_key(SCALAR, ec.SECP256K1()) \
    .public_key().public_numbers()
pub = f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"

CMD = [str(ENGINE), "--target", pub,
       "--start", "0x4000000000000000000000000000000000",
       "--size", "20000000000", "--segment", "20000000000",
       "--threads", "2", "--progress-interval", "0",
       "--table-binary", str(TABLE), "--glv-orbit", "--partner-table"]

CACHE.unlink(missing_ok=True)
for tag in ("cold (full validation)", "warm (cache hit)"):
    t0 = time.monotonic()
    proc = subprocess.run(CMD, text=True, capture_output=True, timeout=120)
    wall = time.monotonic() - t0
    mode = re.search(r"table binary validation: (.+)", proc.stdout)
    build = re.search(r"table build seconds: ([0-9.e+-]+)", proc.stdout)
    ok = "recovery correct: true" in proc.stdout
    print(f"{tag}: wall {wall:.2f} s, build {build.group(1) if build else '?'} s, "
          f"mode='{mode.group(1).strip() if mode else '?'}', recovered={ok}")
    if proc.returncode != 0:
        print(proc.stderr[-300:])
        sys.exit(1)
