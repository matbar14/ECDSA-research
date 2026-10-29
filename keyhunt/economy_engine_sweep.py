#!/usr/bin/env python3
"""Engine-level economy sweep for the (d01,d12) router on a real 3-row bucket.

For each trial the target scalar is start + i0*m + s with s drawn from the
orbit scalars {±lambda^a * j3} of the injected table orbit, so the giant lands
on the genuine 3-row K3 bucket (rows j3, lambda*j3, lambda^2*j3). With
--partner-table every enumerated branch before the match is a full EC check,
which makes "exact candidate checks" directly comparable to the branch
enumeration measured by phase_signature_economy_experiment.py.

Modes compared:
  baseline6   --glv-orbit --partner-table                (six branches, row 1)
  signature   --glv-orbit --partner-table --phase-signature
  canonical   --glv-orbit --partner-table --phase-signature --canonical-k3
"""

import random
import re
import subprocess
import tempfile
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

LAMBDA = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
LAMBDA2 = 0xAC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

START = 1 << 134
SIZE = 1 << 20
J3 = 777
TRIALS = 24


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(binary, target, extra):
    cmd = [binary, "--target", target, "--start", hex(START),
           "--size", str(SIZE), "--segment", str(SIZE), "--threads", "2"] + extra
    result = subprocess.run(cmd, text=True, capture_output=True)
    out = result.stdout

    def counter(name):
        m = re.search(rf"{re.escape(name)}: (\d+)", out)
        return int(m.group(1)) if m else None

    m = re.search(r"recovered private key: 0x([0-9a-fA-F]+)", out)
    recovered = int(m.group(1), 16) if m else None
    return recovered, counter("exact candidate checks"), \
        counter("orbit factors scheduled")


def main():
    rng = random.Random(20261030)
    with tempfile.TemporaryDirectory() as tmpdir:
        rows = [J3, LAMBDA * J3 % N, LAMBDA2 * J3 % N]
        rows += [rng.randrange(1, 1 << 100) for _ in range(6)]
        rows_path = Path(tmpdir) / "rows.txt"
        rows_path.write_text(
            "\n".join([hex(r) for r in rows]) + "\n")
        modes = {
            "baseline6": ["--glv-orbit", "--partner-table",
                          "--table-scalars", str(rows_path)],
            "signature": ["--glv-orbit", "--partner-table", "--phase-signature",
                          "--table-scalars", str(rows_path)],
            "canonical": ["--glv-orbit", "--partner-table", "--phase-signature",
                          "--canonical-k3", "--table-scalars", str(rows_path)],
        }
        stats = {name: {"exact": 0, "factors": 0, "ok": 0}
                 for name in modes}
        trials = 0
        for trial in range(TRIALS):
            a = rng.randrange(3)
            sign = rng.choice((1, -1))
            s = sign * (J3, LAMBDA * J3 % N, LAMBDA2 * J3 % N)[a] % N
            if s == 0:
                continue
            trials += 1
            i0 = rng.randrange(50, 400)
            for name, extra in modes.items():
                m_rows = 7 if name == "canonical" else 9
                target = (START + i0 * m_rows + s) % N
                rec, exact, factors = run(args_binary(), pubkey_hex(target), extra)
                if rec == target:
                    stats[name]["ok"] += 1
                    stats[name]["exact"] += exact or 0
                    stats[name]["factors"] += factors or 0

        print(f"Движковая экономика маршрутизации на 3-строчном bucket "
              f"({trials} целей, {len(modes)} режима):")
        base = stats["baseline6"]
        print(f"  baseline6 : recovery {base['ok']}/{trials}, "
              f"точных EC-проверок на hit {base['exact']/max(1,base['ok']):.3f}, "
              f"итераций ветвей {base['factors']/max(1,base['ok']):.3f}")
        sig = stats["signature"]
        print(f"  signature : recovery {sig['ok']}/{trials}, "
              f"точных EC-проверок на hit {sig['exact']/max(1,sig['ok']):.3f}, "
              f"итераций ветвей {sig['factors']/max(1,sig['ok']):.3f}")
        can = stats["canonical"]
        print(f"  canonical : recovery {can['ok']}/{trials}, "
              f"точных EC-проверок на hit {can['exact']/max(1,can['ok']):.3f}, "
              f"итераций ветвей {can['factors']/max(1,can['ok']):.3f}")
        b, g = base["exact"]/max(1, base["ok"]), sig["exact"]/max(1, sig["ok"])
        if g > 0:
            print(f"  экономия точных проверок baseline/signature: {b/g:.3f}x")


def args_binary():
    return "./eng_bsgs"


if __name__ == "__main__":
    main()
