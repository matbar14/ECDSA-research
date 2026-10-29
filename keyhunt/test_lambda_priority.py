#!/usr/bin/env python3
"""Small end-to-end checks for Lambda ordering and automatic coverage."""

import re
import subprocess
import sys
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from w1 import DomeInvariantAnalyzer, EllipticCurvePoint  # noqa: E402


def public_key(scalar):
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def run(binary, scalar, *options):
    result = subprocess.run(
        [binary, "--target", public_key(scalar), *map(str, options)],
        text=True, capture_output=True, check=False,
    )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]{64})", result.stdout)
    recovered = int(match.group(1), 16) if match else None
    assert result.returncode == 0 and recovered == scalar, (result.returncode, recovered, scalar)
    return result.stdout


def check_lambda_dump(binary):
    scalar = 12345
    output = subprocess.check_output(
        [binary, "--target", public_key(scalar), "--start", "1", "--size", "1", "--lambda-dump"],
        text=True,
    )
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    y = numbers.y
    point = EllipticCurvePoint(numbers.x, y)
    expected = DomeInvariantAnalyzer()
    orbit = expected.analyze_endomorphism_orbit_silent(point)["points"]
    actual = re.findall(r"lambda[0-2]: 0x([0-9a-f]+)", output)
    assert actual == [f"{expected.compute_torus_invariant_lambda(q):064x}" for q in orbit]


def main():
    binary = str(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("eng_bsgs"))
    check_lambda_dump(binary)
    run(binary, 54321, "--start", "1", "--size", "100000", "--segment", "10000",
        "--lambda-priority", "--lambda-sectors", "8")
    output = run(binary, 1249, "--start", "1", "--size", "100", "--segment", "20",
                 "--cover-start", "1000", "--cover-size", "250",
                 "--lambda-priority", "--lambda-sectors", "4", "--glv-orbit")
    assert "automatic source chunk start: 0x4b0 size: 50" in output
    mixed = subprocess.run(
        [binary, "--target", public_key(10), "--start", "1", "--size", "10",
         "--shift", "1", "--cover-start", "1", "--cover-size", "10"],
        text=True, capture_output=True,
    )
    assert mixed.returncode == 2 and "cannot be combined" in mixed.stderr
    print("lambda priority checks: true")


if __name__ == "__main__":
    main()
