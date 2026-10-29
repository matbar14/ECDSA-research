#!/usr/bin/env python3
"""Inspect exact decimal digits and binary bits of Lambda(kG) / p.

For a private scalar k, let a = Lambda(kG) and define the normalized phase

    r_k = a / p,  0 <= r_k < 1.

The decimal and binary expansions of r_k are computed with integer arithmetic:

    decimal digit i = floor(10^i * a / p) mod 10,
    binary bit i    = floor(2^i  * a / p) mod 2.

No floating-point angle and no decimal approximation of pi is used.  The
script reports digit/bit frequencies and scans simple alignments between
phase bits and low bits of the scalar offset.  It is an exploratory screening
test: a deviation in a finite sample is a candidate for a follow-up test, not
a proof of a relation to the private scalar.
"""

import argparse
import math
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import DomeInvariantAnalyzer, EllipticCurvePoint, FIELD_PRIME, GROUP_ORDER


RANGE_START = 1 << 255
SAMPLE_COUNT = 100000
DECIMAL_POSITIONS = 12
BINARY_POSITIONS = 24
OFFSET_BITS = 16
PROGRESS_INTERVAL = 10000


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def phase_decimal_digit(lambda_value, position):
    """Return decimal digit at one-based position after the point of Lambda/p."""
    return (lambda_value * (10 ** position) // FIELD_PRIME) % 10


def phase_binary_bit(lambda_value, position):
    """Return binary bit at one-based position after the point of Lambda/p."""
    return (lambda_value * (1 << position) // FIELD_PRIME) % 2


def run_experiment(
    range_start=RANGE_START,
    sample_count=SAMPLE_COUNT,
    decimal_positions=DECIMAL_POSITIONS,
    binary_positions=BINARY_POSITIONS,
    offset_bits=OFFSET_BITS,
):
    if range_start < 1 or sample_count <= 0 or range_start + sample_count > GROUP_ORDER:
        raise ValueError("sampled range must be contained in the canonical scalar interval")
    if min(decimal_positions, binary_positions, offset_bits) <= 0:
        raise ValueError("all position counts must be positive")

    analyzer = DomeInvariantAnalyzer()
    decimal_counts = [[0] * 10 for _ in range(decimal_positions)]
    binary_ones = [0] * binary_positions
    bit_agreements = [[0] * offset_bits for _ in range(binary_positions)]
    decimal_powers = [10 ** position for position in range(1, decimal_positions + 1)]
    binary_powers = [1 << position for position in range(1, binary_positions + 1)]
    examples = []
    start_time = time.perf_counter()

    for offset in range(sample_count):
        scalar = range_start + offset
        lambda_value = analyzer.compute_torus_invariant_lambda(point_from_scalar(scalar))
        if lambda_value is None:
            raise RuntimeError("Lambda is undefined for a sampled scalar")

        decimal_digits = []
        for index, power in enumerate(decimal_powers):
            digit = (lambda_value * power // FIELD_PRIME) % 10
            decimal_counts[index][digit] += 1
            decimal_digits.append(str(digit))

        binary_digits = []
        for index, power in enumerate(binary_powers):
            bit = (lambda_value * power // FIELD_PRIME) % 2
            binary_ones[index] += bit
            binary_digits.append(str(bit))
            for offset_bit in range(offset_bits):
                bit_agreements[index][offset_bit] += bit == ((offset >> offset_bit) & 1)

        if len(examples) < 3:
            examples.append(
                {
                    "offset": offset,
                    "lambda": lambda_value,
                    "decimal_prefix": "".join(decimal_digits),
                    "binary_prefix": "".join(binary_digits),
                }
            )
        if (offset + 1) % PROGRESS_INTERVAL == 0:
            print(f"  Completed: {offset + 1}/{sample_count}")

    decimal_results = []
    for position, counts in enumerate(decimal_counts, start=1):
        expected = sample_count / 10
        chi_square = math.fsum((count - expected) ** 2 / expected for count in counts)
        decimal_results.append(
            {
                "position": position,
                "counts": tuple(counts),
                "chi_square": chi_square,
            }
        )

    binary_results = []
    for index, one_count in enumerate(binary_ones):
        agreements = bit_agreements[index]
        deviations = [abs(agreement / sample_count - 0.5) for agreement in agreements]
        strongest_offset_bit = max(range(offset_bits), key=lambda bit: deviations[bit])
        binary_results.append(
            {
                "position": index + 1,
                "one_fraction": one_count / sample_count,
                "same_position_agreement": (
                    agreements[index] / sample_count if index < offset_bits else None
                ),
                "strongest_offset_bit": strongest_offset_bit,
                "strongest_agreement": agreements[strongest_offset_bit] / sample_count,
            }
        )

    return {
        "range_start": range_start,
        "range_end": range_start + sample_count - 1,
        "sample_count": sample_count,
        "decimal_results": decimal_results,
        "binary_results": binary_results,
        "examples": examples,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Exact phase digit and bit experiment")
    print(f"  Scalar range start: 0x{report['range_start']:064x}")
    print(f"  Scalar range end:   0x{report['range_end']:064x}")
    print(f"  Samples: {report['sample_count']}")
    print("  First exact normalized-phase prefixes:")
    for example in report["examples"]:
        print(
            f"    offset={example['offset']}"
            f"  Lambda=0x{example['lambda']:064x}"
            f"  r decimal=.{example['decimal_prefix']}"
            f"  r binary=.{example['binary_prefix']}"
        )

    print("  Decimal positions: digit count range and chi-square (9 degrees of freedom)")
    for result in report["decimal_results"]:
        print(
            f"    {result['position']:2d}: {min(result['counts'])}..{max(result['counts'])}"
            f"  chi2={result['chi_square']:.6f}"
        )

    print("  Binary positions: fraction of ones and offset-bit agreement")
    for result in report["binary_results"]:
        same_position = result["same_position_agreement"]
        same_text = f"{same_position:.4%}" if same_position is not None else "n/a"
        print(
            f"    {result['position']:2d}: ones={result['one_fraction']:.4%}"
            f"  same offset bit={same_text}"
            f"  strongest offset bit {result['strongest_offset_bit']}="
            f"{result['strongest_agreement']:.4%}"
        )
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=SAMPLE_COUNT)
    parser.add_argument("--range-start", type=lambda value: int(value, 0), default=RANGE_START)
    parser.add_argument("--decimal-positions", type=int, default=DECIMAL_POSITIONS)
    parser.add_argument("--binary-positions", type=int, default=BINARY_POSITIONS)
    parser.add_argument("--offset-bits", type=int, default=OFFSET_BITS)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    print_report(
        run_experiment(
            range_start=arguments.range_start,
            sample_count=arguments.samples,
            decimal_positions=arguments.decimal_positions,
            binary_positions=arguments.binary_positions,
            offset_bits=arguments.offset_bits,
        )
    )


if __name__ == "__main__":
    main()
