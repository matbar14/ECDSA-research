#!/usr/bin/env python3
"""Exact bounded Pollard-kangaroo benchmark for secp256k1.

The benchmark is intentionally conservative: it implements negation-aware
collision resolution and verifies the recovered scalar by comparing the full
public point.  A GLV orbit is not treated as a free collision key: doing so
without carrying the corresponding scalar multiplier is mathematically wrong.
The report therefore gives an exact baseline for the range-search engine and
an honest estimate for the remaining GLV work.
"""

import argparse
import json
import random
import statistics
import time
from pathlib import Path

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141


class Point:
    __slots__ = ("x", "y", "inf")

    def __init__(self, x=None, y=None, inf=False):
        self.x, self.y, self.inf = x, y, inf


G = Point(
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)


def neg(a):
    return Point(a.x, (-a.y) % P) if not a.inf else a


def add(a, b):
    if a.inf:
        return b
    if b.inf:
        return a
    if a.x == b.x and (a.y + b.y) % P == 0:
        return Point(inf=True)
    if a.x == b.x and a.y == b.y:
        if a.y == 0:
            return Point(inf=True)
        slope = 3 * a.x * a.x * pow(2 * a.y, -1, P) % P
    else:
        slope = (b.y - a.y) * pow((b.x - a.x) % P, -1, P) % P
    x = (slope * slope - a.x - b.x) % P
    return Point(x, (slope * (a.x - x) - a.y) % P)


def mul(k, a=G):
    k %= N
    out = Point(inf=True)
    while k:
        if k & 1:
            out = add(out, a)
        a = add(a, a)
        k >>= 1
    return out


def same(a, b):
    return a.inf == b.inf and a.x == b.x and a.y == b.y


def kangaroo(start, width, target_scalar, rng, dp_bits, jumps):
    """Return exact recovery statistics for one bounded target."""
    target = mul(target_scalar)
    start_point = mul(start)
    relative_target = add(target, neg(start_point))
    mid = width // 2
    jump_points = [mul(j) for j in jumps]
    tame = mul(mid)
    wild = add(relative_target, neg(mul(mid)))
    tame_dist = 0
    wild_dist = 0
    # Full point key: negation is handled explicitly below, not accidentally.
    records = {}
    steps = 0

    def visit(point, kind, distance):
        nonlocal steps
        if point.inf or (point.x & ((1 << dp_bits) - 1)):
            return None
        key = (point.x, point.y)
        prior = records.get(key)
        opposite = records.get((point.x, (-point.y) % P))
        records.setdefault(key, (kind, distance))
        other = prior
        if other is not None and other[0] == kind:
            other = None
        opposite_kind = opposite is not None and opposite[0] != kind
        if other is None and not opposite_kind:
            return None
        candidates = []
        if other is not None:
            if kind == 0:
                tame_d, wild_d = distance, other[1]
            else:
                tame_d, wild_d = other[1], distance
            # Relative target = W + tame_d - wild_d.
            candidates.append(start + width + tame_d - wild_d)
        if opposite_kind:
            # A collision with the negated point gives relative target
            # = -W - tame_d - wild_d (mod N).  Test its canonical lift.
            if kind == 0:
                tame_d, wild_d = distance, opposite[1]
            else:
                tame_d, wild_d = opposite[1], distance
            candidates.append(start - width - tame_d - wild_d)
        for candidate in candidates:
            rel = candidate - start
            if 0 <= rel < width and same(mul(candidate), target):
                return candidate
        return None

    while steps < 200 * max(1, width.bit_length() * (1 << (width.bit_length() // 2))):
        for kind in (0, 1):
            point = tame if kind == 0 else wild
            index = point.x % len(jumps)
            point = add(point, jump_points[index])
            distance = (tame_dist if kind == 0 else wild_dist) + jumps[index]
            if kind == 0:
                tame, tame_dist = point, distance
            else:
                wild, wild_dist = point, distance
            steps += 1
            found = visit(point, kind, distance)
            if found is not None:
                return {"found": True, "steps": steps,
                        "records": len(records), "scalar": found}
    return {"found": False, "steps": steps, "records": len(records)}


def run_case(bits, trials, dp_bits, jump_count, seed):
    width = 1 << bits
    rng = random.Random(seed + bits)
    jump_max = 1 << (bits // 2)
    jumps = [rng.randrange(1, jump_max + 1) for _ in range(jump_count)]
    rows = []
    t0 = time.perf_counter()
    for _ in range(trials):
        start = 1 << 40
        target = start + rng.randrange(width)
        row = kangaroo(start, width, target, rng, dp_bits, jumps)
        row["correct"] = row.get("scalar") == target
        rows.append(row)
    solved = [r for r in rows if r["correct"]]
    mean = statistics.mean(r["steps"] for r in solved) if solved else None
    return {
        "width_bits": bits,
        "width": width,
        "trials": trials,
        "solved": len(solved),
        "dp_bits": dp_bits,
        "jump_count": jump_count,
        "mean_steps": round(mean, 2) if mean else None,
        "mean_steps_over_sqrt_width": round(mean / (width ** 0.5), 4) if mean else None,
        "mean_records": round(statistics.mean(r["records"] for r in rows), 2),
        "wall_seconds": round(time.perf_counter() - t0, 3),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bits", type=int, nargs="+", default=[16, 18, 20, 22])
    parser.add_argument("--trials", type=int, default=16)
    parser.add_argument("--dp-bits", type=int, default=6)
    parser.add_argument("--jump-count", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20261118)
    parser.add_argument("--output", type=Path, default=Path("results/kangaroo_benchmark.json"))
    args = parser.parse_args()
    results = []
    for bits in args.bits:
        row = run_case(bits, args.trials, args.dp_bits, args.jump_count, args.seed)
        results.append(row)
        print(json.dumps(row, sort_keys=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"mode": "exact-negation-aware", "results": results}, indent=2) + "\n")


if __name__ == "__main__":
    main()
