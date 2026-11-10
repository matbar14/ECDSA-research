#!/usr/bin/env python3
"""Test base-centered rotations and zero-sum angle collapse.

Lambda orbit values satisfy a sum identity modulo p. This experiment separates
that algebraic identity from the chosen real-angle lift. It measures the wrap
integer, base-referenced centroid phase, and rotation-invariant enclosing-arc
width, then tests their relation to scalar bands on a holdout sample.
"""

import argparse
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import spearmanr

from w1 import (
    BASE_POINT, DomeInvariantAnalyzer, ENDOMORPHISM_BETA,
    FIELD_PRIME, GROUP_ORDER, EllipticCurvePoint,
)
from w1_circle_uniformity_experiment import point_from_scalar


def lambda_values(point, analyzer):
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(values)


def circular_shape(values):
    turns = sorted(value / FIELD_PRIME for value in values)
    gaps = [turns[1] - turns[0], turns[2] - turns[1], 1.0 - turns[2] + turns[0]]
    largest_gap = max(gaps)
    # Smallest connected arc containing all three vertices.
    return 1.0 - largest_gap, tuple(sorted(gaps))


def centroid_phase(values, reference=0):
    angles = [math.tau * ((value - reference) % FIELD_PRIME) / FIELD_PRIME for value in values]
    x = sum(math.cos(angle) for angle in angles)
    y = sum(math.sin(angle) for angle in angles)
    return (math.atan2(y, x) % math.tau) / math.tau


def row(offset, values, base_lambda, width):
    enclosing_arc, gaps = circular_shape(values)
    wrap_integer = sum(values) // FIELD_PRIME
    base_phase = centroid_phase(values, base_lambda)
    raw_phase = centroid_phase(values)
    # The three possible common rotations that make the lifted sum an integer
    # multiple of one turn differ by 1/3 turn. Their absolute phase is tested.
    zero_sum_rotations = tuple(((-wrap_integer / 3.0) + j / 3.0) % 1.0 for j in range(3))
    return {
        "offset": offset,
        "values": values,
        "enclosing_arc": enclosing_arc,
        "gaps": gaps,
        "wrap_integer": wrap_integer,
        "base_phase": base_phase,
        "raw_phase": raw_phase,
        "zero_sum_rotations": zero_sum_rotations,
        "width": width,
    }


def band_accuracy(rows, feature, width, bands, sectors):
    labels = [item["offset"] * bands // width for item in rows]
    values = []
    for item in rows:
        if feature == "wrap_integer":
            values.append(item[feature])
        elif feature == "enclosing_arc":
            values.append(min(sectors - 1, int(item[feature] * sectors)))
        else:
            values.append(min(sectors - 1, int(item[feature] * sectors)))
    indices = list(range(len(rows)))
    rng = random.Random(20261004)
    rng.shuffle(indices)
    split = len(rows) * 7 // 10
    train, test = indices[:split], indices[split:]
    counts = defaultdict(Counter)
    for index in train:
        counts[values[index]][labels[index]] += 1
    global_counts = Counter(labels[index] for index in train)
    fallback = global_counts.most_common(1)[0][0]
    predictor = {key: counter.most_common(1)[0][0] for key, counter in counts.items()}
    correct = sum(predictor.get(values[index], fallback) == labels[index] for index in test)
    rho, rho_p = spearmanr([item["offset"] / width for item in rows],
                           [item[feature] for item in rows])
    return {
        "bins": len(counts),
        "holdout_accuracy": correct / len(test),
        "random_baseline": 1 / bands,
        "spearman_rho": float(rho),
        "spearman_p_value": float(rho_p),
    }


def run(start, width, samples, mode, bands, sectors, seed):
    if start < 1 or width < 2 or start + width > GROUP_ORDER:
        raise ValueError("invalid interval")
    analyzer = DomeInvariantAnalyzer()
    base_values = lambda_values(BASE_POINT, analyzer)
    base_lambda = base_values[0]
    if mode == "consecutive":
        offsets = list(range(samples))
        effective_width = samples
    else:
        rng = random.Random(seed)
        offsets = sorted({rng.randrange(width) for _ in range(samples)})
        while len(offsets) < samples:
            offsets.append(rng.randrange(width))
        offsets = sorted(offsets[:samples])
        effective_width = width
    rows = []
    for offset in offsets:
        values = lambda_values(point_from_scalar(start + offset), analyzer)
        rows.append(row(offset, values, base_lambda, effective_width))
    wrap_counts = Counter(item["wrap_integer"] for item in rows)
    collapse_thresholds = {
        str(threshold): sum(item["enclosing_arc"] <= threshold for item in rows) / len(rows)
        for threshold in (0.05, 0.10, 0.25, 0.50)
    }
    report = {
        "start": start, "width": width, "effective_width": effective_width,
        "samples": len(rows), "mode": mode, "bands": bands, "sectors": sectors,
        "base_lambda_values": base_values,
        "base_lambda_reference": "Lambda(G), with phi(G) and phi^2(G) boundary values",
        "wrap_integer_counts": dict(wrap_counts),
        "collapse_fraction_by_arc_width": collapse_thresholds,
        "mean_enclosing_arc": sum(item["enclosing_arc"] for item in rows) / len(rows),
        "features": {
            feature: band_accuracy(rows, feature, effective_width, bands, sectors)
            for feature in ("wrap_integer", "enclosing_arc", "base_phase", "raw_phase")
        },
        "sample_rows": rows[:20],
    }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--width", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--samples", type=int, default=30000)
    parser.add_argument("--mode", choices=("random", "consecutive"), default="random")
    parser.add_argument("--bands", type=int, default=16)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20261005)
    parser.add_argument("--json", type=Path, default=Path("results/triangle_rotation_collapse.json"))
    args = parser.parse_args()
    report = run(args.start, args.width, args.samples, args.mode,
                 args.bands, args.sectors, args.seed)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"mode={args.mode}, samples={report['samples']}, effective width={report['effective_width']}")
    print(f"base Lambda orbit={report['base_lambda_values']}")
    print(f"wrap integer counts={report['wrap_integer_counts']}")
    print(f"collapse fractions={report['collapse_fraction_by_arc_width']}")
    print("feature          holdout/base  rho       rho p")
    for name, result in report["features"].items():
        print(f"{name:16s} {result['holdout_accuracy']:.3%}/{result['random_baseline']:.3%} "
              f"{result['spearman_rho']:+.5f} {result['spearman_p_value']:.3g}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
