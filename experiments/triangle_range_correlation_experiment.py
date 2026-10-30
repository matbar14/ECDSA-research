#!/usr/bin/env python3
"""Test whether the complete Lambda triangle locates k in a search window.

The search interval is [start, start + width). Lambda residues are projected
to the same circle by floor(Lambda * width / p). Several triangle features are
quantized, fitted on 70% of the samples, and evaluated on an independent 30%
holdout. The experiment measures correlation only; it never excludes keys.
"""

import argparse
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path

from w1 import ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER, DomeInvariantAnalyzer, EllipticCurvePoint
from w1_circle_uniformity_experiment import point_from_scalar


def lambda_triangle(point, analyzer, width):
    values = []
    for _ in range(3):
        value = analyzer.compute_torus_invariant_lambda(point)
        values.append(value * width // FIELD_PRIME)
        point = EllipticCurvePoint(
            ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
        )
    return tuple(values)


def circular_gaps(positions, width):
    ordered = sorted(positions)
    return (
        ordered[1] - ordered[0],
        ordered[2] - ordered[1],
        width - ordered[2] + ordered[0],
    )


def features(positions, width, sectors):
    q = lambda value: value * sectors // width
    ordered = tuple(q(value) for value in positions)
    sorted_vertices = tuple(sorted(ordered))
    gaps = circular_gaps(positions, width)
    sorted_gaps = tuple(sorted(q(value) for value in gaps))

    angles = [math.tau * value / width for value in positions]
    center_x = sum(math.cos(angle) for angle in angles)
    center_y = sum(math.sin(angle) for angle in angles)
    center_angle = math.atan2(center_y, center_x) % math.tau
    center_sector = int(center_angle / math.tau * sectors) % sectors
    resultant_bucket = min(sectors, int(math.hypot(center_x, center_y) / 3 * sectors))

    return {
        "ordered_vertices": ordered,
        "sorted_vertices": sorted_vertices,
        "sorted_gaps": sorted_gaps,
        "center": (center_sector, resultant_bucket),
    }


def shortest_span(values, fraction, width):
    ordered = sorted(values)
    needed = max(1, math.ceil(len(ordered) * fraction))
    if needed >= len(ordered):
        return 0.0
    best = min(ordered[i + needed - 1] - ordered[i] for i in range(len(ordered) - needed + 1))
    return best / width


def entropy(counter, total):
    return -sum((count / total) * math.log2(count / total)
                for count in counter.values() if count)


def evaluate(rows, width, bands, sectors, seed, interval_fraction):
    rng = random.Random(seed)
    indices = list(range(len(rows)))
    rng.shuffle(indices)
    split = len(rows) * 7 // 10
    training, holdout = indices[:split], indices[split:]
    labels = [offset * bands // width for offset, _ in rows]
    feature_names = ("ordered_vertices", "sorted_vertices", "sorted_gaps", "center")
    output = {}

    global_counts = Counter(labels[i] for i in training)
    default_band = global_counts.most_common(1)[0][0]
    baseline_top = {band for band, _ in global_counts.most_common(max(1, bands // 4))}
    for name in feature_names:
        buckets = defaultdict(list)
        bucket_bands = defaultdict(Counter)
        for index in training:
            key = features(rows[index][1], width, sectors)[name]
            buckets[key].append(rows[index][0])
            bucket_bands[key][labels[index]] += 1
        predictors = {
            key: counts.most_common(1)[0][0]
            for key, counts in bucket_bands.items()
        }
        top_count = max(1, bands // 4)
        top_predictors = {
            key: {band for band, _ in counts.most_common(top_count)}
            for key, counts in bucket_bands.items()
        }
        top1 = sum(predictors.get(
            features(rows[i][1], width, sectors)[name], default_band) == labels[i]
                   for i in holdout) / len(holdout)
        topq = sum(labels[i] in top_predictors.get(
            features(rows[i][1], width, sectors)[name], baseline_top)
                   for i in holdout) / len(holdout)
        spans = [shortest_span(offsets, interval_fraction, width)
                 for offsets in buckets.values()]
        mean_span = sum(spans) / len(spans) if spans else 1.0
        holdout_coverage = 0
        holdout_span = 0.0
        for i in holdout:
            key = features(rows[i][1], width, sectors)[name]
            offsets = buckets.get(key)
            if not offsets:
                holdout_span += 1.0
                continue
            ordered = sorted(offsets)
            needed = max(1, math.ceil(len(ordered) * interval_fraction))
            candidates = [
                (ordered[j + needed - 1] - ordered[j], ordered[j], ordered[j + needed - 1])
                for j in range(len(ordered) - needed + 1)
            ]
            _, low, high = min(candidates)
            holdout_coverage += low <= rows[i][0] <= high
            holdout_span += (high - low) / width
        band_entropy = entropy(global_counts, len(training))
        conditional_entropy = 0.0
        for counts in bucket_bands.values():
            total = sum(counts.values())
            conditional_entropy += total / len(training) * entropy(counts, total)
        output[name] = {
            "training_buckets": len(buckets),
            "mean_bucket_members": len(training) / len(buckets),
            "top1_holdout_accuracy": top1,
            "top1_random_baseline": 1 / bands,
            "top_quarter_holdout_accuracy": topq,
            "top_quarter_random_baseline": len(baseline_top) / bands,
            "mean_training_short_span_fraction": mean_span,
            "holdout_short_span_coverage": holdout_coverage / len(holdout),
            "holdout_short_span_fraction": holdout_span / len(holdout),
            "normalized_mutual_information": (
                (band_entropy - conditional_entropy) / band_entropy
                if band_entropy else 0.0
            ),
        }
    return len(training), len(holdout), output


def run(start, width, samples, bands, sectors, mode, seed, interval_fraction):
    if start < 1 or width < 2 or start + width > GROUP_ORDER:
        raise ValueError("invalid search interval")
    if samples < 100 or samples > width:
        raise ValueError("samples must be between 100 and width")
    analyzer = DomeInvariantAnalyzer()
    if mode == "consecutive":
        offsets = list(range(samples))
    else:
        rng = random.Random(seed)
        offsets = sorted(rng.sample(range(width), samples)) if width < (1 << 63) else sorted({rng.randrange(width) for _ in range(samples)})
        while len(offsets) < samples:
            candidate = rng.randrange(width)
            if candidate not in offsets:
                offsets.append(candidate)
        offsets.sort()
    rows = []
    for offset in offsets:
        rows.append((offset, lambda_triangle(point_from_scalar(start + offset), analyzer, width)))
    train, holdout, result = evaluate(rows, width if mode == "random" else samples,
                                      bands, sectors, seed + 1, interval_fraction)
    return {
        "start": start, "width": width, "samples": samples, "mode": mode,
        "bands": bands, "sectors": sectors, "train_count": train,
        "holdout_count": holdout, "seed": seed, "features": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--width", type=lambda value: int(value, 0), default=1 << 134)
    parser.add_argument("--samples", type=int, default=30000)
    parser.add_argument("--bands", type=int, default=16)
    parser.add_argument("--sectors", type=int, default=8)
    parser.add_argument("--mode", choices=("random", "consecutive"), default="random")
    parser.add_argument("--seed", type=int, default=20261003)
    parser.add_argument("--interval-fraction", type=float, default=0.8)
    parser.add_argument("--json", type=Path, default=Path("results/triangle_range_correlation.json"))
    args = parser.parse_args()
    report = run(args.start, args.width, args.samples, args.bands, args.sectors,
                 args.mode, args.seed, args.interval_fraction)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"mode={args.mode}, interval width={args.width}, samples={args.samples}, "
          f"train={report['train_count']}, holdout={report['holdout_count']}")
    print("feature             buckets  top1/base  top25/base  NMI     span coverage/width")
    for name, row in report["features"].items():
        print(f"{name:19s} {row['training_buckets']:7d}  "
              f"{row['top1_holdout_accuracy']:5.1%}/{row['top1_random_baseline']:5.1%}  "
              f"{row['top_quarter_holdout_accuracy']:5.1%}/{row['top_quarter_random_baseline']:5.1%}  "
              f"{row['normalized_mutual_information']:.4f}  "
              f"{row['holdout_short_span_coverage']:5.1%}/{row['holdout_short_span_fraction']:5.1%}")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
