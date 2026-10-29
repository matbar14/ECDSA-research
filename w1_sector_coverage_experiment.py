#!/usr/bin/env python3
"""Test whether an angular sector from w1.py predicts a scalar range.

For each invariant, train on 70% of distinct known scalars and evaluate on
the other 30%. The angular bins are exact integer partitions of F_p. The
scalar interval is split into equal numeric bands. A useful sector should
predict unseen bands or a short scalar interval with reliable coverage.
"""

import argparse
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import chi2_contingency

from w1 import BASE_POINT, DomeInvariantAnalyzer, FIELD_PRIME, GROUP_ORDER
from w1_circle_uniformity_experiment import NAMES, point_from_scalar, residues


def shortest_interval(values, count):
    ordered = sorted(values)
    width, lo, hi = min(
        (ordered[i + count - 1] - ordered[i], ordered[i], ordered[i + count - 1])
        for i in range(len(ordered) - count + 1)
    )
    return lo, hi, width


def sample_rows(start, count, width, mode, seed):
    analyzer = DomeInvariantAnalyzer()
    if mode == "consecutive":
        point = point_from_scalar(start)
        for offset in range(count):
            yield offset, residues(point, analyzer)
            point = point + BASE_POINT
    else:
        rng = random.Random(seed)
        offsets = set()
        while len(offsets) < count:
            offsets.add(rng.randrange(width))
        for offset in sorted(offsets):
            yield offset, residues(point_from_scalar(start + offset), analyzer)


def analyze(rows, width, sectors, bands, seed, fraction, pairs=False):
    count = len(rows)
    rng = random.Random(seed)
    indices = list(range(count))
    rng.shuffle(indices)
    split = count * 7 // 10
    training, test = indices[:split], indices[split:]
    domain = width
    results = {}
    features = [(name, (column,)) for column, name in enumerate(NAMES)]
    if pairs:
        features.extend((name, columns) for name, columns in (
            ("x3+lambda_0", (0, 4)),
            ("lambda_0+lambda_1", (4, 5)),
            ("radius+area", (1, 2)),
        ))
    for name, columns in features:
        angular = [
            tuple(row[column] * sectors // FIELD_PRIME for column in columns)
            for _, row in rows
        ]
        scalar_bands = [offset * bands // domain for offset, _ in rows]
        by_sector = defaultdict(list)
        band_counts = defaultdict(Counter)
        table = defaultdict(lambda: [0] * bands)
        for i in training:
            sector, band = angular[i], scalar_bands[i]
            by_sector[sector].append(rows[i][0])
            band_counts[sector][band] += 1
            table[sector][band] += 1
        global_bands = Counter(scalar_bands[i] for i in training)
        default_band = global_bands.most_common(1)[0][0]
        predicted_bands = {
            sector: counts.most_common(1)[0][0]
            for sector, counts in band_counts.items()
        }
        top_count = max(1, bands // 4)
        top_bands = {
            sector: {band for band, _ in counts.most_common(top_count)}
            for sector, counts in band_counts.items()
        }
        default_top = {band for band, _ in global_bands.most_common(top_count)}
        intervals = {
            sector: shortest_interval(offsets, max(1, round(len(offsets) * fraction)))
            for sector, offsets in by_sector.items()
        }
        top1_correct = sum(predicted_bands.get(angular[i], default_band) == scalar_bands[i] for i in test)
        top1_baseline = sum(scalar_bands[i] == default_band for i in test)
        topq_correct = sum(scalar_bands[i] in top_bands.get(angular[i], default_top) for i in test)
        topq_baseline = sum(scalar_bands[i] in default_top for i in test)
        interval_hits = sum(
            angular[i] in intervals and intervals[angular[i]][0] <= rows[i][0] <= intervals[angular[i]][1]
            for i in test
        )
        mean_interval_width = sum(
            (intervals[angular[i]][2] / domain) if angular[i] in intervals else 1.0
            for i in test
        ) / len(test)
        observed = list(table.values())
        if len(columns) == 1:
            chi_square, p_value, _, _ = chi2_contingency(observed, correction=False)
            cramers_v = (chi_square / (len(training) * min(len(observed) - 1, bands - 1))) ** 0.5
        else:
            # With S^2 pair cells, many expected cell counts are small;
            # asymptotic chi-square p-values are unreliable here.
            p_value = None
            cramers_v = None
        spans = [max(v) - min(v) + 1 for v in by_sector.values()]
        results[name] = {
            "occupied_train_sectors": len(by_sector),
            "mean_train_sector_span_fraction": sum(spans) / len(spans) / domain,
            "top1_holdout_accuracy": top1_correct / len(test),
            "top1_baseline_accuracy": top1_baseline / len(test),
            "top_quarter_holdout_accuracy": topq_correct / len(test),
            "top_quarter_baseline_accuracy": topq_baseline / len(test),
            "train_short_interval_fraction": fraction,
            "holdout_short_interval_coverage": interval_hits / len(test),
            "mean_short_interval_width_fraction": mean_interval_width,
            "train_band_independence_p_value": p_value,
            "train_cramers_v": cramers_v,
        }
    return len(training), len(test), results


def run(start, count, width, sectors, bands, mode, seed, fraction, pairs=False):
    if start < 1 or count < 100 or start + width > GROUP_ORDER or count > width:
        raise ValueError("require >=100 samples inside a valid group-order interval")
    if sectors < 2 or bands < 2 or count * .7 / (sectors * bands) < 5:
        raise ValueError("need >=5 expected training samples per angular-sector/scalar-band cell")
    if not 0 < fraction < 1:
        raise ValueError("interval fraction must lie between 0 and 1")
    rows = list(sample_rows(start, count, width, mode, seed))
    train_n, test_n, result = analyze(rows, count if mode == "consecutive" else width,
                                      sectors, bands, seed + 1, fraction, pairs)
    return {
        "start": start, "count": count, "width": width, "mode": mode,
        "sectors": sectors, "scalar_bands": bands, "seed": seed,
        "train_count": train_n, "holdout_count": test_n,
        "invariants": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda s: int(s, 0), default=1 << 134)
    parser.add_argument("--count", type=int, default=30000)
    parser.add_argument("--width", type=lambda s: int(s, 0), default=1 << 134)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--bands", type=int, default=16)
    parser.add_argument("--mode", choices=("consecutive", "random"), default="consecutive")
    parser.add_argument("--seed", type=int, default=20260926)
    parser.add_argument("--interval-fraction", type=float, default=0.8)
    parser.add_argument("--pairs", action="store_true", help="also test three pairs of sector coordinates")
    parser.add_argument("--json", type=Path, default=Path("w1_sector_coverage.json"))
    args = parser.parse_args()
    report = run(args.start, args.count, args.width, args.sectors, args.bands,
                 args.mode, args.seed, args.interval_fraction, args.pairs)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"mode={args.mode}, train={report['train_count']}, holdout={report['holdout_count']}, "
          f"angular sectors={args.sectors}, scalar bands={args.bands}")
    print("invariant   span    top1/base     top25/base    interval coverage/width    Cramer's V  p")
    for name, row in report["invariants"].items():
        print(f"{name:11s} {row['mean_train_sector_span_fraction']:6.1%}  "
              f"{row['top1_holdout_accuracy']:5.1%}/{row['top1_baseline_accuracy']:5.1%}  "
              f"{row['top_quarter_holdout_accuracy']:5.1%}/{row['top_quarter_baseline_accuracy']:5.1%}  "
              f"{row['holdout_short_interval_coverage']:5.1%}/{row['mean_short_interval_width_fraction']:5.1%}  "
              f"{row['train_cramers_v']:.4f}  {row['train_band_independence_p_value']:.3g}"
              if row['train_cramers_v'] is not None else
              f"{name:11s} {row['mean_train_sector_span_fraction']:6.1%}  "
              f"{row['top1_holdout_accuracy']:5.1%}/{row['top1_baseline_accuracy']:5.1%}  "
              f"{row['top_quarter_holdout_accuracy']:5.1%}/{row['top_quarter_baseline_accuracy']:5.1%}  "
              f"{row['holdout_short_interval_coverage']:5.1%}/{row['mean_short_interval_width_fraction']:5.1%}  --")
    print(f"report={args.json}")


if __name__ == "__main__":
    main()
