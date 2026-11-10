#!/usr/bin/env python3
"""Test Lambda triangles on a circle whose circumference is the search width.

For a search interval [start, start + width), the scalar coordinate is the
exact offset k-start on a circle of circumference ``width``. A field residue
Lambda in [0,p) is projected to that same coordinate system by
floor(Lambda * width / p). This projection is explicit and reversible only
up to its chosen precision; it does not assert a mathematical relation to k.
"""

import argparse
import json
import math
import random
from pathlib import Path

from scipy.stats import chi2_contingency

from w1 import BASE_POINT, DomeInvariantAnalyzer, FIELD_PRIME, GROUP_ORDER
from w1_circle_uniformity_experiment import point_from_scalar, residues


def sample_rows(start, width, samples, mode, seed):
    if mode == "consecutive":
        offsets = range(samples)
    else:
        rng = random.Random(seed)
        offsets = sorted({rng.randrange(width) for _ in range(samples)})
        while len(offsets) < samples:
            offsets.append(rng.randrange(width))
        offsets = offsets[:samples]
    analyzer = DomeInvariantAnalyzer()
    for offset in offsets:
        point = point_from_scalar(start + offset)
        values = residues(point, analyzer)[4:7]
        positions = tuple(value * width // FIELD_PRIME for value in values)
        yield offset, positions


def draw(rows, width, bands, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 4, figsize=(14, 11), constrained_layout=True)
    palette = plt.get_cmap("tab20", bands)
    selected = [rows[min(len(rows) - 1, i * len(rows) // 12)] for i in range(12)]
    for ax, (offset, positions) in zip(axes.flat, selected):
        for band in range(bands):
            angles = [math.tau * (band + j / 40) / bands for j in range(41)]
            ax.plot([math.cos(t) for t in angles], [math.sin(t) for t in angles],
                    color=palette(band), linewidth=6, alpha=.65)
        angles = [math.tau * position / width for position in positions]
        vertices = [(math.cos(angle), math.sin(angle)) for angle in angles]
        polygon = vertices + vertices[:1]
        ax.plot([x for x, _ in polygon], [y for _, y in polygon], color="#1d2c42", linewidth=1.4)
        ax.scatter([x for x, _ in vertices], [y for _, y in vertices], color="#e02626", s=32, zorder=3)
        scalar_angle = math.tau * offset / width
        ax.scatter([1.17 * math.cos(scalar_angle)], [1.17 * math.sin(scalar_angle)],
                   color=[palette(min(offset * bands // width, bands - 1))], marker="*", s=105, zorder=4)
        ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"offset {offset} / width, band {offset * bands // width}", fontsize=10)
    fig.suptitle("Lambda triangles on a circle whose circumference equals the search width\n"
                 "colored arcs = scalar bands; star = actual k offset; red points = projected Lambda", fontsize=13)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run(start, width, samples, sectors, bands, mode, seed, figure, report_path):
    if start < 1 or width < 2 or start + width > GROUP_ORDER:
        raise ValueError("invalid search interval")
    if samples < 100 or sectors < 2 or bands < 2:
        raise ValueError("samples, sectors, and bands are too small")
    rows = list(sample_rows(start, width, samples, mode, seed))
    matrix = [[0] * bands for _ in range(sectors)]
    same_vertices = 0
    any_same = 0
    for offset, positions in rows:
        band = offset * bands // width
        matching = 0
        for position in positions:
            sector = position * sectors // width
            matrix[sector][band] += 1
            matching += position * bands // width == band
        same_vertices += matching
        any_same += matching > 0
    occupied = [row for row in matrix if sum(row)]
    if len({band for row in occupied for band, value in enumerate(row) if value}) < 2:
        statistic, p_value = None, None
    else:
        statistic, p_value, _, _ = chi2_contingency(occupied, correction=False)
    report = {
        "start": start, "width": width, "samples": samples, "mode": mode,
        "circumference": width, "field_prime": FIELD_PRIME,
        "mapping": "position = floor(Lambda * width / field_prime)",
        "angular_sectors": sectors, "scalar_bands": bands,
        "same_band_vertex_fraction": same_vertices / (3 * len(rows)),
        "same_band_vertex_baseline": 1 / bands,
        "triangle_any_same_band_fraction": any_same / len(rows),
        "sector_band_chi_square_descriptive": statistic,
        "sector_band_p_value_descriptive": p_value,
        "matrix": matrix,
    }
    draw(rows, width, bands, figure)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda s: int(s, 0), default=1 << 134)
    parser.add_argument("--width", type=lambda s: int(s, 0), default=1 << 134)
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--bands", type=int, default=8)
    parser.add_argument("--mode", choices=("consecutive", "random"), default="random")
    parser.add_argument("--seed", type=int, default=20261002)
    parser.add_argument("--figure", type=Path, default=Path("w1_range_circle.png"))
    parser.add_argument("--json", type=Path, default=Path("w1_range_circle.json"))
    args = parser.parse_args()
    report = run(args.start, args.width, args.samples, args.sectors, args.bands,
                 args.mode, args.seed, args.figure, args.json)
    print(f"interval=[{args.start}, {args.start + args.width}), circumference={args.width}")
    print(f"mapping: {report['mapping']}")
    print(f"same-band vertices: {report['same_band_vertex_fraction']:.3%} "
          f"(random baseline {report['same_band_vertex_baseline']:.3%})")
    print(f"triangles with any same-band vertex: {report['triangle_any_same_band_fraction']:.3%}")
    print(f"figure={args.figure}; report={args.json}")


if __name__ == "__main__":
    main()
