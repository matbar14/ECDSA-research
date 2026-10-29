#!/usr/bin/env python3
"""Place each w1.py Lambda-orbit triangle on a circle of length 2^256-1.

The three Lambda values are calculated modulo secp256k1 field prime p.
Their canonical representatives are placed at integer positions on a separate
circle with circumference 2^256-1. A plot colors each triangle by its known
scalar band to inspect whether triangle position locates the band.
"""

import argparse
import json
import math
from pathlib import Path

from scipy.stats import chi2_contingency

from w1 import (
    BASE_POINT, DomeInvariantAnalyzer, EllipticCurvePoint,
    ENDOMORPHISM_BETA, FIELD_PRIME, GROUP_ORDER,
)
from w1_circle_uniformity_experiment import point_from_scalar


CIRCUMFERENCE = (1 << 256) - 1


def triangle(point, analyzer):
    values = []
    for _ in range(3):
        values.append(analyzer.compute_torus_invariant_lambda(point))
        point = EllipticCurvePoint(ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y)
    return tuple(values)


def circular_gap_signature(positions):
    ordered = sorted(positions)
    gaps = (ordered[1] - ordered[0], ordered[2] - ordered[1],
            CIRCUMFERENCE - ordered[2] + ordered[0])
    return tuple(sorted(gaps))


def draw(samples, bands, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 4, figsize=(14, 11), constrained_layout=True)
    palette = plt.get_cmap("tab20", bands)
    selected = [samples[min(len(samples) - 1, i * len(samples) // 12)] for i in range(12)]
    for ax, (offset, values) in zip(axes.flat, selected):
        band = offset * bands // len(samples)
        for scalar_band in range(bands):
            angles = [math.tau * (scalar_band + j / 40) / bands for j in range(41)]
            ax.plot([math.cos(t) for t in angles], [math.sin(t) for t in angles],
                    color=palette(scalar_band), linewidth=6, alpha=.65)
        angles = [math.tau * (position / CIRCUMFERENCE) for position in values]
        vertices = [(math.cos(angle), math.sin(angle)) for angle in angles]
        points = vertices + vertices[:1]
        ax.plot([p[0] for p in points], [p[1] for p in points],
                color="#1d2c42", linewidth=1.4)
        ax.scatter([p[0] for p in vertices], [p[1] for p in vertices],
                   color=["#e02626"], s=32, zorder=3)
        scalar_angle = math.tau * offset / len(samples)
        ax.scatter([1.17 * math.cos(scalar_angle)], [1.17 * math.sin(scalar_angle)],
                   color=[palette(band)], marker="*", s=105, zorder=4)
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"k offset {offset}, band {band}", fontsize=10)
    fig.suptitle("Individual w1 Lambda triangles | circle length 2^256−1\n"
                 "Colored arcs = equal parts of k range; star = k position; red dots = invariant vertices", fontsize=13)
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run(start, count, circle_sectors, bands, figure, json_output):
    if start < 1 or count < 2 or start + count > GROUP_ORDER:
        raise ValueError("invalid scalar interval")
    if circle_sectors < 2 or bands < 2 or count * 3 / (circle_sectors * bands) < 5:
        raise ValueError("need at least five expected vertices per sector/band cell")
    analyzer = DomeInvariantAnalyzer()
    point = point_from_scalar(start)
    samples = []
    matrix = [[0] * bands for _ in range(circle_sectors)]
    same_band_vertex_count = 0
    same_band_triangle_count = 0
    for offset in range(count):
        positions = triangle(point, analyzer)
        samples.append((offset, positions))
        band = offset * bands // count
        matching = sum(position * bands // CIRCUMFERENCE == band for position in positions)
        same_band_vertex_count += matching
        same_band_triangle_count += matching > 0
        for position in positions:
            sector = position * circle_sectors // CIRCUMFERENCE
            matrix[sector][band] += 1
        point = point + BASE_POINT
    occupied = [row for row in matrix if sum(row)]
    statistic, p_value, _, _ = chi2_contingency(occupied, correction=False)
    # The three orbit phases are dependent, so the vertex-table chi-square is
    # descriptive only; measure the same-band event at triangle level too.
    same_band_rate = same_band_vertex_count / (3 * count)
    independent_triangle_count = count
    triangle_rate = same_band_triangle_count / independent_triangle_count
    report = {
        "circumference": CIRCUMFERENCE,
        "field_prime": FIELD_PRIME,
        "start": start,
        "count": count,
        "angular_sectors": circle_sectors,
        "scalar_bands": bands,
        "occupied_angular_sectors": len(occupied),
        "same_k_band_vertex_fraction": same_band_rate,
        "triangles_with_any_vertex_in_same_k_band_fraction": triangle_rate,
        "random_vertex_baseline": 1 / bands,
        "mean_same_band_vertices_per_triangle": same_band_vertex_count / count,
        "vertex_count_by_sector_and_scalar_band": matrix,
        "sector_band_chi_square": statistic,
        "sector_band_p_value_descriptive_only": p_value,
        "sample_triangles": [
            {"offset": offset, "positions": list(positions),
             "circle_sectors": [p * circle_sectors // CIRCUMFERENCE for p in positions],
             "sorted_circular_gaps": list(circular_gap_signature(positions))}
            for offset, positions in samples[:min(20, count)]
        ],
    }
    draw(samples, bands, figure)
    json_output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=lambda s: int(s, 0), default=1 << 134)
    parser.add_argument("--count", type=int, default=10000)
    parser.add_argument("--sectors", type=int, default=32)
    parser.add_argument("--bands", type=int, default=8)
    parser.add_argument("--figure", type=Path, default=Path("w1_triangle_circle.png"))
    parser.add_argument("--json", type=Path, default=Path("w1_triangle_circle.json"))
    args = parser.parse_args()
    report = run(args.start, args.count, args.sectors, args.bands,
                 args.figure, args.json)
    print(f"circumference={CIRCUMFERENCE}")
    print(f"field_prime={FIELD_PRIME}")
    print(f"triangles={args.count}, occupied sectors={report['occupied_angular_sectors']}/{args.sectors}")
    print(f"vertex sector vs scalar band (descriptive; orbit vertices depend): "
          f"chi2={report['sector_band_chi_square']:.3f}, "
          f"p={report['sector_band_p_value_descriptive_only']:.6g}")
    print(f"vertex in same k band: {report['same_k_band_vertex_fraction']:.3%} "
          f"(baseline {report['random_vertex_baseline']:.3%})")
    print(f"triangle with any vertex in same k band: "
          f"{report['triangles_with_any_vertex_in_same_k_band_fraction']:.3%}")
    print(f"figure={args.figure}; report={args.json}")


if __name__ == "__main__":
    main()
