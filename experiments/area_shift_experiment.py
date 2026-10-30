#!/usr/bin/env python3
"""Measure whether Lambda-triangle areas shift across scalar-height bands.

The experiment partitions the exact 256-bit scalar interval [2^255, n) into
equal-width bands.  It samples each band uniformly, computes the Lambda
triangle area S(kG), and compares conditional area distributions.  The report
uses means, standard errors, a Pearson correlation with normalized scalar
height, and the two-sample Kolmogorov-Smirnov distance of each band against
all remaining bands.

This tests a distributional displacement of area rather than attempting to
predict individual bits of k.
"""

import math
import random
import time

from cryptography.hazmat.primitives.asymmetric import ec

from lambda_circle_triangle_experiment import (
    GROUP_ORDER,
    lambda_orbit,
    triangle_area,
    unit_circle_vertices,
)
from w1 import DomeInvariantAnalyzer, EllipticCurvePoint


BAND_COUNT = 8
SAMPLES_PER_BAND = 2000
SEED = 20260927
SVG_OUTPUT = "figures/area_shift_by_scalar_band.svg"


def point_from_scalar(scalar):
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def area_from_scalar(analyzer, scalar):
    values = lambda_orbit(analyzer, point_from_scalar(scalar))
    if values is None:
        return None
    return triangle_area(unit_circle_vertices(values))


def population_standard_deviation(values, mean):
    return math.sqrt(math.fsum((value - mean) ** 2 for value in values) / len(values))


def pearson_correlation(left, right):
    left_mean = math.fsum(left) / len(left)
    right_mean = math.fsum(right) / len(right)
    numerator = math.fsum(
        (left_value - left_mean) * (right_value - right_mean)
        for left_value, right_value in zip(left, right)
    )
    left_norm = math.sqrt(math.fsum((value - left_mean) ** 2 for value in left))
    right_norm = math.sqrt(math.fsum((value - right_mean) ** 2 for value in right))
    return numerator / (left_norm * right_norm) if left_norm and right_norm else 0.0


def ks_distance(first_values, second_values):
    """Return the two-sample empirical Kolmogorov-Smirnov distance."""
    first = sorted(first_values)
    second = sorted(second_values)
    first_index = 0
    second_index = 0
    maximum = 0.0
    while first_index < len(first) and second_index < len(second):
        if first[first_index] <= second[second_index]:
            current = first[first_index]
            while first_index < len(first) and first[first_index] == current:
                first_index += 1
        else:
            current = second[second_index]
            while second_index < len(second) and second[second_index] == current:
                second_index += 1
        maximum = max(maximum, abs(first_index / len(first) - second_index / len(second)))
    return maximum


def run_experiment(
    band_count=BAND_COUNT,
    samples_per_band=SAMPLES_PER_BAND,
    seed=SEED,
):
    if not isinstance(band_count, int) or band_count < 2:
        raise ValueError("band_count must be an integer of at least two")
    if not isinstance(samples_per_band, int) or samples_per_band <= 0:
        raise ValueError("samples_per_band must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    lower = 1 << 255
    width = (GROUP_ORDER - lower) // band_count
    bands = [[] for _ in range(band_count)]
    normalized_heights = []
    all_areas = []
    undefined_count = 0
    start_time = time.perf_counter()

    for band_index in range(band_count):
        band_lower = lower + band_index * width
        band_upper = GROUP_ORDER if band_index == band_count - 1 else band_lower + width
        for _ in range(samples_per_band):
            scalar = generator.randrange(band_lower, band_upper)
            area = area_from_scalar(analyzer, scalar)
            if area is None:
                undefined_count += 1
                continue
            bands[band_index].append(area)
            normalized_heights.append((scalar - lower) / (GROUP_ORDER - lower))
            all_areas.append(area)
        print(f"  Completed band: {band_index + 1}/{band_count}")

    global_mean = math.fsum(all_areas) / len(all_areas)
    global_stddev = population_standard_deviation(all_areas, global_mean)
    band_reports = []
    for band_index, values in enumerate(bands):
        mean = math.fsum(values) / len(values)
        stddev = population_standard_deviation(values, mean)
        other_values = [
            area for other_index, other_band in enumerate(bands)
            if other_index != band_index for area in other_band
        ]
        band_reports.append(
            {
                "index": band_index,
                "count": len(values),
                "mean": mean,
                "stddev": stddev,
                "standard_error": stddev / math.sqrt(len(values)),
                "mean_shift": mean - global_mean,
                "ks_distance": ks_distance(values, other_values),
            }
        )

    return {
        "band_count": band_count,
        "samples_per_band": samples_per_band,
        "seed": seed,
        "processed": len(all_areas),
        "undefined_count": undefined_count,
        "global_mean": global_mean,
        "global_stddev": global_stddev,
        "height_area_correlation": pearson_correlation(normalized_heights, all_areas),
        "bands": band_reports,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def write_svg(report, output_path=SVG_OUTPUT):
    """Render conditional means and standard errors as a compact static plot."""
    width = 820
    height = 350
    margin_left = 72
    margin_right = 34
    margin_top = 42
    margin_bottom = 62
    plot_width = width - margin_left - margin_right
    plot_height = height - margin_top - margin_bottom
    means = [band["mean"] for band in report["bands"]]
    errors = [band["standard_error"] for band in report["bands"]]
    lower = min(mean - error for mean, error in zip(means, errors))
    upper = max(mean + error for mean, error in zip(means, errors))
    padding = max((upper - lower) * 0.15, 1e-5)
    lower -= padding
    upper += padding

    def x_position(index):
        return margin_left + (index + 0.5) * plot_width / report["band_count"]

    def y_position(value):
        return margin_top + (upper - value) * plot_height / (upper - lower)

    global_y = y_position(report["global_mean"])
    marks = []
    labels = []
    for band in report["bands"]:
        x = x_position(band["index"])
        upper_y = y_position(band["mean"] + band["standard_error"])
        lower_y = y_position(band["mean"] - band["standard_error"])
        mean_y = y_position(band["mean"])
        marks.append(
            f'<line x1="{x:.2f}" y1="{upper_y:.2f}" x2="{x:.2f}" y2="{lower_y:.2f}" class="error" />'
            f'<line x1="{x - 5:.2f}" y1="{upper_y:.2f}" x2="{x + 5:.2f}" y2="{upper_y:.2f}" class="error" />'
            f'<line x1="{x - 5:.2f}" y1="{lower_y:.2f}" x2="{x + 5:.2f}" y2="{lower_y:.2f}" class="error" />'
            f'<circle cx="{x:.2f}" cy="{mean_y:.2f}" r="4" class="mean" />'
        )
        labels.append(
            f'<text x="{x:.2f}" y="{height - 30}" text-anchor="middle" class="label">{band["index"]}</text>'
        )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">
  <title id="title">Средняя площадь Lambda-треугольника по полосам скаляра</title>
  <desc id="description">Точки показывают условные средние площади, вертикальные отрезки -- стандартные ошибки. Пунктирная линия -- среднее по всей выборке.</desc>
  <style>
    .title {{ font: 600 16px sans-serif; fill: #18212b; }}
    .label {{ font: 12px sans-serif; fill: #425466; }}
    .axis {{ stroke: #607586; stroke-width: 1; }}
    .reference {{ stroke: #be185d; stroke-width: 1.5; stroke-dasharray: 5 4; }}
    .error {{ stroke: #0f766e; stroke-width: 1.5; }}
    .mean {{ fill: #0f766e; }}
  </style>
  <rect width="{width}" height="{height}" fill="#ffffff" />
  <text x="{margin_left}" y="24" class="title">Средняя площадь по полосам нормированной высоты k</text>
  <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height - margin_bottom}" class="axis" />
  <line x1="{margin_left}" y1="{height - margin_bottom}" x2="{width - margin_right}" y2="{height - margin_bottom}" class="axis" />
  <line x1="{margin_left}" y1="{global_y:.2f}" x2="{width - margin_right}" y2="{global_y:.2f}" class="reference" />
  <text x="{margin_left - 10}" y="{global_y - 6:.2f}" text-anchor="end" class="label">{report['global_mean']:.6f}</text>
  {''.join(marks)}
  {''.join(labels)}
  <text x="{width / 2:.2f}" y="{height - 8}" text-anchor="middle" class="label">полоса k</text>
  <text x="14" y="{height / 2:.2f}" transform="rotate(-90 14 {height / 2:.2f})" text-anchor="middle" class="label">площадь S(kG)</text>
</svg>
'''
    with open(output_path, "w", encoding="utf-8") as output_file:
        output_file.write(svg)


def print_report(report):
    print("Area-shift experiment by scalar-height band")
    print(f"  256-bit scalar bands: {report['band_count']}")
    print(f"  Samples per band: {report['samples_per_band']}")
    print(f"  Processed: {report['processed']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Undefined areas: {report['undefined_count']}")
    print(f"  Global area mean: {report['global_mean']:.12f}")
    print(f"  Global area standard deviation: {report['global_stddev']:.12f}")
    print(f"  Pearson corr(normalized k, area): {report['height_area_correlation']:.8f}")
    print("  band  mean area      mean shift     std. error     KS vs. other bands")
    for band in report["bands"]:
        print(
            f"  {band['index']:4d}  {band['mean']:.12f}"
            f"  {band['mean_shift']:+.12f}"
            f"  {band['standard_error']:.12f}"
            f"  {band['ks_distance']:.8f}"
        )
    print(f"  SVG figure: {SVG_OUTPUT}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    report = run_experiment()
    write_svg(report)
    print_report(report)


if __name__ == "__main__":
    main()
