#!/usr/bin/env python3
"""Independent follow-up for one scalar-height band from area_shift_experiment.

The first stratified screen found the largest KS distance in band 2 of 8.
This script uses fresh samples and compares only that preselected band against
the equal-weight mixture of the other seven bands.  It reports the mean shift
and two-sample Kolmogorov-Smirnov statistic without searching new bands.
"""

import math
import random
import time

from area_shift_experiment import (
    BAND_COUNT,
    GROUP_ORDER,
    area_from_scalar,
    ks_distance,
    population_standard_deviation,
)
from w1 import DomeInvariantAnalyzer


TARGET_BAND = 2
SAMPLES_PER_GROUP = 10000
SEED = 20260928


def band_bounds(index):
    lower = 1 << 255
    width = (GROUP_ORDER - lower) // BAND_COUNT
    band_lower = lower + index * width
    band_upper = GROUP_ORDER if index == BAND_COUNT - 1 else band_lower + width
    return band_lower, band_upper


def ks_asymptotic_pvalue(distance, first_count, second_count):
    """Return the standard asymptotic two-sample KS p-value estimate."""
    effective_count = first_count * second_count / (first_count + second_count)
    root = math.sqrt(effective_count)
    argument = (root + 0.12 + 0.11 / root) * distance
    series = 0.0
    for term_index in range(1, 100):
        term = math.exp(-2.0 * term_index * term_index * argument * argument)
        series += term if term_index % 2 else -term
        if term < 1e-15:
            break
    return max(0.0, min(1.0, 2.0 * series))


def run_experiment(samples_per_group=SAMPLES_PER_GROUP, seed=SEED):
    if not isinstance(samples_per_group, int) or samples_per_group <= 0:
        raise ValueError("samples_per_group must be positive")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    target_lower, target_upper = band_bounds(TARGET_BAND)
    target_areas = []
    other_areas = []
    start_time = time.perf_counter()

    for _ in range(samples_per_group):
        target_scalar = generator.randrange(target_lower, target_upper)
        target_areas.append(area_from_scalar(analyzer, target_scalar))

        other_band = generator.randrange(BAND_COUNT - 1)
        if other_band >= TARGET_BAND:
            other_band += 1
        other_lower, other_upper = band_bounds(other_band)
        other_scalar = generator.randrange(other_lower, other_upper)
        other_areas.append(area_from_scalar(analyzer, other_scalar))

    if any(area is None for area in target_areas + other_areas):
        raise RuntimeError("Lambda area was unexpectedly undefined")

    target_mean = math.fsum(target_areas) / len(target_areas)
    other_mean = math.fsum(other_areas) / len(other_areas)
    target_stddev = population_standard_deviation(target_areas, target_mean)
    other_stddev = population_standard_deviation(other_areas, other_mean)
    standard_error = math.sqrt(
        target_stddev * target_stddev / len(target_areas)
        + other_stddev * other_stddev / len(other_areas)
    )
    distance = ks_distance(target_areas, other_areas)
    return {
        "target_band": TARGET_BAND,
        "samples_per_group": samples_per_group,
        "seed": seed,
        "target_mean": target_mean,
        "other_mean": other_mean,
        "mean_shift": target_mean - other_mean,
        "standard_error": standard_error,
        "mean_z_score": (target_mean - other_mean) / standard_error,
        "ks_distance": distance,
        "ks_pvalue": ks_asymptotic_pvalue(
            distance, len(target_areas), len(other_areas)
        ),
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def print_report(report):
    print("Independent area-shift confirmation")
    print(f"  Target band: {report['target_band']} of {BAND_COUNT}")
    print(f"  Samples per group: {report['samples_per_group']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Target mean: {report['target_mean']:.12f}")
    print(f"  Other-band mean: {report['other_mean']:.12f}")
    print(f"  Mean shift: {report['mean_shift']:+.12f}")
    print(f"  Mean-shift z score: {report['mean_z_score']:.6f}")
    print(f"  KS distance: {report['ks_distance']:.8f}")
    print(f"  KS asymptotic p-value: {report['ks_pvalue']:.8f}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    print_report(run_experiment())


if __name__ == "__main__":
    main()
