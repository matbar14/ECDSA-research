#!/usr/bin/env python3
"""Phase-torus experiment for the Lambda orbit invariant on secp256k1.

The Lambda values are computed with the exact finite-field formula from w1.py.
To draw them in ordinary Euclidean geometry, this experiment explicitly chooses
the additive-character map

    chi(a) = exp(2*pi*i*lift(a)/p),

where lift(a) is the representative in [0, p).  Thus each Lambda value is a
point on the unit circle.  Since the Lambda orbit sum is zero modulo p,

    chi(Lambda(P)) chi(Lambda(phi(P))) chi(Lambda(phi^2(P))) = 1.

The ordered triple therefore lies on the two-dimensional phase torus
{(u0,u1,u2) in (S1)^3 : u0*u1*u2 = 1}.  We use the first two phases as torus
coordinates and draw the three unit-circle points as an inscribed triangle.

The map chi is a declared visualization convention, not a canonical Euclidean
embedding of F_p.  Areas are numerical Euclidean measurements of that model.
Exact invariant-triple collisions are checked with field residues.  Numerical
area coincidences are reported separately at a stated decimal precision.
"""

from collections import defaultdict
import html
import math
import random
import time

from cryptography.hazmat.primitives.asymmetric import ec

from w1 import (
    DomeInvariantAnalyzer,
    ENDOMORPHISM_BETA,
    EllipticCurvePoint,
    FIELD_PRIME,
    GROUP_ORDER,
)


SAMPLE_COUNT = 100000
SEED = 20260917
PROGRESS_INTERVAL = 10000
SVG_OUTPUT = "figures/lambda_phase_torus.svg"
AREA_DECIMAL_PLACES = 12

# A conventional torus embedding used only for the report figure.
TORUS_MAJOR_RADIUS = 2.4
TORUS_MINOR_RADIUS = 0.8


def generate_256_bit_scalar(generator):
    """Return a reproducible scalar whose bit length is exactly 256."""
    return generator.randrange(1 << 255, GROUP_ORDER)


def point_from_scalar(scalar):
    """Compute scalar*G through the native OpenSSL-backed implementation."""
    public_numbers = ec.derive_private_key(
        scalar, ec.SECP256K1()
    ).public_key().public_numbers()
    return EllipticCurvePoint(public_numbers.x, public_numbers.y)


def lambda_orbit(analyzer, point):
    """Return the exact ordered Lambda triple for P, phi(P), phi^2(P)."""
    phi_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * point.x % FIELD_PRIME, point.y
    )
    phi2_point = EllipticCurvePoint(
        ENDOMORPHISM_BETA * phi_point.x % FIELD_PRIME, point.y
    )
    values = tuple(
        analyzer.compute_torus_invariant_lambda(orbit_point)
        for orbit_point in (point, phi_point, phi2_point)
    )
    if any(value is None for value in values):
        return None
    return values


def phase(value):
    """Map a residue to its unit-circle angle in radians.

    ``value / FIELD_PRIME`` is intentionally converted only at this final
    drawing step.  All cryptographic computations remain exact integers.
    """
    return math.tau * (value / FIELD_PRIME)


def unit_circle_vertices(values):
    """Return the three complex unit-circle vertices for a Lambda triple."""
    return tuple(complex(math.cos(phase(value)), math.sin(phase(value))) for value in values)


def triangle_area(vertices):
    """Return the unsigned Euclidean area by the shoelace determinant."""
    first, second, third = vertices
    double_area = (
        first.real * (second.imag - third.imag)
        + second.real * (third.imag - first.imag)
        + third.real * (first.imag - second.imag)
    )
    return abs(double_area) / 2.0


def orbit_fingerprint(values):
    """Exact key of a triangle, invariant under cyclic orbit reordering."""
    return tuple(sorted(values))


def torus_point(values):
    """Embed the first two phases in a conventional Euclidean torus in R^3."""
    major_angle = phase(values[0])
    minor_angle = phase(values[1])
    radius = TORUS_MAJOR_RADIUS + TORUS_MINOR_RADIUS * math.cos(minor_angle)
    return (
        radius * math.cos(major_angle),
        radius * math.sin(major_angle),
        TORUS_MINOR_RADIUS * math.sin(minor_angle),
    )


def format_scalar(scalar):
    return f"0x{scalar:064x}"


def run_experiment(sample_count=SAMPLE_COUNT, seed=SEED, fingerprint_path=None):
    """Evaluate areas and exact collision classes for 256-bit source scalars."""
    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    generator = random.Random(seed)
    fingerprints = defaultdict(list)
    areas = []
    undefined_lambda_count = 0
    nonzero_sum_count = 0
    area_buckets = defaultdict(list)
    representative = None
    start_time = time.perf_counter()

    fingerprint_file = (
        open(fingerprint_path, "w", encoding="ascii")
        if fingerprint_path is not None
        else None
    )
    try:
        for sample_index in range(1, sample_count + 1):
            scalar = generate_256_bit_scalar(generator)
            point = point_from_scalar(scalar)
            values = lambda_orbit(analyzer, point)
            if values is None:
                undefined_lambda_count += 1
                continue

            if sum(values) % FIELD_PRIME != 0:
                nonzero_sum_count += 1

            vertices = unit_circle_vertices(values)
            area = triangle_area(vertices)
            fingerprint = orbit_fingerprint(values)
            fingerprints[fingerprint].append(scalar)
            if fingerprint_file is not None:
                fingerprint_file.write(
                    "".join(f"{value:064x}" for value in fingerprint)
                    + f" {area:.17g} {scalar:064x}\n"
                )
            areas.append(area)
            area_buckets[round(area, AREA_DECIMAL_PLACES)].append(
                (scalar, values)
            )

            if representative is None:
                representative = {
                    "scalar": scalar,
                    "values": values,
                    "vertices": vertices,
                    "area": area,
                    "torus_point": torus_point(values),
                }

            if sample_index % PROGRESS_INTERVAL == 0:
                print(f"  Completed: {sample_index}/{sample_count}")
    finally:
        if fingerprint_file is not None:
            fingerprint_file.close()

    repeated_source_classes = {
        fingerprint: scalars
        for fingerprint, scalars in fingerprints.items()
        if len(scalars) > 1
    }
    source_collision_pairs = sum(
        len(scalars) * (len(scalars) - 1) // 2
        for scalars in repeated_source_classes.values()
    )
    repeated_area_buckets = {
        area: entries for area, entries in area_buckets.items() if len(entries) > 1
    }
    sorted_areas = sorted(areas)
    minimum_area_gap = (
        min(right - left for left, right in zip(sorted_areas, sorted_areas[1:]))
        if len(sorted_areas) > 1
        else None
    )

    return {
        "sample_count": sample_count,
        "seed": seed,
        "processed": len(areas),
        "undefined_lambda_count": undefined_lambda_count,
        "nonzero_sum_count": nonzero_sum_count,
        "source_distinct_triangle_classes": len(fingerprints),
        "source_repeated_triangle_classes": len(repeated_source_classes),
        "source_collision_pairs": source_collision_pairs,
        "repeated_area_buckets": repeated_area_buckets,
        "minimum_area_gap": minimum_area_gap,
        "area_min": min(areas),
        "area_max": max(areas),
        "area_mean": math.fsum(areas) / len(areas),
        "representative": representative,
        "elapsed_seconds": time.perf_counter() - start_time,
    }


def svg_point(point, center_x, center_y, radius):
    return center_x + radius * point.real, center_y - radius * point.imag


def project_torus(point):
    """Fixed oblique projection from R^3 to the right-hand SVG panel."""
    x, y, z = point
    return 635 + 58 * (x - 0.62 * y), 185 - 58 * (z + 0.24 * (x + y))


def torus_mesh_paths():
    """Return projected meridian and longitude paths for a compact wireframe."""
    paths = []
    steps = 72
    for fixed_major in range(0, 12):
        major = math.tau * fixed_major / 12
        points = []
        for index in range(steps + 1):
            minor = math.tau * index / steps
            radius = TORUS_MAJOR_RADIUS + TORUS_MINOR_RADIUS * math.cos(minor)
            points.append(project_torus((radius * math.cos(major), radius * math.sin(major), TORUS_MINOR_RADIUS * math.sin(minor))))
        paths.append(points)
    for fixed_minor in range(0, 12):
        minor = math.tau * fixed_minor / 12
        points = []
        for index in range(steps + 1):
            major = math.tau * index / steps
            radius = TORUS_MAJOR_RADIUS + TORUS_MINOR_RADIUS * math.cos(minor)
            points.append(project_torus((radius * math.cos(major), radius * math.sin(major), TORUS_MINOR_RADIUS * math.sin(minor))))
        paths.append(points)
    return paths


def svg_path(points):
    return "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def write_svg(report, output_path=SVG_OUTPUT):
    """Write one representative phase torus and its unit-circle triangle."""
    sample = report["representative"]
    vertices = sample["vertices"]
    circle_points = [svg_point(vertex, 185, 185, 112) for vertex in vertices]
    polygon = " ".join(f"{x:.2f},{y:.2f}" for x, y in circle_points)
    torus_x, torus_y = project_torus(sample["torus_point"])
    mesh = "\n".join(
        f'<path d="{svg_path(path)}" class="mesh" />' for path in torus_mesh_paths()
    )
    labels = "\n".join(
        f'<text x="{x + 8:.2f}" y="{y - 8:.2f}" class="vertex">L{index}</text>'
        for index, (x, y) in enumerate(circle_points)
    )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="840" height="365" viewBox="0 0 840 365" role="img" aria-labelledby="title description">
  <title id="title">Фазовый тор Lambda и вписанный треугольник</title>
  <desc id="description">Треугольник на единичной окружности построен из трёх значений Lambda орбиты одной точки secp256k1. Справа показана точка этой тройки на параметрическом двухмерном фазовом торе.</desc>
  <style>
    .title {{ font: 600 16px sans-serif; fill: #18212b; }}
    .label {{ font: 13px sans-serif; fill: #425466; }}
    .vertex {{ font: 600 13px sans-serif; fill: #9d174d; }}
    .mesh {{ fill: none; stroke: #9fb8c5; stroke-width: .8; }}
    .frame {{ fill: none; stroke: #607586; stroke-width: 1.2; }}
  </style>
  <rect width="840" height="365" fill="#ffffff" />
  <text x="36" y="32" class="title">Единичная окружность: chi(Lambda_i)</text>
  <circle cx="185" cy="185" r="112" class="frame" />
  <line x1="62" y1="185" x2="308" y2="185" stroke="#d5dde3" />
  <line x1="185" y1="62" x2="185" y2="308" stroke="#d5dde3" />
  <polygon points="{polygon}" fill="#f9a8d4" fill-opacity=".28" stroke="#be185d" stroke-width="2" />
  {''.join(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4" fill="#9d174d" />' for x, y in circle_points)}
  {labels}
  <text x="68" y="340" class="label">A = {sample['area']:.12f}</text>
  <text x="492" y="32" class="title">Фазовый тор T²: chi(L0) chi(L1) chi(L2) = 1</text>
  {mesh}
  <circle cx="{torus_x:.2f}" cy="{torus_y:.2f}" r="5" fill="#0f766e" />
  <text x="{torus_x + 10:.2f}" y="{torus_y - 9:.2f}" class="vertex">(L0, L1)</text>
  <text x="492" y="340" class="label">k = {html.escape(format_scalar(sample['scalar']))}</text>
</svg>
'''
    with open(output_path, "w", encoding="utf-8") as output_file:
        output_file.write(svg)


def print_report(report):
    sample = report["representative"]
    print("Lambda phase-torus triangle experiment")
    print(f"  256-bit source scalars: {report['sample_count']}")
    print(f"  Reproducibility seed: {report['seed']}")
    print(f"  Processed: {report['processed']}")
    print(f"  Undefined Lambda triples: {report['undefined_lambda_count']}")
    print(f"  Nonzero Lambda orbit sums: {report['nonzero_sum_count']}")
    print(f"  Distinct exact triangle classes among sources: {report['source_distinct_triangle_classes']}")
    print(f"  Repeated exact triangle classes among sources: {report['source_repeated_triangle_classes']}")
    print(f"  Exact source collision pairs: {report['source_collision_pairs']}")
    print("  Euclidean areas under chi(a) = exp(2*pi*i*a/p):")
    print(f"    minimum: {report['area_min']:.12f}")
    print(f"    maximum: {report['area_max']:.12f}")
    print(f"    mean:    {report['area_mean']:.12f}")
    print(f"  Repeated areas after rounding to {AREA_DECIMAL_PLACES} decimals: {len(report['repeated_area_buckets'])}")
    if report["minimum_area_gap"] is not None:
        print(f"  Minimum numerical area gap: {report['minimum_area_gap']:.16g}")
    print("  Representative sample")
    print(f"    k: {format_scalar(sample['scalar'])}")
    for index, value in enumerate(sample["values"]):
        print(f"    Lambda_{index}: 0x{value:064x}")
    print(f"    triangle area: {sample['area']:.12f}")
    print("    torus point: (" + ", ".join(f"{value:.12f}" for value in sample["torus_point"]) + ")")
    print(f"  SVG figure: {SVG_OUTPUT}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    report = run_experiment()
    write_svg(report)
    print_report(report)


if __name__ == "__main__":
    main()
