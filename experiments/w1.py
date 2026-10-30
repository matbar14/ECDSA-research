#!/usr/bin/env python3
"""Exploratory analysis of endomorphism-orbit invariants on secp256k1.

For every sampled scalar k the script computes P = kG, its endomorphism
orbit, torus coordinates (x^3, y), and Lambda(P). It also reports two
separate assertions:

* the algebraic identity P = kG;
* the user-supplied endpoint hypothesis: d turns are completed at h = k,
  and the torus is declared absent at that endpoint.

The endpoint hypothesis is reported as a model assumption. It is not derived
from the secp256k1 group law or from Lambda(P).
"""

from fractions import Fraction
import secrets
import time


# secp256k1 parameters
FIELD_PRIME = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
GROUP_ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
ENDOMORPHISM_BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE

# Analysis parameters. Change these values for a different experiment.
SAMPLE_COUNT = 8
TORUS_TURNS = 3


class EllipticCurvePoint:
    """Affine secp256k1 point, including the point at infinity."""

    def __init__(self, x=None, y=None, is_infinity=False):
        if is_infinity:
            self.x = None
            self.y = None
            self.is_infinity = True
        else:
            self.x = x % FIELD_PRIME if x is not None else None
            self.y = y % FIELD_PRIME if y is not None else None
            self.is_infinity = False

    def __add__(self, other):
        if not isinstance(other, EllipticCurvePoint):
            return NotImplemented
        if self.is_infinity:
            return other
        if other.is_infinity:
            return self

        if self.x == other.x:
            if self.y != other.y or self.y == 0:
                return EllipticCurvePoint(is_infinity=True)
            slope = (3 * self.x * self.x * pow(2 * self.y, -1, FIELD_PRIME)) % FIELD_PRIME
        else:
            slope = (
                (other.y - self.y) * pow(other.x - self.x, -1, FIELD_PRIME)
            ) % FIELD_PRIME

        x3 = (slope * slope - self.x - other.x) % FIELD_PRIME
        y3 = (slope * (self.x - x3) - self.y) % FIELD_PRIME
        return EllipticCurvePoint(x3, y3)

    def __mul__(self, scalar):
        if not isinstance(scalar, int):
            return NotImplemented
        if scalar < 0:
            raise ValueError("Scalar must be non-negative")

        result = EllipticCurvePoint(is_infinity=True)
        addend = self
        while scalar:
            if scalar & 1:
                result = result + addend
            addend = addend + addend
            scalar >>= 1
        return result

    def __rmul__(self, scalar):
        return self * scalar


BASE_POINT = EllipticCurvePoint(
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)


class DomeInvariantAnalyzer:
    """Computes orbit-based quantities for sampled private scalars."""

    def __init__(self, base_point=BASE_POINT):
        self.base_point = base_point

    @staticmethod
    def points_equal(left, right):
        return (
            isinstance(left, EllipticCurvePoint)
            and isinstance(right, EllipticCurvePoint)
            and left.is_infinity == right.is_infinity
            and left.x == right.x
            and left.y == right.y
        )

    @staticmethod
    def is_valid_curve_point(point):
        return (
            isinstance(point, EllipticCurvePoint)
            and not point.is_infinity
            and point.x is not None
            and point.y is not None
            and (point.y * point.y - point.x * point.x * point.x - 7) % FIELD_PRIME == 0
        )

    def generate_large_random_key(self):
        """Returns a scalar with a bit length in the interval [200, 256]."""
        bit_size = secrets.randbelow(57) + 200
        scalar = secrets.randbits(bit_size) % GROUP_ORDER
        return scalar if scalar else 1

    def verify_height_equals_k(self, height, projection):
        """Checks the algebraic statement P = kG for a canonical scalar k."""
        canonical_height = isinstance(height, int) and 1 <= height < GROUP_ORDER
        valid_projection = self.is_valid_curve_point(projection)
        projection_matches_height = (
            self.points_equal(height * self.base_point, projection)
            if canonical_height and valid_projection
            else False
        )
        return {
            "canonical_height": canonical_height,
            "valid_projection": valid_projection,
            "projection_matches_height": projection_matches_height,
            "valid": canonical_height and valid_projection and projection_matches_height,
        }

    def analyze_dome_structure(self, height):
        projection = height * self.base_point
        height_proof = self.verify_height_equals_k(height, projection)

        print("Dome structure")
        print(f"  Height k: {height}")
        print(f"  Height size: {height.bit_length()} bits")
        print(f"  Projection P.x: 0x{projection.x:064x}")
        print(f"  Projection P.y: 0x{projection.y:064x}")
        print(f"  Identity P = kG: {'valid' if height_proof['valid'] else 'invalid'}")

        return {
            "height": height,
            "height_bits": height.bit_length(),
            "projection": projection,
            "height_proof": height_proof,
        }

    @staticmethod
    def torus_coordinates(point):
        return pow(point.x, 3, FIELD_PRIME), point.y

    def analyze_dome_base_torus(self, projection):
        torus_x, torus_y = self.torus_coordinates(projection)
        radius = (torus_x + torus_y) % FIELD_PRIME
        area = (torus_x * torus_y) % FIELD_PRIME
        curvature = (torus_x - torus_y) % FIELD_PRIME

        print("Torus coordinates")
        print(f"  x^3: 0x{torus_x:064x}")
        print(f"  y:   0x{torus_y:064x}")
        print(f"  Radius:    0x{radius:064x}")
        print(f"  Area:      0x{area:064x}")
        print(f"  Curvature: 0x{curvature:064x}")

        return {
            "coordinates": (torus_x, torus_y),
            "radius": radius,
            "area": area,
            "curvature": curvature,
        }

    def analyze_dome_geometry(self, height, projection, torus_data):
        base_radius = torus_data["radius"]
        base_area = torus_data["area"]
        height_to_radius = (
            height * pow(base_radius, -1, FIELD_PRIME) % FIELD_PRIME
            if base_radius
            else None
        )
        volume = height * base_area % FIELD_PRIME
        slope_x = height * pow(projection.x, -1, FIELD_PRIME) % FIELD_PRIME
        slope_y = height * pow(projection.y, -1, FIELD_PRIME) % FIELD_PRIME

        print("Derived quantities")
        print(f"  Height/radius: {format_hex(height_to_radius)}")
        print(f"  Height*area:   0x{volume:064x}")
        print(f"  Height/P.x:    0x{slope_x:064x}")
        print(f"  Height/P.y:    0x{slope_y:064x}")

        return {
            "height_to_radius": height_to_radius,
            "volume": volume,
            "slope_x": slope_x,
            "slope_y": slope_y,
        }

    def analyze_endomorphism_orbit(self, projection):
        """Computes {P, phi(P), phi^2(P)} for phi(x, y) = (beta*x, y)."""
        phi_projection = EllipticCurvePoint(
            ENDOMORPHISM_BETA * projection.x % FIELD_PRIME, projection.y
        )
        phi2_projection = EllipticCurvePoint(
            ENDOMORPHISM_BETA * ENDOMORPHISM_BETA * projection.x % FIELD_PRIME,
            projection.y,
        )
        orbit = [projection, phi_projection, phi2_projection]
        orbit_torus = [self.torus_coordinates(point) for point in orbit]
        center_x = sum(point[0] for point in orbit_torus) % FIELD_PRIME
        center_y = sum(point[1] for point in orbit_torus) % FIELD_PRIME

        print("Endomorphism orbit")
        for index, point in enumerate(orbit):
            print(f"  phi^{index}(P): (0x{point.x:064x}, 0x{point.y:064x})")
        print(f"  Orbit torus center: (0x{center_x:064x}, 0x{center_y:064x})")

        return {
            "points": orbit,
            "torus_coordinates": orbit_torus,
            "center": (center_x, center_y),
        }

    def analyze_endomorphism_orbit_silent(self, projection):
        """Computes orbit data without producing sample-level output."""
        phi_projection = EllipticCurvePoint(
            ENDOMORPHISM_BETA * projection.x % FIELD_PRIME, projection.y
        )
        phi2_projection = EllipticCurvePoint(
            ENDOMORPHISM_BETA * ENDOMORPHISM_BETA * projection.x % FIELD_PRIME,
            projection.y,
        )
        orbit = [projection, phi_projection, phi2_projection]
        orbit_torus = [self.torus_coordinates(point) for point in orbit]
        return {
            "points": orbit,
            "torus_coordinates": orbit_torus,
            "center": (
                sum(point[0] for point in orbit_torus) % FIELD_PRIME,
                sum(point[1] for point in orbit_torus) % FIELD_PRIME,
            ),
        }

    def compute_torus_invariant_lambda(self, projection):
        """Compute Lambda(P), with zero assigned at a boundary denominator."""
        x_base = self.base_point.x
        x_projection = projection.x
        phi_base_x = ENDOMORPHISM_BETA * x_base % FIELD_PRIME
        phi_projection_x = ENDOMORPHISM_BETA * x_projection % FIELD_PRIME

        numerator = (
            (x_base - phi_base_x) * (x_projection - phi_projection_x)
        ) % FIELD_PRIME
        denominator = (
            (x_base - phi_projection_x) * (x_projection - phi_base_x)
        ) % FIELD_PRIME
        if denominator == 0:
            # Boundary convention: assign the first-radius angle theta = 0.
            return 0
        return numerator * pow(denominator, -1, FIELD_PRIME) % FIELD_PRIME

    def analyze_torus_invariants(self, orbit_data):
        lambda_values = [
            self.compute_torus_invariant_lambda(point) for point in orbit_data["points"]
        ]
        invariant_defined = all(value is not None for value in lambda_values)
        orbit_sum = sum(lambda_values) % FIELD_PRIME if invariant_defined else None

        print("Torus invariants")
        for index, value in enumerate(lambda_values):
            print(f"  Lambda_{index}: {format_hex(value)}")
        print(f"  Sum Lambda_i: {format_hex(orbit_sum)}")
        print(f"  Sum Lambda_i = 0: {'yes' if orbit_sum == 0 else 'no'}")

        return {
            "lambda_values": lambda_values,
            "orbit_sum": orbit_sum,
            "sum_is_zero": orbit_sum == 0,
        }

    @staticmethod
    def analyze_torus_turn_hypothesis(height, turns):
        """Reports the endpoint condition stipulated by the turn hypothesis."""
        if not isinstance(turns, int) or turns <= 0:
            raise ValueError("TORUS_TURNS must be a positive integer")

        final_phase = Fraction(turns * height, height)
        prior_phase = Fraction(turns * (height - 1), height)
        result = {
            "turns": turns,
            "disappearance_height": height,
            "phase_before_endpoint": prior_phase,
            "phase_at_endpoint": final_phase,
            "completed_turns_at_endpoint": final_phase == turns,
            "torus_present_before_endpoint": height > 0,
            "torus_absent_at_endpoint": True,
        }

        print("Torus turn hypothesis")
        print(f"  Claimed turns d: {turns}")
        print(f"  Phase at h = k: {final_phase}")
        print(f"  Completed turns at h = k: {'yes' if result['completed_turns_at_endpoint'] else 'no'}")
        print("  Torus state before h = k: present")
        print("  Torus state at h = k: absent (model assumption)")
        return result

    def analyze_sample(self, sample_id, turns):
        print(f"Sample {sample_id}")
        height = self.generate_large_random_key()
        dome_data = self.analyze_dome_structure(height)
        torus_data = self.analyze_dome_base_torus(dome_data["projection"])
        geometry_data = self.analyze_dome_geometry(height, dome_data["projection"], torus_data)
        orbit_data = self.analyze_endomorphism_orbit(dome_data["projection"])
        invariant_data = self.analyze_torus_invariants(orbit_data)
        turn_data = self.analyze_torus_turn_hypothesis(height, turns)

        return {
            "sample_id": sample_id,
            "dome": dome_data,
            "torus": torus_data,
            "geometry": geometry_data,
            "orbit": orbit_data,
            "invariants": invariant_data,
            "turn_hypothesis": turn_data,
        }


def format_hex(value):
    return "undefined" if value is None else f"0x{value:064x}"


def print_lambda_orbit_identity():
    """Print the symbolic identity used by the orbit-sum experiment."""
    beta_cubed = pow(ENDOMORPHISM_BETA, 3, FIELD_PRIME)
    cyclotomic_factor = (
        ENDOMORPHISM_BETA * ENDOMORPHISM_BETA + ENDOMORPHISM_BETA + 1
    ) % FIELD_PRIME

    print("Lambda orbit identity modulo p")
    print("  z = x(P) / x(G)")
    print("  L(z) = (1-beta)^2*z / ((1-beta*z)*(z-beta))")
    print("  S(z) = L(z) + L(beta*z) + L(beta^2*z)")
    print("  numerator(S) =")
    print("    z*(beta-1)^2*(beta*z-1)^2*(beta^2+beta+1)*Q(beta,z)")
    print(f"  beta^3 mod p: {beta_cubed}")
    print(f"  beta^2 + beta + 1 mod p: {cyclotomic_factor}")
    print("  Therefore S(z) = 0 mod p when the Lambda denominators are nonzero.")


def analyze_patterns(results):
    key_sizes = [result["dome"]["height_bits"] for result in results]
    invariant_sums = [result["invariants"]["sum_is_zero"] for result in results]
    height_claims = [result["dome"]["height_proof"]["valid"] for result in results]
    centers = {result["orbit"]["center"] for result in results}

    print("Analysis summary")
    print(f"  Samples: {len(results)}")
    print(f"  Scalar size range: {min(key_sizes)}-{max(key_sizes)} bits")
    print(f"  Valid P = kG identities: {sum(height_claims)}/{len(results)}")
    print(f"  Lambda orbit sums equal to zero: {sum(invariant_sums)}/{len(results)}")
    print(f"  Distinct orbit centers: {len(centers)}")
    print(f"  Turn hypothesis d: {TORUS_TURNS}")


def run_256_bit_invariant_experiment(sample_count, turns):
    """Run a non-verbose invariant experiment on exactly 256-bit scalars.

    Public points are computed by the native OpenSSL-backed secp256k1
    implementation provided by ``cryptography``. The orbit and Lambda
    calculations remain those implemented in this file.
    """
    from cryptography.hazmat.primitives.asymmetric import ec

    if not isinstance(sample_count, int) or sample_count <= 0:
        raise ValueError("sample_count must be a positive integer")
    if not isinstance(turns, int) or turns <= 0:
        raise ValueError("turns must be a positive integer")

    analyzer = DomeInvariantAnalyzer()
    minimum_256_bit_scalar = 1 << 255
    failures = {
        "invalid_curve_point": 0,
        "undefined_lambda": 0,
        "nonzero_lambda_sum": 0,
        "endpoint_phase_not_d": 0,
    }
    start_time = time.perf_counter()

    for sample_index in range(1, sample_count + 1):
        scalar = minimum_256_bit_scalar + secrets.randbelow(
            GROUP_ORDER - minimum_256_bit_scalar
        )
        public_numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
        projection = EllipticCurvePoint(public_numbers.x, public_numbers.y)

        if not analyzer.is_valid_curve_point(projection):
            failures["invalid_curve_point"] += 1

        orbit_data = analyzer.analyze_endomorphism_orbit_silent(projection)
        lambda_values = [
            analyzer.compute_torus_invariant_lambda(point) for point in orbit_data["points"]
        ]
        if any(value is None for value in lambda_values):
            failures["undefined_lambda"] += 1
        elif sum(lambda_values) % FIELD_PRIME != 0:
            failures["nonzero_lambda_sum"] += 1

        if Fraction(turns * scalar, scalar) != turns:
            failures["endpoint_phase_not_d"] += 1

        if sample_index % 10000 == 0:
            print(f"  Completed: {sample_index}/{sample_count}")

    elapsed_seconds = time.perf_counter() - start_time
    return {
        "sample_count": sample_count,
        "turns": turns,
        "scalar_bit_length": 256,
        "failures": failures,
        "elapsed_seconds": elapsed_seconds,
    }


def print_256_bit_experiment_report(report):
    print("256-bit invariant experiment")
    print(f"  Samples: {report['sample_count']}")
    print(f"  Turn count d: {report['turns']}")
    print(f"  Invalid curve points: {report['failures']['invalid_curve_point']}")
    print(f"  Undefined Lambda values: {report['failures']['undefined_lambda']}")
    print(f"  Nonzero Lambda sums: {report['failures']['nonzero_lambda_sum']}")
    print(f"  Endpoint phases not equal to d: {report['failures']['endpoint_phase_not_d']}")
    print(f"  Elapsed time: {report['elapsed_seconds']:.2f} s")


def main():
    analyzer = DomeInvariantAnalyzer()
    print_lambda_orbit_identity()
    results = [
        analyzer.analyze_sample(sample_id, TORUS_TURNS)
        for sample_id in range(1, SAMPLE_COUNT + 1)
    ]
    analyze_patterns(results)


if __name__ == "__main__":
    main()
