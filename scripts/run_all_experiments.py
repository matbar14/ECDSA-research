#!/usr/bin/env python3
"""Batch smoke-run of ECDSA research experiments with small parameters.

Runs each experiment in a subprocess with reduced sample counts, records
exit codes, wall time and the last output lines, and prints a summary table.
"""
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/home/z/my-project/ECDSA")

# script -> (extra args, timeout seconds)
RUNS = {
    "angle_anchor_experiment.py": (["--range-size", "2000"], 240),
    "angle_global_locality_experiment.py": (["--samples", "2000"], 240),
    "angle_quantization_experiment.py": ([], 240),
    "angle_table_recovery_experiment.py": ([], 240),
    "area_lookup_table_experiment.py": ([], 240),
    "area_ratio_experiment.py": ([], 240),
    "area_multi_ratio_experiment.py": ([], 240),
    "area_angle_signature_experiment.py": ([], 240),
    "area_shift_confirmation.py": ([], 240),
    "orbit_period_experiment.py": ([], 240),
    "invariant_leakage_experiment.py": ([], 240),
    "invariant_candidate_search.py": (["--samples", "200"], 240),
    "cubic_normalized_invariant_experiment.py": (["--samples", "2000"], 240),
    "symmetric_orbit_invariant_experiment.py": ([], 240),
    "private_angle_experiment.py": ([], 240),
    "phase_difference_experiment.py": ([], 240),
    "phase_digit_bit_experiment.py": ([], 240),
    "lambda_quadratic_branch_experiment.py": ([], 240),
    "d_from_height_experiment.py": ([], 240),
    "lambda_circle_triangle_experiment.py": ([], 240),
    "torus_trajectory_experiment.py": (["--range-size", "1000"], 240),
    "marker_encoding_benchmark.py": (["--range-size", "20000"], 300),
    "marker_decoder_search.py": ([], 300),
    "storage_scaling_benchmark.py": (["--range-size", "20000"], 300),
    "segmented_bsgs_benchmark.py": ([], 300),
    "segmented_bsgs_recovery.py": ([], 300),
    "compact_angle_table_benchmark.py": (["--range-size", "20000"], 300),
    "glv_bsgs_benchmark.py": (["--sizes", "10000"], 300),
    "next_step_benchmark.py": (["--sizes", "10000", "--workers", "2"], 300),
    "parallel_segmented_bsgs.py": (["--window-size", "20000", "--workers", "2"], 300),
    "w1_circle_uniformity_experiment.py": (["--count", "2000"], 300),
    "w1_triangle_circle_experiment.py": ([], 300),
    "triangle_range_correlation_experiment.py": (["--samples", "5000"], 300),
    "triangle_rotation_collapse_experiment.py": (["--samples", "5000"], 300),
    "triangle_descending_g_experiment.py": (["--max-scalar", "2000"], 300),
    "random_135bit_recovery_experiment.py": ([], 300),
    "random_segmented_bsgs_test.py": ([], 300),
    "parallel_angle_table_recovery.py": ([], 300),
    "compact_parallel_angle_table.py": (["--range-size", "20000", "--workers", "2"], 300),
}


def main() -> int:
    only = sys.argv[1:] if len(sys.argv) > 1 else None
    results = []
    for script, (extra, timeout) in RUNS.items():
        if only and not any(o in script for o in only):
            continue
        t0 = time.time()
        try:
            proc = subprocess.run(
                [sys.executable, script, *extra],
                cwd=ROOT, text=True, capture_output=True, timeout=timeout,
            )
            dt = time.time() - t0
            status = "OK" if proc.returncode == 0 else f"FAIL({proc.returncode})"
            tail = (proc.stdout + proc.stderr).strip().splitlines()[-3:]
        except subprocess.TimeoutExpired:
            dt = time.time() - t0
            status = f"TIMEOUT(>{timeout}s)"
            tail = []
        results.append((script, status, dt, tail))
        print(f"[{status:>12}] {dt:7.1f}s  {script}", flush=True)
        for line in tail:
            print(f"               | {line}")

    print("\n==== SUMMARY ====")
    ok = sum(1 for _, s, _, _ in results if s == "OK")
    print(f"passed {ok}/{len(results)}")
    for script, status, dt, _ in results:
        if status != "OK":
            print(f"  {status}: {script}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
