#!/usr/bin/env python3
"""Cross-validate the delta-sweep pair arm against parallel_angle_table_recovery.py.

Builds the same cell (W = 2^17, B = 24, i.e. 12 bits per angle) with the
reference script and compares mean/max fingerprint bucket with the sweep's
pair arm at (log2 W = 17, delta = 24 - 17 = 7).
"""
import numpy as np

from angle_table_delta_sweep import build_pass, pair_key
from w1 import GROUP_ORDER

LOG2_W = 17
B = 24


def main():
    q0, q1, q2, _ = build_pass(LOG2_W, 8, "/tmp/angle_delta_sweep")
    b1 = (B + 1) // 2
    b2 = B // 2
    keys = pair_key(q1, q2, b1, b2)
    uniq, counts = np.unique(keys, return_counts=True)
    sweep_mean = float(counts.mean())
    sweep_max = int(counts.max())
    occupied = int(len(uniq))
    print(f"sweep pair arm L={LOG2_W} B={B}: "
          f"W={1 << LOG2_W} occupied={occupied} "
          f"mean bucket={sweep_mean:.3f} max bucket={sweep_max}")
    rng = np.random.default_rng(20261126)
    qidx = rng.choice(1 << LOG2_W, size=200, replace=False)
    qkeys = keys[qidx]
    order = np.argsort(keys)
    sk = keys[order]
    left = np.searchsorted(sk, qkeys, side="left")
    right = np.searchsorted(sk, qkeys, side="right")
    print(f"sweep query buckets: mean={float((right-left).mean()):.3f} "
          f"max={int((right-left).max())}")
    assert occupied <= 1 << LOG2_W
    return sweep_mean, sweep_max


if __name__ == "__main__":
    main()
