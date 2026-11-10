#!/usr/bin/env python3
"""Compare baseline --glv-orbit vs --glv-orbit --phase-signature on a 135-bit window."""
import re
import subprocess
import sys

BIN = "/home/z/my-project/ECDSA/keyhunt/eng_bsgs"
START = 1 << 134


def run(extra):
    cmd = [BIN, "--target", TARGET, "--start", hex(START), "--size", "100000",
           "--segment", "100000", "--threads", "2", *extra]
    r = subprocess.run(cmd, text=True, capture_output=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    out = r.stdout
    def g(name):
        m = re.search(rf"{re.escape(name)}: (.+)", out)
        return m.group(1).strip() if m else "?"
    return {
        "raw": int(g("raw baby rows built")),
        "unique": int(g("unique baby rows built")),
        "dup": int(g("duplicate baby keys")),
        "hits": int(g("invariant hits")),
        "factors": int(g("orbit factors scheduled")),
        "exact": int(g("exact candidate checks")),
        "routes": int(g("signature routes")) if "signature routes" in out else None,
        "anomalies": int(g("signature anomalies (6-branch fallback)")) if "signature anomalies" in out else None,
        "build_s": float(g("table build seconds")),
        "search_s": float(g("giant search seconds")),
        "recovered": g("recovered scalar"),
        "correct": "recovery correct: true" in out,
    }


if __name__ == "__main__":
    TARGET = sys.argv[1]
    base = run(["--glv-orbit"])
    sig = run(["--glv-orbit", "--phase-signature"])
    print("baseline --glv-orbit:")
    for k, v in base.items():
        print(f"  {k}: {v}")
    print("--glv-orbit --phase-signature:")
    for k, v in sig.items():
        print(f"  {k}: {v}")
    print("\ncomparison:")
    print(f"  same recovered scalar: {base['recovered'] == sig['recovered']}")
    print(f"  both correct: {base['correct'] and sig['correct']}")
    print(f"  factors scheduled: baseline {base['factors']} vs signature {sig['factors']}")
    print(f"  exact candidate checks: baseline {base['exact']} vs signature {sig['exact']}")
    print(f"  unique K3 keys: {sig['unique']} of {sig['raw']} rows (duplicates {sig['dup']})")
    print(f"  signature routes: {sig['routes']}, anomalies: {sig['anomalies']}")
    print(f"  table build seconds: baseline {base['build_s']:.4f} vs signature {sig['build_s']:.4f}")
