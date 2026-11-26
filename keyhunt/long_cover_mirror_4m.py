#!/usr/bin/env python3
"""Paired H=4M CBT1 validation for default and mirror giant tilings."""

import argparse
import hashlib
import json
import os
import re
import statistics
import subprocess
import threading
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec


ROOT = Path(__file__).resolve().parents[1]
START = 1 << 134
ROWS = 4_000_000
DEFAULT_STRIDE = ROWS
MIRROR_STRIDE = 2 * ROWS - 1
GIANTS_PER_CHUNK = 1_000_000
CHUNK = GIANTS_PER_CHUNK * DEFAULT_STRIDE
COVER_CHUNKS = 12
COVER_SIZE = CHUNK * COVER_CHUNKS
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

INPUTS = {
    "canonical": ROOT / "tables/representatives_canonical_H4000000.txt",
    "mirror": ROOT / "tables/representatives_mirror_H4000000.txt",
}


def pubkey_hex(scalar):
    public = ec.derive_private_key(scalar % N, ec.SECP256K1()) \
        .public_key().public_numbers()
    return f"{2 | (public.y & 1):02x}{public.x:064x}"


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    os.replace(temp, path)


def parse_checkpoint(path):
    if not path.exists():
        return None
    fields = {}
    for line in path.read_text().splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            fields[key.strip()] = value.strip()
    return fields


def match_int(pattern, text):
    match = re.search(pattern, text)
    return int(match.group(1)) if match else None


def sample_peak_rss(pid, sink, stop):
    status = Path(f"/proc/{pid}/status")
    while not stop.is_set():
        try:
            for line in status.read_text().splitlines():
                if line.startswith("VmHWM:"):
                    sink[0] = max(sink[0], int(line.split()[1]))
                    break
        except (FileNotFoundError, ProcessLookupError):
            return
        stop.wait(0.01)


def run_process(command, cwd, timeout):
    proc = subprocess.Popen(command, cwd=cwd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    peak = [0]
    stop = threading.Event()
    sampler = threading.Thread(target=sample_peak_rss,
                               args=(proc.pid, peak, stop), daemon=True)
    sampler.start()
    started = time.perf_counter()
    killed = False
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        killed = True
        proc.kill()
        stdout, stderr = proc.communicate()
    stop.set()
    sampler.join(timeout=0.5)
    counters = {}
    for name in (
        "giant steps", "giant steps planned", "giant steps/second",
        "search hits", "invariant hits", "exact candidate checks",
        "partner recoveries", "automatic chunks attempted",
    ):
        value = match_int(rf"{re.escape(name)}: (\d+)", stdout)
        if value is None:
            number = re.search(
                rf"{re.escape(name)}: ([0-9]+(?:\.[0-9]+)?)", stdout)
            value = float(number.group(1)) if number else None
        if value is not None:
            counters[name] = value
    recovered = re.search(
        r"recovered private key: 0x([0-9a-fA-F]{64})", stdout)
    resumed = match_int(
        r"chunks completed before this session: (\d+)", stdout)
    return {
        "killed": killed,
        "returncode": None if killed else proc.returncode,
        "wall seconds": round(time.perf_counter() - started, 4),
        "peak rss kb": peak[0],
        "recovered": recovered.group(1).lower() if recovered else None,
        "resumed chunks before": resumed,
        "already recovered": "already recovered" in stdout,
        "counters": counters,
        "stdout tail": stdout[-1200:],
        "stderr tail": stderr[-1200:],
    }


def stride_for(mode):
    return DEFAULT_STRIDE if mode == "canonical" else MIRROR_STRIDE


def is_mirror(mode):
    return mode == "mirror"


def table_path(workdir, mode):
    return workdir / f"{mode}_H4000000.cbt1"


def target_flags(mode, stride):
    flags = ["--glv-orbit", "--partner-table"]
    if is_mirror(mode):
        flags += ["--giant-stride", str(stride)]
    return flags


def cover_command(binary, workdir, mode, stride, table, scalar, threads,
                  progress, checkpoint=None, resume=False):
    command = [
        str(binary), "--target", pubkey_hex(scalar),
        "--start", hex(START), "--size", str(CHUNK),
        "--segment", str(CHUNK), "--threads", str(threads),
        "--progress-interval", str(progress),
        "--cover-start", hex(START), "--cover-size", hex(COVER_SIZE),
        "--table-binary", str(table),
    ] + target_flags(mode, stride)
    if checkpoint is not None:
        command += ["--checkpoint-file", str(checkpoint),
                    "--checkpoint-steps", "1"]
    if resume:
        command.append("--resume")
    return command


def build_table(binary, workdir, mode, input_path, threads, target_scalar):
    table = table_path(workdir, mode)
    if table.exists():
        header = table.read_bytes()[:40]
        if len(header) < 40 or header[:4] != b"CBT1":
            raise RuntimeError(f"existing table is not CBT1: {table}")
        rows = int.from_bytes(header[32:40], "big")
        if rows != ROWS or table.stat().st_size != 256 + 65 * ROWS:
            raise RuntimeError(f"existing table has unexpected shape: {table}")
        return {
            "path": str(table), "bytes": table.stat().st_size,
            "sha256": sha256_file(table), "reused": True,
        }

    command = [
        str(binary), "--target", pubkey_hex(target_scalar),
        "--start", hex(START), "--size", str(CHUNK),
        "--segment", str(CHUNK), "--threads", str(threads),
        "--progress-interval", "0",
        "--table-scalars", str(input_path),
        "--glv-orbit", "--canonical-k3", "--partner-table",
        "--binary-format", "cbt1",
        "--dump-table-binary", str(table),
    ]
    result = run_process(command, binary.parent, timeout=1800)
    if result["killed"] or result["returncode"] != 0:
        raise RuntimeError(
            f"{mode} CBT1 build failed: rc={result['returncode']}\n"
            f"{result['stderr tail']}\n{result['stdout tail']}"
        )
    if result["recovered"] != f"{target_scalar:064x}":
        raise RuntimeError(f"{mode} table build control target mismatch")
    if table.stat().st_size != 256 + 65 * ROWS:
        raise RuntimeError(f"{mode} CBT1 has unexpected size")
    return {
        "path": str(table), "bytes": table.stat().st_size,
        "sha256": sha256_file(table), "reused": False,
        "build wall seconds": result["wall seconds"],
        "build peak rss kb": result["peak rss kb"],
        "build counters": result["counters"],
    }


def mirror_representation(offset):
    residue = offset % MIRROR_STRIDE
    if residue >= ROWS:
        return {"sign": "+", "center": offset // MIRROR_STRIDE,
                "row scalar": residue}
    return {"sign": "-", "center": offset // MIRROR_STRIDE + 1,
            "row scalar": MIRROR_STRIDE - residue}


def run_proof(binary, workdir, mode, stride, table, offset, threads,
              artifact, out_path, label):
    scalar = START + offset
    command = cover_command(binary, workdir, mode, stride, table, scalar,
                            threads, progress=0)
    result = run_process(command, binary.parent, timeout=300)
    expected = f"{scalar:064x}"
    ok = result["returncode"] == 0 and result["recovered"] == expected
    record = {
        "offset": offset, "expected scalar": expected,
        "pass": ok, "run": result,
    }
    if mode == "mirror":
        record["mirror representation"] = mirror_representation(offset)
    artifact["proofs"][f"{mode}:{label}"] = record
    if not ok:
        atomic_json(out_path, artifact)
        raise RuntimeError(
            f"{mode} {label} recovery failed: rc={result['returncode']} "
            f"recovered={result['recovered']}\n{result['stdout tail']}"
        )
    atomic_json(out_path, artifact)
    print(f"[{mode} proof {label}] PASS, "
          f"{result['counters'].get('giant steps')} actual giants")


def run_last_resume(binary, workdir, mode, stride, table, threads,
                    kill_after, artifact):
    state = artifact["resume"].setdefault(mode, {
        "offset": COVER_SIZE - 1, "sessions": [], "done": False,
    })
    if state.get("done"):
        print(f"[{mode} last resume] already PASS")
        return
    cp = workdir / f"cp_{mode}_last.txt"
    if cp.exists() and not state["sessions"]:
        raise RuntimeError(
            f"checkpoint exists without saved session metadata: {cp}; "
            "rename it before starting a fresh chain"
        )
    scalar = START + COVER_SIZE - 1
    while len(state["sessions"]) < 40:
        old_cp = parse_checkpoint(cp)
        use_resume = old_cp is not None
        command = cover_command(
            binary, workdir, mode, stride, table, scalar, threads,
            progress=5, checkpoint=cp, resume=use_resume,
        )
        session = run_process(command, binary.parent, timeout=kill_after)
        cp_after = parse_checkpoint(cp)
        session["checkpoint after"] = None if cp_after is None else {
            "chunks completed": int(cp_after.get("chunks completed", -1)),
            "found": cp_after.get("found"),
            "giant steps": int(cp_after.get("counter giant steps", -1)),
            "plan fingerprint": cp_after.get("plan fingerprint"),
        }
        state["sessions"].append(session)
        atomic_json(Path(artifact["artifact path"]), artifact)

        if session["killed"]:
            if cp_after is None:
                raise RuntimeError(
                    f"{mode} timeout occurred before first checkpoint"
                )
            if cp_after.get("found") not in ("0", "1"):
                raise RuntimeError(
                    f"{mode} interrupted session has invalid found state"
                )
            continue
        if session["returncode"] != 0:
            raise RuntimeError(
                f"{mode} final-target run failed rc={session['returncode']}\n"
                f"{session['stderr tail']}\n{session['stdout tail']}"
            )
        if session["recovered"] != f"{scalar:064x}":
            raise RuntimeError(f"{mode} last target recovered wrong scalar")
        if cp_after is None or cp_after.get("found") != "1":
            raise RuntimeError(f"{mode} final checkpoint is missing found=1")
        if int(cp_after.get("chunks completed", -1)) != COVER_CHUNKS:
            raise RuntimeError(f"{mode} final checkpoint did not finish cover")
        state["done"] = True
        break
    else:
        raise RuntimeError(f"{mode} exceeded the 40-session resume limit")

    sessions = state["sessions"]
    chain_ok = True
    prior = None
    for session in sessions:
        if prior is not None:
            resumed_at = session["resumed chunks before"]
            no_op_resume = (session["already recovered"]
                            and resumed_at is None
                            and prior == COVER_CHUNKS)
            if resumed_at != prior and not no_op_resume:
                chain_ok = False
        cp_info = session.get("checkpoint after")
        if cp_info is not None:
            prior = cp_info["chunks completed"]
    state["killed sessions"] = sum(s["killed"] for s in sessions)
    state["resume chain pass"] = chain_ok
    state["pass"] = (chain_ok and state["killed sessions"] > 0
                      and sessions[-1]["recovered"] == f"{scalar:064x}")
    if not state["pass"]:
        raise RuntimeError(f"{mode} kill/resume chain validation failed")
    atomic_json(Path(artifact["artifact path"]), artifact)
    print(f"[{mode} last resume] PASS, {len(sessions)} sessions, "
          f"{state['killed sessions']} killed")


def expected_miss_steps(mode):
    stride = stride_for(mode)
    override = is_mirror(mode)
    per_chunk = (CHUNK - 1) // stride + 1 + int(override)
    return (COVER_CHUNKS + 1) * per_chunk


def run_miss(binary, workdir, mode, stride, table, threads, repeat,
             artifact):
    state = artifact["miss runs"].setdefault(mode, [])
    if len(state) >= repeat:
        return
    scalar = START + COVER_SIZE + 12_345
    command = cover_command(binary, workdir, mode, stride, table, scalar,
                            threads, progress=5)
    result = run_process(command, binary.parent, timeout=900)
    planned = result["counters"].get("giant steps planned")
    actual = result["counters"].get("giant steps")
    no_partner_recovery = result["counters"].get("partner recoveries", 0) == 0
    boundary_exercised = (
        mode != "mirror" or
        (result["counters"].get("search hits", 0) > 0 and
         result["counters"].get("exact candidate checks", 0) > 0)
    )
    expected_steps = expected_miss_steps(mode)
    ok = (result["returncode"] == 1 and result["recovered"] is None
          and planned == expected_steps and actual == expected_steps
          and no_partner_recovery and boundary_exercised)
    record = {
        "pass": ok, "expected giant steps": expected_steps,
        "no partner recovery": no_partner_recovery,
        "out-of-cover boundary exercised": boundary_exercised,
        "run": result,
    }
    state.append(record)
    atomic_json(Path(artifact["artifact path"]), artifact)
    if not ok:
        raise RuntimeError(
            f"{mode} complete miss failed: rc={result['returncode']}, "
            f"planned={planned}, actual={actual}, expected={expected_steps}\n"
            f"{result['stdout tail']}\n{result['stderr tail']}"
        )
    print(f"[{mode} full miss {len(state)}/{repeat}] PASS, "
          f"{actual:,} giants, {result['wall seconds']:.2f}s")


def summarize(artifact):
    summaries = {}
    for mode in ("canonical", "mirror"):
        runs = artifact["miss runs"].get(mode, [])
        times = [r["run"]["wall seconds"] for r in runs if r.get("pass")]
        rates = [r["run"]["counters"].get("giant steps/second")
                 for r in runs if r.get("pass") and
                 r["run"]["counters"].get("giant steps/second") is not None]
        summaries[mode] = {
            "miss repeats": len(times),
            "median miss wall seconds": statistics.median(times) if times else None,
            "median giant steps per second": statistics.median(rates) if rates else None,
            "planned giant steps per miss": expected_miss_steps(mode),
        }
    if summaries["canonical"]["median miss wall seconds"] and \
            summaries["mirror"]["median miss wall seconds"]:
        summaries["paired wall-time speedup"] = (
            summaries["canonical"]["median miss wall seconds"] /
            summaries["mirror"]["median miss wall seconds"]
        )
        summaries["giant-step reduction factor"] = (
            summaries["canonical"]["planned giant steps per miss"] /
            summaries["mirror"]["planned giant steps per miss"]
        )
    artifact["summary"] = summaries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--workdir", type=Path,
                        default=Path("/tmp/ecdsa_research_mirror4m"))
    parser.add_argument("--out", type=Path,
                        default=ROOT / "results/mirror_tiling_4m_comparison.json")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--kill-after", type=float, default=4.0)
    parser.add_argument("--miss-repeats", type=int, default=2)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.exists():
        raise SystemExit(f"missing engine: {binary}")
    if args.threads < 1 or args.miss_repeats < 1:
        raise SystemExit("threads and miss-repeats must be positive")
    args.workdir.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    config = {
        "rows": ROWS, "default stride": DEFAULT_STRIDE,
        "mirror stride": MIRROR_STRIDE, "chunk scalar width": CHUNK,
        "cover chunks": COVER_CHUNKS, "cover size": COVER_SIZE,
        "threads": args.threads, "kill timeout seconds": args.kill_after,
        "miss repeats": args.miss_repeats,
        "engine sha256": sha256_file(binary),
    }
    input_meta = {mode: {"path": str(path), "bytes": path.stat().st_size,
                         "sha256": sha256_file(path)}
                  for mode, path in INPUTS.items()}
    if args.out.exists():
        artifact = json.loads(args.out.read_text())
        if artifact.get("config") != config or artifact.get("inputs") != input_meta:
            raise SystemExit("existing result artifact has different inputs/config")
    else:
        artifact = {
            "status": "running", "config": config, "inputs": input_meta,
            "binary": str(binary), "workdir": str(args.workdir.resolve()),
            "artifact path": str(args.out.resolve()),
            "tables": {}, "proofs": {}, "resume": {}, "miss runs": {},
        }
        atomic_json(args.out, artifact)

    # Both tables use the same 65-byte CBT1 row layout and 4M rows.
    for mode, input_path in INPUTS.items():
        if mode not in artifact["tables"]:
            artifact["tables"][mode] = build_table(
                binary, args.workdir, mode, input_path, args.threads,
                START + MIRROR_STRIDE - 1_234,
            )
            atomic_json(args.out, artifact)
            print(f"[{mode} table] {artifact['tables'][mode]['bytes']:,} bytes, "
                  f"sha256 {artifact['tables'][mode]['sha256'][:16]}")

    # A positive mirror row in chunk 0 and a negative mirror row in chunk 1.
    offsets = {
        "positive": 100_000 * MIRROR_STRIDE + 6_000_000,
        "negative": 520_000 * MIRROR_STRIDE - 6_000_000,
    }
    if mirror_representation(offsets["positive"])["sign"] != "+":
        raise RuntimeError("positive target construction is invalid")
    if mirror_representation(offsets["negative"])["sign"] != "-":
        raise RuntimeError("negative target construction is invalid")

    for mode in ("canonical", "mirror"):
        stride = stride_for(mode)
        table = table_path(args.workdir, mode)
        for label, offset in offsets.items():
            if artifact["proofs"].get(f"{mode}:{label}", {}).get("pass"):
                continue
            run_proof(binary, args.workdir, mode, stride, table, offset,
                      args.threads, artifact, args.out, label)

    for mode in ("canonical", "mirror"):
        run_last_resume(binary, args.workdir, mode, stride_for(mode),
                        table_path(args.workdir, mode), args.threads,
                        args.kill_after, artifact)

    # Alternate layouts on each repetition to reduce order and thermal bias.
    for _ in range(args.miss_repeats):
        for mode in ("canonical", "mirror"):
            run_miss(binary, args.workdir, mode, stride_for(mode),
                     table_path(args.workdir, mode), args.threads,
                     args.miss_repeats, artifact)

    summarize(artifact)
    artifact["status"] = "complete"
    atomic_json(args.out, artifact)
    print("\nsummary:")
    print(json.dumps(artifact["summary"], indent=2))
    print(f"[artifact] {args.out.resolve()}")


if __name__ == "__main__":
    main()
