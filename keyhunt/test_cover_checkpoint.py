#!/usr/bin/env python3
"""Black-box верификация checkpoint/resume для --cover-start/--cover-size.

Дизайн под убивающую среду (PROGRESS_MEMO.md, следующий рабочий шаг):
длинный прогон должен переживать SIGKILL в середине внешнего chunk и
продолжаться с атомарного checkpoint, а не с начала плана.

Проверки:

  [1] kill-тест: движок убивается SIGKILL после завершения chunk 0
      ("checkpoint saved" в stdout); процесс умирает, target не найден,
      checkpoint-файл цел (атомарная запись tmp+rename).
  [2] resume: тот же план с --resume продолжает с chunk 1, находит
      контрольную цель (засеянную в chunk 1), печатает корректный ключ
      и не переделывает chunk 0 (automatic chunks attempted == 1).
  [3] повторный --resume после успеха: цель уже в checkpoint - движок
      печатает найденный ключ из checkpoint и не ищет заново
      (правило "прогон не считается успешным без >=1 найденной цели").
  [4] negative control: --resume того же checkpoint против ДРУГОГО плана
      (другой --size) отклоняется по fingerprint, код 2.

Харнесс передаёт движку только public key/start/size - private key
используется лишь для сверки внутри тестового процесса.
"""

import argparse
import os
import re
import signal
import subprocess
import time
from pathlib import Path

from cryptography.hazmat.primitives.asymmetric import ec

START = 1 << 134

FAILURES = []


def pubkey_hex(scalar):
    numbers = ec.derive_private_key(scalar, ec.SECP256K1()).public_key().public_numbers()
    return f"{2 | (numbers.y & 1):02x}{numbers.x:064x}"


def parse_harness(stdout):
    match = re.search(r"recovered private key: 0x([0-9a-fA-F]{64})", stdout)
    return int(match.group(1), 16) if match else None


def kill_after_checkpoint(cmd, workdir):
    """Запускает движок и убивает его SIGKILL вскоре после checkpoint chunk 0.

    Первый session живёт под 5-секундным потолком по наибольшему из
    (первый "checkpoint saved", старт chunk 1) - в тестовом chunk ~2s
    giant-работы, поэтому окно kill заведомо шире задержки пула.
    """
    proc = subprocess.Popen(cmd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, cwd=str(workdir))
    sig_saved = re.compile(r"checkpoint saved: chunks (\d+)")
    chunks_at_kill = None
    t0 = time.time()
    while True:
        line = proc.stdout.readline()
        if not line:
            break
        m = sig_saved.search(line)
        if m:
            chunks_at_kill = int(m.group(1))
            break
        if time.time() - t0 > 300:
            proc.send_signal(signal.SIGKILL)
            raise RuntimeError("kill window missed: no checkpoint line in 300s")
    proc.send_signal(signal.SIGKILL)
    stdout, _ = proc.communicate()
    return parse_harness(stdout), chunks_at_kill, proc


def check(name, ok, detail=""):
    print(f"[{name}] {'PASS' if ok else 'FAIL'}{(' - ' + detail) if detail else ''}")
    if not ok:
        FAILURES.append(name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="./eng_bsgs")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--tmp", default=None)
    args = parser.parse_args()

    binary = Path(args.binary).resolve()

    import tempfile
    tmp = Path(args.tmp) if args.tmp else Path(tempfile.mkdtemp(prefix="cover_cp_"))
    cp = tmp / "cover_checkpoint.txt"

    # план: chunk = 1e12 кандидатов, m = 1e6 строк, cover = 3 chunk.
    # контрольная цель - внутри chunk 1 (i=1), чтобы kill в chunk 1
    # гарантированно не находил её, а resume с chunk 1 - находил быстро.
    chunk = 10 ** 12
    target_scalar = START + chunk + 123_456
    target_pub = pubkey_hex(target_scalar)
    base_cmd = [
        str(binary),
        "--target", target_pub,
        "--start", hex(START),
        "--size", str(chunk),
        "--baby-rows", "1000000",
        "--threads", str(args.threads),
        "--cover-start", hex(START),
        "--cover-size", hex(3 * chunk),
        "--checkpoint-file", str(cp),
        "--checkpoint-steps", "1",
    ]

    # ---- [1] kill в составе run 1 ----
    recovered, chunks_at_kill, _ = kill_after_checkpoint(base_cmd, tmp)
    check("kill", recovered is None and chunks_at_kill == 1,
          f"recovered={'yes' if recovered is not None else 'no'}, "
          f"chunks saved={chunks_at_kill}")
    cp_exists = cp.is_file()
    stale_tmp = Path(str(cp) + ".tmp").exists()
    check("checkpoint file atomic", cp_exists and not stale_tmp,
          f"exists={cp_exists}, stale_tmp={stale_tmp}")
    if cp_exists:
        fields = cp.read_text()
        digest = re.search(r"chunks completed: (\d+)", fields)
        found_flag = re.search(r"^found: (\d+)$", fields, re.M)
        check("checkpoint content",
              digest and digest.group(1) == "1" and found_flag and found_flag.group(1) == "0",
              f"chunks_completed={digest.group(1) if digest else '?'}, "
              f"found={found_flag.group(1) if found_flag else '?'}")

    # ---- [2] resume с checkpoint ----
    t0 = time.perf_counter()
    proc_resumed = subprocess.run(base_cmd + ["--resume"], text=True,
                                  capture_output=True, cwd=str(tmp))
    resumed_out = proc_resumed.stdout
    resumed_recovered = parse_harness(resumed_out)
    resumed_correct = resumed_recovered == target_scalar
    resumed_chunks = re.search(r"automatic chunks attempted: (\d+)", resumed_out)
    resumed_marker = "checkpoint resumed: chunks completed before this session: 1" in resumed_out
    check("resume recovery", resumed_correct and proc_resumed.returncode == 0,
          f"exit={proc_resumed.returncode}, "
          f"recovered={'yes' if resumed_correct else 'no'}")
    check("resume continuity", resumed_marker and resumed_chunks
          and resumed_chunks.group(1) == "1",
          f"marker={'yes' if resumed_marker else 'no'}, "
          f"chunks_attempted={resumed_chunks.group(1) if resumed_chunks else '?'}")
    check("resume accounting",
          "time to solution seconds:" in resumed_out
          and "exact candidate checks per hit:" in resumed_out,
          "split counters/tts in output")
    print(f"    resume finished in {time.perf_counter() - t0:.1f}s")

    # ---- [3] повторный resume: цель уже найдена ----
    proc_again = subprocess.run(base_cmd + ["--resume"], text=True,
                                capture_output=True, cwd=str(tmp))
    again_recovered = parse_harness(proc_again.stdout)
    check("resume after success", again_recovered == target_scalar
          and proc_again.returncode == 0
          and "earlier session" in proc_again.stdout,
          f"exit={proc_again.returncode}")

    # ---- [4] negative control: другой план противтого же checkpoint ----
    proc_wrong = subprocess.run(
        [str(binary), "--target", target_pub, "--start", hex(START),
         "--size", str(2 * chunk),  # другой план
         "--baby-rows", "1000000", "--threads", str(args.threads),
         "--cover-start", hex(START), "--cover-size", hex(6 * chunk),
         "--checkpoint-file", str(cp), "--resume"],
        text=True, capture_output=True, cwd=str(tmp))
    check("fingerprint guard", proc_wrong.returncode == 2
          and "belongs to a different plan" in proc_wrong.stderr,
          f"exit={proc_wrong.returncode}")

    # ---- очистка ----
    for f in (cp, Path(str(cp) + ".tmp")):
        if f.exists():
            f.unlink()
    tmp.rmdir()

    if FAILURES:
        print(f"\nFAILURES: {len(FAILURES)}: {', '.join(FAILURES)}")
        raise SystemExit(1)
    print("\nALL PASS")


if __name__ == "__main__":
    main()
