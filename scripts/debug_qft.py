#!/usr/bin/env python3
"""Отладка qft_pow2: |j> -> QFT и косет-состояние против FFT-пути."""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/home/z/my-project/ECDSA/experiments")
from qip_engine import Statevector, coset_state_posterior  # noqa: E402

# 1) QFT против DFT на 3 кубитах
w = 3
for j in (1, 5):
    sv = Statevector(w)
    sv.psi[j] = 1.0
    sv.qft_pow2([0, 1, 2])  # q0=MSB регистра
    amps = sv.psi
    exact = np.array([np.exp(2j * np.pi * j * k / 8) / np.sqrt(8)
                      for k in range(8)])
    err = np.abs(amps - exact).max()
    exact_rev = exact.reshape([2, 2, 2]).transpose(2, 1, 0).ravel()
    err_rev = np.abs(amps - exact_rev).max()
    print(f"j={j}: max|amp-DFT|={err:.4f}, bit-reversed={err_rev:.4f}")
    print("   amps:", np.round(amps, 3))
    print("   dft :", np.round(exact, 3))

# 2) косет n=128: движок (qft_pow2 на готовом косете) против fft2
n = 128
x, z0 = 12, 77
t = np.arange(n)
sv = Statevector(14)
for tv in t:
    a = (z0 - tv * x) % n
    b = tv
    idx = (a << 7) | b  # a на q13..q7, b на q6..q0
    sv.psi[idx] = 1.0 / np.sqrt(n)
sv.qft_pow2(list(range(13, 6, -1)))   # QFT по a-регистру
sv.qft_pow2(list(range(6, -1, -1)))   # QFT по b-регистру
p_gate = np.abs(sv.psi.reshape(n, n)) ** 2
p_fft = np.abs(np.fft.fft2(coset_state_posterior(n, x, z0),
                           norm="ortho")) ** 2

nz_gate = np.argwhere(p_gate > 1e-12)
nz_fft = np.argwhere(p_fft > 1e-12)
print("nonzero gate:", len(nz_gate), "| nonzero fft:", len(nz_fft),
      "| sum:", round(float(p_gate.sum()), 4), round(float(p_fft.sum()), 4))
print("gate on lattice:", all(v == u * x % n for u, v in nz_gate),
      "| fft on lattice:", all(v == u * x % n for u, v in nz_fft))
tg = sorted(((float(p_gate[u, v]), int(u), int(v)) for u, v in nz_gate),
            reverse=True)[:6]
tf = sorted(((float(p_fft[u, v]), int(u), int(v)) for u, v in nz_fft),
            reverse=True)[:6]
print("gate top:", [(round(p_, 4), u, v) for p_, u, v in tg])
print("fft  top:", [(round(p_, 4), u, v) for p_, u, v in tf])
