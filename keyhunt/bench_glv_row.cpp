// Micro-benchmark of the per-row GLV build components, linking the engine's
// own primitives. Measures: EC walk addition (add_safe), K3 = x^3, scalar-side
// orbit rotations (2x ModMulK1order vs incremental ModAdd), the x-only
// Lambda triple (3 cross-ratios = 3 ModInv), and the x-only rotations
// (2x ModMulK1 field multiplications).
#include <cstdio>
#include <chrono>
#include "secp256k1/Int.h"
#include "secp256k1/IntGroup.h"
#include "secp256k1/Point.h"
#include "secp256k1/SECP256k1.h"

static Secp256K1 curve;
static Int beta;
static Int beta_g;
static Int glv_lambda;
static Int glv_lambda2;

// mirror of engine lambda_value_for_x (x-only cross-ratio)
Int lambda_value_for_x(Int x) {
  if (x.IsZero()) return Int(0);
  Int beta_x;
  beta_x.ModMulK1(&beta, &x);
  Int left_g; left_g.ModSub(&curve.G.x, &beta_g);
  Int left_x; left_x.ModSub(&x, &beta_x);
  Int numerator; numerator.ModMulK1(&left_g, &left_x);
  Int right_g; right_g.ModSub(&curve.G.x, &beta_x);
  Int right_x; right_x.ModSub(&x, &beta_g);
  Int denominator; denominator.ModMulK1(&right_g, &right_x);
  if (denominator.IsZero()) return Int(0);
  denominator.ModInv();
  numerator.ModMulK1(&denominator);
  return numerator;
}

int main() {
  curve.Init();
  beta.SetBase16("7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE");
  beta_g.ModMulK1(&beta, &curve.G.x);
  glv_lambda.SetBase16("5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72");
  glv_lambda2.SetBase16("AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE");

  const int N = 200000;
  Point current = curve.G;
  Int lam1; lam1.SetInt32(0);
  Int lam2; lam2.SetInt32(0);
  volatile uint64_t sink = 0;

  auto t0 = std::chrono::steady_clock::now();
  for (int i = 1; i <= N; ++i) {
    current = curve.AddDirect(current, curve.G);
  }
  auto t1 = std::chrono::steady_clock::now();
  printf("EC walk add:            %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);

  Int x = current.x;
  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < N; ++i) {
    Int cube; cube.ModCube(&x);
    sink += cube.GetInt64() & 0xffffff;
  }
  t1 = std::chrono::steady_clock::now();
  printf("K3 = x^3 (ModCube):     %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);

  Int s(&x);
  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < N; ++i) {
    Int r1; r1.Set(&glv_lambda); r1.ModMulK1order(&s);
    Int r2; r2.Set(&glv_lambda2); r2.ModMulK1order(&s);
    sink += r1.GetInt64() & 0xffffff + r2.GetInt64() & 0xffffff;
  }
  t1 = std::chrono::steady_clock::now();
  printf("scalar rotations 2xMMK1order: %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);

  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < N; ++i) {
    lam1.ModAdd(&glv_lambda); lam1.Mod(&curve.order);
    lam2.ModAdd(&glv_lambda2); lam2.Mod(&curve.order);
    sink += lam1.GetInt64() & 0xffffff + lam2.GetInt64() & 0xffffff;
  }
  t1 = std::chrono::steady_clock::now();
  printf("scalar rotations incremental: %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);

  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < N; ++i) {
    Int x1; x1.ModMulK1(&beta, &x);
    Int x2; x2.ModMulK1(&beta, &x1);
    sink += x1.GetInt64() & 0xffffff + x2.GetInt64() & 0xffffff;
  }
  t1 = std::chrono::steady_clock::now();
  printf("x-only rotations 2xMMK1: %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);

  t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < N; ++i) {
    Int x1; x1.ModMulK1(&beta, &x);
    Int x2; x2.ModMulK1(&beta, &x1);
    Int v0 = lambda_value_for_x(x);
    Int v1 = lambda_value_for_x(x1);
    Int v2 = lambda_value_for_x(x2);
    sink += v0.GetInt64() & 0xffffff + v1.GetInt64() & 0xffffff + v2.GetInt64() & 0xffffff;
  }
  t1 = std::chrono::steady_clock::now();
  printf("Lambda triple (3 ModInv): %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);

  // batch-inverted Lambda triple over a 64-row chunk: 1 ModInv + mults
  const int CHUNK = 64;
  t0 = std::chrono::steady_clock::now();
  for (int c = 0; c < N / CHUNK; ++c) {
    Int denom[CHUNK];
    for (int i = 0; i < CHUNK; ++i) {
      Int x1; x1.ModMulK1(&beta, &x);
      Int x2; x2.ModMulK1(&beta, &x1);
      // denominators of the three cross-ratios, as in lambda_value_for_x
      Int beta_x; beta_x.ModMulK1(&beta, &x);
      Int d0; d0.ModSub(&curve.G.x, &beta_x); d0.ModMulK1(&d0, &x);
      // NOTE: simplified denom only for timing shape; correctness path
      // stays in the engine.
      denom[i] = d0;
      sink += d0.GetInt64() & 0xffffff;
    }
    Int acc; acc.SetInt32(1);
    for (int i = 0; i < CHUNK; ++i) { acc.ModMulK1(&denom[i]); }
    acc.ModInv();
    Int inv; inv.Set(&acc);
    for (int i = CHUNK - 1; i >= 0; --i) {
      Int out; out.ModMulK1(&inv, &denom[i]);
      inv.ModMulK1(&denom[i]);
      sink += out.GetInt64() & 0xffffff;
    }
  }
  t1 = std::chrono::steady_clock::now();
  printf("batch-invert denom (64): %.3f us/row\n",
         std::chrono::duration<double, std::micro>(t1 - t0).count() / N);
  return 0;
}
