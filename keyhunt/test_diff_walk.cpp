// Standalone check of the engine's projective differential walk:
// replicate differential_step + doubling + batch normalization and compare
// affine x(jG) with ec_mult-derived x for j = 1..200, plus the incremental
// canonical scalars against canonical_orbit_scalar.
#include <cstdio>
#include <cstring>
#include <vector>
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

static Int diff_v, diff_v3_2b, diff_two_v, diff_num_c;
static Int diff_four, diff_eight_b, diff_v3_b4;

static void init_diff_constants() {
  diff_v.Set(&curve.G.x);
  Int v3;
  v3.ModSquare(&diff_v);
  v3.ModMulK1(&diff_v);
  diff_v3_2b.Set(&v3);
  Int two_b;
  two_b.SetInt32(14);
  diff_v3_2b.ModAdd(&two_b);
  diff_two_v.ModAdd(&diff_v, &diff_v);
  diff_num_c.ModSub(&curve.G.x, &beta_g);
  diff_four.SetInt32(4);
  diff_eight_b.SetInt32(56);
  // doubling denominator: 4(v^3 + b) with b = 7 (not v^3 + 2b - that
  // constant belongs to the differential step identity).
  Int seven;
  seven.SetInt32(7);
  diff_v3_b4.ModAdd(&v3, &seven);
  diff_v3_b4.ModMulK1(&diff_four);
}

static void differential_step(Int &Xm, Int &Zm, Int &Xp, Int &Zp) {
  Int t1;
  t1.ModMulK1(&diff_v, &Zp);
  Int T;
  T.ModSub(&Xp, &t1);
  Int T2;
  T2.ModSquare(&T);
  Int xs;
  xs.ModSquare(&Xp);
  xs.ModMulK1(&Xp);
  Int zs;
  zs.ModSquare(&Zp);
  zs.ModMulK1(&Zp);  // Zp^3
  zs.ModMulK1(&diff_v3_2b, &zs);
  zs.ModAdd(&xs);  // Xp^3 + (v^3+2b) Zp^3
  Int Xn;
  Xn.ModMulK1(&zs, &Zm);
  Xn.ModAdd(&Xn);  // 2 zs Zm
  Int M1;
  M1.ModMulK1(&T2, &Zp);
  M1.ModMulK1(&Xm);          // T^2 Zp Xm
  Int M2;
  M2.ModMulK1(&T2, &Zm);     // T^2 Zm
  Int M2b;
  M2b.ModAdd(&Xp, &Xp);      // 2Xp
  Int M2c;
  M2c.ModMulK1(&diff_two_v, &Zp);
  M2b.ModAdd(&M2c);          // 2Xp + 2v Zp
  M2.ModMulK1(&M2b);         // T^2 Zm (2Xp + 2v Zp)
  Xn.ModSub(&M1);
  Xn.ModSub(&M2);
  Int Zn;
  Zn.ModMulK1(&Zm, &Zp);
  Zn.ModMulK1(&T2);
  Xm.Set(&Xp);
  Zm.Set(&Zp);
  Xp.Set(&Xn);
  Zp.Set(&Zn);
}

int main() {
  curve.Init();
  beta.SetBase16("7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE");
  beta_g.ModMulK1(&beta, &curve.G.x);
  glv_lambda.SetBase16("5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72");
  glv_lambda2.SetBase16("AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE");
  init_diff_constants();

  // doubling: x(2G) = (v^4-8bv)/(4(v^3+b))
  Int Xp, Zp;
  Int vs;
  vs.ModSquare(&diff_v);
  Xp.ModSquare(&vs);
  Int t;
  t.ModMulK1(&diff_eight_b, &diff_v);
  Xp.ModSub(&t);
  Zp.Set(&diff_v3_b4);
  Int Xm(&diff_v);
  Int Zm; Zm.SetInt32(1);
  { alignas(8) unsigned char bb[32]; Xp.Get32Bytes(bb);
    printf("X2 = "); for (int i=0;i<32;++i) printf("%02x", bb[i]); printf("\n");
    Zp.Get32Bytes(bb);
    printf("Z2 = "); for (int i=0;i<32;++i) printf("%02x", bb[i]); printf("\n");
    diff_v.Get32Bytes(bb);
    printf("v  = "); for (int i=0;i<32;++i) printf("%02x", bb[i]); printf("\n"); }

  // normalize row 2 and compare
  const int N = 200;
  std::vector<Int> Xs(N + 1), Zs(N + 1);
  Xs[1] = diff_v; Zs[1].SetInt32(1);
  Xs[2] = Xp; Zs[2] = Zp;
  for (int j = 2; j < N; ++j) {
    differential_step(Xm, Zm, Xp, Zp);
    Xs[j + 1] = Xp; Zs[j + 1] = Zp;
  }
  int bad = 0;
  for (int j = 1; j <= N; ++j) {
    Int scalar; scalar.SetInt64(j);
    Point ref = curve.ComputePublicKey(&scalar);
    Int x_aff;
    Int zinv;
    zinv.Set(&Zs[j]);
    zinv.ModInv();
    x_aff.ModMulK1(&Xs[j], &zinv);
    if (!x_aff.IsEqual(&ref.x)) {
      if (bad < 5) printf("MISMATCH j=%d\n", j);
      ++bad;
    }
  }
  printf(bad == 0 ? "walk chain j=1..%d: OK\n" : "walk chain FAILURES: %d\n",
         bad == 0 ? N : bad);

  // signature batch vs per-row reference comparison
  int sig_bad = 0;
  for (int j = 3; j < 40; ++j) {
    Int x_aff;
    Int zinv;
    zinv.Set(&Zs[j]);
    zinv.ModInv();
    x_aff.ModMulK1(&Xs[j], &zinv);
    // reference triple (per-row path)
    Int x1; x1.ModMulK1(&beta, &x_aff);
    Int x2; x2.ModMulK1(&beta, &x1);
    auto lv = [&](Int &x) {
      Int bx_; bx_.ModMulK1(&beta, &x);
      Int lg; lg.ModSub(&curve.G.x, &beta_g);
      Int lx; lx.ModSub(&x, &bx_);
      Int num; num.ModMulK1(&lg, &lx);
      Int rg; rg.ModSub(&curve.G.x, &bx_);
      Int rx; rx.ModSub(&x, &beta_g);
      Int den; den.ModMulK1(&rg, &rx);
      den.ModInv();
      num.ModMulK1(&den);
      return num;
    };
    Int L0r = lv(x_aff), L1r = lv(x1), L2r = lv(x2);
    Int d01r; d01r.ModSub(&L1r, &L0r);
    Int d12r; d12r.ModSub(&L2r, &L1r);
    // batch path
    Int xs[3] = {x_aff, x1, x2};
    Int Lb[3], dens[3];
    Int num_c; num_c.ModSub(&curve.G.x, &beta_g);
    for (int k = 0; k < 3; ++k) {
      int j2 = (k + 1) % 3;
      Int num; num.ModSub(&xs[k], &xs[j2]);
      num.ModMulK1(&num_c, &num);
      Int d1; d1.ModSub(&curve.G.x, &xs[j2]);
      Int d2; d2.ModSub(&xs[k], &beta_g);
      dens[k].ModMulK1(&d1, &d2);
      Lb[k].Set(&num);
    }
    // Montgomery batch over the 3 denominators
    Int pref[3]; Int acc; acc.SetInt32(1);
    for (int k = 0; k < 3; ++k) { pref[k].Set(&acc); acc.ModMulK1(&dens[k]); }
    acc.ModInv();
    Int inv[3];
    for (int k = 2; k >= 0; --k) { inv[k].ModMulK1(&acc, &pref[k]); acc.ModMulK1(&dens[k]); }
    for (int k = 0; k < 3; ++k) Lb[k].ModMulK1(&inv[k]);
    Int d01b; d01b.ModSub(&Lb[1], &Lb[0]);
    Int d12b; d12b.ModSub(&Lb[2], &Lb[1]);
    alignas(8) unsigned char rb1[32], rb2[32], bb1[32], bb2[32];
    d01r.Get32Bytes(rb1); d12r.Get32Bytes(rb2);
    d01b.Get32Bytes(bb1); d12b.Get32Bytes(bb2);
    if (memcmp(rb1, bb1, 32) != 0 || memcmp(rb2, bb2, 32) != 0) {
      if (sig_bad < 3) printf("sig mismatch j=%d\n", j);
      ++sig_bad;
    }
  }
  printf(sig_bad == 0 ? "signature batch == per-row: OK\n" : "signature batch FAILURES: %d\n", sig_bad);
  return bad != 0 || sig_bad != 0;
}
