#include "int2048.h"

#include <algorithm>
#include <cstdint>
#include <functional>

namespace sjtu {

// Base 10^9 representation
static const long long BASE = 1000000000LL;
static const int BASE_DIGITS = 9;

// ---------- Helpers ----------
static void trim(std::vector<long long> &d) {
  while (d.size() > 1 && d.back() == 0) d.pop_back();
}

static bool isZero(const std::vector<long long> &d) {
  return d.size() == 1 && d[0] == 0;
}

// Compare absolute values: -1 if |a|<|b|, 0 if equal, +1 if |a|>|b|
static int absCmp(const std::vector<long long> &a,
                  const std::vector<long long> &b) {
  if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
  for (int i = (int)a.size() - 1; i >= 0; --i) {
    if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
  }
  return 0;
}

// a += b (unsigned digit vectors)
static void absAdd(std::vector<long long> &a, const std::vector<long long> &b) {
  if (a.size() < b.size()) a.resize(b.size(), 0);
  long long carry = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    long long v = a[i] + carry + (i < b.size() ? b[i] : 0);
    if (v >= BASE) {
      carry = 1;
      v -= BASE;
    } else {
      carry = 0;
    }
    a[i] = v;
  }
  if (carry) a.push_back(carry);
}

// a -= b, assuming |a| >= |b|
static void absSub(std::vector<long long> &a, const std::vector<long long> &b) {
  long long borrow = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    long long v = a[i] - borrow - (i < b.size() ? b[i] : 0);
    if (v < 0) {
      v += BASE;
      borrow = 1;
    } else {
      borrow = 0;
    }
    a[i] = v;
  }
  trim(a);
}

// ---------- NTT multiplication ----------
static const unsigned long long MOD1 = 998244353ULL;   // 119*2^23+1, g=3
static const unsigned long long MOD2 = 985661441ULL;   // 235*2^22+1, g=3
static const unsigned long long MOD3 = 754974721ULL;   // 45*2^24+1,  g=11
static const unsigned long long G = 3;
static const unsigned long long G3 = 11;

static unsigned long long pw(unsigned long long a, unsigned long long b,
                             unsigned long long m) {
  unsigned long long r = 1 % m;
  a %= m;
  while (b) {
    if (b & 1) r = r * a % m;
    a = a * a % m;
    b >>= 1;
  }
  return r;
}

static void ntt(std::vector<unsigned long long> &a, bool inv,
                unsigned long long mod, unsigned long long g) {
  int n = (int)a.size();
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (int len = 2; len <= n; len <<= 1) {
    unsigned long long w = pw(g, (mod - 1) / len, mod);
    if (inv) w = pw(w, mod - 2, mod);
    int half = len >> 1;
    std::vector<unsigned long long> roots(half);
    roots[0] = 1;
    for (int k = 1; k < half; ++k) roots[k] = roots[k - 1] * w % mod;
    for (int i = 0; i < n; i += len) {
      for (int j = 0; j < half; ++j) {
        unsigned long long u = a[i + j];
        unsigned long long v = a[i + j + half] * roots[j] % mod;
        a[i + j] = (u + v) % mod;
        a[i + j + half] = (u + mod - v) % mod;
      }
    }
  }
  if (inv) {
    unsigned long long nInv = pw((unsigned long long)n, mod - 2, mod);
    for (int i = 0; i < n; ++i) a[i] = a[i] * nInv % mod;
  }
}

// Multiply two non-zero absolute digit vectors.
static std::vector<long long> mulAbs(const std::vector<long long> &A,
                                     const std::vector<long long> &B) {
  if (isZero(A) || isZero(B)) return {0};

  size_t na = A.size(), nb = B.size();

  // Schoolbook for small inputs.
  if (na <= 64 || nb <= 64) {
    std::vector<long long> C(na + nb, 0);
    for (size_t i = 0; i < na; ++i) {
      unsigned long long carry = 0;
      unsigned long long ai = (unsigned long long)A[i];
      size_t j = 0;
      for (; j < nb; ++j) {
        unsigned long long cur =
            (unsigned long long)C[i + j] + ai * (unsigned long long)B[j] + carry;
        C[i + j] = (long long)(cur % BASE);
        carry = cur / BASE;
      }
      size_t k = i + nb;
      while (carry) {
        unsigned long long cur = (unsigned long long)C[k] + carry;
        C[k] = (long long)(cur % BASE);
        carry = cur / BASE;
        ++k;
      }
    }
    trim(C);
    return C;
  }

  // NTT with base-10^9 digits directly; use 3 primes and CRT.
  // Coefficients in convolution bounded by: digit^2 * n <= (10^9)^2 * n.
  // For n up to 2^18, this is ~2^77. Three 30-bit primes give ~2^90 modulus via CRT — safe.
  size_t result_len = na + nb;
  size_t n = 1;
  while (n < result_len) n <<= 1;

  std::vector<unsigned long long> fa(n, 0), fb(n, 0);
  for (size_t i = 0; i < na; ++i) fa[i] = (unsigned long long)A[i];
  for (size_t i = 0; i < nb; ++i) fb[i] = (unsigned long long)B[i];

  // NTT mod MOD1
  std::vector<unsigned long long> a1 = fa, b1 = fb;
  ntt(a1, false, MOD1, G);
  ntt(b1, false, MOD1, G);
  for (size_t i = 0; i < n; ++i) a1[i] = a1[i] * b1[i] % MOD1;
  ntt(a1, true, MOD1, G);
  { std::vector<unsigned long long> tmp; tmp.swap(b1); }

  // NTT mod MOD2
  std::vector<unsigned long long> a2 = fa, b2 = fb;
  ntt(a2, false, MOD2, G);
  ntt(b2, false, MOD2, G);
  for (size_t i = 0; i < n; ++i) a2[i] = a2[i] * b2[i] % MOD2;
  ntt(a2, true, MOD2, G);
  { std::vector<unsigned long long> tmp; tmp.swap(b2); }

  // NTT mod MOD3
  std::vector<unsigned long long> a3 = std::move(fa);
  std::vector<unsigned long long> b3 = std::move(fb);
  ntt(a3, false, MOD3, G3);
  ntt(b3, false, MOD3, G3);
  for (size_t i = 0; i < n; ++i) a3[i] = a3[i] * b3[i] % MOD3;
  ntt(a3, true, MOD3, G3);
  { std::vector<unsigned long long> tmp; tmp.swap(b3); }

  // CRT three primes.
  // x mod MOD1 = r1; x mod MOD2 = r2; x mod MOD3 = r3.
  // Step 1: combine r1, r2 => x12 mod MOD1*MOD2 using u64 (MOD1*MOD2 ~ 2^60).
  // Step 2: combine (x12, r3) to get actual value in two 64-bit parts (hi, lo) or as 128-bit.
  unsigned long long inv_m1_m2 = pw(MOD1 % MOD2, MOD2 - 2, MOD2);
  unsigned long long m12_mod_m3 = MOD1 % MOD3 * (MOD2 % MOD3) % MOD3;
  unsigned long long inv_m12_m3 = pw(m12_mod_m3, MOD3 - 2, MOD3);

  // We carry into the digit array incrementally to avoid materializing huge coefficients.
  // Each coefficient is at most ~ (10^9)^2 * n ~ 2^77. We keep running carry as __int128.
  std::vector<long long> C;
  C.reserve(result_len + 2);
  __uint128_t carry = 0;
  const __uint128_t BASE128 = (__uint128_t)BASE;

  for (size_t i = 0; i < n; ++i) {
    unsigned long long r1 = a1[i];
    unsigned long long r2 = a2[i];
    unsigned long long r3 = a3[i];
    // x12 = r1 + MOD1 * k2, where k2 = (r2 - r1) * inv(MOD1) mod MOD2
    unsigned long long k2 = (r2 + MOD2 - r1 % MOD2) % MOD2 * inv_m1_m2 % MOD2;
    // x12 fits in u64 (< MOD1*MOD2 < 2^60)
    unsigned long long x12 = r1 + MOD1 * k2;
    // x = x12 + MOD1 * MOD2 * k3, where k3 = (r3 - x12 mod MOD3) * inv(MOD1*MOD2 mod MOD3) mod MOD3
    unsigned long long x12_mod_m3 = x12 % MOD3;
    unsigned long long k3 =
        (r3 + MOD3 - x12_mod_m3) % MOD3 * inv_m12_m3 % MOD3;
    // x = x12 + (MOD1*MOD2) * k3
    // Compute as 128-bit
    __uint128_t m12 = (__uint128_t)MOD1 * MOD2;
    __uint128_t x = (__uint128_t)x12 + m12 * k3;

    __uint128_t v = x + carry;
    C.push_back((long long)(v % BASE128));
    carry = v / BASE128;
  }
  while (carry) {
    C.push_back((long long)(carry % BASE128));
    carry /= BASE128;
  }
  trim(C);
  return C;
}

// Divide A by a small positive integer d. Returns quotient; remainder in rem.
static std::vector<long long> divSmall(const std::vector<long long> &A,
                                       long long d, long long &rem) {
  std::vector<long long> C(A.size(), 0);
  unsigned long long r = 0;
  for (int i = (int)A.size() - 1; i >= 0; --i) {
    unsigned long long cur = r * (unsigned long long)BASE + (unsigned long long)A[i];
    C[i] = (long long)(cur / (unsigned long long)d);
    r = cur % (unsigned long long)d;
  }
  rem = (long long)r;
  trim(C);
  return C;
}

// Shift-left by k base-10^9 digits (multiply by BASE^k).
static std::vector<long long> shiftLeft(const std::vector<long long> &A, int k) {
  if (isZero(A) || k == 0) return A;
  std::vector<long long> R(A.size() + k, 0);
  for (size_t i = 0; i < A.size(); ++i) R[i + k] = A[i];
  return R;
}

// Return top k digits of A (most-significant) as a new vector (little-endian).
// If A has fewer than k digits, return A itself.
static std::vector<long long> topDigits(const std::vector<long long> &A, int k) {
  int n = (int)A.size();
  if (k >= n) return A;
  std::vector<long long> R(A.begin() + (n - k), A.end());
  trim(R);
  return R;
}

// Schoolbook long division of A by B (absolute values). A.size() and B.size() small.
// Quotient and remainder returned.
static std::vector<long long> divAbsSchool(const std::vector<long long> &A,
                                           const std::vector<long long> &B,
                                           std::vector<long long> &R);

// Compute reciprocal of B with precision d: returns floor(BASE^d / B).
// Uses Newton's iteration.
static std::vector<long long> reciprocal(const std::vector<long long> &B, int d);

// Divide A by B where both are absolute. Returns quotient Q, remainder R satisfies A = Q*B + R, 0<=R<B.
static std::vector<long long> divAbs(const std::vector<long long> &A,
                                     const std::vector<long long> &B,
                                     std::vector<long long> &R) {
  int cmp = absCmp(A, B);
  if (cmp < 0) {
    R = A;
    return {0};
  }
  if (cmp == 0) {
    R = {0};
    return {1};
  }
  int m = (int)B.size();
  if (m == 1) {
    long long rem;
    std::vector<long long> Q = divSmall(A, B[0], rem);
    R = {rem};
    return Q;
  }
  // Threshold for schoolbook
  int n = (int)A.size();
  if (n - m <= 32 || m <= 32) {
    return divAbsSchool(A, B, R);
  }

  // Newton's method: compute reciprocal of B with precision d = n (number of base digits in A).
  int d = n;  // We want quotient length at most n - m + 1.
  std::vector<long long> inv = reciprocal(B, d);
  // Q ~ floor(A * inv / BASE^d)
  std::vector<long long> prod = mulAbs(A, inv);
  // Divide prod by BASE^d (shift right by d digits)
  std::vector<long long> Q;
  if ((int)prod.size() <= d) {
    Q = {0};
  } else {
    Q.assign(prod.begin() + d, prod.end());
    trim(Q);
  }
  // Verify and correct Q. Compute R = A - Q*B.
  auto computeR = [&](const std::vector<long long> &Qv,
                      std::vector<long long> &Rv) -> int {
    // Rv = A - Q*B; allow negative representation.
    std::vector<long long> QB = mulAbs(Qv, B);
    int c = absCmp(A, QB);
    if (c >= 0) {
      Rv = A;
      absSub(Rv, QB);
      return 1;  // non-negative R = A - QB
    } else {
      Rv = QB;
      absSub(Rv, A);
      return -1;  // R is -(QB - A) i.e. negative
    }
  };

  std::vector<long long> Rv;
  int sgn = computeR(Q, Rv);
  // Adjust Q up or down until 0 <= R < B.
  std::vector<long long> one = {1};
  int iters = 0;
  while (iters < 64) {
    if (sgn < 0) {
      // A - Q*B < 0 => Q too large; decrement
      absSub(Q, one);
      if (isZero(Q)) {
        // Q = 0
      }
      sgn = computeR(Q, Rv);
    } else {
      // R >= 0. Check R < B.
      if (absCmp(Rv, B) < 0) break;
      // R >= B => Q too small
      absAdd(Q, one);
      sgn = computeR(Q, Rv);
    }
    ++iters;
  }
  R = Rv;
  if (sgn < 0) R = {0};  // shouldn't happen
  trim(Q);
  return Q;
}

// Knuth's Algorithm D (TAOCP, Vol 2, §4.3.1) in base BASE = 10^9.
static std::vector<long long> divAbsSchool(const std::vector<long long> &A,
                                           const std::vector<long long> &B,
                                           std::vector<long long> &R) {
  int m0 = (int)B.size();
  // Single-digit divisor
  if (m0 == 1) {
    long long rem;
    std::vector<long long> Q = divSmall(A, B[0], rem);
    R = {rem};
    return Q;
  }
  // Normalize so B's top digit is >= BASE/2.
  long long d = BASE / (B.back() + 1);
  std::vector<long long> u, v;
  // u = A * d
  {
    u.assign(A.size() + 1, 0);
    unsigned long long carry = 0;
    for (size_t i = 0; i < A.size(); ++i) {
      unsigned long long cur = (unsigned long long)A[i] * (unsigned long long)d + carry;
      u[i] = (long long)(cur % BASE);
      carry = cur / BASE;
    }
    u[A.size()] = (long long)carry;
  }
  // v = B * d
  {
    v.assign(B.size(), 0);
    unsigned long long carry = 0;
    for (size_t i = 0; i < B.size(); ++i) {
      unsigned long long cur = (unsigned long long)B[i] * (unsigned long long)d + carry;
      v[i] = (long long)(cur % BASE);
      carry = cur / BASE;
    }
    // v.back() should now be >= BASE/2; no overflow by choice of d.
  }
  int n = (int)v.size();
  int m = (int)u.size() - n;  // number of quotient digits
  std::vector<long long> q(m, 0);

  for (int j = m - 1; j >= 0; --j) {
    // Step D3: Calculate qhat
    unsigned long long top = (unsigned long long)u[j + n] * (unsigned long long)BASE +
                             (unsigned long long)u[j + n - 1];
    unsigned long long qhat = top / (unsigned long long)v[n - 1];
    unsigned long long rhat = top % (unsigned long long)v[n - 1];
    if (qhat >= (unsigned long long)BASE) {
      qhat = BASE - 1;
      rhat = top - qhat * (unsigned long long)v[n - 1];
    }
    // Refine
    while (true) {
      unsigned long long lhs = qhat * (unsigned long long)v[n - 2];
      unsigned long long rhs = (unsigned long long)BASE * rhat +
                               (unsigned long long)u[j + n - 2];
      // We need lhs <= rhs; if not, decrement qhat.
      if (lhs > rhs) {
        qhat--;
        rhat += (unsigned long long)v[n - 1];
        if (rhat >= (unsigned long long)BASE) break;
      } else {
        break;
      }
    }

    // Step D4: Multiply and subtract. Compute u[j..j+n] -= qhat * v.
    long long borrow = 0;
    for (int i = 0; i < n; ++i) {
      unsigned long long p = qhat * (unsigned long long)v[i];
      long long plo = (long long)(p % (unsigned long long)BASE);
      long long phi = (long long)(p / (unsigned long long)BASE);
      long long cur = u[j + i] - plo - borrow;
      if (cur < 0) {
        // bring up by BASE
        long long k = (-cur + BASE - 1) / BASE;
        cur += k * BASE;
        borrow = phi + k;
      } else {
        borrow = phi;
      }
      u[j + i] = cur;
    }
    long long top_cur = u[j + n] - borrow;

    // Step D5: Test remainder
    if (top_cur < 0) {
      // Step D6: Add back
      qhat--;
      long long addCarry = 0;
      for (int i = 0; i < n; ++i) {
        long long vv = u[j + i] + v[i] + addCarry;
        if (vv >= BASE) {
          vv -= BASE;
          addCarry = 1;
        } else {
          addCarry = 0;
        }
        u[j + i] = vv;
      }
      top_cur += addCarry;
    }
    u[j + n] = top_cur;
    q[j] = (long long)qhat;
  }

  trim(q);
  // Denormalize remainder
  std::vector<long long> rem(u.begin(), u.begin() + n);
  trim(rem);
  long long dummy;
  R = divSmall(rem, d, dummy);
  trim(R);
  return q;
}

// Compute floor(BASE^d / B) using Newton's iteration on truncated/top-digits of B.
// The idea: compute reciprocal of the leading digits of B at increasing precision.
//
// For a divisor B of length m (base-10^9 digits), we first compute a "small" reciprocal
// using schoolbook Algorithm D, then iteratively refine with:
//   new_x = x * (2 * BASE^p - B_top * x / BASE^p_old) / BASE^(p_old - (p - p_old))
// But we use an equivalent simpler formulation that doubles precision each step.
static std::vector<long long> reciprocal(const std::vector<long long> &B, int d) {
  int m = (int)B.size();
  if (d - m <= 20) {
    std::vector<long long> numer(d + 1, 0);
    numer[d] = 1;
    std::vector<long long> remDummy;
    return divAbsSchool(numer, B, remDummy);
  }

  // Start with base reciprocal computed via schoolbook.
  // We compute x_0 = floor(BASE^(2*m_small) / B_top) where B_top is top digits of B.
  // Then refine precision doubling.
  int base_k = 32;  // initial precision width in base-10^9 digits
  if (base_k > d) base_k = d;

  // x approximates BASE^k / B_{[top k digits]}, we maintain invariant:
  //   x is a k_x-digit (roughly) number such that x ~= floor(BASE^k / B).
  // We iterate: given x with precision k, compute x_new with precision 2k
  // using formula: x_new = 2*x*BASE^k - B_{top 2k+m digits of B} * x^2 shifted appropriately.
  //
  // Cleanest implementation: at precision k, we have inv satisfying
  //   inv = floor(BASE^(k+m) / B)   (approximately)
  // where B is truncated/used in full. When k >= d - m, we're done.

  auto truncB = [&](int len) -> std::vector<long long> {
    if (len >= m) return B;
    std::vector<long long> R(B.end() - len, B.end());
    trim(R);
    return R;
  };

  // Initial x: use schoolbook to get inv at small precision k0.
  int k = base_k;
  std::vector<long long> x;
  {
    // Compute floor(BASE^(k+m) / B_top(k+m)) using schoolbook.
    int Btop_len = std::min(m, k + m);
    std::vector<long long> Btop = truncB(Btop_len);
    std::vector<long long> numer(k + m + 1, 0);
    numer[k + m] = 1;
    std::vector<long long> remDummy;
    x = divAbsSchool(numer, Btop, remDummy);
    // x has size approximately k+1
  }

  // Now iteratively double k until k >= d - m.
  int target_k = d - m;
  while (k < target_k) {
    int new_k = std::min(target_k, 2 * k);
    // Newton step: given x is "reciprocal of B" at precision k (i.e., x ~= BASE^(k+m) / B),
    // compute x' at precision new_k.
    // The formula: x' = 2*x*BASE^(new_k - k) - floor(B * x^2 / BASE^(2k + m))
    //   where the result x' ~= BASE^(new_k + m) / B.
    // But we truncate B to (new_k + m) digits and x^2 correspondingly.

    // Compute x^2
    std::vector<long long> xsq = mulAbs(x, x);
    // Use full B (safer for correctness).
    std::vector<long long> Bxsq = mulAbs(B, xsq);
    // shift right by (2k + m - new_k). Derivation:
    // x ~ BASE^(k+m)/B. y_scaled = x * BASE^(new_k-k) ~ BASE^(new_k+m)/B.
    // Newton: y_new = 2*y_scaled - B*y_scaled^2 / BASE^(new_k+m).
    //       = 2 x BASE^(new_k-k) - B x^2 BASE^(2(new_k-k)) / BASE^(new_k+m)
    //       = 2 x BASE^(new_k-k) - B x^2 / BASE^(2k + m - new_k).
    int shift_right = 2 * k + m - new_k;
    std::vector<long long> Bxsq_shifted;
    if ((int)Bxsq.size() <= shift_right) {
      Bxsq_shifted = {0};
    } else {
      Bxsq_shifted.assign(Bxsq.begin() + shift_right, Bxsq.end());
      trim(Bxsq_shifted);
    }
    // 2*x*BASE^(new_k - k)
    std::vector<long long> twoX = shiftLeft(x, new_k - k);
    {
      long long carry = 0;
      for (size_t i = 0; i < twoX.size(); ++i) {
        long long v = twoX[i] * 2 + carry;
        if (v >= BASE) {
          carry = v / BASE;
          v %= BASE;
        } else {
          carry = 0;
        }
        twoX[i] = v;
      }
      if (carry) twoX.push_back(carry);
    }

    std::vector<long long> xnew;
    if (absCmp(twoX, Bxsq_shifted) >= 0) {
      xnew = twoX;
      absSub(xnew, Bxsq_shifted);
    } else {
      xnew = {0};
    }
    x = xnew;
    k = new_k;
  }

  // Now x ~= BASE^(d) / B (since k + m == d). Actually we had k such that
  // x ~= BASE^(k+m)/B, so k+m = d means k = d - m.
  // If k > d - m, we need to shift right by (k - (d-m)).
  // If k < d - m, shouldn't happen after our loop since we target d-m.
  if (k > target_k) {
    int sh = k - target_k;
    std::vector<long long> R;
    if ((int)x.size() <= sh) {
      x = {0};
    } else {
      R.assign(x.begin() + sh, x.end());
      trim(R);
      x = R;
    }
  }
  // Now x is approximately floor(BASE^d / B). Apply correction.
  std::vector<long long> base_d(d + 1, 0);
  base_d[d] = 1;
  std::vector<long long> one = {1};
  // Make sure B*x <= base_d; otherwise decrement.
  for (int t = 0; t < 8; ++t) {
    std::vector<long long> Bx = mulAbs(B, x);
    if (absCmp(Bx, base_d) > 0) {
      if (isZero(x)) break;
      absSub(x, one);
    } else {
      // check (x+1)*B <= base_d
      std::vector<long long> xp = x;
      absAdd(xp, one);
      std::vector<long long> Bxp = mulAbs(B, xp);
      if (absCmp(Bxp, base_d) <= 0) {
        x = xp;
      } else {
        break;
      }
    }
  }
  return x;
}

// ---------- Constructors ----------
int2048::int2048() : sign(0), digits(1, 0) {}

int2048::int2048(long long x) {
  if (x == 0) {
    sign = 0;
    digits = {0};
    return;
  }
  if (x < 0) {
    sign = 1;
    // handle LLONG_MIN safely
    unsigned long long ux = (unsigned long long)(-(x + 1)) + 1ULL;
    digits.clear();
    while (ux > 0) {
      digits.push_back((long long)(ux % (unsigned long long)BASE));
      ux /= (unsigned long long)BASE;
    }
  } else {
    sign = 0;
    digits.clear();
    unsigned long long ux = (unsigned long long)x;
    while (ux > 0) {
      digits.push_back((long long)(ux % (unsigned long long)BASE));
      ux /= (unsigned long long)BASE;
    }
  }
}

int2048::int2048(const std::string &s) { read(s); }

int2048::int2048(const int2048 &o) : sign(o.sign), digits(o.digits) {}

// ---------- read / print ----------
void int2048::read(const std::string &s) {
  sign = 0;
  digits.clear();
  int start = 0;
  if (!s.empty() && (s[0] == '-' || s[0] == '+')) {
    if (s[0] == '-') sign = 1;
    start = 1;
  }
  // Skip leading zeros to check for effective length (but keep if all zero).
  while (start < (int)s.size() - 1 && s[start] == '0') start++;

  // Parse from the end in chunks of BASE_DIGITS digits.
  for (int i = (int)s.size(); i > start; i -= BASE_DIGITS) {
    int l = std::max(start, i - BASE_DIGITS);
    long long v = 0;
    for (int k = l; k < i; ++k) v = v * 10 + (s[k] - '0');
    digits.push_back(v);
  }
  if (digits.empty()) digits.push_back(0);
  trim(digits);
  if (isZero(digits)) sign = 0;
}

void int2048::print() {
  if (sign == 1 && !isZero(digits)) std::putchar('-');
  std::printf("%lld", digits.back());
  for (int i = (int)digits.size() - 2; i >= 0; --i) {
    std::printf("%09lld", digits[i]);
  }
}

// ---------- add / minus (Integer1) ----------
int2048 &int2048::add(const int2048 &o) {
  if (sign == o.sign) {
    absAdd(digits, o.digits);
  } else {
    int c = absCmp(digits, o.digits);
    if (c == 0) {
      digits = {0};
      sign = 0;
    } else if (c > 0) {
      absSub(digits, o.digits);
    } else {
      std::vector<long long> tmp = o.digits;
      absSub(tmp, digits);
      digits = std::move(tmp);
      sign = o.sign;
    }
  }
  if (isZero(digits)) sign = 0;
  return *this;
}

int2048 add(int2048 a, const int2048 &b) {
  a.add(b);
  return a;
}

int2048 &int2048::minus(const int2048 &o) {
  // a - b = a + (-b)
  int2048 neg = o;
  if (!isZero(neg.digits)) neg.sign ^= 1;
  return add(neg);
}

int2048 minus(int2048 a, const int2048 &b) {
  a.minus(b);
  return a;
}

// ---------- Unary operators ----------
int2048 int2048::operator+() const { return *this; }

int2048 int2048::operator-() const {
  int2048 r = *this;
  if (!isZero(r.digits)) r.sign ^= 1;
  return r;
}

// ---------- Assignment ----------
int2048 &int2048::operator=(const int2048 &o) {
  if (this != &o) {
    sign = o.sign;
    digits = o.digits;
  }
  return *this;
}

// ---------- +=, -= ----------
int2048 &int2048::operator+=(const int2048 &o) { return add(o); }
int2048 operator+(int2048 a, const int2048 &b) {
  a += b;
  return a;
}

int2048 &int2048::operator-=(const int2048 &o) { return minus(o); }
int2048 operator-(int2048 a, const int2048 &b) {
  a -= b;
  return a;
}

// ---------- *= ----------
int2048 &int2048::operator*=(const int2048 &o) {
  if (isZero(digits) || isZero(o.digits)) {
    sign = 0;
    digits = {0};
    return *this;
  }
  digits = mulAbs(digits, o.digits);
  sign ^= o.sign;
  if (isZero(digits)) sign = 0;
  return *this;
}

int2048 operator*(int2048 a, const int2048 &b) {
  a *= b;
  return a;
}

// ---------- /= and %= ----------
// Floor division toward negative infinity.
int2048 &int2048::operator/=(const int2048 &o) {
  // Compute truncated division on absolute values, then adjust for sign.
  int s1 = sign, s2 = o.sign;
  std::vector<long long> R;
  std::vector<long long> Q = divAbs(digits, o.digits, R);
  // Q and R are absolute values. Quotient sign = s1 ^ s2. If signs differ and remainder != 0,
  // truncated toward zero gives -|Q|. Floor wants -|Q| - 1 and remainder adjusted.
  int qSign = s1 ^ s2;
  digits = Q;
  sign = qSign;
  if (isZero(digits)) sign = 0;
  // Floor adjustment: if signs differ and R != 0, subtract 1 from quotient.
  if (s1 != s2 && !isZero(R)) {
    // digits = -|Q| - 1 ->  |digits| should become |Q| + 1, sign stays negative.
    // If Q was zero: digits = {0}, sign = 0; we want -1.
    std::vector<long long> one = {1};
    absAdd(digits, one);
    sign = 1;
  }
  if (isZero(digits)) sign = 0;
  return *this;
}

int2048 operator/(int2048 a, const int2048 &b) {
  a /= b;
  return a;
}

int2048 &int2048::operator%=(const int2048 &o) {
  // r = a - (a / b) * b
  int2048 q = *this;
  q /= o;
  int2048 prod = q;
  prod *= o;
  *this -= prod;
  return *this;
}

int2048 operator%(int2048 a, const int2048 &b) {
  a %= b;
  return a;
}

// ---------- Stream ops ----------
std::istream &operator>>(std::istream &in, int2048 &x) {
  std::string s;
  in >> s;
  x.read(s);
  return in;
}

std::ostream &operator<<(std::ostream &out, const int2048 &x) {
  if (x.sign == 1 && !(x.digits.size() == 1 && x.digits[0] == 0)) out << '-';
  out << x.digits.back();
  char buf[16];
  for (int i = (int)x.digits.size() - 2; i >= 0; --i) {
    std::snprintf(buf, sizeof(buf), "%09lld", x.digits[i]);
    out << buf;
  }
  return out;
}

// ---------- Comparisons ----------
bool operator==(const int2048 &a, const int2048 &b) {
  return a.sign == b.sign && a.digits == b.digits;
}

bool operator!=(const int2048 &a, const int2048 &b) { return !(a == b); }

bool operator<(const int2048 &a, const int2048 &b) {
  bool aZero = (a.digits.size() == 1 && a.digits[0] == 0);
  bool bZero = (b.digits.size() == 1 && b.digits[0] == 0);
  int sa = aZero ? 0 : a.sign;
  int sb = bZero ? 0 : b.sign;
  if (sa != sb) return sa > sb;  // negative (sign=1) is smaller
  int c = absCmp(a.digits, b.digits);
  if (sa == 0) return c < 0;
  return c > 0;
}

bool operator>(const int2048 &a, const int2048 &b) { return b < a; }
bool operator<=(const int2048 &a, const int2048 &b) { return !(b < a); }
bool operator>=(const int2048 &a, const int2048 &b) { return !(a < b); }

}  // namespace sjtu
