#!/usr/bin/env python3
"""
batch_Sj_attack.py -- can batch power sums S_j = sum_{k<=K} k^j
for j=0..J, J~m/3, beat O(J*K) / O(J^2)?

Routes:
  (1) Faulhaber + Bernoulli recurrence (classical O(J^2))
  (2) repeated-multiply per k (O(J K))
  (3) Stirling-2nd: k^j = sum S(j,i) * falling(k,i);
      sum_k falling(k,i) = falling(K+1,i+1)/(i+1)  -- closed!
      So S_j = sum_{i=0}^j S(j,i) * (K+1)^{falling i+1}/(i+1)

  If Stirling-2nd row S(j,*) for all j<=J is available in soft-linear,
  batch S_j is soft-linear.  Stirling-2nd generating function:
      sum_{j} S(j,i) x^j / j! = (e^x - 1)^i / i!   -- EGF, 1/j! not in Z_2
      OGF: sum S(j,i) x^j = x^i / prod_{r=1}^i (1 - r x)

  Computing all S(j,i) for i,j <= J is Theta(J^2) table (or soft-J^2
  with FFT triangle algorithms still ~ J^2 / log?).

  (4) Pull: Harvey Bernoulli is O(J^2 log^{2+eps}); no soft-linear known.
"""
from __future__ import annotations
import time

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def power_sums_naive(K, J, m):
    MOD = mpz(1) << m
    S = [mpz(0)] * (J + 1)
    for k in range(1, K + 1):
        p = mpz(1)
        for j in range(J + 1):
            S[j] = (S[j] + p) % MOD
            p = (p * k) % MOD
    return S


def stirling2_triangle(J, MOD):
    """S2[j][i] = Stirling 2nd kind S(j,i) mod MOD, j,i <= J."""
    S2 = [[mpz(0)] * (J + 1) for _ in range(J + 1)]
    S2[0][0] = mpz(1)
    for j in range(1, J + 1):
        for i in range(1, j + 1):
            S2[j][i] = (S2[j - 1][i - 1] + i * S2[j - 1][i]) % MOD
    return S2


def falling_sum(K, i):
    """sum_{k=1}^K falling(k,i) = falling(K+1,i+1)/(i+1) - [i==0]."""
    acc = mpz(1)
    for t in range(i + 1):
        acc *= (K + 1 - t)
    assert acc % (i + 1) == 0
    out = acc // (i + 1)
    if i == 0:
        out -= 1  # drop k=0 term
    return out


def power_sums_stirling(K, J, m):
    """Via Stirling 2nd + closed falling sums."""
    guard = m + J.bit_length() + 8
    GMOD = mpz(1) << guard
    S2 = stirling2_triangle(J, GMOD)
    fall = [falling_sum(K, i) % GMOD for i in range(J + 1)]
    S = [mpz(0)] * (J + 1)
    for j in range(J + 1):
        acc = mpz(0)
        for i in range(j + 1):
            acc = (acc + S2[j][i] * fall[i]) % GMOD
        S[j] = acc % (1 << m)
    return S


def bernoulli_upto(J, MOD):
    """B_0..B_J as rationals mapped into 2-adics when denom odd part inverted.
    Uses Akiyama-Tanigawa / recursive: (n+1) B_n = -sum binom(n+1,k) B_k for k<n
    Actually: sum_{k=0}^m binom(m+1,k) B_k = 0 for m>0, B_0=1.
    B_m = -1/(m+1) sum_{k=0}^{m-1} binom(m+1,k) B_k
    """
    B = [mpz(0)] * (J + 1)
    B[0] = mpz(1)
    # keep numerator/denom separately for 2-adic embed
    # Store B[k] * some common approach: compute as pairs (num, den)
    nums = [mpz(1)] + [mpz(0)] * J
    dens = [mpz(1)] + [mpz(1)] * J
    for m in range(1, J + 1):
        # B_m = -1/(m+1) * sum_{k=0}^{m-1} C(m+1,k) B_k
        s_num, s_den = mpz(0), mpz(1)
        for k in range(m):
            # add C(m+1,k) * nums[k]/dens[k]
            from math import comb
            c = mpz(comb(m + 1, k))
            # s_num/s_den += c * nums[k]/dens[k]
            s_num = s_num * dens[k] + c * nums[k] * s_den
            s_den = s_den * dens[k]
            g = mpz(math_gcd(s_num, s_den))
            s_num //= g
            s_den //= g
        # B_m = -s_num/s_den/(m+1)
        nums[m] = -s_num
        dens[m] = s_den * (m + 1)
        g = mpz(math_gcd(abs(int(nums[m])), int(dens[m])))
        nums[m] //= g
        dens[m] //= g
    # embed into Z/2^m if dens odd
    out = []
    for m in range(J + 1):
        d = dens[m]
        if int(d) % 2 == 0:
            # strip factors of 2 from den into valuation — for Faulhaber need care
            v = 0
            dd = int(d)
            while dd % 2 == 0:
                dd //= 2
                v += 1
            inv = pow(mpz(dd), -1, MOD)
            val = (nums[m] * inv) % MOD
            # actual B has factor 2^{-v}; store as (val, v) — skip for now
            out.append((val, v))
        else:
            inv = pow(d % MOD, -1, MOD)
            out.append(((nums[m] * inv) % MOD, 0))
    return out


def math_gcd(a, b):
    while b:
        a, b = b, a % b
    return abs(a)


def part_stirling_validate():
    P("=" * 70)
    P("PART A: Stirling-2nd batch power sums == naive")
    for K, J, m in [(30, 20, 64), (50, 25, 80)]:
        a = power_sums_naive(K, J, m)
        b = power_sums_stirling(K, J, m)
        ok = all(a[j] == b[j] for j in range(J + 1))
        P(f"  K={K} J={J} m={m}  match={ok}")
        assert ok


def part_complexity_note():
    P("=" * 70)
    P("PART B: timing Stirling triangle O(J^2) vs naive O(JK)")
    for K, J, m in [(500, 80, 96), (1000, 120, 128)]:
        t0 = time.perf_counter()
        power_sums_naive(K, J, m)
        t1 = time.perf_counter()
        power_sums_stirling(K, J, m)
        t2 = time.perf_counter()
        P(f"  K={K} J={J}  naive={t1-t0:.3f}s  stirling={t2-t1:.3f}s")
    P("""
  DIAGNOSIS: Stirling route replaces O(JK) with O(J^2 + J log K) arithmetic,
  which for J~K~N is still Theta(N^2).  Soft-linear Stirling triangle
  (all S(j,i) for j,i<=N) is NOT known below ~N^2;
  EGF (e^x-1)^i / i! blocked in Z_2 by factorials.

  EQUIVALENCE: soft-linear batch S_j  <=>  soft-linear Bernoulli batch
  (Newton/Faulhaber invertible triangular)  <=>  Open Problem E via
  H = sum (-8)^j S_j.
""")


def part_ogf_stirling():
    P("=" * 70)
    P("PART C: OGF for Stirling-2nd column — still quadratic mass")
    P("  sum_j S(j,i) x^j = x^i / prod_{r=1}^i (1-r x)")
    P("  To get all columns i=1..J truncated to deg J: product tree of")
    P("  (1-r x) factors + series invert — cost soft-O(J^2) ring ops.")
    P("  Same quadratic barrier.")


def main():
    P("batch_Sj_attack -- power-sum batch vs Open Problem E")
    part_stirling_validate()
    part_complexity_note()
    part_ogf_stirling()
    P()
    P("VERDICT: batch S_j does not slip under N^2; glued to Bernoulli mass.")


if __name__ == "__main__":
    main()
