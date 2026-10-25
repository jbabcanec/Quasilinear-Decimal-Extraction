#!/usr/bin/env python3
"""
odd_prod_lift.py -- pull OEIS A395199 + PE592 ideas into the dual-product attack.

Facts used
----------
(1) OEIS A395199: for k>=2,
      (2^{k+1}-1)!!  ≡  (2^k-1)!!^2   (mod  2^{3k-1})
    so odd double-factorials lift by squaring (O(log) big-muls), after a
    short base range.

(2) Dual identity (already validated):
      prod_{k=1}^K (R + mu_k) ≡ P0 * (1 + R H)  (mod R^2),  R=2^m.

(3) Product of ALL odds < R is (R-1)!! ; our P_R is the 1-mod-8 subproduct
    of odds in (R, R+8K].  Filter via four residue classes mod 8.

Goal: soft-linear (in m) evaluation of (2^m-1)!! mod 2^{2m}, then of the
1-mod-8 subproduct over a SHORT window of length Theta(m), then H.
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


def v2(n):
    n = int(n)
    return 999 if n == 0 else (n & -n).bit_length() - 1


def odd_double_fact_naive(k, mod_bits):
    """(2^k - 1)!! = 1*3*...* (2^k - 1)  mod 2^{mod_bits}."""
    MOD = mpz(1) << mod_bits
    p = mpz(1)
    lim = 1 << k
    for j in range(1, lim, 2):
        p = (p * j) % MOD
    return p


def odd_double_fact_lift(n, mod_bits=None):
    """
    Compute (2^n - 1)!! mod 2^{mod_bits} via OEIS squaring lift.
    Base: naive up to k0 ~ mod_bits/3; then square-lift.
    """
    if mod_bits is None:
        mod_bits = 2 * n
    MOD = mpz(1) << mod_bits
    # Need base depth so that 3*k-1 >= mod_bits for the congruence modulus
    # When squaring from k to k+1: congruence mod 2^{3k-1}.
    # Choose k0 so 3*k0 - 1 >= mod_bits, i.e. k0 >= (mod_bits+1)/3
    k0 = max(3, (mod_bits + 1 + 2) // 3)
    if k0 >= n:
        return odd_double_fact_naive(n, mod_bits)

    # naive base at k0
    cur = odd_double_fact_naive(k0, mod_bits)
    # Also need the intermediate (2^{k}-1)!! for k=k0..n-1 via squaring,
    # BUT the OEIS relation is only (2^{k+1}-1)!! ≡ ((2^k-1)!!)^2 mod 2^{3k-1}.
    # For k large enough that 3k-1 >= mod_bits, squaring mod MOD is exact for us.
    for k in range(k0, n):
        # When 3k-1 < mod_bits, squaring alone is NOT enough — need the
        # missing odd factors in (2^k , 2^{k+1}) .
        if 3 * k - 1 >= mod_bits:
            cur = (cur * cur) % MOD
        else:
            # must multiply the new odds in (2^k, 2^{k+1})
            # (2^{k+1}-1)!! = (2^k-1)!! * prod_{j=2^{k-1}}^{2^k - 1} (2j+1)
            #               = (2^k-1)!! * prod_{odd o in (2^k, 2^{k+1})} o
            pr = mpz(1)
            for o in range((1 << k) + 1, (1 << (k + 1)), 2):
                pr = (pr * o) % MOD
            cur = (cur * pr) % MOD
    return cur


def part_A_validate_lift():
    P("=" * 70)
    P("PART A: OEIS squaring lift for (2^n-1)!!")
    for n, mb in [(6, 12), (8, 16), (10, 20), (12, 24), (14, 28)]:
        a = odd_double_fact_lift(n, mb)
        b = odd_double_fact_naive(n, mb)
        P(f"  n={n:2d} mod_bits={mb:2d}  MATCH={a==b}")
        assert a == b


def part_B_timing():
    P()
    P("=" * 70)
    P("PART B: lift timing vs naive (naive dies exponentially)")
    for n in [8, 10, 12, 14, 16, 18]:
        mb = 2 * n
        t0 = time.perf_counter()
        odd_double_fact_lift(n, mb)
        dt_lift = time.perf_counter() - t0
        if n <= 14:
            t0 = time.perf_counter()
            odd_double_fact_naive(n, mb)
            dt_naive = time.perf_counter() - t0
            P(f"  n={n:2d}  lift={dt_lift:.4f}s  naive={dt_naive:.4f}s  "
              f"speedup={dt_naive/max(dt_lift,1e-9):.1f}x")
        else:
            P(f"  n={n:2d}  lift={dt_lift:.4f}s  naive=skipped (2^{n} odds)")


def part_C_window_1mod8():
    """
    Our P_R = prod_{k=1}^K (R+8k+1) is a SHORT window (K=Theta(m)) of
    1-mod-8 numbers just above R, not the full (R-1)!!.

    Full odds in (R, R+8K]: length 4K odds approx; we want 1/4 of them.
    """
    P()
    P("=" * 70)
    P("PART C: short window 1-mod-8 product above R=2^m")
    for m, K in [(16, 8), (20, 12), (24, 16), (32, 20)]:
        R = mpz(1) << m
        MOD = R * R
        # direct
        pr = mpz(1)
        for k in range(1, K + 1):
            pr = (pr * (R + 8 * k + 1)) % MOD
        # 1-mod-8 window: R+9, R+17, ..., R+8K+1  (exclude R+1)
        pr2 = mpz(1)
        for o in range(R + 9, R + 8 * K + 2, 8):
            pr2 = (pr2 * o) % MOD
        ok = pr == pr2
        P(f"  m={m} K={K}  window filter MATCH={ok}")
        assert ok

        # Dual extract H
        P0 = mpz(1)
        for k in range(1, K + 1):
            P0 = (P0 * (8 * k + 1)) % MOD
        u = (pr * pow(P0, -1, MOD)) % MOD
        H = ((u - 1) // R) % R
        Ht = sum(pow(mpz(8 * k + 1), -1, R) for k in range(1, K + 1)) % R
        P(f"           H extract MATCH={H==Ht}")


def part_D_log_faulhaber_short_product():
    """
    PE592-style: prod_{k=1}^K (1+8k) via exp sum log, with Faulhaber —
    that's P0, already Borwein-cracked.  Here: try prod (R+1+8k) =
    prod R(1 + (1+8k)/R) = R^K * prod (1 + mu/R).

    Mod R^2: R^K = 0 for K>=2.  Wrong scaling — factors are R+mu not R(1+...).
    Stick with dual: need sum mu^{-1}.

    Try computing sum_{k=1}^K Log(R+mu_k) via Faulhaber expansion of Log(R(1+mu/R))
    = K Log R + sum Log(1+mu/R) — Log R not 2-adic for R=2^m.
    """
    P()
    P("=" * 70)
    P("PART D: Log(R+mu) expansion (Iwasawa) — structural note")
    P("  Log(R+mu) is NOT Iwasawa-Log of an odd unit in the usual sense")
    P("  when we only have mod R^2; the dual product is the right language.")
    P("  Faulhaber-Log route for P0 is Borwein.  For P_R it restates H.")


def part_E_residue_class_via_four_odd_prods():
    """
    prod_{j<=X, j=1 mod 8} j  =  ?
    Using: for X=R+8K+1, compare four double-factorial-like products.

    Let f_r(X) = prod_{1<=j<=X, j≡r mod 8} j.
    Then f_1(R+8K+1) / f_1(R+1) = our product (approximately).

    Can f_r(2^m) be lifted like (2^m-1)!!?
    """
    P()
    P("=" * 70)
    P("PART E: residue-class products at powers of 2")

    def f_r_naive(X, r, mod_bits):
        MOD = mpz(1) << mod_bits
        p = mpz(1)
        for j in range(r, X + 1, 8):
            p = (p * j) % MOD
        return p

    # Compare f_1(2^m) growth / pattern
    P("  f_1(2^m) mod 2^{2m}:")
    for m in range(4, 12):
        mb = 2 * m
        v = f_r_naive(1 << m, 1, mb)
        P(f"    m={m}  f_1={hex(v)[:28]}...  bits_nonzero={v.bit_length()}")

    # Squaring relation? f_1(2^{m+1}) vs f_1(2^m)^2 * (product of new terms)
    P("  check f_1(2^{m+1}) / related to f_1(2^m):")
    for m in range(4, 10):
        mb = 2 * (m + 1) + 4
        a = f_r_naive(1 << m, 1, mb)
        b = f_r_naive(1 << (m + 1), 1, mb)
        # new terms: 2^m+1, 2^m+9, ..., 2^{m+1}-7
        pr = mpz(1)
        MOD = mpz(1) << mb
        for j in range((1 << m) + 1, (1 << (m + 1)), 8):
            pr = (pr * j) % MOD
        predict = (a * pr) % MOD
        P(f"    m={m}: b==a*new? {b==predict}")


def main():
    P("odd_prod_lift -- OEIS/PE592 pull into dual-product barrier")
    part_A_validate_lift()
    part_B_timing()
    part_C_window_1mod8()
    part_D_log_faulhaber_short_product()
    part_E_residue_class_via_four_odd_prods()
    P()
    P("NEXT: short-window product of K=Theta(m) terms (R+8k+1) still looks")
    P("like Theta(K) muls unless a squaring/FFT structure appears for")
    P("windows of length m sitting at height 2^m.  Full (2^m-1)!! DOES lift.")
    P("Gap = full-odds lift  -/-  short 1-mod-8 window at height 2^m.")


if __name__ == "__main__":
    main()
