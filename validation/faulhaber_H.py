#!/usr/bin/env python3
"""
faulhaber_H.py -- H via geometric series + power sums (no Bernoulli table).

    H = sum_{k=1}^K 1/(8k+1)
      = sum_{k=1}^K sum_{j>=0} (-8k)^j
      = sum_{j>=0} (-8)^j * S_j(K),   S_j = sum_{k=1}^K k^j

2-adically v2(term_j) grows ~ 3j - v2(S_j), truncate at j ~ m/3.

Power sums for K,m = Theta(N) by repeated multiply (no Bernoulli):
  cost O((m/3) * K * M(m)/balance) = O(N^2 log N) — worse than merge,
  but CORRECT and avoids inverse-merge.  Useful as an independent engine
  and as a scaffold for trying soft-linear power-sum techniques.

Also records the OEIS new_block squaring observation for the LONG
(dyadic-length) 1-mod-8 product — wrong regime for digit-K, but real structure.
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


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


def power_sums(K, J, m):
    """S[j] = sum_{k=1}^K k^j mod 2^m for j=0..J."""
    MOD = mpz(1) << m
    S = [mpz(0)] * (J + 1)
    for k in range(1, K + 1):
        p = mpz(1)
        for j in range(J + 1):
            S[j] = (S[j] + p) % MOD
            p = (p * k) % MOD
    return S


def H_faulhaber(K, m):
    J = (m + 20) // 3 + 2
    # terms need precision m + some guard for intermediate (-8)^j
    guard = m + 2 * J + 16
    S = power_sums(K, J, guard)
    MOD = mpz(1) << guard
    acc = mpz(0)
    pow8 = mpz(1)  # (-8)^j
    for j in range(J + 1):
        term = (pow8 * S[j]) % MOD
        acc = (acc + term) % MOD
        pow8 = (pow8 * (-8)) % MOD
    return acc % (1 << m)


def part_validate():
    P("=" * 70)
    P("VALIDATE Faulhaber-geometric H")
    for K, m in [(16, 48), (32, 64), (64, 96), (128, 120)]:
        a = atom_H(K, m)
        b = H_faulhaber(K, m)
        P(f"  K={K:4d} m={m:4d}  MATCH={a==b}")
        assert a == b


def part_timing():
    P()
    P("=" * 70)
    P("TIMING vs direct inverse-sum (same asymptotic class as sequential dual)")
    for K in [100, 200, 400, 800]:
        m = K + 40
        t0 = time.perf_counter()
        atom_H(K, m)
        d0 = time.perf_counter() - t0
        t0 = time.perf_counter()
        H_faulhaber(K, m)
        d1 = time.perf_counter() - t0
        P(f"  K={K:4d} m={m:4d}  direct={d0:.4f}s  faulhaber={d1:.4f}s")


def part_newblock_squaring_lemma():
    P()
    P("=" * 70)
    P("LEMMA (machine): let L=2^{m-3}, R=2^m,")
    P("  B_m = prod_{t=0}^{L-1} (R+8t+1).")
    P("Then v2( B_{m+1} / B_m^2  -  1 ) = 2m-3.")
    P("(Squaring lift for the LONG dyadic 1-mod-8 block; length exponential")
    P(" in m — not the digit regime K=Theta(m).)")

    def new_block(m, mb):
        MOD = mpz(1) << mb
        R = mpz(1) << m
        p = mpz(1)
        for t in range(1 << (m - 3)):
            p = (p * (R + 8 * t + 1)) % MOD
        return p

    for m in range(6, 14):
        mb = 4 * m
        MOD = mpz(1) << mb
        a = new_block(m, mb)
        b = new_block(m + 1, mb)
        ratio = (b * pow((a * a) % MOD, -1, MOD)) % MOD
        c = v2(ratio - 1)
        expect = 2 * m - 3
        P(f"  m={m:2d}  v2(B_{{m+1}}/B_m^2 - 1)={c:3d}  expect={expect}  OK={c==expect}")
        assert c == expect


def main():
    P("faulhaber_H -- geometric engine + new_block lemma")
    part_validate()
    part_timing()
    part_newblock_squaring_lemma()
    P()
    P("STATUS: Faulhaber engine bit-exact; O(N^2) not O(N log^2).")
    P("new_block squaring is real but wrong length regime.")
    P("Still hunting soft-linear short-window product at height 2^m.")


if __name__ == "__main__":
    main()
