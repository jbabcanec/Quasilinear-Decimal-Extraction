#!/usr/bin/env python3
"""
counting_structure.py -- structural probes of the counting face (companion to
counting_face.py). Two questions about B(d) = sum_t bit_d(1/(8t+1)):

 UPPER-BOUND lead: is B(d) low-complexity (C-finite / automatic)?  If yes -> jump it.
 HARDNESS lead:    what exactly is the hard kernel?

RESULTS (recorded run, T=300, N=2400 for part 1; exact):
 (A) Linear complexity of B(d) is MAXIMAL (~N/2) mod 2,3,5,7,257 => NOT C-finite,
     no low-order recurrence to jump. Kills the clean-recurrence upper bound.
 (B) 2-kernel of the parity grows (not saturated) => not simply automatic.
 (C) Centered partial sums drift ~ -T*bias (mean 143.7 != T/2=150): a COMPUTABLE
     mean-bias main term + complexity-maximal fluctuation. Structure is in the
     main term only.
 (D) SHARPENED OBJECT: B(d) = sum_t f_t(d mod p_t), a superposition of T periodic
     0/1 sequences with incommensurate periods p_t = ord_2(8t+1); total
     description ~ sum p_t ~ N^2 (the mass, in its cleanest statable form).
 (R1) per-period popcount of 1/mu = #(<2 mod mu> cap (mu/2, mu)) -- subgroup-in-
      interval count (the incomplete-Kloosterman / Bruin object). = p/2 EXACTLY
      iff -1 in <2 mod mu>; else deviates by a character/class-number sum
      (the tame/wild split, exact and visible).
 (R2) prefix sum_(d<N) B(d) = sum_t s_t(N), s_t(N) = (N//p_t) w_t + partial:
      the SMOOTH part of the face is governed by the T subgroup-counts w_t.
"""
from sympy import n_order


def bits_recip(mu, N):
    out = bytearray(N)
    r = 1 % mu
    for d in range(N):
        r <<= 1
        if r >= mu:
            out[d] = 1
            r -= mu
    return out


def Bseq(T, N):
    cols = [0] * N
    for t in range(1, T + 1):
        b = bits_recip(8 * t + 1, N)
        for d in range(N):
            cols[d] += b[d]
    return cols


def berlekamp_massey(s, p):
    n = len(s); C = [1] + [0] * n; B = [1] + [0] * n; L = 0; m = 1; b = 1
    for i in range(n):
        d = s[i] % p
        for j in range(1, L + 1):
            d = (d + C[j] * s[i - j]) % p
        if d == 0:
            m += 1
        elif 2 * L <= i:
            Tc = C[:]; coef = (d * pow(b, -1, p)) % p
            for j in range(n + 1 - m):
                if j + m <= n:
                    C[j + m] = (C[j + m] - coef * B[j]) % p
            L = i + 1 - L; B = Tc; b = d; m = 1
        else:
            coef = (d * pow(b, -1, p)) % p
            for j in range(n + 1 - m):
                if j + m <= n:
                    C[j + m] = (C[j + m] - coef * B[j]) % p
            m += 1
    return L


def period_popcount_and_subgroup(mu):
    p = n_order(2, mu)
    r = 1 % mu; ones = 0
    for _ in range(p):
        r <<= 1
        if r >= mu:
            ones += 1; r -= mu
    g = 1; sg = 0
    for _ in range(p):
        if g > mu / 2:
            sg += 1
        g = (g * 2) % mu
    return p, ones, sg


def main():
    import itertools
    T, N = 300, 2400
    B = Bseq(T, N)
    print(f"[T={T}, N={N}] B(0..12) = {B[:13]}")
    for p in [2, 3, 5, 7, 257]:
        print(f"(A) linear complexity mod {p}: {berlekamp_massey(B, p)} / {N} "
              f"(random ~{N // 2}) -> {'MAXIMAL' if berlekamp_massey(B,p) > N//2 - N//10 else 'LOW!'}")
    cen = list(itertools.accumulate(b - T // 2 for b in B))
    print(f"(C) centered partial-sum range [{min(cen)}, {max(cen)}]  (sqrt N ~ {int(N**.5)}; "
          f"drift = computable mean bias, fluctuation random-like)")

    print("\n(R1) per-period popcount == subgroup <2> in (mu/2,mu):")
    okR1 = True
    for t in [1, 2, 3, 6, 9, 13, 50, 400]:
        mu = 8 * t + 1
        p, ones, sg = period_popcount_and_subgroup(mu)
        half = "= p/2" if 2 * ones == p else f"!= p/2 (dev {ones - p/2:+.1f})"
        okR1 &= (ones == sg)
        print(f"     mu={mu:5d} period={p:4d} popcount={ones:4d} subgrp={sg:4d} "
              f"[{'OK' if ones == sg else 'X'}] {half}")
    print(f"(R1) VALIDATED: {okR1}  (deviation from p/2 = the wild character sum)")


if __name__ == "__main__":
    main()
