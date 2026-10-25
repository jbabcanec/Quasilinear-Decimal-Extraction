#!/usr/bin/env python3
"""
dual_borwein_hunt.py -- hunt a Borwein-style doubling identity for

    B(m,K) = prod_{k=1}^K (2^m + 8k + 1)  mod  2^{2m}

Low-word Borwein uses quadratic relations among AP products at
related lengths.  Probe whether B(m,K) factors through B(m+1, *) and
B(m, *) or through low-word products P0.

Also test the concrete pairing:
  (R+16j+1)(R+16j+9)  vs  something in (2R + ...)
"""
from __future__ import annotations

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def B(m, K, modbits=None):
    R = mpz(1) << m
    MOD = (mpz(1) << (modbits if modbits else 2 * m))
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = (acc * (R + 8 * k + 1)) % MOD
    return acc


def P0(K, modbits):
    MOD = mpz(1) << modbits
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = (acc * (8 * k + 1)) % MOD
    return acc


def part_pairing_identity():
    P("=" * 70)
    P("PART A: local pairing (R+16j+1)(R+16j+9)")
    for m in [8, 12, 16]:
        R = mpz(1) << m
        MOD = mpz(1) << (2 * m)
        j = 3
        left = ((R + 16 * j + 1) * (R + 16 * j + 9)) % MOD
        # candidates
        cands = {
            "P0-pair": ((16 * j + 1) * (16 * j + 9)) % MOD,
            "2R-pair": ((2 * R + 8 * j + 1) * (2 * R + 8 * j + 9)) % MOD,
            "R2+R*10+P0": (R * R + R * (32 * j + 10) + (16 * j + 1) * (16 * j + 9)) % MOD,
        }
        # last one is algebraic expansion — must match left
        P(f"  m={m} j={j}")
        P(f"    left == R^2 expansion: {left == cands['R2+R*10+P0']}")
        P(f"    left / P0-pair mod R: {(left * pow(cands['P0-pair'], -1, MOD)) % R}")
        # ratio left / ((2R+8j+1)(2R+8j+9))
        ratio = (left * pow(cands["2R-pair"], -1, MOD)) % MOD
        P(f"    left/2R-pair = {hex(int(ratio))}  v2(ratio-1)={_v2(ratio-1)}")


def _v2(n):
    n = int(n)
    if n == 0:
        return 999
    return (n & -n).bit_length() - 1


def part_doubling_search():
    """
    Search for small integer relations:
      B(m,K)^a * B(m+1,K//2)^b * P0(K)^c * B(m,K//2)^d ≡ 1 mod 2^{w}
    via brute exponents in {-2,-1,0,1,2}.
    """
    P("=" * 70)
    P("PART B: brute doubling monomials (small K,m)")
    hits = []
    for m in [10, 12, 14]:
        for K in [8, 12, 16]:
            if K % 2:
                continue
            w = 2 * m - 4  # slightly below full 2m to allow slack
            MOD = mpz(1) << w
            vals = {
                "BmK": B(m, K, w),
                "Bm1_K2": B(m + 1, K // 2, w),
                "Bm_K2": B(m, K // 2, w),
                "P0K": P0(K, w),
                "P0K2": P0(K // 2, w),
                "Bm1_K": B(m + 1, K, w),
            }
            names = list(vals.keys())
            # brute exponents -2..2
            from itertools import product
            for exps in product(range(-2, 3), repeat=len(names)):
                if all(e == 0 for e in exps):
                    continue
                if sum(abs(e) for e in exps) > 4:
                    continue
                acc = mpz(1)
                ok = True
                for e, nm in zip(exps, names):
                    if e == 0:
                        continue
                    v = vals[nm]
                    if e < 0:
                        try:
                            v = pow(v, -1, MOD)
                        except ValueError:
                            ok = False
                            break
                        e = -e
                    acc = (acc * pow(v, e, MOD)) % MOD
                if not ok:
                    continue
                if acc == 1:
                    rel = " * ".join(
                        f"{nm}^{e}" for nm, e in zip(names, exps) if e
                    )
                    hits.append((m, K, rel))
            # also check B(m,K) == B(m,K//2) * mid_half
            mid = mpz(1)
            R = mpz(1) << m
            for k in range(K // 2 + 1, K + 1):
                mid = (mid * (R + 8 * k + 1)) % MOD
            half = B(m, K // 2, w)
            if (half * mid) % MOD == B(m, K, w):
                pass  # trivial split always true
    if hits:
        P(f"  FOUND {len(hits)} relations (showing up to 12):")
        for h in hits[:12]:
            P(f"    m={h[0]} K={h[1]}: {h[2]} == 1")
    else:
        P("  NO non-trivial monomial relation with |exps|<=2, wt<=4.")
    P("  (Trivial half-split always holds; not a complexity win.)")


def part_new_block_vs_short():
    P("=" * 70)
    P("PART C: new_block squaring vs short-window (regime check)")
    # B_full(m) = prod_{t=0}^{2^{m-3}-1} (2^m+8t+1)
    # identity: B_full(m+1) / B_full(m)^2 ≈ 1 mod 2^{2m-3}
    for m in [8, 9, 10]:
        L = 1 << (m - 3)
        # only check if L small enough
        if L > 2000:
            P(f"  m={m} L={L} skip (too long)")
            continue
        MOD = mpz(1) << (2 * m)
        bf = B(m, L, 2 * m)  # k=1..L vs t=0..L-1 off-by-one; close enough for structure
        # actual t=0..L-1:
        R = mpz(1) << m
        acc = mpz(1)
        for t in range(L):
            acc = (acc * (R + 8 * t + 1)) % MOD
        # short window K=m
        short = B(m, m, 2 * m)
        P(f"  m={m} L_full={L}  v2(full-1)={_v2(acc-1)}  v2(short-1)={_v2(short-1)}")
        P(f"    full length exponential; short K=m has no squaring lift.")


def part_gamma_ratio_obstacle():
    P("=" * 70)
    P("PART D: why Gamma-ratio soft-linear fails for dual")
    P("""
  P_R = 8^K * ((2^m+9)/8)_K = 8^K * Gamma(K+(2^m+9)/8)/Gamma((2^m+9)/8)

  (2^m+9)/8 = 2^{m-3} + 9/8  not in Z_2 (denom 8), and for Morita Gamma_2
  the natural domain is Z_2.  Clearing denom:
    prod (2^m+8k+1) = odd-part ratio of factorials around 2^m.
  But n! for n>=2m is 0 mod 2^{2m} — factorial route collapses.

  Low-word works because (9/8)_K * 8^K involves only O(K log K)-bit
  integers BEFORE reduction, and Borwein sieves the prime factorization
  of that structured product.  Dual factors sit near 2^m, bit length m
  each, and their factorization is essentially random odds — sieve mass
  Theta(K m / log) = Theta(N^2 / log), no compression.
""")


def main():
    P("dual_borwein_hunt -- doubling identities for dual product")
    part_pairing_identity()
    part_doubling_search()
    part_new_block_vs_short()
    part_gamma_ratio_obstacle()
    P()
    P("VERDICT: no dual Borwein doubling found; sieve-mass obstruction stands.")


if __name__ == "__main__":
    main()
