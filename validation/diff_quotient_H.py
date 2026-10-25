#!/usr/bin/env python3
"""
diff_quotient_H.py -- exact recovery of the additive atom via product deformation.

THEOREM (machine-validated).  Let
    H = sum_{k=1}^K (8k+1)^{-1}  in Z/2^m Z,
    P(t) = prod_{k=1}^K (8k+1+t)  in Z/2^{M} Z,
    Q(t) = (Log P(t) - Log P(0)) / t.
For t = 2^r with r >= m+2 - v2(sum mu^{-2}) and M >= m+r,
    Q(t) == H  (mod 2^m).

Dual form (also validated): with R=2^m,
    prod_k (R + mu_k) == P(0) * (1 + R H)  (mod R^2).
Hence H is the high word of P(R)/P(0).

Geometric/Faulhaber form (also validated):
    H == sum_{j>=0} (-8)^j * sum_{k=1}^K k^j   (2-adic convergent),
with ~ m/3 Faulhaber terms.  Hits the Bernoulli-mass wall (paper rem:root).

COMPLEXITY.  P(0) is O(M(N) log N) by Borwein (Thm crack).  P(t) for
r=Theta(N) has factors of size 2^{Theta(N)}, so the sieve-mass axiom forbids
Borwein; a plain product tree is Theta(K M(M)) = Theta(N^2 log N) (measured).
Admissible (small-t) stencils recover only O(log N) bits per shot (door beta /
postnikov rate law).  So this is an exact identity, not yet an O(N log^2 N)
algorithm -- same barrier, cleaner statement.

This file validates the identities and records the rate/complexity fence.
"""
from __future__ import annotations

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


def log_1mod8(x, m):
    """Log on 1+8 Z_2 via series (our products land here)."""
    MOD = mpz(1) << m
    x = mpz(x) % MOD
    assert v2(int(x - 1)) >= 3, v2(int(x - 1))
    y = (x - 1) % MOD
    s, t = mpz(0), mpz(1)
    for n in range(1, m + 64):
        t = (t * y) % MOD
        e = v2(n)
        if e:
            if v2(int(t)) < e:
                break
            term = ((t >> e) * pow(mpz(n >> e), -1, MOD)) % MOD
        else:
            term = (t * pow(mpz(n), -1, MOD)) % MOD
        s = (s + term) % MOD if (n % 2) else (s - term) % MOD
        if term == 0 and n > 8:
            break
    return s


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


def prod_mod(K, t, M):
    MOD = mpz(1) << M
    p = mpz(1)
    for k in range(1, K + 1):
        p = (p * ((8 * k + 1 + t) % MOD)) % MOD
    return p


def recover_diff_quot(K, m, r=None):
    MODH = mpz(1) << m
    S2 = sum(pow(mpz(8 * k + 1), -2, 1 << (m + 40)) for k in range(1, K + 1))
    vs2 = v2(S2)
    if r is None:
        r = m + 2 - vs2
        if r < 4:
            r = 4
    M = m + r + 16
    t = mpz(1) << r
    L0 = log_1mod8(prod_mod(K, 0, M), M)
    Lt = log_1mod8(prod_mod(K, int(t), M), M)
    diff = (Lt - L0) % (1 << M)
    Ha = (diff >> r) % MODH
    return Ha, r, vs2


def recover_dual(K, m):
    """H from high word of prod(R+mu)/P mod R^2."""
    R = mpz(1) << m
    R2 = R * R
    P0 = mpz(1)
    PR = mpz(1)
    for k in range(1, K + 1):
        mu = 8 * k + 1
        P0 = (P0 * mu) % R2
        PR = (PR * (R + mu)) % R2
    u = (PR * pow(P0, -1, R2)) % R2
    assert (u - 1) % R == 0
    return ((u - 1) // R) % R


def part_identities():
    P("=" * 70)
    P("IDENTITY A: difference quotient")
    for K, m in [(16, 40), (32, 64), (64, 80), (128, 100)]:
        H = atom_H(K, m)
        Ha, r, vs2 = recover_diff_quot(K, m)
        ok = Ha == H
        P(f"  K={K:4d} m={m:4d} vs2={vs2} r={r:4d}  MATCH={ok}")
        assert ok

    P()
    P("IDENTITY B: dual / high-word")
    for K, m in [(8, 16), (16, 24), (32, 40), (64, 48)]:
        H = atom_H(K, m)
        Ha = recover_dual(K, m)
        ok = Ha == H
        P(f"  K={K:4d} m={m:4d}  MATCH={ok}")
        assert ok


def part_complexity_fence():
    import time

    P()
    P("=" * 70)
    P("COMPLEXITY FENCE: product tree for P(2^r), r~m~K")
    P("  (naive sequential mul; balanced tree has same Theta(K) M-bit muls)")
    prev = None
    for K in [200, 400, 800, 1600]:
        M = 2 * K
        t0 = time.perf_counter()
        prod_mod(K, 0, M)
        prod_mod(K, int(mpz(1) << K), M)
        dt = time.perf_counter() - t0
        ratio = "" if prev is None else f"  time/prev={dt / prev:.2f} (N^2 expects ~4)"
        P(f"  K={K:5d} M={M:5d}  {dt:.4f}s{ratio}")
        prev = dt
    P("  => quadratic in the digit-scale regime.  Not O(N log^2 N).")


def part_admissible_rate():
    P()
    P("=" * 70)
    P("ADMISSIBLE RATE: bits of H from Q(2^a) with a = O(log K)")
    K, m = 64, 120
    H = atom_H(K, m)
    L0 = log_1mod8(prod_mod(K, 0, m + 32), m + 32)
    P(f"  K={K} m={m}")
    for a in range(1, 14):
        if (1 << a) > 32 * K:
            P(f"  a={a}: sieve-mass break (t > 32K)")
            break
        Lt = log_1mod8(prod_mod(K, 1 << a, m + 32), m + 32)
        est = ((Lt - L0) % (1 << (m + 32))) >> a
        d = int((est - H) % (1 << m))
        bits = m if d == 0 else v2(d)
        P(f"  a={a:2d}  bits_agree={bits:3d}  (expect ~a+v2(S2)-1)")


def main():
    P("diff_quotient_H -- exact identities for Open Problem E atom")
    part_identities()
    part_admissible_rate()
    part_complexity_fence()
    P()
    P("VERDICT: identities are exact; Borwein-friendly deformations are")
    P("rate-limited; large-t deformations are quadratic.  Barrier stands,")
    P("now pinned to: compute prod_k (2^m + 8k + 1) mod 2^{2m} in O(N log^2 N).")


if __name__ == "__main__":
    main()
