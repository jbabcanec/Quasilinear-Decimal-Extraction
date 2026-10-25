#!/usr/bin/env python3
"""
attack_dual_product.py -- attack the pinned barrier:

    compute  P_R = prod_{k=1}^K (2^m + 8k + 1)  mod  2^{2m}
    in O(N log^2 N), then H = high_word(P_R * P_0^{-1}).

Routes tried here:
  (1) half-word / limb decomposition with carry tracking
  (2) recursive AP-halving: H_K = H_{K/2} + H^{offset}_{K/2}
  (3) multifactorial reduction to 8-step swing-style primes
  (4) check whether P_R mod 2^{2m} equals a simple function of cracked P_0
      and endpoint factors (LLL)

If any route yields bit-exact H at cost visibly below N log^3, escalate.
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
    R = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, R) for k in range(1, K + 1)) % R


def P0_mod(K, M):
    MOD = mpz(1) << M
    p = mpz(1)
    for k in range(1, K + 1):
        p = (p * (8 * k + 1)) % MOD
    return p


def PR_naive(K, m):
    R = mpz(1) << m
    MOD = R * R
    p = mpz(1)
    for k in range(1, K + 1):
        p = (p * (R + 8 * k + 1)) % MOD
    return p


def H_from_PR(K, m, PR):
    R = mpz(1) << m
    MOD = R * R
    P0 = P0_mod(K, 2 * m)
    u = (PR * pow(P0, -1, MOD)) % MOD
    assert (u - 1) % R == 0
    return ((u - 1) // R) % R


# ---------------------------------------------------------------------------
# Route 1: limb / half-word product (still Theta(K) muls -- baseline)
# ---------------------------------------------------------------------------

def PR_halfword(K, m):
    """Multiply as (low + R) factors with explicit low*low carry into high."""
    R = mpz(1) << m
    MODR = R
    # maintain acc = a0 + R * a1  (a0,a1 < R), representing value mod R^2
    a0, a1 = mpz(1), mpz(0)
    for k in range(1, K + 1):
        mu = 8 * k + 1  # < R for K << 2^{m-3}
        # (a0 + R a1) * (mu + R * 1) = a0*mu + R*(a0 + a1*mu)  + R^2 * ...
        lo = a0 * mu
        carry = lo >> m
        a0_new = lo & (R - 1)
        a1_new = (a0 + a1 * mu + carry) & (R - 1)
        a0, a1 = a0_new, a1_new
    return a0 + R * a1


def route1():
    P("=" * 70)
    P("ROUTE 1: half-word dual multiply (correctness + timing)")
    for K, m in [(32, 48), (64, 80)]:
        assert K < (1 << (m - 3))
        PR = PR_halfword(K, m)
        assert PR == PR_naive(K, m)
        assert H_from_PR(K, m, PR) == atom_H(K, m)
        P(f"  K={K} m={m}  MATCH")
    prev = None
    for K in [200, 400, 800, 1600]:
        m = 2 * K  # digit-scale; K << 2^{m-3} ok
        t0 = time.perf_counter()
        PR_halfword(K, m)
        dt = time.perf_counter() - t0
        ratio = "" if prev is None else f"  x{dt / prev:.2f}"
        P(f"  time K={K}: {dt:.4f}s{ratio}")
        prev = dt
    P("  still ~quadratic.  half-word does not remove the K factor.")


# ---------------------------------------------------------------------------
# Route 2: recursive AP halving -- cost model
# ---------------------------------------------------------------------------

def H_offset(K, m, b):
    """sum_{k=1}^K 1/(8k+b) mod 2^m, b odd."""
    R = mpz(1) << m
    return sum(pow(mpz(8 * k + b), -1, R) for k in range(1, K + 1)) % R


def route2():
    P()
    P("=" * 70)
    P("ROUTE 2: AP halving identity")
    # H_K = H_{K/2} + sum_{k=1}^{K/2} 1/(8k + 4K + 1)
    for K, m in [(32, 64), (64, 96), (128, 120)]:
        assert K % 2 == 0
        H = atom_H(K, m)
        Hl = atom_H(K // 2, m)
        b = 4 * K + 1
        Hr = H_offset(K // 2, m, b)
        ok = (Hl + Hr) % (1 << m) == H
        P(f"  K={K} m={m}  H=H_left+H_offset  MATCH={ok}")
        assert ok
    P("  Recurrence is exact.  Merge is another AP-harmonic of size K/2:")
    P("  T(K)=2T(K/2)+Cost(offset-H) => need Cost(offset-H)=O(M(m) log K)")
    P("  to win.  Offset-H is the SAME problem (door beta).  NO SAVE.")


# ---------------------------------------------------------------------------
# Route 3: 8-step multifactorial / prime-swing sketch
# ---------------------------------------------------------------------------

def route3():
    P()
    P("=" * 70)
    P("ROUTE 3: multifactorial reduction")
    P("  prod_k (R+8k+1) = (R+8K+1)!^{(8)} / (R+1)")
    P("  Swing/prime product needs primes <= R+8K+1 ~ 2^m = 2^{Theta(N)}.")
    P("  Sieve mass exponential.  KILL under sieve-mass axiom.")


# ---------------------------------------------------------------------------
# Route 4: LLL -- is PR a low-height combo of cracked objects?
# ---------------------------------------------------------------------------

def route4(K=36, m=96):
    P()
    P("=" * 70)
    P("ROUTE 4: is P_R in span of cracked P_0 / endpoints / G(admissible)?")
    R = mpz(1) << m
    MOD = R * R
    PR = PR_naive(K, m)
    P0 = P0_mod(K, 2 * m)
    # cracked-ish basket at precision 2m
    basket = [P0, (P0 * P0) % MOD, pow(P0, -1, MOD), mpz(1), R,
              mpz(8 * K + 1), pow(mpz(8 * K + 1), -1, MOD),
              (P0 * (1 + R)) % MOD,  # = P_R if H=1 nonsense
              ]
    # admissible G products
    for a in range(1, 6):
        t = 1 << a
        p = mpz(1)
        for k in range(1, K + 1):
            p = (p * (8 * k + 1 + t)) % MOD
        basket.append(p)
    tags = [f"b{i}" for i in range(len(basket))]
    target = int(PR)
    best = 0
    for i, b in enumerate(basket):
        for s in [1, -1]:
            d = (s * int(b) - target) % int(MOD)
            sc = 2 * m if d == 0 else v2(d)
            best = max(best, sc)
    P(f"  best height-1 agree vs P_R: {best}/{2*m}")
    # specifically P0*(1+R*Htruth) must match -- uses H, not cracked
    H = atom_H(K, m)
    expect = (P0 * (1 + R * H)) % MOD
    P(f"  P0*(1+R H) == P_R? {expect == PR}  (sanity)")
    P(f"  no cracked stand-in for the (1+R H) factor without H.")


def main():
    P("attack_dual_product -- barrier = P_R mod 2^{2m}")
    route1()
    route2()
    route3()
    route4()
    P()
    P("VERDICT: dual product still requires either Theta(K) M-bit muls,")
    P("exponential sieve, or an H-oracle.  Barrier holds.")
    P()
    P("PULLS LANDED (see odd_prod_lift.py, faulhaber_H.py, egf_power_sums.py):")
    P("  * OEIS A395199: (2^n-1)!! lifts by squaring — FULL odds, poly(n).")
    P("  * new_block lemma: B_{m+1}/B_m^2 = 1+O(2^{2m-3}) — LONG 1-mod-8")
    P("    dyadic block (length 2^{m-3}), wrong regime for digit K=Theta(m).")
    P("  * Faulhaber-geometric H bit-exact but O(N^2).")
    P("  * EGF batch power sums: 1/n! not in Z_2 — blocked at p=2.")
    P("Next: OGF/Newton identities for AP power sums; dual Borwein;")
    P("      short-window structure at height 2^m.")


if __name__ == "__main__":
    main()
