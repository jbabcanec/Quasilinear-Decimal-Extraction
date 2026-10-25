#!/usr/bin/env python3
"""
dual_jump.py -- hunt geometric m-jumps for the dual product.

Have: B(m)^2 ≡ P0 * B(m+1)  (mod 2^{2m+})
Want: B(m)^{2^s} ≡ ? * B(m+s) with v2(err) >= 2(m+s)  (full target width)

Or: relation B(2m, K) vs B(m, K) with usable valuation.
"""
from __future__ import annotations
from itertools import product as iproduct

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def v2(n):
    n = int(n)
    if n == 0:
        return 9999
    if n < 0:
        n = -n
    return (n & -n).bit_length() - 1


def B(m, K):
    R = mpz(1) << m
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= R + 8 * k + 1
    return acc


def P0(K):
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= 8 * k + 1
    return acc


def part_chain_power():
    P("=" * 70)
    P("PART A: v2( B(m)^{2^s} - P0^{2^s - 1} B(m+s) ) vs need 2(m+s)")
    for m in [8, 10, 12]:
        for K in [m, 8]:
            for s in [1, 2, 3, 4]:
                # B(m)^{2^s} vs P0^{2^s - 1} * B(m+s)
                # from iterating B(j)^2 = P0 B(j+1) approximately
                lhs = B(m, K) ** (1 << s)
                rhs = (P0(K) ** ((1 << s) - 1)) * B(m + s, K)
                ve = v2(lhs - rhs)
                need = 2 * (m + s)
                P(f"  m={m} K={K} s={s}  v2={ve:4d}  need={need:4d}  "
                  f"margin={ve-need:4d}  full={ve>=need}")


def part_double_m():
    P("=" * 70)
    P("PART B: B(2m) vs polynomials in B(m), P0")
    # Brute small exponents: B(2m)^a * B(m)^b * P0^c == 1 mod 2^w?
    for m in [6, 8, 10]:
        K = m
        w = 2 * m  # compare at low width first
        MOD = mpz(1) << w
        vals = {
            "B2m": B(2 * m, K) % MOD,
            "Bm": B(m, K) % MOD,
            "P0": P0(K) % MOD,
            "B2m_halfK": B(2 * m, K // 2) % MOD if K >= 2 else mpz(1),
            "Bm1": B(m + 1, K) % MOD,
        }
        hits = []
        names = list(vals.keys())
        for exps in iproduct(range(-2, 3), repeat=len(names)):
            if all(e == 0 for e in exps):
                continue
            if sum(abs(e) for e in exps) > 5:
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
            if ok and acc == 1:
                rel = "*".join(f"{n}^{e}" for n, e in zip(names, exps) if e)
                hits.append(rel)
        P(f"  m={m} K={K} w={w}: {len(hits)} hits")
        for h in hits[:8]:
            P(f"    {h}")


def part_explicit_double():
    """
    Algebra: B(2m)/B(m) = prod (2^{2m}+mu)/(2^m+mu)
      = prod (2^m - mu + mu^2/(2^m+mu))
    Not a unit-power of P0.  Check v2(B(2m) - B(m)^2 / something).
    """
    P("=" * 70)
    P("PART C: direct B(2m) vs B(m)^2 / P0^{?} ")
    for m in [6, 8, 10, 12]:
        for K in [8, m]:
            Bm = B(m, K)
            B2 = B(2 * m, K)
            p0 = P0(K)
            for e in range(0, 5):
                # compare B2 * p0^e vs Bm^2 ?
                # or B2 vs Bm^2 / p0^e
                try:
                    # integer: v2(B2 * p0^e - Bm^2)
                    ve = v2(B2 * (p0 ** e) - Bm * Bm)
                except Exception:
                    ve = -1
                P(f"  m={m} K={K} e={e}  v2(B2*P0^e - Bm^2)={ve}  "
                  f"2m={2*m} 4m={4*m}")


def part_inv_power_batch():
    """
    Inverse power sums Q_j = sum mu^{-j}.
    Q_1 = H is the barrier.  Are Q_2, Q_4, ... cheaper?
    Newton from Log: if we had P(t)=prod(mu+t) as series, Log P gives all.
    Dyadic: Q_{2^r} via Freshman?  (a^{-1})^{2^r} = a^{-2^r}, sum still additive.
    """
    P("=" * 70)
    P("PART D: inverse-power sums -- all in the additive cluster")
    for K, m in [(20, 48), (40, 64)]:
        MOD = mpz(1) << m
        Q = []
        for j in range(1, 6):
            s = sum(pow(mpz(8 * k + 1), -j, MOD) for k in range(1, K + 1)) % MOD
            Q.append(s)
        # check Newton-ish: no free Q_2 from P0 alone
        P(f"  K={K} m={m}  Q1..Q5 computed directly (each Theta(K) inversions)")
        P(f"    Q1={hex(int(Q[0]))[:18]}...")
    P("  No soft-linear batch for Q_j known; Q_1 is Op E, Q_j>=2 same tower.")


def main():
    P("dual_jump -- geometric m-jump hunt")
    part_chain_power()
    part_double_m()
    part_explicit_double()
    part_inv_power_batch()
    P()
    P("VERDICT: no full-width geometric jump; additive tower intact.")


if __name__ == "__main__":
    main()
