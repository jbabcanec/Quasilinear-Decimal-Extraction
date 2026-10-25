#!/usr/bin/env python3
"""
diff_approach.py -- DIFFERENT route to two logs.

Stop merging inverses. Write the atom H via Euler-Maclaurin as
  H = (easy Log + endpoint) + (Bernoulli correction),
and try to replace the Bernoulli sum by a 2-adic special function
built from exp/log (Caruso quasi-linear).

If corr = F(phi(u), u, Log(...)) for phi = u/(e^u-1), we win:
  u = 8/(8K+1) has v2(u)=3, inside the p=2 exp disk, and
  exp/log at precision m cost O(M(m) log m) = O(N log^2 N).

UPDATE (see also diff_quotient_H.py, attack_dual_product.py):
  Exact identities found and validated:
    H == (Log P(2^r) - Log P(0)) / 2^r     (r large)
    P(2^m) == P(0) * (1 + 2^m H)  mod 2^{2m}
    H == sum_j (-8)^j * Faulhaber_j(K)       (2-adic)
  Complexity fence: large-r / dual product is Theta(K M(m));
  admissible stencils rate-limited to O(log N) bits (door beta).
  Barrier pinned to: prod_k (2^m + 8k + 1) mod 2^{2m} in O(N log^2 N).
  Dual squaring lemma (validation/dual_lift.py): with R=2^m, mu_k=8k+1,
    B(m,K)^2 - P0(K)*B(m+1,K) has v2 >= 2m,  exact ratio
    prod (1 + R^2/(mu(mu+2R))).  Corrected m->m+1 lift needs sum mu^{-2}
    (additive tower); +1 steps => O(N^2 log N) by rate law (door beta).
  Exp/log basket hunt below: no low-height hit for EM correction.
"""
from __future__ import annotations
import itertools

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


def odd_inv(a, m):
    a = mpz(a)
    assert int(a) % 2 == 1
    return pow(a, -1, 1 << m)


def exp2adic(x, m):
    """exp(x) mod 2^m; requires v2(x) >= 2 for guaranteed convergence at p=2."""
    MOD = mpz(1) << m
    x = mpz(x) % MOD
    assert v2(x) >= 2, "exp disk"
    s, t = mpz(1), mpz(1)
    for n in range(1, m + 16):
        t = (t * x) % MOD
        e = v2(n)
        if v2(int(t)) < e:
            break
        if e:
            t = ((t >> e) * odd_inv(n >> e, m)) % MOD
        else:
            t = (t * odd_inv(n, m)) % MOD
        s = (s + t) % MOD
        if t == 0:
            break
    return s


def log2adic_unit(x, m):
    """Log2(x) for odd x."""
    MOD = mpz(1) << m
    y = (mpz(x) - 1) % MOD
    s, t = mpz(0), mpz(1)
    for n in range(1, m + 16):
        t = (t * y) % MOD
        e = v2(n)
        if e:
            if v2(int(t)) < e:
                break
            term = ((t >> e) * odd_inv(n >> e, m)) % MOD
        else:
            term = (t * odd_inv(n, m)) % MOD
        s = (s + term) % MOD if (n % 2) else (s - term) % MOD
    return s


def div_by_2pow(x, e, m):
    """x / 2^e mod 2^m, requiring v2(x) >= e."""
    x = mpz(x)
    assert v2(x) >= e
    return (x >> e) % (1 << m)


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(odd_inv(8 * k + 1, m) for k in range(1, K + 1)) % MOD


def main_term(K, m):
    """F(K)-F(0)+(f(K)-f(0))/2 with F=Log(8x+1)/8, f=1/(8x+1)."""
    MOD = mpz(1) << m
    Log_end = log2adic_unit(8 * K + 1, m)  # Log(1+8K)
    # Log(1+8K)/8 : Log has v2 >= 3
    FK = div_by_2pow(Log_end, 3, m)
    fK = odd_inv(8 * K + 1, m)
    f0 = mpz(1)
    return (FK + div_by_2pow((fK - f0) % MOD, 1, m)) % MOD


def phi_u(u, m):
    """u/(e^u - 1) mod 2^m (both even; cancel 2-powers)."""
    MOD = mpz(1) << m
    eu = exp2adic(u, m)
    emu1 = (eu - 1) % MOD
    vu, ve = v2(u), v2(emu1)
    e = min(vu, ve)
    return (((u >> e) * odd_inv(emu1 >> e, m)) << (vu - e)) % MOD


def hunt(K=32, m=192):
    P("=" * 70)
    P("DIFF APPROACH: special-function EM, not inverse-merge")
    P(f"K={K} m={m}")
    MOD = mpz(1) << m
    H = atom_H(K, m)
    main = main_term(K, m)
    corr = (H - main) % MOD
    u = (mpz(8) * odd_inv(8 * K + 1, m)) % MOD
    P(f"v2(u)={v2(u)} (need >=2 for exp); corr v2={v2(corr)}")

    phi = phi_u(u, m)
    eu = exp2adic(u, m)
    # basket of exp/log-buildable quantities
    basket_names = ["phi", "u", "u2", "u3", "eu", "1", "fK", "FK", "Log_end"]
    fK = odd_inv(8 * K + 1, m)
    FK = div_by_2pow(log2adic_unit(8 * K + 1, m), 3, m)
    Log_end = log2adic_unit(8 * K + 1, m)
    basket = [
        phi, u, (u * u) % MOD, (u * u * u) % MOD, eu,
        mpz(1), fK, FK, Log_end,
    ]

    # also phi-1, phi*u, Log via eu
    basket_names += ["phi-1", "phi*u", "u/(eu)"]
    basket += [
        (phi - 1) % MOD,
        (phi * u) % MOD,
        (u * odd_inv(eu, m)) % MOD,
    ]

    P(f"basket size {len(basket)}: {basket_names}")
    P("brute height-2 integer combos against corr...")

    best = None
    # height 1 first
    for i, b in enumerate(basket):
        for s in (+1, -1):
            d = (s * b - corr) % MOD
            sc = m if d == 0 else v2(d)
            if sc >= m // 2:
                P(f"  h1 HIT {s}*{basket_names[i]} agree={sc}/{m}")
                best = (1, sc, (s, i))
    # height 2
    if best is None or best[1] < m - 8:
        for i, j in itertools.combinations(range(len(basket)), 2):
            for si, sj in itertools.product((-2, -1, 1, 2), repeat=2):
                acc = (si * basket[i] + sj * basket[j]) % MOD
                d = (acc - corr) % MOD
                sc = m if d == 0 else v2(d)
                if sc >= (best[1] if best else m // 2):
                    if sc >= m - 8:
                        P(f"  h2 HIT {si}*{basket_names[i]} + {sj}*{basket_names[j]} agree={sc}/{m}")
                        best = (2, sc, (si, i, sj, j))
                        if sc >= m - 2:
                            return best
    if best is None:
        P("  no low-height hit. Showing agreement of individuals vs corr:")
        for n, b in zip(basket_names, basket):
            d = (b - corr) % MOD
            P(f"    {n:12s} v2(diff)={0 if d==0 else v2(d)}")
        # Does main+c*phi recover H for some small c?
        P("  scan c such that main + c*phi/8 ~= H:")
        for c in range(-16, 17):
            # phi/8
            if v2(phi) < 3:
                continue
            trial = (main + c * div_by_2pow(phi, 3, m)) % MOD
            d = (trial - H) % MOD
            sc = m if d == 0 else v2(d)
            if sc >= 16:
                P(f"    c={c:3d} agree={sc}/{m}")
    else:
        P(f"BEST {best}")
    return best


def em_bernoulli_reference(K=24, m=128, IMAX=40):
    """Confirm EM with Bernoulli still matches H (sanity)."""
    from fractions import Fraction
    from math import factorial, comb

    def bernoulli(n):
        B = [Fraction(0)] * (n + 1)
        B[0] = Fraction(1)
        for md in range(1, n + 1):
            B[md] = -sum(comb(md + 1, j) * B[j] for j in range(md)) / (md + 1)
        return B

    B = bernoulli(2 * IMAX + 2)

    def to_mod(fr):
        return (fr.numerator * odd_inv(fr.denominator, m)) % (1 << m)

    # use piece formula we verified: piece_i = -B_{2i}/(2i)*8^{2i-1}*((8K+1)^{-2i}-1)
    corr = Fraction(0)
    for i in range(1, IMAX + 1):
        piece = -B[2 * i] / Fraction(2 * i) * Fraction(8) ** (2 * i - 1) * (
            Fraction(1, (8 * K + 1) ** (2 * i)) - 1
        )
        corr += piece

    # main in Q via truncated? use 2-adic main
    main = main_term(K, m)
    H = atom_H(K, m)
    corr_m = to_mod(corr)
    got = (main + corr_m) % (1 << m)
    d = (got - H) % (1 << m)
    sc = m if d == 0 else v2(d)
    P("=" * 70)
    P(f"SANITY EM+Bernoulli: K={K} m={m} IMAX={IMAX} agree={sc}/{m}")
    return sc


if __name__ == "__main__":
    em_bernoulli_reference()
    P()
    hunt(32, 192)
    P()
    hunt(64, 256)
