#!/usr/bin/env python3
"""
bet3_leave_one_out.py -- Live bet 3:
  compress leave-one-out e_{T-1} = sum_i prod_{j!=i} mu_j
  the way Borwein compresses e_T = prod mu_j.

Borwein works because e_T is SMOOTH (prime factors <= 8K+1).
e_{T-1} is a generic integer ~ exp(Theta(K log K)) bits with large prime factors.

Attacks:
  A. Measure smoothness / largest prime factor of e_{T-1} (small K)
  B. Factor e_{T-1} in special form: e_{T-1} = P0 * H_rational = P0 * (num/den)
     but H_rational = e_{T-1}/P0 already
  C. Write e_{T-1} as resultant / determinant / Hankel and look for structure
  D. Segmented sieve for leave-one-out? (doesn't make sense — not a product)
  E. Pairing / duplication: e_{T-1}(2K) in terms of e_{T-1}(K) and cracked stuff
  F. Represent e_{T-1} in factored form via primes p <= 8K with
     exponents = sum_i 1_{p|prod_{j!=i} mu_j} = (#terms not divisible by p)*v_p
     WAIT: e_{T-1} = sum of products, NOT a product.  Cannot take v_p of a sum
     as sum of v_p.  Unless one term dominates p-adically — rare.
  G. Common denominator representation and cancel against P0 for H mod 2^m
"""
from __future__ import annotations
import math

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def mu(k):
    return 8 * k + 1


def P0(K):
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= mu(k)
    return acc


def e_Tm1(K):
    """Leave-one-out numerator."""
    s = mpz(0)
    for i in range(1, K + 1):
        term = mpz(1)
        for j in range(1, K + 1):
            if j != i:
                term *= mu(j)
        s += term
    return s


def factor_trial(n, bound=None):
    """Trial factorization; return list of (p,e)."""
    n = int(n)
    if n < 0:
        n = -n
    if bound is None:
        bound = int(math.isqrt(n)) + 1
    factors = []
    # strip 2
    e = 0
    while n % 2 == 0:
        n //= 2
        e += 1
    if e:
        factors.append((2, e))
    p = 3
    while p * p <= n and p <= bound:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            factors.append((p, e))
        p += 2
    if n > 1:
        factors.append((n, 1))  # cofactor prime or composite
    return factors


def part_A_smoothness():
    P("=" * 70)
    P("BET3-A: smoothness of e_{T-1} vs e_T=P0")
    for K in [10, 15, 20, 25, 30]:
        p0 = P0(K)
        e = e_Tm1(K)
        fp0 = factor_trial(p0, bound=8 * K + 1)
        # for e_{T-1}, trial up to 8K first then report cofactor size
        fe_small = []
        rem = int(e)
        for p in range(2, 8 * K + 2):
            if rem % p == 0:
                ee = 0
                while rem % p == 0:
                    rem //= p
                    ee += 1
                fe_small.append((p, ee))
        P(f"  K={K}")
        P(f"    P0 bits={p0.bit_length()}  max_prime_factor<={8*K+1}  "
          f"fully_smooth_over_8K={rem==1 if False else all(p<=8*K+1 for p,e in factor_trial(p0))}")
        P(f"    eT-1 bits={e.bit_length()}  after stripping p<=8K+1: "
          f"cofactor_bits={rem.bit_length()}  cofactor={rem if rem.bit_length()<80 else 'big'}")
        if rem > 1 and rem.bit_length() <= 60:
            # check primality-ish
            P(f"      cofactor factors: {factor_trial(rem)}")


def part_B_duplication():
    P("=" * 70)
    P("BET3-B: e_{T-1}(2K) vs e_{T-1}(K) and P0")
    for K in [4, 6, 8]:
        eK = e_Tm1(K)
        e2 = e_Tm1(2 * K)
        pK = P0(K)
        p2 = P0(2 * K)
        # search small relations e2 = a*eK*pK + b*pK^2 + ...
        # e2 / (eK * something)
        P(f"  K={K}  bits eK={eK.bit_length()} e2={e2.bit_length()} pK={pK.bit_length()}")
        # structural: terms for i<=K and i>K
        # e2 = sum_{i=1}^{2K} prod_{j!=i} mu_j
        # = sum_{i<=K} (prod_{j<=K,j!=i} mu_j) * (prod_{j=K+1}^{2K} mu_j)
        # + sum_{i>K} (prod_{j<=K} mu_j) * (prod_{j>K,j!=i} mu_j)
        # = (prod_{K+1..2K} mu) * e_{T-1}(K) + P0(K) * e_{T-1}^{shifted}(K)
        shift_P0 = mpz(1)
        for j in range(K + 1, 2 * K + 1):
            shift_P0 *= mu(j)
        # shifted leave-one-out for second half
        e_shift = mpz(0)
        for i in range(K + 1, 2 * K + 1):
            term = mpz(1)
            for j in range(K + 1, 2 * K + 1):
                if j != i:
                    term *= mu(j)
            e_shift += term
        recon = shift_P0 * eK + pK * e_shift
        P(f"    reconstruction identity: {recon==e2}")
        P(f"    still needs e_shift = leave-one-out on SHIFTED AP (same hardness)")


def part_C_determinant():
    P("=" * 70)
    P("BET3-C: e_{T-1} as coefficient / resultant")
    P("  e_{T-1} = [x^{K-1}] prod (x - (-mu)) * (-1)^{...} = elem sym e_{K-1}")
    P("  Computing charpoly of diag(mu) is product tree = Theta(K M) for mod 2^m.")
    P("  No known soft-linear for a SINGLE elementary symmetric e_{K-1}")
    P("  without the others when the roots are an AP — unless via H*P0.")


def part_D_vp_of_sum():
    P("=" * 70)
    P("BET3-D: can v_p(e_{T-1}) be read from AP counting like Borwein?")
    for K in [12, 18]:
        e = e_Tm1(K)
        # pick a medium prime
        for p in [3, 5, 7, 11, 13]:
            # v_p of sum is NOT sum of v_p
            ve = 0
            ee = int(e)
            while ee % p == 0:
                ee //= p
                ve += 1
            # compare to min v_p of terms
            min_term = 99
            for i in range(1, K + 1):
                vt = 0
                for j in range(1, K + 1):
                    if j == i:
                        continue
                    m = mu(j)
                    while m % p == 0:
                        m //= p
                        vt += 1
                min_term = min(min_term, vt)
            P(f"  K={K} p={p}  v_p(e)={ve}  min_v_p(terms)={min_term}  "
              f"gap={ve-min_term}")
    P("  Often v_p(sum) > min v_p(terms) (cancellation).  Counting alone insufficient.")


def part_E_mod_2m_without_full():
    P("=" * 70)
    P("BET3-E: need e_{T-1} mod 2^m * (via H = e_{T-1} * P0^{-1})")
    P("  Equivalently e_{T-1} mod 2^m after fixing P0 inv.")
    P("  Borwein gives P0 mod 2^m.  Soft-linear e_{T-1} mod 2^m IS Open Problem E.")
    for K, m in [(20, 40), (40, 48)]:
        e = e_Tm1(K)
        p0 = P0(K)
        MOD = mpz(1) << m
        H = (e * pow(p0 % MOD, -1, MOD)) % MOD
        Href = sum(pow(mpz(mu(k)), -1, MOD) for k in range(1, K + 1)) % MOD
        P(f"  K={K} m={m}  H_from_e={H==Href}")


def part_F_sieve_exponent_fantasy():
    P("=" * 70)
    P("BET3-F: fantasy Borwein for sums — KILL")
    P("""
  For products, v_p(prod mu_k) = sum_k v_p(mu_k) = AP count, soft-linear assemble.

  For sums S = sum_i a_i with a_i = prod_{j!=i} mu_j:
    v_p(S) is NOT determined by the multiset {v_p(a_i)} alone (cancellation).
    The p-adic digits of S require the actual residues a_i mod p^{e+c},
    i.e. essentially the leave-one-out terms mod p-powers —
    which is the same problem at every odd prime, then CRT
    (Omega(N^2) for exact, or back to 2-adic for mod 2^m).

  Compression that works for products (exponent vectors) fails for sums.
""")


def main():
    P("bet3_leave_one_out -- Borwein-style compression for e_{T-1}")
    part_A_smoothness()
    part_B_duplication()
    part_C_determinant()
    part_D_vp_of_sum()
    part_E_mod_2m_without_full()
    part_F_sieve_exponent_fantasy()
    P()
    P("BET3 VERDICT: e_{T-1} is NOT smooth; cofactor has large primes.")
    P("  Duplication reduces to shifted leave-one-out (same hardness).")
    P("  v_p(sum) has cancellation; Borwein exponent calculus does not transfer.")
    P("  Leave-one-out compression remains open — and looks structurally blocked.")


if __name__ == "__main__":
    main()
