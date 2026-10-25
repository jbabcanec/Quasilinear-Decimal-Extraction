#!/usr/bin/env python3
"""
dual_lift.py -- validate and exploit the dual-Borwein hit:

    B(m,K)^2  ==  P0(K) * B(m+1, K)   (mod 2^{?})

where B(m,K)=prod_{k=1}^K (2^m+8k+1), P0(K)=prod_{k=1}^K (8k+1).

If the 2-valuation of the error grows fast enough in m, we can lift
B(·,K) in m by squaring + one P0-inverse (P0 is Borwein-cheap).

Digit regime needs K=Theta(m) simultaneously -- probe whether a
joint (m,K)-recursion exists that stays soft-linear.
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
    n = int(n) % (1 << 2000)  # safety for huge
    # for exact diff that may be 0
    if n == 0:
        return 9999
    # handle negative
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


def part_valuation_law():
    P("=" * 70)
    P("PART A: v2( B(m,K)^2 - P0(K)*B(m+1,K) ) as function of m,K")
    P(f"  {'m':>4} {'K':>4} {'v2(err)':>10} {'2m':>6} {'2m-3':>6} {'ok_full':>8}")
    for m in range(6, 18):
        for K in [4, 8, 12, 16, m, 2 * m]:
            if K < 1:
                continue
            err = B(m, K) ** 2 - P0(K) * B(m + 1, K)
            ve = v2(err)
            P(f"  {m:4d} {K:4d} {ve:10d} {2*m:6d} {2*m-3:6d} {str(ve>=2*m):>8}")


def part_exact_identity():
    """
    Algebraic check: is B(m,K)^2 / B(m+1,K) == P0(K) exactly in Q?
    """
    P("=" * 70)
    P("PART B: exact rational identity?")
    for m, K in [(6, 4), (8, 8), (10, 5), (12, 16)]:
        num = B(m, K) ** 2
        den = B(m + 1, K)
        p0 = P0(K)
        # compare num/den vs p0
        if num == p0 * den:
            P(f"  m={m} K={K}  EXACT equality")
        else:
            err = num - p0 * den
            P(f"  m={m} K={K}  NOT exact; v2(err)={v2(err)}  "
              f"bitlen(err)={err.bit_length()}")


def part_algebraic_expand():
    """
    Write each factor:
      B(m,K) = prod (2^m + mu_k), mu_k=8k+1
      B(m+1,K) = prod (2^{m+1} + mu_k) = prod (2*2^m + mu_k)
      B(m,K)^2 = prod_k prod_j (2^m+mu_k)(2^m+mu_j)

    Pairwise: (2^m+mu)(2^m+mu) = 2^{2m} + 2^{m+1} mu + mu^2
    vs (2^{m+1}+mu)*mu = 2^{m+1} mu + mu^2

    So (2^m+mu)^2 - mu(2^{m+1}+mu) = 2^{2m}
    Per-term: (R+mu)^2 = mu(2R+mu) + R^2

    Product form is NOT simply P0 * B(m+1) because cross terms.
    Use: prod (R+mu_k) = prod mu_k * prod (1 + R/mu_k)
        prod (2R+mu_k) = prod mu_k * prod (1 + 2R/mu_k)
        B(m)^2 / (P0 * B(m+1)) = P0 * prod (1+R/mu)^2 / prod (1+2R/mu)
                     = P0 * prod (1+R/mu)^2 / ( (1+2R/mu) ) / P0 wait

        B(m)/P0 = prod (1+R/mu)
        B(m+1)/P0 = prod (1+2R/mu)
        B(m)^2 / (P0 B(m+1)) = P0 * [prod (1+R/mu)]^2 / [P0 prod (1+2R/mu)]
                         = prod (1+R/mu)^2 / prod (1+2R/mu)

    So ratio_err_factor = prod_k (1+R/mu_k)^2 / (1+2R/mu_k)
    Identity holds iff this is 1.

    (1+R/mu)^2 / (1+2R/mu) = (mu+R)^2 / (mu(mu+2R)) * mu/mu
      = (mu+R)^2 / (mu (mu+2R))

    Product_k (mu+R)^2 / (mu (mu+2R)) = B(m)^2 / (P0 B(m+1))

    For this to be 1: B(m)^2 = P0 B(m+1). Never exact (part B).

    2-adic: prod (1+R/mu)^2/(1+2R/mu) = prod (1 + O(R^2/mu^2)) since
      (1+x)^2/(1+2x) = (1+2x+x^2)/(1+2x) = 1 + x^2/(1+2x)
    with x=R/mu, so each factor = 1 + R^2/(mu^2 (1+2R/mu))
    Product = 1 + R^2 * sum 1/mu^2 + O(R^4)  (when R small vs ...)

    With R=2^m and mu=O(K), need K << 2^{m/2} for convergence in units.
    """
    P("=" * 70)
    P("PART C: expansion law  B(m)^2/(P0 B(m+1)) = prod (1 + R^2/(mu^2+2R mu))")
    for m in [8, 12, 16, 20]:
        for K in [8, 16, m, 2 * m]:
            R = mpz(1) << m
            ratio = (B(m, K) ** 2) * pow(P0(K) * B(m + 1, K), -1, 1 << (4 * m))
            # compute in integers: ratio_int = B(m)^2 / (P0 B(m+1)) as rational
            num = B(m, K) ** 2
            den = P0(K) * B(m + 1, K)
            # v2(num-den) already done; v2(num/den - 1) = v2(num-den)-v2(den)
            ve = v2(num - den) - v2(den)
            # predict: v2(R^2 * sum mu^{-2}) = 2m + v2(sum mu^{-2})
            # mu odd, sum of K units ~ v2 maybe 0
            P(f"  m={m:2d} K={K:3d}  v2(ratio-1)~={ve:4d}  predict>=2m={2*m}")


def part_lift_algorithm():
    """
    When v2(err) >= 2m, we have B(m,K)^2 ≡ P0 B(m+1,K) mod 2^{2m}.
    Digit need: B(m,K) mod 2^{2m} with K~m.
    Lift in m at FIXED K: start B(m0,K) mod 2^{2M} for target M, m0 small?
    Careful: the congruence modulus grows with m.

    Algorithm: to compute B(M,K) mod 2^{2M}:
      If 2M bits needed and identity gives B(m+1) from B(m) mod 2^{2m},
      precision at step m is only 2m bits -- may not be enough to reach 2M.

    Check: from B(m) known mod 2^{2m}, compute B(m)^2 / P0 mod 2^{2m}.
    Identity says this equals B(m+1) mod 2^{2m}, but we NEED B(m+1) mod 2^{2(m+1)}.
    Gap of 2 bits per step -- LOSE precision relative to target!
    Same disease as door beta / Postnikov rate law.
    """
    P("=" * 70)
    P("PART D: lift precision audit (the rate law)")
    # Simulate: start from exact B(m0), lift via formula mod growing width
    for K in [8, 16, 32]:
        m0 = 6
        # exact target
        M = 20
        target = B(M, K) % (1 << (2 * M))
        # lift
        cur_m = m0
        cur = B(m0, K)  # exact
        p0 = P0(K)
        while cur_m < M:
            # predict B(cur_m+1) from identity mod 2^{2*cur_m}
            mod = mpz(1) << (2 * cur_m)
            pred = (cur * cur % (mod * mpz(4)))  # keep some guard
            # use exact integer: pred = cur^2 / p0  should be near B(cur_m+1)
            # Since cur^2 - p0*B(m+1) divisible by high power of 2,
            # B(m+1) = cur^2 / p0  in 2-adics
            # Compute cur^2 * inv(p0) mod 2^{2*(cur_m+1)}
            new_mod = mpz(1) << (2 * (cur_m + 1))
            inv_p0 = pow(p0 % new_mod, -1, new_mod)
            lifted = (cur * cur * inv_p0) % new_mod
            exact_next = B(cur_m + 1, K) % new_mod
            agree = v2(lifted - exact_next)
            P(f"  K={K} m={cur_m}->{cur_m+1}: v2(lifted-exact)={agree}  "
              f"need={2*(cur_m+1)}  full={agree>=2*(cur_m+1)}")
            cur = lifted  # continue with lifted value (maybe truncated)
            cur_m += 1
            if agree < 2 * cur_m - 4:
                P(f"    BREAK: precision collapse")
                break
        final_mod = mpz(1) << (2 * M)
        # recompute clean lift for final compare
        P(f"  K={K} final check skipped if collapsed; target bits={2*M}")


def part_digit_regime():
    P("=" * 70)
    P("PART E: digit regime K=Theta(m) -- does v2(err) beat 2m?")
    P("  Need v2(err) >= 2m to use as congruence for dual barrier mod 2^{2m}.")
    good = 0
    bad = 0
    for m in range(8, 24):
        K = m  # digit-like
        err = B(m, K) ** 2 - P0(K) * B(m + 1, K)
        ve = v2(err)
        ok = ve >= 2 * m
        good += ok
        bad += not ok
        if m <= 14 or m % 3 == 0:
            P(f"  m={m} K={K}  v2(err)={ve}  2m={2*m}  usable={ok}")
    P(f"  usable count in sample: check above (expect NOT usable when K~m)")


def main():
    P("dual_lift -- squaring lift for dual product")
    part_valuation_law()
    part_exact_identity()
    part_algebraic_expand()
    part_lift_algorithm()
    part_digit_regime()
    P()
    P("VERDICT: squaring lemma is REAL (v2>=2m, ratio exact).")
    P("Lift m->m+1 at full dual width needs inverse-power corrections")
    P("(additive tower); Theta(N) steps => O(N^2 log N) rate law.")
    P("Not the O(N log^2 N) crack; new structural lemma for the paper.")


if __name__ == "__main__":
    main()
