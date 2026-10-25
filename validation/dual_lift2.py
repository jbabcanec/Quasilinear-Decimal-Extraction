#!/usr/bin/env python3
"""
dual_lift2.py -- corrected lift using next term of the expansion.

LEMMA.  For R=2^m, mu_k=8k+1,
  B(m,K)^2 / (P0(K) B(m+1,K)) = prod_k (1 + R^2/(mu_k(mu_k+2R)))
  and each summand has v2 = 2m, so v2(ratio-1) >= 2m.

Expansion mod 2^{2m+c}:
  ratio = 1 + R^2 * sum 1/(mu(mu+2R)) + O(R^4)
        = 1 + R^2 * sum mu^{-2} + O(2^{3m})   (since mu+2R = mu*(1+2R/mu),
          1/(mu+2R)=mu^{-1}(1+2R/mu)^{-1}=mu^{-1}(1-2R/mu+O(R^2/mu^2)),
          so 1/(mu(mu+2R))=mu^{-2}(1-2R/mu+O(R^2)))

Since mu=1 mod 8, mu^{-2}=1 mod 8, sum_{k=1}^K mu^{-2} ≡ K mod 8.

CORRECTED LIFT:
  B(m+1) ~= B(m)^2 / P0 / (1 + R^2 K)  mod 2^{2m+3}
           = B(m)^2 / P0 * (1 - R^2 K)  mod 2^{2m+3}

If this reaches full 2(m+1) bits, we can lift m->m+1 preserving the
dual-barrier modulus ratio -- still Theta(N) steps => O(N^2 log N)
unless we can geometrically jump.  But first verify bit-exact correction.
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


def part_lemma_check():
    P("=" * 70)
    P("PART A: ratio = prod (1 + R^2/(mu(mu+2R))) exact")
    for m, K in [(8, 8), (12, 10), (16, 16)]:
        R = mpz(1) << m
        num = B(m, K) ** 2
        den = P0(K) * B(m + 1, K)
        # compute prod as rational / integer with common structure
        # ratio_num / ratio_den for prod (mu(mu+2R) + R^2) / (mu(mu+2R))
        pn, pd = mpz(1), mpz(1)
        for k in range(1, K + 1):
            mu = 8 * k + 1
            d = mu * (mu + 2 * R)
            pn *= d + R * R
            pd *= d
        # num/den should equal pn/pd
        ok = num * pd == den * pn
        P(f"  m={m} K={K}  lemma_exact={ok}  v2(num-den)={v2(num-den)}")


def part_corrected_lift():
    P("=" * 70)
    P("PART B: corrected lift B(m+1) <- B(m)^2/P0 * (1 - K*R^2)")
    for m in range(6, 22):
        for K in [m, 2 * m, 8, 16]:
            R = mpz(1) << m
            need = 2 * (m + 1)
            MOD = mpz(1) << (need + 4)  # guard
            p0 = P0(K)
            Bm = B(m, K)
            inv_p0 = pow(p0 % MOD, -1, MOD)
            raw = (Bm * Bm * inv_p0) % MOD
            corr = (raw * (1 - K * R * R)) % MOD
            # also try corr2 with sum mu^{-2} exact mod 2^{need}
            s2 = sum(pow(mpz(8 * k + 1), -2, MOD) for k in range(1, K + 1)) % MOD
            corr2 = (raw * (1 - s2 * R * R)) % MOD
            # higher: 1/(mu(mu+2R)) exactly summed
            s_exact = mpz(0)
            for k in range(1, K + 1):
                mu = mpz(8 * k + 1)
                s_exact += pow(mu * (mu + 2 * R), -1, MOD)
            s_exact %= MOD
            corr3 = (raw * (1 - s_exact * R * R)) % MOD
            exact = B(m + 1, K) % MOD
            P(
                f"  m={m:2d} K={K:3d}  raw_v2={v2(raw-exact):4d}  "
                f"K-corr={v2(corr-exact):4d}  "
                f"s2-corr={v2(corr2-exact):4d}  "
                f"sex-corr={v2(corr3-exact):4d}  need={need}"
            )


def part_lift_chain():
    """Chain corrected lifts from m0 to M at fixed K, precision 2M throughout."""
    P("=" * 70)
    P("PART C: chain lift with s2-correction to target precision")
    for K, M in [(8, 24), (16, 28), (32, 32)]:
        m0 = max(4, K.bit_length() + 2)
        MOD = mpz(1) << (2 * M + 8)
        cur = B(m0, K) % MOD
        p0 = P0(K)
        inv_p0 = pow(p0 % MOD, -1, MOD)
        s2 = sum(pow(mpz(8 * k + 1), -2, MOD) for k in range(1, K + 1)) % MOD
        ok_all = True
        for m in range(m0, M):
            R = mpz(1) << m
            raw = (cur * cur * inv_p0) % MOD
            cur = (raw * (1 - s2 * R * R)) % MOD
            exact = B(m + 1, K) % MOD
            # compare at precision 2*(m+1)
            cmp_mod = mpz(1) << (2 * (m + 1))
            if (cur - exact) % cmp_mod != 0:
                # measure agreement
                ag = v2(cur - exact)
                P(f"  K={K} fail at m={m}->{m+1}: agree_v2={ag} need={2*(m+1)}")
                ok_all = False
                break
        if ok_all:
            final_ag = v2(cur - (B(M, K) % MOD))
            P(f"  K={K} M={M} m0={m0}  CHAIN OK  final_v2(diff)={final_ag}  "
              f"need={2*M}")


def part_complexity():
    P("=" * 70)
    P("PART D: complexity even if corrected lift works")
    P("""
  Per +1 step in m: O(M(N)) at target precision 2N.
  Steps: Theta(N).  Total: O(N M(N)) = O(N^2 log N).

  Soft-linear would need geometric jumps (m -> 2m) with a doubling
  identity that preserves valuation margin.  Empirically the +1
  identity's surplus valuation is O(1) bits (often exactly 0--2 above
  2m), not Theta(m).  Doubling jumps are not supported at full width.

  Same rate-law family as door beta / Postnikov / admissible t=2^r.
""")


def main():
    P("dual_lift2 -- corrected dual squaring lift")
    part_lemma_check()
    part_corrected_lift()
    part_lift_chain()
    part_complexity()


if __name__ == "__main__":
    main()
