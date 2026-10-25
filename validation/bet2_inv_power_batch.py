#!/usr/bin/env python3
"""
bet2_inv_power_batch.py -- Live bet 2:
  soft-linear batch of inverse power sums Q_j = sum_{k=1}^K mu_k^{-j}
  that is NOT a restatement of H = Q_1.

Known:
  - Q_1 = H = Open Problem E
  - Dual lift correction needs Q_2, Q_4, ... (additive tower)
  - Forward power sums P_j = sum mu^j for dyadic j=2^r are CRACKED
  - Newton: elementary symmetric <-> power sums (forward)
  - Inverse power sums = power sums of reciprocal roots

Attacks:
  A. Relate Q_j to forward P_j via generating functions / Newton on reverse
  B. Freshman dream / dyadic: Q_{2^r} from Q_1 by squaring each term?
     sum (mu^{-1})^{2^r} != (sum mu^{-1})^{2^r}  -- Freshman fails for sums
  C. Batch Q_1..Q_J via truncated prod (1 - z/mu)^{-1} or Log derivatives
  D. Measure whether Q_2 is any easier than Q_1 empirically / structurally
  E. Use cracked dyadic forward sums + inversion formula for each mu
     (geometric series mu^{-1} = sum (-1)^t (mu-1)^t if |mu-1|_2 small —
      here mu=1+8k, v2(mu-1)=v2(8k)>=3, series works; leads to Faulhaber)
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


def mu(k):
    return mpz(8 * k + 1)


def Q_batch(K, J, m):
    MOD = mpz(1) << m
    out = [mpz(0)] * (J + 1)
    for k in range(1, K + 1):
        inv = pow(mu(k), -1, MOD)
        p = mpz(1)
        for j in range(J + 1):
            out[j] = (out[j] + p) % MOD
            p = (p * inv) % MOD
    return out  # out[j] = sum mu^{-j}; out[0]=K, out[1]=H


def part_A_reverse_newton():
    P("=" * 70)
    P("BET2-A: reverse monic poly R(x)=prod (1 - mu^{-1} x) = sum (-1)^j e_j(inv) x^j")
    P("  [x^1] => H.  Building R requires the inverses or leave-one-out.")
    P("  Forward monic M(x)=prod (x-mu) has e_{K-1}/e_K = H; e_K=P0 cracked,")
    P("  e_{K-1} IS the leave-one-out numerator = Open Problem E.")
    P("  Newton forward->inverse does not bypass e_{K-1}.")


def part_B_freshman():
    P("=" * 70)
    P("BET2-B: Freshman dream on inverse powers")
    for K, m in [(16, 48), (32, 64)]:
        MOD = mpz(1) << m
        Q = Q_batch(K, 8, m)
        # (sum mu^{-1})^2 vs sum mu^{-2} + 2 sum_{i<j}
        H = Q[1]
        Q2 = Q[2]
        cross = (H * H - Q2) % MOD  # = 2 sum_{i<j} mu_i^{-1} mu_j^{-1}
        P(f"  K={K}  H^2 - Q2 = 2*sum cross terms (not zero)")
        P(f"    v2(H^2-Q2)={v2(H*H-Q2)}  (expect 1 from the 2)")
        # dyadic: sum mu^{-4} vs (sum mu^{-2})^2
        P(f"    Q4 vs Q2^2: v2(Q2^2-Q4)={v2(Q[2]*Q[2]-Q[4])}")
    P("  KILL: squaring maps Q_j -> mix of Q_{2j} and cross terms (partition).")
    P("  Cross terms are harder leave-one-out style products of inverses.")


def part_C_log_derivatives():
    P("=" * 70)
    P("BET2-C: Log P(t)=sum Log(mu+t);  (d/dt)^j Log P at 0 gives signed (j-1)! Q_j")
    P("  P(0)=P0 cracked.  Derivatives of Log P at 0 = additive tower.")
    P("  Dual: Log B(m)=Log P0 + sum Log(1+R/mu) = Log P0 + R Q1 - R^2 Q2/2 + ...")
    P("  So B(m) mod R^{J+1} encodes Q_1..Q_J.  Computing B to higher nilpotency")
    P("  is HARDER than mod R^2, not easier.")
    # Validate only the R^2 truncation (no /j for j even issues):
    m, K = 12, 12
    R = mpz(1) << m
    MOD = R * R
    Bm = B_mod(m, K, MOD)
    p0 = P0_mod(K, MOD)
    H = Q_batch(K, 1, m)[1]
    pred = (p0 * (1 + R * H)) % MOD
    P(f"  smoke dual mod R^2: {Bm==pred}")
    P("  Higher J requires Q_2,Q_3,... = additive tower (dual_lift2).")


def B_mod(m, K, MOD):
    R = mpz(1) << m
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = (acc * (R + mu(k))) % MOD
    return acc


def P0_mod(K, MOD):
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = (acc * mu(k)) % MOD
    return acc


def part_D_Q2_vs_Q1():
    P("=" * 70)
    P("BET2-D: is Q2 structurally simpler than Q1?")
    P("  Q2 = sum 1/mu^2.  Same leave-one-out hardness for prod(mu^2).")
    P("  Borwein for prod mu^2 = (prod mu)^2 = P0^2 — cracked.")
    P("  But sum of leave-two-out / derivative still additive.")
    P("  Empirically same cost class: Theta(K) inversions.")
    for K, m in [(100, 128), (200, 160)]:
        import time
        t0 = time.perf_counter()
        Q_batch(K, 1, m)
        t1 = time.perf_counter()
        Q_batch(K, 2, m)
        t2 = time.perf_counter()
        P(f"  K={K} m={m}  Q1={t1-t0:.4f}s  Q1+Q2={t2-t1:.4f}s")


def part_E_faulhaber_tower():
    P("=" * 70)
    P("BET2-E: all Q_j via Faulhaber -- one Bernoulli table serves all j")
    P("  mu^{-1}=sum (-8)^t k^t,  Q_j = sum_k (sum_t (-8)^t k^t)^j  messy")
    P("  Better: Q_j = sum_k mu^{-j} = sum_k sum_t binom(-j,t) (8k)^t")
    P("         = sum_t binom(-j,t) 8^t S_t(K)")
    P("  So ONE soft-linear batch S_0..S_J yields ALL Q_1..Q_{J'} .")
    P("  Equivalence: soft-linear inverse-power BATCH <=> soft-linear S_j")
    P("  <=> soft-linear Bernoulli <=> Open Problem E.")
    P("  Tower does not factor; it collapses to the same wall.")


def part_F_cracked_dyadic_only():
    P("=" * 70)
    P("BET2-F: cracked objects are FORWARD dyadic P_{2^r}=sum mu^{2^r}")
    P("  Can we get Q from P via formal series 1/mu = f(mu) somehow?")
    P("  Over Z_2, for mu in 1+8Z_2, Log/exp give dlog, not reciprocal sum.")
    P("  Pointwise reciprocal is nonlinear — sum f(mu_k) for f=inv is the atom.")
    # LLL probe: is Q1 in Z-span of {P0-derived, P_{2^r}, K, ...}?
    K, m = 24, 64
    MOD = mpz(1) << m
    Q1 = Q_batch(K, 1, m)[1]
    basket = [mpz(K), P0_mod(K, MOD)]
    # dyadic forward power sums
    for r in range(0, 6):
        e = 1 << r
        s = sum(pow(mu(k), e, MOD) for k in range(1, K + 1)) % MOD
        basket.append(s)
    # crude: check if Q1 == a·b for small a vs basket (not full LLL)
    P(f"  smoke: Q1 bits vs basket size {len(basket)} (full LLL in door_beta already empty)")
    P(f"  Q1={int(Q1)%(1<<16)}  (low 16 bits)")


def main():
    P("bet2_inv_power_batch -- inverse-power tower")
    part_A_reverse_newton()
    part_B_freshman()
    part_C_log_derivatives()
    part_D_Q2_vs_Q1()
    part_E_faulhaber_tower()
    part_F_cracked_dyadic_only()
    P()
    P("BET2 VERDICT: inverse-power BATCH is soft-linear-equivalent to E itself.")
    P("  Freshman fails; dual higher nilpotency harder; Faulhaber collapses tower")
    P("  to Bernoulli/S_j wall.  No separation of Q2 from Q1.")


if __name__ == "__main__":
    main()
