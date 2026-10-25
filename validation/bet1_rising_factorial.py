#!/usr/bin/env python3
"""
bet1_rising_factorial.py -- Live bet 1:
  soft-linear short rising factorial at height a = 2^m

    f(K) = (a)_K^{(step 8)} := prod_{t=0}^{K-1} (a + 8t)
         = prod_{k=1}^{K} (2^m + 8k + 1)   [shift]   = B(m,K) up to index

Digit regime: K,m = Theta(N).  Need f(K) mod 2^{2m}.

Attacks in this file:
  A. Nilpotent binomial expansion of each factor around 2^m
  B. Write as P0 * prod (1 + R/mu) and study prod via Log series
     truncated with cracked dyadic power sums of mu^{-1} -- circular for Q1
  C. Polynomial method: (a + 8x)_K as value of rising-factorial poly
     built from Stirling1 * a^k -- cost of Stirling1 row
  D. Recurrence / transfer scheme using dual squaring lemma as a
     *downward* reduction in m at fixed K (precision audit)
  E. Multiplicative decomposition f = oddpart( factorial_ratio )
     and measure 2-valuation collapse
"""
from __future__ import annotations
import math
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
    if n == 0:
        return 9999
    if n < 0:
        n = -n
    return (n & -n).bit_length() - 1


def B(m, K, MOD=None):
    R = mpz(1) << m
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= R + 8 * k + 1
        if MOD is not None:
            acc %= MOD
    return acc


def P0(K, MOD=None):
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= 8 * k + 1
        if MOD is not None:
            acc %= MOD
    return acc


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


# ---------------------------------------------------------------------------
# A. Nilpotent expansion
# ---------------------------------------------------------------------------
def part_A_nilpotent():
    P("=" * 70)
    P("BET1-A: Log B = sum Log(R+mu) = K Log-ish + sum_j (-1)^{j+1} R^j/j Q_j")
    P("  where Q_j = sum mu^{-j}.  Mod R^2 only Q_1=H survives.")
    P("  Soft-linear B mod R^2 <=> soft-linear H.  Tautology.")
    for m, K in [(16, 16), (20, 20)]:
        R = mpz(1) << m
        MOD = R * R
        Bm = B(m, K, MOD)
        p0 = P0(K, MOD)
        H = atom_H(K, m)
        pred = (p0 * (1 + R * H)) % MOD
        P(f"  m={m} K={K}  dual_id={Bm==pred}")


# ---------------------------------------------------------------------------
# B. Stirling1 rising factorial poly at a=R+9
# ---------------------------------------------------------------------------
def stirling1_signed_row(n):
    """Unsigned cycle numbers c(n,k); signed s(n,k)=(-1)^{n-k} c(n,k).
    Recurrence: c(n,k) = c(n-1,k-1) + (n-1)*c(n-1,k)."""
    # return list s[0..n] for fixed n: signed stirling 1st kind s(n,k)
    prev = [mpz(1)]  # s(0,0)=1
    for nn in range(1, n + 1):
        cur = [mpz(0)] * (nn + 1)
        for k in range(nn + 1):
            a = prev[k - 1] if k - 1 >= 0 and k - 1 < len(prev) else mpz(0)
            b = prev[k] if k < len(prev) else mpz(0)
            cur[k] = a - (nn - 1) * b
        prev = cur
    return prev


def part_B_stirling():
    P("=" * 70)
    P("BET1-B: B(m,K) = 8^K * ( (R+9)/8 )_K = sum_k s(K,k) * (R+9)^k * ?")
    P("  Rising (x)_K = sum_k s(K,k) x^k  with x=(R+9)/8 not in Z.")
    P("  Integer: prod_{t=0}^{K-1}(A+8t) with A=R+9 = sum c_j A^j.")
    P("  Coeffs = 8^{K-j} * (signed stirling related). Cost of full row: Theta(K^2).")
    for K in [8, 12]:
        srow = stirling1_signed_row(K)
        # verify (x)_K = sum s(K,k) x^k for x=5
        x = mpz(5)
        lhs = mpz(1)
        for t in range(K):
            lhs *= x + t
        rhs = mpz(0)
        xp = mpz(1)
        # s(K,k) for k=0..K; (x)_K = sum_k |s| with sign
        # standard: x(x+1)...(x+K-1) = sum_k s(K,k) * x^k ? 
        # Actually rising x^{(K)} = sum S1_unsigned * x^k with signs:
        # x(x+1)...(x+n-1) = sum_k s(n,k) x^k  where s = signed stirling 1st
        for k in range(K + 1):
            rhs += srow[k] * (x ** k)
        P(f"  K={K} rising_check (x=5): {lhs==rhs}")
    P("  Evaluating sum_k s(K,k) A^k mod 2^{2m} still needs the Stirling row")
    P("  (Theta(K^2) or soft-K^2) — no compression for short K=Theta(m).")


# ---------------------------------------------------------------------------
# C. Downward dual-square reduction
# ---------------------------------------------------------------------------
def part_C_downward():
    P("=" * 70)
    P("BET1-C: reduce B(m) from B(m-1) via dual squaring (precision margin)")
    # B(m-1)^2 = P0 * B(m) + O(2^{2(m-1)})
    # so B(m) = B(m-1)^2 / P0  mod 2^{2m-2}  -- 2 bits short of 2m
    for m in range(10, 18):
        K = m
        MOD = mpz(1) << (2 * m)
        p0 = P0(K)
        prev = B(m - 1, K)
        pred = (prev * prev * pow(p0 % MOD, -1, MOD)) % MOD
        exact = B(m, K) % MOD
        ag = v2(pred - exact)
        P(f"  m={m} K={K}  agree_v2={ag}  need={2*m}  short={2*m-ag}")


# ---------------------------------------------------------------------------
# D. Factorial ratio / odd-part
# ---------------------------------------------------------------------------
def odd_part(n):
    """Remove factors of 2 from n! roughly: product of odds <= n."""
    acc = mpz(1)
    for i in range(1, n + 1, 2):
        acc *= i
    return acc


def part_D_factorial():
    P("=" * 70)
    P("BET1-D: express B as ratio of odd parts near 2^m")
    for m, K in [(8, 4), (10, 5)]:
        R = mpz(1) << m
        # numbers R+9, R+17, ..., R+8K+1
        # among odds from 1 to R+8K+1, the residue-1-mod-8 ones in range
        Bm = B(m, K)
        # crude: odd_part(R+8K+1) / odd_part(R+1) includes ALL odd residues
        top = R + 8 * K + 1
        bot = R + 1
        # Can't divide factorials mod 2^{2m} -- both eventually 0 for large
        v_top = sum(top // (1 << i) for i in range(1, top.bit_length() + 2))
        P(f"  m={m} K={K}  top={int(top)}  v2((top)!)~={v_top}  "
          f"modbits={2*m}  vanishes={v_top>=2*m}")
    P("  DIAGNOSIS: for m large, factorial at height 2^m has v2 >> 2m.")
    P("  Odd-part isolation is the whole problem (structured vs random odds).")


# ---------------------------------------------------------------------------
# E. Structured product via CRT over odd primes (rational H route)
# ---------------------------------------------------------------------------
def part_E_rational_CRT():
    P("=" * 70)
    P("BET1-E: exact rational H = e_{T-1}/P0 then reduce mod 2^m")
    P("  Numerator bit length ~ K log K = Theta(N^2) bits before reduction.")
    for K in [20, 40, 80]:
        # bit length of P0 and of leave-one-out numerator
        p0 = P0(K)
        # e_{T-1} = P0 * H_rational; compute leave-one-out exactly
        num = mpz(0)
        for i in range(1, K + 1):
            term = mpz(1)
            for j in range(1, K + 1):
                if j != i:
                    term *= 8 * j + 1
            num += term
        P(f"  K={K}  bits(P0)={p0.bit_length()}  bits(num)={num.bit_length()}  "
          f"bits(K*logK)~={int(K*math.log2(max(K,2))+1)}")
    P("  CRT of H mod many odd primes to recover exact rational costs Omega(N^2).")
    P("  No shortcut to H mod 2^m without the 2-adic problem.")


# ---------------------------------------------------------------------------
# F. Fresh attempt: prefix products via dyadic decomposition of K
#    B(m,K) = B(m,K1)*B_shifted(m,K2) when K=K1+K2
#    Shifted: prod_{k=K1+1}^{K1+K2} (R+8k+1) = prod_{j=1}^{K2} (R+8K1+8j+1)
#            = B but with R' = R+8K1 — NOT a pure power of 2.
#    If 8K1 is small, maybe Postnikov-compose with cracked methods?
# ---------------------------------------------------------------------------
def part_F_shift_compose():
    P("=" * 70)
    P("BET1-F: dyadic split K with shifted height R+8K1")
    for m, K in [(16, 16), (20, 20)]:
        K1 = K // 2
        R = mpz(1) << m
        MOD = R * R
        left = B(m, K1, MOD)
        # right = prod_{k=K1+1}^K (R+8k+1)
        right = mpz(1)
        for k in range(K1 + 1, K + 1):
            right = (right * (R + 8 * k + 1)) % MOD
        full = B(m, K, MOD)
        P(f"  m={m} K={K}  split_ok={(left*right)%MOD==full}")
        # Is right equal to some B(m',K')? Height is R+8(K1+1)+... not 2^{m'}
        shift = 8 * K1
        # right = prod_{j=1}^{K-K1} (R+shift+8j+1)
        # Compare to prod (2^m + shift + 8j + 1)
        P(f"    shift={shift}  v2(shift)={v2(shift)}  vs m={m}")
        P(f"    shift not a high power of 2 => dual-square lemma does not apply")
        P(f"    directly; Postnikov rate: admissible shift O(poly N) yields")
        P(f"    only O(log N) bits (door beta).  Here shift=Theta(K)=Theta(m).")


def main():
    P("bet1_rising_factorial -- short rising factorial at height 2^m")
    part_A_nilpotent()
    part_B_stirling()
    part_C_downward()
    part_D_factorial()
    part_E_rational_CRT()
    part_F_shift_compose()
    P()
    P("BET1 VERDICT: no soft-linear rising-factorial algorithm found.")
    P("  Dual id tautological; Stirling row quadratic; downward lift 2 bits short;")
    P("  factorial vanishes; CRT Omega(N^2); shifted split leaves door-beta regime.")


if __name__ == "__main__":
    main()
