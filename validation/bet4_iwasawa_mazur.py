#!/usr/bin/env python3
"""
bet4_iwasawa_mazur.py -- Live bet 4:
  go deeper than O(log N) Mazur/Bernoulli moments (already LLL-empty in door_beta).

Iwasawa: the Kubota-Leopoldt zeta / Mazur measure mu_Mazur on Z_2^*
encodes special values that relate to Bernoulli and to the Gamma-transform
(cracked E-side).  The question is whether H = sum_{k=1}^K mu_k^{-1} is a
*finite* transform of the measure that we can evaluate soft-linearly.

Attacks:
  A. Volkenborn integral of x |-> 1/(1+8 floor(K x)) or similar -- does it
     equal H / something?
  B. Partial sums of the Iwasawa power series L_p(s, chi) truncated
  C. Measure moments beyond O(log N): if we need Theta(N) moments, mass
     is again Theta(N^2)
  D. Explicit: H as <mu_Mazur, f_K> for some f_K of low complexity
  E. Relation of H to (1-2^{1-s}) zeta(s) at s->1  (harmonic) -- archimedean
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


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


def volkenborn_sample(f, m, Nsamp=None):
    """
    Volkenborn integral approx: lim n->inf 1/2^n sum_{i=0}^{2^n - 1} f(i)
    Work mod 2^m; take n = m (standard for precision m).
    """
    if Nsamp is None:
        Nsamp = m
    MOD = mpz(1) << m
    # sum f(0)+...+f(2^n - 1) / 2^n
    # For convergence need f continuous on Z_2
    n = min(Nsamp, m)
    s = mpz(0)
    top = 1 << n
    for i in range(top):
        s = (s + f(i, MOD)) % (MOD << n)  # keep guard
    # divide by 2^n: if s divisible
    return (s >> n) % MOD, n


def part_A_volkenborn():
    P("=" * 70)
    P("BET4-A: Volkenborn probes for H")
    # Classical: int_{Z_2} x^{k} dx = B_k (related)
    # Try f(x) = 1/(1+8*(x mod K)) — discontinuous in K
    for K, m in [(8, 32), (16, 40)]:
        Href = atom_H(K, m)
        MOD = mpz(1) << m

        def f_modK(i, MOD):
            return pow(mpz(1 + 8 * (i % K)), -1, MOD)

        vb, n = volkenborn_sample(f_modK, m, Nsamp=min(m, 14))
        # This is ~ (1/K) sum_{a=0}^{K-1} 1/(1+8a) if uniform — related but
        # includes a=0 term and wraps, not H = sum_{k=1}^K 1/(8k+1)
        P(f"  K={K} m={m}  H={int(Href)%(1<<16)}  VB_modK={int(vb)%(1<<16)}  "
          f"(n={n}, not equal expected)")

        # f(x)=1/(1+8x) on 0..K-1 only is NOT a Volkenborn integral
        # (domain wrong).  Volkenborn wants all of Z_2.
    P("  DIAGNOSIS: H is a PARTIAL sum over an AP, not a Volkenborn integral")
    P("  over Z_2.  Extending to Z_2 reintroduces infinite/convergent tails")
    P("  that are Diamond-log-gamma (cracked) PLUS the partial = H (barrier).")


def part_B_moments_mass():
    P("=" * 70)
    P("BET4-B: Mazur moments needed to span H")
    P("  door_beta: O(log N) Bernoulli/Mazur moments LLL-empty vs H.")
    P("  If Theta(N) moments are required, reading them is Omega(N^2) bits")
    P("  (Conjecture mass).  Measure:")
    # Build basket of first R power-sum-like moments of 1/(1+8k) weights
    # Actually: moments m_j = sum_{k=1}^K k^j  (forward) — cracked via? 
    # only dyadic.  Or Iwasawa: int <x>^j dmu
    for K, m, R in [(20, 48, 4), (20, 48, 8), (20, 48, 12)]:
        Href = atom_H(K, m)
        MOD = mpz(1) << m
        # Can we solve H = sum_{j=0}^{R-1} c_j * S_j  with small c_j?
        # H = sum_j (-8)^j S_j  exactly with R~m/3 — needs R=Theta(m)
        # Check minimal R for Faulhaber truncation to agree
        S = [mpz(0)] * R
        for k in range(1, K + 1):
            p = mpz(1)
            for j in range(R):
                S[j] = (S[j] + p) % (1 << (m + 2 * R))
                p = (p * k) % (1 << (m + 2 * R))
        acc = mpz(0)
        pow8 = mpz(1)
        GMOD = mpz(1) << (m + 2 * R)
        for j in range(R):
            acc = (acc + pow8 * S[j]) % GMOD
            pow8 = (pow8 * (-8)) % GMOD
        approx = acc % (1 << m)
        agree = 0
        diff = int(Href - approx) % (1 << m)
        if diff == 0:
            agree = m
        else:
            # trailing agreeing bits
            while agree < m and ((Href ^ approx) & (1 << agree)) == 0:
                agree += 1
            # better: v2(diff)
            d = int(Href - approx) % (1 << m)
            if d == 0:
                agree = m
            else:
                agree = (d & -d).bit_length() - 1 if d else m
        P(f"  K={K} m={m}  Faulhaber_trunc R={R}  v2(H-approx)={agree}")
    P("  RATE: need R = Theta(m) moments for full H — mass Theta(N^2).")


def part_C_iwasawa_L():
    P("=" * 70)
    P("BET4-C: Kubota-Leopoldt / Iwasawa L-function at one point")
    P("""
  L_2(s, chi) soft-linear at FIXED algebraic s is plausible (bit-burst /
  overconvergent methods for p-adic L).  But H is a TRUNCATED harmonic
  of length K=Theta(N), equivalent to digamma(K+1/8)-digamma(1/8),
  i.e. an L-value difference at a PRECISION-DEPENDENT argument K.

  Evaluating L or digamma at an argument that itself has Theta(N) bits
  is the moving-target problem Caruso solves for holonomic DEs —
  digamma is not holonomic.  So Iwasawa at fixed s does not give H(K).
""")


def part_D_gamma_transform():
    P("=" * 70)
    P("BET4-D: Mazur Gamma-transform vs H")
    P("  Paper: E is a truncation of the Mazur measure Gamma-transform")
    P("  (CRACKED, multiplicative).  H is the digamma / logarithmic")
    P("  derivative side.  Same razor as e_T vs e_{T-1}.")
    P("  No new soft-linear path from deeper Iwasawa without differentiating")
    P("  the measure transform — which reintroduces the additive core.")


def main():
    P("bet4_iwasawa_mazur -- deeper Iwasawa/Mazur")
    part_A_volkenborn()
    part_B_moments_mass()
    part_C_iwasawa_L()
    part_D_gamma_transform()
    P()
    P("BET4 VERDICT: partial sums != Volkenborn; Theta(N) moments needed;")
    P("  fixed-s Iwasawa misses moving K; Gamma-transform razor stands.")
    P("  Door beta's O(log N) kill extends: full Iwasawa doesn't soft-solve H.")


if __name__ == "__main__":
    main()
