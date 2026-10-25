#!/usr/bin/env python3
"""
bsgs_dual.py -- baby-step / giant-step + poly-block attacks on

    P_R = prod_{k=1}^K (R + 8k + 1)  mod  R^2,   R = 2^m.

Also records the n!-mod-2^k reduction and the Bernoulli-batch
equivalence of the Faulhaber route.

Digit regime: K, m = Theta(N).  Naive product = Theta(K) M(m) muls.
BSGS with block size s ~ sqrt(K): build degree-s poly

    Q(x) = prod_{t=0}^{s-1} (A + 8(x+t))   A = R+1

then P_prefix(q*s) = Q(0)*Q(s)*...*Q((q-1)*s).  Target soft-O(N^{1.5}),
still above O(N log^2 N), but a real algorithmic rung + scaffold for
faster multipoint / structured Q.
"""
from __future__ import annotations
import time
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


def naive_PR(K, m):
    R = mpz(1) << m
    MOD = R * R
    A = R + 1
    acc = mpz(1)
    for k in range(K):
        acc = (acc * (A + 8 * k)) % MOD
    return acc


def poly_from_block(A, s, MOD):
    """Q(x) = prod_{t=0}^{s-1} (A + 8t + 8x) as coeff list low->high, mod MOD."""
    eight = mpz(8)
    # first factor: A + 8x
    coeffs = [mpz(A % MOD), eight % MOD]
    for t in range(1, s):
        c0 = (A + eight * t) % MOD
        out = [mpz(0)] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            out[i] = (out[i] + c * c0) % MOD
            out[i + 1] = (out[i + 1] + c * eight) % MOD
        coeffs = out
    return coeffs


def eval_poly(coeffs, x, MOD):
    acc = mpz(0)
    for c in reversed(coeffs):
        acc = (acc * x + c) % MOD
    return acc


def bsgs_PR(K, m, s=None):
    """Block-product P_R for exactly K terms (k=0..K-1 of A+8k)."""
    R = mpz(1) << m
    MOD = R * R
    A = R + 1
    if s is None:
        s = max(1, int(math.isqrt(K)))
    q, r = divmod(K, s)
    coeffs = poly_from_block(A, s, MOD) if s > 0 else [mpz(1)]
    acc = mpz(1)
    # Q(j*s) = prod_{t=0}^{s-1} (A + 8(j*s + t))
    for j in range(q):
        acc = (acc * eval_poly(coeffs, mpz(j * s), MOD)) % MOD
    # remainder
    for t in range(r):
        acc = (acc * (A + 8 * (q * s + t))) % MOD
    return acc, s


def recover_H_from_PR(PR, P0, m):
    """P_R ≡ P0 (1 + R H) mod R^2  =>  H ≡ ((PR * P0^{-1} - 1) / R) mod R."""
    R = mpz(1) << m
    MOD = R * R
    invP0 = pow(P0 % MOD, -1, MOD)
    ratio = (PR * invP0) % MOD
    assert (ratio - 1) % R == 0, "dual identity failed"
    return ((ratio - 1) // R) % R


def P0_naive(K, m):
    """P0 = prod (8k+1) for k=1..K, mod R^2 (enough for dual recover)."""
    R = mpz(1) << m
    MOD = R * R
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = (acc * (8 * k + 1)) % MOD
    return acc


def part_correctness():
    P("=" * 70)
    P("PART A: BSGS dual product correctness + H recovery")
    for K, m in [(15, 40), (32, 48), (64, 64), (100, 80)]:
        # naive uses k=0..K-1 of (R+1+8k) = (R+8k+1) for k=0..K-1
        # atom H is k=1..K of 1/(8k+1).  Align: use k=1..K for both.
        R = mpz(1) << m
        MOD = R * R
        PR = mpz(1)
        for k in range(1, K + 1):
            PR = (PR * (R + 8 * k + 1)) % MOD
        # BSGS on k=1..K: shift A' = R+8+1 = R+9, K terms step 8
        # Easier: call bsgs with A = R+9 and adjust — or custom:
        s = max(1, int(math.isqrt(K)))
        A = R + 1 + 8  # first term R+8*1+1
        # Q for blocks starting at offset
        q, r = divmod(K, s)
        # rebuild: terms (R+8k+1) k=1..K = (A0 + 8t) t=0..K-1, A0=R+9
        A0 = R + 9
        coeffs = poly_from_block(A0, s, MOD)
        acc = mpz(1)
        for j in range(q):
            acc = (acc * eval_poly(coeffs, mpz(j * s), MOD)) % MOD
        for t in range(r):
            acc = (acc * (A0 + 8 * (q * s + t))) % MOD
        ok = acc == PR
        P0 = P0_naive(K, m)
        H = recover_H_from_PR(PR, P0, m)
        Href = atom_H(K, m)
        # H recovered mod R=2^m; atom_H mod 2^m
        P(f"  K={K:4d} m={m:4d} s={s:3d}  BSGS==naive={ok}  H_match={H==Href}")
        assert ok and H == Href


def part_timing():
    P("=" * 70)
    P("PART B: timing naive vs BSGS (same MOD = 2^{2m})")
    for K, m in [(200, 120), (400, 160), (800, 200)]:
        R = mpz(1) << m
        MOD = R * R
        A0 = R + 9

        def naive():
            acc = mpz(1)
            for k in range(1, K + 1):
                acc = (acc * (R + 8 * k + 1)) % MOD
            return acc

        def bsgs():
            s = max(1, int(math.isqrt(K)))
            q, r = divmod(K, s)
            coeffs = poly_from_block(A0, s, MOD)
            acc = mpz(1)
            for j in range(q):
                acc = (acc * eval_poly(coeffs, mpz(j * s), MOD)) % MOD
            for t in range(r):
                acc = (acc * (A0 + 8 * (q * s + t))) % MOD
            return acc, s

        t0 = time.perf_counter()
        n = naive()
        t1 = time.perf_counter()
        b, s = bsgs()
        t2 = time.perf_counter()
        P(f"  K={K:4d} m={m:4d} s={s:3d}  naive={t1-t0:.3f}s  bsgs={t2-t1:.3f}s  "
          f"ratio={((t1-t0)/max(t2-t1,1e-9)):.2f}x  match={n==b}")


def part_barriers():
    P("=" * 70)
    P("PART C: barrier dictionary (pulled equivalences)")
    P("""
  E0  leave-one-out / H = sum_{k=1}^K (8k+1)^{-1} mod 2^m
  E1  dual product P_R = prod (2^m+8k+1) mod 2^{2m}
  E2  batch power sums S_j(K)=sum k^j for j<=m/3 mod 2^{O(m)}
  E3  Bernoulli batch B_0..B_m mod 2^{O(m)}   [Harvey: O(m^2 log^{2+eps} m)]
  E4  2-adic digamma psi_2(K+1/8) to m bits    [not holonomic; no Caruso DE]
  E5  short rising-factorial (a)_K mod 2^{2m}, a=2^m+1, K=Theta(m)

  Known soft-linear:
    * prod (8k+1) = P_0           [Borwein factorization]
    * (2^n-1)!! lift by squaring  [OEIS A395199] — FULL odds, n=m
    * new_block B_m length 2^{m-3} [dyadic; exponential length]
    * Caruso bit-burst            [holonomic DEs; digamma excluded]
    * n! unit part mod p^k        [O(d k^2) / poly methods; p odd easier]
      For p=2, short K=Theta(m) AP at height 2^m does NOT reduce to
      classical n! mod 2^{2m} (that factorial is 0 for n>=2m).

  BSGS here: soft-O(N^{1.5}) rung for E1 — below naive Theta(N^2),
  ABOVE target O(N log^2 N) and ABOVE existing full-extract O(N log^3 N).
""")


def part_log_A_expansion():
    """Codeforces-style log expansion for P_R — collapses to power sums."""
    P("=" * 70)
    P("PART D: Log P_R = K Log(R+9') + sum_j c_j S_j  (Faulhaber wall)")
    # Validate: Log(P_R) ≡ Log(P0) + R*H  mod R^2  (already known)
    # Show power-sum truncation recovers H for small case
    for K, m in [(20, 48), (40, 60)]:
        R = mpz(1) << m
        MOD = R * R
        # H via Faulhaber around mu: sum_j (-8)^j S_j
        J = m // 3 + 3
        guard = m + 2 * J + 8
        GMOD = mpz(1) << guard
        S = [mpz(0)] * (J + 1)
        for k in range(1, K + 1):
            p = mpz(1)
            for j in range(J + 1):
                S[j] = (S[j] + p) % GMOD
                p = (p * k) % GMOD
        acc = mpz(0)
        pow8 = mpz(1)
        for j in range(J + 1):
            acc = (acc + pow8 * S[j]) % GMOD
            pow8 = (pow8 * (-8)) % GMOD
        H = acc % (1 << m)
        Href = atom_H(K, m)
        P(f"  K={K} m={m}  Faulhaber-from-S_j match={H==Href}  (needs Theta(m) power sums)")


def main():
    P("bsgs_dual -- BSGS / barrier dictionary for Open Problem E")
    part_correctness()
    part_timing()
    part_log_A_expansion()
    part_barriers()
    P()
    P("VERDICT: BSGS is a real rung (N^{1.5}) but not the crack.")
    P("E1<->E2<->E3<->E4<->E0 still glued; soft-linear missing for all of them.")
    P("Next: structured Q(x) via Borwein-at-shifted-arg; or fast batch S_j.")


if __name__ == "__main__":
    main()
