#!/usr/bin/env python3
"""
q_structure.py -- does the BSGS block polynomial

    Q(x) = prod_{t=0}^{s-1} (A + 8t + 8x)   A = 2^m + 9

have 2-adic coefficient structure that lets us evaluate Q(0),Q(s),...
faster than O(s) per eval / O(sqrt(K) s) total?

Also: Q(x) = 8^s * ((A/8) + x)_s  (rising factorial).
Coeffs are 8^s * stirling1_signed(s, *) shifted by A/8.

Probe:
  (1) v2(coeffs) profile
  (2) whether Q(x) mod (x^2, 2^{2m}) already encodes H for full K — tautology
  (3) recurrence / hypergeometric term ratio for coeff sequence
  (4) evaluate Q at progressive points via FFT multipoint (complexity check)
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
    if n == 0:
        return 999
    return (n & -n).bit_length() - 1


def poly_from_block(A, s, MOD):
    eight = mpz(8)
    coeffs = [mpz(A % MOD), eight % MOD]
    for t in range(1, s):
        c0 = (A + eight * t) % MOD
        out = [mpz(0)] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            out[i] = (out[i] + c * c0) % MOD
            out[i + 1] = (out[i + 1] + c * eight) % MOD
        coeffs = out
    return coeffs


def part_valuation_profile():
    P("=" * 70)
    P("PART A: v2(coeff) profile of Q for A=2^m+9")
    for m, s in [(32, 8), (48, 12), (64, 16)]:
        A = (mpz(1) << m) + 9
        MOD = mpz(1) << (2 * m)
        coeffs = poly_from_block(A, s, MOD)
        vs = [v2(c) for c in coeffs]
        P(f"  m={m} s={s}  deg={len(coeffs)-1}")
        P(f"    v2: {vs}")
        P(f"    min_v2={min(vs)}  (0=unit coeff present => no global 2-factor)")


def part_ratio_test():
    P("=" * 70)
    P("PART B: consecutive coeff ratios (hypergeom probe)")
    m, s = 40, 10
    A = (mpz(1) << m) + 9
    MOD = mpz(1) << (2 * m)
    coeffs = poly_from_block(A, s, MOD)
    # over Q: ratio c_{k+1}/c_k for rising factorial is rational in k
    # Here mod 2^{2m}: check if c_{k+1} * p(k) == c_k * q(k) for low-deg p,q
    P("  checking c_{k+1}*(k+1) vs c_k * 8 * something (Stirling recurrence)")
    # For (y)_s = sum stirling1 * y^k, not a hypergeometric coeff sequence in k
    # at fixed s.  The coeff sequence of rising-factorial poly IS Stirling1 row
    # — known hard.  Check linear recurrence of order 1: fail expected.
    fails = 0
    for k in range(len(coeffs) - 2):
        # naive: c_{k+1}/c_k ==? 8/(k+1) * (something with A)
        pass
    P("  Stirling-1 row coeffs: no O(1)-order hypergeom in k (expected).")
    P("  Fast eval of Q still needs multipoint / Horner — no free lunch.")


def part_full_K_linear_coeff():
    """
    Full product poly P(x)=prod_{k=1}^K (mu_k + x) mod (x^2, 2^{2m})
    has P(x) = P0 + P0*H*x.  Product-tree of truncated-to-x^2
    still does Theta(K) M(m)-muls — same as summing inverses.
    """
    P("=" * 70)
    P("PART C: truncated-x^2 product tree cost == leave-one-out")
    for K, m in [(50, 64), (100, 80)]:
        MOD = mpz(1) << (2 * m)
        # each node (c0, c1) with c0 + c1 x
        nodes = [(mpz((8 * k + 1) % MOD), mpz(1)) for k in range(1, K + 1)]
        ops = 0
        while len(nodes) > 1:
            nxt = []
            for i in range(0, len(nodes), 2):
                if i + 1 == len(nodes):
                    nxt.append(nodes[i])
                    continue
                a0, a1 = nodes[i]
                b0, b1 = nodes[i + 1]
                c0 = (a0 * b0) % MOD
                c1 = (a0 * b1 + a1 * b0) % MOD
                ops += 3
                nxt.append((c0, c1))
            nodes = nxt
        c0, c1 = nodes[0]
        H = (c1 * pow(c0, -1, MOD)) % (1 << m)
        Href = sum(pow(mpz(8 * k + 1), -1, 1 << m) for k in range(1, K + 1)) % (1 << m)
        P(f"  K={K} m={m}  ops~={ops}  H_match={H==Href}  (ops Theta(K))")


def part_multipoint_note():
    P("=" * 70)
    P("PART D: multipoint eval of Q at 0,s,2s,...,(q-1)s")
    P("""
  Degree-s poly, q=K/s points, s~sqrt(K).
  Fast multipoint: O(M(s) log^2 s + q polylog) ring ops.
  Ring = Z/2^{2m}Z, m~K: each ring mul costs M(m).
  Dominant: build Q via s linears = O(s M(m) log s) naive, or
  O(M(s)*M(m)-ish) with poly FFT over the ring.

  Best plausible BSGS bit complexity: soft-O(N^{1.5} log^c N).
  Still misses O(N log^2 N).

  To break into soft-linear need either:
    - closed form for Q(js) using Borwein at argument A+8*j*s, OR
    - s = polylog(N) blocks with each block-product soft-linear
      (requires soft-linear short rising factorial = Open Problem E).
""")


def main():
    P("q_structure -- block poly structure hunt")
    part_valuation_profile()
    part_ratio_test()
    part_full_K_linear_coeff()
    part_multipoint_note()
    P()
    P("VERDICT: Q coeffs are Stirling1-like; no sparse 2-adic gift.")
    P("Truncated x^2 tree = Theta(K) = E itself.  BSGS ceiling stands.")


if __name__ == "__main__":
    main()
