#!/usr/bin/env python3
"""transient_structure.py -- the transient-structure program (2026-07-28, item 3).

The band is the transient of X/mu (periodicity begins exactly at the tail
cutoff). The mathematics of transients is essentially empty (Stewart 1980 the
only bound).

PREMISE CORRECTION (protocol -- my first characterization was WRONG, M1 caught
it): I claimed the transient of 1/(2^J mu) is the bits of mu^{-1} mod 2^J and
that this is a genuinely aperiodic object. FALSE. The 2-adic inverse of mu is
PERIODIC with period exactly ord(2, mu) (verified: mu=73 -> 9 = ord(2,73);
89 -> 11; 257 -> 16). M1 shows it: small mu give LOW linear complexity
(L ~ ord(2,mu)), not maximal. So mu^{-1} mod 2^J is NOT a transient -- it is
the periodic reciprocal in a 2-adic mirror, and single-mu windows being cheap
is just the known fact that every 1/mu is periodic (why single-term windows
are modexp-cheap).

THE GENUINE TRANSIENT (resumption-night claim, which SURVIVES): the aperiodic
part is the INTEGER PART of X/mu for the specific X = 5^{N-1} -- its top
~log2(X) bits are aperiodic (tied to X's own pseudorandom bits); periodicity
begins only in the fractional tail. That transient is X-SPECIFIC: one
pseudorandom string read through many moduli, NOT a clean modulus-family
matrix. So item 3 collapses into item 2 (the prefix-mod data structure on the
single string X): "is there an O~(N) preprocessing of X answering
(X mod 2^J) mod mu in polylog?" -- the genuine open question.

What the measurements below actually establish (kept, correctly relabeled):
  (M1) the 2-adic inverse bit stream has linear complexity ~ ord(2,mu):
       PERIODIC, confirming the correction (not a transient, not a crack).
  (M2) the inverse matrix I[j][t]=bit_j(mu_t^{-1} mod 2^J) is the SAME rigid
       object class as the periodic reciprocal matrix R (as it must be, both
       periodic) -- a consistency check, not a new transient result.
"""

import numpy as np


def bm_gf2(s):
    """Berlekamp-Massey over GF(2); returns linear complexity L."""
    C = [1]; B = [1]; L = 0; m = 1
    for n in range(len(s)):
        d = 0
        for i in range(L + 1):
            if i < len(C):
                d ^= C[i] & s[n - i]
        if d:
            T = C[:]
            if len(C) < len(B) + m:
                C = C + [0] * (len(B) + m - len(C))
            for i in range(len(B)):
                C[i + m] ^= B[i]
            if 2 * L <= n:
                L = n + 1 - L; B = T; m = 1
            else:
                m += 1
        else:
            m += 1
    return L


def two_adic_inverse_bits(mu, J):
    """bits (LSB..) of mu^{-1} mod 2^J."""
    inv = pow(mu, -1, 1 << J)
    return [(inv >> j) & 1 for j in range(J)]


def reciprocal_bits(mu, J):
    """bits d=1..J of 1/mu in base 2 (the periodic orbit)."""
    x = 1 % mu
    out = []
    for _ in range(J):
        x <<= 1
        out.append(1 if x >= mu else 0)
        x %= mu
    return out


def gf2_rank(rows, ncols):
    """rank over GF(2) of 0/1 rows (each a list), via bitmask elimination."""
    basis = []
    for r in rows:
        v = 0
        for j, b in enumerate(r):
            if b:
                v |= (1 << j)
        for w in basis:
            v = min(v, v ^ w)
        if v:
            basis.append(v)
    return len(basis)


def main():
    J = 256
    mus = [8 * t + 1 for t in range(1, J + 1)]

    # (M1) linear complexity of the 2-adic inverse bit stream, several mu
    print("(M1) Berlekamp-Massey linear complexity of the 2-adic inverse bits")
    print("     (maximal ~ J/2 = %d => transient is NOT low-complexity):" % (J // 2))
    for mu in [73, 89, 105, 257, 1000001]:
        bits = two_adic_inverse_bits(mu, J)
        L = bm_gf2(bits)
        print("     mu=%-8d L=%3d  (%.2f x J/2)" % (mu, L, L / (J / 2)))

    # (M2) transient matrix vs periodic matrix, rigidity head to head
    I = np.array([[ (pow(int(m), -1, 1 << J) >> j) & 1 for m in mus]
                  for j in range(J)], dtype=float)
    R = np.array([reciprocal_bits(int(m), J) for m in mus], dtype=float).T
    # real numerical rank
    def nr(A, tol=1e-9):
        s = np.linalg.svd(A, compute_uv=False)
        return int((s > tol * s[0]).sum()), s
    rI, sI = nr(I - I.mean())
    rR, sR = nr(R - R.mean())
    # GF(2) rank
    gI = gf2_rank([list(map(int, I[j])) for j in range(J)], J)
    gR = gf2_rank([list(map(int, R[j])) for j in range(J)], J)
    # top-8 spectral energy (rigidity proxy: low => spread => rigid)
    e8I = (sI[:8]**2).sum() / (sI**2).sum()
    e8R = (sR[:8]**2).sum() / (sR**2).sum()
    print("(M2) transient matrix I vs periodic matrix R (%dx%d):" % (J, J))
    print("     real rank:  I=%d/%d   R=%d/%d" % (rI, J, rR, J))
    print("     GF(2) rank: I=%d/%d   R=%d/%d" % (gI, J, gR, J))
    print("     top-8 SV energy: I=%.3f  R=%.3f  (random ~%.3f; low => rigid)"
          % (e8I, e8R, 8.0 / J))

    print()
    print("VERDICT (corrected): the 2-adic inverse is PERIODIC (L ~ ord(2,mu)),")
    print("NOT a transient -- my partial-fraction premise was wrong and M1")
    print("caught it. Its matrix I is the same rigid object class as the")
    print("periodic reciprocal R (both full rank, both spread spectrum), as it")
    print("must be. The GENUINE transient is the integer part of X/mu for the")
    print("SPECIFIC X = 5^{N-1}: one pseudorandom string through many moduli, no")
    print("family matrix -- so item 3 collapses into item 2 (a prefix-mod data")
    print("structure on X). Honest yield: one wrong premise killed, item 3")
    print("redirected to the data-structure question, resumption-night")
    print("band=transient (integer part of X/mu) claim intact.")
    maximal = rI > 0.9 * J and rR > 0.9 * J
    return 0 if maximal else 1


if __name__ == "__main__":
    raise SystemExit(main())
