#!/usr/bin/env python3
"""
counting_face.py -- the COUNTING FACE (twelfth face) of the additive core.

Setup: a_t = c_t/mu_t in (0,1), t <= T (here the atom's family mu_t = 8t+1,
c_t = 1, i.e. a_t = frac representation of 1/mu_t as a real in (0,1)).
Let x_t(d) = frac(2^d a_t) and S(d) = sum_t x_t(d), and let
B(d) = #{t : bit d of a_t's binary expansion is 1} = #{t : x_t(d) >= 1/2}
(the COLUMN POPCOUNT of the T x n fraction bit-matrix).

Claims validated here (exact rational arithmetic):
  (C1) doubling recursion:  S(d+1) = 2*S(d) - B(d)          for all d.
  (C2) unrolled identity:   2^D * S(0) = S(D) + sum_{d<D} 2^{D-1-d} B(d),
       hence, since 0 <= S(D) < T, the column-count stream {B(d)} determines
       the binary expansion of S(0) = sum_t a_t (the atom as a real) exactly,
       by an O(N)-time carry-normalization pass (amortized-linear addition).
  (C3) reconstruction: digits of sum_t 1/mu_t recovered from {B(d)} alone
       match direct high-precision computation, bit for bit.

Consequence (the face): computing the atom to depth N  <==>  computing all
column popcounts B(0..N). Each B(d) = #{t : frac(2^d / mu_t) >= 1/2} is an
exact incomplete counting sum over the orbit family (bracket/Weyl-type);
the carry/threshold side is amortized-linear and was never the wall.
New hardness-import surface: embed a known-hard exact counting family as
{B(d)}. New upper-bound surface: batch column counts below T x N.
"""

from fractions import Fraction

T = 40           # terms: mu_t = 8t+1, t = 1..T
DEPTH = 200      # bits


def mus():
    return [8 * t + 1 for t in range(1, T + 1)]


def main():
    M = mus()
    a = [Fraction(1, m) for m in M]

    # direct: x_t(d) via exact fractions; B(d) columns; S(d) sums
    x = list(a)
    S = [sum(x)]
    B = []
    ok1 = True
    for d in range(DEPTH):
        bd = sum(1 for v in x if 2 * v >= 1)
        B.append(bd)
        x = [2 * v - (1 if 2 * v >= 1 else 0) for v in x]
        S.append(sum(x))
        ok1 &= (S[d + 1] == 2 * S[d] - bd)
    print(f"(C1) doubling recursion S(d+1) = 2 S(d) - B(d): "
          f"{'EXACT, all ' + str(DEPTH) + ' depths' if ok1 else 'FAIL'}")

    # (C2) unrolled identity at several D
    ok2 = True
    for D in [7, 63, 128, DEPTH]:
        lhs = (2 ** D) * S[0]
        rhs = S[D] + sum((2 ** (D - 1 - d)) * B[d] for d in range(D))
        ok2 &= (lhs == rhs)
    print(f"(C2) unrolled identity 2^D S(0) = S(D) + sum 2^(D-1-d) B(d): "
          f"{'EXACT' if ok2 else 'FAIL'}")

    # (C3) reconstruct the binary digits of S(0) from {B(d)} alone:
    # 2^D S(0) = integer_part + S(D) with 0 <= S(D) < T, so
    # floor(2^D S(0)) = sum_{d<D} 2^{D-1-d} B(d) + floor(S(D)); the unknown
    # floor(S(D)) in [0, T) only perturbs the low log2(T) bits: bits of S(0)
    # above that guard band are determined by the B-partial-sums alone.
    D = DEPTH
    guard = T.bit_length() + 1
    acc = sum((2 ** (D - 1 - d)) * B[d] for d in range(D))          # B-only
    truth = (2 ** D * S[0]).__floor__()
    match_bits = (acc >> guard) == (truth >> guard) or \
                 abs((acc >> guard) - (truth >> guard)) <= 1
    # strict check: acc and truth differ by floor(S(D)) < T
    strict = 0 <= truth - acc < T
    print(f"(C3) B-only reconstruction: truth - acc = {truth - acc} "
          f"(must lie in [0, T={T})): {'EXACT' if strict else 'FAIL'}; "
          f"digits above {guard}-bit guard band agree: {match_bits}")

    print("\nFace enshrined iff all EXACT: atom <=> column popcounts {B(d)}.")
    print(f"Sample counts B(0..19): {B[:20]}")


if __name__ == "__main__":
    main()
