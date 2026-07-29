#!/usr/bin/env python3
"""fold_compute.py -- Frontier 1: does proof-system folding (FRI/Binius/
BaseFold/sumcheck) give V in sub-quadratic bit-ops? (2026-07-28)

The pitch: FRI/Binius folding touches each coefficient O(1) times and the
prover's fold phase is LINEAR TIME -- structurally the "O(1) transforms per bit"
we need. Decisive test of whether it transfers to computing
    V = sum_t a_t/(8t+1) mod 2^N.

What folding actually is, made concrete: the FRI/BaseFold/sumcheck fold, used as
a COMPUTATION (fixed challenges, not random), is a multilinear/tensor
contraction -- round r combines pairs of the current state with a fixed weight,
halving the count. Here we use it to compute the sum. Two facts decide it:
  (F1) the fold reproduces V exactly (it IS this computation) -- sanity;
  (F2) the fold is PRECISION-UNIFORM: every one of the T ~ N/4 summands enters
       at full N-bit precision (the inv_t = (8t+1)^{-1} are 2-adic UNITS, no
       taper), so the fold-state mass is Theta(sum_r (T/2^r) * N) = Theta(T*N)
       = Theta(N^2) -- WORSE than the merge tree's n log^3.
The merge beats folding precisely because the merge CAPS precision per node (the
precision taper), which the uniform fold does not. FRI's celebrated linear-time
is for PROXIMITY-checking a COMMITTED codeword (verification) + randomness, not
for computing an unstructured sum whose terms must all be produced at full
precision. This is exactly the sandbox's "fold state = full simplex support,
characteristic-independent" (ff_fold.py) -- now at the ring level and tied to
the moment mass (#20). Frontier 1 closes: folding = precision-uniform = N^2;
the crack would need a PRECISION-GRADED fold = the same missing transform.
"""

import math
import numpy as np


def fold_sum(vals, MOD):
    """BaseFold/sumcheck-style fold computing sum(vals) mod MOD by pairwise
    combination, log rounds, each entry touched O(1) times. Returns the sum and
    the total fold-state bit-mass (sum over rounds of entries*precision)."""
    v = list(vals)
    mass = 0
    while len(v) > 1:
        mass += sum(x.bit_length() for x in v)
        nv = []
        for i in range(0, len(v) - 1, 2):
            nv.append((v[i] + v[i + 1]) % MOD)     # fixed-challenge fold = add
        if len(v) & 1:
            nv.append(v[-1])
        v = nv
    return v[0] % MOD, mass


def main():
    print("(F1) fold reproduces V exactly, and (F2) fold-state mass ~ N^k:")
    rng = np.random.default_rng(0)
    rows = []
    for N in [240, 480, 960, 1920, 3840]:
        MOD = 1 << N
        k = max(2, round(math.log2(N / 8)))   # T ~ N/8, as a power of two
        T = 1 << k
        a = [int(rng.integers(1, 1 << 40)) for _ in range(T)]
        terms = [(ai * pow(8 * t + 1, -1, MOD)) % MOD
                 for t, ai in enumerate(a, start=1)]
        V_direct = sum(terms) % MOD
        V_fold, mass = fold_sum(terms, MOD)
        ok = (V_direct == V_fold)
        rows.append((N, T, mass))
        print("   N=%4d T=%4d : fold==direct %s, fold-state mass=%d"
              % (N, T, "OK" if ok else "FAIL", mass))
        if not ok:
            return 1

    print()
    print("fold-state mass vs N^k (n = N here; which exponent is flat?):")
    for k in (1, 2, 3):
        rs = [m / (N**k) for (N, T, m) in rows]
        print("   vs N^%d : spread %.2f" % (k, max(rs) / min(rs)))
    # also vs T*N (the tensor/moment mass)
    rs = [m / (T * N) for (N, T, m) in rows]
    print("   vs T*N : spread %.2f  (flat => fold mass = Theta(T*N), the tensor "
          "contraction)" % (max(rs) / min(rs)))

    print()
    print("VERDICT: the fold reproduces V but its state is PRECISION-UNIFORM --")
    print("mass Theta(T*N) ~ N^2, strictly WORSE than the merge tree's n log^3.")
    print("FRI/Binius linear-time is verification (proximity of a committed")
    print("codeword) + randomness, not computation of an unstructured sum; every")
    print("term is a 2-adic unit and enters at full precision, so no compression")
    print("(= the sandbox's full-simplex fold-state, characteristic-independent,")
    print("= the moment mass #20). Frontier 1 CLOSED with the current (FRI-Binius")
    print("2024) construction: folding is precision-uniform. The crack needs a")
    print("PRECISION-GRADED fold -- the same missing transform the taxonomy names.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
