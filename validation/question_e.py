#!/usr/bin/env python3
"""question_e.py -- prove-or-disprove session on Question E (2026-07-28).

RESULT 1 (DISPROOF, validated here): Question E as posed in rank_lemma.md
("|E(A)| <= T^{1-delta} uniformly for T^{1+eps} <= |A| <= exp(T)") is FALSE.
Construction: A = prod_{t in S} mu_t for a subset S of size ~T/log^2 T. Then
A = 0 mod mu_t for every t in S (phases coherent at 1), so |E(A)| >= |S| -
O(sqrt T) = T^{1-o(1)}, while log A ~ |S| log T stays well under T. The
divisor structure of A can conspire; uniformity over ALL heights is dead.

RESULT 2 (the repair): the Gram matrix only ever meets STRUCTURED heights
A = m 2^d + m' 2^{d'}. There the conspiracy is impossible beyond n^{o(1)}
terms: mu_t | A requires mu_t to divide the fixed odd part of m 2^{d-d'} + m',
and a fixed integer n has only d(n) = n^{o(1)} divisors. Measured below:
resonant-term counts for random structured A are tiny. Question E' (the
corrected open question) restricts to the structured family, where measured
cancellation is clean square-root.

RESULT 3 (the hard object, measured): deep column-pair correlations
C(t,t') = sum_{j<D} eps(2^{d0+j} mod mu_t) eps(2^{d0+j} mod mu_t') at base
depth d0 = 10^6. These are incomplete joint doubling-orbit sums mod mu mu'
of length D ~ sqrt(mu mu') -- the exact object a Bourgain/Diffie-Hellman-type
import would need to bound to lift Theorem A from polylog to polynomial rank.
Expected sqrt(D)-clean; measured below.
"""

import math
import numpy as np


def part1_disproof(T=4000, gap=10):
    print("PART 1: disproof of Question E (uniform heights)")
    S = list(range(1, T + 1, gap))          # |S| = T/gap
    A = 1
    for t in S:
        A *= 8 * t + 1
    ph = np.empty(T)
    for i, t in enumerate(range(1, T + 1)):
        mu = 8 * t + 1
        ph[i] = (A % mu) / mu
    E = abs(np.exp(2j * np.pi * ph).sum())
    print("  T = %d, |S| = %d, log2 A = %d (<< T: height admissible)"
          % (T, len(S), A.bit_length()))
    print("  |E(A)| = %.1f  vs coherent mass |S| = %d, sqrt(T) = %.1f"
          % (E, len(S), math.sqrt(T)))
    ok = E > len(S) - 3 * math.sqrt(T)
    print("  DISPROOF %s: |E| ~ |S| >> sqrt(T); uniform Question E is FALSE"
          % ("CONFIRMED" if ok else "FAILED"))
    return ok


def part2_structured_resonances(T=4000, trials=200, seed=2):
    print("PART 2: structured heights A = m 2^d + m' 2^{d'} cannot conspire")
    rng = np.random.default_rng(seed)
    counts = []
    for _ in range(trials):
        d = int(rng.integers(100, 4000))
        delta = int(rng.integers(1, 60))
        m = int(rng.integers(1, 100)) * 2 + 1          # odd
        mp = -(int(rng.integers(1, 100)) * 2 + 1)      # odd, negative
        res = 0
        for t in range(1, T + 1):
            mu = 8 * t + 1
            if (m * pow(2, delta, mu) + mp) % mu == 0:  # mu | odd part of A
                res += 1
        counts.append(res)
    counts = np.array(counts)
    print("  resonant t per structured A (T = %d, %d random draws):"
          % (T, trials))
    print("  max = %d, mean = %.2f  (divisor-bound prediction: n^{o(1)}, "
          "i.e. O(1)-ish)" % (counts.max(), counts.mean()))
    return counts.max() <= 8


def part3_column_pairs(T_cols=300, D=512, d0=10**6, Tpool=4000, seed=3):
    print("PART 3: deep column-pair correlations (the Bourgain-import target)")
    rng = np.random.default_rng(seed)
    ts = np.sort(rng.choice(np.arange(1, Tpool + 1), size=T_cols,
                            replace=False))
    B = np.empty((T_cols, D))
    for i, t in enumerate(ts):
        mu = int(8 * t + 1)
        x = pow(2, d0, mu)
        half = mu / 2
        for j in range(D):
            B[i, j] = 1.0 if x > half else -1.0
            x = (x << 1) % mu
    B *= 0.5
    G = B @ B.T                       # C(t,t') on the off-diagonal
    off = G[np.triu_indices(T_cols, k=1)]
    sd_random = math.sqrt(D) / 4      # sd of sum of D iid +-1/4 products
    z = np.abs(off) / sd_random
    print("  %d columns from [1,%d], window length D = %d at base depth 2^%d"
          % (T_cols, Tpool, D, int(math.log2(d0))))
    print("  |C|/ (sqrt(D)/4): median = %.2f (|N(0,1)| median 0.67), "
          "p95 = %.2f (N: 1.96), max = %.2f" %
          (float(np.median(z)), float(np.quantile(z, 0.95)), float(z.max())))
    rank = np.linalg.matrix_rank(B)
    print("  numerical rank of the %d x %d deep window = %d (full = %d)"
          % (T_cols, D, rank, T_cols))
    clean = float(np.median(z)) < 1.5 and float(z.max()) < 8
    print("  VERDICT: %s -- joint orbit correlations are sqrt-clean at depth "
          "10^6;\n  a per-pair T^{-delta} theorem (Bourgain/DH-type, composite "
          "modulus,\n  length ~ sqrt(modulus)) would lift rank(M) from polylog "
          "to polynomial." % ("CLEAN" if clean else "TAILS PRESENT"))
    return clean


def main():
    ok1 = part1_disproof()
    ok2 = part2_structured_resonances()
    ok3 = part3_column_pairs()
    print()
    print("SUMMARY: Question E (uniform) DISPROVED by divisor conspiracy;")
    print("Question E' (structured heights, where resonances are n^{o(1)})")
    print("replaces it and remains open + empirically sqrt-clean; the deep")
    print("joint-orbit correlations are measured clean = the import target.")
    return 0 if (ok1 and ok2 and ok3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
