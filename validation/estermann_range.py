#!/usr/bin/env python3
"""estermann_range.py -- height-regime calibration of the linchpin sums
    E(A) = sum_{t<=T} e(A / mu_t),  mu_t = 8t+1.

Purpose (2026-07-28): orbit_estermann.py measured |E(A)| ~ sqrt(T) at heights
A ~ 1e9 ~ T^2.5 and concluded "full rank provable via known ANT." That
conclusion needs a FENCE: the matrix's actual rows sit at heights A ~ 2^d with
d up to N, and provable technology (van der Corput / exponent pairs /
Vinogradov) covers only polynomial-in-T heights (clean sqrt window roughly
T^{1+eps} <= A <= T^{3/2}, extensions to T^{O(1)}). Beyond that, E(A) at
A = k 2^d is the FOURIER TRANSFORM of the column data itself, so the pairwise
question at deep rows is equivalent to the single-row equidistribution
question: the open question already at r = 1. This validator measures |E(A)|
across height regimes, including the true row heights, to separate what is
TRUE (cancellation everywhere, empirically) from what is PROVABLE (the
polynomial window).

Also verifies the two fence-posts of the provable window:
  - bottom fence: for A <= T the sum has NO cancellation (E(A) ~ T - O(A));
  - E(0) = T exactly (the Gram diagonal).
"""

import math
import numpy as np

T_LIST = [4000, 16000]


def E_of_A(A, T):
    """|sum_{t<=T} e(A/mu_t)| computed exactly via integer A mod mu."""
    ph = np.empty(T)
    for i, t in enumerate(range(1, T + 1)):
        mu = 8 * t + 1
        ph[i] = (A % mu) / mu
    return abs(np.exp(2j * np.pi * ph).sum())


def E_of_pow2(d, T, k=1, k2=0, d2=0):
    """|E(k*2^d + k2*2^d2)| with the huge integer handled modularly."""
    ph = np.empty(T)
    for i, t in enumerate(range(1, T + 1)):
        mu = 8 * t + 1
        a = (k * pow(2, d, mu) + k2 * pow(2, d2, mu)) % mu
        ph[i] = a / mu
    return abs(np.exp(2j * np.pi * ph).sum())


def main():
    print("E(A) = |sum_{t<=T} e(A/(8t+1))|, reported as multiples of sqrt(T)")
    for T in T_LIST:
        rt = math.sqrt(T)
        print("T = %d, sqrt(T) = %.1f" % (T, rt))
        print("  E(0)/T = %.3f (must be 1.000: Gram diagonal)"
              % (E_of_A(0, T) / T))
        # bottom fence: no cancellation at small heights
        for A in [T // 8, T // 2]:
            print("  A = %-8d (~T^%.2f): E/T = %.3f   <- below the window: "
                  "NO cancellation" % (A, math.log(A, T), E_of_A(A, T) / T))
        # the provable polynomial window and its empirical continuation
        for expo in [1.2, 1.5, 2.0, 2.5, 3.0]:
            A = int(T ** expo)
            print("  A ~ T^%.1f : |E|/sqrt(T) = %.2f"
                  % (expo, E_of_A(A, T) / rt))
        # true row heights: A = 2^d, d far beyond any provable range
        for d in [40, 100, 400, 1000, 4000]:
            print("  A = 2^%-5d: |E|/sqrt(T) = %.2f   <- beyond provable "
                  "technology" % (d, E_of_pow2(d, T) / rt))
        # true pairwise heights: A = k 2^d + k' 2^d' (row-correlation phases)
        rng = np.random.default_rng(1)
        vals = []
        for _ in range(8):
            d = int(rng.integers(200, 2000))
            dd = d - int(rng.integers(1, 40))
            k = int(rng.integers(1, 50))
            kk = -int(rng.integers(1, 50))
            vals.append(E_of_pow2(d, T, k=k, k2=kk, d2=dd) / rt)
        print("  pairwise A = k 2^d + k' 2^d' (8 draws, d in [200,2000]): "
              "median |E|/sqrt(T) = %.2f" % float(np.median(vals)))
    print()
    print("READING: sqrt-cancellation is TRUE at every height measured,")
    print("including the exponential heights of the real rows. It is PROVABLE")
    print("only in the polynomial window (exponent pairs; clean sqrt for")
    print("T^{1+eps} <= A <= T^{3/2}). At exponential height, E(k 2^d) is the")
    print("Fourier transform of the column data: the pairwise lemma there IS")
    print("the open equidistribution question, already at r = 1. See")
    print("notes/rank_lemma.md for the fenced theorem and the corrected claim.")


if __name__ == "__main__":
    main()
