#!/usr/bin/env python3
"""merge_reuse.py -- ROUND 2 of the transform-reuse attack (2026-08-01).

ROUND 1 (`merge_levels.py`) showed work(l) = mass(l) * log2(b_l) with b_l
doubling per level, so the third log is the GROWTH of the operands, not a
re-transform of reusable data.  Across levels there is nothing to reuse: every
operand is a fresh product.  That leaves exactly two recoverable quantities,
both CONSTANTS, both measured here:

  (R1) WITHIN-MERGE TRANSFORM REUSE.  A merge does three multiplications --
       nL*dR, nR*dL, dL*dR -- over only FOUR distinct operands: dL and dR each
       appear twice.  Naive: 6 forward + 3 inverse = 9 transforms.  Cached:
       4 forward + 3 inverse = 7.

  (R2) WRAPPED (mod 2^cap) MULTIPLICATION.  Every node result is truncated to
       cap = Jmax - shift bits, so a full 2b-bit product is computed and half
       discarded.  A cyclic/negacyclic convolution of length cap computes the
       low cap bits directly (Schoenhage-Strassen multiplies mod 2^k+1
       natively), halving the transform length wherever the cap binds.

Both are constants, so neither touches the exponent.  They are reported because
the 10^12 record run is scheduled and these are implementable there.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decimal_bbp_extract import bbp_pi_terms, cutoff

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def M(b):
    return b * math.log2(b + 2)


def profile(N, G=192):
    K = cutoff(N, G, 'pi')
    terms = bbp_pi_terms(N, K)
    LX = (mpz(5) ** (N - 1)).bit_length()
    band = [(sg, -E, mu) for sg, E, mu in terms if E < 0 and -E < LX + G + 8]
    band.sort(key=lambda t: t[1])
    if not band:
        return None
    Jmax = max(J for _, J, _ in band)
    nodes = [(mpz(sg), mpz(mu), Jmax - J) for sg, J, mu in band]

    st = {'naive': 0.0, 'cached': 0.0, 'wrapped': 0.0, 'both': 0.0,
          'prod_bits': 0, 'kept_bits': 0, 'capped_work': 0.0, 'work': 0.0}

    while len(nodes) > 1:
        nxt = []
        for i in range(0, len(nodes) - 1, 2):
            (nL, dL, sL), (nR, dR, sR) = nodes[i], nodes[i + 1]
            s = sL if sL < sR else sR
            cap = Jmax - s
            bnL, bdL = nL.bit_length(), dL.bit_length()
            bnR, bdR = nR.bit_length(), dR.bit_length()

            # transform-count model: a b1 x b2 multiply needs transforms of
            # length ~ (b1 + b2); charge M(len) per transform.
            for (x, y) in ((bnL, bdR), (bnR, bdL), (bdL, bdR)):
                L = x + y
                st['naive'] += 3 * M(L)                # 2 forward + 1 inverse
                Lw = min(L, cap) if cap < L else L     # wrapped length
                st['wrapped'] += 3 * M(Lw)
                st['prod_bits'] += L
                st['kept_bits'] += min(L, cap)
                st['work'] += M(L)
                if L > cap:
                    st['capped_work'] += M(L)
            # cached: dL and dR are each transformed once instead of twice.
            # 9 transforms -> 7, at the sizes actually used.
            LnLdR, LnRdL, LdLdR = bnL + bdR, bnR + bdL, bdL + bdR
            st['cached'] += (3 * M(LnLdR) + 3 * M(LnRdL) + 3 * M(LdLdR)
                             - M(LdLdR) - M(LdLdR))    # save 2 forwards
            st['both'] += ((3 * M(min(LnLdR, cap)) + 3 * M(min(LnRdL, cap))
                            + 3 * M(min(LdLdR, cap)))
                           - 2 * M(min(LdLdR, cap)))

            n_ = ((nL * dR) << (sL - s)) + ((nR * dL) << (sR - s))
            d_ = dL * dR
            if n_.bit_length() > cap:
                n_ &= (mpz(1) << cap) - 1
            if d_.bit_length() > cap:
                d_ &= (mpz(1) << cap) - 1
            nxt.append((n_, d_, s))
        if len(nodes) & 1:
            nxt.append(nodes[-1])
        nodes = nxt
    return Jmax, len(band), st


def main():
    print("ROUND 2: what is actually recoverable?  (both are CONSTANTS)")
    print("%-9s %-8s %-11s %-11s %-11s %-12s"
          % ("N", "T", "R1 cache", "R2 wrap", "R1+R2", "discarded"))
    for N in (10000, 30000, 100000):
        r = profile(N)
        if r is None:
            continue
        Jmax, T, st = r
        r1 = 1 - st['cached'] / st['naive']
        r2 = 1 - st['wrapped'] / st['naive']
        both = 1 - st['both'] / st['naive']
        disc = 1 - st['kept_bits'] / float(st['prod_bits'])
        print("%-9d %-8d %-11s %-11s %-11s %-12s"
              % (N, T, "%.1f%%" % (100 * r1), "%.1f%%" % (100 * r2),
                 "%.1f%%" % (100 * both), "%.1f%% of product bits"
                 % (100 * disc)))
    print()
    print("R1 (within-merge transform caching): 9 transforms -> 7 per merge.")
    print("R2 (wrapped mod-2^cap multiply): only where the cap binds, i.e. the")
    print("   top of the tree, which is also where the operands are largest.")
    print("NEITHER TOUCHES THE EXPONENT -- work stays Theta(n log^3 n).")


if __name__ == '__main__':
    main()
