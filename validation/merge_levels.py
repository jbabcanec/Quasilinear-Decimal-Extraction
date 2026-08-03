#!/usr/bin/env python3
"""merge_levels.py -- ROUND 1 of the transform-reuse attack (2026-08-01).

merge_mass.py established the TOTALS (data ~ n log^2, work ~ n log^3).  This
splits both by LEVEL of the merge tree, to see which factor is which:

    work(level)  =?=  mass(level) x log2(operand size at that level)

If yes, the third log is exactly the FFT log inside M(b) = b log b, and the
second log is the number of tree levels, each of which re-carries the same
T log N bits.  That decomposition decides what an attack must remove:

    n log^3  =  (T log N bits)  x  (log T levels)  x  (log b per multiply)
                 ^ leaf data      ^ tree depth      ^ fast multiplication

and it says something sharper than "reuse transforms": with FREE (linear)
multiplication the merge would sit at n log^2 = EXACTLY AGM parity, never
better.  So the prize needs one of the first two factors, not the third.
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


def M(b1, b2):
    """quasi-linear multiply proxy, same as merge_mass.py"""
    s = b1 + b2
    return s * math.log2(s + 2)


def build_band(N, G=192):
    K = cutoff(N, G, 'pi')
    terms = bbp_pi_terms(N, K)
    LX = (mpz(5) ** (N - 1)).bit_length()
    band = [(sg, -E, mu) for sg, E, mu in terms
            if E < 0 and -E < LX + G + 8]
    band.sort(key=lambda t: t[1])
    return band


def level_profile(N):
    """Level-synchronous merge; returns per-level (nodes, mass, work, avg_bits)."""
    band = build_band(N)
    if not band:
        return None
    Jmax = max(J for _, J, _ in band)
    nodes = [(mpz(sg), mpz(mu), Jmax - J) for sg, J, mu in band]

    levels = []
    while len(nodes) > 1:
        mass = work = 0
        nxt = []
        for i in range(0, len(nodes) - 1, 2):
            (nL, dL, sL), (nR, dR, sR) = nodes[i], nodes[i + 1]
            s = sL if sL < sR else sR
            cap = Jmax - s
            bnL, bdL = nL.bit_length(), dL.bit_length()
            bnR, bdR = nR.bit_length(), dR.bit_length()
            mass += bnL + bdR + bnR + bdL + bdL + bdR
            work += M(bnL, bdR) + M(bnR, bdL) + M(bdL, bdR)
            n_ = ((nL * dR) << (sL - s)) + ((nR * dL) << (sR - s))
            d_ = dL * dR
            if n_.bit_length() > cap:
                n_ &= (mpz(1) << cap) - 1
            if d_.bit_length() > cap:
                d_ &= (mpz(1) << cap) - 1
            nxt.append((n_, d_, s))
        if len(nodes) & 1:
            nxt.append(nodes[-1])
        npairs = len(nodes) // 2
        levels.append((len(nodes), npairs, mass, work,
                       mass / (6.0 * npairs) if npairs else 0))
        nodes = nxt
    return Jmax, len(band), levels


def main():
    print("ROUND 1: is the third log the FFT log?   work(l) =?= mass(l)*log2(b_l)")
    for N in (10000, 30000, 100000):
        prof = level_profile(N)
        if prof is None:
            continue
        Jmax, T, levels = prof
        n = Jmax
        print()
        print("N=%d   band terms T=%d   Jmax=%d   levels=%d"
              % (N, T, Jmax, len(levels)))
        print("  %-4s %-8s %-12s %-14s %-10s %-10s"
              % ("lvl", "nodes", "mass(bits)", "work", "avg b", "work/mass"))
        tot_mass = tot_work = 0
        for l, (nn, npairs, mass, work, avgb) in enumerate(levels):
            tot_mass += mass
            tot_work += work
            if l < 6 or l >= len(levels) - 4:
                print("  %-4d %-8d %-12d %-14.3e %-10.0f %-10.2f"
                      % (l, nn, mass, work, avgb,
                         work / mass if mass else 0))
            elif l == 6:
                print("  ...")
        # the three factors
        TlogN = T * math.log2(8 * len(str(T)) + 2)   # leaf data proxy, see below
        leaf_bits = sum(mu.bit_length() for _, _, mu in build_band(N))
        print("  TOTAL mass=%.4e  work=%.4e   work/mass=%.2f"
              % (tot_mass, tot_work, tot_work / tot_mass))
        print("  leaf data (sum of denominator bits) = %d" % leaf_bits)
        print("  mass / leaf_data = %.2f   (should be ~ #levels = %d)"
              % (tot_mass / leaf_bits, len(levels)))
        print("  work / (mass*log2 n) = %.3f   (1.00 => third log IS the FFT log)"
              % (tot_work / (tot_mass * math.log2(n))))


if __name__ == '__main__':
    main()
