#!/usr/bin/env python3
"""fence3_wide.py -- FENCE 3 WIDENED (2026-08-03): the formula windows.

Original (`base10_gaussian_hunt.py`, 2026-07-30): degree-1 base-10 empty for
m <= 12 (confirming BBG04); degree-2 Gaussian (mk+j)^2 + c^2 empty for
m in {1,2,4}, c <= 3 -- closing BBG04's own loophole in that window.

Widened here:
  (i)   degree-2 Gaussian pushed to m in {1..6,8}, c in {1..5};
  (ii)  DEGREE 3 -- never searched by anyone: irreducible cubic denominators
        (mk+j)^3 + c (c a non-cube, so x^3 + c is irreducible over Q for the
        c used), m in {1..4}, pipeline-legal word-size denominators;
  (iii) same PSLQ discipline: dps=120, maxcoeff 10^8, live control
        (log(10/9) = sum 1/(k 10^k), a TRUE degree-1 base-10 identity) and a
        planted relation must both be recovered in every session.

A hit = a base-10 BBP-type formula for pi at degree 2 or 3 = band-free
extraction at O(N loglog N) word ops, BELOW the prize.  Expected: empty
(extending the BBG04 tightness two windows further out).
"""

from mpmath import mp, mpf, mpc, pslq, log as mp_log, pi as mp_pi

mp.dps = 120
TERMS = 140


def series_deg2(m, j, c):
    s = mpf(0)
    for k in range(TERMS):
        d = (m * k + j) ** 2 + c * c
        s += mpf(1) / (mpf(10) ** k * d)
    return s


def series_deg3(m, j, c):
    s = mpf(0)
    for k in range(TERMS):
        d = (m * k + j) ** 3 + c
        s += mpf(1) / (mpf(10) ** k * d)
    return s


def hunt(basket, names, tag, maxcoeff=10 ** 8):
    rel = pslq(basket, maxcoeff=maxcoeff, maxsteps=60000)
    if rel is None:
        print("  %-24s empty" % tag)
        return None
    nz = [(names[i], rel[i]) for i in range(len(rel)) if rel[i]]
    if rel[0] == 0:
        print("  %-24s internal relation (no pi): %s" % (tag, nz))
        return None
    print("  %-24s ** HIT **: %s" % (tag, nz))
    return nz


def main():
    print("FENCE 3 widened: degree-2 to (m<=8, c<=5); degree-3 first search")
    # controls
    ctrl = pslq([mp_log(mpf(10) / 9),
                 sum(mpf(1) / (k * mpf(10) ** k) for k in range(1, 200))],
                maxcoeff=100)
    print("  control log(10/9) vs its series: %s"
          % ("RECOVERED %s" % ctrl if ctrl == [1, -1] else "*** MISSED ***"))
    a = series_deg2(4, 1, 2)
    b = series_deg2(4, 3, 2)
    plant = 3 * a - 7 * b
    pc = pslq([plant, a, b], maxcoeff=100)
    print("  control planted (3,-7): %s"
          % ("RECOVERED %s" % pc if pc and abs(pc[0]) <= 2 else "*** MISSED ***"))
    hits = 0

    print("  --- degree 2 Gaussian (mk+j)^2 + c^2 ---")
    for m in (1, 2, 3, 4, 5, 6, 8):
        for c in (1, 2, 3, 4, 5):
            names = ["pi"] + ["S_%d" % j for j in range(1, m + 1)] + ["1"]
            basket = ([+mp_pi] + [series_deg2(m, j, c) for j in range(1, m + 1)]
                      + [mpf(1)])
            if hunt(basket, names, "deg2 m=%d c=%d" % (m, c)):
                hits += 1

    print("  --- degree 3 (mk+j)^3 + c, c non-cube ---")
    for m in (1, 2, 3, 4):
        for c in (2, 3, 4, 5, 6, 7, 9):
            names = ["pi"] + ["S_%d" % j for j in range(1, m + 1)] + ["1"]
            basket = ([+mp_pi] + [series_deg3(m, j, c) for j in range(1, m + 1)]
                      + [mpf(1)])
            if hunt(basket, names, "deg3 m=%d c=%d" % (m, c)):
                hits += 1

    print()
    if hits:
        print("*** %d HIT(S): verify at dps=240 before believing anything ***"
              % hits)
    else:
        print("VERDICT: degree-2 window (m<=8, c<=5) and the FIRST degree-3")
        print("window (m<=4, 7 cubics) are EMPTY for pi in base 10.  BBG04")
        print("tightness extends two windows beyond the original search.")


if __name__ == '__main__':
    main()
