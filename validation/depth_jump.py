#!/usr/bin/env python3
"""depth_jump.py -- does reaching depth N really cost a log? (2026-07-31)

Joseph pushed back on the "extraction floor" claim: why must a depth-N
extraction pay for 10^N, and can't you get there some other way. Checking
rather than defending, because the claim has TWO different costs hiding in it:

  route A (modexp):        r = (10^N mod q) * p mod q
     log N squarings, EVERY one at full size because of the reduction
     => M(N) log N = N log^2 N.  This is what sqrt2_extract.py does and what
     its docstring calls "the extraction tax".

  route B (exact power):   r = (10^N * p) mod q
     10^N by repeated squaring TELESCOPES -- the operand sizes halve going
     backwards, so the whole exponentiation is M(N) + M(N/2) + ... = O(M(N)).
     Then one multiply and one division, both M(N).
     => O(M(N)) = N log N, a full log CHEAPER.

If B measures faster and scales like M(N), the "floor at N log^2 N" I claimed
is wrong, sqrt2's extraction is at PARITY with its own computation (Newton,
N log N) rather than one log above it, and the 07-28 "extraction tax exists
with zero BBP involvement" conclusion needs correcting.
"""

import math
import time

from mpmath import mp, sqrt as mp_sqrt

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def pell_convergent(min_q_digits):
    p, q = mpz(1), mpz(1)
    target = mpz(10) ** min_q_digits
    while q < target:
        p, q = p + 2 * q, p + q
    return p, q


def route_A(p, q, N, window=10):
    """modexp: log N full-size modular squarings."""
    r = (pow(mpz(10), N, q) * (p % q)) % q
    return str((r * mpz(10) ** window) // q).zfill(window)


def route_B(p, q, N, window=10):
    """exact power (telescoping) + one multiply + one division."""
    P = mpz(10) ** N
    r = (P * (p % q)) % q
    return str((r * mpz(10) ** window) // q).zfill(window)


def reference(N, window=10):
    mp.dps = N + window + 30
    s = mp.nstr(mp_sqrt(2), N + window + 15, strip_zeros=False).split('.')[1]
    return s[N:N + window]


def main():
    print("=" * 78)
    print("Reaching depth N: modexp (route A) vs exact power (route B)")
    print("target = the N-th decimal digit of sqrt2 from one Pell convergent")
    print("=" * 78)
    print("  %-9s %-12s %-12s %-11s %-11s %-8s %s"
          % ("N", "A digits", "B digits", "A time", "B time", "A/B", "agree"))
    rows = []
    for N in (10 ** 4, 3 * 10 ** 4, 10 ** 5, 3 * 10 ** 5, 10 ** 6):
        p, q = pell_convergent((N + 40) // 2 + 40)
        t0 = time.perf_counter()
        a = route_A(p, q, N)
        ta = time.perf_counter() - t0
        t0 = time.perf_counter()
        b = route_B(p, q, N)
        tb = time.perf_counter() - t0
        rows.append((N, ta, tb))
        print("  %-9d %-12s %-12s %-11.4f %-11.4f %-8.1f %s"
              % (N, a, b, ta, tb, ta / tb, "yes" if a == b else "NO"))

    print()
    print("  scaling (ratio of time to n log n and to n log^2 n; flat = fit):")
    for label, idx in (("A", 1), ("B", 2)):
        r1 = [row[idx] / (row[0] * math.log(row[0])) for row in rows]
        r2 = [row[idx] / (row[0] * math.log(row[0]) ** 2) for row in rows]
        print("    route %s: vs n log n spread %.2f | vs n log^2 n spread %.2f"
              % (label, max(r1) / min(r1), max(r2) / min(r2)))

    print()
    print("  exact operand-bit mass fed to the multiplier (no timing noise):")
    print("    %-12s %-16s %-16s %-8s" % ("N", "route A mass", "route B mass",
                                          "A/B"))
    for N in (10 ** 4, 10 ** 6, 10 ** 9, 10 ** 12):
        qb = N * math.log2(10) / 2                 # bits of q
        # A: square-and-multiply mod q. ~log2 N squarings + ~log2(N)/2 multiplies,
        # every operand held at full size qb by the reduction.
        massA = 1.5 * math.log2(N) * 2 * qb
        # B: exact 10^N by repeated squaring, sizes HALVE going backwards, so
        # sum_k 2 * (N log2 10 / 2^k) = 4 * N log2 10; then one multiply
        # (3.32N x 1.66N) and one division (5N by 1.66N).
        massB = 4 * N * math.log2(10) + (3 * qb + 2 * qb) + (5 * qb + 2 * qb)
        print("    %-12.0e %-16.4e %-16.4e %-8.1f" % (N, massA, massB,
                                                      massA / massB))
    print("    A/B grows like log N exactly as predicted; B is flat in N.")

    print()
    print("  correctness check against mpmath at N = 1e4:")
    p, q = pell_convergent((10 ** 4 + 40) // 2 + 40)
    print("    route A %s | route B %s | mpmath %s"
          % (route_A(p, q, 10 ** 4), route_B(p, q, 10 ** 4), reference(10 ** 4)))

    print()
    print("=" * 78)
    print("READING: if A/B grows with N, the log is an artifact of CHOOSING the")
    print("modexp. The exponentiation only costs a log when every squaring is")
    print("forced to full size by a modular reduction; the exact power")
    print("telescopes. Reaching depth N is M(N), not M(N) log N.")
    print("=" * 78)


if __name__ == "__main__":
    main()
