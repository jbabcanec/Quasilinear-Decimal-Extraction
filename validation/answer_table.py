#!/usr/bin/env python3
"""answer_table.py -- the X-specific horn of the representation dichotomy
(2026-08-03, second block).

THE DICHOTOMY (derivation, logged tonight): any bi-local representation
principle (the five-adjective spec from the invention-target sweep) yields
the answer table q_t = (X mod 2^{4t}) mod mu_t.  If it works for ARBITRARY
X it solves universal CRT-to-positional conversion (strictly stronger than
the prize: Door 2, middle bits of powers, Graham's number all fall).  If it
is X-SPECIFIC it must exploit structure of 5^N's table that random strings
lack, and producing the table is then equivalent to fast exact
incomplete-Kloosterman column counting for 5^N's phase family (the Bruin
watch item).  No third form exists.

THIS VALIDATOR measures the X-specific horn's foothold: is the answer table
of X = 5^N statistically distinguishable from the table of a random odd X
of the same length?  Batteries: parity balance + serial correlation,
Berlekamp-Massey linear complexity of the parity stream, zlib
compressibility of the packed table, uniformity of q_t/mu_t.  Any
separation = a wire trip (an exploitable 5^N fingerprint).  None =
the X-specific horn offers no visible purchase beyond the generic one.
"""

import random
import zlib

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def table(X):
    LX = int(X).bit_length()
    T = (LX - 8) // 4
    qs, xs = [], []
    for t in range(1, T + 1):
        mu = 8 * t + 1
        q = int((X & ((mpz(1) << (4 * t)) - 1)) % mu)
        qs.append(q)
        xs.append(q / mu)
    return qs, xs


def bm_gf2(bits):
    """Berlekamp-Massey linear complexity over GF(2)."""
    C, B, L, m = [1], [1], 0, 1
    for n in range(len(bits)):
        d = bits[n]
        for i in range(1, L + 1):
            if i < len(C) and C[i]:
                d ^= bits[n - i]
        if d == 0:
            m += 1
        elif 2 * L <= n:
            Tc = C[:]
            if len(C) < len(B) + m:
                C = C + [0] * (len(B) + m - len(C))
            for i, b in enumerate(B):
                if b:
                    C[i + m] ^= 1
            L, B, m = n + 1 - L, Tc, 1
        else:
            if len(C) < len(B) + m:
                C = C + [0] * (len(B) + m - len(C))
            for i, b in enumerate(B):
                if b:
                    C[i + m] ^= 1
            m += 1
    return L


def battery(name, qs, xs):
    T = len(qs)
    bits = [q & 1 for q in qs]
    bal = sum(bits) / T
    ser = sum(1 for i in range(T - 1) if bits[i] == bits[i + 1]) / (T - 1.0)
    Lc = bm_gf2(bits[:3000])
    n_bm = min(T, 3000)
    packed = b"".join(q.to_bytes(3, "little") for q in qs)
    comp = len(zlib.compress(packed, 9)) / float(len(packed))
    mean = sum(xs) / T
    # discrepancy proxy: max deviation of the empirical cdf on a 20-grid
    disc = max(abs(sum(1 for x in xs if x < g / 20.0) / T - g / 20.0)
               for g in range(1, 20))
    print("  %-12s T=%-6d parity %.4f  serial %.4f  BM %d/%d  "
          "zlib %.4f  mean %.4f  disc %.4f"
          % (name, T, bal, ser, Lc, n_bm, comp, mean, disc))
    return (bal, ser, Lc / (n_bm / 2.0), comp, mean, disc)


def main():
    print("the X-specific horn: answer table of 5^N vs random X (same length)")
    print("expect: 5^N indistinguishable => no foothold; separation => WIRE TRIP")
    rng = random.Random(17)
    for Np in (2000, 6000):
        X5 = mpz(5) ** Np
        LX = int(X5).bit_length()
        Xr = mpz(rng.getrandbits(LX - 2) << 1 | 1) | (mpz(1) << (LX - 1))
        print("N'=%d (%d bits):" % (Np, LX))
        r5 = battery("5^N", *table(X5))
        rr = battery("random", *table(Xr))
        gaps = [abs(a - b) for a, b in zip(r5, rr)]
        sep = (gaps[0] > 0.02 or gaps[1] > 0.02 or gaps[2] > 0.05
               or gaps[3] > 0.01 or gaps[5] > 0.02)
        print("  separation: %s"
              % ("*** WIRE TRIP -- investigate ***" if sep else
                 "none (within noise)"))
    print()
    print("VERDICT: no visible 5^N fingerprint in the answer table => the")
    print("X-specific horn offers no purchase beyond the generic one; the")
    print("representation principle, if built, is the universal converter")
    print("or the Kloosterman counter.  The dichotomy stands measured.")


if __name__ == '__main__':
    main()
