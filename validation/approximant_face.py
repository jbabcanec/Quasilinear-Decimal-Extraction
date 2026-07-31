#!/usr/bin/env python3
"""approximant_face.py -- is there a FUNDAMENTALLY different representation?
(2026-07-31)

Joseph: "what about a fundamentally different or new way of representing pi."

`sqrt2_extract.py` (2026-07-28) is the clue nobody followed up. sqrt2 has NO
series, NO BBP, NO merge tree: ONE Pell convergent p/q with q ~ 10^(N/2), then
ONE modexp, and the digit falls out at N log^2 N. So the representation that
beats a series is not a better series -- it is a single RATIONAL APPROXIMANT,
provided p and q are computable without computing the constant's digits.

That gives a universal floor and a test:

  FLOOR.  Any depth-N decimal extraction must apply 10^N to something, which
  is log N squarings at M(N) each = N log^2 N. Nothing can go below it (short
  of beating addition chains). So N log^2 N is not just "AGM parity", it is
  the extraction floor.

  TEST.   A constant attains the floor iff it has a rational approximant
  p/q -> c with |c - p/q| < 10^-N, q of O(N) bits, and p, q computable in
  N log^2 N WITHOUT computing c's digits. Call the ratio

      redundancy = (bits of approximant data) / (bits of answer)

  A constant with redundancy 1 attains the floor. Redundancy log N costs a log.

This script:
  (1) BUILDS a decimal extractor for e from ONE continued-fraction convergent
      plus ONE modexp -- no series, no merge, the sqrt2 method -- and verifies
      it against mpmath. e's CF is patterned ([2;1,2,1,1,4,1,1,6,...]) so its
      convergents are computable without knowing e.
  (2) MEASURES redundancy for sqrt2, e, and pi.
  (3) States exactly what pi is missing.
"""

import math

from mpmath import mp, mpf, e as mp_e, sqrt as mp_sqrt

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


# ---------------------------------------------------------------------------
# (1) e's decimal extraction from ONE convergent (the sqrt2 method)
# ---------------------------------------------------------------------------

def e_cf_terms(n):
    """e = [2; 1,2,1, 1,4,1, 1,6,1, ...]:  a_{3k+2} = 2(k+1), else 1."""
    out = [2]
    i = 1
    while len(out) < n:
        out.append(2 * ((i - 2) // 3 + 1) if i % 3 == 2 else 1)
        i += 1
    return out


def e_convergent(min_q_bits):
    """Smallest CF convergent p/q of e with q >= 2^min_q_bits."""
    pm1, p = 1, 2          # p_{-1}, p_0
    qm1, q = 0, 1          # q_{-1}, q_0
    i = 1
    while q.bit_length() < min_q_bits:
        a = 2 * ((i - 2) // 3 + 1) if i % 3 == 2 else 1
        pm1, p = p, a * p + pm1
        qm1, q = q, a * q + qm1
        i += 1
    return mpz(p), mpz(q), i


def extract_from_convergent(p, q, N, window=10):
    """Decimal digits of p/q at positions N+1..N+window."""
    r = (pow(mpz(10), N, q) * (p % q)) % q
    return str((r * mpz(10) ** window) // q).zfill(window)


def reference_digits(const, N, window=10):
    mp.dps = N + window + 30
    v = mp_e if const == 'e' else mp_sqrt(2)
    s = mp.nstr(v, N + window + 15, strip_zeros=False)
    s = s.split('.')[1]
    return s[N:N + window]


def part1():
    print("=" * 74)
    print("(1) DECIMAL DIGIT EXTRACTION FOR e FROM ONE CF CONVERGENT")
    print("    no series, no BBP, no merge tree: one convergent + one modexp")
    print("=" * 74)
    print("  %-8s %-9s %-12s %-12s %-12s %s"
          % ("N", "CF terms", "q bits", "extracted", "reference", "ok"))
    ok_all = True
    for N in (100, 1000, 5000, 20000):
        guard = 40
        qbits = int((N + guard) * math.log2(10) / 2) + 64
        p, q, nterms = e_convergent(qbits)
        got = extract_from_convergent(p, q, N)
        ref = reference_digits('e', N)
        ok = got == ref
        ok_all &= ok
        print("  %-8d %-9d %-12d %-12s %-12s %s"
              % (N, nterms, q.bit_length(), got, ref, "OK" if ok else "MISMATCH"))
    print("  VERDICT: %s" % ("e extracts from ONE convergent, exactly like sqrt2"
                             if ok_all else "FAILED -- do not trust this face"))
    return ok_all


# ---------------------------------------------------------------------------
# (2) redundancy: approximant data bits per answer bit
# ---------------------------------------------------------------------------

def part2():
    print()
    print("=" * 74)
    print("(2) REDUNDANCY of the best computable rational approximant")
    print("    redundancy = (approximant data bits) / (answer bits)")
    print("    1.0 = attains the N log^2 N extraction floor; log N = one log up")
    print("=" * 74)
    print("  %-10s %-24s %-12s %-12s %-10s"
          % ("N", "constant / approximant", "data bits", "answer bits", "redund"))
    for N in (1000, 10000, 100000):
        ans = N * math.log2(10)                       # bits of answer needed

        # sqrt2: ONE Pell convergent, p,q ~ 10^(N/2); the CF is periodic so the
        # convergent comes from ONE fast exponentiation in Z[sqrt2].
        data = ans                                     # p and q, q ~ 2^(ans/2)
        print("  %-10d %-24s %-12.4e %-12.4e %-10.3f"
              % (N, "sqrt2 / Pell convergent", data, ans, data / ans))

        # e: ONE CF convergent. Measure the REAL numbers.
        qbits = int(ans / 2) + 64
        p, q, nterms = e_convergent(qbits)
        data = p.bit_length() + q.bit_length()
        print("  %-10s %-24s %-12.4e %-12.4e %-10.3f"
              % ("", "e / CF convergent", data, ans, data / ans))

        # pi: the best COMPUTABLE approximant family is binary splitting's P/Q
        # (Chudnovsky). Q = product of T ratio-denominators, T = ans/47.11.
        T = ans / 47.11
        data = 2 * T * (3 * math.log2(max(T, 2)) + math.log2(640320 ** 3 / 24))
        print("  %-10s %-24s %-12.4e %-12.4e %-10.3f"
              % ("", "pi / Chudnovsky P,Q", data, ans, data / ans))

        # pi: its CF convergents ARE tight (q ~ 2^(ans/2)) but the partial
        # quotients are not computable without knowing pi.
        print("  %-10s %-24s %-12.4e %-12.4e %-10.3f  %s"
              % ("", "pi / CF convergent", ans, ans, 1.0,
                 "<- tight but NOT computable"))
        print()


def part3():
    print("=" * 74)
    print("(3) WHAT PI IS MISSING")
    print("=" * 74)
    print("  sqrt2: convergents come from a GROUP LAW. (1+sqrt2)^n by repeated")
    print("         squaring in Z[sqrt2] gives p,q in log N ring mults. The")
    print("         approximant is reachable without the digits.")
    print("  e:     convergents come from a PATTERNED CF (equivalently, from")
    print("         factorial-denominator Pade approximants). Denominator mass")
    print("         T log T = O(N): tight. Reachable without the digits.")
    print("  pi:    its CF convergents are equally tight, but no rule generates")
    print("         them -- you need pi first. Every family that IS generated by")
    print("         a rule (Machin, BBP, Chudnovsky, Pade, Beukers) has")
    print("         denominator mass Theta(N log N): the approximant-numerator")
    print("         invariant. That log IS the third log.")
    print()
    print("  So the fifth paradigm has a one-line specification:")
    print("    ** a rational p/q with |pi - p/q| < 10^-N, q of O(N) bits,")
    print("       computable in N log^2 N WITHOUT computing pi's digits. **")
    print("  Find that and extraction is one modexp, exactly as for sqrt2 and e.")
    print("  It is not a new series. It is a group law, or a patterned CF, or")
    print("  any rule whose approximant data is information-theoretically tight.")


def main():
    ok = part1()
    part2()
    part3()
    print()
    print("=" * 74)
    print("STATUS: %s" % ("e joins sqrt2 at the extraction floor; pi still one "
                          "log above." if ok else "part 1 failed"))
    print("=" * 74)


if __name__ == "__main__":
    main()
