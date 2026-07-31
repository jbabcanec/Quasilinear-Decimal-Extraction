#!/usr/bin/env python3
"""redundancy_screen.py -- a 10-second test for "can this idea reach parity?"
(2026-07-31)

Everything measured over the last two days collapses to ONE number. Let a
representation produce rational approximants p_n/q_n -> pi by a rule. Define

    REDUNDANCY = (data bits) / (accuracy bits)
               = log2(p_n * q_n) / log2(1 / |pi - p_n/q_n|)

Why this is THE number:
  * the data bits are exactly the mass fed to the assembly tree
    (log q_n = sum of the per-step entry sizes), so the tree costs
    M(data) * log n = redundancy * N log^2 N;
  * so redundancy FLAT in n  => N log^2 N = parity;
    redundancy GROWING       => N log^3 N = where we are.
  * it subsumes the g_eff series criterion, the product-face dichotomy and
    the approximant face: all three were this ratio in different clothes.

For a generalized continued fraction with partial numerators a_k, the error is
|Prod a_k| / (q_n q_{n+1}), so

    redundancy = 1 + (sum log a_k) / (accuracy bits).

A SIMPLE continued fraction (all a_k = 1) has redundancy exactly 1. That is
why sqrt2 and e reach the floor and pi does not: pi's simple-CF partial
quotients are not rule-generated, and every rule-generated family for pi is a
GENERALIZED CF whose numerator product sum log a_k is precisely the third log.

This screens every family we can name, including one never scored before (the
Beukers-type integrals behind the irrationality-measure bounds for pi).
"""

import math
from fractions import Fraction

from mpmath import mp, mpf, mpmathify, pi as mp_pi, e as mp_e, sqrt as mp_sqrt

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def score(p, q, err):
    """redundancy = data bits / accuracy bits."""
    if p == 0 or q == 0 or err == 0:
        return None
    data = math.log2(abs(p)) + math.log2(abs(q))
    acc = -math.log2(abs(err))
    return data, acc, data / acc


# ---------------------------------------------------------------------------
# families
# ---------------------------------------------------------------------------

def fam_pell(steps):
    """sqrt2 from (1+sqrt2)^n: CONSTANT matrix, so the product is a POWER."""
    mp.dps = 400
    tgt = mp_sqrt(2)
    p, q = mpz(1), mpz(1)
    out = []
    for i in range(steps):
        p, q = p + 2 * q, p + q
        if i > 4 and i % 8 == 0:
            r = score(p, q, mpf(p) / mpf(q) - tgt)
            if r:
                out.append(r)
    return out


def fam_simple_cf(const, steps):
    """simple CF convergents: a_k = 1 always, so redundancy must be 1."""
    mp.dps = 600
    tgt = {'e': mp_e, 'pi': mp_pi, 'sqrt2': mp_sqrt(2)}[const]
    x = tgt
    pm1, p, qm1, q = mpz(1), mpz(int(x)), mpz(0), mpz(1)
    x = x - int(x)
    out = []
    for i in range(steps):
        if x == 0:
            break
        x = 1 / x
        a = int(x)
        x -= a
        pm1, p = p, mpz(a) * p + pm1
        qm1, q = q, mpz(a) * q + qm1
        if i > 4 and i % 8 == 0:
            r = score(p, q, mpf(p) / mpf(q) - tgt)
            if r:
                out.append(r)
    return out


def fam_gcf(a_of_k, b_of_k, target, steps, sample=8):
    """generalized CF b_0 + a_1/(b_1 + a_2/(b_2 + ...)) by the matrix recurrence,
       with rational a_k, b_k cleared to integers."""
    mp.dps = 600
    pm1, p = Fraction(1), Fraction(b_of_k(0))
    qm1, q = Fraction(0), Fraction(1)
    out = []
    for k in range(1, steps + 1):
        a, b = Fraction(a_of_k(k)), Fraction(b_of_k(k))
        pm1, p = p, b * p + a * pm1
        qm1, q = q, b * q + a * qm1
        if k > 4 and k % sample == 0:
            val = mpf(p.numerator) / mpf(p.denominator) / \
                (mpf(q.numerator) / mpf(q.denominator))
            r = Fraction(p / q)
            sc = score(r.numerator, r.denominator, val - target)
            if sc:
                out.append(sc)
    return out


def fam_chudnovsky(nmax):
    """binary-splitting P/Q for Chudnovsky: the standard computable family."""
    out = []
    for n in (200, 800, 3200, 12800):
        if n > nmax:
            break
        # accuracy 47.11 bits/term; Q = prod of the ratio denominators
        acc = 47.11 * n
        data = 2 * n * (3 * math.log2(max(n, 2)) + math.log2(640320 ** 3 / 24))
        out.append((data, acc, data / acc))
    return out


def fam_beukers(nmax=11):
    """NEVER SCORED: the Beukers-type integrals behind pi's irrationality
       measure.  I_n = int_0^1 x^{4n}(1-x)^{4n}/(1+x^2) dx = R_n + b*pi/4,
       giving pi ~ -4R_n/b with error 4 I_n / |b|.  n=1 is the 22/7 classic."""
    mp.dps = 400
    out = []
    for n in range(1, nmax + 1):
        # P(x) = x^{4n} (1-x)^{4n}, coefficients ascending
        cf = [0] * (4 * n) + [math.comb(4 * n, j) * (-1) ** j
                              for j in range(4 * n + 1)]
        # divide by (1 + x^2): synthetic, descending
        P = list(cf)
        deg = len(P) - 1
        Q = [0] * (deg - 1)
        for d in range(deg, 1, -1):
            c = P[d]
            if c:
                Q[d - 2] = c
                P[d] -= c
                P[d - 2] -= c
        a, b = P[1], P[0]
        if a != 0:
            print("    n=%d: ln2 term present (a=%d), skipped" % (n, a))
            continue
        R = sum(Fraction(Q[j], j + 1) for j in range(len(Q)))     # int_0^1 Q
        approx = Fraction(-4, 1) * R / b
        mp.dps = 60 + 3 * n
        I = mp.quad(lambda x: x ** (4 * n) * (1 - x) ** (4 * n) / (1 + x * x),
                    [0, 1])
        err = 4 * I / abs(b)
        if err == 0:
            continue
        out.append(score(approx.numerator, approx.denominator, err))
    return out


def show(name, rows, note=""):
    if not rows:
        print("  %-34s (no data)" % name)
        return
    reds = [r[2] for r in rows]
    trend = reds[-1] / reds[0]
    verdict = "FLAT -> parity-capable" if trend < 1.15 else \
              "GROWING -> costs a log"
    print("  %-34s %-28s %-22s %s"
          % (name, " ".join("%.2f" % r for r in reds[:5]), verdict, note))


def main():
    print("=" * 100)
    print("REDUNDANCY SCREEN: data bits per accuracy bit, across increasing n")
    print("flat = the assembly tree costs M(N) log n = N log^2 N = PARITY")
    print("=" * 100)
    print("  %-34s %-28s %-22s %s"
          % ("family", "redundancy (increasing n)", "verdict", "note"))

    show("sqrt2: Pell (constant matrix)", fam_pell(200),
         "matrix POWER: telescopes to M(N)")
    show("sqrt2: simple CF", fam_simple_cf('sqrt2', 200))
    show("e: simple CF (patterned)", fam_simple_cf('e', 200),
         "rule-generated (Euler)")
    show("pi: simple CF", fam_simple_cf('pi', 200),
         "NOT rule-generated: needs pi")

    mp.dps = 600
    show("pi: Brouncker GCF (-> 4/pi)", fam_gcf(lambda k: (2 * k - 1) ** 2,
                                                lambda k: 2 if k else 1,
                                                4 / mp_pi, 400, 40),
         "a_k = (2k-1)^2 unbounded")
    z = Fraction(1, 5)
    show("arctan(1/5): Gauss GCF", fam_gcf(
        lambda k: z if k == 1 else (k - 1) ** 2 * z * z,
        lambda k: (2 * k - 1) if k else 0,
        mp.atan(mpf(1) / 5), 120, 8),
         "Machin engine; a_k ~ k^2")
    show("pi: Chudnovsky binary splitting", fam_chudnovsky(12800),
         "the practical champion")
    print("  scoring the Beukers/irrationality-measure family (never done):")
    show("pi: Beukers integrals", fam_beukers(11),
         "n=1 is the 22/7 classic")

    print()
    print("=" * 100)
    print("READING: redundancy 1.0 is the floor and only SIMPLE continued")
    print("fractions attain it. Every rule-generated family for pi is a")
    print("GENERALIZED CF, and its partial numerators a_k contribute")
    print("sum log a_k of pure overhead -- that sum IS the third log.")
    print()
    print("So the hunt target, in one line:")
    print("  ** a RULE generating 2x2 integer matrices with BOUNDED entries")
    print("     whose product's ratio converges geometrically to pi. **")
    print("  Bounded entries => data bits = Theta(accuracy bits) => redundancy")
    print("  O(1) => assembly M(N) log n = N log^2 N, then one exact power and")
    print("  one division for the digit. sqrt2 has it (a constant matrix).")
    print("  e has it (entries grow but accuracy grows with them).")
    print("=" * 100)


if __name__ == "__main__":
    main()
