#!/usr/bin/env python3
"""lemma_check.py -- exact checks of every stated identity in the paper.
(2026-07-31, from the full claim-by-claim read)

One test per claim, in paper order. Everything is exact integer or Fraction
arithmetic; nothing here uses floating point except where a decimal constant
in the text is being checked to its stated precision.
"""

import math
import random
from fractions import Fraction

OK = []


def claim(name, ok):
    OK.append(ok)
    print("  [%s] %s" % ("ok " if ok else "FAIL", name))


# --- eq:digit ---------------------------------------------------------------
def check_digit_formula():
    # d_N(x) = floor(10 * frac(10^(N-1) x)); test against the digits of 1/7
    x = Fraction(1, 7)                      # 0.142857142857...
    want = "142857142857142857"
    got = "".join(str((10 * ((10 ** (N - 1) * x) % 1)).__floor__())
                  for N in range(1, len(want) + 1))
    claim("eq:digit  d_N = floor(10 frac(10^{N-1} x))", got == want)


# --- lem:normal -------------------------------------------------------------
def check_normal_form():
    # 10^n * 16^-k * (4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6))
    #   = 2^{A+2}X/(8k+1) - 2^{A-1}X/(2k+1) - 2^A X/(8k+5) - 2^{A-1}X/(4k+3)
    ok = True
    for n in (7, 23, 60):
        X = Fraction(5) ** n
        for k in range(0, 12):
            A = n - 4 * k
            lhs = (Fraction(10) ** n / Fraction(16) ** k) * (
                Fraction(4, 8 * k + 1) - Fraction(2, 8 * k + 4)
                - Fraction(1, 8 * k + 5) - Fraction(1, 8 * k + 6))
            rhs = (Fraction(2) ** (A + 2) * X / (8 * k + 1)
                   - Fraction(2) ** (A - 1) * X / (2 * k + 1)
                   - Fraction(2) ** A * X / (8 * k + 5)
                   - Fraction(2) ** (A - 1) * X / (4 * k + 3))
            ok &= (lhs == rhs)
    claim("lem:normal  termwise normal form", ok)


# --- lem:band and lem:h -----------------------------------------------------
def check_band_and_h():
    ok_band = ok_h = True
    rnd = random.Random(11)
    for _ in range(400):
        X = rnd.getrandbits(rnd.randint(20, 300))
        J = rnd.randint(1, 200)
        mu = rnd.randrange(3, 1 << 20) | 1
        H, r = divmod(X, 1 << J)
        h, u = H % mu, Fraction(r, 1 << J)
        ok_band &= (Fraction(X, mu << J) % 1 == Fraction(h + u, mu))
        R, Yt = X % mu, (X % (1 << J)) % mu
        ok_h &= (h == (R - Yt) * pow(pow(2, J, mu), -1, mu) % mu)
    claim("lem:band  frac(X/(mu 2^J)) = (h+u)/mu", ok_band)
    claim("lem:h     h = (R - Ytilde) 2^{-J} mod mu", ok_h)
    # the paper's word-size inverse formula for 2^{-1}
    claim("2^{-J} = ((mu+1)/2)^J mod mu",
          all(pow(2, -J, mu) == pow((mu + 1) // 2, J, mu)
              for mu in (3, 5, 17, 8 * 977 + 1) for J in (1, 5, 64)))


# --- lem:euler and lem:fermat ----------------------------------------------
def check_euler_fermat():
    ok_e = ok_f = True
    for q in [p for p in range(3, 4000) if all(p % d for d in range(2, int(p ** .5) + 1))]:
        ok_f &= (pow(2, 4 * q, q) == 16 % q)
        k0 = next(k for k in range(q) if (8 * k + 1) % q == 0)
        s = (8 * k0 + 1) // q
        if s > 7 or s % 2 == 0:
            ok_e = False
            break
        leg = pow(2, (q - 1) // 2, q)            # Euler criterion value
        ok_e &= (pow(2, 4 * k0, q)
                 == pow(leg, s, q) * pow(2, (s - 1) // 2, q) % q)
    claim("lem:euler   2^{4k0} = (2|q)^s 2^{(s-1)/2},  s<=7 odd", ok_e)
    claim("lem:fermat  2^{4q} = 16 mod q", ok_f)


# --- lem:bezout -------------------------------------------------------------
def check_bezout():
    ok = True
    rnd = random.Random(5)
    for _ in range(300):
        X = rnd.getrandbits(rnd.randint(30, 400))
        J = rnd.randint(1, 160)
        mu = rnd.randrange(3, 1 << 20) | 1
        u = pow(pow(2, J, mu), -1, mu)
        w = pow(mu, -1, 1 << J)
        lhs = Fraction(X, mu << J) % 1
        rhs = (Fraction((X % mu) * u % mu, mu)
               + Fraction(X * w % (1 << J), 1 << J)) % 1
        ok &= (lhs == rhs)
    claim("lem:bezout  odd/dyadic split", ok)


# --- lem:consolidate --------------------------------------------------------
def check_consolidate():
    rnd = random.Random(7)
    ok = True
    for _ in range(60):
        T = rnd.randint(2, 12)
        X = rnd.getrandbits(rnd.randint(40, 300))
        terms = [(rnd.choice((1, -1)), rnd.randint(1, 90),
                  rnd.randrange(3, 1 << 16) | 1) for _ in range(T)]
        Jmax = max(J for _, J, _ in terms)
        direct = sum(sg * Fraction(X * pow(mu, -1, 1 << J) % (1 << J), 1 << J)
                     for sg, J, mu in terms) % 1
        V = sum(sg * (1 << (Jmax - J)) * pow(mu, -1, 1 << Jmax)
                for sg, J, mu in terms) % (1 << Jmax)
        ok &= (direct == Fraction(X * V % (1 << Jmax), 1 << Jmax) % 1)
    claim("lem:consolidate  dyadic halves collapse into V", ok)


# --- assorted stated identities ---------------------------------------------
def check_misc():
    # Strassnitzky-Dase, as the paper writes it in Z[i]
    def mul(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    claim("(2+i)(5+i)(8+i) = 65(1+i)",
          mul(mul((2, 1), (5, 1)), (8, 1)) == (65, 65))
    # and that it really is pi/4
    claim("arctan(1/2)+arctan(1/5)+arctan(1/8) = pi/4",
          abs(math.atan(.5) + math.atan(.2) + math.atan(.125) - math.pi / 4)
          < 1e-15)
    # rising-factorial form of the denominator product
    K = 9
    lhs = math.prod(8 * k + 1 for k in range(K + 1))
    rhs = Fraction(8) ** (K + 1) * math.prod(
        Fraction(1, 8) + j for j in range(K + 1))
    claim("prod(8k+1) = 8^{K+1} (1/8)^{rising K+1}", lhs == rhs)
    # (Z/2^m)* = <-1> x <3>, and the 1-mod-8 units are exactly <9>
    ok = True
    for m in (5, 8, 11):
        M = 1 << m
        gen3 = {pow(3, j, M) for j in range(1 << (m - 2))}
        ok &= (len(gen3) == 1 << (m - 2))
        ok &= ({u for u in range(1, M, 2)} ==
               {s * g % M for s in (1, M - 1) for g in gen3})
        gen9 = {pow(9, j, M) for j in range(1 << (m - 3))}
        ok &= (gen9 == {u for u in range(1, M, 8)})
    claim("(Z/2^m)* = <-1>x<3>, and {u = 1 mod 8} = <9> of order 2^{m-3}", ok)
    # bit length of X, as the paper defines L
    claim("L = ceil(n log2 5) is the bit length of 5^n (n>=1)",
          all((5 ** n).bit_length() == math.ceil(n * math.log2(5))
              for n in range(1, 4000)))


def main():
    print("Exact checks of the paper's stated identities:")
    check_digit_formula()
    check_normal_form()
    check_band_and_h()
    check_euler_fermat()
    check_bezout()
    check_consolidate()
    check_misc()
    print()
    print("ALL PASSED" if all(OK) else "SOME FAILED")
    return 0 if all(OK) else 1


if __name__ == "__main__":
    raise SystemExit(main())
