"""Validate the 2-adic Euler-Maclaurin identity for S = sum_{k=1}^K Log2(1+8k):

  S = F(K) - F(0) + (f(K)-f(0))/2 + sum_{i>=1} B_{2i}/(2i)! * (f^{(2i-1)}(K) - f^{(2i-1)}(0))

with f(x) = Log2(1+8x)  (2-adic log, convergent series),
     F(x) = ((1+8x)(Log2(1+8x) - 1) + 1)/8   (antiderivative, F(0)=0),
     f^{(r)}(x) = (-1)^{r-1} (r-1)! 8^r / (1+8x)^r,
and the Bernoulli series 2-adically CONVERGENT (term valuation ~ 6i).

Everything is computed in exact rational arithmetic; both sides are reduced
mod 2^m at the end (denominators must be odd). A match validates the
equivalence  E <-> 2-adic log-gamma evaluation  claimed in Paper II.
"""
from fractions import Fraction

M = 224
K = 257
JLOG = (M + 60) // 3 + 4    # Log series truncation (term val = 3j - v2(j) > M+40)
IMAX = (M + 80) // 6 + 6    # Bernoulli series truncation

def log2adic(kx8):
    """Log2(1+y) with y = 8k, as exact Fraction of the truncated series."""
    y = Fraction(kx8)
    s = Fraction(0)
    t = Fraction(1)
    for j in range(1, JLOG):
        t *= y
        s += (Fraction(-1) ** (j + 1)) * t / j
    return s

def bernoulli_list(n):
    """B_0..B_n exact (B_1 = -1/2 convention) via the standard recurrence."""
    B = [Fraction(0)] * (n + 1)
    B[0] = Fraction(1)
    from math import comb
    for mdeg in range(1, n + 1):
        acc = Fraction(0)
        for j in range(mdeg):
            acc += comb(mdeg + 1, j) * B[j]
        B[mdeg] = -acc / (mdeg + 1)
    return B

def mod2m(fr, m):
    num, den = fr.numerator, fr.denominator
    assert den % 2 == 1, 'denominator not odd: %s' % den
    return (num * pow(den, -1, 1 << m)) % (1 << m)

# ---- direct side ----
direct = Fraction(0)
for k in range(1, K + 1):
    direct += log2adic(8 * k)

# ---- Euler-Maclaurin side ----
fK = log2adic(8 * K)
f0 = Fraction(0)
FK = (Fraction(1 + 8 * K) * (fK - 1) + 1) / 8
F0 = Fraction(0)
em = FK - F0 + (fK - f0) / 2

B = bernoulli_list(2 * IMAX + 2)
from math import factorial
for i in range(1, IMAX + 1):
    r = 2 * i - 1
    # f^{(r)}(x) = (-1)^{r-1} (r-1)! 8^r / (1+8x)^r
    dK = (Fraction(-1) ** (r - 1)) * factorial(r - 1) * Fraction(8 ** r, (1 + 8 * K) ** r)
    d0 = (Fraction(-1) ** (r - 1)) * factorial(r - 1) * Fraction(8 ** r, 1)
    em += B[2 * i] / factorial(2 * i) * (dK - d0)

lhs = mod2m(direct, M)
rhs = mod2m(em, M)
print('direct mod 2^%d: %x' % (M, lhs))
print('euler-maclaurin: %x' % rhs)
print('MATCH' if lhs == rhs else 'MISMATCH (diff valuation: v2 = %d)'
      % ((lhs - rhs) % (1 << M)).bit_length() if lhs != rhs else 'MATCH')
if lhs != rhs:
    d = (lhs - rhs) % (1 << M)
    v = 0
    while d % 2 == 0 and d:
        d //= 2
        v += 1
    print('agreement to 2^%d of 2^%d' % (v, M))
