"""Validate the H-probe EM identity (digamma tail):

  H(K) = sum_{k=1}^{K} (8k+1)^{-1}  mod 2^m
       = F(K) + f(K)/2 + sum_i B_{2i}/(2i)! (f^{(2i-1)}(K) - f^{(2i-1)}(0))

with f(x) = (8x+1)^{-1} (2-adic), F = antiderivative = Log2(8x+1)/8,
f^{(r)}(x) = (-1)^r r! 8^r (8x+1)^{-(r+1)}.  2-adically convergent
(i-th correction valuation ~ 6i).  If exact: the unweighted inverse-sum
probe is EM-reducible to the DERIVATIVE of the cracked log-gamma tail --
confirming the height-barrier classification (weights innocent,
derivative-structure is the enemy).
"""
from fractions import Fraction
from math import comb, factorial

M, K = 200, 150
MOD = 1 << M
IMAX = (M + 80) // 6 + 6
JLOG = (M + 60) // 3 + 4

def log2adic(y):
    s, t = Fraction(0), Fraction(1)
    for j in range(1, JLOG):
        t *= y
        s += (Fraction(-1) ** (j + 1)) * t / j
    return s

def bern(n):
    B = [Fraction(0)] * (n + 1)
    B[0] = Fraction(1)
    for md in range(1, n + 1):
        acc = Fraction(0)
        for j in range(md):
            acc += comb(md + 1, j) * B[j]
        B[md] = -acc / (md + 1)
    return B

def mod2m(fr):
    assert fr.denominator % 2 == 1
    return (fr.numerator * pow(fr.denominator, -1, MOD)) % MOD

# direct
H = Fraction(0)
for k in range(1, K + 1):
    H += Fraction(1, 8 * k + 1)

# EM side
f = lambda x: Fraction(1, 8 * x + 1)
FK = log2adic(Fraction(8 * K)) / 8          # F(K)-F(0) = Log2(1+8K)/8
em = FK + f(K) / 2 - f(0) / 2 + Fraction(1, 2) * 0
# careful with the constant: EM: sum_{k=1}^{K} f(k) = F(K)-F(0) + (f(K)-f(0))/2 + ...
em = FK + (f(K) - f(0)) / 2
B = bern(2 * IMAX + 2)
for i in range(1, IMAX + 1):
    r = 2 * i - 1
    dK = (Fraction(-1) ** r) * factorial(r) * Fraction(8 ** r) / Fraction((8 * K + 1) ** (r + 1))
    d0 = (Fraction(-1) ** r) * factorial(r) * Fraction(8 ** r)
    em += B[2 * i] / factorial(2 * i) * (dK - d0)

a, b = mod2m(H), mod2m(em)
print('H direct  mod 2^%d: %x' % (M, a))
print('H via EM (digamma tail): %x' % b)
print('MATCH' if a == b else 'MISMATCH')
