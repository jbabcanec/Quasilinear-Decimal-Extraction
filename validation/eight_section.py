"""Validate the 8-section identity in Q_2(zeta_8):

  sum_{t>=0} 16^t/(8t+1)  =  (1/(8y)) * sum_{zeta^8=1} zeta^{-1} * (-Log(1 - zeta*y))

with y = sqrt(2) = zeta8 + zeta8^{-1}, working in K = Q_2[x]/(x^4+1)
(zeta8 = x, ramification index 2). All arithmetic: polynomials of degree < 4
with Fraction coefficients; Log via its series (converges: v(zeta*y) = 1/2).

If EXACT to the working precision: the complete 2-adic sum is 8 ramified
logarithms => quasi-linear; the barrier is purely the TRUNCATION (8th face).
"""
from fractions import Fraction

PREC_TERMS_MAIN = 220   # terms of sum 16^t/(8t+1): valuation 4t
PREC_TERMS_LOG = 1800   # terms of Log series: valuation ~ n/2
CHECK_BITS = 800        # compare to this many 2-adic bits

# ---- K = Q_2[x]/(x^4+1): elements = [c0,c1,c2,c3] Fractions ----
def kmul(a, b):
    c = [Fraction(0)] * 7
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    c[i + j] += ai * bj
    # reduce x^4 = -1
    return [c[0] - c[4], c[1] - c[5], c[2] - c[6], c[3]]

def kadd(a, b):
    return [x + y for x, y in zip(a, b)]

def kscale(a, s):
    return [x * s for x in a]

ZETA = [Fraction(0), Fraction(1), Fraction(0), Fraction(0)]         # zeta8 = x
def kpow(a, e):
    r = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    for _ in range(e):
        r = kmul(r, a)
    return r

# y = sqrt(2) = x + x^{-1} = x - x^3   (since x^{-1} = -x^3 in x^4=-1)
Y = [Fraction(0), Fraction(1), Fraction(0), Fraction(-1)]
# sanity: y^2 should be 2
y2 = kmul(Y, Y)
assert y2 == [Fraction(2), Fraction(0), Fraction(0), Fraction(0)], y2

def klog1m(z, nterms):
    """-Log(1 - z) = sum_{n>=1} z^n / n in K."""
    s = [Fraction(0)] * 4
    p = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    for n in range(1, nterms + 1):
        p = kmul(p, z)
        s = kadd(s, kscale(p, Fraction(1, n)))
    return s

# RHS = (1/(8y)) * sum_zeta zeta^{-1} * (-Log(1 - zeta y))
rhs = [Fraction(0)] * 4
for k in range(8):
    zk = kpow(ZETA, k)            # zeta = zeta8^k
    zk_inv = kpow(ZETA, (8 - k) % 8)
    term = klog1m(kmul(zk, Y), PREC_TERMS_LOG)
    rhs = kadd(rhs, kmul(zk_inv, term))
# divide by 8y: multiply by y/(8*2) since 1/y = y/2
rhs = kscale(kmul(rhs, Y), Fraction(1, 16))

# LHS = complete sum, a plain rational
lhs = Fraction(0)
p = Fraction(1)
for t in range(PREC_TERMS_MAIN):
    lhs += p / (8 * t + 1)
    p *= 16

# compare in K: rhs should be [lhs, 0, 0, 0] up to 2-adic precision
def v2frac(fr):
    if fr == 0:
        return 10 ** 9
    n, d = abs(fr.numerator), fr.denominator
    v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    while d % 2 == 0:
        d //= 2; v -= 1
    return v

diff0 = v2frac(rhs[0] - lhs)
others = [v2frac(c) for c in rhs[1:]]
print('v2(rhs[0] - lhs) = %s  (want >= %d)' % (diff0, CHECK_BITS))
print('v2(rhs[1..3])    = %s  (want all >= %d)' % (others, CHECK_BITS))
ok = diff0 >= CHECK_BITS and all(o >= CHECK_BITS for o in others)
print('8-SECTION IDENTITY:', 'VALIDATED to %d bits' % CHECK_BITS if ok else 'MISMATCH')
