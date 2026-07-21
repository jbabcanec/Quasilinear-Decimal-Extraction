"""Validate the seventh face: huge-exponent power sums linearize to E.

Claim: for u = 1 mod 8,  u^{2^r} = exp(2^r Log2 u)  gives
    sum_{k<=K} (1+8k)^{2^r}  ==  K  +  2^r * sum_{k<=K} Log2(1+8k)   (mod 2^{2r+5})
so computing the elementary power sum  S(K,r) = sum (1+8k)^{2^r} mod 2^{2r}
in Otilde(K+r) is EQUIVALENT to the barrier (E to r bits).
"""
from fractions import Fraction

K, r = 200, 120
M2 = 2 * r + 5
MOD = 1 << M2

# power sum, directly (repeated squaring per term -- the slow way, for truth)
S = 0
for k in range(1, K + 1):
    S = (S + pow(1 + 8 * k, 1 << r, MOD)) % MOD

# Log side, exact rationals
JL = (M2 + 60) // 3 + 4
def log2adic(u):
    y = Fraction(u - 1)
    s_, t_ = Fraction(0), Fraction(1)
    for j in range(1, JL):
        t_ *= y
        s_ += (Fraction(-1) ** (j + 1)) * t_ / j
    return s_

Lsum = Fraction(0)
for k in range(1, K + 1):
    Lsum += log2adic(1 + 8 * k)
Lval = (Lsum.numerator * pow(Lsum.denominator, -1, MOD)) % MOD

lhs = S % MOD
rhs = (K + (Lval << r)) % MOD
d = (lhs - rhs) % MOD
v = 0
dd = d
while dd and dd % 2 == 0:
    dd //= 2
    v += 1
print('power sum vs K + 2^r*LogSum: %s (agree to 2^%d of 2^%d)'
      % ('EXACT' if d == 0 else 'differ', v if d else M2, M2))
# recovery of E-precision: LogSum mod 2^r from the power sum alone
rec = ((S - K) % MOD) >> r
want = Lval % (1 << (r + 5 - 0))
print('recovered LogSum mod 2^%d: %s' % (r, 'MATCH' if (rec - Lval) % (1 << r) == 0 else 'FAIL'))
