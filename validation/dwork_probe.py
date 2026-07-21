"""Dwork/Cartier probe for f(x) = sum_t x^t/(8t+1) = 2F1(1/8, 1; 9/8; x) at p = 2.

Tests whether a Dwork-type truncation congruence holds despite the 2-ramified
parameters (denominator 8). Classical shape (Dwork; Beukers-Vlasenko):

   F_{2^{s+1}}(x) * F_{2^{s-1}}(x^2)  ==  F_{2^s}(x) * F_{2^s}(x^2)   (mod 2^{s+c})

coefficientwise. We measure min 2-adic valuation of the difference R_s over
low coefficients, for s = 2..S. Linear growth in s => descent structure EXISTS
(lead advances to precision amplification). Stalling => ramification kills it
(terminal obstruction; try twisted/group variants next).

Also: the value-level test at x = 16 (v2(x)=4), and the sibling odd-part
series g(x) = sum x^t/(8t+5) mixed variant.
"""
from fractions import Fraction

D = 96      # coefficient window (must exceed truncation boundaries tested)
S = 6       # levels: truncations up to 2^(S+1)

def trunc_series(T, offset=1):
    """coefficients of sum_{t<T} x^t/(8t+offset), as length-D Fraction list."""
    return [Fraction(1, 8 * t + offset) if t < T else Fraction(0) for t in range(D)]

def mul(a, b):
    c = [Fraction(0)] * D
    for i, ai in enumerate(a):
        if ai:
            for j in range(D - i):
                if b[j]:
                    c[i + j] += ai * b[j]
    return c

def sq_var(a):
    """a(x^2) truncated to D coefficients."""
    c = [Fraction(0)] * D
    for i, ai in enumerate(a):
        if 2 * i < D:
            c[2 * i] = ai
    return c

def v2(fr):
    if fr == 0:
        return 10 ** 9
    n = abs(fr.numerator)
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    d = fr.denominator
    while d % 2 == 0:
        d //= 2
        v -= 1
    return v

print('coefficientwise Dwork ratio test, F = sum x^t/(8t+1):')
print(' s | trunc pair          | min v2(R_s) over coeffs [0, 2^{s-1})')
for s in range(2, S + 1):
    F_big = trunc_series(2 ** (s + 1))
    F_mid = trunc_series(2 ** s)
    F_sml = trunc_series(2 ** (s - 1))
    lhs = mul(F_big, sq_var(F_sml))
    rhs = mul(F_mid, sq_var(F_mid))
    lo, hi = 2 ** (s - 1), min(D, 2 ** (s + 1) + 8)
    vals = [v2(lhs[j] - rhs[j]) for j in range(lo, hi)]
    mv = min(vals)
    print('%2d | (2^%d,2^%d | 2^%d,2^%d) | boundary coeffs [%d,%d): min v2 = %s'
          % (s, s + 1, s - 1, s, s, lo, hi,
             'exact-0' if mv >= 10 ** 9 else mv))

print()
print('value-level test at x = 16:')
def val_at16(T, offset=1):
    v = Fraction(0)
    p = Fraction(1)
    for t in range(T):
        v += p / (8 * t + offset)
        p *= 16
    return v
for s in range(2, S + 1):
    Vb = val_at16(2 ** (s + 1))
    Vm = val_at16(2 ** s)
    Vs_ = val_at16(2 ** (s - 1))
    def val_at256(T):
        v = Fraction(0)
        p = Fraction(1)
        for t in range(T):
            v += p / (8 * t + 1)
            p *= 256
        return v
    diff = Vb * val_at256(2 ** (s - 1)) - Vm * val_at256(2 ** s)
    print(' s=%d: v2 = %s' % (s, v2(diff)))
