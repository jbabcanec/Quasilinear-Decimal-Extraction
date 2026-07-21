"""Validate the fifth-face identities (Gamma-transform formulation of E).

(1) Sparse-rational Amice transform of the counting measure on 1+8Z in [0,2^t):
      sum_{u<2^t, u=1 mod 8} (1+T)^u  ==  (1+T)((1+T)^{2^t}-1)/((1+T)^8-1)
    checked as power series mod (T^R, 2^m).
(2) E-extraction consistency: sum dlog_9(u) == Log2(prod u)/Log2(9) mod 2^m,
    with dlog_9(u) = Log2(u)/Log2(9) for u = 1 mod 8.
"""
from fractions import Fraction

t, R, m = 9, 6, 96   # 2^t range, series degree, 2-adic precision
MOD = 1 << m

# ---- (1) power series check ----
def poly_mul(a, b):
    c = [0] * R
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if i + j < R:
                    c[i + j] = (c[i + j] + ai * bj) % MOD
    return c

def poly_pow(base, e):
    r = [1] + [0] * (R - 1)
    while e:
        if e & 1:
            r = poly_mul(r, base)
        base = poly_mul(base, base)
        e >>= 1
    return r

T1 = [1, 1] + [0] * (R - 2)          # (1+T)
lhs = [0] * R
for u in range(1, 1 << t, 8):
    lhs = [(x + y) % MOD for x, y in zip(lhs, poly_pow(T1, u))]

num = poly_mul(T1, [(x - y) % MOD for x, y in zip(poly_pow(T1, 1 << t),
                                                  [1] + [0] * (R - 1))])
den = [(x - y) % MOD for x, y in zip(poly_pow(T1, 8), [1] + [0] * (R - 1))]
# den = 8T + 28T^2 + ... : divide num by den as power series over Z/2^m:
# lowest term of den is 8T -> factor T out of both, then invert (8 + 28T + ...)
num_s = num[1:] + [0]
den_s = den[1:] + [0]
# den_s[0] = 8: not invertible mod 2^m -- but num_s is divisible by 8 too in
# the quotient sense; check via cross-multiplication instead: lhs * den == num
chk = poly_mul(lhs, den)
ok1 = all((chk[i] - num[i]) % MOD == 0 for i in range(R))
print('(1) Amice sparse-rational form:', 'EXACT' if ok1 else 'MISMATCH')

# ---- (2) E two ways ----
JL = (m + 40) // 3 + 4
def log2adic(u):          # Log2(u) for u = 1 mod 8, exact Fraction (truncated)
    y = Fraction(u - 1)
    s_, t_ = Fraction(0), Fraction(1)
    for j in range(1, JL):
        t_ *= y
        s_ += (Fraction(-1) ** (j + 1)) * t_ / j
    return s_

def val2m(fr):
    return (fr.numerator * pow(fr.denominator, -1, MOD)) % MOD

L9 = log2adic(9)
E_direct = Fraction(0)
prod = 1
for u in range(9, 1 << t, 8):
    E_direct += log2adic(u) / L9
    prod = prod * u % (1 << (m + 60))
# careful: Log(prod mod 2^big) needs prod's low bits only up to precision
E_via_prod = log2adic(prod % (1 << (m + 48))) / L9
d = (val2m(E_direct) - val2m(E_via_prod)) % MOD
# expect agreement modulo the group order caveat (2^{m-3}-ish)
v = 0
dd = d
while dd and dd % 2 == 0:
    dd //= 2
    v += 1
print('(2) E direct vs Log(product): agree to 2^%d of 2^%d %s'
      % (v if d else m, m, '(EXACT)' if d == 0 else ''))
