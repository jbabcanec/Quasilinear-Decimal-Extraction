"""Identify the cascade multiplier: does psi(D_{s+1}) = U(x) * D_s(x) for a
FIXED series U independent of s?

d_s(x) := D_s(x) / x^{2^s}  (normalized to start at constant term)
U_s(x) := psi(D_{s+1})(x) / D_s(x)  computed as a power-series quotient at
the aligned support; if v2(U_{s+1} - U_s) grows with s, the multiplier is a
stable 2-adic series => the defect cascade EXISTS with explicit U =>
iterate for amplification.
"""
from fractions import Fraction

D = 200
S = 6

def trunc_series(T):
    return [Fraction(1, 8 * t + 1) if t < T else Fraction(0) for t in range(D)]

def mul(a, b):
    c = [Fraction(0)] * D
    for i, ai in enumerate(a):
        if ai:
            for j in range(D - i):
                if b[j]:
                    c[i + j] += ai * b[j]
    return c

def sq_var(a):
    c = [Fraction(0)] * D
    for i, ai in enumerate(a):
        if 2 * i < D:
            c[2 * i] = ai
    return c

def v2(fr):
    if fr == 0:
        return 10 ** 9
    n, d = abs(fr.numerator), fr.denominator
    v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    while d % 2 == 0:
        d //= 2; v -= 1
    return v

defect = {}
for s in range(2, S + 1):
    F_big = trunc_series(2 ** (s + 1))
    F_mid = trunc_series(2 ** s)
    F_sml = trunc_series(2 ** (s - 1))
    lhs = mul(F_big, sq_var(F_sml))
    rhs = mul(F_mid, sq_var(F_mid))
    defect[s] = [(lhs[j] - rhs[j]) / Fraction(2 ** (s + 2)) for j in range(D)]

def norm(s):
    """d_s = D_s / x^{2^s}, first L coefficients."""
    L = 12
    base = 2 ** s
    return [defect[s][base + j] if base + j < D else Fraction(0) for j in range(L)]

def psi_norm(s):
    """psi(D_{s+1}) / x^{2^s}: psi picks D_{s+1}[2j]; support of D_{s+1} is 2^{s+1}
    so psi-support starts at j = 2^s. Return first L coefficients from there."""
    L = 12
    base = 2 ** s
    return [defect[s + 1][2 * (base + j)] if 2 * (base + j) < D else Fraction(0)
            for j in range(L)]

def series_div(a, b):
    """a/b as power series, b[0] != 0, length L."""
    L = len(a)
    q = [Fraction(0)] * L
    for j in range(L):
        acc = a[j]
        for i in range(j):
            acc -= q[i] * b[j - i]
        q[j] = acc / b[0]
    return q

U = {}
for s in range(3, S):
    a, b = psi_norm(s), norm(s)
    if b[0] == 0:
        print('s=%d: leading coeff zero, skip' % s)
        continue
    U[s] = series_div(a, b)
    lead = [(u.numerator * pow(u.denominator, -1, 256)) % 256 if u.denominator % 2 else '~'
            for u in U[s][:6]]
    print('s=%d: U_s mod 256 (first 6): %s' % (s, lead))

print()
ss = sorted(U)
for i in range(len(ss) - 1):
    s1, s2 = ss[i], ss[i + 1]
    diffs = [v2(U[s2][j] - U[s1][j]) for j in range(8)]
    print('v2(U_%d - U_%d) coeffwise: %s' % (s2, s1, diffs))
