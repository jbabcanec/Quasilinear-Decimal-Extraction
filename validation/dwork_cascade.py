"""Defect cascade probe: do the Dwork defects themselves descend?

D_s(x) := [F_{2^{s+1}}(x) F_{2^{s-1}}(x^2) - F_{2^s}(x) F_{2^s}(x^2)] / 2^{s+2}

Tests: (a) v2 profile of D_s across the boundary window (structure?);
(b) Cartier matching D_{s+1}[2j] vs D_s[j]: v2(D_{s+1}[2j] - D_s[j]) and
    v2(D_{s+1}[2j] - u*D_s[j]) for small odd u -- growth in s means the
    defects satisfy their own congruence => amplification cascade exists.
"""
from fractions import Fraction

D = 160
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

print('(a) v2 profile of D_s at boundary coeffs (first 8 nonvacuous):')
for s in range(2, S + 1):
    lo = 2 ** (s - 1)
    prof = [v2(defect[s][j]) for j in range(lo, min(lo + 8, D))]
    print('  s=%d: coeffs[%d..]: %s' % (s, lo, prof))

print()
print('(b) cascade test AT TRUE SUPPORT (j >= 2^s): v2(D_{s+1}[2j] - u*D_s[j]):')
for s in range(2, S):
    lo = 2 ** s
    hi = min(lo + 12, D // 2)
    best = None
    for u in (1, -1, 3, -3, 5, -5, 7, -7):
        vals = [v2(defect[s + 1][2 * j] - u * defect[s][j]) for j in range(lo, hi)]
        mv = min(vals)
        if best is None or mv > best[1]:
            best = (u, mv)
    raw = min(v2(defect[s + 1][2 * j]) for j in range(lo, hi))
    print('  s=%d->%d: raw min v2(D_{s+1}[2j]) = %s ; best u = %+d : min v2(diff) = %s'
          % (s, s + 1, raw, best[0], best[1]))

print()
print('(c) leading ratios D_{s+1}[2j]/D_s[j] mod 32 (diagnostic):')
for s in range(2, S):
    lo = 2 ** s
    out = []
    for j in range(lo, min(lo + 6, D // 2)):
        a, b = defect[s + 1][2 * j], defect[s][j]
        if b != 0 and v2(a) == v2(b):
            r = a / b
            num, den = r.numerator, r.denominator
            out.append((num * pow(den, -1, 32)) % 32)
        else:
            out.append('~')
    print('  s=%d->%d: %s' % (s, s + 1, out))
