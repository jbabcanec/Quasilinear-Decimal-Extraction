"""Cross-formula tower probe: does the Dwork tower structure depend on the
BBP formula family? Compare denominators 8t+1 (hex BBP, ramified atom,
block bonus +3) vs 10t+1 (Bellard, unramified atom, block bonus +1) vs
6t+1 (control: 6 = 2*3) vs 4t+1 (control: 4 = 2^2).

Measure min v2 of R_s at the boundary for each family; the collapse lemma
generalizes (2(dk+1) + (di+1) = d*j+3 for i+2k=j), so towers should exist
for all d with precision tracking v2(d)-alignment. Prediction: bonus tracks
v2(d) + 2. Data decides.
"""
from fractions import Fraction

D = 96
S = 6

def v2(fr):
    if fr == 0:
        return 10 ** 9
    n, d_ = abs(fr.numerator), fr.denominator
    v = 0
    while n % 2 == 0:
        n //= 2; v += 1
    while d_ % 2 == 0:
        d_ //= 2; v -= 1
    return v

def tower_profile(d):
    def trunc(T):
        return [Fraction(1, d * t + 1) if t < T else Fraction(0) for t in range(D)]
    def mul(a, b):
        c = [Fraction(0)] * D
        for i, ai in enumerate(a):
            if ai:
                for j in range(D - i):
                    if b[j]:
                        c[i + j] += ai * b[j]
        return c
    def sq(a):
        c = [Fraction(0)] * D
        for i, ai in enumerate(a):
            if 2 * i < D:
                c[2 * i] = ai
        return c
    out = []
    for s in range(2, S + 1):
        lhs = mul(trunc(2 ** (s + 1)), sq(trunc(2 ** (s - 1))))
        rhs = mul(trunc(2 ** s), sq(trunc(2 ** s)))
        lo, hi = 2 ** s, min(D, 2 ** (s + 1) + 8)
        mv = min(v2(lhs[j] - rhs[j]) for j in range(lo, hi))
        out.append(mv if mv < 10 ** 8 else 'exact')
    return out

print('family      | min v2(R_s) at boundary, s = 2..%d  | v2(d)' % S)
for d in (8, 10, 6, 4, 16, 12):
    prof = tower_profile(d)
    vd = (d & -d).bit_length() - 1
    print('d = %-2d      | %-28s | %d' % (d, prof, vd))
