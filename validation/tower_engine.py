"""Validate the elementary engine of the Dwork tower:

(1) Collapse lemma (exact algebra, spot-checked): for i+2k=j,
      1/((8i+1)(8k+1)) = (2/(8i+1) + 1/(8k+1)) / (8j+3)
    => R_s[j] * (8j+3) = [difference of SINGLE harmonic sums over the
       P/N index regions] -- verify against the direct convolution.

(2) Complete-block inversion: sum_{v==1 (8), v < 2^L} v^{-1}  ==  sum v
    (mod 2^L) -- the source of the tower's O(s) precision.
"""
from fractions import Fraction

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

# ---- (1) collapsed formula vs direct convolution for R_s[j] ----
def a(T, i):
    return Fraction(1, 8 * i + 1) if 0 <= i < T else Fraction(0)

def Rs_direct(s, j):
    tot = Fraction(0)
    for k in range(j // 2 + 1):
        i = j - 2 * k
        tot += a(2 ** (s + 1), i) * a(2 ** (s - 1), k) - a(2 ** s, i) * a(2 ** s, k)
    return tot

def Rs_collapsed(s, j):
    # P region: i in [2^s, 2^{s+1}), k < 2^{s-1}, i+2k=j  (+)
    # N region: i < 2^s, k in [2^{s-1}, 2^s), i+2k=j      (-)
    tot = Fraction(0)
    for k in range(j // 2 + 1):
        i = j - 2 * k
        inP = (2 ** s <= i < 2 ** (s + 1)) and (k < 2 ** (s - 1))
        inN = (i < 2 ** s) and (2 ** (s - 1) <= k < 2 ** s)
        if inP:
            tot += Fraction(2, 8 * i + 1) + Fraction(1, 8 * k + 1)
        if inN:
            tot -= Fraction(2, 8 * i + 1) + Fraction(1, 8 * k + 1)
    return tot / (8 * j + 3)

ok = True
for s in (3, 4, 5):
    for j in (2 ** s, 2 ** s + 3, 2 ** s + 8, 3 * 2 ** (s - 1)):
        d1, d2 = Rs_direct(s, j), Rs_collapsed(s, j)
        if d1 != d2:
            ok = False
            print('MISMATCH s=%d j=%d' % (s, j))
print('(1) collapse lemma reproduces R_s[j]:', 'EXACT at all tested (s,j)' if ok else 'FAIL')

# ---- (2) complete-block inversion congruence ----
print('(2) block-sum congruence  sum v^{-1} == sum v  (mod 2^L), class 1 mod 8:')
for L in (6, 8, 10, 12):
    M = 1 << L
    sinv = sum(pow(v, -1, M) for v in range(1, M, 8)) % M
    sdir = sum(v for v in range(1, M, 8)) % M
    print('   L=%2d: %s' % (L, 'HOLDS' if sinv == sdir else 'fails (v2 of diff = %d)'
          % v2(Fraction(sinv - sdir))))

# ---- (3) consequence check: derived precision floor s+2 on R_s ----
print('(3) tower precision from block structure: min v2(R_s[j]) at boundary:')
for s in (3, 4, 5):
    mv = min(v2(Rs_direct(s, j)) for j in range(2 ** s, 2 ** s + 16))
    print('   s=%d: min v2 = %d (predicted >= %d)' % (s, mv, s + 2))
