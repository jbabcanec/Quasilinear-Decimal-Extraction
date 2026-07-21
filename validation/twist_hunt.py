"""Salie-twist and reflection-partner hunt on V.

In finite-field land, Salie sums (quadratic-character-twisted Kloosterman)
collapse to closed form while untwisted Kloosterman sums are hard. Test the
analog: do character-twisted versions of V = sum 16^s/(A-8s) mod 2^P relate
algebraically to V, to each other, or to completes? Also test the corrected
reflection partner Vr (denominators B+8s, B = 2*(A mod 8) - A), which no
previous hunt included.

Basis: V; six +/-1-twists of period 4 and 8 in s; Vr; interval log; 1.
Same calibrated LLL harness (controls at height O(1)).
"""

P = 240
WREL = 225
MODP = 1 << P
MODW = 1 << WREL
A = 483
T = 56

def inv(o):
    return pow(o % MODP, -1, MODP)

def vsum(twist):
    acc, p16 = 0, 1
    for s in range(T):
        acc = (acc + twist(s) * p16 * inv(A - 8 * s)) % MODP
        p16 = p16 * 16 % MODP
    return acc

V0 = vsum(lambda s: 1)
tw = {
    'V_alt':  vsum(lambda s: (-1) ** s),
    'V_p4a':  vsum(lambda s: 1 if s % 4 in (0, 1) else -1),
    'V_p4b':  vsum(lambda s: 1 if s % 4 in (0, 3) else -1),
    'V_p8a':  vsum(lambda s: 1 if s % 8 in (0, 1, 2, 3) else -1),
    'V_p8b':  vsum(lambda s: 1 if s % 8 in (0, 1, 4, 5) else -1),
    'V_p8c':  vsum(lambda s: 1 if s % 8 in (0, 3, 4, 7) else -1),
}

r0 = A % 8
B = 2 * r0 - A          # corrected reflection parameter (negative odd)
Vr, p16 = 0, 1
for s in range(T):
    Vr = (Vr + p16 * inv(B + 8 * s)) % MODP
    p16 = p16 * 16 % MODP

LG = 0
w = 8 * T * inv(A) % MODP
wp = 1
for n in range(1, 80):
    wp = wp * w % MODP
    e = (n & -n).bit_length() - 1
    LG = (LG - (wp >> e) * inv(n >> e)) % MODP

# ---- LLL harness ----
def lll(Bin):
    Bm = [list(r) for r in Bin]
    n = len(Bm)
    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))
    d = [0] * (n + 1)
    d[0] = 1
    lam = [[0] * n for _ in range(n)]
    def redi(k, l):
        if 2 * abs(lam[k][l]) > d[l + 1]:
            q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1])
            Bm[k] = [x - q * y for x, y in zip(Bm[k], Bm[l])]
            lam[k][l] -= q * d[l + 1]
            for i in range(l):
                lam[k][i] -= q * lam[l][i]
    def swapi(k, kmax):
        Bm[k], Bm[k - 1] = Bm[k - 1], Bm[k]
        for j in range(k - 1):
            lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
        lm = lam[k][k - 1]
        Bv = (d[k - 1] * d[k + 1] + lm * lm) // d[k]
        for i in range(k + 1, kmax + 1):
            t = lam[i][k]
            lam[i][k] = (d[k + 1] * lam[i][k - 1] - lm * t) // d[k]
            lam[i][k - 1] = (Bv * t + lm * lam[i][k]) // d[k + 1]
        d[k] = Bv
    def incremental(k):
        for j in range(k + 1):
            u = dot(Bm[k], Bm[j])
            for i in range(j):
                u = (d[i + 1] * u - lam[k][i] * lam[j][i]) // d[i]
            if j < k:
                lam[k][j] = u
            else:
                d[k + 1] = u
    kmax = 0
    incremental(0)
    k = 1
    while k < n:
        if k > kmax:
            kmax = k
            incremental(k)
        while True:
            redi(k, k - 1)
            if 4 * (d[k + 1] * d[k - 1] + lam[k][k - 1] ** 2) < 3 * d[k] ** 2:
                swapi(k, kmax)
                k = max(k - 1, 1)
            else:
                for l in range(k - 2, -1, -1):
                    redi(k, l)
                k += 1
                break
    return Bm

def hunt(names, vals, tag, maxcoef=1 << 40):
    dd = len(vals)
    K = MODW
    rows = []
    for i in range(dd):
        r = [0] * dd + [K * (vals[i] % MODW)]
        r[i] = 1
        rows.append(r)
    rows.append([0] * dd + [K * MODW])
    red = lll(rows)
    found = []
    for r in red:
        if r[-1] == 0 and any(r[:dd]):
            c = r[:dd]
            if max(abs(x) for x in c) < maxcoef:
                if sum(ci * vi for ci, vi in zip(c, vals)) % MODW == 0:
                    found.append(c)
    print('[%s] dim %d: %s' % (tag, dd,
          'RELATIONS: ' + '; '.join(
              str([(n_, c_) for n_, c_ in zip(names, c) if c_]) for c in found[:4])
          if found else 'EMPTY (coeffs < 2^40, mod 2^%d)' % WREL), flush=True)
    return found

# control
planted = (5 * V0 - 9 * tw['V_alt'] + 4 * LG) % MODW
hunt(['V', 'V_alt', 'LG', 'planted'], [V0, tw['V_alt'], LG, planted], 'control')

# the twist hunt
names = ['V'] + list(tw.keys()) + ['Vr', 'LG', '1']
vals = [V0] + list(tw.values()) + [Vr, LG, 1]
hunt(names, vals, 'salie-twist+reflection')

# second parameter point
A2, T2 = 417, 48
def vsum2(Aa, Tt, twist):
    acc, p16 = 0, 1
    for s in range(Tt):
        acc = (acc + twist(s) * p16 * inv(Aa - 8 * s)) % MODP
        p16 = p16 * 16 % MODP
    return acc
V0b = vsum2(A2, T2, lambda s: 1)
twb = [vsum2(A2, T2, f) for f in (
    lambda s: (-1) ** s,
    lambda s: 1 if s % 4 in (0, 1) else -1,
    lambda s: 1 if s % 8 in (0, 1, 2, 3) else -1)]
Bb = 2 * (A2 % 8) - A2
Vrb, p16 = 0, 1
for s in range(T2):
    Vrb = (Vrb + p16 * inv(Bb + 8 * s)) % MODP
    p16 = p16 * 16 % MODP
hunt(['V', 'V_alt', 'V_p4a', 'V_p8a', 'Vr', '1'],
     [V0b] + twb + [Vrb, 1], 'salie-twist+reflection@A2')
