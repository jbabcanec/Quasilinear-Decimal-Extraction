"""LLL relation hunt: mine S(K) = sum (1+8k)^{2^r} mod 2^{2r+5} data for exact
integer relations -- the fingerprint any hidden reciprocity/descent must leave.

Components per (r,K):  S(K), S(2K), S(4K),  2^r*D(K), 2^r*D(2K), 2^r*D(4K),
                       B(K), B(2K), B(4K),  K, 1
with P = product, D = dlog9(P), B = boundary term (1+8K)^{2^r}.

Positive controls: the seventh-face identity  S(cK) - 2^r D(cK) - cK = 0 (mod 2^{2r+5})
must be found for c = 1,2,4. Any OTHER short relation that persists across
(r,K) instances is a crack fingerprint. A second lattice excludes D (pure-S hunt).
"""
from fractions import Fraction

def log2adic_mod(u, M):
    """Log2(u) mod 2^M for u = 1 mod 8, exact via truncated series w/ Fractions."""
    y = u - 1
    JL = (M + 40) // 3 + 4
    num, den = 0, 1
    t = 1
    for j in range(1, JL):
        t *= y
        term_num = (-1) ** (j + 1) * t
        num = num * j + term_num * den * (1 if True else 1)
        den *= j
        # keep fraction reduced modestly small: reduce mod nothing; sizes ok for JL small
    # num/den is Log; reduce mod 2^M (den odd after cancellation guaranteed by valuations)
    from math import gcd
    g = gcd(num, den)
    num //= g
    den //= g
    assert den % 2 == 1
    return (num * pow(den, -1, 1 << M)) % (1 << M)

def dlog9(u, M):
    L9 = log2adic_mod(9, M + 8)
    Lu = log2adic_mod(u % (1 << (M + 40)), M + 8)
    v = 3  # v2(Log u) >= 3, v2(Log 9) = 3
    return ((Lu >> v) * pow(L9 >> v, -1, 1 << M)) % (1 << M)

def lll(B, delta=Fraction(99, 100)):
    B = [row[:] for row in B]
    n = len(B)
    def gso(B):
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in B[i]]
            for j in range(i):
                d = sum(a * b for a, b in zip(Bs[j], Bs[j]))
                mu[i][j] = sum(Fraction(a) * b for a, b in zip(B[i], Bs[j])) / d if d else Fraction(0)
                v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu
    Bs, mu = gso(B)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                B[k] = [a - q * b for a, b in zip(B[k], B[j])]
                Bs, mu = gso(B)
        nk = sum(x * x for x in Bs[k])
        nk1 = sum(x * x for x in Bs[k - 1])
        if nk >= (delta - mu[k][k - 1] ** 2) * nk1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gso(B)
            k = max(k - 1, 1)
    return B

def hunt(r, K, include_D=True):
    M = 2 * r + 5
    MOD = 1 << M
    xs, names = [], []
    P = {}
    S = {}
    for c in (1, 2, 4):
        Kc = c * K
        s = 0
        p = 1
        big = 1 << (M + 48)
        for k in range(1, Kc + 1):
            u = 1 + 8 * k
            s = (s + pow(u, 1 << r, MOD)) % MOD
            p = (p * u) % big
        S[c], P[c] = s, p
    for c in (1, 2, 4):
        xs.append(S[c]); names.append('S(%dK)' % c)
    if include_D:
        for c in (1, 2, 4):
            xs.append((log2adic_mod(P[c], M) << r) % MOD); names.append('2^r*LogP(%dK)' % c)
    for c in (1, 2, 4):
        xs.append(pow(1 + 8 * c * K, 1 << r, MOD)); names.append('B(%dK)' % c)
    xs.append(K % MOD); names.append('K')
    xs.append(1); names.append('1')
    n = len(xs)
    rows = [[1 if i == j else 0 for j in range(n)] + [xs[i]] for i in range(n)]
    rows.append([0] * n + [MOD])
    red = lll(rows)
    rels = []
    for v in red:
        if v[-1] == 0 and any(v[:n]) and max(abs(c) for c in v[:n]) <= (1 << 14):
            # normalize sign
            for c in v[:n]:
                if c:
                    if c < 0:
                        v = [-a for a in v]
                    break
            rels.append(tuple(v[:n]))
    return names, sorted(set(rels))

for include_D, tag in ((True, 'WITH dlog components'), (False, 'PURE-S hunt')):
    print('=== %s ===' % tag)
    common = None
    for (r, K) in ((30, 700), (30, 1000), (41, 900)):
        names, rels = hunt(r, K, include_D)
        print(' (r=%d,K=%d): %d relation(s)' % (r, K, len(rels)))
        for rel in rels:
            print('   ', ' + '.join('%d*%s' % (c, names[i]) for i, c in enumerate(rel) if c), '= 0')
        common = set(rels) if common is None else (common & set(rels))
    print(' PERSISTENT across all instances:', len(common) if common else 0)
    for rel in (common or []):
        print('   >>>', ' + '.join('%d*%s' % (c, names[i]) for i, c in enumerate(rel) if c), '= 0')
