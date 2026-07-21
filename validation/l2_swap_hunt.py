"""L2 swap-hunt: does a Lerch-FE-type (argument <-> parameter) reciprocity
exist 2-adically for the incomplete geometric-harmonic sum?

Object (true V geometry): W = sum_{s<T} 16^s/(A-8s) mod 2^P, A odd ~ 2P,
T just before A/8 -- the T-term, slope-8 incomplete sum.

Swap thesis (Landsberg-Schaar / Lerch FE shape): W relates algebraically to
an 8-term, slope-T dual over Q_2(zeta_8), possibly with prefactor 16^T,
plus complete logs. If ANY such relation exists at low height, LLL finds it;
the real-analytic version has a transcendental Fresnel remainder -- the
question is whether 2-adic discreteness kills the remainder.

Basis: W; components of linear dual D1 = sum_k zeta8^k/(A-Tk) and quadratic
dual D2 = sum_k zeta8^{k^2}/(A-Tk); same scaled by 16^T; the interval log
LG = Log2(1 - 8T/A); 1. All mod 2^150. Controls: planted relation.
"""

P = 240
WREL = 225
MODP = 1 << P
MODW = 1 << WREL
A = 483
T = 56

def inv(o):
    return pow(o, -1, MODP)

# ---- the incomplete object ----
Wobj = 0
p16 = 1
for s in range(T):
    Wobj = (Wobj + p16 * inv(A - 8 * s)) % MODP
    p16 = p16 * 16 % MODP

# ---- duals over Q_2(zeta_8): components of sum_k zeta8^{e(k)}/(A - T k) ----
def dual(efun):
    comp = [0, 0, 0, 0]
    for k in range(8):
        e = efun(k) % 8
        c, sg = e % 4, (-1) ** (e // 4)
        comp[c] = (comp[c] + sg * inv(A - T * k)) % MODP
    return comp

D1 = dual(lambda k: k)
D2 = dual(lambda k: k * k)

# ---- interval log: Log2(1 - 8T/A) = -sum (8T/A)^n / n ----
LG = 0
w = 8 * T * inv(A) % MODP
wp = 1
for n in range(1, 80):
    wp = wp * w % MODP
    e = (n & -n).bit_length() - 1
    LG = (LG - (wp >> e) * inv(n >> e)) % MODP

S16T = pow(16, T, MODP)

# ---- LLL harness (all-integer, Cohen 2.6.7) ----
import importlib.util as _il
spec = _il.spec_from_file_location('agm', 'agm_law_hunt_lll.py')
# inline instead: copy of the lll + hunt utilities
def lll(Bin):
    B = [list(r) for r in Bin]
    n = len(B)
    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))
    d = [0] * (n + 1)
    d[0] = 1
    lam = [[0] * n for _ in range(n)]
    def redi(k, l):
        if 2 * abs(lam[k][l]) > d[l + 1]:
            q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1])
            B[k] = [x - q * y for x, y in zip(B[k], B[l])]
            lam[k][l] -= q * d[l + 1]
            for i in range(l):
                lam[k][i] -= q * lam[l][i]
    def swapi(k, kmax):
        B[k], B[k - 1] = B[k - 1], B[k]
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
            u = dot(B[k], B[j])
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
    return B

def hunt(names, vals, tag, maxcoef=1 << 40):
    d = len(vals)
    K = MODW
    rows = []
    for i in range(d):
        r = [0] * d + [K * (vals[i] % MODW)]
        r[i] = 1
        rows.append(r)
    rows.append([0] * d + [K * MODW])
    red = lll(rows)
    found = []
    for r in red:
        if r[-1] == 0 and any(r[:d]):
            c = r[:d]
            if max(abs(x) for x in c) < maxcoef:
                if sum(ci * vi for ci, vi in zip(c, vals)) % MODW == 0:
                    found.append(c)
    print('[%s] dim %d: %s' % (tag, d,
          'RELATIONS: ' + '; '.join(
              str([(n_, c_) for n_, c_ in zip(names, c) if c_]) for c in found[:3])
          if found else 'EMPTY (coeffs < 2^40, mod 2^%d)' % WREL), flush=True)
    return found

# control: planted
planted = (11 * Wobj - 3 * D1[0] + 2 * LG) % MODW
hunt(['W', 'D1c0', 'LG', 'planted'], [Wobj, D1[0], LG, planted], 'control:planted')

# the swap hunt
names = (['W'] +
         ['D1c%d' % i for i in range(4)] + ['D2c%d' % i for i in range(4)] +
         ['16^T*D1c%d' % i for i in range(4)] + ['16^T*D2c%d' % i for i in range(4)] +
         ['LG', '16^T', '1'])
vals = ([Wobj] + D1 + D2 +
        [S16T * x % MODW for x in D1] + [S16T * x % MODW for x in D2] +
        [LG, S16T, 1])
hunt(names, vals, 'L2-swap')

# variant: swap with the COMPLEMENT length (A-8T small side) as dual step
Tc = (A - 8 * T)  # = 33, the residual small parameter
D1c = dual(lambda k: k)  # same
D3 = [0, 0, 0, 0]
for k in range(8):
    e = k % 8
    c, sg = e % 4, (-1) ** (e // 4)
    D3[c] = (D3[c] + sg * inv(A - Tc * k if (A - Tc * k) % 2 else A - Tc * k + 1)) % MODP
names2 = ['W'] + ['D3c%d' % i for i in range(4)] + ['LG', '16^T', '1']
hunt(names2, [Wobj] + D3 + [LG, S16T, 1], 'L2-swap-complement')
