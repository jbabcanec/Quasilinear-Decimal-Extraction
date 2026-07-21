"""AGM-law hunt: does the truncated 2-adic logarithm obey ANY algebraic
relation across dyadic truncation scales, beyond the known functional
equation  T_{2M}(z) + T_{2M}(-z) = T_M(z^2)  (exact, even half only)?

If an AGM/Landen-type depth-doubling iteration exists, the odd halves
O_M(z) = sum_{n<M odd} z^n/n must satisfy some low-height algebraic
(linear or quadratic) relation with objects at other scales / arguments /
the complete logs. LLL hunts for integer relations mod 2^W.

Baby atom: z = 6 in Q_2 (v(z)=1, non-degenerate: L(6) = -Log2(-5) != 0).
All objects are integers mod 2^W (every term z^n/n is 2-adically integral).

Hunts:
  0. harness control: recover the FE (1,1,-1) and a planted relation
  1. linear cross-scale: odd parts O_{2^s}, completes, 1
  2. quadratic (AGM-shape): pairwise products of T's and L's
  3. exp-side: E_s = exp(2 T_{2^s}(6)) family
Verdict: any verified nontrivial relation outside the FE-ideal = crack.
"""
from fractions import Fraction

W = 300          # relation checked mod 2^W
P = W + 48       # working precision
MODW = 1 << W
MODP = 1 << P

def v2int(n):
    return (n & -n).bit_length() - 1 if n else 10 ** 9

def T(M, z, P=P):
    """sum_{1<=n<M} z^n/n mod 2^P for even z (v(z)>=1). Exact 2-adic."""
    mod = 1 << (P + 24)
    zp, acc = 1, 0
    for n in range(1, M):
        zp = zp * z % mod
        e = v2int(n)
        o = n >> e
        acc += (zp >> e) * pow(o, -1, MODP)
    return acc % MODP % (1 << P)

def O(M, z):
    """odd part: sum_{n<M, n odd} z^n/n mod 2^P."""
    mod = 1 << (P + 24)
    zp, acc = 1, 0
    for n in range(1, M):
        zp = zp * z % mod
        if n % 2 == 1:
            acc += zp * pow(n, -1, MODP)
    return acc % MODP

def L(z, vz=1):
    """complete log -Log2(1-z): series to self-truncation."""
    M = (P + 32) // vz
    return T(M, z)

def expser(a):
    """exp(a) mod 2^W for v2(a) >= 2."""
    acc, term = 1, 1
    for k in range(1, W):
        term = term * a % MODP
        e = 0
        kk = k
        while kk % 2 == 0:
            kk //= 2; e += 1
        # divide by k: k = 2^e * kk ; term must keep integrality via v(a^k/k!)
        term = (term >> e) * pow(kk, -1, MODP) % MODP
        acc = (acc + term) % MODP
        if term == 0:
            break
    return acc % MODW

# ---------------- all-integer LLL (Cohen Alg. 2.6.7, no Fractions) ----------------
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

def hunt(names, vals, tag, maxcoef=1 << 60):
    """seek small integer c with sum c_i vals_i == 0 mod 2^W."""
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
                chk = sum(ci * vi for ci, vi in zip(c, vals)) % MODW
                if chk == 0:
                    found.append(c)
    print('[%s] dim %d: %s' % (tag, d,
          'RELATIONS: ' + '; '.join(str([(n, c_) for n, c_ in zip(names, c) if c_])
                                    for c in found[:3])
          if found else 'EMPTY (no relation, coeffs < 2^60, mod 2^300)'),
          flush=True)
    return found

# ---------------- objects ----------------
z = 6
T16, T32, T64, T128, T256 = (T(M, 6) for M in (16, 32, 64, 128, 256))
T16n, T32n = T(16, -6), T(32, -6)
Tq16, Tq32, Tq64 = (T(M, 36) for M in (16, 32, 64))
L6, L36 = L(6), L(36)
O16, O32, O64, O128, O256 = (O(M, 6) for M in (16, 32, 64, 128, 256))

# hunt 0a: FE control  T(32,6) + T(32,-6) - T(16,36) = 0
hunt(['T32(6)', 'T32(-6)', 'T16(36)'], [T32, T32n, Tq16], 'control:FE')
# hunt 0b: planted relation detection at full scale
planted = (7 * O32 - 5 * O64) % MODW
hunt(['O32', 'O64', 'planted'], [O32, O64, planted], 'control:planted')

# hunt 1: linear cross-scale on odd parts + completes + 1
hunt(['O16', 'O32', 'O64', 'O128', 'O256', 'L(6)', 'L(36)', '1'],
     [O16, O32, O64, O128, O256, L6, L36, 1], 'linear-cross-scale')

# hunt 2: quadratic (AGM-shape) -- products among the doubling chain
lin = [('T16', T16), ('T32', T32), ('T64', T64), ('Tq16', Tq16),
       ('Tq32', Tq32), ('L6', L6), ('L36', L36)]
quad = [('T16', T16), ('T32', T32), ('T64', T64), ('L6', L6)]
names, vals = ['1'], [1]
for n1, v1 in lin:
    names.append(n1); vals.append(v1)
for i in range(len(quad)):
    for j in range(i, len(quad)):
        names.append(quad[i][0] + '*' + quad[j][0])
        vals.append(quad[i][1] * quad[j][1] % MODW)
hunt(names, vals, 'quadratic-AGM', maxcoef=1 << 50)

# hunt 3: exp-side family
E = [expser(2 * T(1 << s, 6) % MODP) for s in (4, 5, 6, 7)]
inv25 = pow(25, -1, MODW)
hunt(['E16', 'E32', 'E64', 'E128', '1/25', '1'],
     E + [inv25, 1], 'exp-side')
