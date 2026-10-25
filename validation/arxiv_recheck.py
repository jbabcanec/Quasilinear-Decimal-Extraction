#!/usr/bin/env python3
"""
Independent re-verification of every checkable claim in decimal_extraction_paper.tex,
written directly from the statements in the paper (not from the existing validators).
Run:  python validation/arxiv_recheck.py
Every check prints PASS/FAIL; exit code 0 iff all pass.
"""
from fractions import Fraction
import math, random, sys

random.seed(20260805)
FAILURES = []

def report(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f"  -- {detail}" if detail else ""))
    if not ok:
        FAILURES.append((name, detail))

def v2(x):
    x = int(x)
    if x == 0:
        return None  # infinity
    return (x & -x).bit_length() - 1

def P0(K):
    p = 1
    for k in range(1, K + 1):
        p *= 8 * k + 1
    return p

def B(m, K):
    R = 1 << m
    p = 1
    for k in range(1, K + 1):
        p *= R + 8 * k + 1
    return p

def inv_mod2(a, m):
    """inverse of odd a mod 2^m"""
    return pow(a, -1, 1 << m)

# ----------------------------------------------------------------------
# 1. Dual squaring lemma: exact ratio identity (as rationals)
# ----------------------------------------------------------------------
ok = True
for (m, K) in [(3, 5), (5, 4), (6, 7), (8, 12), (10, 16), (12, 5)]:
    R = 1 << m
    lhs = Fraction(B(m, K) ** 2, P0(K) * B(m + 1, K))
    rhs = Fraction(1)
    for k in range(1, K + 1):
        mu = 8 * k + 1
        rhs *= 1 + Fraction(R * R, mu * (mu + 2 * R))
    if lhs != rhs:
        ok = False
        report("dual-square ratio identity", False, f"(m,K)=({m},{K})")
        break
report("dual-square ratio identity (exact rationals, 6 cases)", ok)

# ----------------------------------------------------------------------
# 2. Dual squaring valuation law: v2(B(m,K)^2 - P0*B(m+1,K)) == 2m + v2(K)
#    Paper claims m>=5, v2(K)<m.  Test m=5..24 x K=1..96, plus high-v2 K,
#    plus m=5 boundary explicitly, plus spot checks at larger m.
# ----------------------------------------------------------------------
ok = True; ntest = 0; worst = None
cases = [(m, K) for m in range(5, 25) for K in range(1, 97)]
cases += [(m, 1 << j) for m in (5, 6, 12, 20, 30) for j in range(0, 11) if j < m]
cases += [(40, 64), (40, 96), (60, 128), (30, 3 * 2 ** 8)]
for (m, K) in cases:
    if v2(K) >= m:
        continue
    val = v2(B(m, K) ** 2 - P0(K) * B(m + 1, K))
    ntest += 1
    if val != 2 * m + v2(K):
        ok = False; worst = (m, K, val)
report("dual-square valuation law = 2m+v2(K)", ok,
       f"{ntest} pairs incl. m=5 and v2(K) up to 10" + ("" if ok else f"; first fail {worst}"))

# where does it break below m=5?  (informational, not a paper claim)
small = []
for m in range(1, 5):
    bad = sum(1 for K in range(1, 33) if v2(K) < m and
              v2(B(m, K) ** 2 - P0(K) * B(m + 1, K)) != 2 * m + v2(K))
    small.append((m, bad))
print(f"       (info: failures for m=1..4 out of eligible K<=32: {small})")

# ----------------------------------------------------------------------
# 3. Jump law: v2(B(m)^(2^s) - P0^(2^s-1) B(m+s)) == 2m + (s-1) + v2(K)
# ----------------------------------------------------------------------
ok = True; ntest = 0
for m in (5, 6, 8, 10, 12):
    for K in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32):
        for s in (1, 2, 3, 4, 5):
            if v2(K) >= m:
                continue
            lhsv = v2(B(m, K) ** (1 << s) - P0(K) ** ((1 << s) - 1) * B(m + s, K))
            ntest += 1
            if lhsv != 2 * m + (s - 1) + v2(K):
                ok = False
                report("jump law", False, f"(m,K,s)=({m},{K},{s}): got {lhsv}")
report("jump valuation law = 2m+(s-1)+v2(K)", ok, f"{ntest} triples")

# full-width jump criterion: works at width 2(m+s) iff s <= v2(K)-1
ok = True
for m in (8, 10):
    for K in (8, 16, 32, 64):
        for s in (1, 2, 3, 4, 5, 6):
            lhs = (2 * m + (s - 1) + v2(K)) >= 2 * (m + s)
            rhs = (s <= v2(K) - 1)
            if lhs != rhs:
                ok = False
report("full-width jump iff s <= v2(K)-1 (algebraic consequence)", ok)

# ----------------------------------------------------------------------
# 4. Faulhaber/von-Staudt step used in the lemma's proof:
#    Q2 = sum mu_k^{-2} (2-adically) satisfies v2(Q2 - K) > v2(K), so v2(Q2)=v2(K).
#    Also C = sum 1/(mu(mu+2R)) == Q2 mod 2^{m+1}.
# ----------------------------------------------------------------------
ok = True
W = 90
for K in [1, 2, 3, 4, 8, 12, 16, 31, 32, 64, 100, 128, 256, 512, 1024, 96, 160]:
    q2 = sum(inv_mod2(8 * k + 1, W) ** 2 for k in range(1, K + 1)) % (1 << W)
    d = (q2 - K) % (1 << W)
    if d != 0 and v2(d) <= v2(K):
        ok = False
        report("v2(Q2-K) > v2(K)", False, f"K={K}: v2={v2(d)}")
report("v2(Q2 - K) > v2(K), hence v2(Q2)=v2(K)  (17 K-values to 1024)", ok)

ok = True
for (m, K) in [(8, 5), (10, 12), (16, 33), (20, 64)]:
    R = 1 << m
    W2 = m + 1
    C = sum(inv_mod2((8 * k + 1) * (8 * k + 1 + 2 * R) % (1 << (W2 + 8)), W2 + 8)
            for k in range(1, K + 1)) % (1 << W2)
    q2 = sum(inv_mod2(8 * k + 1, W2) ** 2 for k in range(1, K + 1)) % (1 << W2)
    if C != q2 % (1 << W2):
        ok = False
report("C == Q2 mod 2^(m+1)", ok)

# ----------------------------------------------------------------------
# 5. Prop E-equiv (ii): B(m,K) == P0 (1 + 2^m H) mod 2^(2m)
# ----------------------------------------------------------------------
ok = True
for (m, K) in [(16, 20), (37, 100), (64, 64), (128, 101), (256, 300)]:
    H = sum(inv_mod2(8 * k + 1, m) for k in range(1, K + 1)) % (1 << m)
    lhs = B(m, K) % (1 << (2 * m))
    rhs = (P0(K) * (1 + (H << m))) % (1 << (2 * m))
    if lhs != rhs:
        ok = False
        report("E-equiv (ii)", False, f"(m,K)=({m},{K})")
report("E-equiv dual product: B == P0(1+2^m H) mod 2^(2m)  (5 sizes)", ok)

# ----------------------------------------------------------------------
# 6. Geometric expansion H = sum_j (-8)^j S_j(K), convergence 3 bits/term
# ----------------------------------------------------------------------
ok = True
m, K = 60, 50
H = sum(inv_mod2(8 * k + 1, m) for k in range(1, K + 1)) % (1 << m)
acc = 0
Jmax = m // 3 + 2
for j in range(Jmax + 1):
    Sj = sum(pow(k, j, 1 << m) for k in range(1, K + 1))
    acc = (acc + pow(-8, j, 1 << m) * Sj) % (1 << m)
    if j >= 1:
        d = (acc - H) % (1 << m)
        if d != 0 and v2(d) < 3 * (j + 1) - 4:  # generous slack; nominal rate 3j
            ok = False
report("geometric expansion H = sum (-8)^j S_j converges at ~3 bits/term", ok,
       f"final error v2 >= {m if (acc-H)%(1<<m)==0 else v2((acc-H)%(1<<m))} at j={Jmax}, m={m}")

# ----------------------------------------------------------------------
# 7. 2-adic Log: isometry + homomorphism on 1-units (Lemma flatten),
#    and Log(B/P0) = sum ((-1)^(j+1)/j) 2^(wj) Q_j (Theorem conservation display)
# ----------------------------------------------------------------------
def log2adic(u, W, guard=16):
    """2-adic Log of u = 1 + t, v2(t)>=2, computed mod 2^W (uses W+guard working bits)."""
    WW = W + guard
    t = (u - 1) % (1 << WW)
    if t == 0:
        return 0
    assert v2(t) >= 2
    acc = 0
    tj = t
    j = 1
    while tj != 0 and j * 2 <= WW + 40:
        term = tj * inv_mod2(j >> v2(j), WW)
        e = v2(j)
        # tj has v2 >= 2j > e, exact shift
        term >>= e
        if j % 2 == 1:
            acc = (acc + term) % (1 << WW)
        else:
            acc = (acc - term) % (1 << WW)
        j += 1
        tj = (tj * t) % (1 << WW)
    return acc % (1 << W)

W = 200
ok = True
for _ in range(60):
    a = (random.getrandbits(W - 4) << 3) | 0  # 8Z2 element
    b = (random.getrandbits(W - 4) << 3) | 0
    u = (1 + a) % (1 << W)
    v = (1 + b) % (1 << W)
    lu, lv = log2adic(u, W - 16), log2adic(v, W - 16)
    if a and v2(lu) != v2(a):
        ok = False; report("flatten isometry", False, f"v2(Log u)={v2(lu)} vs v2(u-1)={v2(a)}")
    luv = log2adic(u * v % (1 << W), W - 16)
    if (lu + lv - luv) % (1 << (W - 20)) != 0:
        ok = False; report("flatten homomorphism", False, "Log(uv)!=Logu+Logv")
report("flatten: Log isometry v2(Log u)=v2(u-1) and homomorphism (60 random pairs, W=200)", ok)

ok = True
for (w, K) in [(8, 12), (12, 30), (16, 20)]:
    W = 6 * w
    b_over_p = B(w, K) * inv_mod2(P0(K), W + 24) % (1 << (W + 24))
    lhs = log2adic(b_over_p, W)
    rhs = 0
    j = 1
    while w * j < W + 24:
        Qj = sum(pow(inv_mod2(8 * k + 1, W + 24), j, 1 << (W + 24)) for k in range(1, K + 1))
        term = (Qj << (w * j)) * inv_mod2(j >> v2(j), W + 24) >> v2(j)
        rhs = (rhs + (term if j % 2 == 1 else -term)) % (1 << (W + 24))
        j += 1
    if (lhs - rhs) % (1 << W) != 0:
        ok = False; report("Log(B/P0) tower", False, f"(w,K)=({w},{K})")
report("Log(B/P0) = sum ((-1)^(j+1)/j) 2^(wj) Q_j  (3 sizes)", ok)

# ----------------------------------------------------------------------
# 8. Scale conservation: numeric adversarial check of the certified bound
#    max_i (p_i - w_i), with worst-case error draws.
# ----------------------------------------------------------------------
ok = True
Wbig = 260
K = 20
profiles = [((8, 30), (12, 40), (16, 44)), ((10, 25), (20, 55)), ((6, 40), (7, 20), (9, 33))]
for prof in profiles:
    bound = max(p - w for (w, p) in prof)
    for trial in range(200):
        aa = [random.randint(-8, 8) for _ in prof]
        c = sum(a * (1 << w) for a, (w, p) in zip(aa, prof))
        if c == 0:
            continue
        # worst-case certified error: some draw with v2(a_i eps_i) exactly min(p_i+v2(a_i))
        efloor = min(p + (v2(a) if a else 10 ** 9) for a, (w, p) in zip(aa, prof))
        cert = efloor - v2(c)
        if cert > bound:
            ok = False
            report("scale conservation bound", False, f"profile {prof}, a={aa}: cert {cert} > bound {bound}")
report("scale conservation: certified precision never exceeds max(p_i - w_i) (600 combos)", ok)

# ----------------------------------------------------------------------
# 9. Band identity, h-recovery, Bezout split, consolidation (core pipeline)
# ----------------------------------------------------------------------
ok = True
for _ in range(300):
    X = random.getrandbits(400) | 1
    mu = random.randrange(3, 1 << 16, 2)
    J = random.randint(1, 380)
    Hfl = X >> J
    h = Hfl % mu
    u = Fraction(X % (1 << J), 1 << J)
    lhs = Fraction(X, mu << J)
    lhs = lhs - (lhs.numerator // lhs.denominator)
    if lhs != Fraction(h + u.numerator / Fraction(u.denominator), mu) and lhs != (h + u) / mu:
        ok = False
report("band identity {X/(mu 2^J)} = (h+u)/mu  (300 random)", ok)

ok = True
for _ in range(300):
    X = random.getrandbits(400) | 1
    mu = random.randrange(3, 1 << 16, 2)
    J = random.randint(1, 380)
    R = X % mu
    Yt = (X % (1 << J)) % mu
    h = ((R - Yt) * pow(pow(2, J, mu), -1, mu)) % mu
    if h != (X >> J) % mu:
        ok = False
report("h-recovery h = (R - Ytilde) 2^(-J) mod mu  (300 random)", ok)

ok = True
for _ in range(300):
    X = random.getrandbits(300) | 1
    mu = random.randrange(3, 1 << 12, 2)
    J = random.randint(2, 280)
    u = pow(pow(2, J, mu), -1, mu)
    w = pow(mu, -1, 1 << J)
    lhs = Fraction(X, mu << J)
    lhs -= lhs.numerator // lhs.denominator
    rhs = Fraction((X % mu) * u % mu, mu) + Fraction(X * w % (1 << J), 1 << J)
    rhs -= rhs.numerator // rhs.denominator
    if lhs != rhs:
        ok = False
report("Bezout split (300 random)", ok)

ok = True
for _ in range(40):
    X = random.getrandbits(200) | 1
    terms = []
    for _t in range(25):
        mu = random.randrange(3, 1 << 10, 2)
        J = random.randint(2, 150)
        sg = random.choice([1, -1])
        terms.append((sg, mu, J))
    Jmax = max(J for _, _, J in terms)
    Vv = 0
    for sg, mu, J in terms:
        Vv = (Vv + sg * ((1 << (Jmax - J)) * pow(mu, -1, 1 << Jmax))) % (1 << Jmax)
    lhs = Fraction(0)
    for sg, mu, J in terms:
        f = Fraction(X * pow(mu, -1, 1 << J) % (1 << J), 1 << J)
        lhs += sg * f
    lhs -= lhs.numerator // lhs.denominator
    rhs = Fraction(X * Vv % (1 << Jmax), 1 << Jmax)
    if (lhs - rhs).numerator % 1 != 0 or (lhs - rhs) != Fraction((lhs - rhs).numerator // (lhs - rhs).denominator):
        # compare mod 1
        d = lhs - rhs
        if d != int(d):
            ok = False
report("dyadic consolidation sum{X mu^-1 2^-J} == {XV/2^Jmax} mod 1  (40 random multisets)", ok)

# ----------------------------------------------------------------------
# 10. Euler starters lemma
# ----------------------------------------------------------------------
def is_prime(n):
    if n < 2: return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0: return n == p
    d = n - 1; r = 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

ok = True
count = 0
q = 3
while q < 20000:
    if is_prime(q):
        # minimal k0 with q | 8k0+1
        k0 = (pow(8, -1, q) * (q - 1) * 1) % q  # 8k0 == -1 mod q -> k0 = -inv8 mod q
        k0 = (-pow(8, -1, q)) % q
        s = (8 * k0 + 1) // q
        assert (8 * k0 + 1) % q == 0
        if not (s <= 7 and s % 2 == 1):
            ok = False; break
        leg = pow(2, (q - 1) // 2, q)  # = (2/q) mod q
        pred = pow(leg, s, q) * pow(2, (s - 1) // 2, q) % q
        if pow(2, 4 * k0, q) != pred:
            ok = False
            report("Euler starters", False, f"q={q}")
        count += 1
    q += 2
report("Euler starters lemma (all odd primes < 20000)", ok, f"{count} primes")

# ----------------------------------------------------------------------
# 11. Power-sum congruence: sum (1+8k)^(2^r) == K + 2^r sum Log(1+8k) mod 2^(2r+5)
# ----------------------------------------------------------------------
ok = True
for (K, r) in [(20, 20), (50, 30), (7, 12)]:
    Wc = 2 * r + 5
    lhs = sum(pow(1 + 8 * k, 1 << r, 1 << Wc) for k in range(1, K + 1)) % (1 << Wc)
    ls = 0
    for k in range(1, K + 1):
        ls = (ls + log2adic(1 + 8 * k, Wc)) % (1 << Wc)
    rhs = (K + (ls << r)) % (1 << Wc)
    if lhs != rhs:
        ok = False
        report("power-sum congruence", False, f"(K,r)=({K},{r})")
report("power-sum congruence mod 2^(2r+5)  (3 sizes)", ok)

# ----------------------------------------------------------------------
# 12. 2-adic Euler-Maclaurin (Lemma em), small case, exact von-Staudt handling
# ----------------------------------------------------------------------
def bernoulli_list(n):
    """B_0..B_n as Fractions, B1 = -1/2 convention"""
    A = [Fraction(0)] * (n + 1)
    B_ = []
    for mm in range(n + 1):
        A[mm] = Fraction(1, mm + 1)
        for j in range(mm, 0, -1):
            A[j - 1] = j * (A[j - 1] - A[j])
        B_.append(A[0])
    return B_

def frac_mod2(fr, W):
    """rational with v2 >= 0 reduced mod 2^W"""
    num, den = fr.numerator, fr.denominator
    e = v2(den)
    assert e == 0 or v2(num) >= e, "negative valuation"
    if e:
        num >>= 0
        # divide num by 2^e exactly if possible; else shift into num's v2
        assert v2(num) >= e
        num = num >> e
        den = den >> e
    return num * inv_mod2(den, W) % (1 << W)

ok = True
K, W = 12, 70
lhs = 0
for k in range(1, K + 1):
    lhs = (lhs + log2adic(1 + 8 * k, W)) % (1 << W)
# F(x) = ((1+8x)(f(x)-1)+1)/8 with f = Log(1+8x)
fK = log2adic(1 + 8 * K, W + 6)
FK = ((1 + 8 * K) * (fK - 1) + 1) % (1 << (W + 6))
assert FK % 8 == 0
FK = (FK // 8) % (1 << W)
half_fK = fK * inv_mod2(1, W)  # (1/2) f(K): f(K) has v2>=3 so exact shift
assert v2(fK % (1 << (W + 6))) >= 1
rhs = (FK + ((fK % (1 << (W + 6))) >> 1)) % (1 << W)
Bs = bernoulli_list(40)
for i in range(1, 15):
    # f^(2i-1)(x) = (-1)^(2i-2) (2i-2)! 8^(2i-1) (1+8x)^(-(2i-1))
    r_ = 2 * i - 1
    coef = Fraction(Bs[2 * i], math.factorial(2 * i)) * math.factorial(r_ - 1) * (8 ** r_)
    # term = coef * ((1+8K)^-r - 1^-r)... f^(r)(0) uses (1+0)^{-r} = 1
    # 2-adic evaluation:
    num, den = coef.numerator, coef.denominator
    e = v2(den)
    valnum = num * (pow(inv_mod2(1 + 8 * K, W + 20), r_, 1 << (W + 20)) - 1)
    if valnum == 0:
        continue
    assert v2(valnum) >= e
    term = (valnum >> e) * inv_mod2(den >> e, W) % (1 << W)
    rhs = (rhs + term) % (1 << W)
ok = (lhs == rhs)
report("2-adic Euler-Maclaurin identity (K=12, 70 bits)", ok,
       "" if ok else f"diff v2 = {v2((lhs-rhs)%(1<<W))}")

# ----------------------------------------------------------------------
# 13. Prop atom: Pi'(0)/Pi(0) = sum 1/(8k+1)  (exact rationals, small K)
# ----------------------------------------------------------------------
ok = True
K = 9
mus = [8 * k + 1 for k in range(1, K + 1)]
Pi0 = 1
for mu in mus: Pi0 *= mu
dPi = sum(Pi0 // mu for mu in mus)
ok = Fraction(dPi, Pi0) == sum(Fraction(1, mu) for mu in mus)
report("Pi'(0)/Pi(0) = sum 1/(8k+1)  (K=9 exact)", ok)

# ----------------------------------------------------------------------
# 14. B(m,K) = P0(K + 2^(m-3)) / P0(2^(m-3))  (offset-segment identity)
# ----------------------------------------------------------------------
ok = True
for (m, K) in [(6, 10), (8, 5), (10, 20)]:
    if B(m, K) * P0(1 << (m - 3)) != P0(K + (1 << (m - 3))):
        ok = False
report("B(m,K) is the 1-mod-8 product segment at offset 2^(m-3)", ok)

# ----------------------------------------------------------------------
# 15. Normal form: numeric check that the four-family decomposition sums to pi
#     (Lemma normal) via high-precision floats
# ----------------------------------------------------------------------
import mpmath as mp
mp.mp.dps = 120
n = 40
tot = mp.mpf(0)
X5 = mp.mpf(5) ** n
for k in range(0, 400):
    A = n - 4 * k
    tot += mp.mpf(2) ** (A + 2) * X5 / (8 * k + 1)
    tot -= mp.mpf(2) ** (A - 1) * X5 / (2 * k + 1)
    tot -= mp.mpf(2) ** (A) * X5 / (8 * k + 5)
    tot -= mp.mpf(2) ** (A - 1) * X5 / (4 * k + 3)
ok = abs(tot - mp.mpf(10) ** n * mp.pi) < mp.mpf(10) ** (n - 100 + 60)
report("normal form sums to 10^n pi  (n=40, 400 terms, 120 dps)", ok,
       f"residual ~ {mp.nstr(abs(tot - mp.mpf(10)**n*mp.pi)/mp.mpf(10)**n, 3)} of pi units")

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILURES:")
    for f in FAILURES:
        print("  ", f)
    sys.exit(1)
print("ALL CHECKS PASS")
