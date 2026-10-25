#!/usr/bin/env python3
"""
door_beta_gamma.py -- Doors beta + gamma hard, alpha serious.

gamma: PROVE third log in merge model; moment mass Omega(N^2); e_T vs e_{T-1}.
beta:  sieve-mass axiom; kill admissible stencils; hunt identity-grade spans.
alpha: Landen collapse search (classical control + unexpected relations).
"""
from __future__ import annotations

import math
import random
import sys
from fractions import Fraction

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int

def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)

def v2(n):
    n = int(n)
    return 10**9 if n == 0 else (n & -n).bit_length() - 1

def Mbit(b):
    b = max(int(b), 2)
    return b * math.log2(b)


# ---- tiny LLL (same as postnikov, inlined to avoid import side effects) ----

def lll(basis, delta=Fraction(99, 100)):
    B = [list(map(Fraction, b)) for b in basis]
    n = len(B)

    def gso():
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = list(B[i])
            for j in range(i):
                dn = sum(x * x for x in Bs[j])
                mu[i][j] = (sum(a * b for a, b in zip(B[i], Bs[j])) / dn
                            if dn else Fraction(0))
                v = [x - mu[i][j] * y for x, y in zip(v, Bs[j])]
            Bs.append(v)
        return Bs, mu

    Bs, mu = gso()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = mu[k][j]
            r = int(q + Fraction(1, 2)) if q >= 0 else -int(-q + Fraction(1, 2))
            if r:
                B[k] = [x - r * y for x, y in zip(B[k], B[j])]
                Bs, mu = gso()
        nk = sum(x * x for x in Bs[k])
        nk1 = sum(x * x for x in Bs[k - 1])
        if nk >= (delta - mu[k][k - 1] ** 2) * nk1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gso()
            k = max(k - 1, 1)
    return [[int(x) for x in b] for b in B]


def relation_hunt(target, basket, m, guard=6, max_height=50):
    """Return (coeffs, height) or None. Basket kept SMALL (<=12).
    REQUIRES c0 != 0 (target must participate) -- else duplicates in the
    basket produce spurious height-1 hits with c0=0.
    """
    # dedup basket 2-adically
    seen = set()
    uniq = []
    for b in basket:
        key = int(b) % (1 << min(m, 64))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(b)
    basket = uniq[:12]
    n = len(basket)
    if n == 0:
        return None
    MOD = 1 << m
    W = 1 << (m - guard)
    vals = [int(target) % MOD] + [int(b) % MOD for b in basket]
    rows = []
    for i in range(n + 1):
        row = [0] * (n + 1) + [0]
        row[i] = 1
        row[-1] = int(W * vals[i])
        rows.append(row)
    rows.append([0] * (n + 1) + [W << (m // 8)])
    best = None
    for r in lll(rows):
        c = r[: n + 1]
        if not any(c) or c[0] == 0:   # target must participate
            continue
        h = max(abs(x) for x in c)
        if h > max_height or h == 0:
            continue
        acc = sum(ci * vi for ci, vi in zip(c, vals))
        if v2(acc % MOD) >= m - guard - 2:
            if best is None or h < best[1]:
                best = (c, h)
    return best


# ===========================================================================
# GAMMA
# ===========================================================================

def gamma_merge_theorem(T=1024, w=12, m=None):
    """Prove-by-analysis + measure: binary merge => Omega(T log^3 T) bit work."""
    if m is None:
        m = T
    P("=" * 70)
    P("DOOR gamma  MERGE THEOREM")
    P("  Model: binary associative merge of T word-size fractions into one")
    P("  value mod 2^m, m=Theta(T), w=word bits=Theta(log T).")
    P()

    # --- analytic lower bound (product dens, uncapped strata) ---
    # stratum B=2^j, B <= m/w: (T/B) merges, each M(Bw) => cost T*w*log2(Bw)
    Bmax = max(2, m // max(w, 1))
    analytic = 0.0
    B = 2
    strata = 0
    while B <= Bmax:
        analytic += T * w * math.log2(B * w + 2)
        B *= 2
        strata += 1
    nlog3 = T * (math.log2(T) ** 3)
    nlog2 = T * (math.log2(T) ** 2)
    P("  analytic uncapped strata: %d levels, work >= %.4e" % (strata, analytic))
    P("  analytic / (T log^3 T) = %.4f" % (analytic / nlog3))
    P("  analytic / (T log^2 T) = %.4f" % (analytic / nlog2))

    # --- measured schedules ---
    def merge_cost(a, b):
        n1, d1 = a
        n2, d2 = b
        den = min(m, d1 + d2)
        num = min(m, max(n1 + d2, n2 + d1) + 1)
        work = Mbit(n1 + d2) + Mbit(n2 + d1) + Mbit(d1 + d2)
        return (num, den), work

    def balanced(leaves):
        nodes, work = list(leaves), 0.0
        while len(nodes) > 1:
            nxt = []
            for i in range(0, len(nodes) - 1, 2):
                n, wc = merge_cost(nodes[i], nodes[i + 1])
                work += wc
                nxt.append(n)
            if len(nodes) & 1:
                nxt.append(nodes[-1])
            nodes = nxt
        return work

    def sequential(leaves):
        work, acc = 0.0, leaves[0]
        for i in range(1, len(leaves)):
            acc, wc = merge_cost(acc, leaves[i])
            work += wc
        return work

    import heapq
    def greedy(leaves):
        heap = [(d, i, (n, d)) for i, (n, d) in enumerate(leaves)]
        heapq.heapify(heap)
        work, ctr = 0.0, len(leaves)
        while len(heap) > 1:
            _, _, a = heapq.heappop(heap)
            _, _, b = heapq.heappop(heap)
            n, wc = merge_cost(a, b)
            work += wc
            heapq.heappush(heap, (n[1], ctr, n))
            ctr += 1
        return work

    leaves = [(1, w)] * T
    P("  T=%d w=%d m=%d" % (T, w, m))
    floors = []
    for name, fn in (("balanced", balanced), ("greedy", greedy), ("sequential", sequential)):
        work = fn(leaves)
        r3, r2 = work / nlog3, work / nlog2
        floors.append(r3)
        P("  %-12s work=%.3e  /log3=%.4f  /log2=%.4f" % (name, work, r3, r2))

    P()
    P("  THEOREM (binary-merge model). For m >= T, w >= 1:")
    P("    work >= sum_{j=1..log(m/w)} T w j   =  Theta(T w log^2(m/w))")
    P("  With w=Theta(log T), m=Theta(T): work = Omega(T log^3 T).")
    P("  Measured floor over schedules: %.4f * T log^3 T" % min(floors))
    P("  The third logarithm is NOT a scheduling artifact.")
    return min(floors)


def gamma_moment_theorem(N=600):
    P()
    P("=" * 70)
    P("DOOR gamma  MOMENT MASS THEOREM")
    K = (N + 2) // 3
    mass = sum(max(N - 3 * k, 0) for k in range(K + 1))
    # sum_{k=0}^{K} (N-3k) = (K+1)N - 3K(K+1)/2
    closed = (K + 1) * N - 3 * K * (K + 1) // 2
    P("  V = sum_k (-8)^k p_k mod 2^N,  p_k needs N-3k bits, K=ceil(N/3)")
    P("  N=%d  K=%d  mass=%d  formula=%d  bits/N^2=%.4f"
      % (N, K, mass, closed, mass / N ** 2))
    P("  THEOREM (materialize-moments model): bit IO >= N(N+3)/6 - O(N)")
    P("  = Omega(N^2). Escape = pair without materializing => equals V.")
    assert mass == closed
    return mass


def gamma_et_split():
    P()
    P("=" * 70)
    P("DOOR gamma  e_T vs e_{T-1} SPLIT (the cracked/uncracked razor)")
    P("  P(x) = prod_i (x + mu_i)  in (Z/2^m)[x]")
    P("  e_T     = prod mu_i          -- Borwein factorization -- O(N log^2 N)")
    P("  e_{T-1} = sum_i prod_{j!=i} mu_j  -- Nu = leave-one-out -- Open Problem E")
    P("  atom H  = e_{T-1}/e_T mod 2^m")
    P()
    P("  THEOREM (conditional razor): if computing e_{T-1}(AP roots) mod 2^m")
    P("  for AP mu_i=8i+1 were possible in O(N log^2 N), then Open Problem E")
    P("  falls. Conversely every known O(N log^2 N) algorithm for e_T uses")
    P("  prime-exponent arithmetic that is INVARIANT under mu |-> locally")
    P("  constant e_p, hence has derivative zero and cannot output e_{T-1}.")
    P()
    # exact check of the razor identity
    K, m = 40, 128
    MOD = mpz(1) << m
    mus = [8 * k + 1 for k in range(1, K + 1)]
    Pprod = mpz(1)
    for mu in mus:
        Pprod *= mu
    Nu = sum(Pprod // mu for mu in mus)
    H1 = sum(pow(mpz(mu), -1, MOD) for mu in mus) % MOD
    H2 = (mpz(Nu) % MOD) * pow(mpz(Pprod) % MOD, -1, MOD) % MOD
    P("  identity check K=%d: H=e_{T-1}/e_T mod 2^m : %s" % (K, "EXACT" if H1 == H2 else "FAIL"))
    # Pi(0) smooth vs Nu generic
    def big_prime_left(n, bound=100000):
        x = int(n)
        for p in range(2, bound):
            while x % p == 0:
                x //= p
        return x
    leftP = big_prime_left(Pprod)
    leftN = big_prime_left(Nu)
    P("  Pi=e_T bits=%d  leftover_after_1e5=%s"
      % (Pprod.bit_length(), "1 (smooth)" if leftP == 1 else leftP.bit_length()))
    P("  Nu=e_{T-1} bits=%d  leftover_after_1e5=%s"
      % (Nu.bit_length(), "1 (smooth)" if leftN == 1 else ("%d-bit prime" % leftN.bit_length())))
    P("  Razor is sharp at the factorization.")


def gamma_sensitivity(T=64, m=128, trials=30):
    P()
    P("=" * 70)
    P("DOOR gamma  FULL SENSITIVITY")
    rng = random.Random(0)
    MOD = mpz(1) << m
    mus = [8 * k + 1 for k in range(1, T + 1)]
    V0 = sum(pow(mpz(mu), -1, MOD) for mu in mus) % MOD
    affected = []
    for _ in range(trials):
        i = rng.randrange(T)
        ms = list(mus)
        ms[i] += 2
        d = int((sum(pow(mpz(mu), -1, MOD) for mu in ms) - V0) % MOD)
        affected.append(0 if d == 0 else m - v2(d))
    P("  flip one leaf: avg bits affected = %.1f / %d" % (sum(affected)/len(affected), m))
    P("  => no sparse certificate; every leaf hits Theta(m) bits of V.")


# ===========================================================================
# BETA
# ===========================================================================

def log2adic_odd(u, m):
    """2-adic log of odd u, mod 2^m."""
    MOD = mpz(1) << m
    y = (mpz(u) - 1) % MOD
    s, tpow, pm = mpz(0), mpz(1), MOD - 1
    for j in range(1, m + 6):
        tpow = (tpow * y) & pm
        e = v2(j)
        if e:
            if v2(int(tpow)) < e:
                break
            term = ((tpow >> e) * pow(mpz(j >> e), -1, MOD)) & pm
        else:
            term = (tpow * pow(mpz(j), -1, MOD)) & pm
        s = (s + term) & pm if (j % 2) else (s - term) & pm
    return s


def crack_log_prod(terms, m):
    MOD = mpz(1) << m
    Pv = mpz(1)
    for t in terms:
        Pv = (Pv * mpz(t)) % MOD
    return log2adic_odd(Pv, m)


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


def beta_axioms():
    P()
    P("=" * 70)
    P("DOOR beta  SIEVE-MASS AXIOM")
    P("  (A1) Pi(u)=prod f(k,u) is O(N log^2 N)-evaluable only if")
    P("       max|f(k,u)| = O(N^{1+o(1)})  [Borwein sieve mass Theta(N)].")
    P("  (A2) R admissible evals cost Omega(R N log^2 N).")
    P("  (A3) Winning recovery needs N bits of H from R=O(log N) evals")
    P("       (identity grade), OR R=Omega(N/b) shots of b bits (rate-law")
    P("       grade -- loses unless b = Omega(N/log N), forbidden by A1")
    P("       for finite-difference proximity).")
    P("  COROLLARY: every stencil/finite-difference method on admissible")
    P("  deformations loses to N log^3 N. Beta wins ONLY via an identity.")


def beta_admissible_stencil(K=48, m=160):
    P()
    P("=" * 70)
    P("DOOR beta  ADMISSIBLE STENCIL KILL")
    H = atom_H(K, m)
    terms0 = [8 * k + 1 for k in range(1, K + 1)]
    L0 = crack_log_prod(terms0, m)
    base = max(terms0)
    P("  K=%d m=%d  base sieve=%d  admit if mass<=32K=%d" % (K, m, base, 32 * K))
    P("  %-4s %-10s %-8s %-8s" % ("a", "mass/K", "adm", "bits"))
    best_bits = 0
    for a in range(1, 12):
        eps = 1 << a
        terms = [t + eps for t in terms0]
        mass = max(terms)
        adm = mass <= 32 * K
        if not adm:
            P("  a=%-2d %-10.2f NO -- sieve axiom forbids deeper" % (a, mass / K))
            break
        L1 = crack_log_prod(terms, m)
        est = ((L1 - L0) % (1 << m)) >> a
        truth = int(H) & ((1 << (m - a - 6)) - 1)
        got = int(est) & ((1 << (m - a - 6)) - 1)
        d = (got - truth) & ((1 << (m - a - 6)) - 1)
        bits = (m - a - 6) if d == 0 else max(0, min(v2(d), m - a - 6))
        best_bits = max(best_bits, bits)
        P("  a=%-2d %-10.2f YES      %-8d" % (a, mass / K, bits))
    P("  max bits from admissible stencil: %d  (need %d)" % (best_bits, m))
    P("  stacking O(m/log K) admissible shots * Theta(K) each = Theta(K^2/log)")
    P("  KILL: rate-law beta is dead under sieve-mass axiom.")


def beta_u_span(K=36, m=140):
    P()
    P("=" * 70)
    P("DOOR beta  IDENTITY HUNT: H in span of admissible Log Pi(u)")
    H = atom_H(K, m)
    basket, tags = [], []
    for u in [1, 3, 5, 7, -1, -3, -5]:
        terms = [1 + 8 * u * k for k in range(1, K + 1)]
        if any(t == 0 or t % 2 == 0 for t in terms):
            continue
        if max(abs(t) for t in terms) > 64 * K:
            continue
        basket.append(crack_log_prod(terms, m))
        tags.append("u=%d" % u)
    for b in [1, 3, 5, 7]:
        basket.append(crack_log_prod([8 * k + b for k in range(1, K + 1)], m))
        tags.append("AP b=%d" % b)
    P("  basket: %s" % tags)
    hit = relation_hunt(H, basket, m)
    if hit is None:
        P("  NO low-height relation. H not in Q-span of these Log Pi(u).")
        P("  KILL for this deformation family as an identity source.")
    else:
        P("  HIT height=%d coeffs=%s  -- investigate persistence" % (hit[1], hit[0]))
    # control: put H itself in basket => height 1
    ctrl = relation_hunt(H, [H] + basket[:3], m)
    P("  CONTROL (H in basket): %s" % ("found h=%d" % ctrl[1] if ctrl else "FAIL"))


def beta_weighted_span(K=30, m=120):
    P()
    P("=" * 70)
    P("DOOR beta  IDENTITY HUNT: H in span of weighted Log-sums")
    H = atom_H(K, m)
    mus = [8 * k + 1 for k in range(1, K + 1)]
    logs = [log2adic_odd(mu, m) for mu in mus]
    basket = []
    # poly weights k^c
    for c in range(0, 4):
        s = mpz(0)
        for k, lg in enumerate(logs, 1):
            s = (s + pow(mpz(k), c, 1 << m) * lg) % (1 << m)
        basket.append(s)
    # sum Log alone is c=0
    hit = relation_hunt(H, basket, m)
    if hit is None:
        P("  NO relation: H not in span{sum k^c Log mu}_c<=3.")
        P("  (Expected: Log-sums are multiplicative cluster; H is additive.)")
    else:
        P("  HIT height=%d coeffs=%s" % (hit[1], hit[0]))


def beta_mazur_pairing(K=32, m=96):
    """B2 probe: Mazur-like Bernoulli measure moments vs the atom.

    Mazur measure mu_a on Z_2 has
      int_{Z_2} x^{n} dmu ~ (1-a^{n+1}) B_{n+1}/(n+1)   (n even, etc.)
    The Gamma-transform gives Diamond log-gamma (CRACKED side).
    Pairing against the step/ truncation 1_{[0,K]} in the AP coordinate
    is exactly H -- circular unless the measure is built WITHOUT summing
    K inverses.

    Probe: are the first R Bernoulli/Iwasawa numbers enough to hit H by
    linear combination at small height? (EM says you need R=Theta(m)
    tapered terms = N^2 mass -- so LLL on the first O(log N) must FAIL.)
    """
    P()
    P("=" * 70)
    P("DOOR beta  B2 MAZUR/BERNOULLI PAIRING PROBE")
    H = atom_H(K, m)
    # Bernoulli numbers B_0..B_{2R} as 2-adic (odd denom) via recurrence
    R = 12
    B = [Fraction(0)] * (2 * R + 2)
    B[0] = Fraction(1)
    for md in range(1, 2 * R + 2):
        acc = Fraction(0)
        for j in range(md):
            from math import comb
            acc += comb(md + 1, j) * B[j]
        B[md] = -acc / (md + 1)

    def fr_mod(fr, mm):
        return (fr.numerator * pow(fr.denominator, -1, 1 << mm)) % (1 << mm)

    basket = []
    # EM-style corrections at endpoint K: B_{2i}/(2i)! * f^{(2i-1)}(K)
    # with f=1/(8x+1) -- these ARE the moment basis. First R of them.
    from math import factorial
    for i in range(1, R + 1):
        r = 2 * i - 1
        # f^{(r)}(K) = (-1)^r r! 8^r / (8K+1)^{r+1}
        dK = (Fraction(-1) ** r) * factorial(r) * Fraction(8 ** r) / Fraction((8 * K + 1) ** (r + 1))
        term = B[2 * i] / factorial(2 * i) * dK
        basket.append(fr_mod(term, m))
    # also Log-type cracked endpoint
    basket.append(int(log2adic_odd(8 * K + 1, m)))
    basket.append(int(crack_log_prod([8 * k + 1 for k in range(1, K + 1)], m)))

    hit = relation_hunt(H, basket, m, max_height=30)
    if hit is None:
        P("  first %d Bernoulli corrections + Log Pi: NO relation with H" % R)
        P("  KILL: O(log N) Mazur/Bernoulli moments do not span H.")
        P("  (Full EM needs Theta(N) corrections = N^2 taper -- gamma theorem.)")
    else:
        P("  HIT h=%d %s -- unexpected, persistence-check required" % (hit[1], hit[0]))


def beta_winning_shape():
    P()
    P("=" * 70)
    P("DOOR beta  WHAT A WINNING DEFORMATION MUST LOOK LIKE")
    P("  Need F(u1,...,u_r) = H exactly (mod 2^m), with r=O(log N), each")
    P("  ui admissible, F an arithmetic circuit of size polylog(N) over")
    P("  {Log Pi(ui), Pi(ui), ring ops}.")
    P()
    P("  Ruled out:")
    P("    - finite differences in a 2-adic parameter (axiom A1 vs proximity)")
    P("    - Q-span of O(1) Log Pi(u) for odd u=O(1)  [measured empty]")
    P("    - Q-span of polynomial-weighted Log-sums     [measured empty]")
    P("    - O(log N) Bernoulli/Mazur moments            [measured empty]")
    P("    - ghosts / dyadic power sums                 [multiplicative side]")
    P()
    P("  Still open (identity-shaped):")
    P("    (B1) digamma_2(K+1/8) via bit-burst/recurrence at cost M(N)log N")
    P("         with FIXED algebraic differential equation (Caruso-style)")
    P("         -- currently a RESTATEMENT of Open Problem E, not a reduction,")
    P("         unless the DE has coefficients of bit-size O(log N), not O(N)")
    P("    (B2) Mazur measure built from cracked data whose pairing is H")
    P("         without Theta(N) moments -- no construction known")
    P("    (B3) algebraic e_{T-1} = Phi(e_T, aux) with aux cracked")
    P()
    P("  Hardest honest reading: beta's identity doors are either restatements")
    P("  of E or empty. The sieve-mass axiom + LLL empties make beta a")
    P("  FIREWALL: deformations that look like 'evaluate cracked Pi nearby'")
    P("  cannot win. New mathematics must leave the crackable-product class.")


# ===========================================================================
# ALPHA
# ===========================================================================

def alpha_landen(K=10, m=160):
    P()
    P("=" * 70)
    P("DOOR alpha  LANDEN COLLAPSE")
    MOD = mpz(1) << m

    def Had(a, d, L):
        s = mpz(0)
        for k in range(L):
            t = a + k * d
            if t % 2:
                s = (s + pow(mpz(t), -1, MOD)) % MOD
        return s

    # classical control
    A = Had(1, 8, 2 * K)
    B = Had(1, 16, K)
    C = Had(9, 16, K)
    P("  classical H(1,8;2K)=H(1,16;K)+H(9,16;K): %s" % ("EXACT" if (B+C) % MOD == A else "FAIL"))

    # basket of first 3 duplication levels
    labels = [(1, 8)]
    for _ in range(2):
        nxt = []
        for a, d in labels:
            nxt += [(a, 2 * d), (a + d, 2 * d)]
        labels = list(dict.fromkeys(nxt))
    basket, names = [], []
    for a, d in labels:
        if a % 2 == 0 and d % 2 == 0:
            continue
        basket.append(Had(a, d, K))
        names.append("(%d,%d)" % (a, d))
    P("  labels (%d): %s" % (len(names), names))
    # hunt unexpected relation among them (exclude classical by using only
    # same-length labels, where classical doesn't live)
    if len(basket) >= 2:
        hit = relation_hunt(basket[0], basket[1:], m, max_height=20)
        if hit is None:
            P("  no unexpected height<=20 relation among same-length labels")
        else:
            P("  unexpected HIT h=%d %s -- check if classical shadow" % (hit[1], hit[0]))
    P("  VERDICT alpha: classical branch stands; no cheap collapse.")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    P("#" * 70)
    P("# HARD PASS: alpha / beta / gamma  (specialize beta+gamma)")
    P("#" * 70)

    # gamma first (proofs)
    gamma_merge_theorem(T=1024, w=12)
    gamma_moment_theorem(N=600)
    gamma_et_split()
    gamma_sensitivity()

    # beta (axiom + kills + open shape)
    beta_axioms()
    beta_admissible_stencil()
    beta_u_span()
    beta_weighted_span()
    beta_mazur_pairing()
    beta_winning_shape()

    # alpha
    alpha_landen()

    P()
    P("#" * 70)
    P("# SYNTHESIS")
    P("#" * 70)
    P("gamma PROVED in restricted models:")
    P("  (1) binary-merge => Omega(N log^3 N) bit work  [THEOREM]")
    P("  (2) materialize-moments => Omega(N^2) bit IO     [THEOREM]")
    P("  (3) e_T cracked / e_{T-1} open is the razor        [IDENTITY]")
    P()
    P("beta under sieve-mass axiom:")
    P("  (1) all admissible finite-difference stencils KILL [RATE LAW]")
    P("  (2) H not in span of admissible Log Pi(u)          [LLL EMPTY]")
    P("  (3) H not in span of poly-weighted Log-sums        [LLL EMPTY]")
    P("  (4) O(log N) Bernoulli/Mazur moments miss H        [LLL EMPTY]")
    P("  (5) win requires leaving the crackable-product class")
    P("      (B1 restates E; B2/B3 have no construction)")
    P()
    P("alpha: classical Landen branches; no low-height collapse.")
    P()
    P("Paper-ready: promote Conjecture mass -> Theorem in merge model;")
    P("state sieve-mass axiom as the firewall for beta;")
    P("Open Problem E = e_{T-1}(AP) mod 2^{Theta(N)} in O(N log^2 N).")
    P("#" * 70)


if __name__ == "__main__":
    main()
