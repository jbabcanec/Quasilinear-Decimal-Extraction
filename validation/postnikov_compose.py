#!/usr/bin/env python3
"""
LEAD L1 (cracks_dossier.md): compose the Borwein crack with the Hiary-Postnikov
block decomposition of the additive core.

The crack (Thm thm:crack) makes every MULTIPLICATIVE object cheap: for any
arithmetic progression with word-size terms,

    Pi_{u<B} (mu_u + t)  mod 2^m      in O(M(N) log N),

hence the one cheap ADDITIVE functional of the terms, G(t) = log_2 Pi(mu_u + t).
Postnikov blocks turn the atom into inverse power sums over a block.  Composing
them means: measure the block interiors A_r with the cheap G(t).

    PART A   block tensor identity  (exact, mod 2^m)
             S = sum_i (-1)^i 2^{(j+3)i} P_i A_{i+1},
             P_i = sum_{v<W} v^i,   A_r = sum_{u<B} mu_u^{-r}
    PART B   measurement identity   (exact, mod 2^m)
             G(t) - G(0) = sum_{r>=1} (-1)^{r-1} t^r A_r / r
    PART C   THE RATE LAW           how many bits of A_1 per crack-evaluation,
             measured against ground truth for stencils at 2-adic depth a
    PART D   module hunt            is the atom in the small-height Z-span of
             crack-cheap log-sums at all?  (LLL, with +/- controls)

If PART C's rate is superconstant in the cost, or PART D finds a relation, the
composition is a crack and the digit goes to O(N log^2 N).  Otherwise the route
is closed WITH A MEASURED LAW.
"""

import random
from fractions import Fraction

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


# ---------------------------------------------------------------- 2-adic tools

def v2(x):
    """2-adic valuation of a nonzero integer."""
    x = int(x)
    if x == 0:
        return None
    return (x & -x).bit_length() - 1


def log2adic(x, m):
    """Iwasawa 2-adic logarithm of an odd integer x, mod 2^m.

    log(x) := (1/2) log(x^2), x^2 in 1 + 8Z_2, so the series in y = x^2 - 1
    converges with v(y^n/n) >= 3n - log2(n).
    """
    MOD = mpz(1) << (m + 32)          # guard bits for the 1/2 and the 1/n
    y = (mpz(x) * mpz(x) - 1) % MOD
    term, total, n = y, mpz(0), 1
    while True:
        e = v2(n)
        # v(y^n / n) = 3n - e  (v(y) >= 3);  stop when it exceeds the guard
        if 3 * n - e > m + 24:
            break
        odd = n >> e
        total = (total + (term // (mpz(1) << e)) * pow(odd, -1, MOD)
                 * (1 if n % 2 else -1)) % MOD
        n += 1
        term = (term * y) % MOD
    half = total // 2                  # exact: v(total) >= 3
    return int(half % (mpz(1) << m))


def sum_log(terms, m):
    """sum of 2-adic logs = log of the product (the crack-cheap functional)."""
    MOD = mpz(1) << m
    prod = mpz(1)
    for x in terms:
        prod = (prod * x) % MOD
    return log2adic(prod, m)


def inv_power_sum(B, r, m, d=8, b=1):
    """A_r = sum_{u<B} (d u + b)^{-r} mod 2^m."""
    MOD = mpz(1) << m
    return int(sum(pow(mpz(d * u + b), -r, MOD) for u in range(B)) % MOD)


def atom(K, m, d=8, b=1):
    """S = sum_{k<K} (d k + b)^{-1} mod 2^m -- the additive core, unweighted."""
    return inv_power_sum(K, 1, m, d, b)


# --------------------------------------------------- PART A: the block tensor

def part_A(j=6, W=64, m=384):
    """S over K = 2^j * W terms  ==  sum_i (-1)^i 2^{(j+3)i} P_i A_{i+1}."""
    B, MOD = 1 << j, mpz(1) << m
    K = B * W
    lhs = atom(K, m)

    rhs, i = mpz(0), 0
    while (j + 3) * i < m:
        P_i = sum(mpz(v) ** i for v in range(W)) % MOD
        A_r = inv_power_sum(B, i + 1, m)
        rhs = (rhs + (-1) ** i * (mpz(1) << ((j + 3) * i)) * P_i * A_r) % MOD
        i += 1
    ok = lhs == int(rhs % MOD)
    print("PART A  block tensor split   K=2^%d*%d=%d  m=%d  i-terms=%d   %s"
          % (j, W, K, m, i, "EXACT" if ok else "*** MISMATCH ***"))
    return ok


# ------------------------------------------- PART B: the measurement identity

def part_B(j=6, m=384, a=4):
    """G(t) - G(0) = sum_r (-1)^{r-1} t^r A_r / r   for t = 2^a  (exact)."""
    B, MOD = 1 << j, mpz(1) << m
    mus = [8 * u + 1 for u in range(B)]
    t = 1 << a

    lhs = (sum_log([mu + t for mu in mus], m) - sum_log(mus, m)) % (1 << m)

    rhs, r = mpz(0), 1
    while True:
        e = v2(r)
        if a * r - e >= m:
            break
        A_r = inv_power_sum(B, r, m)
        term = (mpz(t) ** r // (mpz(1) << e)) * pow(r >> e, -1, MOD) * A_r
        rhs = (rhs + (-1) ** (r - 1) * term) % MOD
        r += 1
    ok = lhs == int(rhs % MOD)
    print("PART B  shift measurement    B=2^%d  t=2^%d  m=%d  r-terms=%d   %s"
          % (j, a, m, r - 1, "EXACT" if ok else "*** MISMATCH ***"))
    return ok


# ------------------------------------------------------ PART C: the rate law

def stencil_weights(nodes, order=1):
    """Rational w with sum_l w_l nodes_l^r = [r == order] for r = 1..len(nodes).

    Vandermonde solve: this is the finite-difference stencil that isolates A_order
    and annihilates every other A_r up to the stencil length.
    """
    L = len(nodes)
    # rows r = 1..L, columns l:  M[r-1][l] = nodes[l]^r
    M = [[Fraction(n) ** r for n in nodes] for r in range(1, L + 1)]
    rhs = [Fraction(1) if r == order else Fraction(0) for r in range(1, L + 1)]
    # gaussian elimination (exact)
    for c in range(L):
        p = next(i for i in range(c, L) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        rhs[c], rhs[p] = rhs[p], rhs[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        rhs[c] = rhs[c] / pv
        for i in range(L):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
                rhs[i] = rhs[i] - f * rhs[c]
    return rhs


def recover_A1(nodes, mus, m):
    """Combine the crack-cheap measurements G(t)-G(0) over the stencil.

    Returns (estimate, bits_lost_to_stencil_denominators).
    """
    w = stencil_weights([-t for t in nodes])   # sign: coefficient of A_r is
    # (-1)^{r-1} t^r / r, so use nodes -> -t and fold the 1/r in below.
    # We instead solve directly on the true coefficient matrix:
    L = len(nodes)
    M = [[Fraction((-1) ** (r - 1)) * Fraction(t) ** r / r for t in nodes]
         for r in range(1, L + 1)]
    rhs = [Fraction(1) if r == 1 else Fraction(0) for r in range(1, L + 1)]
    for c in range(L):
        p = next(i for i in range(c, L) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        rhs[c], rhs[p] = rhs[p], rhs[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        rhs[c] = rhs[c] / pv
        for i in range(L):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
                rhs[i] = rhs[i] - f * rhs[c]
    w = rhs

    # clear denominators: w_l = num_l / (2^e * odd)
    from math import gcd
    den = 1
    for x in w:
        den = den * x.denominator // gcd(den, x.denominator)
    e = v2(den) if den % 2 == 0 else 0
    odd = den >> e
    MOD = mpz(1) << m
    G0 = sum_log(mus, m)
    acc = mpz(0)
    for wl, t in zip(w, nodes):
        c = int(wl * den)                      # integer
        R = (sum_log([mu + t for mu in mus], m) - G0) % (1 << m)
        acc = (acc + c * R) % MOD
    acc = acc * pow(odd, -1, MOD) % MOD        # divide by the odd part
    est = int(acc >> e)                        # divide by 2^e: lose e bits
    return est, e


def part_C(j=6, m=320, depths=(1, 2, 4, 6, 8), Ls=(1, 2, 3, 4, 5, 6)):
    """MEASURED: bits of A_1 recovered per crack-evaluation, vs 2-adic depth a."""
    B = 1 << j
    mus = [8 * u + 1 for u in range(B)]
    truth = inv_power_sum(B, 1, m)
    print("PART C  rate law   B=2^%d  m=%d   (bits of A_1 recovered)" % (j, m))
    print("        stencil t = l*2^a, l = 1..L")
    print("        %-6s" % "a" + "".join("  L=%-6d" % L for L in Ls)
          + "   slope(bits/eval)")
    table = {}
    for a in depths:
        row, bits = [], []
        for L in Ls:
            nodes = [l << a for l in range(1, L + 1)]
            est, lost = recover_A1(nodes, mus, m)
            d = (est - truth) % (1 << (m - lost - 4))
            got = (m - lost - 4) if d == 0 else v2(d)
            got = max(0, min(got, m))
            bits.append(got)
            row.append("  %-8d" % got)
        slope = (bits[-1] - bits[0]) / (Ls[-1] - Ls[0]) if len(Ls) > 1 else 0
        table[a] = (bits, slope)
        print("        a=%-4d" % a + "".join(row) + "   %.2f" % slope)
    return table


def part_C2(j=6, m=320, a=6, Ls=(2, 4, 6, 8)):
    """STENCIL SHOOTOUT: is the a+1 bits/evaluation rate an artefact of the
    stencil?  Any family with superlinear bits(L) would be the crack.
    """
    B = 1 << j
    mus = [8 * u + 1 for u in range(B)]
    truth = inv_power_sum(B, 1, m)

    def yield_bits(nodes):
        est, lost = recover_A1(nodes, mus, m)
        d = (est - truth) % (1 << (m - lost - 4))
        return max(0, min((m - lost - 4) if d == 0 else v2(d), m))

    fams = {
        "arithmetic  l*2^a": lambda L: [l << a for l in range(1, L + 1)],
        "geometric   2^(a+l)": lambda L: [1 << (a + l) for l in range(L)],
        "central     +-l*2^a": lambda L: [((l + 2) // 2) * (1 if l % 2 else -1) << a
                                          for l in range(1, L + 1)],
        "mixed depth l*2^(a+l)": lambda L: [l << (a + l) for l in range(1, L + 1)],
        "random even v2=a": lambda L: [(2 * random.Random(7 + l).randrange(1, 999) + 1)
                                       << a for l in range(1, L + 1)],
    }
    # A crack-evaluation at shift t must sieve to 8B + |t|; charge that mass,
    # in units of the unshifted evaluation (8B).  Bits per EVALUATION is the
    # wrong metric -- deep stencils buy bits by paying an exponential sieve.
    def cost(nodes):
        return sum((8 * B + abs(t)) / float(8 * B) for t in nodes)

    print("PART C2 stencil shootout   B=2^%d  m=%d  a=%d" % (j, m, a))
    print("        cost unit = one unshifted crack-evaluation (sieve to 8B)")
    print("        %-22s" % "family" + "".join("  L=%-5d" % L for L in Ls)
          + "  bits/eval  bits/COST")
    for name, gen in fams.items():
        bits = [yield_bits(gen(L)) for L in Ls]
        cs = [cost(gen(L)) for L in Ls]
        per_eval = (bits[-1] - bits[0]) / float(Ls[-1] - Ls[0])
        per_cost = (bits[-1] - bits[0]) / (cs[-1] - cs[0])
        print("        %-22s" % name + "".join("  %-7d" % b for b in bits)
              + "  %-9.2f  %.3f" % (per_eval, per_cost))
    print("        Deep families (geometric, mixed) beat a+1=%d bits per" % (a + 1))
    print("        EVALUATION, and lose badly per unit COST: their yield is")
    print("        ~c*log2(cost), the shallow family's is ~c*cost.  Both need")
    print("        Theta~(N^2/log N) to reach m = Theta(N) bits.")


def part_C3(js=(5, 6, 8), m=384, L=6, depths=range(1, 17)):
    """THE DEPTH OPTIMUM.  The composed route's whole complexity is set by the
    depth a that maximises bits per unit cost.  Claim: the optimum sits at
    a = log2(8B) -- shallow shifts ride the sieve you already pay for, deeper
    ones grow the sieve exponentially for a linear gain.  If the optimum moved
    faster than log2(8B), the composition would beat N^2.
    """
    print("PART C3 depth optimum   m=%d  L=%d  (central stencil +-l*2^a)" % (m, L))
    print("        %-9s %-11s %-28s %s" % ("block B", "log2(8B)",
                                           "argmax_a bits/cost", "peak bits/cost"))
    for j in js:
        B = 1 << j
        mus = [8 * u + 1 for u in range(B)]
        truth = inv_power_sum(B, 1, m)
        best = (None, -1.0)
        curve = []
        for a in depths:
            nodes = [((l + 2) // 2) * (1 if l % 2 else -1) << a
                     for l in range(1, L + 1)]
            est, lost = recover_A1(nodes, mus, m)
            d = (est - truth) % (1 << (m - lost - 4))
            bits = max(0, min((m - lost - 4) if d == 0 else v2(d), m))
            cost = sum((8 * B + abs(t)) / float(8 * B) for t in nodes)
            r = bits / cost
            curve.append((a, bits, r))
            if r > best[1]:
                best = (a, r)
        print("        %-9d %-11.1f %-28s %.3f"
              % (B, (8 * B).bit_length() - 1, "a = %d" % best[0], best[1]))
        print("           a:bits/cost  " + "  ".join("%d:%.2f" % (a, r)
                                                     for a, _, r in curve[:12]))
    print("        Peak yield grows ~ +1 bit per doubling of the sieve bound")
    print("        (4.95, 5.84, 7.78 at log2(8B) = 8, 9, 11), i.e. ~log2(8B)-3")
    print("        bits per unit cost, and the argmax moves with B.  Depth is")
    print("        free only while the shift hides inside the sieve already")
    print("        being paid for: yield ~ log N bits per Theta~(N) evaluation.")


def part_C_cost(K_bits=(10, 12, 14, 16, 18, 20)):
    """The cost side of the measured rate.

    PART C measures: one crack-evaluation at 2-adic depth a = v2(t) buys a+1
    bits of A_1 (slope 5.00/7.00/9.00 at a = 4/6/8).
    The crack needs the shifted terms mu_u + t to be sieve-factorable, so its
    cost is Theta~(max term) = Theta~(8K + 2^a) bit-mass.  Hence

        bits bought / cost  =  (a+1) / (8K + 2^a),

    maximised at a = log2(8K) -- beyond that the sieve grows exponentially
    while the yield grows linearly.  Deeper shifts are not crackable at all.
    """
    print("PART C' cost side   (measured yield a+1 bits ; sieve mass 8K + 2^a)")
    print("        N = m = Theta(K).  evaluations = m/(a+1), each O~(8K + 2^a)")
    print("        %-8s %-11s %-13s %-13s %s"
          % ("log2 K", "a=log2(8K)", "evals", "total/N^2", "verdict"))
    for kb in K_bits:
        K, a = 1 << kb, kb + 3
        evals = K / (a + 1.0)              # m = Theta(K)
        total = evals * K                  # each evaluation ~ N mass (x polylog)
        print("        %-8d %-11d %-13.1f %-13s %s"
              % (kb, a, evals, "1/%.1f" % (a + 1.0), "N^2/log N"))
    print("        => a+1 bits per Theta~(2^a) work: 2-adic depth bought from")
    print("           the crack costs EXPONENTIALLY in the depth.")
    print("           Optimum a = log2(8K):  composed cost = Theta~(N^2 / log N)")
    print("           * WORSE than the existing N log^3 N pipeline. *")


# ------------------------------------------------- PART D: the module hunt

def lll(basis):
    """Textbook LLL over Q (delta = 0.99). basis: list of integer vectors."""
    B = [list(map(Fraction, b)) for b in basis]
    n = len(B)

    def gso():
        Bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = list(B[i])
            for jj in range(i):
                dn = sum(x * x for x in Bs[jj])
                mu[i][jj] = (sum(a * b for a, b in zip(B[i], Bs[jj])) / dn
                             if dn else Fraction(0))
                v = [x - mu[i][jj] * y for x, y in zip(v, Bs[jj])]
            Bs.append(v)
        return Bs, mu

    Bs, mu = gso()
    k = 1
    while k < n:
        for jj in range(k - 1, -1, -1):
            q = mu[k][jj]
            r = int(q + Fraction(1, 2)) if q >= 0 else -int(-q + Fraction(1, 2))
            if r:
                B[k] = [x - r * y for x, y in zip(B[k], B[jj])]
                Bs, mu = gso()
        nk = sum(x * x for x in Bs[k])
        nk1 = sum(x * x for x in Bs[k - 1])
        if nk >= (Fraction(99, 100) - mu[k][k - 1] ** 2) * nk1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            Bs, mu = gso()
            k = max(k - 1, 1)
    return [[int(x) for x in b] for b in B]


def relation_hunt(target, basket, m):
    """Smallest integer relation c0*target + sum ci*basket_i = 0 mod 2^m."""
    n = len(basket)
    W = mpz(1) << m                                   # forces the congruence
    vals = [target] + list(basket)
    rows = []
    for i in range(n + 1):
        row = [0] * (n + 1)
        row[i] = 1
        row.append(int(W * vals[i]))
        rows.append(row)
    rows.append([0] * (n + 1) + [int(W * (mpz(1) << m))])
    best = None
    for r in lll(rows):
        c = r[:n + 1]
        if any(c) and r[-1] == 0:
            h = max(abs(x) for x in c)
            if best is None or h < best[1]:
                best = (c, h)
    return best


BASKET_SPEC = ((8, 1), (8, 3), (8, 5), (8, 7), (16, 1), (16, 9),
               (4, 1), (4, 3), (24, 1), (2, 1))


def build_basket(K, m):
    return [sum_log([d * k + b for k in range(K)], m) for d, b in BASKET_SPEC]


def hunt_and_verify(make_target, label, K, m, m_hi):
    """Find the relation at precision m, then TEST IT at precision m_hi.

    A genuine relation holds at every precision; an LLL artefact dies.  This
    is the decisive test -- height alone is not, since the noise floor for an
    (n+1)-term basket is 2^(m/(n+1)).
    """
    n = len(BASKET_SPEC)
    floor = m / (n + 1.0)
    got = relation_hunt(make_target(m), build_basket(K, m), m)
    if got is None:
        print("        %-26s none found" % label)
        return
    c, h = got
    hb = (h.bit_length() - 1) if h else 0
    # persistence: does the SAME integer relation hold at higher precision?
    hi_vals = [make_target(m_hi)] + build_basket(K, m_hi)
    resid = sum(ci * vi for ci, vi in zip(c, hi_vals)) % (1 << m_hi)
    holds = (resid == 0)
    surviving = m_hi if holds else (v2(resid) if resid else m_hi)
    print("        %-26s height 2^%-5.1f (floor 2^%.1f)  survives to 2^%-4d %s"
          % (label, hb, floor, surviving,
             "REAL RELATION" if holds else "-> ARTEFACT"))


def part_D(K=512, m=192, m_hi=384):
    """Is the atom in the small-height Z-span of crack-cheap log-sums?"""
    print("PART D  module hunt   K=%d  find at m=%d, verify at m=%d" % (K, m, m_hi))
    print("        basket: %d crack-cheap log-sums log P(dk+b), (d,b) = %s"
          % (len(BASKET_SPEC), str(BASKET_SPEC[:4])[:-1] + ", ...)"))
    hunt_and_verify(lambda mm: atom(K, mm),
                    "TARGET atom S=sum 1/(8k+1)", K, m, m_hi)
    hunt_and_verify(lambda mm: (3 * build_basket(K, mm)[0]
                                - 5 * build_basket(K, mm)[1]
                                + build_basket(K, mm)[4]) % (1 << mm),
                    "control(+) 3L0-5L1+L4", K, m, m_hi)
    hunt_and_verify(lambda mm: sum_log([8 * k + 1 for k in range(K, 2 * K)], mm),
                    "control(+) log P over [K,2K)", K, m, m_hi)
    hunt_and_verify(lambda mm: int(random.Random(11).getrandbits(mm)),
                    "control(-) random residue", K, m, m_hi)


# ------------------------------- PART E: the crack under index-graded weights

def part_E(K_list=(1 << 8, 1 << 10, 1 << 12), decay=1):
    """Does the crack survive the dyadic weights that turn the atom into V?

    V weights term t by 2^(Jmax - J_t) with J_t affine in t.  The multiplicative
    analogue is prod (8k+1)^(w_k), w_k = 2^(decay*(K-k)).  Legendre counting
    still works (e_p = sum of a geometric series over an AP -- closed form), so
    the FACTORISATION survives; what dies is the COMPRESSION:  unweighted
    exponents are O(log K) bits each (total Theta(N)), weighted ones are
    Theta(N) bits each (total Theta(N^2)), and exponent-bit assembly runs one
    product tree per exponent BIT.
    """
    print("PART E  index-graded weights vs the crack's compression")
    print("        %-8s %-14s %-14s %-14s %s"
          % ("K", "term mass", "unweighted", "weighted", "ratio"))
    for K in K_list:
        terms = [8 * k + 1 for k in range(1, K + 1)]
        term_mass = sum(t.bit_length() for t in terms)
        # unweighted exponent mass: e_p = #\{k : p | 8k+1\} etc.
        exps, wexps = {}, {}
        for k, t in enumerate(terms, start=1):
            x, p = t, 3
            w = decay * (K - k)                    # log2 of the dyadic weight
            while x > 1 and p * p <= x:
                while x % p == 0:
                    exps[p] = exps.get(p, 0) + 1
                    wexps[p] = wexps.get(p, 0) + (1 << w)
                    x //= p
                p += 2
            if x > 1:
                exps[x] = exps.get(x, 0) + 1
                wexps[x] = wexps.get(x, 0) + (1 << w)
        um = sum(max(1, e.bit_length()) for e in exps.values())
        wm = sum(max(1, e.bit_length()) for e in wexps.values())
        print("        %-8d %-14d %-14d %-14d %.1fx"
              % (K, term_mass, um, wm, wm / float(um)))
    print("        unweighted exponent mass ~ #primes * log K  (the crack's")
    print("        Theta(N)-bit compressed representation);  weighted mass is")
    print("        Theta(K) bits PER PRIME = Theta(N^2/log N) -- the compression")
    print("        that the crack runs on is destroyed by index-graded weights.")
    print("        => V fails the composition for a SECOND, independent reason.")


if __name__ == '__main__':
    okA = part_A()
    okB = part_B()
    print()
    part_C()
    print()
    part_C2()
    print()
    part_C3()
    print()
    part_C_cost()
    print()
    part_D()
    print()
    part_E()
    print()
    print("IDENTITIES", "PASSED" if (okA and okB) else "FAILED")
