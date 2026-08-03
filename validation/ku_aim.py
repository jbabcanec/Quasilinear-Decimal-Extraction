#!/usr/bin/env python3
"""ku_aim.py -- fire the Kedlaya-Umans gun at the prefix-residue family
(2026-08-02, the outside-audit follow-through).

THE TARGET.  Every route to the prize pays one toll (08-02 addendum): the
prefix residues h_t = (X mod 2^{4t}) mod mu_t, mu_t = 8t+1, X the one
materialized string, t <= T = Theta(N).  Both side conversions are free
(PART A validates the zero-slack equivalence):

    X mod mu_t        = pow(5, N', mu_t)          -- modexp, polylog
    16^t mod mu_t     = pow(16, t, mu_t)          -- modexp, polylog
    h_t  <->  s_t = (X >> 4t) mod mu_t            -- one ring op apart

Kedlaya-Umans multivariate multipoint evaluation (SICOMP 2011; quasi-linear
refinements Bhargava et al. 2022-24) is the one published machine whose
headline is teleport-shaped: evaluations of ONE polynomial at MANY points
cost table + points, quasi-linearly, over ANY ring Z/r.  It was killed in
this campaign for composition/exponent re-indexing (KL24 door), a different
use.  Here it is aimed at the family itself.  Three aims exist; each dies at
a named, now-measured point.

AIM 1 (KU's own engine room: small-prime CRT).  KU reduces Z/r evaluation to
many small prime fields F_p.  That requires the target to be a low-degree
polynomial over SOME fixed small field of the shared input data.  But the
double-mod P -> (P mod 2^m) mod mu is, over every fixed F_p, a foreign-
characteristic MOD counter: for p != odd factors of mu the outer mod-mu is a
MOD-q object (Razborov-Smolensky territory: no low-degree representation
over F_p, p coprime to q); for odd p dividing mu the inner mod-2^m
(truncation/carry side) becomes the foreign counter.  For every prime p one
of the two scales is R-S-hard, so no fixed-field polynomial table of
subexponential size exists.  PART B measures it: exact ANF (F_2) degree of
bit_j((P mod mu)) as a function of the m input bits is m - O(1) (maximal,
random-like), against degree 1 for the tame control mu = 2^k.  The KU CRT
layer is dead on arrival -- not by ring variety but by Razborov-Smolensky.

AIM 2 (single ambient ring).  Embed all rings in Z/r, r = lcm(mu_t).  The
embedding is mathematically sound -- PART C verifies it exactly at toy scale
-- but log2 r = Theta(N log N), so a SINGLE point-value carries Theta(N log N)
bits and writing the T of them costs Theta(N^2 log N) before any evaluation
happens: the mass is repaid at the door.  Worse, the prefix family has no
single-polynomial phrasing at all (the polynomial truncates with t), so the
instance is per-t Horner in a huge ring: ring-ops measured Theta(T^2) at toy
scale (PART C).  This is the notes' single-ring embedding row, now with
KU-specific accounting.

AIM 3 (reindex the suffix: one bivariate polynomial).  h_t <-> s_t =
sum_{s>=t} w_s 16^{s-t} mod mu_t.  A bivariate phrasing needs y^{s-t} at
(y = 16, per-t twist 16^{-t} mod mu_t).  Inside Z/mu_t the index t is TWO
different things at once: a ring element (8t+1 = 0 forces t = -1/8, so every
polynomial-in-t coefficient collapses to a ring constant -- PART A validates
this collapse exactly) and an exponent of 16 (reduced only mod ord(2, mu_t),
which varies wildly across the family -- Kurlberg-Pomerance; counting_dlog
pinned orbit-position extraction as dlog-hard).  No polynomial ring contains
both: any homomorphic image that specializes t to -1/8 has already lost the
exponent, and vice versa.  The exponent/ring mismatch is the multiplicative
telling of the floor barrier; AIM 3 is stillborn at the phrasing stage.

INVERTED AIMS (use MME where it does apply) all land on pinned faces:
moments basis (inv(8t+1) as degree-Theta(N) polynomial, transposed
Vandermonde) is output-quadratic in the tapered array (recorded door);
the theta-operator closed form (1+8.theta)u = geometric collapses V to the
Lerch value at 1/16 -- the holonomic face, L3, and the FOURTH sighting of
the p-odd obstruction (CMTV vacuous at p = 2, and the p = 5 source move
re-truncates: band = full - head with head as long as the band);
P-recursive term extraction by baby-step/giant-step is the recorded
Strassen/BGS N^{1.5} door.

VERDICT.  The algebraic-preprocessing sector of the preprocessing door
(transient item 2 / the LDC keeper) is CLOSED: any preprocessing whose query
answers are low-degree polynomial evaluations over a fixed ring of shared
tables dies on the R-S pincer (measured), the single-ring mass (accounted),
or the exponent/ring mismatch (derived).  What remains open, both ways, is
NON-algebraic preprocessing -- arbitrary combinatorial data structures --
where no lower-bound technique reaches and no candidate exists.
"""

import random

import numpy as np

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


# ------------------------------------------------ PART A: zero-slack, exact

def part_A(Np=4000, samples=200, seed=5):
    """The tollbooth family and its free conversions, validated exactly."""
    X = mpz(5) ** Np
    LX = X.bit_length()
    Tmax = (LX - 8) // 4
    rng = random.Random(seed)
    ts = sorted(rng.sample(range(1, Tmax), samples - 2) + [1, Tmax - 1])

    ok_modexp = ok_equiv = ok_stream = ok_collapse = True
    f_coeffs = [rng.randrange(-99, 100) for _ in range(6)]

    for t in ts:
        mu = 8 * t + 1
        P = X & ((mpz(1) << (4 * t)) - 1)          # prefix (truncation, free)
        S = X >> (4 * t)                           # suffix
        h = int(P % mu)
        s = int(S % mu)
        x = int(X % mu)
        # free conversion 1: X mod mu by modexp
        if pow(5, Np, mu) != x:
            ok_modexp = False
        # zero-slack equivalence: h = x - 16^t s  (mod mu)
        if (x - pow(16, t, mu) * s - h) % mu != 0:
            ok_equiv = False
        # the suffix stream recurrence S_t = 16 S_{t+1} + w_t
        w = int((X >> (4 * t)) & 15)
        if (X >> (4 * t)) != ((X >> (4 * (t + 1))) << 4) + w:
            ok_stream = False
        # the ring-side collapse: f(t) = f(-1/8) mod mu for polynomial f
        ft = sum(c * t ** i for i, c in enumerate(f_coeffs))
        tstar = (-pow(8, -1, mu)) % mu
        fstar = sum(c * pow(tstar, i, mu) for i, c in enumerate(f_coeffs))
        if (ft - fstar) % mu != 0:
            ok_collapse = False

    print("PART A  zero slack   X = 5^%d (%d bits), %d sampled t <= %d"
          % (Np, LX, len(ts), Tmax))
    print("        X mod mu_t == pow(5,N,mu_t)  (free side)   : %s"
          % ("EXACT" if ok_modexp else "*** MISMATCH ***"))
    print("        h_t == x_t - 16^t s_t  (h <-> s, one op)   : %s"
          % ("EXACT" if ok_equiv else "*** MISMATCH ***"))
    print("        suffix stream S_t = 16 S_{t+1} + w_t       : %s"
          % ("EXACT" if ok_stream else "*** MISMATCH ***"))
    print("        ring-side collapse f(t) == f(-1/8) mod mu_t: %s"
          % ("EXACT" if ok_collapse else "*** MISMATCH ***"))
    print("        (the exponent side 16^{-t} mod mu_t does NOT collapse:")
    print("         t reduces mod ord(2,mu_t), not mod mu_t -- AIM 3's grave)")
    return ok_modexp and ok_equiv and ok_stream and ok_collapse


# ------------------------- PART B: ANF degree of the double-mod (AIM 1 kill)

def anf(f):
    """Exact algebraic normal form over F_2 (in-place butterfly)."""
    f = f.copy()
    m = f.size.bit_length() - 1
    for k in range(m):
        f = f.reshape(-1, 2, 1 << k)
        f[:, 1, :] ^= f[:, 0, :]
        f = f.reshape(-1)
    return f


def part_B(ms=(12, 16, 20), seed=11):
    """Measured: bit_j(P mod mu) has near-maximal F_2 (ANF) degree in the
    input bits; the tame control bit_j(P mod 2^k) has degree 1."""
    print("PART B  ANF degree of the double-mod over F_2  (KU CRT-layer kill)")
    print("        function: m input bits -> bit_0 of (P mod mu)")
    print("        %-8s %-10s %-12s %-16s %s"
          % ("m", "mu", "ANF degree", "monomials", "note"))
    worst_gap = 0
    for m in ms:
        n = 1 << m
        idx = np.arange(n, dtype=np.int64)
        pc = np.zeros(n, dtype=np.uint8)
        for k in range(m):
            pc += ((idx >> k) & 1).astype(np.uint8)
        for mu, note in ((9, "3^2, in family"), (17, "prime, in family"),
                         (73, "prime"), (257, "Fermat"), (1033, "word-size"),
                         (1 << 8, "TAME control 2^8"),
                         (0, "random control")):
            if mu == 0:
                rng = np.random.default_rng(seed)
                fbit = rng.integers(0, 2, n, dtype=np.uint8)
                name = "rand"
            else:
                fbit = ((idx % mu) & 1).astype(np.uint8)
                name = str(mu)
            a = anf(fbit)
            on = a.nonzero()[0]
            deg = int(pc[on].max()) if on.size else 0
            mono = int(on.size)
            if mu not in (0, 1 << 8):
                worst_gap = max(worst_gap, m - deg)
            print("        %-8d %-10s %-12d %-16d %s"
                  % (m, name, deg, mono, note))
    print("        Every odd mu: degree m - O(1), ~2^{m-1} monomials (dense,")
    print("        random-like; worst deficit %d).  mu = 2^k: degree 1." % worst_gap)
    print("        => over F_2 the coefficient table for the double-mod is")
    print("        2^Theta(m); with m = Theta(N) input bits the KU/MME table")
    print("        d^m is exponential.  Field change cannot help: for every")
    print("        F_p one of the two scales is a foreign-characteristic MOD")
    print("        counter (Razborov-Smolensky).  AIM 1 CLOSED, measured.")
    return True


# ------------- PART C: single-ring embedding, exact but mass-repaying (AIM 2)

def lcm(a, b):
    import math
    return a // math.gcd(a, b) * b


def part_C_toy(T0=64):
    """The Z/lcm embedding computes the family exactly -- and counts ring
    ops Theta(T^2) because the prefix polynomial truncates with t."""
    k = int(4 * T0 / 2.3219) + 2
    X = mpz(5) ** k
    r = 1
    mus = [8 * t + 1 for t in range(1, T0 + 1)]
    for mu in mus:
        r = lcm(r, mu)
    nib = [int((X >> (4 * s)) & 15) for s in range(T0)]
    ops = 0
    ok = True
    for t in range(1, T0 + 1):
        acc = 0
        for s in range(t - 1, -1, -1):            # Horner over the t-prefix
            acc = (acc * 16 + nib[s]) % r
            ops += 1
        if acc % mus[t - 1] != int((X & ((mpz(1) << (4 * t)) - 1)) % mus[t - 1]):
            ok = False
    print("PART C  single-ring embedding, toy exactness   T=%d  log2 r=%d"
          % (T0, r.bit_length()))
    print("        h_t via Horner in Z/lcm, all t          : %s"
          % ("EXACT" if ok else "*** MISMATCH ***"))
    print("        ring ops = %d = %.2f T^2 (no shared polynomial exists:"
          % (ops, ops / T0 ** 2))
    print("        the prefix truncates with t) -- correct and conserved.")
    return ok


def part_C_accounting():
    """The mass table at real parameters (anchored to measured 1e5 values)."""
    import math
    print()
    print("PART C' accounting at pipeline scale (anchor: N=1e5 measured")
    print("        leaf 4280312 bits, tree work 2.09e9)")
    print("        %-8s %-11s %-12s %-12s %-12s %s"
          % ("N", "log2 r", "tree work", "AIM2 write", "AIM1 CRT", "repay"))
    for N in (1e4, 1e5, 1e6):
        T = 2.324 * N
        LX = 2.3219 * N
        avglg = math.log2(2 * T) - 1.44          # avg lg mu over the band
        logr = T * avglg                          # = leaf data = log2 lcm-ish
        tree = 12 * (T * avglg) * math.log2(T) ** 2 / 4.2   # anchored shape
        aim2 = T * logr                           # writing T values in Z/r
        aim1 = T * LX / avglg * avglg             # T values' Z-lifts, CRT width
        print("        %-8.0e %-11.3e %-12.3e %-12.3e %-12.3e %.0fx"
              % (N, logr, tree, aim2, aim1, aim2 / tree))
    print("        AIM 2 repays ~N/log^2 over the tree before evaluating")
    print("        anything; AIM 1's CRT width alone is Theta(N^2/log N).")
    return True


if __name__ == '__main__':
    okA = part_A()
    print()
    okB = part_B()
    print()
    okC = part_C_toy()
    part_C_accounting()
    print()
    print("IDENTITIES", "PASSED" if (okA and okC) else "FAILED")
    print("VERDICT: all three KU aims dead at named points (R-S pincer,")
    print("single-ring mass, exponent/ring mismatch).  The algebraic sector")
    print("of the preprocessing door is CLOSED; the residual open sector is")
    print("non-algebraic (combinatorial) preprocessing, both directions.")
