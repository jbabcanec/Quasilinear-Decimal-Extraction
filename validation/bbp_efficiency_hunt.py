#!/usr/bin/env python3
"""bbp_efficiency_hunt.py -- how large can g_eff get over ALL BBP-type formulas?
(2026-07-30)

formula_mass.py pinned the third log to ONE formula property:

    g_eff = (answer bits) / (band terms) = w / r
            for a base-2^w degree-1 formula with r denominators per step,

    third log = log2(mu_max) / g_eff .

Parity needs g_eff ~ log2 N.  The hex BBP formula has g_eff = 1.00, Bellard
1.43.  So: is g_eff bounded over the whole BBP-type family, or can a formula
with a big base and few denominators exist?

Two experiments.

(A) THE DYADIC-MACHIN DREAM, KILLED EXACTLY.
    arctan(1/2^a) alone would give g_eff = 2a with ONE word-size denominator
    per term -- unbounded efficiency, instant parity.  We prove no rational
    combination of them is pi/4:
      * odd denominator D: e^{i D pi/4} drags in zeta_8, and the
        (1+i)-adic valuation of the identity is 0 on one side, odd on the
        other -- contradiction (proved here, restated in code comments);
      * even D: the identity forces a Z-relation among the Gaussian
        factorizations of z_a = 2^a + i.  We COMPUTE the exponent matrix
        (which of the two Gaussian primes over each p divides z_a, with
        multiplicity) and check its rank.  Full rank = no relation exists.

(B) THE HUNT.  PSLQ over degree-1 BBP-type formulas
        pi = sum_j a_j * S_j(m, b),   S_j = sum_k b^{-k}/(mk+j),  b = +-2^w,
    for a sweep of (m, w).  Every hit is scored by g_eff = w / #nonzero a_j.
    Known formulas must reappear (hex BBP: m=8, w=4, r=4, g_eff=1) or the
    harness is broken.
"""

import math
from fractions import Fraction

from mpmath import mp, mpf, pslq, pi as mp_pi
from sympy import factorint


# ===========================================================================
# (A) dyadic Machin impossibility
# ===========================================================================

def gaussian_prime_for(p):
    """Canonical Gaussian prime x+yi over p = 1 mod 4, with x odd > 0, y > 0."""
    for x in range(1, int(math.isqrt(p)) + 1):
        y2 = p - x * x
        y = math.isqrt(y2)
        if y * y == y2 and y > 0:
            if x % 2 == 0:
                x, y = y, x
            return x, y
    raise ValueError("no representation for %d" % p)


def dyadic_arctan_rank(A=28, verbose=True):
    """Exponent-difference matrix f[a][p] = v_pi(z_a) - v_pibar(z_a),
    z_a = 2^a + i, norm 4^a + 1.  Returns (rank, nrows)."""
    rows = []
    primes = {}
    for a in range(1, A + 1):
        Na = 4 ** a + 1
        fac = factorint(Na)
        row = {}
        for p, e in fac.items():
            assert p % 4 == 1, (a, p)            # no (1+i) ever divides z_a
            x, y = gaussian_prime_for(p)
            primes.setdefault(p, (x, y))
            # i == -x * y^{-1} (mod p) in Z[i]/(x+yi); pi_p | z_a iff
            # 2^a == x*y^{-1} (mod p), else the conjugate divides it.
            t = (x * pow(y, -1, p)) % p
            row[p] = e if pow(2, a, p) == t else -e
        rows.append(row)
    cols = sorted(primes)
    idx = {p: j for j, p in enumerate(cols)}
    M = [[Fraction(r.get(p, 0)) for p in cols] for r in rows]
    # exact Gaussian elimination over Q
    rank = 0
    ncols = len(cols)
    for c in range(ncols):
        piv = None
        for r in range(rank, len(M)):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        pv = M[rank][c]
        M[rank] = [v / pv for v in M[rank]]
        for r in range(len(M)):
            if r != rank and M[r][c] != 0:
                f = M[r][c]
                M[r] = [vr - f * vk for vr, vk in zip(M[r], M[rank])]
        rank += 1
        if rank == len(M):
            break
    if verbose:
        print("(A) dyadic Machin: z_a = 2^a + i for a = 1..%d" % A)
        print("    distinct rational primes involved: %d" % len(cols))
        print("    rank of the exponent-difference matrix: %d / %d rows"
              % (rank, len(rows)))
        if rank == len(rows):
            print("    => FULL RANK: no nontrivial Z-relation among the z_a")
            print("       up to conjugation, so no combination of arctan(1/2^a)")
            print("       is any rational multiple of pi (even-D case).")
        else:
            print("    => RANK DEFICIENT: a relation may exist -- investigate!")
        print("    odd-D case: pi/4 * D with D odd forces zeta_8, and")
        print("       v_(1+i)(LHS) = 0 while v_(1+i)(RHS) is odd. Contradiction.")
        print("    CONSEQUENCE: the 'one word-size denominator per 2a bits'")
        print("       formula (g_eff = 2a, unbounded) does not exist for pi.")
    return rank, len(rows)


# ===========================================================================
# (B) BBP-type efficiency hunt
# ===========================================================================

def S_values(m, w, sign):
    """S_j = sum_k b^{-k}/(mk+j), j=1..m, b = sign*2^w."""
    terms = int(mp.dps * 3.33 / w) + 6
    out = []
    for j in range(1, m + 1):
        s = mpf(0)
        num = mpf(1)
        for k in range(terms):
            s += num / (m * k + j)
            num /= mpf(sign) * (mpf(2) ** w)
        out.append(s)
    return out


def hunt(mlist, wlist, maxcoeff=10 ** 8, dps=140, verbose=True):
    mp.dps = dps
    hits = []
    for m in mlist:
        for w in wlist:
            for sign in (1, -1):
                vec = [mp_pi] + S_values(m, w, sign)
                rel = pslq(vec, maxcoeff=maxcoeff, maxsteps=20000,
                           tol=mpf(10) ** (-(dps - 25)))
                if not rel or rel[0] == 0:
                    continue
                a = rel[1:]
                r = sum(1 for v in a if v != 0)
                if r == 0:
                    continue
                # verify at higher precision
                mp.dps = dps + 60
                v2 = [mp_pi] + S_values(m, w, sign)
                resid = sum(c * x for c, x in zip(rel, v2))
                mp.dps = dps
                if abs(resid) > mpf(10) ** (-(dps - 20)):
                    continue
                hits.append((w / r, m, w, sign, r, rel))
    hits.sort(reverse=True)
    if verbose:
        print()
        print("(B) BBP-type degree-1 formulas for pi found by PSLQ")
        print("    scored by g_eff = w / (nonzero denominators per step)")
        print("    %-7s %-4s %-4s %-6s %-4s  %s"
              % ("g_eff", "m", "w", "base", "r", "coefficients (pi, S_1..S_m)"))
        for g, m, w, sign, r, rel in hits[:25]:
            base = ("+2^%d" % w) if sign > 0 else ("-2^%d" % w)
            print("    %-7.3f %-4d %-4d %-6s %-4d  %s"
                  % (g, m, w, base, r, rel))
        if not hits:
            print("    (no relations found in this window)")
    return hits


def main():
    print("=" * 74)
    print("Can any BBP-type formula reach g_eff ~ log2 N (parity)?")
    print("=" * 74)
    dyadic_arctan_rank(A=28)

    # sweep: small modulus x large base is where high g_eff would live
    hits = hunt(mlist=list(range(2, 13)), wlist=list(range(1, 41)))
    hits += hunt(mlist=[13, 14, 16, 18, 20, 24], wlist=list(range(1, 25)),
                 verbose=False)
    hits.sort(reverse=True)
    print()
    print("=" * 74)
    if hits:
        best = hits[0]
        print("BEST g_eff FOUND: %.3f  (m=%d, base %s2^%d, r=%d denominators)"
              % (best[0], best[1], "-" if best[3] < 0 else "+", best[2],
                 best[4]))
        print("hex BBP is g_eff = 1.00; Bellard's mixed-modulus formula 1.43.")
        print("Parity needs g_eff ~ log2 N (21 at N=1e6, 40 at N=1e12).")
    print("=" * 74)


if __name__ == "__main__":
    main()
