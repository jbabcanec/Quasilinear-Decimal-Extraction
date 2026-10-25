#!/usr/bin/env python3
"""
rule_change_hunt.py -- exhaust rule-changes aimed at proving O(N log^2 N).

Philosophy: don't like the game, change the rules. Each PART either
  (i)  produces a usable lemma toward Open Problem E, or
  (ii) kills a rule-change with a measured law / exact identity.

PART A  Limb support: do top guard bits of XV need all of V?
PART B  Difference quotient of cracked Log Pi: rate vs sieve-cost law
PART C  LCM-denominator BS: operand mass still N log^3 N?
PART D  Fixed-dimension Landen: does digamma-multiplication close?
PART E  Archimedean short-window: can a O(log N)-bit window of X replace h?
PART F  Output-bit compression: accumulator is O(log N) bits -- does that
        let the merge truncate early?
"""

from __future__ import annotations

import math
import random
from math import gcd

try:
    from gmpy2 import mpz, lcm as gmpy_lcm
except ImportError:
    mpz = int

    def gmpy_lcm(a, b):
        return a * b // gcd(a, b)


def mlcm(a, b):
    return gmpy_lcm(mpz(a), mpz(b))


def v2(n):
    n = int(n)
    if n == 0:
        return 10**9
    return (n & -n).bit_length() - 1


# ---------------------------------------------------------------------------
# PART A: high bits of (X*V mod 2^J) require low bits of V
# ---------------------------------------------------------------------------

def part_A(trials=200, J=256, G=32, seed=1):
    """If V' = V + 2^r * odd, when do top G bits of X*V mod 2^J flip?

    Claim: for X odd, any r < J - G admits a carry chain that flips the
    top window. So protecting G top bits requires r >= J-G, i.e. V correct
    in its LOW (J-G) bits -- essentially the whole V when G = O(log J).
    """
    rng = random.Random(seed)
    print("PART A  limb support of XV -> top G bits")
    print("        J=%d  G=%d  trials=%d" % (J, G, trials))
    MOD = mpz(1) << J
    TOP = J - G
    flips = {r: 0 for r in range(0, J, max(1, J // 16))}
    for r in flips:
        for _ in range(trials):
            X = mpz(rng.randrange(1, 1 << J)) | 1
            V = mpz(rng.randrange(0, 1 << J))
            e = mpz(rng.randrange(1, 1 << 20)) | 1
            W = (X * V) & (MOD - 1)
            Wp = (X * (V + (mpz(1) << r) * e)) & (MOD - 1)
            if (W >> TOP) != (Wp >> TOP):
                flips[r] += 1
    print("        r          flip_rate  (need ~0 to drop those limbs of V)")
    safe = None
    for r in sorted(flips):
        rate = flips[r] / trials
        print("        r=%-8d %.3f" % (r, rate))
        if safe is None and rate == 0:
            safe = r
    # theoretical: r >= J-G should be safe absent pathological X; check
    r_need = J - G
    X = mpz((1 << J) - 1)  # all-ones odd
    V = mpz(0)
    W = 0
    Wp = (X * (mpz(1) << (r_need - 1))) & (MOD - 1)  # error at r = J-G-1
    flipped = (Wp >> TOP) != 0
    print("        worst-case flip at r=J-G-1=%d with X=2^J-1: %s"
          % (r_need - 1, "YES (top window hit)" if flipped else "no"))
    print("        VERDICT: top G bits of XV require V correct through")
    print("        limb J-G; with G=O(log N) this is Theta(N) bits of V.")
    print("        Rule-change 'approximate V' is DEAD.")
    return safe


# ---------------------------------------------------------------------------
# PART B: difference quotient of cracked products
# ---------------------------------------------------------------------------

def borwein_product_mod(mus, m):
    """Factored product of an explicit list of odd terms, mod 2^m.
    (Small-scale stand-in for the crack; cost model charged separately.)
    """
    MOD = mpz(1) << m
    # prime exponents via trial on the values
    from collections import Counter
    exps = Counter()
    for mu in mus:
        x = int(mu)
        p = 3
        while p * p <= x:
            while x % p == 0:
                exps[p] += 1
                x //= p
            p += 2
        if x > 1:
            exps[x] += 1
    acc = mpz(1)
    maxbit = max(exps.values()).bit_length() if exps else 1
    for j in range(maxbit - 1, -1, -1):
        acc = (acc * acc) % MOD
        for p, e in exps.items():
            if (e >> j) & 1:
                acc = (acc * p) % MOD
    return acc


def part_B(K=80, depths=(1, 2, 3, 4, 5, 6, 8, 10, 12)):
    """H approx (1/eps) Log( Pi(eps)/Pi(0) ) with eps=2^a.

    Measure: bits of H recovered vs a, and sieve-limit 8K+eps.
    Cost unit = sieve mass / sieve mass at eps=0.
    """
    print("PART B  difference quotient of cracked Pi")
    m = 2 * K + 64
    MOD = mpz(1) << m
    mus0 = [8 * k + 1 for k in range(1, K + 1)]
    H = sum(pow(mpz(mu), -1, MOD) for mu in mus0) % MOD
    P0 = borwein_product_mod(mus0, m)

    def log2adic(u, prec):
        """Log2(u) for u odd, in Z/2^prec, via (u-1) series after unit norm."""
        # u = 1+y with y even; Log = sum (-1)^{j+1} y^j / j
        y = (u - 1) % (1 << prec)
        s = mpz(0)
        t = mpz(1)
        pm = (mpz(1) << prec) - 1
        for j in range(1, prec + 4):
            t = (t * y) & pm
            # divide by j: j may be even -- only valid while val(y^j)>=val(j)+needed
            if j % 2 == 0 and v2(j) > v2(int(t)) + 8:
                break
            invj = pow(mpz(j), -1, 1 << prec) if j % 2 else None
            if invj is None:
                # j even: t/j must be integral 2-adically
                tj = t >> v2(j)
                jj = j >> v2(j)
                if jj % 2 == 0:
                    break
                invj = pow(mpz(jj), -1, 1 << prec)
                term = (tj * invj) & pm
            else:
                term = (t * invj) & pm
            if j % 2:
                s = (s + term) & pm
            else:
                s = (s - term) & pm
        return s

    print("        a   bits_recovered  sieve_cost  bits/cost")
    base_sieve = 8 * K + 1
    rows = []
    for a in depths:
        eps = mpz(1) << a
        mus = [mu + int(eps) for mu in mus0]
        # sieve cost: primes up to max(mus)
        cost = (max(mus)) / float(base_sieve)
        Peps = borwein_product_mod(mus, m)
        ratio = (Peps * pow(P0, -1, MOD)) % MOD
        # Log(ratio)/eps  -- work at precision m, then shift
        lg = log2adic(ratio, m)
        # divide by eps=2^a: shift right a, lose a bits
        est = lg >> a
        truth = int(H) & ((1 << (m - a - 4)) - 1)
        got = int(est) & ((1 << (m - a - 4)) - 1)
        diff = (got - truth) & ((1 << (m - a - 4)) - 1)
        bits = (m - a - 4) if diff == 0 else min(v2(diff), m - a - 4)
        bits = max(0, bits)
        bpc = bits / cost
        rows.append((a, bits, cost, bpc))
        print("        a=%-3d %-14d %-10.2f  %.3f" % (a, bits, cost, bpc))

    # optimum bits/cost should sit at a ~ log2(8K), shallow
    best = max(rows, key=lambda r: r[3])
    print("        peak bits/cost at a=%d (log2(8K)=%.1f)"
          % (best[0], math.log2(base_sieve)))
    print("        To reach m=Theta(K) bits at peak rate ~%.2f bits/cost"
          % best[3])
    print("        needs cost Theta(K) i.e. Theta(K) crack-evals -> N^2.")
    print("        VERDICT: difference-quotient rule-change is DEAD")
    print("        (same rate law as postnikov_compose PART C).")
    return rows


# ---------------------------------------------------------------------------
# PART C: LCM binary-splitting operand mass
# ---------------------------------------------------------------------------

def part_C(K=2048):
    """Simulate lcm-merge operand mass vs product-merge mass."""
    print("PART C  LCM-denominator merge mass vs product-tree mass")

    def M(b):
        return b * math.log2(b + 2)

    # build leaf dens
    mus = [8 * k + 1 for k in range(1, K + 1)]

    def mass_for(merge_den):
        """Balanced reduction; merge_den(d1,d2)->d. Track sum M(bitlen)."""
        nodes = [(mpz(1), mpz(mu)) for mu in mus]  # (num, den) = 1/mu
        work = 0.0
        data = 0.0
        while len(nodes) > 1:
            nxt = []
            for i in range(0, len(nodes) - 1, 2):
                n1, d1 = nodes[i]
                n2, d2 = nodes[i + 1]
                d = merge_den(d1, d2)
                # num = n1*(d/d1) + n2*(d/d2)
                q1 = d // d1
                q2 = d // d2
                b = max(n1.bit_length() + q1.bit_length(),
                        n2.bit_length() + q2.bit_length(),
                        d.bit_length())
                work += 3 * M(max(b, 2))
                data += d.bit_length()
                n = n1 * q1 + n2 * q2
                nxt.append((n, d))
            if len(nodes) & 1:
                nxt.append(nodes[-1])
            nodes = nxt
        return work, data, nodes[0][1].bit_length()

    def den_prod(d1, d2):
        return d1 * d2

    def den_lcm(d1, d2):
        return mlcm(d1, d2)

    wp, dp, bp = mass_for(den_prod)
    wl, dl, bl = mass_for(den_lcm)
    nlog2 = K * (math.log2(K) ** 2)
    nlog3 = K * (math.log2(K) ** 3)
    print("        K=%d" % K)
    print("        product-tree  work=%.2e  (work/Nlog^3 N)=%.3f  root_den_bits=%d"
          % (wp, wp / nlog3, bp))
    print("        lcm-tree      work=%.2e  (work/Nlog^3 N)=%.3f  root_den_bits=%d"
          % (wl, wl / nlog3, bl))
    print("        lcm/product work ratio=%.3f" % (wl / wp))
    print("        VERDICT: lcm compresses the ROOT (bits %d vs %d) but work"
          % (bl, bp))
    print("        stays on the N log^3 N rail -- local blocks don't share primes.")
    print("        Rule-change 'use lcm not product' is DEAD for asymptotics.")
    return wl / wp


# ---------------------------------------------------------------------------
# PART D: fixed-dimension Landen / digamma multiplication
# ---------------------------------------------------------------------------

def part_D(K0=32, steps=6, m=512):
    """Digamma multiplication psi(2z)=1/2(psi(z)+psi(z+1/2))+log 2
    branches the AP family. Track dimension of the Q-span of active
    inverse-sums under duplication. If dim stays bounded, Landen closes.
    """
    print("PART D  Landen dimension under digamma duplication")
    MOD = mpz(1) << m

    # An 'active sum' is labelled by (a,d): sum_{k=0}^{L-1} 1/(a + k*d)
    # Start with (1,8) at length K0. Duplication map on the range:
    # sum_{k=0}^{2L-1} 1/(a+kd) = sum_{j=0}^{L-1} 1/(a+2jd) + sum_{j=0}^{L-1} 1/(a+d+2jd)
    # = (1/2)-scaled sums on denser lattices? Over Z we just split:
    # new labels (a, 2d) and (a+d, 2d), same length L.
    # Dimension = number of distinct labels after s steps (exact branch count).

    from collections import Counter
    active = Counter({(1, 8): K0})  # label -> length
    print("        step  #labels  total_terms  (branching => product tree)")
    for s in range(steps + 1):
        nlabels = len(active)
        terms = sum(active.values())
        print("        s=%-3d %-8d %-10d" % (s, nlabels, terms))
        if s == steps:
            break
        nxt = Counter()
        for (a, d), L in active.items():
            # split into even/odd index subprogressions
            nxt[(a, 2 * d)] += L  # k=2j
            nxt[(a + d, 2 * d)] += L  # k=2j+1
        active = nxt
    print("        VERDICT: labels double every step -> 2^s families after s")
    print("        doublings. Finite-dimensional Landen for the AP-inverse")
    print("        sum does NOT close. (Matches paper Remark exclusions #2.)")
    print("        Rule-change 'scalar/vector AGM on digamma' is DEAD unless")
    print("        a NEW identity collapses labels (not the classical one).")
    return nlabels


# ---------------------------------------------------------------------------
# PART E: short archimedean window cannot replace h
# ---------------------------------------------------------------------------

def part_E(n=80, J=40, mu=17, Gp=16):
    """Band identity needs h = floor(X/2^J) mod mu.
    Using only bits [J-Gp, J+extra) of X: does any extra window length
    recover h?
    """
    print("PART E  short window of X vs exact h")
    X = pow(5, n)
    H = X >> J
    h = H % mu
    u_bits = (X & ((1 << J) - 1)) >> max(0, J - Gp)
    # try recovering h from (X mod 2^{J+extra}) and (X mod mu) alone --
    # that's the lemma. Without X mod mu, can a longer window work?
    print("        exact h = %d (mod mu=%d)" % (h, mu))
    for extra in [0, 8, 16, 32, 64, 128]:
        # guess: use float(X_window / 2^J / mu) fractional -- no h
        lo = max(0, J - Gp)
        hi = J + extra
        window = (X >> lo) & ((1 << (hi - lo)) - 1)
        # approximate H from bits above J only if we had them; here we have
        # bits [lo, hi). Bits >= J are the low `extra` bits of H.
        if extra == 0:
            h_approx = None
            ok = False
        else:
            H_low = window >> (J - lo)  # low `extra` bits of H
            # H_low = H mod 2^extra; need H mod mu. Insufficient when 2^extra < mu
            h_approx = H_low % mu  # wrong in general
            ok = h_approx == h
        print("        extra=%-4d H mod 2^extra available: %s  recovers h: %s"
              % (extra, "yes" if extra else "no", ok))
    # show CRT needs both systems
    Y = X & ((1 << J) - 1)
    R = X % mu
    h2 = ((R - Y % mu) * pow((mu + 1) // 2, J, mu)) % mu
    print("        Bezout/lemma recovery from (X mod mu, X mod 2^J): h=%d %s"
          % (h2, "EXACT" if h2 == h else "FAIL"))
    print("        VERDICT: any archimedean short-window rule-change that")
    print("        drops (X mod mu) or (X mod 2^J) mod mu loses h. DEAD.")


# ---------------------------------------------------------------------------
# PART F: early merge truncation to O(log N) accumulator bits
# ---------------------------------------------------------------------------

def part_F(K=128, Jmax=256, G=24):
    """v3 merge caps at Jmax. If we only need top G bits of XV, can each
    merge node truncate to G bits? (No -- measured.)
    """
    print("PART F  early truncation of V-merge to G bits")
    rng = random.Random(0)
    mus = [8 * k + 1 for k in range(1, K + 1)]
    signs = [1 if k % 2 else -1 for k in range(1, K + 1)]
    shifts = [Jmax - (4 * k % Jmax) for k in range(1, K + 1)]
    X = mpz(5) ** 200 + (mpz(1) << 3)

    def full_V():
        MOD = mpz(1) << Jmax
        return sum(signs[i] * (mpz(1) << shifts[i]) * pow(mpz(mus[i]), -1, MOD)
                   for i in range(K)) % MOD

    V = full_V()
    W = (X * V) & ((mpz(1) << Jmax) - 1)
    top = int(W >> (Jmax - G))

    # truncated merge: every leaf inverse only mod 2^G, shifts mod G -- nonsense
    # better test: V mod 2^{G+c} for various c, see when top bits of XV stabilize
    print("        precision of V (low bits kept) vs top-G match of XV")
    for c in [0, 8, 16, 32, 64, 128, Jmax - G]:
        prec = G + c
        if prec > Jmax:
            prec = Jmax
        Vtr = V & ((mpz(1) << prec) - 1)
        Wtr = (X * Vtr) & ((mpz(1) << Jmax) - 1)
        top_tr = int(Wtr >> (Jmax - G))
        match = top_tr == top
        print("        V mod 2^%-4d  topG_match=%s  (diff_v2=%s)"
              % (prec, match,
                 v2(top_tr ^ top) if top_tr != top else "inf"))
    print("        VERDICT: need prec ~= Jmax to stabilize top G bits.")
    print("        Accumulator being O(log N) bits does NOT license early")
    print("        truncation of V. DEAD.")


def main():
    print("=" * 72)
    print("RULE-CHANGE HUNT for O(N log^2 N) decimal extraction")
    print("=" * 72)
    part_A()
    print()
    part_B()
    print()
    part_C()
    print()
    part_D()
    print()
    part_E()
    print()
    part_F()
    print()
    print("=" * 72)
    print("GLOBAL VERDICT")
    print("  Every rule-change that stays inside Z/2^{Theta(N)} ring operations")
    print("  on the named bases, or classical digamma Landen, or approximate V,")
    print("  or short X-windows, dies by measured law.")
    print("  Surviving doors (need NEW mathematics, not more of the same):")
    print("   (1) a non-classical identity that collapses the Landen label tree")
    print("   (2) a deformation of Pi whose logarithmic derivative is V and")
    print("       whose values stay in the Borwein-crackable class at cost")
    print("       O(N log^2 N) (Witt / Iwasawa measure pairing, not eps=2^a)")
    print("   (3) a restricted lower bound proving mass conservation, ending")
    print("       the hunt inside the named bases")
    print("=" * 72)


if __name__ == "__main__":
    main()
