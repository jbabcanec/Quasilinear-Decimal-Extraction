#!/usr/bin/env python3
"""merge_kary.py -- ROUND 3 of the transform-reuse attack (2026-08-01).

ROUND 1 pinned work = leaf_data x sum_l log2(b_l), with b_l doubling per level:
the third log is the operand growth over the tree's log T levels.  Since
nothing is reusable across levels (every operand is a fresh product), the only
asymptotic escape is to REMOVE LEVELS -- merge k fractions at once instead of
pairwise.  This round asks what that costs.

Merging k fractions n_i/d_i needs

    numerator  =  sum_i n_i * prod_{j != i} d_j        <-- LEAVE-ONE-OUT SUM
    denominator = prod_j d_j

and at the leaves (n_i = +-1, d_i = mu_i) the numerator is
sum_t +- 2^(Jmax-J_t) prod_{s != t} mu_s -- which is verbatim the numerator
form of V in Open Problem op:E of the paper.

So: collapsing the merge tree to O(1) levels IS Open Problem op:E.  The route
"attack the pipeline's third log by restructuring the merge" is not an
independent attack on the wall; it is the wall.  Verified two ways below.
"""

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from decimal_bbp_extract import bbp_pi_terms, cutoff

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def M(b):
    return b * math.log2(b + 2)


def check_equivalence(K=64, m=512):
    """(A) k-ary merged numerator == the leave-one-out sum == V's numerator."""
    MOD = mpz(1) << m
    mus = [mpz(8 * t + 1) for t in range(1, K + 1)]
    sgn = [mpz(1 if t % 3 else -1) for t in range(1, K + 1)]
    Jmax = 4 * K
    shifts = [Jmax - 4 * t for t in range(1, K + 1)]

    # pairwise tree (the v3 merge, division-free)
    nodes = [(sgn[i], mus[i], shifts[i]) for i in range(K)]
    while len(nodes) > 1:
        nxt = []
        for i in range(0, len(nodes) - 1, 2):
            (nL, dL, sL), (nR, dR, sR) = nodes[i], nodes[i + 1]
            s = sL if sL < sR else sR
            n_ = ((nL * dR) << (sL - s)) + ((nR * dL) << (sR - s))
            nxt.append((n_ % MOD, (dL * dR) % MOD, s))
        if len(nodes) & 1:
            nxt.append(nodes[-1])
        nodes = nxt
    tree_num, tree_den, s0 = nodes[0]

    # direct leave-one-out sum (the one-level k-ary merge)
    P = mpz(1)
    for mu in mus:
        P = P * mu
    direct = mpz(0)
    for i in range(K):
        direct += sgn[i] * (P // mus[i]) * (mpz(1) << (shifts[i] - s0))
    ok1 = (tree_num % MOD) == (direct % MOD)
    ok2 = (tree_den % MOD) == (P % MOD)

    # and it equals V, the additive core, as sum of scaled inverses
    V = sum(sgn[i] * (mpz(1) << (shifts[i] - s0)) * pow(mus[i], -1, MOD)
            for i in range(K)) % MOD
    ok3 = (direct % MOD) * pow(P % MOD, -1, MOD) % MOD == V
    print("(A) k-ary merge numerator == leave-one-out sum : %s"
          % ("EXACT" if ok1 else "MISMATCH"))
    print("    k-ary denominator     == prod mu           : %s"
          % ("EXACT" if ok2 else "MISMATCH"))
    print("    leave-one-out / prod  == V (op:E)          : %s"
          % ("EXACT" if ok3 else "MISMATCH"))
    return ok1 and ok2 and ok3


def cost_of_k_ary(b=64, ks=(2, 4, 8, 16, 32, 64)):
    """(B) Does a k-ary node beat the pairwise subtree it replaces?

    Pairwise subtree over k leaves of b bits: level l has k/2^(l+1) merges of
    3 multiplies at operand size 2^l b.
    Direct k-ary via prefix/suffix products: 3k multiplies, sizes up to k*b.
    """
    print()
    print("(B) cost of ONE k-ary node vs the pairwise subtree it replaces")
    print("    %-6s %-16s %-16s %s" % ("k", "pairwise subtree", "direct k-ary",
                                       "ratio"))
    for k in ks:
        tree = 0.0
        l = 0
        while (1 << l) < k:
            merges = k // (1 << (l + 1))
            sz = (1 << l) * b
            tree += merges * 3 * M(2 * sz)
            l += 1
        # prefix products, suffix products, then one multiply per term
        direct = 0.0
        for i in range(1, k):
            direct += M((i + 1) * b)          # prefix
            direct += M((i + 1) * b)          # suffix
        for i in range(k):
            direct += M(k * b)                # n_i * (prefix*suffix)
        print("    %-6d %-16.3e %-16.3e %.2fx %s"
              % (k, tree, direct, direct / tree,
                 "(k-ary WORSE)" if direct > tree else "(k-ary better)"))
    print("    The pairwise tree IS the efficient leave-one-out algorithm;")
    print("    prefix/suffix collapse costs MORE, not less.  A k-ary node that")
    print("    beat it by a log would, at k = T, be a fast V -- i.e. op:E.")


if __name__ == '__main__':
    print("ROUND 3: removing merge levels == Open Problem op:E")
    ok = check_equivalence()
    cost_of_k_ary()
    print()
    print("EQUIVALENCE", "VERIFIED" if ok else "FAILED")
    print("CONCLUSION: the transform-reuse route is not independent of V.")
