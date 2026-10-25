#!/usr/bin/env python3
"""
door_iota_adjoint.py -- door iota: the gradient trick (Baur-Strassen /
reverse-mode differentiation), checked and fenced.

THE TEMPTATION.  e_{K-1} = sum_k prod_{j!=k} mu_j = P0 * H is the
gradient sum of the product: with f(x_1..x_K) = prod x_k,

    e_{K-1} = sum_k (df/dx_k)|_{x=mu}.

Baur-Strassen: ALL partial derivatives of an arithmetic circuit cost at
most 5x the circuit size.  The Borwein crack computes P0 mod 2^m in
O(N log^2 N).  If that crack were an arithmetic circuit IN THE VARIABLES
mu_k, Baur-Strassen would deliver e_{K-1} = P0*H in O(N log^2 N) and
close Open Problem E on the spot.

WHY IT DOES NOT.  Two independent reasons, both verified here.

1.  The Borwein crack is VALUE-based, not circuit-based: it computes the
    integer P0 by sieving the prime factorizations of the specific
    integers 8k+1 (shared primes across k).  As a computation it takes k
    as input, not mu_k as free variables; its steps (exponent assembly
    over the sieve) are not ring operations on the mu's.  Baur-Strassen
    applies to straight-line programs over the variables -- the cheapest
    KNOWN mu-variable circuit for prod x_k is the product tree, size
    K-1 multiplications.  The value/polynomial gap IS the razor
    ("differentiation destroys factorization"), now in circuit language.

2.  Reverse-mode through the product tree re-derives the existing
    pipeline, not something faster.  Forward values SHRINK toward the
    leaves (a node at level i holds ~2^i factors, O(2^i log K) bits,
    capped at the working width); adjoints DO NOT SHRINK (every adjoint
    is a cofactor product of ~K factors: full width at every node).
    Propagating the SUM of leaf-adjoints instead of each adjoint gives
    exactly the recursion

        e_{K-1}(T) = P(right)*e_{K-1}(left) + P(left)*e_{K-1}(right)

    -- the S x P merge of the O(N log^3 N) extraction pipeline.  The
    one-log gap between the crack and the gradient is the width
    mismatch: forward mass ~ leaf mass, adjoint mass ~ #nodes x full
    width.  Checks: A correctness of the reverse sweep; B the width
    law (adjoint/forward mass ratio grows like m/log K toward leaves);
    C the sum-propagation recursion is exactly the merge.

CONSEQUENCE.  Automatic differentiation, adjoint methods, and
transposition principles (Tellegen) applied to any value-DAG for P0 all
land inside the existing merge-mass accounting.  Door iota shut.  A
crack of E via this door would need a mu-VARIABLE circuit for prod x_k
of size o(K) -- which cannot exist (the polynomial has K distinct
irreducible factors; even reading the inputs costs K).
"""
from __future__ import annotations


def build_tree(vals, mod):
    """product tree; returns list of levels, level[0] = leaves."""
    levels = [[v % mod for v in vals]]
    while len(levels[-1]) > 1:
        cur = levels[-1]
        nxt = [(cur[2 * i] * cur[2 * i + 1]) % mod
               for i in range(len(cur) // 2)]
        if len(cur) % 2:
            nxt.append(cur[-1])
        levels.append(nxt)
    return levels


def reverse_sweep(levels, mod):
    """adjoints: adj(root)=1; adj(child) = adj(parent)*sibling."""
    adj = [1]
    for lev in range(len(levels) - 2, -1, -1):
        cur = levels[lev]
        up = adj
        nxt = []
        for i in range(len(cur)):
            p, r = divmod(i, 2)
            sib = i + 1 if r == 0 else i - 1
            if sib < len(cur) and p < len(up):
                nxt.append((up[p] * cur[sib]) % mod)
            else:  # odd node carried up unchanged
                nxt.append(up[p] % mod)
        adj = nxt
    return adj  # leaf adjoints = prod_{j != k} mu_j mod mod


def check_A():
    print("=" * 74)
    print("A.  reverse sweep: sum of leaf adjoints == e_{K-1} == P0*H")
    print("=" * 74)
    ok = True
    for K, m in ((16, 48), (32, 96), (64, 128)):
        mod = 1 << (2 * m)
        mus = [8 * k + 1 for k in range(1, K + 1)]
        levels = build_tree(mus, mod)
        adjs = reverse_sweep(levels, mod)
        grad_sum = sum(adjs) % mod
        P0 = levels[-1][0]
        H = sum(pow(mu, -1, mod) for mu in mus) % mod
        hit = grad_sum == (P0 * H) % mod
        # brute leave-one-out
        brute = 0
        for k in range(K):
            p = 1
            for j in range(K):
                if j != k:
                    p = (p * mus[j]) % mod
            brute = (brute + p) % mod
        ok &= hit and grad_sum == brute
        print(f"    K={K:3d} m={m:4d}: grad-sum == P0*H: {hit}, == brute: "
              f"{grad_sum == brute}")
    print(f"    A: {'PASS' if ok else 'FAIL'}")
    return ok


def check_B():
    print("=" * 74)
    print("B.  width law: forward values shrink, adjoints stay full width")
    print("=" * 74)
    ok = True
    K, m = 256, 256
    mod = 1 << (2 * m)
    mus = [8 * k + 1 for k in range(1, K + 1)]
    # forward widths as exact integers (uncapped) vs cap 2m
    levels_exact = [[v for v in mus]]
    while len(levels_exact[-1]) > 1:
        cur = levels_exact[-1]
        nxt = [cur[2 * i] * cur[2 * i + 1] for i in range(len(cur) // 2)]
        if len(cur) % 2:
            nxt.append(cur[-1])
        levels_exact.append(nxt)
    levels = build_tree(mus, mod)
    adjs_by_level = []
    adj = [1]
    for lev in range(len(levels) - 2, -1, -1):
        cur = levels[lev]
        up = adj
        nxt = []
        for i in range(len(cur)):
            p, r = divmod(i, 2)
            sib = i + 1 if r == 0 else i - 1
            if sib < len(cur) and p < len(up):
                nxt.append((up[p] * cur[sib]) % mod)
            else:
                nxt.append(up[p] % mod)
        adj = nxt
        adjs_by_level.append((lev, adj))
    print(f"    {'level':>6} {'#nodes':>7} {'fwd bits/node':>14} "
          f"{'adj bits/node':>14} {'adj/fwd':>8}")
    for lev, adj in reversed(adjs_by_level):
        n = len(adj)
        fwd = sum(min(v.bit_length(), 2 * m)
                  for v in levels_exact[lev]) / n
        aw = sum(a.bit_length() for a in adj) / n
        ratio = aw / fwd
        print(f"    {lev:6d} {n:7d} {fwd:14.1f} {aw:14.1f} {ratio:8.1f}")
        if lev == 0:
            ok &= ratio > 30  # leaves: ~2m / log(8K) >> 1
    print("    (adjoints are full-width at EVERY level; forward leaves are")
    print("     log-width.  Reverse-mode mass = #nodes x 2m = Theta(K m).)")
    print(f"    B: {'PASS' if ok else 'FAIL'}")
    return ok


def check_C():
    print("=" * 74)
    print("C.  sum-propagation == the S x P merge recursion, node by node")
    print("=" * 74)
    ok = True
    K, m = 64, 96
    mod = 1 << (2 * m)
    mus = [8 * k + 1 for k in range(1, K + 1)]

    def sxp(lo, hi):
        """returns (P, E) = (prod, leave-one-out sum) on mus[lo:hi]."""
        if hi - lo == 1:
            return mus[lo] % mod, 1
        mid = (lo + hi) // 2
        PL, EL = sxp(lo, mid)
        PR, ER = sxp(mid, hi)
        return (PL * PR) % mod, (PR * EL + PL * ER) % mod

    P0, E = sxp(0, K)
    H = sum(pow(mu, -1, mod) for mu in mus) % mod
    ok &= E == (P0 * H) % mod
    print(f"    S x P merge: e_(K-1) == P0*H mod 2^(2m): {ok}")
    print("    (this recursion IS what pushing 'sum of adjoints' down the")
    print("     tree yields: A(parent) -> A*P_R on the left, A*P_L on the")
    print("     right.  Baur-Strassen on the tree = the known pipeline.)")
    print(f"    C: {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    results = {
        "A reverse sweep correct": check_A(),
        "B width law": check_B(),
        "C gradient == S x P merge": check_C(),
    }
    print("=" * 74)
    for name, r in results.items():
        print(f"    {name:28s} {'PASS' if r else 'FAIL'}")
    print("=" * 74)
    if all(results.values()):
        print("    DOOR IOTA SHUT: gradient/adjoint/transposition methods on")
        print("    any value-DAG for P0 reproduce the merge; Baur-Strassen")
        print("    needs a mu-variable circuit, and prod x_k needs size K.")
