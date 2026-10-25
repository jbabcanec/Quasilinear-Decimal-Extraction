#!/usr/bin/env python3
"""
bet10_window_rigidity.py -- Bet 10: the index-decomposition axis, closed.

MOTIVATION.  Every merge scheme in this project decomposes the index range
{1..K} into a NESTED family of blocks -- a partition refinement, i.e. a tree,
i.e. a cover whose nerve is 0-dimensional.  The merge-mass accounting
(leaf mass x depth) assumes that nesting.  Topology does not: a cover may
overlap, and gluing local data over an overlapping cover is governed by a
Cech H^1 obstruction rather than a sum over levels.  Bet 10 asks whether a
NON-NESTED (signed / overlapping) index decomposition can beat a partition.

FIRST ANSWER (trivial, and it kills the naive form).  H(S) = sum_{k in S}
mu_k^{-1} is finitely ADDITIVE in S.  For an additive functional the Cech
obstruction vanishes identically: inclusion-exclusion over any cover is exact.
So there is never an obstruction, and never a saving either -- you pay for
every intersection.  Cover mass >= partition mass, equality iff partition.
(check F).

SECOND ANSWER (the real content).  The saving would have to come not from
the combinatorics but from the blocks: use blocks whose local sum has a
CLOSED FORM, so their cost is O(1) instead of proportional to their size,
and write the short window as a signed combination of those.  This is the
"virtual cover" version, strictly more general than the union-of-cosets
question the log already noted as dead.

    The only inversion-closed blocks are complete cosets.  1+8Z_2 mod 2^m is
    CYCLIC of order M = 2^{m-3}, so its only subgroups are the filtration
    G_j = 1+2^j Z_2.  Inversion permutes any subgroup, so sum over a coset
    gG_j of mu^{-1} is g^{-1} sum_{h in G_j} h -- an arithmetic series.  In
    index coordinates (mu = 8k+1), a coset of G_j is an AP of step 2^{j-3}
    and length 2^{m-j}: total index span exactly M, ALWAYS.  (checks A,B,C)

THEOREM (window rigidity).  Let M = 2^{m-3} and 0 < K < M/2.  Suppose

    1_[1,K]  =  sum_b c_b 1_{B_b}  +  sum_{k in S} d_k 1_{{k}}

with B_b any inversion-closed blocks (any steps, any offsets, any rational
c_b, any number of terms) and S a finite set of singletons.  Then |S| >= K.

PROOF.  Pass to generating polynomials in Z[z] and reduce mod
Phi_M(z) = z^{M/2} + 1.  A closed block of step s has length M/s >= 2, so
s <= M/2 and s | M/2; its elements pair off as x <-> x + M/2 inside the
block, so its image is 0.  The window has degree K < M/2, so its image is
itself, with K nonzero coefficients.  Folding can only decrease the number of
nonzero coefficients of the singleton part, so |S| >= K.   QED   (check D)

EQUIVALENT FORM (check E).  The Z-span of ALL closed blocks is exactly the
space of M/2-periodic functions on the index line, dimension M/2 -- the
step-(M/2) blocks {a, a+M/2} are already a basis of it, and every other
closed block is a sum of those.  Closed blocks see only the M/2-periodic
part of an index set.  A window of length K < M/2 has none.

CONSEQUENCE.  Complete cosets contribute exactly NOTHING to a short window:
any exact representation still needs K singletons, i.e. K separate
inversions, i.e. the naive algorithm.  No index decomposition -- nested or
overlapping, partition or signed cover -- beats term-by-term.  This upgrades
the log's informal "window << modulus" observation to a theorem, and closes
"a genuinely new decomposition axis" for decompositions OF THE INDEX SET.
The residue in Open Problem E is not a new way to cut up the range.
"""
from __future__ import annotations

import random


def v2(n: int) -> int:
    n = abs(int(n))
    return (n & -n).bit_length() - 1 if n else -1


# --------------------------------------------------- A: the subgroup lattice

def check_A():
    print("=" * 74)
    print("A.  1+8Z_2 mod 2^m is cyclic of order M = 2^(m-3); only subgroups")
    print("    are the filtration G_j = 1+2^j Z_2.  So 'coset' has ONE meaning.")
    print("=" * 74)
    ok = True
    print(f"    {'m':>4} {'M=2^(m-3)':>10} {'ord(9)':>8} {'cyclic?':>8} {'#subgroups':>11}")
    for m in range(6, 17):
        mod = 1 << m
        M = 1 << (m - 3)
        o, x = 1, 9 % mod
        while x != 1:
            x = x * 9 % mod
            o += 1
        # subgroups of a cyclic 2-group of order M: one per divisor, m-2 of them
        nsub = (m - 3) + 1
        good = (o == M)
        ok &= good
        print(f"    {m:>4} {M:>10} {o:>8} {str(good):>8} {nsub:>11}")
    print(f"    -> generator 9 has full order in all cases: {ok}")
    return ok


# ------------------------------------------------------- B: closure of cosets

def check_B():
    print()
    print("=" * 74)
    print("B.  closure: sum of inverses over a complete coset is an arithmetic")
    print("    series (inversion permutes the subgroup).")
    print("=" * 74)
    ok = True
    print(f"    {'m':>4} {'j':>4} {'|coset|':>9}  direct sum == closed form?")
    for m in (10, 12, 14):
        mod = 1 << m
        for j in (4, 6, 8):
            if j >= m:
                continue
            size = 1 << (m - j)
            g = random.choice([x for x in range(1, mod, 8)])
            coset = [(g + (1 << j) * t) % mod for t in range(size)]
            direct = sum(pow(c, -1, mod) for c in coset) % mod
            # closed form: g^{-1} * sum_{h in G_j} h,  G_j = {1 + 2^j u}
            sub = (size + (1 << j) * (size * (size - 1) // 2)) % mod
            closed = pow(g, -1, mod) * sub % mod
            good = direct == closed
            ok &= good
            print(f"    {m:>4} {j:>4} {size:>9}  {good}")
    print(f"    -> closure holds: {ok}")
    return ok


# -------------------------------------------- C: closed blocks are APs, span M

def closed_blocks(m):
    """index sets of all inversion-closed blocks of length >= 2, as (step, offset)."""
    M = 1 << (m - 3)
    out = []
    s = 1
    while s <= M // 2:                 # length M/s >= 2
        for a in range(s):
            out.append((s, a))
        s *= 2
    return out


def check_C():
    print()
    print("=" * 74)
    print("C.  in index coordinates a coset of G_j is an AP of step 2^(j-3) and")
    print("    length 2^(m-j): index span = M ALWAYS, whatever the step.")
    print("=" * 74)
    ok = True
    m = 12
    mod, M = 1 << m, 1 << (m - 3)
    print(f"    m={m}, M={M}")
    print(f"    {'j':>4} {'step':>6} {'length':>8} {'span':>7}  index set is the coset?")
    for j in range(3, m):
        s, L = 1 << (j - 3), 1 << (m - j)
        a = random.randrange(s)
        idx = [a + s * t for t in range(L)]
        mus = sorted((8 * k + 1) % mod for k in idx)
        g = (8 * a + 1) % mod
        coset = sorted((g + (1 << j) * t) % mod for t in range(L))
        good = (mus == coset) and (s * L == M)
        ok &= good
        print(f"    {j:>4} {s:>6} {L:>8} {s*L:>7}  {good}")
    print(f"    -> every closed block spans exactly M indices: {ok}")
    return ok


# ------------------------------------------------------- D: THE THEOREM

def fold(coeffs, half):
    """reduce a Z[z] polynomial (dict deg->coef) mod z^half + 1."""
    out = {}
    for d, c in coeffs.items():
        r, q = d % half, d // half
        out[r] = out.get(r, 0) + (c if q % 2 == 0 else -c)
    return {r: c for r, c in out.items() if c}


def check_D():
    print()
    print("=" * 74)
    print("D.  THEOREM.  mod Phi_M(z) = z^(M/2)+1 every closed block dies and")
    print("    the window survives with all K coefficients.  => |S| >= K.")
    print("=" * 74)
    ok = True
    print(f"    {'m':>4} {'M':>6} {'#blocks':>9} {'all fold to 0?':>15}"
          f" {'K':>5} {'window nonzeros':>16}")
    for m in range(7, 15):
        M = 1 << (m - 3)
        half = M // 2
        blocks = closed_blocks(m)
        allzero = True
        for (s, a) in blocks:
            L = M // s
            poly = {a + s * t: 1 for t in range(L)}
            if fold(poly, half):
                allzero = False
                break
        # digit-regime window: K = Theta(m), certainly < M/2 for m >= 7
        K = min(m, half - 1)
        w = fold({k: 1 for k in range(1, K + 1)}, half)
        good = allzero and len(w) == K
        ok &= good
        print(f"    {m:>4} {M:>6} {len(blocks):>9} {str(allzero):>15}"
              f" {K:>5} {len(w):>16}")
    print()
    print("    Any signed combination of closed blocks folds to 0, so the")
    print("    singleton part alone must reproduce z + ... + z^K after folding.")
    print("    Folding never increases support => at least K singletons.")
    print(f"    -> theorem verified on all blocks, m = 7..14: {ok}")
    return ok


# ------------------------------------- E: span of closed blocks, exactly

def check_E():
    print()
    print("=" * 74)
    print("E.  span(all closed blocks) = {M/2-periodic functions}, dim M/2.")
    print("    The step-(M/2) blocks {a, a+M/2} are already a basis; every")
    print("    other closed block is a sum of those.  Window is not periodic.")
    print("=" * 74)
    ok = True
    print(f"    {'m':>4} {'M':>6} {'every block = sum of basis?':>28} {'window in span?':>17}")
    for m in range(7, 15):
        M = 1 << (m - 3)
        half = M // 2
        good = True
        for (s, a) in closed_blocks(m):
            L = M // s
            elts = [(a + s * t) % M for t in range(L)]   # cyclic index model
            # candidate expansion: sum over residues r of {r, r+half}
            need = sorted(elts)
            built = sorted({r for e in elts for r in (e % half, e % half + half)})
            if need != built:
                good = False
                break
        K = min(m, half - 1)
        win = [1 if 1 <= i <= K else 0 for i in range(M)]
        periodic = all(win[i] == win[i + half] for i in range(half))
        ok &= good and not periodic
        print(f"    {m:>4} {M:>6} {str(good):>28} {str(periodic):>17}")
    print()
    print("    Closed blocks see ONLY the M/2-periodic part of an index set.")
    print("    A window of length K < M/2 has zero periodic content.")
    print(f"    -> span characterisation verified: {ok}")
    return ok


# ------------------------------------ F: the Cech idea is empty for sums

def check_F():
    print()
    print("=" * 74)
    print("F.  the overlapping-cover / Cech idea, for the record: for an")
    print("    ADDITIVE functional the obstruction vanishes identically and")
    print("    inclusion-exclusion strictly increases mass unless it is a")
    print("    partition.  There is nothing to exploit in the combinatorics.")
    print("=" * 74)
    m, K = 14, 40
    mod = 1 << m
    a = {k: pow(8 * k + 1, -1, mod) for k in range(1, K + 1)}
    target = sum(a.values()) % mod
    ok = True

    def window(lo, hi):                       # clamped to [1, K]
        return list(range(max(1, lo), min(K, hi) + 1))

    covers = [
        ("partition (4 blocks)",
         [window(1 + i * K // 4, (i + 1) * K // 4) for i in range(4)]),
        ("overlapping (5 blocks)",
         [window(1 + i * K // 6, i * K // 6 + K // 2) for i in range(5)]),
        ("heavy overlap (8 blocks)",
         [window(1 + i * K // 12, i * K // 12 + 2 * K // 3) for i in range(8)]),
    ]
    print(f"    {'cover':>26} {'#terms':>8} {'mass':>6} {'covers?':>8} {'exact?':>8}")
    for label, cover in covers:
        n = len(cover)
        assert set().union(*map(set, cover)) == set(range(1, K + 1)), label
        tot, mass = 0, 0
        sets = [set(c) for c in cover]
        for bits in range(1, 1 << n):
            inter = None
            for i in range(n):
                if bits >> i & 1:
                    inter = sets[i] if inter is None else inter & sets[i]
            if not inter:
                continue
            sign = -1 if bin(bits).count("1") % 2 == 0 else 1
            tot += sign * sum(a[k] for k in inter)
            mass += len(inter)
        exact = tot % mod == target
        ok &= exact and mass >= K
        print(f"    {label:>26} {n:>8} {mass:>6} {'True':>8} {str(exact):>8}")
    print()
    print(f"    partition mass = K = {K} is the minimum; every cover is exact")
    print("    (H^1 = 0 for an additive functional) and every cover costs more.")
    return ok


if __name__ == "__main__":
    res = [check_A(), check_B(), check_C(), check_D(), check_E(), check_F()]
    print()
    print("=" * 74)
    print(f"VERDICT: {'all checks pass' if all(res) else 'FAILURES PRESENT'}")
    print("Bet 10: index-decomposition axis CLOSED.  Closed blocks span only")
    print("the M/2-periodic functions; a short window needs K singletons no")
    print("matter how the blocks are combined.  The third axis, if it exists,")
    print("is not a way of cutting up the range.")
    print("=" * 74)
