#!/usr/bin/env python3
"""
door_eps_zeta.py -- doors epsilon (completion / 2-adic Fourier) and zeta
(geometric blocks / dlog intervals): the two loopholes bet 10 left open.

Bet 10's window-rigidity theorem covers signed covers by INVERSION-CLOSED
blocks (= cosets) with RATIONAL coefficients.  Two families slip its
hypotheses, and this script measures both.

DOOR ZETA (geometric blocks).  In dlog coordinates the atom is a pure
geometric: with u = dlog_9(1+8x) on Z/M, M = 2^(m-3),

    (8x+1)^{-1} = 9^{-u},   so   H = sum_{u in U(K)} 9^{-u},

U(K) = dlog image of the window [1,K].  Any AP in u-space (any odd step
lam, any offset, any length) is a closed-form block: its H-mass is a
geometric series, cost O(polylog).  These blocks are NOT cosets and NOT
inversion-closed -- the folding proof does not touch them.  Question: does
U(K) compress into few u-APs?

Lower bound instrument.  For a single common step lam, representing 1_U by
AP-indicators of step lam (+ singletons) needs at least runs_lam(U) blocks,
where runs_lam(U) = number of maximal lam-APs inside U -- concretely, runs
of the dilated set lam^{-1}U under consecutive adjacency.  (Each maximal
run forces a block boundary; signed coefficients cannot cancel a jump.
Polynomial or telescoping weights on intervals do not evade this: a
piecewise-closed-form representation of a 0/1 function still needs a piece
boundary at every run boundary, and each boundary evaluation is one term
of naive mass.)

EXACT LAW (derived, then machine-checked here).  Adjacency u, u-r both in
U (r = lam^{-1} mod M) pulls back through the isometry to the value side:

    runs_lam(U) = K - |W cap 9^r W|,   W = 1+8[1,K] mod 2^m,

the overlap of the value-window with its multiplicative translate by
w = 9^r.  A translate w = 1+8t (t >= 1, dlog odd) meets W in ~K/w points,
so the best STRUCTURED base is w = 9 itself: runs_1(U) = K - floor((K-1)/9)
= (8/9)K.  The measured min over ALL bases sits slightly below that,
exactly tracking the extreme-value floor of a RANDOM K-subset scanned over
M/2 dilations (check C prints both) -- i.e. the shaving beyond 8K/9 is
birthday noise, not structure, and the min stays Theta(K).

MIXED BASES (check C2).  A cover may mix steps.  The jump argument does
not bound mixed covers, but geometric escape does: a lam-AP of length L
inside U is a geometric progression v, vw, vw^2, ..., w = 9^{lam^{-1}},
trapped in the value window 1+8[1,K]; since powers of 9 grow, a trapped
orbit has L <= log_9(8K) + O(1) (plus rare modular luck).  So EVERY
geometric block -- any base -- covers O(log K) window elements, and any
mixed signed cover needs >= K/O(log K) blocks.  Best case, door zeta saves
one log factor off naive Theta(K) mass: Theta(N^2/log N), far above even
the BGS rung N^{3/2}, let alone the target.  Door zeta closes by an exact
combinatorial law plus geometric escape, not by equidistribution.

DOOR EPSILON (completion / cyclotomic coefficients).  Fourier-completing
the window with M-th roots of unity z (irrational coefficients -- outside
bet 10's scope) trades H for the completed transforms

    ghat(a) = sum_{x mod M} (8x+1)^{-1} z^{ax}   in  Z/2^m [z]/(z^{M/2}+1).

Archimedean precedent (Gauss/Salie/stationary phase, Hiary's cascades)
says complete sums mod p^m collapse to O(1) closed forms for ODD p.  At
the wild prime p = 2 there is no such collapse: check D measures the
support of ghat(a) in the fold ring and finds every ODD frequency fully
DENSE, indistinguishable from the transform of a random value-function
(coefficients are pairwise differences of odd values, so all even --
parity-forced, no concentration).  The only sparse frequencies are the
2-power-degenerate ones (a = M/4 has support 2): exactly the coset
frequencies bet 10's theorem already covers.  Coherent picture: the
closed-form-accessible part of the spectrum IS the M/2-periodic part, and
the window has none of it.  No stationary phase => completion replaces
one opaque object by M = 2^(m-3) opaque objects.  And even given closed
forms, exact inversion divides by M = 2^(m-3) (the inverse-DFT
denominator; V Vbar = M I), i.e. costs m-3 guard bits -- the same width
doubling as the dual product.  Door epsilon closes on density + count +
denominator, three independent kills.

Frequency-residue reciprocity (noted, not a door): substituting x -> ax
gives ghat_{8k+1}(a) = a * ghat_{8k+a}(1) -- frequency a trades into the
AP residue class.  A symmetry permuting the M objects, not reducing them.

CONSEQUENCE.  The two loopholes of bet 10 are measured shut.  Combined
with bet 10 itself: no decomposition of the index range -- rational or
cyclotomic coefficients, coset or geometric blocks -- beats term-by-term.
"""
from __future__ import annotations

import random

random.seed(20260823)


def v2(n: int) -> int:
    n = abs(int(n))
    return (n & -n).bit_length() - 1 if n else -1


def inv_mod(a: int, mod: int) -> int:
    return pow(a, -1, mod)


def dlog_table(m: int):
    """u -> x(u) = (9^u - 1)/8 on Z/M, M = 2^(m-3), plus inverse table."""
    M = 1 << (m - 3)
    mod = 1 << m
    x_of_u = [0] * M
    t = 1
    for u in range(M):
        x_of_u[u] = ((t - 1) // 8) % M
        t = (t * 9) % mod
    u_of_x = [0] * M
    for u, x in enumerate(x_of_u):
        u_of_x[x] = u
    return M, mod, x_of_u, u_of_x


def cyclic_runs(mask: list, M: int) -> int:
    return sum(1 for u in range(M) if mask[u] and not mask[u - 1])


# ------------------------------------------------ A: dlog coordinates sanity

def check_A():
    print("=" * 74)
    print("A.  dlog coordinates: bijection, value identity, isometry")
    print("=" * 74)
    ok = True
    for m in (8, 10, 12, 14):
        M, mod, x_of_u, u_of_x = dlog_table(m)
        bij = sorted(x_of_u) == list(range(M))
        val = all(
            inv_mod(8 * x_of_u[u] + 1, mod) == pow(inv_mod(9, mod), u, mod)
            for u in random.sample(range(M), min(M, 40))
        )
        iso = all(
            v2((x_of_u[u] - x_of_u[v]) % M) == v2((u - v) % M)
            for u in random.sample(range(M), 20)
            for v in random.sample(range(M), 20)
            if u != v
        )
        # adjacency pullback used by the exact law:  u-1 in U  <->  x in 9[1,K]+1
        adj = all(
            x_of_u[u - 1] == ((x_of_u[u] - 1) * inv_mod(9, M)) % M
            for u in range(1, M)
        )
        ok &= bij and val and iso and adj
        print(f"    m={m:3d}  M={M:5d}  bijection={bij}  value-id={val}"
              f"  isometry={iso}  adjacency={adj}")
    print(f"    A: {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------- B: run count of the dlog window, base lam = 1

def check_B():
    print("=" * 74)
    print("B.  runs(U(K)) = K - floor((K-1)/9) exactly while 9K+1 <= M")
    print("    (each run boundary is one forced block; Theta(K) always)")
    print("=" * 74)
    ok = True
    print(f"    {'m':>4} {'M':>6} {'K':>6} {'runs':>6} {'K-floor((K-1)/9)':>17} "
          f"{'exact?':>7} {'runs/K':>7}")
    for m in (10, 12, 14, 16):
        M, mod, x_of_u, u_of_x = dlog_table(m)
        for K in (m, 4 * m, 16 * m, M // 8, M // 4):
            if not (0 < K < M // 2):
                continue
            mask = [False] * M
            for x in range(1, K + 1):
                mask[u_of_x[x]] = True
            r = cyclic_runs(mask, M)
            pred = K - (K - 1) // 9
            exact = (9 * K + 1 <= M)
            hit = (r == pred) if exact else True
            ok &= hit and r >= K // 3
            print(f"    {m:4d} {M:6d} {K:6d} {r:6d} {pred:17d} "
                  f"{('yes' if r == pred else ('n/a' if not exact else 'NO')):>7} "
                  f"{r / K:7.3f}")
    print(f"    B: {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------- C: base-change scan -- min over ALL geometric-block bases

def check_C():
    print("=" * 74)
    print("C.  min over all odd steps lam of runs_lam(U): the exact law")
    print("    runs_lam(U) = K - |W cap 9^r W| (r = lam^{-1}); best base is 9")
    print("=" * 74)
    ok = True
    # exact-regime cases: 9K+1 <= M, the digit regime K = Theta(m) << M
    for m, K in ((10, 12), (12, 48), (13, 100)):
        M, mod, x_of_u, u_of_x = dlog_table(m)
        U = [u_of_x[x] for x in range(1, K + 1)]
        Uset = set(U)

        best, best_lam, total = None, None, 0
        n_lam = 0
        law_ok = True
        maxchain, maxchain_lam = 0, None
        for lam in range(1, M, 2):
            r = inv_mod(lam, M)
            # runs_lam(U) = #{u in U : u - r not in U}; also longest chain
            runs = 0
            for u in U:
                if (u - r) % M not in Uset:
                    runs += 1
                    L, t = 1, (u + r) % M
                    while t in Uset:
                        L += 1
                        t = (t + r) % M
                    if L > maxchain:
                        maxchain, maxchain_lam = L, lam
            # exact law: K - |W cap 9^r W| with W = 1+8[1,K] mod 2^m
            w = pow(9, r, mod)
            ov = 0
            for x in range(1, K + 1):
                y = ((1 + 8 * x) * w) % mod
                if y % 8 == 1 and 1 <= (y - 1) // 8 <= K:
                    ov += 1
            law_ok &= (runs == K - ov)
            total += runs
            n_lam += 1
            if best is None or runs < best:
                best, best_lam = runs, lam

        # random baseline: same scan on random K-subsets of Z/M
        rbest_all = []
        for _ in range(3):
            S = set(random.sample(range(M), K))
            rb = min(
                sum(1 for u in S if (u - inv_mod(lam, M)) % M not in S)
                for lam in range(1, M, 2)
            )
            rbest_all.append(rb)
        rbest = sum(rbest_all) / len(rbest_all)

        import math
        pred = K - (K - 1) // 9
        esc = math.log(8 * K + 1) / math.log(9)
        print(f"    m={m:3d} K={K:4d}: base-9 value {pred}; min over bases = "
              f"{best} at lam={best_lam}; mean {total / n_lam:.1f}; "
              f"random-set min ~ {rbest:.1f}")
        print(f"      exact law runs = K - overlap: "
              f"{'holds for all lam' if law_ok else 'VIOLATED'}")
        print(f"      C2 longest chain any base: {maxchain} (at lam="
              f"{maxchain_lam}); geometric-escape cap log_9(8K+1) = {esc:.1f}")
        ok &= law_ok and best >= (2 * K) // 3
        ok &= maxchain <= 3 * esc + 3
    print(f"    C: {'PASS' if ok else 'FAIL'}  (min over bases = Theta(K),")
    print("       tracking the random extreme-value floor; longest block O(log K)")
    print("       => any mixed cover needs >= K/O(log K) blocks; door zeta shut)")
    return ok


# --------------------- D: completed transforms are DENSE (no stationary phase)

def fold_transform(vals, a: int, M: int, mod: int):
    """sum_x vals[x] z^(ax) in Z/mod [z]/(z^(M/2)+1); returns coeff list."""
    half = M // 2
    c = [0] * half
    for x in range(M):
        e = (a * x) % M
        if e < half:
            c[e] = (c[e] + vals[x]) % mod
        else:
            c[e - half] = (c[e - half] - vals[x]) % mod
    return c


def check_D():
    print("=" * 74)
    print("D.  door epsilon: ghat(a) dense at all odd frequencies, sparse only")
    print("    at 2-power (coset) frequencies; matches random value-functions")
    print("=" * 74)
    ok = True
    for m in (10, 12):
        M, mod, x_of_u, u_of_x = dlog_table(m)
        g = [inv_mod(8 * x + 1, mod) for x in range(M)]
        grand = [random.randrange(1, mod, 2) for _ in range(M)]
        sample_a = sorted({1, 2, 3, 5, M // 4, M // 2 - 1,
                           random.randrange(1, M, 2), random.randrange(1, M, 2)})
        print(f"    m={m}: M={M}, fold-ring dim {M // 2}")
        print(f"    {'a':>6} {'nnz(g)':>8} {'units(g)':>9} "
              f"{'nnz(rand)':>10} {'units(rand)':>12}")
        for a in sample_a:
            cg = fold_transform(g, a, M, mod)
            cr = fold_transform(grand, a, M, mod)
            nnz_g = sum(1 for c in cg if c)
            uni_g = sum(1 for c in cg if c % 2)
            nnz_r = sum(1 for c in cr if c)
            uni_r = sum(1 for c in cr if c % 2)
            print(f"    {a:6d} {nnz_g:8d} {uni_g:9d} {nnz_r:10d} {uni_r:12d}")
            if a % 2 == 1:  # full-support frequencies
                ok &= nnz_g >= (M // 2) * 9 // 10
        # reciprocity spot check: ghat_{8k+1}(a) = a * ghat_{8k+a}(1), a odd unit
        a = 5 if M > 8 else 3
        left = fold_transform(g, a, M, mod)
        h = [(inv_mod(8 * x + a, mod) * a) % mod for x in range(M)]
        right = fold_transform(h, 1, M, mod)
        rec = (left == right)
        print(f"    reciprocity ghat_(8k+1)(a) = a*ghat_(8k+a)(1) at a={a}: {rec}")
        ok &= rec
    print(f"    D: {'PASS' if ok else 'FAIL'}  (no closed forms to complete into;")
    print("       and exact inversion carries denominator M = 2^(m-3) regardless)")
    return ok


if __name__ == "__main__":
    results = {
        "A dlog coordinates": check_A(),
        "B run-count law": check_B(),
        "C base-change scan": check_C(),
        "D completion density": check_D(),
    }
    print("=" * 74)
    for name, r in results.items():
        print(f"    {name:28s} {'PASS' if r else 'FAIL'}")
    print("=" * 74)
    if all(results.values()):
        print("    DOORS EPSILON AND ZETA: measured shut.")
        print("    Geometric blocks: single base saves the factor 8/9 + noise;")
        print("    mixed bases save at most one log (escape law).  Completion:")
        print("    no stationary phase at p=2, and the M-denominator regardless.")
    else:
        print("    UNEXPECTED STRUCTURE FOUND -- investigate before logging.")
