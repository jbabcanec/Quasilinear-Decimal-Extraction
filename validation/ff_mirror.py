#!/usr/bin/env python3
"""
ff_mirror.py -- the function-field (F_2[x]) mirror of the additive core.

Claim under test (the ELEVENTH FACE / model separation):

  Over F_2[x], the carry-free twin of the atom  H = sum_{k<=K} 1/(8k+1) mod 2^N
  is CLOSED-FORM per dyadic block, hence quasi-linear overall, because dyadic
  index-blocks become F_2-SUBSPACES and the product over a subspace is a
  LINEARIZED (additive) polynomial  l_W(Y) = sum c_i Y^{2^i}:  sparse, with
  constant derivative  l'_W = c_0.  Therefore

      sum_{w in W} 1/(A + w)  =  l'_W(A)/l_W(A)  =  c_0 / l_W(A),

  one inversion per block, log K blocks per arbitrary cutoff K.
  Over Z_2 the same blocks are unit-group cosets only to precision O(log|block|)
  (the complete-block/Dwork O(s) ceiling); F_2[x] has the structure at FULL
  precision.  The separation isolates what the wall is made of.

  PART 1: machine-validate the subspace closed form and the arbitrary-cutoff
          (interval) decomposition, exactly, as power series mod x^P.
  PART 1c: display the valuation collapse (FF blocks cancel to order ~3(2^j-1);
          Z_2 blocks stay near valuation 0) -- the visible fingerprint of carries.
  PART 2: the WEIGHTED mirror (true V-twin)  F = sum_v x^{4 int(v)}/(1 + x^3 pol(v)),
          self-truncating like V.  Hunt a Mahler/Frobenius-linear relation
          sum_i p_i(x) F^{2^i} = p(x)  (exists iff F is algebraic, by Ore/Christol).
          Found relation => weighted mirror is fast too (separation includes weights).
          Empty hunt => geometric-index weights resist even carry-free (evidence).

Exact arithmetic throughout; no floats.  Pure Python, ints as F_2[x] elements.
"""

import random
random.seed(8)

MASK = lambda P: (1 << P) - 1


def gmul(a, b, P=0):
    """Carryless (F_2[x]) product, optionally truncated mod x^P."""
    if a.bit_count() < b.bit_count():
        a, b = b, a
    r = 0
    while b:
        low = b & -b
        r ^= a << (low.bit_length() - 1)
        b ^= low
    if P:
        r &= MASK(P)
    return r


def series_inv(a, P):
    """Inverse of a unit power series (const term 1) mod x^P, Newton in char 2."""
    assert a & 1, "series_inv needs constant term 1"
    v, k = 1, 1
    while k < P:
        k = min(2 * k, P)
        e = gmul(a, v, k) ^ 1          # e = a*v - 1, O(x^{k/2})
        v ^= gmul(v, e, k)             # v <- v + v*e ; then a*v = 1 + e^2
    return v


def subspace_poly(basis, P):
    """Linearized poly of span(basis): coefficients c[i] of Y^{2^i}, mod x^P.
    Recursion: l_{i+1}(Y) = l_i(Y)^2 + l_i(g_i) * l_i(Y).  Returns list c."""
    c = [1]                            # l_0(Y) = Y
    for g in basis:
        lg = 0                         # l_i(g_i) mod x^P
        gg = g & MASK(P)
        for ci in c:
            lg ^= gmul(ci, gg, P)
            gg = gmul(gg, gg, P)       # Frobenius: g^(2^{i+1})
        newc = [gmul(lg, c[0], P)]
        for i in range(1, len(c)):
            newc.append(gmul(c[i - 1], c[i - 1], P) ^ gmul(lg, c[i], P))
        newc.append(gmul(c[-1], c[-1], P))
        c = newc
    return c


def eval_linearized(c, A, P):
    """l(A) = sum c_i A^{2^i} mod x^P."""
    r, Af = 0, A & MASK(P)
    for ci in c:
        r ^= gmul(ci, Af, P)
        Af = gmul(Af, Af, P)
    return r


def pol(n):
    """Bijection N -> F_2[x]: binary digits become coefficients (the int map)."""
    return n


def blocks_of_interval(M):
    """[0, M) as disjoint dyadic blocks (base, j): base + {0..2^j-1}, disjoint supports."""
    out = []
    for b in range(M.bit_length() - 1, -1, -1):
        if (M >> b) & 1:
            base = (M >> (b + 1)) << (b + 1)
            out.append((base, b))
    return out


# ---------------------------------------------------------------- PART 1
def part1(P=2048):
    print("== PART 1: subspace closed form  sum_{w in W} 1/(A+w) = c_0/l_W(A) ==")
    ok = True
    for j in [2, 3, 4, 5, 6, 7]:
        w0 = random.getrandbits(5)
        A = 1 ^ (pol(w0) << 3) ^ (random.getrandbits(4) << (j + 4))
        A |= 1
        basis = [1 << (3 + i) for i in range(j)]        # x^3 * V_j
        direct = 0
        for v in range(1 << j):
            direct ^= series_inv(A ^ (pol(v) << 3), P)
        c = subspace_poly(basis, P)
        closed = gmul(c[0], series_inv(eval_linearized(c, A, P) | 1, P), P) \
            if (eval_linearized(c, A, P) & 1) else None
        # l_W(A) has constant term prod(A const)=1 since every A+w is a unit
        lA = eval_linearized(c, A, P)
        assert lA & 1, "l_W(A) must be a unit"
        closed = gmul(c[0], series_inv(lA, P), P)
        match = direct == closed
        ok &= match
        print(f"  j={j}: 2^{j}={1<<j} terms  vs  1 inversion + {j+1} sparse evals : "
              f"{'EXACT MATCH' if match else 'MISMATCH'}  "
              f"(block valuation v_x = {val(direct)})")
    print("  PART 1:", "PASS" if ok else "FAIL")
    return ok


def val(s):
    return (s & -s).bit_length() - 1 if s else None


def part1b(P=2048, K=0b11010110):
    print(f"== PART 1b: arbitrary cutoff K={K}: interval = {len(blocks_of_interval(K+1))} "
          "closed-form blocks ==")
    direct = 0
    for u in range(K + 1):
        direct ^= series_inv(1 ^ (pol(u) << 3), P)
    closed = 0
    for base, j in blocks_of_interval(K + 1):
        A = 1 ^ (pol(base) << 3)
        basis = [1 << (3 + i) for i in range(j)]
        c = subspace_poly(basis, P)
        lA = eval_linearized(c, A, P)
        assert lA & 1
        closed ^= gmul(c[0], series_inv(lA, P), P)
    match = direct == closed
    print(f"  {K+1} direct inversions == sum of {len(blocks_of_interval(K+1))} "
          f"closed blocks: {'EXACT MATCH' if match else 'MISMATCH'}")
    print("  PART 1b:", "PASS" if match else "FAIL")
    return match


def part1c():
    print("== PART 1c: valuation collapse (FF) vs carry diffusion (Z_2) ==")
    print("   FF block sums S_j = sum_{v<2^j} 1/(1+x^3 v):    v_x(S_j) ~ 3(2^j - 1)")
    P = 2048
    for j in [3, 5, 7]:
        s = 0
        for v in range(1 << j):
            s ^= series_inv(1 ^ (pol(v) << 3), P)
        print(f"     j={j}:  v_x = {val(s)}   [3(2^j-1) = {3*((1<<j)-1)}]")
    print("   Z_2 block sums H_j = sum_{k<2^j} 1/(8k+1) mod 2^2048:")
    M = 1 << 2048
    for j in [3, 5, 7]:
        h = 0
        for k in range(1 << j):
            h = (h + pow(8 * k + 1, -1, M)) % M
        v2 = (h & -h).bit_length() - 1 if h else None
        print(f"     j={j}:  v_2 = {v2}   (stays O(log): no cancellation -- carries diffuse)")
    print("   (contrast displayed; no assertion)")
    return True


# ---------------------------------------------------------------- PART 2
def weighted_F(P, shift=4, dshift=3):
    """F = sum_v x^{shift*int(v)} / (1 + x^dshift * pol(v))  mod x^P (self-truncating)."""
    F = 0
    vmax = P // shift + 1
    for v in range(vmax + 1):
        t = gmul(1 << (shift * v), series_inv(1 ^ (pol(v) << dshift), P), P)
        F ^= t
    return F & MASK(P)


def mahler_hunt(F, P, imax=5, degmax=24, label="F"):
    """Hunt sum_{i<=imax} p_i(x) F^{2^i} + p(x) = 0, deg p_i <= degmax, over F_2.
    Columns: x^e * F^(2^i) and x^e * 1, as bit-vectors of length P.  Nullspace by
    Gaussian elimination on column-combination bitmasks."""
    cols, tags = [], []
    Fi = F
    for i in range(imax + 1):
        for e in range(degmax + 1):
            cols.append(gmul(1 << e, Fi, P))
            tags.append((i, e))
        Fi = gmul(Fi, Fi, P)
    for e in range(degmax + 1):
        cols.append((1 << e) & MASK(P))
        tags.append(("one", e))
    n = len(cols)
    # row-reduce the P x n system: represent each column with an identity tail
    vecs = [(cols[i] << n) | (1 << i) for i in range(n)]
    pivots = {}
    null = []
    for v in vecs:
        while v >> n:
            p = (v >> n).bit_length() - 1
            if p in pivots:
                v ^= pivots[p]
            else:
                pivots[p] = v
                v = 0
                break
        if v:
            null.append(v & MASK(n))
    print(f"  [{label}] hunt: {n} columns x {P} rows -> {len(null)} nullspace vector(s)")
    rels = []
    for nv in null[:4]:
        terms = [tags[i] for i in range(n) if (nv >> i) & 1]
        rels.append(terms)
        print(f"    relation support: {terms[:12]}{' ...' if len(terms) > 12 else ''}")
    return rels


def part2():
    print("== PART 2: weighted mirror  F = sum x^{4 int(v)}/(1 + x^3 pol(v)) ==")
    P = 1024
    # positive control: g = sum x^{2^i} satisfies g^2 + g + x = 0
    g = 0
    for i in range(11):
        if (1 << i) < P:
            g ^= 1 << (1 << i)
    rels = mahler_hunt(g, P, imax=2, degmax=2, label="control+ (sum x^2^i)")
    ctrl_pos = len(rels) > 0
    # negative control: sum x^{i^2} (non-automatic => transcendental, Christol)
    h = 0
    i = 0
    while i * i < P:
        h ^= 1 << (i * i)
        i += 1
    rels = mahler_hunt(h, P, imax=5, degmax=24, label="control- (sum x^i^2)")
    ctrl_neg = len(rels) == 0
    print(f"  controls: positive {'FOUND' if ctrl_pos else 'MISSED'}, "
          f"negative {'clean' if ctrl_neg else 'FALSE POSITIVE'}")
    # the object
    F = weighted_F(P)
    rels = mahler_hunt(F, P, imax=5, degmax=24, label="F weighted mirror")
    if rels:
        print("  candidate relation(s) found -- verifying at doubled precision...")
        P2 = 2048
        F2 = weighted_F(P2)
        verified = []
        for terms in rels:
            acc = 0
            Fi = {0: F2}
            cur = F2
            for i in range(1, 6):
                cur = gmul(cur, cur, P2)
                Fi[i] = cur
            for t in terms:
                i, e = t
                if i == "one":
                    acc ^= (1 << e) & MASK(P2)
                else:
                    acc ^= gmul(1 << e, Fi[i], P2)
            # relation must hold beyond the hunted window (allow tail slack degmax+P)
            holds = (acc & MASK(P2 - 64)) == 0
            verified.append(holds)
            print(f"    relation verified to x^{P2-64}: {holds}")
        print("  PART 2 verdict: ALGEBRAIC (weighted mirror falls carry-free)"
              if any(verified) else
              "  PART 2 verdict: relations spurious -> treat as NO RELATION")
        return any(verified), ctrl_pos and ctrl_neg
    else:
        print("  PART 2 verdict: NO Mahler relation at this scale "
              "(evidence: geometric-index weights resist even carry-free)")
        return False, ctrl_pos and ctrl_neg


if __name__ == "__main__":
    ok1 = part1()
    ok1b = part1b()
    part1c()
    alg, ctrls = part2()
    print()
    print("SUMMARY")
    print(f"  subspace closed form (FF atom, unweighted): {'VALIDATED' if ok1 and ok1b else 'FAILED'}")
    print(f"  controls sane: {ctrls}")
    print(f"  weighted mirror algebraic: {alg}")
    print()
    print("Reading: the carry-free twin of the ATOM (arbitrary-cutoff reciprocal sum)")
    print("is closed-form per block -- the entire Z_2 wall for the unweighted atom is")
    print("carry-arithmetic (torsion-limited block structure), since F_2[x] with")
    print("full-precision additive polynomials collapses it.")
