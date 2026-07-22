#!/usr/bin/env python3
"""
ff_fold.py -- FRI/Binius-style subspace FOLDING of FF-V (sandbox round 4).

Key identity (char 2, additive): for the quotient map q(y) = y(y + G0),
    (A + Y)(A + G0 + Y) = q(A) + q(Y),  q additive => q(Y(b')) is again a
    LINEAR FORM of the remaining bits. So the two children of the bit-0 split
    COMBINE over a common denominator that is the level-(l+1) subspace form:
    the branching becomes a single descendant with a transformed numerator.

Fold (exact, termwise):
    sum_b chi(b) N(b) / (A + Y(b))
  = sum_b' chi'(b') N1(b') / (A1 + Y1(b')),
    A1 = q(A), Y1(b') = sum b'_t q(G_t),
    N1(b') = N(0,b')*(A + G0 + Y'(b')) + lam0*N(1,b')*(A + Y'(b')),
    Y'(b') = sum b'_t G_t  (parent forms).

Numerators live in the MULTILINEAR algebra over the bits (b_t^2 = b_t, char-2
cross terms of squared forms vanish). The experiment: fold i levels, at each
level (a) verify EXACT value match against brute force, (b) record the number
of active multilinear monomials and the minimum x-valuation of their
coefficients, with and without mod-x^P truncation. Bounded/truncatable growth
=> quasi-linear sandbox algorithm. Flat-valuation growth => mass conservation
with the sharpest mechanism yet.
"""

MASK = lambda P: (1 << P) - 1


def gmul(a, b, P=0):
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
    assert a & 1
    v, k = 1, 1
    while k < P:
        k = min(2 * k, P)
        e = gmul(a, v, k) ^ 1
        v ^= gmul(v, e, k)
    return v


def val(s):
    return (s & -s).bit_length() - 1 if s else None


def brute_value(i, P, gens, A, lams, N=None):
    """Direct sum over all b: chi(b) * N(b) / (A + Y(b)); N multilinear dict or None=1."""
    tot = 0
    for b in range(1 << i):
        Y = 0
        chi = 1
        for t in range(i):
            if (b >> t) & 1:
                Y ^= gens[t]
                chi = gmul(chi, lams[t], P)
        num = 1
        if N is not None:
            num = 0
            for mono, coef in N.items():
                if mono & ~b:
                    continue
                num ^= coef
            num &= MASK(P)
        tot ^= gmul(gmul(chi, num, P), series_inv(A ^ Y, P), P)
    return tot & MASK(P)


def fold_once(i, P, gens, A, lams, N):
    """One fold on bit 0. N: dict {bitmask over bits 0..i-1: coeff}. Returns new state."""
    G0, lam0 = gens[0], lams[0]
    # split N into N(b0=0,.) and N(b0=1,.) as multilinear dicts over bits 1..i-1
    N0, N1v = {}, {}
    for mono, c in N.items():
        if mono & 1:
            N1v[mono >> 1] = N1v.get(mono >> 1, 0) ^ c
        else:
            N0[mono >> 1] = N0.get(mono >> 1, 0) ^ c
    # N(1, b') = N0' + N1v' (setting b0=1 adds the b0-monomials)
    Nat1 = dict(N0)
    for m, c in N1v.items():
        Nat1[m] = Nat1.get(m, 0) ^ c
    # multiply-by-form helper: F(b') = base + sum_{t} b'_t * gens[t+1]
    def mul_affine(D, base, formgens, P):
        out = {}
        for m, c in D.items():
            if c == 0:
                continue
            out[m] = out.get(m, 0) ^ gmul(c, base, P)
            for t, g in enumerate(formgens):
                mm = m | (1 << t)          # b'_t * monomial (b^2=b)
                out[mm] = out.get(mm, 0) ^ gmul(c, g, P)
        return {m: c & MASK(P) for m, c in out.items() if c & MASK(P)}
    formgens = gens[1:]
    partA = mul_affine(N0, A ^ G0, formgens, P)               # N(0)*(A+G0+Y')
    partB = mul_affine({m: gmul(c, lam0, P) for m, c in Nat1.items()},
                       A, formgens, P)                        # lam0*N(1)*(A+Y')
    Nnew = dict(partA)
    for m, c in partB.items():
        Nnew[m] = Nnew.get(m, 0) ^ c
    Nnew = {m: c for m, c in Nnew.items() if c}
    newgens = [gmul(g, g ^ G0, P) for g in gens[1:]]          # q(G_t)
    newA = gmul(A, A ^ G0, P)                                 # q(A)
    newlams = lams[1:]
    return i - 1, newgens, newA, newlams, Nnew


def run(i, P):
    gens = [1 << (3 + t) for t in range(i)]
    lams = [1 << (4 * (1 << t)) for t in range(i)]
    A = 1
    N = {0: 1}
    ref = brute_value(i, P, gens, A, lams)      # the true F
    print(f"i={i}, P={P}: direct F computed. Folding:")
    lvl_i, lvl_gens, lvl_A, lvl_lams, lvl_N = i, gens, A, lams, N
    ok_all = True
    for step in range(i):
        lvl_i, lvl_gens, lvl_A, lvl_lams, lvl_N = fold_once(
            lvl_i, P, lvl_gens, lvl_A, lvl_lams, lvl_N)
        cur = brute_value(lvl_i, P, lvl_gens, lvl_A, lvl_lams, lvl_N)
        ok = (cur == ref)
        ok_all &= ok
        monos = len(lvl_N)
        maxdeg = max((m.bit_count() for m in lvl_N), default=0)
        vals = sorted(val(c) for c in lvl_N.values())
        # effective monomials if we truncate coefficients below x^P (already done)
        print(f"  after fold {step+1}: {lvl_i} bits left | {'OK' if ok else 'VALUE MISMATCH'}"
              f" | monomials {monos} (max deg {maxdeg})"
              f" | coeff valuations min {vals[0] if vals else '-'}"
              f" med {vals[len(vals)//2] if vals else '-'}"
              f" max {vals[-1] if vals else '-'}")
    final = gmul(next(iter(lvl_N.values())) if lvl_N else 0,
                 series_inv(lvl_A, P), P) if lvl_N else 0
    print(f"  final closed value == direct F: {final == ref}")
    return ok_all and final == ref


if __name__ == "__main__":
    ok = True
    for (i, P) in [(6, 512), (8, 1024), (10, 2048)]:
        ok &= run(i, P)
        print()
    print("VERDICT:", "fold EXACT at every level" if ok else "FOLD BROKEN — fix algebra")
