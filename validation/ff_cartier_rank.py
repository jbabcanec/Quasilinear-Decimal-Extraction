#!/usr/bin/env python3
"""
ff_cartier_rank.py -- GF(2)-RANK of the 2-kernel (Cartier orbit) of the
weighted function-field mirror F = sum_v x^{4 int(v)}/(1 + x^3 pol(v)).

This targets the Christol/Cartier criterion DIRECTLY: F is algebraic over
F_2(x) iff the F_2-span of its 2-kernel {a[r::2^k]} is finite-dimensional.
The earlier probe counted DISTINCT kernel rows (maximal growth); distinct rows
could in principle still span a small space, so here we Gaussian-eliminate and
report the rank per level, against an algebraic control.

Rank growing without bound (tracking the row-count budget min(2^k, P/2^k))
= the numerical shadow of the transcendence theorem (attack queue #6);
the formal proof route is the pole/birth-filtration argument in
notes/ff_sandbox.md.
"""

P = 1 << 15


def stream_inv(taps_int, deg, length):
    s = bytearray(length)
    s[0] = 1
    hist, tapmask, mask = 1, taps_int >> 1, (1 << deg) - 1
    for e in range(1, length):
        b = (hist & tapmask).bit_count() & 1
        s[e] = b
        hist = ((hist << 1) | b) & mask
    return s


def build_F(P):
    a = bytearray(P)
    v = 0
    while 4 * v < P:
        q = 1 | (v << 3)
        L = P - 4 * v
        s = stream_inv(q, q.bit_length() - 1, L)
        base = 4 * v
        for e in range(L):
            a[base + e] ^= s[e]
        v += 1
    return a


def pack(bits):
    val = 0
    for b in reversed(bits):
        val = (val << 1) | b
    return val


def gf2_rank(rows):
    rank = 0
    pivots = []
    for r in rows:
        for p in pivots:
            if r.bit_length() == p.bit_length():
                r ^= p
        while r:
            placed = False
            for p in pivots:
                if r.bit_length() == p.bit_length():
                    r ^= p
                    placed = True
                    break
            if not placed:
                pivots.append(r)
                rank += 1
                break
    return rank


def rank_profile(a, P, kmax, label):
    print(f"  {label}:")
    out = []
    for k in range(kmax + 1):
        L = P >> k
        if L < 32:
            break
        rows = [pack(list(a[r::1 << k])[:L]) for r in range(1 << k)]
        rk = gf2_rank(rows)
        bound = min(1 << k, L)
        out.append((k, rk, bound))
        print(f"    k={k}: rank {rk} of {1 << k} rows (budget {bound})")
    return out


if __name__ == "__main__":
    a = build_F(P)
    print(f"F built to x^{P}")
    prof = rank_profile(a, P, 9, "F (weighted mirror)")
    ctrl = bytearray(stream_inv(0b1011, 3, P))
    cprof = rank_profile(ctrl, P, 9, "control 1/(1+x+x^3) [algebraic]")
    full = all(rk == bound for _, rk, bound in prof)
    print(f"\nF kernel span at FULL budget rank at every level: {full}")
    print(f"control span bounded: {max(rk for _, rk, _ in cprof)} (must stay O(1))")
    print("VERDICT:", "kernel span dimension saturates its budget -- the Christol"
          " criterion's numerical shadow says TRANSCENDENTAL"
          if full else "span dimension sub-maximal -- inspect before claiming")
