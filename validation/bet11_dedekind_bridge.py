#!/usr/bin/env python3
"""
bet11_dedekind_bridge.py -- Bet 11: the Dedekind / reciprocity-cascade
family, the last live door from the 2026-08-23 external sweep.

WHY THIS FAMILY.  The sweep established that the ONLY objects in the
literature that compute modular-inverse-flavored sums in polylog time are
the Dedekind/cotangent/Vasyunin family: exact reciprocity under the
SL2(Z) action on h/k, Euclidean cascade, O(log k) evaluation
(word-rewriting algorithms; Bettin-Conrey period functions).  p-adic
valued versions exist (Rosen-Snyder / Snyder / Kudo).  Bet 11 asks: can
any member or combination of this family reach the atom

    H = sum_{k=1..K} (8k+1)^{-1} mod 2^m ?

STRUCTURAL CLOSURE (three steps, the actual content):

1.  Every state of the reciprocity cascade is a COMPLETE weighted sum
    sum_{a mod c} w(a/c) with sawtooth/Bernoulli weights.  Reciprocity
    maps complete sums to complete sums; the cascade shrinks the modulus
    and the window TOGETHER.  There is no incomplete member anywhere in
    the family (sweep: "no literature object has that shape").

2.  At the 2-power cusp (modulus 2^t): a complete polynomially-weighted
    sum is in the closed span -- weighted Faulhaber/Eulerian towers over
    complete cosets, killed in section 9 of the log, folded to zero by
    the window-rigidity theorem.  Cascade values at denominator 2^m are
    window-blind.  (Check A verifies the cascade computes them exactly
    in O(log) steps -- they are CHEAP, and cheap here means
    cracked-side, not door-side.)

3.  At odd moduli q = 8k+1: each s_p(h,q) is a rational of height
    O(q^{p+1}) -- an O(log q)-bit object.  A FIXED number of them cannot
    span the Theta(m) fresh bits of H; a Theta(K)-indexed family of them
    is per-k residue data feeding the same merge tree as every other
    weighting (mass law; the section-9 weighted-tower kill again).

WHAT REMAINS FALSIFIABLE: finitely many STRUCTURED couplings -- single
cascade values whose arguments are built from (K, 2^m), e.g.
s(2^m mod q*, q*) at q* = 8K+1, s(q*^{-1} mod 2^m, 2^m), Apostol
variants.  Check C hunts for a low-height Z-relation between H and a
pool of such values, at precisions 2^512 and 2^768 and two window
sizes, by lattice reduction.  The reducer is validated in check B by two
witnesses it MUST find:

  - the true identity  H + 8*(sum k q_k^{-1}) - K = 0  (since
    k q_k^{-1} = (1 - q_k^{-1})/8);
  - a planted fake relation.

EXPECTED (and the honest reading of a negative): the hunt comes back
empty at heights far below the lattice's random floor, uniformly across
the grid => no structured coupling exists, and by 1-3 the family is
closed.  Bet 11 then records: the reciprocity cascade lives entirely on
the complete-sum side of the rigidity boundary; its p-adic
interpolations move in the weight direction (Iwasawa/Mazur), which is
bet 4's dead zone.  The last external family is inside the fence.
"""
from __future__ import annotations

from fractions import Fraction
import math
import random

random.seed(20260823)


# ------------------------------------------------------------ Dedekind sums

def saw(x: Fraction) -> Fraction:
    fx = x - math.floor(x)
    return Fraction(0) if fx == 0 else fx - Fraction(1, 2)


def s_brute(h: int, k: int) -> Fraction:
    return sum((saw(Fraction(a, k)) * saw(Fraction(h * a, k))
                for a in range(1, k)), Fraction(0))


def s_cascade(h: int, k: int):
    """classical Dedekind sum via reciprocity; returns (value, steps)."""
    steps = 0
    sign, acc = 1, Fraction(0)
    h %= k
    while h > 1:
        acc += sign * (Fraction(-1, 4)
                       + Fraction(h * h + k * k + 1, 12 * h * k))
        h, k = k % h, h
        sign = -sign
        steps += 1
    if h == 1:  # s(1,k) closed form
        acc += sign * Fraction((k - 1) * (k - 2), 12 * k)
    return acc, steps


B3 = lambda x: x ** 3 - Fraction(3, 2) * x ** 2 + Fraction(1, 2) * x


def s3_brute(h: int, k: int) -> Fraction:
    tot = Fraction(0)
    for a in range(1, k):
        x = Fraction(h * a, k)
        tot += saw(Fraction(a, k)) * B3(x - math.floor(x))
    return tot


def red2(fr: Fraction, s: int) -> int:
    """2-adic reduction of a rational with odd denominator, mod 2^s."""
    den = fr.denominator
    assert den % 2 == 1, f"even denominator {den}"
    return (fr.numerator * pow(den, -1, 1 << s)) % (1 << s)


# ------------------------------------------- A: cascade exactness + the cusp

def check_A():
    print("=" * 74)
    print("A.  cascade == brute; O(log) steps; 2-power-cusp values are cheap")
    print("=" * 74)
    ok = True
    for h, k in ((5, 7), (7, 12), (31, 97), (89, 144), (100, 253)):
        if math.gcd(h, k) != 1:
            continue
        v, st = s_cascade(h, k)
        ok &= (v == s_brute(h, k))
    print(f"    cascade == brute on coprime pairs: {ok}")
    for h, k in ((5, 8), (11, 30), (89, 144)):
        lhs = s_brute(h, k) + s_brute(k, h)
        rhs = Fraction(-1, 4) + Fraction(h * h + k * k + 1, 12 * h * k)
        ok &= (lhs == rhs)
    print(f"    reciprocity exact on 3 pairs: {ok}")
    print(f"    {'t':>4} {'h':>6} {'match':>6} {'steps':>6}")
    for t in (6, 8, 10):
        k = 1 << t
        for h in (3, k // 4 + 1, k - 3):
            if math.gcd(h, k) != 1:
                continue
            v, st = s_cascade(h, k)
            hit = (v == s_brute(h, k))
            ok &= hit and st <= 2 * t + 2
            print(f"    {t:4d} {h:6d} {str(hit):>6} {st:6d}")
    hbig = pow(105, -1, 1 << 512)
    v, st = s_cascade(hbig, 1 << 512)
    intcheck = (v * 6 * (1 << 512)).denominator == 1
    print(f"    s(105^-1 mod 2^512, 2^512): steps={st}, 6k*s integral={intcheck}")
    ok &= intcheck and st <= 800
    print(f"    A: {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------- lattice reducer

def lll(rows, delta=Fraction(3, 4)):
    """LLL on integer row vectors; exact Fraction GSO with incremental
    swap updates (formulas derived from the GSO invariants; validated by
    check B's witnesses)."""
    b = [list(r) for r in rows]
    n = len(b)
    dim = len(b[0])

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    def full_gso():
        mu = [[Fraction(0)] * n for _ in range(n)]
        B = [Fraction(0)] * n
        bs = []
        for i in range(n):
            v = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = (sum(Fraction(b[i][t]) * bs[j][t]
                                for t in range(dim)) / B[j]
                            if B[j] else Fraction(0))
                v = [v[t] - mu[i][j] * bs[j][t] for t in range(dim)]
            bs.append(v)
            B[i] = sum(x * x for x in v)
        return mu, B

    mu, B = full_gso()
    k = 1
    passes = 0
    while k < n:
        passes += 1
        if passes > 200000:
            print("      (LLL pass cap hit -- partial reduction)")
            break
        # size-reduce b_k against all lower rows
        for l in range(k - 1, -1, -1):
            q = math.floor(mu[k][l] + Fraction(1, 2))
            if q:
                b[k] = [x - q * y for x, y in zip(b[k], b[l])]
                for j in range(l):
                    mu[k][j] -= q * mu[l][j]
                mu[k][l] -= q
        if B[k] < (delta - mu[k][k - 1] ** 2) * B[k - 1]:
            # swap rows k-1, k and update GSO incrementally
            b[k], b[k - 1] = b[k - 1], b[k]
            m0 = mu[k][k - 1]
            Bk1, Bk = B[k - 1], B[k]
            C = Bk + m0 * m0 * Bk1
            for j in range(k - 1):
                mu[k - 1][j], mu[k][j] = mu[k][j], mu[k - 1][j]
            for i in range(k + 1, n):
                t1, t2 = mu[i][k - 1], mu[i][k]
                mu[i][k - 1] = (t2 * Bk + m0 * t1 * Bk1) / C
                mu[i][k] = t1 - m0 * t2
            mu[k][k - 1] = m0 * Bk1 / C
            B[k - 1] = C
            B[k] = Bk1 * Bk / C
            k = max(k - 1, 1)
        else:
            k += 1
    return b


def relation_hunt(vals, s, height_bound, label):
    """vals: residues mod 2^s (vals[0] = target).  Hunt small c with
    sum c_i vals_i == 0 mod 2^s.  Returns (found, (height, c))."""
    n = len(vals)
    W = 1 << 64
    rows = []
    for i in range(n):
        r = [0] * n + [W * (vals[i] % (1 << s))]
        r[i] = 1
        rows.append(r)
    rows.append([0] * n + [W * (1 << s)])
    red = lll(rows)
    small = []
    best = None
    for r in red:
        if any(r[:n]) and r[n] == 0:
            h = max(abs(x) for x in r[:n])
            if h <= height_bound:
                small.append((h, list(r[:n])))
            if best is None or h < best[0]:
                best = (h, list(r[:n]))
    small.sort()
    hstr = f"2^{best[0].bit_length() - 1}" if best else "none"
    print(f"    [{label}] best null-vector height ~ {hstr} "
          f"(bound 2^{height_bound.bit_length() - 1}); "
          f"{len(small)} small null vector(s)")
    for h, c in small:
        chk = sum(ci * vi for ci, vi in zip(c, vals)) % (1 << s)
        print(f"      c = {c}  verifies: {chk == 0}")
    return small, best


# ------------------------------------- B: pipeline validation (two witnesses)

def check_B():
    print("=" * 74)
    print("B.  reducer validation: must find the true identity and a plant")
    print("=" * 74)
    ok = True
    s, K = 256, 12
    mod = 1 << s
    H = sum(pow(8 * k + 1, -1, mod) for k in range(1, K + 1)) % mod
    B7 = sum(k * pow(8 * k + 1, -1, mod) for k in range(1, K + 1)) % mod
    assert (H + 8 * B7 - K) % mod == 0
    small, best = relation_hunt([H, 1, B7], s, 1 << 6, "true identity")
    ok &= bool(small) and best[0] <= 16
    X = random.randrange(mod) | 1
    F = (5 * H + 17 + 3 * X) % mod
    small2, best2 = relation_hunt([H, 1, X, F], s, 1 << 6, "planted fake")
    ok &= bool(small2) and best2[0] <= 32
    print(f"    B: {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------------------ C: the bridge hunt over structured couplings

def pool_for(K: int, s: int):
    mod = 1 << s
    qs = [8 * k + 1 for k in range(1, K + 1)]
    H = sum(pow(q, -1, mod) for q in qs) % mod
    qstar = 8 * K + 1
    P0 = 1
    for q in qs:
        P0 = (P0 * q) % mod
    atoms = {"one": 1, "P0": P0}
    # cusp values at denominator 2^s: integers 6*2^s*s(h,2^s) mod 2^s.
    # NOTE: s(hbar,k) = s(h,k) (inverse-argument symmetry) -- the first
    # run of this hunt REDISCOVERED that theorem from raw residues at
    # height 1, validating the machinery; inverse-argument atoms are
    # excluded as known duplicates.
    for name, h in (("Ds_3", 3),
                    ("Ds_5", 5),
                    ("Ds_q*", qstar),
                    ("Ds_2K1", 2 * K + 1)):
        v, _ = s_cascade(h % mod, mod)
        num = v * 6 * mod
        assert num.denominator == 1
        atoms[name] = num.numerator % mod
    # odd-cusp values at q* and neighbors (2-adic reductions of 2*s(h,q))
    for name, (h, q) in (("Dq*", (pow(2, s, qstar), qstar)),
                         ("Dq*b", (pow(2, s + 1, qstar), qstar)),
                         ("Dq8", (pow(2, s, qstar + 8), qstar + 8)),
                         ("Dq16", (pow(2, s, qstar + 16), qstar + 16))):
        v, _ = s_cascade(h % q, q)
        atoms[name] = red2(2 * v, s)
    # Apostol weight-3 at (2^s mod q*, q*)
    v3 = s3_brute(pow(2, s, qstar) % qstar, qstar)
    atoms["D3q*"] = red2(4 * v3, s)
    return H, atoms


def check_C():
    print("=" * 74)
    print("C.  bridge hunt: c_H*H + sum c_i*atom_i == 0 mod 2^s, small height")
    print("=" * 74)
    ok = True
    for K, s in ((12, 512), (16, 512), (16, 768)):
        H, atoms = pool_for(K, s)
        names = sorted(atoms)
        vals = [H] + [atoms[nm] for nm in names]
        n1 = len(vals) + 1
        floor_bits = (64 + s) // n1
        small, best = relation_hunt(vals, s, 1 << 20, f"K={K}, s={s}")
        print(f"      pool (n={len(vals)}): H, {', '.join(names)}")
        print(f"      random-lattice floor ~ 2^{floor_bits}; "
              f"relation bound 2^20")
        bridge = [sc for sc in small if sc[1][0] != 0]
        internal = [sc for sc in small if sc[1][0] == 0]
        for h, c in internal:
            nz = [(names[i - 1], c[i]) for i in range(1, len(c)) if c[i]]
            print(f"      internal atom identity (no H): {nz} -- "
                  f"investigate/known symmetry")
        if bridge:
            print("      !!! BRIDGE CANDIDATE (involves H) -- check "
                  "stability across grid before believing")
        ok &= not bridge
    print(f"    C: {'PASS (no bridge to H; family closed)' if ok else 'BRIDGE?'}")
    return ok


if __name__ == "__main__":
    results = {
        "A cascade + cusp": check_A(),
        "B reducer validation": check_B(),
        "C bridge hunt": check_C(),
    }
    print("=" * 74)
    for name, r in results.items():
        print(f"    {name:26s} {'PASS' if r else 'FAIL'}")
    print("=" * 74)
    if all(results.values()):
        print("    BET 11 CLOSED: the reciprocity cascade transports only")
        print("    complete-sum data (window-blind at the 2-power cusp by")
        print("    rigidity; O(log q)-bit rationals at odd cusps); no")
        print("    structured coupling to H at any tested precision.")
        print("    The last external family is inside the fence.")
    else:
        print("    ANOMALY -- investigate before logging.")
