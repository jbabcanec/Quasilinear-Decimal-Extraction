#!/usr/bin/env python3
"""
carlitz_control.py -- the function-field control experiment for Open
Problem E: the atom has NO function-field counterpart.

QUESTION (from the 2026-08-23 external sweep).  Over Z the atom is

    H = sum_{k=1..K} (8k+1)^{-1} mod 2^m,

and every kill in the campaign traces back to one mismatch: the window
[1,K] is an ARCHIMEDEAN interval while all closed structure mod 2^m is
2-ADIC (cosets, balls).  Control: move to F_2[T], where the natural
bounded index set "deg a < d" IS an ultrametric ball (an F_2-subspace).
If the analogue collapses there, the hardness of E is pinned to the seam
between the two absolute values of Z -- Z has an archimedean place,
F_2[T] does not.

ANALOGUE.  Place 2 -> place T, 8 = 2^3 -> T^3, window -> V = {a : deg a
< d} (an F_2-vector space of dimension d, |V| = 2^d = "K"):

    H_FF = sum_{a in V} (1 + T^3 a)^{-1}  mod T^mm.

COLLAPSE (classical Carlitz; Goss ch. 1, Thakur).  For an F_q-subspace V,

    l_V(z) = prod_{a in V} (z + a)

is an ADDITIVE (q-linearized) polynomial: l_V(z) = sum_i c_i z^(2^i),
only d+1 monomials, built by d Ore steps

    l_{i+1}(z) = l_i(z)^2 + l_i(T^i) * l_i(z).

In char 2 the formal derivative of z^(2^i) vanishes for i >= 1, so
l_V'(z) = c_0, and the log-derivative identity gives the CLOSED FORM

    sum_{a in V} 1/(z0 + a) = l_V'(z0)/l_V(z0) = c_0 / l_V(z0)

for any z0 of negative valuation (z0 = u T^{-e}: then z0 + a is never
zero and every term is honestly invertible in F_2((T))).  With
z0 = T^{-3}:  H_FF = T^{-3} c_0 / l_V(T^{-3}) -- computable in O(d) ring
multiplications, i.e. polylog in the window size 2^d.  The 2^d-term sum
costs d multiplications.

READING.  The function-field world solves E instantly because there the
window IS a ball, products over balls are 2-sparse additive polynomials,
and the derivative of an additive polynomial is a constant.  This is the
same phenomenon as bet 10's "complete coset sums are trivially closed",
seen from the other side of the mirror: over F_2[T] there IS no other
kind of bounded window.  Over Z the interval is not a ball, no additive
polynomials exist (Frobenius is not Z-linear), and the campaign's
rigidity/conservation theorems say no substitute structure exists.  E is
a phenomenon of the archimedean/ultrametric seam; it cannot even be
STATED in F_2[T].  (Independent echoes from the sweep: Igusa/cell-
decomposition closed forms stop exactly at the definable = cell = coset
boundary [Cluckers, Zuniga-Galindo]; the exact Lerch/Glaisher window
evaluations exist only at windows that are fixed FRACTIONS of the
modulus, i.e. completion-compatible ones.)

Checks: A Ore ladder (sparsity, additivity, product identity); B closed
form vs brute at generic negative-valuation points; C the H_FF atom
itself, with operation counts.
"""
from __future__ import annotations

import random

random.seed(20260823)

# ---------------------------------------------------------------- F2[T] core
# polynomial = python int, bit i = coeff of T^i; all ops carryless


def cmul(a: int, b: int) -> int:
    r = 0
    while b:
        low = b & -b
        r ^= a << (low.bit_length() - 1)
        b ^= low
    return r


def cinv(f: int, mm: int) -> int:
    """inverse of f with f(0)=1 mod T^mm (char-2 Newton)."""
    assert f & 1
    g, t = 1, 1
    while t < mm:
        t = min(2 * t, mm)
        e = (cmul(f, g) ^ 1) & ((1 << t) - 1)
        g = (g ^ cmul(g, e)) & ((1 << t) - 1)
    return g


# ------------------------------------------------------- Laurent in 1/T land
# element = (val, unit): T^val * unit, unit odd, kept to prec bits


class Laur:
    __slots__ = ("val", "u")

    def __init__(self, val: int, u: int, prec: int):
        while u and not (u & 1):
            u >>= 1
            val += 1
        if u:
            self.val, self.u = val, u & ((1 << prec) - 1)
        else:
            self.val, self.u = 10 ** 9, 0

    def is_zero(self):
        return self.u == 0

    def series_mod(self, mm: int):
        """as element of F2[[T]] mod T^mm; None if negative valuation."""
        if self.is_zero():
            return 0
        if self.val < 0:
            return None
        return (self.u << self.val) & ((1 << mm) - 1)


def lmul(x: Laur, y: Laur, prec: int) -> Laur:
    if x.is_zero() or y.is_zero():
        return Laur(0, 0, prec)
    return Laur(x.val + y.val, cmul(x.u, y.u), prec)


def ladd(x: Laur, y: Laur, prec: int) -> Laur:
    if x.is_zero():
        return y
    if y.is_zero():
        return x
    v = min(x.val, y.val)
    u = (x.u << (x.val - v)) ^ (y.u << (y.val - v))
    return Laur(v, u, prec)


def linv(x: Laur, prec: int) -> Laur:
    return Laur(-x.val, cinv(x.u, prec), prec)


# ------------------------------------------------- A: Ore ladder for l_V(z)

def ore_ladder(d: int):
    """coefficients {i: c_i in F2[T]} of l_V(z) = sum c_i z^(2^i),
    V = {deg < d}; returns (coeffs, ring-op count)."""
    coeffs = {0: 1}          # l_0(z) = z
    ops = 0
    for i in range(d):
        s = 0                # s = l_i(T^i)
        for j, c in coeffs.items():
            s ^= c << (i * (1 << j))
            ops += 1
        new = {}
        for j, c in coeffs.items():
            new[j + 1] = new.get(j + 1, 0) ^ cmul(c, c)      # l_i(z)^2
            ops += 1
            new[j] = new.get(j, 0) ^ cmul(s, c)               # s * l_i(z)
            ops += 1
        coeffs = new
    return coeffs, ops


def eval_ladder(coeffs, z0: Laur, prec: int):
    """(l_V(z0), op count) by Frobenius towers on z0."""
    den = Laur(0, 0, prec)
    ops = 0
    zp = z0
    j_prev = 0
    for j in sorted(coeffs):
        for _ in range(j - j_prev):
            zp = lmul(zp, zp, prec)
            ops += 1
        j_prev = j
        den = ladd(den, lmul(Laur(0, coeffs[j], prec), zp, prec), prec)
        ops += 1
    return den, ops


def check_A():
    print("=" * 74)
    print("A.  Ore ladder: l_V is additive and (d+1)-sparse; O(d) ops")
    print("=" * 74)
    ok = True
    for d in (2, 4, 6, 8):
        coeffs, ops = ore_ladder(d)
        sparse = (len(coeffs) == d + 1 and sorted(coeffs) == list(range(d + 1)))
        z0 = random.getrandbits(6) | 1

        def ev(z):
            r = 0
            for j, c in coeffs.items():
                zp = z
                for _ in range(j):
                    zp = cmul(zp, zp)
                r ^= cmul(c, zp)
            return r
        brute = 1
        for a in range(1 << d):
            brute = cmul(brute, z0 ^ a)
        xs, ys = random.getrandbits(8), random.getrandbits(8)
        addv = ev(xs ^ ys) == (ev(xs) ^ ev(ys))
        ok &= sparse and (brute == ev(z0)) and addv
        print(f"    d={d:2d}  sparse(d+1 terms)={sparse}  "
              f"product-id={brute == ev(z0)}  additive={addv}  ring-ops={ops}")
    print(f"    A: {'PASS' if ok else 'FAIL'}")
    return ok


# --------------------------- B: closed form c0/l(z0) vs brute inverse sums

def check_B():
    print("=" * 74)
    print("B.  sum_(a in V) 1/(z0+a) = c0 / l_V(z0) at generic z0 = u T^{-e}")
    print("=" * 74)
    ok = True
    mm = 64
    for d, e in ((3, 1), (5, 2), (7, 1), (9, 4)):
        u = (random.getrandbits(5) << 1) | 1
        prec = (e + 3) * (1 << d) + mm + 128
        # brute: 1/(z0+a) = T^e / (u + T^e a)
        acc = 0
        for a in range(1 << d):
            acc ^= cinv((u ^ (a << e)) & ((1 << (mm + 64)) - 1), mm + 64)
        brute = Laur(e, acc, mm + 64).series_mod(mm)
        coeffs, _ = ore_ladder(d)
        z0 = Laur(-e, u, prec)
        den, _ = eval_ladder(coeffs, z0, prec)
        form = lmul(Laur(0, coeffs[0], prec), linv(den, prec),
                    prec).series_mod(mm)
        hit = (form is not None) and form == brute
        ok &= hit
        print(f"    d={d:2d} e={e}  K=2^d={1 << d:4d} terms  "
              f"closed form == brute: {hit}")
    print(f"    B: {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------------------- C: the atom itself, H_FF, in O(d) ring ops

def check_C():
    print("=" * 74)
    print("C.  H_FF = sum_(deg a<d) (1+T^3 a)^{-1} = T^{-3} c0 / l_V(T^{-3})")
    print("    -- the E-atom analogue collapses to O(d) ring ops")
    print("=" * 74)
    ok = True
    mm = 64
    for d in (3, 5, 7, 9, 11):
        brute = 0
        for a in range(1 << d):
            brute ^= cinv((1 ^ (a << 3)) & ((1 << mm) - 1), mm)
        prec = 6 * (1 << d) + mm + 128
        coeffs, ops = ore_ladder(d)
        z0 = Laur(-3, 1, prec)
        den, ops2 = eval_ladder(coeffs, z0, prec)
        h = lmul(lmul(Laur(0, coeffs[0], prec), z0, prec),
                 linv(den, prec), prec)
        form = h.series_mod(mm)
        hit = (form is not None) and form == brute
        ok &= hit
        print(f"    d={d:2d}  K=2^d={1 << d:5d}  match mod T^{mm}: {hit}   "
              f"ring-ops ~ {ops + ops2 + 2} (vs {1 << d} inversions brute)")
    print(f"    C: {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    results = {
        "A Ore ladder / additivity": check_A(),
        "B ball inverse-sum closed form": check_B(),
        "C H_FF collapses (the control)": check_C(),
    }
    print("=" * 74)
    for name, r in results.items():
        print(f"    {name:34s} {'PASS' if r else 'FAIL'}")
    print("=" * 74)
    if all(results.values()):
        print("    CONTROL RESULT: the E-atom has no function-field counterpart.")
        print("    Over F_2[T] the window is a ball, the ball-product is a")
        print("    (d+1)-term additive polynomial, and the inverse sum is")
        print("    c0/l(z0): polylog.  E lives on the archimedean/2-adic seam.")
