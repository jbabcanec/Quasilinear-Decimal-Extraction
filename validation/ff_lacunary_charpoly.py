#!/usr/bin/env python3
"""
ff_lacunary_charpoly.py -- the sparse Frobenius-lacunary characteristic
polynomial of FF-V's jet sequences (sandbox campaign round 2).

Object: F_i(Y) = sum_{b in {0,1}^i} chi(b)/(Y + w(b)),  chi(b) = prod lambda_t^{b_t},
lambda_t = x^{4*2^t}, w(b) = sum b_t x^{3+t}  (FF-V at Y=1).
Jet sequence at A:  J_n = D^(n) F_i(A) = sum_b chi(b) / (A + w(b))^{n+1}.

CLAIM (derived 2026-07-22, round 2): (J_n) is C-finite with characteristic
polynomial, after clearing by l_W(A) = prod_b (A + w(b)):

    C(S) = S^{2^i} * l_W(A)  +  sum_{t=0..i} c_t * S^{2^i - 2^t},

where l_W(Y) = sum_t c_t Y^{2^t} is the linearized subspace polynomial of
W = span{x^{3+t}}.  Proof: prod_b((A+w(b))S + 1) = S^{2^i} l_W(A + 1/S) and
l_W is additive, so l_W(A + 1/S) = l_W(A) + sum c_t S^{-2^t}.  Sparse: i+2
monomials with Frobenius-gap tap offsets 2^i - 2^t.

Consequence tested here: l_W(A) * J_{n+2^i} = sum_t c_t * J_{n + 2^i - 2^t}.

Significance: the sandbox analog of the Z-side contiguity/BSGS structure; the
order stays 2^i (initial data exponential -- dimension count preserves mass,
consistent with mass conservation), but the taps are Frobenius-lacunary, a
structure with no Z_2 analog. Banked as attack surface for term-jumping
literature (sparse/lacunary recurrence jumping).
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


def subspace_poly(basis, P):
    c = [1]
    for g in basis:
        lg, gg = 0, g & MASK(P)
        for ci in c:
            lg ^= gmul(ci, gg, P)
            gg = gmul(gg, gg, P)
        newc = [gmul(lg, c[0], P)]
        for t in range(1, len(c)):
            newc.append(gmul(c[t - 1], c[t - 1], P) ^ gmul(lg, c[t], P))
        newc.append(gmul(c[-1], c[-1], P))
        c = newc
    return c


def run(i, P, A):
    basis = [1 << (3 + t) for t in range(i)]
    c = subspace_poly(basis, P)              # l_W(Y) = sum c_t Y^{2^t}
    lA, Af = 0, A & MASK(P)
    for ct in c:
        lA ^= gmul(ct, Af, P)
        Af = gmul(Af, Af, P)                 # l_W(A)
    # jet sequence J_n = sum_b chi(b) rho_b^{n+1}, rho_b = 1/(A + w(b))
    nmax = (1 << i) + 6
    J = [0] * (nmax + 1)
    for b in range(1 << i):
        w = 0
        for t in range(i):
            if (b >> t) & 1:
                w ^= 1 << (3 + t)
        chi = 1 << (4 * b)                   # x^{4*int(b)}
        rho = series_inv(A ^ w, P)
        p = rho
        for n in range(nmax + 1):
            J[n] ^= gmul(chi, p, P)
            p = gmul(p, rho, P)
    # test the lacunary recurrence at several n
    ok = True
    for n in range(5):
        lhs = gmul(lA, J[n + (1 << i)], P)
        rhs = 0
        for t in range(i + 1):
            rhs ^= gmul(c[t], J[n + (1 << i) - (1 << t)], P)
        ok &= (lhs == rhs)
    print(f"  i={i} (2^{i}={1<<i} points, recurrence order {1<<i}, "
          f"{i+2} sparse taps): {'EXACT MATCH x5 offsets' if ok else 'MISMATCH'}")
    return ok


if __name__ == "__main__":
    import random
    random.seed(11)
    print("Sparse Frobenius-lacunary characteristic polynomial of jet sequences:")
    allok = True
    for i in [2, 3, 4, 5]:
        A = 1 | (random.getrandbits(6) << 9) | 1
        allok &= run(i, 512, A | 1)
    print("VERDICT:", "VALIDATED -- C(S) = S^{2^i} l(A) + sum c_t S^{2^i-2^t} exact"
          if allok else "FAILED")
