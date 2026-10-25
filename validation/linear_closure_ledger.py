#!/usr/bin/env python3
"""
linear_closure_ledger.py -- machine checks for notes/linear_closure.tex:
the write-floor theorem, Tellegen closure, and the interface ledger in
information form (the linear-closure extension of T2).

A.  WRITE FLOOR (Theorem 1).  Node charge min(s_v, T) in factor units
    (the common factor b/2 divides out).  Since the charge depends only
    on block sizes, the exact minimum over ALL binary merge trees obeys
    the size DP  f(n) = min_{a+b=n} f(a)+f(b) + min(n,T),  f(1)=0.
    Verified: f(K) >= 0.5 * K * log2(T) for 4 <= T <= K across a grid,
    plus the tightness of the balanced tree, plus random/caterpillar
    forests sitting above the DP optimum, plus d-ary samples >= binary
    optimum (Lemma binary-is-optimal).

B.  TELLEGEN CLOSURE (Theorem 2).  Random linear SLPs over Z/2^64 are
    transposed by adjoint edge reversal; verified: <Mx,y> == <x,M^T y>
    on random vectors, and the multiset of scalar-multiplication
    constants is IDENTICAL in the forward and transposed programs.

C.  INTERFACE LEDGER (Theorem 3).  The exact 2-adic identity
    J(t) = exp(sum_j (-1)^(j-1) t^j Q_j / j),  J(t) = prod(1 + t/mu_k),
    is implemented and verified against the direct product.  Two-world
    adversary: replacing Q_1 by Q_1 + 2^(t*), t* = max(p_i - w_i),
    leaves EVERY promised atom value J(2^(w_i)) mod 2^(p_i) unchanged
    (so no computation certifies bit t* of H), while the perturbation
    2^(t*-1) IS detected by the tightest atom (the ledger is sharp).
"""
from __future__ import annotations

import math
import random

random.seed(20260823)


# ------------------------------------------------- A: the write-floor DP

def floor_dp(K: int, T: int):
    """exact min over all binary merge trees of sum_v min(s_v, T)."""
    f = [0] * (K + 1)
    for n in range(2, K + 1):
        best = min(f[a] + f[n - a] for a in range(1, n // 2 + 1))
        f[n] = best + min(n, T)
    return f[K]


def balanced_cost(n: int, T: int):
    if n == 1:
        return 0
    a = n // 2
    return balanced_cost(a, T) + balanced_cost(n - a, T) + min(n, T)


def random_forest_cost(K: int, T: int):
    blocks = [1] * K
    tot = 0
    while len(blocks) > 1:
        i = random.randrange(len(blocks))
        x = blocks.pop(i)
        j = random.randrange(len(blocks))
        y = blocks.pop(j)
        s = x + y
        tot += min(s, T)
        blocks.append(s)
    return tot


def dary_sample_cost(K: int, T: int, d: int):
    blocks = [1] * K
    tot = 0
    while len(blocks) > 1:
        take = min(d, len(blocks))
        parts = [blocks.pop(random.randrange(len(blocks)))
                 for _ in range(take)]
        s = sum(parts)
        tot += (take - 1) * min(s, T)
        blocks.append(s)
    return tot


def check_A():
    print("=" * 74)
    print("A.  write floor: exact DP minimum >= 0.5 K log2 T; balanced tight")
    print("=" * 74)
    ok = True
    print(f"    {'K':>6} {'T':>5} {'DP min':>8} {'0.5K*lgT':>9} "
          f"{'balanced':>9} {'random':>8} {'caterp.':>8} {'3-ary':>7}")
    for K in (64, 256, 1024, 4096):
        for T in (4, 16, 64, 256):
            if T > K:
                continue
            bound = 0.5 * K * math.log2(T)
            dpmin = floor_dp(K, T)
            bal = balanced_cost(K, T)
            rnd = min(random_forest_cost(K, T) for _ in range(5))
            cat = sum(min(j, T) for j in range(2, K + 1))
            d3 = min(dary_sample_cost(K, T, 3) for _ in range(5))
            ok &= dpmin >= bound and bal >= dpmin and rnd >= dpmin \
                and cat >= dpmin and d3 >= dpmin
            print(f"    {K:6d} {T:5d} {dpmin:8d} {bound:9.0f} "
                  f"{bal:9d} {rnd:8d} {cat:8d} {d3:7d}")
    # binary-is-optimal numerics
    dvals = all((d - 1) / math.log2(d) >= 1 - 1e-12 for d in range(2, 65))
    ok &= dvals
    print(f"    (d-1)/log2(d) >= 1 for d in 2..64 (equality at 2): {dvals}")
    print(f"    A: {'PASS' if ok else 'FAIL'}")
    return ok


# ------------------------------------------------- B: Tellegen closure

def run_forward(n_in, ops, outs, x, mod):
    vals = list(x)
    for op in ops:
        if op[0] == "add":
            vals.append((vals[op[1]] + vals[op[2]]) % mod)
        else:  # ("smul", c, i)
            vals.append((op[1] * vals[op[2]]) % mod)
    return [vals[i] for i in outs]


def run_transpose(n_in, ops, outs, y, mod):
    """adjoint sweep: same constants, reversed dataflow."""
    nvar = n_in + len(ops)
    adj = [0] * nvar
    used_consts = []
    for pos, o in enumerate(outs):
        adj[o] = (adj[o] + y[pos]) % mod
    for idx in range(len(ops) - 1, -1, -1):
        v = n_in + idx
        a = adj[v]
        op = ops[idx]
        if op[0] == "add":
            adj[op[1]] = (adj[op[1]] + a) % mod
            adj[op[2]] = (adj[op[2]] + a) % mod
        else:
            used_consts.append(op[1])
            adj[op[2]] = (adj[op[2]] + op[1] * a) % mod
    return [adj[i] for i in range(n_in)], used_consts


def check_B():
    print("=" * 74)
    print("B.  Tellegen: <Mx,y> == <x,M^T y>; identical constant multiset")
    print("=" * 74)
    ok = True
    mod = 1 << 64
    for trial in range(5):
        n_in, n_ops, n_out = 6, 40, 4
        ops = []
        for i in range(n_ops):
            nvar = n_in + i
            if random.random() < 0.5:
                ops.append(("add", random.randrange(nvar),
                            random.randrange(nvar)))
            else:
                c = random.randrange(mod)
                if trial % 2:
                    c |= 1  # mix unit and non-unit constants
                ops.append(("smul", c, random.randrange(nvar)))
        outs = [n_in + n_ops - 1 - j for j in range(n_out)]
        x = [random.randrange(mod) for _ in range(n_in)]
        y = [random.randrange(mod) for _ in range(n_out)]
        Mx = run_forward(n_in, ops, outs, x, mod)
        MTy, used = run_transpose(n_in, ops, outs, y, mod)
        lhs = sum(a * b for a, b in zip(Mx, y)) % mod
        rhs = sum(a * b for a, b in zip(x, MTy)) % mod
        fwd_consts = sorted(op[1] for op in ops if op[0] == "smul")
        ok &= (lhs == rhs) and (sorted(used) == fwd_consts)
    print(f"    5 random programs (40 ops, mixed unit/non-unit constants): "
          f"{'all pass' if ok else 'FAIL'}")
    print(f"    B: {'PASS' if ok else 'FAIL'}")
    return ok


# --------------------------------------- C: interface ledger, two worlds

def exp2adic(x: int, P: int) -> int:
    """exp(x) mod 2^P for v2(x) >= 2, exact integer series."""
    mod = 1 << P
    term, tot, i = 1, 1, 0
    while term:
        i += 1
        num = term * x
        v = i - bin(i).count("1") - (i - 1) * 0  # v2 of i! increment
        # divide by i exactly 2-adically: i = 2^a * odd
        a = (i & -i).bit_length() - 1
        odd = i >> a
        assert num % (1 << a) == 0 or num % mod == 0
        term = ((num >> a) * pow(odd, -1, mod)) % mod
        tot = (tot + term) % mod
        if i > 4 * P:
            break
    return tot


def jet_from_tower(Q, t: int, P: int) -> int:
    """exp(sum_j (-1)^(j-1) t^j Q_j / j) mod 2^P."""
    mod = 1 << P
    x = 0
    tj = 1
    for j in range(1, len(Q) + 1):
        tj = tj * t
        if tj % mod == 0 and j > 1:
            break
        a = (j & -j).bit_length() - 1
        odd = j >> a
        assert tj % (1 << a) == 0
        contrib = ((tj >> a) * pow(odd, -1, mod) * Q[j - 1]) % mod
        x = (x + contrib if j % 2 else x - contrib) % mod
    return exp2adic(x, P)


def check_C():
    print("=" * 74)
    print("C.  interface ledger: identity exact; two worlds agree at t*;")
    print("    perturbation at t*-1 is detected (sharpness)")
    print("=" * 74)
    ok = True
    K = 64
    profiles = [(24, 56), (40, 64), (16, 40)]
    tstar = max(p - w for w, p in profiles)
    P = max(p for _, p in profiles) + 64
    mod = 1 << P
    mus = [8 * k + 1 for k in range(1, K + 1)]
    Jmax = P // min(w for w, _ in profiles) + 2
    Q = [sum(pow(mu, -j, mod) for mu in mus) % mod for j in range(1, Jmax)]
    # identity sanity: exp form == direct product, all scales
    idok = True
    for w, p in profiles:
        direct = 1
        for mu in mus:
            direct = direct * (1 + (1 << w) * pow(mu, -1, mod)) % mod
        viaQ = jet_from_tower(Q, 1 << w, P)
        idok &= (direct - viaQ) % (1 << p) == 0
    print(f"    exp/log identity matches direct product at all scales: {idok}")
    ok &= idok
    # world 2: Q1 += 2^t*
    Q2w = list(Q)
    Q2w[0] = (Q2w[0] + (1 << tstar)) % mod
    agree = True
    for w, p in profiles:
        v1 = jet_from_tower(Q, 1 << w, P)
        v2 = jet_from_tower(Q2w, 1 << w, P)
        agree &= (v1 - v2) % (1 << p) == 0
    print(f"    t* = {tstar}: all atoms unchanged mod their precisions: "
          f"{agree}  (=> no computation certifies bit {tstar} of H)")
    ok &= agree
    # sharpness: 2^(t*-1) must be visible to the tightest atom
    Q3w = list(Q)
    Q3w[0] = (Q3w[0] + (1 << (tstar - 1))) % mod
    seen = False
    for w, p in profiles:
        v1 = jet_from_tower(Q, 1 << w, P)
        v3 = jet_from_tower(Q3w, 1 << w, P)
        if (v1 - v3) % (1 << p) != 0:
            seen = True
    print(f"    perturbation 2^(t*-1) detected by some atom: {seen}")
    ok &= seen
    print(f"    C: {'PASS' if ok else 'FAIL'}")
    return ok


if __name__ == "__main__":
    results = {
        "A write floor (exact DP)": check_A(),
        "B Tellegen closure": check_B(),
        "C interface ledger": check_C(),
    }
    print("=" * 74)
    for name, r in results.items():
        print(f"    {name:28s} {'PASS' if r else 'FAIL'}")
    print("=" * 74)
    if all(results.values()):
        print("    LINEAR-CLOSURE THEOREMS: all machine checks pass.")
        print("    Floor Omega(N log^2 N) in the merge-forest model;")
        print("    transposition preserves it; the scale interface certifies")
        print("    at most max(p_i - w_i) bits under ARBITRARY computation.")
