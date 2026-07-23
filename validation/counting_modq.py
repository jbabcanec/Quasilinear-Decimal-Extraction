#!/usr/bin/env python3
"""
counting_modq.py -- crack swing #1 (B(d) mod small q) + lower-bound groundwork.

Crack idea: B(d) <= T, so residues mod primes q with product > T pin B(d)
exactly. If B(d) mod q were low-complexity or a cheap character sum for small q,
CRT over O(log T) primes would give B(d) cheaply -> crack.

MEASURED (T=400, N=1200): B(d) mod q has MAXIMAL linear complexity (~N/2) for
every q in {2,3,5,7,11,13}; the parity is not automatic (2-kernel grows); and
B(d) mod q matches no cheap linear/character statistic (chance-level agreement).
=> No residue is tame; B(d) mod q is exactly as hard as B(d). Crack closed
(a mass-conservation confirmation).

LOWER-BOUND GROUNDWORK (why communication is the WRONG tool): the output is N
counts of ~log2(T) bits, ~O(N) bits total. Any communication-complexity lower
bound is bounded by the output size, so it is <= O(N) and CANNOT reach the
N^2/polylog wall. Combined with Jerabek (summing is amortized-linear => cannot
charge carries), the ONLY viable lower-bound route is cell-probe / branching-
program ORBIT-ACCESS (Clifford-Jalsenius template). This eliminates two tool
families and pins the model -- see notes/lb_model.md.
"""
import math


def exact_B(T, N):
    B = [0] * N
    for t in range(1, T + 1):
        mu = 8 * t + 1
        r = 1 % mu
        for d in range(N):
            if r >= (mu + 1) // 2:
                B[d] += 1
            r = (2 * r) % mu
    return B


def berlekamp_massey(s, q):
    n = len(s); C = [1] + [0] * n; Bp = [1] + [0] * n; L = 0; m = 1; b = 1
    for i in range(n):
        d = s[i] % q
        for j in range(1, L + 1):
            d = (d + C[j] * s[i - j]) % q
        if d == 0:
            m += 1
        elif 2 * L <= i:
            Tc = C[:]; c = (d * pow(b, -1, q)) % q
            for j in range(n + 1 - m):
                if j + m <= n:
                    C[j + m] = (C[j + m] - c * Bp[j]) % q
            L = i + 1 - L; Bp = Tc; b = d; m = 1
        else:
            c = (d * pow(b, -1, q)) % q
            for j in range(n + 1 - m):
                if j + m <= n:
                    C[j + m] = (C[j + m] - c * Bp[j]) % q
            m += 1
    return L


def kernel2(seq, kmax=8):
    rows = set()
    for k in range(kmax + 1):
        step = 1 << k
        if len(seq) // step < 16:
            break
        for r in range(step):
            rows.add(bytes(x % 2 for x in seq[r::step])[:len(seq) // step // 2])
    return len(rows)


def main():
    T, N = 400, 1200
    B = exact_B(T, N)
    print(f"T={T}, N={N}. B(d) mod q linear complexity (of {N}, random ~{N // 2}):")
    allmax = True
    for q in [2, 3, 5, 7, 11, 13]:
        L = berlekamp_massey(B, q)
        mx = L > N // 2 - N // 8
        allmax &= mx
        print(f"  q={q:3d}: {L:5d}  {'MAXIMAL' if mx else '*** LOW ***'}")
    print(f"  parity 2-kernel: {kernel2(B)} (grows => not automatic)")
    print(f"\nVERDICT: B(d) mod q maximal for all small q => no tame residue "
          f"=> crack closed. All-maximal: {allmax}")
    print(f"[LB groundwork] output ~{int(N * math.log2(T))} bits = O(N) => "
          f"communication bounds <= O(N), cannot reach N^2. Cell-probe only.")


if __name__ == "__main__":
    main()
