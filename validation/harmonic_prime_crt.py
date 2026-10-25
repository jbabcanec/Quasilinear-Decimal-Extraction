"""Prime-power partial fractions for the UNWEIGHTED harmonic atom.

This is a new experimental route, separate from the decimal extractor.
It computes an exact rational sum and then reduces it modulo 2**m.
No large product cofactor is constructed per term. See the companion note.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import time
from fractions import Fraction

from gmpy2 import mpz, invert


def primes_to(bound):
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(bound) + 1):
        if sieve[p]:
            sieve[p*p::p] = b"\x00" * ((bound - p*p)//p + 1)
    return (p for p in range(3, bound + 1, 2) if sieve[p])


def rational_tree(nodes):
    """Unreduced rational sum; leaf denominators are pairwise coprime."""
    if not nodes:
        return mpz(0), mpz(1)
    while len(nodes) > 1:
        nxt = []
        for i in range(0, len(nodes)-1, 2):
            a, b = nodes[i]
            c, d = nodes[i+1]
            nxt.append((a*d + c*b, b*d))
        if len(nodes) & 1:
            nxt.append(nodes[-1])
        nodes = nxt
    return nodes[0]


def harmonic_prime_crt(K, m=None, *, start=9, step=8):
    """sum_{k=0}^{K-1} 1/(start+step*k), exactly, and optionally mod 2^m.

    start must be positive and odd and step a positive power of two.
    Complexity bound in the note assumes start=O(K), step=O(1), m=O(K).
    """
    if K < 1 or start < 1 or start % 2 != 1 or step < 1 or step & (step-1):
        raise ValueError("require K>=1, odd positive start, power-of-two step")
    if step == 1 and K > 1:
        raise ValueError("all denominators must be odd")
    t0 = time.perf_counter()
    bound = start + step*(K-1)
    # Both sums below have <= K+bound terms; the absolute rounding error
    # in their difference is < 1/32, independent of integer boundaries.
    guard = (K + bound).bit_length() + 5
    unit = 1 << guard
    real_h_lower = sum(unit // (start + step*k) for k in range(K))
    real_t_lower = 0
    nodes = []
    incidences = 0
    leaf_bits = 0
    for p in primes_to(bound):
        q = p
        while q <= bound // p:
            q *= p
        first = (-start * pow(step, -1, p)) % p
        c = 0
        for k in range(first, K, p):
            mu = start + step*k
            u, pv = mu, 1
            while u % p == 0:
                u //= p
                pv *= p
            c += (q // pv) * int(invert(u, q))
            c %= q
            incidences += 1
        if c:
            nodes.append((mpz(c), mpz(q)))
            real_t_lower += (c * unit) // q
            leaf_bits += q.bit_length()
    t_sieve = time.perf_counter() - t0
    correction = (real_h_lower - real_t_lower + unit//2) // unit
    t1 = time.perf_counter()
    num, den = rational_tree(nodes)
    num += correction * den
    t_tree = time.perf_counter() - t1
    t1 = time.perf_counter()
    residue = None
    if m is not None:
        if m < 1:
            raise ValueError("m must be positive")
        modulus = mpz(1) << m
        residue = (num * pow(den, -1, modulus)) % modulus
    t_reduce = time.perf_counter() - t1
    stats = dict(K=K, m=m, start=start, step=step,
                 primes=len(nodes), incidences=incidences, leaf_bits=leaf_bits,
                 numerator_bits=num.bit_length(), denominator_bits=den.bit_length(),
                 correction=correction, sieve=t_sieve, tree=t_tree,
                 reduce=t_reduce, total=t_sieve+t_tree+t_reduce)
    return num, den, residue, stats


def harmonic_product_tree(K, m):
    """Reference: ordinary rational binary splitting without reductions."""
    num, den = rational_tree([(mpz(1), mpz(8*k+1)) for k in range(1, K+1)])
    mod = mpz(1) << m
    return (num * pow(den, -1, mod)) % mod


def selftest():
    rng = random.Random(20260906)
    cases = [(k, 9, 8) for k in range(1, 151)]
    cases += [(rng.randrange(1, 90), rng.randrange(1, 100)*2+1,
               1 << rng.randrange(1, 5)) for _ in range(100)]
    for K, start, step in cases:
        m = rng.randrange(1, 400)
        num, den, residue, _ = harmonic_prime_crt(K, m, start=start, step=step)
        ref = sum((Fraction(1, start+step*k) for k in range(K)), Fraction())
        assert num * ref.denominator == den * ref.numerator, (K, start, step)
        mod = 1 << m
        assert residue == (ref.numerator * pow(ref.denominator, -1, mod)) % mod
    for K in [257, 1000, 3000, 10000]:
        _, _, residue, _ = harmonic_prime_crt(K, K)
        assert residue == harmonic_product_tree(K, K)
    # An exact higher-moment identity refutes the unrestricted interface
    # ledger with a free complete moment tower: x^(1+2^(m-3)) == x mod 2^m
    # for x in 1+8 Z_2. It is a model counterexample, not a fast evaluator.
    for m in range(5, 13):
        mod = 1 << m
        j = 1 + (1 << (m-3))
        for K in [1, 9, 17, 100]:
            H = sum(pow(8*k+1, -1, mod) for k in range(1, K+1)) % mod
            Qj = sum(pow(8*k+1, -j, mod) for k in range(1, K+1)) % mod
            assert H == Qj
    print(json.dumps(dict(selftest="PASS", exact_fraction_cases=len(cases),
                          large_tree_cases=4, higher_moment_cases=32)))


def benchmark(sizes):
    for K in sizes:
        _, _, residue, st = harmonic_prime_crt(K, K)
        t0 = time.perf_counter()
        ref = harmonic_product_tree(K, K)
        st["product_tree"] = time.perf_counter() - t0
        st["match"] = bool(residue == ref)
        assert st["match"]
        print(json.dumps(st), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bench", nargs="+", type=int)
    args = ap.parse_args()
    if args.selftest:
        selftest()
    if args.bench:
        benchmark(args.bench)
    if not args.selftest and not args.bench:
        ap.print_help()
