#!/usr/bin/env python3
"""
Search for Machin-like identities with "smooth" denominators.
============================================================

Goal: find integer relations of the form

    sum_i a_i * arctan(1 / b_i) = pi/4

where b_i are "smooth" (typically 2^u * 5^v) so that base-10 digit extraction
might avoid the nasty base-conversion + large-prime bottlenecks (e.g., 239).

Key exact criterion (no floating point):

Let theta_b = arctan(1/b). Then:
    exp(2 i theta_b) = (b + i)/(b - i)

Thus the identity is equivalent to:
    Π_i ((b_i + i)/(b_i - i))^{a_i} = i

Write the left side as Num/Den in Gaussian integers. Then the condition is:
    Num == i * Den   (exact, in Z[i])

This script brute-forces small coefficient searches for 2-term and 3-term
relations over a candidate b-set.

If you find any hits, validate numerically (MPFR/decimal) separately before
building an extractor around them.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from typing import Iterable

import argparse


@dataclass(frozen=True)
class G:
    """Gaussian integer a + b i with Python big ints."""
    a: int
    b: int

    def __mul__(self, other: "G") -> "G":
        return G(self.a * other.a - self.b * other.b, self.a * other.b + self.b * other.a)

    def conj(self) -> "G":
        return G(self.a, -self.b)


def g_pow(z: G, e: int) -> G:
    if e < 0:
        raise ValueError("g_pow expects nonnegative exponent")
    r = G(1, 0)
    base = z
    while e:
        if e & 1:
            r = r * base
        e >>= 1
        if e:
            base = base * base
    return r


def i_mul(z: G) -> G:
    """Multiply by i: i*(a+bi) = -b + a i."""
    return G(-z.b, z.a)


def smooth_bs(max_u: int = 10, max_v: int = 10, extra: Iterable[int] = ()) -> list[int]:
    s = set(extra)
    for u in range(max_u + 1):
        for v in range(max_v + 1):
            s.add((2**u) * (5**v))
    s.discard(1)  # arctan(1) is slow series; keep search focused
    return sorted(s)


def check_relation(terms: list[tuple[int, int]]) -> bool:
    """
    terms: [(a_i, b_i), ...]
    Returns True if sum a_i arctan(1/b_i) == pi/4 exactly.
    """
    num = G(1, 0)
    den = G(1, 0)
    for a, b in terms:
        z = G(b, 1)       # b + i
        zc = G(b, -1)     # b - i
        if a >= 0:
            num = num * g_pow(z, a)
            den = den * g_pow(zc, a)
        else:
            aa = -a
            num = num * g_pow(zc, aa)
            den = den * g_pow(z, aa)
    return num == i_mul(den)


def search_two_term(bs: list[int], coeff_bound: int) -> list[list[tuple[int, int]]]:
    hits: list[list[tuple[int, int]]] = []
    for b1, b2 in combinations(bs, 2):
        for a1 in range(-coeff_bound, coeff_bound + 1):
            if a1 == 0:
                continue
            for a2 in range(-coeff_bound, coeff_bound + 1):
                if a2 == 0:
                    continue
                terms = [(a1, b1), (a2, b2)]
                if check_relation(terms):
                    hits.append(terms)
    return hits


def search_three_term(bs: list[int], coeff_bound: int, b_limit: int = 12) -> list[list[tuple[int, int]]]:
    """
    Brute 3-term search; b_limit caps the number of b values used (combinatorics explode).
    """
    hits: list[list[tuple[int, int]]] = []
    bs2 = bs[:b_limit]
    for b1, b2, b3 in combinations(bs2, 3):
        for a1, a2, a3 in product(range(-coeff_bound, coeff_bound + 1), repeat=3):
            if a1 == 0 or a2 == 0 or a3 == 0:
                continue
            terms = [(a1, b1), (a2, b2), (a3, b3)]
            if check_relation(terms):
                hits.append(terms)
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--deep", action="store_true", help="Run a much larger search (can be slow).")
    args = ap.parse_args()

    bs = smooth_bs(max_u=10, max_v=6, extra=(3, 6, 12, 15, 24, 30, 40, 75))
    print(f"Candidate b count: {len(bs)}", flush=True)
    print("First 30 b values:", bs[:30], flush=True)

    # Quick by default; deep mode increases bounds.
    coeff_bound_2 = 25 if args.deep else 10
    b_take_2 = 24 if args.deep else 14
    print(f"\nSearching 2-term relations with |a| <= {coeff_bound_2} ...", flush=True)
    hits2 = search_two_term(bs[:b_take_2], coeff_bound_2)
    for h in hits2[:20]:
        print("HIT 2-term:", " + ".join([f"{a}*atan(1/{b})" for a, b in h]), "= pi/4", flush=True)
    print(f"2-term hits: {len(hits2)}", flush=True)

    coeff_bound_3 = 10 if args.deep else 6
    b_limit_3 = 14 if args.deep else 10
    print(f"\nSearching 3-term relations with |a| <= {coeff_bound_3} ...", flush=True)
    hits3 = search_three_term(bs, coeff_bound_3, b_limit=b_limit_3)
    for h in hits3[:20]:
        print("HIT 3-term:", " + ".join([f"{a}*atan(1/{b})" for a, b in h]), "= pi/4", flush=True)
    print(f"3-term hits: {len(hits3)}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

