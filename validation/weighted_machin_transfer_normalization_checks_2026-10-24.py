"""Small exact checks for universal normalization of Machin transfer matrices.

Checked 24 October 2026.

This checks algebra and valuation bounds, not a complete extractor or its
running-time exponent. Uses only Python's standard library.
"""
from fractions import Fraction
from math import comb, factorial, gcd, lcm
import json
import random


def primes(n):
    return [p for p in range(2, n + 1)
            if all(p % q for q in range(2, int(p ** 0.5) + 1))]


def valuation(n, p):
    assert n
    out = 0
    while n % p == 0:
        n //= p
        out += 1
    return out


def floor_log(n, p):
    out = 0
    while n >= p:
        n //= p
        out += 1
    return out


def normalizer(s):
    return factorial(s) // lcm(*range(1, s + 1)) if s else 1


def odd_normalizer(s):
    value = normalizer(s)
    return value // (value & -value)


def pqt_merge(left, right):
    lp, lq, lt = left
    rp, rq, rt = right
    return lp*rp, lq*rq, lt*rq+lp*rt


def pqt_leaf(A, a, C, j):
    if j == 0:
        return 1, 1, 1
    divisor = gcd(a, C)
    assert divisor == 2
    p, q = (a//divisor)*j, (C//divisor)*(A+a*j)
    return p, q, p


def raw_pqt(A, a, C, lo, hi):
    value = (1, 1, 0)
    for j in range(lo, hi):
        value = pqt_merge(value, pqt_leaf(A, a, C, j))
    return value


def normalized_pqt(A, a, C, lo, hi, counts):
    if hi == lo + 1:
        return pqt_leaf(A, a, C, lo)
    mid = (lo+hi)//2
    left = normalized_pqt(A, a, C, lo, mid, counts)
    right = normalized_pqt(A, a, C, mid, hi, counts)
    l, r = mid-max(lo, 1), hi-mid
    child_factor = odd_normalizer(l)*odd_normalizer(r)
    assert odd_normalizer(l+r) % child_factor == 0
    factor = odd_normalizer(l+r)//child_factor
    value = pqt_merge(left, right)
    assert all(x % factor == 0 for x in value)
    counts['real_pqt_merges'] += 1
    return tuple(x//factor for x in value)


def check_pqt(A, a, p, lo, hi, counts):
    assert A % 2 == 1 and a >= 2 and (a & (a-1)) == 0
    C = p*p+1
    raw = raw_pqt(A, a, C, lo, hi)
    positive = hi-max(lo, 1)
    factor = odd_normalizer(positive)
    assert all(x % factor == 0 for x in raw)
    reduced = normalized_pqt(A, a, C, lo, hi, counts)
    assert reduced == tuple(x//factor for x in raw)
    # Check the represented series independently with exact rationals.
    product, total = Fraction(1), Fraction(0)
    for j in range(lo, hi):
        ratio = Fraction(a*j, C*(A+a*j)) if j else Fraction(1)
        product *= ratio
        total += product
    assert Fraction(reduced[0], reduced[1]) == product
    assert Fraction(reduced[2], reduced[1]) == total
    counts['real_pqt_intervals'] += 1


def compose(after, before):
    a, b, c, d = after
    x, y, z, w = before
    return a*x, a*y+b*z, c*z, d*w


def leaf(r, D, k):
    return -k, -(r-k+1), -D*(r-k+1), D*k


def raw_transfer(r, D, lo, hi):
    result = (1, 0, 1, 1)
    for k in range(lo + 1, hi + 1):
        result = compose(leaf(r, D, k), result)
    return result


def normalized_transfer(r, D, lo, hi):
    if hi == lo:
        return 1, 0, 1, 1
    if hi == lo + 1:
        return leaf(r, D, hi)
    mid = (lo + hi) // 2
    left = normalized_transfer(r, D, lo, mid)
    right = normalized_transfer(r, D, mid, hi)
    value = compose(right, left)
    divisor = normalizer(hi-lo) // (normalizer(mid-lo)*normalizer(hi-mid))
    assert all(x % divisor == 0 for x in value)
    return tuple(x // divisor for x in value)


def coefficient_row(r, D):
    C = D + 1
    e = [(C**r - 1) // D]
    b = [-1]
    for k in range(1, r + 1):
        b.append((-1)**(k+1) * comb(r, k))
        numerator = b[k] - e[-1]
        assert numerator % D == 0
        e.append(numerator // D)
    assert e[-1] == 0
    return e, b


def check_case(r, p, lo, hi, counts):
    D = p*p
    s = hi-lo
    raw = raw_transfer(r, D, lo, hi)
    factor = normalizer(s)
    assert all(x % factor == 0 for x in raw)
    reduced = normalized_transfer(r, D, lo, hi)
    assert reduced == tuple(x // factor for x in raw)
    counts['matrix_intervals'] += 1
    e, b = coefficient_row(r, D)
    a, off, c, denominator = reduced
    assert a*e[lo] + off*b[lo] == denominator*e[hi]
    assert c*b[lo] == denominator*b[hi]
    for q in primes(r + 3):
        if q == p:
            continue
        loss = valuation(denominator, q)
        bound = floor_log(r, q) + floor_log(s, q)
        assert loss <= bound <= 2*floor_log(r, q)
        counts['valuation_bounds'] += 1
        # One sparse-division propagation with enough logarithmic padding.
        target = q**2
        source = q**(2 + 2*floor_log(r, q))
        wx = (a*(e[lo] % source) + off*(b[lo] % source)) % source
        wy = c*(b[lo] % source) % source
        common = gcd(denominator, source)
        assert wx % common == wy % common == 0
        assert (source // common) % target == 0
        unit = denominator // common
        inverse = pow(unit, -1, target)
        assert (wx // common)*inverse % target == e[hi] % target
        assert (wy // common)*inverse % target == b[hi] % target
        counts['sparse_divisions'] += 1


def main():
    counts = dict(matrix_intervals=0, valuation_bounds=0,
                  sparse_divisions=0, merge_factors=0,
                  real_identities=0, lucas_prefixes=0,
                  real_pqt_intervals=0, real_pqt_merges=0)
    for left in range(34):
        for right in range(34):
            s = left + right
            denominator = normalizer(left)*normalizer(right)
            assert normalizer(s) % denominator == 0
            reconstructed = 1
            for q in primes(s):
                exponent = 0
                z = q
                while z <= s:
                    term = s//z-left//z-right//z+(z<=left)+(z<=right)-1
                    assert term >= 0
                    exponent += term
                    z *= q
                reconstructed *= q**exponent
            assert reconstructed == normalizer(s)//denominator
            counts['merge_factors'] += 1
    for p in (3, 5, 239):
        for r in range(1, 17):
            for lo in range(r):
                for hi in range(lo + 1, r + 1):
                    check_case(r, p, lo, hi, counts)
    rng = random.Random(20261025)
    for _ in range(200):
        r = rng.randrange(17, 97)
        lo = rng.randrange(r)
        hi = rng.randrange(lo + 1, r + 1)
        check_case(r, rng.choice((3, 5, 13, 239)), lo, hi, counts)
    for p in (3, 5, 239):
        C, c = p*p+1, (p*p+1)//2
        for r in range(1, 13):
            e, _ = coefficient_row(r, p*p)
            for A in (1, 3, 17, 101):
                for a in (2, 4, 8):
                    term = total = Fraction(1)
                    for j in range(1, r):
                        term *= Fraction(a*j, C*(A+a*j))
                        total += term
                    partial = C**(r-1)-sum(
                        (Fraction(a*k*e[k], A+a*k) for k in range(1, r)),
                        Fraction(0))
                    assert partial == C**(r-1)*total
                    smaller = c**(r-1)*total
                    denominator = lcm(*(A+a*k for k in range(1, r)))
                    assert denominator % smaller.denominator == 0
                    assert smaller.denominator % 2 == 1
                    counts['real_identities'] += 1
    for _ in range(400):
        p = rng.choice((3, 5, 239))
        D, C = p*p, p*p+1
        r = rng.randrange(1, 180)
        q = rng.choice([q for q in primes(47) if q not in (2, p)])
        k = rng.randrange(r)
        r1, r0 = divmod(r, q)
        h, k0 = divmod(k, q)
        direct = sum(comb(r, j)*pow(D, j, q) for j in range(k+1)) % q
        high = sum(comb(r1, j)*pow(D, j, q) for j in range(h)) % q
        weight = comb(r1, h)*pow(D, h, q) % q
        low = sum(comb(r0, j)*pow(D, j, q)
                  for j in range(min(k0, r0)+1)) % q
        grouped = (pow(C, r0, q)*high + weight*low) % q
        assert direct == grouped
        e, _ = coefficient_row(r, D)
        formula = (-1)**k*pow(pow(D, k+1, q), -1, q)*(pow(C, r, q)-direct) % q
        assert formula == e[k] % q
        counts['lucas_prefixes'] += 1
    for p in (3, 5, 239):
        for a in (2, 4, 8, 16):
            for hi in range(1, 17):
                check_pqt(1, a, p, 0, hi, counts)
                check_pqt(101, a, p, 7, 7+hi, counts)
    for _ in range(240):
        lo = rng.randrange(32)
        hi = lo+rng.randrange(1, 33)
        check_pqt(2*rng.randrange(5000)+1, 1 << rng.randrange(1, 7),
                  rng.choice((3, 5, 13, 239)), lo, hi, counts)
    print(json.dumps(dict(status='PASS', **counts), indent=2))


if __name__ == '__main__':
    main()
