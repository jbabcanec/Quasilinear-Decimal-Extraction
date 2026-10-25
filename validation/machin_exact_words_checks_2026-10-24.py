"""Small exact checks for a separate proposed Machin word-fraction fallback.

This file is a correctness prototype, not a benchmark or a production extractor.
Only coefficient generation uses bounded chunks; Fraction is the independent oracle.
"""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import random


def primes_to(bound):
    sieve = bytearray(b"\1") * (bound + 1)
    if bound >= 0:
        sieve[0] = 0
    if bound >= 1:
        sieve[1] = 0
    for prime in range(2, int(bound**0.5) + 1):
        if sieve[prime]:
            for multiple in range(prime * prime, bound + 1, prime):
                sieve[multiple] = 0
    return [prime for prime in range(3, bound + 1, 2) if sieve[prime]]


def prime_power(prime, bound):
    power = prime
    while power <= bound // prime:
        power *= prime
    return power


def family(p0, n, cutoff):
    return {"p0": p0, "n": n, "K": cutoff,
            "b": 2 if p0 == 5 else 10,
            "s": n if p0 == 5 else 0,
            "c": 16 if p0 == 5 else 4,
            "sign": 1 if p0 == 5 else -1}


def direct_word(f, mu):
    unit, valuation = mu, 0
    while unit % f["p0"] == 0:
        unit //= f["p0"]
        valuation += 1
    if unit == 1:
        return 0, unit, valuation
    residue = (f["c"] * pow(f["b"], f["n"], unit)
               * pow(f["p0"], f["s"] - mu - valuation, unit)) % unit
    signed = f["sign"] * (-1 if ((mu - 1) // 2) & 1 else 1) * residue
    return signed, unit, valuation


def direct_sum(families):
    total = Fraction(0)
    for f in families:
        for k in range(f["K"] + 1):
            signed, unit, _ = direct_word(f, 2 * k + 1)
            total += Fraction(signed, unit)
    return total % 1


def row_coefficient(f, prime, modulus, chunk_length):
    bound = 2 * f["K"] + 1
    if prime == f["p0"] or prime > bound:
        return 0
    sign = f["sign"] * (-1 if ((prime - 1) // 2) & 1 else 1)
    weight = (sign * f["c"] * pow(f["b"], f["n"], modulus)
              * pow(f["p0"], f["s"] - prime, modulus)) % modulus
    step = -pow(f["p0"], -2 * prime, modulus) % modulus
    coefficient, mu = 0, prime
    while mu <= bound:
        chunk, prefix, product = [], [], 1
        while mu <= bound and len(chunk) < chunk_length:
            unit, power = mu, 1
            while unit % prime == 0:
                unit //= prime
                power *= prime
            chunk.append((unit, power, weight))
            prefix.append(product)
            product = product * unit % modulus
            weight = weight * step % modulus
            mu += 2 * prime
        inverse = pow(product, -1, modulus)
        for index in range(len(chunk) - 1, -1, -1):
            unit, power, weight_at_mu = chunk[index]
            inverse_unit = inverse * prefix[index] % modulus
            inverse = inverse * unit % modulus
            coefficient = (coefficient + (modulus // power)
                           * weight_at_mu * inverse_unit) % modulus
    return coefficient


def exact_fraction(families, chunk_length=7):
    bound = max(2 * f["K"] + 1 for f in families)
    leaves = []
    for prime in primes_to(bound):
        modulus = prime_power(prime, bound)
        coefficient = sum(row_coefficient(f, prime, modulus, chunk_length)
                          for f in families) % modulus
        if coefficient:
            leaves.append((coefficient, modulus))
    if not leaves:
        return 0, 1
    while len(leaves) > 1:
        merged = []
        for index in range(0, len(leaves) - 1, 2):
            ln, ld = leaves[index]
            rn, rd = leaves[index + 1]
            merged.append((ln * rd + rn * ld, ld * rd))
        if len(leaves) & 1:
            merged.append(leaves[-1])
        leaves = merged
    numerator, denominator = leaves[0]
    return numerator % denominator, denominator


def main():
    generator = random.Random(20261024)
    cases = []
    for n in range(5):
        for cutoff in (0, 1, 3, 7, 19, 63, 121):
            cases.extend(([family(5, n, cutoff)], [family(239, n, cutoff)]))
    for _ in range(60):
        n = generator.randrange(100)
        cases.append([family(5, n, generator.randrange(601)),
                      family(239, n, generator.randrange(601))])
    projection_checks = 0
    for index, families in enumerate(cases):
        reference = direct_sum(families)
        numerator, denominator = exact_fraction(families, (1, 3, 7, 16)[index % 4])
        assert Fraction(numerator, denominator) == reference, (index, families)
        for guard in (8, 64, 192, 4096):
            assert (numerator << guard) // denominator == (reference.numerator << guard) // reference.denominator
            projection_checks += 1

    # Exercise deep removed p0 powers without constructing a huge whole sum.
    primary_checks = 0
    for p0, mu in ((5, 5**3 * 3**2 * 7), (5, 5**7 * 11),
                   (239, 239**2 * 3), (239, 239 * 5**2 * 7)):
        for n in (0, 1, 17, 43):
            f = family(p0, n, (mu - 1) // 2)
            signed, unit, _ = direct_word(f, mu)
            for prime in (3, 5, 7, 11):
                if prime == p0 or mu % prime:
                    continue
                cofactor, power = mu, 1
                while cofactor % prime == 0:
                    cofactor //= prime
                    power *= prime
                modulus = prime_power(prime, mu)
                original = ((modulus // power) * signed * pow(unit // power, -1, modulus)) % modulus
                weight = (f["sign"] * (-1 if ((mu - 1) // 2) & 1 else 1)
                          * f["c"] * pow(f["b"], n, modulus)
                          * pow(p0, f["s"] - mu, modulus)) % modulus
                cancelled = ((modulus // power) * weight * pow(cofactor, -1, modulus)) % modulus
                assert original == cancelled, (p0, mu, n, prime)
                primary_checks += 1

    source = Path(__file__).resolve()
    report = {"status": "PASS", "whole_fraction_cases": len(cases),
              "large_guard_projection_checks": projection_checks,
              "deep_p0_cancellation_components": primary_checks,
              "maximum_whole_sum_endpoint": max(2 * f["K"] + 1 for row in cases for f in row),
              "maximum_projection_guard": 4096,
              "scope": "Exact word fractions only; not an integrated adaptive extractor or timing claim",
              "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    output = source.with_name("machin_exact_words_validation.json")
    with output.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
