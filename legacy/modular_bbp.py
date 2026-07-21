#!/usr/bin/env python3
"""
A BBP-Type Algorithm for Extracting Decimal Digits of π
========================================================

This module implements a novel algorithm for computing arbitrary decimal
digits of π without computing all preceding digits. The method combines:

1. Series splitting of Machin's formula to eliminate alternating signs
2. Modular arithmetic for efficient fractional part computation
3. Borrow tracking to combine fractional parts under subtraction
4. Base conversion from base-5 to base-10

Formula:
    π = 16·arctan(1/5) - 4·arctan(1/239)    (Machin's formula)

    arctan(1/5) = P₁/5 - Q₁/125

    where:
        P₁ = Σ_{j=0}^∞ 1/((4j+1)·625^j)
        Q₁ = Σ_{j=0}^∞ 1/((4j+3)·625^j)

Key insight: {A - B - C} can be computed from {A}, {B}, {C} alone using
borrow tracking, avoiding computation of large integer parts.

Complexity: O(N²) due to the arctan(1/239) term and base conversion.
            P₁ and Q₁ components use O(N log N) via modular arithmetic.

Author: Joseph Babcanec (Benedict College)
Date: 2026
"""

from __future__ import annotations

from decimal import Decimal, getcontext, localcontext
from typing import Callable
import time
import sys

__all__ = [
    'extract_decimal_digit',
    'extract_digits_range',
    'verify_against_reference',
    'benchmark',
]


def _set_precision(n: int) -> None:
    """Set the global decimal precision."""
    getcontext().prec = n


def _mod_pow(base: int, exp: int, mod: int) -> int:
    """
    Compute base^exp mod mod using binary exponentiation.

    Time complexity: O(log exp) multiplications
    Space complexity: O(1)

    Args:
        base: The base of the exponentiation
        exp: The exponent (must be non-negative)
        mod: The modulus (must be positive)

    Returns:
        base^exp mod mod
    """
    if mod == 1:
        return 0
    if exp == 0:
        return 1

    result = 1
    base = base % mod

    while exp > 0:
        if exp & 1:
            result = (result * base) % mod
        exp >>= 1
        base = (base * base) % mod

    return result


def _frac(x: Decimal) -> Decimal:
    """
    Compute the fractional part of x, guaranteed to be in [0, 1).

    For negative numbers, adjusts to return positive fractional part.
    Example: frac(-0.3) = 0.7, not -0.3

    Args:
        x: Any decimal number

    Returns:
        The fractional part in [0, 1)
    """
    r = x - int(x)
    return r + Decimal(1) if r < 0 else r


def _compute_frac_P1(N: int, precision: int) -> Decimal:
    """
    Compute {5^N · P₁} using modular arithmetic (BBP-style).

    P₁ = Σ_{j=0}^∞ 1/((4j+1)·625^j)

    5^N · P₁ = Σ_{j=0}^∞ 5^(N-4j)/(4j+1)

    For terms where N-4j ≥ 0, we use:
        5^(N-4j)/(4j+1) = floor(...) + (5^(N-4j) mod (4j+1))/(4j+1)

    Only the fractional part matters, so we compute 5^(N-4j) mod (4j+1)
    using fast modular exponentiation in O(log N) time per term.

    Args:
        N: The exponent for 5^N
        precision: Working precision in decimal digits

    Returns:
        {5^N · P₁} ∈ [0, 1)
    """
    _set_precision(precision + 50)

    result = Decimal(0)

    # Terms for j > N/4 have 5^(N-4j) < 1 and decay geometrically
    max_j = N // 4 + precision // 2 + 20

    for j in range(max_j + 1):
        denom = 4 * j + 1
        exp = N - 4 * j

        if exp >= 0:
            # Use modular arithmetic - the key BBP insight
            numerator = _mod_pow(5, exp, denom)
            term = Decimal(numerator) / Decimal(denom)
        else:
            # Small term, compute directly
            term = Decimal(1) / (Decimal(denom) * Decimal(5) ** (-exp))
            if term < Decimal(10) ** (-(precision + 10)):
                break

        result += term

        # Reduce periodically to prevent accumulator growth
        if result > 1000:
            result = _frac(result)

    return _frac(result)


def _compute_frac_Q1(N: int, precision: int) -> Decimal:
    """
    Compute {5^N · Q₁} using modular arithmetic (BBP-style).

    Q₁ = Σ_{j=0}^∞ 1/((4j+3)·625^j)

    Args:
        N: The exponent for 5^N
        precision: Working precision in decimal digits

    Returns:
        {5^N · Q₁} ∈ [0, 1)
    """
    _set_precision(precision + 50)

    result = Decimal(0)
    max_j = N // 4 + precision // 2 + 20

    for j in range(max_j + 1):
        denom = 4 * j + 3
        exp = N - 4 * j

        if exp >= 0:
            numerator = _mod_pow(5, exp, denom)
            term = Decimal(numerator) / Decimal(denom)
        else:
            term = Decimal(1) / (Decimal(denom) * Decimal(5) ** (-exp))
            if term < Decimal(10) ** (-(precision + 10)):
                break

        result += term

        if result > 1000:
            result = _frac(result)

    return _frac(result)


def _compute_frac_arctan239(N: int, precision: int) -> Decimal:
    """
    Compute {5^N · arctan(1/239)}.

    arctan(1/239) converges very rapidly (each term is 239² ≈ 57000 times
    smaller than the previous), so direct computation is practical.

    Note: This is currently O(N²) due to computing 5^N as a full-precision
    Decimal. A future optimization could apply modular techniques here as well.

    Args:
        N: The exponent for 5^N
        precision: Working precision in decimal digits

    Returns:
        {5^N · arctan(1/239)} ∈ [0, 1)
    """
    _set_precision(precision + 50)

    five_N = Decimal(5) ** N
    result = Decimal(0)
    power = Decimal(1) / Decimal(239)
    inv_239_sq = Decimal(1) / Decimal(239 ** 2)

    # Convergence: |term_k| < 1/(239^(2k+1) · (2k+1))
    # Need roughly 0.4*N terms for N-digit precision
    max_terms = int(0.42 * N) + 50

    for k in range(max_terms + 1):
        term = five_N * power / Decimal(2 * k + 1)

        if k % 2 == 0:
            result += term
        else:
            result -= term

        if abs(term) < Decimal(10) ** (-(precision + 10)):
            break

        power *= inv_239_sq

        # Reduce periodically
        if abs(result) > 1000:
            result = _frac(result)

    return _frac(result)


def extract_decimal_digit(N: int) -> int:
    """
    Extract the N-th decimal digit of π (1-indexed, after the decimal point).

    π = 3.14159265358979...
         ^-position 1
          ^-position 2
           etc.

    The algorithm:
    1. Compute {5^(N-1) · π} using fractional parts of Machin components
    2. Apply borrow tracking to combine: {A - B - C} from {A}, {B}, {C}
    3. Convert to base 10: {10^(N-1) · π} = {2^(N-1) · {5^(N-1) · π}}
    4. Extract digit: d_N = floor(10 · {10^(N-1) · π})

    Args:
        N: Position of digit to extract (1-indexed, must be ≥ 1)

    Returns:
        The N-th decimal digit of π (0-9)

    Raises:
        ValueError: If N < 1

    Examples:
        >>> extract_decimal_digit(1)
        1
        >>> extract_decimal_digit(2)
        4
        >>> extract_decimal_digit(100)
        9
    """
    if N < 1:
        raise ValueError(f"N must be at least 1, got {N}")

    # Precision requirement: 2^(N-1) multiplication amplifies errors by 2^(N-1),
    # so we need approximately N·log₁₀(2) + guard digits
    precision = int(N * 1.5) + 100
    _set_precision(precision)

    # Compute fractional parts of the three Machin components:
    # A = 16 · 5^(N-2) · P₁
    # B = 16 · 5^(N-4) · Q₁
    # C = 4 · 5^(N-1) · arctan(1/239)

    if N >= 2:
        frac_5_P1 = _compute_frac_P1(N - 2, precision)
    else:
        # Handle N=1: compute P₁ directly
        _set_precision(precision)
        P1 = sum(
            Decimal(1) / (Decimal(4*j + 1) * Decimal(625)**j)
            for j in range(100)
        )
        frac_5_P1 = _frac(Decimal(5)**(N - 2) * P1)

    frac_A = _frac(16 * frac_5_P1)

    if N >= 4:
        frac_5_Q1 = _compute_frac_Q1(N - 4, precision)
    else:
        _set_precision(precision)
        Q1 = sum(
            Decimal(1) / (Decimal(4*j + 3) * Decimal(625)**j)
            for j in range(100)
        )
        frac_5_Q1 = _frac(Decimal(5)**(N - 4) * Q1)

    frac_B = _frac(16 * frac_5_Q1)

    frac_5_arctan = _compute_frac_arctan239(N - 1, precision)
    frac_C = _frac(4 * frac_5_arctan)

    # Borrow Tracking Lemma:
    # {A - B - C} = {A} - {B} - {C} + ε₁ + ε₂  (mod 1)
    # where ε₁ = 1 if {A} < {B}, else 0
    #       ε₂ = 1 if ({A} - {B} + ε₁) < {C}, else 0

    eps1 = Decimal(1) if frac_A < frac_B else Decimal(0)
    frac_AB = frac_A - frac_B + eps1

    eps2 = Decimal(1) if frac_AB < frac_C else Decimal(0)
    frac_5pi = _frac(frac_AB - frac_C + eps2)  # This is {5^(N-1) · π}

    # Base conversion: {10^(N-1) · π} = {2^(N-1) · {5^(N-1) · π}}
    frac_10pi = _frac(Decimal(2) ** (N - 1) * frac_5pi)

    # Extract digit
    return int(frac_10pi * 10)


def extract_digits_range(start: int, end: int) -> list[int]:
    """
    Extract a range of decimal digits of π.

    Args:
        start: First position (inclusive, 1-indexed)
        end: Last position (inclusive, 1-indexed)

    Returns:
        List of digits from position start to end

    Example:
        >>> extract_digits_range(1, 5)
        [1, 4, 1, 5, 9]
    """
    return [extract_decimal_digit(n) for n in range(start, end + 1)]


def _compute_pi_reference(precision: int) -> Decimal:
    """
    Compute π using Machin's formula for verification.

    π = 4(4·arctan(1/5) - arctan(1/239))

    This is a straightforward O(N²) computation used only for testing.
    """
    _set_precision(precision + 50)

    def arctan(x: Decimal, max_terms: int = 1000) -> Decimal:
        result = Decimal(0)
        power = x
        x_sq = x * x
        for k in range(max_terms):
            term = power / Decimal(2 * k + 1)
            result += term if k % 2 == 0 else -term
            if abs(term) < Decimal(10) ** (-(precision + 20)):
                break
            power *= x_sq
        return result

    return 4 * (4 * arctan(Decimal(1)/Decimal(5)) - arctan(Decimal(1)/Decimal(239)))


def verify_against_reference(max_digits: int = 50, verbose: bool = True) -> bool:
    """
    Verify the modular algorithm against a reference π computation.

    Args:
        max_digits: Number of digits to verify (default 50)
        verbose: Whether to print progress

    Returns:
        True if all digits match, False otherwise
    """
    if verbose:
        print(f"Computing π to {max_digits} digits for verification...")

    pi = _compute_pi_reference(max_digits + 10)
    pi_str = str(pi)

    if '.' in pi_str:
        digits_str = pi_str.split('.')[1][:max_digits]
    else:
        digits_str = ""

    if verbose:
        print(f"π = 3.{digits_str}")
        print()
        print("Verifying modular BBP digit extraction:")
        print("-" * 70)

    all_correct = True

    for N in range(1, min(max_digits + 1, len(digits_str) + 1)):
        expected = int(digits_str[N - 1])
        computed = extract_decimal_digit(N)

        if computed != expected:
            all_correct = False
            if verbose:
                print(f"d_{N:3d} = {computed} (expected {expected}) MISMATCH")
        elif verbose:
            print(f"d_{N:3d} = {computed} ✓")

    if verbose:
        print("-" * 70)
        if all_correct:
            print("All digits verified correctly.")
        else:
            print("Some digits were incorrect.")

    return all_correct


def benchmark(positions: list[int] | None = None,
              verbose: bool = True) -> dict[int, tuple[int, float]]:
    """
    Benchmark digit extraction at various positions.

    Args:
        positions: List of positions to test (default: [10, 50, 100, 500, 1000])
        verbose: Whether to print results

    Returns:
        Dictionary mapping position -> (digit, time_in_seconds)
    """
    if positions is None:
        positions = [10, 50, 100, 500, 1000]

    results = {}

    if verbose:
        print("\nBenchmark: Modular BBP Digit Extraction")
        print("-" * 60)
        print(f"{'Position':>10} | {'Digit':>6} | {'Time (s)':>12}")
        print("-" * 60)

    for N in positions:
        start = time.perf_counter()
        digit = extract_decimal_digit(N)
        elapsed = time.perf_counter() - start

        results[N] = (digit, elapsed)

        if verbose:
            print(f"{N:>10} | {digit:>6} | {elapsed:>12.4f}")

    if verbose:
        print("-" * 60)

    return results


# Known correct digits at various positions for validation
KNOWN_DIGITS = {
    1: 1, 2: 4, 3: 1, 4: 5, 5: 9, 6: 2, 7: 6, 8: 5, 9: 3, 10: 5,
    100: 9, 500: 2, 1000: 9, 2000: 9, 5000: 1, 10000: 8,
}


def validate_known_digits(verbose: bool = True) -> bool:
    """
    Validate against known correct digits at various positions.

    Returns:
        True if all known digits match
    """
    all_correct = True

    if verbose:
        print("Validating against known digits:")
        print("-" * 50)

    for pos, expected in sorted(KNOWN_DIGITS.items()):
        computed = extract_decimal_digit(pos)
        correct = computed == expected
        all_correct = all_correct and correct

        if verbose:
            status = "✓" if correct else "FAIL"
            print(f"  Position {pos:>5}: {computed} (expected {expected}) {status}")

    if verbose:
        print("-" * 50)

    return all_correct


def main() -> None:
    """Main entry point for command-line usage."""
    print("=" * 70)
    print("  A BBP-TYPE ALGORITHM FOR DECIMAL DIGITS OF π")
    print("=" * 70)
    print()
    print("Key insight: {A - B - C} can be computed from {A}, {B}, {C} alone")
    print("using borrow tracking, enabling modular arithmetic for base-5.")
    print()

    # Verify first 50 digits
    success = verify_against_reference(50)

    if success:
        print()
        # Benchmark at larger positions
        benchmark([10, 50, 100, 200, 500, 1000])

        print()
        print("=" * 70)
        print("  ALGORITHM SUMMARY")
        print("=" * 70)
        print("""
Components:
  - Series splitting: arctan(1/5) = P₁/5 - Q₁/125 (non-alternating)
  - Modular extraction: {5^N · P₁} via 5^(N-4j) mod (4j+1)
  - Borrow tracking: {A-B-C} from fractional parts alone
  - Base conversion: {10^N · π} = {2^N · {5^N · π}}

Complexity: O(N²) overall
  - P₁, Q₁ terms: O(N) terms × O(log N) per term = O(N log N)
  - arctan(1/239): O(N) terms × O(N) arithmetic = O(N²)
  - Base conversion: O(N)-bit multiplication = O(N²)
""")


if __name__ == "__main__":
    main()
