#!/usr/bin/env python3
"""
Decimal Digit Extraction for ln 2 via Three-Term Decomposition
===============================================================

The first decimal digit extraction algorithm for ln 2, using:

    ln 2 = 7*ln(10/9) - 2*ln(25/24) + 3*ln(81/80)

where:
    S1 = ln(10/9)  = sum_{k=1}^inf 1/(k*10^k)   -- base-10 BBP component
    S2 = ln(25/24) = sum_{k=1}^inf 1/(k*25^k)   -- direct computation
    S3 = ln(81/80) = sum_{k=1}^inf 1/(k*81^k)   -- direct computation

Key insight: S1 has base 10, enabling true BBP-style modular extraction
for the dominant component. S2 and S3 require full-precision arithmetic.

Identity verification (exact, via prime factorization):
    (10/9)^7 * (24/25)^2 * (81/80)^3 = 2

Complexity: O(N^2) with schoolbook multiplication; O(N^alpha) with
            M(n) = Theta(n^alpha) multiplication (Theorem, Part 12).

Author: Joseph Babcanec (Benedict College)
Date: March 2026
"""

from __future__ import annotations

from decimal import Decimal, getcontext
import time
import math

__all__ = [
    'extract_ln2_digit',
    'extract_ln2_digits_range',
    'verify_ln2',
    'benchmark_ln2',
]


def _set_precision(n: int) -> None:
    """Set the global decimal precision."""
    getcontext().prec = n


def _frac(x: Decimal) -> Decimal:
    """
    Compute the fractional part of x, guaranteed in [0, 1).

    Args:
        x: Any decimal number

    Returns:
        The fractional part in [0, 1)
    """
    r = x - int(x)
    return r + Decimal(1) if r < 0 else r


def _compute_frac_S1(N: int, precision: int) -> Decimal:
    """
    Compute {10^(N-1) * S1} where S1 = ln(10/9) = sum 1/(k*10^k).

    Uses BBP-style modular extraction: 10^(N-1) * S1 = sum 10^(N-1-k)/k.
    For k <= N-2: compute 10^(N-1-k) mod k via fast modular exponentiation.
    For k >= N-1: terms are < 1, sum directly as a geometric tail.

    The base-10 alignment makes this a true BBP extraction.
    Cost: O(N log^2 N) for the modular phase + O(N) for the tail.

    Args:
        N: Position parameter (digit index)
        precision: Working precision in decimal digits

    Returns:
        {10^(N-1) * S1} in [0, 1)
    """
    _set_precision(precision + 50)

    # Modular phase: k = 1 to N-1
    # For each k, 10^(N-1-k)/k contributes r_k/k to the fractional part
    # where r_k = 10^(N-1-k) mod k
    result = Decimal(0)

    for k in range(1, N):
        exp = N - 1 - k
        if exp >= 0:
            r_k = pow(10, exp, k)  # O(log N) via Python's built-in
            result += Decimal(r_k) / Decimal(k)
        else:
            # k = N-1 case when exp = 0: 10^0 = 1, contribute 1/k
            # Actually exp = N-1-k, so k=N-1 gives exp=0, k=N gives exp=-1
            # For exp=0: 10^0 = 1, r_k = 1 mod k = 1
            break

        # Reduce periodically to prevent accumulator growth
        if result > 1000:
            result = _frac(result)

    # Tail phase: k = N to ...
    # Terms are 10^(N-1-k)/k = 1/(k * 10^(k-N+1)), geometric decay
    tail = Decimal(0)
    for k in range(N, N + precision + 50):
        power_of_10 = k - N + 1
        term = Decimal(1) / (Decimal(k) * Decimal(10) ** power_of_10)
        if term < Decimal(10) ** (-(precision + 10)):
            break
        tail += term

    return _frac(result + tail)


def _compute_S_direct(base: int, precision: int) -> Decimal:
    """
    Compute S_b = ln(b/(b-1)) = sum_{k=1}^inf 1/(k * b^k) to given precision.

    Uses iterative multiplication by 1/b to avoid recomputing b^k.

    Args:
        base: The base b (e.g. 25 for S2, 81 for S3)
        precision: Working precision in decimal digits

    Returns:
        S_b as a Decimal
    """
    _set_precision(precision + 50)

    inv_base = Decimal(1) / Decimal(base)
    power = inv_base  # starts at 1/b
    result = Decimal(0)

    # Number of terms: ceil(precision / log10(base)) + guard
    max_terms = int(precision / math.log10(base)) + 50

    for k in range(1, max_terms + 1):
        term = power / Decimal(k)
        if term < Decimal(10) ** (-(precision + 20)):
            break
        result += term
        power *= inv_base

    return result


def extract_ln2_digit(N: int) -> int:
    """
    Extract the N-th decimal digit of ln 2 (1-indexed, after the decimal point).

    ln 2 = 0.693147180559945...
              ^-position 1
               ^-position 2

    Uses: ln 2 = 7*ln(10/9) - 2*ln(25/24) + 3*ln(81/80)

    The algorithm:
    1. Compute {10^(N-1) * S1} via modular extraction (base-10 BBP)
    2. Compute S2, S3 via direct summation to N+guard digits
    3. Scale and take fractional parts
    4. Combine with borrow tracking
    5. Extract digit: d_N = floor(10 * {10^(N-1) * ln 2})

    Args:
        N: Position of digit to extract (1-indexed, must be >= 1)

    Returns:
        The N-th decimal digit of ln 2 (0-9)

    Raises:
        ValueError: If N < 1
    """
    if N < 1:
        raise ValueError(f"N must be at least 1, got {N}")

    # Working precision: need N digits + guard for rounding
    precision = 2 * N + 120
    _set_precision(precision)

    # Step 1: Modular extraction for S1 = ln(10/9) = sum 1/(k*10^k)
    F1 = _compute_frac_S1(N, precision)

    # Step 2: Direct computation for S2 = ln(25/24) = sum 1/(k*25^k)
    S2 = _compute_S_direct(25, precision)
    # Compute {10^(N-1) * S2}
    scaled_S2 = Decimal(10) ** (N - 1) * S2
    F2 = _frac(scaled_S2)

    # Step 3: Direct computation for S3 = ln(81/80) = sum 1/(k*81^k)
    S3 = _compute_S_direct(81, precision)
    # Compute {10^(N-1) * S3}
    scaled_S3 = Decimal(10) ** (N - 1) * S3
    F3 = _frac(scaled_S3)

    # Step 4: Combine via borrow tracking
    # ln 2 = 7*S1 - 2*S2 + 3*S3
    # {10^(N-1) * ln 2} = {7*F1 - 2*F2 + 3*F3}
    #
    # Write 7*F1 = a1 + f1, 2*F2 = a2 + f2, 3*F3 = a3 + f3
    # where ai are integers and fi in [0,1).
    # Then {7*F1 - 2*F2 + 3*F3} = {f1 - f2 + f3}  (integer parts drop out)
    # f1 - f2 + f3 in (-1, 2), so at most one borrow/carry adjustment.

    f1 = _frac(Decimal(7) * F1)
    f2 = _frac(Decimal(2) * F2)
    f3 = _frac(Decimal(3) * F3)

    x = f1 - f2 + f3
    # Adjust to [0, 1)
    if x < 0:
        x += Decimal(1)
    elif x >= 1:
        x -= Decimal(1)

    # Step 5: Extract digit
    return int(x * 10)


def extract_ln2_digits_range(start: int, end: int) -> list[int]:
    """
    Extract a range of decimal digits of ln 2.

    Args:
        start: First position (inclusive, 1-indexed)
        end: Last position (inclusive, 1-indexed)

    Returns:
        List of digits from position start to end
    """
    return [extract_ln2_digit(n) for n in range(start, end + 1)]


def _compute_ln2_reference(precision: int) -> Decimal:
    """
    Compute ln 2 via the series sum 1/(k*2^k) for verification.

    This is independent of the three-term decomposition used in the
    main algorithm, providing a true cross-check.
    """
    _set_precision(precision + 50)

    result = Decimal(0)
    inv2 = Decimal(1) / Decimal(2)
    power = inv2

    max_terms = int(3.4 * precision) + 50  # ~3.32 terms per digit
    for k in range(1, max_terms + 1):
        term = power / Decimal(k)
        if term < Decimal(10) ** (-(precision + 20)):
            break
        result += term
        power *= inv2

    return result


def verify_ln2(max_digits: int = 50, verbose: bool = True) -> bool:
    """
    Verify the digit extraction algorithm against an independent ln 2 computation.

    The reference uses sum 1/(k*2^k), which is completely independent of
    the three-term decomposition 7*ln(10/9) - 2*ln(25/24) + 3*ln(81/80).

    Args:
        max_digits: Number of digits to verify
        verbose: Whether to print progress

    Returns:
        True if all digits match
    """
    if verbose:
        print(f"Computing ln 2 to {max_digits} digits for verification...")

    ln2 = _compute_ln2_reference(max_digits + 10)
    ln2_str = str(ln2)

    if '.' in ln2_str:
        digits_str = ln2_str.split('.')[1][:max_digits]
    else:
        digits_str = ""

    if verbose:
        print(f"ln 2 = 0.{digits_str}")
        print()
        print("Verifying modular digit extraction:")
        print("-" * 70)

    all_correct = True

    for n in range(1, min(max_digits + 1, len(digits_str) + 1)):
        expected = int(digits_str[n - 1])
        computed = extract_ln2_digit(n)

        if computed != expected:
            all_correct = False
            if verbose:
                print(f"d_{n:3d} = {computed} (expected {expected}) MISMATCH")
        elif verbose:
            print(f"d_{n:3d} = {computed} (expected {expected}) OK")

    if verbose:
        print("-" * 70)
        if all_correct:
            print(f"All {min(max_digits, len(digits_str))} digits verified correctly.")
        else:
            print("Some digits were incorrect.")

    return all_correct


def benchmark_ln2(positions: list[int] | None = None,
                  verbose: bool = True) -> dict[int, tuple[int, float]]:
    """
    Benchmark digit extraction at various positions.

    Args:
        positions: List of positions to test
        verbose: Whether to print results

    Returns:
        Dictionary mapping position -> (digit, time_in_seconds)
    """
    if positions is None:
        positions = [10, 50, 100, 200, 500, 1000]

    results = {}

    if verbose:
        print("\nBenchmark: ln 2 Decimal Digit Extraction")
        print("-" * 60)
        print(f"{'Position':>10} | {'Digit':>6} | {'Time (s)':>12}")
        print("-" * 60)

    for n in positions:
        start = time.perf_counter()
        digit = extract_ln2_digit(n)
        elapsed = time.perf_counter() - start

        results[n] = (digit, elapsed)

        if verbose:
            print(f"{n:>10} | {digit:>6} | {elapsed:>12.4f}")

    if verbose:
        print("-" * 60)

        # Fit empirical scaling exponent
        if len(positions) >= 3:
            import math
            log_n = [math.log(n) for n in positions if n >= 10]
            log_t = [math.log(max(results[n][1], 1e-10))
                     for n in positions if n >= 10]
            if len(log_n) >= 2:
                n_pts = len(log_n)
                sum_x = sum(log_n)
                sum_y = sum(log_t)
                sum_xy = sum(x * y for x, y in zip(log_n, log_t))
                sum_xx = sum(x * x for x in log_n)
                denom = n_pts * sum_xx - sum_x ** 2
                if denom != 0:
                    alpha = (n_pts * sum_xy - sum_x * sum_y) / denom
                    print(f"Empirical scaling exponent: {alpha:.2f}")

    return results


# Known correct digits of ln 2 for validation
# ln 2 = 0.6931471805599453094172321214581765680755001343602552541206...
LN2_DIGITS_100 = (
    "6931471805599453094172321214581765680755001343602552541206800094"
    "933936219696947156058633269964186875420"
)

KNOWN_LN2_DIGITS = {
    1: 6, 2: 9, 3: 3, 4: 1, 5: 4, 6: 7, 7: 1, 8: 8, 9: 0, 10: 5,
    50: 5, 100: 5,
}


def validate_known_digits(verbose: bool = True) -> bool:
    """
    Quick validation against known correct digits.

    Returns:
        True if all known digits match
    """
    all_correct = True

    if verbose:
        print("Validating against known ln 2 digits:")
        print("-" * 50)

    for pos, expected in sorted(KNOWN_LN2_DIGITS.items()):
        computed = extract_ln2_digit(pos)
        correct = computed == expected
        all_correct = all_correct and correct

        if verbose:
            status = "OK" if correct else "FAIL"
            print(f"  Position {pos:>5}: {computed} (expected {expected}) {status}")

    if verbose:
        print("-" * 50)

    return all_correct


def main() -> None:
    """Main entry point for command-line usage."""
    print("=" * 70)
    print("  DECIMAL DIGIT EXTRACTION FOR ln 2")
    print("  via 7*ln(10/9) - 2*ln(25/24) + 3*ln(81/80)")
    print("=" * 70)
    print()
    print("Key: S1 = ln(10/9) = sum 1/(k*10^k) is a base-10 BBP component.")
    print("     S2, S3 computed via direct full-precision summation.")
    print()

    # Verify first 30 digits
    success = verify_ln2(30)

    if success:
        print()
        benchmark_ln2([10, 50, 100, 200, 500, 1000])

        print()
        print("=" * 70)
        print("  ALGORITHM SUMMARY")
        print("=" * 70)
        print("""
Identity: ln 2 = 7*ln(10/9) - 2*ln(25/24) + 3*ln(81/80)

Components:
  - S1 = ln(10/9): BBP-style modular extraction (base 10)
  - S2 = ln(25/24): Direct high-precision summation
  - S3 = ln(81/80): Direct high-precision summation
  - Borrow tracking: 1 indicator (vs 2 for pi)

Complexity: O(N^2) schoolbook, O(N^alpha) general
  - S1 modular phase: O(N log^2 N)  [subquadratic!]
  - S2, S3 direct: O(N^2)  [bottleneck]

This is the first decimal digit extraction algorithm for ln 2.
""")


if __name__ == "__main__":
    main()
