#!/usr/bin/env python3
"""product_face.py -- can pi be extracted from a PRODUCT instead of a sum?
(2026-07-31)

Joseph: "what about a product series."  This is the right question to ask,
because products are the ONE thing this campaign actually cracked: Borwein
factorization computes prod(8k+1) mod 2^m in O(N log^2 N) = parity, precisely
because a product's information compresses into an exponent vector (Theta(N)
bits) while its term list carries Theta(N log N).  Sums have no factorization.
So if pi came out of a product with word-size factors, the crack would apply
and the third log would be gone.

Four measurements, which together close the door:

(A) THE BRIDGE THAT LOOKS REAL.  pi*sqrt2 = Gamma(1/4)Gamma(3/4), and Euler's
    limit turns that into a ratio of FOUR AP-products, every one of them
    Borwein-crackable:
        pi*sqrt2 = lim_n  (n!)^2 * n * 16^(n+1) / (prod(4k+1) * prod(4k+3)).
    Verified here.  Then measured: it converges like 1/n, so 2^N factors.

(B) THE DICHOTOMY.  A product with word-size factors is a fixed rational
    function of k of degree j.  Matching j leading coefficients makes it
    converge like n^-j, so T = 2^(N/j) factors.  Word-size factors force
    j = O(1) and T exponential.  Measured convergence exponents for the
    classical products.

(C) THE OTHER BRANCH.  Geometric convergence forces factors 1 + c_k 2^(-4k),
    whose numerators 2^(4k)+c_k carry Theta(k) bits, not word-size.  Taking
    logs lands back on a BBP-type sum -- with MORE terms, by the divisor
    count.  Measured against the real band.

(D) THE ACCELERATION.  Stirling/Euler-Maclaurin fixes (A)'s 1/n convergence,
    at the cost of a tapered Bernoulli array.  Measured: Theta(N^2).
"""

import math

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int

from mpmath import mp, mpf, sqrt as mp_sqrt, pi as mp_pi


# ---------------------------------------------------------------------------
# (A) pi*sqrt(2) as a ratio of four Borwein-crackable AP-products
# ---------------------------------------------------------------------------

def gamma_quarter_product(n):
    """(n!)^2 * n * 16^(n+1) / (prod_{k=0..n}(4k+1) * prod_{k=0..n}(4k+3))."""
    fact = mpz(1)
    for k in range(2, n + 1):
        fact *= k
    p1 = mpz(1)
    p3 = mpz(1)
    for k in range(0, n + 1):
        p1 *= 4 * k + 1
        p3 *= 4 * k + 3
    num = fact * fact * n * mpz(16) ** (n + 1)
    return mpf(num) / mpf(p1 * p3)


def partA():
    mp.dps = 60
    target = mp_pi * mp_sqrt(2)
    print("=" * 74)
    print("(A) pi*sqrt2 = Gamma(1/4)Gamma(3/4) as a ratio of FOUR AP-products")
    print("    (n!)^2 n 16^(n+1) / (prod(4k+1) prod(4k+3))  ->  pi*sqrt2")
    print("    every factor family is Borwein-crackable at N log^2 N")
    print("=" * 74)
    print("  %-8s %-24s %-14s %-12s" % ("n", "value", "rel error", "error * n"))
    prev = None
    for n in (10, 100, 1000, 4000, 16000):
        v = gamma_quarter_product(n)
        err = abs(v - target) / target
        print("  %-8d %-24s %-14.4e %-12.6f"
              % (n, mp.nstr(v, 18), float(err), float(err * n)))
        prev = err
    print("  target  %s" % mp.nstr(target, 18))
    print("  READING: error*n is CONSTANT => convergence is Theta(1/n).")
    print("  To reach 2^-N you need n ~ 2^N factors. The identity is exact and")
    print("  every product in it is crackable; the CONVERGENCE is what kills it.")


# ---------------------------------------------------------------------------
# (B) convergence exponent of products with word-size rational factors
# ---------------------------------------------------------------------------

def measure_exponent(factor, limit, n0=400, reps=5):
    """Measured j in error ~ n^-j by doubling n against an EXACT limit."""
    mp.dps = 50
    errs = []
    n = n0
    p = mpf(1)
    k = 1
    for _ in range(reps):
        while k <= n:
            p *= factor(k)
            k += 1
        errs.append(abs(p - limit) / abs(limit))
        n *= 2
    js = [math.log(float(errs[i] / errs[i + 1])) / math.log(2)
          for i in range(len(errs) - 1) if errs[i + 1] != 0]
    return js[-1] if js else float('nan')


def shift_family(a, d):
    """prod_k prod_i (k+a_i)/(k+d_i), with its EXACT limit
       prod_i Gamma(1+d_i)/Gamma(1+a_i) (valid when sum a = sum d)."""
    from mpmath import gamma as mp_gamma

    def f(k):
        num = mpf(1)
        den = mpf(1)
        for x in a:
            num *= (k + x)
        for x in d:
            den *= (k + x)
        return num / den
    lim = mpf(1)
    for x in d:
        lim *= mp_gamma(1 + x)
    for x in a:
        lim /= mp_gamma(1 + x)
    return f, lim


def partB():
    mp.dps = 50
    print()
    print("=" * 74)
    print("(B) how fast can a product with WORD-SIZE factors converge?")
    print("    (exact Gamma-closed-form limits, so the exponent is unbiased)")
    print("=" * 74)
    f2, l2 = shift_family([0, 3], [1, 2])            # sums 3 = 3
    f3, l3 = shift_family([1, 5], [2, 4])            # sums 6 = 6
    f4, l4 = shift_family([1, 5, 6], [2, 3, 7])      # sums 12, squares 62 both
    tests = [
        ("Wallis 4k^2/(4k^2-1) -> pi/2  [1 moment]",
         lambda k: mpf(4 * k * k) / (4 * k * k - 1), mp_pi / 2),
        ("(k)(k+3)/((k+1)(k+2))         [1 moment]", f2, l2),
        ("(k+1)(k+5)/((k+2)(k+4))       [1 moment]", f3, l3),
        ("(k+1)(k+5)(k+6)/((k+2)(k+3)(k+7)) [2 mom]", f4, l4),
        ("Gamma-limit ratio for pi*sqrt2 (part A)", None, None),
    ]
    print("  %-46s %s" % ("product", "measured j in error ~ n^-j"))
    for name, f, lim in tests:
        if f is None:
            e = []
            for n in (500, 1000, 2000, 4000):
                e.append(abs(gamma_quarter_product(n) - mp_pi * mp_sqrt(2)))
            j = math.log(float(e[-2] / e[-1])) / math.log(2)
            print("  %-46s %.3f" % (name, j))
            continue
        print("  %-46s %.3f" % (name, measure_exponent(f, lim)))
    print()
    print("  A degree-j factor rule matches j moments and converges like n^-j,")
    print("  so T = 2^(N/j) factors are needed. Word-size factors => j = O(1)")
    print("  => T exponential. Getting j = Theta(N/log N) needs factors that")
    print("  are degree-N/log N polynomials, i.e. Theta(N)-bit coefficients:")
    print("  the mass comes straight back. At the extreme T=1, j=N it IS just")
    print("  'pi ~ p/q with a CF convergent', which requires knowing pi.")


# ---------------------------------------------------------------------------
# (C) the geometrically-convergent branch: factors 1 + c_k 2^(-4k)
# ---------------------------------------------------------------------------

def partC(N=10 ** 6):
    Jmax = int(2.3219 * N)
    T_band = Jmax                       # measured: hex BBP band, 1 term/bit
    # log(prod (1 + c_k 2^-4k)) = sum_{k,j} (-1)^{j+1} c_k^j 2^{-4kj} / j
    # term count = #{(k,j) : 4kj <= Jmax} = divisor summatory D(Jmax/4)
    M = Jmax // 4
    terms = sum(M // j for j in range(1, M + 1)) if M < 10 ** 6 else \
        int(M * (math.log(M) + 2 * 0.5772156649 - 1))
    print()
    print("=" * 74)
    print("(C) the geometrically-convergent product, N = %.0e" % N)
    print("=" * 74)
    print("  a factor 1 + c 2^(-4k) has numerator 2^(4k)+c : %d bits at k, NOT"
          % (4 * (Jmax // 4)))
    print("  word-size. So the Borwein crack does not apply to it at all.")
    print("  Taking logs lands on a BBP-type sum. Term counts:")
    print("    hex BBP band                       %12d terms" % T_band)
    print("    log-expansion of the same product  %12d terms  (%.2fx MORE)"
          % (terms, terms / T_band))
    print("  The blowup is the divisor count sum_j Jmax/(4j) ~ (Jmax/4) ln Jmax.")
    print("  So the product route is a factor ~%.1f WORSE than the sum it came"
          % (terms / T_band))
    print("  from, and identical in shape. Products and sums are the same")
    print("  object under log/exp once the factors approach 1.")


# ---------------------------------------------------------------------------
# (D) accelerating (A) with Stirling: the Bernoulli taper
# ---------------------------------------------------------------------------

def partD():
    print()
    print("=" * 74)
    print("(D) fixing (A)'s 1/n convergence with Euler-Maclaurin / Stirling")
    print("=" * 74)
    print("  %-10s %-10s %-12s %-16s %-10s %-10s" %
          ("N", "n needed", "B-terms m", "tapered mass", "mass/N^2",
           "x log2 n"))
    for N in (10 ** 4, 10 ** 5, 10 ** 6, 10 ** 9):
        # optimal Stirling truncation: error ~ e^{-2 pi n}, need 2^-N
        n = N * math.log(2) / (2 * math.pi)
        m = math.pi * n                      # optimal term count
        ln = math.log2(max(n, 2))
        # B_{2j} needed to N - (2j-1)*log2(n) bits; taper dies at j = N/(2 ln)
        jmax = min(m, N / (2 * ln))
        mass = sum(N - 2 * j * ln for j in range(1, int(jmax) + 1)) \
            if jmax < 10 ** 6 else jmax * N - ln * jmax * jmax
        print("  %-10.0e %-10.3e %-12.3e %-16.4e %-10.4f %-10.4f"
              % (N, n, m, mass, mass / N ** 2, mass * ln / N ** 2))
    print("  Bernoulli mass is Theta(N^2 / log N) (last column flat) -- the")
    print("  campaign's already-recorded")
    print("  'materialized tapered Bernoulli array: output size Omega(N^2)'.")
    print("  This is exactly where the product route dies.")


def main():
    partA()
    partB()
    partC()
    partD()
    print()
    print("=" * 74)
    print("VERDICT: products are crackable but never converge; convergent")
    print("products are sums in disguise, with more terms. The one bridge")
    print("(Gamma(1/4)Gamma(3/4) = four crackable AP-products) is exact and")
    print("dies on convergence, and its acceleration dies on Bernoulli mass.")
    print("=" * 74)


if __name__ == "__main__":
    main()
