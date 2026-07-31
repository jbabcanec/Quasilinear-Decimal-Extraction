#!/usr/bin/env python3
"""sqrt2_extract.py -- the start-over calibration (2026-07-28).

*** CORRECTED 2026-07-31 (see validation/depth_jump.py). The conclusion below
*** that extraction costs "a log MORE than full computation" is WRONG, and the
*** phrase "the extraction tax" should not be used for sqrt2. It is an artifact
*** of choosing a MODEXP (10^N mod q), where the reduction pins every one of
*** the log N squarings at full size. Computing the exact power 10^N instead
*** TELESCOPES (operand sizes halve going backwards), so the whole depth-N jump
*** is O(M(N)) = N log N: one multiply and one division after it. Measured
*** operand-mass ratio modexp/exact-power = 2.0, 3.0, 4.5, 6.0 at
*** N = 1e4, 1e6, 1e9, 1e12 -- a clean log. So sqrt2 extraction is at PARITY
*** with its own computation (Newton, N log N), NOT one log above it. pi is
*** the only constant measured so far with a genuine extraction gap.
*** The code below is correct; only its complexity conclusion was wrong.

Question (Joseph): what if we abandon BBP entirely? Control experiment: digit
extraction for sqrt(2), the maximally structured non-BBP constant (periodic
continued fraction, single-convergent approximation, zero series involved).

Method: Pell convergents p/q of sqrt2 (p^2 - 2q^2 = +-1) give
|sqrt2 - p/q| ~ 1/(2.83 q^2), so with q ~ 10^{(N+G)/2} one rational suffices:
  frac(10^N sqrt2) = (10^N p mod q)/q + tiny known correction.
Cost: one modexp 10^N mod q with q an N/2-DIGIT modulus = Theta(log N)
squarings at M(N) each = N log^2 N bit ops, versus N log N to compute ALL N
digits by Newton. CONCLUSION the experiment documents: even for an algebraic
number with no BBP and no series, single-digit extraction by the best known
route costs a log MORE than full computation -- the "extraction tax" is not
a BBP artifact. And sqrt2 shows exactly what structure buys: a quadratic
irrational needs ONE convergent (one modexp, no sum); e (patterned CF,
hypergeometric convergents) needs the merge tree; pi likewise. The tree
appears precisely when infinitely many independent terms must be combined;
BBP is not the cause of the wall, it is the reason pi's binary digits are
EASY compared to generic transcendentals.
"""

from mpmath import mp, mpf, sqrt


def pell_convergent(min_q_digits):
    """Smallest Pell convergent with q >= 10^min_q_digits (p/q -> sqrt2)."""
    p, q = 1, 1                      # (1+sqrt2)^1: p+q*sqrt2
    target = 10 ** min_q_digits
    while q < target:
        p, q = p + 2 * q, p + q      # multiply by (1+sqrt2)
    return p, q


def sqrt2_digits_extract(N, window=10, guard=25):
    """Digits of sqrt2 at decimal positions N+1..N+window via one convergent."""
    p, q = pell_convergent((N + guard) // 2 + guard)
    r = (pow(10, N, q) * (p % q)) % q
    # window of r/q needs only the leading digits of r and q, plus the known
    # correction 10^N (sqrt2 - p/q) = sign/(q(p+q sqrt2)) * 10^N, tiny by
    # construction; fold both at small fixed precision.
    mp.dps = window + guard
    frac = mpf(r) / mpf(q)           # mpf rounds; leading digits suffice
    sign = 1 if p * p - 2 * q * q < 0 else -1
    corr = sign * mpf(10) ** N / (mpf(q) * (mpf(p) + mpf(q) * sqrt(2)))
    val = (frac + corr) % 1
    return ('%0.*f' % (window + 5, val))[2:2 + window]   # truncate, not round


def main():
    for N in (100, 1000, 5000, 20000):
        got = sqrt2_digits_extract(N)
        mp.dps = N + 40
        ref = mp.nstr(sqrt(2), N + 20, strip_zeros=False)
        ref_digits = ref.split('.')[1][N:N + 10]
        ok = got == ref_digits
        print("N=%-6d extract=%s reference=%s %s"
              % (N, got, ref_digits, "OK" if ok else "MISMATCH"))
        if not ok:
            return 1
    print("VERDICT: sqrt2 digit extraction works from ONE convergent + one")
    print("modexp (no series, no fraction sum, no BBP): ~N log^2 N bit ops,")
    print("one log ABOVE full computation (Newton, N log N). The extraction")
    print("tax exists with zero BBP involvement; the merge tree appears only")
    print("when infinitely many independent terms must combine (e, pi).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
