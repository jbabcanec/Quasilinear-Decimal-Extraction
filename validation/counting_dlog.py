#!/usr/bin/env python3
"""
counting_dlog.py -- the GENERALIZED counting face contains DISCRETE LOG.

Generalized orbit-interval counting problem (the natural generalization of the
column popcounts B(d) of the counting face):
    C(D, mu, c, I) = #{ 0 <= d < D : (c * 2^d) mod mu  in  I },
for an interval/set I (mod mu). The counting face's B(d) is the instance
c = 1, I = [mu/2, mu), read at all prefixes D, summed over the modulus family.

REDUCTION (validated here): given a counting oracle for C, the discrete log of
any y in <2> mod mu is recovered with O(log ord) oracle calls: with I = {y}
(singleton), f(D) = C(D, mu, 1, {y}) is a 0/1 step function on one period,
jumping from 0 to 1 at D = dlog_2(y) + 1; binary search on D finds the jump.

CONSEQUENCE: if the GENERAL problem (arbitrary singleton/interval, arbitrary
modulus mu) were solvable in polylog(mu) per query, discrete log in <2> mod mu
would fall in polylog(mu). So exact orbit-interval counting is DLOG-HARD at
cryptographic modulus sizes -- a genuine hardness anchor for the general
counting face, replacing the falsified Girstmair-popcount anchor with a
correct, self-validated one.

HONEST CAVEAT (why this does NOT prove OUR pipeline instance hard): the
extraction pipeline's moduli are mu_t = 8t+1 with t <= T ~ N, i.e.
POLYNOMIALLY SMALL (mu = poly(N)). At poly-size moduli discrete log is
poly(N)-easy (even trivial: O(mu) baby steps), so the dlog wall is vacuous for
our instance. Our hardness is NOT dlog-hardness; it is the MASS of computing
T ~ N such counts jointly to N depths (Theta(N^2/polylog), the total orbit mass
Sigma ord(2,mu_t) -- see fusion_mass.py). The two hardnesses are different
axes: dlog-hardness is per-query at large mu (crypto sizes); our wall is
aggregate at small mu (the T x N grid). This script proves the FORMER exactly
and states the gap to the latter explicitly -- so nobody mistakes the general
hardness for a proof about the pipeline.
"""

from sympy import n_order


def counting_oracle(D, mu, c, target):
    """C(D, mu, c, {target}) = #{0<=d<D : c*2^d mod mu == target}. Reference impl."""
    cnt = 0
    v = c % mu
    for _ in range(D):
        if v == target:
            cnt += 1
        v = (2 * v) % mu
    return cnt


def dlog_via_counting(mu, y):
    """Recover d with 2^d == y (mod mu) using only the counting oracle, by
    binary search on the step function f(D) = C(D, mu, 1, {y})."""
    p = n_order(2, mu)                       # period (search range [0, p))
    # f(D) = #{d < D : 2^d == y} is 0 for D <= dlog, 1 for D > dlog (within one
    # period). Find least D with f(D) = 1; then dlog = D - 1.
    lo, hi = 0, p
    calls = 0
    while lo < hi:
        mid = (lo + hi) // 2
        calls += 1
        if counting_oracle(mid, mu, 1, y) >= 1:
            hi = mid
        else:
            lo = mid + 1
    return lo - 1, calls


def main():
    print("Reduction: discrete log in <2> mod mu from a counting oracle only.")
    print("  mu     y      recovered_d   check 2^d==y   oracle_calls   period")
    ok = True
    trials = [(101, 5), (257, 10), (1009, 3), (2003, 7), (4099, 42), (7919, 100)]
    for mu, base_exp in trials:
        p = n_order(2, mu)
        d_true = base_exp % p
        y = pow(2, d_true, mu)
        d_rec, calls = dlog_via_counting(mu, y)
        good = (pow(2, d_rec, mu) == y)
        ok &= good
        print(f"{mu:6d}  {y:6d}   {d_rec:8d}      {'OK' if good else 'X':5s}"
              f"       {calls:6d}       {p}")
    print(f"\nReduction correct on all trials: {ok}")
    print(f"Oracle calls per dlog = O(log period) (binary search) -- confirmed above.")
    print("\nVERDICT: general orbit-interval counting is DLOG-HARD (this reduction).")
    print("CAVEAT: pipeline moduli mu_t = 8t+1 are poly-small => dlog vacuous there;")
    print("        our wall is the AGGREGATE T x N orbit mass, not per-query dlog.")


if __name__ == "__main__":
    main()
