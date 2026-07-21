#!/usr/bin/env python3
"""
Validates Proposition (the additive core is an isolated harmonic number) of
decimal_extraction_paper.tex.

Stripping the dyadic weights 2^(Jmax-Jt) from V leaves the numerator Pi'(0) of
Pi(t)=prod_{k<=K}(8k+1+t), and Pi'(0)/Pi(0) = sum_{k<=K} 1/(8k+1) mod 2^m -- a
2-adic digamma value / (via n! H_n = |s(n+1,2)|) an unsigned Stirling number of
the first kind.  Also exhibits WHY the multiplicative crack does not cross to
the additive side: Pi(0)=prod(8k+1) is SMOOTH (Borwein's factorization applies)
while its t-derivative Pi'(0) is a GENERIC integer with a large prime factor --
differentiation destroys the factorization.
"""


def smooth_part(n, bound=100000):
    """Strip prime factors < bound; leftover > 1 means a large prime remains."""
    x = n
    for p in range(2, bound):
        while x % p == 0:
            x //= p
    return n.bit_length(), x


def _desc(left):
    return ("smooth (<1e5)" if left == 1
            else "has a large prime factor (%d bits)" % left.bit_length())


def check(K=48, m=256):
    mod = 1 << m
    mus = [8 * k + 1 for k in range(K + 1)]
    # Pi(t)=prod(8k+1+t) truncated at t^2:  c0 = Pi(0) = P,  c1 = Pi'(0) = numerator
    c0, c1 = 1, 0
    for u in mus:
        c0, c1 = c0 * u, c1 * u + c0
    P, Nu = c0, c1
    assert Nu == sum(P // u for u in mus), "Pi'(0) != leave-one-out sum"
    # harmonic sum mod 2^m, two ways
    H1 = sum(pow(u, -1, mod) for u in mus) % mod
    H2 = (Nu % mod) * pow(P % mod, -1, mod) % mod
    assert H1 == H2, "harmonic numerator identity FAILED"
    bP, leftP = smooth_part(P)
    bN, leftN = smooth_part(Nu)
    print("K=%4d: P =Pi(0)  %5d bits, %s (largest prime <= 8K+1=%d)  COMPACT"
          % (K, bP, _desc(leftP), 8 * K + 1))
    print("        Nu=Pi'(0) %5d bits, %s  GENERIC" % (bN, _desc(leftN)))
    print("        Nu/P = sum_k 1/(8k+1) mod 2^%d : VERIFIED" % m)
    return H1 == H2


if __name__ == '__main__':
    ok = check(48, 256) and check(200, 512)
    print("\nATOM CHECK", "PASSED" if ok else "FAILED")
    print("Pi(0) factors (Borwein crack); Pi'(0) does not -- d/dt mixes the primes,")
    print("which is exactly why the multiplicative faces fall and the additive core does not.")
