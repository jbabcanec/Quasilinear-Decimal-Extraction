"""Validate the Borwein-factorization crack:

  Pi_{k=1}^{K} (8k+1) mod 2^m   via   prime factorization + exponent-bit assembly

vs the direct term product. If correct and faster-scaling, the E-face of the
barrier drops to O(N log^2 N): E = Log2(product)/Log2(9).

Method: e_p = sum_j #{1<=k<=K : p^j | 8k+1} by Legendre-type counting on the AP;
then product = prod_j ( prod_{p: bit j of e_p set} p )^{2^j} mod 2^m,
with each inner product built by a truncated product tree over DISTINCT primes
(total prime mass Theta(K) bits, vs Theta(K log K) for the term list).
"""
import time

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b'\x00\x00'
    i = 2
    while i * i <= n:
        if sieve[i]:
            sieve[i * i::i] = b'\x00' * len(sieve[i * i::i])
        i += 1
    return [i for i in range(3, n + 1) if sieve[i]]  # odd primes only (2 never divides 8k+1)


def count_ap(K, q):
    """#{1<=k<=K : q | 8k+1} for odd q."""
    k0 = (-pow(8, -1, q)) % q
    if k0 == 0:
        k0 = q
    return 0 if k0 > K else (K - k0) // q + 1


def product_tree_mod(vals, MOD):
    level = [mpz(v) for v in vals]
    while len(level) > 1:
        nxt = [(level[i] * level[i + 1]) % MOD for i in range(0, len(level) - 1, 2)]
        if len(level) & 1:
            nxt.append(level[-1])
        level = nxt
    return level[0] if level else mpz(1)


def borwein_product(K, m):
    MOD = mpz(1) << m
    ps = primes_upto(8 * K + 1)
    # exponents via Legendre counting
    exps = {}
    for p in ps:
        e, q = 0, p
        while q <= 8 * K + 1:
            c = count_ap(K, q)
            if c == 0:
                break
            e += c
            q *= p
        if e:
            exps[p] = e
    # exponent-bit scheduling
    maxbit = max(e.bit_length() for e in exps.values())
    acc = mpz(1)
    for j in range(maxbit - 1, -1, -1):
        acc = (acc * acc) % MOD
        Lj = product_tree_mod([p for p, e in exps.items() if (e >> j) & 1], MOD)
        acc = (acc * Lj) % MOD
    return acc % MOD


def direct_product(K, m):
    MOD = mpz(1) << m
    return product_tree_mod([8 * k + 1 for k in range(1, K + 1)], MOD)


for K in (1000, 10000, 60000):
    m = int(2.4 * K) + 64
    t0 = time.perf_counter()
    a = borwein_product(K, m)
    tb = time.perf_counter() - t0
    t0 = time.perf_counter()
    b = direct_product(K, m)
    td = time.perf_counter() - t0
    print('K=%-7d m=%-7d %s  borwein %.3fs  direct %.3fs  (ratio %.2fx)'
          % (K, m, 'MATCH' if a == b else '*** MISMATCH ***', tb, td, td / tb if tb > 0 else 0))
