"""Validation: sieve-batched exponentiation along a modulus progression.

Claim: {2^{4k} mod (8k+1) : k <= K} is computable in O(K loglog K) mulmods via
  - segmented factorization of the AP by sieve,
  - per-prime-power geometric chains (2^{4k} mod q^e steps by *2^{4 q^e}),
  - CRT recombination,
instead of ~K * log(exponent) mulmods for per-term modexp.

This validates correctness against direct pow() and counts operations.
"""
import sys
from math import gcd

K = 20000
M = 8 * K + 1

# smallest-prime-factor sieve up to M
spf = list(range(M + 1))
i = 2
while i * i <= M:
    if spf[i] == i:
        for j in range(i * i, M + 1, i):
            if spf[j] == j:
                spf[j] = i
    i += 1

mulmods_batched = 0

# chain caches: value[(q^e)] -> dict k -> 2^{4k} mod q^e along k ≡ k0 (mod q^e)
chain_val = {}

def build_chain(P, is_prime):
    """geometric chain of 2^{4k} mod P over k ≡ k0 (mod P), k <= K."""
    global mulmods_batched
    inv8 = pow(8, -1, P)
    k0 = (-inv8) % P            # 8k+1 ≡ 0 (mod P)
    vals = {}
    if k0 > K:
        return vals
    if is_prime:
        # Euler starter, O(1): 8*k0+1 = P*s with s < 8 odd, so
        # 2^{4k0} = 2^{(P*s-1)/2} = (2|P)^s * 2^{(s-1)/2}  (mod P)
        s = (8 * k0 + 1) // P
        leg = 1 if P % 8 in (1, 7) else P - 1      # (2|P) for odd prime P
        g = (pow(leg, s, P) * (1 << ((s - 1) >> 1))) % P
        mulmods_batched += 4
    else:
        g = pow(2, 4 * k0, P)   # rare prime-power case: real modexp
        mulmods_batched += 2 * (4 * k0).bit_length() if k0 else 1
    # chain constant 2^{4P} mod P: for prime P, Fermat gives
    # 2^{4P} = (2^{P-1})^4 * 2^4 ≡ 16 (mod P)  -- O(1), no modexp.
    if is_prime:
        c = 16 % P
        mulmods_batched += 1
    else:
        c = pow(2, 4 * P, P)
        mulmods_batched += 2 * (4 * P).bit_length()
    k = k0
    while k <= K:
        vals[k] = g
        g = (g * c) % P         # one mulmod per incidence
        mulmods_batched += 1
        k += P
    return vals

def factor_mu(mu):
    f = {}
    while mu > 1:
        q = spf[mu]
        e = 0
        while mu % q == 0:
            mu //= q
            e += 1
        f[q] = e
    return f

bad = 0
checked = 0
for k in range(1, K + 1):
    mu = 8 * k + 1
    # batched value via CRT over prime powers
    r, mod = 0, 1
    for q, e in factor_mu(mu).items():
        P = q ** e
        if P not in chain_val:
            chain_val[P] = build_chain(P, e == 1)
        vP = chain_val[P][k]
        # CRT combine (r mod mod, vP mod P)
        inv = pow(mod % P, -1, P)
        r = r + mod * (((vP - r) * inv) % P)
        mod *= P
        mulmods_batched += 3
    want = pow(2, 4 * k, mu)
    checked += 1
    if r % mu != want:
        bad += 1
        if bad < 5:
            print('MISMATCH k=%d mu=%d got=%d want=%d' % (k, mu, r % mu, want))

# baseline cost: per-term modexp, ~2 mulmods per exponent bit
mulmods_direct = sum(2 * (4 * k).bit_length() for k in range(1, K + 1))

# Euler bonus check: prime mu => value 1
euler_bad = 0
for k in range(1, K + 1):
    mu = 8 * k + 1
    if spf[mu] == mu and pow(2, 4 * k, mu) != 1:
        euler_bad += 1

print('checked %d terms: %s' % (checked, 'ALL CORRECT' if bad == 0 else '%d BAD' % bad))
print('Euler criterion (prime mu -> 2^{4k} == 1): %s'
      % ('CONFIRMED' if euler_bad == 0 else '%d violations' % euler_bad))
print('mulmods: batched %d vs direct %d  -> %.1fx fewer'
      % (mulmods_batched, mulmods_direct, mulmods_direct / mulmods_batched))
print('incidence count Sum omega(mu): %d  (= %.2f per term; loglog(8K) ~ %.2f)'
      % (sum(len(factor_mu(8 * j + 1)) for j in range(1, K + 1)),
         sum(len(factor_mu(8 * j + 1)) for j in range(1, K + 1)) / K,
         __import__('math').log(__import__('math').log(8 * K))))
sys.exit(1 if bad else 0)
