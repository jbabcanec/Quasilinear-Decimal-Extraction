"""Validate the block-anchor reduction (the pairing-kernel derivation):

  sum_{tau<b} 2^{-4 tau} (mu0 + 8 tau)^{-1}
    ==  mu0^{-1} * sum_j z^j W_j,   z = -8 mu0^{-1},  W_j = sum_{tau<b} 2^{-4 tau} tau^j

2-adically mod 2^m (scaled by 2^{4(b-1)} to clear negative powers). If exact,
per-family V = sum_i w_i W(z_i) with ONE fixed small rational W -- the
tapered-pairing formulation is sound.
"""
m = 512
MOD = 1 << m
b = 40
mu0 = 8 * 12345 + 1

# LHS: scaled exact sum  sum_tau 2^{4(b-1-tau)} * inverse(mu0+8tau)
lhs = 0
for tau in range(b):
    lhs = (lhs + (1 << (4 * (b - 1 - tau))) * pow(mu0 + 8 * tau, -1, MOD)) % MOD

# RHS: mu0^{-1} sum_j z^j * Ws_j,  Ws_j = sum_tau 2^{4(b-1-tau)} tau^j  (scaled W)
inv0 = pow(mu0, -1, MOD)
z = (-8 * inv0) % MOD
rhs = 0
zj = 1
j = 0
while True:
    Wj = sum((1 << (4 * (b - 1 - tau))) * pow(tau, j, MOD) for tau in range(b)) % MOD
    rhs = (rhs + zj * Wj) % MOD
    j += 1
    zj = (zj * z) % MOD
    if 3 * j > m + 8:   # v2(z^j) >= 3j
        break
rhs = (rhs * inv0) % MOD

print('block-anchor identity:', 'EXACT' if lhs == rhs else 'MISMATCH',
      '(%d j-terms, m=%d)' % (j, m))
