"""Validate the Bezout band-term split and the dyadic consolidation.

Per-term identity:  frac(X/(mu 2^J)) = frac( (R*u mod mu)/mu + (X*w mod 2^J)/2^J )
with R = X mod mu, u = 2^{-J} mod mu, w = mu^{-1} mod 2^J.

Consolidation: sum over terms of the dyadic parts equals
frac( X * V / 2^{Jmax} ),  V = sum_t s_t * 2^{Jmax-J_t} * (mu_t^{-1} mod 2^{Jmax})  mod 2^{Jmax}.
"""
import random
from fractions import Fraction

rnd = random.Random(7)
N = 601
n = N - 1
X = 5 ** n

terms = []
for _ in range(300):
    k = rnd.randrange(N // 4 + 1, int(0.8 * N))
    J = 4 * k - (n + 2)
    if J < 1:
        continue
    mu = 8 * k + 1
    s = rnd.choice((1, -1))
    terms.append((s, J, mu))

# per-term identity check (exact rationals)
bad = 0
for s, J, mu in terms:
    lhs = Fraction(X, mu << J) % 1
    R = X % mu
    u = pow(2, -J, mu)
    w = pow(mu, -1, 1 << J)
    rhs = (Fraction((R * u) % mu, mu) + Fraction((X * w) % (1 << J), 1 << J)) % 1
    if lhs != rhs:
        bad += 1
print('per-term Bezout split: %s (%d terms)' % ('ALL EXACT' if bad == 0 else '%d BAD' % bad, len(terms)))

# consolidation check: sum of dyadic parts vs single X*V/2^Jmax, mod 1, to 2^-160
Jmax = max(J for _, J, _ in terms)
V = 0
direct = Fraction(0)
for s, J, mu in terms:
    w_full = pow(mu, -1, 1 << Jmax)
    V = (V + s * ((1 << (Jmax - J)) * w_full)) % (1 << Jmax)
    direct += s * Fraction((X * pow(mu, -1, 1 << J)) % (1 << J), 1 << J)
direct %= 1
consolidated = Fraction((X * V) % (1 << Jmax), 1 << Jmax)
err = abs(direct - consolidated) % 1
err = min(err, 1 - err)
print('consolidation |direct - X*V/2^Jmax| mod 1 = %s (should be 0)' % ('0 EXACT' if err == 0 else float(err)))
