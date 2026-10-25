#!/usr/bin/env python3
"""
egf_power_sums.py -- batch power sums via EGF closed form in 2-adic FS.

    sum_{n>=0} S_n t^n / n!  =  sum_{k=1}^K e^{k t}
                              =  e^t (e^{K t} - 1) / (e^t - 1)

If this truncated series is computable in O(M(m) log m) over Z/2^m,
all S_0..S_{J} (J~m/3) fall out, and

    H = sum_j (-8)^j S_j

is then O(M(m) log m) after the Borel-weight j! is applied — TWO LOGS.

Obstacles at p=2: exp series, division by (e^t-1) ~ t, factorial vals.
"""
from __future__ import annotations
import time

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def v2(n):
    n = int(n)
    return 999 if n == 0 else (n & -n).bit_length() - 1


def odd_inv(a, m):
    return pow(mpz(a), -1, 1 << m)


class FS:
    """Truncated formal power series over Z/2^m Z, length L (coeffs 0..L-1)."""

    __slots__ = ("c", "m", "L", "MOD")

    def __init__(self, coeffs, m, L=None):
        self.m = m
        self.MOD = mpz(1) << m
        if L is None:
            L = len(coeffs)
        self.L = L
        self.c = [mpz(coeffs[i]) % self.MOD if i < len(coeffs) else mpz(0)
                  for i in range(L)]

    def copy(self):
        return FS(self.c, self.m, self.L)

    def __add__(self, other):
        assert self.m == other.m and self.L == other.L
        return FS([(self.c[i] + other.c[i]) % self.MOD for i in range(self.L)],
                  self.m, self.L)

    def __sub__(self, other):
        assert self.m == other.m and self.L == other.L
        return FS([(self.c[i] - other.c[i]) % self.MOD for i in range(self.L)],
                  self.m, self.L)

    def __mul__(self, other):
        """Naive O(L^2); swap to NTT later if validated."""
        assert self.m == other.m and self.L == other.L
        out = [mpz(0)] * self.L
        for i in range(self.L):
            if self.c[i] == 0:
                continue
            for j in range(self.L - i):
                out[i + j] = (out[i + j] + self.c[i] * other.c[j]) % self.MOD
        return FS(out, self.m, self.L)

    def scale(self, s):
        s = mpz(s) % self.MOD
        return FS([(a * s) % self.MOD for a in self.c], self.m, self.L)

    def shift(self, k):
        """Multiply by t^k."""
        out = [mpz(0)] * self.L
        for i in range(self.L - k):
            out[i + k] = self.c[i]
        return FS(out, self.m, self.L)


def fs_exp(x: FS) -> FS:
    """exp(x) for x with x(0)=0 and v2(coeffs) high enough, or general via
    sum x^n/n! truncated.  For p=2 need care on /n!."""
    m, L, MOD = x.m, x.L, x.MOD
    # exp = sum_{n=0}^{L-1} x^n / n!
    acc = FS([1] + [0] * (L - 1), m, L)  # 1
    if all(v == 0 for v in x.c):
        return acc
    pw = FS([1] + [0] * (L - 1), m, L)
    fact = mpz(1)
    for n in range(1, L):
        pw = pw * x
        # / n!  = / (fact * n) with fact=(n-1)!
        # compute inv of n! tracking 2-valuation
        # n! = 1*2*...*n
        # For coeff-wise: each coeff of pw / n!
        e = 0
        odd = 1
        # v2(n!) and odd part
        # better: maintain fact_odd and v2fact
        # simpler approach: divide by n iteratively
        # pw_n = pw_{n-1} * x / n
        # Restart: iterative Horner-like
        break
    # Iterative: e = 1; term = 1; for n=1..: term *= x/n; e += term
    e = FS([1] + [0] * (L - 1), m, L)
    term = FS([1] + [0] * (L - 1), m, L)
    for n in range(1, L):
        term = term * x
        # divide term by n
        vn = v2(n)
        oddn = n >> vn
        inv_odd = odd_inv(oddn, m)
        newc = []
        ok = True
        for coeff in term.c:
            vc = v2(coeff) if coeff else 999
            if vc < vn:
                # can't divide exactly in Z_2 — series may still work if
                # those coeffs should be zero at this precision
                newc.append(mpz(0))
                ok = False
            else:
                newc.append((((coeff >> vn) % MOD) * inv_odd) % MOD)
        term = FS(newc, m, L)
        e = e + term
    return e


def fs_inv(x: FS) -> FS:
    """Inverse of series with x(0) odd unit, Newton."""
    m, L, MOD = x.m, x.L, x.MOD
    assert x.c[0] % 2 == 1
    # Newton: y <- y (2 - x y), start y0 = inv(x0)
    y = FS([odd_inv(x.c[0], m)] + [0] * (L - 1), m, L)
    n = 1
    while n < L:
        n = min(L, 2 * n)
        # truncate mul to n for speed later; for now full
        xy = x * y
        two_minus = FS([(mpz(2) if i == 0 else mpz(0)) - xy.c[i]
                        for i in range(L)], m, L)
        # fix mod
        two_minus = FS([c % MOD for c in two_minus.c], m, L)
        y = y * two_minus
    return y


def egf_power_sum_series(K, m, L):
    """G(t) = sum_{k=1}^K exp(k t) as FS length L over Z/2^m.
    = exp(t) * (exp(K t) - 1) * inv(exp(t) - 1)
    """
    # t as series
    t = FS([0, 1] + [0] * (L - 2), m, L)
    # Kt
    Kt = t.scale(K)
    et = fs_exp(t)
    eKt = fs_exp(Kt)
    # et - 1
    etm1 = FS([(et.c[0] - 1) % et.MOD] + et.c[1:], m, L)
    # (et - 1)/t = 1 + t/2! + t^2/3! + ... has constant 1 (odd) after /t?
    # Better: inv(et-1) carefully — et-1 has zero constant, val_t = 1.
    # Write et-1 = t * u(t) with u(0)=1, inv(et-1) = inv(t) * inv(u) — not FS.
    # So compute (eKt - 1) / (et - 1) as series using the identity
    # = sum_{j=0}^{K-1} exp(j t)
    # For general K this is still a long sum — BAD.
    #
    # Use: (eKt-1)/(et-1) = quot of series after canceling t:
    # eKt-1 = Kt * uK, et-1 = t * u1, ratio = K * uK * inv(u1)
    u1 = FS(etm1.c[1:] + [0], m, L)  # (et-1)/t , length L — drops last
    # fix length
    u1 = FS([etm1.c[i + 1] if i + 1 < L else 0 for i in range(L)], m, L)
    eKtm1 = FS([(eKt.c[0] - 1) % eKt.MOD] + eKt.c[1:], m, L)
    uK = FS([eKtm1.c[i + 1] if i + 1 < L else 0 for i in range(L)], m, L)
    # ratio_without_K = uK * inv(u1); full = K * that? 
    # (eKt-1)/(et-1) = (Kt uK)/(t u1) = K uK/u1
    inv_u1 = fs_inv(u1)
    ratio = uK * inv_u1
    ratio = ratio.scale(K)
    # G = et * ratio? 
    # sum_{k=1}^K e^{k t} = e^t (e^{Kt}-1)/(e^t-1) = et * (eKt-1)/(et-1)
    # = et * K * uK/u1
    G = et * ratio
    return G


def extract_S_from_egf(G: FS, J, m):
    """S_n = n! * [t^n] G(t), for n=0..J."""
    S = []
    fact = mpz(1)
    MOD = mpz(1) << m
    for n in range(J + 1):
        # S_n = n! * G.c[n]
        # compute n! * G[n] with 2-adic care
        if n >= 1:
            # fact = n!
            pass
        # maintain fact
        if n == 0:
            fact = mpz(1)
        else:
            fact = (fact * n)  # exact integer factorial — may be huge
        # We only need mod 2^m, so fact mod 2^m * (accounting if v2(fact)>=m -> 0)
        if v2(fact) >= m:
            S.append(mpz(0))
        else:
            # fact = 2^e * odd, G[n] needs to absorb
            e = v2(fact)
            odd = fact >> e
            gn = G.c[n] if n < G.L else mpz(0)
            # S = fact * gn = 2^e * odd * gn. Need v2(gn)>=0 always.
            S.append((odd * gn * (mpz(1) << e)) % MOD)
            # wait that's fact*gn if no overflow; use
            S[-1] = (mpz(fact % (1 << (m + e + 2))) * gn) % MOD
            # simpler:
            S[-1] = (mpz(int(fact)) * gn) % MOD if fact.bit_length() < 2 * m + 10 else mpz(0)
    return S


def part_validate_exp():
    P("=" * 70)
    P("PART A: fs_exp smoke test")
    m, L = 64, 20
    # exp(8) via series for unit... use x = 8 as constant series? 
    # exp of series with only t^1 coeff = 8: exp(8t)
    x = FS([0, 8] + [0] * (L - 2), m, L)
    e = fs_exp(x)
    # [t^0]=1, [t^1]=8, [t^2]=8^2/2=32, [t^3]=8^3/6= 512/6 not integer in Q but
    # 8^3/6 = 512/6 = 256/3 — not in Z_2! 
    # So exp(8t) has non-2-adic-integer coeffs. Problem!
    P(f"  exp(8t) coeffs: {[hex(c) for c in e.c[:6]]}")
    P("  NOTE: 8^n/n! not always in Z_2 — EGF over Q_2 needs valuation tracking.")


def part_diagnose_vals():
    P()
    P("=" * 70)
    P("PART B: v2(8^n / n!) and v2(k^n / n!) for EGF integrality")
    for n in range(0, 25):
        # v2(8^n / n!) = 3n - v2(n!)
        vf = 0
        pp = n
        while pp:
            pp //= 2
            vf += pp
        P(f"  n={n:2d}  v2(8^n/n!)={3*n - vf:3d}  v2(1^n/n!)={0-vf:3d}")


def main():
    P("egf_power_sums -- EGF attack on batch power sums")
    part_validate_exp()
    part_diagnose_vals()
    P()
    P("DIAGNOSIS: exp(t) as EGF for sum e^{k t} has coeffs 1/n! not in Z_2.")
    P("exp(8t) has 8^n/n! with v2=3n-v2(n!) >= n > 0 — OK in Z_2 for n>=1.")
    P("But sum_k exp(k t) needs exp(t) with 1/n! — NOT a Z_2 series.")
    P("Must work in Q_2 with bounded denom, or use OGF not EGF.")
    P()
    P("OGF pivot: sum_n S_n z^n = sum_k 1/(1 - k z) for |kz|<1,")
    P("  = (1/z) (psi(1/z) - psi(1/z - K)) — digamma again.")
    P("  2-adic: sum_k inv(1 - k z) for z=-8: sum inv(1+8k)=H tautology.")


if __name__ == "__main__":
    main()
