#!/usr/bin/env python3
"""
Decimal digit extraction for pi (and ln 2) from BINARY BBP formulas
===================================================================

Method ("three-phase split"). Write, using the classic hex BBP formula

    pi = sum_k 16^-k [ 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6) ],

the shifted value as a sum of signed elementary terms

    10^(N-1) * pi = sum  +/- 2^E * 5^(N-1) / mu        (mu odd),

folding all powers of two (numerator coefficients 4, 2, 1 and even parts of
the denominators 8k+4, 8k+6) into the exponent E = E(k, j).

Phase 1 (easy, E >= 0): the numerator is an integer, so mod 1 the term is
    (2^E * 5^(N-1) mod mu) / mu  --  two word-sized modexps per term.

Phase 2 (band, E < 0): term = X / (mu * 2^J) with X = 5^(N-1), J = -E > 0.
    Exact identity (h in [0, mu), u in [0, 1)):
        frac(X / (mu * 2^J)) = (h + u) / mu
        h = floor(X / 2^J) mod mu
          = ( (X mod mu) - ((X mod 2^J) mod mu) ) * 2^(-J)  mod mu
        u = frac(X / 2^J)   -- a ~Gp-bit window into X's bit-string at depth J
    Every band term shares the SINGLE precomputed bit-string of X = 5^(N-1).
    Per-term cost: word-sized modexps + an O(1) window read + (v1 only) one
    big-prefix mod; the prefix mods are batchable to polylog total via a
    prefix-remainder tree (see research log Part 16).

Phase 3 (tail, J >= bits(X)): floor part is 0; window read only.

No digits of pi before position N are ever computed.
Auxiliary memory: the bit-string of 5^(N-1), ~0.2903 * N bytes.

Usage:
  python decimal_bbp_extract.py --pos 100000              # pi digits at N
  python decimal_bbp_extract.py --const ln2 --pos 10000   # ln 2 digits at N
  python decimal_bbp_extract.py --selftest                # verify vs reference
  python decimal_bbp_extract.py --selftest --big          # add N = 1e5
  python decimal_bbp_extract.py --audit                   # exact per-term audit
  python decimal_bbp_extract.py --bench                   # scaling benchmark
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from math import isqrt

L2_10 = math.log2(10)

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0x7fffffff)

_P10 = {}


def _pow10(e):
    v = _P10.get(e)
    if v is None:
        v = _P10[e] = 10 ** e
    return v


def int_to_dec(x, ndigits):
    """Exactly ndigits decimal chars of non-negative x (zero-padded), via
    divide-and-conquer -- CPython's builtin conversion is quadratic."""
    if ndigits <= 4000:
        return format(x, '0%dd' % ndigits)
    half = ndigits >> 1
    hi, lo = divmod(x, _pow10(half))
    return int_to_dec(hi, ndigits - half) + int_to_dec(lo, half)


# ----------------------------------------------------------------------------
# Term generators: yield (sign, E, mu) with term = sign * 2^E * 5^(N-1) / mu
# ----------------------------------------------------------------------------

def bbp_pi_terms(N, K):
    """Hex BBP formula, decimal-shifted. A = N-1-4k."""
    for k in range(K + 1):
        A = N - 1 - 4 * k
        yield (1, A + 2, 8 * k + 1)   # +4/(8k+1)
        yield (-1, A - 1, 2 * k + 1)  # -2/(8k+4) = -1/(2(2k+1)) -> 2^(A-1)/(2k+1)
        yield (-1, A, 8 * k + 5)      # -1/(8k+5)
        yield (-1, A - 1, 4 * k + 3)  # -1/(8k+6) = -1/(2(4k+3))
    return


def ln2_terms(N, K):
    """ln 2 = sum_{k>=1} 1/(k 2^k), decimal-shifted. k = 2^s * mu."""
    for k in range(1, K + 1):
        s = (k & -k).bit_length() - 1
        yield (1, N - 1 - k - s, k >> s)
    return


def cutoff(N, G, const):
    if const == 'pi':
        return int(((N - 1) * L2_10 + G + 24) / 4) + 2
    return int((N - 1) * L2_10 + G + 24) + 2


# ----------------------------------------------------------------------------
# Core extractor
# ----------------------------------------------------------------------------

def extract(N, const='pi', G=192, Gp=144, window=12, stats=None):
    """Return `window` decimal digits of `const` at positions N..N+window-1
    (1-indexed after the decimal point), plus a stats dict."""
    K = cutoff(N, G, const)
    terms = bbp_pi_terms(N, K) if const == 'pi' else ln2_terms(N, K)

    t0 = time.perf_counter()
    n5 = N - 1
    X = pow(5, n5)
    LX = X.bit_length()
    Xb = X.to_bytes((LX + 7) >> 3, 'little')
    t_pre = time.perf_counter() - t0

    MOD = 1 << G
    GPMASK = (1 << Gp) - 1
    acc = 0
    n_easy = n_band = n_tail = 0

    t0 = time.perf_counter()
    for sign, E, mu in terms:
        if E >= 0:
            n_easy += 1
            if mu == 1:
                continue
            r = (pow(2, E, mu) * pow(5, n5, mu)) % mu
            if r == 0:
                continue
            f = (r << G) // mu
        else:
            J = -E
            if J >= LX + G + 8:
                continue  # term below accumulator resolution
            if J >= LX:
                n_tail += 1
                h = 0
                uG = (X >> (J - Gp)) if J >= Gp else (X << (Gp - J))
            else:
                n_band += 1
                nby = (J + 7) >> 3
                Ylow = int.from_bytes(Xb[:nby], 'little') & ((1 << J) - 1)
                if mu == 1:
                    h = 0
                else:
                    R = pow(5, n5, mu)
                    h = ((R - Ylow % mu) * pow((mu + 1) >> 1, J, mu)) % mu
                if J >= Gp:
                    o = J - Gp
                    lo = o >> 3
                    chunk = int.from_bytes(Xb[lo:lo + (Gp >> 3) + 3], 'little') >> (o & 7)
                    uG = chunk & GPMASK
                else:
                    uG = Ylow << (Gp - J)
            f = ((h << G) + (uG << (G - Gp))) // mu
        acc = (acc + sign * f) % MOD
    t_loop = time.perf_counter() - t0

    digs = format((acc * 10 ** window) >> G, '0%dd' % window)
    st = dict(pre=t_pre, loop=t_loop, mem_bytes=len(Xb), K=K,
              easy=n_easy, band=n_band, tail=n_tail)
    if stats is not None:
        stats.update(st)
    return digs, st


# ----------------------------------------------------------------------------
# v2: batched band reduction via one accumulating remainder tree
# (paper Theorem 7 / research log 16.4-slim; Costa-Gerbicz-Harvey style)
# ----------------------------------------------------------------------------

def _bits_at(Xb, pos, nbits):
    """Bits [pos, pos+nbits) of the little-endian byte string Xb, zero-filled
    outside [0, 8*len(Xb)); pos may be negative (low zeros)."""
    if nbits <= 0:
        return 0
    if pos < 0:
        return _bits_at(Xb, 0, nbits + pos) << (-pos)
    lo = pos >> 3
    if lo >= len(Xb):
        return 0
    chunk = int.from_bytes(Xb[lo:lo + ((nbits + 7) >> 3) + 1], 'little')
    return (chunk >> (pos & 7)) & ((1 << nbits) - 1)


def batched_prefix_mods(Xb, Js, mus):
    """For sorted J_0 <= J_1 <= ... and word moduli mus, return lists
    (Y, S) with Y[t] = (X mod 2^{J_t}) mod mu_t and S[t] = 2^{J_t} mod mu_t,
    where X is given by its little-endian bytes Xb.

    One accumulating remainder tree over products of the 2x2 upper-triangular
    matrices [[1, c_t], [0, 2^{d_t}]], d_t = J_{t+1}-J_t, c_t = bits
    [J_t, J_{t+1}) of X: the pair v_t = (X mod 2^{J_t}, 2^{J_t}) satisfies
    v_{t+1} = A_t v_t, and leaf t reduces the running product mod mu_t.
    Cost O(M(L + T*w) log T); workspace = the two trees (block/forest variant
    for production runs is future work)."""
    T = len(mus)
    if T == 0:
        return [], []
    Ys = [0] * T
    Ss = [0] * T

    tree = {}

    def build(lo, hi):
        if hi - lo == 1:
            if lo < T - 1:
                d = Js[lo + 1] - Js[lo]
                mat = (_bits_at(Xb, Js[lo], d), 1 << d)
            else:
                mat = (0, 1)  # identity; A_{T-1} is never used in a prefix
            node = (mat, mus[lo])
        else:
            mid = (lo + hi) >> 1
            (bE, dE), ML = build(lo, mid)   # earlier block
            (bL, dL), MR = build(mid, hi)   # later block
            node = ((bE + bL * dE, dL * dE), ML * MR)
        tree[(lo, hi)] = node
        return node

    build(0, T)

    def down(lo, hi, p, s):
        # invariant: (p, s) === (X mod 2^{J_lo}, 2^{J_lo})  (mod Mprod[lo,hi))
        if hi - lo == 1:
            mu = mus[lo]
            Ys[lo] = p % mu
            Ss[lo] = s % mu
            return
        mid = (lo + hi) >> 1
        (b, d), ML = tree[(lo, mid)]
        _, MR = tree[(mid, hi)]
        down(lo, mid, p % ML, s % ML)
        down(mid, hi, (p + b * s) % MR, (d * s) % MR)

    down(0, T, _bits_at(Xb, 0, Js[0]), 1 << Js[0])
    return Ys, Ss


def extract_v2(N, const='pi', G=192, Gp=144, window=12):
    """Same output as extract(), but all band/tail prefix residues come from
    one batched accumulating remainder tree instead of per-term reductions."""
    K = cutoff(N, G, const)
    terms = bbp_pi_terms(N, K) if const == 'pi' else ln2_terms(N, K)

    n5 = N - 1
    X = pow(5, n5)
    LX = X.bit_length()
    Xb = X.to_bytes((LX + 7) >> 3, 'little')

    MOD = 1 << G
    acc = 0
    band = []  # (J, mu, sign)
    for sign, E, mu in terms:
        if E >= 0:
            if mu == 1:
                continue
            r = (pow(2, E, mu) * pow(5, n5, mu)) % mu
            if r:
                acc = (acc + sign * ((r << G) // mu)) % MOD
        else:
            J = -E
            if J < LX + G + 8:
                band.append((J, mu, sign))

    band.sort(key=lambda t: t[0])
    Js = [t[0] for t in band]
    mus = [t[1] for t in band]
    Ys, Ss = batched_prefix_mods(Xb, Js, mus)

    for (J, mu, sign), Yt, St in zip(band, Ys, Ss):
        if mu == 1:
            h = 0
        else:
            R = pow(5, n5, mu)
            h = ((R - Yt) * pow(St, -1, mu)) % mu
        o = J - Gp
        if o >= 0:
            uG = _bits_at(Xb, o, Gp)
        else:
            uG = _bits_at(Xb, 0, J) << (-o)
        f = ((h << G) + (uG << (G - Gp))) // mu
        acc = (acc + sign * f) % MOD

    return format((acc * 10 ** window) >> G, '0%dd' % window)


# ----------------------------------------------------------------------------
# v3: the refined pipeline (paper Sec. 4) -- Bezout split + division-free
# 2-adic consolidation. All band-term dyadic parts merge into ONE integer V
# mod 2^Jmax via truncated multiplications only (no big divisions), then one
# Newton-Hensel inversion and one product X*V finish the job.
# Uses gmpy2 (GMP) transparently when available for quasi-linear multiplies.
# ----------------------------------------------------------------------------

try:
    from gmpy2 import mpz as _mpz
    HAVE_GMPY2 = True
except ImportError:
    _mpz = int
    HAVE_GMPY2 = False


def extract_v3(N, const='pi', G=192, window=12, stats=None):
    """Same digits as extract(); band terms via the Bezout/2-adic normal form."""
    K = cutoff(N, G, const)
    terms = bbp_pi_terms(N, K) if const == 'pi' else ln2_terms(N, K)

    t0 = time.perf_counter()
    n5 = N - 1
    X = _mpz(5) ** n5 if HAVE_GMPY2 else pow(5, n5)
    LX = X.bit_length()
    t_pre = time.perf_counter() - t0

    MOD = 1 << G
    acc = 0
    band = []
    t0 = time.perf_counter()
    for sign, E, mu in terms:
        if E >= 0:
            if mu == 1:
                continue
            r = (pow(2, E, mu) * pow(5, n5, mu)) % mu
            if r:
                acc = (acc + sign * ((r << G) // mu)) % MOD
        else:
            J = -E
            if J < LX + G + 8:
                band.append((sign, J, mu))
                if mu > 1:
                    # odd half of the Bezout split: word arithmetic
                    R = pow(5, n5, mu)
                    u = pow((mu + 1) >> 1, J, mu)   # 2^{-J} mod mu
                    ro = (R * u) % mu
                    if ro:
                        acc = (acc + sign * ((ro << G) // mu)) % MOD
    t_word = time.perf_counter() - t0

    # dyadic halves -> single V mod 2^Jmax, division-free balanced merge.
    # Nodes are (num, den, shift): value = 2^shift * num/den, with num, den
    # capped at their own needed precision 2^(Jmax-shift). Leaves are word
    # size; shifts align only at merges, so memory stays O(N log N) bits.
    t0 = time.perf_counter()
    Jmax = max(J for _, J, _ in band) if band else 1

    def merge(a, b):
        nL, dL, sL = a
        nR, dR, sR = b
        s = sL if sL < sR else sR
        cap = Jmax - s
        n_ = ((nL * dR) << (sL - s)) + ((nR * dL) << (sR - s))
        d_ = dL * dR
        # mask lazily: building a cap-bit mask at every merge is quadratic churn
        if n_.bit_length() > cap:
            n_ &= (_mpz(1) << cap) - 1
        if d_.bit_length() > cap:
            d_ &= (_mpz(1) << cap) - 1
        return (n_, d_, s)

    def reduce_nodes(nodes):
        while len(nodes) > 1:
            nxt = [merge(nodes[i], nodes[i + 1])
                   for i in range(0, len(nodes) - 1, 2)]
            if len(nodes) & 1:
                nxt.append(nodes[-1])
            nodes = nxt
        return nodes[0]

    if band:
        band.sort(key=lambda t: t[1])
        BLK = 1 << 15
        roots = []
        for lo in range(0, len(band), BLK):
            leaves = [(_mpz(sg), _mpz(mu), Jmax - J)
                      for sg, J, mu in band[lo:lo + BLK]]
            roots.append(reduce_nodes(leaves))
        n_root, d_root, s_root = reduce_nodes(roots)
        cap = Jmax - s_root
        capm = (_mpz(1) << cap) - 1
        # Newton-Hensel inverse of odd d_root mod 2^cap
        inv = _mpz(1)
        prec = 1
        while prec < cap:
            prec = min(2 * prec, cap)
            pm = (_mpz(1) << prec) - 1
            inv = (inv * (2 - d_root * inv)) & pm
        V = (((n_root * inv) & capm) << s_root) & ((_mpz(1) << Jmax) - 1)
        W = (X * V) & ((_mpz(1) << Jmax) - 1)   # X*V mod 2^Jmax
        dy = int(W >> (Jmax - G)) if Jmax >= G else int(W << (G - Jmax))
        acc = (acc + dy) % MOD
    t_v = time.perf_counter() - t0

    st = dict(pre=t_pre, word=t_word, vphase=t_v,
              total=t_pre + t_word + t_v, gmp=HAVE_GMPY2,
              mem_bytes=(LX + 7) // 8)
    if stats is not None:
        stats.update(st)
    return format((acc * 10 ** window) >> G, '0%dd' % window), st


# ----------------------------------------------------------------------------
# v4: v3 + sieve-batched word phase (paper Thm. "sieve-batched exponentiation").
# Every word-phase value is v_k = 2^{E(k)} * 5^n mod mu(k) with E linear in k
# and mu on an AP. For each small prime q: the hits k ≡ k0 (mod q) form an AP
# along which 2^{E} mod q is GEOMETRIC with ratio 2^{-4q} ≡ 16^{-1} (Fermat);
# 5^n mod q is cached once. Large-prime cofactors use Fermat-reduced
# exponents. Streaming CRT reassembles each term; k is processed in chunks.
# ----------------------------------------------------------------------------

def _small_primes(limit):
    s = bytearray(b'\x01') * (limit + 1)
    s[0:2] = b'\x00\x00'
    i = 2
    while i * i <= limit:
        if s[i]:
            s[i * i::i] = b'\x00' * len(s[i * i::i])
        i += 1
    return [i for i in range(3, limit + 1) if s[i]]   # odd primes only


def extract_v4(N, const='pi', G=192, window=12, stats=None):
    """v3 pipeline with the sieve-batched word phase. pi only."""
    assert const == 'pi'
    K = cutoff(N, G, const)
    t0 = time.perf_counter()
    n5 = N - 1
    X = _mpz(5) ** n5 if HAVE_GMPY2 else pow(5, n5)
    LX = X.bit_length()
    t_pre = time.perf_counter() - t0

    MOD = 1 << G
    acc = 0
    band = []          # (sign, J, mu) for the V phase
    Jcut = LX + G + 8

    # families: mu = a*k+b, E = e0 - 4k, sign
    # e0 follows bbp_pi_terms exactly: A = n5 - 4k; E = A+2, A-1, A, A-1.
    FAMS = ((1, 8, 1, n5 + 2), (-1, 2, 1, n5 - 1), (-1, 8, 5, n5),
            (-1, 4, 3, n5 - 1))

    t0 = time.perf_counter()
    Mstar = 8 * K + 6
    root = int(Mstar ** 0.5) + 1
    primes = _small_primes(root)

    CH = 1 << 19       # k-chunk size
    for sign, a, b, e0 in FAMS:
        # per-prime persistent chain state: q -> [next_k, g]  (g = 2^{E(next_k)} mod q)
        # plus caches inv16[q], pow5[q]; big: large-prime cache
        # c -> [last_k, last_g2, pow5, inv16]
        big = {}
        chain = {}
        inv16 = {}
        pow5c = {}
        for q in primes:
            if a % q == 0:
                continue
            k0 = (-b * pow(a, -1, q)) % q
            if k0 > K:
                continue
            chain[q] = [k0, pow(2, (e0 - 4 * k0) % (q - 1), q)]
            inv16[q] = pow(16, -1, q)
            # for q=5 no Euler reduction (5 not coprime); pow handles n5 < 1
            pow5c[q] = pow(5, n5 % (q - 1), q) if q != 5 else pow(5, n5, 5)
        for lo in range(0, K + 1, CH):
            hi = min(lo + CH, K + 1)
            nk = hi - lo
            rem = [a * k + b for k in range(lo, hi)]
            rr = [0] * nk          # CRT residue so far
            mm = [1] * nk          # CRT modulus so far
            for q, st in chain.items():
                k, g = st
                if k >= hi:
                    continue
                i16 = inv16[q]
                p5 = pow5c[q]
                while k < hi:
                    i = k - lo
                    mu = rem[i]
                    e = 1
                    mu //= q
                    while mu % q == 0:
                        mu //= q
                        e += 1
                    if e == 1:
                        v = (g * p5) % q
                        P = q
                    else:
                        P = q ** e
                        lam = P // q * (q - 1)
                        p5 = pow(5, n5, P) if q == 5 else pow(5, n5 % lam, P)
                        v = (pow(2, (e0 - 4 * k) % lam, P) * p5) % P
                    rem[i] = mu
                    # CRT combine (rr[i] mod mm[i]) with (v mod P)
                    m_ = mm[i]
                    rr[i] = rr[i] + m_ * (((v - rr[i]) * pow(m_ % P, -1, P)) % P)
                    mm[i] = m_ * P
                    g = (g * i16) % q      # exponent drops by 4q along the AP
                    k += q
                st[0], st[1] = k, g
            for i in range(nk):
                k = lo + i
                E = e0 - 4 * k
                if E < 0:
                    J = -E
                    if J >= Jcut:
                        continue
                    band.append((sign, J, a * k + b))
                mu_full = a * k + b
                if mu_full == 1:
                    continue
                c = rem[i]
                if c > 1:      # large prime cofactor
                    st_ = big.get(c)
                    if st_ is None:
                        # first occurrence: one modexp for each part
                        g2 = pow(2, E % (c - 1), c)
                        big[c] = [k, g2, pow(5, n5 % (c - 1), c),
                                  pow(16, -1, c)]
                        v = (g2 * big[c][2]) % c
                    else:
                        # later occurrence: exponent dropped by 4q per step of q
                        t = (k - st_[0]) // c
                        g2 = (st_[1] * pow(st_[3], t, c)) % c
                        st_[0], st_[1] = k, g2
                        v = (g2 * st_[2]) % c
                    m_ = mm[i]
                    rr[i] = rr[i] + m_ * (((v - rr[i]) * pow(m_ % c, -1, c)) % c)
                    mm[i] = m_ * c
                r = rr[i] % mu_full
                if r:
                    acc = (acc + sign * ((r << G) // mu_full)) % MOD
    t_word = time.perf_counter() - t0

    # ---- V phase: identical to extract_v3 ----
    t0 = time.perf_counter()
    Jmax = max(J for _, J, _ in band) if band else 1

    def merge(a_, b_):
        nL, dL, sL = a_
        nR, dR, sR = b_
        s = sL if sL < sR else sR
        cap = Jmax - s
        n_ = ((nL * dR) << (sL - s)) + ((nR * dL) << (sR - s))
        d_ = dL * dR
        if n_.bit_length() > cap:
            n_ &= (_mpz(1) << cap) - 1
        if d_.bit_length() > cap:
            d_ &= (_mpz(1) << cap) - 1
        return (n_, d_, s)

    def reduce_nodes(nodes):
        while len(nodes) > 1:
            nxt = [merge(nodes[i], nodes[i + 1])
                   for i in range(0, len(nodes) - 1, 2)]
            if len(nodes) & 1:
                nxt.append(nodes[-1])
            nodes = nxt
        return nodes[0]

    if band:
        band.sort(key=lambda t: t[1])
        BLK = 1 << 15
        roots = []
        for lo in range(0, len(band), BLK):
            leaves = [(_mpz(sg), _mpz(mu), Jmax - J)
                      for sg, J, mu in band[lo:lo + BLK]]
            roots.append(reduce_nodes(leaves))
        n_root, d_root, s_root = reduce_nodes(roots)
        cap = Jmax - s_root
        capm = (_mpz(1) << cap) - 1
        inv = _mpz(1)
        prec = 1
        while prec < cap:
            prec = min(2 * prec, cap)
            pm = (_mpz(1) << prec) - 1
            inv = (inv * (2 - d_root * inv)) & pm
        V = (((n_root * inv) & capm) << s_root) & ((_mpz(1) << Jmax) - 1)
        W = (X * V) & ((_mpz(1) << Jmax) - 1)
        dy = int(W >> (Jmax - G)) if Jmax >= G else int(W << (G - Jmax))
        acc = (acc + dy) % MOD
    t_v = time.perf_counter() - t0

    st = dict(pre=t_pre, word=t_word, vphase=t_v,
              total=t_pre + t_word + t_v, gmp=HAVE_GMPY2)
    if stats is not None:
        stats.update(st)
    return format((acc * 10 ** window) >> G, '0%dd' % window), st


# ----------------------------------------------------------------------------
# Independent cross-check: Dase 1844 identity pi/4 = atan(1/2)+atan(1/5)+atan(1/8)
# (2-5-smooth arguments; the 1/5 series band uses exact big-int fallback here)
# ----------------------------------------------------------------------------

def extract_dase(N, G=192, window=12):
    """Same digits via a completely different formula; slow exact band math.
    For cross-validation at small/medium N only."""
    MOD = 1 << G
    acc = 0
    n5 = N - 1
    # b = 2: 4*(-1)^k 2^(N-2k) 5^(N-1) / (2k+1); b = 8: exponent N-2-6k
    for (e_step, e0) in ((2, 0), (6, 2)):  # (per-k decrement, offset): E = N - e0 - e_step*k
        k = 0
        while True:
            E = N - e0 - e_step * k
            mu = 2 * k + 1
            mag_bits = -(E + n5 * math.log2(5)) + math.log2(mu)
            if mag_bits > G + 8:
                break
            if E >= 0:
                r = (pow(2, E, mu) * pow(5, n5, mu)) % mu
                f = (r << G) // mu
            else:
                den = mu << (-E)
                f = ((pow(5, n5) % den) << G) // den
            acc = (acc + (1 if k % 2 == 0 else -1) * f) % MOD
            k += 1
    # b = 5: 4*(-1)^k 2^(N+1) 5^(N-2-2k) / (2k+1)
    k = 0
    while True:
        F5 = N - 2 - 2 * k
        mu = 2 * k + 1
        mag_bits = -((N + 1) + F5 * math.log2(5)) + math.log2(mu)
        if mag_bits > G + 8:
            break
        if F5 >= 0:
            r = (pow(2, N + 1, mu) * pow(5, F5, mu)) % mu
            f = (r << G) // mu
        else:
            den = mu * 5 ** (-F5)
            f = ((pow(2, N + 1) % den) << G) // den
        acc = (acc + (1 if k % 2 == 0 else -1) * f) % MOD
        k += 1
    return format((acc * 10 ** window) >> G, '0%dd' % window)


# ----------------------------------------------------------------------------
# Reference digit generators (integer-only binary splitting)
# ----------------------------------------------------------------------------

def pi_reference(p):
    """First p decimal digits of pi after the point, via Chudnovsky."""
    n = p // 14 + 4
    C24 = 640320 ** 3 // 24

    def bs(a, b):
        if b - a == 1:
            if a == 0:
                P = Q = 1
            else:
                P = (6 * a - 5) * (2 * a - 1) * (6 * a - 1)
                Q = a * a * a * C24
            T = P * (13591409 + 545140134 * a)
            return P, Q, (-T if a & 1 else T)
        m = (a + b) // 2
        P1, Q1, T1 = bs(a, m)
        P2, Q2, T2 = bs(m, b)
        return P1 * P2, Q1 * Q2, T1 * Q2 + P1 * T2

    _, Q, T = bs(0, n)
    prec = p + 30
    sq = isqrt(10005 * 10 ** (2 * prec))
    val = 426880 * sq * Q // T
    s = int_to_dec(val, prec + 1)
    assert s[:6] == '314159', 'chudnovsky reference broken: ' + s[:20]
    return s[1:p + 1]


def ln2_reference(p):
    """First p decimal digits of ln 2 after the point, via 2*atanh(1/3)."""
    T = int(p / (2 * math.log10(3))) + 10

    def bs(a, b):  # sum_{k in [a,b)} 1/((2k+1) 9^(k-a)) as (num, den)
        if b - a == 1:
            return 1, 2 * a + 1
        m = (a + b) // 2
        n1, d1 = bs(a, m)
        n2, d2 = bs(m, b)
        nine = 9 ** (m - a)
        return n1 * d2 * nine + n2 * d1, d1 * d2 * nine

    n_, d_ = bs(0, T)
    val = 2 * n_ * 10 ** (p + 10) // (3 * d_)
    s = int_to_dec(val, p + 10)
    assert s[:6] == '693147', 'ln2 reference broken: ' + s[:20]
    return s[:p]


_REF_CACHE = {}


def ref_digits(const, p):
    key = (const, p)
    for (c, q), v in _REF_CACHE.items():
        if c == const and q >= p:
            return v[:p]
    v = pi_reference(p) if const == 'pi' else ln2_reference(p)
    _REF_CACHE[key] = v
    return v


# ----------------------------------------------------------------------------
# Test / audit / bench drivers
# ----------------------------------------------------------------------------

def run_selftest(big=False, huge=False):
    ok = True

    def check(const, N, ref, tag, m=8):
        nonlocal ok
        got, st = extract(N, const=const)
        want = ref[N - 1:N - 1 + m]
        good = got[:m] == want
        ok = ok and good
        print('  %-4s N=%-8d got=%s want=%s  %s   (%.3fs, %.1f MB)'
              % (const, N, got[:m], want, 'OK ' if good else '*** FAIL ***',
                 st['pre'] + st['loop'], st['mem_bytes'] / 1e6))
        return good

    print('[1] pi positions 1..60 (continuous, 8-digit windows)')
    ref = ref_digits('pi', 300)
    bad = 0
    for N in range(1, 61):
        got, _ = extract(N)
        if got[:8] != ref[N - 1:N + 7]:
            print('   FAIL at N=%d: got %s want %s' % (N, got[:8], ref[N - 1:N + 7]))
            bad += 1
            ok = False
    print('   %d/60 correct' % (60 - bad))

    print('[2] pi spot checks')
    ref = ref_digits('pi', 10100)
    for N in (100, 1000, 5000, 10000):
        check('pi', N, ref, 'spot')

    print('[3] Dase 1844 cross-check (independent formula, same digits)')
    for N in (1, 47, 100, 500, 1000):
        d1, _ = extract(N)
        d2 = extract_dase(N)
        good = d1[:10] == d2[:10]
        ok = ok and good
        print('   N=%-6d bbp16=%s dase=%s  %s' % (N, d1[:10], d2[:10],
                                                  'OK ' if good else '*** FAIL ***'))

    print('[4] ln2 positions 1..40 + spots')
    ref = ref_digits('ln2', 10100)
    bad = 0
    for N in range(1, 41):
        got, _ = extract(N, const='ln2')
        if got[:8] != ref[N - 1:N + 7]:
            print('   FAIL at N=%d: got %s want %s' % (N, got[:8], ref[N - 1:N + 7]))
            bad += 1
            ok = False
    print('   %d/40 correct' % (40 - bad))
    for N in (100, 1000, 10000):
        check('ln2', N, ref, 'spot')

    if big:
        print('[5] N = 100000 (pi + ln2)')
        ref = ref_digits('pi', 100100)
        check('pi', 100000, ref, 'big')
        ref = ref_digits('ln2', 100100)
        check('ln2', 100000, ref, 'big')

    if huge:
        print('[6] N = 1000000 (pi)  [reference gen takes a while]')
        t0 = time.perf_counter()
        ref = ref_digits('pi', 1000100)
        print('   reference generated in %.1fs' % (time.perf_counter() - t0))
        check('pi', 1000000, ref, 'huge')

    print('SELFTEST', 'PASSED' if ok else 'FAILED')
    return ok


def run_audit(N=997, samples=400):
    """Exact big-rational audit of individual band/tail terms."""
    import random
    rnd = random.Random(12345)
    K = cutoff(N, 192, 'pi')
    n5 = N - 1
    X = pow(5, n5)
    LX = X.bit_length()
    Xb = X.to_bytes((LX + 7) >> 3, 'little')
    G, Gp = 192, 144
    worst = 0
    cnt = 0
    all_terms = [t for t in bbp_pi_terms(N, K) if t[1] < 0]
    rnd.shuffle(all_terms)
    for sign, E, mu in all_terms[:samples]:
        J = -E
        if J >= LX + G + 8:
            continue
        # exact: frac(X/(mu 2^J)) in G-bit fixed point
        den = mu << J
        fe = ((X % den) << G) // den
        # windowed path (mirrors extract())
        if J >= LX:
            h = 0
            uG = (X >> (J - Gp)) if J >= Gp else (X << (Gp - J))
        else:
            Ylow = int.from_bytes(Xb[:(J + 7) >> 3], 'little') & ((1 << J) - 1)
            if mu == 1:
                h = 0
            else:
                R = pow(5, n5, mu)
                h = ((R - Ylow % mu) * pow((mu + 1) >> 1, J, mu)) % mu
            if J >= Gp:
                o = J - Gp
                chunk = int.from_bytes(Xb[o >> 3:(o >> 3) + (Gp >> 3) + 3], 'little') >> (o & 7)
                uG = chunk & ((1 << Gp) - 1)
            else:
                uG = Ylow << (Gp - J)
        fm = ((h << G) + (uG << (G - Gp))) // mu
        worst = max(worst, abs(fe - fm))
        cnt += 1
    print('audited %d band/tail terms at N=%d: worst fixed-point deviation = %d ulp(2^-%d)'
          % (cnt, N, worst, G))
    lim = 1 << (G - Gp + 2)
    print('bound  (should be < %d): %s' % (lim, 'OK' if worst < lim else '*** FAIL ***'))
    return worst < lim


def run_bench(big=False):
    sizes = [2000, 5000, 10000, 20000, 50000, 100000]
    if big:
        sizes += [200000, 500000, 1000000]
    rows = []
    for N in sizes:
        digs, st = extract(N)
        t = st['pre'] + st['loop']
        rows.append((N, t, st))
        print('N=%-8d %10.3fs  (pre %.3fs, loop %.3fs)  mem=%.2f MB  terms:%d easy/%d band/%d tail  -> %s'
              % (N, t, st['pre'], st['loop'], st['mem_bytes'] / 1e6,
                 st['easy'], st['band'], st['tail'], digs))
    if len(rows) >= 3:
        xs = [math.log(r[0]) for r in rows[1:]]
        ys = [math.log(r[1]) for r in rows[1:]]
        n = len(xs)
        sx, sy = sum(xs), sum(ys)
        sxx = sum(x * x for x in xs)
        sxy = sum(x * y for x, y in zip(xs, ys))
        slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        print('empirical exponent (total time): N^%.3f' % slope)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--pos', type=int, help='extract digits at this position')
    ap.add_argument('--const', choices=('pi', 'ln2'), default='pi')
    ap.add_argument('--window', type=int, default=12)
    ap.add_argument('--engine', choices=('v1', 'v2', 'v3', 'v4'), default='v3',
                    help='v3 (default): division-free 2-adic pipeline, N^1.08 measured; '
                         'v1: naive band (reference); v2: accumulating tree; '
                         'v4: sieve-batched word phase (correctness oracle for the C port)')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--big', action='store_true', help='include N=1e5 tests / larger bench')
    ap.add_argument('--huge', action='store_true', help='include N=1e6 test')
    ap.add_argument('--audit', action='store_true')
    ap.add_argument('--bench', action='store_true')
    ap.add_argument('--dump-x5', nargs=2, metavar=('N', 'FILE'),
                    help='write little-endian bytes of 5^(N-1) to FILE (for the C++ driver)')
    ap.add_argument('--v2-check', action='store_true',
                    help='validate the batched accumulating-tree path (extract_v2) against v1')
    ap.add_argument('--v3-check', action='store_true',
                    help='validate the division-free 2-adic pipeline (extract_v3) against v1')
    ap.add_argument('--beatquad', action='store_true',
                    help='benchmark extract_v3 scaling to demonstrate the subquadratic slope')
    args = ap.parse_args()

    if args.v3_check:
        ok = True
        t0 = time.perf_counter()
        for const in ('pi', 'ln2'):
            for N in list(range(1, 41)) + [100, 997, 2500, 10000]:
                a, _ = extract(N, const=const)
                b, _ = extract_v3(N, const=const)
                if a != b:
                    print('MISMATCH %s N=%d v1=%s v3=%s' % (const, N, a, b))
                    ok = False
        print('V3-CHECK %s (gmpy2=%s, %.2fs)'
              % ('PASSED' if ok else 'FAILED', HAVE_GMPY2, time.perf_counter() - t0))
        sys.exit(0 if ok else 1)

    if args.beatquad:
        sizes = [10000, 30000, 100000, 300000, 1000000]
        if args.big:
            sizes += [3000000, 10000000]
        print('v3 pipeline (gmpy2=%s):' % HAVE_GMPY2)
        pts = []
        for N in sizes:
            digs, st = extract_v3(N)
            pts.append((N, st['total']))
            print('N=%-9d total %9.2fs  (pre %.2f, word %.2f, V %.2f)  -> %s'
                  % (N, st['total'], st['pre'], st['word'], st['vphase'], digs),
                  flush=True)
        xs = [math.log(p[0]) for p in pts]
        ys = [math.log(p[1]) for p in pts]
        nn = len(xs)
        sx, sy = sum(xs), sum(ys)
        sxx = sum(x * x for x in xs)
        sxy = sum(x * y for x, y in zip(xs, ys))
        slope = (nn * sxy - sx * sy) / (nn * sxx - sx * sx)
        print('v3 fitted exponent: N^%.3f' % slope)
        print('v1 reference (same machine, recorded): 0.12s @1e4, 6.4s @1e5, '
              '571s @1e6 (N^1.95 tail); C++ v1: 32.3s @1e6, 3421s @1e7 (N^2.03)')
        return

    if args.v2_check:
        ok = True
        t0 = time.perf_counter()
        for const in ('pi', 'ln2'):
            for N in list(range(1, 41)) + [100, 997, 2500, 10000]:
                a, _ = extract(N, const=const)
                b = extract_v2(N, const=const)
                if a != b:
                    print('MISMATCH %s N=%d v1=%s v2=%s' % (const, N, a, b))
                    ok = False
        print('V2-CHECK %s (%d positions x2 constants, %.2fs)'
              % ('PASSED' if ok else 'FAILED', 44, time.perf_counter() - t0))
        sys.exit(0 if ok else 1)

    if args.dump_x5:
        N, path = int(args.dump_x5[0]), args.dump_x5[1]
        t0 = time.perf_counter()
        X = pow(5, N - 1)
        b = X.to_bytes((X.bit_length() + 7) >> 3, 'little')
        with open(path, 'wb') as fh:
            fh.write(b)
        print('wrote %d bytes (%d bits) of 5^(N-1), N=%d, in %.3fs'
              % (len(b), X.bit_length(), N, time.perf_counter() - t0))
        return

    sys.setrecursionlimit(100000)

    if args.selftest:
        sys.exit(0 if run_selftest(big=args.big, huge=args.huge) else 1)
    if args.audit:
        sys.exit(0 if run_audit() else 1)
    if args.bench:
        run_bench(big=args.big)
        return
    if args.pos:
        if args.engine == 'v1':
            digs, st = extract(args.pos, const=args.const, window=args.window)
            tot = st['pre'] + st['loop']
        elif args.engine == 'v2':
            digs = extract_v2(args.pos, const=args.const, window=args.window)
            st, tot = {}, float('nan')
        elif args.engine == 'v4':
            digs, st = extract_v4(args.pos, window=args.window)
            tot = st['total']
        else:
            digs, st = extract_v3(args.pos, const=args.const, window=args.window)
            tot = st['total']
        print('%s digits at position %d..%d:  %s   [engine %s]' %
              (args.const, args.pos, args.pos + args.window - 1, digs, args.engine))
        if st:
            print('time %.3fs  (%s)' %
                  (tot, ', '.join('%s %.3fs' % (k, v) for k, v in st.items()
                                  if isinstance(v, float))))
        return
    ap.print_help()


if __name__ == '__main__':
    main()
