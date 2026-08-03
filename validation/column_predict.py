#!/usr/bin/env python3
"""column_predict.py -- can the column stream be PREDICTED? (2026-08-03,
Joseph's question: can the fractions be predicted or controlled?)

"Controlled" is the formula axis (g_eff, priced, ceiling 1.4286, swept
through degree 3) and the telescoping route is blocked by transcendence
itself.  "Predicted" has one sense never measured: every prior test of the
column stream B(d) = #{t : 2^d mod mu_t in top half} was LINEAR or spectral
(Berlekamp-Massey, mod-q complexity, DFT, SVD, zlib).  A NONLINEAR learned
predictor -- Markov/context models of increasing order, with and without
phase features -- was never run.  If the stream has finite-context
structure invisible to linear tests, this sees it; an accuracy wedge above
noise would be compressible structure = a lever on Summit 1.

Streams tested (train first half, test second half):
  wild:   b(d) = B(d) mod 2 and c(d) = [B(d) > T/2] for the real family
          mu_t = 8t+1, t <= T
  tame:   same stream for the family {odd divisors > 1 of 2^12 - 1}
          (orbit periods divide 12 => the stream is 12-periodic; any
          context model with k >= 12 must approach 100% -- probe control)
  random: seeded coin flips (must sit at 50%)

Predictors: pure k-gram context models, k in {4, 8, 12, 16}, plus an
order-8 model with phase feature d mod 8.  Unseen context -> global
majority.  Expected on the wild stream: 50% everywhere (the 15th
independent confirmation, this time against nonlinear prediction); any
value many sigma above 50% = WIRE TRIP.

RESULT + DIAGNOSIS (2026-08-03, same session).  Parity: coin-flip at every
order (0.498-0.515, <= 1.6 sigma) -- the exact content is unpredictable
also nonlinearly.  Median stream tripped at ~0.75 (+28 sigma); diagnosed
BENIGN in minutes: base rate P(c=1) = 0.256 (window sign bias -24.6 of
sd 57.7 -- the frozen-phase/finite-window drift of long orbits; "always
guess low" alone gives 0.744) plus lag-1 deviation autocorrelation 0.689
decaying to noise by lag ~16 ~ log T -- the SAME locality constant
kakeya_hop measured for carries, now seen from a second independent probe.
No exactness lever: +-1 column errors are fatal and sign bias gives none.
"""

import random

import numpy as np


def column_bits(mus, D):
    mus = np.array(mus, dtype=np.int64)
    r = np.ones(len(mus), dtype=np.int64)
    b = np.empty(D, dtype=np.uint8)
    c = np.empty(D, dtype=np.uint8)
    T = len(mus)
    for d in range(D):
        r = (r << 1) % mus                    # r = 2^{d+1} mod mu
        cnt = int((2 * r > mus).sum())
        b[d] = cnt & 1
        c[d] = 1 if cnt > T // 2 else 0
    return b, c


def kgram_accuracy(bits, k, phase=0):
    n = len(bits)
    half = n // 2
    table = {}
    ctx = 0
    mask = (1 << k) - 1
    for i in range(k, half):
        ctx = 0
        for j in range(k):
            ctx = (ctx << 1) | int(bits[i - k + j])
        key = (ctx, i % phase) if phase else ctx
        s, t = table.get(key, (0, 0))
        table[key] = (s + int(bits[i]), t + 1)
    glob = int(bits[:half].sum()) * 2 >= half
    hit = tot = 0
    for i in range(half + k, n):
        ctx = 0
        for j in range(k):
            ctx = (ctx << 1) | int(bits[i - k + j])
        key = (ctx, i % phase) if phase else ctx
        s, t = table.get(key, (0, 0))
        pred = (s * 2 >= t) if t else glob
        hit += int(pred == bool(bits[i]))
        tot += 1
    return hit / tot, tot


def battery(name, bits):
    row = [name]
    for k, ph in ((4, 0), (8, 0), (12, 0), (16, 0), (8, 8)):
        acc, tot = kgram_accuracy(bits, k, ph)
        sigma = (acc - 0.5) / (0.5 / tot ** 0.5)
        row.append("%.3f(%+.1fs)" % (acc, sigma))
    print("  %-12s %s" % (row[0], "  ".join(row[1:])))
    return row


def main():
    T, D = 2048, 6144
    print("nonlinear prediction of the column stream  T=%d D=%d" % (T, D))
    print("  %-12s %-13s %-13s %-13s %-13s %s"
          % ("stream", "k=4", "k=8", "k=12", "k=16", "k=8+phase"))

    b, c = column_bits([8 * t + 1 for t in range(1, T + 1)], D)
    battery("wild parity", b)
    battery("wild median", c)

    tame_mus = [m for m in range(3, 4096, 2) if 4095 % m == 0]
    tb, tc = column_bits(tame_mus, D)
    battery("TAME parity", tb)

    rng = random.Random(23)
    rb = np.array([rng.getrandbits(1) for _ in range(D)], dtype=np.uint8)
    battery("random ctrl", rb)

    print()
    print("expect: tame -> ~1.000 at k>=12 (probe validated); random ~0.500;")
    print("wild ~0.500 everywhere = the fractions' aggregate is unpredictable")
    print("also in the learned/nonlinear sense.  >> 0.5 by many sigma = WIRE.")


if __name__ == '__main__':
    main()
