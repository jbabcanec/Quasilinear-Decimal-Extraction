#!/usr/bin/env python3
"""merge_mass.py -- is the third log REAL or a loose analysis? (2026-07-28)

Joseph: "log^3 is arbitrary, push harder." Fair. The provable bound on the v3
merge is O(N log^3 N), but v3 MEASURES N^1.079 wall-clock -- more consistent
with N log^2 than N log^3. This instruments the ACTUAL v3 division-free merge
and counts the exact total operand bit-mass fed to the multiplier, then checks
which curve n log^k n it tracks (k=1,2,3). If it tracks n log^2, the third log
is an artifact of loose analysis and the extractor is EMPIRICALLY at parity
with AGM -- a major honest reframe. If n log^3, the bound is tight.
"""

import math
from decimal_bbp_extract import bbp_pi_terms, cutoff
try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def band_mass(N, G=192):
    K = cutoff(N, G, 'pi')
    terms = bbp_pi_terms(N, K)
    n5 = N - 1
    LX = (mpz(5) ** n5).bit_length()
    band = []
    for sign, E, mu in terms:
        if E < 0:
            J = -E
            if J < LX + G + 8:
                band.append((sign, J, mu))
    if not band:
        return None
    band.sort(key=lambda t: t[1])
    Jmax = max(J for _, J, _ in band)
    stats = {'mbits': 0, 'nmul': 0, 'work': 0.0}

    def _M(b1, b2):
        s = b1 + b2
        return s * math.log2(s + 2)          # quasilinear-multiply proxy

    def merge(a, b, depth=0):
        nL, dL, sL = a
        nR, dR, sR = b
        s = sL if sL < sR else sR
        cap = Jmax - s
        m = (nL.bit_length() + dR.bit_length()
             + nR.bit_length() + dL.bit_length()
             + dL.bit_length() + dR.bit_length())
        stats['mbits'] += m
        stats['work'] += (_M(nL.bit_length(), dR.bit_length())
                          + _M(nR.bit_length(), dL.bit_length())
                          + _M(dL.bit_length(), dR.bit_length()))
        stats['nmul'] += 3
        n_ = ((nL * dR) << (sL - s)) + ((nR * dL) << (sR - s))
        d_ = dL * dR
        if n_.bit_length() > cap:
            n_ &= (mpz(1) << cap) - 1
        if d_.bit_length() > cap:
            d_ &= (mpz(1) << cap) - 1
        return (n_, d_, s)

    def reduce_nodes(nodes):
        while len(nodes) > 1:
            nxt = [merge(nodes[i], nodes[i + 1])
                   for i in range(0, len(nodes) - 1, 2)]
            if len(nodes) & 1:
                nxt.append(nodes[-1])
            nodes = nxt
        return nodes[0]

    BLK = 1 << 15
    roots = []
    for lo in range(0, len(band), BLK):
        leaves = [(mpz(sg), mpz(mu), Jmax - J)
                  for sg, J, mu in band[lo:lo + BLK]]
        roots.append(reduce_nodes(leaves))
    reduce_nodes(roots)
    return len(band), Jmax, stats['mbits'], stats['work']


def main():
    print("Instrumenting the real v3 merge. n = Jmax (answer-object bits).")
    print("MASS = total operand bits; WORK = sum M(b) ~ b log b (bit-ops proxy).")
    print("%-9s %10s %14s | mass fit    %14s | work fit" %
          ("N", "Jmax", "operand_bits", "M(b)-work"))
    rows = []
    for N in [1000, 3000, 10000, 30000, 100000, 300000, 1000000]:
        r = band_mass(N)
        if r is None:
            continue
        Tb, Jmax, mbits, work = r
        rows.append((N, Jmax, mbits, work))
        print("%-9d %10d %14d | %14.0f" % (N, Jmax, mbits, work))
    print()
    print("Which n log^k n does each track? (spread max/min of the ratio; "
          "1.00 = perfect)")
    for label, idx in [("MASS", 2), ("WORK", 3)]:
        print("  %s:" % label)
        for k in (1, 2, 3):
            rs = [row[idx] / (row[1] * math.log(row[1])**k) for row in rows]
            print("     vs n log^%d n : spread %.2f%s"
                  % (k, max(rs) / min(rs),
                     "  <== flattest" if False else ""))
        # find flattest
        best = min((1, 2, 3),
                   key=lambda k: (lambda rs: max(rs) / min(rs))(
                       [row[idx] / (row[1] * math.log(row[1])**k) for row in rows]))
        print("     => best fit: n log^%d n" % best)
    print()
    print("READING: MASS is the data size (~n log^2). WORK is the bit-operations.")
    print("If WORK also tracks n log^2, v3's band IS at N log^2 = AGM parity")
    print("(the third log is a loose analysis). If WORK tracks n log^3, the extra")
    print("log is real and comes from re-transforming (FFT) at every merge level")
    print("-- which points the attack at FFT-reuse / relaxed multiplication.")


if __name__ == "__main__":
    main()
