#!/usr/bin/env python3
"""formula_mass.py -- what does a FORMULA have to look like to kill the third log?
(2026-07-30)

Joseph: "other forms of BBP formulas or other formulas altogether to get below
log^3."  This makes that question quantitative instead of a matter of taste.

merge_mass.py measured the v3 band merge: operand MASS ~ n log^2 n, bit-WORK
~ n log^3 n.  The mass is not mysterious.  A band term is a fraction
2^s / mu with mu a word-size odd number, and the merge tree carries the
DENOMINATOR PRODUCT: at every one of the log T levels the live data is the
whole denominator list, sum_t log2(mu_t) bits.  So

    MASS  ~  D * log T,      D = sum over band terms of log2(mu_t)
    WORK  ~  M(D) * log T  =  D log D log T .

The answer object is Jmax ~ 2.32 N bits.  Therefore

    D = Theta(Jmax)          =>  WORK = Theta(n log^2 n)   PARITY
    D = Theta(Jmax log N)    =>  WORK = Theta(n log^3 n)   where we are.

D = (band term count) x (denominator bits), so the formula controls it through
exactly one number:

    g_eff = Jmax / (number of band terms) = ANSWER BITS PER BAND TERM.

    D / Jmax = log2(mu_max) / g_eff.

Parity needs g_eff >= log2(mu_max) ~ log2 N.  The hex BBP formula as
implemented has g_eff = 1.00 (four denominators per four bits: one band term
per bit of answer), so D/Jmax = log2 N and the third log is exactly that.

This script:
  (1) builds a bit-length accounting model of the real v3 merge and CHECKS it
      against merge_mass.py's instrumented numbers;
  (2) sweeps synthetic BBP-type formulas (base 2^w, r denominators per step,
      so g_eff = w/r) and fits WORK to n log^k n -- testing the prediction that
      g_eff fixed gives log^3 and g_eff ~ log n gives log^2;
  (3) censuses the real formula candidates (hex BBP, Bellard 2^10, and the
      same-denominator collision collapse) at matched N.
"""

import math

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


# ---------------------------------------------------------------------------
# (1) bit-length accounting model of the v3 division-free merge
#     (mirrors decimal_bbp_extract.extract_v3 / merge_mass.band_mass exactly,
#      but tracks bit lengths only, so it scales past what mpz can run)
# ---------------------------------------------------------------------------

def _M(b1, b2):
    s = b1 + b2
    return s * math.log2(s + 2)


def sim_merge(band, Jmax, BLK=1 << 15):
    """band: list of (J, mu_bits) sorted by J ascending. Returns (mass, work)."""
    stats = {'mass': 0, 'work': 0.0}

    def merge(a, b):
        nL, dL, sL = a
        nR, dR, sR = b
        s = sL if sL < sR else sR
        cap = Jmax - s
        stats['mass'] += nL + dR + nR + dL + dL + dR
        stats['work'] += _M(nL, dR) + _M(nR, dL) + _M(dL, dR)
        n_ = max(nL + dR + (sL - s), nR + dL + (sR - s)) + 1
        d_ = dL + dR
        return (min(n_, cap), min(d_, cap), s)

    def reduce_nodes(nodes):
        while len(nodes) > 1:
            nxt = [merge(nodes[i], nodes[i + 1])
                   for i in range(0, len(nodes) - 1, 2)]
            if len(nodes) & 1:
                nxt.append(nodes[-1])
            nodes = nxt
        return nodes[0]

    roots = []
    for lo in range(0, len(band), BLK):
        leaves = [(1, mb, Jmax - J) for J, mb in band[lo:lo + BLK]]
        roots.append(reduce_nodes(leaves))
    reduce_nodes(roots)
    return stats['mass'], stats['work']


# ---------------------------------------------------------------------------
# real hex-BBP band (for validation and for the census)
# ---------------------------------------------------------------------------

def real_bbp_band(N, G=192):
    """Returns (band, Jmax) with band = sorted [(J, mu)] -- raw mu."""
    from decimal_bbp_extract import bbp_pi_terms, cutoff
    K = cutoff(N, G, 'pi')
    LX = (mpz(5) ** (N - 1)).bit_length()
    band = []
    for sign, E, mu in bbp_pi_terms(N, K):
        if E < 0:
            J = -E
            if J < LX + G + 8:
                band.append((J, mu))
    band.sort()
    Jmax = max(J for J, _ in band)
    return band, Jmax


def bits(band):
    return [(J, mu.bit_length()) for J, mu in band]


# ---------------------------------------------------------------------------
# (2) synthetic BBP-type formula: base 2^w, r denominator families per step
#     step k contributes r band terms at J ~ w*k, denominators ~ 8k+odd.
#     g_eff = w / r answer bits per band term.
# ---------------------------------------------------------------------------

def synthetic_band(Jmax, w, r, mscale=8):
    band = []
    K = int(Jmax / w) + 1
    for k in range(K + 1):
        for i in range(r):
            J = w * k + (w * i) // r
            if 0 < J <= Jmax:
                mu = mscale * k + 2 * i + 1
                band.append((J, max(1, mu).bit_length()))
    band.sort()
    return band


def fit(rows, idx, ns_idx=1):
    """rows: tuples; returns dict k -> spread of rows[idx]/(n log^k n)."""
    out = {}
    for k in (1, 2, 3):
        rs = [row[idx] / (row[ns_idx] * math.log(row[ns_idx]) ** k)
              for row in rows]
        out[k] = max(rs) / min(rs)
    return out


def best_fit(rows, idx, ns_idx=1):
    f = fit(rows, idx, ns_idx)
    return min(f, key=f.get), f


# ---------------------------------------------------------------------------

def part1_validate():
    print("=" * 74)
    print("(1) accounting model vs the instrumented real merge (merge_mass.py)")
    print("=" * 74)
    real = {  # from merge_mass.py, gmpy2 instrumented
        1000:    (2519, 688968, 5778886),
        3000:    (7163, 2586830, 23926329),
        10000:   (23415, 11133370, 113613071),
        30000:   (69855, 41820041, 464147308),
        100000:  (232390, 173396541, 2086555301),
    }
    print("%-8s %9s %14s %8s %16s %8s" %
          ("N", "Jmax", "mass(sim)", "vs real", "work(sim)", "vs real"))
    for N in sorted(real):
        band, Jmax = real_bbp_band(N)
        mass, work = sim_merge(bits(band), Jmax)
        rJ, rmass, rwork = real[N]
        print("%-8d %9d %14d %7.3fx %16.0f %7.3fx" %
              (N, Jmax, mass, mass / rmass, work, work / rwork))
    print("  (ratios ~1 => the bit-length model reproduces the real merge)")


def part2_sweep():
    print()
    print("=" * 74)
    print("(2) synthetic formula sweep: does g_eff = answer-bits-per-band-term")
    print("    control the exponent?   band = base 2^w, r denominators per step")
    print("=" * 74)
    JM = [(10 ** e, int(2.3219 * 10 ** e)) for e in (3, 4, 5, 6)]

    def run(label, wr_of_Jmax):
        rows = []
        for N, Jmax in JM:
            w, r = wr_of_Jmax(Jmax)
            band = synthetic_band(Jmax, w, r)
            mass, work = sim_merge(band, Jmax)
            D = sum(mb for _, mb in band)
            rows.append((N, Jmax, mass, work, len(band), D, w / r))
        kM, fM = best_fit(rows, 2)
        kW, fW = best_fit(rows, 3)
        last = rows[-1]
        print("  %-26s g_eff=%6.2f  D/Jmax=%6.2f   MASS~n log^%d  WORK~n log^%d"
              % (label, last[6], last[5] / last[1], kM, kW))
        print("      %-22s work spreads  vs n log n %.2f | n log^2 n %.2f | "
              "n log^3 n %.2f" % ("", fW[1], fW[2], fW[3]))
        print("      %-22s work at N=1e6: %.4e   (slope loglog %.3f)"
              % ("", last[3],
                 math.log(rows[-1][3] / rows[0][3]) /
                 math.log(rows[-1][1] / rows[0][1])))
        return rows

    print("  fixed formulas (g_eff constant):")
    run("hex BBP shape  w=4,r=4", lambda J: (4, 4))
    run("Bellard shape  w=10,r=7", lambda J: (10, 7))
    run("hypothetical   w=32,r=4", lambda J: (32, 4))
    print("  N-dependent formulas (g_eff grows with N):")
    run("g_eff = log2(Jmax)", lambda J: (int(math.log2(J)) * 4, 4))
    run("g_eff = 3*log2(Jmax)", lambda J: (int(3 * math.log2(J)) * 4, 4))
    print()
    print("  PREDICTION: g_eff fixed -> WORK ~ n log^3 n (the wall).")
    print("              g_eff ~ log n -> WORK ~ n log^2 n (AGM parity).")


def part3_census():
    print()
    print("=" * 74)
    print("(3) census of REAL formula candidates for the pipeline")
    print("=" * 74)
    print("%-8s %-22s %10s %10s %8s %8s %12s" %
          ("N", "formula", "bandterms", "distinct", "g_eff", "D/Jmax", "work"))
    for N in (10 ** 4, 10 ** 5, 10 ** 6):
        band, Jmax = real_bbp_band(N)
        from decimal_bbp_extract import bbp_pi_terms, cutoff
        K = cutoff(N, 192, 'pi')
        LX = (mpz(5) ** (N - 1)).bit_length()
        mus = [mu for sign, E, mu in bbp_pi_terms(N, K)
               if E < 0 and -E < LX + 192 + 8]
        D = sum(math.log2(mu) for mu in mus)
        mass, work = sim_merge(bits(band), Jmax)
        print("%-8d %-22s %10d %10d %8.2f %8.2f %12.3e" %
              (N, "hex BBP (v3, now)", len(band), len(set(mus)),
               Jmax / len(band), D / Jmax, work))

        # same-denominator collapse: terms sharing mu fuse into one leaf
        first = {}
        coll = []
        for J, mu in band:
            if mu in first:
                continue
            first[mu] = J
            coll.append((J, mu.bit_length()))
        coll.sort()
        Dc = sum(math.log2(mu) for mu in set(mus))
        massc, workc = sim_merge(coll, Jmax)
        print("%-8s %-22s %10d %10d %8.2f %8.2f %12.3e  (%.2fx work)" %
              ("", "  + collapse dup mu", len(coll), len(coll),
               Jmax / len(coll), Dc / Jmax, workc, workc / work))

        # Bellard base-2^10: 7 denominators per 10 bits
        bel = synthetic_band(Jmax, 10, 7)
        Db = sum(mb for _, mb in bel)
        massb, workb = sim_merge(bel, Jmax)
        print("%-8s %-22s %10d %10d %8.2f %8.2f %12.3e  (%.2fx work)" %
              ("", "Bellard 2^10 shape", len(bel), len(bel),
               Jmax / len(bel), Db / Jmax, workb, workb / work))
        print()


def part4_computation_side():
    """The SAME criterion explains the known complexities of pi/e computation.

    Binary splitting of a term list costs M(D) * log T with D = total
    ratio-denominator bits.  D = Theta(n) gives n log^2 n, D = Theta(n log n)
    gives n log^3 n.  This is why e is n log^2 (Brent 1976) while every
    hypergeometric pi series -- Chudnovsky included -- is n log^3, and why AGM
    (no term list at all) is the only n log^2 route to pi.
    """
    print()
    print("=" * 74)
    print("(4) the same criterion on the COMPUTATION side (sanity check")
    print("    against textbook complexities)")
    print("=" * 74)
    print("%-28s %10s %12s %10s %8s" %
          ("series (n bits of output)", "terms T", "D bits", "D/n", "=> cost"))
    for n in (10 ** 6, 10 ** 9, 10 ** 12):
        print("  n = %.0e" % n)
        # e = sum 1/k! : T log2 T = n, ratio denominator = k
        T = 2
        while T * math.log2(T) < n:
            T = int(T * 1.05) + 1
        D = sum(math.log2(k) for k in range(2, 2000)) if T < 2000 else \
            T * math.log2(T) - T * 1.4427
        print("    %-26s %10.3e %12.3e %10.2f %8s"
              % ("e = sum 1/k! (factorial)", T, D, D / n, "n log^2"))
        # Chudnovsky: 47.11 bits/term, ratio denominator ~ k^3 * 640320^3/24
        T = n / 47.11
        D = T * (3 * math.log2(T) + math.log2(640320 ** 3 / 24))
        print("    %-26s %10.3e %12.3e %10.2f %8s"
              % ("Chudnovsky (geometric)", T, D, D / n, "n log^3"))
        # our extraction band: 1 band term per answer bit, mu ~ 2.86 * Jmax
        Jmax = 2.3219 * n
        T = Jmax
        D = T * math.log2(2.86 * Jmax)
        print("    %-26s %10.3e %12.3e %10.2f %8s"
              % ("hex BBP band (v3)", T, D, D / n, "n log^3"))
        print("    %-26s %10s %12s %10s %8s"
              % ("AGM (no term list)", "-", "-", "-", "n log^2"))


def main():
    part1_validate()
    part2_sweep()
    part3_census()
    part4_computation_side()
    print("=" * 74)
    print("READING")
    print("  The third log is D/Jmax = log2(mu_max)/g_eff, and g_eff is a")
    print("  property of the FORMULA alone: answer bits per band term.")
    print("  Every fixed BBP-type formula has g_eff = O(1), hence one extra")
    print("  log.  Parity requires g_eff ~ log2 N, i.e. a formula whose base")
    print("  2^w GROWS with N while the denominators stay word-size.")


if __name__ == "__main__":
    main()
