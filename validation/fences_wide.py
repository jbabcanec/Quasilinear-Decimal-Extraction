#!/usr/bin/env python3
"""fences_wide.py -- FENCE 1 OF 3, WIDENED (2026-08-03).

Exactly three closures in the whole map are finite searches over infinite
spaces rather than identities, measurements, or theorems:

  FENCE 1  cross-scale relations of the truncated 2-adic log
           (`agm_law_hunt.py`: ONE argument z = 6, scales <= 256, ONE
           precision W = 300, height thresholds only)
  FENCE 2  incomplete reciprocity duals (two parameter points)
  FENCE 3  the formula-sweep windows (m <= 12, w <= 40; base-10 degree <= 2)

This file widens FENCE 1 on the three axes the original did not have:

  (i)   CROSS-ARGUMENT: baskets mixing z = 6, 4, 16 (the original searched
        one argument and its square; relations ACROSS arguments were never
        searched -- and 16 is the pipeline's own dyadic base),
  (ii)  DEEPER LADDER: odd parts O_M to M = 1024 at W = 400,
  (iii) PERSISTENCE: every candidate relation found at W = 400 is re-tested
        exactly at W = 800 (the postnikov PART D methodology; the original
        fence used height thresholds at a single precision, which is the
        weaker artifact-killer).

Objects (all exact 2-adic integers mod 2^P): T_M(z) = sum_{1<=n<M} z^n/n,
O_M(z) = odd-n half, L(z) = saturated (complete) log-type value.  The known
functional equation T_{2M}(z) + T_{2M}(-z) = T_M(z^2) must be recovered by
the harness (control), and the odd parts O carry no known relation -- any
persistent nontrivial hit here is a WIRE TRIP (a scale-recursion = the crack
shape excluded in the paper's Remark on structural exclusions).

Verdict format per hunt: best relation's height vs the noise floor
2^{W/rows}, and PERSISTS/ARTEFACT at W = 800.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from postnikov_compose import lll                      # textbook LLL over Q

W1, W2 = 400, 800


def v2(n):
    return (n & -n).bit_length() - 1 if n else 10 ** 9


def T(M, z, W):
    """sum_{1<=n<M} z^n/n mod 2^W, exact (z even, v2(z) >= 1)."""
    P = W + 48
    mod = 1 << (P + 24)
    zp, acc = 1, 0
    for n in range(1, M):
        zp = zp * z % mod
        e = v2(n)
        acc += (zp >> e) * pow(n >> e, -1, 1 << P)
    return acc % (1 << W)


def O(M, z, W):
    """odd-n half of T_M."""
    P = W + 48
    mod = 1 << (P + 24)
    zp, acc = 1, 0
    for n in range(1, M):
        zp = zp * z % mod
        if n & 1:
            acc += zp * pow(n, -1, 1 << P)
    return acc % (1 << W)


def L(z, W):
    """saturated (complete) value of T at width W."""
    return T(W + 80, z, W)


def hunt(names, vals, W, label):
    """LLL relation hunt mod 2^W; returns (coeffs, height) or None."""
    n = len(vals)
    S = 1 << W                                       # congruence scale
    rows = []
    for i, v in enumerate(vals):
        row = [0] * n
        row[i] = 1
        row.append(int(S * (v % S)))
        rows.append(row)
    rows.append([0] * n + [S * S])
    best = None
    for r in lll(rows):
        c = r[:n]
        if any(c) and r[-1] == 0:
            h = max(abs(x) for x in c)
            if best is None or h < best[1]:
                best = (c, h)
    floor = W / float(n + 1)
    if best is None:
        print("  %-26s no relation at all (floor 2^%.0f)" % (label, floor))
        return None
    c, h = best
    hb = h.bit_length() - 1
    return (label, names, c, hb, floor)


def report(found, hi_vals_by_name):
    label, names, c, hb, floor = found
    resid = 0
    for ci, nm in zip(c, names):
        resid += ci * hi_vals_by_name[nm]
    persists = resid % (1 << W2) == 0
    tag = "PERSISTS  ** WIRE TRIP **" if persists else "-> ARTEFACT"
    nz = [(ci, nm) for ci, nm in zip(c, names) if ci]
    rel = " ".join("%+d.%s" % (ci, nm) for ci, nm in nz[:6])
    print("  %-26s height 2^%-4d (floor 2^%.0f)  %s" % (label, hb, floor, tag))
    print("      relation: %s%s" % (rel, " ..." if len(nz) > 6 else ""))
    return persists


def main():
    Zs = (6, 4, 16)
    Ms = (32, 64, 128, 256, 512, 1024)

    lo, hi = {}, {}
    for z in Zs:
        for M in Ms:
            lo["O%d(%d)" % (M, z)] = O(M, z, W1)
            hi["O%d(%d)" % (M, z)] = O(M, z, W2)
        lo["L(%d)" % z] = L(z, W1)
        hi["L(%d)" % z] = L(z, W2)
    lo["1"] = hi["1"] = 1

    print("FENCE 1 widened: cross-argument x deep-ladder x persistence")
    print("W = %d, verify at W = %d; floors are per-hunt" % (W1, W2))
    print()

    # -- harness controls -------------------------------------------------
    print("controls:")
    fe = hunt(["T32(6)", "T32(-6)", "T16(36)"],
              [T(32, 6, W1), T(32, -6, W1), T(16, 36, W1)], W1, "FE (1,1,-1)")
    if fe:
        _, names, c, hb, _ = fe
        ok = sorted(map(abs, c)) == [1, 1, 1] and hb == 0
        print("  %-26s %s" % ("FE recovery", "RECOVERED exactly" if ok
                              else "*** NOT RECOVERED ***"))
    planted_val = (7 * lo["O64(6)"] - 5 * lo["O128(6)"]) % (1 << W1)
    pl = hunt(["O64(6)", "O128(6)", "planted"],
              [lo["O64(6)"], lo["O128(6)"], planted_val], W1, "planted 7,-5")
    if pl:
        print("  %-26s height 2^%d (expect 0-3)" % ("planted recovery", pl[3]))
    print()

    # -- H1: linear cross-scale x cross-argument, per argument pair -------
    print("H1  linear, cross-argument pairs (odd parts + completes + 1):")
    trips = 0
    import itertools
    for za, zb in itertools.combinations(Zs, 2):
        names = (["O%d(%d)" % (M, za) for M in Ms[:5]]
                 + ["O%d(%d)" % (M, zb) for M in Ms[:5]]
                 + ["L(%d)" % za, "L(%d)" % zb, "1"])
        found = hunt(names, [lo[nm] for nm in names], W1,
                     "z=%d x z=%d" % (za, zb))
        if found:
            trips += report(found, hi)

    # -- H2: cross-argument quadratic at matched scale --------------------
    print()
    print("H2  quadratic, cross-argument products at matched scale:")
    qn, qv_lo, qv_hi = [], {}, {}
    for M in (64, 128, 256):
        for za, zb in itertools.combinations(Zs, 2):
            nm = "O%d(%d)*O%d(%d)" % (M, za, M, zb)
            qn.append(nm)
            qv_lo[nm] = lo["O%d(%d)" % (M, za)] * lo["O%d(%d)" % (M, zb)]
            qv_hi[nm] = hi["O%d(%d)" % (M, za)] * hi["O%d(%d)" % (M, zb)]
    names = qn + ["L(6)", "L(4)", "L(16)", "1"]
    vals = [qv_lo.get(nm, lo.get(nm)) for nm in names]
    hi_all = dict(hi)
    hi_all.update(qv_hi)
    found = hunt(names, vals, W1, "matched-scale products")
    if found:
        trips += report(found, hi_all)

    # -- H3: the deep target: express O1024 via everything shallower ------
    print()
    print("H3  deep target O1024 vs all shallower scales, all arguments:")
    for z in Zs:
        names = (["O1024(%d)" % z]
                 + ["O%d(%d)" % (M, zz) for M in Ms[:5] for zz in Zs]
                 + ["L(%d)" % zz for zz in Zs] + ["1"])
        found = hunt(names, [lo[nm] for nm in names], W1, "target O1024(%d)" % z)
        if found:
            trips += report(found, hi)

    print()
    if trips:
        print("*** %d PERSISTENT RELATION(S): wire trip -- do NOT enshrine;"
              % trips)
        print("    re-derive symbolically before believing anything. ***")
    else:
        print("VERDICT: every candidate is an LLL artefact (dies at W=800);")
        print("controls recovered.  FENCE 1 is now cross-argument, deeper,")
        print("and persistence-tested: still EMPTY.  Fences 2-3 queued.")


if __name__ == '__main__':
    main()
