#!/usr/bin/env python3
"""fence2_wide.py -- FENCE 2 WIDENED (2026-08-03): the incomplete-reciprocity
dual search (`l2_swap_hunt.py`: ONE point (A,T)=(483,56), P=240, linear +
quadratic zeta8-duals, single precision) widened on three axes:
  (i)   EIGHT parameter points (A,T) instead of one/two,
  (ii)  a CUBIC dual zeta8^{k^3} the original never searched,
  (iii) PERSISTENCE: every candidate relation found at W1=224 is re-verified
        exactly at W2=448 (recomputed values; artefacts cannot survive).
Object: W(A,T) = sum_{s<T} 16^s/(A-8s) mod 2^P -- the T-term slope-8
incomplete sum (true V geometry).  Swap thesis: W relates at low height to
8-term slope-T duals over Q_2(zeta8) + complete logs + prefactor 16^T.
Any PERSISTENT relation = the 2-adic incomplete Landsberg-Schaar = a crack
shape.  Expected: planted control recovered; all candidates die at W2.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from postnikov_compose import lll

PTS = [(389,), (483,), (517,), (611,), (733,), (807,), (901,), (1013,)]


def values(A, T, P):
    MOD = 1 << P
    inv = lambda o: pow(o, -1, MOD)
    Wv, p16 = 0, 1
    for s in range(T):
        Wv = (Wv + p16 * inv(A - 8 * s)) % MOD
        p16 = p16 * 16 % MOD

    def dual(efun):
        comp = [0, 0, 0, 0]
        for k in range(8):
            e = efun(k) % 8
            c, sg = e % 4, (-1) ** (e // 4)
            comp[c] = (comp[c] + sg * inv(A - T * k)) % MOD
        return comp

    D1, D2, D3 = dual(lambda k: k), dual(lambda k: k * k), dual(lambda k: k ** 3)
    LG, w, wp = 0, 8 * T * inv(A) % MOD, 1
    for n in range(1, P // 2 + 8):
        wp = wp * w % MOD
        e = (n & -n).bit_length() - 1
        LG = (LG - (wp >> e) * inv(n >> e)) % MOD
    return Wv, D1, D2, D3, LG, pow(16, T, MOD)


def hunt(names, vals, W):
    n = len(vals)
    S = 1 << W
    rows = []
    for i, v in enumerate(vals):
        r = [0] * n
        r[i] = 1
        r.append(int(S * (v % S)))
        rows.append(r)
    rows.append([0] * n + [S * S])
    best = None
    for r in lll(rows):
        c = r[:n]
        if any(c) and r[-1] == 0:
            h = max(abs(x) for x in c)
            if best is None or h < best[1]:
                best = (c, h)
    return best


def main():
    W1, W2 = 224, 448
    trips = 0
    print("FENCE 2 widened: 8 points x {lin,quad,cubic} duals x persistence")
    for (A,) in PTS:
        T = ((A // 8) - 2) & ~1                      # even: A - Tk stays odd
        Wv, D1, D2, D3, LG, S16 = values(A, T, W1 + 16)
        Wh, E1, E2, E3, LGh, S16h = values(A, T, W2 + 16)
        for tag, Dl, Dh in (("lin+quad", D1 + D2, E1 + E2),
                            ("quad+cubic", D2 + D3, E2 + E3)):
            # drop structural zeros (zeta8^{k^j} support misses classes) and
            # dedupe components shared between duals -- both would otherwise
            # gift height-1 "relations" that mean nothing
            seen, keep = set(), []
            for i, (v, vh) in enumerate(zip(Dl, Dh)):
                if v % (1 << W1) != 0 and v not in seen:
                    seen.add(v)
                    keep.append((i, v, vh))
            names = ["W"] + ["d%d" % i for i, _, _ in keep] + ["LG", "1"]
            vals = [Wv] + [v for _, v, _ in keep] + [LG, 1]
            hi = [Wh] + [vh for _, _, vh in keep] + [LGh, 1]
            if 4 * T < W1:                       # else 16^T = 0 mod 2^W1
                names.insert(-1, "16^T")
                vals.insert(-1, S16)
                hi.insert(-1, S16h)
            got = hunt(names, vals, W1)
            if got is None:
                print("  A=%-5d T=%-4d %-10s no relation" % (A, T, tag))
                continue
            c, h = got
            resid = sum(ci * vi for ci, vi in zip(c, hi)) % (1 << W2)
            persists = resid == 0
            # only W-involving persistent relations are reciprocity-shaped;
            # cW = 0 persistents are point-specific unit-fraction identities
            # among the dual denominators (measured: A=389)
            if persists and c[0] != 0:
                trips += 1
                verdict = "PERSISTS with cW!=0 ** WIRE TRIP **"
            elif persists:
                verdict = "persists, cW=0 -> point coincidence"
            else:
                verdict = "-> ARTEFACT"
            print("  A=%-5d T=%-4d %-10s height 2^%-4d floor 2^%-3.0f %s"
                  % (A, T, tag, h.bit_length() - 1, W1 / 12.0, verdict))
    # planted control at the first point
    A, T = PTS[0][0], PTS[0][0] // 8 - 2
    Wv, D1, D2, D3, LG, S16 = values(A, T, W1 + 16)
    plant = (3 * Wv - 7 * D1[0] + 2 * LG) % (1 << W1)
    got = hunt(["W", "d0", "LG", "plant"], [Wv, D1[0], LG, plant], W1)
    ok = got and got[1] <= 8
    print("  control planted (3,-7,2): %s"
          % ("RECOVERED h=%d" % got[1] if ok else "*** MISSED ***"))
    print()
    if trips:
        print("*** %d persistent relation(s): investigate before believing ***"
              % trips)
    else:
        print("VERDICT: fence 2 widened (8 points, cubic dual, persistence)")
        print("and still EMPTY: no 2-adic incomplete Landsberg-Schaar at")
        print("height < floor across the widened window.")


if __name__ == '__main__':
    main()
