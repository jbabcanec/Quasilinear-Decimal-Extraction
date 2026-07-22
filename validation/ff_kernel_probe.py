#!/usr/bin/env python3
"""
ff_kernel_probe.py -- 2-kernel growth of the weighted function-field mirror
F = sum_v x^{4 int(v)} / (1 + x^3 pol(v))  over F_2[[x]]  (the true V-twin;
companion to ff_mirror.py Part 2).

Christol: F algebraic over F_2(x)  <=>  its coefficient sequence is 2-automatic
<=>  the 2-kernel { (a_{2^k m + r})_m } is FINITE.  We generate a_m for m < P
via LFSR streams (each 1/(1+x^3 v) is a linear recurring sequence; no big-int
work), then count distinct kernel rows per level.

RESULT (recorded run, P = 2^13, 2026-07-22):
    F (weighted mirror):   k=0..8 -> 1,2,4,8,16,32,64,128,256 distinct rows
                           == the MAXIMUM 2^k at every level (fully saturated)
    control 1/(1+x+x^3):   saturates at 7 rows from k=3 on (algebraic, as must)
DEEPER RUN (P = 2^15, same day): saturation continues -- 2^k exactly through
    k=9 (512/512); k=10 gives 1013/1024, the shortfall being a compare-length
    artifact (16-byte prefixes), not structure.
Maximal kernel growth through the testable range = strong evidence F is
TRANSCENDENTAL over F_2(x), consistent with the empty Mahler/Frobenius-relation
hunt in ff_mirror.py.  Reading: the index-geometric weight (the positional map
int: F_2[x] -> N in the exponent) destroys algebraicity even carry-free, while
degree-graded weights stay closed-form.  V carries two independent obstructions:
carries (values) and index-grading (weights).
Open exercise (candidate THEOREM, attack queue): prove transcendence by
exhibiting infinitely many distinct kernel rows from the incommensurate LFSR
periods of the streams 1/(1+x^3 v).
"""

P = 1 << 13


def stream_inv(taps_int, deg, length):
    """Coefficient stream of 1/q(x): s_e = XOR_{i>=1} q_i s_{e-i}, s_0 = 1."""
    s = bytearray(length)
    s[0] = 1
    hist = 1
    tapmask = taps_int >> 1
    for e in range(1, length):
        b = bin(hist & tapmask).count("1") & 1
        s[e] = b
        hist = ((hist << 1) | b) & ((1 << deg) - 1) if deg else 0
    return s


def build_F(P):
    a = bytearray(P)
    v = 0
    while 4 * v < P:
        q = 1 | (v << 3)
        s = stream_inv(q, q.bit_length() - 1, P - 4 * v)
        base = 4 * v
        for e in range(P - 4 * v):
            a[base + e] ^= s[e]
        v += 1
    return a


def kernel_growth(a, P, kmax=8, label=""):
    print(f"  {label}: kernel sizes by level k:")
    sizes = []
    for k in range(kmax + 1):
        L = P >> k
        if L < 32:
            break
        rows = set()
        for r in range(1 << k):
            row = bytes(a[(m << k) + r] for m in range(min(L, (P - r) >> k)))
            rows.add(row[: L // 2])
        sizes.append(len(rows))
        print(f"    k={k}: {len(rows)} distinct rows (max {1 << k})")
    return sizes


if __name__ == "__main__":
    a = build_F(P)
    print(f"F built to precision x^{P}; density {sum(a)}/{P}")
    sF = kernel_growth(a, P, 8, "F (weighted mirror)")
    ctrl = bytearray(stream_inv(0b1011, 3, P))
    sC = kernel_growth(ctrl, P, 8, "control 1/(1+x+x^3) [algebraic]")
    saturated = all(sF[k] == 1 << k for k in range(len(sF)))
    print(f"\nF kernel fully saturated (2^k at every level): {saturated}")
    print(f"control saturates finitely: {max(sC) < 16}")
    print("VERDICT:", "maximal kernel growth -- evidence-transcendental"
          if saturated else "kernel growth sub-maximal -- investigate")
