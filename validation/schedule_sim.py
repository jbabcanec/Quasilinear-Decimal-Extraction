#!/usr/bin/env python3
"""schedule_sim.py -- is the third log a SCHEDULING artifact or INTRINSIC?
(2026-07-28, the push-harder decisive test.)

merge_mass.py measured: the v3 merge does n log^2 DATA but n log^3 WORK.
Joseph's intuition: the extra log is arbitrary. Decisive test: simulate the
merge's bit-WORK (sum M(b), b log b) under many DIFFERENT correct tree
schedules -- faithful size bookkeeping only, no real arithmetic, so we can
try thousands of leaves and many schedules cheaply. If some schedule's work
drops toward n log^2, the log is a removable scheduling artifact (a LEAD). If
every schedule sits at n log^3, the triangle is intrinsic to fraction-merging.

Node model (faithful to extract_v3): a node has (num_bits, den_bits, Jhi)
where Jhi = max J in its subtree; cap = Jhi (v3: cap = Jmax - min_shift = Jhi).
Merge of L, R (Jn = max(Jhi_L, Jhi_R), cap = Jn):
  den = min(den_L + den_R, cap)
  num = min(max(num_L + den_R + (Jn - Jhi_L),
                num_R + den_L + (Jn - Jhi_R)) + 1, cap)
  work += M(num_L, den_R) + M(num_R, den_L) + M(den_L, den_R),  M(x,y)=(x+y)lg(x+y)
"""

import heapq
import math
from decimal_bbp_extract import bbp_pi_terms, cutoff


def real_Js(N, G=192):
    K = cutoff(N, G, 'pi')
    LX = ((N - 1) * math.log2(5))
    Js = []
    for sign, E, mu in bbp_pi_terms(N, K):
        if E < 0:
            J = -E
            if J < LX + G + 8:
                Js.append((J, max(1, int(math.log2(mu)) + 1)))
    return Js


def M(x, y):
    s = x + y
    return s * math.log2(s + 2)


def merge(a, b, work):
    numA, denA, JA = a
    numB, denB, JB = b
    Jn = JA if JA > JB else JB
    cap = Jn
    work[0] += M(numA, denB) + M(numB, denA) + M(denA, denB)
    den = min(denA + denB, cap)
    num = min(max(numA + denB + (Jn - JA), numB + denB + (Jn - JB)) + 1, cap)
    return (num, den, Jn)


def sched_balanced(leaves):
    """v3: leaves already in an order; balanced pairwise reduction."""
    work = [0.0]
    nodes = leaves
    while len(nodes) > 1:
        nxt = [merge(nodes[i], nodes[i + 1], work)
               for i in range(0, len(nodes) - 1, 2)]
        if len(nodes) & 1:
            nxt.append(nodes[-1])
        nodes = nxt
    return work[0]


def sched_greedy_minJ(leaves):
    """precision-greedy: always merge the two nodes of smallest Jhi (keep caps
    small as long as possible)."""
    work = [0.0]
    heap = [(J, i, (num, den, J)) for i, (num, den, J) in enumerate(leaves)]
    heapq.heapify(heap)
    ctr = len(leaves)
    while len(heap) > 1:
        _, _, a = heapq.heappop(heap)
        _, _, b = heapq.heappop(heap)
        n = merge(a, b, work)
        heapq.heappush(heap, (n[2], ctr, n))
        ctr += 1
    return work[0]


def sched_greedy_minsize(leaves):
    """merge the two nodes whose product would be cheapest (smallest den+num)."""
    work = [0.0]
    heap = [(num + den, i, (num, den, J)) for i, (num, den, J) in enumerate(leaves)]
    heapq.heapify(heap)
    ctr = len(leaves)
    while len(heap) > 1:
        _, _, a = heapq.heappop(heap)
        _, _, b = heapq.heappop(heap)
        n = merge(a, b, work)
        heapq.heappush(heap, (n[0] + n[1], ctr, n))
        ctr += 1
    return work[0]


def fit(rows, idx):
    best = None
    for k in (1, 2, 3):
        rs = [row[idx] / (row[0] * math.log(row[0])**k) for row in rows]
        spread = max(rs) / min(rs)
        if best is None or spread < best[1]:
            best = (k, spread)
    return best


def main():
    print("Simulated merge WORK (sum M(b)) under different schedules.")
    print("%-9s %8s %14s %14s %14s" %
          ("N", "T", "balanced(v3)", "greedy-minJ", "greedy-minsize"))
    rowsB, rowsGJ, rowsGS = [], [], []
    for N in [3000, 10000, 30000, 100000, 300000]:
        Js = real_Js(N)
        if not Js:
            continue
        T = len(Js)
        n = max(J for J, _ in Js)
        leavesJ = sorted(Js)                       # by J (v3)
        leaves = [(1, d, J) for J, d in leavesJ]
        wB = sched_balanced(leaves)
        wGJ = sched_greedy_minJ(leaves)
        wGS = sched_greedy_minsize(leaves)
        rowsB.append((n, wB)); rowsGJ.append((n, wGJ)); rowsGS.append((n, wGS))
        print("%-9d %8d %14.3e %14.3e %14.3e" % (N, T, wB, wGJ, wGS))
    print()
    for name, rows in [("balanced(v3)", rowsB), ("greedy-minJ", rowsGJ),
                       ("greedy-minsize", rowsGS)]:
        k, spread = fit(rows, 1)
        print("  %-16s best fit n log^%d n (spread %.2f)" % (name, k, spread))
    print()
    print("READING: if every schedule fits n log^3, the third log is INTRINSIC")
    print("to balanced fraction-merging (the precision triangle). If a schedule")
    print("reaches n log^2, the log is a removable scheduling artifact = a LEAD.")


if __name__ == "__main__":
    main()
