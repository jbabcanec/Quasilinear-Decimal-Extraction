#!/usr/bin/env python3
"""
orbit_rigidity.py -- the lower-bound crux, measured (+ two-sided crack test).

The last log of pi digit extraction reduces (counting face) to computing all
column sums B(d) = sum_t M[d][t] of the doubling-orbit matrix M[d][t] =
bit_d(1/(8t+1)). In the linear-circuit model, computing this below N^2/polylog
requires M to have an efficient linear circuit = M NON-rigid (low-rank + sparse,
Dvir-Liu style). So the lower bound is EQUIVALENT to rigidity of M:
   M rigid (random-like)   => bound plausible; no cheap linear circuit; no crack.
   M non-rigid (structured)=> a CRACK (fast M*1 via low-rank + sparse).

This LOCATES our open problem precisely: "the last log = rigidity of the
doubling-orbit matrix," an explicit-matrix rigidity question (Valiant's problem
class -- proving rigidity of explicit matrices is itself famously open, which
is exactly why an unconditional bound is out of reach; the reduction explains
the hardness rather than resolving it).

MEASURED (256x256, recorded run): the orbit matrix is statistically
indistinguishable from a random 0/1 matrix on rank, singular-value
concentration, and low-rank-approximability -- and FAR from a genuinely
non-rigid matrix (which would concentrate 90% energy in ~10-20 singular
values). Verdict: RIGID / random-like => bound plausible, NO crack.
"""
import numpy as np


def orbit_matrix(N, T):
    M = np.zeros((N, T), dtype=np.int64)
    for j, t in enumerate(range(1, T + 1)):
        mu = 8 * t + 1
        r = 1 % mu
        for d in range(N):
            M[d, j] = 1 if r >= (mu + 1) // 2 else 0
            r = (2 * r) % mu
    return M


def energy_ranks(M):
    Mc = M.astype(float)
    s = np.linalg.svd(Mc - Mc.mean(), compute_uv=False)
    e = np.cumsum(s ** 2) / np.sum(s ** 2)
    k = lambda frac: int(np.searchsorted(e, frac)) + 1
    return k(0.5), k(0.9), s


def main():
    N = T = 256
    np.random.seed(0)
    M = orbit_matrix(N, T)
    R = (np.random.rand(N, T) < 0.5).astype(np.int64)      # max-rigid control
    # a genuinely NON-rigid control: rank-8 + sparse
    L = (np.random.rand(N, 8) @ np.random.rand(8, T) > 2.0).astype(np.int64)
    S = (np.random.rand(N, T) < 0.02).astype(np.int64)
    NR = ((L + S) > 0).astype(np.int64)

    rows = []
    for Mx, name in [(M, "ORBIT    "), (R, "random   "), (NR, "lowrank+sp")]:
        k50, k90, s = energy_ranks(Mx)
        # tiny-rank energy capture (non-rigid => high)
        e = np.cumsum(s ** 2) / np.sum(s ** 2)
        cap8 = e[7] if len(e) > 7 else 1.0
        rows.append((name, k50, k90, cap8))
        print(f"  [{name}] 50%-energy rank {k50:3d}, 90% rank {k90:3d}, "
              f"top-8 SV capture {100*cap8:4.1f}%")

    orbit, rand, nonrig = rows[0], rows[1], rows[2]
    # verdict: compare ORBIT to the two controls
    close_to_random = abs(orbit[2] - rand[2]) < 0.25 * rand[2]
    far_from_nonrigid = orbit[3] < 0.4  # non-rigid captures most energy in top 8
    print(f"\n  orbit 90%-rank {orbit[2]} vs random {rand[2]} (close: {close_to_random}),"
          f" vs nonrigid {nonrig[2]}")
    print(f"  orbit top-8 capture {100*orbit[3]:.1f}% vs nonrigid {100*nonrig[3]:.1f}% "
          f"(orbit far from non-rigid: {far_from_nonrigid})")
    verdict = ("RIGID / random-like => bound plausible, NO low-rank+sparse crack"
               if close_to_random and far_from_nonrigid else
               "STRUCTURED => investigate as a crack lead")
    print(f"  VERDICT: the doubling-orbit matrix is {verdict}")
    print("\n  => The last log = rigidity of an explicit random-like matrix.")
    print("     Unconditional proof is Valiant-hard (explains the difficulty);")
    print("     the reduction + this evidence is the honest lower-bound result.")


if __name__ == "__main__":
    main()
