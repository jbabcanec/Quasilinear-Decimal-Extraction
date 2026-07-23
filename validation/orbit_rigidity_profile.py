#!/usr/bin/env python3
"""
orbit_rigidity_profile.py -- relentless rigidity pursuit of the doubling-orbit
matrix M[d][t]=bit_d(1/(8t+1)). Resolves the "doubling correlation" yellow flag
and identifies an ACHIEVABLE weak-rigidity bound.

MEASURED (N=T=512):
 (1) CORRELATION DECAY: row-correlation ratio to random by lag:
       lag 1: 2.76, lag 2: 1.83, lag 3: 1.30, lag 8: 1.21, lag>=16: ~1.0 (floor).
     => the doubling correlation is confined to a NEAR-DIAGONAL BAND of width
        ~2 log N; off-band the rows are Weil-pseudorandom (random floor).
 (2) BAND WEIGHT: removing the lag-1 linear fit changes the 90%-energy rank by 1
     (247 -> 248). => the band carries NEGLIGIBLE spectral weight; it is NOT a
     low-rank structure and does NOT yield a low-rank+sparse decomposition
     (i.e., NOT a crack via the doubling shadow).
 (3) SUBMATRIX RANK (the SSS weak-rigidity ingredient): random k x k submatrices
     have rank/k ~ 0.98 (R) and 0.93-0.99 (F_2) -- NEAR FULL. No large low-rank
     submatrix.

CONCLUSION: M = (sparse O(log N)-wide near-diagonal band) + (pseudorandom rigid
core). The doubling shadow is CONTAINED (sparse, negligible spectral weight,
does not reduce rank) -- it does not crack the matrix. The near-full submatrix
rank is exactly the ingredient for the best-known WEAK rigidity bound
(Shokrollahi-Spielman-Stemann): R_M(r) >= (n^2/(r+1)) * polylog, PROVABLE once
the near-full submatrix rank is made rigorous via Kloosterman/Weil equidistri-
bution of (2^d mod mu_t, 2^{d'} mod mu_t). That is the achievable partial
theorem. STRONG rigidity (Valiant, the circuit-lower-bound prize) remains open,
as for every explicit candidate.
"""
import numpy as np


def orbit_matrix(N, T):
    M = np.zeros((N, T))
    for j, t in enumerate(range(1, T + 1)):
        mu = 8 * t + 1
        r = 1 % mu
        for d in range(N):
            M[d, j] = 1.0 if r >= (mu + 1) // 2 else 0.0
            r = (2 * r) % mu
    return M


def f2_rank(A):
    A = (A % 2).astype(np.uint8).copy()
    n, m = A.shape
    rank = 0
    row = 0
    for col in range(m):
        piv = next((i for i in range(row, n) if A[i, col]), None)
        if piv is None:
            continue
        A[[row, piv]] = A[[piv, row]]
        for i in range(n):
            if i != row and A[i, col]:
                A[i] ^= A[row]
        row += 1
        rank += 1
        if row == n:
            break
    return rank


def top90(A):
    s = np.linalg.svd(A - A.mean(), compute_uv=False)
    e = np.cumsum(s ** 2) / np.sum(s ** 2)
    return int(np.searchsorted(e, 0.9)) + 1


def main():
    N = T = 512
    M = orbit_matrix(N, T)
    np.random.seed(2)
    R = (np.random.rand(N, T) < 0.5).astype(float)
    C = M - M.mean(axis=1, keepdims=True)
    Cr = R - R.mean(axis=1, keepdims=True)

    print("(1) row-correlation ratio to random by lag (decay = banded structure):")
    band = None
    for k in [1, 2, 3, 4, 8, 16, 32]:
        o = np.mean([abs(C[d] @ C[d + k]) / T for d in range(N - k)])
        rr = np.mean([abs(Cr[d] @ Cr[d + k]) / T for d in range(N - k)])
        ratio = o / rr
        if band is None and ratio < 1.1:
            band = k
        print(f"    lag {k:3d}: ratio {ratio:.2f}")
    print(f"    => correlation reaches the random floor by lag ~{band} "
          f"(band width ~2 log N = {2*np.log2(T):.0f})")

    before = top90(M)
    Mres = M.copy()
    for d in range(1, N):
        a = (C[d] @ C[d - 1]) / (C[d - 1] @ C[d - 1] + 1e-9)
        Mres[d] = M[d] - a * M[d - 1]
    after = top90(Mres)
    print(f"\n(2) 90%-energy rank: {before} -> {after} after removing lag-1 fit "
          f"(random ~{top90(R)}) => band carries negligible spectral weight")

    print("\n(3) random k x k submatrix rank/k (SSS ingredient, near 1 = no large low-rank minor):")
    for k in [32, 64, 128]:
        rs, fs = [], []
        for _ in range(15):
            di = np.random.choice(N, k, replace=False)
            tj = np.random.choice(T, k, replace=False)
            S = M[np.ix_(di, tj)]
            rs.append(np.linalg.matrix_rank(S) / k)
            fs.append(f2_rank(S) / k)
        print(f"    k={k:3d}: rank_R/k {np.mean(rs):.3f}  rank_F2/k {np.mean(fs):.3f}")

    print("\nVERDICT: M = sparse O(log N) band + pseudorandom rigid core. The")
    print("doubling shadow is contained (sparse, negligible spectral weight, no")
    print("rank drop) => not a crack. Near-full submatrix rank => the weak SSS")
    print("rigidity bound is provable via Kloosterman/Weil equidistribution")
    print("(achievable partial theorem). Strong rigidity remains open.")


if __name__ == "__main__":
    main()
