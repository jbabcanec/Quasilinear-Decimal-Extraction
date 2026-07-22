#!/usr/bin/env python3
"""
excellent_lift_p5.py -- empirical test of the Beukers-Vlasenko excellent-lift
amplification (DC-III Conjecture 7.5 analog) in a RAMIFIED case: our family
F(t) = sum_n t^n/(5n+1) at p = 5 (5 divides the parameter denominator, so no
B-V theorem applies; the classical tower nevertheless holds with measured law
min v_5(D_s) = s, `dwork_p5_probe.py`).

Object: D_s(sigma) = F_{5^{s+1}}(t) * F_{5^{s-1}}(t^sigma)
                   - F_{5^s}(t)   * F_{5^s}(t^sigma),
with Frobenius lift t^sigma = t^5 * (1 + 5*u(t)). Classical lift u = 0 gives
min val exactly s. B-V excellent lifts (where their theory applies) should
give ~2s. Protocol here:
  (1) solve u mod 5 from s=1 (push val(D_1) from 1 to >= 2),
  (2) refine u mod 25 from s=2 (push val(D_2) toward >= 4),
  (3) FREEZE sigma and measure s=3: the conjecture-analog PREDICTS >= 6
      against classical 3. The s=3 number is the finding, whatever it is.
Exact arithmetic mod 5^L throughout.
"""

P, L = 5, 12
MOD = P ** L


def v5(c):
    c %= MOD
    if c == 0:
        return L
    v = 0
    while c % P == 0:
        c //= P
        v += 1
    return v


def mul(a, b, T):
    out = [0] * T
    for i, ai in enumerate(a):
        if ai:
            hi = min(T - i, len(b))
            for j in range(hi):
                if b[j]:
                    out[i + j] = (out[i + j] + ai * b[j]) % MOD
    return out


def coeffs(N):
    return [pow(5 * n + 1, -1, MOD) for n in range(N)]


def make_sigma(u, T):
    """t^sigma = t^5 * (1 + 5 u(t)) as series list mod t^T."""
    base = [0] * T
    if 5 < T:
        base[5] = 1
    for k, uk in enumerate(u):
        if 5 + k < T and uk:
            base[5 + k] = (base[5 + k] + 5 * uk) % MOD
    return base


def compose_F(N, ts, T):
    """F_N(t^sigma) mod t^T; stop when n*val(ts) >= T."""
    vts = next((i for i, c in enumerate(ts) if c % MOD), T)
    c = coeffs(N)
    acc = [0] * T
    cur = [0] * T
    cur[0] = 1
    for n in range(N):
        if n * vts >= T:
            break
        for i in range(T):
            if cur[i]:
                acc[i] = (acc[i] + c[n] * cur[i]) % MOD
        cur = mul(cur, ts, T)
    return acc


def direct_F(N, T):
    c = coeffs(min(N, T))
    return [c[n] if n < N else 0 for n in range(T)]


def D_s(s, u, T):
    ts = make_sigma(u, T)
    A = mul(direct_F(P ** (s + 1), T), compose_F(P ** (s - 1), ts, T), T)
    B = mul(direct_F(P ** s, T), compose_F(P ** s, ts, T), T)
    return [(A[i] - B[i]) % MOD for i in range(T)]


def minval(D, lo):
    vals = [v5(c) for c in D[lo:]]
    return min(vals) if vals else L


def solve_level(s, u, T, K, target_val):
    """Adjust u by multiples of 5^(target_val-2) t^k to push val(D_s) to target."""
    step = P ** (target_val - 2)
    scale = P ** (target_val - 1)
    base = D_s(s, u, T)
    cur = minval(base, P ** s)
    if cur >= target_val:
        return u, True, cur, True
    if cur < target_val - 1:
        return u, False, cur, False   # too coarse to fix at this step scale
    bvec = [(-(base[i] // scale)) % P for i in range(T)]
    # sensitivities: direction u + step * t^k
    Svecs = []
    for k in range(K):
        du = list(u) + [0] * max(0, K - len(u))
        du[k] = (du[k] + step) % MOD
        Dk = D_s(s, du, T)
        Svecs.append([((Dk[i] - base[i]) // scale) % P for i in range(T)])
    # solve sum x_k S_k = bvec over F_5 (Gaussian elimination on columns)
    rows = list(range(T))
    x = [0] * K
    piv = []
    M = [[Svecs[k][i] for k in range(K)] for i in range(T)]
    rhs = [bvec[i] for i in range(T)]
    r = 0
    for c in range(K):
        pr = None
        for i in range(r, T):
            if M[i][c] % P:
                pr = i
                break
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        rhs[r], rhs[pr] = rhs[pr], rhs[r]
        inv = pow(M[r][c], -1, P)
        M[r] = [(m * inv) % P for m in M[r]]
        rhs[r] = (rhs[r] * inv) % P
        for i in range(T):
            if i != r and M[i][c] % P:
                f = M[i][c]
                M[i] = [(M[i][j] - f * M[r][j]) % P for j in range(K)]
                rhs[i] = (rhs[i] - f * rhs[r]) % P
        piv.append((r, c))
        r += 1
    # consistency: any nonzero rhs row with zero M-row => unsolvable
    solvable = all(any(M[i][j] % P for j in range(K)) or rhs[i] % P == 0
                   for i in range(T))
    for (rr, cc) in piv:
        x[cc] = rhs[rr]
    unew = list(u) + [0] * (K - len(u))
    for k in range(K):
        unew[k] = (unew[k] + step * x[k]) % MOD
    achieved = minval(D_s(s, unew, T), P ** s)
    return unew, achieved >= target_val, achieved, solvable


if __name__ == "__main__":
    print("Excellent-lift hunt for the 5-ramified family F = sum t^n/(5n+1):")
    u = []
    # baseline
    for s, T in [(1, 45), (2, 90)]:
        print(f"  classical val(D_{s}) = {minval(D_s(s, [], T), P**s)} (law: {s})")
    # level 1: solve u mod 5 to push D_1 val 1 -> 2
    u, ok1, ach1, solv1 = solve_level(1, [], 45, 30, 2)
    print(f"  s=1 solve (u mod 5): solvable={solv1}, achieved val(D_1) = {ach1}")
    # level 2: refine u mod 25 to push D_2 -> 4 (via 3 first)
    for tgt in [3, 4]:
        u, ok2, ach2, solv2 = solve_level(2, u, 90, 30, tgt)
        print(f"  s=2 refine (target {tgt}): solvable={solv2}, achieved val(D_2) = {ach2}")
        if not ok2:
            break
    # s=3 with frozen sigma: THE PREDICTION TEST (classical 3; conjecture-analog ~6)
    T3 = 220
    d3 = minval(D_s(3, u, T3), P ** 3)
    print(f"  s=3 FROZEN-sigma measurement: val(D_3) = {d3} "
          f"(classical law 3; excellent-lift analog predicts ~6)")
    print("FINDING:", "amplification persists in the ramified case"
          if d3 >= 5 else
          "amplification does NOT persist in the ramified case"
          if d3 <= 4 else "ambiguous")
