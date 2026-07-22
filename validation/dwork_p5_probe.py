#!/usr/bin/env python3
"""
dwork_p5_probe.py -- the p = 5 Dwork tower for the Gaussian-pipeline family
(source-level move: re-plumb through Zudilin's Z[i] series => additive core at
p = 5, where excellent Frobenius lifts are NOT vacuous).

p = 2 baseline (validation/dwork_probe.py etc., 2026-07-21): for f(x) =
sum x^t/(d t+1), the Dwork ratio
    D_s = F_{p^{s+1}}(x) * F_{p^{s-1}}(x^p)  -  F_{p^s}(x) * F_{p^s}(x^p)
has min coefficient valuation obeying the cross-formula law
    min v_2 = s + v_2(d) - 1   (slope 1, offset from d; measured d = 4..16).

THIS PROBE: same object at p = 5, families d = 5, 10, 25 (all dt+1 prime to 5).
Predicted analog law (if the tower is characteristic-uniform):
    min v_5(D_s) = s + v_5(d) - 1.
Anything steeper than slope 1 at p = 5 would be the first sign of amplification
beyond the p = 2 ceiling; slope 1 = the wall relocates, as mass conservation
predicts. Either outcome is new data (odd-p tower for this family: unmeasured).

Exact arithmetic mod 5^L, coefficientwise valuations.
"""

P = 5
L = 18
MOD = P ** L


def v5(c):
    if c % MOD == 0:
        return L  # at-or-above tracking precision
    c %= MOD
    v = 0
    while c % P == 0:
        c //= P
        v += 1
    return v


def F_trunc(n, d):
    """Coefficients of F_n(x) = sum_{t<n} x^t/(d t+1) mod 5^L."""
    return [pow(d * t + 1, -1, MOD) for t in range(n)]


def spread(coeffs):
    """g(x) -> g(x^P)."""
    out = [0] * ((len(coeffs) - 1) * P + 1)
    for t, c in enumerate(coeffs):
        out[t * P] = c
    return out


def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj:
                    out[i + j] = (out[i + j] + ai * bj) % MOD
    return out


def dwork_ratio_minval(s, d):
    A = polymul(F_trunc(P ** (s + 1), d), spread(F_trunc(P ** (s - 1), d)))
    B = polymul(F_trunc(P ** s, d), spread(F_trunc(P ** s, d)))
    n = max(len(A), len(B))
    A += [0] * (n - len(A))
    B += [0] * (n - len(B))
    vals = [v5(A[t] - B[t]) for t in range(n)]
    # index below p^s must vanish identically (exact agreement region), per p=2 pattern
    low_exact = all(vals[t] >= L for t in range(min(P ** s, n)))
    tail = [v for t, v in enumerate(vals) if t >= P ** s and v < L]
    return (min(tail) if tail else L), low_exact


def main():
    print(f"p = {P}, precision 5^{L}. Dwork ratio min valuations (predicted s + v5(d) - 1):")
    for d in [5, 10, 25]:
        vd = 0
        dd = d
        while dd % P == 0:
            dd //= P
            vd += 1
        row = []
        for s in [1, 2, 3]:
            mv, low = dwork_ratio_minval(s, d)
            row.append((s, mv, low))
        pred = [s + vd - 1 for s, _, _ in row]
        meas = [mv for _, mv, _ in row]
        lows = all(low for _, _, low in row)
        print(f"  d={d} (v5(d)={vd}): measured {meas}, predicted {pred}, "
              f"below-p^s exact-zero region: {lows}   "
              f"{'LAW HOLDS' if meas == pred else 'DEVIATION — investigate'}")


if __name__ == "__main__":
    main()
