#!/usr/bin/env python3
"""
fusion_mass.py -- BUILD the remainder-tree x geometric-multipoint fusion for the
counting face and MEASURE its intermediate mass. (The concrete upper-bound
target named by the 2026-07-23 sweeps: fuse amortization-over-moduli with
amortization-over-exponent for the grid {2^d mod mu_t}, aggregate output B(d).)

Column-sum generating function:
   F(x) = sum_d B(d) x^d = sum_t G_t(x)/(1 - x^{p_t})  mod x^N,
   p_t = ord(2, mu_t),  G_t(x) = sum_{j<p_t} [2^j mod mu_t > mu_t/2] x^j.

Two amortization levers and their measured cost:
 (L1) DENOMINATOR sharing (remainder-tree flavor): terms sharing a period p
      combine over one (1 - x^p); cost = sum over DISTINCT periods.
 (L2) NUMERATOR / orbit mass (geometric-multipoint flavor): each G_t needs the
      orbit threshold pattern, sum_t p_t total.

RESULT (recorded run, T = 100..1600, N = 4T; correctness verified T<=800):
   sum_t p_t          ~ T^1.876   (numerator / total orbit mass)
   sum distinct p     ~ T^1.856   (denominator mass -- sharing saves a CONSTANT
                                    ~1.2x, NOT an exponent)
   # distinct periods ~ T^0.914   (nearly linear count, but each period ~T)
=> MASS CONSERVED at Theta(T^2/polylog): the fusion computes B(d) correctly but
does NOT break the wall. The wall is pinned to the total multiplicative-orbit
mass  sum_t ord(2, mu_t) = Theta(T^2/polylog)  (Kurlberg-Pomerance: average
order of 2 is mu/polylog). The remainder-tree axis gives O~(T)/column x N
columns; the geometric axis gives sum p_t; combining them does NOT
cross-amortize (doubling touches every residue every step -- no column-to-column
locality). The "unclaimed gap between the amortizations" is hereby claimed,
built, measured, and closed: it conserves mass (13th independent confirmation).
The residual free structure: residues 2^d mod mu_t are a rational 1/(1-2x) mod
mu_t (free); ALL cost is the THRESHOLD pattern G_t = the incomplete-Kloosterman
data. Consistent with counting_face.py (thresholds were never the wall to
RESOLVE, but they ARE the wall to COMPRESS).
"""
import math
from sympy import n_order


def measure(T, N, verify):
    mus = [8 * t + 1 for t in range(1, T + 1)]
    pts = [n_order(2, mu) for mu in mus]
    sum_p = sum(pts)
    distinct = {}
    for p in pts:
        distinct[p] = distinct.get(p, 0) + 1
    sum_distinct = sum(distinct)
    ok = None
    if verify:
        B = [0] * N
        for mu in mus:
            for d in range(N):
                if pow(2, d, mu) > mu / 2:
                    B[d] += 1
        F = [0] * N
        for mu, p in zip(mus, pts):
            G = [1 if pow(2, j, mu) > mu / 2 else 0 for j in range(p)]
            for d in range(N):
                F[d] += G[d % p]
        ok = (F == B)
    return sum_p, sum_distinct, len(distinct), ok


def fit(xs, ys):
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(xs)
    sx, sy = sum(lx), sum(ly)
    sxx = sum(a * a for a in lx)
    sxy = sum(a * b for a, b in zip(lx, ly))
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def main():
    print("T      sum_p      sum_distinct  #distinct  correct")
    data = []
    for T in [100, 200, 400, 800, 1600]:
        sp, sd, nd, ok = measure(T, 4 * T, verify=(T <= 800))
        data.append((T, sp, sd, nd))
        print(f"{T:5d}  {sp:10d}  {sd:10d}  {nd:8d}   {ok if ok is not None else 'skip'}")
    Ts = [d[0] for d in data]
    e_p = fit(Ts, [d[1] for d in data])
    e_d = fit(Ts, [d[2] for d in data])
    e_n = fit(Ts, [d[3] for d in data])
    print(f"\n  sum_p        ~ T^{e_p:.3f}   (orbit/numerator mass; ~2 => wall)")
    print(f"  sum_distinct ~ T^{e_d:.3f}   (denominator mass; sharing = constant only)")
    print(f"  #distinct    ~ T^{e_n:.3f}")
    print(f"  naive T*N baseline ~ T^2.")
    verdict = "MASS CONSERVED (fusion does not break the wall)" if e_p > 1.6 \
        else "SUB-QUADRATIC -- INVESTIGATE"
    print(f"  VERDICT: {verdict}")


if __name__ == "__main__":
    main()
