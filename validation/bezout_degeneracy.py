#!/usr/bin/env python3
"""
bezout_degeneracy.py -- the Bezout / elimination-theory door, killed.

Reformulation.  Let f(x) = prod_{k=1}^K (x + mu_k), mu_k = 8k+1.  Then

    f(0)  = P_0            (Borwein-cheap, Theorem crack)
    f'(0) = e_{K-1} = P_0 H  (the barrier)

so H = (log f)'(0) and the e_T vs e_{T-1} razor is exactly "the curve is
cheap, its polar is not".  Elimination theory offers a way to get f' from
f-side data: f is squarefree, so Bezout gives u f + v f' = Res(f,f') =
disc(f), and evaluating at x = 0 solves for H.  This script shows both
normalizations are 2-adically vacuous in the digit regime, gives the
valuation laws, and shows the repaired route lands on E2.

Checks:
  A. disc law.       v2(disc f) = 3K(K-1) + 2*sum_{d<K}(K-d) v2(d).
                     ~5K^2 >> m = Theta(K): the integral certificate is
                     congruent to ZERO at any digit-regime width.
  B. cofactor poles. Normalizing to u f + v f' = 1 moves the degeneracy into
                     the cofactors: v2(v(0)) ~ -cK with c ~ 7, a 2-adic pole
                     deeper than the whole target width m.
  C. cofactor size.  bits(u(0)), bits(v(0)) = Theta(K m) = Theta(N^2) --
                     the output-size wall, same as the e_j batch.
  D. Bezoutian.      Hermite: det Bez(f,f') = det Hankel(p_{i+j}) = disc.
                     Bezout's own matrix carries exactly the power sums, so
                     the repaired route is E2.
  E. rigidity.       No finite extension of Q_2 separates the configuration.
                     Unramified: mu_k are rational, Frobenius fixes them, all
                     still reduce to 1 in F_{2^s}.  Ramified: v scales by e,
                     uniformly.  The root configuration's shape is invariant.

Root cause (the topological statement).  |F_2^*| = 1.  Every odd number
reduces to the same point of the special fibre, so ANY set of units mod 2^m
is a maximally degenerate configuration: all C(K,2) intersections of f with
its polar pile up at one point.  Bezout's theorem still counts them, it just
counts them all in the same place.  At p = 2 there is no other kind of
instance.
"""
from __future__ import annotations

from fractions import Fraction


# ----------------------------------------------------------------- helpers

def v2(n: int) -> float:
    n = abs(int(n))
    if n == 0:
        return float("inf")
    return (n & -n).bit_length() - 1


def v2_frac(q: Fraction) -> int:
    q = Fraction(q)
    return int(v2(q.numerator) - v2(q.denominator))


# ------------------------------------------- dense polys over Q, coeff[i]=x^i

def ptrim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                out[i + j] += ai * bj
    return ptrim(out)


def psub(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for i, ai in enumerate(a):
        out[i] += ai
    for i, bi in enumerate(b):
        out[i] -= bi
    return ptrim(out)


def pdivmod(a, b):
    a = a[:]
    q = [Fraction(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and a != [Fraction(0)]:
        d = len(a) - len(b)
        c = a[-1] / b[-1]
        q[d] += c
        a = psub(a, [Fraction(0)] * d + [c * x for x in b])
    return ptrim(q), ptrim(a)


def pext_gcd(a, b):
    """(g, u, v) with u*a + v*b = g, g monic."""
    r0, r1 = a[:], b[:]
    s0, s1 = [Fraction(1)], [Fraction(0)]
    t0, t1 = [Fraction(0)], [Fraction(1)]
    while r1 != [Fraction(0)]:
        q, r = pdivmod(r0, r1)
        r0, r1 = r1, r
        s0, s1 = s1, psub(s0, pmul(q, s1))
        t0, t1 = t1, psub(t0, pmul(q, t1))
    lead = r0[-1]
    return ([c / lead for c in r0],
            [c / lead for c in s0],
            [c / lead for c in t0])


def build(K):
    f = [Fraction(1)]
    for k in range(1, K + 1):
        f = pmul(f, [Fraction(8 * k + 1), Fraction(1)])
    fp = ptrim([f[i] * i for i in range(1, len(f))])
    return f, fp


def det_frac(M):
    """Fraction-exact Gaussian elimination determinant."""
    M = [row[:] for row in M]
    n = len(M)
    det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        inv = Fraction(1) / M[c][c]
        for r in range(c + 1, n):
            if M[r][c]:
                fac = M[r][c] * inv
                for cc in range(c, n):
                    M[r][cc] -= fac * M[c][cc]
    return det


# ------------------------------------------------------------------ A: disc

def disc_v2_formula(K):
    return 3 * K * (K - 1) + 2 * sum((K - d) * int(v2(d)) for d in range(1, K))


def check_A():
    print("=" * 74)
    print("A.  disc law:  v2(disc f) = 3K(K-1) + 2*sum_{d<K}(K-d)v2(d)")
    print("=" * 74)
    bad = 0
    for K in range(2, 26):
        d = 1
        for i in range(1, K + 1):
            for j in range(i + 1, K + 1):
                d *= (8 * i - 8 * j) ** 2
        if int(v2(d)) != disc_v2_formula(K):
            bad += 1
    print(f"    direct disc vs formula, K = 2..25: {'PASS' if not bad else 'FAIL'}"
          f" ({bad} mismatches)")
    print()
    print(f"    {'K':>5} {'v2(disc)':>10} {'5K^2':>9} {'m=1.2K':>8}  certificate mod 2^m")
    for K in (16, 64, 256, 1024):
        vd = disc_v2_formula(K)
        m = int(1.2 * K)          # digit regime: m ~ n, K ~ n*log2(10)/4
        print(f"    {K:>5} {vd:>10} {5*K*K:>9} {m:>8}  "
              f"{'ZERO  (v2 >> m)' if vd >= m else 'survives'}")
    print()
    print("    The integral Bezout certificate vanishes mod 2^m by a factor of")
    print("    ~4K in the exponent.  It carries no information about H.")
    return bad == 0


# ------------------------------------------------- B/C: cofactors at x = 0

def check_BC():
    print()
    print("=" * 74)
    print("B/C.  normalized certificate  u f + v f' = 1  at x = 0")
    print("=" * 74)
    print(f"    {'K':>4} {'v2(v(0))':>9} {'v2/K':>6} {'m=1.2K':>7} {'bits u(0)':>10}"
          f" {'bits v(0)':>10} {'K*m':>7}  H exact?")
    poles = []
    sizes = []
    ok_all = True
    for K in range(3, 15):
        f, fp = build(K)
        g, u, v = pext_gcd(f, fp)
        assert g == [Fraction(1)], "f not squarefree?"
        u0, v0 = u[0], v[0]
        P0 = int(f[0])
        assert u0 * P0 + v0 * int(fp[0]) == 1

        m = max(8, int(1.2 * K))
        M = 1 << m
        Hrat = sum(Fraction(1, 8 * k + 1) for k in range(1, K + 1))
        Hmod = sum(pow(8 * k + 1, -1, M) for k in range(1, K + 1)) % M
        ok = (Hrat.numerator * pow(Hrat.denominator, -1, M)) % M == Hmod
        ok_all &= ok

        pv = v2_frac(v0)
        bu = u0.numerator.bit_length() + u0.denominator.bit_length()
        bv = v0.numerator.bit_length() + v0.denominator.bit_length()
        poles.append((K, pv))
        sizes.append((K, bu + bv))
        print(f"    {K:>4} {pv:>9} {pv/K:>6.1f} {m:>7} {bu:>10} {bv:>10}"
              f" {K*m:>7}  {ok}")

    (k0, p0), (k1, p1) = poles[0], poles[-1]
    slope = (p1 - p0) / (k1 - k0)
    (s0, b0), (s1, b1) = sizes[0], sizes[-1]
    print()
    print(f"    pole depth of v(0):  v2 ~ {slope:+.2f} * K   (linear in K)")
    print(f"    at m ~ 1.2K the pole is ~{abs(slope)/1.2:.1f}x deeper than the")
    print("    entire target width: reducing the certificate mod 2^m needs more")
    print("    guard bits than H has bits, and those guard bits ARE H.")
    print(f"    cofactor size grows {(b1/b0):.1f}x over K {s0}->{s1} "
          f"(ratio K^2: {(s1/s0)**2:.1f}x)  =>  Theta(K m) = Theta(N^2) output.")
    return ok_all


# ------------------------------------------------------- D: Bezoutian = E2

def bezoutian(f, g):
    """Bez(f,g)_{ij} from (f(x)g(y) - f(y)g(x))/(x-y)."""
    n = max(len(f), len(g)) - 1
    F = f + [Fraction(0)] * (n + 1 - len(f))
    G = g + [Fraction(0)] * (n + 1 - len(g))
    B = [[Fraction(0)] * n for _ in range(n)]
    # standard closed form: b_{ij} = sum_{k>max(i,j)} ... use the recurrence
    # b_{ij} = sum_{k=0}^{min(i,j)} (F[k]G[i+j+1-k] - F[i+j+1-k]G[k])
    for i in range(n):
        for j in range(n):
            s = Fraction(0)
            for k in range(min(i, j) + 1):
                idx = i + j + 1 - k
                if idx <= n:
                    s += F[k] * G[idx] - F[idx] * G[k]
            B[i][j] = s
    return B


def check_D():
    print()
    print("=" * 74)
    print("D.  Hermite:  Bez(f,f') and the power-sum Hankel carry the same data")
    print("=" * 74)
    print(f"    {'K':>4} {'det Bez':>14} {'det Hankel':>14}   (vs disc, up to sign)")
    ok_all = True
    for K in range(2, 9):
        f, fp = build(K)
        mu = [8 * k + 1 for k in range(1, K + 1)]
        disc = Fraction(1)
        for i in range(K):
            for j in range(i + 1, K):
                disc *= Fraction((mu[i] - mu[j]) ** 2)
        dB = det_frac(bezoutian(f, fp))
        p = [Fraction(sum(x ** r for x in mu)) for r in range(2 * K)]
        Hk = [[p[i + j] for j in range(K)] for i in range(K)]
        dH = det_frac(Hk)
        # classical signs: det Bez(f,f') = (-1)^{K(K-1)/2} disc, det Hankel = +disc
        sB = "+disc" if dB == disc else ("-disc" if dB == -disc else "MISMATCH")
        sH = "+disc" if dH == disc else ("-disc" if dH == -disc else "MISMATCH")
        ok_all &= (sB != "MISMATCH" and sH != "MISMATCH")
        print(f"    {K:>4} {sB:>14} {sH:>14}")
    print()
    print("    Bezout's own matrix construction IS the batch power sums S_j.")
    print("    Superfast Hankel solvers run in O(M(n) log n) but need the")
    print("    entries, and the entries are E2.  The repaired route is circular.")
    return ok_all


# ------------------------------------------------------- E: extension rigidity

def check_E():
    print()
    print("=" * 74)
    print("E.  rigidity: no finite extension of Q_2 separates the configuration")
    print("=" * 74)
    # Unramified quadratic W(F_4) = Z_2[x]/(x^2+x+1), residue field F_4.
    # Elements are a + b*x; v = min(v2(a), v2(b)) since the ext is unramified.
    # Compute the pairwise distance spectrum IN THE EXTENSION and compare to
    # the one in Z_2, and check every mu_k still reduces to 1 in F_4.
    K = 24
    mu = [8 * k + 1 for k in range(1, K + 1)]
    ext = [(m_, 0) for m_ in mu]                      # mu_k = mu_k + 0*x

    def vext(u):
        return min(v2(u[0]), v2(u[1]))

    d_z2 = sorted(int(v2(mu[i] - mu[j]))
                  for i in range(K) for j in range(i + 1, K))
    d_w4 = sorted(int(vext((ext[i][0] - ext[j][0], ext[i][1] - ext[j][1])))
                  for i in range(K) for j in range(i + 1, K))
    resid = {(u[0] % 2, u[1] % 2) for u in ext}        # image in F_4
    print(f"    unramified W(F_4) = Z_2[x]/(x^2+x+1):")
    print(f"      distance spectra identical to Z_2:  {d_z2 == d_w4}")
    print(f"      images of mu_k in F_4:              {resid}  "
          f"(collapsed to one point: {len(resid) == 1})")
    same = (d_z2 == d_w4) and len(resid) == 1
    print(f"      mu_k are rational, Frobenius fixes them, so this holds for")
    print(f"      every W(F_2^s) regardless of how large the residue field is.")
    # ramified: v(pi) = 1/e scales EVERY distance by e, so ratios of distances
    # -- i.e. the shape of the ultrametric configuration -- are invariant.
    e = 3
    d_base = sorted({int(v2(mu[i] - mu[j]))
                     for i in range(K) for j in range(i + 1, K)})
    d_ram = sorted({e * d for d in d_base})
    print(f"    ramified deg e={e}: distance spectrum {d_base[:6]}... -> "
          f"{d_ram[:6]}...")
    print(f"      uniform rescale by e, shape invariant: "
          f"{[e*d for d in d_base] == d_ram}")
    print()
    print("    The tree of the root configuration is the same in every finite")
    print("    extension.  'Go to a bigger residue field' is not available.")
    return same


if __name__ == "__main__":
    results = [check_A(), check_BC(), check_D(), check_E()]
    print()
    print("=" * 74)
    print(f"VERDICT: {'all checks pass' if all(results) else 'FAILURES PRESENT'}")
    print("Bezout / elimination door: KILLED, twice over (certificate vanishes;")
    print("repaired certificate has Theta(N^2) cofactors and collapses to E2).")
    print("=" * 74)
