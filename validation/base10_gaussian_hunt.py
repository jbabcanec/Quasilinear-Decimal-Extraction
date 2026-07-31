#!/usr/bin/env python3
"""base10_gaussian_hunt.py -- the no-band dream, and the one loophole in BBG04.
(2026-07-30)

Why base 10 is the biggest prize on the formula axis:

  With a base-2^w formula the extraction splits into a word phase (terms whose
  dyadic shift is still an integer -- one modexp each, cheap) and the BAND
  (terms landing inside the bit-string of X = 5^(N-1)).  The band is
  LX = 2.32 N bits wide, carries ~LX/g_eff terms, and IS the whole war.

  With a base-10^w formula there is no band at all: 10^N / 10^{wk} is an
  integer for every wk <= N, so every term is word phase, and terms past
  wk > N die faster than the window.  Extraction would cost one modexp per
  term = O(N loglog N) word operations ~ N log N bit ops, BELOW parity --
  better than the prize the campaign is chasing.

Borwein-Borwein-Galway (Canad. J. Math 56 (2004)) proved no MACHIN-TYPE
(degree 1, rational coefficients) BBP formula for pi exists in a non-binary
base.  That leaves an explicit loophole: DEGREE 2 with an irreducible
quadratic denominator, i.e. Gaussian denominators (mk+j)^2 + c^2.  These are
still word-size integers, so they are still pipeline-legal -- the word phase
just runs in Z[i]/(mu).  Nobody has searched that space for base 10.

Experiments (all with controls):
  C0  harness control: log(10/9) from the m=1 base-10 degree-1 basis.
  C1  harness control: pi from the m=8 base-16 degree-1 basis (BBP).
  C2  pi from degree-1 bases at b = 5, 10, 100 -- must be EMPTY (BBG04).
  C3  pi from Gaussian degree-2 bases at b = 10, 100 -- the open loophole.
"""

from mpmath import mp, mpf, pslq, log as mp_log, pi as mp_pi


def deg1_basis(m, b, terms):
    """S_j = sum_k b^{-k} / (mk+j), j = 1..m."""
    out = []
    for j in range(1, m + 1):
        s = mpf(0)
        num = mpf(1)
        for k in range(terms):
            s += num / (m * k + j)
            num /= mpf(b)
        out.append(s)
    return out, ["S_%d" % j for j in range(1, m + 1)]


def gauss_basis(m, b, cset, terms):
    """Re/Im parts of sum_k b^{-k}/((mk+j) - c i):
       A = sum b^{-k}(mk+j)/((mk+j)^2+c^2),  B = sum b^{-k} c/((mk+j)^2+c^2).
       Denominators (mk+j)^2+c^2 are word-size integers -> pipeline-legal."""
    out, names = [], []
    for c in cset:
        for j in range(1, m + 1):
            a = mpf(0)
            bb = mpf(0)
            num = mpf(1)
            for k in range(terms):
                d = mpf((m * k + j) ** 2 + c * c)
                a += num * (m * k + j) / d
                bb += num * c / d
                num /= mpf(b)
            out.append(a)
            names.append("A_%d,%d" % (j, c))
            out.append(bb)
            names.append("B_%d,%d" % (j, c))
    return out, names


def find_relation(target, basis, names, maxcoeff=10 ** 6, dps=220, tries=6):
    """PSLQ for target in span(basis). Internal relations that do not involve
    the target are stripped by dropping their heaviest member and retrying."""
    mp.dps = dps
    vec = [target] + list(basis)
    nm = ["TARGET"] + list(names)
    for _ in range(tries):
        rel = pslq(vec, maxcoeff=maxcoeff, maxsteps=30000,
                   tol=mpf(10) ** (-(dps - 30)))
        if not rel:
            return None
        if rel[0] == 0:
            # internal relation among the basis: drop its heaviest member
            k = max(range(1, len(rel)), key=lambda i: abs(rel[i]))
            del vec[k]
            del nm[k]
            continue
        # verify at higher precision by rebuilding is caller's job; here we
        # re-check the residual against the same (already converged) values
        resid = sum(c * x for c, x in zip(rel, vec))
        if abs(resid) > mpf(10) ** (-(dps - 25)):
            return None
        return [(nm[i], rel[i]) for i in range(len(rel)) if rel[i] != 0]
    return None


def report(label, target, basis, names, **kw):
    hit = find_relation(target, basis, names, **kw)
    if hit:
        r = sum(1 for n, c in hit if n != "TARGET")
        print("  %-46s FOUND  r=%d  %s"
              % (label, r, "  ".join("%s*%s" % (c, n) for n, c in hit)))
    else:
        print("  %-46s empty" % label)
    return hit


def main():
    mp.dps = 260
    print("=" * 78)
    print("base-10 / Gaussian-denominator hunt: can the BAND be removed entirely?")
    print("=" * 78)

    print("controls (the harness must find these):")
    bas, nm = deg1_basis(1, 10, 300)
    report("C0  log(10/9) from base-10, m=1, degree 1", mp_log(mpf(10) / 9),
           bas, nm)
    bas, nm = deg1_basis(8, 16, 260)
    report("C1  pi from base-16, m=8, degree 1 (BBP)", mp_pi, bas, nm)

    print()
    print("C2  degree-1 formulas for pi in non-binary bases (BBG04 says none):")
    for b in (5, 10, 20, 50, 100, 1000):
        for m in (1, 2, 4, 6, 8, 10, 12):
            terms = int(260 * 3.33 / max(1.0, __import__('math').log10(b))) + 10
            bas, nm = deg1_basis(m, b, min(terms, 1400))
            report("    b=%-5d m=%-3d" % (b, m), mp_pi, bas, nm,
                   maxcoeff=10 ** 5, dps=200)

    print()
    print("C3  degree-2 GAUSSIAN denominators (mk+j)^2+c^2 -- BBG04's loophole:")
    for b in (10, 100):
        for m in (1, 2, 4):
            for cset in ([1], [1, 2], [1, 2, 3]):
                terms = int(260 * 3.33 / __import__('math').log10(b)) + 10
                bas, nm = gauss_basis(m, b, cset, min(terms, 1000))
                d1, n1 = deg1_basis(m, b, min(terms, 1000))
                report("    b=%-4d m=%-2d c=%-9s" % (b, m, cset), mp_pi,
                       bas + d1, nm + n1, maxcoeff=10 ** 4, dps=230)

    print()
    print("=" * 78)
    print("READING: a hit in C3 would mean band-free decimal extraction at")
    print("O(N loglog N) word ops -- below the parity prize. An empty C3 makes")
    print("BBG04's degree-1 exclusion effectively degree-2-tight for base 10.")


if __name__ == "__main__":
    main()
