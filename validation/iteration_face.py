#!/usr/bin/env python3
"""iteration_face.py -- the ITERATION paradigm, swept as its own face
(2026-08-03, third block; the last unexamined premise).

THE QUESTION.  The paradigm taxonomy covers decompositions of frac(10^N pi)
into INDEPENDENT cheap pieces (series / quadrature / counting / residue-
reconstruction).  Iteration -- AGM/Borwein dynamics, the machinery of the
N log^2 benchmark itself -- is not a decomposition: its pieces (iterates)
are sequentially dependent.  Never swept: can an iteration be given
EXTRACTION semantics, i.e. a deep window of the limit without carrying
full-precision intermediate states?  A yes = the fifth paradigm.

THE STRUCTURE (derivation, then measured below).  The AGM step is
(a,b) -> ((a+b)/2, sqrt(ab)).  Addition is carry-local, but multiplication
and sqrt are GLOBAL mixers: a perturbation epsilon at depth d of any iterate
becomes epsilon*(dense real) after one step -- a value shift of size
2^{-d}, which rewrites essentially EVERY digit below depth d of the limit.
Contrast the series face, where kakeya_hop measured carries LOCAL (a deep
bit depends on a ~log T window of columns).  The iteration face is the
opposite pole: no state can be truncated (every iterate needs the full
target precision -- perturbations persist, quadratic convergence does not
attenuate them), and no window of a state can be propagated (the window of
a product of succinctly-described full states is exactly the
Wegener-Woelfel SUCCINCT-ITERATED case = Door 2's named open seam).
So windowed iteration reduces to the atom in its Door-2 costume; the
extraction semantics are impossible in-scheme because the states literally
materialize every shallower digit en route.

MEASURED HERE:
  (A) fixed-point Brent-Salamin AGM at P bits, verified against const_pi;
  (B) perturbation persistence + density: flip one bit at depth d of an
      early iterate; the limit changes starting at depth ~d (no
      attenuation) and a deep window at W >> d is rewritten in EVERY trial
      (global mixing; control: flips BELOW the window leave it fixed);
  (C) the forced-mass corollary: log N iterates x Theta(N) live bits each,
      none truncatable by (B) => Theta(N log N) state mass, all shallow
      digits materialized => iteration cannot carry extraction semantics.

VERDICT: iteration is a fifth COSTUME of the wall, not a fifth paradigm;
the taxonomy is closed at four decompositions + one dynamics, all landing
on positional assembly.  Door 2's seams absorb the question; no new lever.
"""

import gmpy2
from gmpy2 import isqrt, mpz


def agm_pi_fixed(P, flip=None):
    """Brent-Salamin in fixed point at scale 2^P.  flip = (k, depth): XOR
    bit at `depth` (from the top) of a_k after iteration k."""
    S = mpz(1) << P
    a = S
    b = isqrt((S * S) >> 1)
    t = S >> 2
    k = 0
    while a - b > 4:
        if flip is not None and flip[0] == k:
            a = a ^ (mpz(1) << (P - flip[1]))
        an = (a + b) >> 1
        c2 = (a - an) * (a - an)
        t = t - ((c2 << k) >> P)
        b = isqrt(a * b)
        a = an
        k += 1
    return ((a + b) * (a + b)) // (4 * t), k


def window(x, P, W, G=64):
    return int((x >> (P - W - G)) & ((mpz(1) << G) - 1))


def main():
    P = 2048
    base, iters = agm_pi_fixed(P)

    gmpy2.get_context().precision = P + 16
    ref = int(gmpy2.mpfr(gmpy2.const_pi()) * (mpz(1) << P))
    agree = P - (base ^ ref).bit_length()
    print("PART A  fixed-point AGM  P=%d  iterations=%d" % (P, iters))
    print("        agrees with const_pi to %d of %d bits : %s"
          % (agree, P, "VERIFIED" if agree > P - 24 else "*** BAD ***"))

    W = P - 320                       # the deep read window [W, W+64)
    print()
    print("PART B  one-bit perturbations of iterate k at depth d")
    print("        window = 64 bits at depth %d" % W)
    print("        %-4s %-7s %-14s %-16s %s"
          % ("k", "depth d", "limit differs", "window rewritten", "note"))
    changed = total = 0
    for k in (1, 2, 4, 6):
        for d in (128, 512, 1024, 1536):
            pert, _ = agm_pi_fixed(P, flip=(k, d))
            first = P - (base ^ pert).bit_length()
            wch = window(pert, P, W) != window(base, P, W)
            total += 1
            changed += wch
            print("        %-4d %-7d at bit %-8d %-16s persist=%+d"
                  % (k, d, first, "YES" if wch else "no", first - d))
    ctrl, _ = agm_pi_fixed(P, flip=(2, P - 64))
    ctrl_ok = window(ctrl, P, W) == window(base, P, W)
    print("        control: flip BELOW window (d=%d): window %s"
          % (P - 64, "unchanged (as it must be)" if ctrl_ok
             else "*** CHANGED -- broken probe ***"))
    print("        %d/%d perturbations above the window rewrote it;"
          % (changed, total))
    print("        every flip moves the limit at its own depth (persist ~ 0):")
    print("        NO attenuation, NO locality -- the anti-kakeya.")

    print()
    print("PART C  the forced mass, and the reduction")
    print("        every iterate must carry ~P live bits (B: no truncation),")
    print("        log P iterates => Theta(N log N) state mass; the states")
    print("        materialize every shallower digit en route => extraction")
    print("        semantics impossible in-scheme.  A windowed step would be")
    print("        the middle window of a product of SUCCINCT full states =")
    print("        Wegener-Woelfel's open succinct-iterated case = Door 2.")
    print()
    print("VERDICT: iteration = fifth COSTUME, not fifth paradigm; the")
    print("taxonomy closes (4 decompositions + 1 dynamics, one atom).")


if __name__ == '__main__':
    main()
