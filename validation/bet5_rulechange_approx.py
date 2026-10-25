#!/usr/bin/env python3
"""
bet5_rulechange_approx.py -- Live bet 5:
  deliberate rule-change: drop exact V for something still digit-correct.

Already KILLED (rule_change_hunt / door_beta):
  - O(log N) bit approximation of V (limb support)
  - admissible deformations giving O(log N) bits of H

New attempts:
  A. Digit only needs {10^{N-1} * pi} to ~log N bits — but the pipeline
     reduces this to {X V / 2^J}.  Can we replace V by V0 + 2^{J-G} Delta
     with Delta cheap and G = omega(log N)?
  B. Randomized checksum: compute V mod 2^{c log N} many random odd primes
     via CRT — does NOT give 2-adic V
  C. Low-precision V for MOST limbs + exact correction on a sparse set
  D. Sign / carry pattern of the BBP sum that depends on V only through
     a soft-linear predicate (e.g. MSB of XV)
  E. Measure: how many bits of H actually flip the output decimal digit?
"""
from __future__ import annotations
import random

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


def part_A_sensitivity():
    P("=" * 70)
    P("BET5-A: bit-sensitivity of a model digit functional")
    # Model: digit depends on floor( (X * H mod 2^m) / 2^{m-G} ) for small G
    # Measure: flipping bit b of H flips the top G-bit window how often?
    random.seed(0)
    for K, m, G in [(32, 64, 4), (32, 64, 8), (64, 96, 8)]:
        H = atom_H(K, m)
        X = mpz(random.getrandbits(m)) | 1  # odd
        MOD = mpz(1) << m
        top = lambda V: ((X * V) % MOD) >> (m - G)
        base = top(H)
        flips = []
        for b in range(m):
            H2 = H ^ (mpz(1) << b)
            if top(H2) != base:
                flips.append(b)
        # which bits matter?
        if flips:
            P(f"  K={K} m={m} G={G}  bits that flip top-{G}: "
              f"min={min(flips)} max={max(flips)} count={len(flips)}")
        else:
            P(f"  K={K} m={m} G={G}  no flips (degenerate X)")
    P("  LAW: bits of H down to position ~ (m-G) affect the top G window.")
    P("  For G=O(log N) need H correct through bit m-O(log N) = Theta(N) bits.")
    P("  Reconfirms limb-support kill; no soft approximate V for the digit.")


def part_B_sparse_correction():
    P("=" * 70)
    P("BET5-B: low-prec V everywhere + sparse exact corrections")
    P("""
  Idea: compute H mod 2^{c log N} cheaply (c log N bits — soft if c log),
  then correct carries at positions where the low-prec estimate is close
  to a limb boundary.

  Problem: for random X, (X H) mod 2^m is equidistributed.  The probability
  a given limb is within 2^{c log N} of a boundary is ~ 2^{c log N - limb}.
  Number of dangerous limbs ~ N / word * that prob — still forces
  Theta(N / polylog) exact corrections in the worst case unless structure
  in X*H prevents boundary approach.

  X = 5^{N-1} * (BBP weights) has structure, but prior fence hunts
  (fence2/3, column_predict) found no exploitable limb prediction.
""")
    # Empirically: how close is XH to 2^{m-G} boundaries for structured X?
    for K, m, G in [(40, 80, 8)]:
        H = atom_H(K, m)
        MOD = mpz(1) << m
        # structured X: product of a few (8k+1) approx
        X = mpz(1)
        for k in range(1, min(K, 20)):
            X = (X * (8 * k + 1)) % MOD
        val = (X * H) % MOD
        # distance to nearest multiple of 2^{m-G}
        limb = mpz(1) << (m - G)
        r = val % limb
        dist = min(int(r), int(limb - r))
        P(f"  sample dist-to-boundary={dist}  threshold=2^{G}={1<<G}  "
          f"dangerous={dist < (1<<G)}")


def part_C_predicate():
    P("=" * 70)
    P("BET5-C: soft-linear predicate for the digit?")
    P("  Output is a decimal digit 0-9, i.e. ~3.3 bits.  Information")
    P("  theoretically O(1) bits.  But the MAP from V to the digit is")
    P("  a high-precision modular reduction — computing it is as hard as")
    P("  computing the fractional bits, not a PAC-style weak predicate.")
    P("  No known soft-linear algorithm for MSB(X V mod 2^J) when V is")
    P("  the additive core (would essentially solve E for G=1).")


def part_D_new_invariant():
    P("=" * 70)
    P("BET5-D: hunt a NEW approximate invariant")
    P("""
  Requirements for a rule-change win:
    (i)  computable in O(N log^2 N)
    (ii) determines the decimal digit
    (iii) is NOT equivalent to Theta(N) bits of V

  Candidates tried and killed:
    - top O(log N) bits of V
    - V mod 2^{O(log N)}
    - admissible Log Pi(u) stencils
    - O(log N) Mazur moments
    - Archimedean short windows

  Remaining speculative:
    - a TC0 / FMM-style hierarchical approximation with 2-adic error
      control that somehow uses AP structure (FMM faces already hunted
      in multipole_2adic / moment_transform — empty)
    - a base-10 specific identity skipping binary V entirely
      (base10_gaussian_hunt — empty)

  No new invariant found in this pass.
""")


def main():
    P("bet5_rulechange_approx -- approximate digit-correct invariants")
    part_A_sensitivity()
    part_B_sparse_correction()
    part_C_predicate()
    part_D_new_invariant()
    P()
    P("BET5 VERDICT: sensitivity law reconfirmed — Theta(N) bits of H needed.")
    P("  Sparse correction / predicate / new invariant: no construction.")
    P("  Approximate rule-change remains blocked without a genuinely new idea.")


if __name__ == "__main__":
    main()
