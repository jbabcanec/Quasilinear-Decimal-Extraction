#!/usr/bin/env python3
"""
stress_dual_claims.py -- stress-test the load-bearing claims of the
Open Problem E campaign (dual squaring lemma, lift rate law, bet verdicts).

Tests:
  A. Lemma at scale: v2(B(m,K)^2 - P0*B(m+1,K)) >= 2m for m up to 1000
     (previously machine-checked only to m=21).
  B. SHARP law conjecture: v2(err) = 2m + v2(K) exactly.
     Reason: err/(P0*B(m+1)) = ratio-1 = 2^{2m}*C + O(2^{4m}),
     C = sum 1/(mu(mu+2R)) === Q2 mod 2^{m+1}, and Q2 === K mod 16
     (mu^{-2} === 1 mod 16), deeper terms all have v2 > v2(K).
  C. Jump margin law: v2(B(m)^{2^s} - P0^{2^s-1} B(m+s)) = 2m + (s-1) + v2(K).
     Corollary: full-width jump iff s <= v2(K) - 1 (not a flat "s>=3 fails").
  D. Repaired one-step lift: the bet1-C "1-2 bits short" lift is repairable
     FOR FREE with correction factor (1 + 2^{2(m-1)}*Ccheap), Ccheap from
     O(1) Faulhaber terms.  Surplus bits ~ depth of Ccheap.
  E. Chained lift rate law: with O(1)-term cheap corrections, final precision
     of the chain m0 -> M grows like M + const, NOT 2M.  With exact
     corrections (Theta(K) work per level) chain is exact -- the N^2 wall.
  F. Wider monomial relation hunt (meet-in-the-middle, 7 atoms, exps +-4):
     at width 2m-4 we must re-find the lemma family (sanity);
     at width 3m there must be NO relation (else a geometric jump exists).
     Previous hunt (dual_jump part B): 5 atoms, exps +-2, width 2m only.
  G. e_{T-1} smoothness at larger K (to 64), stripping to 10^6 and
     Miller-Rabin on the cofactor (bet3-A only went to K=30, strip 8K+1).
  H. Paper sanity: Prop E-equiv items (ii) and (iii) at larger sizes.
"""
from __future__ import annotations
import math

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def v2m(x, W):
    """v2 of x computed mod 2^W; returns W if x === 0 (meaning >= W)."""
    x = int(x) % (1 << W)
    if x == 0:
        return W
    return (x & -x).bit_length() - 1


def B(m, K, MOD=None):
    R = mpz(1) << m
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= R + 8 * k + 1
        if MOD is not None:
            acc %= MOD
    return acc


def P0(K, MOD=None):
    acc = mpz(1)
    for k in range(1, K + 1):
        acc *= 8 * k + 1
        if MOD is not None:
            acc %= MOD
    return acc


def atom_H(K, m):
    MOD = mpz(1) << m
    return sum(pow(mpz(8 * k + 1), -1, MOD) for k in range(1, K + 1)) % MOD


# ---------------------------------------------------------------------------
# A. Lemma at scale
# ---------------------------------------------------------------------------
def part_A_scale():
    P("=" * 72)
    P("A: lemma v2(err) >= 2m at scale (log only had m<=21)")
    for m in [50, 100, 200, 500, 1000]:
        for K in [m, 2 * m]:
            W = 2 * m + 64
            MOD = mpz(1) << W
            err = (B(m, K, MOD) ** 2 - P0(K, MOD) * B(m + 1, K, MOD)) % MOD
            ve = v2m(err, W)
            pred = 2 * m + v2m(K, 60)
            P(f"  m={m:5d} K={K:5d}  v2(err)={ve:6d}  2m={2*m:6d}  "
              f"pred(2m+v2K)={pred:6d}  lemma_ok={ve >= 2*m}  sharp={ve == pred}")


# ---------------------------------------------------------------------------
# B. Sharp law sweep
# ---------------------------------------------------------------------------
def part_B_sharp():
    P("=" * 72)
    P("B: sharp law v2(err) == 2m + v2(K), sweep m=6..40, K=1..128")
    bad = []
    total = 0
    for m in range(6, 41):
        W = 2 * m + 40
        MOD = mpz(1) << W
        p0 = mpz(1)
        bm = mpz(1)
        bm1 = mpz(1)
        R = mpz(1) << m
        R2 = mpz(2) << m
        for K in range(1, 129):
            mu = 8 * K + 1
            p0 = p0 * mu % MOD
            bm = bm * (R + mu) % MOD
            bm1 = bm1 * (R2 + mu) % MOD
            err = (bm * bm - p0 * bm1) % MOD
            ve = v2m(err, W)
            pred = 2 * m + v2m(K, 30)
            total += 1
            if ve != pred:
                bad.append((m, K, ve, pred))
    P(f"  checked {total} pairs; exceptions: {len(bad)}")
    for m, K, ve, pred in bad[:20]:
        P(f"    EXC m={m} K={K}  v2={ve}  pred={pred}")
    P("  high-v2(K) spot: m=30, K=2^j")
    m = 30
    W = 2 * m + 40
    MOD = mpz(1) << W
    for j in range(1, 13):
        K = 1 << j
        err = (B(m, K, MOD) ** 2 - P0(K, MOD) * B(m + 1, K, MOD)) % MOD
        ve = v2m(err, W)
        P(f"    K=2^{j:<2d}  v2(err)={ve:4d}  pred={2*m+j:4d}  ok={ve == 2*m+j}")


# ---------------------------------------------------------------------------
# C. Jump margin law
# ---------------------------------------------------------------------------
def part_C_jump():
    P("=" * 72)
    P("C: jump law v2(B(m)^{2^s} - P0^{2^s-1}B(m+s)) == 2m+(s-1)+v2(K)")
    P("   => full-width (need 2(m+s)) iff s <= v2(K)-1")
    bad = 0
    for m in [8, 10, 12]:
        for K in [8, 12, 16, 32, 48, 64]:
            for s in range(1, 6):
                need = 2 * (m + s)
                W = need + 40
                MOD = mpz(1) << W
                lhs = pow(B(m, K, MOD), 1 << s, MOD)
                rhs = pow(P0(K, MOD), (1 << s) - 1, MOD) * B(m + s, K, MOD) % MOD
                ve = v2m(lhs - rhs, W)
                pred = 2 * m + (s - 1) + v2m(K, 30)
                full = ve >= need
                pfull = s <= v2m(K, 30) - 1
                tag = "" if (ve == pred and full == pfull) else "  <-- EXC"
                if tag:
                    bad += 1
                if s <= 3 or tag:
                    P(f"  m={m} K={K:3d} s={s}  v2={ve:3d} pred={pred:3d} "
                      f"need={need:3d} full={full} pred_full={pfull}{tag}")
    P(f"  exceptions: {bad}")


# ---------------------------------------------------------------------------
# D. Repaired one-step lift
# ---------------------------------------------------------------------------
def S_j(K, j):
    return sum(mpz(k) ** j for k in range(1, K + 1))


def C_cheap(K, R, depth_terms):
    """Faulhaber expansion of C = sum 1/(mu(mu+2R)) === Q2 mod 2^{m+1}.
    Q2 = sum mu^{-2} = sum_j (j+1)(-8)^j S_j(K).  depth_terms = # of j terms."""
    acc = mpz(0)
    for j in range(depth_terms):
        acc += (j + 1) * (mpz(-8) ** j) * S_j(K, j)
    return acc


def part_D_repair():
    P("=" * 72)
    P("D: one-step lift m-1 -> m repaired with cheap Faulhaber correction")
    P("   naive short = 2 - v2K-ish (bet1-C); corrected surplus ~ 3*depth")
    for m, K in [(16, 16), (20, 20), (24, 17), (24, 32)]:
        W = 2 * m + 30
        MOD = mpz(1) << W
        p0 = P0(K, MOD)
        prev = B(m - 1, K, MOD)
        exact = B(m, K, MOD)
        inv_p0 = pow(p0, -1, MOD)
        naive = prev * prev % MOD * inv_p0 % MOD
        v_naive = v2m(naive - exact, W)
        row = f"  m={m} K={K:3d}  naive_v2={v_naive:3d} (need {2*m})"
        for depth in [1, 2, 4, 8]:
            Cc = C_cheap(K, mpz(1) << (m - 1), depth) % MOD
            corr = (1 + (mpz(1) << (2 * (m - 1))) * Cc) % MOD
            fixed = naive * pow(corr, -1, MOD) % MOD
            v_fix = v2m(fixed - exact, W)
            row += f"  d{depth}:{v_fix}"
        P(row + f"   (W={W} means >=)")
    P("  reading: depth-1 (C~=K, free) already covers need=2m;")
    P("  each extra term buys ~3 more bits; Theta(m) surplus = Theta(m) terms")
    P("  = Faulhaber/Bernoulli mass = E3 wall.  Chain still rate-limited (E).")


# ---------------------------------------------------------------------------
# E. Chained lift rate law
# ---------------------------------------------------------------------------
def part_E_chain():
    P("=" * 72)
    P("E: chain m0->M, cheap corrections only: final precision ~ M+const, not 2M")
    m0 = 8
    for M in [20, 26, 32, 40]:
        for K, depth in [(31, 1), (31, 4), (32, 1), (32, 4)]:
            W = 2 * M + 20
            MOD = mpz(1) << W
            p0 = P0(K, MOD)
            inv_p0 = pow(p0, -1, MOD)
            cur = B(m0, K, MOD)
            for mm in range(m0, M):
                Cc = C_cheap(K, mpz(1) << mm, depth) % MOD
                corr = (1 + (mpz(1) << (2 * mm)) * Cc) % MOD
                cur = cur * cur % MOD * inv_p0 % MOD * pow(corr, -1, MOD) % MOD
            exact = B(M, K, MOD)
            vf = v2m(cur - exact, W)
            P(f"  m0={m0} M={M} K={K} depth={depth}:  final_v2={vf:3d}  "
              f"2M={2*M:3d}  M+m0={M+m0:3d}")
    P("  exact-correction control (Theta(K) invs per level == the N^2 cost):")
    M, K = 26, 31
    W = 2 * M + 20
    MOD = mpz(1) << W
    p0 = P0(K, MOD)
    inv_p0 = pow(p0, -1, MOD)
    cur = B(m0, K, MOD)
    for mm in range(m0, M):
        R = mpz(1) << mm
        Cx = sum(pow((8 * k + 1) * (8 * k + 1 + 2 * R) % MOD, -1, MOD)
                 for k in range(1, K + 1)) % MOD
        corr = (1 + (mpz(1) << (2 * mm)) * Cx) % MOD
        cur = cur * cur % MOD * inv_p0 % MOD * pow(corr, -1, MOD) % MOD
    P(f"  exact-C chain m0={m0}->M={M}, K={K}: final_v2={v2m(cur - B(M,K,MOD), W)}"
      f"  (>= 2M+ means exact)")


# ---------------------------------------------------------------------------
# F. Wider monomial relation hunt (MITM)
# ---------------------------------------------------------------------------
def part_F_relation_hunt():
    P("=" * 72)
    P("F: MITM relation hunt, 7 atoms, exponents in [-4,4]")
    P("   prev hunt (dual_jump B): 5 atoms, [-2,2], width 2m.  Now width 3m.")
    EXP = 4
    for m, K in [(24, 24), (28, 20)]:
        for w, expect in [(2 * m - 4, "lemma family MUST appear"),
                          (3 * m, "MUST be empty")]:
            MOD = mpz(1) << w
            atoms = {
                "P0": P0(K) % MOD,
                "Bm": B(m, K) % MOD,
                "Bm1": B(m + 1, K) % MOD,
                "Bm2": B(m + 2, K) % MOD,
                "B2m": B(2 * m, K) % MOD,
                "Bm_2K": B(m, 2 * K) % MOD,
                "Bm_hK": B(m, K // 2) % MOD,
            }
            names = list(atoms.keys())
            # precompute powers e in [-EXP, EXP]
            pw = {}
            for nm in names:
                v = atoms[nm]
                iv = pow(v, -1, MOD)
                d = {0: mpz(1)}
                for e in range(1, EXP + 1):
                    d[e] = d[e - 1] * v % MOD
                    d[-e] = d[-(e - 1)] * iv % MOD
                pw[nm] = d
            left_names, right_names = names[:4], names[4:]
            table = {}
            rng = range(-EXP, EXP + 1)
            for e0 in rng:
                a0 = pw[left_names[0]][e0]
                for e1 in rng:
                    a1 = a0 * pw[left_names[1]][e1] % MOD
                    for e2 in rng:
                        a2 = a1 * pw[left_names[2]][e2] % MOD
                        for e3 in rng:
                            key = int(a2 * pw[left_names[3]][e3] % MOD)
                            table.setdefault(key, []).append((e0, e1, e2, e3))
            hits = []
            for f0 in rng:
                b0 = pw[right_names[0]][-f0]
                for f1 in rng:
                    b1 = b0 * pw[right_names[1]][-f1] % MOD
                    for f2 in rng:
                        key = int(b1 * pw[right_names[2]][-f2] % MOD)
                        if key in table:
                            for le in table[key]:
                                vec = le + (f0, f1, f2)
                                if any(vec):
                                    hits.append(vec)
            P(f"  (m,K)=({m},{K}) w={w} [{expect}]: {len(hits)} nontrivial hits")
            # print up to 6, primitive-looking ones first
            shown = 0
            for vec in sorted(hits, key=lambda v: sum(abs(x) for x in v)):
                if shown >= 6:
                    break
                rel = "*".join(f"{n}^{e}" for n, e in zip(names, vec) if e)
                P(f"    {rel}")
                shown += 1


# ---------------------------------------------------------------------------
# G. e_{T-1} smoothness, larger K
# ---------------------------------------------------------------------------
def is_probable_prime(n, rounds=20):
    n = int(n)
    if n < 2:
        return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        if n % p == 0:
            return n == p
    d = n - 1
    r = 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def part_G_smooth():
    P("=" * 72)
    P("G: e_{T-1} smoothness, K to 64, strip to 1e6, Miller-Rabin cofactor")
    for K in [30, 40, 50, 64]:
        p0 = P0(K)
        e = mpz(0)
        for i in range(1, K + 1):
            e += p0 // (8 * i + 1)
        rem = int(e)
        for p in range(2, 8 * K + 2):
            while rem % p == 0:
                rem //= p
        bits_8K = rem.bit_length()
        p = 8 * K + 2
        while p < 1000000 and rem > 1:
            while rem % p == 0:
                rem //= p
            p += 1
        P(f"  K={K}: bits(e)={e.bit_length()}  cofactor after 8K+1: {bits_8K} bits"
          f"  after 1e6: {rem.bit_length()} bits"
          f"  prime={is_probable_prime(rem) if rem > 1 else '-'}")
    P("  (cofactor staying huge after 1e6 => not poly(K)-smooth either;")
    P("   Borwein-style exponent assembly stays dead beyond the 8K bound.)")


# ---------------------------------------------------------------------------
# H. Paper sanity at larger sizes
# ---------------------------------------------------------------------------
def part_H_paper():
    P("=" * 72)
    P("H: Prop E-equiv (ii) B===P0(1+2^m H) mod 2^{2m}; (iii) H=sum(-8)^j S_j")
    for m, K in [(64, 64), (128, 101), (256, 300)]:
        MOD = mpz(1) << (2 * m)
        ok = B(m, K, MOD) == P0(K, MOD) * (1 + (mpz(1) << m) * atom_H(K, m)) % MOD
        P(f"  (ii)  m={m} K={K}: {ok}")
    for m, K in [(48, 50), (96, 200)]:
        MOD = mpz(1) << m
        H = atom_H(K, m)
        acc = mpz(0)
        last = None
        for j in range(m // 3 + 8):
            acc = (acc + (mpz(-8) ** j % MOD) * (S_j(K, j) % MOD)) % MOD
            if acc == H and last is None:
                last = j
        P(f"  (iii) m={m} K={K}: converged_at_j={last} (m/3={m//3})  "
          f"final_match={acc == H}")


def main():
    part_A_scale()
    part_B_sharp()
    part_C_jump()
    part_D_repair()
    part_E_chain()
    part_F_relation_hunt()
    part_G_smooth()
    part_H_paper()


if __name__ == "__main__":
    main()
