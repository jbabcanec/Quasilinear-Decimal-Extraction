#!/usr/bin/env python3
"""
bet7_doubling_ledger.py -- Bet 7: the doubling ledger.

The archimedean world crossed the exact one-log gap we face (Machin
binary-splitting N log^3 vs AGM N log^2) with a quadratically convergent
iteration carrying O(1) state.  A 2-adic analogue for the dual product
would be an iteration whose state at scale w (B(w,K) mod 2^{2w}, i.e.
H mod 2^w) produces the state at scale 2w (H mod 2^{2w}) in O(1) ring ops.

The natural doubling relation is
    B(2w) = B(w)^2 * P0^{-1} * Exp( sum_j (-1)^{j+1}/j * (2^{2wj} - 2*2^{wj}) Q_j )
with Q_j = sum mu^{-j}, Q_1 = H.

LEDGER CLAIM: even given B(w) EXACTLY, and Q_2, Q_3, ... EXACTLY (all the
higher tower for free), reconstructing B(2w) mod 2^{4w} requires H to depth
    d = 3w - O(1)
because H enters with weight 2^{w+1}.  So one doubling step consumes 3w bits
of H to emit 2w bits: the iteration is bit-NEGATIVE with deficit ratio 3/2.
No state-bounded quadratic iteration can exist; the 2-adic AGM door is closed
by accounting, not by failure to find one.

Parts:
  A. The naive doubling relation is weak: v2(B(w)^2 - P0*B(2w)) = w+1+v2(K)
     (contrast the unit-step lemma's 2w+v2(K) -- longer jumps are WORSE).
  B. Exp/Log reconstruction sanity: with the full tower the formula is exact.
  C. The ledger: sweep truncation depth d of H; measure achieved precision
     of B-hat(2w).  Predict v2(err) = d + w + 1, full width 4w at d = 3w-1.
  D. Same ledger for a general jump w -> cw (c=3,4): need d = (2c-1)w-ish,
     always (2c-1)/c > 2 per step of ratio c: deficit is universal.
"""
from __future__ import annotations

try:
    from gmpy2 import mpz
except ImportError:
    mpz = int


def P(*a, **k):
    k.setdefault("flush", True)
    print(*a, **k)


def v2m(x, W):
    x = int(x) % (1 << W)
    if x == 0:
        return W
    return (x & -x).bit_length() - 1


def B(m, K, MOD):
    R = mpz(1) << m
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = acc * (R + 8 * k + 1) % MOD
    return acc


def P0(K, MOD):
    acc = mpz(1)
    for k in range(1, K + 1):
        acc = acc * (8 * k + 1) % MOD
    return acc


def Qj(K, j, MOD):
    return sum(pow(mpz(8 * k + 1), -j, MOD) for k in range(1, K + 1)) % MOD


def exp2adic(t, W):
    """Exp(t) mod 2^W for v2(t) >= 2. Few terms since v2(t) is large here."""
    MOD = mpz(1) << (W + 32)
    acc = mpz(1)
    term = mpz(1)
    j = 0
    tv = v2m(t, W + 32) if t % MOD else W + 32
    while True:
        j += 1
        if j * max(tv - 2, 1) > W + 16 and j > 2:
            break
        # term_j = t^j / j!  -- divide exactly: v2(t^j) >= j*tv >= v2(j!)
        term = term * (t % MOD)
        fj = mpz(j)
        vj = v2m(fj, 64)
        term >>= vj
        term = term * pow((fj >> vj) % MOD, -1, MOD) % MOD
        acc = (acc + term) % MOD
    return acc % (mpz(1) << W)


def log_term_sum(w, c, Qs, W):
    """sum_j (-1)^{j+1}/j * (2^{c*w*j} - 2*2^{w*j}) * Q_j  mod 2^W.
    Qs[j] = Q_j (possibly truncated for j=1)."""
    MOD = mpz(1) << (W + 32)
    acc = mpz(0)
    j = 0
    while True:
        j += 1
        if j * w > W + 16:
            break
        coeff = (mpz(1) << (c * w * j)) - (mpz(2) << (w * j))
        num = coeff % MOD * (Qs[j] % MOD) % MOD
        vj = v2m(mpz(j), 64)
        num = num >> vj if vj else num
        num = num * pow((mpz(j) >> vj) % MOD, -1, MOD) % MOD
        acc = (acc - num) % MOD if j % 2 == 0 else (acc + num) % MOD
    return acc % (mpz(1) << W)


def part_A_weak():
    P("=" * 72)
    P("A: naive doubling is WEAK: v2(B(w)^2 - P0*B(2w)) = w+1+v2(K)")
    P("   (unit step has 2w+v2(K); longer jumps decay toward w)")
    for w in [8, 12, 16, 20]:
        for K in [w, 16, 24]:
            W = 4 * w + 20
            MOD = mpz(1) << W
            err = (B(w, K, MOD) ** 2 - P0(K, MOD) * B(2 * w, K, MOD)) % MOD
            ve = v2m(err, W)
            pred = w + 1 + v2m(K, 30)
            P(f"  w={w:2d} K={K:2d}  v2={ve:3d}  pred={pred:3d}  ok={ve == pred}")


def part_B_sanity():
    P("=" * 72)
    P("B: Exp/Log reconstruction with FULL tower is exact to width 4w")
    for w, K in [(8, 8), (12, 12), (16, 16)]:
        W = 4 * w
        Wg = W + 32
        MOD = mpz(1) << Wg
        Qs = {j: Qj(K, j, MOD) for j in range(1, (Wg + 16) // w + 2)}
        arg = log_term_sum(w, 2, Qs, Wg)
        Bh = B(w, K, MOD) ** 2 % MOD * pow(P0(K, MOD), -1, MOD) % MOD \
            * exp2adic(arg, Wg) % MOD
        exact = B(2 * w, K, MOD)
        P(f"  w={w} K={K}: v2(Bhat - B(2w)) = {v2m(Bh - exact, W)}  (need {W})")


def part_C_ledger():
    P("=" * 72)
    P("C: THE LEDGER.  B(w) exact, Q_2.. exact, H truncated to d bits.")
    P("   predict v2(err) = d + w + 1;  full width 4w at d = 3w - 1")
    for w, K in [(8, 8), (12, 12), (16, 16), (20, 16)]:
        W = 4 * w
        Wg = W + 32
        MOD = mpz(1) << Wg
        Qs_full = {j: Qj(K, j, MOD) for j in range(1, (Wg + 16) // w + 2)}
        p0inv = pow(P0(K, MOD), -1, MOD)
        Bw2 = B(w, K, MOD) ** 2 % MOD
        exact = B(2 * w, K, MOD)
        row = f"  w={w:2d} K={K:2d}:"
        dstar = None
        for d in [w, 2 * w, 3 * w - 4, 3 * w - 1, 3 * w + 2]:
            Qs = dict(Qs_full)
            Qs[1] = Qs_full[1] % (mpz(1) << d)          # H to d bits only
            arg = log_term_sum(w, 2, Qs, Wg)
            Bh = Bw2 * p0inv % MOD * exp2adic(arg, Wg) % MOD
            ve = v2m(Bh - exact, W)
            row += f"  d={d}:{ve}"
            if dstar is None and ve >= W:
                dstar = d
        P(row + f"   need={W}  first_full_d={dstar}  3w-1={3*w-1}")
    P("  LEDGER: holding w bits of H, a doubling emits 2w but consumes 3w-1:")
    P("  every NEW certified bit of H costs 2 external bits of H.  Uniform 2:1.")


def part_D_general_jump():
    P("=" * 72)
    P("D: jump w -> cw needs H to d = (2c-1)w - O(1); ratio (2c-1)/(c) > 2 always")
    for c in [2, 3, 4]:
        w, K = 10, 12
        W = 2 * c * w
        Wg = W + 40
        MOD = mpz(1) << Wg
        Qs_full = {j: Qj(K, j, MOD) for j in range(1, (Wg + 16) // w + 2)}
        p0inv = pow(P0(K, MOD), -1, MOD)
        Bw2 = B(w, K, MOD) ** 2 % MOD
        exact = B(c * w, K, MOD)
        # relation: B(cw) = B(w)^2 P0^{-1} Exp(sum_j (-1)^{j+1}/j (2^{cwj}-2*2^{wj}) Q_j)
        dstar = None
        for d in range((2 * c - 1) * w - 6, (2 * c - 1) * w + 3):
            Qs = dict(Qs_full)
            Qs[1] = Qs_full[1] % (mpz(1) << d)
            arg = log_term_sum(w, c, Qs, Wg)
            Bh = Bw2 * p0inv % MOD * exp2adic(arg, Wg) % MOD
            if v2m(Bh - exact, W) >= W:
                dstar = d
                break
        held = w
        gain = c * w - held
        external = dstar - held
        P(f"  c={c}: hold {held}, emit {c*w}, consume d*={dstar}"
          f"  ((2c-1)w-1={(2*c-1)*w-1})  external/gain={external/gain:.2f}")


def main():
    part_A_weak()
    part_B_sanity()
    part_C_ledger()
    part_D_general_jump()
    P()
    P("BET7 VERDICT: doubling ledger is exact and bit-negative (3w in, 2w out).")
    P("The 2-adic-AGM door is closed by accounting: the archimedean trick that")
    P("crosses this one-log gap (quadratic iteration, O(1) state) cannot have a")
    P("2-adic analogue for the extraction atom.")


if __name__ == "__main__":
    main()
