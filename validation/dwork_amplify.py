#!/usr/bin/env python3
"""
Dwork precision-amplification probe for the additive core.

Object: f(x) = 2F1(1/8, 1; 9/8; x) = sum_{k>=0} x^k / (8k+1),  p = 2.
Truncations: f_s(x) = sum_{k < 2^s} A_k x^k,  A_k = (8k+1)^{-1} in Z_2.

THE QUESTION.  Dwork's congruence relates consecutive truncation levels through
Frobenius x -> x^p.  Define the classical residual
        D_s(x) = f_{s+1}(x) f_{s-1}(x^p) - f_s(x) f_s(x^p).
Its minimum 2-adic coefficient valuation  val(s) = min_j v_2([x^j] D_s)  is EXACTLY
the number of 2-adic digits of the level-(s+1) truncation that the lower levels
predict for free.  We also track the Dwork-unit ratio rho_s = f_{s+1}/f_s(x^p)
and its stabilization v_2(rho_s - rho_{s-1}).

If val(s) = Theta(s): to reach N terms (depth s = log2 N) we get only O(log N)
digits -> the barrier is real.  If val(s) grows super-linearly (e.g. ~2^s):
Frobenius descent amplifies precision -> the crack.
"""
import sys

p = 2
W = 1 << 12                     # 2-adic working precision (bits); plenty
MOD = 1 << W
INF = W                         # valuation cap


def v2(c):
    c %= MOD
    if c == 0:
        return INF
    return (c & -c).bit_length() - 1


def polymul(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if not a:
            continue
        for j, b in enumerate(B):
            if b:
                C[i + j] = (C[i + j] + a * b) % MOD
    return C


def polysub(A, B):
    n = max(len(A), len(B))
    A = A + [0] * (n - len(A))
    B = B + [0] * (n - len(B))
    return [(A[i] - B[i]) % MOD for i in range(n)]


def frob(A):                    # A(x^p)
    C = [0] * ((len(A) - 1) * p + 1)
    for i, a in enumerate(A):
        C[i * p] = a
    return C


def inv_series(A, D):           # 1/A mod x^D, A[0] a 2-adic unit
    inv = [0] * D
    inv[0] = pow(A[0] % MOD, -1, MOD)
    for n in range(1, D):
        s = 0
        for j in range(1, n + 1):
            if j < len(A) and A[j]:
                s = (s + A[j] * inv[n - j]) % MOD
        inv[n] = (-inv[0] * s) % MOD
    return inv


def series_mul_trunc(A, B, D):
    C = [0] * D
    for i in range(min(len(A), D)):
        if not A[i]:
            continue
        for j in range(min(len(B), D - i)):
            if B[j]:
                C[i + j] = (C[i + j] + A[i] * B[j]) % MOD
    return C


def A_of(k):
    return pow(8 * k + 1, -1, MOD)


def main():
    S = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    # truncations f_0 .. f_{S+1}
    f = []
    for l in range(S + 2):
        f.append([A_of(k) for k in range(1 << l)])

    print("f(x) = 2F1(1/8,1;9/8;x) = sum x^k/(8k+1),  p=2,  precision 2^%d" % W)
    print()
    print(" s | terms=2^(s+1) | val(s)=min v2 of Dwork residual | rho stabilization")
    print("---+---------------+---------------------------------+------------------")
    quad_vals = []
    for s in range(1, S):
        # quadratic Dwork residual  D_s = f_{s+1} f_{s-1}(x^p) - f_s f_s(x^p)
        left = polymul(f[s + 1], frob(f[s - 1]))
        right = polymul(f[s], frob(f[s]))
        Ds = polysub(left, right)
        val = min(v2(c) for c in Ds)
        quad_vals.append((s, val))

        # Dwork-unit ratio stabilization rho_s = f_{s+1} / f_s(x^p), compare to rho_{s-1}
        Dout = 1 << s
        rho_s = series_mul_trunc(f[s + 1], inv_series(frob(f[s]), Dout), Dout)
        rho_sm = series_mul_trunc(f[s], inv_series(frob(f[s - 1]), Dout), Dout)
        diff = polysub(rho_s, rho_sm)
        rstab = min(v2(c) for c in diff[:Dout])

        print(" %d | %13d | %31d | %d"
              % (s, 1 << (s + 1), val, rstab))

    # slope of val(s) vs s (linear fit) and vs 2^s
    print()
    xs = [s for s, _ in quad_vals if quad_vals]
    ys = [v for _, v in quad_vals]
    if len(xs) >= 3:
        n = len(xs)
        sx, sy = sum(xs), sum(ys)
        sxx = sum(x * x for x in xs)
        sxy = sum(x * y for x, y in zip(xs, ys))
        slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
        inter = (sy - slope * sx) / n
        print("linear fit  val(s) ~ %.3f * s + %.3f" % (slope, inter))
        # ratios val(s)/s and val(s)/2^s to see the growth law
        print("val(s)/s   :", ["%.2f" % (v / s) for s, v in quad_vals])
        print("val(s+1)-val(s):", [ys[i+1]-ys[i] for i in range(len(ys)-1)])
        print()
        print("READING: constant first-difference => val(s) linear in s => O(log N) digits")
        print("         at depth log N => BARRIER.  Growing first-difference => amplification.")


if __name__ == '__main__':
    main()
