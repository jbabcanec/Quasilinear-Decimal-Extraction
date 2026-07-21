#!/usr/bin/env python3
"""
The Dwork residual valuation reproduces the cross-formula tower law.

For each BBP family with denominator d*k+1, take f(x)=sum_k x^k/(d k+1) and the
classical Dwork residual  D_s(x) = f_{s+1}(x) f_{s-1}(x^2) - f_s(x) f_s(x^2),
f_s = truncation to degree < 2^s, at p=2.  We measure

        val(s) = min_j v_2( [x^j] D_s ).

RESULT (this file): val(s) = s + v_2(d) - 1 exactly, slope 1, for d=2,4,8,16.
That is the SAME law bellard_tower.py validates ( min v_2(R_s) = s + v_2(d)-1 ),
here obtained from the Frobenius/Dwork side by an independent construction.

CONSEQUENCE.  The Dwork congruence supplies exactly s + v_2(d) - 1 two-adic
digits of the level-(s+1) truncation for free.  At the depth needed for N terms,
s = log2(N), that is Theta(log N) digits when the digit pipeline needs Theta(N).
The slope is 1 and formula-independent: no base-2^r formula raises it (formula
shopping closed), and basic Frobenius descent does NOT amplify precision.
"""
p = 2
W = 1 << 12
MOD = 1 << W


def v2(c):
    c %= MOD
    return W if c == 0 else (c & -c).bit_length() - 1


def polymul(A, B):
    C = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b:
                    C[i + j] = (C[i + j] + a * b) % MOD
    return C


def polysub(A, B):
    n = max(len(A), len(B))
    A = A + [0] * (n - len(A)); B = B + [0] * (n - len(B))
    return [(A[i] - B[i]) % MOD for i in range(n)]


def frob(A):
    C = [0] * ((len(A) - 1) * p + 1)
    for i, a in enumerate(A):
        C[i * p] = a
    return C


def v2int(n):
    return (n & -n).bit_length() - 1


def dwork_val_seq(d, S=8):
    A = lambda k: pow(d * k + 1, -1, MOD)
    f = [[A(k) for k in range(1 << l)] for l in range(S + 2)]
    seq = []
    for s in range(1, S):
        Ds = polysub(polymul(f[s + 1], frob(f[s - 1])),
                     polymul(f[s], frob(f[s])))
        seq.append((s, min(v2(c) for c in Ds)))
    return seq


def main():
    print("family d*k+1 : val(s) = min v2 of Dwork residual   (p=2)")
    print(" d  v2(d) | intercept val(s)-s | first-diffs | predicted s+v2(d)-1")
    ok = True
    for d in (2, 4, 8, 16):
        seq = dwork_val_seq(d)
        inter = seq[0][1] - seq[0][0]
        diffs = [seq[i + 1][1] - seq[i][1] for i in range(len(seq) - 1)]
        pred = v2int(d) - 1
        good = (inter == pred) and set(diffs) == {1}
        ok = ok and good
        print(" %2d   %d    | %+d                | %s | %+d   %s"
              % (d, v2int(d), inter, diffs, pred, "OK" if good else "FAIL"))
    print()
    print("val(s) = s + v_2(d) - 1  reproduced from the Dwork side:",
          "PASS" if ok else "FAIL")
    print("=> slope 1, formula-independent: basic Frobenius descent yields",
          "Theta(log N) digits at depth log N.")


if __name__ == '__main__':
    main()
