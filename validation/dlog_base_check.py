#!/usr/bin/env python3
"""dlog_base_check.py -- is E = Log_2(P)/Log_2(3) or /Log_2(9)?  (2026-07-31)

The paper defines (Sec. 6.2)

    prod_{k<=K}(1+8k) = +-3^E  (mod 2^m),   E = sum_k dlog_3(1+8k),

and states dlog_3(u) = log_2(u)/log_2(3). But the proof sketch of the crack
theorem writes E = Log_2(prod)/Log_2(9). Those differ by a factor 2, since
Log_2(9) = 2 Log_2(3). This decides it by exact 2-adic arithmetic.

Note every 1+8k is 1 mod 8, hence a square in the unit group: the units that
are 1 mod 8 are exactly <9> = <3^2>, so dlog_3 of each factor is EVEN and E is
even. That is consistent with either reading, so the test is whether
3^E = P or 9^E = P.
"""


def v2(n):
    return (n & -n).bit_length() - 1 if n else 10 ** 9


def log2adic(u, m):
    """2-adic Log of a unit u = 1 (mod 4), as an integer mod 2^m."""
    M = 1 << m
    x = (u - 1) % M
    vx = v2(x)
    assert vx >= 2, "need u = 1 mod 4"
    total, n = 0, 1
    while True:
        a = v2(n)
        b = n >> a
        xn = pow(x, n, 1 << (m + a + 2))
        t = (xn >> a) * pow(b, -1, M) % M
        total = (total + t) % M if n % 2 else (total - t) % M
        if vx * (n + 1) - v2(n + 1) >= m:
            break
        n += 1
    return total


def main():
    m, K = 96, 40
    M = 1 << m

    P = 1
    for k in range(K + 1):
        P = P * (8 * k + 1) % M
    assert P % 8 == 1

    LP = log2adic(P, m)
    L9 = log2adic(9, m)
    L3 = log2adic((-3) % M, m)          # Log(3) = Log(-1) + Log(-3) = Log(-3)

    print("K = %d, m = %d" % (K, m))
    print("  v2(Log P) = %d,  v2(Log 3) = %d,  v2(Log 9) = %d"
          % (v2(LP), v2(L3), v2(L9)))
    print("  Log 9 == 2 * Log 3 ?  %s" % (L9 == 2 * L3 % M))

    # E candidates, taken as exact 2-adic quotients
    def quot(a, b, prec):
        va, vb = v2(a), v2(b)
        assert va >= vb, "quotient not integral"
        return (a >> vb) * pow(b >> vb, -1, 1 << prec) % (1 << prec)

    prec = m - 4
    E3 = quot(LP, L3, prec)
    E9 = quot(LP, L9, prec)
    print("  E via /Log 3 = %d" % E3)
    print("  E via /Log 9 = %d   (= E3/2 ? %s)" % (E9, E9 == E3 // 2))

    Mt = 1 << prec
    ok3 = pow(3, E3, Mt) == P % Mt
    ok9 = pow(3, E9, Mt) == P % Mt
    print()
    print("  3^(Log P / Log 3) == P ?   %s" % ok3)
    print("  3^(Log P / Log 9) == P ?   %s" % ok9)
    print("  9^(Log P / Log 9) == P ?   %s" % (pow(9, E9, Mt) == P % Mt))
    print()
    if ok3 and not ok9:
        print("VERDICT: with E = sum dlog_3 as the paper defines it, the correct")
        print("normalization is E = Log_2(P)/Log_2(3). The proof sketch's")
        print("'/Log_2 9' computes E/2 and is off by a factor of two.")
    elif ok9 and not ok3:
        print("VERDICT: '/Log_2 9' is right and the definition needs dlog_9.")
    else:
        print("VERDICT: inconclusive, investigate.")


if __name__ == "__main__":
    main()
