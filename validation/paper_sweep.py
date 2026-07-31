#!/usr/bin/env python3
"""paper_sweep.py -- mechanical audit of the manuscript (2026-07-31)

Written after a full claim-by-claim read turned up a factor-of-two error, a
transposed matrix, a bound that contradicted the theorem it cited, and a
paragraph silently corrupted by shell escape collapsing. Each of those has a
mechanical signature; this checks for all of them so the class does not
recur.
"""

import io
import re
import sys

TEX = r"C:\Users\josep\dropbox\babcanec works\mathematics\bbp algos\decimal_extraction_paper.tex"

BS = chr(92)          # backslash, kept out of literals to survive any shell


def check(name, ok, detail=""):
    print("  [%s] %s%s" % ("ok " if ok else "FAIL", name,
                           ("  -- " + detail) if detail and not ok else ""))
    return ok


def main():
    s = io.open(TEX, encoding="utf-8", newline="").read()
    allok = True

    # 1. control characters from collapsed escapes (\t, \r, \f, \v, \b, \a)
    stray = [i for i, c in enumerate(s)
             if c in "\t\x0b\x0c\x08\x07"
             or (c == "\r" and s[i:i + 2] != "\r\n")]
    allok &= check("no stray control chars (collapsed escapes)", not stray,
                   "%d found" % len(stray))
    for i in stray[:5]:
        print("        near: %r" % s[max(0, i - 50):i + 15])

    # 2. big operators must carry an index
    for op in ("prod", "sum"):
        bare = re.findall(BS + BS + op + r"(?![_^{a-zA-Z])", s)
        allok &= check("no bare " + BS + op, not bare, "%d bare" % len(bare))

    # 3. the 2-adic logarithm must not be written with the binary-log macro.
    #    \log_2 is legitimate for bit counts; flag it only next to 2-adic words.
    for m in re.finditer(BS + BS + r"log_2", s):
        ctx = s[max(0, m.start() - 260):m.start() + 260]
        if "2$-adic" in ctx or "2-adic" in ctx:
            allok &= check("bare " + BS + "log_2 near 2-adic text", False,
                           repr(s[m.start() - 60:m.start() + 60]))

    # 4. discrete-log base used consistently (code + validators use base 9)
    b3 = len(re.findall(r"dlog}_3|dlog_3", s))
    b9 = len(re.findall(r"dlog}_9|dlog_9", s))
    l3 = len(re.findall(r"Log}_2 3|Log}_2 9", s))
    allok &= check("dlog base 9 throughout (matches borwein_crack.py)",
                   b3 == 0 and b9 > 0, "base-3 sites: %d" % b3)
    print("        (dlog_9 sites %d, Log_2 base sites %d)" % (b9, l3))

    # 5. orbit-matrix orientation must agree between the two subsections
    tn = "T{" + BS + "times} N" in s or r"T\times N" in s
    nt = r"N\times T" in s
    allok &= check("orbit matrix orientation stated once (T x N)",
                   tn and not nt, "T x N: %s, N x T: %s" % (tn, nt))

    # 6. arithmetic constants that appear as decimals in the text
    import math
    checks = [
        ("0.2902 bytes per digit of 5^n", 0.2902, math.log2(5) / 8, 1e-4),
        ("3.32 N terms  = 4(K+1), K ~ n log2(10)/4", 3.32, math.log2(10), 5e-3),
        ("N^1.077 interpreted exponent plausible", 1.077,
         math.log(111.6 / 0.06) / math.log(10 ** 7 / 10 ** 4), 0.02),
    ]
    for name, stated, actual, tol in checks:
        allok &= check(name, abs(stated - actual) <= tol,
                       "stated %.4f vs %.4f" % (stated, actual))

    # 7. every \ref and \cite target exists
    labels = set(re.findall(r"label\{([^}]+)\}", s))
    refs = set(re.findall(r"ref\{([^}]+)\}", s))
    bibs = set(re.findall(r"bibitem\{([^}]+)\}", s))
    cites = set()
    for c in re.findall(r"cite\{([^}]+)\}", s):
        cites.update(x.strip() for x in c.split(","))
    allok &= check("all " + BS + "ref targets defined", refs <= labels,
                   str(sorted(refs - labels)))
    allok &= check("all " + BS + "cite keys defined", cites <= bibs,
                   str(sorted(cites - bibs)))
    unused = bibs - cites
    print("        (%d bib entries, %d uncited: %s)"
          % (len(bibs), len(unused), sorted(unused) if unused else "none"))

    print()
    print("SWEEP", "PASSED" if allok else "FAILED")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
