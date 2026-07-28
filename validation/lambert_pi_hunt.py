#!/usr/bin/env python3
"""lambert_pi_hunt.py -- PSLQ hunt: is pi in the rational span of natively
window-cheap constants (Lambert/divisor and theta values at DYADIC arguments)?

Why this basis: for q = 2^-r, the series Sum d(n) q^n, Sum sigma_1(n) q^n,
theta_3(q) = 1 + 2 Sum q^{n^2}, and Euler's phi(q) = Prod(1-q^n) (pentagonal,
lacunary) are all extractable at position N in O(N loglog N) word ops given the
shared string X = 5^(N-1): coefficients come from a sieve, and every term is a
word times a WINDOW READ of X at an arithmetic (or quadratic/pentagonal)
position. No fraction sum, no merge tree, no atom. A rational relation
pi = Sum r_i * (these) would therefore give decimal digit extraction in
O(N log N) bits -- below the campaign prize -- through a door no theorem closes
(BBG 2004 covers Machin-type BBP only; Nesterenko's independence theorem at
algebraic q constrains the E_2/E_4/E_6 values among THEMSELVES, and is silent
on pi).

Expectation (recorded before running): EMPTY. The modular analysis says pi
lives in the quasimodular ANOMALY E_2(-1/tau) = tau^2 E_2(tau) + 12 tau/(2 pi i),
which is only visible from a PAIR of nomes q = e^{2 pi i tau}, q' = e^{-2 pi i/tau};
if q is dyadic the partner q' is transcendental and never dyadic, so the
pi-carrying combination is never expressible inside the dyadic-cheap world.
The hunt tests for an ACCIDENTAL relation anyway, because both prior cracks of
this campaign entered through representations absent from every list.

Protocol: two height-1 positive controls must be found (harness validity);
main vectors reported with the height bound (maxcoeff) and working precision,
so an empty result is a calibrated exclusion: "no relation of height <= H at
precision 10^-P." Any hit is re-verified at doubled precision before being
believed, and even then is only a LEAD (the word stays reserved).
"""

import time
from mpmath import mp, mpf, pi, log, sqrt, pslq

DPS = 400          # hunt precision (digits)
VERIFY_DPS = 800   # re-verification precision for any hit
MAXCOEFF = 10**10  # height bound for the exclusion statement
MAXSTEPS = 500000


def lambert_sigma(a, q):
    """Sum_{n>=1} sigma_a(n) q^n via Lambert form Sum_k k^a q^k/(1-q^k), exact
    to working precision (terms decay like q^k)."""
    K = int(mp.dps * mp.log(10) / -mp.log(q)) + 10
    s = mpf(0)
    for k in range(1, K + 1):
        qk = q ** k
        s += (mpf(k) ** a) * qk / (1 - qk)
    return s


def theta3(q):
    K = int(sqrt(mp.dps * mp.log(10) / -mp.log(q))) + 5
    return 1 + 2 * sum(q ** (n * n) for n in range(1, K + 1))


def theta4(q):
    K = int(sqrt(mp.dps * mp.log(10) / -mp.log(q))) + 5
    return 1 + 2 * sum((-1) ** n * q ** (n * n) for n in range(1, K + 1))


def euler_phi(q):
    K = int(mp.dps * mp.log(10) / -mp.log(q)) + 10
    p = mpf(1)
    for n in range(1, K + 1):
        p *= (1 - q ** n)
    return p


def harmonic_lambert(q):
    """Sum_{k>=1} q^k / k = -log(1-q); computed as the SERIES (control)."""
    K = int(mp.dps * mp.log(10) / -mp.log(q)) + 10
    return sum(q ** k / k for k in range(1, K + 1))


def hunt(name, vec, labels):
    t0 = time.time()
    rel = pslq(vec, maxcoeff=MAXCOEFF, maxsteps=MAXSTEPS)
    dt = time.time() - t0
    if rel is None:
        print("  %-14s EMPTY  (no relation, height <= %.0e, %d digits, %.1fs)"
              % (name, MAXCOEFF, mp.dps, dt))
        return None
    print("  %-14s HIT %s  (%.1fs)  -- re-verifying at %d digits"
          % (name, rel, dt, VERIFY_DPS))
    for c, l in zip(rel, labels):
        if c:
            print("      %+d * %s" % (c, l))
    return rel


def main():
    mp.dps = DPS
    q2, q4, q16 = mpf(1) / 2, mpf(1) / 4, mpf(1) / 16

    print("lambert_pi_hunt: PSLQ at %d digits, height bound %.0e" % (DPS, MAXCOEFF))

    # ---- positive controls (must both HIT at height 1) ----
    print("controls:")
    c1 = pslq([log(2), harmonic_lambert(q2)], maxcoeff=100, maxsteps=100000)
    print("  ln2 = Lambert-harmonic(1/2):", c1, "(expect [1, -1] or [-1, 1])")
    c2 = pslq([theta3(q2) * theta4(q2), theta4(q4) ** 2], maxcoeff=100,
              maxsteps=100000)
    print("  Jacobi t3(q)t4(q) = t4(q^2)^2:", c2, "(expect [1, -1] or [-1, 1])")
    ok_controls = (c1 is not None and sorted(map(abs, c1)) == [1, 1]
                   and c2 is not None and sorted(map(abs, c2)) == [1, 1])
    if not ok_controls:
        print("CONTROLS FAILED -- harness invalid, aborting")
        raise SystemExit(1)
    print("  controls PASSED (harness finds height-1 truths)")

    # ---- the dyadic-cheap basis ----
    basis = [
        (mpf(1), "1"),
        (log(2), "ln2"),
        (lambert_sigma(0, q2), "L_d(1/2)"),
        (lambert_sigma(0, q4), "L_d(1/4)"),
        (lambert_sigma(0, q16), "L_d(1/16)"),
        (lambert_sigma(1, q2), "L_s1(1/2)"),
        (lambert_sigma(1, q4), "L_s1(1/4)"),
        (lambert_sigma(3, q2), "L_s3(1/2)"),
        (theta3(q2), "t3(1/2)"),
        (theta3(q4), "t3(1/4)"),
        (theta4(q2), "t4(1/2)"),
        (euler_phi(q2), "phi(1/2)"),
    ]
    bvals = [b[0] for b in basis]
    blabs = [b[1] for b in basis]

    print("main hunts (target + %d basis constants):" % len(basis))
    hits = []
    for tname, tval in [("pi", +pi), ("1/pi", 1 / pi), ("pi^2", pi ** 2)]:
        r = hunt(tname, [tval] + bvals, [tname] + blabs)
        if r is not None:
            hits.append((tname, r))

    # multiplicative sector: monomial relations among logs of positives
    logvec = [log(pi), log(2), log(lambert_sigma(0, q2)),
              log(lambert_sigma(1, q2)), log(theta3(q2)), log(theta4(q2)),
              log(theta3(q4)), log(euler_phi(q2)), log(lambert_sigma(0, q4))]
    loglab = ["log pi", "ln2", "log L_d(1/2)", "log L_s1(1/2)", "log t3(1/2)",
              "log t4(1/2)", "log t3(1/4)", "log phi(1/2)", "log L_d(1/4)"]
    r = hunt("monomial", logvec, loglab)
    if r is not None:
        hits.append(("monomial", r))

    if not hits:
        print("VERDICT: EMPTY -- pi has no relation of height <= %.0e over the"
              % MAXCOEFF)
        print("dyadic Lambert/theta/Euler basis at %d digits (controls passed)."
              % DPS)
        print("The dyadic-cheap door is CLOSED at this height; the quasimodular-")
        print("anomaly analysis stands (pi lives across a modular transform whose")
        print("partner nome is never dyadic).")
        return 0

    # any hit: re-verify at doubled precision
    print("re-verifying hits at %d digits..." % VERIFY_DPS)
    mp.dps = VERIFY_DPS
    # (recompute everything; a surviving relation must satisfy to ~10^-VERIFY_DPS)
    q2, q4, q16 = mpf(1) / 2, mpf(1) / 4, mpf(1) / 16
    bvals2 = [mpf(1), log(2), lambert_sigma(0, q2), lambert_sigma(0, q4),
              lambert_sigma(0, q16), lambert_sigma(1, q2), lambert_sigma(1, q4),
              lambert_sigma(3, q2), theta3(q2), theta3(q4), theta4(q2),
              euler_phi(q2)]
    tmap = {"pi": +pi, "1/pi": 1 / pi, "pi^2": pi ** 2}
    for tname, rel in hits:
        if tname == "monomial":
            print("  monomial hit: re-verify manually (logs)")
            continue
        v = rel[0] * tmap[tname] + sum(c * b for c, b in zip(rel[1:], bvals2))
        print("  %s: residual = %s (must be ~1e-%d to survive)"
              % (tname, mp.nstr(v, 8), VERIFY_DPS - 20))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
