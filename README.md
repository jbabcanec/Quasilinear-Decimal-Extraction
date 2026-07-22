# Quasi-linear decimal digit extraction for π

Compute the N-th **decimal** digit of π **without computing any earlier digit**, in
quasi-linear time: O(N log³N) bit operations, measured at N^1.077 across positions
10⁴–10⁷. To our knowledge this is the first implemented and verified subquadratic
decimal digit-extraction algorithm.

Every published decimal extractor before this (Plouffe 1996, Bellard 1997, Gourdon
2003, Plouffe 2022) runs in essentially quadratic time. The method multiplies the
classic *binary* BBP formula termwise by 10^{N−1} and handles the terms that defeated
prior attempts (including Zudilin's 2024 note, arXiv:2409.10097) through one shared,
read-only object: the binary expansion of 5^{N−1}, about 0.29·N bytes.

## Results

Single-threaded Python over GMP, default engine:

| N | time | note |
|---|---|---|
| 10⁴ | 0.06 s | |
| 10⁶ | 8.9 s | |
| 10⁷ | 111.6 s | digits `725915133612`, verified against public corpora |

The fitted exponent is N^1.077 over three decades with no upward bend. Digits are
verified against an independent Chudnovsky computation (through 10⁶) and against the
pi.delivery and angio.net record corpora (at 10⁷).

## Quickstart

```bash
pip install gmpy2          # optional but recommended (GMP integers)

# the millionth decimal digit of pi (the v3 engine is the default)
python decimal_bbp_extract.py --pos 1000000

# ln 2 also works; engines v1/v2/v3/v4 are selectable
python decimal_bbp_extract.py --const ln2 --pos 100000

# verification and benchmarks
python decimal_bbp_extract.py --selftest   # digits vs an independent Chudnovsky reference
python decimal_bbp_extract.py --v3-check   # the division-free pipeline vs the reference
python decimal_bbp_extract.py --beatquad   # scaling benchmark
```

## Contents

- **`decimal_extraction_paper.tex` / `.pdf`** — the paper: the algorithm, the
  O(N log³N) theorem, the sieve-batched exponentiation and 2-adic normal form, the
  generalizations (any base-2^r BBP constant, any output base), and an analysis of the
  one remaining logarithmic factor as an open problem.
- **`decimal_bbp_extract.py`** — the reference implementation, with four engines
  (v1 naive band reduction; v2 accumulating-remainder-tree batching; v3 the
  division-free 2-adic pipeline, used by default; v4 sieve-batched word phase), plus
  self-tests against integer-Chudnovsky references and scaling benchmarks.
- **`three_phase_extract.cpp`** — a dependency-free C++/OpenMP implementation of the
  band identity, used as an independent cross-check.
- **`validation/`** — standalone exact-arithmetic validators for every identity in
  the paper.
- **`legacy/`** — earlier approaches that this work supersedes.

## The open problem

The algorithm is one logarithmic factor above the O(N log²N) cost of computing all N
digits by the arithmetic–geometric mean. Section 6 of the paper isolates that factor
as a single problem: computing a consolidated 2-adic integer V, equivalent to an
isolated harmonic number modulo 2^N. The multiplicative part of V is computable in
O(N log²N); the additive part has no known fast algorithm. Solving it would give
decimal digit extraction in O(N log²N).

## Reference

J. Babčanec, *Extracting decimal digits of π in quasi-linear time from binary BBP
formulas*, 2026.
