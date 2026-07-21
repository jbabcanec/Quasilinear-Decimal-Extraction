# Quasi-Linear Decimal Digit Extraction for π

The N-th **decimal** digit of π, computed **without computing any earlier digit**, in
quasi-linear time — the first subquadratic decimal digit-extraction algorithm in any
model, measured at **N^1.077** across positions 10⁴–10⁷.

Every published decimal extractor before this (Plouffe 1996, Bellard 1997, Gourdon
2003, Plouffe 2022) is quadratic-time or worse. The method here multiplies the classic
*binary* BBP formula termwise by 10^{N−1} and resolves the terms that defeated all
prior attempts (including Zudilin's 2024 note, arXiv:2409.10097) through one shared
read-only object: the binary expansion of 5^{N−1}, ~0.29·N bytes.

## Headline numbers (single-threaded Python + GMP, `--engine v3`)

| N | time | note |
|---|---|---|
| 10⁴ | 0.06 s | |
| 10⁶ | 8.9 s | vs 571 s (naive Python), 32.3 s (naive C++, 16 threads) |
| 10⁷ | 111.6 s | digits `725915133612`, verified against public corpora |

Fitted exponent N^1.077 over three decades, no bend. Digits verified against
independent Chudnovsky computation (≤10⁶) and the pi.delivery / angio.net record
corpora (10⁷).

## Quickstart

```bash
pip install gmpy2          # optional but strongly recommended (GMP integers)

# millionth digit of pi (v3 engine is the default)
python decimal_bbp_extract.py --pos 1000000

# ln 2 works too; engines v1/v2/v4 selectable
python decimal_bbp_extract.py --const ln2 --pos 100000
python decimal_bbp_extract.py --pos 1000000 --engine v1

# verification batteries and benchmarks
python decimal_bbp_extract.py --selftest       # digits vs Chudnovsky reference
python decimal_bbp_extract.py --v3-check       # v3 vs v1 battery
python decimal_bbp_extract.py --beatquad       # scaling benchmark (fits exponent)
```

## Contents

- `decimal_extraction_paper.tex` / `.pdf` — **the paper**: the algorithm, the
  O(N log³N) theorem, the sieve-batched exponentiation and 2-adic normal form, the
  generalizations (any base-2^r BBP constant; any output base), and the barrier to
  going faster, stated in five machine-validated equivalent forms with an explicit
  open problem.
- `decimal_bbp_extract.py` — all four engines:
  - **v3 (default)**: Bézout word phase + division-free 2-adic V merge — the
    measured-N^1.08 pipeline. *Use this one until we break it.*
  - v1: naive band reduction (reference implementation of the basic three-phase split)
  - v2: accumulating-remainder-tree batching
  - v4: sieve-batched word phase — validated at 10⁷; slower than v3 under CPython
    (interpreter tax, quantified in the log); serves as the correctness oracle and
    structural template for the C port
  - plus exact per-term auditors, self-tests vs integer-Chudnovsky references, and
    scaling benchmarks
- `three_phase_extract.cpp` — dependency-free C++/OpenMP port of v1 (verified to 10⁷)
- `validation/` — standalone exact validators for every identity in the paper:
  Bézout split & dyadic consolidation, sieve-batched chains (Euler/Fermat),
  2-adic Euler–Maclaurin, and the Γ-transform faces
- `refs/` — archived copies of the load-bearing sources (Gourdon 2003 incl. the
  unpublished Theorem 2, Plouffe 2022, Zudilin 2024, Bailey's survey, citation data)
- `research_log.md` — the complete research record, including every failed route
- `legacy/` — the March 2026 work this supersedes (Machin series-splitting extractor,
  Guillera analysis, ln 2 three-term decomposition, Gourdon comparison builds, and the
  original paper)

## Status / roadmap

- Proven: O(N log³N) bit operations, Θ(N)-bit read-only auxiliary space; all
  non-tree phases O(N loglog N) word operations.
- Measured: N^1.077 (v3). Next speedup: C+GMP port of the v4 structure (~3–7×
  word phase + threading; v4 is the oracle).
- Open: the last logarithm — one precisely-stated problem with five equivalent
  faces (see paper §6); solving it gives O(N log²N).

## Reference

J. Babčanec, *Extracting decimal digits of π in quasi-linear time from binary BBP
formulas*, manuscript, July 2026 (this repository).
