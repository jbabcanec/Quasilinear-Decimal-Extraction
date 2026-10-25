# External comparison baselines

These sources support Section 8 of the manuscript. The experiment compares
fixed-precision numerical outputs. It does not time the complete certified
algorithm, and it is not a contest against every available implementation.

## Bellard

The experiment uses Fabrice Bellard's faster Gosper-based `pi1.c`, dated
February 26, 1997, from <https://bellard.org/pi/pi1.c>. It is unchanged and
compiled with the author's `HAS_LONG_LONG` option:

```sh
gcc -O3 -march=native -std=c17 -DHAS_LONG_LONG bellard_pi1.c -lm -o bellard_pi1.exe
```

The original attribution remains in the downloaded file. Its SHA-256 is
recorded in `provenance.json`. This third-party source is not redistributed
inside the manuscript ZIP: `compare_external.py --build` downloads the
exact source and checks its hash before compiling. The manuscript itself
compiles without any network access or ancillary execution.

Bellard outputs nine digits with the original double accumulator and
truncation parameter based on N+20. The benchmark does not substitute the
slower original `pi.c`, use floating modular multiplication, or add a
reference calculation to the timed program.

## Gourdon

`gourdon_window.cpp` is an adapter of this repository's
`legacy/gourdon_digit_extract.cpp`, implementing Algorithms 1 and 2 of
Xavier Gourdon's 2003 manuscript:
<https://plouffe.fr/simon/articles/nthdecimaldigit.pdf>.

This is the manuscript author's reimplementation. It is not Gourdon's
original `pidec` executable, and timings must never be labeled as `pidec`.
It retains the repository implementation's modular binomial recurrence,
symmetry optimization, parameter selection, default n0=12, and cached
prime sieve. Its fractional accumulator uses 64 bits.

The adapter prints nine digits from the same accumulated fractional part,
with exact integer scaling. It removes the original CLI's second
extraction and optional MPFR check. It does not add nine separate digit
calculations. Unused helpers are retained to make provenance easier to
audit. The original and adapter hashes and modifications are recorded.

The cached sieve means that this executable does not enforce the paper's
O(log^2 N) memory bound. Its measured process working set is reported
separately. It has not received the same degree of tuning as a production
implementation; ratios describe these sources, not every implementation
of Gourdon's method.

## Full-prefix MPFR control (context only)

`mpfr_prefix.cpp` constructs pi through decimal position N+8, with 192
additional precision bits, forms the scaled integer prefix, and returns
the final nine digits. It uses GMP and MPFR; it is not y-cruncher.
This program is outside the no-prefix decimal-extraction category.
Its measurements are retained in a separate panel as context and are
excluded from the extraction ranking.

## Reproduction

From the ancillary directory (or repository root), with GCC, GMP, MPFR,
Python, gmpy2, and psutil available:

```sh
python -B validation/compare_external.py --build
python -B validation/compare_external.py --check
python -B validation/compare_external.py --bench 10000 100000 --repeats 3 --output external_repeat.json
```

The Windows harness prepends `C:\msys64\ucrt64\bin` to PATH.
Compilation is completed before measurement. The five engines run in
fresh processes, sequentially, with rotating order. The timed window is
process wall time, including startup and printing. Peak process working
set is sampled every 1 ms. No warm-up executions are part of the benchmark;
the adapter correctness checks are recorded separately.

The baseline/control engines produce nine digits; our unchanged v5 and
weighted engines produce twelve. All twelve digits from our engines are
verified, and the first nine from every engine must agree. Native code
uses GCC 16.1.0 with -O3 -march=native; weighted Python uses the installed
GMP through gmpy2. Precision policies and implementation languages differ,
so this is an implementation comparison, not an experiment isolating the
effect of the weighted-block algorithm alone.

Verification uses explicit flooring of an MPFR reference at two guard
precisions and occurs outside all timed runs. Since the prefix control
also uses MPFR, its verification is a precision/adapter check, not an
independent-library check. None of this numerical agreement supplies the
paper's exact output certificate.
