# Decimal extraction of pi

Research on computing a fixed-width base-10 window of pi without computing
its preceding digits. Position 1 is the first digit after the decimal point.

The manuscript, **Weighted block reconstruction for certified decimal
extraction of pi**, presents an O(N log^(5/2) N) bit-complexity construction.
The new weighted Machin implementation combines that reconstruction idea
with Machin's smaller prime-power cores and an integer output certificate.
Its core has a separate O(N log^(5/2) N) derivation, and it now beats the
earlier log-cubed Machin implementation in matched local tests. The bound
applies to a pass with G = O(log N); a complete worst-case bound for the
native adaptive guard retries is still outstanding. No external speed
record is established.

The 24 October optimized weighted prototype is about three times faster
than the original. In matched local tests it beats the optimized certified
Euler control by about 8% at one million and 11% at ten million. Machin
remains faster than that BBP-based version. That program certifies successful fixed-guard outputs;
it is not the complete adaptive implementation of the theorem. See the
[new benchmark and build instructions](validation/README.md#weighted-optimization-24-october-2026).

The latest [weighted Machin comparison](validation/README.md#weighted-machin-24-october-2026)
uses exact LCM assembly and nested Newton inversion. Matched three-run
medians at ten million were 23.018 s for the previous weighted
program, 17.983 s for the new version, and 24.750 s for the
faster of the original and improved log-cubed controls. The new version
takes 21.9% less time than the previous weighted program and
27.3% less than the faster log-cubed control in this batch. All 36
timed outputs passed independent directed-rounding MPFR checks.
The [LCM derivation](validation/weighted_machin_lcm_derivation_2026-10-24.md)
preserves the one-pass log^(5/2) bound and linear stored data; it does not
establish log^2. Earlier source snapshots and benchmark reports are preserved.

The [parallel weighted build](validation/README.md#parallel-weighted-machin-and-direct-reconstruction-24-october-2026)
now runs the two independent prime cores and word ranges concurrently. In
the latest matched three-run batch, eight threads took 10.114 s at ten
million, compared with 16.615 s for weighted serial and 23.895 s for the
improved log-cubed control. All 24 timed windows passed MPFR verification.
This is a wall-time improvement with the existing log^(5/2) work bound.
The new exact direct-H reconstruction removes recursive corrections but
was slower (24.840 s); it remains experimental and establishes no log^2 bound.

See [RESEARCH_STATUS.md](RESEARCH_STATUS.md) for the current claims,
limitations, and the practical success criterion recorded on 7 September 2026.

The [Lean formalization](formalization/README.md), added 27 September 2026,
checks reconstruction, certification, and cost-balancing components. The
complete extractor and its bit-complexity theorem are not yet formally verified.

## Paper

- [Manuscript PDF](decimal_extraction_paper.pdf) and [LaTeX source](decimal_extraction_paper.tex)
- [Source package](output/arxiv/decimal_extraction_source.zip)
- [Submission notes](output/arxiv/SUBMISSION.md) and [metadata](output/arxiv/metadata.txt)
- [Build verification](output/arxiv/verification.json)

The paper contains mathematical definitions, proofs, complexity comparisons,
and measured performance. File mappings and commands belong in the
[reproduction guide](validation/README.md). Nothing has been submitted to arXiv.

## Native extraction

Build with a C++17 compiler and GMP. Python and MPFR are not needed by this
extractor. On the Windows workstation used for the measurements:

```powershell
$env:PATH = 'C:/msys64/ucrt64/bin;' + $env:PATH
g++ -O3 -march=native -std=c++17 decimal_extract_weighted_machin.cpp -lgmpxx -lgmp -o extract_weighted_machin.exe
./extract_weighted_machin.exe --selftest
./extract_weighted_machin.exe 1000000 --timing
```

For the faster multicore build on this workstation:

```powershell
g++ -O3 -march=native -std=c++17 -pthread decimal_extract_weighted_machin_parallel.cpp -lgmpxx -lgmp -o extract_weighted_machin_parallel.exe
./extract_weighted_machin_parallel.exe 1000000 --threads=8 --timing
```

The parallel program defaults to at most eight logical threads;
`--threads=1` uses its original serial path. Its `total` timing is elapsed
wall time. Per-phase times sum overlapping workers and should not be added
to estimate elapsed time.

The last command returns `130927562832`. The executable returns twelve
digits only after its integer interval test passes. The new source is
self-contained. It defaults to 192 guard bits and retries an insufficient
guard; `--guard=N` changes the starting precision. `--block=N` is an
experimental partition override, outside the default-schedule cost claim.
The program enforces input and precision limits and reports an error if
it cannot certify a result. The previous `decimal_extract_machin.cpp`
remains available, including its `--euler` comparison mode.

## Measurements

The 7 September Machin batch used sequential fresh processes, rotating method order,
three runs per method and position, and independent MPFR interval checks
outside the timers. Median wall times for Machin and certified Euler were:

- Position 100,000: 0.108 and 0.231 seconds.
- Position 1,000,000: 1.256 and 2.756 seconds.
- Position 10,000,000: 19.086 and 38.195 seconds.

At ten million, peak working sets were 138.1 MiB and 75.1 MiB respectively.
All 27 timed outputs, including the archived fixed-guard Euler control,
passed the independent checks. The [raw report](validation/machin_benchmark_2026-09-07.json)
records individual samples, guards, source hashes, and validation results.
These measurements establish a local implementation comparison through
ten million. They do not establish a crossover for weighted blocks or Lucas,
an extrapolation to larger positions, or a world record.

## Research and reproduction

- [Current research status](RESEARCH_STATUS.md): proof, program, and performance claims.
- [Reproduction guide](validation/README.md): source mappings, builds, datasets, and controls.
- [Chronological research log](OPEN_PROBLEM_E_ATTACK_LOG.md): retained experiments and corrections.
- `notes/`: local proof drafts, audits, and manuscript backups; excluded from the source package.
- `legacy/`: earlier implementations retained for reference.
- `tmp/`: disposable build and rendering output; excluded from version control.
- `output/pdf/` and `output/arxiv/`: current deliverables.

Rebuild the manuscript and source package with:

```powershell
pwsh -File validation/package_manuscript.ps1
```

This command compiles an isolated copy of the submission sources and checks
the package. PDF visual inspection is recorded separately after rendering.
The previous manuscript, README, package, and audit are preserved locally in
`notes/manuscript_backups/before_consolidation_2026-09-07/`.
