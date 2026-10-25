# Reproduction guide

Run commands from the project root. Source names are retained so recorded
hashes and historical measurements remain traceable. Algorithm names in
the paper do not depend on those filenames.

The current fastest implementation and matched comparison are documented
under [Weighted Machin](#weighted-machin-24-october-2026).

## Earlier certified Machin and Euler comparison

Requirements: C++17, GMP, MPFR for the independent verifier, and PowerShell 7
for the benchmark harness. The extractor itself uses GMP only. The recorded
Windows build used GCC 16.1.0, GMP 6.3.0, and MPFR 4.2.2.

```powershell
$env:PATH = 'C:/msys64/ucrt64/bin;' + $env:PATH
g++ -O3 -march=native -std=c++17 decimal_extract_machin.cpp -lgmpxx -lgmp -o extract_machin.exe
g++ -O3 -march=native -std=c++17 decimal_extract_shared_residues.cpp -lgmpxx -lgmp -o extract_shared_residues.exe
g++ -O3 -march=native -std=c++17 validation/mpfr_window_interval.cpp -lmpfr -lgmpxx -lgmp -o validation/mpfr_window_interval.exe
./extract_machin.exe --selftest
./extract_machin.exe 1000000 --timing
./extract_machin.exe 1000000 --euler --timing
pwsh -File validation/benchmark_machin.ps1 -OutputPath validation/machin_benchmark_local.json
```

Use a new report filename when rerunning. The default report path in the
historical harness names the original dataset, so the explicit `-OutputPath`
above prevents replacing it. The default benchmark covers positions
100,000, 1,000,000, and 10,000,000, with three fresh processes per method
and position, rotating order, and no discarded samples. Memory is peak
working set sampled every 25 ms; short processes may be under-sampled.
All returned outputs are checked by directed MPFR intervals outside timers.
The reference program computes a prefix and is used only for verification.

[machin_benchmark_2026-09-07.json](machin_benchmark_2026-09-07.json) contains
27 timed runs, 224 smoke comparisons, 10 forced guard retries, 2,445 residue
checks, and 108 exact prime-power core checks. It hashes the extractor,
its shared source dependency, the verifier, and the harness. The certified
modes start at 128 guard bits; the archived Euler control uses a fixed 192.

## Source mapping

- `decimal_extract_weighted.cpp`: weighted reconstruction, exact core and
  fixed-guard decimal output; the main manuscript's prototype.
- `decimal_extract_v3.cpp`: direct binary splitting in the earlier tables.
- `decimal_extract_v5.cpp`: earlier 2-adic Euler summation in those tables.
- `decimal_extract_shared_residues.cpp`: shared-residue Euler optimization;
  the archived fixed-guard control in the latest comparison.
- `decimal_extract_machin.cpp`: certified Machin prime-power cores and the
  matched certified Euler mode; includes the shared-residue source.
- `decimal_extract_lucas.cpp`: current Lucas experiment; includes the
  weighted source. This experiment is outside the submitted manuscript.
- `decimal_extract_lucas_measured_2026-09-07.cpp`: preserved source of the
  original Lucas measurements, before the leaf-operation repair.
- `baselines/gourdon_window.cpp`: documented Gourdon low-memory
  reimplementation; `baselines/README.md` records its limitations.
- `baselines/mpfr_prefix.cpp` and `mpfr_window_interval.cpp`: independent
  prefix controls, excluded from no-prefix performance rankings.

## Earlier measurements and validation

The paper labels earlier benchmark batches separately from the latest
certified comparison. Do not combine their times into a matched speed ratio.

- [external_comparison_2026-09-07.json](external_comparison_2026-09-07.json):
  historical Bellard, Gourdon, Euler, and MPFR batch; baseline provenance
  is in [baselines/provenance.json](baselines/provenance.json).
- [weighted_cpp_performance_2026-09-07.json](weighted_cpp_performance_2026-09-07.json):
  nine weighted timings, measured in a later batch on the same workstation.
- [cpp_euler_final_benchmark_2026-09-06.json](cpp_euler_final_benchmark_2026-09-06.json):
  earlier direct-splitting and Euler samples, including the larger positions.
- [weighted_cpp_checks_2026-09-07.json](weighted_cpp_checks_2026-09-07.json):
  large-core and decimal validation for the weighted implementation.
- [shared_residue_benchmark_2026-09-07.json](shared_residue_benchmark_2026-09-07.json):
  exploratory shared-residue batch; its larger trials encountered concurrent
  CPU work and are not the latest matched comparison.
- [lucas_benchmark_raw_2026-09-07.json](lucas_benchmark_raw_2026-09-07.json):
  original Lucas samples; [normalized report](lucas_benchmark_2026-09-07.json)
  retains the same measurements. Its original source hash refers to the
  preserved measured snapshot, not the repaired current source.
- [lucas_leaf_audit_2026-09-07.json](lucas_leaf_audit_2026-09-07.json):
  checks and single-run comparisons after repairing leaf operations.

The weighted native self-test can be reproduced with:

```powershell
g++ -O3 -march=native -std=c++17 decimal_extract_weighted.cpp -lgmpxx -lgmp -o extract_weighted.exe
./extract_weighted.exe --selftest
./extract_weighted.exe 1000000 --timing
```

Earlier mathematical reference programs and benchmark drivers are retained
for provenance and independent checks. They use `gmpy2==2.3.1` and
`psutil==5.9.6`; they are optional for native extraction. Their commands and
the earlier experimental narrative are preserved in the local README backup
under `notes/manuscript_backups/before_consolidation_2026-09-07/`.

## Weighted optimization, 24 October 2026

`decimal_extract_weighted_optimized.cpp` preserves the weighted-block
construction and adds shared BBP word residues, contiguous prime incidence
groups, batch inversion, reciprocal word reduction, packed small leaves,
a tighter certified real-series cutoff, and one final modular inverse.
It includes `decimal_extract_weighted_words.hpp`. The earlier weighted
source remains unchanged; the benchmark also preserves its exact snapshot
as `validation/decimal_extract_weighted_baseline_2026-10-24.cpp`.

The optimized program prints a window only after an integer interval
certificate passes. The default guard is 192 bits. `--guard=N` changes it;
an insufficient guard produces an error. This is a bounded fixed-guard
implementation, not the manuscript's complete adaptive exact-word
procedure. `--block=N` is an explicit experimental partition override;
the default retains the manuscript's `ceil(sqrt(log2 N))` schedule.

Build and run a fresh matched comparison without overwriting reports:

```powershell
pwsh -File validation/run_weighted_optimization.ps1 -OutputPath validation/weighted_optimization_local.json
```

The wrapper builds the preserved weighted source, the optimized source,
and certified Euler and Machin controls. The control copies under
`validation/weighted_controls_common/` receive the same applicable word
arithmetic, accumulation, and Newton-inversion improvements. All three
certified methods start with 192 guard bits. The controls retain their
existing retry behavior; guard and pass counts are recorded. Runs use
sequential fresh processes, rotating order, three repetitions, source and
executable hashes, and independent MPFR interval checks outside the timers.
No simultaneous performance work should run during this comparison.

Independent arithmetic checks, including exact GMP rational oracles,
32-bit boundary cases, negative odd denominators, and varied block sizes:

```powershell
$env:PATH = 'C:/msys64/ucrt64/bin;' + $env:PATH
g++ -O3 -std=c++17 validation/weighted_optimized_arithmetic_checks.cpp -lgmpxx -lgmp -o tmp/weighted_optimized_arithmetic_checks.exe
./tmp/weighted_optimized_arithmetic_checks.exe
```

These tests and timings do not prove an asymptotic exponent. The optimized
program's complete adaptive bit-complexity theorem remains outside their
scope.

## Weighted Machin, 24 October 2026

`decimal_extract_weighted_machin.cpp` transfers weighted CRT reconstruction
to Machin's prime-power cores. The signed base is -p^2; a local p-power
clears nonunit denominators, and a certified positive Euler series recovers
the missing signed integer. The [block derivation](weighted_machin_derivation_2026-10-24.md)
includes blocks starting at denominator one.

The current [LCM assembly and inverse derivation](weighted_machin_lcm_derivation_2026-10-24.md)
replaces repeated denominator factors with metadata-based LCM merges and a
balanced streaming stack. Nested Newton lifting reuses moduli and operand
residues. The same inversion and binary scaling optimizations are applied
to the improved log-cubed control.

The default schedule supports a one-pass O(N log^(5/2) N) bound for
G=O(log N), with O(N) stored data plus arithmetic workspace. It does not
establish log^2 or the complete adaptive worst-case bound. Experimental
`--block=N` overrides are outside the default-schedule cost claim.

Reproduce the current comparison with a fresh report path:

```powershell
pwsh -File validation/run_weighted_machin_lcm.ps1 -OutputPath validation/weighted_machin_lcm_local.json
```

The wrapper builds four engines with identical flags and initial guard 192:
the original log-cubed control, its improved arithmetic variant, the previous
weighted source snapshot, and current weighted LCM. It checks their native
selftests, both weighted cores, and independent arithmetic before timing.
Three repetitions at 100,000, one million, and ten million use sequential
fresh processes and rotating order. Independent directed MPFR verification
runs outside timers; source and executable hashes are checked before/after.

The [final report](weighted_machin_lcm_final_2026-10-24.json) contains 36
verified outputs:

| Position | Original log-cubed | Improved log-cubed | Previous weighted | LCM weighted |
| ---: | ---: | ---: | ---: | ---: |
| 100,000 | 0.115628 s | 0.101732 s | 0.092878 s | 0.079033 s |
| 1,000,000 | 2.408116 s | 2.307414 s | 2.238210 s | 1.867780 s |
| 10,000,000 | 29.792373 s | 24.749944 s | 23.018332 s | 17.983396 s |

At ten million the new version uses
21.9% less time than the previous weighted source and
27.3% less than the faster log-cubed control. Compare times within
this batch: absolute host runtime varied across the day's experiments.

Reproduce the [additional checks](weighted_machin_lcm_checks_2026-10-24.json):

```powershell
$env:PATH = 'C:/msys64/ucrt64/bin;' + $env:PATH
g++ -O3 -march=native -std=c++17 decimal_extract_weighted_machin.cpp -lgmpxx -lgmp -o extract_weighted_machin.exe
g++ -O3 -march=native -std=c++17 validation/weighted_machin_arithmetic_checks.cpp -lgmpxx -lgmp -o tmp/weighted_machin_arithmetic_checks_final.exe
g++ -O3 -march=native -std=c++17 validation/weighted_machin_position_core_checks.cpp -lgmpxx -lgmp -o tmp/weighted_machin_position_core_checks.exe
python validation/validate_weighted_machin.py --mpfr tmp/weighted_machin_lcm_final/mpfr_window_interval.exe --output validation/weighted_machin_lcm_checks_local.json
```

Checks cover 612 inverses against GMP, 72 exact signed blocks including
denominators near UINT32_MAX, 540 exact signed scaling cases, 14 complete
prime-power cores through one million, 43 decimal windows, three forced
guard retries, partition overrides, and rejected CLI inputs. Native tests
add 2,445 word residues, 108 direct cores, 180 exact blocks, 216 geometric
core comparisons, 5,489 denominator invariants, and 2,489 exact LCM merges.

The [pre-change profile](weighted_machin_assembly_profile_2026-10-24.json)
and `decimal_extract_weighted_machin_profile_2026-10-24.cpp` preserve
the assembly investigation. Build that source with the same flags and GMP
libraries, then run `10000000 --timing` to reproduce the instrumentation.
Its sample is separate from the final matched timing batch.

Earlier [45-run results](weighted_machin_final_2026-10-24.json) and
[checks](weighted_machin_checks_2026-10-24.json) remain unchanged.
Their weighted source is preserved byte-for-byte in
`decimal_extract_weighted_machin_pre_lcm_2026-10-24.cpp`.
The still-earlier four-engine exploration and half-width Newton source
remain in `weighted_machin_2026-10-24.json` and
`decimal_extract_weighted_machin_half_newton_2026-10-24.cpp`.
Older wrappers remain historical comparison tools; the LCM wrapper above
is the current reproduction entry point.

The exact-word fallback remains prospective, with separate Python checks:
130 fraction identities, 520 projections through 4096 guard bits, and 24
deep prime-cancellation cases. It and its guard schedule are not integrated.

## Parallel weighted Machin and direct reconstruction, 24 October 2026

`decimal_extract_weighted_machin_parallel.cpp` preserves the selected serial
weighted arithmetic and certificate. A bounded worker group runs the two
prime cores and disjoint word ranges concurrently. Each job owns its GMP
integers and counters; exact signed term floors are combined after joining.
Invariant diagnostics use the serial path. The serial production source
and earlier benchmark records remain unchanged.

The [matched report](weighted_machin_parallel_reconstruction_2026-10-24.json)
contains three fresh-process repetitions per engine and position, rotating
engine order, guard 192, and independent directed-rounding MPFR verification
outside the timers. All 24 outputs passed, with unchanged tracked sources
and executables. Medians of the extractor's elapsed `total` time were:

| Position | Improved log-cubed, serial | Weighted, serial | Weighted, 8 threads | Direct-H Lucas, serial |
| ---: | ---: | ---: | ---: | ---: |
| 1,000,000 | 1.858599 s | 1.532399 s | 0.654414 s | 2.109336 s |
| 10,000,000 | 23.895055 s | 16.614774 s | 10.114286 s | 24.840227 s |

At ten million the parallel version is 1.643 times as fast as weighted
serial and 2.363 times as fast as the improved log-cubed control. This
machine reports 16 logical processors; the timed parallel version used
eight. These gains use additional processors. Fixed-thread total work
retains the existing one-pass O(N log^(5/2) N) scope and O(N) stored data;
neither a log^2 bound nor an adaptive worst-case theorem follows. Host
runtime varied: weighted serial at one million ranged from 1.107 to
1.641 seconds in this batch. Compare within the stored batch, not with
earlier absolute times. Parallel phase fields sum overlapping worker
elapsed times; use `total` or process wall time for comparisons.

The [additional interval checks](weighted_machin_parallel_checks_2026-10-24.json)
passed 43 decimal windows, three forced guard retries, four block overrides,
and rejected CLI inputs. The [native exact-check record](weighted_machin_parallel_native_checks_2026-10-24.json)
separately records 147 serial/parallel accumulator comparisons across guards
1 through 257 and thread counts 2, 3, and 8; four matching adaptive
certificates; and the serial fallback for invariant diagnostics.

Reproduce with a fresh report path:

```powershell
pwsh -File validation/run_weighted_machin_parallel.ps1 -OutputPath validation/weighted_machin_parallel_local.json
```

The runner builds all controls and checkers with the same optimization and
thread-link flags, validates before timing, and saves a build manifest with
compiler/runtime information, commands, hashes, and exact-check records.
It uses the byte-preserved pre-reconstruction weighted baseline. The
recorded experiment above used the same optimization flags; `-pthread`
was additionally used to link the new parallel executable.
The runner itself passed an end-to-end smoke run at position 1,000;
its [build and exact-check manifest](weighted_machin_parallel_build_checks_2026-10-24.json)
records all eight builds and the full native checks, including the final
13-adic native batching and Newton inversion.

The experimental `decimal_extract_weighted_machin_lucas.cpp` makes a
different mathematical change. It reconstructs H directly, using
1 <= H < C/(C-1) and the exceptional 13-primary fraction, so the integer
correction needs no child reconstruction. The first recursive prototype
had 17 recursive invocations at one million; the direct variant needs
only the two initial solves. Independent tests passed 594 fixed-prime
projections and 378 complete H comparisons. Its global-precision 13-adic
kernel remains costly, and the matched table shows that it loses to the
selected serial engine. A fixed cutoff retains the recursive backend
outside the finite direct-H regime. The
[derivation and remaining cost barrier](machin_lucas_recursive_design_2026-10-24.md)
make no claim of a new global exponent. This prototype is not selected
for the fast build.

## Manuscript build

```powershell
pwsh -File validation/package_manuscript.ps1
```

This creates a source-only archive with an inline bibliography, compiles
that archive in a fresh directory with shell escape disabled, checks source
hashes and the new timing table against the stored report, and writes the
current PDF and verification report. It uses PDFLaTeX, pdfinfo, pdftotext,
and pdffonts. Tool paths can be passed as script parameters. It does not
rerun benchmarks, infer a complexity exponent from timings, submit to arXiv,
or certify the mathematics. Rendered-page inspection is a separate check.
