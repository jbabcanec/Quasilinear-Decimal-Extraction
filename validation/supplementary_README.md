# Supplementary reproduction material

For *Weighted block reconstruction for certified decimal extraction of pi*,
Joseph Babcanec, September 2026.

The manuscript's main file is one directory above this supplementary folder.
Its inline bibliography needs no BibTeX step. From the archive root, run
PDFLaTeX twice with shell escape disabled:

```sh
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error decimal_extraction_paper.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error decimal_extraction_paper.tex
```

Run the following extraction commands from this `anc/` directory.

## Certified Machin and Euler comparison

The extractor requires C++17 and GMP. The independent verifier additionally
requires MPFR; the benchmark harness uses PowerShell 7. On the measured
Windows workstation, add `C:/msys64/ucrt64/bin` to PATH.

```sh
g++ -O3 -march=native -std=c++17 decimal_extract_machin.cpp -lgmpxx -lgmp -o extract_machin.exe
g++ -O3 -march=native -std=c++17 decimal_extract_shared_residues.cpp -lgmpxx -lgmp -o extract_shared_residues.exe
g++ -O3 -march=native -std=c++17 validation/mpfr_window_interval.cpp -lmpfr -lgmpxx -lgmp -o validation/mpfr_window_interval.exe
./extract_machin.exe --selftest
./extract_machin.exe 1000000 --timing
./extract_machin.exe 1000000 --euler --timing
pwsh -File validation/benchmark_machin.ps1 -OutputPath validation/machin_benchmark_local.json
```

The extractor includes the shared-residue source; keep both files together.
It returns twelve digits only after an integer interval certificate passes.
Position 1 is the first decimal digit. Position one million returns
`130927562832`. MPFR is used only by the external verifier, outside timers.
It computes a prefix and is excluded from the extraction comparison.

The stored `validation/machin_benchmark_2026-09-07.json` contains all 27 timed
runs, the 224 smoke checks, ten forced guard retries, 2,445 residue checks,
and 108 exact modular-core checks. It records source hashes. Both certified
modes begin with 128 guard bits; the archived Euler mode uses a fixed 192.
Runs were fresh processes, sequential and in rotating order, with no
discarded samples. Peak working set was sampled every 25 ms. Short processes
may be under-sampled. Always use a new output filename when rerunning the
harness to preserve the historical dataset.

The Machin attempt has the conservative bound O(N log^3 N) for
G = O(log N), with O(N log log N) stored bits plus arithmetic workspace.
Large adaptive guards can introduce quadratic word-accumulation work.
The program therefore does not establish the main theorem's smaller bound
or an unconditional O(N log^3 N) bound for its full adaptive procedure.

## Earlier implementations and tables

The source mapping is:

- Weighted reconstruction: `decimal_extract_weighted.cpp`.
- Direct binary splitting: `decimal_extract_v3.cpp`.
- Earlier 2-adic Euler summation: `decimal_extract_v5.cpp`.
- Shared-residue Euler control: `decimal_extract_shared_residues.cpp`.
- Certified Machin and Euler: `decimal_extract_machin.cpp`.

The earlier three programs use a fixed 192-bit guard. The weighted program
implements the exact weighted core but does not implement the complete
exact word aggregation and adaptive output procedure in the main theorem.

```sh
g++ -O3 -march=native -std=c++17 decimal_extract_weighted.cpp -lgmpxx -lgmp -o extract_weighted.exe
g++ -O3 -march=native -std=c++17 decimal_extract_v3.cpp -lgmpxx -lgmp -o extract_v3.exe
g++ -O3 -march=native -std=c++17 decimal_extract_v5.cpp -lgmpxx -lgmp -o extract_v5.exe
./extract_weighted.exe --selftest
```

The earlier reports are:

- `validation/external_comparison_2026-09-07.json`: historical Bellard,
  Gourdon, Euler, and MPFR comparison, including all samples.
- `validation/weighted_cpp_performance_2026-09-07.json`: nine weighted
  timings from a subsequent batch on the same workstation.
- `validation/cpp_euler_final_benchmark_2026-09-06.json`: the earlier
  direct-splitting and Euler larger-position samples.
- `validation/weighted_cpp_checks_2026-09-07.json`: complete large-core and
  decimal checks, supplementing the 300-block and 120-core native self-test.

The paper retains these batches separately. The subsequent optimized
Euler control is a different implementation from the earlier Euler source.
Do not combine their timings to form a matched speed ratio.

`validation/baselines/README.md` documents the Bellard source and the Gourdon
reimplementation, including its memory and precision policies. Bellard's
source is not redistributed here; the historical comparison driver obtains
the author's source with a recorded hash check. No download is needed for
the manuscript build or for the native Machin and weighted self-tests.

## Optional reference checks

The reference implementations and earlier benchmark drivers are retained
for reproduction. They require the versions in `requirements.txt`:

```sh
python -m pip install -r requirements.txt
python -B decimal_bbp_extract.py --selftest
python -B validation/validate_weighted_cpp.py
```

`selftest_output.json` is the preserved reference self-test report from the
earlier package. Source files, raw reports, and baseline provenance have
been retained without rewriting historical paths or measurements.
Python is not required by any native extractor or by the latest benchmark
harness. Numerical validation does not establish the complexity exponent,
independent verification of the proof, or an external speed record.
