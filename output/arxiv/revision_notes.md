# Manuscript revision

## Current consolidation (September 7, 2026)

This entry supersedes the status summaries below, which describe earlier
revisions. The main theorem and the proof sections are unchanged. The
numerical section now includes a self-contained description of the Machin
prime-power computation, its integer certificate, and a third performance
table with all three controls at three positions. The historical comparison
tables retain their original samples and are explicitly labeled as earlier
batches. File paths, development versions, and reference-language details
remain outside the article.

The text records that the weighted prototype has not demonstrated the
desired practical advantage. The latest Machin program is faster than its
matched certified Euler control through ten million, with higher memory
use; its restricted-guard O(N log^3 N) estimate is not the main theorem's
bound or a complete unrestricted-retry bound. No world record or observed
weighted/Lucas crossover is asserted. The Lucas proof drafts and prototypes
remain outside this manuscript and its source archive.

The repository now has one current research-status note, a concise README,
and separate reproduction instructions. Older audits and the attack log
are marked as historical. The previous manuscript and package were copied
to a local backup before editing. The source archive retains all earlier
code and raw measurement payloads, updates the ancillary guide and article,
and adds the five files needed for the latest certified comparison.

The PowerShell packaging script checks every displayed performance row,
measured-source hashes, theorem preservation, references, embedded fonts,
and an isolated three-pass TeX build. Current page count, file hashes, and
rendering-review status are in `verification.json`. This is an editorial
and reproduction check, not independent mathematical review or submission.

## Historical revision entries

## Academic presentation revision (September 7, 2026)

The manuscript now places the theoretical comparison in Section 8 and
the numerical experiments in Section 9. Repository paths, command-line
instructions, development version labels, and Python performance
comparisons have been removed from the article. They remain in the
supplementary reproduction material where appropriate.

The experimental methods are named weighted reconstruction, direct binary
splitting, and 2-adic Euler summation. Each is defined in the text. The
Euler identity is included to explain the comparison method; it is the
identity already implemented and checked in the archived source.

Two native C/C++ performance tables replace the three development-oriented
tables. Wall time is used consistently. The historical baseline comparison
retains three-run medians, full ranges, output widths, and memory. The
larger-position table distinguishes the weighted three-run median from
single-run direct-splitting and Euler observations. Full-prefix MPFR
results remain in a separate contextual panel.

The main theorem and mathematical Sections 2 through 7 are unchanged.
The fixed-guard implementation limitation is retained once in the
experimental methods, and numerical validation is reported concisely.
Measured code and data are unchanged. Supplementary instructions map
the new method names to the original source and dataset identifiers.

The preceding manuscript and package are preserved under
`notes/manuscript_backups/before_academic_edit_2026-09-07/`.

## No-prefix comparison scope (September 7, 2026)

The abstract and introduction now define the comparison task as base-10
digit extraction without first computing the prefix, with linear stored
arithmetic data allowed. The theorem and mathematical proofs are unchanged.

The complexity table groups eligible decimal extractors separately from
other bases, full-prefix algorithms, and contextual formulas. The external
performance table has a no-prefix panel and a separate MPFR prefix-control
panel. Prefix controls are excluded from extraction rankings; their faster
timings do not adjudicate a no-prefix speed claim.

The results explicitly state that C++ v5 has the lowest median wall time
among the tested no-prefix programs at both benchmark positions. At
100,000 its ratios are 299.2 against Bellard and 51.1 against our Gourdon
adapter. Native weighted ratios remain 175.8 and 30.0. All measurements,
source versions, output widths, and validation limitations are retained.

The O(N log^(5/2) N) theorem applies to the weighted construction, not v5.
Gourdon's earlier linear-memory discussion remains relevant prior work;
no fastest-known no-prefix bound or world record is asserted. The README,
ancillary instructions, submission notes, metadata, and archive are aligned
with this scope. Earlier entries below describe preceding revisions.

The preceding source, PDF, README, package, and verification records are
preserved under `notes/manuscript_backups/before_no_prefix_scope_2026-09-07/`.

## Performance-table update

Section 8 now contains two measured performance tables. The first gives
three-run medians, full timing ranges, peak process working set, and
block/group counts for the Python weighted-block prototype. Its nine
fresh runs are checked against MPFR. The second retains the four
recorded C++ v3/v5 samples, including phase costs and memory.

The paper states the hardware, software, sample counts, timing boundaries,
memory units, and verification method. It distinguishes the older C++
engines from the weighted-block theorem and excludes projected runs from
the measurement tables. Raw data and the additional C++ sources and
measurement scripts are included in the updated source archive.

The preceding 14-page version is preserved in
`notes/manuscript_backups/before_performance_tables_2026-09-06/`.

The paper has been rewritten around the weighted-block algorithm and
its certified decimal-output extension. The original source, PDF, and
README are preserved under
`notes/manuscript_backups/pre_arxiv_revision_2026-09-06/`.

## Mathematical changes

- Theorem 1.1 states the certified fixed-width decimal result.
- The least common denominator bound is proved by prime valuations,
  including uniform control of short final blocks.
- Prime-power partial fractions recover a block modulo the integers.
  A transformed real series, an explicit remainder inequality, and
  exact rounding recover its missing integer correction.
- The weighted-core proof accounts for factor enumeration, small-modulus
  powers and inverses, sequential record sorting, every tree level,
  consecutive-group assembly, and arithmetic workspace.
- The two local tree costs are balanced explicitly to obtain the
  O(N log^(5/2) N) bound.
- Exact aggregation of word contributions permits linear guard precision
  without N separate full-width accumulator updates.
- An integer decimal-cell test handles both carries and wraparound.
  The published inequality mu(pi) < 8 bounds first-pass precision for
  all sufficiently large positions; explicit doubling handles the
  finitely many exceptions.
- The algorithm uses actual denominator bounds on repeated passes,
  so small-input termination does not depend on asymptotic constants.

## Claims and references

The rewrite removes unsupported priority and speed-record claims,
claims that all preceding decimal algorithms are quadratic,
overbroad exclusions of alternative methods, equivalences between
unweighted and weighted tasks, and blanket generalizations to all
BBP constants and output bases. These are not retained as theorems.

The introduction discusses Gourdon's explicit linear-memory regime
and the faster full-prefix AGM bound. The contribution is the
specific constructive algorithm and its analysis.

The bibliography contains twelve directly relevant sources after the
external-comparison revision.
Primary sources were checked for the BBP formula, fast multiplication,
computer-arithmetic primitives, hypergeometric reconstruction, historical
extraction methods, AGM comparison, and the irrationality bound.
Barry Gergel's author initial is corrected to B.

## Presentation and reproduction

The source is one self-contained LaTeX article with numbered lemmas,
theorems, equations, and internal cross-references. It uses standard
packages, embedded vector fonts, a fixed date, and an inline bibliography.
The algorithm summary is kept together, and the references use a
consistent compact layout.

The full certified algorithm is explicitly distinguished from the
fixed-guard prototype. The source archive includes the prototype's local
imports, pinned dependency, recorded millionth-position diagnostic,
and output from a fresh reproduction check.

Build warnings, source-archive contents, isolated compilation, font
embedding, PDF page inspection, and package self-test results are
recorded in `verification.json`. These are editorial and reproducibility
checks. No formal proof-assistant verification, independent referee
review, or arXiv submission is claimed.

## External comparison revision (September 7, 2026)

The external-comparison revision initially had 18 pages, three performance tables, and a
complexity table covering ten methods or memory regimes. The complexity
table distinguishes bit operations, word operations, and the historical
Plouffe count, as well as output base, prefix computation, and stored
data. It includes Gourdon's intermediate and linear-memory regimes,
the AGM bound, Chudnovsky/Ramanujan binary splitting, and Plouffe's
2022 formulas without assigning the latter an unsupported total cost.

Thirty new timed runs compare Bellard's unchanged, faster Gosper C
implementation, our Gourdon adapter, C++ v5, weighted Python, and a
full-prefix MPFR control at 10,000 and 100,000. There are three
repetitions of every engine/position combination. Every output passed
the separately timed numerical reference check.

At position 100,000, median wall-time ratios favor C++ v5 by 299.2
against Bellard and 51.1 against the Gourdon adapter. The weighted
Python ratios are 46.5 and 7.9. The full-prefix control is faster
than both local extractors. The paper reports memory and different
output widths, and limits these speedup claims to the tested programs.
The Gourdon entry is explicitly a reimplementation, not original pidec.

The archive includes the experiment, data, adapters, build commands,
and provenance. Bellard's original source is fetched with a hash check
by the build command; it is not redistributed inside the archive.

## Native weighted implementation

The current manuscript is 19 pages. The principal weighted extractor is
now native C++/GMP, with Python retained as a reference. The native code
also reuses modular weights along a prime's incidence progression using
the exact ratio 16^(-p) modulo its local prime power. The full certified
algorithm remains specified in the paper; native output uses G=192.

The native source passed 300 exact rational block checks, 120 complete-core
checks, five comparisons of large complete residues against the separate
Euler construction, and 45 MPFR decimal checks through position one million.
Nine fresh runs give computation medians 0.053, 0.640, and 9.507 seconds at
10,000, 100,000, and one million. The Python median at one million was
28.720 seconds. Native peak memory is 17.4 MiB there.

Table 1 compares native and Python computation medians. Table 4 adds the
native wall times as a separate measurement batch. At 100,000, native
weighted wall time is 0.654 seconds: ratios 175.8 against Bellard and
30.0 against the Gourdon adapter. The earlier C++ v5 remains faster at
tested overlapping sizes. No crossover or overall speed record is claimed.
