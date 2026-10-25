# Research status

Updated 24 October 2026; the new weighted comparison is dated 24 October.
Earlier numerical measurements retain their original 7 September scope.
This file supersedes status statements in older
research logs, audits, and benchmark summaries. Measurements remain attached
to the exact sources and batches that produced them.

## Practical success criterion

The project requires an analytic improvement to produce a verified runtime
advantage in a corresponding implementation at relevant large positions.
A smaller exponent alone does not satisfy that practical objective.
The new weighted Machin construction meets the measured one-pass objective:
its derived half-log core improvement now beats the earlier Machin program.
The native complete adaptive worst-case bound remains unestablished.
The earlier optimized weighted BBP implementation beats matched Euler
controls but remains slower than both Machin variants. The measured Lucas
prototypes have not met the objective. Earlier
statements implying that their theoretical bounds already established the
fastest program, an observed crossover, or a world record are withdrawn.

For mathematical accuracy, a failed runtime comparison refutes a claimed
measured advantage; it does not alone refute numerical correctness or a
big-O upper bound. Two upper bounds need not describe tight growth rates.
Even tight bounds with different exponents do not locate a finite crossover
without constants and lower-order terms. A complexity claim about a program
also requires that its actual operations implement the analyzed algorithm.
We will not treat an unobserved crossover as a demonstrated benefit.

## Lean verification, 27 September 2026

The new [formalization](formalization/README.md) proves finite modular
reconstruction, exact correction rounding, integer certificate soundness,
binary-splitting identities, and the asymptotic balance of the displayed
cost expression. It connects reconstruction to exact returned rational
blocks and the BBP tail to decimal certification, with the remaining
analytic and representation hypotheses explicit. Ordinary irrationality
of pi is taken from mathlib's proof and gives qualitative eventual passing.

This is not a formal proof of the complete extractor or its bit complexity.
Open obligations include the BBP and transformed-series identities,
denominator-size estimates, the published quantitative irrationality bound,
and refinement of an executable into a justified bit-cost model. The native
weighted prototype remains unchanged. Formal checking supplies neither a
new benchmark advantage nor an independent human review of the full proof.

## What is established in the manuscript

The weighted-block manuscript gives a derivation of O(N log^(5/2) N) bit
time for certified fixed-width decimal extraction, with O(N) stored bits
plus arithmetic workspace. Its main theorem uses exact word aggregation
and a stopping argument based on an irrationality bound for pi. The native
weighted prototype implements the exact modular core but retains a
fixed-guard output path. Its timings are not timings of the complete
certified theorem. The proof has not received independent human review.

The numerical section includes the faster Machin comparison separately.
Classical Machin identities, modular fraction decomposition, and binary
splitting are not claimed as new. No first-ever quasi-linear extraction
claim is made; the paper retains Gourdon's earlier linear-memory discussion.

## Current methods

**Weighted blocks.** The manuscript's bound is O(N log^(5/2) N). The original
native prototype was slower than direct splitting and Euler in September.
The separate optimized implementation, added 24 October, is about three
times faster than that original and beats the matched Euler control at the
three tested positions. It certifies successful fixed-guard outputs and
rejects failures; it does not implement the complete adaptive exact-word
procedure. The new weighted Machin transfer, described below, is faster.

**Weighted Machin blocks, added later on 24 October.** The separate
`decimal_extract_weighted_machin.cpp` transfers local CRT reconstruction
and certified real integer correction to the signed base `-p^2`, for
p = 5 and 239. At one million the number of core denominators falls from
2,322,125 in weighted BBP to 425,622, and the real-series count falls from
863,686 to 227,893. The geometric Machin control has the same 425,622
denominators; the weighted transfer reduces its reconstruction work.

The [new derivation](validation/weighted_machin_derivation_2026-10-24.md)
covers denominators near one using an aggregate bound on real-series
lengths. An integer count schedule replaces reliance on floating scheduling
for that bound. The original implementation used local-size modulus checks
and bounded forests. The current LCM assembly, audited in the
[follow-up derivation](validation/weighted_machin_lcm_derivation_2026-10-24.md),
instead keeps exact interval denominators in one streaming balanced tree.
Together these support O(N log^(5/2) N) bit work for the default core and
for a pass with G = O(log N), with O(N) stored data plus arithmetic
workspace, in the manuscript's multiplication model. This is an audited
derivation, not a completed Lean proof or independent human review.

Successful output carries the existing integer interval certificate.
The executable's guard doubling, fixed input limits, and ordinary word
rounding still prevent assigning it the complete adaptive worst-case
theorem. A separately checked exact-word aggregation identity and a
direct jump to linear guard precision are documented as prospective
components; they are not part of the measured executable.

## Weighted Machin LCM follow-up, 24 October 2026

The current implementation replaces repeated denominator products with exact
LCM merges from retained prime-power metadata. One balanced streaming stack
replaces capped forests and their intermediate modular reductions. Every
local integer correction remains present. Balanced Newton targets, reused
prime-power moduli, nested operand residues, and binary multiplier stripping
also improve the ordinary Machin control; those shared fixes are included
in the new control instead of attributing their entire gain to weighted blocks.

The final matched batch uses three repetitions, rotating sequential fresh
processes, identical compiler flags, and initial guard 192:

| Position | Original log-cubed | Improved log-cubed | Previous weighted | LCM weighted |
| ---: | ---: | ---: | ---: | ---: |
| 100,000 | 0.115628 s | 0.101732 s | 0.092878 s | 0.079033 s |
| 1,000,000 | 2.408116 s | 2.307414 s | 2.238210 s | 1.867780 s |
| 10,000,000 | 29.792373 s | 24.749944 s | 23.018332 s | 17.983396 s |

At ten million, the new implementation uses
21.9% less time than the previous weighted implementation and
27.3% less than the faster measured log-cubed control.
Assembly falls from 13.159 s to 8.289 s (separate phase medians).
All 36 outputs match directed MPFR intervals, every certificate succeeds on
the first pass, and source/executable hashes remain unchanged during timing.
Host runtime varied substantially between batches; these ratios use only
this final batch. These finite measurements do not establish an empirical
complexity exponent or a crossover outside the measured positions.

The [pre-change profile](validation/weighted_machin_assembly_profile_2026-10-24.json)
isolated inversion at 5.771 s of 12.914 s assembly in its own ten-million
sample. Generic GCD normalization and direct GMP inversion were rejected
after slower exploratory comparisons. The retained factor metadata avoids
generic GCD work while preserving exact rational addition.

The [source audit and derivation](validation/weighted_machin_lcm_derivation_2026-10-24.md)
establish the same one-pass O(N log^(5/2) N) bound with O(N) stored data
in the stated multiplication model. Nested inversion now has O(M(N)) work.
The stack diagnostic counts retained numerator/denominator bits and factor
vector capacity; it is not peak RSS or an arithmetic-workspace measurement.
Log^2 and a complete adaptive worst-case bound remain open for this program.
The note records the separate, unimplemented Lucas transfer toward a
candidate O(N log^2 N log log N) construction and its unresolved requirements.

See the [36-run report](validation/weighted_machin_lcm_final_2026-10-24.json),
[additional validation](validation/weighted_machin_lcm_checks_2026-10-24.json),
and [reproduction guide](validation/README.md#weighted-machin-24-october-2026).
The previous weighted source is preserved byte-for-byte as
`validation/decimal_extract_weighted_machin_pre_lcm_2026-10-24.cpp`.
The manuscript and Lean files have not been updated to claim this follow-up
as a new verified exponent.

## Weighted Machin measurements, 24 October 2026 (earlier batch)

That earlier batch used five repetitions per method and position, identical
compiler flags and initial G=192, sequential fresh processes, and rotating
order. Both the previous fastest log-cubed control and its general-prime
half-width Newton variant are measured; the comparison uses the faster
control median. That batch's weighted source retained the original
general-prime Newton formula. Median compute times were:

| Position | Previous Machin | Half-width Newton control | Weighted Machin | Less time vs faster control |
| ---: | ---: | ---: | ---: | ---: |
| 100,000 | 0.072690 s | 0.073415 s | 0.063078 s | 13.2% |
| 1,000,000 | 1.823595 s | 1.848537 s | 1.518385 s | 16.7% |
| 10,000,000 | 20.257954 s | 20.574741 s | 16.111973 s | 20.5% |

All 45 timed outputs matched independent directed-rounding MPFR intervals;
all certificates passed on the first attempt, and all measured source and
executable hashes remained unchanged. Absolute runtime varied between
batches; the ratios above use this batch only. At ten million, every weighted
sample was faster than every sample from either control (weighted range
15.651-18.184 s; previous-control range 19.775-24.869 s).

Additional validation checks 612 inverses against GMP, 72 boundary-range
signed blocks, 14 complete position cores through one million, 43 decimal
windows, three insufficient-guard retries, and four partition overrides.
The native selftests and separate prospective exact-word checks also pass.

At that stage, the measured bottleneck was rational assembly: 9.285 s of the
16.112 s median total at ten million, followed by real correction at 2.659 s
and word residues at 2.692 s. Phase medians are taken separately and need
not sum exactly to the median total.

See [the final raw benchmark](validation/weighted_machin_final_2026-10-24.json),
[additional checks](validation/weighted_machin_checks_2026-10-24.json), and
[reproduction instructions](validation/README.md#weighted-machin-24-october-2026).
The earlier four-engine exploration is retained separately, including the
preserved source with half-width Newton inversion. No historical timing is
used to compute the displayed speed ratios.

## Weighted optimization measurements, 24 October 2026

The new batch builds all engines with the same compiler and flags, rotates
method order in sequential fresh processes, uses three repetitions at each
position, and checks every window by independent MPFR intervals outside
timers. The certified Euler and Machin control copies receive applicable
common word-arithmetic, accumulation, and Newton-inversion improvements.
All certified methods start with 192 guard bits and pass on the first
attempt. No historical timings are mixed into the following comparisons.

Median compute times for original weighted, optimized weighted, matched
Euler, and matched Machin respectively:

- Position 100,000: 0.466, 0.130, 0.134, and 0.072 seconds.
- Position 1,000,000: 5.920, 1.832, 1.981, and 1.086 seconds.
- Position 10,000,000: 84.755, 28.186, 31.237, and 17.504 seconds.

Thus weighted has a measured advantage over this Euler control (about 8%
and 11% higher throughput at one and ten million). This does not establish
an asymptotic exponent, an external record, or a win over Machin. The
manuscript, original weighted source, and Lean formalization are unchanged.

The bottlenecks were duplicate modular exponentiation, per-term allocations
and hash lookups, repeated word inversions, and excessive exact real-series
work. Optimizations preserve weighted reconstruction and its default block
schedule. At one million the real-series term count falls from 1,258,432 to
863,686 while retaining the exact integer remainder check. Large real,
rational, and assembly trees now dominate the remaining runtime.

All 36 timed outputs matched independent MPFR intervals. Additional checks
cover 43 decimal windows, the complete 2,322,125-bit residue at one million,
word-range boundary arithmetic, and rejection of an insufficient guard.
See [the raw benchmark](validation/weighted_optimization_2026-10-24.json),
[additional checks](validation/weighted_optimized_checks_2026-10-24.json),
and [reproduction instructions](validation/README.md).

## Other current methods

**Lucas and sparse prime-power batches.** The subsequent proof draft targets
O(N log^2 N log log N), with a later parameter refinement. This is a candidate
analytic result outside the manuscript. A complete correspondence between
the native implementation and the abstract cost model is not established.
The original measured prototype was slower than its controls at one million
and ten million. Repairing two leaf operations did not reverse that result.
The newest parameter schedule has not been implemented and measured.

**Earlier Machin prime-power cores.** This was the fastest program before
the weighted Machin transfer and certifies every returned decimal window. It reduces the leading
number of core terms from approximately 1.29917 N for the compared Euler
construction to 0.42556 N. The conservative time bound for an attempt with
G = O(log N) guard bits is still O(N log^3 N). Stored arithmetic data have
the conservative bound O(N log log N), plus primitive workspace. If adaptive
certification needs a much larger guard, the separate word rounding can
introduce quadratic work. The native program lacks the exact-word fallback
needed to remove that cost; the restricted-guard estimate is not an
unconditional complexity bound for the complete adaptive program.

**Euler controls.** The repository retains both the earlier fixed-guard
implementation and a matched certified mode in the Machin executable.
The latest certified comparison uses the same style of stopping test and
the same initial guard in both methods. Older Euler timings refer to older
source versions and separate batches.

In the 7 September batch at ten million, Machin took a three-run median of 19.086 seconds versus
38.195 seconds for the certified Euler control, using 138.1 versus 75.1 MiB
peak working set. This supports a roughly twofold local speed improvement
with higher memory use. It does not establish an external record or a
smaller logarithmic exponent.

## Evidence and remaining work

The [reproduction guide](validation/README.md) identifies measured source
snapshots, complete samples, and validation. Full-prefix MPFR computations
are independent controls outside the no-prefix decimal-extraction ranking.
Historical external comparisons concern Bellard's original implementation
and the documented Gourdon low-memory reimplementation, not every possible
implementation or Gourdon's linear-memory regime.

The unresolved combined objective is a complete, audited implementation of
a smaller-bound algorithm that also wins matched certified benchmarks at
the intended scale. No runtime at one billion, one trillion, or one
quadrillion is measured by the present reports. No crossover or run-time
estimate for those positions is asserted here.

Local derivations and audits:

- [Weighted-block proof](notes/weighted_block_crt_2026-09-06.md)
- [Lucas proof draft](notes/lucas_sparse_prime_power_batch_2026-09-07.md)
- [Lucas implementation audit](notes/lucas_proof_implementation_audit_2026-09-07.md)
- [Lucas parameter refinement](notes/lucas_smaller_parameters_2026-09-07.md)
- [Machin construction and certificate](notes/machin_prime_core_extraction_2026-09-07.md)

This cleanup preserves the theorem text and historical evidence. It updates
the current interpretation and presentation; it is not a new independent
proof audit.
