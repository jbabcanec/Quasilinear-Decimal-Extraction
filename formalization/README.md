# Lean formalization of weighted BBP components

27 September 2026. **Partial formalization: the complete extractor and its
`O(N log^(5/2) N)` bit-complexity theorem are not yet formally verified.**
This directory contains actual Lean proofs, not theorem stubs. Missing
mathematics remains visible as hypotheses or as obligations described below;
it is not replaced by new axioms.

The [27 September verification](verification-2026-09-27.json) passed the
full build and axiom audit for 56 explicitly named theorems in seven modules.
The audit reports 182 theorem declarations and 258 total declarations after
including Lean-generated equations and helpers.

## Reproduce the checks

The project pins Lean **4.33.0** and mathlib **v4.33.0**, commit
`db584cd6d46c92f209a44c0f1c829460d327499d`. With elan installed, run from
this directory:

```powershell
lake update
lake exe cache get
pwsh -File verify.ps1
```

The first two commands prepare a fresh checkout. Once dependencies exist,
only the verification command is needed. The direct checks are:

```powershell
lake build
lake env lean Audit.lean
```

`verify.ps1` checks the installed mathlib revision, builds the complete
library, audits every public `BBP` declaration's transitive axioms, and
records tool versions and source hashes in `build-verification.json`.
It rejects changes to source files during verification. Full command output
is retained under `.lake/verification-*.log`. The dated verification report
saved beside this README records the checked source snapshot; subsequent
edits require a new run.

Only `propext`, `Classical.choice`, and `Quot.sound` are allowed by the axiom
audit. These are the ordinary logical foundations used by mathlib. The audit
rejects `sorryAx`, custom axioms, and `Lean.ofReduceBool` (native evaluation
used as a proof). A passing axiom audit does **not** discharge hypotheses in
a theorem's type or prove that the statement captures the entire algorithm.
Generated equations and helper theorems are included in the audit's counts.

## What is checked

- **Reconstruction.lean:** rational and shifted merges; finite CRT
  reconstruction; local modular-inverse algebra. The strongest result,
  `partial_fractions_of_local_inverses`, derives integrality of the difference
  between the weighted sum and its prime fractions from explicit denominator
  factorizations, inverse equations, and coefficient congruences.
- **Certificate.lean:** exact nearest-integer recovery, including negative
  integers; the paper's error budget; decimal certificate soundness and
  exclusion of wraparound; positive decimal-grid distance using mathlib's
  actual proof of `irrational_pi`; eventual passing for all sufficiently
  large guard exponents whenever a valid approximation is supplied.
- **Block.lean:** connects CRT integrality to rounding. Under the explicit
  modular data and approximation budget, `weighted_block_crt_reconstruction`
  proves that the returned corrected numerator and denominator represent the
  intended finite weighted sum exactly.
- **Series.lean:** executable integer splitting-tree operations and their
  correctness for arbitrary tree shapes; specialization to the hypergeometric
  factors; ratio contraction and finite coefficient bounds; BBP term bounds
  and the exact integer cutoff's geometric remainder estimate. The BBP tail
  theorem explicitly assumes `HasSum bbpTerm Real.pi`.
- **Guard.lean:** proves the strict integer guard inequality, scales an
  exponent-eight rational-approximation inequality to decimal boundaries,
  and proves the resulting error margin. The quantitative approximation
  inequality is a hypothesis, not a formal proof of the published result.
- **CostBalance.lean:** proves that choosing `u = ceil(sqrt L)` balances the
  four stated costs, including the logarithmic overhead. It proves a genuine
  `Asymptotics.IsBigO` theorem for that expression, then specializes to
  `N * (log N)^(5/2)`. It does not assign this expression to a program's
  execution cost. Natural logarithms in this module change the asymptotic
  constants from the manuscript's base-two convention. A further theorem
  covers natural positions and the manuscript's discrete
  `L = max(1, ceil(log2 N))`, with `u = ceil(sqrt L)`.
- **Pipeline.lean:** `certified_bbp_window` connects the proved tail bound
  and integer certificate to the requested digits. It assumes the BBP
  identity and an exact rational representation of the finite shifted sum.
  `initial_pass_certified_bbp_window` additionally proves that the stated
  initial guard passes when the quantitative approximation hypothesis holds
  at the specified denominator.

Read the theorem signatures before interpreting these results. For example,
an exact approximation supplied as a premise does not establish that the
current native program computes it, nor does qualitative eventual passing
bound the cost of an executable loop.

## Obligations still open

1. Prove the BBP identity for `Real.pi` and the transformed real-series
   identity, with the required extension and truncation error estimates.
   The geometric BBP tail estimate is already checked *conditional on the
   identity*; the finite correction coefficient bounds are already checked.
2. Formalize the short-progression lcm/factorial divisibility and uniform
   denominator-size estimate, the prime-power valuation facts, and the
   prime-incidence bounds used in the routing cost.
3. Formalize the published quantitative irrationality bound that supplies
   the exponent-eight hypothesis eventually. Mathlib's ordinary irrationality
   proof supplies qualitative termination, not the linear worst-case guard.
4. Specify an executable with exact word aggregation, block generation,
   certified rounding, forest assembly, and the adaptive output test; prove
   it supplies every premise of the integrated block and pipeline theorems.
   The existing fixed-guard C++ weighted prototype is unchanged.
5. Give the executable a formal sequential bit-cost and storage model,
   justify the chosen fast-arithmetic implementations in that model, and
   prove the denominator, tree, routing, assembly, and exceptional-input
   bounds. Then connect those bounds to `balanced_cost_isBigO` and prove
   the intended no-prefix dataflow property.

The last two items are essential: correctness of a formula and a bound on
an abstract cost expression do not constitute an algorithm-complexity proof.
The formalization also establishes no empirical crossover or speed record.

[PROOF_AUDIT.md](PROOF_AUDIT.md) gives the detailed dependency analysis and
source references. It is a prose audit, not an additional Lean proof.
