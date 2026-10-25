# Audit and dependency plan for the weighted-block theorem

Date: 27 September 2026. Scope: the current `decimal_extraction_paper.tex`,
especially lines 166–972, and the two weighted-block derivation/audit notes.
This is a mathematical and formalization audit, not a Lean verification.
It does not certify the manuscript or any executable.

## Finding

The manuscript supplies a plausible constructive argument for certified,
fixed-width decimal extraction in `O(N log^(5/2) N)` bit operations. This
audit found no counterexample to the principal algebraic identities or the
displayed exponent balance. That is weaker than either an independently
reviewed proof or a complete kernel-checked theorem.

The hardest remaining obligations are the connection between an executable
algorithm and a bit-cost model, the fast arithmetic algorithms in that model,
the analytic series identities and bounds, and a quantitative irrationality
bound for pi. Proving the final inequality after assuming the manuscript's
cost expression would verify the balance of those costs only. It would not
verify the claimed complexity of an algorithm.

The existing native weighted implementation uses a fixed guard and omits the
exact word aggregation and adaptive certificate of the complete theorem.
`RESEARCH_STATUS.md:26` and `decimal_extraction_paper.tex:1122` state this
explicitly. Its numerical tests and benchmarks cannot discharge the missing
formal obligations. Historical `notes/paper_audit.md` concerns an earlier
manuscript and expressly disclaims verification of this result. The September
7 consolidation review was editorial and preserved the proof text.

## Exact claim to formalize

Fix an output width `w >= 1`. For each position `N >= 1`, the mathematical
output is

`floor (10^w * frac (10^(N-1) * Real.pi))`.

A complete theorem should identify a particular deterministic program, prove
that it terminates and returns that integer, and bound its execution cost by
`C * N * (log N)^(5/2)` for sufficiently large `N`. Constants and thresholds
may depend on the fixed width. Use `max 1 (ceil (log2 N))` or an equivalent
positive logarithm in intermediate estimates to avoid artificial failures
at `N = 1`.

The space claim is `O(N) + S_arith(C1*N)` bits. The stronger `O(N)` total
space conclusion is conditional on linear workspace for the actual chosen
arithmetic primitives. A multiplication time theorem alone does not prove
that workspace condition.

No-prefix extraction is an additional operational property of this program,
not a consequence of output correctness or complexity. Specify its dataflow:
it constructs modular residues and the fractional part of a shifted sum,
without invoking a prefix calculation or forming that sum's integer part.
Do not define the property as “no intermediate integer ever equals the
prefix”: incidental numerical equality on small inputs would be irrelevant.

## Mathematical proof obligations, in dependency order

### 1. Integer conventions, finite sums, and modular inverses

Use positive natural moduli, integer numerators, and Euclidean remainders.
Specify empty products, empty prime sets, and the modulus-one convention.
For odd positive denominators, prove all inverses modulo `2^m` exist. For
prime-power partial fractions, prove the denominator with its full p-part
removed is coprime to the corresponding prime power. Natural subtraction
must not silently replace signed exponents: introduce `J = -E` only under
`E < 0`, and prove every remaining exponent nonnegative.

An executable rational representation should store an unreduced numerator
and positive denominator. Using Lean's normalized `Rat` representation for
execution can insert gcd reductions absent from the analyzed algorithm.
`Rat` remains a useful specification semantics.

### 2. Short progression denominator bound

Manuscript lines 241–286. For odd `A`, power-of-two `a >= 2`, and `t >= 1`,
prove `odd(t!)` divides `prod (A+a*j)` and that the lcm of the progression
divides the product divided by `odd(t!)`, times `odd(lcm(1,...,t))`.

The valuation argument needs a formal count of one invertible residue class
in an interval of `t` consecutive indices, Legendre's factorial valuation,
valuation of a finite lcm, and oddness of all relevant integers. Explicitly
prove the integer quotient exists before dividing.

Then prove the uniform logarithmic size bound
`log D <= t*log(C*N/t) + O(t)`, including arbitrarily short final blocks.
The finite divisibility theorem alone does not establish this bound.
The manuscript uses a factorial lower bound (Stirling is stronger than
necessary) and `log(lcm(1,...,x)) = O(x)`.

### 3. Elementary prime estimates

Manuscript lines 224–237. Prove the required Chebyshev estimate, the lcm
bound, and `sum_{p<=x} 1/p = O(log log(x+2))`. The proposed route uses the
central binomial coefficient, dyadic intervals, the prime-power identity
for the lcm, and partial summation. No prime number theorem is required.

Count prime-power incidences and large residual prime factors explicitly.
It is insufficient to formalize the prime estimates as unproved hypotheses
and then describe the final result as unconditional.

### 4. Partial fractions modulo the integers

Manuscript lines 290–324. Given distinct primes with adequate powers `q_p`,
show that the specified `c_p` satisfy
`sum_j w_j/mu_j - sum_p c_p/q_p` being an integer.

The proof multiplies by `D* = prod q_p`, checks divisibility by every `q_p`,
and combines pairwise coprime divisors. Signed weights are allowed. Extra
primes have zero coefficients. Modulo-the-integers reconstruction is not
yet reconstruction modulo a power of two: the integer correction must not
be omitted.

### 5. Real correction identity and remainder certificate

Manuscript lines 372–470. Prove the infinite real-series identity, its
extension error at most `1/(15*(A+a*t))`, alternating terms, and first
omitted term estimate. The manuscript's integral derivation requires actual
justification of summation/integration interchange and of the beta integral.
An equivalent finite identity with a bounded remainder would also suffice
and might avoid unnecessary special-function infrastructure.

Prove the binary-splitting merge invariant for an arbitrary consecutive
interval, not just its root. At the root, `U/Q = sum_{r<R} h_r` and
`P/Q = h_(R-1)`, with `Q > 0` and `R >= 1`. The integer inequality

`2^(4*t+6) * abs(P) * a * R < 225*A*(A+a*R)*Q`

is exactly the first-omitted-term test for truncation error below `1/64`.
This audit rechecked that factor `225` and the power of two.

Keep two errors distinct in formal statements: the truncated rational
approximation's error and the subsequent flooring error. The integer test
certifies the former; by itself it does not certify the sum of both errors
below the same `1/64`. The final block argument correctly allows both.

For sufficiently large `N`, prove that
`R = ceil((8*ell+16)/u)` already passes. For all other inputs, prove that
doubling `R` terminates using the uniform ratio bound below `1/15`.
An asymptotic inequality alone does not make the small-input algorithm total.

### 6. Exact integer rounding and reconstruction

Manuscript lines 472–515. With the stipulated `g`, bound the sum of flooring
errors by the actual number of fractions, and show

`abs ((B-C0)/2^g - (S-T)) < 1/15 + 1/64 + 1/256 + 1/256 < 1/2`.

Since `S-T` is an integer, nearest-integer rounding recovers it, including
negative corrections. Prove the signed floor formula with Euclidean integer
division; truncation toward zero would be wrong for negative inputs.
Then prove the returned numerator/denominator represents the exact block.

### 7. Weighted modular core and assembly

Manuscript lines 519–669. Prove alignment of each block with its outside
dyadic shift, correctness of the rational merge, and correctness of root
reduction using the odd denominator's inverse modulo `2^(m-s)`.

The linear storage estimate relies on grouping consecutive blocks of one
family. Formalize the group-flush invariant and bound the last-block
overshoot. At every tree level, sum the spans of disjoint consecutive
intervals; a bound on each span separately is not enough. Signs remain
short signed integers until reduction. Storing every block at once would
retain `O(N*u)` bits and fail the advertised linear data bound.

The theorem quantifies over uniformly bounded families and fixed constants
`c,C,a`. Pointwise `Theta(N)` hypotheses without uniform constants are not
sufficient for a uniform complexity conclusion.

### 8. Exact word aggregation

Manuscript lines 679–747. Apply partial fractions with global prime powers,
prove exact recovery of the small integer correction, and bound the complete
rational numerator and denominator by `O(N)` bits. Residues must be
recomputed in prime order as described, with their cost counted. Separately
rounding `O(N)` terms into a `Theta(N)`-bit accumulator is a different
algorithm and can introduce quadratic work.

Prove that every global maximal prime power exceeds `sqrt(B0)`. Together
with the total bit mass this bounds the number of fraction records and
their metadata. Omitting zero fractions changes the denominator to a divisor
of the global product and must be covered in correctness statements.

### 9. Shifted BBP identity and exact fractional consolidation

Manuscript lines 65–70 and 752–833. Prove the BBP series for `Real.pi`, the
four normalized families, the explicit integer cutoff, and its tail bound.
A theorem about an abstract number that is merely assumed to equal the BBP
series is an intermediate conditional result, not yet a theorem about pi.

Prove the Bezout decomposition for negative exponents, including `mu=1`;
lift inverses from `2^J` to the common modulus; and show that `U/Q` is exactly
the fractional part of the finite shifted sum. Preserve all negative-exponent
terms through the cutoff. The actual guarded algorithm can have a much
larger cutoff than the fixed-guard prototype.

### 10. Decimal certificate, termination, and guard growth

Manuscript lines 835–940. A useful independent theorem is certificate
soundness for any real target within `2^(-G)` of a rational approximation.
Its two strict inequalities prove that both numbers lie in the same open
decimal cell, including cells adjacent to fractional wraparound.

Ordinary irrationality of pi and convergence of the approximations prove
eventual termination for each fixed input. They do not provide the uniform
`G = O(N)` bound. The latter needs, at minimum, the following separately
proved quantitative theorem:

`exists q0, forall q >= q0, forall p : Int, abs(pi - p/q) > q^(-8)`.

Here division and powers are interpreted in the reals and `q > 0`.
If the source theorem only quantifies over reduced fractions, formally
derive the all-denominators version, also handling the finitely many small
reduced denominators. The cited Zeilberger–Zudilin result implies this bound
in ordinary mathematics. A citation does not create a Lean theorem.

From it, show the distance to each decimal boundary exceeds
`10^(-7*n-8*w)`. Check the elementary integer/exponential inequality for
`G=32*(N+w+1)`, then prove the first attempt passes for sufficiently large
`N`. The unknown threshold is used only in analysis and must not be an
oracle consumed by the executable.

The finite exceptional positions need a total adaptive implementation.
On retries, use the actual largest denominator and actual precision for
prime enumeration and rounding guards, as required at manuscript lines
964–972. A hardcoded `C*N` limit from first-pass asymptotics would not justify
retry correctness. Finite terminating exceptions can be absorbed into
asymptotic constants without furnishing an explicit practical threshold.

## Required algorithmic and complexity infrastructure

Choose a precise deterministic bit-machine model or an equivalent model
with a proved simulation. Formalize binary encodings, integer sizes,
sequential records, tree-level storage, and primitive step costs. The
manuscript explicitly charges movement and sorting; unit-cost arithmetic
on arbitrary Lean `Nat` or `Int` does not match that model.

The primitive dependency set includes multiplication in `O(b log(b+2))`,
division with remainder in the corresponding bound, odd inversion modulo
`2^b`, and general modular inversion with its extra logarithm. Include
sign handling, shifts, addition, comparisons, record generation and copying.
The Harvey–van der Hoeven theorem is a substantial external result. Unless
an appropriate machine-level formal theorem is located or independently
formalized, expose it as an explicit assumption and label the resulting
complexity theorem conditional. The same applies to reductions proving
division/inversion costs. Lean's native big integers or GMP do not inherit
the cited theorem merely by having a multiplication operation.

Next implement and charge each described procedure: trial-division base
prime enumeration; segmented record sieving; modular powering; prime-power
incidence generation; merge sorting and scans; rational and real trees;
consecutive block groups; exact word aggregation; final modular arithmetic;
and both adaptive loops. Prove that the accumulated cost is bounded by

`N*u*L^2 + N*L^3/u + N*L*u^2 + N*L^2*(log L)^2`

plus the stated sublinear enumeration costs. Only then apply the exponent
balance. Formalization of `L,u,ell`, ceiling bounds, uniformity in block
length, and negligibility of the enumeration terms is still required.

A natural-valued majorant such as `N * L^2 * ceilSqrt(L)` is convenient
for an executable cost bound. Prove its equivalence to the real-valued
`N*(log N)^(5/2)` scale separately. Arithmetic correctness and asymptotic
cost must refer to the same algorithm, not two independently specified
functions with an unproved correspondence.

Track maximum live encoded bits rather than just final output sizes. A
functional implementation may retain persistent old tree levels or list
spines; either prove a storage discipline for it or refine to an explicitly
managed machine. Keep the arithmetic workspace term in the final theorem.

## Suggested deliverables and completion criteria

1. A kernel-checked foundation of finite algebra: modular consolidation,
   signed rounding, rational merges, CRT partial fractions, LCM divisibility,
   and the abstract decimal certificate. State their scope precisely.
2. A kernel-checked analytic layer: the series identity, remainder bounds,
   BBP formula, irrationality termination, and the quantitative guard theorem.
   Isolate any unformalized pi-measure dependency in the theorem signature.
3. An executable algorithm with totality and output-correctness proofs.
   Finite examples supplement, but do not replace, those proofs.
4. A costed refinement with all enumeration, copying, arithmetic and storage
   obligations. Publish conditional complexity results with named primitive
   and analytic hypotheses while dependencies remain open.
5. An unconditional end-to-end theorem only after those hypotheses have
   kernel-checked proofs in the chosen model, plus the operational no-prefix
   property and a documented trust boundary for any extracted executable.

For each milestone pin Lean and mathlib versions, compile from a clean
checkout, and report the exact theorem name and its full type. Check the
transitive axiom dependencies (`#print axioms`) and scan for `sorry`,
`admit`, `sorryAx`, and project-specific axioms. Standard Lean foundational
axioms should be listed; extra trusted evaluation mechanisms must be
disclosed if used. A successful build with assumptions is a successful
conditional formalization, not evidence that those assumptions were proved.

Do not label a proof of only the balance inequality, finite tests, or a
certificate checker as the full `O(N log^(5/2) N)` theorem. Do not label a
theorem about an abstract certified rational approximation as a theorem
about this weighted-block implementation without the missing refinement.
