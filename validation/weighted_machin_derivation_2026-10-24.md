# Weighted Machin blocks: exactness and an aggregate half-log core bound

24 October 2026. This note audits the proposed transfer of weighted-block
reconstruction to the two prime-power cores in Machin's identity. It records
an algebraic derivation and a complexity argument for the core, including
denominators near one. It is not a formal verification, a claim that a
benchmark establishes an exponent, or a complete complexity proof for the
adaptive decimal extractor.

The implementation mapping below refers to `decimal_extract_weighted_machin.cpp` as
read on this date. Native arithmetic remains restricted to its checked word
range. The asymptotic discussion uses the same bit-arithmetic model as the
manuscript, with multiplication cost `M(b) = O(b log(b+2))` and separately
accounted arithmetic workspace.

## 1. The exact prime-power core

Fix the constant prime `p` (5 or 239), a truncation index `K`, and a cancelled
power `s` (`N-1` for the 5-core, zero for the 239-core). Put

\[
 A_0=2K+1,\qquad h=\lfloor\log_p A_0\rfloor,
 \qquad m=A_0+h-s.
\]

There is no core if `m <= 0`. Otherwise let

\[
 R=\min(K+1,\lceil m/2\rceil),\qquad B=-p^2.
\]

Up to the outside numerator and sign already present in the Machin
extractor, the required residue is

\[
 V=p^h\sum_{r=0}^{R-1}\frac{B^r}{A_0-2r}\pmod {p^m}.
\]

Every summand is a p-adic integer: if `v = v_p(A_0-2r)`, then `v <= h`,
and cancellation leaves the p-unit denominator `(A_0-2r)/p^v` and the
factor `p^(2r+h-v)`. Terms with `r >= ceil(m/2)` vanish modulo `p^m`, since
`2r+h-v >= 2r >= m`. Thus this omission does not use an approximation.

For a block `lo <= r < hi`, write `t=hi-lo` and

\[
 A=A_0-2(hi-1),\quad
 h_i=\lfloor\log_p(A+2(t-1))\rfloor,\quad P=p^{h_i},
\]
\[
 F=P\sum_{j=0}^{t-1}\frac{B^{t-1-j}}{A+2j}.
\]

The block's contribution is exactly

\[
 (-1)^{lo}p^{\,2lo+h-h_i}F.
\]

Here `h_i <= h`, and `F` has a p-unit denominator after cancellation. A
block with outside exponent `2lo+h-h_i >= m` vanishes and must be skipped
before requesting precision `m-(2lo+h-h_i)`.

## 2. Why the CRT reconstruction omits the prime p

The multiplier `P` cancels every p-power occurring in a block denominator.
For another prime `q`, let `q^e` be its local maximum power and let
`q^v || mu_j`, where `mu_j=A+2j`. The local coefficient is

\[
 c_q\equiv\sum_{q\mid\mu_j}
 \frac{q^e}{q^v}\,P B^{t-1-j}
 (\mu_j/q^v)^{-1}\pmod {q^e}.
\]

The inverse exists because only the q-part was removed. The ordinary CRT
partial-fraction argument then gives

\[
 F-\sum_{q\ne p}\frac{c_q}{q^e}\in\mathbb Z.
\]

There is no missing fractional contribution from `p`: it was cancelled
term by term by `P`. Signed powers of `B` are ordinary integer weights;
they do not change this argument.

## 3. Recovering the integer correction

Let `C=p^2+1` and define

\[
 S_B(A)=\sum_{j\ge0}\frac{B^{-j-1}}{A+2j}
       =-\frac1{CA}\sum_{r\ge0}H_r,
\qquad
 H_r=\frac{2^r r!}{C^r\prod_{i=1}^r(A+2i)}.
\]

The equality follows by expanding
`1/(B-x^2) = -(1/C)/(1-(1-x^2)/C)` under the integral on `[0,1]`.
All transformed terms are positive, and

\[
 H_{r+1}/H_r=\frac{2(r+1)}{C(A+2(r+1))}<1/C.
\]

The finite block satisfies the exact identity

\[
 F=P\{B^tS_B(A)-S_B(A+2t)\}.
\]

Its omitted second term obeys

\[
 |P S_B(A+2t)|
 \le\frac{P}{(p^2-1)(A+2t)}<\frac1{p^2-1},
\]

because `P <= A+2(t-1)`. This remains uniform at `A=1`.

For a split tree truncated before index `r=q`, let `U/Q` be the sum of the
retained transformed terms and `P_tree/Q=H_(q-1)`. The remaining positive
tail is at most `H_q C/(C-1)`. Therefore the exact integer test

\[
 64\,P\,p^{2t}P_{\rm tree}(2q)
 < C(C-1)A(A+2q)Q
\]

proves an absolute error below `1/64` in the approximation to
`P B^t S_B(A)`. The sign of that quantity is `(-1)^(t+1)`.

Use signed floor division for its fixed-point approximation. With the
same local guard convention as the BBP blocks, each of the real flooring
error and the accumulated CRT-fraction flooring error is below `1/256`.
The total distance from the true integer correction is consequently less
than

\[
 \frac1{p^2-1}+\frac1{64}+\frac1{256}+\frac1{256}<\frac12.
\]

Nearest-integer recovery is exact, including a negative correction. The
inequality also holds for `p=3`, which the prototype uses in its tests.

## 4. Uniform blocks suffice even when A approaches one

For one extraction pass with `G=O(log N)`, both Machin core sizes and
their maximum denominators are `Theta(N)`. Set

\[
 L=\max(1,\lceil\log_2 N\rceil),\quad
 u=\lceil\sqrt L\rceil,\quad
 \ell=\lceil N/2^u\rceil,
\]

up to constant rescalings for the individual core. Let `b=Theta(2^u)`
denote the number of blocks. It is unnecessary to assume `A >= cN` for
every block, or to impose a dyadic partition. An aggregate estimate is
enough.

Order the uniform blocks by increasing starting denominator `A`. Except
for a possible endpoint block, these starts are spaced by `2 ell`.

* For the first `ceil(sqrt(b))` blocks, the bound `H_q < C^-q` gives a
  sufficient count `q=O(ell+L)`. Their combined count is
  `O(N/sqrt(b)+sqrt(b)L)=o(N/u)`.
* Every remaining block has `A >= c ell sqrt(b)` for a fixed positive
  constant `c`. Since a crude sufficient count is `q_0=O(ell+L)=O(ell)`,
  every ratio through `q_0` is at most `2q_0/(CA)=O(1/sqrt(b))`.
  Thus `q=O((ell+L)/u)` suffices. Summing these counts gives
  `O(N/u+bL/u)=O(N/u)`.

Consequently `sum_i q_i=O(N/u)`, even though an individual block close
to one may need `Theta(ell)` transformed terms. Its larger local tree
has `O(ell L)=o(N)` bits, which is still compatible with linear total
data storage.

### An integer schedule that implements this bound

The current implementation uses the following deterministic upper bound before
optionally trying a smaller floating-point estimate. For a block of
length `t`, set

\[
 T=\operatorname{bits}(P)+2t\operatorname{bits}(p)+12,
 \qquad q_0=\left\lceil T/\lfloor\log_2 C\rfloor\right\rceil.
\]

Let

\[
 d=\left\lfloor\log_2
      \left\lfloor\frac{CA}{2q_0}\right\rfloor\right\rfloor.
\]

If the inner integer is zero or `d=0`, retain `q_0`. Otherwise use
`min(q_0,ceil(T/d))`. The original count works because `H_q<C^-q`.
The reduced count is at most `q_0`, so
`H_q <= (2q_0/(CA))^q <= 2^(-dq)`. Both bounds are strong enough for the
integer remainder test above. The same formula with `2` replaced by
the progression step covers the helper's general parameter.

The optional log-gamma estimate only decreases this certified upper
schedule. If its proposal fails the exact remainder test, doubling
reaches an accepted count no larger than twice the deterministic
schedule. Thus floating-point accuracy does not establish either the
returned rational's correctness or the aggregate term-count bound.

## 5. Core work and data bounds

The short-progression LCM estimate only needs an upper bound on the
denominators, not a lower bound on `A`. Every returned block denominator
therefore has `O(ell u)` bits; deleting its p-part can only shorten it.
The corrected numerator has that size plus `O(ell+L)` bits.

The rational CRT trees cost `O(NuL^2)` in total. The real trees cost

\[
 O\left(L^3\sum_i q_i\right)=O(NL^3/u).
\]

Factoring, sorting, and local modular arithmetic have the same bounds
as the manuscript's progression construction, with the fixed constants
`p`, `B`, and `C` changing no exponent. Their aggregate cost fits
`O(N L^2 (log L)^2)+N^(3/4+o(1))`.

All returned block data total `O(Nu)` bits. Process consecutive blocks
in forests whose input mass is a fixed constant times `N`, rather than
retaining every block. There are `O(u)` forests. Store only consecutive
tree levels and a constant number of accumulated roots. At each level,
the exponent spans of disjoint consecutive intervals sum to `O(N)`;
the changes in local valuation are monotone and telescope. Hence each
forest has `O(N)` data and `O(u)` levels. Forest work totals
`O(u^2 M(N))=O(NLu^2)`.

Multiplication by `p^k` has bit cost and must be charged. Generate powers
only when needed, with lifetime bounded by the current operation or
forest. Before constructing `p^precision` merely to trim an integer,
the exact preliminary test

\[
 \operatorname{bits}(|x|)\le
 \text{precision}\,\lfloor\log_2 p\rfloor
\]

proves that no reduction is needed. If this test fails, the requested
modulus has only `O(bits(x))` bits, so its construction is chargeable to
the node already being processed. The same gate can inspect a node's
numerator and denominator together.

More explicitly, if reduction is requested for precision `q`, the gate
implies `max(bits(n),bits(d)) > q floor(log2 p)`. The modulus therefore
has at most a fixed constant times that input size. Both its construction
and the ensuing division fit the node's multiplication-scale work. A
reduction performed by this rule does not increase the numerator or
denominator bit length, even for a negative numerator. The ordinary
unreduced product/shift bounds can consequently be used inductively for
all later levels. Generic prime-power reduction does not introduce a
full-width operation at every small leaf.

After each forest is reduced to a root, it may be combined with a
single running rational root, with numerator and denominator capped at
the applicable p-adic precision. This adds `O(u M(N))` work and preserves
linear storage. One final p-adic inverse evaluates the accumulated root.

Combining the terms, for `u=Theta(sqrt(L))` the exact core construction
has the bound

\[
 O\left(NuL^2+NL^3/u+NLu^2+NL^2(\log L)^2\right)
 =O(NL^{5/2}),
\]

with `O(N)` tree/record data plus the separately accounted arithmetic
workspace. Keeping integer corrections separate through balanced
merges is compatible with these bounds. Adding a full-width correction
to a global residue once per block is not the justified implementation.

## 6. Implementation mapping and scope

| Requirement | Production source inspected |
| --- | --- |
| Local p cancellation and omitted CRT p-part | `weighted_machin_block`, local power and `if(p==prime) continue` |
| Exact signed correction | `machin_real_approx`, integer remainder comparison, signed `floor_div` |
| Non-floating count upper bound | Integer `target_bits`, `geometric_bits`, and `ratio_bits` schedule |
| Small-A blocks | Uniform reverse-r traversal; covered by the aggregate argument above |
| Bounded resident block data | `mass_limit`, `pending`, `flush`, and one `accumulated` root |
| No unconditional giant leaf moduli | `trim_pair` performs the integer bit-length gate first |
| No unbounded power cache | Shift powers and trimming moduli have local lifetimes |
| Vanishing block handling | `shift >= m` is skipped before precision arithmetic |

Storing all block roots and caching arbitrary
large powers, including moduli built before checking whether reduction
is needed, would not meet the stated storage/work argument;
successful numerical comparisons alone do not repair that discrepancy.

The Machin word loop and interval certificate are separate from this
core analysis. At `G=O(log N)`, one pass of the ordinary word arithmetic
fits below the core bound. The program's adaptive guard retries need
their own worst-case accounting, as do exact word aggregation and any
claim of a complete extractor in the manuscript's model. A successful
integer interval certificate proves the returned digits; it does not
by itself prove that all inputs finish within `O(N log^(5/2) N)` time.

## 7. Prospective exact-word fallback lemma

This section describes a proposed extension, not an implemented mode of
the inspected program. Its purpose is to remove the `N*G` cost of
rounding all word fractions individually when the guard is linear in
the requested position. The current implementation and measurements
retain the one-pass scope above.

Fix one Machin family with base prime `p0`, coefficient `c`, numerator
base `b`, cancelled exponent `s`, and outside sign `sigma`. Thus
`(p0,c,b,s,sigma)` is `(5,16,2,n,+1)` or `(239,4,10,0,-1)`, where
`n=N-1`. Write

\[
 \mu=2k+1=p_0^t u,\qquad p_0\nmid u,
\]
\[
 r_k\equiv c b^n p_0^{s-\mu-t}\pmod u,
 \qquad 0\le r_k<u.
\]

The word contribution is `sigma*(-1)^k*r_k/u`; the case `u=1` contributes
zero modulo integers. Let `B0=2K+1`. For each odd prime `l != p0`, use
the global power

\[
 Q_l=l^{\lfloor\log_l B_0\rfloor}.
\]

At an incident denominator with `v=v_l(mu)`, its CRT coefficient can be
computed directly as

\[
 \frac{Q_l}{l^v}\,\sigma(-1)^k c b^n
 p_0^{s-\mu}\,(\mu/l^v)^{-1}\pmod {Q_l}.
\]

To verify the formula, first apply the usual CRT coefficient to `r_k/u`.
Replacing `r_k` by its defining modular expression is valid after
multiplication by `Q_l/l^v`. Moreover

\[
 (u/l^v)^{-1}=p_0^t(\mu/l^v)^{-1}\pmod {Q_l},
\]

which cancels the exponent `-t`. All displayed inverses are of units:
`l != p0` and `mu/l^v` is not divisible by `l`.

For fixed `l`, enumerate `mu=l,3l,5l,... <= B0`. Successive indices
`k` differ by the odd number `l`, so the modular weight advances by
the constant factor

\[
 -p_0^{-2l}\pmod {Q_l}.
\]

The cofactor inverses can be computed individually, or in bounded
chunks using batch inversion. A chunk of `O(L^2)` word records requires
only `O(L^3)` bits; it need not be retained after updating the current
coefficient. Storing all term incidences is unnecessary.

Let `c_l` denote the completed coefficient and form

\[
 T=\sum_{l\ne p_0}c_l/Q_l.
\]

Then `T` equals the word sum modulo the integers. No integer-correction
recovery is needed when only the fractional part is required.

The product of the `Q_l` divides the odd part of
`lcm(1,...,B0)`, so it has `O(B0)` bits. Also every `Q_l > sqrt(B0)`:
for `l <= sqrt(B0)`, maximality gives `Q_l > B0/l`, and the larger primes
already exceed that threshold. Hence the completed fractions and their
metadata occupy `O(B0)` bits. A balanced rational tree computes their
sum in `O(B0 L^2)` bit operations and `O(B0)` tree data. The two families
can then be combined by one rational addition; shared prime factors
between their denominators do not invalidate that addition or the
linear bit-length bound.

For a sequential bit-cost implementation, use the manuscript's
segmented prime enumeration: generate composite marks in segments of
length `ceil(sqrt(B0))`, sort those records, and stream the unmarked
primes. This costs `O(B0 L^2 log L)` work and
`O(sqrt(B0)*polylog(B0))` resident bits. The total number of denominator
incidences and successive valuation divisions is `O(B0 log L)`.
Even independently regenerating the modular weight and inverse at
each incidence gives the conservative bound

\[
 O(B_0L^2(\log L)^2),
\]

below the weighted core bound. The recurrence and chunked inverses may
improve this arithmetic count, but no lower unit-cost RAM estimate is
needed for the claim here. Prime generation, record routing, metadata,
and large-integer arithmetic have all been included.

After reducing the resulting numerator modulo its denominator, floor
the single exact word fraction at precision `G`. This requires one
division with `O(B0+G)`-bit operands, rather than `Theta(B0)` separate
G-bit divisions. Together with at most two core roundings and the
combined analytic tail, the interval certificate can use `E=4` ulps.
For `G=Theta(N)` and `B0=Theta(N)`, the word fallback has linear data
storage and the stated subdominant bit cost.

### A necessary guard-scheduling condition

An exact word fallback alone does not establish the complete adaptive
bound. Recomputing an `O(N L^(5/2))` core for every doubling from
`G=O(log N)` to `G=Theta(N)` permits an extra `Theta(log N)` factor.

A prospective complete construction can use at most a fixed number of
cheap pilot attempts and, upon failure, jump directly to

\[
 G=32(n+w+2),
\]

where `w` is the fixed output width. Use the exact word fallback on
that attempt. The manuscript's cited result `mu(pi)<8` implies, for
all sufficiently large `N`, distance at least
`delta_N=10^(-7n-8w)` from `10^n*pi` to a width-w decimal boundary.
The displayed guard satisfies `8*2^-G < delta_N`. Since `E<=4`, the
certification interval must then fit inside one decimal cell.

The first fallback consequently succeeds for all sufficiently large
positions. For the finite exceptional positions, doubling the fallback
guard and recomputing terminates by irrationality. This is the same
finite-exception convention used in the manuscript; it does not
require knowing the irrationality result's implicit threshold.

For a fallback with `G` a fixed constant times `N`, the core has
`K,m,R=Theta(N)` with different constants. The uniform-block aggregate
argument still applies. These ingredients therefore supply a route to
an end-to-end bound without an `N*G` word-rounding cost or a logarithmic
number of full core recomputations. That route still requires actual
implementation, validation, and a separate complete theorem statement;
it is not a property being assigned to the current native executable.

The current native limits also need attention before that extension:
the proposed linear guard can exceed the command-line guard limit, and
its larger cutoffs can exceed the checked 32-bit word moduli. The
prospective construction requires corresponding representation changes.

The separate `machin_exact_words_checks_2026-10-24.py` validates this
word-fraction identity in 130 exact whole-fraction cases, 520 projections
through 4096 guard bits, and 24 deep prime-cancellation cases. Its saved
JSON report covers that component only, not the unimplemented adaptive
extractor.
