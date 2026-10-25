# Recursive Lucas reconstruction for the positive Machin correction

24 October 2026. Executable design and cost mapping, with explicit scope.

This design replaces the weighted Machin real-correction stage. It does
not replace the word stage or compute a prefix of pi. Its exact recursive
solver is implemented in `machin_lucas_reconstruction.hpp`, using the
shared coefficient kernel in `machin_lucas_coefficients.hpp`. The native
prototype and its integration must pass exact tests and matched timings
before adoption. The cost mapping below specifies a sequential schedule;
it is not automatically a proof that every native container access has
that schedule.

The recursive construction can have many finite-size reconstruction
levels, particularly for p0=239. The current experimental header also
contains the direct-H, fixed-prime variant described in Section 8.
That variant removes those levels, but its globally capped fixed-prime
kernel does not yet have the stronger asymptotic cost bound. The
production weighted Machin engine retains its previously audited
log^(5/2) one-pass core bound. This experimental record neither replaces
that production bound nor claims O(N log^2 N).

## 1. Exact auxiliary problem and coefficients

Fix p0=5 or 239, D=p0^2, C=D+1, odd A>0, and fixed power-of-two a>=2.
Let

    h_j(A)=a^j j!/[C^j product_(i=1)^j(A+a i)],
    H_r(A)=sum_(j=0)^(r-1) h_j(A),
    F_r(A)=C^(r-1) H_r(A),
    E_r(x)=[C^r-(1-x)^r]/[D+x]=sum_(k=0)^(r-1)e_k x^k.

The polynomial quotient is integral. Partial fractions give

    F_r(A)=C^(r-1)-a sum_(k=1)^(r-1) k e_k/(A+a k).       (1)

This identity is exact. In particular, F_r has denominator dividing
L_r(A)=lcm(A+a,...,A+a(r-1)). Every prime dividing this progression,
including p0 and the prime 13 dividing C, belongs in its CRT registry.

For each such prime q, let q^E be its maximum power in the progression.
At an occurrence mu=A+a k=q^v u, use

    c_q += q^(E-v) (-a k e_k) u^(-1) mod q^E.             (2)

Only e_k modulo q^v is needed. Take the final coefficient in [0,q^E)
and put T_r=sum_q c_q/q^E. The standard CRT comparison proves

    Z_r=F_r-T_r is an integer.                           (3)

This is the integer that must be recovered; fractional CRT residues
alone do not determine the desired exact real correction.

## 2. Coefficient generation, including every nonunit case

Set b_k=(-1)^(k+1)binomial(r,k), including b_0=-1. Then

    e_0=(C^r-1)/D,
    D e_k+e_(k-1)=b_k.                                  (4)

For a small query prime q!=p0, generate a full binomial row modulo the
largest needed q-power using valuations and unit parts. Apply (4)
forwards using the inverse of D. A binomial coefficient becoming zero
does not discard its unit state. Precompute the short list of q-powers
instead of exponentiating at each row entry.

For q=p0, D is a nonunit. Generate the coefficient row backwards:

    e_(r-1)=(-1)^(r+1),
    e_(k-1)=b_k-D e_k.                                  (5)

This branch is mandatory even when a finite small-prime threshold is
below 239. No exceptional branch is needed at q=13: D is a unit there,
and no formula below inverts C.

For larger primes use the prefix

    P_k=sum_(j=0)^k binomial(r,j) D^j,
    e_k=(-1)^k D^(-k-1)(C^r-P_k) mod q.                 (6)

Writing r=r1*q+r0 and k=h*q+k0, Lucas' theorem gives

    P_k=C^r0 * sum_(j<h)binomial(r1,j)D^j
          +binomial(r1,h)D^h
             *sum_(j<=min(k0,r0))binomial(r0,j)D^j       (mod q).  (7)

All occurrences q|A+a k in a block have the same k0. Thus only one low
prefix is needed per block and prime; the high indices h are consecutive.
For q<=sqrt(r), generate low/high rows shared across blocks. For
q>sqrt(r), obtain the one low coefficient e_k0 from the shared value
tree, then recover its prefix with (6). Here k0<q. If q>r, there is at
most one occurrence and its tree answer is the required coefficient.
The multiplier (-1)^k D^(-k-1) changes by -D^(-1) between successive
occurrences spaced q apart. Missing this minus sign is an implementation
error. Formula (7) is a prime-modulus formula, not a prime-power formula.

At occurrences q^2|A+a k with q above the small-prime cutoff, obtain
e_k modulo the actual q^v through the sparse exact-division tree. The
ordinary low request may still be useful for the other occurrences.

### Shared exact transfer and sparse padding

The transfer from (e_(k-1),b_(k-1)) to (e_k,b_k) has numerator

    X_k = [ -k,  -(r-k+1) ]
          [  0, -D(r-k+1) ]

and denominator D*k. Normalize an interval of s transfers by
G_s=s!/lcm(1,...,s). The proof in
`weighted_machin_reconstruction_derivation_2026-10-24.md` establishes
that G_s divides all three numerator entries and the denominator, and
that G_(l+r)/(G_l G_r) is an integer O(l+r)-bit merge factor.

For q!=p0 the normalized denominator has valuation at most
2 floor(log_q r). This permits lambda_q=2 floor(log_q r) extra digits
per remaining tree level. A sparse source modulus includes the actual
multiplicity of repeated prime queries. To enter a right child, form
the relevant inflated source S, target M, and normalized transfer U/Q.
Compute g=gcd(Q mod S,S), W=U V mod S, divide W and Q exactly by g,
and multiply by (Q/g)^(-1) modulo M. The padding ensures M|S/g and
gcd(Q/g,M)=1. Reduce Q/g before the inverse; a repeated large-operand
gcd or inverse is not charged as a small-modulus operation.

Ordinary queries have k0<q and q!=p0. Every applied transfer index is
at most k0, so its raw denominator D^k0*k0! is a unit modulo q.
Normalizing cannot add prime factors. Thus those queries need no
sparse padding. Leaves return the incoming state before their own
transfer; query index k must therefore select coefficient e_k.

## 3. Exact recovery from a shorter positive Euler sum

Use a common r for all blocks in a batch. For each block compute

    d_A=floor(log2(C A/(a r))),  d=min_A d_A,
    c0=ceil(log2 C),
    s=ceil((c0*(r-1)+18)/d).                             (8)

The logarithms in (8) are integer comparisons, not floating estimates.
If d<=0, s>=r, or r is below the chosen terminal threshold, evaluate
H_r directly by exact binary splitting. Otherwise recursively obtain
H_s for the same blocks and reconstruct using (3).

For j<=r, h_j/h_(j-1)<=a r/(C A)<=2^-d. Independently all positive
ratios are less than 1/C. Therefore

    0<=C^(r-1)(H_r-H_s)
       <C^(r-1) 2^(-d s) C/(C-1)<2^-17.                (9)

In particular, the exact rational number C^(r-1)H_s-T_r is within
2^-17 of the integer Z_r. Rounding it to its nearest integer recovers
Z_r uniquely; no separation hypothesis about a fractional part is
needed. Then return

    H_r=(T_r+Z_r)/C^(r-1).                              (10)

The implementation uses exact rational arithmetic for that rounding.
A fixed-point version with enough guard bits for the number of prime
fractions is also valid, but is not required for the prototype. The
prototype gcd-reduces the resulting rational; that cost is included
in the reconstruction budget.

This is a strictly shrinking recursion whenever invoked. A fallback
terminates all finite cases, including A close to 1. The asymptotic
halving argument requires a sufficiently large parameter regime; for
p0=239 the constant c0=16 makes a literal halving threshold much more
conservative than for the older 15-series. Correctness does not require
halving, and small finite cases must not be described as if they do.

## 4. Root truncation and small-A handling

Let a weighted outer block have length t, local multiplier P, and
real scale P*p0^(2t)/(C A). For a common maximum t and P, choose an
integer target T>=log2(P*p0^(2t))+12. Start at

    r0=ceil(T/floor(log2 C)).

For any previously certified upper count r_old, the ratios of all
terms through a smaller count are bounded using

    d=max(floor(log2 C),floor(log2(C Amin/(a r_old)))),
    r_new=min(r_old,ceil(T/d)).                         (11)

Repeat while r decreases. Every step is safe because its ratio bound
was established up to the previous larger count. At the final count,
h_r<=2^-T; the geometric tail and the factor 1/(C A) give absolute
scaled error below 2^-12. This is stronger than the required outer
integer-recovery margin. The other omitted outer endpoint retains
its existing bound and must still be included there.

Set L=ceil(log2 N), H=log2(L+2), B=Theta(L^2/H), and
ell=Theta(NH/L^2). The 239 family reaches A=1. Keep the original
exact/direct weighted method on blocks with A below a constant
multiple of N/L, including a boundary block. These blocks span
O(N/L+ell)=O(N/L) terms. Their crude real work O((N/L)L^3)=O(NL^2)
fits the prospective target, and their raw largest real product has
O((N/L)L)=O(N) bits. No pi-prefix computation is introduced.

For the remaining batch, Amin=Omega(N/L). At the initial r0=O(ell),
the useful ratio in (11) is Omega(L/H), so one improvement already
gives a common root R=O(ell/H)=O(N/L^2). Its scaled tail is proved by
(11), not by assuming that all A are Theta(N). The child lengths in
(8) eventually at least halve. Consequently

    sum_j r_j=O(R),
    sum_j r_j[1+log(N/r_j)]=O(RH).                     (12)

The finite constants can use direct fallbacks. A timing experiment
must report when those fallbacks bypass the new algorithm.

## 5. Full cost and stored-data schedule

This mapping uses M(x)=O(x log x), general inverse/gcd cost
O(M(x)log x), fixed p0,a, and the existing segmented factorization,
short-progression LCM, and sequential sorting lemmas. It is a mapping
of those stated lemmas, not an independent replacement for them.

Use the smaller-parameter prescription

    Y=Theta(L^3/H^2), terminal=Theta(sqrt N),
    forest root quota q0=Theta(N/L).

Fixed positive rescalings of these parameters preserve the orders;
a fixed terminal threshold instead of Theta(sqrt N) does not justify
the metadata estimate below. Let B bound the number of regular blocks.

1. **Shared coefficient trees.** The normalized matrix tree at length
   r has O(r log^2 r) retained bits and O(r log^3 r) construction work.
   Summing (12) gives O(RL^2)=O(N) bits and O(RL^3)=O(NL) work.
   The tree is built once per length, not separately for each block.

2. **Small-prime rows.** A full row for one q<=Y costs O(r L H^2)
   with valuation/unit updates and cached q-powers. Summing the shared
   rows over lengths and primes gives O(RY L H^2)=O(NL^2).
   The forced p0 row adds only a fixed number. Shared medium-prime low
   and high row lengths total O(RH), below this budget.

3. **Ordinary requests.** One request per block and large prime has
   total modulus mass sum W_j=O(BRH)=O(N), by the progression LCM
   bound. A nonterminal r>c sqrt N has large primes q>sqrt r with
   Omega(L)-bit lengths, so O(L)-bit query metadata is chargeable to
   that mass. A bounded forest costs O(W L^2+rL^3), including value
   reductions. Even charging a complete shared-value traversal to
   each forest gives O(RL^3(1+N/q0))=O(NL^2). Retained modulus trees
   have O(q0 L)=O(N) bits.

4. **Higher-power requests.** The number of occurrences with q>Y
   and q^2|A+a k is O(BR/Y)+sqrt(N)L^O(1). Each query's inflated
   modulus has O(L^2) bits, including the doubled lambda from G_s.
   Hence total inflated mass Q=O(NH/L)+sqrt(N)L^O(1). Sparse
   traversal costs O(Q L^3)=O(NL^2H); repeated value traversal adds
   O(RL^3(1+Q/q0))=O(NLH) plus sublinear terms. Again retain only one
   O(N)-bit forest. Repeated primes in product moduli do not require
   coprimality, but do require the stated precision multiplicities.

5. **Incidence arithmetic and routing.** Across recursive lengths,
   BR=O(N/H), and all factor incidences number O(BRH)=O(N), up to
   sqrt(N)L^O(1) boundary overhead. Charge O(L M(L))=O(L^2 H) per
   occurrence; this includes fixed-many word-modulus exponentiations
   and all power stripping, not unit-cost modular arithmetic. Sort
   the O(L)-bit records and stream the rows and block states, for
   O(NL^2H) total work. Local factor records and per-block maxima
   have the earlier segmented-factorization bounds. Storing every
   exact coefficient instead would use Theta(r^2) bits and is invalid.

6. **Rational reconstruction.** The coprime prime-power leaf mass of
   T_r(A) is O(r[1+log(N/r)]). All prime-fraction trees therefore cost
   O(BRH L^2)=O(NL^2). Returned H_r rationals have the same order of
   size, including C^(r-1). Canonical gcd reduction fits this budget.
   In (8), the child rational has O(r+L) bits since Amin>=cN/L and
   d is comparable to 1+log(N/r) in the regular regime. The scale,
   rounding, and final division fit the same arithmetic accounting.
   Retaining parent fractions while recursing uses O(BRH)=O(N) bits.

7. **Terminals and original blocks.** Direct terminal splitting costs
   B sqrt(N)L^3=sqrt(N)L^O(1), asymptotically below the target. Original
   weighted block CRT trees retain their O(NuL^2)=O(NL^2H) bound for
   u=Theta(H). LCM assembly is no larger than the earlier bounds.
   The small-A direct part was charged in Section 4.

These give a complete abstract construction for the weighted core with
O(NL^2H) work and O(N) stored bits plus arithmetic workspace, conditional
on the explicitly reused factorization/routing lemmas. For the actual
decimal extractor, the ordinary word stage, fixed one-pass guard
schedule, and output certificate must be included with their own
proved costs. This note does not assert a full adaptive theorem.

## 6. What could invalidate a claimed improvement

- Running a length-r coefficient recurrence separately for every block
  and prime loses the sharing that pays for the construction.
- Routing every prime occurrence through the padded tree gives inflated
  mass O(BR L^2), instead of restricting it to sparse higher powers;
  its traversal can exceed the target by extra powers of L.
- Treating p0 as an ordinary prime attempts to invert D=p0^2 and is
  incorrect. Applying Lucas directly modulo q^v for v>1 is also
  incorrect.
- Setting Amin=1 for a shared root prevents the required shortening
  and can restore O(NL^3) work. The small-A splice is part of the
  algorithm, not an optional proof convenience.
- Keeping all query-product trees, or using root quota Theta(N), can
  retain O(NL) bits. Retain one Theta(N/L)-mass forest at a time and
  free its old GMP allocations when advancing.
- Using base=O(1) in the asymptotic proof makes q>sqrt r insufficient
  to charge large-query metadata to Omega(L)-bit prime lengths.
- A native vector lookup is not by itself a sequential routing proof.
  The abstract schedule sorts low-index requests, streams full rows
  with block states, and sorts sparse exceptions to merge them into
  incidence scans. The prototype may use random indexed vectors for
  finite performance; that correspondence must be audited separately.
- Removing all original block boundaries leaves a global integer carry
  after CRT aggregation. Its recovery is a further problem. The
  present construction retains an O(log log N) factor and does not
  prove strict O(N log^2 N).

This is a concrete recursive algorithm, with exact nonunit handling
and a shrinking auxiliary problem. Its parameterized abstract bound
is stronger than log^(5/2), but neither that bound nor the new local
lemmas predict a practical speedup. The full integrated prototype must
be measured against the current weighted Machin implementation.

## 7. Audit of the integrated one-pass schedule

The first integrated `decimal_extract_weighted_machin_lucas.cpp` uses
B rounded to a power of two above L^2/H, excludes block starts below
Amax/L from the shared real batch, uses base=max(64,sqrt(Amax)/16),
Y=max(8,ceil(L^3/(8H^2))), and quota=max(4096,Amax/L). Here Amax and
the number of retained core terms are each Theta(N) for the fixed
Machin families under G=O(log N). Rounding and fixed rescalings do not
alter the abstract orders. The excluded intervals are contiguous;
one boundary block adds at most ell, which is O(N/L) in the stated
asymptotic regime. The regular batch's minimum A and root count
therefore satisfy Section 4, for fixed C=26 or 57122.

The existing ordinary word loop need not be replaced by exact global
word aggregation for this one-pass scope. There are O(N) terms, and
each simultaneous powering loop uses O(L) modular multiplications on
O(L)-bit quantities. Under fast bit arithmetic this is O(N L M(L))
=O(NL^2H). Stripping powers of the fixed Machin prime costs no more.
For G=O(L), each guard-scaled division has O(L)-bit operands, and the
signed accumulator has O(G+L)=O(L) bits. Its rounding, accumulation
and final reduction are within that bound. Powers 10^(N-1), the final
prime-power moduli, LCM assembly and final inverses have the previously
audited large-integer costs. The successful integer output certificate
is unchanged. Arbitrarily many larger-guard retries are not covered.

There remain implementation-to-model gaps, distinct from correctness:

- Native residue moduli and progression endpoints are explicitly
  bounded by UINT32_MAX. The multiply/remainder primitives and ordinary
  signed-word Euclidean inverse implement that range, not arbitrary
  L-bit words. An unbounded theorem requires an explicit replacement
  by the stated fast modular arithmetic algorithms. In particular,
  the small-row O(r L H^2) estimate assumes fast inversion or a growing
  inversion batch size; a fixed native batch of 256 is sufficient for
  the supported 32-bit range but does not itself establish an
  unbounded fast-inverse theorem.
- The prototype retains indexed vectors for full-row lookups, group
  destinations and query fan-out. Its stored-data bounds are compatible
  with Section 5, but the full sequential record movement there has
  not been proved identical to these native accesses. The explicit
  sorted/streamed schedule is the algorithm covered by that accounting.
- Cutoff and square-root candidates use native floating arithmetic
  followed by exact correction. Exact tests establish their resulting
  values. In an unbounded extension, a bounded number of correction
  steps cannot be inferred from fixed-precision floating arithmetic;
  use integer selection or charge a provably bounded search. For
  example, binary search with O(L) fixed-base exact powers has
  O(M(N)L) cost, within the prospective target.
- A fixed finite streamed exact coefficient backend is selected when
  r<=16384. It does not retain the dense exact row and falls back to
  the normalized value tree above the cutoff. That preserves the
  asymptotic construction as a finite exception. Letting that cutoff
  grow without accounting for the exact row's quadratic bit work would
  not preserve the proof.

Thus the source supports an implemented and tested exact algorithm,
and an explicit candidate schedule with a stronger abstract bound.
It is not appropriate to label the finite executable itself as having
a newly established unbounded theorem on the basis of these tests.

## 8. Removing the real recursion by reconstructing H directly

For p0=5 or 239, C=2*13^e, where e=1 or 4. After canceling the fixed
factor 2 in each Euler ratio, write

    h_j=(a/2)^j j! / [13^(e j) product_(i=1)^j(A+a i)].

For each prefix of length j and each prime power 13^k, the count of
progression multiples differs from floor(j/13^k) by zero or one.
Consequently, with V=floor(log_13(A+a(r-1))),

    0 <= v13(product_(i=1)^j(A+a i))-v13(j!) <= V,
    denominator_13(h_j) <= 13^(e j+V).

Hence J=e(r-1)+V is a safe common exponent for the 13-primary part of
H_r. At every other odd prime q, simply multiply the coefficient in
(2) by the inverse of C^(r-1) modulo its local maximum q-power. The
coefficient precision lemma still applies. The existing progression
CRT groups for q=13 are replaced by one fraction

    c13/13^J,  c13 = 13^J H_r mod 13^J.

The expression defining c13 is evaluated in the 13-adic sense: all
remaining denominator factors are units and are inverted modulo 13^J.
After canonical cancellation of powers of 13, the denominator can be
smaller; the represented fraction is unchanged modulo the integers.

Let T be the sum of this primary fraction and every other prime
fraction of H_r. Their CRT identities prove H_r-T is an integer.
The positive real ratios are all below 1/C, so

    1 <= H_r < C/(C-1),
    |H_r-1| < 1/(C-1) < 1/2.

Therefore round 1-T to the nearest integer Z and return H_r=T+Z.
No shorter real sum, large C scaling, or recursive integer recovery
is required. The prime fractions use pairwise coprime denominators;
their product has O(r[1+log(N/r)]) bits, including the O(r+L)-bit
13-power. A final generic gcd is unnecessary for this size bound.
The prototype keeps recursive handling for C/2 that is not a pure
power of 13, including its p0=3 test fixtures.

### Exact global-cap valuation semiring

`machin_lucas_adic.hpp` computes the fixed-prime component without
forming the full exact rational. Store each raw PQT entry as a unit
coefficient times a known power of 13:

    P=13^vp * p, Q=13^vq * q, T=13^vt * t.

The word "unit" is exact for p and q; t may retain extra powers of
13. All coefficients are kept modulo 13^J. Leaves strip their known
13 factors after the fixed leaf gcd cancellation. Merge by

    vp=vpL+vpR, vq=vqL+vqR,
    vt=min(vtL+vqR, vpL+vtR),
    p=pL*pR, q=qL*qR,
    t=13^(vtL+vqR-vt)*tL*qR
        +13^(vpL+vtR-vt)*pL*tR             (mod 13^J).

Every alignment exponent is nonnegative. Thus coefficient reduction
modulo 13^J is preserved by every later merge. No exact division of
an already truncated coefficient is performed. At the root, the
stored pole K=vq-vt is at most J by the prefix valuation bound above,
even when t retains additional 13 factors. If K<0 the fraction is
13-integral and contributes zero; otherwise

    c13=13^(J-K)*t*q^(-1) mod 13^J.

This proves the global-cap algorithm's exactness. It applies to all
fixed power-of-two steps a used by the solver, because gcd(a,C)=2 and
the remaining numerator factor a/2 is a 13-adic unit.

### The unresolved cost of the fixed-prime evaluation

The global cap bounds each coefficient by O(r+L) bits. It does not
prove that a subtree of s leaves may be capped to O(s+L) bits. In an
affine reverse-prefix merge R=R_left+13^v U_left R_right, the unshifted
left child still needs the full parent precision. In the reversed
orientation the corresponding issue moves to the other child. Unit
denominators have arbitrarily high 13-adic digits; discarding them is
not justified by cancellation of the valuations alone.

Using only the established bounds, the lower value-tree levels can
still have O(rL) aggregate bit mass and O(rL^3) arithmetic work per
block. At B=Theta(L^2/H), that yields the conservative O(NL^3/H)
bound for this new stage, rather than O(NL^2H). It can still improve
finite timings by eliminating repeated levels, but cannot inherit the
recursive route's stronger bound without a new fixed-prime algorithm.
The implementation selects this direct-primary path only for r<=16384.
Above this fixed finite cutoff it uses the recursive algorithm; a child
can enter the primary path after becoming small enough. This preserves
the conditional asymptotic fallback. In an unbounded model a fixed-r
primary call still processes O(L)-bit origins, so it is not literally
constant cost; its O(L)-sized arithmetic per block is below the stated
O(NL^2H) budget. Raising the cutoff as a function of N would require
charging the new fixed-prime work explicitly.

On 24 October the independent native checker passed 378 complete H_r
identities: 126 direct-primary cases for p0=5,239 and 252 recursive or
p0=3 fallback cases. These include r through 257, A near 1 and near
UINT32_MAX, and steps 2,4,8,16. A separate fixed-prime oracle passed
594 comparisons with exact rational splitting through r=500. The
coefficient checker also passed 12699 coefficient checks, 2241 backward
prime-row checks, 3738 ordinary and 7266 sparse queries, 3084 shared-leaf
queries, and 12288 weighted-prefix checks. These counts establish
correctness evidence for the checked implementations, not a runtime
claim or a proof of the unresolved unbounded cost.

### Native leaf batching and nested inverse audit

The wrapper rejects C<2, nonpositive A or r, an even A, a step that is
not a power of two at least 2, and progression endpoints or a*r beyond
UINT32_MAX. The selected base must actually be prime; production uses
13. After dividing a and C by their gcd, the remaining C must be a
positive power of that prime. Thus neither stripping loop receives
zero, and the allocated precision is positive for r>1. The input range
also makes its exponent and valuation arithmetic fit u64.

For a positive leaf, pj/qj=a*j/[C(A+a*j)]<1/C. A positive interval
therefore has T/Q<1/(C-1); one containing index zero has T/Q<C/(C-1)
<=2. With largest bounding every qj, batching s leaves has P,Q at most
largest^s and T<2*largest^s. The source's respective bounds
largest<2^62,2^41,2^31 for s<=2,3,4 therefore put T below 2^125,
2^124,2^125. All native products and sums fit the unsigned 128-bit
temporary. The fixed direct-primary integer recovery still uses
C=26 or 57122, for which C/(C-1)<3/2; the weaker C>=2 wrapper guard
alone is not a proof of that nearest-integer recovery.

For inversion modulo 13^J, retain the ceil-halving precision list and
construct ascending moduli by squaring, dividing exactly by 13 when
the target exponent is odd. Reduce the input once at the largest
modulus, then descend through these nested moduli. Starting with an
inverse modulo 13, the update x <- x*(2-q*x) modulo the next modulus
doubles its valid exponent, covering every next target. The retained
moduli, residues and inverses form geometric bit-size series. Under
the usual regularity assumptions on M, this is O(M(J)) work and O(J)
bits, including all reductions. It avoids repeatedly reducing the
full-size denominator independently at each precision. This is a
read-only proof audit of the native batching and inverse changes;
their regression test and benchmark reports belong to the integration
run rather than being inferred from the earlier test counts above.

## 9. Exact separation of the outer integer corrections

Write each original weighted block exactly as F_b=T_b+Z_b, where
T_b=n_b/d_b is its sum of prime fractions and Z_b is an integer. The
current assembly first forms n_b+Z_b*d_b. Instead, for the existing
nonnegative outer shift s_b and sign epsilon_b, one may accumulate

    U=sum_b epsilon_b*p0^s_b*T_b,
    V=sum_b epsilon_b*p0^s_b*Z_b,
    result=U+V modulo the final p0-power.

This is an exact separation, not a second approximation. Each local
denominator remains coprime to p0. A balanced integer tree represents
V relative to the minimum shift in its interval, restoring differences
with exact powers of p0. That permits one final combination of V with
the assembled rational U. It avoids the individual products Z_b*d_b.

At length ell, the local LCM denominator has
O(ell[1+log(N/ell)]) bits, whereas Z_b has O(ell+L) bits: the geometric
weighted block has logarithm O(ell+L), and T_b is bounded by its number
of prime fractions. Thus the integer side can have smaller local
operands. A balanced integer assembly has O(N) total operand mass
per level, up to negligible metadata terms, and O(M(N)log B) work.

This does not remove the asymptotic height of the rational assembly.
At a node spanning s adjacent blocks, the p0-shift range is O(s*ell),
and its LCM denominator can have O(s*ell[1+log(N/(s*ell))]) bits. Both
terms remain in U even when the local Z_b*d_b products are absent.
Near the root the geometric shift span alone is Theta(N). Summing
the rational merge sizes still permits the previous O(NL(log B)^2)
bound. Separation can reduce constants and local operand sizes; it
does not justify dropping a logarithmic factor from that bound.

The auxiliary denominators A+a*k for 1<=k<r coincide with a prefix of
an original block whenever r<=ell and the same A,a are used. Under
that condition, a future pipeline can reuse that block's factorization
to create its auxiliary coefficient requests. Factor a short extension
separately if r>ell; do not assume the prefix relation in every finite
case. To maintain the space bound, stream T_b into the existing global
LCM stack immediately, retain only bounded-error fixed-point sums for
its later integer recovery, and assemble Z_b after the shared H_r
answers arrive. Holding all original exact block fractions until those
answers arrive can retain O(N log B) bits. This pipeline has not been
implemented in the audited prototype.

## 10. A precise next target for fixed-prime acceleration

There is a useful reverse representation. Put alpha=A/a in the
13-adic integers and

    R_r(z)=sum_(t=0)^(r-1) a_t*z^t,
    a_t=binomial(r-1+alpha,t)/binomial(r-1,t).

Then H_r=h_(r-1)*R_r(C). Its coefficients satisfy
(r-t)*a_t=(r+alpha-t)*a_(t-1). Direct coefficient comparison gives the
exact inhomogeneous differential equation

    z*(z-1)*R'_r+[r-(r+alpha-1)*z]*R_r
      = r-alpha*a_(r-1)*z^r.                            (13)

The homogeneous solution is z^r*(1-z)^(alpha-1). The z^r forcing term
is essential: the formal Gauss expression has bottom parameter 1-r,
a nonpositive integer, so its ordinary infinite hypergeometric series
is singular at the very next coefficient. Although C is 13-adically
small, the valuation of that forcing term is not uniformly above the
required J=e(r-1)+floor(log_13(maximum)). It cannot simply be omitted.

The primary-source algorithms of Caruso, Mezzarobba, Takayama and
Vaccon, [Fast evaluation of some p-adic transcendental functions](https://arxiv.org/abs/2106.09315), give fast binary-splitting and
digit-burst evaluation under explicit coefficient-height and local
solution hypotheses. Their fixed-parameter hypergeometric conclusions
are not a drop-in bound here: r and A grow, the forward argument 1/C
is outside the 13-adic unit disk, and the reverse equation is resonant.
Applying that framework requires charging the growing parameter height
and supplying the appropriate finite solution data. This paragraph
describes the relevant scope of that source, not a theorem obtained
from it for this extractor.

A concrete research target is therefore an exact evaluator of (13),
or an equivalent 13-step transfer, at precision O(r+L) with aggregate
cost below O(rL^3) for these varying A,r. A 13-step product isolates
one divisible factor from each consecutive numerator and denominator
progression, but leaves twelve nondivisible factors on each side.
Their unit product depends on the block index and A to full precision;
replacing it by a constant Frobenius factor is not valid. Any proposed
compression must bound the work to evaluate these unit corrections.

Likewise, dense polynomial multipoint evaluation in A is not yet a
solution: degree-r product coefficients can require O(r^2 log r) bits.
Reducing first modulo a degree-B polynomial of block origins limits
root storage, but the aggregate coefficient height at intermediate
product levels still needs a separate proof. A small root object alone
does not bound the whole tree. These are explicit remaining tasks,
not an implemented O(N log^2 N) algorithm.
