# Exact reconstruction reductions for weighted Machin extraction

24 October 2026. Exact algorithmic reduction and implementation accounting.
The real-PQT and coefficient-transfer reductions below are exact experimental
constructions. They are not adopted in the production real-series kernel and
do not by themselves prove an exponent smaller than the existing weighted bound.
The real-tree variants below are experimental mathematical constructions;
none was selected for production because the measured overhead outweighed
its savings. The half-inverse reconstruction experiment is a separate
change. Refer to the final matched benchmark report for adoption and
timing results.
The algebra allows any positive odd A and power-of-two a>=2. The stated
asymptotic costs keep p and a fixed, as in the production Machin cores
with p=5 or 239 and a=2, and use M(x)=O(x log x).

## 1. Cancel the fixed factor two before splitting

The positive real-series ratio for index j>0 is

    a*j / [C(A+a*j)],    C=p^2+1.

For the tested odd primes and power-of-two a>=2, C is exactly twice an
odd integer. Write c=C/2 and a0=a/2. The same ratio is

    a0*j / [c(A+a*j)].

The index-zero leaf is exceptional: it is (P,Q,T)=(1,1,1), and contains
no positive ratio factor. Its positive-leaf count is zero. A positive
leaf has (P,Q,T)=(a0*j,c(A+a*j),a0*j).

A split interval's raw triple represents its product P/Q and the sum
of its prefix products T/Q. The combine rule remains

    P = P_left P_right,
    Q = Q_left Q_right,
    T = T_left Q_right + P_left T_right.

## 2. A divisor depending only on the number of positive leaves

For s>=1, define

    F_s = odd(s!) / odd(lcm(1,...,s)),    F_0=1.

Here odd denotes the odd part. Each F_s is an integer. It divides all
three entries of every raw PQT triple with s consecutive positive
indices, independently of A, the starting index, and c.

To prove this, fix an odd prime q. A product of h consecutive indices
has q-valuation at least v_q(h!). The same is true of h consecutive
progression denominators A+a*j because a is invertible modulo q.
Factors a0 and c can only add valuation.

Each summand in T consists of h numerator factors followed by s-h
denominator factors, for some h. Its valuation is at least

    v_q(h!) + v_q((s-h)!)
      = v_q(s!) - v_q(binomial(s,h))
      >= v_q(s!) - floor(log_q s)
      = v_q(F_s).

The binomial valuation inequality follows directly from Legendre's
formula: each floor difference is zero or one, and vanishes above s.
The product P and denominator Q have at least v_q(s!) valuation, so
they satisfy the same bound. Every term in T is divisible by F_s, and
therefore so is T. A leading index-zero leaf introduces an additional
constant prefix term, represented by the full denominator Q, so the
same argument applies with the actual positive-leaf count.

No factorization of A+a*j is needed to find this divisor.

## 3. Normalize at every merge without constructing factorials

Suppose the children have already been divided by F_l and F_r. The
parent should be divided by the additional factor

    H(l,r) = F_(l+r)/(F_l F_r).

This is always an integer. For a prime power z=q^e<=s=l+r, its
contribution to the q-valuation of H is

    floor(s/z)-floor(l/z)-floor(r/z)
       + [z<=l] + [z<=r] - 1.

Each contribution is nonnegative. If z is at most both children, it is
the ordinary zero-or-one floor carry plus one. If z is at most exactly
one child, it is that carry. If z exceeds both children, then s<2z and
floor(s/z)=1, so the contribution is zero. This also handles an empty
child by treating F_0=1.

Compute the ordinary PQT merge of the normalized children and divide
P,Q,T exactly by H(l,r). The result is precisely the raw parent triple
divided by F_s. Induction proves that P/Q and T/Q are unchanged at
every node. The current small word batches contain at most four
positive leaves, and F_0 through F_4 are all one, so those batches need
no additional normalization.

The existing exact remainder comparison uses only P/Q, and the real
fixed-point quotient uses only T/Q. Their values and their acceptance
conditions are therefore unchanged. This is representation reduction,
not a new rounding rule or a floating-point approximation.

## 4. Construct and cache H, not F

F_s itself has Theta(s log s) bits and should not be cached. In contrast,

    H(l,r) = odd(binomial(s,l))
              * odd(lcm(1,...,l)) * odd(lcm(1,...,r))
              / odd(lcm(1,...,s))

has O(s) bits. This follows from binomial(s,l)<=2^s and the linear
bit-length bound for lcm(1,...,s).

A direct implementation scans odd primes q<=s, computes the nonnegative
exponent from Section 3, and multiplies the resulting prime powers in a
balanced product tree. This avoids large factorials and big-integer gcd.
Alternatively, balanced products of maximal prime powers supply the
three odd LCMs in the displayed formula. Serial multiplication into a
growing integer is not the justified construction.

A midpoint tree for one requested count R has at most two ordinary
interval lengths per depth, and at most one extra size pair from the
leftmost path containing the index-zero leaf. The sum of all those
distinct pair sizes l+r over the depths is O(R). Caching H for all
distinct pairs encountered in that call therefore requires O(R) integer
bits, not O(R log R).

Across the default weighted blocks, the earlier real-series analysis
gives sum R=O(N/u), where u=Theta(sqrt(log N)). Retaining all resulting
normalizers for one core or extraction pass consequently still uses
O(N/u) integer bits, even without cross-block reuse. The key/descriptor
count contributes O(B L^2) bits, where B=O(2^u), L=ceil(log2 N); this
is o(N). Including all failed floating proposals conservatively changes
that metadata bound to O(B L^3), also o(N), while geometric count doubling
keeps the total normalizer payload within the deterministic count bound.
Caches should have a core/pass lifetime, rather than persist unboundedly
across arbitrary requests.

For a deterministic lookup implementation, a balanced dictionary takes
O(L log(BL^2)) bit work per lookup, including comparison of the word-size
keys. There are O(sum R) actual tree merges, so lookup work is below the
existing core bound. Hashing should not silently be treated as a
worst-case constant-time operation in the proof.

The normalizer arithmetic is also subdominant. Legendre/floor
exponents plus balanced products give the conservative bound
O(R L^2) per call over its distinct size pairs, with the same bound
amortized across cached calls. Summing gives O(NL^2/u). A sieve to the
largest real count and its prime list fit O(N) data in the default
regime. For the manuscript's sequential model, conservatively allowing
an additional log L for prime enumeration and record routing gives
O(NL^2 log L/u), still below O(NL^(5/2)). This is separate from the
ordinary word/CRT stage and its earlier polylogarithmic-in-L costs;
improving the real representation does not remove those factors.

## 5. What the smaller operands do and do not prove

After normalization, a node with s positive leaves has

    bits(Q), bits(T) = O(s[1+log(N/s)]),
    bits(P) = O(s[1+log(R/s)]),

where R bounds its indices. The P bound uses the same factorial
cancellation on its numerator product. The T bound follows because
each positive ratio is less than 1/C, so the prefix sum is uniformly
bounded; the optional leading unit term changes only a constant.

At the root, Q and T shrink from O(RL) bits to
O(R[1+log(N/R)]) bits, and P shrinks to O(R). These are meaningful
reductions in the expensive upper-level products and final division.

The middle levels still matter. Writing d=log(N/R), the aggregate
denominator mass at depth j from the root is O(R[d+j]), while the
logarithm in multiplication cost is roughly log(R)-j. Summing the
corresponding costs can still give

    O(R (log R)^3 + d R (log R)^2) = O(RL^3).

Ignoring smaller logarithms, (d+j)(log R-j) is largest near
j=(log R-d)/2 when that lies in the tree. This corresponds to subtrees
with about sqrt(N) positive leaves, provided R>sqrt(N). A broad range
of middle levels therefore survives the root reduction in the
asymptotic regime. When the denominators in a particular block are
much smaller than N, replace N by their actual upper bound in this
description. Native multiplication thresholds can move the measured
cost distribution, so this analysis is not a prediction of timings.

Exact division by H fits the same node arithmetic bound. Thus this
normalization preserves the existing O(NL^(5/2)) core bound and linear
stored data. It does not prove O(NL^2) or O(NL^2 log L). Practical speed
must be measured against the current LCM/nested-inverse implementation.

The numerator products P are independent of A, which offers another
possible constant reduction if several blocks share a truncation count.
Caching an entire P tree, however, has larger mass than caching H and
needs a separate bounded-memory schedule. That extra optimization is
not assumed in the bounds above.

## 6. Read-only prototype correspondence

Inspected `tmp/reconstruction_next/real_factorial_reduced.cpp` with
`REAL_EXPERIMENT=7`. Its `RealFactorialReductionCache` uses a balanced
`std::map` of ordered child-size pairs and the exact exponent formula
of Section 3. Prime powers are multiplied through balanced in-place
levels. It stores H values only, never factorial F values. Its sieve
can be regenerated on a larger request; charging that request's size
to the new real call preserves the aggregate setup bound.

`factorial_reduced_real_bs` uses
`mid-max(lo,1)` and `hi-mid` as its positive-leaf counts. The recursive
branch has interval length greater than four, so these quantities are
nonnegative; its leftmost path accounts correctly for index zero.
At base intervals it calls the word-batched kernel with divisor
gcd(a,C). For odd p and power-of-two a, this gcd is exactly two;
the progression A+a*j itself is not altered. All three combined entries
are divided by H, and optional invariant checks verify exact divisibility.

The prime-power loop prevents overflow before increasing its power.
Its floor difference is zero or one, so the unsigned intermediate
subtractions cannot underflow; the following signed sum is nonnegative
by the lemma. The resulting prime powers have O(log s) bits individually
and are constructed as GMP integers. The core's RAII scope resets the
cache on entry and exit, including exception unwinding. Active variant
7 does not create the older, larger factorial-product cache retained in
the experimental file for other variants.

This read-only audit found no mismatch between variant 7 and the exact
normalization lemma. It did not run CPU tests while timing runs were in
progress. The prototype's test and benchmark results are separate
evidence and are not an exponent proof.

## 7. Reusing the normalized numerator tree

The separate `tmp/reconstruction_next/real_normalized_p_cache.cpp`
prototype uses a fixed midpoint tree with power-of-two span. Its cached
entry for an interval with s positive indices is

    product(positive indices in that interval) / F_s.

This is independent of A and C. The factor (a/2)^s is restored as a
binary shift when the entry is used. A real-series call constructs only
Q and T, using the cached left P in the cross term. It divides Q and T
by the same H as before. A partially occupied final interval uses its
actual positive counts; Section 3 applies to these unequal child sizes
as well. The root P is a prefix query in the cached tree, needed for the
unchanged exact remainder test.

The interval ending at index zero has no positive factors. The source's
`positive_count` and prefix construction correctly handle that exception.
Its batches of at most four indices fit its checked 128-bit native
product range; the constructor restricts the span to at most 2^32.
The cached tree is independent of the actual denominator values, while
the progression step is preserved in every Q leaf.

Let Rmax bound the largest real count in one core. A complete cached P
tree has O(Rmax log^2 Rmax) integer bits, by summing the interval bound
O(s[1+log(Rmax/s)]) over its levels. Its descriptor array has
O(Rmax log Rmax) bits. The default block schedule gives
Rmax=O(N/2^u), u=Theta(sqrt L), so both are o(N). A larger request
rebuilds a power-of-two span; geometric span growth bounds cumulative
setup by the largest tree, and at most the old and new trees coexist.
Core-scoped ownership prevents retention across an unbounded sequence
of extractions. The setup work is conservatively O(Rmax L^3), still
within the earlier total real-tree bound. These statements concern the
default schedule, not arbitrary external requests with Rmax=Theta(N).

The read-only source audit found its QT merge, prefix P, restored shifts,
and positive-count normalizers consistent with these identities. The
implementation agent separately reported exact tests of 612 triples,
513 prefix shapes and 462847 exact divisions, plus self-tests and core
checks. Those are finite checks, not a new complexity theorem.

The factorial-normalized variant alone was slower in the exploratory
large-position timing run. Caching was tested separately but also was
not selected for production. No speed claim for these variants follows
from their smaller operand bounds; their normalizer and cache overhead
must be counted. Exact timing comparisons belong to the final matched
benchmark report.

## 8. A shared coefficient transfer without large gcds

This section establishes a new local construction for the prospective
Lucas route in `weighted_machin_lcm_derivation_2026-10-24.md`. It is not
implemented in the production extractor and is not a full global
O(NL^2 log L) proof.

Set D=p^2 and C=D+1. Define

    E_r(x) = [C^r-(1-x)^r]/[D+x] = sum_(k=0)^(r-1) e_k x^k,
    b_k = (-1)^(k+1) binomial(r,k), including b_0=-1.

The numerator vanishes at x=-D, so E_r has integer coefficients. They
satisfy

    D e_0=C^r-1,
    D e_k+e_(k-1)=b_k,
    b_k=-(r-k+1)b_(k-1)/k.

Consequently the transfer from (e_(k-1),b_(k-1)) to (e_k,b_k) has the
integer matrix numerator

    X_k = [ -k,  -(r-k+1)  ]
          [  0, -D(r-k+1)  ]

and common denominator D*k. For an interval (l,t], let s=t-l, with
0<=l<t<=r. Its raw numerator is X_t ... X_(l+1), and its denominator is

    D^s product_(k=l+1)^t k.

Define the full, rather than odd, normalizer

    G_s=s!/lcm(1,...,s),    G_0=1.

It divides both diagonal entries, the off-diagonal entry, and the
common denominator of every such raw transfer. For the diagonals and
denominator, this follows because each contains a product of s
consecutive positive integers, whose q-valuation is at least v_q(s!)
for every prime q. Powers of D can only increase valuation.

For the off-diagonal entry, expand the product according to its unique
upper-right transition at h, l<h<=t. Its integer factors, apart from
powers of D and a sign, are a product of t-h consecutive indices from
the upper diagonal and a product of h-l consecutive descending indices
from the lower diagonal and the transition itself. For every prime q,
the valuation of this term is therefore at least

    v_q((t-h)!)+v_q((h-l)!)
      =v_q(s!)-v_q(binomial(s,t-h))
      >=v_q(s!)-floor(log_q s)=v_q(G_s).

Every summand is divisible by G_s. This proof includes q=2; the odd-part
restriction needed by the progression denominator is absent here.

At a merge, divide all three matrix entries and the common denominator
by G_(l+r)/(G_l G_r), where l,r here denote child lengths. This factor
is an integer by exactly the prime-power carry argument of Section 3,
now including prime 2. Its bit length is O(l+r). It can be constructed
and cached from its child lengths without factoring matrix entries or
calling a large gcd. Induction gives the raw transfer divided by G_s
at every node. Exact transfer values are unchanged.

### Retained denominator valuation

The resulting denominator is explicitly

    Q_(l,t) = D^s product_(k=l+1)^t k / G_s.

For a query prime q different from p, D is a unit, so

    v_q(Q_(l,t))
      = sum_(e>=1) (floor(t/q^e)-floor(l/q^e)-floor(s/q^e))
          +floor(log_q s)
      <=floor(log_q r)+floor(log_q s)
      <=2 floor(log_q r).

Each floor difference is zero or one, and vanishes above t<=r. This
bound applies to the deliberately incompletely reduced denominator;
no hidden gcd is assumed. The earlier fully reduced transfer had the
smaller one-log bound, so this simplification uses a larger constant
in the sparse precision padding.

For example, a query modulo q^E can carry source data modulo
q^(E+lambda_q), with lambda_q=2 floor(log_q r). If W is the matrix
numerator applied to that source residue, g=gcd(Q,q^(E+lambda_q))
divides W, the remaining modulus is still a multiple of q^E, and

    (W/g) * (Q/g)^(-1) mod q^E

is the exact propagated residue. This describes one propagation; the
full batched query schedule must retain the appropriate budgets through
its tree. Gcds or prime-power valuations at that sparse division stage
have not been eliminated. The reduction removes large integer gcds
from the shared exact matrix-value tree.

### Size and cost of the shared exact tree

For fixed D and an interval of length s, all normalized entries and its
denominator have O(s[1+log(r/s)]) bits: bound each raw factor by a
constant times r, divide by s!, and restore the O(s)-bit LCM. The at
most s off-diagonal summands add only O(log s) bits. A balanced tree
therefore has O(r log^2 r) retained bits and can be built in
O(r log^3 r) arithmetic work under M(x)=O(x log x). The normalizer setup
fits this bound. These are conservative full-tree bounds, not merely
root-size estimates.

In the prospective parameter regime r=Theta(N/L^2), this shared tree
fits O(N) stored data and its construction costs O(NL). Query routing,
prime-power coefficient work, block reconstruction, recursive scaling,
and the output certificate still need their separate global accounting.

## 9. Explicit formulas for the remaining Lucas implementation

The new transfer lemma provides an exact shared-value tree. The
following identities specify how its coefficient queries connect to
the Machin real correction; they do not substitute for the missing
global implementation or cost proof.

Write

    P_k = sum_(j=0)^k binomial(r,j) D^j.

Unrolling the coefficient recurrence gives, modulo any prime q!=p,

    e_k = (-1)^k D^(-k-1) (C^r-P_k).

The sign and positive D^j inside P_k are significant. For an odd prime
q, write r=r_1*q+r_0 and k=h*q+k_0, with 0<=r_0,k_0<q. Lucas' theorem
and D^q=D modulo q give

    P_k = C^r_0 * sum_(j=0)^(h-1) binomial(r_1,j) D^j
           + binomial(r_1,h) D^h
             * sum_(j=0)^min(k_0,r_0) binomial(r_0,j) D^j  (mod q).

When q>sqrt(r), the high row has fewer than sqrt(r) entries. A query
for e_(k_0) in the shared full-r tree also gives the low prefix, since
the h=0 instance reads

    P_low=C^r-(-1)^k_0 D^(k_0+1)e_(k_0)  (mod q).

The outside multiplier (-1)^k D^(-k-1) advances by -D^(-1) between
successive incidence indices spaced q apart. Prime powers need the
sparse prime-power machinery; the displayed Lucas formulas alone are
only prime-modulus statements.

The fixed prime p is exceptional because D is not invertible there.
It must be explicitly assigned to the small-prime branch, even if a
finite-size threshold would otherwise classify p=239 as large. Use

    e_(r-1)=(-1)^(r+1),
    e_(k-1)=b_k-D*e_k

backwards, with a unit-part binomial row to obtain b_k modulo the needed
power of p. This avoids division by p. Prime 2 is not a denominator
query for the odd progression A+a*k.

### A smaller auxiliary scaling

For odd A and power-of-two a>=2, let

    H_r(A)=sum_(j=0)^(r-1) a^j*j! /
             [C^j product_(i=1)^j(A+a*i)].

The exact partial-fraction identity is

    C^(r-1) H_r(A)
      = C^(r-1)-a sum_(k=1)^(r-1) k*e_k/(A+a*k).

After leafwise cancellation of 2, H_r has an odd denominator. Its
denominator divides c^(r-1) times the progression LCM, where c=C/2.
Therefore c^(r-1)H_r has an odd denominator dividing that progression
LCM. It can be reconstructed using the same odd prime queries and
coefficients

    -a*k*e_k*2^(-(r-1)) mod q,

with constant term c^(r-1). There is no omitted 2-adic denominator:
the preceding denominator-divisibility argument proves its absence.
If the larger C^(r-1) scaling is used first, its numerator against an
odd common denominator is divisible by 2^(r-1), so the change can also
be made by an exact binary shift. This removes a constant factor in
the auxiliary bit length; it is not an exponent improvement.

The 239-core's small-A range can potentially be isolated as described
in the earlier note: O(N/L) such terms cost O(NL^2) by direct exact
splitting, leaving A=Omega(N/L). Parameters
ell=Theta(N log L/L^2), r=Theta(N/L^2) then give A/r=Omega(L).
An exact splice, complete tail inequalities, and recursive costs remain
required. No full-prefix arctangent/AGM replacement is adopted here;
it would change the intended no-prefix computational scope.

## 10. Exact half-inverse reconstruction

The separate `tmp/reconstruction_next/direct_solve.cpp` experiment
solves D*W=T modulo p^m using an inverse only to half the target
precision. Assume p does not divide D and m>=2. Put

    P=p^ceil(m/2), Q=p^floor(m/2), b=D^(-1) mod P,
    x=T*b mod P,          0<=x<P,
    e=(T-D*x)/P,
    y=e*b mod Q,          0<=y<Q,
    W=x+P*y.

The division defining e is exact; e can be negative. Since Q divides
P, b is also an inverse modulo Q. Hence

    D*W = T-P*e+P*D*y = T (mod P*Q).

The stated ranges give 0<=W<P*Q, so W is already canonical. For m=1,
the existing ordinary inverse is used. Signed exact GMP division and
canonical modulo implement this identity without a sign assumption.

At the audited call site, T includes the multiplier Y modulo the
reduced modulus p^(m-shift). After solving, multiplication by p^shift
restores the outer scale, and the final family sign is applied once.
Thus folding Y into T changes neither the core residue nor its output
certificate. Empty accumulators remain zero.

In the exact LCM core, T and D have O(N) bits. The nested half-precision
inverse and a constant number of products, reductions, shifts and
exact divisions cost O(M(N)) and use O(N) stored data under the same
regularity assumptions as the full nested inverse. This changes a
finite constant, not the asymptotic exponent or adaptive-guard scope.
The implementation agent reported 2448 independent exact cases; this
note's source audit is additional algebraic evidence, not a benchmark.

## 11. Falsifiable checks and remaining barrier

`weighted_machin_transfer_normalization_checks_2026-10-24.py` uses only
Python's standard library. Its exact arithmetic checks passed on
24 October 2026 with the following counts:

| Check | Cases |
|---|---:|
| Normalized transfer intervals and coefficient action | 2648 |
| Retained denominator valuation bounds | 17055 |
| Padded sparse exact divisions | 17055 |
| Merge normalizer integrality and prime exponents | 1156 |
| Real partial fractions and smaller odd scaling | 432 |
| Lucas prefixes and coefficient recovery | 400 |
| Odd-factorial real PQT intervals | 624 |
| Exact real PQT merge normalizers | 6701 |

These checks include exhaustive small intervals and deterministic
random larger ones, p=3,5,13,239 where applicable, prime 2 in matrix
normalizers, and query primes dividing matrix denominators. They do
not test an entire Lucas-based extractor, nor do finite tests prove
the global running-time estimate. The real PQT cases use odd A,
power-of-two a, and C=p^2+1; they include intervals containing index
zero and intervals of only positive indices. They compare each
normalized triple against the raw triple divided by F_s, check every
merge division, and independently compare the represented product
and prefix sum with exact rational arithmetic.

The established new ingredients are exact cancellation, its source
correspondence, bounded shared-transfer denominators without large
gcds, and the half-inverse identity. The current weighted real tree
still has its conservative O(NL^3/u) bound, and local CRT reconstruction
still costs O(NuL^2). Removing the balancing parameter u requires the
shared coefficient queries and recursively scaled reconstruction to
replace those stages with a fully accounted schedule. Neither a
smaller final denominator nor a small root matrix establishes that
replacement. O(NL^2 log L) remains a prospective target, strict O(NL^2)
is not established, and practical adoption requires matched timings.
