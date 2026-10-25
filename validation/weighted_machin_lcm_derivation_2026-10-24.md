# Exact LCM assembly for weighted Machin: derivation and source audit

24 October 2026. Mathematical design notes; no timing claims. This extends
the existing one-pass audit, with fixed prime p, maximum denominator O(N),
guard G=O(log N), L=ceil(log2 N), u=Theta(sqrt(L)), and block length
ell=Theta(N/2^u). Multiplication is charged by M(b)=O(b log(b+2)).
The implemented construction is `decimal_extract_weighted_machin.cpp`.
This note supersedes the forest-assembly accounting in
`weighted_machin_derivation_2026-10-24.md`; that earlier note's exact
weighted-block and real-correction lemmas remain the starting point.

## 1. Exact LCM assembly supports a single streaming balanced tree

For a consecutive interval I of s original core terms, represent its
contribution as p^shift*n/d, where d is a positive p-unit and shift is the
smallest outside exponent of its constituent blocks. Do not reduce d
modulo a p-power while using its prime factor metadata.

Every block denominator divides the p-free LCM of its original odd
denominators. Hence an exact LCM merge preserves

    d | lcm(mu_j / p^v_p(mu_j) : j in I).

The short-progression LCM bound gives

    bits(d) = O(s [1 + log(N/s)]).

Here is a direct proof for the odd interval, rather than an appeal to a
lower bound on its starting point. Put mu_j=A+2j (0<=j<s), let P be their
product, and let D be their LCM. For an odd prime l let c_e count the
mu_j divisible by l^e, and let h be their largest l-valuation. Since 2
is a unit modulo l^e, c_e>=floor(s/l^e); in particular every l^e<=s
occurs. Thus

    v_l(P)-v_l(D) = sum_(e=1)^h (c_e-1)
       >= v_l(s!) - floor(log_l s).

Writing odd(x) for the odd part of an integer, this proves

    D | P * odd(lcm(1,...,s)) / odd(s!).

The expression on the right is an integer: v_l(P)>=v_l(s!) for every
odd l. Its binary logarithm is at most

    s log2(max mu) - log2(s!) + v_2(s!)
       + log2(lcm(1,...,s))
      = O(s[1+log(N/s)]),

using Stirling's elementary factorial estimate, v_2(s!)<=s, and the
standard log(lcm(1,...,s))=O(s) bound already used in the manuscript.
Removing the p-part, or any zero-coefficient prime fractions, only
decreases the denominator.

The exact weighted sum, after removing its smallest outside power, has
absolute value at most exp(O(s+L)). Multiplying by d therefore gives

    bits(n) = O(s [1 + log(N/s)] + L).

This estimate does not require cancelling gcd(n,d). In the descending-r
traversal, outside exponents decrease; the difference between extremes
is 2s plus a difference of local p-valuations bounded by O(L).

For adjacent children, let g=gcd(d_left,d_right), computed from their
sorted prime-power metadata. Set a=d_left/g, b=d_right/g, and
shift=min(shift_left,shift_right). The exact merge is

    n = p^(shift_left-shift) n_left b
        + p^(shift_right-shift) n_right a,
    d = a d_right.

The integer corrections already present in each n are retained exactly.
This is ordinary rational addition with a common LCM denominator; it
does not lose the integer carries that a global prime-fraction registry
would discard.

### Streaming storage

Insert consecutive blocks into a binary-counter stack. Keep at most one
node of each dyadic block count; merge adjacent equal-count nodes on a
carry. At any time the interval sizes are geometrically separated.
Writing b=Theta(2^u), their aggregate numerator/denominator size is at most

    sum_(k=0)^log2(b) O(ell 2^k [1+log(N/(ell 2^k))] + L)
      = O(N).

The merge temporaries change only the constant. A final partial block
changes each interval length by at most ell, which does not change this
bound. Fold the remaining stack from small to large, with the older node
as the left operand and the already combined newer suffix as the right
operand. This leaves large operations until the end. A high-to-low fold
is also exact, but may repeat full-width operations O(u) times.

The final result needs one p-power reduction and one p-adic inverse.
There is no mathematical need for intermediate capped forests or
intermediate p-power moduli in this representation.

### Metadata also has linear storage

Keep exactly the nonzero CRT coefficients' (prime, local maximum power)
entries, because the block's current denominator is the product of those
powers. The inspected implementation concatenates increasing stripped
primes at most sqrt(maximum) with sorted residual primes greater than
sqrt(maximum), so its combined metadata is already globally sorted.

For an original block of ell terms, every maximal prime power q in its
factor list exceeds sqrt(ell): if its prime exceeds sqrt(ell), this is
immediate; otherwise the progression contains a multiple of the largest
power of that prime at most ell, which is greater than ell/prime.
Consequently every retained q has Theta(L) bits when u=o(L). Each
O(L)-bit metadata record is chargeable to its denominator factor. A union
at a later node only increases that prime's power. Thus factor metadata
obeys the same O(N) stack bound. For the possible last partial block of
length t, the same charge applies if t>=sqrt(ell), since log(t)=Theta(L).
If t<sqrt(ell), its crude O(t L^2) metadata bound is o(N). These exceptional
records can follow only one ancestor at a time; their work over O(u)
levels remains sublinear times a fixed power of L. This handles partial
blocks without assuming they have constant size. The bound uses the
stated default asymptotic parameter regime.

### Work, including the product of common factors

Intersection of sorted factor lists does not materialize its integer
product for free. There is a useful additional size bound. If q^e divides
both child denominators, it divides one original denominator in each
child. These distinct denominators differ by a nonzero integer of
magnitude less than 2s. Hence

    g | lcm(1,...,2s),     bits(g)=O(s).

A balanced product of the common prime powers therefore costs
O(M(s) log(s+2))=O(s L^2). Exact divisions and rational additions at the
node cost O(M(s[1+log(N/s)]+L)). At a dyadic level all s sum to O(N).
Summing O(u) levels bounds assembly by

    O(N u L^2 + N L u^2).

Factor-list comparisons and copies add O(N u^2) bit traffic under the
metadata charge above. The last stack fold fits these bounds. In
particular, the new scheme preserves the current O(N L^(5/2)) total
bound and O(N) stored data, while replacing repeated p-power reductions
with exact LCM intersections. It is an implementable finite-constant
improvement, not yet a smaller final exponent: the local CRT trees still
cost O(N u L^2).

Computing a generic big-integer gcd at every merge is a different
implementation. Its conservative cost includes another logarithm of
the current operand size; the factor-list scheme avoids charging that
gcd to the larger short-interval denominator.

## 2. Read-only source audit

Inspected the selected production `decimal_extract_weighted_machin.cpp`
on this date. A read-only diff against the audited LCM prototype showed
only its descriptive header, nested Newton inversion, and final
`scale_residue` helper plus the helper's two call sites. The LCM block
construction, factor metadata, streaming stack, and merge formulas are
unchanged. The integrated inverse matches the separately audited helper.
This is an internal proof/code correspondence check, not a formal
verification. Arithmetic tests and measurements are reported separately;
this audit did not run CPU tests during the final timing window.

The code retains a factor entry if and only if its local CRT coefficient
is nonzero. Its unchanged rational tree uses exactly those factors once
each, so their product equals the stored block denominator. Each LCM
merge uses minima of common powers for g and maxima for the merged
support. Its exact divisions and shifted numerator formula implement
Section 1. No numerator or denominator is reduced before final core
evaluation. The final stack fold proceeds from low to high level, with
the older interval on the left and newer suffix on the right.

### Vector capacity and merge temporaries

`support.reserve(left.size()+right.size())` is safe for the linear-space
argument: the union has at least the larger child count, so this reserved
capacity is at most twice its actual record count. A temporary common
factor vector reserves at most the smaller child count; even if the
intersection is empty, this is bounded by the children already being
processed. Factor-product levels are overwritten in place and resized;
only the old and new level coexist, not a retained full product tree.
Their object-array capacity is chargeable to the common factor list.

The binary-counter stack holds O(u) node descriptors. In the installed
GMP C++ header (`C:/msys64/ucrt64/include/gmpxx.h`, inspected with C++17
enabled), `mpz_class`'s move constructor transfers the raw mpz structure
and calls `mpz_init` on the source. A moved-from stack node therefore does
not retain its old limb allocation. Move assignment swaps structures;
at the relevant assignments the source's destructor frees the replaced
allocation. The default-allocator vectors likewise transfer ownership.
Passing stack nodes by value with `std::move` into `merge_pair` empties
the inactive stack slots before their old payloads are destroyed. Thus
inactive slots do not accumulate historical O(N)-bit node allocations.

At a merge there are only a constant number of child, parent, quotient,
shift, and numerator temporaries. Their lengths have the node bound of
Section 1. Add these to the other occupied stack levels to retain O(N)
stored data. The `maximum_pending_bits` diagnostic measures retained
stack values and factor-vector capacity; it is not a measurement of
peak RSS, merge temporaries, or GMP's arithmetic workspace.

### Integrated nested Newton inversion

Let the target exponent be m. Starting from m and repeatedly replacing
q by ceil(q/2), then reversing the list, gives successive precisions
q_new <= 2*q_old. With b an inverse modulo p^q_old, the update

    b <- b (2-a*b) mod p^q_new

is exact because 1-a*b_new=(1-a*b)^2 before reduction. Replacing a by its
residue at the new modulus preserves this identity. The initial inverse
modulo p exists because the LCM denominator omits p. This schedule avoids
an unnecessarily large penultimate precision for targets just above a
power of two; it changes constants, not the proof's exponent.

The production implementation reuses both moduli and operand residues.
Successive exponents are 2q or 2q-1, so squaring the preceding modulus
and, in the latter case, dividing exactly by p gives the next modulus.
It reduces a once at the largest modulus, descends through nested operand
residues, and then lifts upwards. Modulus divisibility
p^q_small | p^q_large makes that preprocessing exact.

The single initial reduction costs O(M(max(bits(a),m))); every later
reduction has bounded operand-size ratio. Squarings, exact division by
the fixed word p, residue reductions, and Newton updates have geometric
operand sizes. Consequently the entire integrated inverse costs
O(M(N)) with O(N) stored data in the audited core regime. The precision
schedule contains only O(L) word entries. This includes the reductions
of a: independently rereading the full denominator at every precision,
as in the earlier prototype, would only give the more conservative
O(NL^2) inverse bound.

### Removing powers of two from the final multiplier

`scale_residue` first takes y=Y mod modulus in the canonical range.
If y=0 it returns zero. Otherwise it writes y=2^t*z with z odd, multiplies
the current residue by z, shifts by t, and reduces. The result is exactly
the original product modulo the modulus, including signed Y. This
avoids multiplying binary zero limbs when the multiplier has a large
power of two. Its scan, shift, product, and final reduction involve
O(N)-bit operands and fit O(M(N)) work and O(N) stored data. It changes
neither the core identity nor the output certificate.

## 3. Scope

This assembly change preserves exact core residues, O(N) stored data,
and the existing one-pass O(N log^(5/2) N) bound for default G=O(log N).
It does not establish a full-adaptive extractor bound or O(N log^2 N).
Local prime-fraction reconstruction and the real correction still have
their earlier bounds O(NuL^2) and O(NL^3/u).

Choosing one-term blocks and assembling LCMs alone does not prove an
O(NL^2) algorithm: at a dyadic level whose intervals have length N/2^d,
the aggregate denominator mass may be O(Nd). Summing the corresponding
multiplication work over Theta(L) levels can still give O(NL^3).
Globally grouping coprime prime-power fractions would reduce that mass,
but their missing integer carry must then be recovered exactly.

The following prospective extension is separate from the implemented
LCM construction and its measurements.

## 4. Future work: transferring the Lucas route toward log squared

**Status: not implemented here, not a complete proof, and no performance
claim.** The existing `notes/lucas_smaller_parameters_2026-09-07.md` and
`notes/lucas_sparse_prime_power_batch_2026-09-07.md` describe a candidate
O(NL^2 log L) construction. Its larger machinery has not established a
practical advantage for the Machin extractor. The identities below show
how to begin transferring that route; they do not assign its bound to
the current production program.

The Machin positive Euler series has a useful partial-fraction
polynomial. For C=p^2+1, define

    H_r(A) = sum_(j=0)^(r-1) a^j j! /
             [C^j product_(i=1)^j (A+a*i)],
    E_r(x) = [C^r-(1-x)^r]/[C-1+x] = sum e_k x^k.

The numerator vanishes at x=1-C, so E_r is an integer polynomial.
Integrating the finite geometric sum in (1-x^a)/C gives

    C^(r-1) H_r(A)
      = C^(r-1) - a sum_(k=1)^(r-1) k*e_k/(A+a*k).

This replaces a large exact product denominator in the real tree by a
short-progression rational reconstruction problem with shared
coefficients across blocks. Writing D=C-1=p^2, those coefficients obey

    D e_0 = C^r-1,
    D e_k + e_(k-1) = (-1)^(k+1) binomial(r,k).

Let b_k=(-1)^(k+1)binomial(r,k). The transfer from k-1 to k for the
vector (e_k,b_k) has integer matrix numerator

    [ -k,  -(r-k+1)   ]
    [  0,  -D(r-k+1)  ]

and common denominator D*k. Unrolling a segment [l,t] shows that its
reduced denominator divides D^(t-l)*binomial(r,l), by the same exact
binomial-ratio argument used in the earlier Lucas note. At a queried
prime q other than p, its lost precision is therefore at most
floor(log_q r). This is the bounded-loss input needed by the existing
sparse prime-power division construction. The sign and fixed constant
change, but the coefficient recurrence's denominator structure survives.

The prime p is an explicit exception: D is not a unit there. It cannot
be silently omitted from the auxiliary real rational reconstruction.
Since p is fixed, generate its row backwards using

    e_(r-1)=(-1)^(r+1),
    e_(k-1)=(-1)^(k+1) binomial(r,k) - D e_k.

This uses no division by p. A unit-part binomial row supplies residues
across multiples of p. Shared fixed-prime row work fits the earlier
small-prime preprocessing budget; it need not be regenerated per block.

The 239-core also reaches A=1, unlike the lower-bound hypothesis in the
old Lucas notes. One way to isolate that issue is to handle the first
O(N/L) small-denominator terms by ordinary exact binary splitting. Its
O((N/L)L^3)=O(NL^2) work and O(N) raw product bits fit the prospective
target. For the remaining terms, A>=cN/L. With

    ell=Theta(N log L/L^2),   r=Theta(N/L^2),

one has A/r=Omega(L), so the positive transformed tail gains
Theta(log L) bits per term. A sufficiently large fixed constant in r
therefore allows the same orders of block and recursive parameters,
after their stopping inequalities are rechecked. The small-A piece
would require an explicit splice into the resulting exact core.

These identities give a concrete route to adapting the earlier
candidate O(NL^2 log L) argument. The coefficient batches, sparse
division invariants, small-A splice, sequential routing, space bounds,
and output scope still require an implementation or a complete mapping
to independently checked lemmas. There is no O(NL^2) theorem here.
Nor does an asymptotic candidate establish that its extra machinery
will outperform the measured weighted Machin implementation.


