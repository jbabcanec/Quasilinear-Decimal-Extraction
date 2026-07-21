# Cracks Dossier — the N log³N → N log²N campaign
**Compiled 2026-07-21 from three deep literature crawls + in-session validation.**
Companion to `decimal_extraction_paper.tex` §6 and `research_log.md` Part 16.

Status line: **one crack validated and banked tonight** (multiplicative faces → O(N log²N));
the barrier now stands on the additive core alone. Below: the crack, every open lead
ranked, every closed door with its reason, and the adjacent-field cost table.

---

## 0. THE VALIDATED CRACK (now Theorem in the paper)

**Borwein-factorization for AP-products.** P. Borwein 1985 computes n! below the
product-tree bound by factoring first (Legendre counts) and assembling Π p^{e_p} by
exponent-bit scheduling. The prerequisite — cheap factorization of the product —
holds for Π(8k+1): e_p = Σ_j #{k ≤ K : p^j | 8k+1} is an O(1) count per prime power.
Deep reason it wins: the term list carries Θ(N log N) bits; the distinct-prime
representation carries Θ(N) — multiplicities ride free as exponents.

- **Result:** Π_{k≤K}(8k+1) mod 2^m in O(M(N) log N) = O(N log²N) — hence E = dlog(Π),
  the power sums Σ(1+8k)^{2^r} mod 2^{2r}, and Diamond log-gamma values at valuation-(−3)
  arguments, all at O(N log²N). First recorded sub-product-tree evaluation of any of these;
  no p-adic Gamma routine exists in FLINT/Arb at all.
- **Validated:** `scratchpad/borwein_crack.py` — exact agreement with direct product at
  6 sizes through K = 3·10⁶ (7.2M bits); time ratio monotone toward the factored method
  (0.30→0.75), the log-factor signature; crossover extrapolates to K ~ 10⁷–10⁸ unoptimized.
- **What it does NOT give:** V (the digit's additive core). Sums have no prime factorization.
  The seven faces split: multiplicative cluster {products, E, power sums, log-gamma,
  truncated-KL} falls; additive core {V, merge depth} stands at N log³.
- **Methodological note:** the crack entered through precisely the one representation
  absent from the Mass Conservation Conjecture's basis list — as the conjecture predicted.

## 1. OPEN LEADS, RANKED

### L1. Hiary–Postnikov block decomposition of the additive core (arXiv:1205.4687)
Hiary computes incomplete character sums mod p^k — *incomplete sums of exponentials of
p-adic dlogs* — by Postnikov's formula: Taylor-expand the dlog on blocks, reduce to
quadratic exponential sums, evaluate those in polylog by theta-descent. **Only published
work whose object matches ours.** Transfer facts established: 2-adically the block-Taylor
step is EXACT (valuations make truncation an identity), so E-type sums decompose exactly
into Faulhaber polynomial sums per block. Honest accounting still lands ≈ Õ(N²) because
the cost migrates to Bernoulli-polynomial data — BUT this decomposition has never been
combined with (a) the Borwein-factorization crack, (b) tapered precision, (c) the
sieve-batched chains. Composition of cracks unexplored. **Priority: attack V with blocks
whose interior sums route through the now-cheap multiplicative objects.**

### L2. Incomplete Landsberg–Schaar reciprocity with algebraic remainder
The existence proof that "incomplete + arbitrary cutoff + exact + log time" is achievable:
**Dedekind sums** (Σ⌊(ak+c)/b⌋) fall to O(log b) via reciprocity + Euclidean descent
because remainders stay polynomial in the ring (Tranbarger–Yang IJNT 2024). One degree up,
Fedotov–Klopp (arXiv:0909.3079) give an exact identity for quadratic exponential sums —
Gauss-map descent, O(log N) levels — but the remainder is a Fresnel-type contour integral,
transcendental. Complete-sum Landsberg–Schaar IS exact. **The wanted object: an incomplete
L–S reciprocity for 2-adic quadratic sums whose remainder stays in ℤ/2^m.** Nobody has
ruled it out. This is the "crack looks like a reciprocity law" thesis.

### L3. Holonomic completion detour (Coleman + CMTV)
CMTV (arXiv:2106.09315) proves quasi-linear-precision 2-adic bit-burst for the holonomic
class. Coleman (Invent. Math. 1982) writes COMPLETE L_p values as finite combinations of
p-adic polylogarithm values (holonomic-evaluable!). Our E/V are incomplete segments.
**Wanted: a completion identity trading the truncation for finitely many holonomic
evaluations** — then CMTV bit-burst gives Õ(N). No such identity known; the complement
tail looks as hard as the segment; ranked third as the only visible route to true
quasi-linearity for the additive core.

### L4. Taper × AP-points × transposed evaluation (unexplored interaction)
Crawler-1's observation: the coefficient taper (c_t needed mod 2^{m−3t}), evaluation
points on an AP, and Tellegen-transposed multipoint evaluation have never been studied
TOGETHER. Naive accounting: constant factors only; but Harvey's Prop. 2 (arXiv:1209.0533)
shows the identical taper arises inside Voronoi congruences — a hybrid (odd auxiliary
primes for bulk + 2-adic reconstruction of a single residue, avoiding the full array)
is unexamined in the literature. Related open object: "incomplete Wilson congruence"
(arXiv:1310.2691) — BGS-based, open conjectures, nothing below √-time.

### L5. Sub-√K for factorial-type targets at KNOWN-factorization modulus
The Shamir/Lipton shield (fast factorials factor the modulus) evaporates at modulus 2^m.
No complexity-theoretic obstruction exists to beating Õ(√K) ring ops for K! mod 2^m —
it is merely unachieved. Any progress transfers instantly.

### L6. q-microscope as a depth-organizing language (weak)
Guo–Zudilin creative microscoping converts congruence depth into vanishing order at
cyclotomic points of q-identities. Truncation points must align with cyclotomic
structure — our arbitrary K does not — but it is the only generating-function handle
on DEPTH in the literature. Watch, don't attack.

## 2. CLOSED DOORS (reason + source)

| Door | Why closed | Source |
|---|---|---|
| D-algebraic bit-burst (general) | PTIME reals are poly-ODE values; time hierarchy forbids universal quasi-linear evaluation | Bournez–Graça–Pouly JACM 2017 |
| Numeric composition route | vdH "Fast composition of numeric power series": Õ(n²) at n-bit precision | vdH 2008 |
| KL/Kinoshita–Li composition | composition cannot re-index exponents; ring-op-counted, blind to 2-adic smallness | arXiv:2404.05177 |
| Power projection / transposed composition | same bound, same blindness | KL24 Thm 6 |
| Modular composition (2026 frontier) | n^{1.343} arithmetic, fields, no exponent re-indexing | arXiv:2601.17422 |
| Chirp/Bluestein re-indexing over ℤ/2^m | unit group has 2-torsion only — no roots of unity to chirp with | crawler-2 sweep |
| Materializing tapered Bernoulli array | output size Ω(N²) unconditionally | Σ(m−6i) count |
| Multisection (Buhler–Harvey) | decimates index at uniform precision; never cuts precision | arXiv:0912.2121 |
| Harvey n^{4/3} single-Bernoulli | gain is amortization across many primes; at fixed p=2 degenerates to quasi-quadratic (his own remark) | arXiv:1209.0533 |
| Pollack–Stevens / Lauder–Vonk | O(m²)-class in precision; quadratic floor in the moment data structure | LV final, PS 2011 |
| Coleman integration transfer | needs finite-dim cohomology + holonomic connection; log-gamma has neither | Best arXiv:1806.03393 |
| Strassen/BGS at precision 2^m | Õ(√K) RING ops → Θ(N)-bit elements → ≥ N^{1.5}; worse than our tree | BGS 2007 |
| Fractional-range congruences at arbitrary K | torsion-locked: cutoffs must sit at rational points of the period | Lehmer 1938, Pan arXiv:0905.0941 |
| Supercongruence depth | record is modulus p^{3r} on range p^r (Guo BLMS 2025) — exponentially shallower than modulus 2^{Θ(K)} | Guo 2025 |
| Fermat-quotient amortization along primes | per-p everywhere; amortization exists only for prefix-structured objects (Wilson) | Dorais–Klyve, CGH |
| Hiary theta-descent DIRECT transfer | polylog(1/ε) approximation ≠ exact mod 2^m; complex values don't determine residues; incomplete L–S remainder transcendental | Hiary 0711.5002, F–K 0909.3079 |
| p-adic fast quadrature | Volkenborn canonical by translation-invariance; EM = its acceleration; Bernoullis = its monomial integrals | paper Remark (root) |
| Formula shopping | geometric denominators ⇒ rational sum; transcendence forces AP moduli | paper Remark (root) |
| Quantum | per-element dlog classically easy; exact summation defeats amplitude estimation | paper Remark (root) |
| Sieve/von Mangoldt restructurings | exact rearrangements conserve the hard core | session validation |
| Iwasawa doubling / ψ-recursion | operator = Taylor shift; state vector mass Θ(N²/s); tapers don't bite (t = log N only) | session analysis |
| Newton/Abel transform of V | closed Δ-coefficients but transformed sum is hypergeometric → splitting again | session validation |
| lcm-thinning | mid-tree mass survives; constants only (~5×) | session analysis |
| CMTV beyond holonomic | their class is exactly p-adic linear ODEs; gamma absent by design | arXiv:2106.09315 |

## 3. WHAT ADJACENT FIELDS PAY (context table)

| Object | Best known (precision m) | Source |
|---|---|---|
| single B_n (exact) | n^{4/3+o(1)} | Harvey 1209.0533 |
| all B_k mod one prime p | Õ(p) | Buhler–Harvey |
| all B_k mod 2^m tapered | Ω(N²) (output size) | — |
| Kubota–Leopoldt / p-adic L | O(m² log m) | Lauder–Vonk |
| Morita Γ_p at a point | no algorithm in standard libraries | FLINT docs |
| p-adic log/exp/polylog/hypergeom | Õ(m) | CMTV 2021, FLINT |
| (p−1)! mod p³ (Wilson) | p^{1/2+ε} per prime; Õ(X) for all p ≤ X | CGH 2014 |
| **Π(8k+1) mod 2^m** | **O(N log²N) — this project, tonight** | Thm (crack) |
| **E, power sums, G₂ values** | **O(N log²N) — this project, tonight** | Thm (crack) |
| **V (additive core) → the digit** | O(M(N log N) log N) = N log³N — this project | Thm (batch) |

## 4. SOURCES (consolidated from three crawls)
Hiary 0711.5002 / 0711.5005 / 1205.4687; Fedotov–Klopp 0909.3079; Tranbarger–Yang IJNT 2024;
Pan 0905.0941; Z.-H. Sun (JNT 128, CMJ); Tauraso 1701.00729 / 0905.3327; Guo BLMS 2025;
Kalinin–Zottor 2602.00206; Goodman 2305.05522; Dorais–Klyve JIS 14; CGH 1209.3436;
Hart–Harvey–Ong 1605.02398; Costa–Kedlaya–Roe ANTS 2020; Harvey 0807.1347 / 1209.0533;
Buhler–Harvey 0912.2121; Lauder–Vonk; Roblot 1110.0246; Pollack–Stevens;
Balakrishnan–Tuitman 1710.01673; Best 1806.03393; CMTV 2106.09315; Besser–de Jeu;
Borwein JA 1985; Farach-Colton–Tsai 1504.05240; Luschny FastFactorial; KL 2404.05177;
Neiger et al. 2601.17422; BGP ICALP16/JACM17; Pouly–Graça 1409.0451; vdH D-algebraic +
numeric composition; BGS SIAM 2007; Bostan–Neiger–Yurkevich 2302.04299; incomplete Wilson
1310.2691; Hesse–Allender–Barrington JCSS 2002; plus full URL lists in the three crawler
reports (session transcripts).
