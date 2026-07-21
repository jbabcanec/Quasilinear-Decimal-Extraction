# arXiv submission — metadata & checklist

## What to upload
Upload **`decimal_extraction_paper_arxiv.tar.gz`** (it contains only `decimal_extraction_paper.tex`).
The paper is fully self-contained — manual `thebibliography` (no BibTeX, no `.bbl`/`.bib`), no figures, no `\input`. Packages used, all on arXiv's TeXLive: `geometry, amsmath, amssymb, amsthm, booktabs, url, hyperref`.
**Verified:** compiles clean from a fresh directory with `pdflatex` ×2 → 17 pp, 0 undefined refs/cites, no overfull-box problems.

## Categories
- **Primary:** `math.NT` (Number Theory)
- **Cross-list:** `cs.CC` (Computational Complexity) — for the barrier / #P-neighborhood analysis
- Alternative cross-list, if you'd rather target an algorithms audience: `cs.SC` (Symbolic Computation) or `cs.DS`.

## Classifications (optional fields)
- **MSC 2020:** 11Y16 (number-theoretic algorithms; complexity) [primary]; 11A63 (radix representation, digital problems); 68Q25 (analysis of algorithms and problem complexity); 11K16 (radix expansions).
- **ACM class:** F.2.1 (Numerical Algorithms and Problems).

## Title
Extracting decimal digits of π in quasi-linear time from binary BBP formulas

## Authors
Joseph Babčanec (Benedict College). The arXiv form accepts UTF-8; if it balks, use "Joseph Babcanec".

## Comments field (suggested)
`17 pages. Reference implementation (Python + C++) and exact-arithmetic validators accompany the source. AI-assistance disclosure in the acknowledgments.`

## Abstract (plain text — paste into the arXiv abstract box)
The Bailey-Borwein-Plouffe (BBP) formula computes isolated binary or hexadecimal digits of pi without computing earlier digits, but no analogous decimal formula exists: Borwein, Borwein and Galway proved that pi admits no Machin-type BBP arctangent formula in any base that is not a power of 2, and every published decimal digit-extraction algorithm runs in essentially quadratic time or worse. We show that the N-th decimal digit of pi can nevertheless be computed from the ordinary hexadecimal BBP formula in quasi-linear time -- O(N log^3 N) bit operations -- at the price of Theta(N) bits of auxiliary space, chiefly the read-only binary expansion of the single integer 5^(N-1): after multiplying the BBP series termwise by 10^(N-1), every term takes the form ±2^E 5^(N-1)/mu with mu odd; terms with E>=0 reduce to word-size modular exponentiation, while the Theta(N) terms with E<0 satisfy an exact identity that expresses their fractional parts through residues and bit-windows of the shared integer 5^(N-1), batchable with accumulating remainder trees. The method applies verbatim to every constant with a BBP formula in a base 2^r (e.g. ln 2), yields the digits of pi in any base B in quasi-linear time, and resolves the obstruction described by Zudilin (arXiv:2409.10097), whose base-5 (originally base-10) attempt fails precisely on the E<0 terms. Implementations in Python and C++ are verified at decimal positions through 10^7 against independent computations and public digit corpora. We then reduce every phase but one to O(N log log N) word operations via sieve-batched exponentiation and a Bezout split, consolidating the remainder into a single 2-adic integer V. The one remaining logarithm is analyzed as a barrier with several machine-validated equivalent formulations; a classical factorization technique of Borwein splits off the multiplicative faces (products, discrete-logarithm sums, power sums, log-gamma values), all falling to O(N log^2 N), while the additive core -- an isolated harmonic number modulo 2^N, equivalently a 2-adic digamma value -- remains at the product-tree bound and is isolated as the precise open problem, shown equivalent to single-argument factorial-mod-prime-power and neighbored by the #P-hard / Kloosterman "wild" side of the tame/wild divide. No factoring-hardness argument protects it. We state the conjectural N log N floor precisely, with the sublinear-time and polylogarithmic-space questions that remain open.

## Pre-submission checklist
- [x] Compiles clean from a clean directory (`pdflatex` ×2), 17 pp
- [x] 37 references — all cited, 0 orphans, 0 dangling `\ref`
- [x] Every citation checked against its primary source (audit: `paper_audit.md`)
- [x] Two prior citation errors fixed (BBG04 author order; ZZ15 venue)
- [x] Main-theorem space/time claim honest (remainder-forest Proposition)
- [x] AI-assistance disclosure present in Acknowledgments
- [x] No `\write18`, no shell-escape, no external files

## Note on the endorsement/moderation angle
This is a strong, verifiable, honestly-scoped result (first implemented+verified N^{1+o(1)} decimal extractor; the barrier is stated as open, not overclaimed). If you need an endorsement for `math.NT`, the reproducible code + public-corpus verification of the 10^7 digit make the result easy for an endorser to check in minutes.
