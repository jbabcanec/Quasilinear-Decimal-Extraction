# Contextualization review of the decimal-extraction paper

Reviewed: `decimal_extraction_paper.tex` at commit a4a5c6f (July 2026 draft, 17 pp.).
Literature check performed September 2026.

This is a referee-style pass on one question only: is the work positioned
correctly relative to what exists? It does not touch the mathematics.
Section 1 is the verdict, Section 2 lists the claims a knowledgeable referee will
push back on with exact replacement text, Section 3 lists the missing
citations with ready-to-paste `\bibitem` entries, Section 4 is the
literature check itself.

## 1. Verdict

The paper is already far better contextualized than most arXiv submissions in
this area: it cites the whole decimal ladder (Plouffe 1996, Bellard 1997,
Gourdon 2003, Plouffe 2022, Zudilin 2024), states the BBG04 obstruction
correctly as "no *Machin-type* formula", has an honest "what is, and is not,
claimed" paragraph, and benchmarks against the AGM rather than against the
weaker extractors. Nothing in the framing is wrong in substance.

Three things will still draw a referee's pen:

1. **One word is an overclaim.** The abstract (and the README) say "first
   implemented and verified *subquadratic* decimal extractor." Gourdon's
   2003 note describes its own Theorem 1 as achieving "sub-quadratic time"
   (his bound is N² log log N / log² N, which is o(N²)), and it was implemented
   as `pidec`. The paper's own intro item (ii) makes the correct claim,
   "N^{1+o(1)} (softly linear)"; the abstract and README should say the same.
2. **The comparison table has no "compute all digits" row.** The text
   correctly names the AGM as the real benchmark, but the table in
   Section 3 only lists extractors, so a skimming reader sees "this paper:
   O(N log³N)" beating everything. Add the AGM / Chudnovsky row.
3. **Three factual statements in the framing carry no citation** and are the
   kind a referee asks for: "binary splitting of the Chudnovsky series",
   "distributed binary-BBP projects", and "a practical need noted by record
   holders". All three are true and easy to source (Section 3).

## 2. Claims to adjust, with replacement text

### 2.1 "Subquadratic" (abstract, l. 51; README)

Current:

> to our knowledge this is the first implemented and verified subquadratic
> decimal extractor.

Replace with:

> to our knowledge this is the first implemented and verified decimal
> extractor of quasi-linear complexity $N^{1+o(1)}$; the strongest previous
> method is quadratic up to polylogarithmic factors.

Reason: Gourdon (2003), first page: "sub-quadratic time"; the bound
N² log log N / log² N is formally o(N²). Saying "subquadratic" invites a
one-line rejection of the novelty claim by anyone who has read his note.
The README has been updated accordingly on this branch; the abstract still
needs the edit (requires a LaTeX rebuild, not available in this environment).

### 2.2 Comparison table (Section 3, "Comparison")

Add one row above "This paper":

```latex
Full computation (AGM~\cite{Brent76,Salamin76}; Chudnovsky binary
splitting~\cite{ChCh89,HP98}) & $O(N\log^{2}N)$ & $\Theta(N)$ bits & any
& implemented (records)\\
\midrule
```

and add to the footnote: "The full-computation row produces every digit and
is the benchmark named in the introduction; the extractors above it produce
one digit."

### 2.3 "N^{1.8} for the strongest published method" (abstract, l. 50)

The abstract presents N^{1.8} as a property of Gourdon's method. Section 7
makes clear it is *your* fit of *your* re-implementation. Referees dislike
an abstract that reads as a measured property of someone else's code.
Replace:

> at a measured exponent of $N^{1.077}$ against $N^{1.8}$ for the strongest
> published method

with:

> at a measured exponent of $N^{1.077}$ over $N\in[10^{4},10^{7}]$, against
> $N^{1.8}$ for a straightforward implementation of Gourdon's method on the
> same machine

### 2.4 Uncited framing statements (Introduction)

| Statement (line) | Cite |
|---|---|
| "softly linear time by binary splitting of the Chudnovsky series" (l. 133) | `ChCh89`, `HP98` |
| "in the style of distributed binary-BBP projects" (l. 147) | `Percival00`, `BBMW13` |
| "a practical need noted by record holders" (Open problems, item 3) | `Yee`, `BBMW13` |
| "O(N log N)-type time" for BBP hex digits (l. 78) | `Bailey06` |

For item 3 of the open problems, the concrete fact to state is: every
decimal record since 2010 (y-cruncher, through the 314-trillion-digit run of
November 2025) is verified by converting the *end* of the decimal expansion
to hexadecimal and spot-checking with a binary BBP formula (Bellard's
variant). No record has ever been spot-checked at an isolated *decimal*
position, because no practical primitive existed. That is the one-sentence
motivation for the algorithm and it is currently missing from the paper.

### 2.5 The Zudilin paragraph (Introduction)

Accurate as written. Two small points. (a) The note is 3 pages and its
abstract says only that it "discusses how to compute (promptly)" base-5
digits; quoting "obvious flaw" is fair since that phrase is his. (b) As of
September 2026 there is still no v3 and no follow-up paper by Zudilin or
anyone else on the base-5 series (Section 4), so "the note has no published
successor" stands. Consider adding "as of [month year]" so the sentence does
not silently go stale.

### 2.6 The Gourdon-tradeoff paragraph

The paragraph about Gourdon's asserted Theorem 2 (memory m, time
N² log³N log log N / (m log²(N/m))) is the most carefully written priority
statement in the paper and should stay exactly as it is. One addition would
strengthen it: note that at m = N his formula gives O(N log³N log log N /
log²(1)), which is *undefined* (log(N/m) = 0), so the "endpoint would
likewise be softly linear" reading is an extrapolation, not something the
formula states. That is in your favour and worth a footnote.

### 2.7 Scope of the BBG04 obstruction

The abstract's "no decimal analogue exists" is immediately qualified by
"Machin-type", which is correct. A referee from the BBP community will
nevertheless want one sentence acknowledging that non-binary BBP formulas
*do* exist for other constants (Broadhurst's base-3 formulas for π² and
ζ(3); Bailey's compendium lists them), so that "π has none in base 10" is
seen as specific to π and to degree-1 formulas, not a general phenomenon.
Suggested insertion after the BBG04 sentence in the introduction:

> (Non-binary BBP formulas are not impossible in general: Broadhurst found
> base-3 formulas for $\pi^{2}$ and $\zeta(3)$~\cite{Broadhurst98}, and
> they have been used for large computations~\cite{BBMW13}; the obstruction
> is specific to $\pi$ and to Machin-type formulas.)

## 3. Missing references (ready to paste)

```latex
\bibitem{ChCh89}
D.~V.~Chudnovsky, G.~V.~Chudnovsky,
\emph{The computation of classical constants},
Proc.\ Nat.\ Acad.\ Sci.\ USA \textbf{86} (1989), 8178--8182.

\bibitem{HP98}
B.~Haible, T.~Papanikolaou,
\emph{Fast multiprecision evaluation of series of rational numbers},
in: ANTS-III, LNCS \textbf{1423}, Springer, 1998, 338--350.

\bibitem{Percival00}
C.~Percival,
\emph{PiHex: A distributed effort to calculate $\pi$}, 1998--2000,
\url{https://www.cecm.sfu.ca/projects/pihex/}.
(Quadrillionth binary digit of $\pi$, September 2000.)

\bibitem{BBMW13}
D.~H.~Bailey, J.~M.~Borwein, A.~Mattingly, G.~Wightwick,
\emph{The computation of previously inaccessible digits of $\pi^{2}$ and
Catalan's constant},
Notices Amer.\ Math.\ Soc.\ \textbf{60} (2013), 844--854.

\bibitem{Broadhurst98}
D.~J.~Broadhurst,
\emph{Polylogarithmic ladders, hypergeometric series and the ten millionth
digits of $\zeta(3)$ and $\zeta(5)$}, arXiv:math/9803067 (1998).

\bibitem{Bailey06}
D.~H.~Bailey,
\emph{The BBP algorithm for pi}, manuscript, 2006,
\url{https://www.davidhbailey.com/dhbpapers/bbp-alg.pdf}.

\bibitem{Yee}
A.~J.~Yee,
\emph{y-cruncher: a multi-threaded pi program}, and the record log
\url{https://www.numberworld.org/y-cruncher/records.html}
(decimal records verified by binary BBP spot checks).

\bibitem{MR07}
R.~J.~McIntosh, E.~L.~Roettger,
\emph{A search for Fibonacci--Wieferich and Wolstenholme primes},
Math.\ Comp.\ \textbf{76} (2007), 2087--2094.
```

`MR07` belongs in Section 6 next to `CGH14`: Wolstenholme-prime searches
are the other large computation of harmonic numbers modulo prime powers,
and they too go through products, not through any direct fast evaluation of
the sum. It supports the "no known fast algorithm for the additive part"
statement with a second independent line of evidence.

Optional but harmless: Adamchik--Wagon, *A simple formula for $\pi$*, Amer.
Math.\ Monthly 104 (1997) 852--855, for the general BBP-formula framework.

## 4. Literature check (September 2026)

Searched for: successors to arXiv:2409.10097; any 2025--2026 preprint on
decimal (or non-binary) digit extraction of π; any claim of quasi-linear or
subquadratic decimal extraction; the current record and how it was verified.

Findings:

- **Zudilin, arXiv:2409.10097.** Still at v2 (17 Sep 2024), 3 pages, title
  "in base 5" (v1: "in base 10"). No later version, no follow-up by him or
  by others found. The paper's statement stands.
- **No competing decimal extractor.** Nothing after Plouffe 2022 was found
  that claims decimal digit extraction below quadratic. The only 2026 item
  returned is a practitioner's implementation guide (Vestergaard, "Practical
  implementation of π algorithms"), which surveys known methods.
- **Records.** Current decimal record: 314 trillion digits (November 2025,
  y-cruncher, single 110-day run). Verification: hexadecimal BBP spot check
  of the trailing digits, performed multiple times. This is the practice the
  open problem (3) should cite explicitly.
- **Gourdon 2003.** The archived PDF describes Theorem 1 as sub-quadratic
  and reports an implementation, `pidec`, roughly 5--10 times faster than
  Bellard's at 5·10³--10⁶. This is the source for the "subquadratic"
  correction in 2.1.

Caveat: web search coverage of arXiv listings is incomplete; a final
pre-submission check of math.NT and cs.DS/cs.SC listings for
"digit extraction" from 2025 onward is recommended.

## 5. What was changed on this branch

- `README.md`: "subquadratic" claim replaced by "quasi-linear (N^{1+o(1)})",
  with the reason stated in 2.1.
- This file.

Not changed (needs a LaTeX rebuild, unavailable here): the abstract wording
(2.1, 2.3), the table row (2.2), the citations (2.4, 3), the insertions
(2.5--2.7). Each is a copy-paste from this document.
